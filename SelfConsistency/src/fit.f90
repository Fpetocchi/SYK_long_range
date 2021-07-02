module fit

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !


   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface reconstruct_lastCoeffs
      module procedure :: reconstruct_lastCoeffs_F
      module procedure :: reconstruct_lastCoeffs_B
   end interface reconstruct_lastCoeffs

   !---------------------------------------------------------------------------!
   !PURPOSE: Module custom types
   !---------------------------------------------------------------------------!
   type AndersonParam
      integer                               :: Norb
      real(8),allocatable                   :: Epsk(:,:,:)
      real(8),allocatable                   :: Vk(:,:,:)
      real(8),allocatable                   :: Eloc(:,:)
      logical                               :: status=.false.
   end type AndersonParam

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif
   !Global parameters
   !Fit variables
   integer,private                          :: Nfit
   integer,private                          :: MaxMom
   integer,private                          :: Nfreq
   !
   complex(8),private                       :: ContinuityConstraint(2)
   !
   integer,parameter,private                :: cg_niter=400
   real(8),parameter,private                :: cg_Ftol=1e-8
   real(8),parameter,private                :: hwband=10d0
   real(8),parameter,private                :: noisefact=0.1
   !Anderson Parameters variables
   type(AndersonParam),private              :: AndPram
   !Moments
   real(8),allocatable,private              :: Moments(:,:,:)
   !Shared
   real(8),allocatable,private              :: wmats(:)
   complex(8),allocatable,private           :: Component(:)
   real(8),private                          :: df_eps=tiny(1d0)

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   public :: fit_moments
   public :: fit_delta
   public :: G_Moments
   public :: S_Moments
   public :: W_Moments

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Rebuild the last Even and Odd coefficients so as to impse continuity
   !betwee the fitted tail and the original function.
   !Here I'm assuming that all the odd coefficients are negative
   !even if M3 should be positive. However the fit is unable to simultaneously
   !find M1>M3 and preserve continuitiy. So I'm keeping negative M3.
   !---------------------------------------------------------------------------!
   subroutine reconstruct_lastCoeffs_F(Coeff,ConstrF)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(inout),allocatable     :: Coeff(:)
      complex(8),intent(in)                 :: ConstrF(2)
      integer                               :: imoment,lastMom
      integer                               :: EvenExp,OddExp
      logical                               :: lastisEven
      real(8)                               :: EvenCoef,OddCoef
      !
      if(verbose)write(*,"(A)") "---- reconstruct_lastCoeffs_F"
      lastMom = size(Coeff)-1
      !
      if(mod(lastMom,2).eq.0)then
         EvenExp = lastMom
         OddExp = lastMom - 1
         lastisEven=.true.
      else
         OddExp = lastMom
         EvenExp = lastMom - 1
         lastisEven=.false.
      endif
      !
      EvenCoef = dreal(ConstrF(1))*(ConstrF(2)**EvenExp)
      do imoment=0,EvenExp-2,2
         EvenCoef = EvenCoef - Coeff(imoment)*(ConstrF(2)**(EvenExp-imoment))
      enddo
      !
      OddCoef = -dimag(ConstrF(1))*(ConstrF(2)**OddExp)
      do imoment=1,OddExp-2,2
         OddCoef = OddCoef - abs(Coeff(imoment))*(ConstrF(2)**(OddExp-imoment))
      enddo
      !
      if(lastisEven)then
         Coeff(lastMom) = EvenCoef
         Coeff(lastMom-1) = abs(OddCoef)
      else
         Coeff(lastMom) = abs(OddCoef)
         Coeff(lastMom-1) = EvenCoef
      endif
      !
   end subroutine reconstruct_lastCoeffs_F
   !
   subroutine reconstruct_lastCoeffs_B(Coeff,ConstrF)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(inout),allocatable     :: Coeff(:)
      real(8),intent(in)                    :: ConstrF(2)
      integer                               :: imoment,lastMom
      integer                               :: EvenExp
      real(8)                               :: EvenCoef
      !
      if(verbose)write(*,"(A)") "---- reconstruct_lastCoeffs_B"
      !
      lastMom = size(Coeff)
      EvenExp = 2*lastMom
      !
      EvenCoef = -ConstrF(1)*(ConstrF(2)**EvenExp)
      do imoment=0,lastMom-1
         EvenCoef = EvenCoef + Coeff(imoment)*(ConstrF(2)**(EvenExp-2*imoment))
      enddo
      !
      Coeff(lastMom) = EvenCoef
      !
   end subroutine reconstruct_lastCoeffs_B


   !=============================== FIT DELTA =================================!


   !---------------------------------------------------------------------------!
   !PURPOSE: Setup/initilize the Anderson parameters.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine setupAndPrams(Norb,dirpath,paramFile)
      !
      use parameters
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: Norb
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      !
      character(len=255)                    :: path
      integer                               :: unit,ierr
      integer                               :: Nh,Nfit_read
      integer                               :: ibath,iorb
      real(8)                               :: de
      real(8),allocatable                   :: ReadLine(:)
      logical                               :: filexists
      real(8)                               :: rnd
      !
      !
      if(verbose)write(*,"(A)") "---- setupAndPrams"
      if(AndPram%status) write(*,"(A)") "Warning: Anderson parameters already initilized."
      AndPram%Norb=Norb
      !
      !
      if(.not.allocated(AndPram%Eloc))allocate(AndPram%Eloc(Norb,Nspin))
      if(.not.allocated(AndPram%Epsk))allocate(AndPram%Epsk(Norb,Nfit,Nspin))
      if(.not.allocated(AndPram%Vk))allocate(AndPram%Vk(Norb,Nfit,Nspin))
      !
      path = reg(dirpath)//reg(paramFile)
      call inquireFile(reg(path),filexists,hardstop=.false.,verb=verbose)
      !
      if(filexists)then
         if(verbose)write(*,"(A)") "     Checking the number of Anderson Parameters in "//reg(path)
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind",iostat=ierr)
         read(unit,*) Nfit_read
         close(unit)
         !
         if(Nfit_read.ne.Nfit)then
            if(verbose)write(*,"(A)") "     Number of Parameters is different. Reinitializing."
            filexists = .false.
         elseif(ierr.ne.0)then
            if(verbose)write(*,"(A)") "     Error in opening file. Reinitializing."
            filexists = .false.
         endif
         !
      endif
      !
      if(filexists)then
         !
         if(verbose)write(*,"(A)") "     Reading Anderson Parameters from "//reg(paramFile)
         allocate(ReadLine(4))
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind")
         read(unit,*) Nfit_read
         do iorb=1,Norb
            read(unit,*) AndPram%Eloc(iorb,:)
            do ibath=1,Nfit
               ReadLine=0d0
               read(unit,*) ReadLine
               AndPram%Epsk(iorb,ibath,1) = ReadLine(1) !Epsk_up
               AndPram%Epsk(iorb,ibath,2) = ReadLine(3) !Epsk_dw
               AndPram%Vk(iorb,ibath,1)   = ReadLine(2) !Vk_up
               AndPram%Vk(iorb,ibath,2)   = ReadLine(4) !Vk_dw
            enddo
            read(unit,*)
         enddo
         deallocate(ReadLine)
         close(unit)
         !
      elseif((.not.filexists).and.(.not.AndPram%status))then
         !
         if(verbose)write(*,"(A)") "     Initializing Anderson Parameters."
         !local energy
         call random_number(rnd)
         AndPram%Eloc=(rnd-0.5)*noisefact
         !bath energies
         Nh=Nfit/2
         if(mod(Nfit,2)==0)then
            de=hwband/max(Nh-1,1)
            call random_number(rnd)
            AndPram%Epsk(:,Nh,:)  = -1.d-3 + (rnd-0.5)*noisefact
            call random_number(rnd)
            AndPram%Epsk(:,Nh+1,:)=  1.d-3 + (rnd-0.5)*noisefact
            do ibath=1,Nh-1
               call random_number(rnd)
               AndPram%Epsk(:,ibath,:)        = -hwband + (ibath-1)*de + (rnd-0.5)*noisefact
               call random_number(rnd)
               AndPram%Epsk(:,Nfit-ibath+1,:) = +hwband - (ibath-1)*de + (rnd-0.5)*noisefact
            enddo
         elseif(mod(Nfit,2)/=0)then
            de=hwband/Nh
            call random_number(rnd)
            AndPram%Epsk(:,Nh+1,:)= 1.d-3 + (rnd-0.5)*noisefact
            do ibath=1,Nh
               call random_number(rnd)
               AndPram%Epsk(:,ibath,:)        = -hwband + (ibath-1)*de + (rnd-0.5)*noisefact
               call random_number(rnd)
               AndPram%Epsk(:,Nfit-ibath+1,:) = +hwband - (ibath-1)*de + (rnd-0.5)*noisefact
            enddo
         endif
         !bath hybridizations
         do ibath=1,Nfit
            call random_number(rnd)
            AndPram%Vk(:,ibath,:) = max(0.1d0,1.d0/sqrt(dble(Nfit))) + (rnd-0.5)*noisefact
         enddo
         !
      endif
      !
      AndPram%status = .true.
      !
   end subroutine setupAndPrams


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file the Delta Anderson parameters.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine dump_AndPrams(dirpath,paramFile)
      !
      use parameters
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      !
      character(len=255)                    :: path
      integer                               :: unit
      integer                               :: ibath,iorb
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- dump_AndPrams"
      if(.not.AndPram%status) stop "dump_AndPrams: Anderson parameters not initilized."
      !
      !
      call inquireDir(reg(dirpath),filexists,verb=verbose)
      path = reg(dirpath)//reg(paramFile)
      if(verbose)write(*,"(A)") "     Dump "//reg(path)//" (readable)"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",action="write",position="rewind")
      write(unit,"(I5,A)") Nfit," Number of bath levels."
      do iorb=1,AndPram%Norb
         write(unit,"(2E20.12)") AndPram%Eloc(iorb,1),AndPram%Eloc(iorb,2)
         do ibath=1,Nfit
            write(unit,"(4E20.12)") AndPram%Epsk(iorb,ibath,1),AndPram%Vk(iorb,ibath,1),AndPram%Epsk(iorb,ibath,2),AndPram%Vk(iorb,ibath,2)
         enddo
         write(unit,*)
      enddo
      close(unit)
      !
   end subroutine dump_AndPrams


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the Anderson Hybridization function of a given orbital
   !         (iorb) from the parameters list [Epsk,Vk,Eo].
   !         The frequency mesh is not global as the routines are made available
   !         in the main.
   !TEST ON:
   !---------------------------------------------------------------------------!
   function DeltaAnderson(AndParaVec,wm) result(Delta)
      use parameters
      implicit none
      real(8),dimension(:)                  :: AndParaVec
      real(8),dimension(:)                  :: wm
      complex(8),dimension(size(wm))        :: Delta
      integer                               :: ibath,io,iw
      integer                               :: Nbath
      real(8),allocatable                   :: Ek(:),Vk(:)
      !
      !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
      !
      Delta=czero
      Nbath=int(size(AndParaVec)/2)
      allocate(Ek(Nbath));Ek=0d0
      allocate(Vk(Nbath));Vk=0d0
      !
      do ibath=1,Nbath
         io = ibath
         Ek(ibath) = AndParaVec(io)
      enddo
      do ibath=1,Nbath
         io = Nfit + ibath
         Vk(ibath) = AndParaVec(io)
      enddo
      !
      do iw=1,size(wm)
         Delta(iw) = sum( Vk(:)*Vk(:)/( dcmplx(0d0,wm(iw)) - Ek(:)) )
      enddo
      !
   end function DeltaAnderson
   !
   function EoDeltaAnderson(AndParaVec,wm) result(Delta)
      use parameters
      implicit none
      real(8),dimension(:)                  :: AndParaVec
      real(8),dimension(:)                  :: wm
      complex(8),dimension(size(wm))        :: Delta
      integer                               :: ibath,io,iw
      integer                               :: Nbath
      real(8),allocatable                   :: Ek(:),Vk(:)
      !
      !\Delta_{aa} = Eloc_{aa} + \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
      !
      Delta=czero
      Nbath=int((size(AndParaVec)-1)/2)
      allocate(Ek(Nbath));Ek=0d0
      allocate(Vk(Nbath));Vk=0d0
      !
      do ibath=1,Nbath
         io = ibath
         Ek(ibath) = AndParaVec(io)
      enddo
      do ibath=1,Nbath
         io = Nfit + ibath
         Vk(ibath) = AndParaVec(io)
      enddo
      !
      do iw=1,size(wm)
         Delta(iw) = AndParaVec(size(AndParaVec)) + sum( Vk(:)*Vk(:)/( dcmplx(0d0,wm(iw)) - Ek(:)) )
      enddo
      !
   end function EoDeltaAnderson


   !---------------------------------------------------------------------------!
   !PURPOSE: Distance between input/output functions minimized by the module.
   !         Here Nfreq is global and defined in the callable fit subroutine.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine chi2_Delta(Npara,AndParaVec,chi2)
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: AndParaVec
      real(8)                               :: chi2
      complex(8),dimension(Nfreq)           :: Delta
      !
      Delta = DeltaAnderson(AndParaVec,wmats)
      !
      chi2=sum(abs(Component(:)-Delta(:))**2)
      !
   end subroutine chi2_Delta
   !
   subroutine chi2_ShiftedDelta(Npara,AndParaVec,chi2)
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: AndParaVec
      real(8)                               :: chi2
      complex(8),dimension(Nfreq)           :: ShiftedDelta
      !
      ShiftedDelta = EoDeltaAnderson(AndParaVec,wmats)
      !
      chi2=sum(abs(Component(:)-ShiftedDelta(:))**2)
      !
   end subroutine chi2_ShiftedDelta


   !---------------------------------------------------------------------------!
   !PURPOSE: Fit a given function [Norn,Nfreq,Nspin] assuming an hybridization
   !         like functional form with an addtional shift.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine fit_Delta(funct,Beta,Nb,dirpath,paramFile,FitMode,Eloc,filename,Wlimit,coef01)
      !
      use parameters
      use utils_misc
      use file_io
      implicit none
      !
      complex(8),intent(in)                 :: funct(:,:,:)
      real(8),intent(in)                    :: Beta
      integer,intent(in)                    :: Nb
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      character(len=*),intent(in)           :: FitMode
      real(8),intent(inout)                 :: Eloc(:,:)
      character(len=*),intent(in),optional  :: filename
      integer,intent(in),optional           :: Wlimit
      real(8),intent(inout),optional        :: coef01(:,:)
      !
      integer                               :: Norb,Niter,Wstart
      integer                               :: iorb,ispin,ifit
      real(8)                               :: chi
      real(8),allocatable                   :: ParamVec(:)
      complex(8),allocatable                :: funct_print(:)
      !
      !
      if(verbose)write(*,"(A)") "---- fit_Delta"
      !
      !
      Nfit = Nb
      Norb = size(funct,dim=1)
      Wstart = 1
      if(present(Wlimit))Wstart = Wlimit
      Nfreq = size(funct,dim=2) - Wstart + 1
      !
      call assert_shape(Eloc,[Norb,Nspin],"fit_Delta","Eloc")
      if(present(coef01))call assert_shape(coef01,[Norb,Nspin],"fit_Delta","coef01")
      allocate(Component(Nfreq));Component=czero
      !
      call setupAndPrams(Norb,dirpath,paramFile)
      call dump_AndPrams(dirpath,"used."//paramFile)
      !
      select case(reg(FitMode))
         case default
            !
            stop "Available modes for moment fitting: Standard, Shifted."
            !
         case("Standard")
            !
            write(*,"(A)")"     Fitting Delta(iw)."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Nfreq) + (Wstart-1)*2d0*pi/Beta
            allocate(ParamVec(2*Nfit));ParamVec=0d0
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,Wstart:Wstart+Nfreq-1,ispin)
                  ParamVec(1:Nfit) = AndPram%Epsk(iorb,:,ispin)
                  ParamVec(1+Nfit:2*Nfit) = AndPram%Vk(iorb,:,ispin)
                  !
                  call fit_wrapper(chi2_Delta,ParamVec,chi,Niter)
                  !
                  write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                  write(*,"(A,I,2(A,F))")"     Iterations: ",Niter," Chi^2: ",chi,", Wstart: ",wmats(1)
                  !
                  AndPram%Epsk(iorb,:,ispin) = ParamVec(1:Nfit)
                  AndPram%Vk(iorb,:,ispin) = ParamVec(1+Nfit:2*Nfit)
                  !
               enddo
            enddo
            !
         case("Shifted")
            !
            write(*,"(A)")"     Fitting Eloc+Delta(iw)."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Nfreq) + (Wstart-1)*2d0*pi/Beta
            allocate(ParamVec(2*Nfit+1));ParamVec=0d0
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,Wstart:Wstart+Nfreq-1,ispin)
                  ParamVec(1:Nfit) = AndPram%Epsk(iorb,:,ispin)
                  ParamVec(1+Nfit:2*Nfit) = AndPram%Vk(iorb,:,ispin)
                  ParamVec(2*Nfit+1) = AndPram%Eloc(iorb,ispin)
                  !
                  call fit_wrapper(chi2_ShiftedDelta,ParamVec,chi,Niter)
                  !
                  write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                  write(*,"(A,I,A,F)")"     Iterations: ",Niter," Chi^2: ",chi
                  !
                  AndPram%Eloc(iorb,ispin) = ParamVec(2*Nfit+1)
                  AndPram%Epsk(iorb,:,ispin) = ParamVec(1:Nfit)
                  AndPram%Vk(iorb,:,ispin) = ParamVec(1+Nfit:2*Nfit)
                  !
               enddo
            enddo
            !
         end select
         !
         deallocate(Component,wmats)
         !
         Eloc = AndPram%Eloc
         if(present(coef01))then
            coef01=0d0
            do ispin=1,Nspin
               do iorb=1,Norb
                  do ifit=1,Nfit
                     coef01(iorb,ispin) = coef01(iorb,ispin) + AndPram%Vk(iorb,ifit,ispin)**2
                  enddo
               enddo
            enddo
         endif
         !
         call dump_AndPrams(dirpath,paramFile)
         !
         if(present(filename))then
            !
            allocate(wmats(size(funct,dim=2)));wmats=0d0
            wmats=FermionicFreqMesh(Beta,size(funct,dim=2))
            !
            allocate(funct_print(size(funct,dim=2)))
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  funct_print=czero
                  funct_print = DeltaAnderson(ParamVec(1:2*Nfit),wmats)
                  call dump_FermionicField(funct_print,reg(dirpath)//"fits/",reg(filename)//"_o"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
                  !
               enddo
            enddo
            !
            deallocate(funct_print,wmats)
            !
         endif
         !
         deallocate(ParamVec)
         if(AndPram%status.eq..true.) then
           if(allocated(AndPram%Eloc))deallocate(AndPram%Eloc)
           if(allocated(AndPram%Epsk))deallocate(AndPram%Epsk)
           if(allocated(AndPram%Vk))deallocate(AndPram%Vk)
           AndPram%status=.false.
         endif
         !
      end subroutine fit_Delta



   !============================== FIT MOMENTS ================================!



   !---------------------------------------------------------------------------!
   !PURPOSE: Setup/initilize the Green's function moments.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine setupMoments(Norb,dirpath,paramFile,refresh)
      !
      use parameters
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: Norb
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      logical,intent(in)                    :: refresh
      !
      character(len=255)                    :: path
      integer                               :: unit,ierr
      integer                               :: Nfit_read
      integer                               :: imoment,ispin,iorb
      real(8),allocatable                   :: ReadLine(:)
      logical                               :: filexists
      real(8)                               :: rnd
      !
      !
      if(verbose)write(*,"(A)") "---- setupMoments"
      !
      !
      if(allocated(Moments)) deallocate(Moments)
      allocate(Moments(Norb,Nspin,0:Nfit))
      !
      path = reg(dirpath)//reg(paramFile)
      call inquireFile(reg(path),filexists,hardstop=.false.,verb=verbose)
      !
      if(filexists)then
         if(verbose)write(*,"(A)") "     Checking the number of coefficients in "//reg(path)
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind",iostat=ierr)
         read(unit,"(1I5)") Nfit_read
         close(unit)
         !
         if(Nfit_read.ne.Nfit)then
            if(verbose)write(*,"(A)") "     Number of coefficients is different. Reinitializing."
            filexists = .false.
         elseif(ierr.ne.0)then
            if(verbose)write(*,"(A)") "     Error in opening file. Reinitializing."
            filexists = .false.
         endif
         !
      endif
      !
      if(filexists.and.(.not.refresh))then
         !
         if(verbose)write(*,"(A)") "     Reading Moments from "//reg(paramFile)
         allocate(ReadLine(Nspin*Norb))
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind")
         read(unit,*) Nfit_read
         do imoment=0,MaxMom
            ReadLine=0d0
            read(unit,*) ReadLine
            Moments(:,1,imoment) = ReadLine(1:Norb)
            Moments(:,2,imoment) = ReadLine(1+Norb:2*Norb)
         enddo
         deallocate(ReadLine)
         close(unit)
         !
      else
         !
         if(verbose)write(*,"(A)") "     Initializing Moments."
         Moments = 0d0
         do imoment=0,MaxMom
            do iorb=1,Norb
               do ispin=1,Nspin
                  call random_number(rnd)
                  Moments(iorb,ispin,imoment) = 1d0 - (rnd*noisefact)
                  if((imoment.gt.0).and.(mod(imoment,2).eq.0)) Moments(iorb,ispin,imoment) = -abs(Moments(iorb,ispin,imoment))
               enddo
            enddo
         enddo
         !
      endif
      !
   end subroutine setupMoments


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file the Moment list.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine dump_Moments(dirpath,paramFile)
      !
      use parameters
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      !
      character(len=255)                    :: path
      integer                               :: unit
      integer                               :: imoment,iorb,Norb
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- dump_Moments"
      if(.not.allocated(Moments)) stop "dump_Moments: Moments not allocated."
      Norb=size(Moments,dim=1)
      !
      !
      call inquireDir(reg(dirpath),filexists,verb=verbose)
      path = reg(dirpath)//reg(paramFile)
      if(verbose)write(*,"(A)") "     Dump "//reg(path)//" (readable)"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",action="write",position="rewind")
      write(unit,"(1I5,A)") Nfit," Number of coefficients."
      do imoment=0,MaxMom
         write(unit,"(999E20.12)") (Moments(iorb,1,imoment),iorb=1,Norb),(Moments(iorb,2,imoment),iorb=1,Norb)
      enddo
      close(unit)
      !
   end subroutine dump_Moments


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates the generic functional form via a moment representation.
   !         Global module variables are not used as the routines are available
   !         in the main.
   !TEST ON:
   !---------------------------------------------------------------------------!
   function G_Moments(ParamVec,wm,Gconstr) result(Gf)
      use parameters
      implicit none
      real(8),dimension(:)                  :: ParamVec
      real(8),dimension(:)                  :: wm
      complex(8),optional                   :: Gconstr(2)
      complex(8),dimension(size(wm))        :: Gf
      real(8),dimension(size(wm))           :: ReGf
      real(8),dimension(size(wm))           :: ImGf
      real(8),allocatable                   :: Coeff(:)
      integer                               :: iw,exp
      integer                               :: imoment,lastMom
      !
      if(verbose)write(*,"(A)") "---- G_Moments"
      !
      if(present(Gconstr))then
         !
         !ParamVec contains from exp 2 and without the last two
         lastMom = size(ParamVec) + 2 + 2 - 1      ! +(fixed M0&M1) + (2 reconstucted) - (starting from 0)
         allocate(Coeff(0:lastMom));Coeff=0d0
         !
         Coeff(0) = 0d0
         Coeff(1) = -1d0
         Coeff(2:lastMom-2) = ParamVec
         call reconstruct_lastCoeffs(Coeff,Gconstr)
         !
      else
         !
         !ParamVec contains also the last two recontructed elements
         lastMom = size(ParamVec) + 2 - 1      ! +(fixed M0&M1) - (starting from 0)
         allocate(Coeff(0:lastMom));Coeff=0d0
         !
         Coeff(0) = 0d0
         Coeff(1) = -1d0
         Coeff(2:lastMom) = ParamVec
         !
      endif
      !
      Gf=czero
      ReGf=0d0
      ImGf=0d0
      do iw=1,size(wm)
         !
         do imoment=0,lastMom
            !
            !the index is equal to the exponent
            exp = imoment
            !
            !Here I'm assuming that all the odd coefficients are negative
            !even if M3 should be positive. However the fit is unable to simultaneously
            !find M1>M3 and preserve continuitiy. So I'm keeping negative M3.
            if(mod(exp,2).eq.0)then
               ReGf(iw) = ReGf(iw) + Coeff(imoment)/(wm(iw)**exp)
            else
               ImGf(iw) = ImGf(iw) - abs(Coeff(imoment))/(wm(iw)**exp)
            endif
            !
         enddo
         !
         Gf(iw) = dcmplx( ReGf(iw) , ImGf(iw) )
         !
      enddo
      !
   end function G_Moments
   !
   function S_Moments(ParamVec,wm,Sconstr) result(Sigma)
      use parameters
      implicit none
      real(8),dimension(:)                  :: ParamVec
      real(8),dimension(:)                  :: wm
      complex(8),optional                   :: Sconstr(2)
      complex(8),dimension(size(wm))        :: Sigma
      real(8),dimension(size(wm))           :: ReSigma
      real(8),dimension(size(wm))           :: ImSigma
      real(8),allocatable                   :: Coeff(:)
      integer                               :: iw,exp
      integer                               :: imoment,lastMom
      !
      if(verbose)write(*,"(A)") "---- S_Moments"
      !
      if(present(Sconstr))then
         !
         !ParamVec contains from exp 0 and without the last two
         lastMom = size(ParamVec) + 2 - 1      ! + (2 reconstucted) - (starting from 0)
         allocate(Coeff(0:lastMom));Coeff=0d0
         !
         Coeff(0:lastMom-2) = ParamVec
         call reconstruct_lastCoeffs(Coeff,Sconstr)
         !
      else
         !
         !ParamVec contains also the last two recontructed elements
         lastMom = size(ParamVec) - 1      ! - (starting from 0)
         allocate(Coeff(0:lastMom));Coeff=0d0
         !
         Coeff(0:lastMom) = ParamVec
         !
      endif
      !
      Sigma=czero
      ReSigma=0d0
      ImSigma=0d0
      do iw=1,size(wm)
         !
         do imoment=0,lastMom
            !
            !the index is equal to the exponent
            exp = imoment
            !
            !Here I'm assuming that all the odd coefficients are negative
            !even if M3 should be positive. However the fit is unable to simultaneously
            !find M1>M3 and preserve continuitiy. So I'm keeping negative M3.
            if(mod(exp,2).eq.0)then
               ReSigma(iw) = ReSigma(iw) + Coeff(imoment)/(wm(iw)**exp)
            else
               ImSigma(iw) = ImSigma(iw) - abs(Coeff(imoment))/(wm(iw)**exp)
            endif
            !
         enddo
         !
         Sigma(iw) = dcmplx( ReSigma(iw) , ImSigma(iw) )
         !
      enddo
      !
   end function S_Moments
   !
   function W_Moments(ParamVec,wm,Wconstr) result(W)
      use parameters
      implicit none
      real(8),dimension(:)                  :: ParamVec
      real(8),dimension(:)                  :: wm
      complex(8),optional                   :: Wconstr(2)
      complex(8),dimension(size(wm))        :: W
      real(8),allocatable                   :: Coeff(:)
      integer                               :: iw,exp
      integer                               :: imoment,lastMom
      !
      if(verbose)write(*,"(A)") "---- W_Moments"
      !
      if(present(Wconstr))then
         !
         !ParamVec contains from exp 0 and without the last one
         lastMom = size(ParamVec) + 1 - 1      ! + (1 reconstucted) - (starting from 0)
         allocate(Coeff(0:lastMom));Coeff=0d0
         !
         Coeff(0:lastMom-1) = ParamVec
         call reconstruct_lastCoeffs(Coeff,dreal(Wconstr))
         !
      else
         !
         !ParamVec contains also the last one recontructed elements
         lastMom = size(ParamVec) - 1      ! - (starting from 0)
         allocate(Coeff(0:lastMom));Coeff=0d0
         !
         Coeff(0:lastMom) = ParamVec
         !
      endif
      !
      W=czero
      do iw=1,size(wm)
         !
         do imoment=0,lastMom
            !
            !the index is equal to the even exponent
            exp = 2*imoment
            !
            W(iw) = W(iw) + Coeff(imoment)/(wm(iw)**exp)
            !
         enddo
         !
      enddo
      !
   end function W_Moments


   !---------------------------------------------------------------------------!
   !PURPOSE: Distance between input/output functions minimized by the module.
   !         Here Nfreq is global and defined in  the callable fit subroutine.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine chi2_G_Moments(Npara,MomentVec,chi2)
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: MomentVec
      real(8)                               :: chi2
      complex(8),dimension(Nfreq)           :: Gf
      !
      Gf = G_Moments(MomentVec,wmats,Gconstr=ContinuityConstraint)
      !
      chi2=sum(abs(Component(:)-Gf(:))**2)
      !
   end subroutine chi2_G_Moments
   !
   subroutine chi2_S_Moments(Npara,MomentVec,chi2)
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: MomentVec
      real(8)                               :: chi2
      complex(8),dimension(Nfreq)           :: Sigma
      !
      Sigma = S_Moments(MomentVec,wmats,Sconstr=ContinuityConstraint)
      !
      chi2=sum(abs(Component(:)-Sigma(:))**2)
      !
   end subroutine chi2_S_Moments
   !
   subroutine chi2_W_Moments(Npara,MomentVec,chi2)
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: MomentVec
      real(8)                               :: chi2
      complex(8),dimension(Nfreq)           :: W
      !
      W = W_Moments(MomentVec,wmats,Wconstr=ContinuityConstraint)
      !
      chi2=sum(abs(Component(:)-W(:))**2)
      !
   end subroutine chi2_W_Moments


   !---------------------------------------------------------------------------!
   !PURPOSE: Fit a given function [Norn,Nfreq,Nspin] assuming a generic
   !         functional moment formulation.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine fit_moments(funct,Beta,dirpath,paramFile,FitMode,MomentsOut,filename,Wlimit,verb,refresh)
      !
      use parameters
      use utils_misc
      use file_io
      implicit none
      !
      complex(8),intent(in)                 :: funct(:,:,:)
      real(8),intent(in)                    :: Beta
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      character(len=*),intent(in)           :: FitMode
      real(8),intent(inout),allocatable     :: MomentsOut(:,:,:)
      character(len=*),intent(in),optional  :: filename
      integer,intent(in),optional           :: Wlimit
      logical,intent(in),optional           :: verb
      logical,intent(in),optional           :: refresh
      !
      integer                               :: Norb,Niter,Wstart
      integer                               :: iorb,ispin
      real(8)                               :: chi
      real(8),allocatable                   :: ParamVec(:),MomentsRebuilt(:)
      complex(8),allocatable                :: funct_print(:)
      logical                               :: verb_,refresh_
      !
      !
      if(verbose)write(*,"(A)") "---- fit_moments"
      !
      !
      Nfit = size(MomentsOut,dim=3)
      if(Nfit.lt.5) stop "fit_moments: the fit will not work with less than five moments (min order is 4)."
      if(.not.allocated(MomentsOut)) stop "fit_moments: output moment container not allocated."
      MaxMom = Nfit - 1
      !
      Norb = size(funct,dim=1)
      Wstart = 1
      if(present(Wlimit))Wstart = Wlimit
      Nfreq = size(funct,dim=2) - Wstart + 1
      !
      verb_=.true.
      if(present(verb))verb_=verb
      refresh_=.false.
      if(present(refresh))refresh_=refresh
      !
      allocate(Component(Nfreq));Component=czero
      !
      call setupMoments(Norb,dirpath,paramFile,refresh_)
      allocate(MomentsRebuilt(0:MaxMom));MomentsRebuilt=0d0
      !
      select case(reg(FitMode))
         case default
            !
            stop "Available modes for moment fitting: Green, Sigma, Boson."
            !
         case("Green")
            !
            if(verb_)write(*,"(A)")"     Fitting moments [2:"//str(MaxMom)//"]."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Nfreq) + (Wstart-1)*2d0*pi/Beta
            allocate(ParamVec(Nfit-4));ParamVec=0d0                             !starts from exponent -2 and interanlly reconstruct the last two
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Moments(iorb,ispin,0) = 0d0
                  Moments(iorb,ispin,1) = 1d0
                  Component = funct(iorb,Wstart:Wstart+Nfreq-1,ispin)
                  ParamVec = Moments(iorb,ispin,2:MaxMom-2)                     !starts from exponent -2 and interanlly reconstruct the last two
                  ContinuityConstraint = [Component(1),dcmplx(wmats(1),0d0)]
                  !
                  call fit_wrapper(chi2_G_Moments,ParamVec,chi,Niter)
                  !
                  !This sucks but in order to start from 0 in the subr I need to pass an allocatable which cannot be sliced
                  MomentsRebuilt(0) = 0d0
                  MomentsRebuilt(1) = 1d0
                  MomentsRebuilt(2:MaxMom-2)=ParamVec
                  call reconstruct_lastCoeffs(MomentsRebuilt,ContinuityConstraint)
                  !
                  if(verb_)then
                     write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                     write(*,"(A,I,2(A,F))")"     Iterations: ",Niter," Chi^2: ",chi,", Wstart: ",wmats(1)
                     write(*,"(A,100E16.6)")"     Moments before fit: ",Moments(iorb,ispin,0:MaxMom)
                     write(*,"(A,100E16.6)")"     Moments after fit:  ",MomentsRebuilt(0:MaxMom)
                  endif
                  !
                  Moments(iorb,ispin,0:MaxMom)=MomentsRebuilt(0:MaxMom)
                  !
               enddo
            enddo
            !
         case("Sigma")
            !
            if(verb_)write(*,"(A)")"     Fitting moments [0:"//str(MaxMom)//"]."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Nfreq) + (Wstart-1)*2*pi/Beta
            allocate(ParamVec(Nfit-2));ParamVec=0d0                             !starts from exponent 0 and interanlly reconstruct the last two
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,Wstart:Wstart+Nfreq-1,ispin)
                  ParamVec = Moments(iorb,ispin,0:MaxMom-2)                     !starts from exponent 0 and interanlly reconstruct the last two
                  ContinuityConstraint = [Component(1),dcmplx(wmats(1),0d0)]
                  !
                  call fit_wrapper(chi2_S_Moments,ParamVec,chi,Niter)
                  !
                  !This sucks but in order to start from 0 in the subr I need to pass an allocatable which cannot be sliced
                  MomentsRebuilt(0:MaxMom-2)=ParamVec
                  call reconstruct_lastCoeffs(MomentsRebuilt,ContinuityConstraint)
                  !
                  if(verb_)then
                     write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                     write(*,"(A,I,2(A,F))")"     Iterations: ",Niter," Chi^2: ",chi,", Wstart: ",wmats(1)
                     write(*,"(A,100E16.6)")"     Moments before fit: ",Moments(iorb,ispin,0:MaxMom)
                     write(*,"(A,100E16.6)")"     Moments after fit:  ",MomentsRebuilt(0:MaxMom)
                  endif
                  !
                  Moments(iorb,ispin,0:MaxMom)=MomentsRebuilt(0:MaxMom)
                  !
               enddo
            enddo
            !
         case("Boson")
            !
            if(verb_)write(*,"(A)")"     Fitting even moments [0:"//str(MaxMom)//"]."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = BosonicFreqMesh(Beta,Nfreq) + (Wstart-1)*2d0*pi/Beta
            allocate(ParamVec(Nfit-2));ParamVec=0d0                             !starts from exponent 0 and interanlly reconstruct the last two
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,Wstart:Wstart+Nfreq-1,ispin)
                  ParamVec = Moments(iorb,ispin,0:MaxMom-2)                     !starts from exponent 0 and interanlly reconstruct the last two
                  ContinuityConstraint = [Component(1),dcmplx(wmats(1),0d0)]
                  !
                  call fit_wrapper(chi2_W_Moments,ParamVec,chi,Niter)
                  !
                  !This sucks but in order to start from 0 in the subr I need to pass an allocatable which cannot be sliced
                  MomentsRebuilt(0:MaxMom-2)=ParamVec
                  call reconstruct_lastCoeffs(MomentsRebuilt,dreal(ContinuityConstraint))
                  !
                  if(verb_)then
                     write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                     write(*,"(A,I,2(A,F))")"     Iterations: ",Niter," Chi^2: ",chi,", Wstart: ",wmats(1)
                     write(*,"(A,100E16.6)")"     Moments before fit: ",Moments(iorb,ispin,0:MaxMom)
                     write(*,"(A,100E16.6)")"     Moments after fit:  ",MomentsRebuilt(0:MaxMom)
                  endif
                  !
                  Moments(iorb,ispin,0:MaxMom)=MomentsRebuilt(0:MaxMom)
                  !
               enddo
            enddo
            !
      end select
      !
      deallocate(Component,ParamVec,MomentsRebuilt,wmats)
      MomentsOut = Moments
      call dump_Moments(dirpath,paramFile)
      !
      if(present(filename))then
         !
         allocate(wmats(size(funct,dim=2)));wmats=0d0
         select case(reg(FitMode))
            case("Green","Sigma")
               wmats = FermionicFreqMesh(Beta,size(funct,dim=2))
            case("Boson")
               wmats = BosonicFreqMesh(Beta,size(funct,dim=2))
         end select
         !
         allocate(funct_print(size(funct,dim=2)))
         do ispin=1,Nspin
            do iorb=1,Norb
               funct_print=czero
               select case(reg(FitMode))
                  case("Green")
                     funct_print = G_Moments(Moments(iorb,ispin,2:Nfit),wmats)
                  case("Sigma")
                     funct_print = S_Moments(Moments(iorb,ispin,:),wmats)
                  case("Boson")
                     funct_print = W_Moments(Moments(iorb,ispin,:),wmats)
               end select
               !
               call dump_FermionicField(funct_print,reg(dirpath)//"fits/",reg(filename)//"_o"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
               !
            enddo
         enddo
         !
         deallocate(funct_print,wmats)
         !
      endif
      !
      if (allocated(Moments)) deallocate(Moments)
      !
   end subroutine fit_moments



   !===========================================================================!



   !---------------------------------------------------------------------------!
   !PURPOSE: Wrapper around the specific minimization routine used
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine fit_wrapper(funct,param,err,outit)
      !
      implicit none
      !
      external                              :: funct
      real(8),allocatable,intent(inout)     :: param(:)
      real(8),intent(out)                   :: err
      integer,intent(out)                   :: outit
      !
      integer                               :: Npara
      integer                               :: mode,iprint,iexit
      real(8)                               :: hh,dfn
      real(8),allocatable,dimension(:)      :: g,hess,w,xprmt
      !
      !
      if(verbose)write(*,"(A)") "---- fit_wrapper"
      !
      !
      Npara=size(param)
      !
      mode=1       !unity Hessian
      dfn=-.5      !frist iteration reduction of energy
      hh=1d-5
      !
      iexit=0
      iprint=0
      if(verbose)iprint=1
      !
      allocate(g(Npara),hess(Npara*Npara),w(1000*Npara),xprmt(Npara))
      g=0.0d0;hess=0.0d0;w=0.0d0
      xprmt=abs(param)+1.d-14
      !
      err=1e5
      outit=0
      !
      call minimize(funct,Npara,param,err,g,hess,w,dfn,xprmt,hh,cg_Ftol,mode,cg_niter,iprint,iexit,outit)
      !
      deallocate(g,hess,w,xprmt)
      !
   end subroutine fit_wrapper


end module fit
