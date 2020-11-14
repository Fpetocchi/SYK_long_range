module fit

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: container for density lookup parameters
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
   integer,private                          :: Nfreq
   integer,parameter,private                :: cg_niter=400
   real(8),parameter,private                :: cg_Ftol=1e-8
   real(8),parameter,private                :: hwband=3d0
   real(8),parameter,private                :: noisefact=0.01
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
            do ibath=2,Nh-1
               call random_number(rnd)
               AndPram%Epsk(:,ibath,:)         = -hwband + (ibath-1)*de + (rnd-0.5)*noisefact
               call random_number(rnd)
               AndPram%Epsk(:,Nfit-ibath+1,:) = +hwband - (ibath-1)*de + (rnd-0.5)*noisefact
            enddo
         elseif(mod(Nfit,2)/=0)then
            de=hwband/Nh
            call random_number(rnd)
            AndPram%Epsk(:,Nh+1,:)= 0.0d0 + (rnd-0.5)*noisefact
            do ibath=2,Nh
               call random_number(rnd)
               AndPram%Epsk(:,ibath,:)         = -hwband + (ibath-1)*de + (rnd-0.5)*noisefact
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
      if(.not.AndPram%status) stop "Anderson parameters not initilized."
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
         Delta(iw) = AndParaVec(2*Nfit+1) + sum( Vk(:)*Vk(:)/( dcmplx(0d0,wm(iw)) - Ek(:)) )
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
   subroutine fit_Delta(funct,Beta,Nb,dirpath,paramFile,FitMode,Eloc)
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
      !
      integer                               :: Norb,Niter
      integer                               :: iorb,ispin
      real(8)                               :: chi
      real(8),allocatable                   :: ParamVec(:)
      complex(8),allocatable                :: funct_print(:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- fit_Delta"
      !
      !
      Nfit = Nb
      Norb = size(funct,dim=1)
      Nfreq = size(funct,dim=2)
      call assert_shape(Eloc,[Norb,Nspin],"fit_Delta","Eloc")
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
            wmats = BosonicFreqMesh(Beta,Nfreq)
            allocate(ParamVec(2*Nfit));ParamVec=0d0
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  !Setting component to fit and parameter first guess
                  Component = funct(iorb,:,ispin)
                  ParamVec(1:Nfit) = AndPram%Epsk(iorb,:,ispin)
                  ParamVec(1+Nfit:2*Nfit) = AndPram%Vk(iorb,:,ispin)
                  !
                  call fit_wrapper(chi2_Delta,ParamVec,chi,Niter)
                  !
                  write(*,"(3(A,I3))")"Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                  write(*,"(A,I,A,F)")"Iterations: ",Niter," Chi^2: ",chi
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
            wmats = BosonicFreqMesh(Beta,Nfreq)
            allocate(ParamVec(2*Nfit+1));ParamVec=0d0
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,:,ispin)
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
         deallocate(Component)
         Eloc = AndPram%Eloc
         call dump_AndPrams(dirpath,paramFile)
         allocate(funct_print(Nfreq))
         do ispin=1,Nspin
            do iorb=1,Norb
               !
               funct_print=czero
               funct_print = DeltaAnderson(ParamVec(1:2*Nfit),wmats)
               call dump_FermionicField(funct_print,reg(dirpath)//"fits/","D_And_w_o"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
               !
            enddo
         enddo
         deallocate(funct_print,ParamVec,wmats)
         !
      end subroutine fit_Delta



   !============================== FIT MOMENTS ================================!



   !---------------------------------------------------------------------------!
   !PURPOSE: Setup/initilize the Green's function moments.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine setupMoments(Norb,dirpath,paramFile)
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
      allocate(Moments(Norb,Nfit,Nspin))
      !
      path = reg(dirpath)//reg(paramFile)
      call inquireFile(reg(path),filexists,hardstop=.false.,verb=verbose)
      !
      if(filexists)then
         if(verbose)write(*,"(A)") "     Checking the number of coefficients in "//reg(path)
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind",iostat=ierr)
         read(unit,*) Nfit_read
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
      if(filexists)then
         !
         if(verbose)write(*,"(A)") "     Reading Moments from "//reg(paramFile)
         allocate(ReadLine(Nspin*Norb))
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind")
         read(unit,*) Nfit_read
         do imoment=1,Nfit
            ReadLine=0d0
            read(unit,*) ReadLine
            Moments(:,imoment,1) = ReadLine(1:Norb)
            Moments(:,imoment,2) = ReadLine(1+Norb:2*Norb)
         enddo
         deallocate(ReadLine)
         close(unit)
         !
      else
         !
         if(verbose)write(*,"(A)") "     Initializing Moments."
         Moments = 0d0
         do imoment=1,Nfit
            do iorb=1,Norb
               do ispin=1,Nspin
                  call random_number(rnd)
                  if(imoment.le.2)then
                     Moments(iorb,imoment,ispin) = 1d0 - (rnd*noisefact)*(imoment-1)!(imoment-1)*noisefact!
                  else
                     !odd moments
                     if(mod(imoment,2).eq.0) Moments(iorb,imoment,ispin) = Moments(iorb,2,ispin) - (imoment/(10d0*Nfit)+(rnd*noisefact))!imoment*noisefact !
                     if(mod(imoment,2).eq.1) Moments(iorb,imoment,ispin) = Moments(iorb,2,ispin) - (imoment/(10d0*Nfit)+(rnd*noisefact))!imoment*noisefact !
                  endif
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
      if(.not.allocated(Moments)) stop "Moments not allocated."
      Norb=size(Moments,dim=1)
      !
      !
      call inquireDir(reg(dirpath),filexists,verb=verbose)
      path = reg(dirpath)//reg(paramFile)
      if(verbose)write(*,"(A)") "     Dump "//reg(path)//" (readable)"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",action="write",position="rewind")
      write(unit,"(I5,A)") Nfit," Number of coefficients."
      do imoment=1,Nfit
         write(unit,"(999E20.12)") (Moments(iorb,imoment,1),iorb=1,Norb),(Moments(iorb,imoment,2),iorb=1,Norb)
      enddo
      close(unit)
      !
   end subroutine dump_Moments


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates the generic functional form via a moment representation.
   !         The frequency mesh is not global as the routines are made available
   !         in the main.
   !TEST ON:
   !---------------------------------------------------------------------------!
   function G_Moments(MomentVec,wm) result(Gf)
      use parameters
      implicit none
      real(8),dimension(:)                  :: MomentVec
      real(8),dimension(:)                  :: wm
      complex(8),dimension(size(wm))        :: Gf
      integer                               :: imoment,iw,exp
      !
      Gf=czero
      !
      do iw=1,size(wm)
         !
         !Moments 0,1
         Gf(iw) = 0d0 + 1d0 / dcmplx(0d0,wm(iw))
         !
         do imoment=1,size(MomentVec)
            !exp = 2*imoment + 1  !Moments 3,5,7,9..
            exp = imoment + 1     !Moments 2,3,4,5..
            Gf(iw) = Gf(iw) + MomentVec(imoment)/(dcmplx(0d0,wm(iw))**exp)
         enddo
         !
      enddo
      !
   end function G_Moments
   !
   function S_Moments(MomentVec,wm) result(Sigma)
      use parameters
      implicit none
      real(8),dimension(:)                  :: MomentVec
      real(8),dimension(:)                  :: wm
      complex(8),dimension(size(wm))        :: Sigma
      integer                               :: imoment,iw,exp
      !
      Sigma=czero
      !
      do iw=1,size(wm)
         !
         !Moments 0,1
         Sigma(iw) = MomentVec(1) + MomentVec(2) / dcmplx(0d0,wm(iw))
         !
         do imoment=3,size(MomentVec)
            !exp = 2*imoment - 3 !Moments 3,5,7,9..
            exp = imoment        !Moments 3,4,5,6..
            Sigma(iw) = Sigma(iw) + MomentVec(imoment)/(dcmplx(0d0,wm(iw))**exp)
         enddo
         !
      enddo
      !
   end function S_Moments
   !
   function W_Moments(MomentVec,wm) result(W)
      use parameters
      implicit none
      real(8),dimension(:)                  :: MomentVec
      real(8),dimension(:)                  :: wm
      complex(8),dimension(size(wm))        :: W
      integer                               :: imoment,iw,exp
      !
      W=czero
      !
      do iw=1,size(wm)
         !
         !Moments 0
         W(iw) = MomentVec(1)
         !
         do imoment=2,size(MomentVec)
            exp = 2*imoment - 2 !Moments 2,4,6,8..
            W(iw) = W(iw) + MomentVec(imoment)/(dcmplx(0d0,wm(iw))**exp)
         enddo
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
      Gf = G_Moments(MomentVec,wmats)
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
      Sigma = S_Moments(MomentVec,wmats)
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
      W = W_Moments(MomentVec,wmats)
      !
      chi2=sum(abs(Component(:)-W(:))**2)
      !
   end subroutine chi2_W_Moments


   !---------------------------------------------------------------------------!
   !PURPOSE: Fit a given function [Norn,Nfreq,Nspin] assuming a generic
   !         functional moment formulation.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine fit_moments(funct,Beta,Nb,dirpath,paramFile,FitMode,MomentsOut,printpath,filename)
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
      real(8),intent(inout)                 :: MomentsOut(:,:,:)
      character(len=255),intent(in),optional:: printpath
      character(len=255),intent(in),optional:: filename
      !
      integer                               :: Norb,Niter
      integer                               :: iorb,ispin
      real(8)                               :: chi
      real(8),allocatable                   :: ParamVec(:)
      complex(8),allocatable                :: funct_print(:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- fit_moments"
      !
      !
      Nfit = Nb
      Norb = size(funct,dim=1)
      Nfreq = size(funct,dim=2)
      call assert_shape(MomentsOut,[Norb,Nfit,Nspin],"fit_moments","MomentsOut")
      allocate(Component(Nfreq));Component=czero
      !
      call setupMoments(Norb,dirpath,paramFile)
      !
      select case(reg(FitMode))
         case default
            !
            stop "Available modes for moment fitting: Green, Sigma, Boson."
            !
         case("Green")
            !
            write(*,"(A)")"     Fitting odd moments of the Green's function."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Nfreq)
            allocate(ParamVec(Nfit-2));ParamVec=0d0
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Moments(iorb,1,ispin) = 0d0
                  Moments(iorb,2,ispin) = 1d0
                  Component = funct(iorb,:,ispin)
                  ParamVec = Moments(iorb,3:Nfit,ispin)
                  !
                  call fit_wrapper(chi2_G_Moments,ParamVec,chi,Niter)
                  !
                  write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                  write(*,"(A,I,A,F)")"     Iterations: ",Niter," Chi^2: ",chi
                  write(*,"(A,100F12.6)")"     Moments before fit: ",Moments(iorb,:,ispin)
                  write(*,"(A,100F12.6)")"     Moments after fit:  ",0d0,1d0,ParamVec
                  !
                  Moments(iorb,3:Nfit,ispin) = ParamVec
                  !
               enddo
            enddo
            !
         case("Sigma")
            !
            write(*,"(A)")"     Fitting Generic odd moments."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Nfreq)
            allocate(ParamVec(Nfit));ParamVec=0d0
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,:,ispin)
                  ParamVec = Moments(iorb,:,ispin)
                  !
                  call fit_wrapper(chi2_S_Moments,ParamVec,chi,Niter)
                  !
                  write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                  write(*,"(A,I,A,F)")"     Iterations: ",Niter," Chi^2: ",chi
                  write(*,"(A,100F12.6)")"     Moments before fit: ",Moments(iorb,:,ispin)
                  write(*,"(A,100F12.6)")"     Moments after fit:  ",ParamVec
                  !
                  Moments(iorb,:,ispin) = ParamVec
                  !
               enddo
            enddo
            !
         case("Boson")
            !
            write(*,"(A)")"     Fitting Generic even moments."
            !
            allocate(wmats(Nfreq));wmats=0d0
            wmats = BosonicFreqMesh(Beta,Nfreq)
            allocate(ParamVec(Nfit));ParamVec=0d0
            !
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,:,ispin)
                  ParamVec = Moments(iorb,:,ispin)
                  !
                  call fit_wrapper(chi2_W_Moments,ParamVec,chi,Niter)
                  !
                  write(*,"(3(A,I3))")"     Results for orb: ",iorb," spin: ",ispin," Npara: ",size(ParamVec)
                  write(*,"(A,I,A,F)")"     Iterations: ",Niter," Chi^2: ",chi
                  write(*,"(A,100F12.6)")"     First 5 moments before fit: ",Moments(iorb,0:4,ispin)
                  write(*,"(A,100F12.6)")"     First 5 moments after fit:  ",ParamVec(1:5)
                  !
                  Moments(iorb,:,ispin) = ParamVec
                  !
               enddo
            enddo
            !
      end select
      !
      deallocate(Component)
      MomentsOut = Moments
      call dump_Moments(dirpath,paramFile)
      allocate(funct_print(Nfreq))
      do ispin=1,Nspin
         do iorb=1,Norb
            !
            funct_print=czero
            select case(reg(FitMode))
               case("Green")
                  funct_print = G_Moments(Moments(iorb,3:Nfit,ispin),wmats)
               case("Sigma")
                  funct_print = S_Moments(Moments(iorb,:,ispin),wmats)
               case("Boson")
                  funct_print = W_Moments(Moments(iorb,:,ispin),wmats)
            end select
            !
            if(present(printpath).and.present(filename))then
               call dump_FermionicField(funct_print,reg(printpath),reg(filename)//"_w_o"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
            else
               call dump_FermionicField(funct_print,reg(dirpath)//"fits/","G_Mom_w_o"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
               if(present(printpath))write(*,"(A)")"Warning: Missing filename, printing default."
               if(present(filename))write(*,"(A)")"Warning: Missing printpath, printing default."
            endif
            !
         enddo
      enddo
      deallocate(funct_print,ParamVec,wmats)
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
