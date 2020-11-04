module fit

   !private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: container for density lookup parameters
   !---------------------------------------------------------------------------!
   type AndersonParam
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

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Setup/initilize the Anderson parameters.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine setupAndPrams(Norb,dirpath,paramFile)
      !
      use parameters
      use utils_misc
      use input_vars
      implicit none
      !
      integer,intent(in)                    :: Norb
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      !
      character(len=255)                    :: path
      integer                               :: unit,Nh
      integer                               :: ibath,iorb
      real(8)                               :: de
      real(8),allocatable                   :: ReadLine(:)
      logical                               :: filexists
      real(8)                               :: rnd
      !
      !
      if(verbose)write(*,"(A)") "---- setupAndPrams"
      if(AndPram%status) write(*,"(A)") "Warning: Anderson parameters already initilized."
      !
      !
      if(.not.allocated(AndPram%Eloc))allocate(AndPram%Eloc(Norb,Nspin))
      if(.not.allocated(AndPram%Epsk))allocate(AndPram%Epsk(Norb,Nbath,Nspin))
      if(.not.allocated(AndPram%Vk))allocate(AndPram%Vk(Norb,Nbath,Nspin))
      !
      path = reg(dirpath)//reg(paramFile)
      call inquireFile(reg(path),filexists,hardstop=.false.)
      !
      if(filexists)then
         !
         if(verbose)write(*,"(A)") "     Reading Anderson Parameters from "//reg(paramFile)
         allocate(ReadLine(4))
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind")
         do iorb=1,Norb
            read(unit,*) AndPram%Eloc(iorb,:)
            do ibath=1,Nbath
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
         AndPram%Eloc=0d0
         !bath energies
         Nh=Nbath/2
         if(mod(Nbath,2)==0)then
            de=hwband/max(Nh-1,1)
            call random_number(rnd)
            AndPram%Epsk(:,Nh,:)  = -1.d-3 + (rnd-0.5)*noisefact
            call random_number(rnd)
            AndPram%Epsk(:,Nh+1,:)=  1.d-3 + (rnd-0.5)*noisefact
            do ibath=2,Nh-1
               call random_number(rnd)
               AndPram%Epsk(:,ibath,:)         = -hwband + (ibath-1)*de + (rnd-0.5)*noisefact
               call random_number(rnd)
               AndPram%Epsk(:,Nbath-ibath+1,:) = +hwband - (ibath-1)*de + (rnd-0.5)*noisefact
            enddo
         elseif(mod(Nbath,2)/=0)then
            de=hwband/Nh
            call random_number(rnd)
            AndPram%Epsk(:,:,Nh+1)= 0.0d0 + (rnd-0.5)*noisefact
            do ibath=2,Nh
               call random_number(rnd)
               AndPram%Epsk(:,ibath,:)         = -hwband + (ibath-1)*de + (rnd-0.5)*noisefact
               call random_number(rnd)
               AndPram%Epsk(:,Nbath-ibath+1,:) = +hwband - (ibath-1)*de + (rnd-0.5)*noisefact
            enddo
         endif
         !bath hybridizations
         do ibath=1,Nbath
            call random_number(rnd)
            AndPram%Vk(:,ibath,:) = max(0.1d0,1.d0/sqrt(dble(Nbath))) + (rnd-0.5)*noisefact
         enddo
         !
      endif
      !
      AndPram%status = .true.
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(Beta,Nmats)
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
      use input_vars
      implicit none
      !
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      !
      character(len=255)                    :: path
      integer                               :: unit
      integer                               :: ibath,iorb,Norb
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- dump_AndPrams"
      if(.not.AndPram%status) stop "Anderson parameters not initilized."
      !
      !
      call inquireDir(reg(dirpath),filexists)
      path = reg(dirpath)//reg(paramFile)
      if(verbose)write(*,"(A)") "     Dump "//reg(path)//" (readable)"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",action="write",position="rewind")
      do iorb=1,Norb
         write(unit,"(2E20.12)") AndPram%Eloc(iorb,:)
         do ibath=1,Nbath
            write(unit,"(4E20.12)") AndPram%Epsk(iorb,ibath,1),AndPram%Vk(iorb,ibath,1),AndPram%Epsk(iorb,ibath,2),AndPram%Vk(iorb,ibath,2)
         enddo
         write(unit,*)
      enddo
      close(unit)
      !
   end subroutine dump_AndPrams


   !---------------------------------------------------------------------------!
   !PURPOSE: Setup/initilize the Green's function moments.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine setupMoments(Norb,dirpath,paramFile)
      !
      use parameters
      use utils_misc
      use input_vars
      implicit none
      !
      integer,intent(in)                    :: Norb
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      !
      character(len=255)                    :: path
      integer                               :: unit,imoment
      integer                               :: ibath,ispin,iorb
      real(8),allocatable                   :: ReadLine(:)
      logical                               :: filexists
      real(8)                               :: rnd
      !
      !
      if(verbose)write(*,"(A)") "---- setupMoments"
      !
      !
      if(allocated(Moments)) deallocate(Moments)
      allocate(Moments(Norb,Nbath,Nspin))
      !
      path = reg(dirpath)//reg(paramFile)
      call inquireFile(reg(path),filexists,hardstop=.false.)
      !
      if(filexists)then
         !
         if(verbose)write(*,"(A)") "     Reading Moments from "//reg(paramFile)
         allocate(ReadLine(Nspin*Norb))
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="old",action="read",position="rewind")
         do ibath=1,Nbath
            ReadLine=0d0
            read(unit,*) ReadLine
            Moments(:,ibath,1) = ReadLine(1:Norb)
            Moments(:,ibath,2) = ReadLine(1+Norb:2*Norb)
         enddo
         deallocate(ReadLine)
         close(unit)
         !
      else
         !
         if(verbose)write(*,"(A)") "     Initializing Moments."
         Moments = 0d0
         do imoment=1,Nbath
            do iorb=1,Norb
               do ispin=1,Nspin
                  call random_number(rnd)
                  Moments(iorb,ibath,ispin) = 1d0 - (rnd-0.5)*noisefact
               enddo
            enddo
         enddo
         !
      endif
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(Beta,Nmats)
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
      use input_vars
      implicit none
      !
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      !
      character(len=255)                    :: path
      integer                               :: unit
      integer                               :: ibath,iorb,Norb
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- dump_Moments"
      if(.not.allocated(Moments)) stop "Moments not allocated."
      !
      !
      call inquireDir(reg(dirpath),filexists)
      path = reg(dirpath)//reg(paramFile)
      if(verbose)write(*,"(A)") "     Dump "//reg(path)//" (readable)"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",action="write",position="rewind")
      do ibath=1,Nbath
         write(unit,"(999E20.12)") (Moments(iorb,ibath,1),iorb=1,Norb),(Moments(iorb,ibath,2),iorb=1,Norb)
      enddo
      close(unit)
      !
   end subroutine dump_Moments


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the Anderson Hybridization function of a given orbital
   !         (iorb) from the parameters list [Epsk,Vk,Eo].
   !TEST ON:
   !---------------------------------------------------------------------------!
   function DeltaAnderson(AndParaVec) result(Delta)
      use input_vars
      use parameters
      implicit none
      real(8),dimension(:)                  :: AndParaVec
      complex(8),dimension(Nmats)           :: Delta
      integer                               :: ibath,io,iw
      real(8),dimension(Nbath)              :: Ek,Vk
      !
      !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
      !
      Delta=czero
      !
      do ibath=1,Nbath
         io = ibath
         Ek(ibath) = AndParaVec(io)
      enddo
      do ibath=1,Nbath
         io = Nbath + ibath
         Vk(ibath) = AndParaVec(io)
      enddo
      !
      do iw=1,Nmats
         Delta(iw) = sum( Vk(:)*Vk(:)/( dcmplx(0d0,wmats(iw)) - Ek(:)) )
      enddo
      !
   end function DeltaAnderson
   !
   function EoDeltaAnderson(AndParaVec) result(Delta)
      use input_vars
      use parameters
      implicit none
      real(8),dimension(:)                  :: AndParaVec
      complex(8),dimension(Nmats)           :: Delta
      integer                               :: ibath,io,iw
      real(8),dimension(Nbath)              :: Ek,Vk
      !
      !\Delta_{aa} = Eloc_{aa} + \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
      !
      Delta=czero
      !
      do ibath=1,Nbath
         io = ibath
         Ek(ibath) = AndParaVec(io)
      enddo
      do ibath=1,Nbath
         io = Nbath + ibath
         Vk(ibath) = AndParaVec(io)
      enddo
      !
      do iw=1,Nmats
         Delta(iw) = AndParaVec(2*Nbath+1) + sum( Vk(:)*Vk(:)/( dcmplx(0d0,wmats(iw)) - Ek(:)) )
      enddo
      !
   end function EoDeltaAnderson


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates the generic functional form via a moment representation.
   !         For Green's function the first moment is known and equal to 1d0.
   !TEST ON:
   !---------------------------------------------------------------------------!
   function GfMoments(MomentVec) result(Gf)
      use input_vars
      use parameters
      implicit none
      real(8),dimension(:)                  :: MomentVec
      complex(8),dimension(Nmats)           :: Gf
      integer                               :: imoment,iw
      !
      Gf=czero
      !
      do iw=1,Nmats
         Gf(iw) = 1d0 / dcmplx(0d0,wmats(iw))
         do imoment=1,size(MomentVec)
            Gf(iw) = Gf(iw) + MomentVec(imoment)/(dcmplx(0d0,wmats(iw))**(imoment))
         enddo
      enddo
      !
   end function GfMoments
   !
   function SigmaMoments(MomentVec) result(Sigma)
      use input_vars
      use parameters
      implicit none
      real(8),dimension(:)                  :: MomentVec
      complex(8),dimension(Nmats)           :: Sigma
      integer                               :: imoment,iw
      !
      Sigma=czero
      !
      do iw=1,Nmats
         do imoment=1,size(MomentVec)
            Sigma(iw) = Sigma(iw) + MomentVec(imoment)/(dcmplx(0d0,wmats(iw))**(imoment))
         enddo
      enddo
      !
   end function SigmaMoments


   !---------------------------------------------------------------------------!
   !PURPOSE: Distance between input/output functions minimized by the module.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine chi2_Delta(Npara,AndParaVec,chi2)
      use input_vars
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: AndParaVec
      real(8)                               :: chi2
      complex(8),dimension(Nmats)           :: Delta
      !
      Delta = DeltaAnderson(AndParaVec)
      !
      chi2=sum(abs(Component(:)-Delta(:))**2)
      !
   end subroutine chi2_Delta
   !
   subroutine chi2_ShiftedDelta(Npara,AndParaVec,chi2)
      use input_vars
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: AndParaVec
      real(8)                               :: chi2
      complex(8),dimension(Nmats)           :: ShiftedDelta
      !
      ShiftedDelta = EoDeltaAnderson(AndParaVec)
      !
      chi2=sum(abs(Component(:)-ShiftedDelta(:))**2)
      !
   end subroutine chi2_ShiftedDelta
   !
   subroutine chi2_GfMoments(Npara,MomentVec,chi2)
      use input_vars
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: MomentVec
      real(8)                               :: chi2
      complex(8),dimension(Nmats)           :: Gf
      !
      Gf = GfMoments(MomentVec)
      !
      chi2=sum(abs(Component(:)-Gf(:))**2)
      !
   end subroutine chi2_GfMoments
   !
   subroutine chi2_SigmaMoments(Npara,MomentVec,chi2)
      use input_vars
      use parameters
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: MomentVec
      real(8)                               :: chi2
      complex(8),dimension(Nmats)           :: Sigma
      !
      Sigma = SigmaMoments(MomentVec)
      !
      chi2=sum(abs(Component(:)-Sigma(:))**2)
      !
   end subroutine chi2_SigmaMoments


   !---------------------------------------------------------------------------!
   !PURPOSE: Fit a given function [Norn,Nmats,Nspin] assuming a generic
   !         functional moment formulation.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine fit_moments(funct,dirpath,paramFile,FitMode,MomentsOut)
      !
      use parameters
      use utils_misc
      use input_vars
      use file_io
      implicit none
      !
      complex(8),intent(in)                 :: funct(:,:,:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      character(len=*),intent(in)           :: FitMode
      real(8),intent(inout)                 :: MomentsOut(:,:,:)
      !
      integer                               :: Norb,iorb,ispin
      real(8)                               :: chi
      real(8),allocatable                   :: ParamVec(:)
      complex(8),allocatable                :: funct_print(:)
      !
      !Minimize stuff
      integer                               :: Niter,Npara
      integer                               :: mode,iprint,iexit
      real(8)                               :: hh,dfn
      real(8),allocatable,dimension(:)      :: g,hess,w,xprmt
      !
      !
      if(verbose)write(*,"(A)") "---- fit_moments"
      !
      !
      hh=1d-5
      iexit=0
      mode=1  !unity Hessian
      dfn=-.5 !frist iteration reduction of energy
      !
      if(size(funct,dim=2).ne.Nmats) stop "Wrong Nmats in function to fit."
      Norb = size(funct,dim=1)
      allocate(Component(Nmats));Component=czero
      !
      call setupMoments(Norb,dirpath,paramFile)
      !
      select case(reg(FitMode))
         case default
            !
            stop "Available modes for moment fitting: Green, Generic."
            !
         case("Green")
            !
            write(*,"(A)")new_line("A")//"Fitting moments of the Green's function."
            allocate(ParamVec(Nbath-1));ParamVec=0d0
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Moments(iorb,1,ispin) = 1d0
                  Component = funct(iorb,:,ispin)
                  ParamVec = Moments(iorb,2:Nbath,ispin)
                  !
                  Npara=size(ParamVec)
                  allocate(g(Npara),hess(Npara*Npara),w(100*Npara),xprmt(Npara))
                  g=0.0d0;w=0.0d0
                  xprmt=abs(ParamVec)+1.d-15
                  call minimize(chi2_GfMoments,Npara,ParamVec,chi,g,hess,w,dfn,xprmt,hh,cg_Ftol,mode,cg_niter,iprint,iexit,Niter)
                  deallocate(g,hess,w,xprmt)
                  !
                  write(*,"(3(A,I3))")"Results for orb: ",iorb," spin: ",ispin," Npara: ",Npara
                  write(*,"(A,I,A,F)")"Iterations: ",Niter," Chi^2: ",chi
                  write(*,"(A,100F12.6)")"Moments before fit: ",Moments(iorb,:,ispin)
                  write(*,"(A,100F12.6)")"Moments after fit:  ",1d0,ParamVec
                  !
                  Moments(iorb,1,ispin) = 1d0
                  Moments(iorb,2:Nbath,ispin) = ParamVec
                  !
               enddo
            enddo
            !
         case("Generic")
            !
            write(*,"(A)")new_line("A")//"Fitting Generic moments."
            allocate(ParamVec(Nbath));ParamVec=0d0
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,:,ispin)
                  ParamVec = Moments(iorb,:,ispin)
                  !
                  Npara=size(ParamVec)
                  allocate(g(Npara),hess(Npara*Npara),w(100*Npara),xprmt(Npara))
                  g=0.0d0;w=0.0d0
                  xprmt=abs(ParamVec)+1.d-15
                  call minimize(chi2_SigmaMoments,Npara,ParamVec,chi,g,hess,w,dfn,xprmt,hh,cg_Ftol,mode,cg_niter,iprint,iexit,Niter)
                  deallocate(g,hess,w,xprmt)
                  !
                  write(*,"(3(A,I3))")"Results for orb: ",iorb," spin: ",ispin," Npara: ",Npara
                  write(*,"(A,I,A,F)")"Iterations: ",Niter," Chi^2: ",chi
                  write(*,"(A,100F12.6)")"Moments before fit: ",Moments(iorb,:,ispin)
                  write(*,"(A,100F12.6)")"Moments after fit:  ",ParamVec
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
      allocate(funct_print(Nmats))
      do ispin=1,Nspin
         do iorb=1,Norb
            !
            funct_print=czero
            select case(reg(FitMode))
               case("Green")
                  funct_print = GfMoments(Moments(iorb,:,ispin))
                  call dump_FermionicField(funct_print,reg(dirpath)//"/fits/","GreenMom_"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
               case("Generic")
                  funct_print = SigmaMoments(Moments(iorb,:,ispin))
                  call dump_FermionicField(funct_print,reg(dirpath)//"/fits/","GenerMom_"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
            end select
            !
         enddo
      enddo
      deallocate(funct_print,ParamVec,wmats)
      !
   end subroutine fit_moments


   !---------------------------------------------------------------------------!
   !PURPOSE: Fit a given function [Norn,Nmats,Nspin] assuming an hybridization
   !         like functional form with an addtional shift.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine fit_Delta(funct,dirpath,paramFile,FitMode,Eloc)
      !
      use parameters
      use utils_misc
      use input_vars
      use file_io
      implicit none
      !
      complex(8),intent(in)                 :: funct(:,:,:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: paramFile
      character(len=*),intent(in)           :: FitMode
      real(8),intent(inout)                 :: Eloc(:,:)
      !
      integer                               :: Norb,iorb,ispin
      real(8)                               :: chi
      real(8),allocatable                   :: ParamVec(:)
      complex(8),allocatable                :: funct_print(:)
      !
      !Minimize stuff
      integer                               :: Niter,Npara
      integer                               :: mode,iprint,iexit
      real(8)                               :: hh,dfn
      real(8),allocatable,dimension(:)      :: g,hess,w,xprmt
      !
      !
      if(verbose)write(*,"(A)") "---- fit_Delta"
      !
      !
      hh=1d-5
      iexit=0
      mode=1  !unity Hessian
      dfn=-.5 !frist iteration reduction of energy
      !
      if(size(funct,dim=2).ne.Nmats) stop "Wrong Nmats in function to fit."
      Norb = size(funct,dim=1)
      allocate(Component(Nmats));Component=czero
      !
      call setupAndPrams(Norb,dirpath,paramFile)
      !
      select case(reg(FitMode))
         case default
            !
            stop "Available modes for moment fitting: Standard, Shifted."
            !
         case("Standard")
            !
            write(*,"(A)")new_line("A")//"Fitting Delta(iw)."
            allocate(ParamVec(2*Nbath));ParamVec=0d0
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,:,ispin)
                  ParamVec(1:Nbath) = AndPram%Epsk(iorb,:,ispin)
                  ParamVec(1+Nbath:2*Nbath) = AndPram%Vk(iorb,:,ispin)
                  !
                  Npara=size(ParamVec)
                  allocate(g(Npara),hess(Npara*Npara),w(100*Npara),xprmt(Npara))
                  g=0.0d0;w=0.0d0
                  xprmt=abs(ParamVec)+1.d-15
                  call minimize(chi2_Delta,Npara,ParamVec,chi,g,hess,w,dfn,xprmt,hh,cg_Ftol,mode,cg_niter,iprint,iexit,Niter)
                  deallocate(g,hess,w,xprmt)
                  !
                  !
                  write(*,"(3(A,I3))")"Results for orb: ",iorb," spin: ",ispin," Npara: ",Npara
                  write(*,"(A,I,A,F)")"Iterations: ",Niter," Chi^2: ",chi
                  !
                  AndPram%Epsk(iorb,:,ispin) = ParamVec(1:Nbath)
                  AndPram%Vk(iorb,:,ispin) = ParamVec(1+Nbath:2*Nbath)
                  !
               enddo
            enddo
            !
         case("Shifted")
            !
            write(*,"(A)")new_line("A")//"Fitting Eloc+Delta(iw)."
            allocate(ParamVec(2*Nbath+1));ParamVec=0d0
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  Component = funct(iorb,:,ispin)
                  ParamVec(1:Nbath) = AndPram%Epsk(iorb,:,ispin)
                  ParamVec(1+Nbath:2*Nbath) = AndPram%Vk(iorb,:,ispin)
                  ParamVec(2*Nbath+1) = AndPram%Eloc(iorb,ispin)
                  !
                  Npara=size(ParamVec)
                  allocate(g(Npara),hess(Npara*Npara),w(100*Npara),xprmt(Npara))
                  g=0.0d0;w=0.0d0
                  xprmt=abs(ParamVec)+1.d-15
                  call minimize(chi2_ShiftedDelta,Npara,ParamVec,chi,g,hess,w,dfn,xprmt,hh,cg_Ftol,mode,cg_niter,iprint,iexit,Niter)
                  deallocate(g,hess,w,xprmt)
                  !
                  !
                  write(*,"(3(A,I3))")"Results for orb: ",iorb," spin: ",ispin," Npara: ",Npara
                  write(*,"(A,I,A,F)")"Iterations: ",Niter," Chi^2: ",chi
                  !
                  AndPram%Eloc(iorb,ispin) = ParamVec(2*Nbath+1)
                  AndPram%Epsk(iorb,:,ispin) = ParamVec(1:Nbath)
                  AndPram%Vk(iorb,:,ispin) = ParamVec(1+Nbath:2*Nbath)
                  !
               enddo
            enddo
            !
         end select
         !
         deallocate(Component)
         Eloc = AndPram%Eloc
         call dump_AndPrams(dirpath,paramFile)
         allocate(funct_print(Nmats))
         do ispin=1,Nspin
            do iorb=1,Norb
               !
               funct_print=czero
               funct_print = DeltaAnderson(ParamVec(1:2*Nbath))
               call dump_FermionicField(funct_print,reg(dirpath)//"/fits/","DmatsAnd_"//str(iorb)//"_s"//str(ispin)//".DAT",wmats)
               !
            enddo
         enddo
         deallocate(funct_print,ParamVec,wmats)
         !
      end subroutine fit_Delta


      !========+=========+=========+=========+=========+=========+=========+=$
      ! PROGRAM: minimize
      ! TYPE   : subroutine ugly as hell
      ! PURPOSE: conjugent gradient search
      ! I/O    :
      ! VERSION: 30-Sep-95
      ! COMMENT: This is a most reliable conjugent gradient routine! It has
      !          served us well for many years, and is capable to cope with
      !          a very large number of variables. Unfortunately, we don't
      !          know who wrote this routine (original name: 'va10a'), and
      !          we find it very obscure. Don't worry, it works just fine.
      !========+=========+=========+=========+=========+=========+=========+=$
            subroutine minimize (funct, n, x, f, g, h, w, dfn, xm, hh, eps, mode, maxfn, iprint, iexit, itn)
            implicit double precision (a-h,o-z)
            dimension x(*), g(*), h(*), w(*), xm(*)
            !external funct
            data zero, half, one, two /0.0d0, 0.5d0, 1.0d0, 2.0d0/
            if (iprint .ne. 0) write (6,1000)
       1000 format (' entry into minimize')
            np = n + 1
            n1 = n - 1
            nn=(n*np)/2
            is = n
            iu = n
            iv = n + n
            ib = iv + n
            idiff = 1
            iexit = 0
            if (mode .eq. 3) go to 15
            if (mode .eq. 2) go to 10
            ij = nn + 1
            do 5 i = 1, n
            do 6 j = 1, i
            ij = ij - 1
         6  h(ij) = zero
         5  h(ij) = one
            go to 15
        10  continue
            ij = 1
            do 11 i = 2, n
            z = h(ij)
            if (z .le. zero) return
            ij = ij + 1
            i1 = ij
            do 11 j = i, n
            zz = h(ij)
            h(ij) = h(ij) / z
            jk = ij
            ik = i1
            do 12 k = i, j
            jk = jk + np - k
            h(jk) = h(jk) - h(ik) * zz
            ik = ik + 1
        12  continue
            ij = ij + 1
        11  continue
            if (h(ij) .le. zero) return
        15  continue
            ij = np
            dmin = h(1)
            do 16 i = 2, n
            if (h(ij) .ge. dmin) go to 16
            dmin = h(ij)
        16  ij = ij + np - i
            if (dmin .le. zero) return
            z = f
            itn = 0
            call funct (n, x, f)
            ifn = 1
            df = dfn
            if (dfn .eq. zero) df = f - z
            if (dfn .lt. zero) df = abs (df * f)
            if (df .le. zero) df = one
        17  continue
            do 19 i = 1, n
            w(i) = x(i)
        19  continue
            link = 1
            if (idiff - 1) 100, 100, 110
        18  continue
            if (ifn .ge. maxfn) go to 90
        20  continue
            if (iprint .eq. 0) go to 21
            if (mod (itn, iprint) .ne. 0) go to 21
             write (6,1001) itn, ifn
      1001  format (1x,'itn = ',i5,' ifn = ',i5)
            write (6,1002) f
      1002  format (1x,'f = ',e15.7)
            if (iprint .lt. 0) go to 21
            write (6,1003) (x(i), i = 1, n)
      !***
      !***
      1003  format (1x,'x = ',4f15.7 / (5x, 4f15.7))
            write (6,1004) (g(i), i = 1, n)
      1004  format (1x,'g = ',4f15.7 / (5x, 4f15.7))
      !**
      !***
        21  continue
            itn = itn + 1
            w(1) = -g(1)
            do 22 i = 2, n
            ij = i
            i1 = i - 1
            z = -g(i)
            do 23 j = 1, i1
            z = z - h(ij) * w(j)
            ij = ij + n - j
        23  continue
        22  w(i) = z
            w(is+n) = w(n) / h(nn)
            ij = nn
            do 25 i = 1, n1
            ij = ij - 1
            z = zero
            do 26 j = 1, i
            z = z + h(ij) * w(is+np-j)
            ij = ij - 1
        26  continue
        25  w(is+n-i) = w(n-i) / h(ij) - z
            z = zero
            gs0 = zero
            do 29 i = 1, n
            if (z * xm(i) .ge. abs (w(is+i))) go to 28
            z = abs (w(is+i)) / xm(i)
        28  gs0 = gs0 + g(i) * w(is+i)
        29  continue
            aeps = eps / z
            iexit = 2
            if (gs0 .ge. zero) go to 92
            alpha = -two * df / gs0
            if (alpha .gt. one) alpha = one
            ff = f
            tot = zero
            int = 0
            iexit = 1
        30  continue
            if (ifn .ge. maxfn) go to 90
            do 31 i = 1, n
            w(i) = x(i) + alpha * w(is+i)
        31  continue
            call funct (n, w, f1)
            ifn = ifn + 1
            if (f1 .ge. f) go to 40
            f2 = f
            tot = tot + alpha
        32  continue
            do 33 i = 1, n
            x(i) = w(i)
        33  continue
            f = f1
            if (int - 1) 35, 49, 50
        35  continue
            if (ifn .ge. maxfn) go to 90
            do 34 i = 1, n
            w(i) = x(i) + alpha * w(is+i)
        34  continue
            call funct (n, w, f1)
            ifn = ifn + 1
            if (f1 .ge. f) go to 50
            if ((f1 + f2 .ge. f + f) .and. (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
            tot = tot + alpha
            alpha = two * alpha
            go to 32
        40  continue
            if (alpha .lt. aeps) go to 92
            if (ifn .ge. maxfn) go to 90
            alpha = half * alpha
            do 41 i = 1, n
            w(i) = x(i) + alpha * w(is+i)
        41  continue
            call funct (n, w, f2)
            ifn = ifn + 1
            if (f2 .ge. f) go to 45
            tot = tot + alpha
            f = f2
            do 42 i = 1, n
            x(i) = w(i)
        42  continue
            go to 49
        45  continue
            z = 0.1d0
            if (f1 + f .gt. f2 + f2) z = one + half * (f - f1) / (f + f1 - f2 - f2)
            if (z .lt. 0.1d0) z = 0.1d0
            alpha = z * alpha
            int = 1
            go to 30
        49  continue
            if (tot .lt. aeps) go to 92
        50  continue
            alpha = tot
            do 56 i = 1, n
            w(i) = x(i)
            w(ib+i) = g(i)
        56  continue
            link = 2
            if (idiff - 1) 100, 100, 110
        54  continue
            if (ifn .ge. maxfn) go to 90
            gys = zero
            do 55 i = 1, n
            w(i) = w(ib+i)
            gys = gys + g(i) * w(is+i)
        55  continue
            df = ff - f
            dgs = gys - gs0
            if (dgs .le. zero) go to 20
            link = 1
            if (dgs + alpha * gs0 .gt. zero) go to 52
            do 51 i = 1, n
            w(iu + i) = g(i) - w(i)
        51  continue
            sig = one / (alpha * dgs)
            go to 70
        52  continue
            zz = alpha / (dgs - alpha * gs0)
            z = dgs * zz - one
            do 53 i = 1, n
            w(iu+i) = z * w(i) + g(i)
        53  continue
            sig = one / (zz * dgs * dgs)
            go to 70
        60  continue
            link = 2
            do 61 i = 1, n
            w(iu+i) = w(i)
        61  continue
            if (dgs + alpha * gs0 .gt. zero) go to 62
            sig = one / gs0
            go to 70
        62  continue
            sig = -zz
        70  continue
            w(iv+1) = w(iu+1)
            do 71 i = 2, n
            ij = i
            i1 = i - 1
            z = w(iu+i)
            do 72 j = 1, i1
            z = z - h(ij) * w(iv+j)
            ij = ij + n - j
        72  continue
            w(iv+i) = z
        71  continue
            ij = 1
            do 75 i = 1, n
            z = h(ij) + sig * w(iv+i) * w(iv+i)
            if (z .le. zero) z = dmin
            if (z .lt. dmin) dmin = z
            h(ij) = z
            w(ib+i) = w(iv+i) * sig / z
            sig = sig - w(ib+i) * w(ib+i) * z
            ij = ij + np - i
        75  continue
            ij = 1
            do 80 i = 1, n1
            ij = ij + 1
            i1 = i + 1
            do 80 j = i1, n
            w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
            h(ij) = h(ij) + w(ib+i) * w(iu+j)
            ij = ij + 1
        80  continue
            go to (60, 20), link
        90  continue
            iexit = 3
            go to 94
        92  continue
            if (idiff .eq. 2) go to 94
            idiff = 2
            go to 17
        94  continue
            if (iprint .eq. 0) return
            write (6,1005) itn, ifn, iexit
      1005  format (1x,'itn = ',i5, ' ifn = ',i5,' iexit = ',i5)
            write (6,1002) f
            write (6,1003) (x(i), i = 1, n)
            write (6,1004) (g(i), i = 1, n)
            return
       100  continue
            do 101 i = 1, n
            z = hh * xm(i)
            w(i) = w(i) + z
            call funct (n, w, f1)
            g(i) = (f1 - f) / z
            w(i) = w(i) - z
       101  continue
            ifn = ifn + n
            go to (18, 54), link
       110  continue
            do 111 i = 1, n
            z = hh * xm(i)
            w(i) = w(i) + z
            call funct (n, w, f1)
            w(i) = w(i) - z - z
            call funct (n, w, f2)
            g(i) = (f1 - f2) / (two * z)
            w(i) = w(i) + z
       111  continue
            ifn = ifn + n + n
            go to (18, 54), link

            end subroutine minimize








end module fit
