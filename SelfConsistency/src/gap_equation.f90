module gap_equation

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   !interface name
   !   module procedure name_a                                                  ! Description
   !   module procedure name_b                                                  ! Description
   !end interface name

   !---------------------------------------------------------------------------!
   !PURPOSE: Module custom types
   !---------------------------------------------------------------------------!
   !type mytype
   !   integer                               :: var1
   !   real(8),allocatable                   :: var2
   !   logical                               :: var3
   !end type mytype

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),allocatable,private              :: omega(:)                         ! Phonon energy on logarithmic grid
   real(8),allocatable,private              :: a2F(:)                           ! alpha^2*F(\Omega) function
   character(len=12),private                :: Phonons_grid="logarithmic"
   logical,private                          :: Phonons_stored=.false.
   logical,private                          :: calc_Phonons=.false.
   !
   logical,private                          :: initialized=.false.
   !
   logical,private                          :: calc_Int_static=.false.
   logical,private                          :: calc_Int_full=.false.
   !
   integer,private                          :: Nkpt3_orig(3)
   real(8),allocatable,private              :: kpt_orig(:,:)
   complex(8),allocatable,private           :: Hk_orig(:,:,:)
   !
   logical,private                          :: interpolate_Hk=.false.
   logical,private                          :: interpolate_Wk=.false.
   integer,private                          :: Nkpt3_Hk(3),Nkpt3_Wk(3)
   real(8),allocatable,private              :: kpt_Hk(:,:),kpt_Wk(:,:)
   !
   real(8),allocatable,private              :: Egrid(:)
   complex(8),allocatable,private           :: Zk(:,:,:)
   !
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   public :: read_a2F
   public :: calc_Zph_e
   public :: calc_Kph_e

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read phonons and initialize only mesh and energy grids
   !---------------------------------------------------------------------------!
   subroutine Initialize_inputs(pathINPUT,mode_ph,mode_el,Nreal,wrealMax,Nkpt3_intp_Hk,Nkpt3_intp_Wk)
      !
      use parameters
      use utils_misc
      use crystal, only : calc_irredBZ
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=*),intent(in)           :: mode_ph
      character(len=*),intent(in)           :: mode_el
      integer,intent(in)                    :: Nreal
      real(8),intent(in)                    :: wrealMax
      integer,intent(in),optional           :: Nkpt3_intp_Hk(3)
      integer,intent(in),optional           :: Nkpt3_intp_Wk(3)
      !
      integer                               :: Ngrid,Nkpti_dum
      integer,allocatable                   :: kptp_dum(:)
      integer,allocatable                   :: pkpt_dum(:,:,:)
      real(8),allocatable                   :: nkstar_dum(:)
      !
      !
      if(verbose)write(*,"(A)") "---- Initialize_inputs"
      !
      !
      !setting up global flags
      select case(reg(mode_ph))
         case default
            !
            stop "Available phonon modes in gap equation: Elk, QEspresso."
            !
         case("None")
            !
            write(*,"(A)")"     Phononic Kernel and renormalization not included in the gap equation."
            !
         case("Elk","QEspresso")
            !
            calc_phonons=.true.
            call read_a2F(reg(pathINPUT),reg(mode_ph))
            !
      end select
      !
      select case(reg(mode_el))
         case default
            !
            stop "Available electronic modes in gap equation: static, static+dynamic."
            !
         case("None")
            !
            write(*,"(A)")"     Electronic Kernel not included in the gap equation."
            !
         case("static")
            !
            calc_Int_static=.true.
            !
         case("static+dynamic")
            !
            calc_Int_full=.true.
            !
      end select
      !
      !compute the interpolated K-grid for Hk
      Nkpt3_Hk = Nkpt3_orig
      if(present(Nkpt3_intp_Hk))then
         !
         if(any(Nkpt3_intp_Hk.eq.0))then
            !
            write(*,"(A)")"     Interpolation grid for Hk provided but one dimension has Nk=0. Interpolation skipped."
            Nkpt3_Hk = Nkpt3_orig
            kpt_Hk = kpt_orig
            !
         elseif(all(Nkpt3_intp_Hk.eq.Nkpt3_orig))then
            !
            write(*,"(A)")"     Interpolation grid for Hk provided but equal to the original. Interpolation skipped."
            Nkpt3_Hk = Nkpt3_orig
            kpt_Hk = kpt_orig
            !
         else
            !
            !compute interpolated kpt for Hk
            Nkpt3_Hk = Nkpt3_intp_Hk
            call calc_irredBZ(reg(pathINPUT),Nkpt3_Hk,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Hk)
            deallocate(kptp_dum,pkpt_dum,nkstar_dum)
            interpolate_Hk=.true.
            !
         endif
         !
      else
         !
         Nkpt3_Hk = Nkpt3_orig
         kpt_Hk = kpt_orig
         !
      endif
      !
      !compute the interpolated K-grid for Wk
      Nkpt3_Wk = Nkpt3_orig
      if(present(Nkpt3_intp_Wk))then
         !
         if(any(Nkpt3_intp_Wk.eq.0))then
            !
            write(*,"(A)")"     Interpolation grid for Wk provided but one dimension has Nk=0. Interpolation skipped."
            Nkpt3_Wk = Nkpt3_orig
            kpt_Wk = kpt_orig
            !
         elseif(all(Nkpt3_intp_Wk.eq.Nkpt3_orig))then
            !
            write(*,"(A)")"     Interpolation grid for Wk provided but equal to the original. Interpolation skipped."
            Nkpt3_Wk = Nkpt3_orig
            kpt_Wk = kpt_orig
            !
         else
            !
            !compute interpolated kpt for Wk
            Nkpt3_Wk = Nkpt3_intp_Wk
            call calc_irredBZ(reg(pathINPUT),Nkpt3_Wk,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Wk)
            deallocate(kptp_dum,pkpt_dum,nkstar_dum)
            interpolate_Wk=.true.
            !
         endif
         !
      else
         !
         Nkpt3_Hk = Nkpt3_orig
         kpt_Wk = kpt_orig
         !
      endif
      !
      !compute the energy axis on a logarithmic grid
      Ngrid=Nreal
      if(mod(Ngrid,2).eq.0)Ngrid=Ngrid+1
      if(mod(Ngrid-1,4).ne.0)Ngrid=Ngrid+mod(Ngrid-1,4)
      allocate(Egrid(Ngrid));Egrid=0d0
      Egrid = denspace(2d0*wrealMax,Ngrid,center=.true.)
      !
      initialized=.true.
      !
   end subroutine Initialize_inputs


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate the Hk and produce the DoS
   !---------------------------------------------------------------------------!
   subroutine calc_DoS(pathOUTPUT,Egrid,DoS,weights)
      !
      use parameters
      use utils_misc
      use linalg, only : eigh
      use crystal, only : wannierinterpolation, tetrahedron_integration
      implicit none
      !
      character(len=*),intent(in)           :: pathOUTPUT
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(out)                   :: DoS(:)
      real(8),intent(out)                   :: weights(:,:,:)
      !
      integer                               :: ik
      integer                               :: Ngrid,Norb,Nkpt_orig,Nkpt
      complex(8),allocatable                :: Hk(:,:,:)
      real(8),allocatable                   :: Ek(:,:)
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- calc_DoS"
      !
      !
      if(.not.initialized)stop "calc_DoS: input meshes not initialized. call Initialize_inputs."
      !
      !Various checks
      Ngrid = size(Egrid)
      Norb = size(Hk_orig,dim=1)
      Nkpt_orig = product(Nkpt3_orig)
      Nkpt = product(Nkpt3_Hk)
      call assert_shape(DoS,[Ngrid],"calc_DoS","DoS")
      call assert_shape(Hk_orig,[Norb,Norb,Nkpt_orig],"calc_DoS","Hk_orig")
      call assert_shape(kpt_orig,[3,Nkpt_orig],"calc_DoS","kpt_orig")
      call assert_shape(kpt_Hk,[3,Nkpt],"calc_DoS","kpt_Hk")
      call assert_shape(weights,[Ngrid,Norb,Nkpt],"calc_DoS","weights")
      !
      !Interpolate Hk to new K-grid
      if(interpolate_Hk)then
         if(Nkpt.eq.Nkpt_orig)stop "calc_DoS: something is wrong with the K-point dimension (interpolation)."
         allocate(Hk(Norb,Norb,Nkpt));Hk=czero
         call cpu_time(start)
         call wannierinterpolation(Nkpt3_orig,kpt_orig,kpt_Hk,Hk_orig,Hk)
         call cpu_time(finish)
         write(*,"(A,F)") "     Gap equation: H(fullBZ) --> H(fullBZ_new) cpu timing:", finish-start
      else
         if(Nkpt.ne.Nkpt_orig)stop "calc_DoS: something is wrong with the K-point dimension."
         Hk = Hk_orig
      endif
      !
      !Diagonalize Hk to the LDA basis
      allocate(Ek(Norb,Nkpt));Ek=0d0
      allocate(Zk(Norb,Norb,Nkpt));Zk=czero
      !
      Zk = Hk
      do ik=1,Nkpt
         call eigh(Zk(:,:,ik),Ek(:,ik))
      enddo
      deallocate(Hk)
      !
      !Compute the DoS from with tetrahedron integration
      call cpu_time(start)
      call tetrahedron_integration(reg(pathOUTPUT),Ek,Nkpt3_Hk,kpt_Hk,Egrid,weights_out=weights,DoS_out=DoS)
      call cpu_time(finish)
      write(*,"(A,F)") "     Gap equation: tetrahedron integration cpu timing:", finish-start
      deallocate(Ek)
      !
   end subroutine calc_DoS



   !---------------------------------------------------------------------------!
   !PURPOSE: Read the output of phonon calculations with  Elk or QuantumEspresso
   !---------------------------------------------------------------------------!
   subroutine read_a2F(pathINPUT,mode)
      !
      use parameters
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=*),intent(in)           :: mode
      !
      integer                               :: ih,Header,unit
      integer                               :: ierr,Nlines
      integer                               :: iomega,Nomega
      real(8)                               :: dwf,dwb
      real(8)                               :: ConversionFactor
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_a2F"
      !
      !
      select case(reg(mode))
         case default
            !
            stop "Available phonon inputs: Elk, QEspresso."
            !
         case("Elk")
            !
            ConversionFactor = H2eV
            Header = 1
            !
         case("QEspresso")
            !
            ConversionFactor = Ry2H*H2eV
            Header = 5
            !
      end select
      !
      !reading the number of phonon frequecies form file depending on the format
      call inquireFile(reg(pathINPUT)//"ALPHA2F.DAT",filexists)
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"ALPHA2F.DAT",form="formatted",action="read",position="rewind")
      do ih=1,Header
         read(unit,*)
      enddo
      ierr=0
      Nlines=0
      do while (ierr.eq.0)
         Nlines = Nlines + 1
         read(unit,*,iostat=ierr)
      enddo
      close(unit)
      Nomega = Nlines - 1
      write(*,"(A)") "      a2F(omega) is read. The number of phononic frequency points is: "//str(Nomega)
      !
      allocate(omega(Nomega));omega=0d0
      allocate(a2F(Nomega));a2F=0d0
      !
      !reading phonons form file
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"ALPHA2F.DAT",form="formatted",action="read",position="rewind")
      do ih=1,Header
         read(unit,*)
      enddo
      do iomega=1,Nomega
        read(unit,"(2G18.10)") omega(iomega),a2F(iomega)
        if(omega(iomega).lt.0d0) a2F(iomega)=0d0
      enddo
      close(unit)
      !
      a2F = a2F * ConversionFactor
      !
      iomegaloop:do iomega=2,Nomega-1
         dwf=abs(omega(iomega)-omega(iomega+1))
         dwb=abs(omega(iomega)-omega(iomega-1))
         if(dwf.eq.dwb)then
            Phonons_grid="uniform"
         else
            Phonons_grid="logarithmic"
            exit iomegaloop
         endif
      enddo iomegaloop
      !
      write(*,"(A)") "      a2F(omega) is read. The frequency mesh is: "//reg(Phonons_grid)
      Phonons_stored=.true.
      !
   end subroutine read_a2F


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic renormalization factor averaged on an energy
   !         grid. Minimal use of shared variables in order to be able to try
   !         different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Zph_e(beta,Egrid,DoS,Zph_e,mode,printZpath)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(in)                    :: DoS(:)
      real(8),intent(out)                   :: Zph_e(:)
      character(len=*),intent(in),optional  :: mode
      character(len=*),intent(in),optional  :: printZpath
      !
      integer                               :: Efermi_ndx,unit
      integer                               :: iE,iE1,iE2,Ngrid
      integer                               :: iomega,Nomega
      real(8)                               :: E1,E2,dE,dw
      real(8),allocatable                   :: a2F_tmp(:),a2F_int(:)
      character(len=12)                     :: mode_used
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Zph_e"
      !
      !
      if(.not.Phonons_stored)stop "calc_Zph_e: a2F(omega) is not stored. call read_a2F."
      !
      Nomega = size(omega)
      Ngrid = size(Egrid)
      call assert_shape(Egrid,[Ngrid],"calc_Zph_e","Egrid")
      call assert_shape(DoS,[Ngrid],"calc_Zph_e","DoS")
      call assert_shape(Zph_e,[Ngrid],"calc_Zph_e","Zph_e")
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      !
      mode_used="symrenorm"
      if(present(mode))mode_used=reg(mode)
      !
      select case(reg(mode_used))
         case default
            !
            stop "Available E->0 liumits for Zph_e: symrenorm, asym, sym."
            !
         case("symrenorm")
            !
            write (*,"(A)") "     Zph_e: Renormalized term assuming symmetrical DoS around Ef."
            write (*,"(A)") "            See PhysRevB 72 024545 eqs. 79-81 for details."
            !
         case("asym")
            !
            write (*,"(A)") "     Zph_e: Used for asymmetrical band structures."
            write (*,"(A)") "            Divergence for E->0 smoothed numerically."
            !
         case("sym")
            !
            write (*,"(A)") "     Zph_e: Full term assuming symmetrical DOS around Ef."
            write (*,"(A)") "            See PhysRevB 72 024545 eqs. 77-78 for details."
            !
      end select
      !
      allocate(a2F_tmp(Nomega));a2F_tmp=0d0
      allocate(a2F_int(Ngrid));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ngrid,Nomega,a2F,omega,Efermi_ndx,beta,Egrid,DoS,Zph_e,mode_used),&
      !$OMP PRIVATE(iE1,iE2,E1,E2,a2F_tmp,a2F_int,iomega,dw,dE)
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         !
         a2F_int=0d0
         do iE2=1,Ngrid
            !
            E1=Egrid(iE1)
            E2=Egrid(iE2)
            !
            !multiply a2F by the proper J function
            a2F_tmp=0d0
            do iomega=1,Nomega
               !
               if(reg(mode_used).eq."symrenorm")then
                  !
                  a2F_tmp(iomega) = a2F(iomega) * ( J(E1,E2,omega(iomega),beta) + J(E1,-E2,omega(iomega),beta) )
                  !
               elseif(reg(mode_used).eq."asym")then
                  !
                  a2F_tmp(iomega) = a2F(iomega) * ( 2d0*Jasym(E1,E2,omega(iomega),beta) - Iasym(E1,E2,omega(iomega),beta) )
                  !
               elseif (reg(mode_used).eq."sym") then
                  !
                  a2F_tmp(iomega) = a2F(iomega) * ( Iprime(E1,E2,omega(iomega),beta) + Iprime(E1,-E2,omega(iomega),beta) )
                  !
               endif
               !
            enddo
            !
            !Integral over phononic frequency - same scheme regargless from Phonons_grid type
            do iomega=2,Nomega
               dw = abs(omega(iomega)-omega(iomega-1))
               a2F_int(iE2) = a2F_int(iE2) + ( a2F_tmp(iomega-1)+a2F_tmp(iomega) ) * (dw/2d0)
            enddo
            !
         enddo !iE2
         !
         !Integral over E2 - same scheme regargless from Energy_grid type
         do iE2=2,Ngrid
            !
            dE = abs(Egrid(iE2)-Egrid(iE2-1))/tanh(beta/2d0*Egrid(iE1))
            if(reg(mode_used).eq."asym")then
               Zph_e(iE1) = Zph_e(iE1) + ( a2F_int(iE2-1) + a2F_int(iE2) ) * (dE/2d0) * (DoS(iE2)/DoS(Efermi_ndx))
            else
               Zph_e(iE1) = Zph_e(iE1) + ( a2F_int(iE2-1) + a2F_int(iE2) ) * (dE/2d0)
            endif
            !
         enddo
         !
         !extra minus compared to PhysRevB.72.024545 where they have dropped it
         !compare to PhysRevB.88.014514 instead where they have it
         Zph_e(iE1) = -Zph_e(iE1)
         !
      enddo !iE1
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(a2F_tmp,a2F_int)
      !
      if(present(printZpath))then
         unit = free_unit()
         open(unit,file=reg(printZpath)//"Zph_e.DAT",form="formatted",status="unknown",position="rewind",action="write")
         do iE=1,Ngrid
            write(unit,"(2F20.10)")Egrid(iE),Zph_e(iE)
         enddo
         close(unit)
      endif
      !
   end subroutine calc_Zph_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic kernel averaged on an energy grid. Minimal use
   !         of shared variables in order to be able to try different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Kph_e(beta,Egrid,DoS,Kph_e,printKpath,printmode)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(in)                    :: DoS(:)
      real(8),intent(out)                   :: Kph_e(:,:)
      character(len=*),intent(in),optional  :: printKpath
      character(len=*),intent(in),optional  :: printmode
      !
      integer                               :: Efermi_ndx,unit
      integer                               :: iE,iE1,iE2,Ngrid
      integer                               :: iomega,Nomega
      real(8)                               :: E1,E2,dw
      real(8)                               :: a2F_int,DoS_Fermi
      real(8),allocatable                   :: a2F_tmp(:)
      character(len=12)                     :: printmode_used
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Kph_e"
      !
      !
      if(.not.Phonons_stored)stop "calc_Kph_e: a2F(omega) is not stored. call read_a2F."
      !
      Nomega = size(omega)
      Ngrid = size(Egrid)
      call assert_shape(Egrid,[Ngrid],"calc_Kph_e","Egrid")
      call assert_shape(DoS,[Ngrid],"calc_Kph_e","DoS")
      call assert_shape(Kph_e,[Ngrid,Ngrid],"calc_Kph_e","Kph_e")
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      DoS_Fermi = DoS(Efermi_ndx)
      write(*,"(A,F)") "     calc_Kph_e: DoS at the Fermi level:",DoS_Fermi
      !
      allocate(a2F_tmp(Nomega));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ngrid,Nomega,a2F,omega,DoS_Fermi,beta,Egrid,Kph_e,Phonons_grid),&
      !$OMP PRIVATE(iE1,iE2,E1,E2,a2F_tmp,a2F_int,iomega,dw)
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         do iE2=1,Ngrid
            !
            E1=Egrid(iE1)
            E2=Egrid(iE2)
            !
            !multiply a2F by the proper J function
            a2F_tmp=0d0
            do iomega=1,Nomega
              a2F_tmp(iomega) = a2F(iomega) * ( I(E1,E2,omega(iomega),beta) - I(E1,-E2,omega(iomega),beta) )
            enddo
            !
            !Integral over phononic frequency - same scheme regargless from Phonons_grid type
            do iomega=2,Nomega
               dw = abs(omega(iomega)-omega(iomega-1))
               a2F_int = a2F_int + ( a2F_tmp(iomega-1)+a2F_tmp(iomega) ) * (dw/2d0)
            enddo
            !
            Kph_e(ie1,ie2) = (2d0/(tanh(beta/2d0*Egrid(iE1))*tanh(beta/2d0*Egrid(iE2)))) * a2F_int / DoS_Fermi
            !
         enddo !iE2
      enddo !iE1
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(a2F_tmp)
      !
      if(present(printKpath))then
         printmode_used="E0"
         if(present(printmode))printmode_used=reg(printmode)
         select case(reg(printmode_used))
            case default
               !
               stop "Available print modes of phononic Kernel: E0, diag, surf, all."
               !
            case("E0")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kph_e_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),Kph_e(iE,Efermi_ndx)
               enddo
               close(unit)
               !
            case("diag")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kph_e_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),Kph_e(iE,iE)
               enddo
               close(unit)
               !
            case("surf")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kph_e_surf.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),Kph_e(iE1,iE2)
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               !
            case("all")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kph_e_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),Kph_e(iE,Efermi_ndx)
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kph_e_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),Kph_e(iE,iE)
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kph_e_surf.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),Kph_e(iE1,iE2)
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               !
               !
         end select
         !
      endif
      !
   end subroutine calc_Kph_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the smoothing function for asymmetric DoS
   !---------------------------------------------------------------------------!
   double precision function psmooth(E,beta)
       double precision,intent(in) :: E,beta
       psmooth = tanh(5d2*beta*E)**4
   end function psmooth


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the thermal functions J, Jasym, I, Iasym, dI/de
   !---------------------------------------------------------------------------!
   double precision function J(E1,E2,omega,beta)
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,beta
      !
      J = Jtilda(omega) - Jtilda(-omega)
      !
      contains
      !
      double precision function Jtilda(omega)
         use utils_misc
         implicit none
         real(8),intent(in)                :: omega
         real(8)                           :: term1,term2
         !
         term1 = ( fermidirac(E2,beta) - fermidirac(E1-omega,beta) ) / ( E1 - E2 - omega )
         term2 = beta * fermidirac(E1-omega,beta) * fermidirac(-E1+omega,beta)
         Jtilda = -(( fermidirac(E1,beta) + boseeinstein(omega,beta) ) / ( E1 - E2 - omega )) * ( term1 - term2 )
         !
      end function Jtilda
      !
   end function J
   !
   double precision function Jasym(E1,E2,omega,beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,beta
      real(8)                              :: term1,term2
      !
      term1 = ( fermidirac(E2,beta) - fermidirac(E1-omega,beta) ) / ( E1 - E2 - omega )
      term1 = term1 -beta*fermidirac(E1-omega,beta) * fermidirac(-E1+omega,beta) * E1 / ( E2 - omega )
      term1 = -term1 * ( fermidirac(E1,beta) + boseeinstein(omega,beta) ) / ( E1 - E2 - omega ) * psmooth(E2-omega,beta)
      !
      term2 = ( fermidirac(E2,beta) - fermidirac(E1+omega,beta) ) / ( E1-E2+omega )
      term2 = term2 -beta*fermidirac(E1+omega,beta) * fermidirac(-E1-omega,beta) * E1 / ( E2 + omega )
      term2 = -term2 * ( fermidirac(E1,beta) + boseeinstein(-omega,beta) ) / ( E1 - E2 + omega ) * psmooth(E2+omega,beta)
      !
      Jasym = term1 - term2
      !
   end function Jasym
   !
   double precision function I(E1,E2,omega,beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,beta
      real(8)                              :: term1,term2
      real(8)                              :: bE1,bE2,bo
      !
      bE1 = beta * E1
      bE2 = beta * E2
      bo = beta * omega
      !
      if(abs(bo).lt.1d-5) write(*,"(A)"), "     Warning in I(E1,E2,omega): beta*omega <1e-5. Carefully check the result for Kph."
      !
      !term 1
      if((E1.ge.0d0).and.((E2+omega).ge.0d0))then
         !
         term1 = ( exp(-be2-bo)-exp(-be1) ) / ( (1d0+exp(-be1)) * (1d0+exp(-be2)) * (1d0-exp(-bo)) * ( E1-E2-omega ) )
         !
      elseif((E1.ge.0d0).and.((E2+omega).lt.0d0))then
         !
         term1 = fermidirac(E2,beta) * boseeinstein(omega,beta) * (1d0-exp(be2+bo-be1)) / ( (1d0+exp(-be1)) * ( E1-E2-omega ) )
         !
      elseif((E1.lt.0d0).and.((E2+omega).ge.0))then
         !
         term1 = fermidirac(E1,beta) * (exp(be1-be2-bo)-1d0) / ( (1d0+exp(-be2)) * (1d0-exp(-bo)) * ( E1-E2-omega ) )
         !
      elseif((E1.lt.0).and.((E2+omega).lt.0))then
         !
         term1 = fermidirac(E1,beta) * fermidirac(E2,beta) * boseeinstein(omega,beta) * ( exp(be1) - exp(be2+bo) ) / ( E1-E2-omega )
         !
      else
         write(*,"(4(A,F))") "E1",E1,"E2",E2,"omega",omega,"beta",beta
         stop "Error in I(E1,E2,omega) - first term."
      endif
      !
      !term 2
      if((E2.ge.0d0).and.((E1+omega).ge.0d0))then
         !
         term2 = 1d0/((1d0+exp(-be1))*(1d0+exp(-be2))*(1d0-exp(-bo)))*(exp(-be1-bo)-exp(-be2))/ (E1-E2+omega)
         !
      elseif((E2.ge.0d0).and.((E1+omega).lt.0d0)) then
         !
         term2 = fermidirac(E1,beta) * boseeinstein(omega,beta) * (1d0-exp(be1+bo-be2)) / ( (1d0+exp(-be2)) * (E1-E2+omega))
         !
      elseif((E2.lt.0d0).and.((E1+omega).ge.0)) then
         !
         term2 = fermidirac(E2,beta) * (exp(be2-be1-bo)-1d0) / ( (1d0+exp(-be1)) * (1d0-exp(-bo)) * (E1-E2+omega) )
         !
      elseif((E2.lt.0).and.((E1+omega).lt.0)) then
         !
         term2 = fermidirac(E1,beta) * fermidirac(E2,beta) * boseeinstein(omega,beta) * (exp(be2)-exp(be1+bo)) / (E1-E2+omega)
         !
      else
         write(*,"(4(A,F))") "E1",E1,"E2",E2,"omega",omega,"beta",beta
         stop "Error in I(E1,E2,omega) - second term."
      endif
      !
      I = term1 - term2
      !
   end function I
   !
   double precision function Iasym(E1,E2,omega,beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,beta
      real(8)                              :: term1,term2,term3,term4
      !
      term1 = ( fermidirac(E1,beta)+boseeinstein(omega,beta)  ) * &
              ( fermidirac(E2,beta)-fermidirac(E1-omega,beta) ) / (E1-E2-omega) * psmooth(E2+omega,beta) / (E2+omega)
      !
      term2 = ( fermidirac(E1,beta)+boseeinstein(-omega,beta) ) * &
              ( fermidirac(E2,beta)-fermidirac(E1+omega,beta) ) / (E1-E2+omega) * psmooth(E2-omega,beta) / (E2-omega)
      !
      term3 = ( fermidirac(E1,beta)+boseeinstein(omega,beta)  ) * &
              ( fermidirac(-E2,beta)-fermidirac(E1-omega,beta)) / (E1+E2-omega) * psmooth(-E2+omega,beta)/ (-E2+omega)
      !
      term4 = ( fermidirac(E1,beta)+boseeinstein(-omega,beta) ) * &
              ( fermidirac(-E2,beta)-fermidirac(E1+omega,beta)) / (E1+E2+omega) * psmooth(-E2-omega,beta)/ (-E2-omega)
      !
      Iasym = term1 - term2 - term3 + term4
      !
   end function Iasym
   !
   double precision function Iprime(E1,E2,omega,beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,beta
      real(8)                              :: term1,term2,term3,term4
      !
      term1 = fermidirac(E2,beta)  * boseeinstein(omega,beta)  / (E1-E2-omega) * &
      ( -beta*fermidirac(-E1,beta) * fermidirac(-E1,beta) + beta*fermidirac(-E1,beta) - fermidirac(-E1,beta)/(E1-E2-omega) )
      !
      term2 = fermidirac(-E2,beta) * boseeinstein(-omega,beta) / (E1-E2-omega) * &
      ( +beta*fermidirac(E1,beta)  * fermidirac(E1,beta)  - beta*fermidirac(E1,beta)  - fermidirac(E1,beta)/(E1-E2-omega)  )
      !
      term3 = fermidirac(-E2,beta) * boseeinstein(omega,beta) / (E1-E2+omega) * &
      ( +beta*fermidirac(E1,beta)  * fermidirac(E1,beta)  - beta*fermidirac(E1,beta)  - fermidirac(E1,beta)/(E1-E2+omega)  )
      !
      term4 = fermidirac(E2,beta)  * boseeinstein(-omega,beta)  / (E1-E2+omega) * &
      ( -beta*fermidirac(-E1,beta) * fermidirac(-E1,beta) + beta*fermidirac(-E1,beta) - fermidirac(-E1,beta)/(E1-E2+omega) )
      !
      Iprime = term1 + term2 - term3 - term4
      !
   end function Iprime


end module gap_equation
