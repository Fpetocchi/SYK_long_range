module gap_equation

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),allocatable,private              :: omega(:)                         !Phonon energy on logarithmic grid
   real(8),allocatable,private              :: a2F(:)                           !alpha^2*F(\Omega) function
   character(len=12),private                :: Phonons_grid="logarithmic"
   !
   !
   integer,private                          :: Nkpt3_orig(3)
   real(8),allocatable,private              :: kpt_orig(:,:)
   complex(8),allocatable,private           :: Hk_orig(:,:,:)
   !
   integer,private                          :: Nkpt3_Hk(3),Nkpt3_Wk(3)
   real(8),allocatable,private              :: kpt_Hk(:,:),kpt_Wk(:,:)
   !
   logical,private                          :: interpolate_Wk=.false.
   complex(8),allocatable,private           :: Zk(:,:,:)
   !
   real(8),allocatable,target,private       :: weights_1(:,:,:),DoS_1(:)
   real(8),allocatable,target,private       :: weights_2(:,:,:),DoS_2(:)
   !
   complex(8),allocatable,private           :: Wk(:,:,:,:)
   !
   logical,private                          :: initialized=.false.
   logical,private                          :: Phonons_stored=.false.
   logical,private                          :: Wk_stored=.false.
   !
   !public
   real(8),allocatable,public,protected     :: Egrid(:)
   real(8),pointer,public,protected         :: weights_Hk(:,:,:),DoS_Hk(:)
   real(8),pointer,public,protected         :: weights_Wk(:,:,:),DoS_Wk(:)
   logical,public,protected                 :: calc_phonons=.false.
   logical,public,protected                 :: calc_Int_static=.false.
   logical,public,protected                 :: calc_Int_full=.false.
   !
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   public :: Initialize_inputs
   public :: store_Wk
   public :: calc_Kel_stat_e
   public :: calc_Kel_dyn_e
   public :: calc_Zph_e
   public :: calc_Kph_e

   !===========================================================================!

contains


   !======================== SETTING UP GLOBAL VARS ===========================!


   !---------------------------------------------------------------------------!
   !PURPOSE: Read phonons and initialize only mesh and energy grids
   !---------------------------------------------------------------------------!
   subroutine Initialize_inputs(pathINPUT,mode_ph,mode_el,Nkpt3,kpt,Hk,Nreal,wrealMax,Nkpt3_intp_Hk,Nkpt3_intp_Wk)
      !
      use parameters
      use utils_misc
      use crystal, only : calc_irredBZ, wannierinterpolation, tetrahedron_integration
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=*),intent(in)           :: mode_ph
      character(len=*),intent(in)           :: mode_el
      integer,intent(in)                    :: Nreal
      real(8),intent(in)                    :: wrealMax
      integer,intent(in)                    :: Nkpt3(3)
      real(8),intent(in)                    :: kpt(:,:)
      complex(8),intent(in)                 :: Hk(:,:,:)
      integer,intent(in),optional           :: Nkpt3_intp_Hk(3)
      integer,intent(in),optional           :: Nkpt3_intp_Wk(3)
      !
      integer                               :: Norb,Ngrid,Nkpti_dum
      integer,allocatable                   :: kptp_dum(:)
      integer,allocatable                   :: pkpt_dum(:,:,:)
      real(8),allocatable                   :: nkstar_dum(:)
      real(8),allocatable                   :: Ek(:,:)
      complex(8),allocatable                :: Hk_intp(:,:,:)
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- Initialize_inputs"
      !
      !
      !setting up shared original K-mesh and Hk
      Nkpt3_orig = Nkpt3
      kpt_orig = kpt
      Hk_orig = Hk
      !
      Norb = size(Hk,dim=1)
      !
      !compute the energy axis on a logarithmic grid
      Ngrid=Nreal
      if(mod(Ngrid,2).eq.0)Ngrid=Ngrid+1
      if(mod(Ngrid-1,4).ne.0)Ngrid=Ngrid+mod(Ngrid-1,4)
      allocate(Egrid(Ngrid));Egrid=0d0
      Egrid = denspace(2d0*wrealMax,Ngrid,center=.true.)
      !
      !setting up global flags
      select case(reg(mode_ph))
         case default
            !
            stop "Available phonon modes in gap equation: None, Elk, QEspresso."
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
      !
      !compute the interpolated K-grid for Hk and the corresponding weights
      write(*,"(A)")"     Weights calculation for electronic energy integrals."
      Nkpt3_Hk = Nkpt3_orig
      if(present(Nkpt3_intp_Hk))then
         !
         if(any(Nkpt3_intp_Hk.eq.0))then
            !
            write(*,"(A)")"     Interpolation grid for Hk provided but one dimension has Nk=0. Interpolation skipped."
            Nkpt3_Hk = Nkpt3_orig
            kpt_Hk = kpt_orig
            call calc_dispersion(Hk_orig,Ek)
            !
         elseif(all(Nkpt3_intp_Hk.eq.Nkpt3_orig))then
            !
            write(*,"(A)")"     Interpolation grid for Hk provided but equal to the original. Interpolation skipped."
            Nkpt3_Hk = Nkpt3_orig
            kpt_Hk = kpt_orig
            call calc_dispersion(Hk_orig,Ek)
            !
         else
            !
            !compute interpolated kpt for Hk
            Nkpt3_Hk = Nkpt3_intp_Hk
            call calc_irredBZ(reg(pathINPUT),Nkpt3_Hk,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Hk)
            deallocate(kptp_dum,pkpt_dum,nkstar_dum)
            !
            allocate(Hk_intp(Norb,Norb,product(Nkpt3_Hk)));Hk_intp=czero
            call cpu_time(start)
            call wannierinterpolation(Nkpt3_orig,kpt_orig,kpt_Hk,Hk_orig,Hk_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     interpolation to ["//str(Nkpt3_Hk(1))//","//str(Nkpt3_Hk(2))//","//str(Nkpt3_Hk(3))//"] K-grid cpu timing:", finish-start
            call calc_dispersion(Hk_intp,Ek)
            deallocate(Hk_intp)
            !
         endif
         !
      else
         !
         Nkpt3_Hk = Nkpt3_orig
         kpt_Hk = kpt_orig
         call calc_dispersion(Hk_orig,Ek)
         !
      endif
      !
      allocate(weights_1(Ngrid,Norb,product(Nkpt3_Hk)));weights_1=0d0
      allocate(DoS_1(Ngrid));DoS_1=0d0
      !
      call cpu_time(start)
      call tetrahedron_integration(reg(pathINPUT),Ek,Nkpt3_Hk,kpt_Hk,Egrid,weights_out=weights_1,DoS_out=DoS_1)
      call cpu_time(finish)
      write(*,"(A,F)") "     tetrahedron integration cpu timing:", finish-start
      deallocate(Ek)
      !
      DoS_Hk => DoS_1
      weights_Hk => weights_1
      !
      !
      !compute the interpolated K-grid for Wk and the corresponding weights
      if(calc_Int_static.or.calc_Int_full)then
         !
         write(*,"(A)")"     Weights calculation for bosonic energy integrals."
         Nkpt3_Wk = Nkpt3_orig
         if(present(Nkpt3_intp_Wk))then
            !
            if(any(Nkpt3_intp_Wk.eq.0))then
               !
               write(*,"(A)")"     Interpolation grid for Wk provided but one dimension has Nk=0. Interpolation skipped."
               Nkpt3_Wk = Nkpt3_orig
               kpt_Wk = kpt_orig
               call calc_dispersion(Hk_orig,Ek,Z=Zk)
               !
            elseif(all(Nkpt3_intp_Wk.eq.Nkpt3_orig))then
               !
               write(*,"(A)")"     Interpolation grid for Wk provided but equal to the original. Interpolation skipped."
               Nkpt3_Wk = Nkpt3_orig
               kpt_Wk = kpt_orig
               call calc_dispersion(Hk_orig,Ek,Z=Zk)
               !
            else
               !
               !compute interpolated kpt for Hk
               Nkpt3_Wk = Nkpt3_intp_Wk
               call calc_irredBZ(reg(pathINPUT),Nkpt3_Wk,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Wk)
               deallocate(kptp_dum,pkpt_dum,nkstar_dum)
               !
               allocate(Hk_intp(Norb,Norb,product(Nkpt3_Wk)));Hk_intp=czero
               call cpu_time(start)
               call wannierinterpolation(Nkpt3_orig,kpt_orig,kpt_Wk,Hk_orig,Hk_intp)
               call cpu_time(finish)
               write(*,"(A,F)") "     interpolation to ["//str(Nkpt3_Wk(1))//","//str(Nkpt3_Wk(2))//","//str(Nkpt3_Wk(3))//"] K-grid cpu timing:", finish-start
               call calc_dispersion(Hk_intp,Ek,Z=Zk)
               deallocate(Hk_intp)
               !
            endif
            !
         else
            !
            Nkpt3_Wk = Nkpt3_orig
            kpt_Wk = kpt_orig
            call calc_dispersion(Hk_orig,Ek,Z=Zk)
            !
         endif
         !
         if(all(Nkpt3_Hk.eq.Nkpt3_Wk))then
            !
            write(*,"(A)") "     tetrahedron integration skipped:"
            deallocate(Ek)
            !
            DoS_Wk => DoS_1
            weights_Wk => weights_1
            !
         else
            !
            allocate(weights_2(Ngrid,Norb,product(Nkpt3_Wk)));weights_2=0d0
            allocate(DoS_2(Ngrid));DoS_2=0d0
            !
            call cpu_time(start)
            call tetrahedron_integration(reg(pathINPUT),Ek,Nkpt3_Wk,kpt_Wk,Egrid,weights_out=weights_2,DoS_out=DoS_2)
            call cpu_time(finish)
            write(*,"(A,F)") "     tetrahedron integration cpu timing:", finish-start
            deallocate(Ek)
            !
            DoS_Wk => DoS_2
            weights_Wk => weights_2
            !
         endif
         !
      endif
      !
      initialized=.true.
      !
      !
      contains
      !
      subroutine calc_dispersion(H,E,Z)
         use linalg, only : eigh
         implicit none
         complex(8),intent(in)              :: H(:,:,:)
         real(8),allocatable,intent(out)    :: E(:,:)
         complex(8),allocatable,intent(out),optional :: Z(:,:,:)
         integer                            :: ik
         complex(8),allocatable             :: Z_(:,:,:)
         !
         if(allocated(E))deallocate(E)
         allocate(E(size(H,dim=1),size(H,dim=3))); E=0d0
         !
         Z_ = H
         do ik=1,size(H,dim=3)
            call eigh(Z_(:,:,ik),Ek(:,ik))
         enddo
         !
         if(present(Z))then
            if(allocated(Z))deallocate(Z)
            Z = Z_
         endif
         !
      end subroutine calc_dispersion
      !
      !
   end subroutine Initialize_inputs


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
      write(*,"(A)") "---- read_a2F"
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
   !PURPOSE: Store the interpolated and rotated NaNb components of the fully
   !         screened interaction
   !---------------------------------------------------------------------------!
   subroutine store_Wk(Wk_orig,Beta,cutoff)
      !
      use parameters
      use utils_misc
      use linalg, only : tensor_transform
      use crystal, only : wannierinterpolation
      implicit none
      !
      complex(8),intent(in)                 :: Wk_orig(:,:,:,:)
      real(8),intent(in)                    :: Beta
      real(8),intent(in)                    :: cutoff
      !
      integer                               :: Ngrid,Norb,Nbp
      integer                               :: iorb,jorb,ib1,ib2
      integer                               :: iw,wndx,Nmats
      integer                               :: ik,Nkpt_orig,Nkpt
      real(8),allocatable                   :: wmats(:)
      complex(8),allocatable                :: Wk_w(:,:,:)                      ![Nbp,Nbp,cutoff]
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- store_Wk"
      !
      !
      if(.not.initialized)stop "store_Wk: input meshes not initialized. call Initialize_inputs."
      !
      !Various checks
      Ngrid = size(Egrid)
      Nbp = size(Wk_orig,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      Nmats = size(Wk_orig,dim=3)
      Nkpt_orig = product(Nkpt3_orig)
      Nkpt = product(Nkpt3_Wk)
      call assert_shape(Wk_orig,[Nbp,Nbp,Nmats,Nkpt_orig],"store_Wk","Wk_orig")
      call assert_shape(kpt_orig,[3,Nkpt_orig],"store_Wk","kpt_orig")
      call assert_shape(kpt_Wk,[3,Nkpt],"store_Wk","kpt_Hk")
      !
      if(interpolate_Wk)then
         if(Nkpt.eq.Nkpt_orig)stop "store_Wk: something is wrong with the K-point dimension (interpolation)."
      else
         if(Nkpt.ne.Nkpt_orig)stop "store_Wk: something is wrong with the K-point dimension."
      endif
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      allocate(wmats(Nmats));wmats=BosonicFreqMesh(Beta,Nmats)
      wndx = Nmats
      if(cutoff.lt.wmats(Nmats))then
         wndx = minloc(abs(wmats-cutoff),dim=1)
         write(*,"(A,F)") "     interaction frequency cutted at iw_["//str(wndx)//"]=",wmats(wndx)
      endif
      deallocate(wmats)
      !
      !Interpolate Hk to new K-grid - this is done separately for each frequency
      !so we dont need to store another full Wk before the orbital rotation
      if(allocated(Wk))deallocate(Wk)
      allocate(Wk(Norb,Norb,wndx,Nkpt));Wk=czero
      allocate(Wk_w(Nbp,Nbp,Nkpt));Wk_w=czero
      call cpu_time(start)
      do iw=1,wndx
         !
         if(interpolate_Wk)then
            call wannierinterpolation(Nkpt3_orig,kpt_orig,kpt_Wk,Wk_orig(:,:,iw,:),Wk_w)
         else
            Wk_w = Wk_orig(:,:,iw,:)
         endif
         !
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(ik,iorb,jorb,ib1,ib2)
         !$OMP DO
         do ik=1,Nkpt
            !
            call tensor_transform(Wk_w(:,:,ik),PhysicalUelements%Full_Map,Zk(:,:,ik),onlyNaNb=.true.)
            !
            do iorb=1,Norb
               do jorb=1,Norb
                  !
                  ib1 = iorb + Norb*(iorb-1)
                  ib2 = jorb + Norb*(jorb-1)
                  !
                  Wk(iorb,jorb,iw,ik) = Wk_w(ib1,ib2,ik)
                  !
               enddo
            enddo
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
      enddo
      deallocate(Wk_w)
      call cpu_time(finish)
      write(*,"(A,F)") "     Wk preparation cpu timing:", finish-start
      !
      Wk_stored=.true.
      !
   end subroutine store_Wk



   !========================== SETTING UP KERNELS =============================!



   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the electronic kernel averaged on an energy grid considering
   !         only the static limit of the fully screened interaction. Minimal use
   !         of shared variables in order to be able to try different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Kel_stat_e(Egrid,weights,Kel_stat_e,printKpath,printmode)
      !
      use parameters
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(in)                    :: weights(:,:,:)
      complex(8),intent(out)                :: Kel_stat_e(:,:)
      character(len=*),intent(in),optional  :: printKpath
      character(len=*),intent(in),optional  :: printmode
      !
      integer                               :: Efermi_ndx,unit
      integer                               :: iE,iE1,iE2,Ngrid

      integer                               :: iorb,jorb,Norb
      integer                               :: ik,Nkpt,Nmats
      real(8)                               :: DoSnorm,dE
      real(8),allocatable                   :: DoS(:)
      complex(8),allocatable                :: Wk_pvt(:,:,:)
      character(len=12)                     :: printmode_used
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//"---- calc_Kel_stat_e"
      !
      !
      if(.not.Wk_stored)stop "calc_Kel_stat_e: fully screened interaction not stored. call store_Wk."
      !
      !Various checks
      Ngrid = size(Egrid)
      Norb = size(Wk,dim=1)
      Nmats = size(Wk,dim=3)
      Nkpt = product(Nkpt3_Wk)
      call assert_shape(Wk,[Norb,Norb,Nmats,Nkpt],"calc_Kel_stat_e","Wk")
      call assert_shape(weights,[Ngrid,Norb,Nkpt],"calc_weights","weights")
      call assert_shape(Kel_stat_e,[Ngrid,Ngrid],"calc_Kph_e","Kel_stat_e")
      Kel_stat_e=czero
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      if(Egrid(Efermi_ndx).ne.0d0) stop "calc_Kel_dyn_e: the energy grid requires the E=0 point."
      !
      DoSnorm=0
      allocate(DoS(Ngrid));DoS=0d0
      do iE=1,Ngrid
         do iorb=1,Norb
            DoS(iE) = DoS(iE) + sum(weights(iE,iorb,:))
         enddo
      enddo
      do iE=2,Ngrid
         dE = abs(Egrid(iE)-Egrid(iE-1))
         DoSnorm = DoSnorm + ( DoS(iE-1)+DoS(iE) ) * (dE/2d0)
         if(iE.eq.Efermi_ndx)then
            write(*,"(A,F)") "     total occupation from DoS:",DoSnorm
            write(*,"(A,F)") "     DoS at the Fermi level:",DoS(iE)
         endif
      enddo
      write(*,"(A,F)") "     DoS normalization:",DoSnorm
      !
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE1,iE2,ik,iorb,jorb,Wk_pvt)
      !
      allocate(Wk_pvt(Norb,Norb,Nkpt));Wk_pvt=czero
      Wk_pvt = Wk(:,:,1,:)
      !
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         if(DoS(iE1).lt.1d-10)cycle
         do iE2=1,Ngrid
            if(DoS(iE2).lt.1d-10)cycle
            !
            do iorb=1,Norb
               do jorb=1,Norb
                  do ik=1,Nkpt
                     !
                     Kel_stat_e(iE1,iE2) = Kel_stat_e(iE1,iE2) + (weights(iE1,iorb,ik)/DoS(iE1)) * (weights(iE2,jorb,ik)/DoS(iE2)) * Wk_pvt(iorb,jorb,ik)
                     !
                  enddo
               enddo
            enddo
            !
         enddo
      enddo
      !$OMP END DO
      deallocate(Wk_pvt)
      !$OMP END PARALLEL
      deallocate(DoS)
      call cpu_time(finish)
      write(*,"(A,F)") "     Calculation of static electronic Kernel cpu timing:", finish-start
      !
      if(present(printKpath).and.(reg(printmode).ne.'None'))then
         printmode_used="E0"
         if(present(printmode))printmode_used=reg(printmode)
         select case(reg(printmode_used))
            case default
               !
               stop "Available print modes of phononic Kernel: E0, 0E, diag, surf, all."
               !
            case("E0")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_stat_e(iE,Efermi_ndx)),dimag(Kel_stat_e(iE,Efermi_ndx))
               enddo
               close(unit)
               !
            case("0E")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_stat_e(Efermi_ndx,iE)),dimag(Kel_stat_e(Efermi_ndx,iE))
               enddo
               close(unit)
               !
            case("diag")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_stat_e(iE,iE)),dimag(Kel_stat_e(iE,iE))
               enddo
               close(unit)
               !
            case("surf")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dreal(Kel_stat_e(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dimag(Kel_stat_e(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               !
            case("all")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_stat_e(iE,Efermi_ndx)),dimag(Kel_stat_e(iE,Efermi_ndx))
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_stat_e(Efermi_ndx,iE)),dimag(Kel_stat_e(Efermi_ndx,iE))
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_stat_e(iE,iE)),dimag(Kel_stat_e(iE,iE))
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dreal(Kel_stat_e(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_stat_e_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dimag(Kel_stat_e(iE1,iE2))
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
   end subroutine calc_Kel_stat_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the electronic kernel averaged on an energy grid considering
   !         only the static limit of the fully screened interaction. Minimal use
   !         of shared variables in order to be able to try different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Kel_dyn_e(Beta,Egrid,weights,MatsStep,Kel_dyn_e,printKpath,printmode)
      !
      use parameters
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(in)                    :: weights(:,:,:)
      real(8),intent(in)                    :: MatsStep
      complex(8),intent(out)                :: Kel_dyn_e(:,:)
      character(len=*),intent(in),optional  :: printKpath
      character(len=*),intent(in),optional  :: printmode
      !
      integer                               :: Efermi_ndx,unit
      integer                               :: iE,iE1,iE2,ix,Ngrid

      integer                               :: iorb,jorb,Norb
      integer                               :: ik,Nkpt,Nmats
      integer                               :: wndx_a,wndx_b
      real(8)                               :: DoSnorm,dE,E1,E2,dx,fact
      real(8)                               :: s_p,s_m,w_p,w_m
      complex(8)                            :: m,q
      complex(8)                            :: Kel_dyn_x,Kel_dyn_e_p,Kel_dyn_e_m
      real(8),allocatable                   :: DoS(:),wmats(:),xgrid(:)
      complex(8),allocatable                :: Wk_pvt(:,:,:,:)
      complex(8),allocatable                :: Kel_dyn_w(:)
      logical                               :: cond1,cond2
      character(len=12)                     :: printmode_used
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//"---- calc_Kel_dyn_e"
      !
      !
      if(.not.Wk_stored)stop "calc_Kel_dyn_e: fully screened interaction not stored. call store_Wk."
      !
      !Various checks
      Ngrid = size(Egrid)
      Norb = size(Wk,dim=1)
      Nmats = size(Wk,dim=3)
      Nkpt = product(Nkpt3_Wk)
      call assert_shape(Wk,[Norb,Norb,Nmats,Nkpt],"calc_Kel_dyn_e","Wk")
      call assert_shape(weights,[Ngrid,Norb,Nkpt],"calc_weights","weights")
      call assert_shape(Kel_dyn_e,[Ngrid,Ngrid],"calc_Kph_e","Kel_dyn_e")
      Kel_dyn_e=czero
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      if(Egrid(Efermi_ndx).ne.0d0) stop "calc_Kel_dyn_e: the energy grid requires the E=0 point."
      !
      if(MatsStep.le.0d0) stop "calc_Kel_dyn_e: MatsStep must be >0."
      write(*,"(A,F)") "     Bosonic frequency step:",MatsStep
      !
      DoSnorm=0
      allocate(DoS(Ngrid));DoS=0d0
      do iE=1,Ngrid
         do iorb=1,Norb
            DoS(iE) = DoS(iE) + sum(weights(iE,iorb,:))
         enddo
      enddo
      do iE=2,Ngrid
         dE = abs(Egrid(iE)-Egrid(iE-1))
         DoSnorm = DoSnorm + ( DoS(iE-1)+DoS(iE) ) * (dE/2d0)
         if(iE.eq.Efermi_ndx)then
            write(*,"(A,F)") "     total occupation from DoS:",DoSnorm
            write(*,"(A,F)") "     DoS at the Fermi level:",DoS(iE)
         endif
      enddo
      write(*,"(A,F)") "     DoS normalization:",DoSnorm
      !
      allocate(xgrid(Ngrid));xgrid=linspace(-1d0,+1d0,Ngrid)
      dx=abs(xgrid(2)-xgrid(1))
      allocate(wmats(Nmats));wmats=BosonicFreqMesh(Beta,Nmats)
      !
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(PRIVATE),&
      !$OMP SHARED(Wk,Ngrid,DoS,Norb,Nmats,Nkpt,weights,wmats,xgrid,dx,MatsStep,Kel_dyn_e)
      !
      allocate(Kel_dyn_w(Nmats));Kel_dyn_w=czero
      allocate(Wk_pvt(Norb,Norb,Nmats,Nkpt));Wk_pvt=czero
      Wk_pvt = Wk
      !
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         if(DoS(iE1).lt.1d-10)cycle
         do iE2=1,Ngrid
            if(DoS(iE2).lt.1d-10)cycle
            !
            E1=Egrid(iE1)
            E2=Egrid(iE2)
            if((E1.eq.0d0).and.(E2.eq.0d0))cycle
            !
            do iorb=1,Norb
               do jorb=1,Norb
                  do ik=1,Nkpt
                     !
                     Kel_dyn_w = Kel_dyn_w + (weights(iE1,iorb,ik)/DoS(iE1)) * (weights(iE2,jorb,ik)/DoS(iE2)) * (Wk_pvt(iorb,jorb,:,ik)-Wk_pvt(iorb,jorb,1,ik))
                     !
                  enddo
               enddo
            enddo
            !
            !integral over auxiliary variable x for E1+E2
            s_p = sign(1d0,E1+E2)
            Kel_dyn_e_p=czero
            do ix=2,Ngrid-1
               !
               w_p = abs(E1+E2) * (1d0+xgrid(ix))/(1d0-xgrid(ix))
               !
               if(w_p.gt.wmats(Nmats))then
                  Kel_dyn_e_p = Kel_dyn_e_p - 0.5d0 * s_p * (2d0/pi) * Kel_dyn_x * dx / (1d0+xgrid(ix-1)**2)
                  exit
               endif
               !
               wndx_a = floor(w_p/MatsStep)+1
               wndx_b = ceiling(w_p/MatsStep)+1
               !
               Kel_dyn_x = Kel_dyn_w(wndx_a) * ( 1d0 - (w_p-wmats(wndx_a))/(wmats(wndx_b)-wmats(wndx_a)) ) + &
                           Kel_dyn_w(wndx_b) * (w_p-wmats(wndx_a))/(wmats(wndx_b)-wmats(wndx_a))
               !
               fact=1d0
               if(ix.eq.2) fact=0.5d0
               !
               Kel_dyn_e_p = Kel_dyn_e_p + fact * s_p * (2d0/pi) * Kel_dyn_x * dx / (1d0+xgrid(ix)**2)
               !
            enddo
            !
            !adding the tail to w-->inf limit of the interaction
            Kel_dyn_e_p = Kel_dyn_e_p + Kel_dyn_w(Nmats) * (pi/2d0 - atan2(wmats(Nmats),(E1+E2)) )
            !
            !adding the fermi function
            Kel_dyn_e_p = Kel_dyn_e_p * (fermidirac(-E1,Beta)-fermidirac(E2,Beta))
            !
            !integral over auxiliary variable x for E1-E2
            s_m = sign(1d0,E1-E2)
            Kel_dyn_e_m=czero
            do ix=2,Ngrid-1
               !
               w_m = abs(E1-E2) * (1d0+xgrid(ix))/(1d0-xgrid(ix))
               !
               if(w_m.gt.wmats(Nmats))then
                  Kel_dyn_e_m = Kel_dyn_e_m - 0.5d0 * s_m * (2d0/pi) * Kel_dyn_x * dx / (1d0+xgrid(ix-1)**2)
                  exit
               endif
               !
               wndx_a = floor(w_m/(2*pi/Beta))+1
               wndx_b = ceiling(w_m/(2*pi/Beta))+1
               !
               Kel_dyn_x = Kel_dyn_w(wndx_a) * ( 1d0 - (w_m-wmats(wndx_a))/(wmats(wndx_b)-wmats(wndx_a)) ) + &
                           Kel_dyn_w(wndx_b) * (w_m-wmats(wndx_a))/(wmats(wndx_b)-wmats(wndx_a))
               !
               fact=1d0
               if(ix.eq.2) fact=0.5d0
               !
               Kel_dyn_e_m = Kel_dyn_e_m + fact * s_m * (2d0/pi) * Kel_dyn_x * dx / (1d0+xgrid(ix)**2)
               !
            enddo
            !
            !adding the tail to w-->inf limit of the interaction
            Kel_dyn_e_m = Kel_dyn_e_m + s_m * Kel_dyn_w(Nmats) * (pi/2d0 - atan2(wmats(Nmats),(E1-E2)) )
            !
            !adding the fermi function
            Kel_dyn_e_m = Kel_dyn_e_m * (fermidirac(E1,Beta)-fermidirac(E2,Beta))
            !
            cond1 = .not.((E1.eq.0d0).and.(E2.ne.0d0))
            cond2 = .not.((E1.ne.0d0).and.(E2.eq.0d0))
            if(cond1.and.cond2)then
               Kel_dyn_e(iE1,iE2) = ( Kel_dyn_e_m + Kel_dyn_e_p ) * (1d0/(tanh(Beta/2d0*E1)*tanh(Beta/2d0*E2)))
            endif
            !
         enddo
      enddo
      !$OMP END DO
      deallocate(Wk_pvt,Kel_dyn_w)
      !$OMP END PARALLEL
      deallocate(DoS,wmats,xgrid)
      !
      !Filling the lines
      do iE=1,Ngrid
         !
         m = ( Kel_dyn_e(iE,Efermi_ndx-1)-Kel_dyn_e(iE,Efermi_ndx-2) ) / ( Egrid(Efermi_ndx-1)-Egrid(Efermi_ndx-2) )
         q = Kel_dyn_e(iE,Efermi_ndx-2)
         Kel_dyn_e(iE,Efermi_ndx) = m * ( 0d0-Egrid(Efermi_ndx-2) ) + q
         !
         m = ( Kel_dyn_e(Efermi_ndx-1,iE)-Kel_dyn_e(Efermi_ndx-2,iE) ) / ( Egrid(Efermi_ndx-1)-Egrid(Efermi_ndx-2) )
         q = Kel_dyn_e(Efermi_ndx-2,iE)
         Kel_dyn_e(Efermi_ndx,iE) = m * ( 0d0-Egrid(Efermi_ndx-2) ) + q
         !
      enddo
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     Calculation of static electronic Kernel cpu timing:", finish-start
      !
      if(present(printKpath).and.(reg(printmode).ne.'None'))then
         printmode_used="E0"
         if(present(printmode))printmode_used=reg(printmode)
         select case(reg(printmode_used))
            case default
               !
               stop "Available print modes of phononic Kernel: E0, 0E, diag, surf, all."
               !
            case("E0")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_dyn_e(iE,Efermi_ndx)),dimag(Kel_dyn_e(iE,Efermi_ndx))
               enddo
               close(unit)
               !
            case("0E")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_dyn_e(Efermi_ndx,iE)),dimag(Kel_dyn_e(Efermi_ndx,iE))
               enddo
               close(unit)
               !
            case("diag")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_dyn_e(iE,iE)),dimag(Kel_dyn_e(iE,iE))
               enddo
               close(unit)
               !
            case("surf")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dreal(Kel_dyn_e(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dimag(Kel_dyn_e(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               !
            case("all")
               !
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_dyn_e(iE,Efermi_ndx)),dimag(Kel_dyn_e(iE,Efermi_ndx))
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_dyn_e(Efermi_ndx,iE)),dimag(Kel_dyn_e(Efermi_ndx,iE))
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  write(unit,"(2F20.10)")Egrid(iE),dreal(Kel_dyn_e(iE,iE)),dimag(Kel_dyn_e(iE,iE))
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dreal(Kel_dyn_e(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
               unit = free_unit()
               open(unit,file=reg(printKpath)//"Kel_dyn_e_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     write(unit,"(2F20.10)")Egrid(iE1),Egrid(iE2),dimag(Kel_dyn_e(iE1,iE2))
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
   end subroutine calc_Kel_dyn_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic renormalization factor averaged on an energy
   !         grid. Minimal use of shared variables in order to be able to try
   !         different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Zph_e(Beta,Egrid,DoS,Zph_e,mode,printZpath)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
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
      call assert_shape(DoS,[Ngrid],"calc_Zph_e","DoS")
      call assert_shape(Zph_e,[Ngrid],"calc_Zph_e","Zph_e")
      Zph_e=0d0
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
      !$OMP SHARED(Ngrid,Nomega,a2F,omega,Efermi_ndx,Beta,Egrid,DoS,Zph_e,mode_used),&
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
                  a2F_tmp(iomega) = a2F(iomega) * ( J(E1,E2,omega(iomega),Beta) + J(E1,-E2,omega(iomega),Beta) )
                  !
               elseif(reg(mode_used).eq."asym")then
                  !
                  a2F_tmp(iomega) = a2F(iomega) * ( 2d0*Jasym(E1,E2,omega(iomega),Beta) - Iasym(E1,E2,omega(iomega),Beta) )
                  !
               elseif (reg(mode_used).eq."sym") then
                  !
                  a2F_tmp(iomega) = a2F(iomega) * ( Iprime(E1,E2,omega(iomega),Beta) + Iprime(E1,-E2,omega(iomega),Beta) )
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
            dE = abs(Egrid(iE2)-Egrid(iE2-1))/tanh(Beta/2d0*E1)
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
   subroutine calc_Kph_e(Beta,Egrid,DoS,Kph_e,printKpath,printmode)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
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
      call assert_shape(DoS,[Ngrid],"calc_Kph_e","DoS")
      call assert_shape(Kph_e,[Ngrid,Ngrid],"calc_Kph_e","Kph_e")
      Kph_e=0d0
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      DoS_Fermi = DoS(Efermi_ndx)
      write(*,"(A,F)") "     calc_Kph_e: DoS at the Fermi level:",DoS_Fermi
      !
      allocate(a2F_tmp(Nomega));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ngrid,Nomega,a2F,omega,DoS_Fermi,Beta,Egrid,Kph_e,Phonons_grid),&
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
              a2F_tmp(iomega) = a2F(iomega) * ( I(E1,E2,omega(iomega),Beta) - I(E1,-E2,omega(iomega),Beta) )
            enddo
            !
            !Integral over phononic frequency - same scheme regargless from Phonons_grid type
            do iomega=2,Nomega
               dw = abs(omega(iomega)-omega(iomega-1))
               a2F_int = a2F_int + ( a2F_tmp(iomega-1)+a2F_tmp(iomega) ) * (dw/2d0)
            enddo
            !
            Kph_e(ie1,ie2) = (2d0/(tanh(Beta/2d0*E1)*tanh(Beta/2d0*E2))) * a2F_int / DoS_Fermi
            !
         enddo !iE2
      enddo !iE1
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(a2F_tmp)
      !
      if(present(printKpath).and.(reg(printmode).ne.'None'))then
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
   double precision function psmooth(E,Beta)
       double precision,intent(in) :: E,Beta
       psmooth = tanh(5d2*Beta*E)**4
   end function psmooth


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the thermal functions J, Jasym, I, Iasym, dI/de
   !---------------------------------------------------------------------------!
   double precision function J(E1,E2,omega,Beta)
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,Beta
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
         term1 = ( fermidirac(E2,Beta) - fermidirac(E1-omega,Beta) ) / ( E1 - E2 - omega )
         term2 = Beta * fermidirac(E1-omega,Beta) * fermidirac(-E1+omega,Beta)
         Jtilda = -(( fermidirac(E1,Beta) + boseeinstein(omega,Beta) ) / ( E1 - E2 - omega )) * ( term1 - term2 )
         !
      end function Jtilda
      !
   end function J
   !
   double precision function Jasym(E1,E2,omega,Beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,Beta
      real(8)                              :: term1,term2
      !
      term1 = ( fermidirac(E2,Beta) - fermidirac(E1-omega,Beta) ) / ( E1 - E2 - omega )
      term1 = term1 -Beta*fermidirac(E1-omega,Beta) * fermidirac(-E1+omega,Beta) * E1 / ( E2 - omega )
      term1 = -term1 * ( fermidirac(E1,Beta) + boseeinstein(omega,Beta) ) / ( E1 - E2 - omega ) * psmooth(E2-omega,Beta)
      !
      term2 = ( fermidirac(E2,Beta) - fermidirac(E1+omega,Beta) ) / ( E1-E2+omega )
      term2 = term2 -Beta*fermidirac(E1+omega,Beta) * fermidirac(-E1-omega,Beta) * E1 / ( E2 + omega )
      term2 = -term2 * ( fermidirac(E1,Beta) + boseeinstein(-omega,Beta) ) / ( E1 - E2 + omega ) * psmooth(E2+omega,Beta)
      !
      Jasym = term1 - term2
      !
   end function Jasym
   !
   double precision function I(E1,E2,omega,Beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,Beta
      real(8)                              :: term1,term2
      real(8)                              :: bE1,bE2,bo
      !
      bE1 = Beta * E1
      bE2 = Beta * E2
      bo = Beta * omega
      !
      if(abs(bo).lt.1d-5) write(*,"(A)"), "     Warning in I(E1,E2,omega): Beta*omega <1e-5. Carefully check the result for Kph."
      !
      !term 1
      if((E1.ge.0d0).and.((E2+omega).ge.0d0))then
         !
         term1 = ( exp(-be2-bo)-exp(-be1) ) / ( (1d0+exp(-be1)) * (1d0+exp(-be2)) * (1d0-exp(-bo)) * ( E1-E2-omega ) )
         !
      elseif((E1.ge.0d0).and.((E2+omega).lt.0d0))then
         !
         term1 = fermidirac(E2,Beta) * boseeinstein(omega,Beta) * (1d0-exp(be2+bo-be1)) / ( (1d0+exp(-be1)) * ( E1-E2-omega ) )
         !
      elseif((E1.lt.0d0).and.((E2+omega).ge.0))then
         !
         term1 = fermidirac(E1,Beta) * (exp(be1-be2-bo)-1d0) / ( (1d0+exp(-be2)) * (1d0-exp(-bo)) * ( E1-E2-omega ) )
         !
      elseif((E1.lt.0).and.((E2+omega).lt.0))then
         !
         term1 = fermidirac(E1,Beta) * fermidirac(E2,Beta) * boseeinstein(omega,Beta) * ( exp(be1) - exp(be2+bo) ) / ( E1-E2-omega )
         !
      else
         write(*,"(4(A,F))") "E1",E1,"E2",E2,"omega",omega,"Beta",Beta
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
         term2 = fermidirac(E1,Beta) * boseeinstein(omega,Beta) * (1d0-exp(be1+bo-be2)) / ( (1d0+exp(-be2)) * (E1-E2+omega))
         !
      elseif((E2.lt.0d0).and.((E1+omega).ge.0)) then
         !
         term2 = fermidirac(E2,Beta) * (exp(be2-be1-bo)-1d0) / ( (1d0+exp(-be1)) * (1d0-exp(-bo)) * (E1-E2+omega) )
         !
      elseif((E2.lt.0).and.((E1+omega).lt.0)) then
         !
         term2 = fermidirac(E1,Beta) * fermidirac(E2,Beta) * boseeinstein(omega,Beta) * (exp(be2)-exp(be1+bo)) / (E1-E2+omega)
         !
      else
         write(*,"(4(A,F))") "E1",E1,"E2",E2,"omega",omega,"Beta",Beta
         stop "Error in I(E1,E2,omega) - second term."
      endif
      !
      I = term1 - term2
      !
   end function I
   !
   double precision function Iasym(E1,E2,omega,Beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,Beta
      real(8)                              :: term1,term2,term3,term4
      !
      term1 = ( fermidirac(E1,Beta)+boseeinstein(omega,Beta)  ) * &
              ( fermidirac(E2,Beta)-fermidirac(E1-omega,Beta) ) / (E1-E2-omega) * psmooth(E2+omega,Beta) / (E2+omega)
      !
      term2 = ( fermidirac(E1,Beta)+boseeinstein(-omega,Beta) ) * &
              ( fermidirac(E2,Beta)-fermidirac(E1+omega,Beta) ) / (E1-E2+omega) * psmooth(E2-omega,Beta) / (E2-omega)
      !
      term3 = ( fermidirac(E1,Beta)+boseeinstein(omega,Beta)  ) * &
              ( fermidirac(-E2,Beta)-fermidirac(E1-omega,Beta)) / (E1+E2-omega) * psmooth(-E2+omega,Beta)/ (-E2+omega)
      !
      term4 = ( fermidirac(E1,Beta)+boseeinstein(-omega,Beta) ) * &
              ( fermidirac(-E2,Beta)-fermidirac(E1+omega,Beta)) / (E1+E2+omega) * psmooth(-E2-omega,Beta)/ (-E2-omega)
      !
      Iasym = term1 - term2 - term3 + term4
      !
   end function Iasym
   !
   double precision function Iprime(E1,E2,omega,Beta)
      use utils_misc
      implicit none
      real(8),intent(in)                   :: E1,E2,omega,Beta
      real(8)                              :: term1,term2,term3,term4
      !
      term1 = fermidirac(E2,Beta)  * boseeinstein(omega,Beta)  / (E1-E2-omega) * &
      ( -Beta*fermidirac(-E1,Beta) * fermidirac(-E1,Beta) + Beta*fermidirac(-E1,Beta) - fermidirac(-E1,Beta)/(E1-E2-omega) )
      !
      term2 = fermidirac(-E2,Beta) * boseeinstein(-omega,Beta) / (E1-E2-omega) * &
      ( +Beta*fermidirac(E1,Beta)  * fermidirac(E1,Beta)  - Beta*fermidirac(E1,Beta)  - fermidirac(E1,Beta)/(E1-E2-omega)  )
      !
      term3 = fermidirac(-E2,Beta) * boseeinstein(omega,Beta) / (E1-E2+omega) * &
      ( +Beta*fermidirac(E1,Beta)  * fermidirac(E1,Beta)  - Beta*fermidirac(E1,Beta)  - fermidirac(E1,Beta)/(E1-E2+omega)  )
      !
      term4 = fermidirac(E2,Beta)  * boseeinstein(-omega,Beta)  / (E1-E2+omega) * &
      ( -Beta*fermidirac(-E1,Beta) * fermidirac(-E1,Beta) + Beta*fermidirac(-E1,Beta) - fermidirac(-E1,Beta)/(E1-E2+omega) )
      !
      Iprime = term1 + term2 - term3 - term4
      !
   end function Iprime


end module gap_equation
