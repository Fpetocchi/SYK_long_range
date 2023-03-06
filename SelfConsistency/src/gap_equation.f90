module gap_equation

   implicit none
   private

   !===========================================================================!

   ! COMMENTS: All the energy scales, included Beta, are converted to the DFT
   !           input energy grid: Hartree for Elk and Rydberg for QuantumEspresso
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface interpFermi
      module procedure interpFermi_Z_d
      module procedure interpFermi_Z_z
      module procedure interpFermi_K_d
      module procedure interpFermi_K_z
   end interface interpFermi

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),allocatable,private              :: omega(:)                         !Phonon energy on logarithmic grid
   real(8),allocatable,private              :: a2F(:)                           !alpha^2*F(\Omega) function
   !
   real(8),private                          :: eV2DFTgrid=1d0                   !this converts eV to the DFT grid
   real(8),private                          :: DFTgrid2eV=1d0
   !
   integer,private                          :: Nkpt3_orig(3)
   real(8),allocatable,private              :: kpt_orig(:,:)
   complex(8),allocatable,private           :: Hk_orig(:,:,:)
   !
   integer,private                          :: Nkpt3_Hk(3),Nkpt3_Wk(3)
   real(8),allocatable,private              :: kpt_Hk(:,:),kpt_Wk(:,:)
   integer,allocatable,private              :: kptdif_Wk(:,:)
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
   logical,private                          :: DFT_DoS_stored=.false.
   logical,private                          :: Wk_stored=.false.
   !
   !public
   real(8),allocatable,public,protected     :: Egrid(:)
   real(8),allocatable,public,protected     :: DoS_DFT(:)
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
      type(SCDFT),intent(in)                :: Inputs
      type(Lattice),intent(in)              :: Lttc


      integer,intent(in)                    :: Nreal
      real(8),intent(in)                    :: wrealMax
      complex(8),intent(in)                 :: Hk_input(:,:,:)
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
      !setting up global flags
      select case(reg(Inputs%mode_ph))
         case default
            !
            stop "Available phonon modes in gap equation: None, Elk, QEspresso."
            !
         case("None")
            !
            write(*,"(A)")"     Phononic Kernel and renormalization not included in the gap equation."
            Ngrid=Nreal
            if(mod(Ngrid,2).eq.0)Ngrid=Ngrid+1
            if(mod(Ngrid-1,4).ne.0)Ngrid=Ngrid+mod(Ngrid-1,4)
            allocate(Egrid(Ngrid));Egrid=0d0
            Egrid=denspace(2d0*wrealMax*eV2DFTgrid,Ngrid,center=.true.)
            !
         case("Elk","QEspresso")
            !
            calc_phonons=.true.
            call read_a2F(reg(pathINPUT),reg(mode_ph))
            call read_DoS_DFT(reg(pathINPUT))
            !
      end select
      !
      select case(reg(Inputs%mode_el))
         case default
            !
            stop "Available electronic modes in gap equation: static, static+dynamic, None."
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
      !setting up shared original K-mesh and Hk
      Norb = size(Hk_input,dim=1)
      !
      !compute the interpolated K-grid for Hk and the corresponding weights
      write(*,"(A)")"     Weights calculation for electronic energy integrals."
      Interpolate_Hk = (.not.(any(Inputs%Nkpt3_intp_Hk.eq.0))) .and. &
                       (.not.(all(Inputs%Nkpt3_intp_Hk.eq.Lttc%Nkpt3)))
      !
      if(Interpolate_Hk)then
         !
         Nkpt3_Hk = Inputs%Nkpt3_intp_Hk
         allocate(Hk_intp(Norb,Norb,product(Nkpt3_Hk)));Hk_intp=czero
         allocate(kpt_Hk(3,product(Nkpt3_Hk)));kpt_Hk=0d0
         !
         call calc_irredBZ(reg(pathINPUT),Nkpt3_Hk,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Hk)
         deallocate(kptp_dum,pkpt_dum,nkstar_dum)
         !
         call cpu_time(start)
         call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,kpt_Hk,(Hk_input*eV2DFTgrid),Hk_intp)
         call cpu_time(finish)
         write(*,"(A,F)") "     interpolation to new ["//str(Nkpt3_Hk(1))//","//str(Nkpt3_Hk(2))//","//str(Nkpt3_Hk(3))//"] K-grid cpu timing:", finish-start
         !
      else
         !
         Nkpt3_Hk = Lttc%Nkpt3
         allocate(Hk_intp(Norb,Norb,product(Nkpt3_Hk)));Hk_intp=czero
         allocate(kpt_Hk(3,product(Nkpt3_Hk)));kpt_Hk=0d0
         !
         kpt_Hk = Lttc%kpt
         Hk_intp = Hk_input * eV2DFTgrid
         !
         write(*,"(A)")"     Interpolation of Hk skipped."
         !
      endif
      !
      allocate(weights_1(Ngrid,Norb,product(Nkpt3_Hk)));weights_1=0d0
      allocate(DoS_1(Ngrid));DoS_1=0d0
      !
      call cpu_time(start)
      call tetrahedron_integration(reg(pathINPUT),Hk_intp,Nkpt3_Hk,kpt_Hk,Egrid,weights_out=weights_1,DoS_out=DoS_1)
      call cpu_time(finish)
      write(*,"(A,F)") "     tetrahedron integration cpu timing:", finish-start
      deallocate(Hk_intp)
      !
      DoS_Hk => DoS_1
      weights_Hk => weights_1
      !
      !compute the interpolated K-grid for Wk and the corresponding weights
      if(calc_Int_static.or.calc_Int_full)then
         !
         write(*,"(A)")"     Weights calculation for bosonic energy integrals."
         Interpolate_Wk = (.not.(any(Inputs%Nkpt3_intp_Wk.eq.0))) .and. &
                          (.not.(all(Inputs%Nkpt3_intp_Wk.eq.Lttc%Nkpt3)))
         !
         if(Interpolate_Wk)then
            !
            Nkpt3_Wk = Inputs%Nkpt3_intp_Wk
            allocate(Hk_intp(Norb,Norb,product(Nkpt3_Wk)));Hk_intp=czero
            allocate(kpt_Wk(3,product(Nkpt3_Wk)));kpt_Wk=0d0
            !
            if(all(Nkpt3_Wk.eq.Nkpt3_Hk))then
               !
               kpt_Wk = kpt_Hk
               write(*,"(A)")"     Interpolation of Hk skipped."
               !
               DoS_Wk => DoS_1
               weights_Wk => weights_1
               write(*,"(A)") "     Tetrahedron integration skipped."
               !
            else
               !
               call calc_irredBZ(reg(pathINPUT),Nkpt3_Wk,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Wk)
               deallocate(kptp_dum,pkpt_dum,nkstar_dum)
               !
               allocate(Hk_intp(Norb,Norb,product(Nkpt3_Wk)));Hk_intp=czero
               call cpu_time(start)
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,kpt_Wk,(Hk_input*eV2DFTgrid),Hk_intp)
               call cpu_time(finish)
               write(*,"(A,F)") "     interpolation to ["//str(Nkpt3_Wk(1))//","//str(Nkpt3_Wk(2))//","//str(Nkpt3_Wk(3))//"] K-grid cpu timing:", finish-start
               !
               allocate(weights_2(Ngrid,Norb,product(Nkpt3_Wk)));weights_2=0d0
               allocate(DoS_2(Ngrid));DoS_2=0d0
               !
               call cpu_time(start)
               call tetrahedron_integration(reg(pathINPUT),Hk_intp,Nkpt3_Wk,kpt_Wk,Egrid,weights_out=weights_2,DoS_out=DoS_2)
               call cpu_time(finish)
               write(*,"(A,F)") "     tetrahedron integration cpu timing:", finish-start
               deallocate(Hk_intp)
               !
               DoS_Wk => DoS_2
               weights_Wk => weights_2
               !
            endif
            !
         else
            !
            Nkpt3_Wk = Lttc%Nkpt3
            allocate(Hk_intp(Norb,Norb,product(Nkpt3_Wk)));Hk_intp=czero
            allocate(kpt_Wk(3,product(Nkpt3_Wk)));kpt_Wk=0d0
            !
            kpt_Wk = Lttc%kpt
            Hk_intp = Hk_input * eV2DFTgrid
            !
            write(*,"(A)")"     Interpolation of Hk skipped."

            kpt_Wk = Lttc%kpt
            Hk_intp = Hk_input * eV2DFTgrid
            !
            DoS_Wk => DoS_1
            weights_Wk => weights_1
            !
         endif


         !keep all the results from Hk
         DoS_Wk => DoS_1
         weights_Wk => weights_1
         call calc_rotation(Hk_intp,Zk)



         if(any(Inputs%Nkpt3_intp_Wk.eq.0))then
            !
            write(*,"(A)")"     Interpolation grid for Hk provided but one dimension has Nk=0. Interpolation skipped."
            Nkpt3_Wk = Lttc%Nkpt3
            kpt_Wk = Lttc%kpt
            Hk_intp = Hk_input * eV2DFTgrid
            !
         elseif(all(Inputs%Nkpt3_intp_Wk.eq.Lttc%Nkpt3))then
            !
            write(*,"(A)")"     Interpolation grid for Hk provided but equal to the original. Interpolation skipped."
            Nkpt3_Wk = Lttc%Nkpt3
            kpt_Wk = Lttc%kpt
            Hk_intp = Hk_input * eV2DFTgrid
            !
         else
            !
            !compute interpolated kpt for Hk
            Nkpt3_Wk = Inputs%Nkpt3_intp_Wk
            call calc_irredBZ(reg(pathINPUT),Nkpt3_Wk,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Wk)
            deallocate(kptp_dum,pkpt_dum,nkstar_dum)
            !
            allocate(Hk_intp(Norb,Norb,product(Nkpt3_Wk)));Hk_intp=czero
            call cpu_time(start)
            call wannierinterpolation(Nkpt3_orig,kpt_orig,kpt_Wk,Hk_orig,Hk_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     interpolation to ["//str(Nkpt3_Wk(1))//","//str(Nkpt3_Wk(2))//","//str(Nkpt3_Wk(3))//"] K-grid cpu timing:", finish-start

            interpolate_Wk=.true.
            !
         endif
         call calc_dispersion(Hk_intp,Ek,Z=Zk)





         Nkpt3_Wk = Nkpt3_orig
         if(present(Nkpt3_intp_Wk))then
            !
            if(any(Nkpt3_intp_Wk.eq.0))then
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
               interpolate_Wk=.true.
               !
            endif
            !
         else
            !
            Nkpt3_Wk = Nkpt3_orig
            kpt_Wk = kpt_orig
            Hk_intp = Hk_orig* SPECIFIC DFT GRID CONVERSION FACTOR EVERYWHERE
            call calc_dispersion(Hk_intp,Ek,Z=Zk)
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
            call tetrahedron_integration(reg(pathINPUT),Hk_intp,Nkpt3_Wk,kpt_Wk,Egrid,weights_out=weights_2,DoS_out=DoS_2)
            call cpu_time(finish)
            write(*,"(A,F)") "     tetrahedron integration cpu timing:", finish-start
            deallocate(Hk_intp,Ek)
            !
            DoS_Wk => DoS_2
            weights_Wk => weights_2
            !
         endif
         !
      endif
      deallocate(Hk_intp)
      !
      initialized=.true.
      !
      !
      !
      contains
      !
      !
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
      integer                               :: ih,unit
      integer                               :: Header,Footer
      integer                               :: ierr,Nlines
      integer                               :: iomega,Nomega
      real(8)                               :: dwf,dwb
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
            eV2DFTgrid = eV2H
            DFTgrid2eV = H2eV
            Header = 0
            Footer = 0
            !
         case("QEspresso")
            !
            eV2DFTgrid = eV2Ry
            DFTgrid2eV = Ry2eV
            Header = 5
            Footer = 1
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
         read(unit,*,iostat=ierr)
         if(ierr.eq.0)Nlines = Nlines + 1
      enddo
      close(unit)
      Nomega = Nlines - Footer
      write(*,"(A)") "     The number of phononic frequency points in a2F(w) is: "//str(Nomega)
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
      if(abs(omega(1)-omega(2)).eq.abs(omega(3)-omega(4)))then
         write(*,"(A)") "     The phononic frequency grid is uniform."
      else
         write(*,"(A)") "     The phononic frequency grid is logarithmic."
      endif
      !
      Phonons_stored=.true.
      !
   end subroutine read_a2F
   !
   subroutine read_DoS_DFT(pathINPUT)
      !
      use parameters
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      !
      integer                               :: ih,unit
      integer                               :: ierr,Nlines
      integer                               :: iE,Ngrid
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_DoS_DFT"
      !
      !
      !reading the number of phonon frequecies form file depending on the format
      call inquireFile(reg(pathINPUT)//"DoS_DFT.DAT",filexists)
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"DoS_DFT.DAT",form="formatted",action="read",position="rewind")
      ierr=0
      Nlines=0
      do while (ierr.eq.0)
         read(unit,*,iostat=ierr)
         if(ierr.eq.0)Nlines = Nlines + 1
      enddo
      close(unit)
      Ngrid = Nlines
      write(*,"(A)") "     The number of points in DFT DoS is: "//str(Ngrid)
      !
      allocate(Egrid(Ngrid));Egrid=0d0
      allocate(DoS_DFT(Ngrid));DoS_DFT=0d0
      !
      !reading phonons form file
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"DoS_DFT.DAT",form="formatted",action="read",position="rewind")
      do iE=1,Ngrid
        read(unit,"(2G18.10)") Egrid(iE),DoS_DFT(iE)
      enddo
      close(unit)
      !
      if(abs(Egrid(1)-Egrid(2)).eq.abs(Egrid(3)-Egrid(4)))then
         write(*,"(A)") "     The DFT DoS energy grid is uniform."
      else
         write(*,"(A)") "     The DFT DoS energy grid is logarithmic."
      endif
      !
      DFT_DoS_stored=.true.
      !
   end subroutine read_DoS_DFT


   !---------------------------------------------------------------------------!
   !PURPOSE: Store the interpolated and rotated NaNb components of the fully
   !         screened interaction
   !---------------------------------------------------------------------------!
   subroutine store_Wk(Wk_orig,Beta,cutoff,pathOUTPUT)
      !
      use parameters
      use utils_misc
      use utils_fields
      use linalg, only : tensor_transform
      use crystal, only : fill_ksumkdiff, wannierinterpolation
      implicit none
      !
      complex(8),intent(in)                 :: Wk_orig(:,:,:,:)
      real(8),intent(in)                    :: Beta
      real(8),intent(in)                    :: cutoff
      character(len=*),intent(in)           :: pathOUTPUT
      !
      integer                               :: Ngrid,Norb,Nbp
      integer                               :: iorb,jorb,ib1,ib2
      integer                               :: iw,wndx,Nmats
      integer                               :: ik1,ik2,iq,Nkpt_orig,Nkpt
      integer,allocatable                   :: kptsum(:,:)
      real(8),allocatable                   :: wmats(:)
      complex(8),allocatable                :: Wk_w(:,:,:),Wk_wq(:,:)
      type(physicalU)                       :: PhysicalUelements
      type(FermionicField)                  :: Wk_io
      real                                  :: start,finish
      logical                               :: Wkexists
      !
      !
      write(*,"(A)") new_line("A")//"---- calc_Wk"
      !
      !
      if(.not.initialized)stop "store_Wk: input meshes not initialized. Call Initialize_inputs."
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
      call fill_ksumkdiff(kpt_Wk,kptsum,kptdif_Wk)
      deallocate(kptsum)
      !
      allocate(wmats(Nmats));wmats=BosonicFreqMesh(Beta,Nmats)
      wndx = Nmats
      if(cutoff.lt.wmats(Nmats))then
         wndx = minloc(abs(wmats-cutoff),dim=1)
         write(*,"(A,F)") "     Interaction frequency cut at iw_["//str(wndx)//"]=",wmats(wndx)*DFTgrid2eV
      endif
      deallocate(wmats)
      !
      !check if Wk is already printed
      call inquireFile(reg(pathOUTPUT)//"Wk_gap_w_k_s1.DAT",Wkexists,hardstop=.false.,verb=verbose)
      !
      if(allocated(Wk))deallocate(Wk)
      allocate(Wk(Norb,Norb,wndx,Nkpt));Wk=czero
      if(Wkexists)then
         !
         call AllocateFermionicField(Wk_io,Norb,wndx,Nkpt=Nkpt)
         call read_FermionicField(Wk_io,reg(pathOUTPUT),"Wk_gap_w",kpt_Wk)
         Wk = Wk_io%wks(:,:,:,:,1)
         call DeallocateField(Wk_io)
         !
      else
         !
         call init_Uelements(Norb,PhysicalUelements)
         !
         !Interpolate Wk to new K-grid - this is done separately for each frequency
         !so we dont need to store another full Wk before the orbital rotation
         allocate(Wk_w(Nbp,Nbp,Nkpt));Wk_w=czero
         call cpu_time(start)
         do iw=1,wndx
            !
            Wk_w=czero
            if(interpolate_Wk)then
               call wannierinterpolation(Nkpt3_orig,kpt_orig,kpt_Wk,Wk_orig(:,:,iw,:),Wk_w)
            else
               Wk_w = Wk_orig(:,:,iw,:)
            endif
            !
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(ik1,ik2,iq,iorb,jorb,ib1,ib2,Wk_wq)
            allocate(Wk_wq(Nbp,Nbp));Wk_wq=czero
            !$OMP DO
            do ik1=1,Nkpt
               do ik2=1,Nkpt
                  !
                  !rotations corresponding to q = (k - k')
                  iq = kptdif_Wk(ik1,ik2)
                  Wk_wq = Wk_w(:,:,iq)
                  call tensor_transform("NN",Wk_wq,PhysicalUelements%Full_Map,Zk(:,:,ik1),Zk(:,:,ik2),Zk(:,:,ik2),Zk(:,:,ik1))
                  !
                  !pick out W_abba in the new basis
                  do iorb=1,Norb
                     do jorb=1,Norb
                        !
                        call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                        Wk(iorb,jorb,iw,iq) = Wk_wq(ib1,ib2)
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
            !$OMP END DO
            deallocate(Wk_wq)
            !$OMP END PARALLEL
            !
         enddo
         deallocate(Wk_w,Wk_wq)
         call cpu_time(finish)
         write(*,"(A,F)") "     Wk preparation cpu timing:", finish-start
         !
         Wk = Wk * eV2DFTgrid
         !
         !this looks orbitally like a Green's function container
         call AllocateFermionicField(Wk_io,Norb,wndx,Nkpt=Nkpt)
         Wk_io%wks(:,:,:,:,1) = Wk
         call dump_FermionicField(Wk_io,reg(pathOUTPUT),"Wk_gap_w",.true.,kpt_Wk,.true.)
         call DeallocateField(Wk_io)
         !
      endif
      !
      Wk_stored=.true.
      !
   end subroutine store_Wk



   !================================ KERNELS ==================================!



   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the electronic kernel averaged on an energy grid considering
   !         only the static limit of the fully screened interaction. Minimal use
   !         of shared variables in order to be able to try different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Kel_stat_e(Beta,Egrid,weights,Kel_stat_e,printKpath,printmode)
      !
      use parameters
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(in)                    :: weights(:,:,:)
      complex(8),intent(out)                :: Kel_stat_e(:,:)
      character(len=*),intent(in),optional  :: printKpath
      character(len=*),intent(in),optional  :: printmode
      !
      integer                               :: Efermi_ndx
      integer                               :: iE,iE1,iE2,Ngrid
      integer                               :: iorb,jorb,Norb
      integer                               :: ik1,ik2,iq,Nkpt,Nmats
      real(8)                               :: DoSnorm,dE,Temp
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
      if(.not.allocated(kptdif_Wk)) stop"calc_Kel_stat_e: kptdif_Wk not allocated"
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
      !$OMP PRIVATE(iE1,iE2,ik1,ik2,iq,iorb,jorb,Wk_pvt)
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
                  do ik1=1,Nkpt
                     do ik2=1,Nkpt
                        !
                        iq = kptdif_Wk(ik1,ik2)
                        !
                        Kel_stat_e(iE1,iE2) = Kel_stat_e(iE1,iE2) + (weights(iE1,iorb,ik1)/DoS(iE1)) * (weights(iE2,jorb,ik2)/DoS(iE2)) * Wk_pvt(iorb,jorb,iq)
                        !
                     enddo
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
      if(present(printKpath).and.(reg(printmode).ne."None"))then
         !
         Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
         printmode_used="E0"
         if(present(printmode))printmode_used=reg(printmode)
         call print_Kernel("electronic",printmode_used,reg(printKpath),"Kel_stat",Temp,Egrid,Kel_stat_e)
         !
      endif
      !
   end subroutine calc_Kel_stat_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the electronic kernel averaged on an energy grid considering
   !         only the static limit of the fully screened interaction. Minimal use
   !         of shared variables in order to be able to try different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Kel_dyn_e(Beta,Egrid,weights,Kel_dyn_e,printKpath,printmode)
      !
      use parameters
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(in)                    :: weights(:,:,:)
      complex(8),intent(out)                :: Kel_dyn_e(:,:)
      character(len=*),intent(in),optional  :: printKpath
      character(len=*),intent(in),optional  :: printmode
      !
      integer                               :: Efermi_ndx
      integer                               :: iorb,jorb,Norb
      integer                               :: ik1,ik2,iq,Nkpt,Nmats
      integer                               :: iE,iE1,iE2,iy,Ngrid
      integer                               :: wndx_a,wndx_b
      real(8)                               :: DoSnorm,Temp,MatsStep
      real(8)                               :: wmax,ymax
      real(8)                               :: E1,E2,DE,DE_m,DE_p
      real(8)                               :: dy,dy_m,dy_p
      real(8)                               :: y_i,y_j
      real(8)                               :: wm,ReW_wm_intp,ImW_wm_intp
      complex(8)                            :: W_wm_i,W_wm_j,Int_i,Int_j
      complex(8)                            :: Kel_dyn_e_p,Kel_dyn_e_m
      real(8),allocatable                   :: DoS(:),wmats(:)
      real(8),allocatable                   :: ygrid_m(:),ygrid_p(:)
      complex(8),allocatable                :: Wee(:)
      complex(8),allocatable                :: Wk_pvt(:,:,:,:)
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
      if(.not.allocated(kptdif_Wk)) stop"calc_Kel_dyn_e: kptdif_Wk not allocated"
      Kel_dyn_e=czero
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
      allocate(wmats(Nmats));wmats=BosonicFreqMesh(Beta,Nmats)
      wmax = wmats(Nmats)
      !
      MatsStep = abs(wmats(2)-wmats(1))
      write(*,"(A,F)") "     Bosonic frequency step:",MatsStep
      !
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(PRIVATE),&
      !$OMP SHARED(Wk,Ngrid,DoS,Norb,Nmats,Nkpt,weights,wmats,MatsStep)
      !
      allocate(Wee(Nmats));Wee=0d0
      allocate(Wk_pvt(Norb,Norb,Nmats,Nkpt));Wk_pvt=czero
      Wk_pvt = Wk
      !
      allocate(ygrid_m(Ngrid));ygrid_m=0d0
      allocate(ygrid_p(Ngrid));ygrid_p=0d0
      !
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         if(DoS(iE1).lt.1d-10)cycle
         do iE2=1,Ngrid
            if(DoS(iE2).lt.1d-10)cycle
            !
            E1=Egrid(iE1)
            E2=Egrid(iE2)
            !
            !avoid the divergence coming from 1/tanh(0) = 1/0
            if((E1.eq.0d0).or.(E2.eq.0d0))cycle
            !avoid the Fermi level since it is interpolated later on
            if((iE1.eq.Efermi_ndx).or.(iE2.eq.Efermi_ndx))cycle
            !
            !bring interaction on the energy gird
            do iorb=1,Norb
               do jorb=1,Norb
                  do ik1=1,Nkpt
                     do ik2=1,Nkpt
                        !
                        iq = kptdif_Wk(ik1,ik2)
                        !
                        Wee = Wee + (weights(iE1,iorb,ik1)/DoS(iE1)) * (weights(iE2,jorb,ik2)/DoS(iE2)) * (Wk_pvt(iorb,jorb,:,iq)-Wk_pvt(iorb,jorb,1,iq))
                        !
                     enddo
                  enddo
               enddo
            enddo
            !
            !integral over auxiliary variable y for DE_m = E1-E2 and DE_p = E1+E2
            DE_m = E1 - E2
            ymax = (wmax-DE_m) / (wmax+DE_m)
            if(ymax.eq.1d0) stop "calc_Kel_dyn_e: divergence will occur for ymax(E1-E2) = 1."
            ygrid_m = linspace(-1d0,ymax,Ngrid)
            dy_m = abs(ygrid_m(2)-ygrid_m(1))
            !
            DE_p = E1 - E2
            ymax = (wmax-DE_p) / (wmax+DE_p)
            if(ymax.eq.1d0) stop "calc_Kel_dyn_e: divergence will occur for ymax(E1+E2) = 1."
            ygrid_p = linspace(-1d0,ymax,Ngrid)
            dy_p = abs(ygrid_p(2)-ygrid_p(1))
            !
            !one loop for both auxiliary y variables
            Kel_dyn_e_m=czero
            Kel_dyn_e_p=czero
            do iy=2,Ngrid
               !
               !-------------------- first term of the sum ---------------------
               DE = DE_m
               dy = dy_m
               y_i = ygrid_m(iy)
               y_j = ygrid_m(iy-1)
               !
               !continous frequency correspnding to iy
               wm = abs(DE) * (1+y_i)/(1-y_i)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               ReW_wm_intp = linear_interp_2y( [wmats(wndx_a),dreal(Wee(wndx_a))] , [wmats(wndx_b),dreal(Wee(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats(wndx_a),dimag(Wee(wndx_a))] , [wmats(wndx_b),dimag(Wee(wndx_b))] , wm )
               W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !continous frequency correspnding to iy-1
               wm = abs(DE) * (1+y_j)/(1-y_j)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               ReW_wm_intp = linear_interp_2y( [wmats(wndx_a),dreal(Wee(wndx_a))] , [wmats(wndx_b),dreal(Wee(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats(wndx_a),dimag(Wee(wndx_a))] , [wmats(wndx_b),dimag(Wee(wndx_b))] , wm )
               W_wm_j = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !integrand for iy and iy-1
               Int_i = W_wm_i / ( 1d0 + y_i**2 )
               Int_j = W_wm_j / ( 1d0 + y_j**2 )
               !
               !trapezoidal integration
               Kel_dyn_e_m = Kel_dyn_e_m + (Int_i+Int_j)*(dy/2d0)
               !
               !
               !------------------- second term of the sum ---------------------
               DE = DE_p
               dy = dy_p
               y_i = ygrid_p(iy)
               y_j = ygrid_p(iy-1)
               !
               !continous frequency correspnding to iy
               wm = abs(DE) * (1+y_i)/(1-y_i)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               ReW_wm_intp = linear_interp_2y( [wmats(wndx_a),dreal(Wee(wndx_a))] , [wmats(wndx_b),dreal(Wee(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats(wndx_a),dimag(Wee(wndx_a))] , [wmats(wndx_b),dimag(Wee(wndx_b))] , wm )
               W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !continous frequency correspnding to iy-1
               wm = abs(DE) * (1+y_j)/(1-y_j)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               ReW_wm_intp = linear_interp_2y( [wmats(wndx_a),dreal(Wee(wndx_a))] , [wmats(wndx_b),dreal(Wee(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats(wndx_a),dimag(Wee(wndx_a))] , [wmats(wndx_b),dimag(Wee(wndx_b))] , wm )
               W_wm_j = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !integrand for iy and iy-1
               Int_i = W_wm_i / ( 1d0 + y_i**2 )
               Int_j = W_wm_j / ( 1d0 + y_j**2 )
               !
               !trapezoidal integration
               Kel_dyn_e_p = Kel_dyn_e_p + (Int_i+Int_j)*(dy/2d0)
               !
            enddo
            !
            !adding the tail to w-->inf limit of the interaction
            Kel_dyn_e_m = Kel_dyn_e_m + Wee(Nmats) * ( pi/2d0 - atan2(wmats(Nmats),DE_m) )
            Kel_dyn_e_p = Kel_dyn_e_p + Wee(Nmats) * ( pi/2d0 - atan2(wmats(Nmats),DE_p) )
            !
            !adding the fermi function differences
            Kel_dyn_e_m = Kel_dyn_e_m * ( fermidirac(+E1,Beta) - fermidirac(+E2,Beta) )*(4d0/pi)
            Kel_dyn_e_p = Kel_dyn_e_p * ( fermidirac(-E1,Beta) - fermidirac(+E2,Beta) )*(4d0/pi)
            !
            !adding the tanh in front
            Kel_dyn_e(iE1,iE2) = ( Kel_dyn_e_m + Kel_dyn_e_p ) / ( tanh(Beta/2d0*E1) * tanh(Beta/2d0*E2) )
            !
         enddo
      enddo
      !$OMP END DO
      deallocate(Wk_pvt,Wee)
      !$OMP END PARALLEL
      deallocate(DoS,wmats,ygrid_m,ygrid_p)
      !
      !Filling the Fermi lines
      call interpFermi(Kel_dyn_e,Egrid,Egrid,Efermi_ndx,Efermi_ndx)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     Calculation of static electronic Kernel cpu timing:", finish-start
      !
      if(present(printKpath).and.(reg(printmode).ne."None"))then
         !
         Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
         printmode_used="E0"
         if(present(printmode))printmode_used=reg(printmode)
         call print_Kernel("electronic",printmode_used,reg(printKpath),"Kel_dyn",Temp,Egrid,Kel_dyn_e)
         !
      endif
      !
   end subroutine calc_Kel_dyn_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic renormalization factor averaged on an energy
   !         grid. Minimal use of shared variables in order to be able to try
   !         different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Zph_e(Beta,Egrid,Zph_e,mode,printZpath)
      !
      use parameters, only : K2eV
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
      real(8),intent(in)                    :: Egrid(:)
      real(8),intent(out)                   :: Zph_e(:)
      character(len=*),intent(in),optional  :: mode
      character(len=*),intent(in),optional  :: printZpath
      !
      integer                               :: Efermi_ndx,unit
      integer                               :: iE,iE1,iE2,Ngrid
      integer                               :: iomega,Nomega
      real(8)                               :: E1,E2,dE,dw,Temp
      real(8),allocatable                   :: a2F_tmp(:),a2F_int(:)
      character(len=12)                     :: mode_used
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Zph_e"
      !
      !
      if(.not.Phonons_stored)stop "calc_Zph_e: a2F(omega) is not stored. call read_a2F."
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      Ngrid = size(Egrid)
      Nomega = size(omega)
      call assert_shape(Zph_e,[Ngrid],"calc_Zph_e","Zph_e")
      !
      mode_used="symrenorm"
      if(present(mode))mode_used=reg(mode)
      !
      if(reg(mode_used).eq."asym")then
         if(.not.allocated(DoS_DFT))stop "calc_Zph_e: DoS_DFT is not stored. call read_a2F."
         if(.not.all(Egrid.eq.Egrid_DFT))stop "calc_Zph_e: mode is asym but the provided Egrid is different than Egrid_DFT."
         if(DoS_DFT(Efermi_ndx).ne.DoS0_DFT)stop "calc_Zph_e: mismatch between DoS_DFT(0) and DoS0_DFT."
      endif
      !
      select case(reg(mode_used))
         case default
            !
            stop "Available E->0 liumits for Zph_e: symrenorm, asym, sym."
            !
         case("symrenorm")
            !
            if(verbose)write (*,"(A)") "     Zph_e: Renormalized term assuming symmetrical DoS around Ef. See PhysRevB 72 024545 eqs. 79-81 for details."
            !
         case("asym")
            !
            if(verbose)write (*,"(A)") "     Zph_e: Used for asymmetrical band structures. Divergence for E->0 smoothed numerically."
            !
         case("sym")
            !
            if(verbose)write (*,"(A)") "     Zph_e: Full term assuming symmetrical DOS around Ef. See PhysRevB 72 024545 eqs. 77-78 for details."
            !
      end select
      !
      Zph_e=0d0
      allocate(a2F_tmp(Nomega));a2F_tmp=0d0
      allocate(a2F_int(Ngrid));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
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
            dE = abs(Egrid(iE2)-Egrid(iE2-1))
            if(reg(mode_used).eq."asym")then
               Zph_e(iE1) = Zph_e(iE1) + ( a2F_int(iE2-1) + a2F_int(iE2) ) * (dE/2d0) * (DoS_DFT(iE2)/DoS0_DFT)
            else
               Zph_e(iE1) = Zph_e(iE1) + ( a2F_int(iE2-1) + a2F_int(iE2) ) * (dE/2d0)
            endif
            !
         enddo
         !
         !extra minus compared to PhysRevB.72.024545 where they have dropped it
         !compare to PhysRevB.88.014514 instead where they have it
         if(E1.ne.0d0)Zph_e(iE1) = -Zph_e(iE1)/tanh(Beta/2d0*E1)
         !
      enddo !iE1
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(a2F_tmp,a2F_int)
      !
      !Filling the Fermi line
      call interpFermi(Zph_e,Egrid,Efermi_ndx)
      !
      write(*,"(A,1E20.10)")"     lambda(Zph): ",Zph_e(Efermi_ndx)
      !
      if(present(printZpath))then
         !
         Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
         unit = free_unit()
         open(unit,file=reg(printZpath)//"Zph_e_T"//str(Temp,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
         do iE=1,Ngrid
            write(unit,"(2F20.10)")Egrid(iE),Zph_e(iE)
         enddo
         close(unit)
         !
      endif
      !
   end subroutine calc_Zph_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic kernel averaged on an energy grid. Minimal use
   !         of shared variables in order to be able to try different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Kph_e(Beta,Egrid,DoS,Kph_e,printKpath,printmode)
      !
      use parameters, only : K2eV
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
      integer                               :: Efermi_ndx1,Efermi_ndx2
      integer                               :: iE1,iE2,Ngrid
      integer                               :: iomega,Nomega
      real(8)                               :: E1,E2,dw,Temp
      real(8)                               :: a2F_int
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
      Efermi_ndx1 = minloc(abs(Egrid),dim=1)
      Efermi_ndx2 = minloc(abs(Egrid_DFT),dim=1)
      if(verbose)write(*,"(A,F)") "     calc_Kph_e: DoS_DFT at the Fermi level:",DoS0_DFT
      !
      allocate(a2F_tmp(Nomega));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE1,iE2,E1,E2,a2F_tmp,a2F_int,iomega,dw)
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         do iE2=1,Ngrid
            !
            E1=Egrid(iE1)
            E2=Egrid_DFT(iE2)
            !
            !multiply a2F by the proper J function
            a2F_tmp=0d0
            do iomega=1,Nomega
              a2F_tmp(iomega) = a2F(iomega) * ( I(E1,E2,omega(iomega),Beta) - I(E1,-E2,omega(iomega),Beta) )
            enddo
            !
            !Integral over phononic frequency - same scheme regargless from Phonons_grid type
            a2F_int=0d0
            do iomega=2,Nomega
               dw = abs(omega(iomega)-omega(iomega-1))
               a2F_int = a2F_int + ( a2F_tmp(iomega-1)+a2F_tmp(iomega) ) * (dw/2d0)
            enddo
            !
            if((E1.ne.0d0).and.(E2.ne.0d0))Kph_e(ie1,ie2) = (2d0/(tanh(Beta/2d0*E1)*tanh(Beta/2d0*E2))) * a2F_int / DoS0_DFT
            !
         enddo !iE2
      enddo !iE1
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(a2F_tmp)
      !
      !Filling the Fermi lines
      call interpFermi(Kph_e,Egrid,Egrid_DFT,Efermi_ndx1,Efermi_ndx2)
      !
      write(*,"(2(A,1E20.10))")"     lambda(Kph): ",-Kph_e(Efermi_ndx1,Efermi_ndx2)*DoS0_DFT,"    DoS_DFT(E=0): ",DoS0_DFT
      !
      if(present(printKpath).and.(reg(printmode).ne."None"))then
         !
         Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
         printmode_used="E0"
         if(present(printmode))printmode_used=reg(printmode)
         call print_Kernel("phononic",printmode_used,reg(printKpath),"Kph",Temp,Egrid,Egrid_DFT,dcmplx(Kph_e,0d0))
         !
      endif
      !
   end subroutine calc_Kph_e



   !========================== AUXILIARY FUNCTIONS ============================!



   !---------------------------------------------------------------------------!
   !PURPOSE: printing Kernel wrapper
   !---------------------------------------------------------------------------!
   subroutine print_Kernel(Kerneltype,printmode,printpath,filename,T,Egrid1,Egrid2,Kernel)
      !
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: Kerneltype
      character(len=*),intent(inout)        :: printmode
      character(len=*),intent(in)           :: printpath
      character(len=*),intent(in)           :: filename
      real(8),intent(in)                    :: T
      real(8),intent(in)                    :: Egrid1(:),Egrid2(:)
      complex(8),intent(in)                 :: Kernel(:,:)
      !
      integer                               :: Ngrid1,Ngrid2,iE,iE1,iE2
      integer                               :: Efermi_ndx1,Efermi_ndx2,unit
      logical                               :: ReK,ImK
      !
      if((reg(Kerneltype).ne."electronic").and.(reg(Kerneltype).ne."phononic"))then
         stop "print_Kernel: available Kerneltype are only electronic or phononic."
      endif
      if((reg(Kerneltype).eq."phononic").and.(reg(printmode).eq."0E"))then
         write(*,"(A)") "     Available print modes of phononic Kernel doen not include 0E, ssetting to E0."
         printmode="E0"
      endif
      !
      Ngrid1 = size(Egrid1)
      Ngrid2 = size(Egrid2)
      call assert_shape(Kernel,[Ngrid1,Ngrid2],"print_Kernel","Kernel")
      Efermi_ndx1 = minloc(abs(Egrid1),dim=1)
      Efermi_ndx2 = minloc(abs(Egrid2),dim=1)
      !
      ImK=.false.
      if(reg(Kerneltype).eq."phononic")ImK=.true.
      ReK = .not.ImK
      !
      write(*,"(A)") "     Printing "//reg(Kerneltype)//" Kernel with mode "//reg(printmode)//" in "//reg(printpath)
      !
      select case(reg(printmode))
         case default
            !
            if(reg(Kerneltype).eq."phononic")stop "Available print modes for phononic Kernel: E0, surf, all."
            if(reg(Kerneltype).eq."electronic")stop "Available print modes for electronic Kernel: E0, 0E, surf, all."
            !
         case("E0")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_e_E0_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid1
               if(ReK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2))
               if(ImK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2)),dimag(Kernel(iE,Efermi_ndx2))
            enddo
            close(unit)
            !
         case("0E")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_e_0E_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid2
               if(ReK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE))
               if(ImK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE)),dimag(Kernel(Efermi_ndx1,iE))
            enddo
            close(unit)
            !
         case("surf")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_e_surf_R_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid1
               do iE2=1,Ngrid2
                  write(unit,"(3F20.10)")Egrid(iE1),Egrid(iE2),dreal(Kernel(iE1,iE2))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(ImK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_e_surf_I_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid1
                  do iE2=1,Ngrid2
                     write(unit,"(3F20.10)")Egrid(iE1),Egrid(iE2),dimag(Kernel(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
            endif
            !
         case("all")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_e_E0_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid1
               if(ReK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2))
               if(ImK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2)),dimag(Kernel(iE,Efermi_ndx2))
            enddo
            close(unit)
            if(reg(Kerneltype).eq."electronic")then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_e_0E_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid2
                  if(ReK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE))
                  if(ImK) write(unit,"(3F20.10)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE)),dimag(Kernel(Efermi_ndx1,iE))
               enddo
               close(unit)
            endif
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_e_surf_R_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid1
               do iE2=1,Ngrid2
                  write(unit,"(3F20.10)")Egrid(iE1),Egrid(iE2),dreal(Kernel(iE1,iE2))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(ImK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_e_surf_I_T"//str(T,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid1
                  do iE2=1,Ngrid2
                     write(unit,"(3F20.10)")Egrid(iE1),Egrid(iE2),dimag(Kernel(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
            endif
            !
            !
      end select
      !
   end subroutine print_Kernel


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
      term1 = term1 -Beta*fermidirac(E1-omega,Beta) * fermidirac(-E1+omega,Beta) * E1 / ( E2 + omega )
      term1 = -term1 * ( fermidirac(E1,Beta) + boseeinstein(omega,Beta) ) / ( E1 - E2 - omega ) * psmooth(E2+omega,Beta)
      !
      term2 = ( fermidirac(E2,Beta) - fermidirac(E1+omega,Beta) ) / ( E1-E2+omega )
      term2 = term2 -Beta*fermidirac(E1+omega,Beta) * fermidirac(-E1-omega,Beta) * E1 / ( E2 - omega )
      term2 = -term2 * ( fermidirac(E1,Beta) + boseeinstein(-omega,Beta) ) / ( E1 - E2 + omega ) * psmooth(E2-omega,Beta)
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
         term1 = 1d0/((1d0+exp(-bE1))*(1d0+exp(-bE2))*(1d0-exp(-bo)))*(exp(-bE2-bo)-exp(-bE1))/(E1-E2-omega)
         !
      elseif((E1.ge.0d0).and.((E2+omega).lt.0d0))then
         !
         term1 = 1d0/(1d0+exp(-bE1))*fermidirac(E2,Beta)*boseeinstein(omega,Beta)*(1d0-exp(bE2+bo-bE1))/(E1-E2-omega)
         !
      elseif((E1.lt.0d0).and.((E2+omega).ge.0d0))then
         !
         term1 = fermidirac(E1,Beta)/((1d0+exp(-bE2))*(1d0-exp(-bo)))*(exp(bE1-bE2-bo)-1d0)/(E1-E2-omega)
         !
      elseif((E1.lt.0d0).and.((E2+omega).lt.0d0))then
         !
         term1 = fermidirac(E1,Beta)*fermidirac(E2,Beta)*boseeinstein(omega,Beta)*(exp(bE1)-exp(bE2+bo))/(E1-E2-omega)
         !
      else
         write(*,"(4(A,F))") "E1",E1,"E2",E2,"omega",omega,"Beta",Beta
         stop "Error in I(E1,E2,omega) - first term."
      endif
      !
      !term 2
      if((E2.ge.0d0).and.((E1+omega).ge.0d0))then
         !
         term2 = 1d0/((1d0+exp(-bE1))*(1d0+exp(-bE2))*(1d0-exp(-bo)))*(exp(-bE1-bo)-exp(-bE2))/(E1-E2+omega)
         !
      elseif((E2.ge.0d0).and.((E1+omega).lt.0d0)) then
         !
         term2 = fermidirac(E1,Beta)/(1d0+exp(-bE2))*boseeinstein(omega,Beta)*(1d0-exp(bE1+bo-bE2))/(E1-E2+omega)
         !
      elseif((E2.lt.0d0).and.((E1+omega).ge.0d0)) then
         !
         term2 = 1d0/(1d0+exp(-bE1))*fermidirac(E2,Beta)/(1d0-exp(-bo))*(exp(bE2-bE1-bo)-1d0)/(E1-E2+omega)
         !
      elseif((E2.lt.0d0).and.((E1+omega).lt.0d0)) then
         !
         term2 = fermidirac(E1,Beta)*fermidirac(E2,Beta)*boseeinstein(omega,Beta)*(exp(bE2)-exp(bE1+bo))/(E1-E2+omega)
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
      term1 = fermidirac(E2,Beta)  * boseeinstein(omega,Beta) / (E1-E2-omega) *  &
      ( -Beta*fermidirac(-E1,Beta) * fermidirac(-E1,Beta) + Beta*fermidirac(-E1,Beta) - fermidirac(-E1,Beta)/(E1-E2-omega) )
      !
      term2 = fermidirac(-E2,Beta) * boseeinstein(-omega,Beta) / (E1-E2-omega) * &
      ( +Beta*fermidirac(E1,Beta)  * fermidirac(E1,Beta)  - Beta*fermidirac(E1,Beta)  - fermidirac(E1,Beta)/(E1-E2-omega)  )
      !
      term3 = fermidirac(-E2,Beta) * boseeinstein(omega,Beta) / (E1-E2+omega) *  &
      ( +Beta*fermidirac(E1,Beta)  * fermidirac(E1,Beta)  - Beta*fermidirac(E1,Beta)  - fermidirac(E1,Beta)/(E1-E2+omega)  )
      !
      term4 = fermidirac(E2,Beta)  * boseeinstein(-omega,Beta) / (E1-E2+omega) * &
      ( -Beta*fermidirac(-E1,Beta) * fermidirac(-E1,Beta) + Beta*fermidirac(-E1,Beta) - fermidirac(-E1,Beta)/(E1-E2+omega) )
      !
      Iprime = term1 + term2 - term3 - term4
      !
   end function Iprime


   !---------------------------------------------------------------------------!
   !PURPOSE: linear fit of the energy dependent Kernel along Fermi levels lines
   !---------------------------------------------------------------------------!
   subroutine interpFermi_Z_d(Z,E,iEo,shift)
      use utils_misc
      implicit none
      real(8),intent(inout)           :: Z(:)
      real(8),intent(in)              :: E(:)
      integer,intent(in)              :: iEo
      integer,intent(in),optional     :: shift
      !
      integer                         :: iE,NE,shift_
      real(8)                         :: yA,yB,xA,xB
      real(8)                         :: y_N,y_S
      !
      !
      if(verbose)write(*,"(A)") "---- interpFermi_Z_d"
      !
      !
      NE = size(E)
      if(E(iEo).ne.0d0) stop "interpFermi_Z_d: provided index does not correspond to Fermi level."
      call assert_shape(Z,[NE],"interpFermi_Z_d","Z")
      !
      shift_=0
      if(present(shift))shift_=shift
      do iE=1,NE
         !
         !linear fit of the (E,0) line - south
         yA = Z(iEo-1-shift_) ; xA = Egrid(iEo-1-shift_)
         yB = Z(iEo-2-shift_) ; xB = Egrid(iEo-2-shift_)
         y_S = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !linear fit of the (E,0) line - north
         yA = Z(iEo+1+shift_) ; xA = Egrid(iEo+1+shift_)
         yB = Z(iEo+2+shift_) ; xB = Egrid(iEo+2+shift_)
         y_N = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         Z(iEo) = ( y_S + y_N )/2d0
         !
      enddo
      !
   end subroutine interpFermi_Z_d
   !
   subroutine interpFermi_Z_z(Z,E,iEo,shift)
      use utils_misc
      implicit none
      complex(8),intent(inout)        :: Z(:)
      real(8),intent(in)              :: E(:)
      integer,intent(in)              :: iEo
      integer,intent(in),optional     :: shift
      !
      integer                         :: iE,NE,shift_
      complex(8)                      :: yA,yB
      real(8)                         :: xA,xB
      real(8)                         :: Rey_N,Rey_S
      real(8)                         :: Imy_N,Imy_S
      !
      !
      if(verbose)write(*,"(A)") "---- interpFermi_Z_z"
      !
      !
      NE = size(E)
      if(E(iEo).ne.0d0) stop "interpFermi_Z_z: provided index does not correspond to Fermi level."
      call assert_shape(Z,[NE],"interpFermi_Z_z","Z")
      !
      shift_=0
      if(present(shift))shift_=shift
      do iE=1,NE
         !
         !linear fit of the (E,0) line - south
         yA = Z(iEo-1-shift_) ; xA = Egrid(iEo-1-shift_)
         yB = Z(iEo-2-shift_) ; xB = Egrid(iEo-2-shift_)
         Rey_S = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_S = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !linear fit of the (E,0) line - north
         yA = Z(iEo+1+shift_) ; xA = Egrid(iEo+1+shift_)
         yB = Z(iEo+2+shift_) ; xB = Egrid(iEo+2+shift_)
         Rey_N = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_N = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         Z(iEo) = dcmplx(( Rey_S + Rey_N ),( Imy_S + Imy_N ))/2d0
         !
      enddo
      !
   end subroutine interpFermi_Z_z
   !
   subroutine interpFermi_K_d(K,E1,E2,iE1o,iE2o,shift)
      use utils_misc
      implicit none
      real(8),intent(inout)           :: K(:,:)
      real(8),intent(in)              :: E1(:),E2(:)
      integer,intent(in)              :: iE1o,iE2o
      integer,intent(in),optional     :: shift
      !
      integer                         :: iE,NE,shift_
      real(8)                         :: yA,yB,xA,xB
      real(8)                         :: y_N,y_S,y_E,y_W
      !
      !
      if(verbose)write(*,"(A)") "---- interpFermi_K_d"
      !
      !
      NE1 = size(E1)
      NE2 = size(E2)
      if(abs(E1(iE1o)).gt.1d-9) stop "interpFermi_K_d: provided index iE1o does not correspond to Fermi level."
      if(abs(E2(iE2o)).gt.1d-9) stop "interpFermi_K_d: provided index iE2o does not correspond to Fermi level."
      call assert_shape(K,[NE1,NE2],"interpFermi_K_d","K")
      !
      shift_=0
      if(present(shift))shift_=shift
      !
      !linear fit of the (E1,0)=(iE,iE2o) column
      do iE=1,NE1
         !
         if(iE.eq.iE1o)cycle
         !
         !at each row iE - east direction
         xA = E2(iE2o+1+shift_); yA = K(iE,iE2o+1+shift_)
         xB = E2(iE2o+2+shift_); yB = K(iE,iE2o+2+shift_)
         y_E = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !at each row iE - west direction
         xA = E2(iE2o-1-shift_); yA = K(iE,iE2o-1-shift_)
         xB = E2(iE2o-2-shift_); yB = K(iE,iE2o-2-shift_)
         y_W = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !the iE,iE2o column is the east/west average
         K(iE,iE2o) = ( y_E + y_W )/2d0
         !
      enddo
      !
      !linear fit of the (0,E2)=(iE1o,iE) row
      do iE=1,NE2
         !
         if(iE.eq.iE2o)cycle
         !
         !at each column iE - north direction
         xA = E1(iE1o+1+shift_); yA = K(iE1o+1+shift_,iE)
         xB = E1(iE1o+2+shift_); yB = K(iE1o+2+shift_,iE)
         y_N = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !at each column iE - south direction
         xA = E1(iE1o-1-shift_); yA = K(iE1o-1-shift_,iE)
         xB = E1(iE1o-2-shift_); yB = K(iE1o-2-shift_,iE)
         y_S = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !the iE1o,iE row is the north/south average
         K(iE1o,iE) = ( y_N + y_S )/2d0
         !
      enddo
      !
      !linear fit around the (0,0)=(iE1o,iE2o) point
      !east direction
      xA = E2(iE2o+1+shift_); yA = K(iE1o,iE2o+1+shift_)
      xB = E2(iE2o+2+shift_); yB = K(iE1o,iE2o+2+shift_)
      y_E = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !west direction
      xA = E2(iE2o-1-shift_); yA = K(iE1o,iE2o-1-shift_)
      xB = E2(iE2o-2-shift_); yB = K(iE1o,iE2o-2-shift_)
      y_W = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !north direction
      xA = E1(iE1o+1+shift_); yA = K(iE1o+1+shift_,iE2o)
      xB = E1(iE1o+2+shift_); yB = K(iE1o+2+shift_,iE2o)
      y_N = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !south direction
      xA = E1(iE1o-1-shift_); yA = K(iE1o-1-shift_,iE2o)
      xB = E1(iE1o-2-shift_); yB = K(iE1o-2-shift_,iE2o)
      y_S = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !
      !the iE1o,iE2o point is the east/west/north/south average
      K(iE1o,iE2o) = ( y_S + y_N + y_W + y_E )/4d0
      !
   end subroutine interpFermi_K_d
   !
   subroutine interpFermi_K_z(K,E1,E2,iE1o,iE2o,shift)
      use utils_misc
      implicit none
      complex(8),intent(inout)        :: K(:,:)
      real(8),intent(in)              :: E1(:),E2(:)
      integer,intent(in)              :: iE1o,iE2o
      integer,intent(in),optional     :: shift
      !
      integer                         :: iE,NE,shift_
      complex(8)                      :: yA,yB
      real(8)                         :: xA,xB
      real(8)                         :: Rey_N,Rey_S,Rey_E,Rey_W
      real(8)                         :: Imy_N,Imy_S,Imy_E,Imy_W
      !
      !
      if(verbose)write(*,"(A)") "---- interpFermi_K_z"
      !
      !
      NE1 = size(E1)
      NE2 = size(E2)
      if(abs(E1(iE1o)).gt.1d-9) stop "interpFermi_K_z: provided index iE1o does not correspond to Fermi level."
      if(abs(E2(iE2o)).gt.1d-9) stop "interpFermi_K_z: provided index iE2o does not correspond to Fermi level."
      call assert_shape(K,[NE1,NE2],"interpFermi_K_z","K")
      !
      shift_=0
      if(present(shift))shift_=shift
      !
      !linear fit of the (E1,0)=(iE,iE2o) column
      do iE=1,NE1
         !
         if(iE.eq.iE1o)cycle
         !
         !at each row iE - east direction
         xA = E2(iE2o+1+shift_); yA = K(iE,iE2o+1+shift_)
         xB = E2(iE2o+2+shift_); yB = K(iE,iE2o+2+shift_)
         Rey_E = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_E = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !at each row iE - west direction
         xA = E2(iE2o-1-shift_); yA = K(iE,iE2o-1-shift_)
         xB = E2(iE2o-2-shift_); yB = K(iE,iE2o-2-shift_)
         Rey_W = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_W = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !the iE,iE2o column is the east/west average
         K(iE,iE2o) = dcmplx(( Rey_E + Rey_W ),( Imy_E + Imy_W ))/2d0
         !
      enddo
      !
      !linear fit of the (0,E2)=(iE1o,iE) row
      do iE=1,NE2
         !
         if(iE.eq.iE2o)cycle
         !
         !at each column iE - north direction
         xA = E1(iE1o+1+shift_); yA = K(iE1o+1+shift_,iE)
         xB = E1(iE1o+2+shift_); yB = K(iE1o+2+shift_,iE)
         Rey_N = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_N = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !at each column iE - south direction
         xA = E1(iE1o-1-shift_); yA = K(iE1o-1-shift_,iE)
         xB = E1(iE1o-2-shift_); yB = K(iE1o-2-shift_,iE)
         Rey_S = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_S = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !the iE1o,iE row is the north/south average
         K(iE1o,iE) = dcmplx(( Rey_N + Rey_S ),( Imy_N + Imy_S ))/2d0
         !
      enddo
      !
      !linear fit around the (0,0)=(iE1o,iE2o) point
      !east direction
      xA = E2(iE2o+1+shift_); yA = K(iE1o,iE2o+1+shift_)
      xB = E2(iE2o+2+shift_); yB = K(iE1o,iE2o+2+shift_)
      Rey_E = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_E = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !west direction
      xA = E2(iE2o-1-shift_); yA = K(iE1o,iE2o-1-shift_)
      xB = E2(iE2o-2-shift_); yB = K(iE1o,iE2o-2-shift_)
      Rey_W = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_W = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !north direction
      xA = E1(iE1o+1+shift_); yA = K(iE1o+1+shift_,iE2o)
      xB = E1(iE1o+2+shift_); yB = K(iE1o+2+shift_,iE2o)
      Rey_N = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_N = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !south direction
      xA = E1(iE1o-1-shift_); yA = K(iE1o-1-shift_,iE2o)
      xB = E1(iE1o-2-shift_); yB = K(iE1o-2-shift_,iE2o)
      Rey_S = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_S = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !
      !the iE1o,iE2o point is the east/west/north/south average
      K(iE1o,iE2o) = dcmplx(( Rey_S + Rey_N + Rey_W + Rey_E ),( Imy_S + Imy_N + Imy_W + Imy_E ))/4d0
      !
   end subroutine interpFermi_K_z


end module gap_equation
