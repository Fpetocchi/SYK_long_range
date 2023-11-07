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

   interface get_Kel
      module procedure get_Kel_integral
      module procedure get_Kel_list
   end interface get_Kel

   interface io_Kel
      module procedure io_Kel_d2
      module procedure io_Kel_d3
   end interface io_Kel

   interface calc_energy_averages
      module procedure calc_energy_averages_integral
      module procedure calc_energy_averages_list
   end interface calc_energy_averages

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),allocatable,private              :: omega(:)                         !Phonon energy on logarithmic grid
   real(8),allocatable,private              :: a2F(:)                           !alpha^2*F(\Omega) function
   !
   integer,private                          :: Nkpt3_Model(3)
   real(8),allocatable,private              :: kpt_Model(:,:)
   real(8),allocatable,private              :: weights_Model(:,:,:)
   integer,allocatable,private              :: finite_weights_Model(:,:)
   !
   real(8),allocatable,private              :: Ek_Model(:,:)
   complex(8),allocatable,private           :: Zk_Model(:,:,:)
   complex(8),allocatable,private           :: Kel_stat(:,:)
   complex(8),allocatable,private           :: Kel_dyn_list(:,:,:)
   complex(8),allocatable,private           :: Wee_dyn(:,:,:)
   !
   logical,private                          :: initialized=.false.
   logical,private                          :: Phonons_stored=.false.
   logical,private                          :: DFT_DoS_stored=.false.
   logical,private                          :: Interpolate2Model=.false.
   logical,private                          :: Kernels_stored=.false.
   logical,private                          :: calc_Int_static=.false.
   logical,private                          :: calc_Int_dynamic=.false.
   !
   integer,private                          :: wmax_ndx                         !input screened interaction cutoff index
   real(8),private                          :: wmax                             !input screened interaction cutoff
   real(8),private                          :: MatsStep                         !input screened interaction frequency step
   real(8),allocatable,private              :: wmats_orig(:)                    !input Matsubara gird in DFT units
   !
   !public
   character(len=2),public,protected        :: DFTgrid="eV"
   real(8),public,protected                 :: eV2DFTgrid=1d0                   !this converts eV to the DFT grid
   real(8),public,protected                 :: DFTgrid2eV=1d0
   !
   real(8),allocatable,public,protected     :: Egrid(:)                         !this is the dft one if phonons are present, custom otherwise
   real(8),allocatable,public,protected     :: DoS_DFT(:)
   real(8),allocatable,public,protected     :: DoS_Model(:)
   !
   logical,public,protected                 :: calc_phonons=.false.
   logical,public,protected                 :: calc_Kel=.false.
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
   public :: overload_G0W0
   public :: overload_DMFT
   public :: calc_energy_averages
   public :: get_Kel
   public :: calc_Zph_e
   public :: calc_Kph_e

   !===========================================================================!

contains


   !======================== SETTING UP GLOBAL VARS ===========================!


   !---------------------------------------------------------------------------!
   !PURPOSE: Read phonons and initialize only mesh and energy grids
   !---------------------------------------------------------------------------!
   subroutine Initialize_inputs(pathINPUT,Inputs,Lttc,Hk_input)
      !
      use parameters
      use utils_misc
      use crystal, only : calc_irredBZ, wannierinterpolation, tetrahedron_integration
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      type(SCDFT),intent(in)                :: Inputs
      type(Lattice),intent(in)              :: Lttc
      complex(8),intent(in)                 :: Hk_input(:,:,:)
      !
      integer                               :: Norb,Ngrid,Nkpt,Nkpti_dum
      integer                               :: iorb,iE,ik,Nweights
      real(8)                               :: wrealMax
      integer,allocatable                   :: kptp_dum(:)
      integer,allocatable                   :: pkpt_dum(:,:,:)
      real(8),allocatable                   :: nkstar_dum(:)
      complex(8),allocatable                :: Hk_intp(:,:,:)
      logical                               :: keep_weight
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- Initialize_inputs"
      !
      !
      if(.not.Inputs%status)stop "Initialize_inputs: input container not properly initialized."
      !
      wrealMax = Inputs%wrealMax
      Ngrid = Inputs%Ngrid
      if(mod(Ngrid,2).eq.0)Ngrid=Ngrid+1
      if(mod(Ngrid-1,4).ne.0)Ngrid=Ngrid+mod(Ngrid-1,4)
      allocate(Egrid(Ngrid));Egrid=0d0
      Egrid=denspace(wrealMax,Ngrid,center=.true.,expfact=Inputs%expfact)
      if(Egrid(minloc(abs(Egrid),dim=1)).ne.0d0) stop "Initialize_inputs: the energy grid requires the E=0 point."
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
            !
         case("Elk","QEspresso")
            !
            if(reg(Inputs%mode_ph).eq."Elk")then
               DFTgrid="Hr"
               eV2DFTgrid = eV2H
               DFTgrid2eV = H2eV
            elseif(reg(Inputs%mode_ph).eq."QEspresso")then
               DFTgrid="Ry"
               eV2DFTgrid = eV2Ry
               DFTgrid2eV = Ry2eV
            endif
            !
            Egrid = Egrid*eV2DFTgrid
            !
            calc_phonons=.true.
            call read_a2F_DFT(reg(pathINPUT),reg(Inputs%mode_ph))
            call read_DoS_DFT(reg(pathINPUT))
            !
      end select
      !
      select case(reg(Inputs%mode_el))
         case default
            !
            stop "Available electronic modes in gap equation: static, dynamic, static+dynamic, None."
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
            calc_Int_static=.true.
            calc_Int_dynamic=.true.
            !
      end select
      calc_Kel = calc_Int_static .or. calc_Int_dynamic
      write(*,"(A)")"     eV->DFTgrid conversion factor: "//str(eV2DFTgrid,5)
      !
      !setting up shared original K-mesh and Hk
      Norb = size(Hk_input,dim=1)
      !
      !compute the interpolated K-grid for Hk and the corresponding weights
      write(*,"(A)")new_line("A")//"     Weights calculation for electronic energy integrals."
      !
      Interpolate2Model = (.not.(any(Inputs%Nkpt3_Model.eq.0))) .and. (.not.(all(Inputs%Nkpt3_Model.eq.Lttc%Nkpt3)))
      if(Interpolate2Model)then
         !
         Nkpt3_Model = Inputs%Nkpt3_Model
         Nkpt = product(Nkpt3_Model)
         allocate(Hk_intp(Norb,Norb,Nkpt));Hk_intp=czero
         allocate(kpt_Model(3,Nkpt));kpt_Model=0d0
         !
         call calc_irredBZ(reg(pathINPUT),Nkpt3_Model,Nkpti_dum,kptp_dum,pkpt_dum,nkstar_dum,kpt_out=kpt_Model)
         deallocate(kptp_dum,pkpt_dum,nkstar_dum)
         !
         call cpu_time(start)
         call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,kpt_Model,(Hk_input*eV2DFTgrid),Hk_intp)
         call cpu_time(finish)
         write(*,"(A,F)") "     Interpolation to new ["//str(Nkpt3_Model(1))//","//str(Nkpt3_Model(2))//","//str(Nkpt3_Model(3))//"] K-grid cpu timing:", finish-start
         !
      else
         !
         Nkpt3_Model = Lttc%Nkpt3
         Nkpt = product(Nkpt3_Model)
         allocate(Hk_intp(Norb,Norb,Nkpt));Hk_intp=czero
         allocate(kpt_Model(3,Nkpt));kpt_Model=0d0
         !
         kpt_Model = Lttc%kpt
         Hk_intp = Hk_input * eV2DFTgrid
         write(*,"(A)")"     Interpolation skipped."
         !
      endif
      !
      allocate(Zk_Model(Norb,Norb,Nkpt));Zk_Model=czero
      allocate(Ek_Model(Norb,Nkpt));Ek_Model=0d0
      call calc_rotation(Hk_intp,Zk_Model,Ek_Model)
      !
      allocate(weights_Model(Ngrid,Norb,Nkpt));weights_Model=0d0
      allocate(DoS_Model(Ngrid));DoS_Model=0d0
      !
      call cpu_time(start)
      call tetrahedron_integration(reg(pathINPUT),Hk_intp,Nkpt3_Model,kpt_Model,Egrid,weights_out=weights_Model,DoS_out=DoS_Model)!,pathOUTPUT=reg(pathINPUT))
      call cpu_time(finish)
      write(*,"(A,F)") "     Tetrahedron integration cpu timing:", finish-start
      deallocate(Hk_intp)
      !
      call cpu_time(start)
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(Inputs%DoSthresh) )
               if(keep_weight) Nweights = Nweights + 1
            enddo
         enddo
      enddo
      !
      allocate(finite_weights_Model(Nweights,3));finite_weights_Model=0
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(Inputs%DoSthresh) )
               if(keep_weight)then
                  Nweights = Nweights + 1
                  finite_weights_Model(Nweights,1) = iE
                  finite_weights_Model(Nweights,2) = iorb
                  finite_weights_Model(Nweights,3) = ik
               endif
            enddo
         enddo
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Reduction of DoS integration points from "//str(Ngrid*Norb*Nkpt)//" to "//str(Nweights)
      write(*,"(A,F)") "     Cpu timing:", finish-start
      !
      initialized=.true.
      !
      !
      !
      contains
      !
      !
      !
      subroutine calc_rotation(H,Z,E)
         !
         use linalg, only : eigh
         use utils_misc, only : assert_shape, F2Bindex
         implicit none
         complex(8),intent(in)              :: H(:,:,:)
         complex(8),intent(out)             :: Z(:,:,:)
         real(8),intent(out)                :: E(:,:)
         !
         integer                            :: i,Ndim,Nk
         complex(8),allocatable             :: Z_(:,:,:)
         !
         Ndim = size(H,dim=1)
         Nk = size(H,dim=3)
         call assert_shape(H,[Ndim,Ndim,Nk],"calc_rotation","H")
         call assert_shape(Z,[Ndim,Ndim,Nk],"calc_rotation","Z")
         call assert_shape(E,[Ndim,Nk],"calc_rotation","E")
         !
         allocate(Z_(Ndim,Ndim,Nk)); Z_=0d0
         !
         Z_ = H
         do i=1,Nk
            call eigh(Z_(:,:,i),E(:,i))
         enddo
         !
         Z = Z_
         deallocate(Z_)
         !
      end subroutine calc_rotation
      !
      !
      !
   end subroutine Initialize_inputs


   !---------------------------------------------------------------------------!
   !PURPOSE: Overload the DoS and rotation with the results computed via the
   !         G0W0 self-energy. Works only if all the required data are present
   !---------------------------------------------------------------------------!
   subroutine overload_G0W0(pathINPUT,pathINPUTtr,Lttc,DoSthresh)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use crystal
      use file_io
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=*),intent(in)           :: pathINPUTtr
      type(Lattice),intent(in)              :: Lttc
      real(8),intent(in)                    :: DoSthresh
      !
      type(FermionicField)                  :: S_G0W0
      complex(8),allocatable                :: Smat(:,:,:,:),SmatE(:,:,:,:)
      complex(8),allocatable                :: Uwan(:,:,:,:),Vxc(:,:,:,:)
      complex(8),allocatable                :: Zk(:,:,:)
      real(8),allocatable                   :: wreal_read(:)
      real(8)                               :: ReS,ImS,mu,eta,Norm
      integer                               :: ik,iw,iE,iorb,ispin
      integer                               :: Norb,Nkpt,Ngrid,Nweights
      integer                               :: Nreal_read = 5000
      logical                               :: paramagnet=.true.
      logical                               :: keep_weight,UWAN_exist,Vxc_exist,G0W0_exist
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- overload_G0W0"
      !
      !
      if(.not.initialized) stop "overload_G0W0: gap equation module not properly initialized. Call Initialize_inputs first."
      !
      call inquireFile(reg(pathINPUT)//"UWAN_used_k_s1.DAT",UWAN_exist,hardstop=.false.)
      call inquireFile(reg(pathINPUT)//"Vxc_k_s1.DAT",Vxc_exist,hardstop=.false.)
      call inquireFile(reg(pathINPUTtr)//"SGoWo_wr_k_s1.DAT",G0W0_exist,hardstop=.false.)
      !
      if((.not.UWAN_exist).or.(.not.Vxc_exist).or.(.not.G0W0_exist).or.Interpolate2Model)then
         if(.not.UWAN_exist)   write(*,"(A)")"     UWAN_used_k_s1.DAT not found. Ignoring call to subroutine."
         if(.not.Vxc_exist)    write(*,"(A)")"     Vxc_k_s1.DAT not found. Ignoring call to subroutine."
         if(.not.G0W0_exist)   write(*,"(A)")"     SGoWo_wr_k_s1.DAT not found. Ignoring call to subroutine."
         if(Interpolate2Model) write(*,"(A)")"     Cannot read G0W0 quantities in another k-grid. Ignoring call to subroutine."
         return
      endif
      !
      mu = Lttc%mu
      write(*,"(A,F)")"     Using mu=",mu
      !
      Norb = size(Zk_Model,dim=1)
      Nkpt = size(Zk_Model,dim=3)
      Ngrid = size(Egrid)
      !
      !read rotation matrix
      allocate(Uwan(Norb,Norb,Nkpt,1));Uwan=czero
      call read_matrix(Uwan(:,:,:,1),reg(pathINPUT)//"UWAN_used_k_s1.DAT")
      allocate(Zk(Norb,Norb,Nkpt));Zk=czero
      do ik=1,Nkpt
         Zk(:,:,ik) = dag(Uwan(:,:,ik,1))
      enddo
      deallocate(Uwan)
      write(*,"(A)")"     Rotation matrix is stored."
      !
      !read Vxc
      allocate(Vxc(Norb,Norb,Nkpt,1));Uwan=czero
      call read_matrix(Vxc(:,:,:,1),reg(pathINPUT)//"Vxc_k_s1.DAT")
      write(*,"(A)")"     Vxc is stored."
      !
      !read G0W0 self-energy
      allocate(wreal_read(Nreal_read));wreal_read=0d0
      call AllocateFermionicField(S_G0W0,Norb,Nreal_read,Nkpt=Nkpt)
      call read_FermionicField(S_G0W0,reg(pathINPUTtr),"SGoWo_wr",Lttc%kpt,axis=wreal_read)
      do iw=Nreal_read,1,-1
         if(abs(wreal_read(iw)).ne.0)exit
      enddo
      Nreal_read = iw-1
      wreal_read = wreal_read(1:Nreal_read)
      write(*,"(A)") "     Input G0W0 self-energy."
      write(*,"(A,1I6)") "     Frequency axis updated to: ",Nreal_read
      write(*,"(A)") "     Frequency boundaries "//str(wreal_read(1),4)//" - "//str(wreal_read(Nreal_read),4)
      !
      !rotate G0W0-Vxc input to LDA basis
      allocate(Smat(Norb,Nreal_read,Nkpt,Nspin));Smat=czero
      do ispin=1,Nspin
         !
         do ik=1,Nkpt
            do iw=1,Nreal_read
               Smat(:,iw,ik,ispin) = diagonal(rotate((S_G0W0%wks(:,:,iw,ik,ispin)-Vxc(:,:,ik,ispin)),Zk(:,:,ik)))
            enddo
         enddo
         if(paramagnet)then
            Smat(:,:,:,Nspin) = Smat(:,:,:,1)
            exit
         endif
         !
      enddo
      call DeallocateField(S_G0W0)
      deallocate(Zk)
      !
      !rescale to DFT units
      wreal_read = wreal_read * eV2DFTgrid
      Smat = Smat * eV2DFTgrid
      !
      !interpolate to logarithmic energy grid
      call cpu_time(start)
      allocate(SmatE(Norb,Ngrid,Nkpt,Nspin));SmatE=czero
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE,ik,iorb,ispin,ReS,ImS)
      !$OMP DO
      do iE=1,Ngrid
         do ik=1,Nkpt
            do ispin=1,Nspin
               !
               do iorb=1,Norb
                  ReS = cubic_interp( wreal_read, dreal(Smat(iorb,1:Nreal_read,ik,ispin)), Egrid(iE) )
                  ImS = cubic_interp( wreal_read, dimag(Smat(iorb,1:Nreal_read,ik,ispin)), Egrid(iE) )
                  if(ImS.gt.0d0)ImS=0d0
                  SmatE(iorb,iE,ik,ispin) = dcmplx(ReS,ImS)
               enddo
               if(paramagnet)then
                  SmatE(:,iE,ik,Nspin) = SmatE(:,iE,ik,1)
                  cycle
               endif
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(wreal_read,Smat)
      !
      !the new weights are replaced the G0W0 spectral functions in the LDA basis. For now spinless.
      eta=0d0
      weights_Model=0d0
      DoS_Model=0d0
      do ik=1,Nkpt
         do iorb=1,Norb
            !
            do iE=1,Ngrid
               weights_Model(iE,iorb,ik) = -dimag( 1d0 / ( dcmplx(Egrid(iE)+mu,eta) - Ek_Model(iorb,ik) - SmatE(iorb,iE,ik,1) ) )
            enddo
            !
            !normalize specific weights component
            Norm=0d0
            do iE=2,Ngrid
               Norm = Norm + (weights_Model(iE,iorb,ik)+weights_Model(iE-1,iorb,ik))*(Egrid(iE)-Egrid(iE-1))/2d0
            enddo
            weights_Model(:,iorb,ik) = weights_Model(:,iorb,ik)/Norm
            !
         enddo
      enddo
      deallocate(SmatE)
      !
      !sum to get the DoS (k average, orbital sum) normalized to the particle density
      do iE=1,Ngrid
         do ik=1,Nkpt
            do iorb=1,Norb
               DoS_Model(iE) = DoS_Model(iE) + weights_Model(iE,iorb,ik)/Nkpt
            enddo
         enddo
      enddo
      !
      !reallocate the weight above threshold
      call cpu_time(start)
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(DoSthresh) )
               if(keep_weight) Nweights = Nweights + 1
            enddo
         enddo
      enddo
      !
      if(allocated(finite_weights_Model))deallocate(finite_weights_Model)
      allocate(finite_weights_Model(Nweights,3));finite_weights_Model=0
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(DoSthresh) )
               if(keep_weight)then
                  Nweights = Nweights + 1
                  finite_weights_Model(Nweights,1) = iE
                  finite_weights_Model(Nweights,2) = iorb
                  finite_weights_Model(Nweights,3) = ik
               endif
            enddo
         enddo
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Reduction of DoS integration points from "//str(Ngrid*Norb*Nkpt)//" to "//str(Nweights)
      write(*,"(A,F)") "     Cpu timing:", finish-start
      !
   end subroutine overload_G0W0


   !---------------------------------------------------------------------------!
   !PURPOSE: Overload the DoS and rotation with the results computed via the
   !         G0W0 self-energy. Works only if all the required data are present
   !---------------------------------------------------------------------------!
   subroutine overload_DMFT(path2MaxEnt,Lttc,DoSthresh)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use crystal
      use file_io
      implicit none
      !
      character(len=*),intent(in)           :: path2MaxEnt
      type(Lattice),intent(in)              :: Lttc
      real(8),intent(in)                    :: DoSthresh
      !
      real(8),allocatable                   :: wreal_read(:),ImG_read(:,:,:)
      real(8)                               :: dw
      integer                               :: ik,iw,iE,iorb,unit
      integer                               :: Norb,Nkpt,Ngrid,Nweights
      integer                               :: Nreal_read,Nreal_old,ierr
      logical                               :: filexists,keep_weight
      character(len=256)                    :: path
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- overload_DMFT"
      !
      !
      if(.not.initialized) stop "overload_DMFT: gap equation module not properly initialized. Call Initialize_inputs first."
      if(Interpolate2Model) write(*,"(A)")"     Cannot read DMFT quantities in another k-grid. Ignoring call to subroutine."
      !
      Norb = size(Zk_Model,dim=1)
      Nkpt = size(Zk_Model,dim=3)
      Ngrid = size(Egrid)
      !
      !replace the rotation matrix with what is used to get the diagonal spectra
      do ik=1,Nkpt
         Zk_Model(:,:,ik) = Lttc%Zk(:,:,ik)
      enddo
      write(*,"(A)")"     Rotation matrix is stored."
      !
      do ik=1,Nkpt
         !
         path = reg(path2MaxEnt)//"MaxEnt_Gk_full_s1/Gk_t_k"//str(ik)//".DAT_dos.dat"
         !
         call inquireFile(reg(path),filexists,hardstop=.false.,verb=.true.)
         if(.not.filexists)then
            write(*,"(A,1I5)") "     Some K-points are missing in the MaxEnt folder. Ignoring call to subroutine."
            return
         endif
         !
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
         !
         Nreal_read=0
         ierr=0
         do while (ierr.eq.0)
            Nreal_read = Nreal_read + 1
            read(unit,*,iostat=ierr)
         enddo
         close(unit)
         !
         !MaxEnt parameters written i the last line
         Nreal_read = Nreal_read - 2
         !
         !write(*,"(A,1I5)") "     The file "//reg(path)//" contains "//str(Nreal_read)//" real frequencies."
         if(Nreal_read.ne.Nreal_old)then
            write(*,"(A,1I5)") "     Real frequency mesh is not consistent among K-points. Ignoring call to subroutine."
            return
         endif
         Nreal_old=Nreal_read
         !
      enddo
      !
      allocate(wreal_read(Nreal_read));wreal_read=0d0
      allocate(ImG_read(Norb,Nreal_read,Nkpt));ImG_read=0d0
      do ik=1,Nkpt
         !
         path = reg(path2MaxEnt)//"MaxEnt_Gk_full_s1/Gk_t_k"//str(ik)//".DAT_dos.dat"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
         do iw=1,Nreal_read
            read(unit,*) wreal_read(iw),(ImG_read(iorb,iw,ik),iorb=1,Norb)
         enddo
         close(unit)
         !
         !Fix normalization
         dw = abs(wreal_read(10)-wreal_read(9))
         do iorb=1,Norb
            ImG_read(iorb,:,ik) = ImG_read(iorb,:,ik) / abs(sum(ImG_read(iorb,:,ik))*dw)
         enddo
         !
      enddo
      write(*,"(A)") "     MaxEnt output on Green's function is read and normalized."
      !
      !rescale to DFT units
      wreal_read = wreal_read * eV2DFTgrid
      ImG_read = ImG_read / eV2DFTgrid
      !
      !interpolate to logarithmic energy grid
      weights_Model=0d0
      DoS_Model=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE,ik,iorb)
      !$OMP DO
      do iE=1,Ngrid
         do ik=1,Nkpt
            do iorb=1,Norb
               weights_Model(iE,iorb,ik) = cubic_interp( wreal_read, ImG_read(iorb,:,ik), Egrid(iE) )
            enddo
         enddo
         DoS_Model(iE) = DoS_Model(iE) + weights_Model(iE,iorb,ik)/Nkpt 
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(wreal_read,ImG_read)
      !
      !recalculate the weight above threshold
      call cpu_time(start)
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(DoSthresh) )
               if(keep_weight) Nweights = Nweights + 1
            enddo
         enddo
      enddo
      !
      if(allocated(finite_weights_Model))deallocate(finite_weights_Model)
      allocate(finite_weights_Model(Nweights,3));finite_weights_Model=0
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(DoSthresh) )
               if(keep_weight)then
                  Nweights = Nweights + 1
                  finite_weights_Model(Nweights,1) = iE
                  finite_weights_Model(Nweights,2) = iorb
                  finite_weights_Model(Nweights,3) = ik
               endif
            enddo
         enddo
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Reduction of DoS integration points from "//str(Ngrid*Norb*Nkpt)//" to "//str(Nweights)
      write(*,"(A,F)") "     Cpu timing:", finish-start
      !
   end subroutine overload_DMFT


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the output of phonon calculations with  Elk or QuantumEspresso
   !---------------------------------------------------------------------------!
   subroutine read_a2F_DFT(pathINPUT,mode)
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
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_a2F_DFT"
      !
      !
      select case(reg(mode))
         case default
            !
            stop "Available phonon inputs: Elk, QEspresso."
            !
         case("Elk")
            !
            Header = 0
            Footer = 0
            !
         case("QEspresso")
            !
            Header = 5
            Footer = 1
            !
      end select
      !
      !reading the number of phonon frequecies form file depending on the format
      call inquireFile(reg(pathINPUT)//"a2F_DFT.DAT",filexists)
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"a2F_DFT.DAT",form="formatted",action="read",position="rewind")
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
      open(unit,file=reg(pathINPUT)//"a2F_DFT.DAT",form="formatted",action="read",position="rewind")
      do ih=1,Header
         read(unit,*)
      enddo
      do iomega=1,Nomega
        read(unit,*) omega(iomega),a2F(iomega)
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
   end subroutine read_a2F_DFT
   !
   subroutine read_DoS_DFT(pathINPUT)
      !
      use parameters
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      !
      integer                               :: unit,ierr,Nlines
      integer                               :: iE,Ngrid_read,Ngrid
      real(8)                               :: Efermi,Edum,Rdum,Ddum,DoS0_DFT
      real(8),allocatable                   :: Egrid_read(:),DoS_DFT_read(:)
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
      read(unit,*,iostat=ierr) Efermi
      do while (ierr.eq.0)
         read(unit,*,iostat=ierr)Edum,Ddum,Rdum
         if(ierr.eq.0)Nlines = Nlines + 1
      enddo
      close(unit)
      Ngrid_read = Nlines
      write(*,"(A)") "     The number of points in DFT DoS is: "//str(Ngrid_read)
      !
      allocate(Egrid_read(Ngrid_read));Egrid_read=0d0
      allocate(DoS_DFT_read(Ngrid_read));DoS_DFT_read=0d0
      !
      !reading phonons form file
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"DoS_DFT.DAT",form="formatted",action="read",position="rewind")
      read(unit,*) Efermi
      do iE=1,Ngrid_read
         read(unit,*)Egrid_read(iE),DoS_DFT_read(iE),Rdum
      enddo
      close(unit)
      !
      !we need the dos per spin
      DoS_DFT_read = DoS_DFT_read/2d0
      !
      DoS0_DFT = DoS_DFT_read(minloc(abs(Egrid_read - Efermi),dim=1))
      write(*,"(A,F10.5)") "     Fermi Energy (DFT-eV): ",Efermi
      write(*,"(A,F10.5)") "     DoS at Fermi (DFT-1/eV): ",DoS0_DFT
      !
      Egrid_read = (Egrid_read - Efermi) * eV2DFTgrid
      DoS_DFT_read = DoS_DFT_read / eV2DFTgrid
      DoS0_DFT = DoS_DFT_read(minloc(abs(Egrid_read),dim=1))
      write(*,"(A,F10.5)") "     Fermi Energy (DFT-"//DFTgrid//"): ",Egrid_read(minloc(abs(Egrid_read),dim=1))
      write(*,"(A,F10.5)") "     DoS at Fermi (DFT-1/"//DFTgrid//"): ",DoS0_DFT
      !
      !interpolation to the custom grid
      Ngrid = size(Egrid)
      allocate(DoS_DFT(Ngrid));DoS_DFT=0d0
      do iE=1,Ngrid
         DoS_DFT(iE) = cubic_interp( Egrid_read, DoS_DFT_read, Egrid(iE) )
      enddo
      where(DoS_DFT.lt.0d0)DoS_DFT=0d0
      deallocate(Egrid_read,DoS_DFT_read)
      !
      write(*,"(A,F10.5)") "     Fermi Energy (interp-"//DFTgrid//"): ",Egrid(minloc(abs(Egrid),dim=1))
      write(*,"(A,F10.5)") "     DoS at Fermi (interp-1/"//DFTgrid//"): ",DoS_DFT(minloc(abs(Egrid),dim=1))
      DoS_DFT(minloc(abs(Egrid),dim=1)) = DoS0_DFT
      !
      DFT_DoS_stored=.true.
      !
   end subroutine read_DoS_DFT


   !---------------------------------------------------------------------------!
   !PURPOSE: Store the interpolated and rotated NaNb components of the fully
   !         screened interaction. This is the only subroutine where beta is in 1/eV
   !---------------------------------------------------------------------------!
   subroutine calc_energy_averages_integral(Wk_orig,Lttc,Beta,cutoff,pathOUTPUT,printmode)
      !
      use parameters
      use utils_misc
      use utils_fields
      use file_io
      use omp_lib
      use linalg, only : tensor_transform
      use crystal, only : fill_ksumkdiff, wannierinterpolation
      implicit none
      !
      complex(8),intent(in),target          :: Wk_orig(:,:,:,:)
      type(Lattice),intent(in)              :: Lttc
      real(8),intent(in)                    :: Beta
      real(8),intent(in)                    :: cutoff
      character(len=*),intent(in)           :: pathOUTPUT
      character(len=*),intent(in)           :: printmode
      !
      integer                               :: iw,Nmats
      integer                               :: Ngrid,Norb,Nbp
      integer                               :: iorb,jorb,ib1,ib2,a,b,c,d
      integer                               :: Wk_dim,row,col,ndx,ndx1,ndx2
      integer                               :: ik1,ik2,iq,Nkpt
      integer                               :: iweig,jweig,iE1,iE2
      integer                               :: ithread,Nthread
      real(8)                               :: DosWeights
      integer,allocatable                   :: kptsum(:,:),kptdif(:,:),map(:,:)
      complex(8),allocatable,target         :: Wk_interp(:,:,:,:)
      complex(8),pointer                    :: Wk_used(:,:,:,:)
      complex(8),allocatable                :: Wk_full(:,:)
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      logical                               :: Kstat_exists,Kdyn_exists
      logical                               :: Wrot_exists,calcKernels
      !
      !
      write(*,"(A)") new_line("A")//"---- calc_energy_averages"
      !
      !
      if(.not.initialized)stop "calc_energy_averages: input meshes not initialized. Call Initialize_inputs."
      !
      !Various checks
      Ngrid = size(Egrid)
      Nbp = size(Wk_orig,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      Nmats = size(Wk_orig,dim=3)
      Nkpt = size(kpt_Model,dim=2)
      call assert_shape(Wk_orig,[Nbp,Nbp,Nmats,Lttc%Nkpt],"calc_energy_averages","Wk_orig")
      call assert_shape(kpt_Model,[3,Nkpt],"calc_energy_averages","kpt_Model")
      !
      if(Interpolate2Model)then
         if(Nkpt.eq.Lttc%Nkpt)stop "calc_energy_averages: something is wrong with the K-point dimension (interpolation)."
      else
         if(Nkpt.ne.Lttc%Nkpt)stop "calc_energy_averages: something is wrong with the K-point dimension."
      endif
      !
      allocate(wmats_orig(Nmats));wmats_orig=BosonicFreqMesh(Beta,Nmats)
      wmax_ndx = minloc(abs(wmats_orig-cutoff),dim=1)
      write(*,"(A)") "     Interaction frequency cut at iw_["//str(wmax_ndx)//"]="//str(wmats_orig(wmax_ndx),5)//" eV -> "//str(wmats_orig(wmax_ndx)*eV2DFTgrid,5)//" "//DFTgrid
      write(*,"(A)") "     Bosonic frequency step="//str(abs(wmats_orig(2)-wmats_orig(1)),5)//" eV -> "//str(abs(wmats_orig(2)-wmats_orig(1))*eV2DFTgrid,5)//" "//DFTgrid
      wmats_orig = wmats_orig * eV2DFTgrid
      wmax = wmats_orig(wmax_ndx)
      MatsStep = abs(wmats_orig(2)-wmats_orig(1))
      !
      !check if any of the needed kernels is already printed
      calcKernels=.false.
      !
      if(calc_Int_static)then
         !
         allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Kel_stat.DAT",Kstat_exists,hardstop=.false.,verb=verbose)
         if(Kstat_exists)then
            call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","read")
         else
            calcKernels=.true.
            deallocate(Kel_stat)
         endif
         !
      endif
      !
      if(calc_Int_dynamic)then
         !
         allocate(Wee_dyn(wmax_ndx,Ngrid,Ngrid));Wee_dyn=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Wee_w.DAT",Kdyn_exists,hardstop=.false.,verb=verbose)
         if(Kdyn_exists)then
            call io_Kel(Wee_dyn,reg(pathOUTPUT)//"Wee_w.DAT","read")
         else
            calcKernels=.true.
            deallocate(Wee_dyn)
         endif
         !
      endif
      !
      if(calcKernels)then
         !
         write(*,"(A)") "     One of the required Kernels is missing. Starting operations on screened interacion in Wannier basis."
         !
         if(Interpolate2Model)then
            !
            allocate(Wk_interp(Nbp,Nbp,wmax_ndx,size(kpt_Model,dim=2)));Wk_interp=czero
            call cpu_time(start)
            do iw=1,wmax_ndx
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,kpt_Model,Wk_orig(:,:,iw,:),Wk_interp(:,:,iw,:))
            enddo
            call cpu_time(finish)
            write(*,"(A,F)") "     Interpolating Wk to the new k grid cpu timing:", finish-start
            !
            Wk_used => Wk_interp
            !
         else
            !
            write(*,"(A)") "     Interpolation skipped, linking to original Wk."
            Wk_used => Wk_orig
            !
         endif
         !
         call init_Uelements(Norb,PhysicalUelements)
         !
         !rotation done externally and stored
         Wk_dim = (Nkpt*Norb) * (Nkpt*Norb+1)/2
         allocate(Wk_full(wmax_ndx,Wk_dim));Wk_full=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Wrot.DAT",Wrot_exists,hardstop=.false.,verb=verbose)
         if(Wrot_exists)then
            !
            call cpu_time(start)
            call io_Kel(Wk_full,reg(pathOUTPUT)//"Wrot.DAT","read")
            call cpu_time(finish)
            write(*,"(A,F)") "     Read interaction in band basis cpu timing:", finish-start
            !
         else
            !
            !map between upper triangular and row,col
            call cpu_time(start)
            allocate(map(2,Wk_dim));map=0
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(row,col,ndx)
            !$OMP DO
            do row=1,Nkpt*Norb
               do col=row,Nkpt*Norb
                  !
                  ! upper triangular map (fixed)
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  map(1,ndx) = row
                  map(2,ndx) = col
                  !
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            write(*,"(A,F)") "     Upper triangular Wmap(dim="//str(Wk_dim)//") stored cpu timing:", finish-start
            !
            call fill_ksumkdiff(kpt_Model,kptsum,kptdif)
            deallocate(kptsum)
            !
            !rotation of the interaction W_(ab)(cd)(q) --> W_(i.k1,j.k2)(j.k2,i.k1)
            !W_(i.k1,j.k2)(j.k2,i.k1) = sum_abcd Zdag_ia(k1) * Z_bj(k2) * Zdag_jc(k2) * Z_di(k1) W_(ab)(cd)(q=k1-k2)
            !NOTE: only the upper triangular (c>=r) part of W_(i.k1,j.k2)(j.k2,i.k1) is computed
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(ndx,row,col,iorb,jorb,ik1,ik2,iq),&
            !$OMP PRIVATE(a,b,c,d,ib1,ib2,ithread)
            Nthread = omp_get_num_threads()
            !$OMP DO
            do ndx=1,Wk_dim
               !
               ithread = omp_get_thread_num()
               print *, "thread", ithread, " / ", Nthread, " ndx: ", ndx, " over: ", Wk_dim
               !
               row = map(1,ndx)
               col = map(2,ndx)
               !
               !row = ik1 + (iorb-1)*Nkpt
               iorb = floor((row-0.01)/Nkpt)+1
               ik1  = row - (iorb-1)*Nkpt
               !col = ik2 + (jorb-1)*Nkpt
               jorb = floor((col-0.01)/Nkpt)+1
               ik2  = col - (jorb-1)*Nkpt
               !
               iq = kptdif(ik1,ik2)
               !
               do ib1=1,Nbp
                  !
                  !diagonal elements
                  a = PhysicalUelements%Full_Map(ib1,ib1,1)
                  b = PhysicalUelements%Full_Map(ib1,ib1,2)
                  c = PhysicalUelements%Full_Map(ib1,ib1,3)
                  d = PhysicalUelements%Full_Map(ib1,ib1,4)
                  !
                  Wk_full(:,ndx) = Wk_full(:,ndx)                                                                      &
                                 + Wk_used(ib1,ib1,1:wmax_ndx,iq) * conjg(Zk_Model(a,iorb,ik1)) * Zk_Model(b,jorb,ik2) &
                                                                  * conjg(Zk_Model(c,jorb,ik2)) * Zk_Model(d,iorb,ik1)
                  !
                  do ib2=ib1+1,Nbp
                     !
                     !off-diagonal elements
                     a = PhysicalUelements%Full_Map(ib1,ib2,1)
                     b = PhysicalUelements%Full_Map(ib1,ib2,2)
                     c = PhysicalUelements%Full_Map(ib1,ib2,3)
                     d = PhysicalUelements%Full_Map(ib1,ib2,4)
                     !
                     Wk_full(:,ndx) = Wk_full(:,ndx)                                                                       &
                                    + Wk_used(ib1,ib2,1:wmax_ndx,iq) * conjg(Zk_Model(a,iorb,ik1)) * Zk_Model(b,jorb,ik2)  &
                                                                     * conjg(Zk_Model(c,jorb,ik2)) * Zk_Model(d,iorb,ik1)  &
                                    + Wk_used(ib2,ib1,1:wmax_ndx,iq) * conjg(Zk_Model(c,iorb,ik1)) * Zk_Model(d,jorb,ik2)  &
                                                                     * conjg(Zk_Model(a,jorb,ik2)) * Zk_Model(b,iorb,ik1)
                     !
                  enddo
               enddo
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            deallocate(map,kptdif)
            !
            !print the interaction in band basis depending explicitly of the two kpoints
            call io_Kel(Wk_full,reg(pathOUTPUT)//"Wrot.DAT","write")
            !
            write(*,"(A,F)") "     Rotation of the interaction to band basis cpu timing:", finish-start
            !
         endif
         if(allocated(Wk_interp))deallocate(Wk_interp)
         nullify(Wk_used)
         !
         !Kernels calculation
         call cpu_time(start)
         if(.not.Kstat_exists)then
            allocate(Kel_stat(Ngrid,Ngrid))
            Kel_stat=czero
         endif
         if(.not.Kdyn_exists)then
            allocate(Wee_dyn(wmax_ndx,Ngrid,Ngrid))
            Wee_dyn=czero
         endif
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(iweig,jweig,iE1,iE2,DosWeights),&
         !$OMP PRIVATE(ndx,ndx1,ndx2,row,col),&
         !$OMP PRIVATE(iorb,jorb,ik1,ik2,ithread)
         Nthread = omp_get_num_threads()
         !$OMP DO
         do jweig=1,size(finite_weights_Model,dim=1)
            !
            iE2 = finite_weights_Model(jweig,1)
            jorb = finite_weights_Model(jweig,2)
            ik2 = finite_weights_Model(jweig,3)
            !
            ithread = omp_get_thread_num()
            print *, "thread", ithread, " / ", Nthread, " jweig: ", jweig, " over: ", size(finite_weights_Model,dim=1)
            !
            do iweig=1,size(finite_weights_Model,dim=1)
               !
               iE1 = finite_weights_Model(iweig,1)
               iorb = finite_weights_Model(iweig,2)
               ik1 = finite_weights_Model(iweig,3)
               !
               DosWeights = (weights_Model(iE1,iorb,ik1)/DoS_Model(iE1)) * (weights_Model(iE2,jorb,ik2)/DoS_Model(iE2))
               !
               !product basis map, the indexes spanned by iweig,jweig cover all the possible
               !(ik1,iorb) pairs, so the whole Wk_full matrix, both LT and UT.
               ndx1 = ik1 + (iorb-1)*Nkpt
               ndx2 = ik2 + (jorb-1)*Nkpt
               !
               if(ndx2.ge.ndx1)then
                  !
                  !I'm looking for an element in the UT. ndx gives me the position
                  row = ndx1 !this is the row
                  col = ndx2 !this is the col
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  if(.not.Kstat_exists) Kel_stat(iE1,iE2) = Kel_stat(iE1,iE2) + Wk_full(1,ndx) * eV2DFTgrid * DosWeights
                  if(.not.Kdyn_exists) Wee_dyn(:,iE1,iE2) = Wee_dyn(:,iE1,iE2) + (Wk_full(:,ndx)-Wk_full(1,ndx)) * eV2DFTgrid * DosWeights
                  !
               else
                  !
                  !I'm looking for an element in the LT. I look via ndx his complex conjg in the UT
                  row = ndx2 !this is the row
                  col = ndx1 !this is the col
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  if(.not.Kstat_exists) Kel_stat(iE1,iE2) = Kel_stat(iE1,iE2) + conjg(Wk_full(1,ndx)) * eV2DFTgrid * DosWeights
                  if(.not.Kdyn_exists) Wee_dyn(:,iE1,iE2) = Wee_dyn(:,iE1,iE2) + conjg(Wk_full(:,ndx)-Wk_full(1,ndx)) * eV2DFTgrid * DosWeights
                  !
               endif
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(Wk_full)
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of static electronic Kernel cpu timing:", finish-start
         !
         if(.not.Kstat_exists) call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","write")
         if(.not.Kdyn_exists) call io_Kel(Wee_dyn,reg(pathOUTPUT)//"Wee_w.DAT","write")
         !
      endif
      deallocate(weights_Model,finite_weights_Model)
      !
      if(calc_Int_static.and.(reg(printmode).ne."None"))then
         call print_Kernel("electronic",reg(printmode),reg(pathOUTPUT),"Kel_stat",Egrid,Egrid,Kel_stat)
      endif
      !
      Kernels_stored = .true.
      !
   end subroutine calc_energy_averages_integral
   !
   subroutine calc_energy_averages_list(Inputs,Wk_orig,Lttc,Beta,pathOUTPUT,printmode)
      !
      use parameters
      use utils_misc
      use utils_fields
      use file_io
      use omp_lib
      use linalg, only : tensor_transform
      use crystal, only : fill_ksumkdiff, wannierinterpolation
      implicit none
      !
      type(SCDFT),intent(in)                :: Inputs
      complex(8),intent(in),target          :: Wk_orig(:,:,:,:)
      type(Lattice),intent(in)              :: Lttc
      real(8),intent(in)                    :: Beta
      character(len=*),intent(in)           :: pathOUTPUT
      character(len=*),intent(in)           :: printmode
      !
      integer                               :: iw,Nmats
      integer                               :: Ngrid,Norb,Nbp
      integer                               :: iorb,jorb,ib1,ib2,a,b,c,d
      integer                               :: Wk_dim,row,col,ndx,ndx1,ndx2
      integer                               :: ik1,ik2,iq,Nkpt
      integer                               :: iweig,jweig,iE1,iE2,iT,Efermi_ndx
      integer                               :: ithread,Nthread
      integer                               :: Ngrid_y=100
      real(8)                               :: E1,E2,dT,tanhs!DoS0
      real(8)                               :: DE_m,DE_p,nF_m,nF_p
      real(8)                               :: DosWeights,cutoff
      integer,allocatable                   :: kptsum(:,:),kptdif(:,:),map(:,:)
      real(8),allocatable                   :: Beta_DFT_list(:)
      complex(8),allocatable,target         :: Wk_interp(:,:,:,:)
      complex(8),pointer                    :: Wk_used(:,:,:,:)
      complex(8),allocatable                :: Wk_full(:,:)
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      logical                               :: Kstat_exists,Kdyn_exists
      logical                               :: Wrot_exists,calcKernels
      !
      !
      write(*,"(A)") new_line("A")//"---- calc_energy_averages_list"
      !
      !
      if(.not.initialized)stop "calc_energy_averages_list: input meshes not initialized. Call Initialize_inputs."
      !
      !Various checks
      Ngrid = size(Egrid)
      Nbp = size(Wk_orig,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      Nmats = size(Wk_orig,dim=3)
      Nkpt = size(kpt_Model,dim=2)
      cutoff = Inputs%Wk_cutoff
      call assert_shape(Wk_orig,[Nbp,Nbp,Nmats,Lttc%Nkpt],"calc_energy_averages_list","Wk_orig")
      call assert_shape(kpt_Model,[3,Nkpt],"calc_energy_averages_list","kpt_Model")
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      !DoS0 = DoS_Model(Efermi_ndx)
      !write(*,"(A,F10.5)") new_line("A")//"     get_Kel_integral: Model DoS at the Fermi level:",DoS0
      if(Egrid(Efermi_ndx).ne.0d0) stop "get_Kel_integral: the energy grid requires the E=0 point."
      !
      if(Interpolate2Model)then
         if(Nkpt.eq.Lttc%Nkpt)stop "calc_energy_averages_list: something is wrong with the K-point dimension (interpolation)."
      else
         if(Nkpt.ne.Lttc%Nkpt)stop "calc_energy_averages_list: something is wrong with the K-point dimension."
      endif
      !
      allocate(wmats_orig(Nmats));wmats_orig=BosonicFreqMesh(Beta,Nmats)
      wmax_ndx = minloc(abs(wmats_orig-cutoff),dim=1)
      write(*,"(A)") "     Interaction frequency cut at iw_["//str(wmax_ndx)//"]="//str(wmats_orig(wmax_ndx),5)//" eV -> "//str(wmats_orig(wmax_ndx)*eV2DFTgrid,5)//" "//DFTgrid
      write(*,"(A)") "     Bosonic frequency step="//str(abs(wmats_orig(2)-wmats_orig(1)),5)//" eV -> "//str(abs(wmats_orig(2)-wmats_orig(1))*eV2DFTgrid,5)//" "//DFTgrid
      wmats_orig = wmats_orig * eV2DFTgrid
      wmax = wmats_orig(wmax_ndx)
      MatsStep = abs(wmats_orig(2)-wmats_orig(1))
      !
      !check if any of the needed kernels is already printed
      calcKernels=.false.
      !
      if(calc_Int_static)then
         !
         allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Kel_stat.DAT",Kstat_exists,hardstop=.false.,verb=verbose)
         if(Kstat_exists)then
            call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","read")
         else
            calcKernels=.true.
            deallocate(Kel_stat)
         endif
         !
      endif
      !
      if(calc_Int_dynamic)then
         !
         allocate(Kel_dyn_list(Inputs%Tsteps,Ngrid,Ngrid));Kel_dyn_list=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Kel_dyn_list.DAT",Kdyn_exists,hardstop=.false.,verb=verbose)
         if(Kdyn_exists)then
            call io_Kel(Kel_dyn_list,reg(pathOUTPUT)//"Kel_dyn_list.DAT","read")
         else
            calcKernels=.true.
            deallocate(Kel_dyn_list)
            !
            allocate(Beta_DFT_list(Inputs%Tsteps));Beta_DFT_list=0d0
            do iT=1,Inputs%Tsteps
               dT=0d0
               if(Inputs%Tsteps.gt.1) dT = (iT-1)*abs(Inputs%Tbounds(2)-Inputs%Tbounds(1))/dble(Inputs%Tsteps-1)
               Beta_DFT_list(iT) = 1d0 / ((Inputs%Tbounds(1) + dT)*K2eV*eV2DFTgrid)
            enddo
         endif
         !
      endif
      !
      if(calcKernels)then
         !
         write(*,"(A)") "     One of the required Kernels is missing. Starting operations on screened interacion in Wannier basis."
         !
         if(Interpolate2Model)then
            !
            allocate(Wk_interp(Nbp,Nbp,wmax_ndx,size(kpt_Model,dim=2)));Wk_interp=czero
            call cpu_time(start)
            do iw=1,wmax_ndx
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,kpt_Model,Wk_orig(:,:,iw,:),Wk_interp(:,:,iw,:))
            enddo
            call cpu_time(finish)
            write(*,"(A,F)") "     Interpolating Wk to the new k grid cpu timing:", finish-start
            !
            Wk_used => Wk_interp
            !
         else
            !
            write(*,"(A)") "     Interpolation skipped, linking to original Wk."
            Wk_used => Wk_orig
            !
         endif
         !
         call init_Uelements(Norb,PhysicalUelements)
         !
         !rotation done externally and stored
         Wk_dim = (Nkpt*Norb) * (Nkpt*Norb+1)/2
         allocate(Wk_full(wmax_ndx,Wk_dim));Wk_full=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Wrot.DAT",Wrot_exists,hardstop=.false.,verb=verbose)
         if(Wrot_exists)then
            !
            call cpu_time(start)
            call io_Kel(Wk_full,reg(pathOUTPUT)//"Wrot.DAT","read")
            call cpu_time(finish)
            write(*,"(A,F)") "     Read interaction in band basis cpu timing:", finish-start
            !
         else
            !
            !map between upper triangular and row,col
            call cpu_time(start)
            allocate(map(2,Wk_dim));map=0
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(row,col,ndx)
            !$OMP DO
            do row=1,Nkpt*Norb
               do col=row,Nkpt*Norb
                  !
                  ! upper triangular map (fixed)
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  map(1,ndx) = row
                  map(2,ndx) = col
                  !
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            write(*,"(A,F)") "     Upper triangular Wmap(dim="//str(Wk_dim)//") stored cpu timing:", finish-start
            !
            call fill_ksumkdiff(kpt_Model,kptsum,kptdif)
            deallocate(kptsum)
            !
            !rotation of the interaction W_(ab)(cd)(q) --> W_(i.k1,j.k2)(j.k2,i.k1)
            !W_(i.k1,j.k2)(j.k2,i.k1) = sum_abcd Zdag_ia(k1) * Z_bj(k2) * Zdag_jc(k2) * Z_di(k1) W_(ab)(cd)(q=k1-k2)
            !NOTE: only the upper triangular (c>=r) part of W_(i.k1,j.k2)(j.k2,i.k1) is computed
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(ndx,row,col,iorb,jorb,ik1,ik2,iq),&
            !$OMP PRIVATE(a,b,c,d,ib1,ib2,ithread)
            Nthread = omp_get_num_threads()
            !$OMP DO
            do ndx=1,Wk_dim
               !
               ithread = omp_get_thread_num()
               print *, "thread", ithread, " / ", Nthread, " ndx: ", ndx, " over: ", Wk_dim
               !
               row = map(1,ndx)
               col = map(2,ndx)
               !
               !row = ik1 + (iorb-1)*Nkpt
               iorb = floor((row-0.01)/Nkpt)+1
               ik1  = row - (iorb-1)*Nkpt
               !col = ik2 + (jorb-1)*Nkpt
               jorb = floor((col-0.01)/Nkpt)+1
               ik2  = col - (jorb-1)*Nkpt
               !
               iq = kptdif(ik1,ik2)
               !
               do ib1=1,Nbp
                  !
                  !diagonal elements
                  a = PhysicalUelements%Full_Map(ib1,ib1,1)
                  b = PhysicalUelements%Full_Map(ib1,ib1,2)
                  c = PhysicalUelements%Full_Map(ib1,ib1,3)
                  d = PhysicalUelements%Full_Map(ib1,ib1,4)
                  !
                  Wk_full(:,ndx) = Wk_full(:,ndx)                                                                      &
                                 + Wk_used(ib1,ib1,1:wmax_ndx,iq) * conjg(Zk_Model(a,iorb,ik1)) * Zk_Model(b,jorb,ik2) &
                                                                  * conjg(Zk_Model(c,jorb,ik2)) * Zk_Model(d,iorb,ik1)
                  !
                  do ib2=ib1+1,Nbp
                     !
                     !off-diagonal elements
                     a = PhysicalUelements%Full_Map(ib1,ib2,1)
                     b = PhysicalUelements%Full_Map(ib1,ib2,2)
                     c = PhysicalUelements%Full_Map(ib1,ib2,3)
                     d = PhysicalUelements%Full_Map(ib1,ib2,4)
                     !
                     Wk_full(:,ndx) = Wk_full(:,ndx)                                                                       &
                                    + Wk_used(ib1,ib2,1:wmax_ndx,iq) * conjg(Zk_Model(a,iorb,ik1)) * Zk_Model(b,jorb,ik2)  &
                                                                     * conjg(Zk_Model(c,jorb,ik2)) * Zk_Model(d,iorb,ik1)  &
                                    + Wk_used(ib2,ib1,1:wmax_ndx,iq) * conjg(Zk_Model(c,iorb,ik1)) * Zk_Model(d,jorb,ik2)  &
                                                                     * conjg(Zk_Model(a,jorb,ik2)) * Zk_Model(b,iorb,ik1)
                     !
                  enddo
               enddo
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            deallocate(map,kptdif)
            !
            !print the interaction in band basis depending explicitly of the two kpoints
            call io_Kel(Wk_full,reg(pathOUTPUT)//"Wrot.DAT","write")
            !
            write(*,"(A,F)") "     Rotation of the interaction to band basis cpu timing:", finish-start
            !
         endif
         if(allocated(Wk_interp))deallocate(Wk_interp)
         nullify(Wk_used)
         !
         !Kernels calculation
         call cpu_time(start)
         if(.not.Kstat_exists)then
            allocate(Kel_stat(Ngrid,Ngrid))
            Kel_stat=czero
         endif
         if(.not.Kdyn_exists)then
            allocate(Kel_dyn_list(Inputs%Tsteps,Ngrid,Ngrid))
            Kel_dyn_list=czero
         endif
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(iweig,jweig,iE1,iE2,DosWeights),&
         !$OMP PRIVATE(ndx,ndx1,ndx2,row,col),&
         !$OMP PRIVATE(iorb,jorb,ik1,ik2,ithread),&
         !$OMP PRIVATE(E1,E2,DE_m,DE_p,nF_m,nF_p,iT,tanhs)
         Nthread = omp_get_num_threads()
         !$OMP DO
         do jweig=1,size(finite_weights_Model,dim=1)
            !
            iE2 = finite_weights_Model(jweig,1)
            jorb = finite_weights_Model(jweig,2)
            ik2 = finite_weights_Model(jweig,3)
            !
            ithread = omp_get_thread_num()
            print *, "thread", ithread, " / ", Nthread, " jweig: ", jweig, " over: ", size(finite_weights_Model,dim=1)
            !
            do iweig=1,size(finite_weights_Model,dim=1)
               !
               iE1 = finite_weights_Model(iweig,1)
               iorb = finite_weights_Model(iweig,2)
               ik1 = finite_weights_Model(iweig,3)
               !
               DosWeights = (weights_Model(iE1,iorb,ik1)/DoS_Model(iE1)) * (weights_Model(iE2,jorb,ik2)/DoS_Model(iE2))
               !
               E1 = Ek_Model(iorb,ik1)
               E2 = Ek_Model(jorb,ik2)
               !
               DE_m = E1-E2
               nF_m = fermidirac(+E1,Beta_DFT_list(iT))-fermidirac(+E2,Beta_DFT_list(iT))
               DE_p = E1+E2
               nF_p = fermidirac(-E1,Beta_DFT_list(iT))-fermidirac(+E2,Beta_DFT_list(iT))
               !
               !product basis map, the indexes spanned by iweig,jweig cover all the possible
               !(ik1,iorb) pairs, so the whole Wk_full matrix, both LT and UT.
               ndx1 = ik1 + (iorb-1)*Nkpt
               ndx2 = ik2 + (jorb-1)*Nkpt
               !
               if(ndx2.ge.ndx1)then
                  !
                  !I'm looking for an element in the UT. ndx gives me the position
                  row = ndx1 !this is the row
                  col = ndx2 !this is the col
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  !static Kernel
                  if(.not.Kstat_exists) Kel_stat(iE1,iE2) = Kel_stat(iE1,iE2) + Wk_full(1,ndx) * eV2DFTgrid * DosWeights
                  !
                  !dynamic Kernel stored for all the required temeperatures
                  if(.not.Kdyn_exists)then
                     !
                     do iT=1,Inputs%Tsteps
                        !
                        tanhs = 0d0
                        if((E1*E2).ne.0d0) tanhs = 1d0/(tanh(Beta_DFT_list(iT)*E1/2d0)*tanh(Beta_DFT_list(iT)*E2/2d0))
                        !
                        Kel_dyn_list(iT,iE1,iE2) = Kel_dyn_list(iT,iE1,iE2) + DosWeights * eV2DFTgrid * (2d0/pi) * tanhs *  &
                        (                                                                                                                                                             &
                           nF_m * ( aux_integral(DE_m,(Wk_full(:,ndx)-Wk_full(1,ndx))) + (Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_m)) ) +                         &
                           nF_p * ( aux_integral(DE_p,(Wk_full(:,ndx)-Wk_full(1,ndx))) + (Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_p)) )                           &
                        )
                        !
                     enddo
                     !
                  endif
                  !
               else
                  !
                  !I'm looking for an element in the LT. I look via ndx his complex conjg in the UT
                  row = ndx2 !this is the row
                  col = ndx1 !this is the col
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  !static Kernel
                  if(.not.Kstat_exists) Kel_stat(iE1,iE2) = Kel_stat(iE1,iE2) + conjg(Wk_full(1,ndx)) * eV2DFTgrid * DosWeights
                  !
                  !dynamic Kernel stored for all the required temeperatures
                  if(.not.Kdyn_exists)then
                     !
                     do iT=1,Inputs%Tsteps
                        !
                        tanhs = 0d0
                        if((E1*E2).ne.0d0) tanhs = 1d0/(tanh(Beta_DFT_list(iT)*E1/2d0)*tanh(Beta_DFT_list(iT)*E2/2d0))
                        !
                        Kel_dyn_list(iT,iE1,iE2) = Kel_dyn_list(iT,iE1,iE2) + DosWeights * eV2DFTgrid * (2d0/pi) * tanhs *  &
                        (                                                                                                                                                             &
                           nF_m * ( aux_integral(DE_m,conjg(Wk_full(:,ndx)-Wk_full(1,ndx))) + conjg(Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_m)) ) +               &
                           nF_p * ( aux_integral(DE_p,conjg(Wk_full(:,ndx)-Wk_full(1,ndx))) + conjg(Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_p)) )                 &
                        )
                        !
                     enddo
                     !
                  endif
                  !
               endif
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(Wk_full)
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of static electronic Kernel cpu timing:", finish-start
         !
         !Filling the Fermi lines
         do iT=1,Inputs%Tsteps
            call interpFermi(Kel_dyn_list(iT,:,:),Egrid,Egrid,Efermi_ndx,Efermi_ndx)
         enddo
         !
         if(.not.Kstat_exists) call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","write")
         if(.not.Kdyn_exists) call io_Kel(Kel_dyn_list,reg(pathOUTPUT)//"Kel_dyn_list.DAT","write")
         !
      endif
      deallocate(weights_Model,finite_weights_Model)
      !
      if(calc_Int_static.and.(reg(printmode).ne."None"))then
         call print_Kernel("electronic",reg(printmode),reg(pathOUTPUT),"Kel_stat",Egrid,Egrid,Kel_stat)
      endif
      !
      Kernels_stored = .true.
      !
      !
      !
   contains
      !
      !
      !
      function ygrid(stop,num,ndx) result(yval)
         implicit none
         real(8),intent(in)                    :: stop
         integer,intent(in)                    :: num,ndx
         real(8)                               :: yval
         real(8)                               :: start,step
         !
         if(num.lt.0)stop "ygrid: N<0, abort."
         if(ndx.le.0)stop "ygrid: ndx<=0, abort."
         start = -1d0
         step = (stop-start)/(dble(num)+1d0)
         yval = start + dble(ndx)*step
         !
      end function ygrid
      !
      function aux_integral(DE,W) result(Integral)
         implicit none
         real(8),intent(in)                    :: DE
         complex(8),intent(in)                 :: W(:)
         complex(8)                            :: Integral
         integer                               :: iy,wndx_a,wndx_b
         real(8)                               :: ymax,dy,y_i,y_j,wm
         real(8)                               :: ReW_wm_intp,ImW_wm_intp
         complex(8)                            :: W_wm_i,W_wm_j,Int_i,Int_j
         !
         ymax = (wmax-abs(DE)) / (wmax+abs(DE))
         dy = abs(ygrid(ymax,Ngrid_y,2)-ygrid(ymax,Ngrid_y,1))
         !
         Integral=czero
         if((DE.ne.0d0).and.(ymax.ne.1d0))then
            !
            do iy=2,Ngrid_y
               !
               !continous frequency correspnding to iy
               y_i = ygrid(ymax,Ngrid_y,iy)
               wm = abs(DE) * (1+y_i)/(1-y_i)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               if(wndx_b.gt.wmax_ndx) stop"aux_integral (DE)_i: the frequency index is beyond the cutoff."
               ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(W(wndx_a))] , [wmats_orig(wndx_b),dreal(W(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(W(wndx_a))] , [wmats_orig(wndx_b),dimag(W(wndx_b))] , wm )
               W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !continous frequency correspnding to iy-1
               y_j = ygrid(ymax,Ngrid_y,iy-1)
               wm = abs(DE) * (1+y_j)/(1-y_j)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               if(wndx_b.gt.wmax_ndx) stop"aux_integral (DE)_j: the frequency index is beyond the cutoff."
               ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(W(wndx_a))] , [wmats_orig(wndx_b),dreal(W(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(W(wndx_a))] , [wmats_orig(wndx_b),dimag(W(wndx_b))] , wm )
               W_wm_j = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !integrand for iy and iy-1
               Int_i = W_wm_i / ( 1d0 + y_i**2 )
               Int_j = W_wm_j / ( 1d0 + y_j**2 )
               !
               !trapezoidal integration
               Integral = Integral + sign(1d0,DE)*(Int_i+Int_j)*(dy/2d0)
               !
            enddo
            !
         endif
         !
      end function aux_integral
      !
      !
      !
   end subroutine calc_energy_averages_list



   !================================ KERNELS ==================================!


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !---------------------------------------------------------------------------!
   subroutine get_Kel_integral(Kel,Beta,printmode,printKpath)
      !
      use parameters
      use utils_misc
      implicit none
      !
      complex(8),intent(out)                :: Kel(:,:)
      real(8),intent(in)                    :: Beta
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printKpath
      !
      integer                               :: Efermi_ndx
      integer                               :: iE1,iE2,iy,Ngrid
      integer                               :: wndx_a,wndx_b
      real(8)                               :: Temp,E1,E2,DE!DoS0
      real(8)                               :: ymax,dy,y_i,y_j
      real(8)                               :: wm,ReW_wm_intp,ImW_wm_intp
      complex(8)                            :: W_wm_i,W_wm_j,Int_i,Int_j
      complex(8)                            :: Kel_dyn_e_p,Kel_dyn_e_m
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- get_Kel_integral"
      !
      !
      if(.not.calc_Kel)stop "get_Kel_integral: inputs not initialized. call Initialize_inputs."
      if(.not.Kernels_stored)stop "get_Kel_integral: fully screened interaction not stored. call calc_energy_averages."
      if(calc_Int_static.and.(.not.allocated(Kel_stat)))stop "get_Kel_integral: strange Kel_stat should be allocated."
      if(calc_Int_dynamic.and.(.not.allocated(Wee_dyn)))stop "get_Kel_integral: strange Wee_dyn should be allocated."
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      !DoS0 = DoS_Model(Efermi_ndx)
      !write(*,"(A,F10.5)") new_line("A")//"     get_Kel_integral: Model DoS at the Fermi level:",DoS0
      if(Egrid(Efermi_ndx).ne.0d0) stop "get_Kel_integral: the energy grid requires the E=0 point."
      !
      Ngrid = size(Egrid)
      call assert_shape(Kel,[Ngrid,Ngrid],"get_Kel_integral","Kel")
      !
      Kel=czero
      if(calc_Int_dynamic)then
         !
         call cpu_time(start)
         !$OMP PARALLEL DEFAULT(PRIVATE),&
         !$OMP SHARED(Ngrid,Egrid,Efermi_ndx,Beta,Kel,Wee_dyn),&
         !$OMP SHARED(wmats_orig,MatsStep,wmax,wmax_ndx)
         !$OMP DO
         do iE1=1,Ngrid
            if(iE1.eq.Efermi_ndx)cycle
            do iE2=1,Ngrid
               if(iE2.eq.Efermi_ndx)cycle
               !
               E1=Egrid(iE1)
               E2=Egrid(iE2)
               !
               !one loop for both auxiliary y variables
               Kel_dyn_e_m=czero
               Kel_dyn_e_p=czero
               do iy=2,Ngrid
                  !
                  !-------------------- first term of the sum ---------------------
                  DE = E1 - E2
                  ymax = (wmax-abs(DE)) / (wmax+abs(DE))
                  !
                  if((DE.ne.0d0).and.(ymax.ne.1d0))then
                     !
                     dy = abs(ygrid(ymax,Ngrid,2)-ygrid(ymax,Ngrid,1))
                     !
                     !continous frequency correspnding to iy
                     y_i = ygrid(ymax,Ngrid,iy)
                     wm = abs(DE) * (1+y_i)/(1-y_i)
                     !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                     wndx_a = floor(wm/MatsStep) + 1
                     wndx_b = wndx_a + 1
                     if(wndx_b.gt.wmax_ndx) stop"get_Kel_integral (E1-E2)_i: the frequency index is beyond the cutoff."
                     ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
                     !
                     !continous frequency correspnding to iy-1
                     y_j = ygrid(ymax,Ngrid,iy-1)
                     wm = abs(DE) * (1+y_j)/(1-y_j)
                     !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                     wndx_a = floor(wm/MatsStep) + 1
                     wndx_b = wndx_a + 1
                     if(wndx_b.gt.wmax_ndx) stop"get_Kel_integral (E1-E2)_j: the frequency index is beyond the cutoff."
                     ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     W_wm_j = dcmplx(ReW_wm_intp,ImW_wm_intp)
                     !
                     !integrand for iy and iy-1
                     Int_i = W_wm_i / ( 1d0 + y_i**2 )
                     Int_j = W_wm_j / ( 1d0 + y_j**2 )
                     !
                     !trapezoidal integration
                     Kel_dyn_e_m = Kel_dyn_e_m + sign(1d0,DE)*(Int_i+Int_j)*(dy/2d0)
                     !
                  endif
                  !
                  !
                  !------------------- second term of the sum ---------------------
                  DE = E1 + E2
                  ymax = (wmax-abs(DE)) / (wmax+abs(DE))
                  !
                  if((DE.ne.0d0).and.(ymax.ne.1d0))then
                     !
                     dy = abs(ygrid(ymax,Ngrid,2)-ygrid(ymax,Ngrid,1))
                     !
                     !continous frequency correspnding to iy
                     y_i = ygrid(ymax,Ngrid,iy)
                     wm = abs(DE) * (1+y_i)/(1-y_i)
                     !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                     wndx_a = floor(wm/MatsStep) + 1
                     wndx_b = wndx_a + 1
                     if(wndx_b.gt.wmax_ndx) stop"get_Kel_integral (E1+E2)_i: the frequency index is beyond the cutoff."
                     ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
                     !
                     !continous frequency correspnding to iy-1
                     y_j = ygrid(ymax,Ngrid,iy-1)
                     wm = abs(DE) * (1+y_j)/(1-y_j)
                     !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                     wndx_a = floor(wm/MatsStep) + 1
                     wndx_b = wndx_a + 1
                     if(wndx_b.gt.wmax_ndx) stop"get_Kel_integral (E1+E2)_j: the frequency index is beyond the cutoff."
                     ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,iE1,iE2))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,iE1,iE2))] , wm )
                     W_wm_j = dcmplx(ReW_wm_intp,ImW_wm_intp)
                     !
                     !integrand for iy and iy-1
                     Int_i = W_wm_i / ( 1d0 + y_i**2 )
                     Int_j = W_wm_j / ( 1d0 + y_j**2 )
                     !
                     !trapezoidal integration
                     Kel_dyn_e_p = Kel_dyn_e_p + sign(1d0,DE)*(Int_i+Int_j)*(dy/2d0)
                     !
                  endif
                  !
               enddo
               !
               !adding the tail to w-->inf limit of the interaction
               Kel_dyn_e_m = Kel_dyn_e_m + Wee_dyn(wmax_ndx,iE1,iE2) * ( pi/2d0 - atan2(wmax,(E1-E2)) )
               Kel_dyn_e_p = Kel_dyn_e_p + Wee_dyn(wmax_ndx,iE1,iE2) * ( pi/2d0 - atan2(wmax,(E1+E2)) )
               !
               !adding the fermi function differences
               Kel_dyn_e_m = Kel_dyn_e_m * ( fermidirac(+E1,Beta) - fermidirac(+E2,Beta) )*(2d0/pi)
               Kel_dyn_e_p = Kel_dyn_e_p * ( fermidirac(-E1,Beta) - fermidirac(+E2,Beta) )*(2d0/pi)
               !
               !adding the tanh in front
               Kel(iE1,iE2) = ( Kel_dyn_e_m + Kel_dyn_e_p ) / ( tanh(Beta/2d0*E1) * tanh(Beta/2d0*E2) )
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         !Filling the Fermi lines
         call interpFermi(Kel,Egrid,Egrid,Efermi_ndx,Efermi_ndx)
         !
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of dynamic electronic Kernel cpu timing:", finish-start
         !
         if(reg(printmode).ne."None")then
            Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
            call print_Kernel("electronic",reg(printmode),reg(printKpath),"Kel_dyn_T"//str(Temp,2),Egrid,Egrid,Kel)
         endif
         !
      endif
      !
      if(calc_Int_static)then
         !
         Kel = Kel + Kel_stat
         !
      endif
      !
      !
      !
   contains
      !
      !
      !
      function ygrid(stop,num,ndx) result(yval)
         implicit none
         real(8),intent(in)                    :: stop
         integer,intent(in)                    :: num,ndx
         real(8)                               :: yval
         real(8)                               :: start,step
         !
         if(num.lt.0)stop "ygrid: N<0, abort."
         if(ndx.le.0)stop "ygrid: ndx<=0, abort."
         start = -1d0
         step = (stop-start)/(dble(num)+1d0)
         yval = start + dble(ndx)*step
         !
      end function ygrid
      !
      !
      !
   end subroutine get_Kel_integral
   !
   subroutine get_Kel_list(Kel,iT,Beta,printmode,printKpath)
      !
      use parameters
      use utils_misc
      implicit none
      !
      complex(8),intent(out)                :: Kel(:,:)
      integer,intent(in)                    :: iT
      real(8),intent(in)                    :: Beta
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printKpath
      integer                               :: Ngrid
      real(8)                               :: Temp
      !
      !
      if(verbose)write(*,"(A)") "---- get_Kel_list"
      !
      !
      if(.not.calc_Kel)stop "get_Kel_list: inputs not initialized. call Initialize_inputs."
      if(.not.Kernels_stored)stop "get_Kel_list: fully screened interaction not stored. call calc_energy_averages."
      if(calc_Int_static.and.(.not.allocated(Kel_stat)))stop "get_Kel_list: strange Kel_stat should be allocated."
      if(calc_Int_dynamic.and.(.not.allocated(Kel_dyn_list)))stop "get_Kel_list: strange Kel_dyn_list should be allocated."
      !
      Ngrid = size(Egrid)
      call assert_shape(Kel,[Ngrid,Ngrid],"get_Kel_list","Kel")
      !
      Kel=czero
      if(calc_Int_dynamic)then
         !
         Kel = Kel_dyn_list(iT,:,:)
         !
         if(reg(printmode).ne."None")then
            Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
            call print_Kernel("electronic",reg(printmode),reg(printKpath),"Kel_dyn_T"//str(Temp,2),Egrid,Egrid,Kel)
         endif
         !
      endif
      !
      if(calc_Int_static) Kel = Kel + Kel_stat
      !
   end subroutine get_Kel_list


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic renormalization factor averaged on an energy
   !         grid. This works only if the DFT DOS and energy grid are provided.
   !---------------------------------------------------------------------------!
   subroutine calc_Zph_e(Beta,Zph_e,mode,printZpath)
      !
      use parameters, only : K2eV
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
      real(8),intent(out)                   :: Zph_e(:)
      character(len=*),intent(in)           :: mode
      character(len=*),intent(in)           :: printZpath
      !
      integer                               :: Efermi_ndx,unit
      integer                               :: iE,iE1,iE2,Ngrid
      integer                               :: iomega,Nomega
      real(8)                               :: E1,E2,dE,dw,Temp,DoS0_DFT
      real(8),allocatable                   :: a2F_tmp(:),a2F_int(:)
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Zph_e"
      !
      !
      if(.not.Phonons_stored)stop "calc_Zph_e: a2F(omega) is not stored. call read_a2F."
      if(.not.DFT_DoS_stored)stop "calc_Zph_e: DoS_DFT is not stored. call read_DoS_DFT."
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      DoS0_DFT = DoS_DFT(Efermi_ndx)
      write(*,"(A,F10.5)") new_line("A")//"     calc_Zph_e: DoS_DFT at the Fermi level:",DoS0_DFT
      !
      Ngrid = size(Egrid)
      Nomega = size(omega)
      call assert_shape(Zph_e,[Ngrid],"calc_Zph_e","Zph_e")
      !
      select case(reg(mode))
         case default
            !
            stop "Available E->0 liumits for Zph_e: symrenorm, asym, sym."
            !
         case("symrenorm")
            !
            if(verbose)write (*,"(A)") "     Zph_e: Renormalized term assuming symmetrical DoS around Ef. See PhysRevB 72.024545 eqs. 79-81 for details."
            !
         case("asym")
            !
            if(verbose)write (*,"(A)") "     Zph_e: Used for asymmetrical band structures. Divergence for E->0 smoothed numerically."
            !
         case("sym")
            !
            if(verbose)write (*,"(A)") "     Zph_e: Full term assuming symmetrical DOS around Ef. See PhysRevB 72.024545 eqs. 77-78 for details."
            !
      end select
      !
      call cpu_time(start)
      Zph_e=0d0
      allocate(a2F_tmp(Nomega));a2F_tmp=0d0
      allocate(a2F_int(Ngrid));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE1,iE2,E1,E2,a2F_tmp,a2F_int,iomega,dw,dE)
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         !
         if(iE1.eq.Efermi_ndx)cycle
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
               if(reg(mode).eq."symrenorm")then
                  !
                  a2F_tmp(iomega) = a2F(iomega) * ( J(E1,E2,omega(iomega),Beta) + J(E1,-E2,omega(iomega),Beta) )
                  !
               elseif(reg(mode).eq."asym")then
                  !
                  a2F_tmp(iomega) = a2F(iomega) * ( 2d0*Jasym(E1,E2,omega(iomega),Beta) - Iasym(E1,E2,omega(iomega),Beta) )
                  !
               elseif (reg(mode).eq."sym") then
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
            if((reg(mode).eq."asym").or.(reg(mode).eq."symrenorm"))then
               Zph_e(iE1) = Zph_e(iE1) + ( a2F_int(iE2-1) + a2F_int(iE2) ) * (dE/2d0) * (DoS_DFT(iE2)/DoS0_DFT)
            else
               Zph_e(iE1) = Zph_e(iE1) + ( a2F_int(iE2-1) + a2F_int(iE2) ) * (dE/2d0)
            endif
            !
         enddo
         !
         !extra minus compared to PhysRevB.72.024545 where they have dropped it compared to PhysRevB.88.014514 where they have it
         if(E1.ne.0d0) Zph_e(iE1) = -Zph_e(iE1)/tanh(Beta/2d0*E1)
         !
      enddo !iE1
      !$OMP END DO
      !$OMP BARRIER
      !$OMP END PARALLEL
      deallocate(a2F_tmp,a2F_int)
      !
      !Filling the Fermi line
      call interpFermi(Zph_e,Egrid,Efermi_ndx)
      !
      write(*,"(A,F10.5)")"     lambda(Zph): ",Zph_e(Efermi_ndx)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     Calculation of phononic Z cpu timing:", finish-start
      !
      !Print Z phonon
      call createDir(reg(printZpath),verb=verbose)
      Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
      unit = free_unit()
      open(unit,file=reg(printZpath)//"Zph_e_T"//str(Temp,2)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
      do iE=1,Ngrid
         write(unit,"(2F20.10)")Egrid(iE),Zph_e(iE)
      enddo
      close(unit)
      !
   end subroutine calc_Zph_e


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic kernel averaged on an energy grid. Minimal use
   !         of shared variables in order to be able to try different grid/meshes.
   !---------------------------------------------------------------------------!
   subroutine calc_Kph_e(Beta,Kph_e,printmode,printKpath)
      !
      use parameters, only : K2eV
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Beta
      real(8),intent(out)                   :: Kph_e(:,:)
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printKpath
      !
      integer                               :: Efermi_ndx,iE1,iE2,Ngrid
      integer                               :: iomega,Nomega
      real(8)                               :: E1,E2,dw,Temp,DoS0_DFT
      real(8)                               :: a2F_int
      real(8),allocatable                   :: a2F_tmp(:)
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Kph_e"
      !
      !
      if(.not.Phonons_stored)stop "calc_Kph_e: a2F(omega) is not stored. call read_a2F."
      if(.not.DFT_DoS_stored)stop "calc_Kph_e: DoS_DFT is not stored. call read_DoS_DFT."
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      DoS0_DFT = DoS_DFT(Efermi_ndx)
      write(*,"(A,F10.5)") new_line("A")//"     calc_Kph_e: DoS_DFT at the Fermi level:",DoS0_DFT
      !
      Ngrid = size(Egrid)
      Nomega = size(omega)
      call assert_shape(Kph_e,[Ngrid,Ngrid],"calc_Kph_e","Kph_e")
      !
      call cpu_time(start)
      Kph_e=0d0
      allocate(a2F_tmp(Nomega));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE1,iE2,E1,E2,a2F_tmp,a2F_int,iomega,dw)
      !$OMP DO SCHEDULE(DYNAMIC)
      do iE1=1,Ngrid
         if(iE1.eq.Efermi_ndx)cycle
         do iE2=1,Ngrid
            if(iE2.eq.Efermi_ndx)cycle
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
            !Integral over phononic frequency
            a2F_int=0d0
            do iomega=2,Nomega
               dw = abs(omega(iomega)-omega(iomega-1))
               a2F_int = a2F_int + ( a2F_tmp(iomega-1)+a2F_tmp(iomega) ) * (dw/2d0)
            enddo
            !
            if((E1.ne.0d0).and.(E2.ne.0d0)) Kph_e(ie1,ie2) = (2d0/(tanh(Beta/2d0*E1)*tanh(Beta/2d0*E2))) * a2F_int / DoS0_DFT
            !
         enddo !iE2
      enddo !iE1
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(a2F_tmp)
      !
      !Filling the Fermi lines
      call interpFermi(Kph_e,Egrid,Egrid,Efermi_ndx,Efermi_ndx)
      !
      write(*,"(A,F10.5)")"     lambda(Kph): ",-Kph_e(Efermi_ndx,Efermi_ndx)*DoS0_DFT
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     Calculation of phononic Kernel cpu timing:", finish-start
      !
      !Print Kernel
      if(reg(printmode).ne."None")then
         Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
         call print_Kernel("phononic",reg(printmode),reg(printKpath),"Kph_T"//str(Temp,2),Egrid,Egrid,dcmplx(Kph_e,0d0))
      endif
      !
   end subroutine calc_Kph_e



   !========================== AUXILIARY FUNCTIONS ============================!



   !---------------------------------------------------------------------------!
   !PURPOSE: printing Kernel wrapper
   !---------------------------------------------------------------------------!
   subroutine print_Kernel(Kerneltype,printmode,printpath,filename,Egrid1,Egrid2,Kernel)
      !
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: Kerneltype
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printpath
      character(len=*),intent(in)           :: filename
      real(8),intent(in)                    :: Egrid1(:),Egrid2(:)
      complex(8),intent(in)                 :: Kernel(:,:)
      !
      integer                               :: Ngrid1,Ngrid2,iE,iE1,iE2
      integer                               :: Efermi_ndx1,Efermi_ndx2,unit
      logical                               :: RealK,CmplxK
      !
      !
      if(verbose)write(*,"(A)") "---- print_Kernel"
      !
      !
      if((reg(Kerneltype).ne."electronic").and.(reg(Kerneltype).ne."phononic"))then
         stop "print_Kernel: available Kerneltype are only electronic or phononic."
      endif
      !
      Ngrid1 = size(Egrid1)
      Ngrid2 = size(Egrid2)
      call assert_shape(Kernel,[Ngrid1,Ngrid2],"print_Kernel","Kernel")
      Efermi_ndx1 = minloc(abs(Egrid1),dim=1)
      Efermi_ndx2 = minloc(abs(Egrid2),dim=1)
      !
      CmplxK=.true.
      if(reg(Kerneltype).eq."phononic")CmplxK=.false.
      RealK = .not.CmplxK
      !
      call createDir(reg(printpath),verb=verbose)
      write(*,"(A)") "     Printing "//reg(Kerneltype)//" Kernel with mode "//reg(printmode)//" in "//reg(printpath)
      !
      select case(reg(printmode))
         case default
            !
            if(reg(Kerneltype).eq."phononic")stop "Available print modes for phononic Kernel: E0, diag, surf, all."
            if(reg(Kerneltype).eq."electronic")stop "Available print modes for electronic Kernel: E0, 0E, diag, surf, all."
            !
         case("E0")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid1
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2)),dimag(Kernel(iE,Efermi_ndx2))
            enddo
            close(unit)
            !
         case("0E")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid2
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE)),dimag(Kernel(Efermi_ndx1,iE))
            enddo
            close(unit)
            !
         case("diag")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid2
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,iE))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,iE)),dimag(Kernel(iE,iE))
            enddo
            close(unit)
            !
         case("surf")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid1
               do iE2=1,Ngrid2
                  write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dreal(Kernel(iE1,iE2))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(CmplxK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid1
                  do iE2=1,Ngrid2
                     write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dimag(Kernel(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
            endif
            !
         case("all")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_E0.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid1
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,Efermi_ndx2)),dimag(Kernel(iE,Efermi_ndx2))
            enddo
            close(unit)
            if(reg(Kerneltype).eq."electronic")then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid2
                  if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE))
                  if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Efermi_ndx1,iE)),dimag(Kernel(Efermi_ndx1,iE))
               enddo
               close(unit)
            endif
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid2
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,iE))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(iE,iE)),dimag(Kernel(iE,iE))
            enddo
            close(unit)
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid1
               do iE2=1,Ngrid2
                  write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dreal(Kernel(iE1,iE2))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(CmplxK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid1
                  do iE2=1,Ngrid2
                     write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dimag(Kernel(iE1,iE2))
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
   !PURPOSE: print/read Kernels in binary format
   !---------------------------------------------------------------------------!
   subroutine io_Kel_d2(K,filepath,mode)
      !
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)           :: K(:,:)
      character(len=*),intent(in)        :: filepath
      character(len=*),intent(in)        :: mode
      integer                            :: unit
      integer                            :: D1,D2,D1_,D2_
      integer                            :: iD1,iD2,iD1_,iD2_
      real(8)                            :: ReW,ImW
      !
      !
      if(verbose)write(*,"(A)") "---- io_Kel_d2"
      !
      !
      D1 = size(K,dim=1)
      D2 = size(K,dim=2)
      !
      select case(reg(mode))
         case default
            !
            stop "io_Kel_d2: Available modes are: read, write."
            !
         case("read")
            !
            write(*,"(A)") "     Read "//reg(filepath)
            unit = free_unit()
            open(unit,file=reg(filepath),form="unformatted",status="unknown",position="rewind",action="read")
            !
            read(unit) D1_,D2_
            if(D1_.ne.D1) stop "io_Kel_d2(read): Kel_stat from file has wrong 1st dimension."
            if(D2_.ne.D2) stop "io_Kel_d2(read): Kel_stat from file has wrong 2nd dimension."
            !
            do iD2=1,D2
               read(unit) iD2_
               if(iD2_.ne.iD2) stop "io_Kel_d2(read): wrong iD2 index."
               do iD1=1,D1
                  read(unit) iD1_,ReW,ImW
                  if(iD1_.ne.iD1) stop "io_Kel_d2(read): wrong iD1 index."
                  K(iD1,iD2) = dcmplx(ReW,ImW)
               enddo
            enddo
            close(unit)
            !
         case("write")
            !
            write(*,"(A)") "     Dump "//reg(filepath)//" (binary)"
            unit = free_unit()
            open(unit,file=reg(filepath),form="unformatted",status="unknown",position="rewind",action="write")
            !
            write(unit) D1,D2
            do iD2=1,D2
               write(unit) iD2
               do iD1=1,D1
                  write(unit) iD1,dreal(K(iD1,iD2)),dimag(K(iD1,iD2))
               enddo
            enddo
            close(unit)
            !
      end select
      !
   end subroutine io_Kel_d2
   !
   subroutine io_Kel_d3(K,filepath,mode)
      !
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)           :: K(:,:,:)
      character(len=*),intent(in)        :: filepath
      character(len=*),intent(in)        :: mode
      integer                            :: unit
      integer                            :: D1,D2,D3,D1_,D2_,D3_
      integer                            :: iD1,iD2,iD3,iD1_,iD2_,iD3_
      real(8)                            :: ReW,ImW
      !
      !
      if(verbose)write(*,"(A)") "---- io_Kel_d3"
      !
      !
      D1 = size(K,dim=1)
      D2 = size(K,dim=2)
      D3 = size(K,dim=3)
      !
      select case(reg(mode))
         case default
            !
            stop "io_Kel_d3: Available modes are: read, write."
            !
         case("read")
            !
            write(*,"(A)") "     Read "//reg(filepath)
            unit = free_unit()
            open(unit,file=reg(filepath),form="unformatted",status="unknown",position="rewind",action="read")
            !
            read(unit) D1_,D2_,D3_
            if(D1_.ne.D1) stop "io_Kel_d3(read): Kel_stat from file has wrong 1st dimension."
            if(D2_.ne.D2) stop "io_Kel_d3(read): Kel_stat from file has wrong 2nd dimension."
            if(D3_.ne.D3) stop "io_Kel_d3(read): Kel_stat from file has wrong 3rd dimension."
            !
            do iD3=1,D3
               read(unit) iD3_
               if(iD3_.ne.iD3) stop "io_Kel_d3(read): wrong iD3 index."
               do iD2=1,D2
                  read(unit) iD2_
                  if(iD2_.ne.iD2) stop "io_Kel_d3(read): wrong iD2 index."
                  do iD1=1,D1
                     read(unit) iD1_,ReW,ImW
                     if(iD1_.ne.iD1) stop "io_Kel_d3(read): wrong iD1 index."
                     K(iD1,iD2,iD3) = dcmplx(ReW,ImW)
                  enddo
               enddo
            enddo
            close(unit)
            !
         case("write")
            !
            write(*,"(A)") "     Dump "//reg(filepath)//" (binary)"
            unit = free_unit()
            open(unit,file=reg(filepath),form="unformatted",status="unknown",position="rewind",action="write")
            !
            write(unit) D1,D2,D3
            do iD3=1,D3
               write(unit) iD3
               do iD2=1,D2
                  write(unit) iD2
                  do iD1=1,D1
                     write(unit) iD1,dreal(K(iD1,iD2,iD3)),dimag(K(iD1,iD2,iD3))
                  enddo
               enddo
            enddo
            close(unit)
            !
      end select
      !
   end subroutine io_Kel_d3


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
      integer                         :: iE,NE1,NE2,shift_
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
      integer                         :: iE,NE1,NE2,shift_
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
