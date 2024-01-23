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
   interface interpFermi_Z
      module procedure interpFermi_Z_d
      module procedure interpFermi_Z_z
   end interface interpFermi_Z

   interface interpFermi_K
      module procedure interpFermi_K_mat_d
      module procedure interpFermi_K_mat_z
      module procedure interpFermi_K_tri_d
      module procedure interpFermi_K_tri_z
   end interface interpFermi_K

   interface io_Kel
      module procedure io_Kel_d1
      module procedure io_Kel_d2
      module procedure io_Kel_d3
   end interface io_Kel

   interface print_Kernel
      module procedure print_Kernel_d1
      module procedure print_Kernel_d2
   end interface print_Kernel


   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),allocatable,private              :: omega(:)                         !Phonon energy on logarithmic grid
   real(8),allocatable,private              :: a2F(:)                           !alpha^2*F(\Omega) function
   !
   integer,private                          :: Nkpt3_Model(3)
   real(8),allocatable,private              :: kpt_Model(:,:)
   real(8),allocatable,private              :: nkstar_irred(:)
   real(8),allocatable,private              :: weights_Model(:,:,:)
   integer,allocatable,private              :: finite_weights_Model(:,:)
   integer,allocatable,private              :: finite_weights_Model_irred(:,:)
   !
   real(8),allocatable,private              :: Ek_Model(:,:)
   complex(8),allocatable,private           :: Zk_Model(:,:,:)
   complex(8),allocatable,private           :: Kel_stat(:)
   complex(8),allocatable,private           :: Wee_dyn(:,:)
   !
   logical,private                          :: BZirred=.true.
   integer,private                          :: Nkpt_irred
   !
   logical,private                          :: initialized=.false.
   logical,private                          :: Phonons_stored=.false.
   logical,private                          :: DFT_DoS_stored=.false.
   logical,private                          :: Interpolate2Model=.false.
   logical,private                          :: Kernels_stored=.false.
   logical,private                          :: calc_Int_static=.false.
   logical,private                          :: calc_Int_dynamic=.false.
   !
   integer,private                          :: Ngrid_aux=0                      !grid of the auxiliary variable for frequency integration
   integer,private                          :: wmax_ndx                         !input screened interaction cutoff index
   real(8),private                          :: wmax                             !input screened interaction cutoff
   real(8),private                          :: MatsStep                         !input screened interaction frequency step
   real(8),allocatable,private              :: wmats_orig(:)                    !input Matsubara gird in DFT units
   !
   !public
   character(len=2),public,protected        :: DFTgrid="eV"
   real(8),public,protected                 :: eV2DFTgrid=1d0                   !this converts eV to the DFT grid
   real(8),public,protected                 :: DFTgrid2eV=1d0                   !this converts the DFT grid to eV
   !
   real(8),allocatable,public,protected     :: Egrid_Phonons(:)                 !logarithmic gird
   real(8),allocatable,public,protected     :: Egrid_Model(:)                   !linear grid
   real(8),allocatable,public,protected     :: DoS_Model(:)
   real(8),allocatable,public,protected     :: DoS_DFT(:)
   integer,allocatable,public,protected     :: bilinear_map(:,:,:)
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
      integer                               :: Norb,Nkpt,Ngrid
      integer,allocatable                   :: kptp_dum(:),pkpt_dum(:,:,:)
      complex(8),allocatable                :: Hk_intp(:,:,:)
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- Initialize_inputs"
      !
      !
      if(.not.Inputs%status)stop "Initialize_inputs: input container not properly initialized."
      !
      Ngrid = Inputs%Ngrid
      !
      allocate(Egrid_Phonons(Ngrid));Egrid_Phonons=0d0
      Egrid_Phonons = denspace(2*abs(Inputs%Ebounds(1)),Ngrid,center=.true.,expfact=Inputs%expfact)
      if(Egrid_Phonons(minloc(abs(Egrid_Phonons),dim=1)).ne.0d0) stop "Initialize_inputs: the energy grid (Egrid_Phonons) requires the E=0 point."
      !
      allocate(Egrid_Model(Ngrid));Egrid_Model=0d0
      Egrid_Model = linspace(-abs(Inputs%Ebounds(1))+eps,+abs(Inputs%Ebounds(1))-eps,Ngrid,istart=.true.,iend=.true.)
      if(Egrid_Model(minloc(abs(Egrid_Model),dim=1)).ne.0d0) stop "Initialize_inputs: the energy grid (Egrid_Model) requires the E=0 point."
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
            Egrid_Phonons = Egrid_Model
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
            Egrid_Phonons = Egrid_Phonons*eV2DFTgrid
            Egrid_Model = Egrid_Model*eV2DFTgrid
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
            Ngrid_aux = Inputs%Ngrid_aux
            !
      end select
      calc_Kel = calc_Int_static .or. calc_Int_dynamic
      write(*,"(A)")"     eV->"//DFTgrid//" conversion factor: "//str(eV2DFTgrid,5)
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
         call calc_irredBZ(reg(pathINPUT),Nkpt3_Model,Nkpt_irred,kptp_dum,pkpt_dum,nkstar_irred,kpt_out=kpt_Model)
         deallocate(kptp_dum,pkpt_dum,nkstar_irred)
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
      call tetrahedron_integration(reg(pathINPUT),Hk_intp,Nkpt3_Model,kpt_Model,Egrid_Model,weights_out=weights_Model,DoS_out=DoS_Model)!,pathOUTPUT=reg(pathINPUT))
      call cpu_time(finish)
      write(*,"(A,F)") "     Tetrahedron integration cpu timing:", finish-start
      deallocate(Hk_intp)
      !
      write(*,"(A,F)")"     Smallest weight (Hk/Hk_renorm):",minval(weights_Model)
      !
      !get the info on irreducible BZ
      if(BZirred)then
         !
         call calc_irredBZ(reg(pathINPUT),Nkpt3_Model,Nkpt_irred,kptp_dum,pkpt_dum,nkstar_irred)
         deallocate(kptp_dum,pkpt_dum)
         !
      else
         !
         allocate(nkstar_irred(Nkpt));nkstar_irred=1d0
         Nkpt_irred = Nkpt
         !
      endif
      !
      !reallocate the weight above threshold
      call set_finite_weights(Inputs%DoSthresh,Norb,Ngrid,Nkpt)
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
   subroutine overload_G0W0(pathINPUT,pathINPUTtr,Lttc,eta,DoSthresh,beta)
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
      real(8),intent(in)                    :: eta
      real(8),intent(in)                    :: DoSthresh
      real(8),intent(in),optional           :: beta
      !
      type(FermionicField)                  :: S_G0W0
      complex(8),allocatable                :: Uwan(:,:,:,:),Vxc(:,:,:,:),Zk(:,:,:)
      complex(8),allocatable                :: Smat(:,:,:,:),Smat_interp(:,:,:,:,:)
      complex(8),allocatable                :: Gmat(:,:,:,:,:),invGf(:,:)
      real(8),allocatable                   :: Egrid_eV(:),wreal_read(:)
      real(8),allocatable                   :: Adiag(:,:,:,:),Aloc(:,:)
      real(8)                               :: ReS,ImS,DoSnorm,dE,Emin,Emax
      integer                               :: ik,iw,iE,iorb,ispin
      integer                               :: Norb,Nkpt,Ngrid
      integer                               :: Nreal_read = 5000
      logical                               :: paramagnet=.true. !this is default for G0W0
      logical                               :: UWAN_exist,Vxc_exist,G0W0_exist
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- overload_G0W0"
      !
      !
      if(.not.initialized) stop "overload_G0W0: gap equation module not properly initialized. Call Initialize_inputs first."
      if(size(Zk_Model,dim=1).ne.Lttc%Norb) stop "overload_G0W0: orbital dimension mismatch."
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
      write(*,"(A)")"     Using mu="//str(Lttc%mu,4)//"[eV] -> "//str(Lttc%mu * eV2DFTgrid,4)//"["//DFTgrid//"]"
      !
      Norb = size(Zk_Model,dim=1)
      Nkpt = size(Zk_Model,dim=3)
      Ngrid = size(Egrid_Model)
      !
      allocate(Egrid_eV(Ngrid));Egrid_eV = Egrid_Model*DFTgrid2eV
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
      allocate(Vxc(Norb,Norb,Nkpt,1));Vxc=czero
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
      !rotate G0W0 input to LDA basis
      allocate(Smat(Norb,Nreal_read,Nkpt,Nspin));Smat=czero
      do ispin=1,Nspin
         !
         do ik=1,Nkpt
            do iw=1,Nreal_read
               Smat(:,iw,ik,ispin) = diagonal(rotate(S_G0W0%wks(:,:,iw,ik,ispin),Zk(:,:,ik)))
            enddo
         enddo
         if(paramagnet)then
            Smat(:,:,:,Nspin) = Smat(:,:,:,1)
            exit
         endif
         !
      enddo
      call DeallocateFermionicField(S_G0W0)
      !
      !interpolate to logarithmic real-frequency mesh and rotate back to Wannier basis
      call cpu_time(start)
      allocate(Smat_interp(Norb,Norb,Ngrid,Nkpt,Nspin));Smat_interp=czero
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE,ik,iorb,ispin,ReS,ImS)
      !$OMP DO
      do iE=1,Ngrid
         do ik=1,Nkpt
            do ispin=1,Nspin
               !
               do iorb=1,Norb
                  ReS = cubic_interp( wreal_read(1:Nreal_read), dreal(Smat(iorb,:,ik,ispin)), Egrid_eV(iE) )
                  ImS = cubic_interp( wreal_read(1:Nreal_read), dimag(Smat(iorb,:,ik,ispin)), Egrid_eV(iE) )
                  if(ImS.gt.0d0)ImS=0d0
                  Smat_interp(iorb,iorb,iE,ik,ispin) = dcmplx(ReS,ImS)
               enddo
               Smat_interp(:,:,iE,ik,ispin) = rotate(Smat_interp(:,:,iE,ik,ispin),dag(Zk(:,:,ik))) - Vxc(:,:,ik,ispin)
               if(paramagnet)then
                  Smat_interp(:,:,iE,ik,Nspin) = Smat_interp(:,:,iE,ik,1)
                  cycle
               endif
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(wreal_read,Smat,Vxc,Zk)
      !
      !compute G0W0 Green's function in [eV]
      allocate(Gmat(Norb,Norb,Ngrid,Nkpt,Nspin));Gmat=czero
      allocate(invGf(Norb,Norb));invGf=czero
      do ispin=1,Nspin
         !
         do ik=1,Nkpt
            do iE=1,Ngrid
               invGf = zeye(Norb)*dcmplx(Egrid_eV(iE)+Lttc%mu,eta) - Lttc%Hk(:,:,ik) - Smat_interp(:,:,iE,ik,ispin)
               call inv(invGf)
               Gmat(:,:,iE,ik,ispin) = invGf
            enddo
         enddo
         if(paramagnet)then
            Gmat(:,:,:,:,Nspin) = Gmat(:,:,:,:,1)
            exit
         endif
         !
      enddo
      deallocate(invGf,Smat_interp)
      !
      !extract spectral function in the diagonal basis
      allocate(Adiag(Ngrid,Norb,Nkpt,Nspin));Adiag=czero
      do ispin=1,Nspin
         !
         do ik=1,Nkpt
            do iw=1,Ngrid
               Adiag(iw,:,ik,ispin) = -dimag(diagonal(rotate(Gmat(:,:,iw,ik,ispin),Lttc%Zk(:,:,ik))))
            enddo
         enddo
         if(paramagnet)then
            Adiag(:,:,:,Nspin) = Adiag(:,:,:,1)
            exit
         endif
         !
      enddo
      deallocate(Gmat)
      !
      !cutoff the spectra with a smooth function
      if(present(beta))then
         !
         Emin = Egrid_eV(minloc(abs(Egrid_eV - 0.9*Egrid_eV(1)),dim=1))
         Emax = Egrid_eV(minloc(abs(Egrid_eV - 0.9*Egrid_eV(Ngrid)),dim=1))
         write(*,"(A)")"     Smoothing the spectra with beta="//str(beta,2)//" eV^-1 between "//str(Emin,4)//" eV and "//str(Emax,4)//" eV."
         do iE=1,Ngrid
            Adiag(iE,:,:,:) = Adiag(iE,:,:,:) * fermidirac(Egrid_eV(iE),Emin,-beta) * fermidirac(Egrid_eV(iE),Emax,beta)
         enddo
         !
      endif
      !
      !normalize to density
      do ispin=1,Nspin
         do ik=1,Nkpt
            do iorb=1,Norb
               DoSnorm=0d0
               do iE=2,Ngrid
                  dE = abs(Egrid_eV(iE)-Egrid_eV(iE-1))
                  DoSnorm = DoSnorm + (Adiag(iE,iorb,ik,ispin)+Adiag(iE-1,iorb,ik,ispin)) * (dE/2d0)
               enddo
               Adiag(:,iorb,ik,ispin) = Adiag(:,iorb,ik,ispin)/(DoSnorm*Nkpt)
            enddo
         enddo
      enddo
      !
      !print local spectra used in LDA basis
      allocate(Aloc(Ngrid,Norb));Aloc=0d0
      do iorb=1,Norb
         !
         do ik=1,Nkpt
            Aloc(:,iorb) = Aloc(:,iorb) + Adiag(:,iorb,ik,1)
         enddo
         !This is just a check that the DoS for each band is normalized to 1
         DoSnorm=0d0
         do iE=2,Ngrid
            dE = abs(Egrid_eV(iE)-Egrid_eV(iE-1))
            DoSnorm = DoSnorm + ( Aloc(iE-1,iorb)+Aloc(iE,iorb) ) * (dE/2d0)
         enddo
         write(*,"(A,F)") "     Spectral function normalization orbital #"//str(iorb)//":", DoSnorm
         call dump_Field_component(Aloc(:,iorb),reg(pathINPUTtr)//"G0W0plots/Aloc_Gap_Equation_Renorm_G0W0/","Aw_G0W0_o"//str(iorb)//"_s1.DAT",Egrid_eV)
         !
      enddo
      deallocate(Aloc,Egrid_eV)
      !
      !convert diagonal spectral function to DFT grid and link to weights_Model and DoS_Model
      weights_Model = Adiag(:,:,:,1) / eV2DFTgrid
      deallocate(Adiag)
      !
      DoS_Model=0d0
      do iE=1,Ngrid
         do ik=1,Nkpt
            do iorb=1,Norb
               DoS_Model(iE) = DoS_Model(iE) + weights_Model(iE,iorb,ik)
            enddo
         enddo
      enddo
      write(*,"(A,F)")"     Smallest weight (G0W0):",minval(weights_Model)
      call cpu_time(start)
      write(*,"(A,F)") "     Cpu timing:", finish-start
      !
      !reallocate the weight above threshold
      call set_finite_weights(DoSthresh,Norb,Ngrid,Nkpt)
      !
   end subroutine overload_G0W0


   !---------------------------------------------------------------------------!
   !PURPOSE: Overload the DoS and rotation with the results computed via the
   !         G0W0 self-energy. Works only if all the required data are present
   !---------------------------------------------------------------------------!
   subroutine overload_DMFT(path2MaxEnt,pathINPUTtr,Lttc,DoSthresh,beta)
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
      character(len=*),intent(in)           :: pathINPUTtr
      type(Lattice),intent(in)              :: Lttc
      real(8),intent(in)                    :: DoSthresh
      real(8),intent(in),optional           :: beta
      !
      real(8),allocatable                   :: ImG_read(:,:,:)
      real(8),allocatable                   :: Egrid_eV(:),wreal_read(:)
      real(8),allocatable                   :: Adiag(:,:,:,:),Aloc(:,:)
      real(8)                               :: DoSnorm,dE,Emin,Emax
      integer                               :: ik,iw,iE,iorb,ispin
      integer                               :: Norb,Nkpt,Ngrid
      integer                               :: Nreal_old,Nreal_read
      integer                               :: unit,ierr
      character(len=256)                    :: path
      logical                               :: paramagnet=.true. !this is default not implemented otherwise
      logical                               :: DMFT_exists
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- overload_DMFT"
      !
      !
      if(.not.initialized) stop "overload_DMFT: gap equation module not properly initialized. Call Initialize_inputs first."
      if(size(Zk_Model,dim=1).ne.Lttc%Norb) stop "overload_DMFT: orbital dimension mismatch."
      if(Interpolate2Model) write(*,"(A)")"     Cannot read DMFT quantities in another k-grid. Ignoring call to subroutine."
      !
      Norb = size(Zk_Model,dim=1)
      Nkpt = size(Zk_Model,dim=3)
      Ngrid = size(Egrid_Model)
      !
      allocate(Egrid_eV(Ngrid));Egrid_eV = Egrid_Model*DFTgrid2eV
      allocate(Adiag(Ngrid,Norb,Nkpt,Nspin));Adiag=czero
      !
      !read MaxEnt data and convert to local spectra
      do ispin=1,Nspin
         !
         !check that all the required k-points are present
         do ik=1,Nkpt
            !
            path = reg(path2MaxEnt)//"MaxEnt_Gk_full_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT_dos.dat"
            !
            call inquireFile(reg(path),DMFT_exists,hardstop=.false.,verb=.true.)
            if(.not.DMFT_exists)then
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
            if((ik.gt.1).and.(Nreal_read.ne.Nreal_old))then
               write(*,"(A,1I5)") "     Real frequency mesh is not consistent among K-points. Ignoring call to subroutine."
               return
            endif
            Nreal_old=Nreal_read
            !
         enddo
         !
         !read each k-point
         allocate(wreal_read(Nreal_read));wreal_read=0d0
         allocate(ImG_read(Norb,Nreal_read,Nkpt));ImG_read=0d0
         do ik=1,Nkpt
            !
            path = reg(path2MaxEnt)//"MaxEnt_Gk_full_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT_dos.dat"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            do iw=1,Nreal_read
               read(unit,*) wreal_read(iw),(ImG_read(iorb,iw,ik),iorb=1,Norb)
            enddo
            close(unit)
            !
         enddo
         write(*,"(A)") "     MaxEnt output on Green's function is read."
         !
         !interpolate to the Model energy grid
         do ik=1,Nkpt
            do iorb=1,Norb
               Adiag(:,iorb,ik,ispin) = cubic_interp(wreal_read, ImG_read(iorb,:,ik), Egrid_eV)
               where(Adiag(:,iorb,ik,ispin).lt.0d0)Adiag(:,iorb,ik,ispin)=0d0
            enddo
         enddo
         deallocate(wreal_read,ImG_read)
         !
         if(paramagnet)then
            Adiag(:,:,:,Nspin) = Adiag(:,:,:,1)
            exit
         endif
         !
      enddo
      !
      !cutoff the spectra with a smooth function
      if(present(beta))then
         !
         Emin = Egrid_eV(minloc(abs(Egrid_eV - 0.9*Egrid_eV(1)),dim=1))
         Emax = Egrid_eV(minloc(abs(Egrid_eV - 0.9*Egrid_eV(Ngrid)),dim=1))
         write(*,"(A)")"     Smoothing the spectra with beta="//str(beta,2)//" eV^-1 between "//str(Emin,4)//" eV and "//str(Emax,4)//" eV."
         do iE=1,Ngrid
            Adiag(iE,:,:,:) = Adiag(iE,:,:,:) * fermidirac(Egrid_eV(iE),Emin,-beta) * fermidirac(Egrid_eV(iE),Emax,beta)
         enddo
         !
      endif
      !
      !normalize to density
      do ispin=1,Nspin
         do ik=1,Nkpt
            do iorb=1,Norb
               DoSnorm=0d0
               do iE=2,Ngrid
                  dE = abs(Egrid_eV(iE)-Egrid_eV(iE-1))
                  DoSnorm = DoSnorm + (Adiag(iE,iorb,ik,ispin)+Adiag(iE-1,iorb,ik,ispin)) * (dE/2d0)
               enddo
               Adiag(:,iorb,ik,ispin) = Adiag(:,iorb,ik,ispin)/(DoSnorm*Nkpt)
            enddo
         enddo
      enddo
      !
      !print local spectra used in LDA basis
      allocate(Aloc(Ngrid,Norb));Aloc=0d0
      do iorb=1,Norb
         !
         do ik=1,Nkpt
            Aloc(:,iorb) = Aloc(:,iorb) + Adiag(:,iorb,ik,1)
         enddo
         !This is just a check that the DoS for each band is normalized to 1
         DoSnorm=0d0
         do iE=2,Ngrid
            dE = abs(Egrid_eV(iE)-Egrid_eV(iE-1))
            DoSnorm = DoSnorm + ( Aloc(iE-1,iorb)+Aloc(iE,iorb) ) * (dE/2d0)
         enddo
         write(*,"(A,F)") "     Spectral function normalization orbital #"//str(iorb)//":", DoSnorm
         call dump_Field_component(Aloc(:,iorb),reg(pathINPUTtr)//"G0W0plots/Aloc_Gap_Equation_Renorm_DMFT/","Aw_DMFT_o"//str(iorb)//"_s1.DAT",Egrid_eV)
         !
      enddo
      deallocate(Aloc,Egrid_eV)
      !
      !convert diagonal spectral function to DFT grid and link to weights_Model and DoS_Model
      weights_Model = Adiag(:,:,:,1) / eV2DFTgrid
      deallocate(Adiag)
      !
      DoS_Model=0d0
      do iE=1,Ngrid
         do ik=1,Nkpt
            do iorb=1,Norb
               DoS_Model(iE) = DoS_Model(iE) + weights_Model(iE,iorb,ik)
            enddo
         enddo
      enddo
      write(*,"(A,F)")"     Smallest weight (DMFT):",minval(weights_Model)
      call cpu_time(start)
      write(*,"(A,F)") "     Cpu timing:", finish-start
      !
      !reallocate the weight above threshold
      call set_finite_weights(DoSthresh,Norb,Ngrid,Nkpt)
      !
   end subroutine overload_DMFT


   !---------------------------------------------------------------------------!
   !PURPOSE: Given the weights allocate the lists linking to those above 
   !         threshold
   !---------------------------------------------------------------------------!
   subroutine set_finite_weights(threshold,Norb,Ngrid,Nkpt)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: threshold
      integer,intent(in)                    :: Norb
      integer,intent(in)                    :: Ngrid
      integer,intent(in)                    :: Nkpt
      !
      integer                               :: iE,iorb,ik,Nweights
      logical                               :: keep_weight,recompute_DoS_Model
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- set_finite_weights"
      !
      !
      if(BZirred.and.(Nkpt.eq.Nkpt_irred)) stop "set_finite_weights: BZirred is true but Nkpt.eq.Nkpt_irred"
      !
      call cpu_time(start)
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(threshold) )
               if(keep_weight) Nweights = Nweights + 1
            enddo
         enddo
      enddo
      !
      recompute_DoS_Model=.false.
      if(Ngrid*Norb*Nkpt.ne.Nweights)then
         DoS_Model=0d0
         recompute_DoS_Model=.true.
      endif
      !
      if(allocated(finite_weights_Model))deallocate(finite_weights_Model)
      allocate(finite_weights_Model(Nweights,3));finite_weights_Model=0
      Nweights=0
      do iE=1,Ngrid
         do iorb=1,Norb
            do ik=1,Nkpt
               keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(threshold) )
               if(keep_weight)then
                  Nweights = Nweights + 1
                  finite_weights_Model(Nweights,1) = iE
                  finite_weights_Model(Nweights,2) = iorb
                  finite_weights_Model(Nweights,3) = ik
                  if(recompute_DoS_Model) DoS_Model(iE) = DoS_Model(iE) + weights_Model(iE,iorb,ik)
               endif
            enddo
         enddo
      enddo
      call cpu_time(finish)
      if(recompute_DoS_Model)then
         write(*,"(A,F)") "     Reduction of DoS integration points (full BZ) from "//str(Ngrid*Norb*Nkpt)//" to "//str(Nweights)//" Cpu timing:", finish-start
      else
         write(*,"(A,F)") "     DoS integration points: "//str(Nweights)//" Cpu timing:", finish-start
      endif
      !
      if(BZirred)then
         !
         call cpu_time(start)
         Nweights=0
         do iE=1,Ngrid
            do iorb=1,Norb
               do ik=1,Nkpt_irred
                  keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(threshold) )
                  if(keep_weight) Nweights = Nweights + 1
               enddo
            enddo
         enddo
         !
         recompute_DoS_Model = Ngrid*Norb*Nkpt_irred.ne.Nweights
         !
         if(allocated(finite_weights_Model_irred))deallocate(finite_weights_Model_irred)
         allocate(finite_weights_Model_irred(Nweights,3));finite_weights_Model_irred=0
         Nweights=0
         do iE=1,Ngrid
            do iorb=1,Norb
               do ik=1,Nkpt_irred
                  keep_weight = ( abs(weights_Model(iE,iorb,ik)).gt.abs(threshold) )
                  if(keep_weight)then
                     Nweights = Nweights + 1
                     finite_weights_Model_irred(Nweights,1) = iE
                     finite_weights_Model_irred(Nweights,2) = iorb
                     finite_weights_Model_irred(Nweights,3) = ik
                  endif
               enddo
            enddo
         enddo
         call cpu_time(finish)
         if(recompute_DoS_Model)then
            write(*,"(A,F)") "     Reduction of DoS integration points (irred BZ) from "//str(Ngrid*Norb*Nkpt_irred)//" to "//str(Nweights)//" Cpu timing:", finish-start
         else
            write(*,"(A,F)") "     DoS integration points (irred BZ): "//str(Nweights)//" Cpu timing:", finish-start
         endif
         !
      else
         !
         if(allocated(finite_weights_Model_irred))deallocate(finite_weights_Model_irred)
         finite_weights_Model_irred = finite_weights_Model
         !
      endif
      !
   end subroutine set_finite_weights


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
      Ngrid = size(Egrid_Phonons)
      allocate(DoS_DFT(Ngrid));DoS_DFT=0d0
      do iE=1,Ngrid
         DoS_DFT(iE) = cubic_interp( Egrid_read, DoS_DFT_read, Egrid_Phonons(iE) )
      enddo
      where(DoS_DFT.lt.0d0)DoS_DFT=0d0
      deallocate(Egrid_read,DoS_DFT_read)
      !
      write(*,"(A,F10.5)") "     Fermi Energy (interp-"//DFTgrid//"): ",Egrid_Phonons(minloc(abs(Egrid_Phonons),dim=1))
      write(*,"(A,F10.5)") "     DoS at Fermi (interp-1/"//DFTgrid//"): ",DoS_DFT(minloc(abs(Egrid_Phonons),dim=1))
      DoS_DFT(minloc(abs(Egrid_Phonons),dim=1)) = DoS0_DFT
      !
      DFT_DoS_stored=.true.
      !
   end subroutine read_DoS_DFT


   !---------------------------------------------------------------------------!
   !PURPOSE: Store the interpolated and rotated NaNb components of the fully
   !         screened interaction. This is the only subroutine where beta is in 1/eV
   !---------------------------------------------------------------------------!
   subroutine calc_energy_averages(Wk_orig,Lttc,Beta,cutoff,pathOUTPUT,printmode)
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
      integer                               :: Wk_dim,Wndx,row,col,ndx1,ndx2
      integer                               :: ik1,ik2,iq,Nkpt
      integer                               :: EE_dim,Endx,iweig,jweig,iE1,iE2
      integer                               :: ithread,Nthread
      real(8)                               :: DosWeights,Elog
      integer,allocatable                   :: kptsum(:,:),kptdif(:,:)
      integer,allocatable                   :: map(:,:),linear_map(:,:)
      complex(8),allocatable,target         :: Wk_interp(:,:,:,:)
      complex(8),pointer                    :: Wk_used(:,:,:,:)
      complex(8),allocatable                :: Wk_full(:,:)
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      logical                               :: Kstat_exists,Kdyn_exists
      logical                               :: calc_Kel_stat,calc_Wee_dyn
      logical                               :: Wrot_exists
      !
      !
      write(*,"(A)") new_line("A")//"---- calc_energy_averages"
      !
      !
      if(.not.initialized)stop "calc_energy_averages: input meshes not initialized. Call Initialize_inputs."
      !
      !Various checks
      Ngrid = size(Egrid_Model)
      Nbp = size(Wk_orig,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      Nmats = size(Wk_orig,dim=3)
      Nkpt = size(kpt_Model,dim=2)
      EE_dim = Ngrid * (Ngrid+1)/2
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
      calc_Kel_stat=.false.
      calc_Wee_dyn=.false.
      !
      if(calc_Int_static)then
         !
         allocate(Kel_stat(EE_dim));Kel_stat=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Kel_stat.DAT",Kstat_exists,hardstop=.false.,verb=verbose)
         if(Kstat_exists)then
            call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","read")
         else
            calc_Kel_stat=.true.
            deallocate(Kel_stat)
         endif
         !
      endif
      !
      if(calc_Int_dynamic)then
         !
         allocate(Wee_dyn(wmax_ndx,EE_dim));Wee_dyn=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Wee_w.DAT",Kdyn_exists,hardstop=.false.,verb=verbose)
         if(Kdyn_exists)then
            call io_Kel(Wee_dyn,reg(pathOUTPUT)//"Wee_w.DAT","read")
         else
            calc_Wee_dyn=.true.
            deallocate(Wee_dyn)
         endif
         !
      endif
      !
      if(calc_Kel_stat.or.calc_Wee_dyn)then
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
            !$OMP PRIVATE(row,col,Wndx)
            !$OMP DO
            do row=1,Nkpt*Norb
               do col=row,Nkpt*Norb
                  !
                  !upper triangular map (fixed)
                  Wndx = rc2ut(row,col,Nkpt*Norb) !(Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  map(1,Wndx) = row
                  map(2,Wndx) = col
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
            !$OMP PRIVATE(Wndx,row,col,iorb,jorb,ik1,ik2,iq),&
            !$OMP PRIVATE(a,b,c,d,ib1,ib2,ithread)
            Nthread = omp_get_num_threads()
            !$OMP DO
            do Wndx=1,Wk_dim
               !
               ithread = omp_get_thread_num()
               print *, "thread", ithread, " / ", Nthread, " Wndx: ", Wndx, " over: ", Wk_dim
               !
               row = map(1,Wndx)
               col = map(2,Wndx)
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
                  Wk_full(:,Wndx) = Wk_full(:,Wndx)                                                                      &
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
                     Wk_full(:,Wndx) = Wk_full(:,Wndx)                                                                       &
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
         if(calc_Kel_stat)then
            allocate(Kel_stat(EE_dim))
            Kel_stat=czero
         endif
         if(calc_Wee_dyn)then
            allocate(Wee_dyn(wmax_ndx,EE_dim))
            Wee_dyn=czero
         endif
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(iweig,jweig,iE1,iE2,DosWeights),&
         !$OMP PRIVATE(Endx,Wndx,ndx1,ndx2),&
         !$OMP PRIVATE(iorb,jorb,ik1,ik2,ithread)
         Nthread = omp_get_num_threads()
         !$OMP DO
         do jweig=1,size(finite_weights_Model,dim=1)
            !
            iE2 = finite_weights_Model(jweig,1) !E_col
            jorb = finite_weights_Model(jweig,2)
            ik2 = finite_weights_Model(jweig,3)
            !
            ithread = omp_get_thread_num()
            print *, "thread", ithread, " / ", Nthread, " jweig: ", jweig, " over: ", size(finite_weights_Model,dim=1)
            !
            do iweig=1,size(finite_weights_Model_irred,dim=1)
               !
               iE1 = finite_weights_Model_irred(iweig,1) !E_row
               iorb = finite_weights_Model_irred(iweig,2)
               ik1 = finite_weights_Model_irred(iweig,3)
               !
               !Storing the UT (E1,E2) matrix
               if(iE1.lt.iE2) cycle
               if(DoS_Model(iE1).eq.0d0) cycle
               if(DoS_Model(iE2).eq.0d0) cycle
               !upper triangular map (fixed)
               Endx = rc2ut(iE1,iE2,Ngrid) !Endx = Ngrid*(iE1-1) - (iE1-1)*iE1/2 + iE2
               !
               DosWeights = (weights_Model(iE1,iorb,ik1)/DoS_Model(iE1)) * (weights_Model(iE2,jorb,ik2)/DoS_Model(iE2)) * nkstar_irred(ik1)
               !
               !product basis map, the indexes spanned by iweig,jweig cover all the possible
               !(ik1,iorb) pairs, so the whole Wk_full matrix, both LT and UT.
               ndx1 = ik1 + (iorb-1)*Nkpt
               ndx2 = ik2 + (jorb-1)*Nkpt
               !
               if(ndx2.ge.ndx1)then
                  !
                  !I'm looking for an element in the UT. Wndx gives me the position
                  !row = ndx1 !this is the row
                  !col = ndx2 !this is the col
                  !Wndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  Wndx = rc2ut(ndx1,ndx2,Nkpt*Norb) 
                  !
                  if(calc_Kel_stat) Kel_stat(Endx) = Kel_stat(Endx) + Wk_full(1,Wndx) * eV2DFTgrid * DosWeights
                  if(calc_Wee_dyn)  Wee_dyn(:,Endx) = Wee_dyn(:,Endx) + (Wk_full(:,Wndx)-Wk_full(1,Wndx)) * eV2DFTgrid * DosWeights
                  !
               else
                  !
                  !I'm looking for an element in the LT. I look via Wndx his complex conjg in the UT
                  !row = ndx2 !this is the row
                  !col = ndx1 !this is the col
                  !Wndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  Wndx = rc2ut(ndx1,ndx2,Nkpt*Norb)
                  !
                  if(calc_Kel_stat) Kel_stat(Endx) = Kel_stat(Endx) + conjg(Wk_full(1,Wndx)) * eV2DFTgrid * DosWeights
                  if(calc_Wee_dyn)  Wee_dyn(:,Endx) = Wee_dyn(:,Endx) + conjg(Wk_full(:,Wndx)-Wk_full(1,Wndx)) * eV2DFTgrid * DosWeights
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
         if(calc_Kel_stat) call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","write")
         if(calc_Wee_dyn) call io_Kel(Wee_dyn,reg(pathOUTPUT)//"Wee_w.DAT","write")
         !
      endif
      deallocate(weights_Model,finite_weights_Model)
      !
      if(calc_Int_static.and.(reg(printmode).ne."None"))then
         call print_Kernel("electronic",reg(printmode),reg(pathOUTPUT),"Kel_stat",Egrid_Model,Kel_stat)
      endif
      !
      !set up the bilinear map 
      if(calc_phonons)then
         !
         allocate(linear_map(Ngrid,2));linear_map=0
         do iE1=1,Ngrid!index of the phonon
            !
            Elog = Egrid_Phonons(iE1)
            Endx = minloc(abs(Egrid_Model-Elog),dim=1)
            !
            if(Elog.le.Egrid_Model(Endx))then
               !
               !Elog is approaching from below or is on top
               linear_map(iE1,1) = Endx-1
               linear_map(iE1,2) = Endx
               !
            else
               !
               !Elog is approaching from above
               linear_map(iE1,1) = Endx
               linear_map(iE1,2) = Endx+1
               !
            endif
            !
         enddo
         !
         allocate(bilinear_map(Ngrid,Ngrid,4));bilinear_map=0
         do iE1=1,Ngrid
            do iE2=1,Ngrid
               bilinear_map(iE1,iE2,1:2) = linear_map(iE1,:)
               bilinear_map(iE1,iE2,3:4) = linear_map(iE2,:)
            enddo
         enddo
         deallocate(linear_map)
         !
      endif
      !
      Kernels_stored = .true.
      !
   end subroutine calc_energy_averages


   !================================ KERNELS ==================================!


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !---------------------------------------------------------------------------!
   subroutine get_Kel(Kel,Beta,printmode,printKpath)
      !
      use parameters
      use utils_misc
      implicit none
      !
      complex(8),intent(out)                :: Kel(:)
      real(8),intent(in)                    :: Beta
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printKpath
      !
      integer                               :: Efermi_ndx
      integer                               :: iE1,iE2,iy,Ngrid
      integer                               :: EE_dim,Endx
      integer                               :: wndx_a,wndx_b
      integer,allocatable                   :: map(:,:)
      real(8)                               :: Temp,E1,E2,DE
      real(8)                               :: ymax,dy,y_i,y_j
      real(8)                               :: wm,ReW_wm_intp,ImW_wm_intp
      complex(8)                            :: W_wm_i,W_wm_j,Int_i,Int_j
      complex(8)                            :: Kel_dyn_e_p,Kel_dyn_e_m
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- get_Kel"
      !
      !
      if(.not.calc_Kel)stop "get_Kel: inputs not initialized. call Initialize_inputs."
      if(.not.Kernels_stored)stop "get_Kel: fully screened interaction not stored. call calc_energy_averages."
      if(calc_Int_static.and.(.not.allocated(Kel_stat)))stop "get_Kel: strange Kel_stat should be allocated."
      if(calc_Int_dynamic.and.(.not.allocated(Wee_dyn)))stop "get_Kel: strange Wee_dyn should be allocated."
      !
      Efermi_ndx = minloc(abs(Egrid_Model),dim=1)
      if(Egrid_Model(Efermi_ndx).ne.0d0) stop "get_Kel: the energy grid requires the E=0 point."
      !
      Ngrid = size(Egrid_Model)
      EE_dim = Ngrid * (Ngrid+1)/2
      call assert_shape(Kel,[EE_dim],"get_Kel","Kel")
      !
      Kel=czero
      if(calc_Int_dynamic)then
         !
         !upper triangular map (fixed)
         allocate(map(2,EE_dim));map=0
         do iE1=1,Ngrid
            do iE2=iE1,Ngrid
               Endx = rc2ut(iE1,iE2,Ngrid) ! Endx = Ngrid*(iE1-1) - (iE1-1)*iE1/2 + iE2
               map(1,Endx) = iE1
               map(2,Endx) = iE2
            enddo
         enddo
         !
         call cpu_time(start)
         !$OMP PARALLEL DEFAULT(PRIVATE),&
         !$OMP SHARED(Ngrid,Ngrid_aux,Egrid_Model,Efermi_ndx,Beta,Kel,Wee_dyn),&
         !$OMP SHARED(wmats_orig,MatsStep,wmax,wmax_ndx,map,EE_dim)
         !$OMP DO
         do Endx=1,EE_dim
            !
            iE1 = map(1,Endx)
            iE2 = map(2,Endx)
            !
            if(iE1.eq.Efermi_ndx)cycle
            if(iE2.eq.Efermi_ndx)cycle
            !
            E1 = Egrid_Model(iE1)
            E2 = Egrid_Model(iE2)
            !
            !one loop for both auxiliary y variables
            Kel_dyn_e_m=czero
            Kel_dyn_e_p=czero
            do iy=2,Ngrid_aux
               !
               !-------------------- first term of the sum ---------------------
               DE = E1 - E2
               ymax = (wmax-abs(DE)) / (wmax+abs(DE))
               !
               if((DE.ne.0d0).and.(ymax.ne.1d0))then
                  !
                  dy = abs(ygrid(ymax,Ngrid_aux,2)-ygrid(ymax,Ngrid_aux,1))
                  !
                  !continous frequency correspnding to iy
                  y_i = ygrid(ymax,Ngrid_aux,iy)
                  wm = abs(DE) * (1+y_i)/(1-y_i)
                  !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                  wndx_a = floor(wm/MatsStep) + 1
                  wndx_b = wndx_a + 1
                  if(wndx_b.gt.wmax_ndx) stop"get_Kel (E1-E2)_i: the frequency index is beyond the cutoff."
                  ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,Endx))] , wm )
                  ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,Endx))] , wm )
                  W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
                  !
                  !continous frequency correspnding to iy-1
                  y_j = ygrid(ymax,Ngrid_aux,iy-1)
                  wm = abs(DE) * (1+y_j)/(1-y_j)
                  !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                  wndx_a = floor(wm/MatsStep) + 1
                  wndx_b = wndx_a + 1
                  if(wndx_b.gt.wmax_ndx) stop"get_Kel (E1-E2)_j: the frequency index is beyond the cutoff."
                  ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,Endx))] , wm )
                  ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,Endx))] , wm )
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
                  dy = abs(ygrid(ymax,Ngrid_aux,2)-ygrid(ymax,Ngrid_aux,1))
                  !
                  !continous frequency correspnding to iy
                  y_i = ygrid(ymax,Ngrid_aux,iy)
                  wm = abs(DE) * (1+y_i)/(1-y_i)
                  !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                  wndx_a = floor(wm/MatsStep) + 1
                  wndx_b = wndx_a + 1
                  if(wndx_b.gt.wmax_ndx) stop"get_Kel (E1+E2)_i: the frequency index is beyond the cutoff."
                  ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,Endx))] , wm )
                  ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,Endx))] , wm )
                  W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
                  !
                  !continous frequency correspnding to iy-1
                  y_j = ygrid(ymax,Ngrid_aux,iy-1)
                  wm = abs(DE) * (1+y_j)/(1-y_j)
                  !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
                  wndx_a = floor(wm/MatsStep) + 1
                  wndx_b = wndx_a + 1
                  if(wndx_b.gt.wmax_ndx) stop"get_Kel (E1+E2)_j: the frequency index is beyond the cutoff."
                  ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dreal(Wee_dyn(wndx_b,Endx))] , wm )
                  ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(Wee_dyn(wndx_a,Endx))] , [wmats_orig(wndx_b),dimag(Wee_dyn(wndx_b,Endx))] , wm )
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
            Kel_dyn_e_m = Kel_dyn_e_m + Wee_dyn(wmax_ndx,Endx) * ( pi/2d0 - atan2(wmax,(E1-E2)) )
            Kel_dyn_e_p = Kel_dyn_e_p + Wee_dyn(wmax_ndx,Endx) * ( pi/2d0 - atan2(wmax,(E1+E2)) )
            !
            !adding the fermi function differences
            Kel_dyn_e_m = Kel_dyn_e_m * ( fermidirac(+E1,Beta) - fermidirac(+E2,Beta) )*(2d0/pi)
            Kel_dyn_e_p = Kel_dyn_e_p * ( fermidirac(-E1,Beta) - fermidirac(+E2,Beta) )*(2d0/pi)
            !
            !adding the tanh in front
            Kel(Endx) = ( Kel_dyn_e_m + Kel_dyn_e_p ) / ( tanh(Beta/2d0*E1) * tanh(Beta/2d0*E2) )
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(map)
         !
         !Filling the Fermi lines
         call interpFermi_K(Kel,Egrid_Model,Efermi_ndx)
         !
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of dynamic electronic Kernel cpu timing:", finish-start
         !
         if(reg(printmode).ne."None")then
            Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
            call print_Kernel("electronic",reg(printmode),reg(printKpath),"Kel_dyn_T"//str(Temp,2),Egrid_Model,Kel)
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
   end subroutine get_Kel


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
      Efermi_ndx = minloc(abs(Egrid_Phonons),dim=1)
      DoS0_DFT = DoS_DFT(Efermi_ndx)
      write(*,"(A,F10.5)") new_line("A")//"     calc_Zph_e: DoS_DFT at the Fermi level:",DoS0_DFT
      !
      Ngrid = size(Egrid_Phonons)
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
            E1=Egrid_Phonons(iE1)
            E2=Egrid_Phonons(iE2)
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
            dE = abs(Egrid_Phonons(iE2)-Egrid_Phonons(iE2-1))
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
      call interpFermi_Z(Zph_e,Egrid_Phonons,Efermi_ndx)
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
         write(unit,"(2F20.10)")Egrid_Phonons(iE),Zph_e(iE)
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
      real(8),intent(out)                   :: Kph_e(:)
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printKpath
      !
      integer                               :: Efermi_ndx,iE1,iE2,Ngrid
      integer                               :: iomega,Nomega,EE_dim,Endx
      real(8)                               :: E1,E2,dw,Temp,DoS0_DFT
      real(8)                               :: a2F_int
      real(8),allocatable                   :: a2F_tmp(:)
      integer,allocatable                   :: map(:,:)
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Kph_e"
      !
      !
      if(.not.Phonons_stored)stop "calc_Kph_e: a2F(omega) is not stored. call read_a2F."
      if(.not.DFT_DoS_stored)stop "calc_Kph_e: DoS_DFT is not stored. call read_DoS_DFT."
      !
      Efermi_ndx = minloc(abs(Egrid_Phonons),dim=1)
      DoS0_DFT = DoS_DFT(Efermi_ndx)
      write(*,"(A,F10.5)") new_line("A")//"     calc_Kph_e: DoS_DFT at the Fermi level:",DoS0_DFT
      !
      Ngrid = size(Egrid_Phonons)
      EE_dim = Ngrid * (Ngrid+1)/2
      call assert_shape(Kph_e,[EE_dim],"calc_Kph_e","Kph_e")
      !
      Nomega = size(omega)
      !
      !upper triangular map (fixed)
      allocate(map(2,EE_dim));map=0
      do iE1=1,Ngrid
         do iE2=iE1,Ngrid
            Endx = rc2ut(iE1,iE2,Ngrid) ! Endx = Ngrid*(iE1-1) - (iE1-1)*iE1/2 + iE2
            map(1,Endx) = iE1
            map(2,Endx) = iE2
         enddo
      enddo
      !
      call cpu_time(start)
      Kph_e=0d0
      allocate(a2F_tmp(Nomega));a2F_int=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE1,iE2,E1,E2,a2F_tmp,a2F_int,iomega,dw,Endx)
      !$OMP DO SCHEDULE(DYNAMIC)
      do Endx=1,EE_dim
         !
         iE1 = map(1,Endx)
         iE2 = map(2,Endx)
         !
         if(iE1.eq.Efermi_ndx)cycle
         if(iE2.eq.Efermi_ndx)cycle
         !
         E1 = Egrid_Phonons(iE1)
         E2 = Egrid_Phonons(iE2)
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
         if((E1.ne.0d0).and.(E2.ne.0d0)) Kph_e(Endx) = (2d0/(tanh(Beta/2d0*E1)*tanh(Beta/2d0*E2))) * a2F_int / DoS0_DFT
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(a2F_tmp,map)
      !
      !Filling the Fermi lines
      call interpFermi_K(Kph_e,Egrid_Phonons,Efermi_ndx)
      !
      Endx = rc2ut(Efermi_ndx,Efermi_ndx,Ngrid)
      write(*,"(A,F10.5)")"     lambda(Kph): ",-Kph_e(Endx)*DoS0_DFT
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     Calculation of phononic Kernel cpu timing:", finish-start
      !
      !Print Kernel
      if(reg(printmode).ne."None")then
         Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
         call print_Kernel("phononic",reg(printmode),reg(printKpath),"Kph_T"//str(Temp,2),Egrid_Phonons,dcmplx(Kph_e,0d0))
      endif
      !
   end subroutine calc_Kph_e



   !========================== AUXILIARY FUNCTIONS ============================!



   !---------------------------------------------------------------------------!
   !PURPOSE: printing Kernels in readable format
   !---------------------------------------------------------------------------!
   subroutine print_Kernel_d1(Kerneltype,printmode,printpath,filename,Egrid,Kernel)
      !
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: Kerneltype
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printpath
      character(len=*),intent(in)           :: filename
      real(8),intent(in)                    :: Egrid(:)
      complex(8),intent(in)                 :: Kernel(:)
      !
      integer                               :: Ngrid,iE,iE1,iE2,Endx
      integer                               :: Efermi_ndx,unit
      logical                               :: RealK,CmplxK
      !
      !
      if(verbose)write(*,"(A)") "---- print_Kernel_d1"
      !
      !
      if((reg(Kerneltype).ne."electronic").and.(reg(Kerneltype).ne."phononic"))then
         stop "print_Kernel_d1: available Kerneltype are only electronic or phononic."
      endif
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      Ngrid = size(Egrid)
      call assert_shape(Kernel,[Ngrid*(Ngrid+1)/2],"print_Kernel_d1","Kernel")
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
            do iE=1,Ngrid
               Endx = rc2ut(iE,Efermi_ndx,Ngrid)
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx)),dimag(Kernel(Endx))
            enddo
            close(unit)
            !
         case("0E")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid
               Endx = rc2ut(Efermi_ndx,iE,Ngrid)
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx)),dimag(Kernel(Endx))
            enddo
            close(unit)
            !
         case("diag")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid
               Endx = rc2ut(iE,iE,Ngrid)
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx)),dimag(Kernel(Endx))
            enddo
            close(unit)
            !
         case("surf")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid
               do iE2=1,Ngrid
                  Endx = rc2ut(iE1,iE2,Ngrid)
                  write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dreal(Kernel(Endx))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(CmplxK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     Endx = rc2ut(iE1,iE2,Ngrid)
                     write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dimag(Kernel(Endx))
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
            do iE=1,Ngrid
               Endx = rc2ut(iE,Efermi_ndx,Ngrid)
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx)),dimag(Kernel(Endx))
            enddo
            close(unit)
            if(reg(Kerneltype).eq."electronic")then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid
                  Endx = rc2ut(Efermi_ndx,iE,Ngrid)
                  if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx))
                  if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx)),dimag(Kernel(Endx))
               enddo
               close(unit)
            endif
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid
               Endx = rc2ut(iE,iE,Ngrid)
               if(RealK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx))
               if(CmplxK) write(unit,"(3E20.12)")Egrid(iE),dreal(Kernel(Endx)),dimag(Kernel(Endx))
            enddo
            close(unit)
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid
               do iE2=1,Ngrid
                  Endx = rc2ut(iE1,iE2,Ngrid)
                  write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dreal(Kernel(Endx))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(CmplxK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid
                  do iE2=1,Ngrid
                     Endx = rc2ut(iE1,iE2,Ngrid)
                     write(unit,"(3E20.12)")Egrid(iE1),Egrid(iE2),dimag(Kernel(Endx))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
            endif
            !
            !
      end select
      !
   end subroutine print_Kernel_d1
   !
   subroutine print_Kernel_d2(Kerneltype,printmode,printpath,filename,Egrid1,Egrid2,Kernel)
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
      if(verbose)write(*,"(A)") "---- print_Kernel_d2"
      !
      !
      if((reg(Kerneltype).ne."electronic").and.(reg(Kerneltype).ne."phononic"))then
         stop "print_Kernel_d2: available Kerneltype are only electronic or phononic."
      endif
      !
      Efermi_ndx1 = minloc(abs(Egrid1),dim=1)
      Efermi_ndx2 = minloc(abs(Egrid2),dim=1)
      !
      Ngrid1 = size(Egrid1)
      Ngrid2 = size(Egrid2)
      call assert_shape(Kernel,[Ngrid1,Ngrid2],"print_Kernel_d2","Kernel")
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
               if(RealK) write(unit,"(3E20.12)")Egrid1(iE),dreal(Kernel(iE,Efermi_ndx2))
               if(CmplxK) write(unit,"(3E20.12)")Egrid1(iE),dreal(Kernel(iE,Efermi_ndx2)),dimag(Kernel(iE,Efermi_ndx2))
            enddo
            close(unit)
            !
         case("0E")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid2
               if(RealK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(Efermi_ndx1,iE))
               if(CmplxK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(Efermi_ndx1,iE)),dimag(Kernel(Efermi_ndx1,iE))
            enddo
            close(unit)
            !
         case("diag")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid2
               if(RealK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(iE,iE))
               if(CmplxK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(iE,iE)),dimag(Kernel(iE,iE))
            enddo
            close(unit)
            !
         case("surf")
            !
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid1
               do iE2=1,Ngrid2
                  write(unit,"(3E20.12)")Egrid1(iE1),Egrid2(iE2),dreal(Kernel(iE1,iE2))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(CmplxK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid1
                  do iE2=1,Ngrid2
                     write(unit,"(3E20.12)")Egrid1(iE1),Egrid2(iE2),dimag(Kernel(iE1,iE2))
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
               if(RealK) write(unit,"(3E20.12)")Egrid1(iE),dreal(Kernel(iE,Efermi_ndx2))
               if(CmplxK) write(unit,"(3E20.12)")Egrid1(iE),dreal(Kernel(iE,Efermi_ndx2)),dimag(Kernel(iE,Efermi_ndx2))
            enddo
            close(unit)
            if(reg(Kerneltype).eq."electronic")then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_0E.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE=1,Ngrid2
                  if(RealK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(Efermi_ndx1,iE))
                  if(CmplxK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(Efermi_ndx1,iE)),dimag(Kernel(Efermi_ndx1,iE))
               enddo
               close(unit)
            endif
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_diag.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE=1,Ngrid2
               if(RealK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(iE,iE))
               if(CmplxK) write(unit,"(3E20.12)")Egrid2(iE),dreal(Kernel(iE,iE)),dimag(Kernel(iE,iE))
            enddo
            close(unit)
            unit = free_unit()
            open(unit,file=reg(printpath)//reg(filename)//"_surf_R.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iE1=1,Ngrid1
               do iE2=1,Ngrid2
                  write(unit,"(3E20.12)")Egrid1(iE1),Egrid2(iE2),dreal(Kernel(iE1,iE2))
               enddo
               write(unit,*)
            enddo
            close(unit)
            if(CmplxK) then
               unit = free_unit()
               open(unit,file=reg(printpath)//reg(filename)//"_surf_I.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iE1=1,Ngrid1
                  do iE2=1,Ngrid2
                     write(unit,"(3E20.12)")Egrid1(iE1),Egrid2(iE2),dimag(Kernel(iE1,iE2))
                  enddo
                  write(unit,*)
               enddo
               close(unit)
            endif
            !
            !
      end select
      !
   end subroutine print_Kernel_d2


   !---------------------------------------------------------------------------!
   !PURPOSE: print/read Kernels in binary format
   !---------------------------------------------------------------------------!
   subroutine io_Kel_d1(K,filepath,mode)
      !
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)           :: K(:)
      character(len=*),intent(in)        :: filepath
      character(len=*),intent(in)        :: mode
      integer                            :: unit
      integer                            :: D1,D1_
      integer                            :: iD1,iD1_
      real(8)                            :: ReW,ImW
      !
      !
      if(verbose)write(*,"(A)") "---- io_Kel_d1"
      !
      !
      D1 = size(K)
      !
      select case(reg(mode))
         case default
            !
            stop "io_Kel_d1: Available modes are: read, write."
            !
         case("read")
            !
            write(*,"(A)") "     Read "//reg(filepath)
            unit = free_unit()
            open(unit,file=reg(filepath),form="unformatted",status="unknown",position="rewind",action="read")
            !
            read(unit) D1_
            if(D1_.ne.D1) stop "io_Kel_d1(read): Kel_stat from file has wrong 1st dimension."
            !
            do iD1=1,D1
               read(unit) iD1_,ReW,ImW
               if(iD1_.ne.iD1) stop "io_Kel_d1(read): wrong iD1 index."
               K(iD1) = dcmplx(ReW,ImW)
            enddo
            close(unit)
            !
         case("write")
            !
            write(*,"(A)") "     Dump "//reg(filepath)//" (binary)"
            unit = free_unit()
            open(unit,file=reg(filepath),form="unformatted",status="unknown",position="rewind",action="write")
            !
            write(unit) D1
            do iD1=1,D1
               write(unit) iD1,dreal(K(iD1)),dimag(K(iD1))
            enddo
            close(unit)
            !
      end select
      !
   end subroutine io_Kel_d1
   !
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
      implicit none
      double precision,intent(in)           :: E,Beta
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
         yA = Z(iEo-1-shift_) ; xA = E(iEo-1-shift_)
         yB = Z(iEo-2-shift_) ; xB = E(iEo-2-shift_)
         y_S = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !linear fit of the (E,0) line - north
         yA = Z(iEo+1+shift_) ; xA = E(iEo+1+shift_)
         yB = Z(iEo+2+shift_) ; xB = E(iEo+2+shift_)
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
         yA = Z(iEo-1-shift_) ; xA = E(iEo-1-shift_)
         yB = Z(iEo-2-shift_) ; xB = E(iEo-2-shift_)
         Rey_S = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_S = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !linear fit of the (E,0) line - north
         yA = Z(iEo+1+shift_) ; xA = E(iEo+1+shift_)
         yB = Z(iEo+2+shift_) ; xB = E(iEo+2+shift_)
         Rey_N = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_N = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         Z(iEo) = dcmplx(( Rey_S + Rey_N ),( Imy_S + Imy_N ))/2d0
         !
      enddo
      !
   end subroutine interpFermi_Z_z
   !
   subroutine interpFermi_K_mat_d(K,E1,E2,iE1o,iE2o,shift)
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
      if(verbose)write(*,"(A)") "---- interpFermi_K_mat_d"
      !
      !
      NE1 = size(E1)
      NE2 = size(E2)
      if(abs(E1(iE1o)).gt.1d-9) stop "interpFermi_K_mat_d: provided index iE1o does not correspond to Fermi level."
      if(abs(E2(iE2o)).gt.1d-9) stop "interpFermi_K_mat_d: provided index iE2o does not correspond to Fermi level."
      call assert_shape(K,[NE1,NE2],"interpFermi_K_mat_d","K")
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
   end subroutine interpFermi_K_mat_d
   !
   subroutine interpFermi_K_mat_z(K,E1,E2,iE1o,iE2o,shift)
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
      if(verbose)write(*,"(A)") "---- interpFermi_K_mat_z"
      !
      !
      NE1 = size(E1)
      NE2 = size(E2)
      if(abs(E1(iE1o)).gt.1d-9) stop "interpFermi_K_mat_z: provided index iE1o does not correspond to Fermi level."
      if(abs(E2(iE2o)).gt.1d-9) stop "interpFermi_K_mat_z: provided index iE2o does not correspond to Fermi level."
      call assert_shape(K,[NE1,NE2],"interpFermi_K_mat_z","K")
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
   end subroutine interpFermi_K_mat_z
   !
   subroutine interpFermi_K_tri_d(K,E,iEo,shift)
      use utils_misc
      implicit none
      real(8),intent(inout)           :: K(:)
      real(8),intent(in)              :: E(:)
      integer,intent(in)              :: iEo
      integer,intent(in),optional     :: shift
      !
      integer                         :: iE,NE,Kdim,shift_
      integer                         :: Endx0,Endx1,Endx2
      real(8)                         :: yA,yB,xA,xB
      real(8)                         :: y_N,y_S,y_E,y_W
      !
      !
      if(verbose)write(*,"(A)") "---- interpFermi_K_tri_d"
      !
      !
      NE = size(E)
      if(abs(E(iEo)).gt.1d-9) stop "interpFermi_K_tri_d: provided index iE1o does not correspond to Fermi level."
      Kdim = NE*(NE+1)/2
      call assert_shape(K,[Kdim],"interpFermi_K_tri_d","K")
      !
      shift_=0
      if(present(shift))shift_=shift
      !
      !linear fit of the (E1,0)=(iE,iEo) column
      do iE=1,NE
         !
         if(iE.eq.iEo)cycle
         !
         !at each row iE - east direction
         Endx1 = rc2ut(iE,iEo+1+shift_,NE)
         Endx2 = rc2ut(iE,iEo+2+shift_,NE)
         !
         xA = E(iEo+1+shift_); yA = K(Endx1)
         xB = E(iEo+2+shift_); yB = K(Endx2)
         y_E = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !at each row iE - west direction
         Endx1 = rc2ut(iE,iEo-1+shift_,NE)
         Endx2 = rc2ut(iE,iEo-2+shift_,NE)
         !
         xA = E(iEo-1-shift_); yA = K(Endx1)
         xB = E(iEo-2-shift_); yB = K(Endx2)
         y_W = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !the iE,iEo column is the east/west average
         Endx0 = rc2ut(iE,iEo,NE)
         K(Endx0) = ( y_E + y_W )/2d0
         !
      enddo
      !
      !linear fit of the (0,E)=(iEo,iE) row
      do iE=1,NE
         !
         if(iE.eq.iEo)cycle
         !
         !at each column iE - north direction
         Endx1 = rc2ut(iEo+1+shift_,iE,NE)
         Endx2 = rc2ut(iEo+2+shift_,iE,NE)
         xA = E(iEo+1+shift_); yA = K(Endx1)
         xB = E(iEo+2+shift_); yB = K(Endx2)
         y_N = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !at each column iE - south direction
         Endx1 = rc2ut(iEo-1+shift_,iE,NE)
         Endx2 = rc2ut(iEo-2+shift_,iE,NE)
         xA = E(iEo-1-shift_); yA = K(Endx1)
         xB = E(iEo-2-shift_); yB = K(Endx2)
         y_S = linear_interp_2y([xA,yA],[xB,yB],0d0)
         !
         !the iEo,iE row is the north/south average
         Endx0 = rc2ut(iEo,iE,NE)
         K(Endx0) = ( y_N + y_S )/2d0
         !
      enddo
      !
      !linear fit around the (0,0)=(iEo,iEo) point
      !east direction
      Endx1 = rc2ut(iEo,iEo+1+shift_,NE)
      Endx2 = rc2ut(iEo,iEo+2+shift_,NE)
      xA = E(iEo+1+shift_); yA = K(Endx1)
      xB = E(iEo+2+shift_); yB = K(Endx2)
      y_E = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !west direction
      Endx1 = rc2ut(iEo,iEo-1-shift_,NE)
      Endx2 = rc2ut(iEo,iEo-2-shift_,NE)
      xA = E(iEo-1-shift_); yA = K(Endx1)
      xB = E(iEo-2-shift_); yB = K(Endx2)
      y_W = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !north direction
      Endx1 = rc2ut(iEo+1+shift_,iEo,NE)
      Endx2 = rc2ut(iEo+2+shift_,iEo,NE)
      xA = E(iEo+1+shift_); yA = K(Endx1)
      xB = E(iEo+2+shift_); yB = K(Endx2)
      y_N = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !south direction
      Endx1 = rc2ut(iEo-1+shift_,iEo,NE)
      Endx2 = rc2ut(iEo-2+shift_,iEo,NE)
      xA = E(iEo-1-shift_); yA = K(Endx1)
      xB = E(iEo-2-shift_); yB = K(Endx2)
      y_S = linear_interp_2y([xA,yA],[xB,yB],0d0)
      !
      !the iEo,iEo point is the east/west/north/south average
      Endx0 = rc2ut(iEo,iEo,NE)
      K(Endx0) = ( y_S + y_N + y_W + y_E )/4d0
      !
   end subroutine interpFermi_K_tri_d
   !
   subroutine interpFermi_K_tri_z(K,E,iEo,shift)
      use utils_misc
      implicit none
      complex(8),intent(inout)        :: K(:)
      real(8),intent(in)              :: E(:)
      integer,intent(in)              :: iEo
      integer,intent(in),optional     :: shift
      !
      integer                         :: iE,NE,Kdim,shift_
      integer                         :: Endx0,Endx1,Endx2
      complex(8)                      :: yA,yB
      real(8)                         :: xA,xB
      real(8)                         :: Rey_N,Rey_S,Rey_E,Rey_W
      real(8)                         :: Imy_N,Imy_S,Imy_E,Imy_W
      !
      !
      if(verbose)write(*,"(A)") "---- interpFermi_K_tri_z"
      !
      !
      NE = size(E)
      if(abs(E(iEo)).gt.1d-9) stop "interpFermi_K_tri_z: provided index iE1o does not correspond to Fermi level."
      Kdim = NE*(NE+1)/2
      call assert_shape(K,[Kdim],"interpFermi_K_tri_z","K")
      !
      shift_=0
      if(present(shift))shift_=shift
      !
      !linear fit of the (E1,0)=(iE,iEo) column
      do iE=1,NE
         !
         if(iE.eq.iEo)cycle
         !
         !at each row iE - east direction
         Endx1 = rc2ut(iE,iEo+1+shift_,NE)
         Endx2 = rc2ut(iE,iEo+2+shift_,NE)
         !
         xA = E(iEo+1+shift_); yA = K(Endx1)
         xB = E(iEo+2+shift_); yB = K(Endx2)
         Rey_E = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_E = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !at each row iE - west direction
         Endx1 = rc2ut(iE,iEo-1+shift_,NE)
         Endx2 = rc2ut(iE,iEo-2+shift_,NE)
         !
         xA = E(iEo-1-shift_); yA = K(Endx1)
         xB = E(iEo-2-shift_); yB = K(Endx2)
         Rey_W = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_W = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !the iE,iEo column is the east/west average
         Endx0 = rc2ut(iE,iEo,NE)
         K(Endx0) = dcmplx(( Rey_E + Rey_W ),( Imy_E + Imy_W ))/2d0
         !
      enddo
      !
      !linear fit of the (0,E)=(iEo,iE) row
      do iE=1,NE
         !
         if(iE.eq.iEo)cycle
         !
         !at each column iE - north direction
         Endx1 = rc2ut(iEo+1+shift_,iE,NE)
         Endx2 = rc2ut(iEo+2+shift_,iE,NE)
         xA = E(iEo+1+shift_); yA = K(Endx1)
         xB = E(iEo+2+shift_); yB = K(Endx2)
         Rey_N = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_N = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !at each column iE - south direction
         Endx1 = rc2ut(iEo-1+shift_,iE,NE)
         Endx2 = rc2ut(iEo-2+shift_,iE,NE)
         xA = E(iEo-1-shift_); yA = K(Endx1)
         xB = E(iEo-2-shift_); yB = K(Endx2)
         Rey_S = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
         Imy_S = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
         !
         !the iEo,iE row is the north/south average
         Endx0 = rc2ut(iEo,iE,NE)
         K(Endx0) = dcmplx(( Rey_N + Rey_S ),( Imy_N + Imy_S ))/2d0
         !
      enddo
      !
      !linear fit around the (0,0)=(iEo,iEo) point
      !east direction
      Endx1 = rc2ut(iEo,iEo+1+shift_,NE)
      Endx2 = rc2ut(iEo,iEo+2+shift_,NE)
      xA = E(iEo+1+shift_); yA = K(Endx1)
      xB = E(iEo+2+shift_); yB = K(Endx2)
      Rey_E = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_E = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !west direction
      Endx1 = rc2ut(iEo,iEo-1-shift_,NE)
      Endx2 = rc2ut(iEo,iEo-2-shift_,NE)
      xA = E(iEo-1-shift_); yA = K(Endx1)
      xB = E(iEo-2-shift_); yB = K(Endx2)
      Rey_W = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_W = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !north direction
      Endx1 = rc2ut(iEo+1+shift_,iEo,NE)
      Endx2 = rc2ut(iEo+2+shift_,iEo,NE)
      xA = E(iEo+1+shift_); yA = K(Endx1)
      xB = E(iEo+2+shift_); yB = K(Endx2)
      Rey_N = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_N = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !south direction
      Endx1 = rc2ut(iEo-1+shift_,iEo,NE)
      Endx2 = rc2ut(iEo-2+shift_,iEo,NE)
      xA = E(iEo-1-shift_); yA = K(Endx1)
      xB = E(iEo-2-shift_); yB = K(Endx2)
      Rey_S = linear_interp_2y([xA,dreal(yA)],[xB,dreal(yB)],0d0)
      Imy_S = linear_interp_2y([xA,dimag(yA)],[xB,dimag(yB)],0d0)
      !
      !the iEo,iEo point is the east/west/north/south average
      Endx0 = rc2ut(iEo,iEo,NE)
      K(Endx0) = dcmplx(( Rey_S + Rey_N + Rey_W + Rey_E ),( Imy_S + Imy_N + Imy_W + Imy_E ))/4d0
      !
   end subroutine interpFermi_K_tri_z


end module gap_equation