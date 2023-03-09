module parameters

   implicit none

   !===========================================================================!

   ! COMMENTS:
   !
   !

   integer,parameter                        :: size_i=4
   integer,parameter                        :: size_d=8
   integer,parameter                        :: size_dc=16
   !
   real(8),parameter                        :: pi=3.14159265358979323846d0
   complex(8),parameter                     :: img=dcmplx(0.d0,1.d0)
   complex(8),parameter                     :: cone=dcmplx(1.d0,0.d0)
   complex(8),parameter                     :: czero=dcmplx(0.d0,0.d0)
   !
   real(8),parameter                        :: bohr=0.5291772108d0
   real(8),parameter                        :: K2eV=8.617333262d-5
   !
   real(8),parameter                        :: H2eV=27.2113831243217d0
   real(8),parameter                        :: Ry2H=0.5d0
   real(8),parameter                        :: Ry2eV=Ry2H*H2eV
   !
   real(8),parameter                        :: eV2H=1d0/H2eV
   real(8),parameter                        :: H2Ry=1d0/Ry2H
   real(8),parameter                        :: eV2Ry=1d0/Ry2eV
   !
   integer,parameter                        :: spatial_directions=4 !100,010,001,110
   integer,parameter                        :: Nspin=2
   !
   real(8),parameter                        :: eps=1e-9


   !---------------------------------------------------------------------------!
   !PURPOSE: container for site-orbitals arrangement
   !---------------------------------------------------------------------------!
   type LocalOrbitals
      integer                               :: Norb
      integer                               :: Nflavor
      character(len=2)                      :: Name
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: CrystalField(:)
      real(8),allocatable                   :: tailFit(:,:)
      !local rotations
      complex(8),allocatable                :: Op(:,:)
      complex(8),allocatable                :: Rot(:,:)
      complex(8),allocatable                :: RotDag(:,:)
      real(8),allocatable                   :: Eig(:)
      !local observables
      real(8),allocatable                   :: rho_Flav(:)                      ![Norb*Nspin]
      real(8),allocatable                   :: rho_OrbSpin(:,:,:)               ![Norb,Norb,Nspin]
      real(8),allocatable                   :: Docc(:,:)
   end type LocalOrbitals


   !---------------------------------------------------------------------------!
   !PURPOSE: container for lattice data
   !---------------------------------------------------------------------------!
   type Lattice
      !generic
      integer                               :: Nsite
      integer                               :: Nkpt3(3)
      integer                               :: Nkpt_irred=0
      integer                               :: iq_gamma=-1
      complex(8),allocatable                :: Hloc(:,:)                        ![Norb,Norb]
      integer,allocatable                   :: kptPos(:)                        ![Nkpt]
      integer,allocatable                   :: kptsum(:,:)                      ![Nkpt,Nkpt]
      integer,allocatable                   :: kptdif(:,:)                      ![Nkpt,Nkpt]
      integer                               :: Norb=0
      real(8)                               :: mu=0d0
      real(8)                               :: density=0d0
      real(8)                               :: D_lower=0d0
      real(8)                               :: D_upper=0d0
      logical                               :: UseDisentangledBS=.false.
      logical                               :: status=.false.
      !full BZ
      integer                               :: Nkpt=0
      real(8),allocatable                   :: kpt(:,:)
      complex(8),allocatable                :: Hk(:,:,:)                        ![Norb,Norb,Nkpt]
      complex(8),allocatable                :: Zk(:,:,:)                        ![Norb,Norb,Nkpt]
      real(8),allocatable                   :: Ek(:,:)                          ![Norb,Nkpt]
      !path along high-symmetry points
      integer                               :: Nkpt_path=0
      real(8),allocatable                   :: kptpath(:,:)                     ![3,Nkpt_path]
      real(8),allocatable                   :: Kpathaxis(:)
      real(8),allocatable                   :: KpathaxisPoints(:)
      complex(8),allocatable                :: Hk_path(:,:,:)                   ![Norb,Norb,Nkpt_path]
      complex(8),allocatable                :: Zk_path(:,:,:)                   ![Norb,Norb,Nkpt_path]
      real(8),allocatable                   :: Ek_path(:,:)                     ![Norb,Nkpt_path]
      logical                               :: pathStored=.false.
      !path along planar sheet on kx,ky
      integer                               :: Nkpt_plane=0
      real(8),allocatable                   :: kptPlane(:,:)                    ![3,Nkpt_path]
      complex(8),allocatable                :: Hk_plane(:,:,:)                  ![Norb,Norb,Nkpt_plane]
      complex(8),allocatable                :: Zk_plane(:,:,:)                  ![Norb,Norb,Nkpt_plane]
      real(8),allocatable                   :: Ek_plane(:,:)                    ![Norb,Nkpt_plane]
      logical                               :: planeStored=.false.
      !hamiltonian in the full BZ with G0W0 self-energy correction
      complex(8),allocatable                :: Hk_qp(:,:,:,:)                   ![Norb,Norb,Nkpt,Nspin]
      real(8),allocatable                   :: Ek_qp(:,:,:)                     ![Norb,Nkpt,Nspin]
   end type Lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: k-dependent and local matrices structures - axis independent
   !---------------------------------------------------------------------------!
   type FermionicField
      complex(8),allocatable                :: wks(:,:,:,:,:)                   ![Norb,Norb,Npoints,Nkpt,Nspin]
      complex(8),allocatable                :: ws(:,:,:,:)                      ![Norb,Norb,Npoints,Nspin]
      complex(8),allocatable                :: N_ks(:,:,:,:)                    ![Norb,Norb,Nkpt,Nspin]
      complex(8),allocatable                :: N_s(:,:,:)                       ![Norb,Norb,Nspin]
      integer                               :: Norb=0
      integer                               :: Npoints=0
      integer                               :: Nkpt=0
      integer                               :: Nsite=0
      real(8)                               :: Beta=0d0
      real(8)                               :: mu=0d0
      logical                               :: local_filled=.false.
      logical                               :: status=.false.
   end type FermionicField


   !---------------------------------------------------------------------------!
   !PURPOSE: k-dependent and local tensor structures - axis independent
   !---------------------------------------------------------------------------!
   type BosonicField
      complex(8),allocatable                :: screened(:,:,:,:)                ![Nbp,Nbp,Npoints,Nkpt]
      complex(8),allocatable                :: bare(:,:,:)                      ![Nbp,Nbp,Nkpt]
      complex(8),allocatable                :: screened_local(:,:,:)            ![Nbp,Nbp,Npoints]
      complex(8),allocatable                :: bare_local(:,:)                  ![Nbp,Nbp]
      integer                               :: Nbp=0
      integer                               :: Npoints=0
      integer                               :: Nkpt=0
      integer                               :: Nsite=0
      integer                               :: iq_gamma=-1
      real(8)                               :: Beta=0d0
      logical                               :: local_filled=.false.
      logical                               :: status=.false.
   end type BosonicField


   !---------------------------------------------------------------------------!
   !PURPOSE: container for density lookup parameters
   !---------------------------------------------------------------------------!
   type musearch
      real(8)                               :: mu=0d0
      real(8)                               :: TargetDensity=0d0
      real(8)                               :: densityRelErr=0d0
      real(8)                               :: muStep=0d0
      integer                               :: muIter=0
      real(8)                               :: muTime=0d0
      logical                               :: local=.false.
      integer,allocatable                   :: orbs(:)                          !orbital restriction (has to be added to the input)
      !I'm writing this as integers due to some mismatch on how to write boolean between fortrann and c++
      integer                               :: mu_scan=0
   end type musearch


   !---------------------------------------------------------------------------!
   !PURPOSE: container for symmetrization variables
   !---------------------------------------------------------------------------!
   type Equivalent
      real(8)                               :: hseed=0d0
      integer                               :: Nset=0                           !set read from the input
      integer                               :: Ntotset=0                        !set actually used
      integer,allocatable                   :: SetNorb(:)
      integer,allocatable                   :: SetOrbs(:,:)
      logical                               :: Gfoffdiag=.true.
      logical                               :: O=.false.
      !I'm writing this as integers due to some mismatch on how to write boolean between fortrann and c++
      integer                               :: para=1
   end type Equivalent


   !---------------------------------------------------------------------------!
   !PURPOSE: container to store the physical interaction elements
   !---------------------------------------------------------------------------!
   type physicalU
      !Interaction with size Nspin*Norb
      logical,allocatable                   :: Flav_Uloc(:,:)
      logical,allocatable                   :: Flav_U1st(:,:)
      logical,allocatable                   :: Flav_U2nd(:,:)
      logical,allocatable                   :: Flav_All(:,:)
      integer,allocatable                   :: Flav_Map(:,:,:)
      integer                               :: Flav_Size
      !Interaction with size Norb^2
      logical,allocatable                   :: Full_Uaa(:,:)
      logical,allocatable                   :: Full_Uab(:,:)
      logical,allocatable                   :: Full_Jsf(:,:)
      logical,allocatable                   :: Full_Jph(:,:)
      logical,allocatable                   :: Full_Imp(:,:)
      logical,allocatable                   :: Full_All(:,:)
      integer,allocatable                   :: Full_Map(:,:,:)
      integer                               :: Full_Size
      !
      logical                               :: status=.false.
   end type physicalU


   !---------------------------------------------------------------------------!
   !PURPOSE: container for QMC Solver variables
   !---------------------------------------------------------------------------!
   type QMC
      integer,allocatable                   :: Time(:)
      real(8)                               :: TargetDensity=0d0
      integer                               :: Nimp
      integer                               :: NtauF
      integer                               :: NtauF_D
      integer                               :: NtauB
      integer                               :: NtauB_K
      integer                               :: NtauF_in
      integer                               :: NtauB_in
      integer                               :: Norder=0
      integer                               :: Nmeas=0
      integer                               :: Gexp=0
      integer                               :: Ntherm=0
      integer                               :: Nshift=0
      integer                               :: Nswap=0
      integer                               :: N_nnt=0
      integer                               :: PrintTime=0
      integer                               :: binlength=0
      integer                               :: binstart=0
      !I'm writing these as integers due to some mismatch on how to write boolean between fortrann and c++
      integer                               :: retarded=0
      integer                               :: mu_scan=0
      integer                               :: removeUhalf=0
      integer                               :: Imprvd_F=0
      integer                               :: Imprvd_B=0
      integer                               :: full_ntOrbSym=0
      integer                               :: tau_uniform_D=1
      integer                               :: tau_uniform_K=1
   end type QMC

   !---------------------------------------------------------------------------!
   !PURPOSE: container for the match beta variables
   !---------------------------------------------------------------------------!
   type OldBeta
      logical                               :: status=.false.
      real(8)                               :: Beta_old=0d0
      real(8)                               :: Beta_new=0d0
      character(len=256)                    :: Path
      real(8)                               :: wmatsMax=0d0
      integer                               :: Nmats_old=0
      integer                               :: Nmats_new=0
   end type OldBeta

   !---------------------------------------------------------------------------!
   !PURPOSE: container for the gap equation variables
   !---------------------------------------------------------------------------!
   type SCDFT
      real(8)                               :: Tbounds(2)=0d0
      integer                               :: Tsteps=0
      integer                               :: loops=0
      real(8)                               :: DeltaErr=0d0
      real(8)                               :: DeltaInit=0d0
      real(8)                               :: DeltaMix=0d0
      logical                               :: HkRenorm
      character(len=255)                    :: mode_ph                          !Elk or QEspresso
      character(len=255)                    :: mode_Zph="symrenorm"
      character(len=255)                    :: mode_el                          !static or static+dynamic
      integer                               :: Nkpt3_Model(3)=0
      real(8)                               :: wstep=0d0
      real(8)                               :: Wk_cutoff=0d0
      character(len=255)                    :: printmode_ph
      character(len=255)                    :: printmode_el
      logical                               :: calc_Tc=.false.
      logical                               :: status=.false.
   end type SCDFT

   !---------------------------------------------------------------------------!
   !PURPOSE: container for the heterostructure variables
   !---------------------------------------------------------------------------!
   type Heterostructures
      integer                               :: Nslab=1
      integer                               :: Nlayer=1
      integer                               :: Norb=0
      integer                               :: Explicit(2)=[1,1]
      integer                               :: tzIndex(2)=[0,0]
      integer                               :: tzRange=0
      real(8),allocatable                   :: tz(:,:,:)                        ![Norb,Nlayer,tzRange]
      complex(8),allocatable                :: tkz(:,:,:,:)                     ![Norb,Norb,Nkpt,Nlayer]
      complex(8),allocatable                :: tkz_path(:,:,:,:)                ![Norb,Norb,Nkpt,Nlayer]
      complex(8),allocatable                :: tkz_plane(:,:,:,:)               ![Norb,Norb,Nkpt,Nlayer]
      complex(8),allocatable                :: P_L(:,:,:,:)                     ![Norb,Norb,Npoints,Nspin]
      complex(8),allocatable                :: P_R(:,:,:,:)                     ![Norb,Norb,Npoints,Nspin]
      logical                               :: status=.false.
      logical                               :: offDiagEk=.false.
   end type Heterostructures


end module parameters
