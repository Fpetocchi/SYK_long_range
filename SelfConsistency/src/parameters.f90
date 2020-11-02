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
   real(8),parameter                        :: H2eV=27.2113831243217d0
   real(8),parameter                        :: bohr=0.5291772108d0
   !
   integer,parameter                        :: spatial_directions=4 !100,010,001,110
   integer,parameter                        :: Nspin=2
   !
   real(8),parameter                        :: eps=1e-9
   !
   integer,parameter                        :: Nbath=10
   integer,parameter                        :: cg_niter=200
   real(8),parameter                        :: cg_Ftol=1e-6
   real(8),parameter                        :: hwband=3d0
   real(8),parameter                        :: noisefact=0.01


   !---------------------------------------------------------------------------!
   !PURPOSE: container for lattice data
   !---------------------------------------------------------------------------!
   type Lattice
      real(8),allocatable                   :: kpt(:,:)
      integer                               :: Nkpt3(3)
      complex(8),allocatable                :: Hk(:,:,:)                        ![Norb,Norb,Nkpt]
      complex(8),allocatable                :: Hloc(:,:)                        ![Norb,Norb]
      complex(8),allocatable                :: Zk(:,:,:)                        ![Norb,Norb,Nkpt]
      real(8),allocatable                   :: Ek(:,:)                          ![Norb,Nkpt]
      integer,allocatable                   :: kptPos(:)                        ![Nkpt]
      integer,allocatable                   :: kptsum(:,:)                      ![Nkpt,Nkpt]
      integer,allocatable                   :: kptdif(:,:)                      ![Nkpt,Nkpt]
      integer,allocatable                   :: kprint(:)
      integer,allocatable                   :: small_ik(:,:)                    ![12,2]
      integer                               :: Nkpt=0
      integer                               :: Nkpt_irred=0
      integer                               :: iq_gamma=-1
      integer                               :: Norb=0
      real(8)                               :: mu=0d0
      real(8)                               :: density=0d0
      logical                               :: UseDisentangledBS=.false.
      logical                               :: status=.false.
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
      real(8)                               :: TargetDensity=0d0
      real(8)                               :: densityRelErr=0d0
      real(8)                               :: muStep=0d0
      integer                               :: muIter=0
      real(8)                               :: muTime=0d0
      integer,allocatable                   :: orbs(:)                          !orbital restriction (has to be added to the input)
      logical                               :: quickloops=.false.
   end type musearch


end module parameters
