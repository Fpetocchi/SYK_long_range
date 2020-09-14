module global_vars

   use parameters
   implicit none


   ! COMMENTS:
   !
   !


   !---------------------------------------------------------------------------!
   !PURPOSE: Parallelization integers
   !---------------------------------------------------------------------------!
   type MpiEnv
      integer                               :: MpiComm
      integer                               :: MpiErr
      integer                               :: MpiRank
      integer                               :: MpiSize
      logical                               :: status=.false.
   end type MpiEnv
   !
   integer                                  :: Nthread


   !---------------------------------------------------------------------------!
   !PURPOSE: Frequency and time Meshes
   !---------------------------------------------------------------------------!
   integer                                  :: Nw_F,Nw_B
   integer                                  :: Ntau_F,Ntau_B
   real(8)                                  :: Beta
   real(8),allocatable                      :: tau_F(:),tau_B(:)
   real(8),allocatable                      :: wm_F(:),wm_B(:)


   !---------------------------------------------------------------------------!
   !PURPOSE: K-points
   !---------------------------------------------------------------------------!
   integer                                  :: Nkpt
   integer                                  :: Nkpt3(3)
   real(8),allocatable                      :: kpt(:,:)                         ![(kx,ky,kz),ik]
   integer                                  :: Nkpt_irred
   integer,allocatable                      :: kptsum(:,:),kptdiff(:,:)


   !---------------------------------------------------------------------------!
   !PURPOSE: Orbital spaces
   !---------------------------------------------------------------------------!
   integer                                  :: Nsite
   integer                                  :: NorbGW
   integer                                  :: NpbGW


   !---------------------------------------------------------------------------!
   !PURPOSE: Hamiltonian Parameters
   !---------------------------------------------------------------------------!
   real(8)                                  :: Uaa
   real(8)                                  :: Uab
   real(8)                                  :: Jh,Jsf,Jph
   real(8)                                  :: ChemPot


   !---------------------------------------------------------------------------!
   !PURPOSE: Observables
   !---------------------------------------------------------------------------!
   real(8)                                  :: DensTarget
   real(8)                                  :: DensTot
   real(8),allocatable                      :: DensOrb(:,:)


   !---------------------------------------------------------------------------!
   !PURPOSE: Fermionic fields - only on Matsubara
   !---------------------------------------------------------------------------!
   type(FermionicField)                     :: Gmats
   type(FermionicField)                     :: SigmaGW
   type(FermionicField)                     :: SigmaDMFT


   !---------------------------------------------------------------------------!
   !PURPOSE: Bosonic fields - only on Matsubara
   !---------------------------------------------------------------------------!
   type(BosonicField)                       :: Ucrpa
   type(BosonicField)                       :: PiGG
   type(BosonicField)                       :: Pi_DMFT


   !---------------------------------------------------------------------------!
   !PURPOSE: Hamiltonians, Rotations and LDA data
   !---------------------------------------------------------------------------!
   real(8),allocatable                      :: Ek(:,:)
   complex(8),allocatable                   :: Hk(:,:,:)
   complex(8),allocatable                   :: HartreeShift(:,:)
   complex(8),allocatable                   :: Vxc(:,:,:)


   !---------------------------------------------------------------------------!
   !PURPOSE: paths
   !---------------------------------------------------------------------------!
   character(len=256)                       :: pathINPUT="GWinput/"


   !---------------------------------------------------------------------------!
   !PURPOSE: logical Flags
   !---------------------------------------------------------------------------!
   logical                                  :: UseXepsKorder
   logical                                  :: UfullStructure








end module global_vars
