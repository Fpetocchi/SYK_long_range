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


   !---------------------------------------------------------------------------!
   !PURPOSE: container for lattice data
   !---------------------------------------------------------------------------!
   type lattice
      real(8),allocatable                   :: kpt(:,:)
      integer                               :: Nkpt3(3)
      complex(8),allocatable                :: Hk(:,:,:)                        ![Norb,Norb,Nkpt]
      complex(8),allocatable                :: Zk(:,:,:)                        ![Norb,Norb,Nkpt]
      real(8),allocatable                   :: Ek(:,:)                          ![Norb,Nkpt]
      integer,allocatable                   :: kptsum(:,:)
      integer,allocatable                   :: kptdif(:,:)
      integer                               :: small_ik(12,2)
      integer                               :: Nkpt=0
      integer                               :: Norb=0
      logical                               :: status=.false.
   end type lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: k-dependent and local matrices structures - axis independent
   !---------------------------------------------------------------------------!
   type FermionicField
      complex(8),allocatable                :: wk(:,:,:,:,:)                    ![Norb,Norb,Npoints,Nkpt,Nspin]
      complex(8),allocatable                :: w(:,:,:,:)                       ![Norb,Norb,Npoints,Nspin]
      integer                               :: Norb=0
      integer                               :: Npoints=0
      integer                               :: Nkpt=0
      integer                               :: Nsite=1
      real(8)                               :: Beta=0d0
      real(8)                               :: mu=0d0
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
      integer                               :: Nsite=1
      real(8)                               :: Beta=0d0
      logical                               :: status=.false.
   end type BosonicField


end module parameters
