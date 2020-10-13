module input_vars

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
   integer                                  :: Nthread=28


   !---------------------------------------------------------------------------!
   !PURPOSE: Double counting types
   !---------------------------------------------------------------------------!
   character(len=10)                        :: VH_type="Ubare"
   logical                                  :: HandleGammaPoint


   !---------------------------------------------------------------------------!
   !PURPOSE: internal imaginary time mesh
   !---------------------------------------------------------------------------!
   logical                                  :: tau_uniform=.false.
   integer                                  :: NtauF=200
   integer                                  :: NtauB=1000
   integer                                  :: Nmats
   integer                                  :: Nreal=2000
   real(8)                                  :: wrealMax=10
   real(8)                                  :: eta=0.04
   real(8)                                  :: wmatsMax=100
   real(8)                                  :: Beta=15


   !---------------------------------------------------------------------------!
   !PURPOSE: paths. Directories must end with "/"
   !---------------------------------------------------------------------------!
   character(len=256)                       :: pathINPUT="InputFiles/"
   integer                                  :: LOGfile=6


   !---------------------------------------------------------------------------!
   !PURPOSE: density lookup
   !---------------------------------------------------------------------------!
   real(8)                                  :: densityPercErr=0.01


   !---------------------------------------------------------------------------!
   !PURPOSE: logical Flags
   !---------------------------------------------------------------------------!
   logical                                  :: UseXepsKorder=.true.
   logical                                  :: paramagneticSPEX=.true.
   logical                                  :: UfullStructure=.true.


   !---------------------------------------------------------------------------!
   !PURPOSE: Site and Orbital space
   !---------------------------------------------------------------------------!
   integer                                  :: Nsite=1
   integer,allocatable                      :: SiteOrbs(:,:)
   character(len=2),allocatable             :: SiteName(:)


   !---------------------------------------------------------------------------!
   !PURPOSE: K-points
   !---------------------------------------------------------------------------!
   integer                                  :: Nkpt3(3)=[8,8,8]


   !---------------------------------------------------------------------------!
   !PURPOSE: The most important variable
   !---------------------------------------------------------------------------!
   character(len=10)                        :: CalculationType="GW+EDMFT"




end module input_vars
