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
   integer                                  :: Nthread


   !---------------------------------------------------------------------------!
   !PURPOSE: Double counting types
   !---------------------------------------------------------------------------!
   character                                :: VH_type="Ubare"


   !---------------------------------------------------------------------------!
   !PURPOSE: K-points
   !---------------------------------------------------------------------------!
   integer                                  :: Nkpt3(3)


   !---------------------------------------------------------------------------!
   !PURPOSE: internal imaginary time mesh
   !---------------------------------------------------------------------------!
   logical                                  :: tau_uniform
   integer                                  :: Ntau


   !---------------------------------------------------------------------------!
   !PURPOSE: paths
   !---------------------------------------------------------------------------!
   character(len=256)                       :: pathINPUT="GWinput/"


   !---------------------------------------------------------------------------!
   !PURPOSE: density lookup
   !---------------------------------------------------------------------------!
   real(8)                                  :: densityPercErr=0.01


   !---------------------------------------------------------------------------!
   !PURPOSE: logical Flags
   !---------------------------------------------------------------------------!
   logical                                  :: UseXepsKorder
   logical                                  :: paramagneticSPEX
   logical                                  :: UfullStructure



end module input_vars
