module post_processing

   use gap_equation, only: store_Wk4gap
   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface dump_MaxEnt
      module procedure :: dump_MaxEnt_Gfunct                                    ![Array(Norb,Np,Nspin),Beta,mode,dirpath,filename,iorb(optinal)]
      module procedure :: dump_MaxEnt_Gfield                                    ![FermionicField,mode,dirpath,filename,orbset(:,:)]
      module procedure :: dump_MaxEnt_Wfunct                                    ![Array(Np),Beta,mode,dirpath,filename,iorb(2)]
      module procedure :: dump_MaxEnt_Wfield                                    ![BosonicField,mode,dirpath,filename,orbset(:,:)]
   end interface dump_MaxEnt

   interface interpolate2Beta
      module procedure :: interpolate2Beta_Fermionic                            ![FermionicField,OldBeta,offDiag]
      module procedure :: interpolate2Beta_Bosonic                              ![BosonicField,OldBeta,offDiag]
   end interface interpolate2Beta

   interface interpolate2kpath
      module procedure :: interpolate2kpath_Fermionic
      module procedure :: interpolate2kpath_Bosonic
   end interface interpolate2kpath

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !functions
   public :: pade
   !subroutines
   public :: dump_MaxEnt
   public :: remove_CDW
   public :: interpolate2Beta
   public :: interpolate2kpath
   public :: calc_Tc
   public :: store_Wk4gap

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE:
   ! - remove_CDW(W,mode,site)
   !---------------------------------------------------------------------------!
   include "post_processing/remove_CDW.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !---------------------------------------------------------------------------!
   include "post_processing/dump_MaxEnt.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !---------------------------------------------------------------------------!
   include "post_processing/calc_pade.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate to a new frequency mesh a Fermionic field
   !---------------------------------------------------------------------------!
   include "post_processing/interpolate_beta.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate to a user provided K-point path a Fermionic field
   !---------------------------------------------------------------------------!
   include "post_processing/interpolate_Kpath.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE: Estimate Tc by solving a gap equation
   !---------------------------------------------------------------------------!
   include "post_processing/gap_equation_interface.f90"


end module post_processing
