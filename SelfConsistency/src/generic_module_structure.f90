module module_name

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   !interface name
   !   module procedure name_a                                                  ! Description
   !   module procedure name_b                                                  ! Description
   !end interface name

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
   !variables
   !public ::
   !public ::
   !subroutines
   !public ::
   !public ::
   !functions
   !public ::
   !public ::

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: subroutine1_name does this
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine subroutine1_name()
      !
      !use mod_a
      !use mod_a
      implicit none
      !
      !input var 1
      !input var 2
      !
      !internal var 1
      !internal var 2
      !
      !
      write(*,"(A)") "--- subroutine1_name ---"
      !
      !
      ! Checks on the input

      !
   end subroutine subroutine1_name


   !---------------------------------------------------------------------------!
   !PURPOSE: subroutine2_name does this
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine subroutine2_name()
      !
      !use mod_a
      !use mod_a
      implicit none
      !
      !input var 1
      !input var 2
      !
      !internal var 1
      !internal var 2
      !
      !
      write(*,"(A)") "--- subroutine2_name ---"
      !
      !
      ! Checks on the input

      !
   end subroutine subroutine2_name


end module module_name
