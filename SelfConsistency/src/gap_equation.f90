module gap_equation

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
   !PURPOSE: Module custom types
   !---------------------------------------------------------------------------!
   !type mytype
   !   integer                               :: var1
   !   real(8),allocatable                   :: var2
   !   logical                               :: var3
   !end type mytype

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),allocatable,private              :: omega(:)                         ! Phonon energy on logarithmic grid
   real(8),allocatable,private              :: a2F(:)                           ! alpha^2*F(\Omega) function
   logical,private                          :: Phonons_stored=.false.
   !
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   public :: calc_Tc

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the output of phonon calculations with  Elk or QuantumEspresso
   !---------------------------------------------------------------------------!
   subroutine read_a2F(pathINPUT,mode)
      !
      use parameters
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=*),intent(in)           :: mode
      !
      integer                               :: iomega,Nomega
      real(8)                               :: ConversionFactor
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_a2F"
      !
      !
      select case(reg(mode))
         case default
            !
            stop "Available phonon inputs: Elk, QEspresso."
            !
         case("Elk")
            !
            ConversionFactor = H2eV
            !
         case("QEspresso")
            !
            ConversionFactor = Ry2H*H2eV
            !
      end select
      !
      !reading phonons form file
      call inquireFile(reg(pathINPUT)//"ALPHA2F.DAT",filexists)
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"ALPHA2F.DAT",form="formatted",action="read",position="rewind")
      read(unit,'("# nw alpha2F: ",I8)') Nomega
      !
      allocate(omega(Nomega));omega=0d0
      allocate(a2F(Nomega));a2F=0d0
      !
      read(unit,*)
      do iomega=1,Nomega
        read(unit,'(2G18.10)') omega(iomega),a2F(iomega)
      enddo
      close(unit)
      !
      a2F = a2F * ConversionFactor
      !
      Phonons_stored=.true.
      !
   end subroutine read_a2F


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the phononic renormalization factor on an energy grid
   !---------------------------------------------------------------------------!
   subroutine calc_Zph_e(beta,Egrid,DoS,Zph_e,mode)
      !
      use utils_misc
      implicit none
      !
      !input var 1
      !input var 2
      !
      !internal var 1
      !internal var 2
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Zph_e"
      !
      !
      ! Checks on the input

      !
   end subroutine calc_Zph_e


end module gap_equation
