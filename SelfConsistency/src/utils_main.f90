module utils_main

   use module_container
   implicit none
   public

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
   type(Lattice)                            :: Crystal
   type(FermionicField)                     :: Glat
   type(FermionicField)                     :: Gimp
   type(FermionicField)                     :: SigmaFull
   type(FermionicField)                     :: SigmaG0W0
   type(FermionicField)                     :: SigmaGW_C,SigmaGW_X
   type(FermionicField)                     :: SigmaGW_Cdc,SigmaGW_Xdc
   type(FermionicField)                     :: SigmaDMFT
   !
   type(BosonicField)                       :: Wlat
   type(BosonicField)                       :: Urpa
   type(BosonicField)                       :: PiGG
   type(BosonicField)                       :: PiEDMFT
   type(BosonicField)                       :: curlyU

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: printHeader
   public :: initialize_DataStructure
   public :: initialize_Lattice
   public :: initialize_Interactions

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: prints the header
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine printHeader()
      !
      implicit none
      character(len=80)                     :: line
      character(:),allocatable              :: header
      integer                               :: i,Lsx,Ldx,Lcalc
      !
      header="*"
      line="*"
      do i=1,79
         line = trim(line)//"*"
      enddo
      Lcalc=len("Calculation type: "//trim(CalculationType))
      Lsx=int((78-Lcalc)/2)-1
      Ldx=78-Lcalc-Lsx
      do i=1,Lsx
         header = header//" "
      enddo
      header = header//"Calculation type: "//trim(CalculationType)
      do i=1,Ldx
         header = header//" "
      enddo
      header = header//"*"
      write(LOGfile,"(A)") line
      write(LOGfile,"(A)") header
      write(LOGfile,"(A)") line//new_line("A")//new_line("A")
      !
   end subroutine printHeader


   !---------------------------------------------------------------------------!
   !PURPOSE: looks for the current iteration number
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine initialize_DataStructure(ItStart)
      !
      implicit none
      integer,intent(out)                   :: ItStart
      character(len=256)                    :: Itpath
      integer                               :: iter,Itfirst
      logical                               :: Itexist
      !
      do iter=0,1000
         Itpath = "Iteration/"//str(iter)
         call inquireDir(reg(Itpath),Itexist,hardstop=.false.)
         if(.not.Itexist)then
            Itfirst = iter
            Exit
         endif
      enddo
      Itpath = "Iteration/"//str(Itfirst)
      !
      if(Itfirst.eq.0)then
         write(LOGfile,"(A)") "Brand new calculation. Initializing "//reg(Itpath)
      else
         write(LOGfile,"(A)") "Last iteration: "//str(Itfirst-1)//". Initializing "//reg(Itpath)
      endif
      call createDir(reg(Itpath))
      !
      ItStart = Itfirst
      !
   end subroutine initialize_DataStructure


   !---------------------------------------------------------------------------!
   !PURPOSE: Initialize Lattice. I could have used the AllocateLattice in
   !         utils_fields but then als useless attributes would have been allocated
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine initialize_Lattice(Lttc,ItStart)
      !
      implicit none
      type(Lattice),intent(out)             :: Lttc
      integer,intent(in)                    :: ItStart
      !
      !
      write(LOGfile,"(A)") "--- initialize_Lattice ---"
      !
      !
      select case(CalculationType)
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0","scGW")
            !
            call read_Hk(pathINPUT,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc)
            !
            Lttc%Norb = size(Lttc%Hk,dim=1)
            Lttc%Nkpt = size(Lttc%Hk,dim=3)
            !
            allocate(Lttc%kptPos(Lttc%Nkpt));Lttc%kptPos=0
            call read_xeps(pathINPUT,Lttc%kpt,Nkpt3,UseXepsKorder, &
            Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,paramagneticSPEX)
            !
            allocate(Lttc%kptsum(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptsum=0
            allocate(Lttc%kptdif(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptdif=0
            call fill_ksumkdiff(Lttc%kpt,Lttc%kptsum,Lttc%kptdif,Nkpt3)
            !
            allocate(Lttc%small_ik(12,2));Lttc%small_ik=0
            call fill_smallk(Lttc%kpt,Lttc%small_ik)
            !
            Lttc%status=.true.
            !
         case("DMFT+statU","DMFT+dynU")
            !
            call read_Hk(pathINPUT,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc)
            !
            Lttc%Norb = size(Lttc%Hk,dim=1)
            Lttc%Nkpt = size(Lttc%Hk,dim=3)
            !
            Lttc%status=.true.
            !
         case("EDMFT","GW+EDMFT")
            !
            call read_Hk(pathINPUT,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc)
            !
            Lttc%Norb = size(Lttc%Hk,dim=1)
            Lttc%Nkpt = size(Lttc%Hk,dim=3)
            !
            allocate(Lttc%kptPos(Lttc%Nkpt));Lttc%kptPos=0
            call read_xeps(pathINPUT,Lttc%kpt,Nkpt3,UseXepsKorder, &
            Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,paramagneticSPEX)
            !
            allocate(Lttc%kptsum(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptsum=0
            allocate(Lttc%kptdif(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptdif=0
            call fill_ksumkdiff(Lttc%kpt,Lttc%kptsum,Lttc%kptdif,Nkpt3)
            !
            allocate(Lttc%small_ik(12,2));Lttc%small_ik=0
            call fill_smallk(Lttc%kpt,Lttc%small_ik)
            !
            Lttc%status=.true.
            !
      end select
      !
      if(ItStart.eq.0)call calc_Glda(0d0,Beta,Lttc)
      !
   end subroutine initialize_Lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Initialize the Fields depending on the starting iteration
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine initialize_Fields(ItStart)
      !
      implicit none
      type(Lattice),intent(out)             :: Lttc
      !
      !
      write(LOGfile,"(A)") "--- initialize_Fields ---"
      !
      !
      select case(CalculationType)
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0","scGW")
            !
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
            !
         case("DMFT+statU")
            !
            call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)




            allocate(Uloc(Nspin*Crystal%Norb,Nspin*Crystal%Norb));Uloc=0d0
            !
         case("DMFT+dynU")
            !
            call AllocateBosonicField(Urpa,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
            !
         case("EDMFT")
            !
            call AllocateBosonicField(Urpa,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
            !
         case("GW+EDMFT")
            !
            call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
            call AllocateBosonicField(Urpa,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
            call AllocateBosonicField(PiGG,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
            !
      end select
      !
   end subroutine initialize_Interactions




end module utils_main
