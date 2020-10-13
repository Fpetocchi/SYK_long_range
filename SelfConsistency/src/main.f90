program test
   !
   use module_container
   implicit none
   integer                                  :: CalculationStart
   integer                                  :: ItStart
   !
   type(Lattice)                            :: Crystal
   !
   type(FermionicField)                     :: Glat
   type(FermionicField)                     :: Slat
   !
   type(BosonicField)                       :: Wlat
   type(BosonicField)                       :: Ulat
   type(BosonicField)                       :: Plat
   !
   call tick(CalculationStart)
   !
   write(LOGfile,"(A)") new_line("A")//"Reading InputFile"//new_line("A")
   !call readInputFile(pathINPUT)
   Nmats = int(Beta*wmatsMax/(2d0*pi))
   !
   write(LOGfile,"(A)") "Setting Nthread"//new_line("A")
   call omp_set_num_threads(Nthread)
   !
   call printHeader()
   !
   call initialize_DataStructure(ItStart)
   !
   call initialize_Lattice(Crystal)
   !
   !
   select case(CalculationType)
      case default
         !
         stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
         !
      case("G0W0","scGW")
         !
         write(LOGfile,"(A)")
         !
      case("DMFT+statU")
         !
         write(LOGfile,"(A)")
         !
      case("DMFT+dynU")
         !
         write(LOGfile,"(A)")
         !
      case("EDMFT")
         !
         write(LOGfile,"(A)")
         !
      case("GW+EDMFT")
         !
         call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
         call AllocateFermionicField(Slat,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
         call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
         call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
         call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%Nkpt,Nsite)
         !
   end select









   !






   write(LOGfile,"(A,1F10.6)") "Self-Consistency finished. Total timing (s): ",tock(CalculationStart)
   !
end program test
