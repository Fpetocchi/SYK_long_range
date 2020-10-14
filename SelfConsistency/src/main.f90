program test
   !
   use module_container
   use utils_main
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: ItStart
   !
   !
   call tick(TimeStart)
   !
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
   call initialize_Lattice(Crystal,ItStart)
   !
   call initialize_Fields(ItStart)
   !










   !






   write(LOGfile,"(A,1F10.6)") "Self-Consistency finished. Total timing (s): ",tock(TimeStart)
   !
end program test
