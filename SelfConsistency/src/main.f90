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
   call read_InputFile("input.in")
   !
   write(LOGfile,"(A,1I4)") "Setting Nthread:",Nthread
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
   !Polarization
   if(calc_PiGG)then
      !
      call calc_Pi(PiGG,Glat,Crystal)
      call dump_BosonicField(PiGG,reg(pathDATA)//str(ItStart)//"/","PiGG.DAT")
      !
      if(merge_Pi) then
         call MergeFields(PiGG,PiEDMFT,alphaPi,SiteOrbs)
         call DeallocateBosonicField(PiEDMFT)
         call dump_BosonicField(PiGG,reg(pathDATA)//str(ItStart)//"/","PiGG_merged.DAT")
      endif
      !
   endif
   !
   !Fully screened interaction
   if(calc_W)then
      !
      if(calc_Wfull)  call calc_W_full(Wlat,Ulat,PiGG,Crystal)
      if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,PiEDMFT,Crystal)
      call dump_BosonicField(Wlat,reg(pathDATA)//str(ItStart)//"/","Wloc.DAT")
      !
      call DeallocateBosonicField(Ulat)
      !
   endif



   !














   !






   write(LOGfile,"(A,1F10.6)") "Self-Consistency finished. Total timing (s): ",tock(TimeStart)
   !
end program test
