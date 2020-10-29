program test
   !
   use module_container
   use utils_main
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: Iteration,ItStart,Itend
   !
   !
   call tick(TimeStart)
   call read_InputFile("input.in")
   write(LOGfile,"(A,1I4)") "Setting Nthread:",Nthread
   call omp_set_num_threads(Nthread)
   call printHeader()
   call initialize_DataStructure(ItStart,Itend)
   call initialize_Lattice(Crystal,ItStart)
   call initialize_Fields(ItStart)
   !



   do Iteration=ItStart,Itend,1
      !
      call printHeader(Iteration)
      !
      !Polarization
      if(calc_PiGG)then
         !
         call calc_Pi(PiGG,Glat,Crystal)
         call dump_BosonicField(PiGG,reg(pathDATA)//str(Iteration)//"/","PiGG.DAT")
         !
         if(merge_Pi.and.solve_DMFT) then !a bit redundant since there is no merge wihtout DMFT
            call MergeFields(PiGG,PiEDMFT,alphaPi,SiteOrbs)
            call DeallocateBosonicField(PiEDMFT)
            call dump_BosonicField(PiGG,reg(pathDATA)//str(Iteration)//"/","PiGG_merged.DAT")
         endif
         !
      endif
      !
      !Fully screened interaction
      if(calc_W)then
         !
         if(calc_Wfull)  call calc_W_full(Wlat,Ulat,PiGG,Crystal)
         if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,PiEDMFT,Crystal)
         call dump_BosonicField(Wlat,reg(pathDATA)//str(Iteration)//"/","Wloc.DAT")
         !
      endif
      !
      !K-dependent self-energy
      if(calc_Sigma)then
         !
         !Hartree shift between G0W0 and LDA
         allocate(VH(Crystal%Norb,Crystal%Norb));VH=czero
         call calc_VH(densityLDA,Glat,Ulat,VH)
         if(solve_DMFT)call DeallocateBosonicField(Ulat)
         !
         !G0W0 contribution and Vexchange readed from SPEX
         allocate(Vxc(Crystal%Norb,Crystal%Norb,Crystal%Nkpt,Nspin));Vxc=czero
         call AllocateFermionicField(SigmaG0W0,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_Sigma_spex(SigmaG0W0,Crystal,verbose,Vxc=Vxc)
         !
         !scGW
         call AllocateFermionicField(SigmaGW_C,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call AllocateFermionicField(SigmaGW_X,Crystal%Norb,0,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call calc_sigmaGW(SigmaGW_C,SigmaGW_X,Glat,Wlat,Crystal)
         !
         !Dc between G0W0 and scGW
         if(Iteration.eq.0)then
            !
            call AllocateFermionicField(SigmaG0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call join_SigmaCX(SigmaG0W0dc,SigmaGW_C,SigmaGW_X)
            call dump_FermionicField(SigmaG0W0dc,reg(pathDATA)//"0/","SigmaG0W0",.true.,Crystal%kpt)
            call DeallocateFermionicField(SigmaG0W0dc)
            !
         elseif(Iteration.gt.0)then
            !
            call AllocateFermionicField(SigmaG0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call read_FermionicField(SigmaG0W0dc,reg(pathDATA)//"0/","SigmaG0W0",kpt=Crystal%kpt)
            !
            call AllocateFermionicField(SigmaGW,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call join_SigmaCX(SigmaGW,SigmaGW_C,SigmaGW_X)
            !
         endif
         call DeallocateFermionicField(SigmaGW_C)
         call DeallocateFermionicField(SigmaGW_X)
         !
         !Merge GW and EDMFT
         if((Iteration.gt.0).and.merge_Sigma.and.solve_DMFT)then !a bit redundant since there is no merge wihtout DMFT
            !
            !First compute the Dc
            call AllocateFermionicField(SigmaGWdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            call AllocateFermionicField(SigmaGW_Cdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            call AllocateFermionicField(SigmaGW_Xdc,Crystal%Norb,0,Nsite=Nsite,Beta=Beta)
            call calc_sigmaGWdc(SigmaGW_Cdc,SigmaGW_Xdc,Glat,Wlat)
            call join_SigmaCX(SigmaGWdc,SigmaGW_Cdc,SigmaGW_Xdc)
            call DeallocateFermionicField(SigmaGW_Cdc)
            call DeallocateFermionicField(SigmaGW_Xdc)
            !
            !Then replace local projection with scGW and EDMFT
            call MergeFields(SigmaGW,SigmaGWdc,SigmaDMFT,alphaSigma,SiteOrbs,DC_type)
            call DeallocateFermionicField(SigmaGWdc)
            call DeallocateFermionicField(SigmaDMFT)
            !
         endif
         !
      endif
      !
      !Put together all the contributions to the full self-energy
      call calc_SigmaFull(Iteration)
      !
      !Compute the Full Green's function
      call calc_Gmats(Glat,Crystal,SigmaFull)
      !
   enddo







   write(LOGfile,"(A,F)") "Self-Consistency finished. Total timing (s): ",tock(TimeStart)
   !
end program test
