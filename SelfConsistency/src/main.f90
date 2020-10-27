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
   calc_PiGG=.false.!TEST
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
   calc_W=.false.!TEST
   if(calc_W)then
      !
      if(calc_Wfull)  call calc_W_full(Wlat,Ulat,PiGG,Crystal)
      if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,PiEDMFT,Crystal)
      call dump_BosonicField(Wlat,reg(pathDATA)//str(ItStart)//"/","Wloc.DAT")
      !
   endif
   !
   !K-dependent self-energy
   if(calc_Sigma)then
      !
      !Hartree shift between G0W0 and LDA
      allocate(VH(Crystal%Norb,Crystal%Norb));VH=czero
      call calc_VH(densityLDA,Glat,Ulat,VH)
      call DeallocateBosonicField(Ulat)
      !
      !G0W0 contribution and Vexchange readed from SPEX
      allocate(Vxc(Crystal%Norb,Crystal%Norb,Crystal%Nkpt,Nspin));Vxc=czero
      call AllocateFermionicField(SigmaG0W0,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
      call read_Sigma_spex(SigmaG0W0,Crystal,verbose,Vxc=Vxc)
      !
      !scGW
      call AllocateFermionicField(SigmaGW_C,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
      call AllocateFermionicField(SigmaGW_X,Crystal%Norb,0,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
      call calc_sigmaGW(SigmaGW_C,SigmaGW_X,Glat,Wlat,Crystal)!call calc_sigmaGW(SigmaGW_C,SigmaGW_X,Glat,Ulat,Crystal)!TEST
      !
      !Dc between G0W0 and scGW
      if(ItStart.eq.0)then
         !
         call AllocateFermionicField(SigmaG0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call join_SigmaCX(SigmaG0W0dc,SigmaGW_C,SigmaGW_X)
         call dump_FermionicField(SigmaG0W0dc,reg(pathDATA)//"0/","SigmaG0W0",.true.,Crystal%kpt)
         call DeallocateFermionicField(SigmaG0W0dc)
         !
      elseif(ItStart.gt.0)then
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
      if(merge_Sigma.and.(ItStart.gt.0))then
         !
         !First compute the Dc
         call AllocateFermionicField(SigmaGWdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
         call AllocateFermionicField(SigmaGW_Cdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
         call AllocateFermionicField(SigmaGW_Xdc,Crystal%Norb,0,Nsite=Nsite,Beta=Beta)
         call calc_sigmaGWdc(SigmaGW_Cdc,SigmaGW_Xdc,Glat,Wlat)!call calc_sigmaGWdc(SigmaGW_Cdc,SigmaGW_Xdc,Glat,Ulat)!TEST
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
   call calc_SigmaFull(ItStart)
   !
   write(LOGfile,"(A,F)") "Self-Consistency finished. Total timing (s): ",tock(TimeStart)
   !
end program test
