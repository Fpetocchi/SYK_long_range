program test
   !
   ! COMMENTS:
   ! tengo collect_QMC_results(ItStart) separato perche forse in un futuro riscrivero il solver
   ! e quindi potro semplicemnte sponstare collect_QMC_results(ItStart) in fondo e tenere
   ! initialize_Fields(ItStart) invariato, Se metto dyson dentro initialize_Fields
   ! dopo se cambio solver mi tocca riscriverlo.
   ! quindi tengo separato collect_QMC_results che crea e scrive le self-energie e Pimp
   ! oi dealloco e rileggo


   use module_container
   use utils_main
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: isite
   integer                                  :: Iteration,ItStart,Itend
   !
   !
   !
   !---------------------------------------------------------------------------!
   !      READING INPUT FILE, INITIALIZING OMP, AND CHECK DATA STRUCTURE       !
   !---------------------------------------------------------------------------!
   call tick(TimeStart)
   call read_InputFile("input.in")
   write(LOGfile,"(A,1I4)") "Setting Nthread:",Nthread
   call omp_set_num_threads(Nthread)
   call printHeader()
   call initialize_DataStructure(ItStart,Itend)
   call initialize_Lattice(Crystal,ItStart)
   !
   !
   !
   !---------------------------------------------------------------------------!
   !       COLLECTING RESULTS FROM THE SOLVER AND SOLVING DYSON EQUATION       !
   !---------------------------------------------------------------------------!
   if(solve_DMFT)call collect_QMC_results()
   !
   !
   !
   !---------------------------------------------------------------------------!
   !SOLVING THE LATTICE PROBLEM AND PRODUCING INPUTS FOR NEXT IMPURITY SOLUTION!
   !---------------------------------------------------------------------------!
   call initialize_Fields(ItStart)
   !
   do Iteration=ItStart,Itend,1
      !
      call printHeader(Iteration)
      !
      !K-dependent Polarization
      if(calc_Plat)then
         !
         call calc_Pi(Plat,Glat,Crystal)
         call dump_BosonicField(Plat,reg(ItFolder),"Plat.DAT")
         !
         if(merge_Pi.and.solve_DMFT) then !a bit redundant since there is no merge wihtout DMFT
            call MergeFields(Plat,PiEDMFT,alphaPi,SiteOrbs)
            call DeallocateBosonicField(PiEDMFT)
            call dump_BosonicField(Plat,reg(ItFolder),"PiGG_merged.DAT")
         endif
         !
      endif
      !
      !Fully screened interaction
      if(calc_W)then
         !
         if(calc_Wfull)  call calc_W_full(Wlat,Ulat,Plat,Crystal)
         if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,PiEDMFT,Crystal)
         call dump_BosonicField(Wlat,reg(ItFolder),"Wloc.DAT")
         !
      endif
      !
      !K-dependent self-energy
      if(calc_Sigmak)then
         !
         !Hartree shift between G0W0 and LDA
         allocate(VH(Crystal%Norb,Crystal%Norb));VH=czero
         call calc_VH(densityLDA,Glat,Ulat,VH)
         if(solve_DMFT.and.bosonicSC)call DeallocateBosonicField(Ulat)
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
      call join_SigmaFull(Iteration)
      !
      !Compute the Full Green's function and set the density
      call calc_Gmats(Glat,Crystal,SigmaFull)
      if(look4dens%TargetDensity.ne.0d0)call set_density(Glat,Crystal,look4dens)
      !
      !Print Gf: local readable and k-dep binfmt
      call dump_FermionicField(Glat,1,reg(ItFolder),"Gloc_up.DAT")
      call dump_FermionicField(Glat,2,reg(ItFolder),"Gloc_dn.DAT")
      call dump_FermionicField(Glat,reg(ItFolder),"Gloc",.true.,Crystal%kpt)
      !
      !Matching the lattice and impurity problems
      if(solve_DMFT)then
         !
         do isite=1,Nsite
            !
            !Extract the hybridization functions and local energies (always diagonal)
            call calc_Delta(isite)
            !
            !Compute local effective interaction (loop on the sites is internal)
            call calc_Interaction(isite,ExpandImpurity)
            !
            if(ExpandImpurity)exit
            !
         enddo
         !
      endif
      !
   enddo
   !
   write(LOGfile,"(A,F)") "Self-Consistency finished. Total timing (s): ",tock(TimeStart)
   !
end program test
