program SelfConsistency
   !
   ! COMMENTS:
   ! tengo collect_QMC_results(ItStart) separato perche forse in un futuro riscrivero il solver
   ! e quindi potro semplicemnte sponstare collect_QMC_results(ItStart) in fondo e tenere
   ! initialize_Fields(ItStart) invariato, Se metto dyson dentro initialize_Fields
   ! dopo se cambio solver mi tocca riscriverlo.
   ! quindi tengo separato collect_QMC_results che crea e scrive le self-energie e Pimp
   ! oi dealloco e rileggo
   ! tutta la roba che alloco qui e' solo temporaea e principalemtne dovuta la fatto che thengo seprate C e X componenti di sigmaGW


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
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   call tick(TimeStart)
   call read_InputFile("input.in")
   write(*,"(A,1I4)") "Setting Nthread:",Nthread
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
   if(solve_DMFT.and.(ItStart.gt.0))call collect_QMC_results()
   !
   !
   !
   !---------------------------------------------------------------------------!
   !SOLVING THE LATTICE PROBLEM AND PRODUCING INPUTS FOR NEXT IMPURITY SOLUTION!
   !---------------------------------------------------------------------------!
   call initialize_Fields(ItStart)
   if(solve_DMFT.and.(ItStart.gt.0))call show_Densities(ItStart-1)
   !
   do Iteration=ItStart,Itend,1
      !
      call printHeader(Iteration)
      !
      !K-dependent Polarization - only G0W0,scGW,GW+EDMFT
      !calc_Plat=.false. !TEST
      if(calc_Plat)then
         !
         call calc_Pi(Plat,Glat,Crystal)
         call dump_BosonicField(Plat,reg(ItFolder),"Plat_w.DAT")
         !
         if(merge_Pi.and.solve_DMFT) then !a bit redundant since there is no merge wihtout DMFT
            call MergeFields(Plat,P_EDMFT,alphaPi,SiteOrbs)
            !call DeallocateBosonicField(P_EDMFT)
            call dump_BosonicField(Plat,reg(ItFolder),"Plat_merged_w.DAT")
         endif
         !
      endif
      !
      !Fully screened interaction - only G0W0,scGW,GW+EDMFT,EDMFT
      !calc_W=.false. !TEST
      if(calc_W)then
         !
         if(calc_Wfull)  call calc_W_full(Wlat,Ulat,Plat,Crystal)
         if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,P_EDMFT,Crystal)
         call dump_BosonicField(Wlat,reg(ItFolder),"Wlat_w.DAT")
         call dump_MaxEnt(Wlat,"mats",reg(ItFolder)//"Convergence/","Wlat",EqvGWndx%SetOrbs)
         !
      endif
      call DeallocateBosonicField(Plat)
      !
      !K-dependent self-energy - only G0W0,scGW,GW+EDMFT
      !calc_Sigmak=.false. !TEST
      if(calc_Sigmak)then
         !
         !Hartree shift between G0W0 and LDA
         allocate(VH(Crystal%Norb,Crystal%Norb));VH=czero
         call calc_VH(densityLDA,Glat,Ulat,VH)
         if(solve_DMFT.and.bosonicSC.and.(.not.Ustart))call DeallocateBosonicField(Ulat)
         !
         !G0W0 contribution and Vexchange readed from SPEX
         allocate(Vxc(Crystal%Norb,Crystal%Norb,Crystal%Nkpt,Nspin));Vxc=czero
         call AllocateFermionicField(S_G0W0,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_Sigma_spex(S_G0W0,Crystal,verbose,Vxc=Vxc,doAC=Sigma_AC)
         !
         !scGW
         call AllocateFermionicField(S_GW_C,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call AllocateFermionicField(S_GW_X,Crystal%Norb,0,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call calc_sigmaGW(S_GW_C,S_GW_X,Glat,Wlat,Crystal)
         !
         !Dc between G0W0 and scGW
         if(Iteration.eq.0)then
            !
            call AllocateFermionicField(S_G0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call join_SigmaCX(S_G0W0dc,S_GW_C,S_GW_X)
            call dump_FermionicField(S_G0W0dc,reg(pathDATA)//"0/","SGoWo",.true.,Crystal%kpt)
            call DeallocateFermionicField(S_G0W0dc)
            !
            !Dump G0W0 local self-energy
            call dump_FermionicField(S_G0W0,reg(ItFolder),"Slat_w")
            !
         elseif(Iteration.gt.0)then
            !
            call AllocateFermionicField(S_G0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call read_FermionicField(S_G0W0dc,reg(pathDATA)//"0/","SGoWo",kpt=Crystal%kpt)
            !
            call AllocateFermionicField(S_GW,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call join_SigmaCX(S_GW,S_GW_C,S_GW_X)
            !
            !Dump scGW local self-energy
            call dump_FermionicField(S_GW,reg(ItFolder),"Slat_w")
            !
         endif
         call DeallocateFermionicField(S_GW_C)
         call DeallocateFermionicField(S_GW_X)
         !
         !Merge GW and EDMFT
         if((Iteration.gt.0).and.merge_Sigma.and.solve_DMFT)then !a bit redundant since there is no merge wihtout DMFT
            !
            !First compute the local Dc
            call AllocateFermionicField(S_GWdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            call AllocateFermionicField(S_GW_Cdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            call AllocateFermionicField(S_GW_Xdc,Crystal%Norb,0,Nsite=Nsite,Beta=Beta)
            call calc_sigmaGWdc(S_GW_Cdc,S_GW_Xdc,Glat,Wlat)
            call join_SigmaCX(S_GWdc,S_GW_Cdc,S_GW_Xdc)
            call DeallocateFermionicField(S_GW_Cdc)
            call DeallocateFermionicField(S_GW_Xdc)
            !
            !Then replace local projection of scGW with EDMFT
            call MergeFields(S_GW,S_GWdc,S_DMFT,alphaSigma,SiteOrbs,DC_type,RotateHloc)
            call DeallocateFermionicField(S_GWdc)
            !
            call dump_FermionicField(S_GW,reg(ItFolder),"Slat_merged_w")
            !
         endif
         !
      endif
      !
      !Put together all the contributions to the full self-energy
      !call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta) !TEST
      call join_SigmaFull(Iteration)
      !
      !Compute the Full Green's function and set the density
      call calc_Gmats(Glat,Crystal,S_Full)
      call DeallocateFermionicField(S_Full)
      if(look4dens%TargetDensity.ne.0d0)call set_density(Glat,Crystal,look4dens)
      !
      !Print Gf: local readable and k-dep binfmt
      call dump_FermionicField(Glat,reg(ItFolder),"Glat_w")
      call dump_FermionicField(Glat,reg(ItFolder),"Glat_w",.true.,Crystal%kpt)
      call dump_MaxEnt(Glat,"mats",reg(ItFolder)//"Convergence/","Glat",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call dump_MaxEnt(Glat,"mats2itau",reg(ItFolder)//"Convergence/","Glat",EqvGWndx%SetOrbs)
      !
      !!Print lattice density
      call dump_Matrix(Glat%N_s(:,:,1),reg(ItFolder)//"Nlat_s1.DAT")
      call dump_Matrix(Glat%N_s(:,:,2),reg(ItFolder)//"Nlat_s2.DAT")
      densityGW=Glat%N_s
      !
      !Matching the lattice and impurity problems
      if(solve_DMFT)then
         !
         do isite=1,Nsite
            !
            !Extract the hybridization functions and local energies (always diagonal)
            call calc_Delta(isite,Iteration)
            !
            !Compute local effective interaction
            call calc_Interaction(isite,Iteration,ExpandImpurity)
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      endif
      !
      if(.not.solve_DMFT)call show_Densities(Iteration)
      !
   enddo !Iteration
   !
   call DeallocateAllFields()
   !
   write(*,"(A,F)") new_line("A")//new_line("A")//"Self-Consistency finished. Total timing (s): ",tock(TimeStart)
   !
   if(ItStart.ne.LastIteration)then
      call execute_command_line(" cp used.input.in "//reg(ItFolder))
      if(solve_DMFT)call execute_command_line(" touch doSolver ")
   endif
   !
end program SelfConsistency
