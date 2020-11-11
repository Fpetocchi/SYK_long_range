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
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
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
   if(solve_DMFT.and.(ItStart.gt.0))call collect_QMC_results()
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
         call dump_BosonicField(Plat,reg(ItFolder),"Plat_w.DAT")
         !
         if(merge_Pi.and.solve_DMFT) then !a bit redundant since there is no merge wihtout DMFT
            call MergeFields(Plat,PiEDMFT,alphaPi,SiteOrbs)
            call DeallocateBosonicField(PiEDMFT)
            if(verbose)call dump_BosonicField(Plat,reg(ItFolder),"Plat_w_merged.DAT")
         endif
         !
      endif
      !
      !Fully screened interaction
      if(calc_W)then
         !
         if(calc_Wfull)  call calc_W_full(Wlat,Ulat,Plat,Crystal)
         if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,PiEDMFT,Crystal)
         call dump_BosonicField(Wlat,reg(ItFolder),"Wlat_w.DAT")
         !
      endif
      !
      !K-dependent self-energy
      if(calc_Sigmak)then
         !
         !Hartree shift between G0W0 and LDA
         allocate(VH(Crystal%Norb,Crystal%Norb));VH=czero
         call calc_VH(densityLDA,Glat,Ulat,VH)
         if(solve_DMFT.and.bosonicSC.and.(.not.Ustart))call DeallocateBosonicField(Ulat)
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
            call dump_FermionicField(SigmaG0W0dc,reg(pathDATA)//"0/","SGoWo",.true.,Crystal%kpt)
            call DeallocateFermionicField(SigmaG0W0dc)
            !
         elseif(Iteration.gt.0)then
            !
            call AllocateFermionicField(SigmaG0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call read_FermionicField(SigmaG0W0dc,reg(pathDATA)//"0/","SGoWo",kpt=Crystal%kpt)
            !
            call AllocateFermionicField(SigmaGW,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call join_SigmaCX(SigmaGW,SigmaGW_C,SigmaGW_X)
            !
         endif
         call DeallocateFermionicField(SigmaGW_C)
         call DeallocateFermionicField(SigmaGW_X)
         !
         !Dump scGW local self-energy
         call dump_FermionicField(Glat,1,reg(ItFolder),"Slat_w_up.DAT")
         call dump_FermionicField(Glat,2,reg(ItFolder),"Slat_w_dw.DAT")
         !
         !Merge GW and EDMFT
         if((Iteration.gt.0).and.merge_Sigma.and.solve_DMFT)then !a bit redundant since there is no merge wihtout DMFT
            !
            !First compute the local Dc
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
            if(verbose)then
               call dump_FermionicField(SigmaGW,1,reg(ItFolder),"Slat_w_up_merged.DAT")
               call dump_FermionicField(SigmaGW,2,reg(ItFolder),"Slat_w_dw_merged.DAT")
            endif
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
      call dump_FermionicField(Glat,1,reg(ItFolder),"Glat_w_up.DAT")
      call dump_FermionicField(Glat,2,reg(ItFolder),"Glat_w_dn.DAT")
      call dump_FermionicField(Glat,reg(ItFolder),"Glat_w",.true.,Crystal%kpt)
      !
      !!Print lattice density
      call dump_Matrix(Glat%N_s(:,:,1),reg(ItFolder)//"Nlat_up.DAT")
      call dump_Matrix(Glat%N_s(:,:,2),reg(ItFolder)//"Nlat_dw.DAT")
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
            if(ExpandImpurity)exit
            !
         enddo
         !
      endif
      !
   enddo
   !
   call DeallocateAllFields()
   call execute_command_line(" cp used.input.in "//reg(ItFolder))
   call execute_command_line(" cp report "//reg(ItFolder))
   call execute_command_line(" cp err "//reg(ItFolder))
   write(LOGfile,"(A,F)") "Self-Consistency finished. Total timing (s): ",tock(TimeStart)
   !
end program test
