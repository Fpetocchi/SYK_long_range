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
   character(len=20)                        :: InputFile="input.in.gap"
   !
   !
   !
   !
   !---------------------------------------------------------------------------!
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   call tick(TimeStart)
   call read_InputFile(reg(InputFile))
   write(*,"(A,1I4)") "Setting Nthread:",Nthread
   call printHeader()
   call initialize_DataStructure(ItStart,Itend)
   call initialize_Lattice(Crystal,ItStart)
   !
   !
   !call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
   !call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
   !call interpolateG2Path(S_Full,Crystal,reg(structure),Nkpt_path,reg(ItFolder))
   !stop

   call calc_Tc(reg(ItFolder),gap_equation,Crystal)
   stop

   !
   !
   !---------------------------------------------------------------------------!
   !       COLLECTING RESULTS FROM THE SOLVER AND SOLVING DYSON EQUATION       !
   !---------------------------------------------------------------------------!
   if(ItStart.gt.0)then
      if(Beta_Match%status)then
         call interpolate_from_oldBeta()
      else
         if(solve_DMFT) call collect_QMC_results()
      endif
   endif
   !
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
      !
      call printHeader(Iteration)
      !
      !
      !Check if needed fields are already present
      call inquireFile(reg(ItFolder)//"Wlat_w.DAT",Wlat_exists,hardstop=.false.,verb=verbose)
      call inquireFile(reg(ItFolder)//"Sfull_w_k_s1.DAT",S_Full_exists,hardstop=.false.,verb=verbose) !add spin2 ?
      !
      if(S_Full_exists)then
         !
         write(*,"(A)") new_line("A")//new_line("A")//"---- skipping S_Full calculation."
         calc_Sigmak=.false.
         call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
         if(print_path) call interpolateG2Path(S_Full,Crystal,reg(structure),Nkpt_path,reg(ItFolder))
         !
         if(Wlat_exists.and.(.not.gap_equation%status))then !(*)
            write(*,"(A)") new_line("A")//new_line("A")//"---- skipping Plat and Wlat calculations."
            calc_Pk=.false.
            calc_W=.false.
            call read_BosonicField(Wlat,reg(ItFolder),"Wlat_w.DAT")
         endif
         !
      endif
      !
      !
      !K-dependent Polarization - only G0W0,scGW,GW+EDMFT
      if(calc_Pk)then
         !
         call calc_Pi(Plat,Glat,Crystal)
         call dump_BosonicField(Plat,reg(ItFolder),"Plat_w.DAT")
         call dump_MaxEnt(Plat,"mats",reg(ItFolder)//"Convergence/","Plat",EqvGWndx%SetOrbs)
         !
         if(merge_P)then
            call MergeFields(Plat,P_EDMFT,alphaPi,SiteOrbs,RotateHloc)
            call dump_BosonicField(Plat,reg(ItFolder),"Plat_merged_w.DAT")
         elseif(calc_Pguess)then
            P_EDMFT%screened_local = Plat%screened_local*alphaPi
         endif
         !
      endif
      !
      !
      !Fully screened interaction - only G0W0,scGW,GW+EDMFT,EDMFT
      if(calc_W)then
         !
         if(calc_Wfull)  call calc_W_full(Wlat,Ulat,Plat,Crystal)
         if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,P_EDMFT,Crystal)
         call dump_BosonicField(Wlat,reg(ItFolder),"Wlat_w.DAT")
         call dump_MaxEnt(Wlat,"mats",reg(ItFolder)//"Convergence/","Wlat",EqvGWndx%SetOrbs)
         !
      endif
      !
      !
      ! Causality correction on curlyU
      if(causal_U) call calc_causality_curlyU_correction()
      call DeallocateBosonicField(Plat)
      !
      !
      !Matching the lattice and impurity problems: Bosons
      if(solve_DMFT)then
         !
         do isite=1,Nsite
            !
            !Compute local effective interaction
            call calc_Interaction(isite,Iteration,ExpandImpurity)
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      endif
      !
      !
      !K-dependent self-energy - only G0W0,scGW,GW+EDMFT
      if(calc_Sigmak)then
         !
         !Hartree shift between G0W0 and LDA
         allocate(VH(Crystal%Norb,Crystal%Norb));VH=czero
         if(.not.Hmodel)call calc_VH(VH,densityLDA,densityDMFT,Ulat) !call calc_VH(VH,densityLDA,Glat,Ulat)
         call dump_Matrix(VH,reg(ItFolder)//"VH.DAT")
         if(.not.VH_use)then
            VH=czero
            write(*,"(A)")"     VH not used."
         endif
         if(solve_DMFT.and.bosonicSC.and.(.not.Ustart))call DeallocateBosonicField(Ulat)
         !
         !read from SPEX G0W0 self-energy and Vexchange
         allocate(Vxc(Crystal%Norb,Crystal%Norb,Crystal%Nkpt,Nspin));Vxc=czero
         call AllocateFermionicField(S_G0W0,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         if(.not.Hmodel)then
            if(Vxc_in)then
               call read_Sigma_spex(S_G0W0,Crystal,verbose,doAC=Sigma_AC,pathOUTPUT=reg(pathINPUTtr))
            else
               call read_Sigma_spex(S_G0W0,Crystal,verbose,doAC=Sigma_AC,pathOUTPUT=reg(pathINPUTtr),Vxc_out=Vxc)
            endif
         endif
         !
         !scGW
         if(Iteration.eq.0)then
            !
            !Compute the Dc between G0W0 and scGW self-energies
            call AllocateFermionicField(S_G0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call calc_sigmaGW(S_G0W0dc,Glat,Wlat,Crystal)
            !
            if(Hmodel)then
               !Use directly the DC since G0W0 is absent for model calculations
               call duplicate(S_G0W0,S_G0W0dc)
            else
               !Store the Dc between G0W0 and scGW self-energies and use G0W0 as self-energy for the first iteration because the DC is equal to scGW
               call dump_FermionicField(S_G0W0dc,reg(ItFolder),"SGoWo_dc_w",.true.,Crystal%kpt)
               call dump_FermionicField(S_G0W0dc,reg(ItFolder),"SGoWo_dc_w")
               call DeallocateFermionicField(S_G0W0dc)
            endif
            !
            call dump_FermionicField(S_G0W0,reg(ItFolder),"Slat_w")
            !
         elseif(Iteration.gt.0)then
            !
            !Read the Dc between G0W0 and scGW if present
            call AllocateFermionicField(S_G0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(.not.Hmodel)call read_FermionicField(S_G0W0dc,reg(pathDATA)//"0/","SGoWo_dc_w",kpt=Crystal%kpt)
            !
            !Compute the scGW self-energy
            call AllocateFermionicField(S_GW,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call calc_sigmaGW(S_GW,Glat,Wlat,Crystal)
            call dump_FermionicField(S_GW,reg(ItFolder),"Slat_w")
            call dump_MaxEnt(S_GW,"mats",reg(ItFolder)//"Convergence/","Slat",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
            !
         endif
         !
         !Merge GW and EDMFT
         if(merge_Sigma)then
            !
            !Compute the local Dc if Tier1 = Tier2
            if(reg(DC_type).eq."GlocWloc")then
               call AllocateFermionicField(S_GWdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call calc_sigmaGWdc(S_GWdc,Glat,Wlat)
            endif
            !
            !Then replace local projection of scGW with EDMFT, S_GWdc is empty for DC_type=Sloc or for model calculations
            call MergeFields(S_GW,S_GWdc,S_DMFT,[alphaSigma,HartreeFact],SiteOrbs,DC_type,RotateHloc)
            call DeallocateFermionicField(S_GWdc)
            call dump_FermionicField(S_GW,reg(ItFolder),"Slat_merged_w")
            !
         endif
         !
      endif
      !
      !
      !Initial Guess for the impurity self-energy only in the 0th iteration
      if(calc_Sguess) call calc_SigmaGuess()
      !
      !
      !Put together all the contributions to the full self-energy and deallocate all the components
      if(.not.S_Full_exists) call join_SigmaFull(Iteration)
      !
      !
      !Solve the Gap equation
      if(gap_equation%status)then
         if(gap_equation%HkRenorm)then
            call calc_Tc(reg(ItFolder),gap_equation,Crystal,Wlat=Wlat,Sfull=S_Full)
         else
            call calc_Tc(reg(ItFolder),gap_equation,Crystal,Wlat=Wlat)
         endif
      endif
      call DeallocateBosonicField(Wlat)
      !
      !
      !Compute the Full Green's function and set the density
      call calc_Gmats(Glat,Crystal,S_Full)
      if(look4dens%TargetDensity.ne.0d0)then
         call set_density(Glat,Crystal,look4dens)
      else
         write(*,"(A,F)")"     Lattice density:",trace(Glat%N_s(:,:,1)+Glat%N_s(:,:,2))
      endif
      !
      !
      !Update the full self-energy, print and compute Glat along the path
      S_Full%mu = Glat%mu
      if(dump_Sigmak)call dump_FermionicField(S_Full,reg(ItFolder),"Sfull_w",.true.,Crystal%kpt)
      if(print_path)call interpolateG2Path(S_Full,Crystal,reg(structure),Nkpt_path,reg(ItFolder))
      call DeallocateFermionicField(S_Full)
      !
      !
      ! Causality correction on Delta
      if(causal_D) call calc_causality_Delta_correction()
      !
      !
      !Print Gf: local readable and k-dep binfmt
      if(dump_Gk)call dump_FermionicField(Glat,reg(ItFolder),"Glat_w",.true.,Crystal%kpt)
      call dump_FermionicField(Glat,reg(ItFolder),"Glat_w")
      call dump_MaxEnt(Glat,"mats",reg(ItFolder)//"Convergence/","Glat",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call dump_MaxEnt(Glat,"mats2itau",reg(ItFolder)//"Convergence/","Glat",EqvGWndx%SetOrbs)
      !
      !
      !Print lattice density
      call dump_Matrix(Glat%N_s(:,:,1),reg(ItFolder)//"Nlat_s1.DAT")
      call dump_Matrix(Glat%N_s(:,:,2),reg(ItFolder)//"Nlat_s2.DAT")
      densityGW=Glat%N_s
      !
      !
      !The local problem must give the same density in the same subset
      if(MultiTier)then
         write(*,*)
         write(*,"(A,1I3)") "     N_READ_IMP updated from "//str(Solver%TargetDensity,4)//" to "//str(get_Tier_occupation(densityGW,SiteOrbs),4)
         Solver%TargetDensity = get_Tier_occupation(densityGW,SiteOrbs)
         call save_InputFile(reg(InputFile))
      endif
      !
      !
      !Matching the lattice and impurity problems: Fermions
      if(solve_DMFT)then
         !
         do isite=1,Nsite
            !
            !Extract the hybridization functions and local energies (always diagonal)
            call calc_Delta(isite,Iteration)
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
      call execute_command_line(" cp used."//reg(InputFile)//" "//reg(ItFolder))
      if(solve_DMFT)call execute_command_line(" touch doSolver ")
   endif
   !
end program SelfConsistency



!
!
!
!------------------------------------------------------------------------------!
!                                   COMMENTS                                   !
!------------------------------------------------------------------------------!
!
!(*): the S_Full_exists flag is present because I need the full K dependent Wlat,
!which is never stored, in ordeer to compute S_Full. If the gap equation is solved
!then the calculation of the full k-dependent Wlat is required.
!
!(**): a bit redundant since there is no merge wihtout DMFT
!
!(***): NON MI DEVO PREOCCUPARE SE AGGIUNGO A S_DMFT ORBITALI CHE STANNO DENTRO
!A TIER-2 TANTO QUANDO MI FACCIO IL DYSON PER CALCOLARE DELTA IMP2LOC MI ESTRAE DA
!S_DMFT SOLO QUELLO CHE PARE A ME POI, DOPO calc_Delta, S_DMFT NON  È  MAI PIU USATA NEL  CODICE
!For the case where the orbital subspaces for the EDMFT and GW calculations are
!the same, DC reduces to the local projection (k sum) of the full GW self-energy.
!If the orbital subspace for the EDMFT calculation is smaller than that of the GW
!calculation the difference between DC and the local projection of the full GW self-energy is...

!call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
!call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
!call interpolateG2Path(S_Full,Crystal,reg(structure),Nkpt_path,reg(ItFolder))























!
