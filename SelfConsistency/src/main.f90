program SelfConsistency
   !
   use module_container
   use utils_main
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: isite
   integer                                  :: Iteration,ItStart,Itend
   !
#ifdef _akw
   character(len=20)                        :: InputFile="input.in.akw"
#elif defined _gap
   character(len=20)                        :: InputFile="input.in.gap"
#else
   character(len=20)                        :: InputFile="input.in"
#endif
   !
   !
   !
   !
   !---------------------------------------------------------------------------!
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   call tick(TimeStart)
   call read_InputFile(reg(InputFile))
   call printHeader()
   call initialize_DataStructure(ItStart,Itend)
   call initialize_Lattice(Crystal,ItStart)
   !
   !
#ifdef _akw
   !
   call inquireFile(reg(ItFolder)//"Sfull_w_k_s1.DAT",S_Full_exists,hardstop=.true.,verb=verbose)
   !
   call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
   call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
   call dump_MaxEnt(S_Full,"mats",reg(ItFolder)//"Convergence/","Sful",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
   if(interp_G) call interpolate2kpath(S_Full,Crystal,reg(MaxEnt_K))
   !
   ! get self-energy at Gamma
   S_Full%ws = S_Full%wks(:,:,:,1,:)
   call dump_FermionicField(S_Full,reg(ItFolder),"Sfull_w_Gamma",paramagnet)
   call dump_MaxEnt(S_Full,"mats",reg(ItFolder)//"Convergence/","Sful_Gamma",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
   call DeallocateFermionicField(S_Full)
   call execute_command_line(" touch doSolver ")
   stop
   !
#elif defined _gap
   !
   call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
   call calc_Tc(reg(ItFolder),gap_equation,Crystal,Wlat)
   stop
   !
#endif
   !
   !
   !
   !
   !
   !---------------------------------------------------------------------------!
   !       COLLECTING RESULTS FROM THE SOLVER AND SOLVING DYSON EQUATION       !
   !---------------------------------------------------------------------------!
   if(ItStart.gt.0)then
      if(Beta_Match%status)then
         call interpolate_from_oldBeta()
      else
         if(collect_QMC) call collect_QMC_results()
      endif
   endif
   !
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
   !
   !
   !
   do Iteration=ItStart,Itend,1
      !
      !
      call printHeader(Iteration)
      call update_DataStructure(Iteration)
      !
      !
      !Check if needed fields are already present
      call inquireFile(reg(ItFolder)//"Sfull_w_k_s1.DAT",S_Full_exists,hardstop=.false.,verb=verbose)
      if(S_Full_exists)then
         write(*,"(A)") new_line("A")//new_line("A")//"---- skipping S_Full calculation."
         calc_Sigmak=.false.
         call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
      endif
      !
      !
      !K-dependent Polarization - only G0W0,scGW,GW+EDMFT
      if(calc_Pk)then
         !
         if(Uscreen)then
            call calc_PiGG(Plat,Glat,Crystal)
            call dump_BosonicField(Plat,reg(ItFolder),"Plat_w.DAT")
            call dump_MaxEnt(Plat,"mats",reg(ItFolder)//"Convergence/","Plat",EqvGWndx%SetOrbs)
         endif
         !
         if(merge_P)then
            call AllocateBosonicField(P_GGdc,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            call calc_PiGGdc(P_GGdc,Glat)
            call dump_BosonicField(P_GGdc,reg(ItFolder),"Plat_dc_w.DAT")
            call MergeFields(Plat,P_EDMFT,P_GGdc,alphaPi,RotateHloc,LocalOrbs)
            call dump_BosonicField(Plat,reg(ItFolder),"Plat_merged_w.DAT")
            call DeallocateBosonicField(P_GGdc)
         elseif(calc_Pguess)then
            P_EDMFT%screened_local = dreal(Plat%screened_local)*alphaPi
         endif
         !
         if(interp_Chi) then
            call calc_chi(Chi,Ulat,Plat,Crystal)!,pathPk=reg(MaxEnt_K)//"/MaxEnt_Chik_path_t")
            call interpolate2kpath(Chi,Crystal,reg(MaxEnt_K),name="C",mode="Trace_NaNa")
            call DeallocateBosonicField(Chi)
         endif
         !
      endif
      !
      !
      !Fully screened interaction - only G0W0,scGW,GW+EDMFT,EDMFT
      if(calc_Wk)then
         !
         if(calc_Wfull)  call calc_W_full(Wlat,Ulat,Plat,Einv,Crystal)
         if(calc_Wedmft) call calc_W_edmft(Wlat,Ulat,P_EDMFT,Einv,Crystal,alpha=alphaPi)
         call dump_BosonicField(Wlat,reg(ItFolder),"Wlat_w.DAT")
         call dump_MaxEnt(Wlat,"mats",reg(ItFolder)//"Convergence/","Wlat",EqvGWndx%SetOrbs)
         !
         if(interp_W) call interpolate2kpath(Wlat,Crystal,reg(MaxEnt_K),name="W",mode="Trace_NaNa")
         if(interp_E) call interpolate2kpath(Einv,Crystal,reg(MaxEnt_K),name="E",mode="Trace_NaNa")
         !if(interp_E) call interpolate2kpath(Einv,Crystal,reg(MaxEnt_K),name="E",mode="Loss",invert=.true.)
         !call execute_command_line(" touch doSolver ")
         !stop
         call DeallocateBosonicField(Einv)
         !
         !Solve the Gap equation
         if(gap_equation%status)call calc_Tc(reg(ItFolder),gap_equation,Crystal,Wlat)
         !
      endif
      !
      !Causality correction on curlyU
      if(causal_U) call calc_causality_curlyU_correction(reg(causal_U_type))
      if(solve_DMFT) call DeallocateBosonicField(Plat)
      !
      !
      !Matching the lattice and impurity problems: Bosons
      if(solve_DMFT)then
         !
         !Compute local effective interaction
         do isite=1,Solver%Nimp
            call calc_Interaction(isite,Iteration,ExpandImpurity)
         enddo
         call DeallocateBosonicField(P_EDMFT)
         !
      endif
      !
      !
      !K-dependent self-energy - only G0W0,scGW,GW+EDMFT
      if(calc_Sigmak)then
         !
         !Hartree shift between G0W0 and scGW
         if(addTierIII)then
            !
            if(Iteration.eq.0)VN_type="None"
            if(MultiTier)VN_type="Nlat"
            !
            call calc_VH(VH_Nlat,densityLDA,densityGW,Ulat)
            call calc_VH(VH_Nimp,densityLDA,densityDMFT,Ulat,local=.true.)
            select case(reg(VN_type))
               case default
                  stop "Wrong entry for VN_TYPE. Available: Nlat, Nimp, None."
               case("Nlat")
                  VH = VH_Nlat
                  call dump_Matrix(VH_Nlat,reg(ItFolder),"VH_Nlat_used.DAT")
                  call dump_Matrix(VH_Nimp,reg(ItFolder),"VH_Nimp.DAT")
               case("Nimp")
                  VH = VH_Nimp
                  call dump_Matrix(VH_Nlat,reg(ItFolder),"VH_Nlat.DAT")
                  call dump_Matrix(VH_Nimp,reg(ItFolder),"VH_Nimp_used.DAT")
               case("None")
                  VH = czero
                  call dump_Matrix(VH_Nlat,reg(ItFolder),"VH_Nlat.DAT")
                  call dump_Matrix(VH_Nimp,reg(ItFolder),"VH_Nimp.DAT")
                  write(*,"(A)")"     VH not used."
            end select
            deallocate(VH_Nlat,VH_Nimp)
            !
            call calc_VH(Hartree_lat,densityGW)
            call dump_Matrix(Hartree_lat,reg(ItFolder),"Hartree_lat.DAT")
            deallocate(Hartree_lat)
            !
         endif
         !
         if(solve_DMFT.and.bosonicSC.and.(.not.Ustart))call DeallocateBosonicField(Ulat)
         !
         !read from SPEX G0W0 self-energy, double counting and Vexchange
         call AllocateFermionicField(S_G0W0,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         if(addTierIII) then
            !
            !G0W0 self-energy: used as self-energy in the 0th iteration or in model calculations
            call read_Sigma_spex(SpexVersion,S_G0W0,Crystal,verbose,RecomputeG0W0,Vxc)
            !
            !G0W0 double counting: this is used only for ab-initio calculations
            call AllocateFermionicField(S_G0W0dc,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            !read if exists otherwise compute a new one
            call inquireFile(reg(pathINPUTtr)//"SGoWo_dc_w_k_s1.DAT",S_G0W0dc_exist,hardstop=.false.,verb=verbose)
            if(S_G0W0dc_exist)then
               call read_FermionicField(S_G0W0dc,reg(pathINPUTtr),"SGoWo_dc_w",Crystal%kpt)
            else
               if(calc_S_G0W0dc)then
                  if(Iteration.gt.0) stop "The internal G0W0dc calculation must be performed in the 0th iteration."
                  write(*,"(A,F)")"     Computing dc between G0W0 and scGW."
                  call calc_sigmaGW(S_G0W0dc,Glat,Wlat,Crystal)
                  call dump_FermionicField(S_G0W0dc,reg(pathINPUTtr),"SGoWo_dc_w",.true.,Crystal%kpt,paramagnet)
                  call dump_FermionicField(S_G0W0dc,reg(pathINPUTtr),"SGoWo_dc_w",paramagnet)
               elseif(spex_S_G0W0dc)then
                  write(*,"(A,F)")"     Reading dc between G0W0 and scGW from SPEX_VERSION: "//reg(SpexVersion)
                  call read_Sigma_spex(SpexVersion,S_G0W0dc,Crystal,verbose,RecomputeG0W0,Vxc,DC=.true.)
               endif
            endif
            !
         endif
         !
         !scGW
         if(Iteration.eq.0)then
            !
            !Use directly the GW formula since G0W0 is absent for model calculations
            if(.not.addTierIII) call calc_sigmaGW(S_G0W0,Glat,Wlat,Crystal)
            !
            call dump_FermionicField(S_G0W0,reg(ItFolder),"Slat_w",paramagnet)
            call dump_MaxEnt(S_G0W0,"mats",reg(ItFolder)//"Convergence/","Slat",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
            !
         elseif(Iteration.gt.0)then
            !
            !Compute the scGW self-energy
            call AllocateFermionicField(S_GW,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call calc_sigmaGW(S_GW,Glat,Wlat,Crystal)
            call dump_FermionicField(S_GW,reg(ItFolder),"Slat_w",paramagnet)
            call dump_MaxEnt(S_GW,"mats",reg(ItFolder)//"Convergence/","Slat",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
            !
         endif
         !
         !Merge GW and EDMFT
         if(merge_Sigma)then
            call AllocateFermionicField(S_GWdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            select case(reg(DC_type_GW))
               case default
                  stop "Wrong entry for DC_TYPE_GW. Available: GlocWloc, Sloc."
               case("GlocWloc")
                  call calc_sigmaGWdc_GlocWloc(S_GWdc,Glat,Wlat)
               case("Slocs")
                  call calc_sigmaGWdc(S_GWdc,Glat,Wlat)
            end select
            call dump_FermionicField(S_GWdc,reg(ItFolder),"Slat_dc_w",paramagnet)
            call MergeFields(S_GW,S_DMFT,S_GWdc,alphaSigma,RotateHloc,LocalOrbs)
            call dump_FermionicField(S_GW,reg(ItFolder),"Slat_merged_w",paramagnet)
            call DeallocateFermionicField(S_GWdc)
         endif
         !
      endif
      if(solve_DMFT) call DeallocateBosonicField(Wlat)
      !
      !
      !Initial Guess for the impurity self-energy only in the 0th iteration
      if(calc_Sguess) call calc_SigmaGuess()
      !
      !
      !Put together all the contributions to the full self-energy and deallocate all non-local components: S_G0W0, S_G0W0dc, S_GW
      if(.not.S_Full_exists) call join_SigmaFull(Iteration)
      !
      !
      !Compute the Full Green's function and set the density
      if(mu_scan)then
         call calc_Gmats(Glat,Crystal,S_Full)
         if(Hetero%status)then
            call set_density(Glat,Crystal,look4dens,S_Full)
         else
            call set_density(Glat,Crystal,look4dens)
         endif
      else
         write(*,*)
         Glat%mu = look4dens%mu
         call calc_Gmats(Glat,Crystal,S_Full)
         write(*,"(A,F)")"     Chemical potential:",Glat%mu
         write(*,"(A,F)")"     Lattice density:",trace(Glat%N_s(:,:,1)+Glat%N_s(:,:,2))
      endif
      call dump_Matrix(Glat%N_s,reg(ItFolder),"Nlat",paramagnet)
      !
      densityGW = Glat%N_s
      Crystal%mu = Glat%mu
      if(.not.S_Full_exists) S_Full%mu = Glat%mu !if it exists I won't change mu
      !
      !
      !Print G0W0 bandstructure - Crystal%mu is used by default
      if(dump_G0W0_bands.and.Uwan_stored)call print_G0W0_dispersion(Crystal,Vxc)
      !
      !
      !Total energy calculation
      Ek = calc_Ek(Glat,Crystal)
      Ep = calc_Ep(Glat,S_Full)
      write(*,"(A,F)")"     Kinetic energy [eV]:",trace(Ek)
      write(*,"(A,F)")"     Potential energy [eV]:",trace(Ep)
      call dump_Matrix(Ek,reg(ItFolder),"Ek.DAT")
      call dump_Matrix(Ep,reg(ItFolder),"Ep.DAT")
      call check_QP_poles(Crystal,S_Full)
      !
      !
      !Print Gf: local readable and k-dep binfmt
      call dump_FermionicField(Glat,reg(ItFolder),"Glat_w",paramagnet)
      if(dump_Gk)call dump_FermionicField(Glat,reg(ItFolder),"Glat_w",.true.,Crystal%kpt,paramagnet)
      call dump_MaxEnt(Glat,"mats",reg(ItFolder)//"Convergence/","Glat",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call dump_MaxEnt(Glat,"mats2itau",reg(ItFolder)//"Convergence/","Glat",EqvGWndx%SetOrbs)
      !
      !
      !Print Potentials if present
      if(Hetero%status)call print_potentials()
      !
      !
      !Print the new full self-energy: local readable, k-dep binfmt (optional) and along path
      call dump_FermionicField(S_Full,reg(ItFolder),"Sfull_w",paramagnet)
      if(dump_Sigmak)call dump_FermionicField(S_Full,reg(ItFolder),"Sfull_w",.true.,Crystal%kpt,paramagnet)
      call dump_MaxEnt(S_Full,"mats",reg(ItFolder)//"Convergence/","Sful",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      if(interp_G)call interpolate2kpath(S_Full,Crystal,reg(MaxEnt_K))
      !
      !Print the new full self-energy at Gamma point
      S_Full%ws = S_Full%wks(:,:,:,1,:)
      call dump_FermionicField(S_Full,reg(ItFolder),"Sfull_w_Gamma",paramagnet)
      call dump_MaxEnt(S_Full,"mats",reg(ItFolder)//"Convergence/","Sful_Gamma",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call DeallocateFermionicField(S_Full)
      !
      !
      !Matching the lattice and impurity problems: Fermions
      if(solve_DMFT)then
         !
         !The local problem must give the same density in the same subset
         if(MultiTier)then
            write(*,*)
            write(*,"(A,1I3)") "     N_READ_IMP updated from "//str(Solver%TargetDensity,4)//" to "//str(get_Tier_occupation(densityGW,LocalOrbs),4)
            Solver%TargetDensity = get_Tier_occupation(densityGW,LocalOrbs)
            call save_InputFile(reg(InputFile))
         endif
         !
         ! Causality correction on Delta
         if(causal_D) call calc_causality_Delta_correction()
         !
         !Extract the hybridization functions and local energies (always diagonal)
         do isite=1,Solver%Nimp
            call calc_Delta(isite,Iteration)
         enddo
         !
      endif
      !
      if(.not.solve_DMFT)call show_Densities(Iteration)
      !
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
