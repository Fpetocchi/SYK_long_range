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
#ifdef _akw
   character(len=20)                        :: InputFile="input.in.akw"
#else
   character(len=20)                        :: InputFile="input.in"
#endif

#ifdef _gap
   character(len=20)                        :: InputFile="input.in.gap"
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
   call calc_Glda(0d0,Beta,Crystal)
   call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
   call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
   if(print_path) call interpolate2kpath(S_Full,Crystal,reg(ItFolder))
   !
   ! get self-energy at Gamma
   call AllocateFermionicField(Slat_Gamma,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
   Slat_Gamma%mu = S_Full%mu
   Slat_Gamma%ws = S_Full%wks(:,:,:,1,:)
   call dump_FermionicField(Slat_Gamma,reg(ItFolder),"Slat_Gamma_w",paramagnet)
   call dump_MaxEnt(Slat_Gamma,"mats",reg(ItFolder)//"Convergence/","Slat_Gamma",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
   call dump_MaxEnt(Slat_Gamma,"mats2itau",reg(ItFolder)//"Convergence/","Slat_Gamma",EqvGWndx%SetOrbs)
   !
   call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
   call read_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w")
   call dump_MaxEnt(S_DMFT,"mats",reg(ItFolder)//"Convergence/","Simp",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
   !
   stop
#elif defined _gap
   call calc_Tc(reg(ItFolder),gap_equation,Crystal)
   stop
#endif
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
         if(print_path) call interpolate2kpath(S_Full,Crystal,reg(ItFolder))
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
         call calc_PiGG(Plat,Glat,Crystal)
         call dump_BosonicField(Plat,reg(ItFolder),"Plat_w.DAT")
         call dump_MaxEnt(Plat,"mats",reg(ItFolder)//"Convergence/","Plat",EqvGWndx%SetOrbs)
         !
         if(merge_P)then
            if(reg(DC_type_P).eq."GlocGloc")then
               call AllocateBosonicField(P_GGdc,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               call calc_PiGGdc(P_GGdc,Glat)
               call MergeFields(Plat,P_EDMFT,alphaPi,LocalOrbs,RotateHloc,PiGG_DC=P_GGdc)
               call DeallocateBosonicField(P_GGdc)
            elseif(reg(DC_type_P).eq."Ploc")then
               call MergeFields(Plat,P_EDMFT,alphaPi,LocalOrbs,RotateHloc)
            endif
            call dump_BosonicField(Plat,reg(ItFolder),"Plat_merged_w.DAT")
         elseif(calc_Pguess)then
            P_EDMFT%screened_local = dreal(Plat%screened_local)*alphaPi
         endif
         !
         if(dump_Chik) then
            call calc_chi(Chi,Ulat,Plat,Crystal)!,pathPk=reg(MaxEnt_K)//"/MaxEnt_Chik_path_t")
            call interpolate2kpath(Chi,Crystal,reg(ItFolder),"C")
            call DeallocateBosonicField(Chi)
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
         if(dump_Wk) call interpolate2kpath(Wlat,Crystal,reg(ItFolder),"W")
         call dump_BosonicField(Wlat,reg(ItFolder),"Wlat_w.DAT")
         call dump_MaxEnt(Wlat,"mats",reg(ItFolder)//"Convergence/","Wlat",EqvGWndx%SetOrbs)
         !
      endif
      !
      ! Causality correction on curlyU
      if(causal_U) call calc_causality_curlyU_correction(reg(causal_U_type))
      if(solve_DMFT) call DeallocateBosonicField(Plat)
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
            call calc_VH(VH_Nlat,densityLDA,densityGW,Ulat)
            call calc_VH(VH_Nimp,densityLDA,densityDMFT,Ulat)
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
         endif
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
                  write(*,"(A,F)")"     Computing dc between G0W0 and scGW."
                  call calc_sigmaGW(S_G0W0dc,Glat,Wlat,Crystal)!,LDAoffdiag=.false.) I believed the scGWdc should have had OD terms removed but its not working
                  call dump_FermionicField(S_G0W0dc,reg(pathINPUTtr),"SGoWo_dc_w",.true.,Crystal%kpt,paramagnet)
                  call dump_FermionicField(S_G0W0dc,reg(pathINPUTtr),"SGoWo_dc_w",paramagnet)
               elseif(spex_S_G0W0dc)then
                  write(*,"(A,F)")"     Reading dc between G0W0 and scGW from SPEX_VERSION: "//reg(SpexVersion)
                  call read_Sigma_spex(SpexVersion,S_G0W0dc,Crystal,verbose,RecomputeG0W0,Vxc,DC=.true.)
               endif
            endif
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
            !
            if(reg(DC_type_S).eq."GlocWloc")then
               call AllocateFermionicField(S_GWdc,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call calc_sigmaGWdc(S_GWdc,Glat,Wlat)
               call MergeFields(S_GW,S_DMFT,[alphaSigma,HartreeFact],LocalOrbs,SigmaGW_DC=S_GWdc)
               call DeallocateFermionicField(S_GWdc)
            elseif(reg(DC_type_S).eq."Sloc")then
               call MergeFields(S_GW,S_DMFT,[alphaSigma,HartreeFact],LocalOrbs)
            endif
            call dump_FermionicField(S_GW,reg(ItFolder),"Slat_merged_w",paramagnet)
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
      !Put together all the contributions to the full self-energy
      !and deallocate all non-local components: S_G0W0, S_G0W0dc, S_GW
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
      if(solve_DMFT) call DeallocateBosonicField(Wlat)
      !
      !
      !Compute the Full Green's function and set the density
      call calc_Gmats(Glat,Crystal,S_Full)
      if(look4dens%TargetDensity.ne.0d0)then
         if(Hetero%status)then
            call set_density(Glat,Crystal,look4dens,S_Full)
         else
            call set_density(Glat,Crystal,look4dens)
         endif
      else
         write(*,"(A,F)")"     Chemical potential:",Glat%mu
         write(*,"(A,F)")"     Lattice density:",trace(Glat%N_s(:,:,1)+Glat%N_s(:,:,2))
      endif
      call dump_Matrix(Glat%N_s,reg(ItFolder),"Nlat",paramagnet)
      !
      densityGW = Glat%N_s
      S_Full%mu = Glat%mu
      !
      !
      !Total energy calculation
      Ek = calc_Ek(Glat,Crystal)
      Ep = calc_Ep(Glat,S_Full)
      write(*,"(A,F)")"     Kinetic energy [eV]:",trace(Ek)
      write(*,"(A,F)")"     Potential energy [eV]:",trace(Ep)
      call dump_Matrix(Ek,reg(ItFolder),"Ek.DAT")
      call dump_Matrix(Ep,reg(ItFolder),"Ep.DAT")
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
      !Print full self-energy: local readable, k-dep binfmt (optional) and along path
      if(calc_Sigmak)then
         call dump_FermionicField(S_Full,reg(ItFolder),"Sfull_w",paramagnet)
         if(dump_Sigmak)call dump_FermionicField(S_Full,reg(ItFolder),"Sfull_w",.true.,Crystal%kpt,paramagnet)
         call dump_MaxEnt(S_Full,"mats",reg(ItFolder)//"Convergence/","Sful",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
         if(print_path)call interpolate2kpath(S_Full,Crystal,reg(ItFolder))
      endif
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





















!
