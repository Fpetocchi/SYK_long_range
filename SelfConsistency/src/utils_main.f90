module utils_main

   use module_container
   implicit none
   public

   !===========================================================================!

   ! COMMENTS:
   ! this is not a library so Im putt much less controls on dimensions here
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   !interface name
   !   module procedure name_a                                                  ! Description
   !   module procedure name_b                                                  ! Description
   !end interface name

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
#ifdef _verb
   logical                                  :: verbose=.true.
#else
   logical                                  :: verbose=.false.
#endif
   !
   type(Lattice)                            :: Crystal
   !
   type(FermionicField)                     :: Glat
   type(FermionicField)                     :: G_DMFT
   !
   type(FermionicField)                     :: S_Full
   type(FermionicField)                     :: S_G0W0,S_G0W0dc
   type(FermionicField)                     :: S_GW,S_GW_C,S_GW_X
   type(FermionicField)                     :: S_GWdc,S_GW_Cdc,S_GW_Xdc
   type(FermionicField)                     :: S_DMFT
   !
   type(FermionicField)                     :: Slat_Gamma
   !
   type(BosonicField)                       :: Wlat
   type(BosonicField)                       :: W_EDMFT
   !
   type(BosonicField)                       :: Ulat
   type(BosonicField)                       :: Plat,P_GGdc
   type(BosonicField)                       :: Chi
   type(BosonicField)                       :: P_EDMFT
   type(BosonicField)                       :: C_EDMFT
   type(BosonicField)                       :: curlyU_EDMFT
   !
   type(FermionicField)                     :: Delta_correction
   type(BosonicField)                       :: curlyU_correction
   !
   real(8),allocatable                      :: Ek(:,:)
   real(8),allocatable                      :: Ep(:,:)
   !
   complex(8),allocatable                   :: densityLDA(:,:,:)
   complex(8),allocatable                   :: densityGW(:,:,:)
   complex(8),allocatable                   :: densityDMFT(:,:,:)
   !
   complex(8),allocatable                   :: Vxc(:,:,:,:)
   complex(8),allocatable                   :: VH(:,:)
   complex(8),allocatable                   :: VH_Nlat(:,:)
   complex(8),allocatable                   :: VH_Nimp(:,:)
   !
   real(8),allocatable                      :: Umat(:,:)
   real(8),allocatable                      :: Kfunct(:,:)
   !
   logical                                  :: Wlat_exists=.false.
   logical                                  :: S_Full_exists=.false.
   !
   logical                                  :: calc_Pk=.false.
   logical                                  :: merge_P=.false.
   logical                                  :: calc_W=.false.
   logical                                  :: calc_Wfull=.false.
   logical                                  :: calc_Wedmft=.false.
   logical                                  :: calc_Sigmak=.false.
   logical                                  :: merge_Sigma=.false.
   !
   logical                                  :: S_G0W0dc_exist=.false.
   logical                                  :: calc_S_G0W0dc=.false.
   logical                                  :: spex_S_G0W0dc=.false.
   !
   character(len=255)                       :: HartreeType="GW" ! "GW" "DMFT"
   real(8)                                  :: HartreeFact=1d0
   logical                                  :: update_curlyG=.true.
   integer                                  :: SigmaMaxMom=3
   !
   logical                                  :: MultiTier=.false.
   logical                                  :: print_path=.false.
   !
   character(len=255)                       :: ItFolder,PrevItFolder,MixItFolder
   character(len=256)                       :: MaxEnt_K

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: printHeader
   public :: initialize_DataStructure
   public :: initialize_Lattice
   public :: initialize_Fields

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: prints the header
   !---------------------------------------------------------------------------!
   subroutine printHeader(Iteration)
      !
      implicit none
      integer,intent(in),optional           :: Iteration
      integer,parameter                     :: bannerLen=100
      character(len=bannerLen)              :: line
      character(:),allocatable              :: header
      integer                               :: i,Lsx,Ldx,Lcalc
      !
      if(present(Iteration))then
         !
         header="<"
         Lcalc=len(" Iteration #: "//str(Iteration,3)//" ")
         Lsx=int(((bannerLen-2)-Lcalc)/2)-1
         Ldx=(bannerLen-2)-Lcalc-Lsx
         do i=1,Lsx
            header = header//"="
         enddo
         header = header//" Iteration #: "//str(Iteration,3)//" "
         do i=1,Ldx
            header = header//"="
         enddo
         header = header//">"
         write(*,"(A)") new_line("A")//new_line("A")//new_line("A")//header//new_line("A")
         !
      else
         !
         header="*"
         line="*"
         do i=1,(bannerLen-1)
            line = trim(line)//"*"
         enddo
         !
         Lcalc=len("Calculation type: "//trim(CalculationType))
         !
         Lsx=int(((bannerLen-2)-Lcalc)/2)-1
         Ldx=(bannerLen-2)-Lcalc-Lsx
         !
         do i=1,Lsx
            header = header//" "
         enddo
         header = header//"Calculation type: "//trim(CalculationType)
         do i=1,Ldx
            header = header//" "
         enddo
         header = header//"*"
         write(*,"(A)") line
         write(*,"(A)") header
         write(*,"(A)") line//new_line("A")//new_line("A")
         !
      endif
      !
   end subroutine printHeader


   !---------------------------------------------------------------------------!
   !PURPOSE: looks for the current iteration number
   !---------------------------------------------------------------------------!
   subroutine initialize_DataStructure(ItStart,Itend)
      !
      implicit none
      integer,intent(out)                   :: ItStart
      integer,intent(out)                   :: Itend
      character(len=256)                    :: Itpath
      integer                               :: isite,ispin
      logical                               :: PrvItexist,ZeroItexist
      logical                               :: JulichDCexists,LundDCexists,DCcond
      !
      !Few general checks
      if(ExpandImpurity.and.AFMselfcons) stop "AFM self-consistency and expansion to real space not implemented."
      if(ExpandImpurity.and.(Nsite.eq.1)) stop "Cannot expand a single site."
      if(AFMselfcons.and.(Nsite.ne.2)) stop "AFM self-consistency is implemented only for lattices with 2 sites."
      if(RotateUloc.and.(.not.RotateHloc)) stop "Rotate the Bosonic impurity problem without rotating the Ferminic one is not allowed."
      !
      if(.not.RemoveHartree)HartreeFact=0d0
      !
      calc_Sguess = calc_Sguess .and. (FirstIteration.eq.0) .and. solve_DMFT
      calc_Pguess = calc_Pguess .and. (FirstIteration.eq.0) .and. solve_DMFT
      !
      Ustart = Ustart .and. (FirstIteration.eq.0) .and. (.not.calc_Pguess)
      !
      causal_D = causal_D .and. ((FirstIteration.ne.0).or.calc_Sguess)
      causal_U = causal_U .and. ((FirstIteration.ne.0).or.calc_Pguess)
      !
      if(Hmodel.or.Umodel)addTierIII=.false.
      !
      if(addTierIII)then
         !
         calc_S_G0W0dc=.true.
         spex_S_G0W0dc=.false.
         !
         call inquireFile(reg(pathINPUT)//"Sigma_dc_imag/spex.out",JulichDCexists,hardstop=.false.,verb=.false.)
         call inquireDir(reg(pathINPUT)//"Sigma_dc_real_01",LundDCexists,hardstop=.false.,verb=.false.)
         !
         DCcond = ( JulichDCexists .and. (reg(SpexVersion).eq."Julich") ) .or. ( LundDCexists .and. (reg(SpexVersion).eq."Lund") )
         if(DCcond)then
            calc_S_G0W0dc=.false.
            spex_S_G0W0dc=.true.
         endif
         !
      endif
      !
      print_path = (reg(path_funct).ne."None") .and. (reg(structure).ne."None")
      !
      call inquireDir(reg(pathDATA)//str(FirstIteration-1),PrvItexist,hardstop=.false.,verb=.false.)
      call inquireDir(reg(pathDATA)//str(0),ZeroItexist,hardstop=.false.,verb=.false.)
      !
      Itpath = reg(pathDATA)//str(FirstIteration)
      ItStart = FirstIteration
      !
      select case(reg(CalculationType))
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0")
            !
            if(ItStart.ne.0) stop "CalculationType is G0W0 but the starting iteration is not 0."
            call createDir(reg(Itpath),verb=verbose)
            do isite=1,Nsite
               call createDir(reg(Itpath)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/fits",verb=verbose)
               if(ExpandImpurity.or.AFMselfcons)exit
            enddo
            !
            if(Uspex)then
               call createDir(reg(pathINPUTtr)//"VW_imag",verb=verbose)
               if(verbose)then
                  call createDir(reg(pathINPUTtr)//"VW_imag_readable",verb=verbose)
                  call createDir(reg(pathINPUTtr)//"VW_real_readable",verb=verbose)
               endif
            endif
            !
            Itend = 0
            !
         case("scGW")
            !
            if(ItStart.ne.0)then
               if(.not.PrvItexist)then
                  write(*,"(A)") "     Previous iteration: "//str(FirstIteration-1)//". Not found, exiting."
                  stop
               endif
            endif
            !
            if(Uspex)then
               call createDir(reg(pathINPUTtr)//"VW_imag",verb=verbose)
               if(verbose)then
                  call createDir(reg(pathINPUTtr)//"VW_imag_readable",verb=verbose)
                  call createDir(reg(pathINPUTtr)//"VW_real_readable",verb=verbose)
               endif
            endif
            !
            call createDir(reg(Itpath),verb=verbose)
            !
            Itend = LastIteration
            !
         case("DMFT+statU","DMFT+dynU","EDMFT","GW+EDMFT")
            !
            if(ItStart.ne.0)then
               if(.not.PrvItexist)then
                  write(*,"(A)") "     Previous iteration: "//str(FirstIteration-1)//". Not found, exiting."
                  stop
               endif
            endif
            !
            if(Uspex)then
               call createDir(reg(pathINPUTtr)//"VW_imag",verb=verbose)
               if(verbose)then
                  call createDir(reg(pathINPUTtr)//"VW_imag_readable",verb=verbose)
                  call createDir(reg(pathINPUTtr)//"VW_real_readable",verb=verbose)
               endif
            endif
            !
            if(ItStart.eq.LastIteration) then
               !
               Itend = ItStart-1 !avoid the lattice in the last it
               !
            else
               !
               call createDir(reg(Itpath),verb=verbose)
               do isite=1,Nsite
                  call createDir(reg(Itpath)//"/Solver_"//reg(LocalOrbs(isite)%Name)//"/fits",verb=verbose)
                  if(ExpandImpurity.or.AFMselfcons)exit
               enddo
               Itend = ItStart
               !
            endif
            !
      end select
      !
      !For memory demanding applications one has to creat the folders at the beginning
      call createDir(reg(Itpath)//"/Convergence",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Glat",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Gimp",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Slat",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Simp",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Cimp",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Ulat",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Wlat",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Wimp",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/curlyUimp",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Plat",verb=verbose)
      call createDir(reg(Itpath)//"/Convergence/Pimp",verb=verbose)
      !
      !Print info to user
      if(FirstIteration.eq.0)then
         write(*,"(A)") "     Brand new calculation. Initializing "//reg(Itpath)
      else
          write(*,"(A)") "     Initializing "//reg(Itpath)
       endif
      !
      !These are used throughout the whole calculation
      ItFolder = reg(pathDATA)//str(ItStart)//"/"
      PrevItFolder = reg(pathDATA)//str(ItStart-1)//"/"
      MixItFolder = reg(pathDATA)//str(ItStart-Mixing_period)//"/"
      !
      MaxEnt_K = reg(ItFolder)//"K_resolved/"
      !
      !Creates folders for the K-resolved data
      if(print_path)then
         !
         call createDir(reg(MaxEnt_K)//"/Sigma_vars",verb=verbose)
         !
         do ispin=1,Nspin
            if(reg(path_funct).eq."G")then
               !
               call createDir(reg(MaxEnt_K)//"/MaxEnt_Gk_full_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on G in the full BZ
               call createDir(reg(MaxEnt_K)//"/MaxEnt_Gk_path_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on G along the K-path
               if(FermiSurf)call createDir(reg(MaxEnt_K)//"/MaxEnt_Gk_plane_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on G in the full BZ to get Fermi surface
               !
            elseif(reg(path_funct).eq."S")then
               !
               if((reg(CalculationType).eq."G0W0").or.(reg(CalculationType).eq."scGW").or.(reg(CalculationType).eq."GW+EDMFT"))then
                  !
                  call createDir(reg(MaxEnt_K)//"/MaxEnt_Sk_path_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on S along the K-path
                  call createDir(reg(MaxEnt_K)//"/Gk_path_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S along the K-path
                  call createDir(reg(MaxEnt_K)//"/Sk_path_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S along the K-path
                  !
                  call createDir(reg(MaxEnt_K)//"/MaxEnt_Sk_full_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on S in the full BZ to get Gloc
                  call createDir(reg(MaxEnt_K)//"/Gk_full_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the full BZ to get Gloc
                  call createDir(reg(MaxEnt_K)//"/Sk_full_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the full BZ to get Gloc
                  !
                  if(FermiSurf)then
                     call createDir(reg(MaxEnt_K)//"/MaxEnt_Sk_plane_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on S in the {kx,ky} plane to get the Fermi surface
                     call createDir(reg(MaxEnt_K)//"/Gk_plane_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the {kx,ky} plane to get the Fermi surface
                     call createDir(reg(MaxEnt_K)//"/Sk_plane_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the {kx,ky} plane to get the Fermi surface
                  endif
                  !
               endif
               !
            elseif(reg(path_funct).eq."GS")then
               !
               call createDir(reg(MaxEnt_K)//"/MaxEnt_Gk_path_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on G along the K-path
               if(FermiSurf)call createDir(reg(MaxEnt_K)//"/MaxEnt_Gk_plane_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on G in the full BZ to get Fermi surface
               !
               if((reg(CalculationType).eq."G0W0").or.(reg(CalculationType).eq."scGW").or.(reg(CalculationType).eq."GW+EDMFT"))then
                  !
                  call createDir(reg(MaxEnt_K)//"/MaxEnt_Sk_path_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on S along the K-path
                  call createDir(reg(MaxEnt_K)//"/Gk_path_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S along the K-path
                  call createDir(reg(MaxEnt_K)//"/Sk_path_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S along the K-path
                  !
                  call createDir(reg(MaxEnt_K)//"/MaxEnt_Sk_full_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on S in the full BZ to get Gloc
                  call createDir(reg(MaxEnt_K)//"/Gk_full_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the full BZ to get Gloc
                  call createDir(reg(MaxEnt_K)//"/Sk_full_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the full BZ to get Gloc
                  !
                  if(FermiSurf)then
                     call createDir(reg(MaxEnt_K)//"/MaxEnt_Sk_plane_t_s"//str(ispin),verb=verbose)  ! This is for MaxEnt on S in the {kx,ky} plane to get the Fermi surface
                     call createDir(reg(MaxEnt_K)//"/Gk_plane_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the {kx,ky} plane to get the Fermi surface
                     call createDir(reg(MaxEnt_K)//"/Sk_plane_wm_s"//str(ispin),verb=verbose)        ! This is for MaxEnt on S in the {kx,ky} plane to get the Fermi surface
                  endif
                  !
               endif
               !
            endif
            if(paramagnet)exit
         enddo
      endif
      !
      if((reg(path_funct).ne."G") .and. (reg(path_funct).ne."S") .and. (reg(path_funct).ne."GS") .and. (reg(path_funct).ne."None"))then
         write(*,"(A)")"     Invalid content of PATH_FUNCT. The calculation on the path will be avoided."
         print_path=.false.
      elseif(reg(path_funct).eq."None")then
         print_path=.false.
      endif
      !
      write(*,"(A)") "     HartreeFact: "//str(HartreeFact,2)
      !
      !Save changes to the inputfile
      call save_InputFile("input.in")
      !
   end subroutine initialize_DataStructure


   !---------------------------------------------------------------------------!
   !PURPOSE: Initialize Lattice. I could have used the AllocateLattice in
   !         utils_fields but then als useless attributes would have been allocated
   !---------------------------------------------------------------------------!
   subroutine initialize_Lattice(Lttc,ItStart)
      !
      implicit none
      type(Lattice),intent(out)             :: Lttc
      integer,intent(in)                    :: ItStart
      integer                               :: isite,iadd,iset,iorb,ik
      integer                               :: iq_gamma_Hk,iq_gamma_XEPS
      integer                               :: shift
      logical                               :: present
      integer,allocatable                   :: oldSetNorb(:),oldSetOrbs(:,:)
      real(8),allocatable                   :: Egrid(:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- initialize_Lattice"
      !
      !
      if(Hmodel)then
         call set_lattice(LatticeVec,ucVec)
      else
         call read_lattice(reg(pathINPUT))
      endif
      !
      if(readHk)then
         call read_Hk(Lttc%Hk,Lttc%kpt,pathINPUT)
      else
         call build_Hk(Lttc%Hk,Lttc%kpt,hopping,Nkpt3,readHr,Hetero,pathOUTPUT=reg(pathINPUT))
      endif
      !
      !Filling the Lattice attributes
      Lttc%Nkpt3 = Nkpt3
      Lttc%Nkpt = product(Nkpt3)
      Lttc%Norb = size(Lttc%Hk,dim=1)
      Lttc%Nsite = Nsite
      !
      call fill_ksumkdiff(Lttc%kpt,Lttc%kptsum,Lttc%kptdif,Nkpt3)
      !
      if(Hmodel)then
         Lttc%Hk = Lttc%Hk * alphaHk
      else
         Lttc%Hk = Lttc%Hk * alphaHk * H2eV
      endif
      !
      if(allocated(Lttc%Hloc))deallocate(Lttc%Hloc)
      allocate(Lttc%Hloc(Lttc%Norb,Lttc%Norb));Lttc%Hloc=czero
      Lttc%Hloc = sum(Lttc%Hk,dim=3)/Lttc%Nkpt
      !
      if(allocated(Lttc%Zk))deallocate(Lttc%Zk)
      if(allocated(Lttc%Ek))deallocate(Lttc%Ek)
      allocate(Lttc%Zk(Lttc%Norb,Lttc%Norb,Lttc%Nkpt));Lttc%Zk=czero
      allocate(Lttc%Ek(Lttc%Norb,Lttc%Nkpt));Lttc%Ek=0d0
      do ik=1,Lttc%Nkpt
         Lttc%Ek(:,ik) = 0d0
         Lttc%Zk(:,:,ik) = Lttc%Hk(:,:,ik)
         call eigh(Lttc%Zk(:,:,ik),Lttc%Ek(:,ik))
      enddo
      !
      !additional data neede for ab-initio calculations
      iq_gamma_Hk = 1
      if(Lttc%Nkpt.gt.1) iq_gamma_Hk = find_vec([0d0,0d0,0d0],Lttc%kpt,eps)
      write(*,"(A)")"     Gamma point index: "//str(iq_gamma_Hk)
      !
      if(.not.Hmodel)then
         allocate(Lttc%kptPos(Lttc%Nkpt));Lttc%kptPos=0
         call read_xeps(pathINPUT,Lttc%kpt,Nkpt3,UseXepsKorder, &
         Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,iq_gamma_XEPS,paramagneticSPEX)
         XEPSisread=.true.
         if(iq_gamma_Hk.eq.iq_gamma_XEPS)then
            Lttc%iq_gamma = iq_gamma_Hk
         else
            stop "Index of the Gamma point is not the same in H(k) and XEPS."
         endif
      else
         Lttc%Nkpt_irred = Lttc%Nkpt
         Lttc%iq_gamma = iq_gamma_Hk
      endif
      !
      Lttc%status=.true.
      !
      if(Lttc%Nkpt.ne.(Nkpt3(1)*Nkpt3(2)*Nkpt3(3)))stop "Total number of K-points does not match with number of K-points per dimension."
      if(XEPSisread)write(*,"(A,1I4)")"     Number of irreducible K-points: ",Lttc%Nkpt_irred
      !
      !Estimate of the bandwidth plus
      if(wrealMax.eq.0d0)wrealMax = 1.5*maxval(abs(Lttc%Ek))
      !
      !Store the local Hamiltonian
      call dump_Matrix(Lttc%Hloc,reg(pathINPUT),"Hloc.DAT")
      !
      !allocate Ek and Ep matrices
      if(allocated(Ek))deallocate(Ek)
      if(allocated(Ep))deallocate(Ep)
      allocate(Ek(Lttc%Norb,Lttc%Norb));Ek=0d0
      allocate(Ep(Lttc%Norb,Lttc%Norb));Ep=0d0
      !
      !print some info
      write(*,"(A)") new_line("A")//new_line("A")//"---- site and orbital structure"
      write(*,"(A,1I3)") "     Number of inequivalent sites: ",Nsite
      write(*,"(A,1I3)") "     Number of solved impurities : ",Solver%Nimp
      do isite=1,Nsite
         write(*,"(2(A,I3),A,10I3)")"     site: ",isite,", orbital space: ",LocalOrbs(isite)%Norb,", orbital indexes: ",LocalOrbs(isite)%Orbs
      enddo
      if(sum(LocalOrbs(:)%Norb).ne.Lttc%Norb)then
         !
         MultiTier = .true.
         if(Hetero%status) stop "MultiTier construction and Heterostructured setup are not allowed together."
         if(RotateHloc.or.RotateUloc) stop "MultiTier construction and Rotations of the local space are not allowed together."
         !
      endif
      !
      !
      !Store the local rotation of each site and add to the input list the Impurity equivalent orbitals
      if(RotateHloc)then
         !
         write(*,"(A)") new_line("A")//new_line("A")//"---- Rotations of the local LDA Hamiltonian (used)"
         call build_rotations("Hloc",LatticeOp=Lttc%Hloc)
         call update_ImpEqvOrbs()
         !
      else
         !
         call update_ImpEqvOrbs()
         !
      endif
      !
      !
      !Symmetrization list:
      !if(EqvGWndx%Nset.eq.0) --> Nothing to symmetrize
      !if((EqvGWndx%Nset.ne.0).and.(.not.ExpandImpurity)) --> Singularly include in the list all the orbitals not included in the user provided list (if any)
      !if((EqvGWndx%Nset.ne.0).and.ExpandImpurity) ---------> Expand the list like the orbitals are expanded
      if(sym_mode.eq.0)EqvGWndx%Nset=0
      if(EqvGWndx%Nset.eq.0)then
         !
         EqvGWndx%O=.false.
         EqvGWndx%Ntotset=EqvGWndx%Nset
         EqvGWndx%Gfoffdiag=.false.
         !
      elseif((EqvGWndx%Nset.ne.0).and.(.not.ExpandImpurity))then
         !
         EqvGWndx%O = sym_mode.le.2
         EqvGWndx%Ntotset = EqvGWndx%Nset
         EqvGWndx%Gfoffdiag = sym_mode.le.2
         !
         if(sum(EqvGWndx%SetNorb).lt.Lttc%Norb)then
            !
            EqvGWndx%Ntotset = EqvGWndx%Nset + (Lttc%Norb-sum(EqvGWndx%SetNorb))
            !
            !reshape of SetNorb
            allocate(oldSetNorb(size(EqvGWndx%SetNorb)))
            oldSetNorb=EqvGWndx%SetNorb
            deallocate(EqvGWndx%SetNorb)
            allocate(EqvGWndx%SetNorb(EqvGWndx%Ntotset));EqvGWndx%SetNorb=0
            EqvGWndx%SetNorb(1:EqvGWndx%Nset) = oldSetNorb
            EqvGWndx%SetNorb(1+EqvGWndx%Nset:EqvGWndx%Ntotset) = 1
            deallocate(oldSetNorb)
            !reshape of SetOrbs
            allocate(oldSetOrbs(EqvGWndx%Nset,size(EqvGWndx%SetOrbs,dim=2)))
            oldSetOrbs=EqvGWndx%SetOrbs
            deallocate(EqvGWndx%SetOrbs)
            allocate(EqvGWndx%SetOrbs(EqvGWndx%Ntotset,size(EqvGWndx%SetOrbs,dim=2)));EqvGWndx%SetOrbs=0
            EqvGWndx%SetOrbs(1:EqvGWndx%Nset,:)=oldSetOrbs
            !refilling of SetOrbs
            iadd=0
            do iorb=1,Lttc%Norb
               present=.false.
               do iset=1,EqvGWndx%Nset
                  present = present.or.any(oldSetOrbs(iset,:).eq.iorb)
               enddo
               if(.not.present)then
                  iadd=iadd+1
                  EqvGWndx%SetOrbs(EqvGWndx%Nset+iadd,1)=iorb
               endif
            enddo
            deallocate(oldSetOrbs)
            !
         endif
         !
      elseif((EqvGWndx%Nset.ne.0).and.ExpandImpurity)then
         !
         if(EqvGWndx%Nset.ne.1) stop "Orbital symmetrization is not implemented for more than one set if the impurity has to be expanded."
         !
         EqvGWndx%O = sym_mode.le.2
         EqvGWndx%Ntotset = EqvGWndx%Nset + (Nsite-1)
         EqvGWndx%Gfoffdiag = .false.
         !
         !reshape of SetNorb
         allocate(oldSetNorb(size(EqvGWndx%SetNorb)))
         oldSetNorb=EqvGWndx%SetNorb
         deallocate(EqvGWndx%SetNorb)
         allocate(EqvGWndx%SetNorb(EqvGWndx%Ntotset));EqvGWndx%SetNorb=0
         EqvGWndx%SetNorb(1:EqvGWndx%Nset) = oldSetNorb
         EqvGWndx%SetNorb(1+EqvGWndx%Nset:EqvGWndx%Ntotset) = oldSetNorb(1)
         deallocate(oldSetNorb)
         !reshape of SetOrbs
         allocate(oldSetOrbs(1,size(EqvGWndx%SetOrbs,dim=2)))
         oldSetOrbs(1,:)=EqvGWndx%SetOrbs(1,:)
         deallocate(EqvGWndx%SetOrbs)
         allocate(EqvGWndx%SetOrbs(EqvGWndx%Ntotset,size(EqvGWndx%SetOrbs,dim=2)));EqvGWndx%SetOrbs=0
         EqvGWndx%SetOrbs(1,:)=oldSetOrbs(1,:)
         !refilling of SetOrbs
         do isite=2,Nsite
            ! only two possible arrangements
            if(abs(oldSetOrbs(1,2)-oldSetOrbs(1,1)).eq.1)then
               shift = size(oldSetOrbs(1,:))*(isite-1)
            elseif(abs(oldSetOrbs(1,2)-oldSetOrbs(1,1)).eq.Nsite)then
               shift = isite-1
            endif
            EqvGWndx%SetOrbs(isite,:) = oldSetOrbs(1,:) + shift
         enddo
         deallocate(oldSetOrbs)
         !
      endif
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- symmetrized orbital sets"
      write(*,"(A,1I3)") "     Provided sets: ",EqvGWndx%Nset
      write(*,"(A,1I3)") "     Expanded sets: ",EqvGWndx%Ntotset
      do iset=1,size(EqvGWndx%SetOrbs,dim=1)
         write(*,"(2(A,I3),A,10I3)")"     set: ",iset,", number of orbitals: ",EqvGWndx%SetNorb(iset),", indexes: ",EqvGWndx%SetOrbs(iset,1:EqvGWndx%SetNorb(iset))
      enddo
      write(*,"(A,L)") "     Symmetrizing off-diagonal: ",EqvGWndx%Gfoffdiag
      if(sym_mode.eq.1)write(*,"(A)") "     Only lattice quantities will be symmetrized."
      if(sym_mode.eq.2)write(*,"(A)") "     Both lattice and impurity quantities will be symmetrized."
      if(sym_mode.eq.3)write(*,"(A)") "     Only impurity quantities will be symmetrized."
      !
      !
      !Dump some LDA results
      if(reg(structure).eq."User") call set_UserPath(UserPath)
      if(ItStart.eq.0)then
         !
         call calc_Glda(0d0,Beta,Lttc)
         if(Hetero%status)call print_potentials(pathINPUT)
         !
         allocate(Egrid(Nreal));Egrid=0d0
         Egrid = linspace(-wrealMax,+wrealMax,Nreal)
         call tetrahedron_integration(reg(pathINPUT),Lttc%Hk,Lttc%Nkpt3,Lttc%kpt,Egrid,fact_intp=2,pathOUTPUT=reg(pathINPUT))
         deallocate(Egrid)
         !
         if(reg(structure).ne."None")call interpolateHk2Path(Lttc,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),doplane=.true.,hetero=Hetero)
         !
      endif
      !
   end subroutine initialize_Lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute or read the site-dependent rotations with respect of the
   !         LatticeOp. Fixed rotations are stored in pathINPUT, iteration
   !         dependent rotations are stored in the iteration folder.
   !---------------------------------------------------------------------------!
   subroutine build_rotations(Opname,LatticeOp,read)
      !
      implicit none
      !
      character(len=*),intent(in)           :: Opname
      complex(8),intent(in),optional        :: LatticeOp(:,:)
      logical,intent(in),optional           :: read
      !
      character(len=255)                    :: Folder
      integer                               :: isite,Norb
      real(8),allocatable                   :: EigR(:,:)
      logical                               :: storeIt,read_,build_
      !
      if(.not.allocated(LocalOrbs)) stop "build_rotations: LocalOrbs not properly initialized."
      !
      read_=.false.
      if(present(read))read_=read
      build_=.false.
      if(present(LatticeOp))build_=.true.
      !
      if(read_.and.build_)then
         write(*,"(A)") "     Lattice operator for rotation contruction provided but not used. Rotations will be read from pathINPUT."
         build_=.false.
      elseif((.not.(read_)).and.(.not.(build_)))then
         stop "build_rotations: lattice operator for rotation contruction not provided."
      endif
      !
      select case(reg(Opname))
         case default
            stop "Available rotation Operators are: Hloc, Hren, Nloc."
         case("Hloc","Hren")
            Folder=reg(pathINPUT)
            storeIt=.false.
         case("Nloc")
            storeIt=.true.
      end select
      !
      do isite=1,Nsite
         if(allocated(LocalOrbs(isite)%Op))    deallocate(LocalOrbs(isite)%Op)
         if(allocated(LocalOrbs(isite)%Rot))   deallocate(LocalOrbs(isite)%Rot)
         if(allocated(LocalOrbs(isite)%RotDag))deallocate(LocalOrbs(isite)%RotDag)
         if(allocated(LocalOrbs(isite)%Eig))   deallocate(LocalOrbs(isite)%Eig)
      enddo
      !
      do isite=1,Nsite
         !
         if(storeIt) Folder = reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/"
         !
         Norb = LocalOrbs(isite)%Norb
         if( ExpandImpurity .and. (Norb.ne.LocalOrbs(1)%Norb) ) stop "The orbital space is not the same among the different impurities."
         !
         allocate(LocalOrbs(isite)%Op(Norb,Norb));LocalOrbs(isite)%Op=czero
         allocate(LocalOrbs(isite)%Eig(Norb));LocalOrbs(isite)%Eig=0d0
         allocate(LocalOrbs(isite)%Rot(Norb,Norb));LocalOrbs(isite)%Rot=czero
         allocate(LocalOrbs(isite)%RotDag(Norb,Norb));LocalOrbs(isite)%RotDag=czero
         !
         if(build_)then
            !
            !Extract then local Hamiltonian for each site
            call loc2imp(LocalOrbs(isite)%Op,LatticeOp,LocalOrbs(isite)%Orbs)
            !
            !Always real rotations assumed for now
            LocalOrbs(isite)%Rot = dreal(LocalOrbs(isite)%Op)
            !
            !Rotate
            call eigh(LocalOrbs(isite)%Rot,LocalOrbs(isite)%Eig)
            LocalOrbs(isite)%RotDag = conjg(transpose(LocalOrbs(isite)%Rot))
            !
            !Save
            call dump_Matrix(LocalOrbs(isite)%Op,reg(Folder),reg(Opname)//"Site_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            call dump_Matrix(LocalOrbs(isite)%Rot,reg(Folder),reg(Opname)//"Rot_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            call dump_Matrix(LocalOrbs(isite)%RotDag,reg(Folder),reg(Opname)//"RotDag_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            call dump_Matrix(diag(LocalOrbs(isite)%Eig),reg(Folder),reg(Opname)//"Eig_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            !
         elseif(read_)then
            !
            call read_Matrix(LocalOrbs(isite)%Op,reg(Folder)//reg(Opname)//"Site_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            call read_Matrix(LocalOrbs(isite)%Rot,reg(Folder)//reg(Opname)//"Rot_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            call read_Matrix(LocalOrbs(isite)%RotDag,reg(Folder)//reg(Opname)//"RotDag_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            allocate(EigR(Norb,Norb));EigR=0d0
            call read_Matrix(EigR,reg(Folder)//reg(Opname)//"Eig_"//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//".DAT")
            LocalOrbs(isite)%Eig = diagonal(EigR)
            deallocate(EigR)
            !
         endif
         !
         !SO(3)check
         write(*,"(A,"//str(Norb)//"I3)") "     Orbs: ",LocalOrbs(isite)%Orbs
         write(*,"(A,2F)") "     det(Rot) of "//reg(LocalOrbs(isite)%Name)//"_"//str(isite)//" :",det(LocalOrbs(isite)%Rot)
         !
      enddo
      !
   end subroutine build_rotations


   !---------------------------------------------------------------------------!
   !PURPOSE: Initialize the site-dependent equivalent sets
   !---------------------------------------------------------------------------!
   subroutine update_ImpEqvOrbs()
      !
      implicit none
      integer                               :: iset,isite,Nimp
      integer                               :: set_orb,set_ndx,unit
      character(len=255)                    :: file
      logical                               :: contained
      !
      Nimp = Nsite
      if(ExpandImpurity.or.AFMselfcons) Nimp=1
      !
      if(allocated(EqvImpndx))deallocate(EqvImpndx)
      allocate(EqvImpndx(Nimp))
      !
      if(RotateHloc)then
         !
         !loop over sites
         do isite=1,Nsite
            !
            !look for pattern given by the diagonal lattice operator
            call get_pattern(EqvImpndx(isite)%SetOrbs,LocalOrbs(isite)%Eig,1e4*eps)
            !
            if(allocated(EqvImpndx(isite)%SetOrbs))then
               !
               EqvImpndx(isite)%Nset = size(EqvImpndx(isite)%SetOrbs,dim=1)
               allocate(EqvImpndx(isite)%SetNorb(EqvImpndx(isite)%Nset))
               do iset=1,EqvImpndx(isite)%Nset
                  EqvImpndx(isite)%SetNorb(iset) = size( pack( EqvImpndx(isite)%SetOrbs(iset,:), EqvImpndx(isite)%SetOrbs(iset,:).gt.0 ) )
               enddo
               EqvImpndx(isite)%O=.true.
               !
            else
               !
               EqvImpndx(isite)%Nset = 0
               EqvImpndx(isite)%O=.false.
               !
            endif
            !
            !This might cause some problem when only a subset of impurity orbitals are symmetrized
            EqvImpndx(isite)%Ntotset = EqvImpndx(isite)%Nset
            EqvImpndx(isite)%para = EqvGWndx%para
            EqvImpndx(isite)%Gfoffdiag = .true. ! beacuse of curlyU
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      else
         !
         !loop over sites
         do isite=1,Nsite
            !
            !dimensional initialization
            EqvImpndx(isite) = EqvGWndx
            EqvImpndx(isite)%Nset = 0
            EqvImpndx(isite)%SetNorb = 0
            EqvImpndx(isite)%SetOrbs = 0
            EqvImpndx(isite)%para = EqvGWndx%para
            !
            !loop over lattice sets
            do iset=1,EqvGWndx%Nset
               !
               !scroll elements in the set
               contained=.true.
               do set_ndx=1,EqvGWndx%SetNorb(iset)
                  set_orb = EqvGWndx%SetOrbs(iset,set_ndx)
                  contained = contained.and.any( LocalOrbs(isite)%Orbs .eq. set_orb )
               enddo
               !
               !the site contains all the indexes within the set
               if(contained.and.(EqvGWndx%SetNorb(iset).gt.1))then
                  EqvImpndx(isite)%Nset = EqvImpndx(isite)%Nset+1
                  EqvImpndx(isite)%SetNorb(EqvImpndx(isite)%Nset) = EqvGWndx%SetNorb(iset)
                  EqvImpndx(isite)%SetOrbs(EqvImpndx(isite)%Nset,:) = EqvGWndx%SetOrbs(iset,:) - (LocalOrbs(isite)%Orbs(1)-1)
                  EqvImpndx(isite)%O=.true.
               endif
               !
            enddo
            !
            !This might cause some problem when only a subset of impurity orbitals are symmetrized
            EqvImpndx(isite)%Ntotset = EqvImpndx(isite)%Nset
            EqvImpndx(isite)%para = EqvGWndx%para
            EqvImpndx(isite)%Gfoffdiag = .true. ! beacuse of curlyU
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      endif
      !
      !write lists
      if(sym_mode.gt.1)then
         !
         if(FirstIteration.ne.LastIteration)then
            do isite=1,Nsite
               if(EqvImpndx(isite)%Nset.gt.0)then
                  file = reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/Eqv.DAT"
                  unit = free_unit()
                  open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
                  do iset=1,EqvImpndx(isite)%Nset
                     write(unit,"("//str(EqvImpndx(isite)%SetNorb(iset))//"I4)")EqvImpndx(isite)%SetOrbs(iset,1:EqvImpndx(isite)%SetNorb(iset))-1
                  enddo
                  close(unit)
               endif
               if(ExpandImpurity.or.AFMselfcons)exit
            enddo
         endif
         !
         call add_separator("Solver symmetrizations")
         do isite=1,Nsite
            call append_to_input_list(EqvImpndx(isite)%Nset,"EQV_"//str(isite)//"_SETS","Number of sets of locally equivalent orbitals in site number "//str(isite)) !User cannot set the EQV_IMP_* fields as they are deduced from EQV_*, EXPAND and ROTATE_F.
            if(EqvImpndx(isite)%Nset.gt.0)then
               do iset=1,EqvImpndx(isite)%Nset
                  call append_to_input_list(EqvImpndx(isite)%SetNorb(iset),"EQV_"//str(isite)//"_NORB_"//str(iset),"Number of equivalent orbitals in the set number "//str(iset)//" in site number "//str(isite))
               enddo
               do iset=1,EqvImpndx(isite)%Nset
                  call append_to_input_list(EqvImpndx(isite)%SetOrbs(iset,1:EqvImpndx(isite)%SetNorb(iset)),"EQV_"//str(isite)//"_ORBS_"//str(iset),"Orbital indexes of equivalent set number "//str(iset)//" in site number "//str(isite))
               enddo
            endif
            if(ExpandImpurity.or.AFMselfcons)exit
         enddo
         !
      endif
      !
      !update input file
      call save_InputFile("input.in")
      !
   end subroutine update_ImpEqvOrbs


   !---------------------------------------------------------------------------!
   !PURPOSE: Initialize the Fields depending on the starting iteration
   !---------------------------------------------------------------------------!
   subroutine initialize_Fields(ItStart)
      !
      implicit none
      integer,intent(in)                    :: ItStart
      logical                               :: filexists
      integer                               :: unit,idum,ib1,ik
      integer                               :: iorb,isite,ispin
      character(len=255)                    :: file
      real(8)                               :: muQMC
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- initialize_Fields"
      !
      !
      if(.not.Crystal%status) stop "Crystal is not properly initialized."
      !
      select case(reg(CalculationType))
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0","scGW")
            !
            !Unscreened interaction
            if(Uspex)then
               call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call read_U_spex(Ulat,save2readable=verbose,kpt=Crystal%kpt,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.gt.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph,Hetero,LocalOnly=.false.)
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
               else
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal,Hetero)
               endif
            endif
            !
            !Fully screened interaction
            call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            !Polarization
            call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            !
            if(dump_Chik)call AllocateBosonicField(Chi,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            !
            !Lattice Gf
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(Glat,reg(PrevItFolder),"Glat_w",Crystal%kpt)
               if(look4dens%TargetDensity.eq.0d0) Glat%mu = look4dens%mu
            else
               Glat%mu = look4dens%mu
               if((look4dens%TargetDensity.ne.0d0).and.Hmodel) call set_density(Glat%mu,Beta,Crystal,look4dens)
               call calc_Gmats(Glat,Crystal)
            endif
            call calc_density(Glat,Crystal,Glat%N_ks)
            call calc_density(Glat,Glat%N_s)
            !
            !Logical Flags
            calc_Pk = .true.
            calc_Wfull = .true.
            calc_Sigmak = .true.
            !
         case("DMFT+statU")
            !
            !Hubbard interaction
            allocate(Umat(Crystal%Norb**2,Crystal%Norb**2));Umat=0d0
            if(Uspex)then
               call read_U_spex(Umat)
            elseif(Umodel)then
               call inquireFile(reg(pathINPUT)//"Umat_model.DAT",filexists,hardstop=.false.,verb=verbose)
               if(filexists)then
                  call read_Matrix(Umat,reg(pathINPUT)//"Umat_model.DAT")
               else
                  call build_Umat(Umat,Uaa,Uab,J)
                  call dump_Matrix(Umat,reg(pathINPUT),"Umat_model.DAT")
               endif
            endif
            !
            !Impurity Self-energy
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0) call read_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w")
            !
            !Lattice local Gf
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(Glat,reg(PrevItFolder),"Glat_w")
               if(look4dens%TargetDensity.eq.0d0) Glat%mu = look4dens%mu
            else
               Glat%mu = look4dens%mu
               if((look4dens%TargetDensity.ne.0d0).and.Hmodel) call set_density(Glat%mu,Beta,Crystal,look4dens)
               call calc_Gmats(Glat,Crystal)
            endif
            call calc_density(Glat,Glat%N_s)
            !
         case("DMFT+dynU")
            !
            !Unscreened interaction
            if(Uspex)then
               call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
               call read_U_spex(Ulat,save2readable=verbose,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.gt.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
                  call inquireFile(reg(pathINPUTtr)//"Uloc_mats.DAT",filexists,hardstop=.false.,verb=verbose)
                  if(filexists)then
                     call read_BosonicField(Ulat,reg(pathINPUTtr),"Uloc_mats.DAT")
                  else
                     call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph,Hetero)
                  endif
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
               else
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal,Hetero,LocalOnly=.true.)
               endif
            endif
            !
            !Impurity Self-energy
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0) call read_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w")
            !
            !Lattice local Gf
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(Glat,reg(PrevItFolder),"Glat_w")
               if(look4dens%TargetDensity.eq.0d0) Glat%mu = look4dens%mu
            else
               Glat%mu = look4dens%mu
               if((look4dens%TargetDensity.ne.0d0).and.Hmodel) call set_density(Glat%mu,Beta,Crystal,look4dens)
               call calc_Gmats(Glat,Crystal)
            endif
            call calc_density(Glat,Glat%N_s)
            !
         case("EDMFT")
            !
            !Unscreened interaction
            if(Uspex)then
               call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call read_U_spex(Ulat,save2readable=verbose,kpt=Crystal%kpt,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.gt.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph,Hetero,LocalOnly=.false.)
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
               else
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal,Hetero)
               endif
            endif
            !
            !Fully screened interaction
            if((ItStart.eq.0).and.(.not.Ustart))then
               call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            else
               call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=0,Nsite=Nsite,Beta=Beta)
            endif
            !
            if(dump_Chik)call AllocateBosonicField(Chi,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            !
            !Polarization
            call AllocateBosonicField(P_EDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            if(ItStart.ne.0)call read_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_w.DAT")
            !
            !Impurity Self-energy
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w")
               call read_Matrix(S_DMFT%N_s,reg(PrevItFolder)//"Hartree_UNimp",paramagnet)
            endif
            !
            !Lattice local Gf
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(Glat,reg(PrevItFolder),"Glat_w")
               if(look4dens%TargetDensity.eq.0d0) Glat%mu = look4dens%mu
            else
               Glat%mu = look4dens%mu
               if((look4dens%TargetDensity.ne.0d0).and.Hmodel) call set_density(Glat%mu,Beta,Crystal,look4dens)
               call calc_Gmats(Glat,Crystal)
            endif
            call calc_density(Glat,Glat%N_s)
            !
            !Logical Flags
            if(ItStart.eq.0)then
               if(.not.Ustart)then
                  calc_Pk = .true.
                  calc_Wfull = .true.
               endif
            else
               calc_Wedmft = .true.
            endif
            !
         case("GW+EDMFT")
            !
            !Unscreened interaction
            if(Uspex)then
               call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call read_U_spex(Ulat,save2readable=verbose,kpt=Crystal%kpt,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.gt.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph,Hetero,LocalOnly=.false.)
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
               else
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal,Hetero)
               endif
            endif
            !
            !Fully screened interaction
            call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            !Polarization
            call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            call AllocateBosonicField(P_EDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            if(ItStart.ne.0)call read_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_w.DAT")
            !
            if(dump_Chik)call AllocateBosonicField(Chi,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            !
            !Impurity Self-energy and Hartree contribution
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w")
               call read_Matrix(S_DMFT%N_s,reg(PrevItFolder)//"Hartree_UNimp",paramagnet)
            endif
            !
            !Lattice Gf
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(Glat,reg(PrevItFolder),"Glat_w",Crystal%kpt)
               if(look4dens%TargetDensity.eq.0d0) Glat%mu = look4dens%mu
            else
               Glat%mu = look4dens%mu
               if((look4dens%TargetDensity.ne.0d0).and.Hmodel) call set_density(Glat%mu,Beta,Crystal,look4dens)
               call calc_Gmats(Glat,Crystal)
            endif
            call calc_density(Glat,Crystal,Glat%N_ks)
            call calc_density(Glat,Glat%N_s)
            !
            !
            !just a sanity check
            do ik=1,Glat%Nkpt
               call check_Hermiticity(Glat%N_ks(:,:,ik,1),eps,enforce=.false.,hardstop=.false.,name="Nlat_k"//str(ik)//"_s1",verb=.true.)
               if(.not.paramagnet)call check_Hermiticity(Glat%N_ks(:,:,ik,Nspin),eps,enforce=.false.,hardstop=.false.,name="Nlat_k"//str(ik)//"_s"//str(Nspin),verb=.true.)
            enddo
            !
            !
            !Logical Flags
            calc_Pk = .true.
            calc_Wfull = .true.
            calc_Sigmak = .true.
            if(ItStart.gt.0)then
               merge_P = .true.
               merge_Sigma = .true.
            endif
            !
      end select
      !
      !
      if(calc_Sigmak)then
         !
         if(allocated(VH))deallocate(VH)
         if(allocated(VH_Nlat))deallocate(VH_Nlat)
         if(allocated(VH_Nimp))deallocate(VH_Nimp)
         allocate(VH(Crystal%Norb,Crystal%Norb));VH=czero
         allocate(VH_Nlat(Crystal%Norb,Crystal%Norb));VH_Nlat=czero
         allocate(VH_Nimp(Crystal%Norb,Crystal%Norb));VH_Nimp=czero
         !
         if(allocated(Vxc))deallocate(Vxc)
         if(.not.Vxc_in)then
            allocate(Vxc(Crystal%Norb,Crystal%Norb,Crystal%Nkpt,Nspin))
            Vxc=czero
         endif
         !
      endif
      !
      calc_W = calc_Wedmft .or. calc_Wfull
      !
      dump_Chik = dump_Chik .and. calc_Pk .and. (reg(structure).ne."None")
      dump_Wk = dump_Wk .and. calc_W .and. (reg(structure).ne."None")
      !
      !
      !Allocate and initialize different density matrices
      allocate(densityLDA(Crystal%Norb,Crystal%Norb,Nspin));densityLDA=czero
      allocate(densityGW(Crystal%Norb,Crystal%Norb,Nspin));densityGW=czero
      allocate(densityDMFT(Crystal%Norb,Crystal%Norb,Nspin));densityDMFT=czero
      do isite=1,Nsite
         allocate(LocalOrbs(isite)%rho_OrbSpin(LocalOrbs(isite)%Norb,LocalOrbs(isite)%Norb,Nspin))
         LocalOrbs(isite)%rho_OrbSpin=0d0
      enddo
      if(ItStart.eq.0)then
         !
         densityLDA = Glat%N_s
         densityGW = Glat%N_s
         densityDMFT=czero
         !
         call dump_Matrix(densityLDA,reg(pathINPUT),"Nlda",paramagnet)
         !
      else
         !
         if(addTierIII)call read_Matrix(densityLDA,reg(pathINPUT)//"Nlda",paramagnet)
         call read_Matrix(densityDMFT,reg(PrevItFolder)//"Nimp",paramagnet)
         densityGW=Glat%N_s
         !
         do isite=1,Nsite
            file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/Nqmc.DAT"
            call inquireFile(reg(file),filexists,hardstop=.false.,verb=verbose)
            if(filexists)then
               unit = free_unit()
               open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
               read(unit,*) muQMC
               do ib1=1,LocalOrbs(isite)%Nflavor
                  iorb = (ib1+mod(ib1,2))/2
                  ispin = abs(mod(ib1,2)-2)
                  read(unit,*) idum,LocalOrbs(isite)%rho_OrbSpin(iorb,iorb,ispin)
               enddo
               close(unit)
            endif
            if(ExpandImpurity.or.AFMselfcons)exit
         enddo
         !
      endif
      !
      write(*,"(A,F)") new_line("A")//"     Lattice chemical potential:  ",Glat%mu
      if(ItStart.gt.0)then
         write(*,"(A,F)") "     Impurity chemical potential: ",muQMC
         write(*,"(A,F)") "     Difference Glat%mu-muQMC:    ",Glat%mu-muQMC
      endif
      !
   end subroutine initialize_Fields


   !---------------------------------------------------------------------------!
   !PURPOSE: Estimate the impurity self-energy for the 0th iteration
   !this is the the local G0W0 contribution if present, the Hartree term otherwise
   !---------------------------------------------------------------------------!
   subroutine calc_SigmaGuess()
      !
      implicit none
      integer                               :: iorb,jorb,korb,lorb
      integer                               :: ispin,ib1,ib2,iw
      real(8),allocatable                   :: Uinst_0th(:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_SigmaGuess"
      !
      !
      if(.not.solve_DMFT) stop "calc_SigmaGuess: no guess needed if DMFT is not performed."
      if(FirstIteration.ne.0) stop "calc_SigmaGuess: this subroutine works only in the 0th iteration."
      if((reg(CalculationType).eq."GW+EDMFT").and.(.not.S_G0W0%status).and.addTierIII)then
         call AllocateFermionicField(S_G0W0,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_Sigma_spex(SpexVersion,S_G0W0,Crystal,verbose,RecomputeG0W0,Vxc)
      endif
      !
      !This has to be fixed because it works only for calculations with Wlat stored
      allocate(Uinst_0th(Crystal%Norb**2,Crystal%Norb**2));Uinst_0th=0d0
      if(Ustart)then
         Uinst_0th = dreal(Ulat%screened_local(:,:,1))
      else
         Uinst_0th = dreal(Wlat%screened_local(:,:,1))
      endif
      !
      !Initial inpurity self-energy guess
      !This differs from the one done later because the both
      !S_DMFT and Nlda can contain (real) off-diag elements.
      do ispin=1,Nspin
         do iorb=1,S_DMFT%Norb
            do jorb=1,S_DMFT%Norb
               do korb=1,S_DMFT%Norb
                  do lorb=1,S_DMFT%Norb
                     !
                     call F2Bindex(S_DMFT%Norb,[iorb,jorb],[korb,lorb],ib1,ib2)
                     !
                     S_DMFT%N_s(iorb,jorb,ispin) = S_DMFT%N_s(iorb,jorb,ispin) + Uinst_0th(ib1,ib2)*dreal(Glat%N_s(korb,lorb,ispin))
                     !
                  enddo
              enddo
            enddo
         enddo
      enddo
      call dump_Matrix(S_DMFT%N_s,reg(ItFolder),"Hartree0",paramagnet)
      deallocate(Uinst_0th)
      !
      !local projection of G0W0 self-energy otherwise just Hartree
      if(reg(CalculationType).eq."GW+EDMFT")then
         S_DMFT%ws = S_G0W0%ws
      else
         do iw=1,S_DMFT%Npoints
            S_DMFT%ws(:,:,iw,:) = S_DMFT%N_s
         enddo
      endif
      !
   end subroutine calc_SigmaGuess


   !---------------------------------------------------------------------------!
   !PURPOSE: Join the all the component of the self-energy.
   !---------------------------------------------------------------------------!
   subroutine join_SigmaFull(Iteration)
      !
      implicit none
      integer,intent(in)                    :: Iteration
      integer                               :: ik,iw,ispin,iorb,jorb
      real(8),allocatable                   :: Z_qpsc(:,:,:)
      complex(8),allocatable                :: Vxc_loc(:,:,:)
      type(FermionicField)                  :: S_EMB

      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- join_SigmaFull"
      !
      !
      call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta,mu=Glat%mu)
      !
      select case(reg(CalculationType))
         case default
            !
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
            !
         case("G0W0","scGW","GW+EDMFT")
            !
            !
            if(.not.S_G0W0%status) stop "join_SigmaFull: S_G0W0 not properly initialized."
            if((Iteration.gt.0).and.(.not.S_GW%status)) stop "join_SigmaFull: S_GW not properly initialized."
            if((Iteration.gt.0).and.addTierIII.and.(.not.S_G0W0dc%status)) stop "join_SigmaFull: S_G0W0dc not properly initialized."
            if(.not.allocated(VH)) stop "join_SigmaFull: VH not allocated."
            if(.not.allocated(Vxc))then !stop "join_SigmaFull: Vxc not allocated."
               if(addTierIII)write(*,"(A)")"     Allocating empty Vxc ( enclosed in SigmaGoWo(k,iw) )"
               allocate(Vxc(Crystal%Norb,Crystal%Norb,Crystal%Nkpt,Nspin))
               Vxc=czero
               allocate(Vxc_loc(S_Full%Norb,S_Full%Norb,Nspin));Vxc_loc=czero
               call read_Matrix(Vxc_loc,reg(pathINPUT)//"Vxc",paramagnet)
            else
               allocate(Vxc_loc(S_Full%Norb,S_Full%Norb,Nspin));Vxc_loc=czero
               do ik=1,S_Full%Nkpt
                  Vxc_loc = Vxc_loc + Vxc(:,:,ik,:)/S_Full%Nkpt
               enddo
               if(verbose)call dump_matrix(Vxc_loc,reg(ItFolder),"Vxc",paramagnet)
            endif
            !
            !
            !Check that the difference between G0W0 and scDGdc is locally causal
            if(addTierIII.and.(Iteration.gt.0))then
               S_G0W0%wks = S_G0W0%wks - S_G0W0dc%wks
               call check_S_G0W0()
            endif
            !
            !
            do ispin=1,Nspin
               !
               do ik=1,S_Full%Nkpt
                  do iw=1,S_Full%Npoints
                     !
                     if(Iteration.eq.0)then
                        !
                        !Keep only single-shot GW: G0W0 - Vxc + VH
                        if(addTierIII)then
                           S_Full%wks(:,:,iw,ik,ispin) = S_G0W0%wks(:,:,iw,ik,ispin) - Vxc(:,:,ik,ispin) + VH(:,:)
                        else
                           S_Full%wks(:,:,iw,ik,ispin) = S_G0W0%wks(:,:,iw,ik,ispin)
                        endif
                        !
                     else
                        !
                        !Put together all the terms: G0W0(already with removed DC) + scGW(already containing Simp)
                        if(addTierIII)then
                           S_Full%wks(:,:,iw,ik,ispin) = S_G0W0%wks(:,:,iw,ik,ispin) + S_GW%wks(:,:,iw,ik,ispin) - Vxc(:,:,ik,ispin) + VH(:,:)
                        else
                           S_Full%wks(:,:,iw,ik,ispin) = S_GW%wks(:,:,iw,ik,ispin)
                        endif
                        !
                     endif
                     !
                  enddo
               enddo
               !
               if(paramagnet)then
                  S_Full%wks(:,:,:,:,Nspin) = S_Full%wks(:,:,:,:,1)
                  exit
               endif
               !
            enddo
            call FermionicKsum(S_Full)
            !
            do ispin=1,Nspin
               S_Full%N_s(:,:,ispin) = Vxc_loc(:,:,ispin) - VH
            enddo
            !
            deallocate(VH,Vxc,Vxc_loc)
            call DeallocateFermionicField(S_G0W0)
            call DeallocateFermionicField(S_G0W0dc)
            if(.not.causal_D)call DeallocateFermionicField(S_GW)
            !
            !
         case("DMFT+statU","DMFT+dynU","EDMFT")
            !
            !
            if(.not.S_DMFT%status) stop "join_SigmaFull: S_DMFT not properly initialized."
            !
            !Put together all the terms
            do ispin=1,Nspin
               do ik=1,S_Full%Nkpt
                  do iw=1,S_Full%Npoints
                     S_Full%wks(:,:,iw,ik,ispin) = S_DMFT%ws(:,:,iw,ispin) - HartreeFact*S_DMFT%N_s(:,:,int(Nspin/ispin))
                  enddo
               enddo
               S_Full%N_s(:,:,ispin) = HartreeFact*S_DMFT%N_s(:,:,int(Nspin/ispin))
               if(paramagnet)then
                  S_Full%wks(:,:,:,:,Nspin) = S_Full%wks(:,:,:,:,1)
                  exit
               endif
            enddo
            call FermionicKsum(S_Full)
            !
            !
      end select
      !
      !
      !Add constant embedding self-energy
      select case(reg(Embedding))
         case default
            !
            stop "Available Embedding self-energies types are: loc (filename: Semb_w_s[1,2].DAT), nonloc (filename: Semb_w_k_s[1,2].DAT)."
            !
         case("None")
            !
            if(verbose)write(*,"(A)") "     No embedding self-energy included."
            !
         case("nonloc")
            !
            write(*,"(A)") "     Adding constant non-local embedding self-energy."
            call AllocateFermionicField(S_EMB,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            call read_FermionicField(S_EMB,reg(pathINPUTtr),"Semb_w",Crystal%kpt)
            !
            !Add non-local embedding
            S_Full%wks = S_Full%wks + S_EMB%wks
            call DeallocateFermionicField(S_EMB)
            call FermionicKsum(S_Full)
            !
         case("loc")
            !
            write(*,"(A)") "     Adding constant local embedding self-energy."
            call AllocateFermionicField(S_EMB,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            call read_FermionicField(S_EMB,reg(pathINPUTtr),"Semb_w")
            !
            !Add local embedding to all k-points
            do ik=1,S_Full%Nkpt
               S_Full%wks(:,:,:,ik,:) = S_Full%wks(:,:,:,ik,:) + S_EMB%ws
            enddo
            call DeallocateFermionicField(S_EMB)
            call FermionicKsum(S_Full)
            !
      end select
      !
      !
      !Compute local quasiparticle weigth in the Wannier basis for the full self-energy
      allocate(Z_qpsc(S_Full%Norb,S_Full%Norb,Nspin));Z_qpsc=0d0
      do ispin=1,Nspin
         do iorb=1,S_Full%Norb
            do jorb=1,S_Full%Norb
               Z_qpsc(iorb,jorb,ispin) = 1d0 / (1d0 + abs(dimag(S_Full%ws(iorb,jorb,1,ispin)))*S_Full%Beta/pi)
            enddo
         enddo
      enddo
      call dump_Matrix(Z_qpsc,reg(ItFolder),"Z_qpsc",paramagnet)
      deallocate(Z_qpsc)
      !
      !
      !
   contains
      !
      !
      !
      subroutine check_S_G0W0()
         !
         implicit none
         integer                               :: iw,ik,iorb,ispin,isite
         real(8),allocatable                   :: wmats_orig(:)
         real(8)                               :: ImS_1,ImS_2
         logical                               :: causal_G0W0_loc
         type(FermionicField)                  :: S_G0W0_DMFT
         type(FermionicField)                  :: S_G0W0_imp
         !
         call FermionicKsum(S_G0W0)
         !
         !Compute the G0W0 contribution to the local self-energy with removed DC
         call AllocateFermionicField(S_G0W0_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta,mu=Glat%mu)
         do isite=1,Nsite
            !
            !Extract the local G0W0 self-energy for each site
            call AllocateFermionicField(S_G0W0_imp,LocalOrbs(isite)%Norb,Nmats,Beta=Beta)
            call loc2imp(S_G0W0_imp,S_G0W0,LocalOrbs(isite)%Orbs)
            !
            !Put it into an object that contains only the site indexes
            call imp2loc(S_G0W0_DMFT,S_G0W0_imp,isite,LocalOrbs,.false.,.false.,.false.,name="S_G0W0_imp")
            !
            call DeallocateField(S_G0W0_imp)
            !
         enddo
         !
         !Check for the causality in the G0W0 contribution to the local self-energy
         !at all frequencies since it's not done in check_S_G0W0
         causal_G0W0_loc = GoWoDC_loc
         if(causal_G0W0_loc)then
            causaloop: do ispin=1,Nspin
               do iorb=1,S_G0W0_DMFT%Norb
                  do iw=1,S_G0W0_DMFT%Npoints
                     ImS_1 = dimag(S_G0W0_DMFT%ws(iorb,iorb,iw,ispin))
                     if(ImS_1.gt.0d0)then
                        write(*,"(A)")"     Warning: the local G0W0 self-energy has been found non-causal at iw="//str(iw)//" iorb="//str(iorb)//" ispin="//str(ispin)
                        causal_G0W0_loc=.false.
                        exit causaloop
                     endif
                     if(iw.le.10)then
                        ImS_2 = dimag(S_G0W0_DMFT%ws(iorb,iorb,iw+1,ispin))
                        if(ImS_2.gt.ImS_1) write(*,"(A)")"     Warning: the local G0W0 self-energy seems not to scale as a Fermi-liquid. If Delta(tau) is non-causal try to set G0W0DC_LOC=F."
                     endif
                  enddo
               enddo
            enddo causaloop
         endif
         !
         !From the S_G0W0^{SPEX}_{ij} + S_G0W0^{SPEX}_{i} - S_G0W0^{DC}_{ij} - S_G0W0^{DC}_{i}
         !we remove [ S_G0W0^{SPEX}_{i} - S_G0W0^{DC}_{i} ]
         !here, if Vxc is inside S_G0W0, also the local contribution from Vxc is removed
         if((.not.GoWoDC_loc).or.(.not.causal_G0W0_loc))then
            write(*,"(A)")"     Local G0W0-scGW_DC self-energy removed."
            do ik=1,S_G0W0%Nkpt
               S_G0W0%wks(:,:,:,ik,:) = S_G0W0%wks(:,:,:,ik,:) - S_G0W0_DMFT%ws
            enddo
         else
            write(*,"(A)")"     Local G0W0-scGW_DC kept."
         endif
         call DeallocateField(S_G0W0_DMFT)
         !
      endsubroutine check_S_G0W0
      !
      !
      !
   end subroutine join_SigmaFull


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the correction to the self-consistency equations
   !         for non-local calculations. See arxiv:2011.05311
   !---------------------------------------------------------------------------!
   subroutine calc_causality_Delta_correction()
      !
      implicit none
      integer                               :: ik,iw,ispin
      complex(8),allocatable                :: GS(:,:),SG(:,:),SGS(:,:)
      complex(8),allocatable                :: invG(:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_causality_Delta_correction"
      !
      !
      if(.not.S_GW%status) stop "calc_causality_Delta_correction: S_GW not properly initialized."
      if(.not.Glat%status) stop "calc_causality_Delta_correction: Glat not properly initialized."
      !
      call AllocateFermionicField(Delta_correction,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
      allocate(GS(S_GW%Norb,S_GW%Norb));GS=czero
      allocate(SG(S_GW%Norb,S_GW%Norb));SG=czero
      allocate(SGS(S_GW%Norb,S_GW%Norb));SGS=czero
      allocate(invG(S_GW%Norb,S_GW%Norb));invG=czero
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(S_GW,Glat,Delta_correction),&
      !$OMP PRIVATE(ispin,iw,ik,GS,SG,SGS,invG)
      !$OMP DO
      do ispin=1,Nspin
         do iw=1,S_GW%Npoints
            !
            GS=czero;SG=czero;SGS=czero
            !
            do ik=1,S_GW%Nkpt
               !
               GS = GS + matmul(Glat%wks(:,:,iw,ik,ispin),S_GW%wks(:,:,iw,ik,ispin))/S_GW%Nkpt
               SG = SG + matmul(S_GW%wks(:,:,iw,ik,ispin),Glat%wks(:,:,iw,ik,ispin))/S_GW%Nkpt
               SGS = SGS + matmul(S_GW%wks(:,:,iw,ik,ispin),matmul(Glat%wks(:,:,iw,ik,ispin),S_GW%wks(:,:,iw,ik,ispin)))/S_GW%Nkpt
               !
            enddo
            !
            invG=czero
            invG = Glat%ws(:,:,iw,ispin)
            call inv(invG)
            !
            ! Put toghether all the pieces. arxiv:2011.05311 Eq.7
            Delta_correction%ws(:,:,iw,ispin) = + SGS                             &
                                            - matmul(SG,matmul(invG,GS))      &
                                            + 2d0*S_GW%ws(:,:,iw,ispin)       & !<-- this is S_DMFT + all the local shifts
                                            - matmul(SG,invG)                 &
                                            - matmul(invG,GS)
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(GS,SG,SGS,invG)
      !
      call DeallocateFermionicField(S_GW)
      !
      call symmetrize_GW(Delta_correction,EqvGWndx)
      call dump_FermionicField(Delta_correction,reg(ItFolder),"Delta_correction_w",paramagnet)
      !
   end subroutine calc_causality_Delta_correction
   !
   subroutine calc_causality_curlyU_correction(mode)
      !
      implicit none
      character(len=*),intent(in)           :: mode
      integer                               :: iq,iw,Norb
      complex(8),allocatable                :: P(:,:),W(:,:)
      complex(8),allocatable                :: WP(:,:),PW(:,:),PWP(:,:)
      complex(8),allocatable                :: W_P(:,:),P_W(:,:)
      complex(8),allocatable                :: invWa(:,:),invWb(:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_causality_curlyU_correction"
      !
      !
      if(.not.Plat%status) stop "calc_causality_curlyU_correction: Plat not properly initialized."
      if(.not.Wlat%status) stop "calc_causality_curlyU_correction: Wlat not properly initialized."
      !
      call AllocateBosonicField(curlyU_correction,Crystal%Norb,Nmats,Crystal%iq_gamma,Beta=Beta,no_bare=.true.)
      !
      select case(reg(mode))
         case default
            !
            stop "calc_causality_curlyU_correction: Available types are: curlyU, Ploc."
            !
         case("curlyU")
            !
            Norb = int(sqrt(dble(Plat%Nbp)))
            allocate(P(Plat%Nbp,Plat%Nbp));P=czero
            allocate(W(Plat%Nbp,Plat%Nbp));W=czero
            allocate(WP(Plat%Nbp,Plat%Nbp));WP=czero
            allocate(PW(Plat%Nbp,Plat%Nbp));PW=czero
            allocate(PWP(Plat%Nbp,Plat%Nbp));PWP=czero
            allocate(W_P(Plat%Nbp,Plat%Nbp));W_P=czero
            allocate(P_W(Plat%Nbp,Plat%Nbp));P_W=czero
            allocate(invWa(Plat%Nbp,Plat%Nbp));invWa=czero
            allocate(invWb(Plat%Nbp,Plat%Nbp));invWb=czero
            !
            !$OMP PARALLEL DEFAULT(PRIVATE),&
            !$OMP SHARED(Plat,Wlat,curlyU_correction)
            !$OMP DO
            do iw=1,Plat%Npoints
               !
               WP=czero;PW=czero;PWP=czero
               W_P=czero;P_W=czero
               invWa=czero;invWb=czero
               !
               do iq=1,Plat%Nkpt
                  !
                  P = Plat%screened(:,:,iw,iq)
                  W = Wlat%screened(:,:,iw,iq)
                  !
                  WP = WP + matmul(W,P)/Plat%Nkpt
                  PW = PW + matmul(P,W)/Plat%Nkpt
                  PWP = PWP + matmul(P,matmul(W,P))/Plat%Nkpt
                  !
               enddo
               !
               P = Plat%screened_local(:,:,iw)
               W = Wlat%screened_local(:,:,iw)
               !
               W_P = matmul(W,P)
               P_W = matmul(P,W)
               !
               invWa = P + PWP
               call inv(invWa)
               !
               invWb = P + matmul(P,matmul(W,P))
               call inv(invWb)
               !
               ! Put toghether all the pieces. arxiv:2011.05311 Eq.22
               curlyU_correction%screened_local(:,:,iw) = matmul(PW,matmul(invWa,WP)) - matmul(P_W,matmul(invWb,W_P))
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(P,W,WP,PW,PWP,W_P,P_W,invWa,invWb)
            !
            call isReal(curlyU_correction)
            !
            do iw=1,curlyU_correction%Npoints
               call check_Symmetry(curlyU_correction%screened_local(:,:,iw),1e7*eps,enforce=.true.,hardstop=.false.,name="curlyU_correction_w"//str(iw),verb=.true.)
            enddo
            !
         case("Ploc")
            !
            allocate(WP(Plat%Nbp,Plat%Nbp));WP=czero
            allocate(PW(Plat%Nbp,Plat%Nbp));PW=czero
            allocate(PWP(Plat%Nbp,Plat%Nbp));PWP=czero
            allocate(invWa(Plat%Nbp,Plat%Nbp));invWa=czero
            !
            !$OMP PARALLEL DEFAULT(PRIVATE),&
            !$OMP SHARED(Plat,Wlat,curlyU_correction)
            !$OMP DO
            do iw=1,Plat%Npoints
               !
               WP=czero;PW=czero;PWP=czero;invWa=czero
               !
               do iq=1,Plat%Nkpt
                  !
                  WP = WP + matmul(Wlat%screened(:,:,iw,iq),Plat%screened(:,:,iw,iq))/Plat%Nkpt
                  PW = PW + matmul(dag(Plat%screened(:,:,iw,iq)),Wlat%screened(:,:,iw,iq))/Plat%Nkpt
                  PWP = PWP + matmul(dag(Plat%screened(:,:,iw,iq)),matmul(Wlat%screened(:,:,iw,iq),Plat%screened(:,:,iw,iq)))/Plat%Nkpt
                  !
               enddo
               !
               !This inversion is not possible for multi-orbital systems
               invWa = Wlat%screened_local(:,:,iw)
               call inv(invWa)
               !
               ! Put toghether all the pieces. Jiyu equation.
               curlyU_correction%screened_local(:,:,iw) = + PWP - matmul(PW,matmul(invWa,WP))   &
                                                          - matmul(PW,invWa) - matmul(invWa,WP) &
                                                          + dag(Plat%screened_local(:,:,iw)) + Plat%screened_local(:,:,iw)
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(WP,PW,PWP,invWa)
            !
      end select
      !
      call symmetrize_GW(curlyU_correction,EqvGWndx)
      call dump_BosonicField(curlyU_correction,reg(ItFolder),"curlyU_correction_w.DAT")
      !
   end subroutine calc_causality_curlyU_correction


   !---------------------------------------------------------------------------!
   !PURPOSE: compute, dump and fit the diagonal hybridization function for each
   !         spin-orbital flavor
   !---------------------------------------------------------------------------!
   subroutine calc_Delta(isite,Iteration)
      !
      implicit none
      integer,intent(in)                    :: isite
      integer,intent(in)                    :: Iteration
      !
      type(FermionicField)                  :: Gloc
      type(FermionicField)                  :: SigmaImp
      type(FermionicField)                  :: DeltaCorr
      type(FermionicField)                  :: FermiPrint
      type(FermionicField)                  :: DeltaOld,DeltaSym
      integer                               :: Norb,Norb_MultiTier,unit
      integer                               :: ispin,iorb,iw
      integer                               :: itau,ndx,wndx
      real(8),allocatable                   :: wmats(:),tau(:),Moments(:,:,:)
      real(8),allocatable                   :: Eloc(:,:),ElocOld(:,:),Eloc_s(:,:,:),PrintLine(:),coef01(:,:)
      real(8)                               :: CrystalField,taup
      complex(8),allocatable                :: zeta(:,:,:),invG(:,:)
      complex(8),allocatable                :: Dfit(:,:,:),Dmats(:,:,:),Ditau(:,:,:)
      complex(8),allocatable                :: invCurlyG(:,:,:),invCurlyG_MultiTier(:,:,:)
      real(8),allocatable                   :: tauF(:),ReadLine(:),DitauOld(:,:,:)
      character(len=255)                    :: file,MomDir
      logical                               :: filexists
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Delta of "//reg(LocalOrbs(isite)%Name)
      !
      !
      if(.not.S_DMFT%status) stop "calc_Delta: S_DMFT not properly initialized."
      if(.not.Glat%status) stop "calc_Delta: Glat not properly initialized."
      if(causal_D.and.(.not.Delta_correction%status)) stop "calc_Delta: requested causality correction but Delta_correction not properly initialized."
      if(.not.allocated(LocalOrbs)) stop "calc_Delta: LocalOrbs not properly initialized."
      !
      Norb = LocalOrbs(isite)%Norb
      !
      allocate(Eloc(Norb,Nspin));Eloc=0d0
      allocate(coef01(Norb,Nspin));coef01=0d0
      allocate(invCurlyG(Norb,Nmats,Nspin));invCurlyG=czero
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(Beta,Nmats)
      allocate(zeta(Norb,Nmats,Nspin))
      do iorb=1,Norb
         do iw=1,Nmats
            zeta(iorb,iw,:) = dcmplx( Glat%mu , wmats(iw) )
         enddo
      enddo
      !
      if(MultiTier.and.(Nsite.eq.1))then
         !
         !The Dyson equation occurs inside all the local space. Rotations not allowed here.
         Norb_MultiTier = Glat%Norb
         allocate(invCurlyG_MultiTier(Norb_MultiTier,Nmats,Nspin));invCurlyG_MultiTier=czero
         allocate(invG(Norb_MultiTier,Norb_MultiTier));invG=czero
         do ispin=1,Nspin
            do iw=1,Nmats
               !
               invG = Glat%ws(:,:,iw,ispin)
               call inv(invG)
               !
               do iorb=1,Norb_MultiTier
                  !
                  invCurlyG_MultiTier(iorb,iw,ispin) = invG(iorb,iorb) + S_DMFT%ws(iorb,iorb,iw,ispin)
                  if(causal_D) invCurlyG(iorb,iw,ispin) = invG(iorb,iorb) + S_DMFT%ws(iorb,iorb,iw,ispin) - Delta_correction%ws(iorb,iorb,iw,ispin)
                  !
               enddo
               !
            enddo
         enddo
         deallocate(invG)
         !
         !LocalOrbs are extracted
         do iorb=1,Norb
            invCurlyG(iorb,:,:) = invCurlyG_MultiTier(LocalOrbs(isite)%Orbs(iorb),:,:)
         enddo
         deallocate(invCurlyG_MultiTier)
         !
      else
         !
         !The Dyson equation occurs directly on LocalOrbs.
         call AllocateFermionicField(SigmaImp,Norb,Nmats,Beta=Beta)
         call AllocateFermionicField(Gloc,Norb,Nmats,Beta=Beta)
         if(causal_D)call AllocateFermionicField(DeltaCorr,Norb,Nmats,Beta=Beta)
         !
         !Extract and rotate from local (non-diagonal) to imp (diagonal) the given sites
         call clear_attributes(Gloc)
         call clear_attributes(SigmaImp)
         if(causal_D)call clear_attributes(DeltaCorr)
         if(RotateHloc)then
            !
            call loc2imp(Gloc,Glat,LocalOrbs(isite)%Orbs,U=LocalOrbs(isite)%Rot)
            call loc2imp(SigmaImp,S_DMFT,LocalOrbs(isite)%Orbs,U=LocalOrbs(isite)%Rot)
            if(causal_D)call loc2imp(DeltaCorr,Delta_correction,LocalOrbs(isite)%Orbs,U=LocalOrbs(isite)%Rot)
            !
         else
            !
            call loc2imp(Gloc,Glat,LocalOrbs(isite)%Orbs)
            call loc2imp(SigmaImp,S_DMFT,LocalOrbs(isite)%Orbs)
            if(causal_D)call loc2imp(DeltaCorr,Delta_correction,LocalOrbs(isite)%Orbs)
            !
         endif
         !
         !Compute the fermionic Weiss field invCurlyG.
         allocate(invG(Norb,Norb));invG=czero
         do ispin=1,Nspin
            do iw=1,Nmats
               !
               invG = Gloc%ws(:,:,iw,ispin)
               call inv(invG)
               !
               do iorb=1,Norb
                  !
                  if(causal_D) then
                     !self-consistency is only on SigmaXC: invCurlyG(iorb,iw,ispin) = invG(iorb,iorb) + ( SigmaImp%ws(iorb,iorb,iw,ispin) - SigmaImp%N_s(iorb,iorb,ispin) ) - DeltaCorr%ws(iorb,iorb,iw,ispin)
                     invCurlyG(iorb,iw,ispin) = invG(iorb,iorb) + SigmaImp%ws(iorb,iorb,iw,ispin) - DeltaCorr%ws(iorb,iorb,iw,ispin)
                  else
                     !self-consistency is only on SigmaXC: invCurlyG(iorb,iw,ispin) = invG(iorb,iorb) + ( SigmaImp%ws(iorb,iorb,iw,ispin) - SigmaImp%N_s(iorb,iorb,ispin) )
                     invCurlyG(iorb,iw,ispin) = invG(iorb,iorb) + SigmaImp%ws(iorb,iorb,iw,ispin)
                  endif
                  !
               enddo
               !
            enddo
         enddo
         deallocate(invG)
         !Recalculate and replace the Hartree term with the lattice densities would go here
         !
         !Print the impurity fields used to compute delta
         call dump_Matrix(Gloc%N_s,reg(ItFolder),"Solver_"//reg(LocalOrbs(isite)%Name)//"/Nloc_"//reg(LocalOrbs(isite)%Name),paramagnet)
         call dump_FermionicField(Gloc,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Glat_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
         if(causal_D)call dump_FermionicField(DeltaCorr,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Delta_correction_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
         call dump_FermionicField(SigmaImp,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Simp_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
         call DeallocateFermionicField(SigmaImp)
         call DeallocateFermionicField(Gloc)
         !
      endif
      !
      !Extract the local energy
      allocate(Dfit(Norb,Nmats,Nspin));Dfit=czero
      Dfit = zeta - invCurlyG
      !
      select case(reg(DeltaFit))
         case default
            !
            stop "Available modes for Delta fitting: Inf, Analytic, Moments."
            !
         case("Analytic")
            !
            file = "DeltaPara_"//reg(LocalOrbs(isite)%Name)//".DAT"
            MomDir = reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/"
            wndx = minloc(abs(wmats-0.85*wmatsMax),dim=1)
            call fit_Delta(Dfit,Beta,Nfit,reg(MomDir),reg(file),"Shifted",Eloc,filename="DeltaAnd",Wlimit=wndx,coef01=coef01)
            !
         case("Moments","Inf")
            !
            file = "DeltaMom_"//reg(LocalOrbs(isite)%Name)//".DAT"
            MomDir = reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/"
            !
            allocate(Moments(Norb,Nspin,0:min(SigmaMaxMom,Nfit)));Moments=0d0
            wndx = minloc(abs(wmats-0.85*wmatsMax),dim=1)
            call fit_moments(Dfit,Beta,reg(MomDir),reg(file),"Sigma",Moments,filename="DeltaMom",Wlimit=wndx)
            !
            if(reg(reg(DeltaFit)).eq."Inf") Eloc = dreal(Dfit(:,Nmats,:))
            if(reg(reg(DeltaFit)).eq."Moments") Eloc = Moments(:,:,0)
            coef01 = Moments(:,:,1)
            !
            deallocate(Moments)
            !
      end select
      deallocate(Dfit)
      !
      !Compute Delta on matsubara
      allocate(Dmats(Norb,Nmats,Nspin));Dmats=czero
      do ispin=1,Nspin
         do iorb=1,Norb
            !
            do iw=1,Nmats
               Dmats(iorb,iw,ispin) = dcmplx( Glat%mu , wmats(iw) ) - Eloc(iorb,ispin) - invCurlyG(iorb,iw,ispin)
            enddo
            !
         enddo
      enddo
      deallocate(invCurlyG)
      !
      !Mixing local energy and optionally Delta(iw)
      if((Mixing_Delta.gt.0d0).and.(Iteration.gt.0))then
         !
         if(.not.Mixing_Delta_tau)then
            write(*,"(A)")"     Mixing Delta(iw) with "//str(Mixing_Delta,3)//" of old solution."
            if(.not.Mixing_Delta_tau)then
               write(*,"(A)")"     Mixing Delta(iw)."
               call AllocateFermionicField(DeltaOld,Norb,Nmats,Beta=Beta)
               call read_FermionicField(DeltaOld,reg(MixItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Delta_"//reg(LocalOrbs(isite)%Name)//"_w")
               do ispin=1,Nspin
                  do iw=1,Nmats
                     do iorb=1,Norb
                        Dmats(iorb,iw,ispin) = (1d0-Mixing_Delta)*Dmats(iorb,iw,ispin) + Mixing_Delta*DeltaOld%ws(iorb,iorb,iw,ispin)
                     enddo
                  enddo
               enddo
               call DeallocateFermionicField(DeltaOld)
            endif
         endif
         !
         write(*,"(A)")"     Mixing Eloc with "//str(Mixing_Delta,3)//" of old solution."
         allocate(ElocOld(Norb,Nspin));ElocOld=0d0
         file = reg(MixItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/Eloc.DAT"
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="read")
         read(unit,*)
         do iorb=1,Norb
            read(unit,"(2E20.12)") (ElocOld(iorb,ispin),ispin=1,Nspin)
         enddo
         close(unit)
         Eloc = (1d0-Mixing_Delta)*Eloc + Mixing_Delta*ElocOld
         deallocate(ElocOld)
         !
      endif
      !
      !Symmetrizations - If sets are defined Delta and Eloc are always symmetrized
      if(sym_mode.gt.0)then
         !
         call AllocateFermionicField(DeltaSym,Norb,Nmats,Beta=Beta)
         do iorb=1,Norb
            DeltaSym%ws(iorb,iorb,:,:) = Dmats(iorb,:,:)
         enddo
         if(RotateHloc)then
            call symmetrize_imp(DeltaSym,LocalOrbs(isite)%Eig)
            if(causal_D)call symmetrize_imp(DeltaCorr,LocalOrbs(isite)%Eig)
         else
            call symmetrize_GW(DeltaSym,EqvImpndx(isite))
            if(causal_D)call symmetrize_GW(DeltaCorr,EqvImpndx(isite))
         endif
         do iorb=1,Norb
            Dmats(iorb,:,:) = DeltaSym%ws(iorb,iorb,:,:)
         enddo
         call DeallocateFermionicField(DeltaSym)
         !
         allocate(Eloc_s(Norb,Norb,Nspin));Eloc_s=0d0
         do ispin=1,Nspin
            Eloc_s(:,:,ispin) = diag(Eloc(:,ispin))
         enddo
         if(RotateHloc)then
            call symmetrize_imp(Eloc_s,LocalOrbs(isite)%Eig)
         else
            call symmetrize_GW(Eloc_s,EqvImpndx(isite))
         endif
         do ispin=1,Nspin
            Eloc(:,ispin) = diagonal(Eloc_s(:,:,ispin))
         enddo
         deallocate(Eloc_s)
         !
      endif
      !
      if(EqvGWndx%para.eq.1)then
         Dmats(:,:,1) = (Dmats(:,:,1) + Dmats(:,:,Nspin))/2d0
         Dmats(:,:,Nspin) = Dmats(:,:,1)
      endif
      !
      !Fourier transform to the tau axis
      allocate(Ditau(Norb,Solver%NtauF,Nspin));Ditau=0d0
      do ispin=1,Nspin
         !add correction that allows to use the tail correction in the FT
         do iorb=1,Norb
            if(reg(CalculationType).ne."GW+EDMFT") coef01(iorb,ispin) = abs(dimag(Dmats(iorb,Nmats,ispin)))*wmats(Nmats)
            write(*,"(A,E10.3)")"     First coeff orb: "//str(iorb)//" spin: "//str(ispin),coef01(iorb,ispin)
            Dmats(iorb,:,ispin) = Dmats(iorb,:,ispin)/coef01(iorb,ispin)
         enddo
         !FT
         call Fmats2itau_vec(Beta,Dmats(:,:,ispin),Ditau(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
         !remove correction
         do iorb=1,Norb
            Dmats(iorb,:,ispin) = Dmats(iorb,:,ispin)*coef01(iorb,ispin)
            Ditau(iorb,:,ispin) = Ditau(iorb,:,ispin)*coef01(iorb,ispin)
         enddo
      enddo
      !
      !Check on Delta(tau) causality
      do ispin=1,Nspin
         do iorb=1,Norb
            do itau=1,Solver%NtauF
               if(dreal(Ditau(iorb,itau,ispin)).gt.0d0)then
                  write(*,"(A,E10.3)")"     Warning: Removing non-causality from Delta(tau) at orb: "//str(iorb)//" spin: "//str(ispin)//" itau: "//str(itau)
                  if(dreal(Ditau(iorb,int(Solver%NtauF/2),ispin)).lt.0d0)then
                     Ditau(iorb,itau,ispin) = Ditau(iorb,int(Solver%NtauF/2),ispin)
                  else
                     Ditau(iorb,itau,ispin) = czero
                  endif
               endif
            enddo
         enddo
      enddo
      !
      !Mixing Delta(tau)
      if((Mixing_Delta.gt.0d0).and.(Iteration.gt.0).and.Mixing_Delta_tau)then
         !
         write(*,"(A)")"     Mixing Delta(tau) with "//str(Mixing_Delta,3)//" of old solution."
         !
         !reading old Delta for each site
         allocate(tauF(Solver%NtauF_in));tauF = linspace(0d0,Beta,Solver%NtauF_in)
         allocate(DitauOld(Norb,Solver%NtauF_in,Nspin));DitauOld=czero
         allocate(ReadLine(LocalOrbs(isite)%Nflavor))
         file = reg(MixItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/Delta_t.DAT"
         call inquireFile(reg(file),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
         do itau=1,Solver%NtauF_in
            ReadLine=0d0
            read(unit,*) taup,ReadLine
            if(abs(taup-tauF(itau)).gt.eps) stop "Reading Delta_t.DAT from previous iteration: Impurity fermionic tau mesh does not coincide."
            !
            ndx=1
            do iorb=1,Norb
               do ispin=1,Nspin
                  DitauOld(iorb,itau,ispin) = dcmplx(ReadLine(ndx),0d0)
                  ndx=ndx+1
               enddo
            enddo
            !
         enddo
         deallocate(ReadLine,tauF)
         !
         Ditau = (1d0-Mixing_Delta)*Ditau + Mixing_Delta*DitauOld
         !
         deallocate(DitauOld)
         !
         !Recompute Delta(iw) because CurlyG(iw) needs to be updated
         Dmats=czero
         do ispin=1,Nspin
            call Fitau2mats_vec(Beta,Ditau(:,:,ispin),Dmats(:,:,ispin),tau_uniform=.true.)
         enddo
         !
      endif
      !
      !Write Eloc and chemical potential
      file = reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/Eloc.DAT"
      unit = free_unit()
      open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(1E20.12)") Glat%mu
      CrystalField=0d0
      do iorb=1,Norb
         if(addCF) CrystalField = LocalOrbs(isite)%CrystalField(iorb)
         write(unit,"(2E20.12)")  Eloc(iorb,1)-EqvGWndx%hseed+CrystalField, Eloc(iorb,Nspin)+EqvGWndx%hseed+CrystalField
      enddo
      close(unit)
      !
      !Write Delta(tau)
      allocate(tau(Solver%NtauF));tau=0d0
      tau = linspace(0d0,Beta,Solver%NtauF)
      allocate(PrintLine(LocalOrbs(isite)%Nflavor))
      file = reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/Delta_t.DAT"
      unit = free_unit()
      open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
      do itau=1,Solver%NtauF
         ndx=1
         PrintLine=0d0
         do iorb=1,Norb
            do ispin=1,Nspin
               PrintLine(ndx) = real(Ditau(iorb,itau,ispin))
               ndx=ndx+1
            enddo
         enddo
         write(unit,"(2000E20.12)") tau(itau),PrintLine
      enddo
      close(unit)
      deallocate(PrintLine)
      !
      !fields that are going to be needed in the following iterations
      call AllocateFermionicField(FermiPrint,Norb,Nmats,Beta=Beta)
      !
      !CurlyG(iw) - recomputed from the eventually mixed Delta
      call clear_attributes(FermiPrint)
      do ispin=1,Nspin
         do iw=1,Nmats
            do iorb=1,Norb
               FermiPrint%ws(iorb,iorb,iw,ispin) = 1d0/(dcmplx( Glat%mu , wmats(iw) ) - Eloc(iorb,ispin) - Dmats(iorb,iw,ispin))
            enddo
         enddo
      enddo
      FermiPrint%mu=Glat%mu
      call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","G0_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
      !
      !Delta(iw)
      call clear_attributes(FermiPrint)
      do iorb=1,Norb
         FermiPrint%ws(iorb,iorb,:,:) = Dmats(iorb,:,:)
      enddo
      FermiPrint%mu=Glat%mu
      call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Delta_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
      !
      !Delta(iw) without causality correction
      if(causal_D)then
         call clear_attributes(FermiPrint)
         do ispin=1,Nspin
            do iorb=1,Norb
               FermiPrint%ws(iorb,iorb,:,ispin) = Dmats(iorb,:,ispin) - DeltaCorr%ws(iorb,iorb,:,ispin)
            enddo
         enddo
         call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Delta_notCorr_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
         call DeallocateFermionicField(DeltaCorr)
      endif
      !
      !test of the fit and FT procedures
      if(verbose)then
         !
         Dmats=czero
         do ispin=1,Nspin
            call Fitau2mats_vec(Beta,Ditau(:,:,ispin),Dmats(:,:,ispin),tau_uniform=.true.)
         enddo
         !
         call clear_attributes(FermiPrint)
         do ispin=1,Nspin
            do iorb=1,Norb
               FermiPrint%ws(iorb,iorb,:,ispin) = Dmats(iorb,:,ispin)
            enddo
         enddo
         call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","testFT_Delta_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
         !
         call clear_attributes(FermiPrint)
         do ispin=1,Nspin
            do iorb=1,Norb
               FermiPrint%ws(iorb,iorb,:,ispin) = Eloc(iorb,ispin) + Dmats(iorb,:,ispin)
            enddo
         enddo
         call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/fits/","test_Eo+Delta_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
         !
      endif
      call DeallocateFermionicField(FermiPrint)
      !
      deallocate(Eloc,zeta,Dmats,Ditau,tau,wmats)
      !
   end subroutine calc_Delta


   !---------------------------------------------------------------------------!
   !PURPOSE: compute and dump the istantanous and retarded interactions. Since
   !         if(ExpandImpurity) I have to be sure that all the local interactions
   !         are the same, here I'm extracting the tensor for all the sites and
   !         using the isite input just to print out the wanted one.
   !---------------------------------------------------------------------------!
   subroutine calc_Interaction(isite,Iteration,checkInvariance)
      !
      implicit none
      integer,intent(in)                    :: isite
      integer,intent(in)                    :: Iteration
      logical,intent(in)                    :: checkInvariance
      !
      type(BosonicField)                    :: Wloc
      type(BosonicField)                    :: Pimp
      type(BosonicField)                    :: curlyU,curlyUold
      type(BosonicField)                    :: curlyUcorr
      type(physicalU)                       :: PhysicalUelements
      integer                               :: Norb,Norb_MultiTier,Nbp,unit
      integer                               :: ib1,ib2,itau,iw
      integer                               :: isitecheck
      real(8),allocatable                   :: Uinst(:,:),Ucheck(:,:)
      real(8),allocatable                   :: Kfunct(:,:,:),Kpfunct(:,:,:)
      real(8),allocatable                   :: ScreeningMat(:,:)
      character(len=255)                    :: file
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Interaction of "//reg(LocalOrbs(isite)%Name)
      !
      !
      if(.not.P_EDMFT%status) stop "calc_Interaction: P_EDMFT not properly initialized."
      if(.not.Wlat%status) stop "calc_Interaction: Wlat not properly initialized."
      if(causal_U.and.(.not.curlyU_correction%status)) stop "calc_Interaction: requested causality correction but curlyU_correction not properly initialized."
      if(.not.allocated(LocalOrbs)) stop "calc_Interaction: LocalOrbs not properly initialized."
      !
      Norb = LocalOrbs(isite)%Norb
      Nbp = Norb**2
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      allocate(Uinst(LocalOrbs(isite)%Nflavor,LocalOrbs(isite)%Nflavor));Uinst=0d0
      call AllocateBosonicField(curlyU,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
      !
      !get interaction
      select case(reg(CalculationType))
         case default
            !
            stop "If you got so far somethig is wrong."
            !
         case("DMFT+statU","DMFT+dynU")
            !
            call loc2imp(curlyU,Ulat,LocalOrbs(isite)%Orbs)
            if(RotateUloc) call TransformBosonicField(curlyU,LocalOrbs(isite)%Rot,PhysicalUelements%Full_Map)
            !
         case("EDMFT","GW+EDMFT")
            !
            if(Ustart)then
               !
               write(*,"(A)") "     Using local Ucrpa as effective interaction."
               call loc2imp(curlyU,Ulat,LocalOrbs(isite)%Orbs)
               if(RotateUloc) call TransformBosonicField(curlyU,LocalOrbs(isite)%Rot,PhysicalUelements%Full_Map)
               !
            else
               !
               if(MultiTier.and.(Nsite.eq.1))then
                  !
                  !The Dyson equation occurs inside all the local space. Rotations not allowed here.
                  write(*,"(A)") "     Computing the local effective interaction (MultiTier)."
                  Norb_MultiTier = int(sqrt(dble(Wlat%Nbp)))
                  call AllocateBosonicField(curlyU_EDMFT,Norb_MultiTier,Nmats,Crystal%iq_gamma,Beta=Beta)
                  !
                  if(causal_U)then
                     call calc_curlyU(curlyU_EDMFT,Wlat,P_EDMFT,curlyUcorr=curlyU_correction,mode=reg(causal_U_type))
                  else
                     call calc_curlyU(curlyU_EDMFT,Wlat,P_EDMFT)
                  endif
                  !
                  !LocalOrbs are extracted
                  call loc2imp(curlyU,curlyU_EDMFT,LocalOrbs(isite)%Orbs)
                  call DeallocateBosonicField(curlyU_EDMFT)
                  !
               else
                  !
                  !The Dyson equation occurs directly on LocalOrbs.
                  write(*,"(A)") "     Computing the local effective interaction."
                  call AllocateBosonicField(Pimp,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
                  call AllocateBosonicField(Wloc,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
                  !
                  call loc2imp(Pimp,P_EDMFT,LocalOrbs(isite)%Orbs)
                  call loc2imp(Wloc,Wlat,LocalOrbs(isite)%Orbs)
                  !
                  if(causal_U)then
                     call AllocateBosonicField(curlyUcorr,Norb,Nmats,Crystal%iq_gamma,Beta=Beta,no_bare=.true.)
                     call loc2imp(curlyUcorr,curlyU_correction,LocalOrbs(isite)%Orbs)
                     call calc_curlyU(curlyU,Wloc,Pimp,curlyUcorr=curlyUcorr,mode=reg(causal_U_type))
                  else
                     call calc_curlyU(curlyU,Wloc,Pimp)
                  endif
                  !
                  if(RotateUloc)then
                     call TransformBosonicField(curlyU,LocalOrbs(isite)%Rot,PhysicalUelements%Full_Map)
                     if(causal_U)then
                        call TransformBosonicField(curlyUcorr,LocalOrbs(isite)%Rot,PhysicalUelements%Full_Map)
                        call dump_BosonicField(curlyUcorr,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","curlyU_correction_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
                     endif
                  endif
                  !
                  call DeallocateBosonicField(curlyUcorr)
                  call DeallocateBosonicField(Pimp)
                  call DeallocateBosonicField(Wloc)
                  !
               endif
               !
            endif
            !
            !Mixing curlyU
            if((Mixing_curlyU.gt.0d0).and.(Iteration.gt.0))then
               write(*,"(A)")"     Mixing curlyU with "//str(Mixing_curlyU,3)//" of old solution."
               call dump_BosonicField(curlyU,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","curlyU_noMix_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
               call AllocateBosonicField(curlyUold,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
               call read_BosonicField(curlyUold,reg(MixItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","curlyU_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
               curlyU%bare_local = (1d0-Mixing_curlyU)*curlyU%bare_local + Mixing_curlyU*curlyUold%bare_local
               do iw=1,Nmats
                  curlyU%screened_local(:,:,iw) = (1d0-Mixing_curlyU)*curlyU%screened_local(:,:,iw) + Mixing_curlyU*curlyUold%screened_local(:,:,iw)
               enddo
               call DeallocateBosonicField(curlyUold)
            endif
            !
      end select
      !
      !Symmetrizations - If sets are defined curlyU is always symmetrized
      if(sym_mode.gt.0) call symmetrize_GW(curlyU,EqvImpndx(isite))
      !
      !Extract istantaneous interaction and screening function
      select case(reg(CalculationType))
         case("DMFT+statU")
            call calc_QMCinteractions(curlyU,Uinst)
         case("DMFT+dynU","EDMFT","GW+EDMFT")
            allocate(Kfunct(LocalOrbs(isite)%Nflavor,LocalOrbs(isite)%Nflavor,Solver%NtauB));Kfunct=0d0
            allocate(Kpfunct(LocalOrbs(isite)%Nflavor,LocalOrbs(isite)%Nflavor,Solver%NtauB));Kpfunct=0d0
            allocate(ScreeningMat(LocalOrbs(isite)%Nflavor,LocalOrbs(isite)%Nflavor));ScreeningMat=0d0
            call calc_QMCinteractions(curlyU,Uinst,Kfunct=Kfunct,Kpfunct=Kpfunct,Screening=ScreeningMat)
      end select
      !
      !Print curlyU in the solver basis
      call dump_BosonicField(curlyU,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","curlyU_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
      !
      !Istantaneous interaction
      call dump_Matrix(Uinst,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Umat.DAT")
      deallocate(Uinst)
      !
      !Print data for retarded interactions
      if(allocated(Kfunct))then
         call print_K(Kfunct,"K_t")
         deallocate(Kfunct)
      endif
      if(allocated(Kpfunct))then
         call print_K(Kpfunct,"Kp_t")
         deallocate(Kpfunct)
      endif
      !
      if(allocated(ScreeningMat))then
         call dump_Matrix(ScreeningMat,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Screening.DAT")
         deallocate(ScreeningMat)
      endif
      !
      !Check if the screened effective local interaction at iw=0 is the same for all the sites. Same Orbital dimension assumed.
      if(checkInvariance)then
         !
         do isitecheck=1,Nsite
            !
            Norb = LocalOrbs(isitecheck)%Norb
            !
            write(*,"(A)")"     Checking site "//reg(LocalOrbs(1)%Name)//"_"//str(isitecheck)
            write(*,"(A,10I3)")"     Norb: "//str(Norb)//" Orbitals: ",LocalOrbs(isitecheck)%Orbs
            !
            if(isitecheck.eq.1)allocate(Uinst(LocalOrbs(isitecheck)%Nflavor,LocalOrbs(isitecheck)%Nflavor))
            if(isitecheck.ne.1)allocate(Ucheck(LocalOrbs(isitecheck)%Nflavor,LocalOrbs(isitecheck)%Nflavor))
            !
            call clear_attributes(curlyU)
            !
            if(Ustart)then
               call loc2imp(curlyU,Ulat,LocalOrbs(isitecheck)%Orbs)
              !call isReal(curlyU)
            else
               call AllocateBosonicField(Pimp,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
               call AllocateBosonicField(Wloc,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
               call loc2imp(Pimp,P_EDMFT,LocalOrbs(isitecheck)%Orbs)
               call loc2imp(Wloc,Wlat,LocalOrbs(isitecheck)%Orbs)
               !
               call calc_curlyU(curlyU,Wloc,Pimp,sym=.false.)
               !
               call DeallocateBosonicField(Pimp)
               call DeallocateBosonicField(Wloc)
               !
            endif
            !
            if(isitecheck.eq.1)then
               call calc_QMCinteractions(curlyU,Uinst,sym=.false.)
               call dump_Matrix(Uinst,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isitecheck)%Name)//"/","Umat_noSym_"//str(isitecheck)//".DAT")
            else
               call calc_QMCinteractions(curlyU,Ucheck,sym=.false.)
               call dump_Matrix(Ucheck,reg(ItFolder)//"Solver_"//reg(LocalOrbs(isitecheck)%Name)//"/","Umat_noSym_"//str(isitecheck)//".DAT")
               !
               do ib1=1,LocalOrbs(isitecheck)%Nflavor
                  do ib2=1,LocalOrbs(isitecheck)%Nflavor
                     if(abs(Ucheck(ib1,ib2)-Uinst(ib1,ib2)).gt.1e-3)then
                        write(*,"(A,F,A,F)")"     Warning: Element["//str(ib1)//"]["//str(ib2)//"] is different:",Ucheck(ib1,ib2)," instead of: ",Uinst(ib1,ib2)
                     endif
                  enddo
               enddo
               deallocate(Ucheck)
               !
            endif
            !
         enddo !isitecheck
         deallocate(Uinst)
         !
      endif
      !
      call DeallocateBosonicField(curlyU)
      !
      !
      !
   contains
      !
      !
      !
      subroutine print_K(K,filename)
         implicit none
         real(8),intent(in)                 :: K(:,:,:)
         character(len=*),intent(in)        :: filename
         integer                            :: ndx
         real(8),allocatable                :: tau(:),PrintLine(:)
         !
         allocate(tau(Solver%NtauB));tau=0d0
         tau = linspace(0d0,Beta,Solver%NtauB)
         !
         file = reg(ItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/"//reg(filename)//".DAT"
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
         allocate(PrintLine(size(K,dim=1)*(size(K,dim=1)+1)/2));PrintLine=0d0
         do itau=1,Solver%NtauB
            ndx=1
            !print diagonal and LT
            do ib1=1,size(K,dim=1)
               do ib2=1,ib1
                  !
                  PrintLine(ndx) = K(ib1,ib2,itau)
                  if(Kdiag) PrintLine(ndx) = K(ib1,ib1,itau)
                  ndx=ndx+1
                  !
               enddo
            enddo
            write(unit,"(999E20.12)") tau(itau),PrintLine
         enddo
         deallocate(tau,PrintLine)
         close(unit)
         !
      end subroutine print_K
      !
      !
      !
   end subroutine calc_Interaction


   !---------------------------------------------------------------------------!
   !PURPOSE: fetch the impurity fields for each site solved by the QMC solver
   !         and plugs into a container with the same dimension of the lattice
   !         to facilitate the merge. If(ExpandImpurity) the impurity Gf is
   !         rotated back from diagonal to Wannier basis.
   !---------------------------------------------------------------------------!
   subroutine collect_QMC_results()
      !
      implicit none
      !
      integer                               :: iorb,jorb,ispin,jspin
      integer                               :: ib1,ib2,isite,idum
      integer                               :: unit,ndx,itau,iw,wndx,wndxOpt
      real(8)                               :: taup,muQMC
      real(8),allocatable                   :: tauF(:),tauB(:),wmats(:)
      real(8),allocatable                   :: ReadLine(:)
      real(8),allocatable                   :: Moments(:,:,:)
      character(len=255)                    :: file,MomDir
      logical                               :: filexists,fitSigmaTail
      !Impurity Green's function
      type(FermionicField),allocatable      :: Gimp(:),Fimp(:)
      complex(8),allocatable                :: Gitau(:,:,:),Fitau(:,:,:)
      complex(8),allocatable                :: Gmats(:,:),Fmats(:,:)
      !Impurity self-energy and fermionic Dyson equation
      type(FermionicField),allocatable      :: Simp(:)
      type(FermionicField)                  :: G0imp
      type(BosonicField)                    :: curlyU
      complex(8),allocatable                :: Sfit(:,:,:),SmatsTail(:)
      real(8),allocatable                   :: n0(:,:),Uinst(:,:)
      real(8)                               :: A(2),B(2)
      !Impurity susceptibilities
      real(8),allocatable                   :: nnt(:,:,:)
      complex(8),allocatable                :: NNitau(:,:,:,:,:)
      complex(8),allocatable                :: ChiMitau(:),ChiMmats(:)
      type(BosonicField)                    :: ChiCitau,ChiCmats
      !Impurity polarization and bosonic Dyson equation
      type(BosonicField)                    :: Pimp
      type(BosonicField)                    :: Wimp
      real(8),allocatable                   :: CDW(:,:),Integral(:,:),Z_dmft(:,:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- collect_QMC_results"
      !
      !
      if(.not.allocated(LocalOrbs)) stop "collect_QMC_results: LocalOrbs not properly initialized."
      !
      !Impurity observables initialization
      if(ExpandImpurity)then
         allocate(Gimp(1))
         allocate(Simp(1))
         if(Dyson_Imprvd_F)allocate(Fimp(1))
      else
         allocate(Gimp(Nsite))
         allocate(Simp(Nsite))
         if(Dyson_Imprvd_F)allocate(Fimp(Nsite))
      endif
      do isite=1,Nsite
         !
         write(*,"(A)") "     Setting up site: "//reg(LocalOrbs(isite)%Name)
         !
         if(allocated(LocalOrbs(isite)%rho_Flav))deallocate(LocalOrbs(isite)%rho_Flav)
         allocate(LocalOrbs(isite)%rho_Flav(LocalOrbs(isite)%Nflavor))
         LocalOrbs(isite)%rho_Flav=0d0
         !
         if(allocated(LocalOrbs(isite)%rho_OrbSpin))deallocate(LocalOrbs(isite)%rho_OrbSpin)
         allocate(LocalOrbs(isite)%rho_OrbSpin(LocalOrbs(isite)%Norb,LocalOrbs(isite)%Norb,Nspin))
         LocalOrbs(isite)%rho_OrbSpin=0d0
         !
         if(allocated(LocalOrbs(isite)%Docc))deallocate(LocalOrbs(isite)%Docc)
         allocate(LocalOrbs(isite)%Docc(LocalOrbs(isite)%Norb,LocalOrbs(isite)%Norb))
         LocalOrbs(isite)%Docc=0d0
         !
         call AllocateFermionicField(Gimp(isite),LocalOrbs(isite)%Norb,Nmats,Beta=Beta)
         call clear_attributes(Gimp(isite))
         !
         call AllocateFermionicField(Simp(isite),LocalOrbs(isite)%Norb,Nmats,Beta=Beta)
         call clear_attributes(Simp(isite))
         !
         if(Dyson_Imprvd_F)then
            call AllocateFermionicField(Fimp(isite),LocalOrbs(isite)%Norb,Nmats,Beta=Beta)
            call clear_attributes(Fimp(isite))
         endif
         !
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      fitSigmaTail = (ReplaceTail_Simp.ne.0d0).and.(ReplaceTail_Simp.lt.wmatsMax)
      !
      !
      ! COLLECT IMPURITY OCCUPATION --------------------------------------------
      allocate(densityDMFT(Crystal%Norb,Crystal%Norb,Nspin));densityDMFT=czero
      do isite=1,Nsite
         !
         write(*,"(A)") new_line("A")//"     Collecting occupation of site "//reg(LocalOrbs(isite)%Name)
         !
         !Read the impurity occupation
         file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/Nqmc.DAT"
         call inquireFile(reg(file),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
         read(unit,*) muQMC
         do ib1=1,LocalOrbs(isite)%Nflavor
            !
            read(unit,*) idum,LocalOrbs(isite)%rho_Flav(ib1)
            !
            iorb = (ib1+mod(ib1,2))/2
            ispin = abs(mod(ib1,2)-2)
            LocalOrbs(isite)%rho_OrbSpin(iorb,iorb,ispin) = LocalOrbs(isite)%rho_Flav(ib1)
            !
         enddo
         close(unit)
         !
         !Insert or Expand to the Lattice basis (this is just for printing)
         call imp2loc(densityDMFT,dcmplx(LocalOrbs(isite)%rho_OrbSpin,0d0),isite,LocalOrbs,ExpandImpurity,AFMselfcons,RotateHloc,name="Nimp")
         !
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      !symmetrize and print
      if(verbose)call dump_Matrix(densityDMFT,reg(PrevItFolder),"Nimp_noSym",paramagnet)
      call symmetrize_GW(densityDMFT,EqvGWndx)
      !
      call dump_Matrix(densityDMFT,reg(PrevItFolder),"Nimp",paramagnet)
      deallocate(densityDMFT)
      !
      !
      !
      !
      ! COLLECT IMPURITY GF ----------------------------------------------------
      call AllocateFermionicField(G_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
      do isite=1,Nsite
         !
         write(*,"(A)") new_line("A")//"     Collecting the impurity Green's function of site "//reg(LocalOrbs(isite)%Name)
         !
         !Read the impurity Green's function
         allocate(tauF(Solver%NtauF_in));tauF = linspace(0d0,Beta,Solver%NtauF_in)
         allocate(Gitau(LocalOrbs(isite)%Norb,Solver%NtauF_in,Nspin));Gitau=czero
         allocate(ReadLine(LocalOrbs(isite)%Nflavor))
         file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/Gimp_t.DAT"
         call inquireFile(reg(file),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
         do itau=1,Solver%NtauF_in
            ReadLine=0d0
            read(unit,*) taup,ReadLine
            if(abs(taup-tauF(itau)).gt.eps) stop "Reading Gimp_t.DAT from previous iteration: Impurity fermionic tau mesh does not coincide."
            !
            ndx=1
            do iorb=1,LocalOrbs(isite)%Norb
               do ispin=1,Nspin
                  Gitau(iorb,itau,ispin) = dcmplx(ReadLine(ndx),0d0)
                  ndx=ndx+1
               enddo
            enddo
            !
         enddo
         close(unit)
         deallocate(ReadLine,tauF)
         call dump_MaxEnt(Gitau,"itau",reg(PrevItFolder)//"Convergence/","Gqmc_"//reg(LocalOrbs(isite)%Name))
         !
         !FT to the matsubara axis
         allocate(Gmats(LocalOrbs(isite)%Norb,Nmats));Gmats=czero
         do ispin=1,Nspin
            call Fitau2mats_vec(Beta,Gitau(:,:,ispin),Gmats,tau_uniform=.true.)
            do iw=1,Nmats
               Gimp(isite)%ws(:,:,iw,ispin) = diag(Gmats(:,iw))
            enddo
            if(paramagnet)then
               Gimp(isite)%ws(:,:,:,Nspin) = Gimp(isite)%ws(:,:,:,1)
               exit
            endif
         enddo
         deallocate(Gitau,Gmats)
         call dump_FermionicField(Gimp(isite),reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Gimp_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
         !
         !Insert or Expand to the Lattice basis
         call imp2loc(G_DMFT,Gimp(isite),isite,LocalOrbs,ExpandImpurity,AFMselfcons,RotateHloc,name="Gimp")
         !
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      !symmetrize and print
      if(verbose)call dump_FermionicField(G_DMFT,reg(PrevItFolder),"Gimp_noSym_w",paramagnet)
      call symmetrize_GW(G_DMFT,EqvGWndx)
      !
      call dump_FermionicField(G_DMFT,reg(PrevItFolder),"Gimp_w",paramagnet)
      call dump_MaxEnt(G_DMFT,"mats",reg(PrevItFolder)//"Convergence/","Gimp",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call dump_MaxEnt(G_DMFT,"mats2itau",reg(PrevItFolder)//"Convergence/","Gimp",EqvGWndx%SetOrbs)
      call DeallocateFermionicField(G_DMFT)
      !
      !
      !
      !
      ! FERMIONIC DYSON EQUATION -----------------------------------------------
      call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
      !
      ! Check if the fermionic improved estimators are present
      if(Dyson_Imprvd_F)then
         do isite=1,Nsite
            file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/Fimp_S_t.DAT"
            call inquireFile(reg(file),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists)then
               write(*,"(A)") "     The static improved esitmators for site "//reg(LocalOrbs(isite)%Name)//" is missing. Switching to standard Dyson equation."
               Dyson_Imprvd_F=.false.
               deallocate(Fimp)
               exit
            endif
         enddo
      endif
      !
      ! Solve the Dyson equation for each site
      if(Dyson_Imprvd_F)then
         !
         ! Collect static improved estimators
         do isite=1,Nsite
            !
            write(*,"(A)") new_line("A")//"     Collecting the static impurity improved estimator of site "//reg(LocalOrbs(isite)%Name)
            !
            !Read the impurity Green's function
            allocate(tauF(Solver%NtauF_in));tauF = linspace(0d0,Beta,Solver%NtauF_in)
            allocate(Fitau(LocalOrbs(isite)%Norb,Solver%NtauF_in,Nspin));Fitau=czero
            allocate(ReadLine(LocalOrbs(isite)%Nflavor))
            file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/Fimp_S_t.DAT"
            call inquireFile(reg(file),filexists,verb=verbose)
            unit = free_unit()
            open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
            do itau=1,Solver%NtauF_in
               ReadLine=0d0
               read(unit,*) taup,ReadLine
               if(abs(taup-tauF(itau)).gt.eps) stop "Reading Fimp_S_t.DAT from previous iteration: Impurity fermionic tau mesh does not coincide."
               !
               ndx=1
               do iorb=1,LocalOrbs(isite)%Norb
                  do ispin=1,Nspin
                     Fitau(iorb,itau,ispin) = dcmplx(ReadLine(ndx),0d0)
                     ndx=ndx+1
                  enddo
               enddo
               !
            enddo
            close(unit)
            deallocate(ReadLine,tauF)
            call dump_MaxEnt(Fitau,"itau",reg(PrevItFolder)//"Convergence/","Fqmc_S_"//reg(LocalOrbs(isite)%Name))
            !
            !FT to the matsubara axis
            allocate(Fmats(LocalOrbs(isite)%Norb,Nmats));Fmats=czero
            do ispin=1,Nspin
               call Fitau2mats_vec(Beta,Fitau(:,:,ispin),Fmats,tau_uniform=.true.)
               do iw=1,Nmats
                  Fimp(isite)%ws(:,:,iw,ispin) = diag(Fmats(:,iw))
               enddo
               if(paramagnet)then
                  Fimp(isite)%ws(:,:,:,Nspin) = Fimp(isite)%ws(:,:,:,1)
                  exit
               endif
            enddo
            deallocate(Fitau,Fmats)
            call dump_FermionicField(Fimp(isite),reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Fimp_S_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
         ! Collect retarded improved estimators
         do isite=1,Nsite
            !
            write(*,"(A)") new_line("A")//"     Collecting the retarded impurity improved estimator of site "//reg(LocalOrbs(isite)%Name)
            !
            !Read the impurity Green's function
            allocate(tauF(Solver%NtauF_in));tauF = linspace(0d0,Beta,Solver%NtauF_in)
            allocate(Fitau(LocalOrbs(isite)%Norb,Solver%NtauF_in,Nspin));Fitau=czero
            allocate(ReadLine(LocalOrbs(isite)%Nflavor))
            file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/Fimp_R_t.DAT"
            call inquireFile(reg(file),filexists,verb=verbose,hardstop=.false.)
            if(.not.filexists)then
               write(*,"(A)")"     File "//reg(file)//" not found, using only static improved estimator."
               cycle
            endif
            unit = free_unit()
            open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
            do itau=1,Solver%NtauF_in
               ReadLine=0d0
               read(unit,*) taup,ReadLine
               if(abs(taup-tauF(itau)).gt.eps) stop "Reading Fimp_R_t.DAT from previous iteration: Impurity fermionic tau mesh does not coincide."
               !
               ndx=1
               do iorb=1,LocalOrbs(isite)%Norb
                  do ispin=1,Nspin
                     Fitau(iorb,itau,ispin) = dcmplx(ReadLine(ndx),0d0)
                     ndx=ndx+1
                  enddo
               enddo
               !
            enddo
            close(unit)
            deallocate(ReadLine)
            !
            !Linear interpolation of the tau boundaries
            do iorb=1,LocalOrbs(isite)%Norb
               do ispin=1,Nspin
                  if(tauF(1).eq.0d0) then
                     A = [tauF(2),dreal(Fitau(iorb,2,ispin))]
                     B = [tauF(3),dreal(Fitau(iorb,3,ispin))]
                     Fitau(iorb,1,ispin) = dcmplx(linear_interp(A,B,0d0),0d0)
                  endif
                  if(tauF(Solver%NtauF_in).eq.Beta) then
                     A = [tauF(Solver%NtauF_in-1),dreal(Fitau(iorb,Solver%NtauF_in-1,ispin))]
                     B = [tauF(Solver%NtauF_in-2),dreal(Fitau(iorb,Solver%NtauF_in-2,ispin))]
                     Fitau(iorb,Solver%NtauF_in,ispin) = dcmplx(linear_interp(A,B,Beta),0d0)
                  endif
               enddo
            enddo
            deallocate(tauF)
            call dump_MaxEnt(Fitau,"itau",reg(PrevItFolder)//"Convergence/","Fqmc_R_"//reg(LocalOrbs(isite)%Name))
            !
            !FT to the matsubara axis and add to the existing field
            allocate(Fmats(LocalOrbs(isite)%Norb,Nmats));Fmats=czero
            do ispin=1,Nspin
               call Fitau2mats_vec(Beta,Fitau(:,:,ispin),Fmats,tau_uniform=.true.)
               do iw=1,Nmats
                  Fimp(isite)%ws(:,:,iw,ispin) = Fimp(isite)%ws(:,:,iw,ispin) + diag(Fmats(:,iw))
               enddo
               if(paramagnet)then
                  Fimp(isite)%ws(:,:,:,Nspin) = Fimp(isite)%ws(:,:,:,1)
                  exit
               endif
            enddo
            deallocate(Fitau,Fmats)
            call dump_FermionicField(Fimp(isite),reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Fimp_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
         !Fermionic Dyson equation in the solver basis (always diagonal)
         do isite=1,Nsite
            !
            write(*,"(A)") new_line("A")//"     Solving fermionic Dyson of site "//reg(LocalOrbs(isite)%Name)//" using improved estimator."
            do ispin=1,Nspin
               do iorb=1,LocalOrbs(isite)%Norb
                  Simp(isite)%ws(iorb,iorb,:,ispin) = Fimp(isite)%ws(iorb,iorb,:,ispin) / Gimp(isite)%ws(iorb,iorb,:,ispin)
               enddo
               if(paramagnet)then
                  Simp(isite)%ws(:,:,:,Nspin) = Simp(isite)%ws(:,:,:,1)
                  exit
               endif
            enddo
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      else
         !
         do isite=1,Nsite
            !
            write(*,"(A)") new_line("A")//"     Collecting curlyG of site "//reg(LocalOrbs(isite)%Name)
            !
            !Read curlyG
            call AllocateFermionicField(G0imp,LocalOrbs(isite)%Norb,Nmats,Beta=Beta)
            call read_FermionicField(G0imp,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","G0_"//reg(LocalOrbs(isite)%Name)//"_w")
            !
            !Adjust with the chemical potential if the solver has changed it
            if((abs(G0imp%mu-muQMC).gt.1e-5).and.update_curlyG)then
               write(*,"(A)") new_line("A")//"     Updating the chemical potential of curlyG from "//str(G0imp%mu)//" to "//str(muQMC)
               do ispin=1,Nspin
                  do iw=1,Nmats
                     do iorb=1,LocalOrbs(isite)%Norb
                        G0imp%ws(iorb,iorb,iw,ispin) = 1d0/(1d0/G0imp%ws(iorb,iorb,iw,ispin) - G0imp%mu + muQMC)
                     enddo
                  enddo
                  if(paramagnet)then
                     G0imp%ws(:,:,:,Nspin) = G0imp%ws(:,:,:,1)
                     exit
                  endif
               enddo
               G0imp%mu=muQMC
               call dump_FermionicField(G0imp,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","G0_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
            endif
            !
            !Fermionic Dyson equation in the solver basis (always diagonal)
            write(*,"(A)") new_line("A")//"     Solving fermionic Dyson of site "//reg(LocalOrbs(isite)%Name)
            do ispin=1,Nspin
               do iorb=1,LocalOrbs(isite)%Norb
                  Simp(isite)%ws(iorb,iorb,:,ispin) = 1d0/G0imp%ws(iorb,iorb,:,ispin) - 1d0/Gimp(isite)%ws(iorb,iorb,:,ispin)
               enddo
               if(paramagnet)then
                  Simp(isite)%ws(:,:,:,Nspin) = Simp(isite)%ws(:,:,:,1)
                  exit
               endif
            enddo
            call DeallocateFermionicField(G0imp)
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      endif
      !
      !Self-energy manipulations
      do isite=1,Nsite
         !
         !Fit the impurity self-energy tail
         if(fitSigmaTail)then
            !
            call dump_FermionicField(Simp(isite),reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","S_noFit_"//reg(LocalOrbs(isite)%Name)//"_w",paramagnet)
            !
            write(*,"(A)") new_line("A")//"     Fitting self-energy moments of site "//reg(LocalOrbs(isite)%Name)
            !
            !define the frequency index from which substitute the tail
            allocate(wmats(Nmats));wmats=FermionicFreqMesh(Beta,Nmats)
            wndx = minloc(abs(wmats-ReplaceTail_Simp),dim=1)
            wndxOpt = wndx
            do ispin=1,Nspin
               do iorb=1,LocalOrbs(isite)%Norb
                  wndxOpt = min(wndxOpt,minloc(abs(Simp(isite)%ws(iorb,iorb,:,ispin)-0.3d0*minval(dimag(Simp(isite)%ws(iorb,iorb,:,ispin)))),dim=1))
               enddo
            enddo
            !
            !perform the fit on the diagonal
            file = "SimpMom_"//reg(LocalOrbs(isite)%Name)//".DAT"
            MomDir = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/"
            !
            allocate(Moments(LocalOrbs(isite)%Norb,Nspin,0:min(SigmaMaxMom,Nfit)));Moments=0d0
            allocate(Sfit(LocalOrbs(isite)%Norb,Nmats,Nspin))
            do ispin=1,Nspin
               do iorb=1,LocalOrbs(isite)%Norb
                  Sfit(iorb,:,ispin) = Simp(isite)%ws(iorb,iorb,:,ispin)
               enddo
            enddo
            call fit_moments(Sfit,Beta,reg(MomDir),reg(file),"Sigma",Moments,filename="Simp",Wlimit=wndx)
            !
            allocate(SmatsTail(Nmats));SmatsTail=czero
            write(*,"(A,F)") new_line("A")//"     Replacing Sigma tail starting from iw_["//str(wndx)//"]=",wmats(wndx)
            if(abs(wmats(wndx)-wmats(wndxOpt)).gt.1d0)write(*,"(A,F)") "     Optimal would be iw_["//str(wndxOpt)//"]=",wmats(wndxOpt)
            do ispin=1,Nspin
               do iorb=1,LocalOrbs(isite)%Norb
                  SmatsTail = S_Moments(Moments(iorb,ispin,:),wmats)
                  Simp(isite)%ws(iorb,iorb,wndx:Nmats,ispin) = SmatsTail(wndx:Nmats)
               enddo
               if(paramagnet)then
                  Simp(isite)%ws(:,:,:,Nspin) = Simp(isite)%ws(:,:,:,1)
                  exit
               endif
            enddo
            deallocate(wmats,Sfit,Moments,SmatsTail)
            !
         endif
         !
         !Fill up the N_s attribute that correspond to the Hartree term of the
         !self-energy always diagonal in the solver basis
         call AllocateBosonicField(curlyU,LocalOrbs(isite)%Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
         call read_BosonicField(curlyU,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","curlyU_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
         select case(reg(HartreeType))
            case default
               stop "collect_QMC_results: Available HartreeT types: GW or DMFT."
            case("GW")
               !
               Simp(isite)%N_s = czero
               do ispin=1,Nspin
                  do iorb=1,LocalOrbs(isite)%Norb
                     do jorb=1,LocalOrbs(isite)%Norb
                        call F2Bindex(LocalOrbs(isite)%Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                        Simp(isite)%N_s(iorb,iorb,ispin) = Simp(isite)%N_s(iorb,iorb,ispin) + curlyU%screened_local(ib1,ib2,1)*LocalOrbs(isite)%rho_OrbSpin(jorb,jorb,ispin)
                        !
                     enddo
                  enddo
               enddo
               call DeallocateBosonicField(curlyU)
               !
               !The magnetization will be given only by the self-energy beyond Hartree
               Simp(isite)%N_s(:,:,1) = (Simp(isite)%N_s(:,:,1)+Simp(isite)%N_s(:,:,Nspin))
               Simp(isite)%N_s(:,:,Nspin) = Simp(isite)%N_s(:,:,1)
               !
            case("DMFT")
               !
               allocate(Uinst(LocalOrbs(isite)%Nflavor,LocalOrbs(isite)%Nflavor));Uinst=0d0
               call calc_QMCinteractions(curlyU,Uinst)
               call DeallocateBosonicField(curlyU)
               !
               Simp(isite)%N_s = czero
               do ib1=1,LocalOrbs(isite)%Nflavor
                  iorb = (ib1+mod(ib1,2))/2
                  ispin = abs(mod(ib1,2)-2)
                  do ib2=1,LocalOrbs(isite)%Nflavor
                     if(ib1.eq.ib2) cycle !Simp(isite)%N_s(iorb,iorb,ispin) = Simp(isite)%N_s(iorb,iorb,ispin) + Uinst(2*iorb-1,2*iorb)*LocalOrbs(isite)%rho_Flav(ib2)
                     Simp(isite)%N_s(iorb,iorb,ispin) = Simp(isite)%N_s(iorb,iorb,ispin) + Uinst(ib1,ib2)*LocalOrbs(isite)%rho_Flav(ib2)
                  enddo
               enddo
               deallocate(Uinst)
               !
               !The magnetization will be given only by the self-energy beyond Hartree
               Simp(isite)%N_s(:,:,1) = (Simp(isite)%N_s(:,:,1)+Simp(isite)%N_s(:,:,Nspin))/2d0
               Simp(isite)%N_s(:,:,Nspin) = Simp(isite)%N_s(:,:,1)
               !
         end select
         !
         !Hartree term in the Solver basis
         call dump_Matrix(Simp(isite)%N_s,reg(PrevItFolder),"Solver_"//reg(LocalOrbs(isite)%Name)//"/Hartree_UNimp_"//reg(LocalOrbs(isite)%Name),paramagnet)
         !
         !Insert or Expand to the Lattice basis
         if(RotateHloc.and.(sym_mode.gt.1))call symmetrize_imp(Simp(isite),LocalOrbs(isite)%Eig)
         call imp2loc(S_DMFT,Simp(isite),isite,LocalOrbs,ExpandImpurity,AFMselfcons,RotateHloc,name="Simp")
         !
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      !symmetrize and print
      if(verbose)call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_noSym_w",paramagnet)
      call symmetrize_GW(S_DMFT,EqvGWndx)
      !
      call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w",paramagnet)
      call dump_MaxEnt(S_DMFT,"mats",reg(PrevItFolder)//"Convergence/","Simp",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call dump_Matrix(S_DMFT%N_s,reg(PrevItFolder),"Hartree_UNimp",paramagnet)
      !
      !Compute local EDMFT quasiparticle weigth in the Wannier basis
      allocate(Z_dmft(S_DMFT%Norb,S_DMFT%Norb,Nspin));Z_dmft=0d0
      do ispin=1,Nspin
         do iorb=1,S_DMFT%Norb
            do jorb=1,S_DMFT%Norb
               Z_dmft(iorb,jorb,ispin) = 1d0 / (1d0 + abs(dimag(S_DMFT%ws(iorb,jorb,1,ispin)))*S_DMFT%Beta/pi)
            enddo
         enddo
      enddo
      call dump_Matrix(Z_dmft,reg(PrevItFolder),"Z_dmft",paramagnet)
      deallocate(Z_dmft)
      call DeallocateFermionicField(S_DMFT)
      !
      !
      !
      !
      ! COLLECT IMPURITY CHARGE SUSCEPTIBILITY AND BOSONIC DYSON ---------------
      if(bosonicSC)then
         !
         call AllocateBosonicField(C_EDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
         call AllocateBosonicField(curlyU_EDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
         call AllocateBosonicField(P_EDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
         call AllocateBosonicField(W_EDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
         !
         do isite=1,Nsite
            !
            write(*,"(A)") new_line("A")//"     Collecting the impurity susceptibilities of site "//reg(LocalOrbs(isite)%Name)
            !
            !Read the impurity N(tau)N(0)
            allocate(tauB(Solver%NtauB_in));tauB=0d0
            tauB = linspace(0d0,Beta,Solver%NtauB_in)
            allocate(nnt(LocalOrbs(isite)%Nflavor,LocalOrbs(isite)%Nflavor,Solver%NtauB_in));nnt=0d0
            allocate(ReadLine(LocalOrbs(isite)%Nflavor*(LocalOrbs(isite)%Nflavor+1)/2))
            file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/nn_t.DAT"
            call inquireFile(reg(file),filexists,verb=verbose)
            unit = free_unit()
            open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
            do itau=1,Solver%NtauB_in
               !
               ReadLine=0d0
               read(unit,*) taup,ReadLine
               if(abs(taup-tauB(itau)).gt.eps) stop "Reading nn_t.DAT from previous iteration: Impurity bosonic tau mesh does not coincide."
               !
               !this is cumbersome but better be safe
               ndx=1
               do ib1=1,LocalOrbs(isite)%Nflavor
                  do ib2=1,ib1
                     nnt(ib1,ib2,itau) = ReadLine(ndx)
                     if(ib1.ne.ib2)nnt(ib2,ib1,itau) = ReadLine(ndx)
                     ndx=ndx+1
                  enddo
               enddo
               !
            enddo
            close(unit)
            deallocate(ReadLine)
            !
            !Reshape N(tau)N(0)for easier handling
            allocate(NNitau(LocalOrbs(isite)%Norb,LocalOrbs(isite)%Norb,Nspin,Nspin,Solver%NtauB_in));NNitau=czero
            do ib1=1,LocalOrbs(isite)%Nflavor
               do ib2=1,LocalOrbs(isite)%Nflavor
                  !
                  iorb = (ib1+mod(ib1,2))/2
                  jorb = (ib2+mod(ib2,2))/2
                  ispin = abs(mod(ib1,2)-2)
                  jspin = abs(mod(ib2,2)-2)
                  !
                  call halfbeta_symm(nnt(ib1,ib2,:))
                  NNitau(iorb,jorb,ispin,jspin,:) = dcmplx(nnt(ib1,ib2,:),0d0)
                  !
               enddo
            enddo
            deallocate(nnt)
            !
            !
            !Spin susceptibility------------------------------------------------
            !ChiM(tau) = Sum_ab <S(tau)S(0)> in the solver basis
            allocate(ChiMitau(Solver%NtauB_in));ChiMitau=czero
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,LocalOrbs(isite)%Norb
                     do jorb=1,LocalOrbs(isite)%Norb
                        !
                        ChiMitau = ChiMitau + NNitau(iorb,jorb,ispin,jspin,:)*(-1d0)**(ispin-jspin)/4d0
                        !
                     enddo
                  enddo
               enddo
            enddo
            call dump_BosonicField(ChiMitau,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","ChiM_"//reg(LocalOrbs(isite)%Name)//"_t.DAT",tauB)
            !
            !FT to get ChiM(iw)
            allocate(wmats(Nmats));wmats=BosonicFreqMesh(Beta,Nmats)
            allocate(ChiMmats(Nmats));ChiMmats=czero
            call Bitau2mats(Beta,ChiMitau,ChiMmats,tau_uniform=.true.)
            call dump_BosonicField(ChiMmats,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","ChiM_"//reg(LocalOrbs(isite)%Name)//"_w.DAT",wmats)
            deallocate(ChiMitau,ChiMmats,wmats)
            !
            !
            !Charge susceptibility----------------------------------------------
            !TEST>>>
            !Collect n(tau=0) instead of the density as it is obtained with the same statistics
            allocate(n0(LocalOrbs(isite)%Norb,Nspin));n0=0d0
            allocate(ReadLine(LocalOrbs(isite)%Nflavor))
            file = reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/resultsQMC/n_t.DAT"
            call inquireFile(reg(file),filexists,verb=verbose)
            unit = free_unit()
            open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
            ReadLine=0d0
            read(unit,*) taup,ReadLine
            if(abs(taup).gt.0d0) stop "Reading n_t.DAT from previous iteration: first tau point is not zero."
            ndx=1
            do iorb=1,LocalOrbs(isite)%Norb
               do ispin=1,Nspin
                  n0(iorb,ispin) = ReadLine(ndx)
                  ndx=ndx+1
               enddo
            enddo
            close(unit)
            deallocate(ReadLine)
            !>>>TEST
            !
            !ChiC(tau) = Sum_s1s2 <Na(tau)Nb(0)> - <Na><Nb> in the solver basis
            call AllocateBosonicField(ChiCitau,LocalOrbs(isite)%Norb,Solver%NtauB_in,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            do iorb=1,LocalOrbs(isite)%Norb
               do jorb=1,LocalOrbs(isite)%Norb
                  !
                  call F2Bindex(LocalOrbs(isite)%Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                  !
                  do ispin=1,Nspin
                     do jspin=1,Nspin
                        !
                        if(Test_flag_4)then
                           ChiCitau%screened_local(ib1,ib2,:) = ChiCitau%screened_local(ib1,ib2,:) + ( NNitau(iorb,jorb,ispin,jspin,:) - n0(iorb,ispin) * n0(jorb,jspin) )
                           !( NNitau(iorb,jorb,ispin,jspin,:) - LocalOrbs(isite)%rho_OrbSpin(iorb,iorb,ispin) * LocalOrbs(isite)%rho_OrbSpin(jorb,jorb,jspin) )
                        else
                           ChiCitau%screened_local(ib1,ib2,:) = ChiCitau%screened_local(ib1,ib2,:) + &
                           ( NNitau(iorb,jorb,ispin,jspin,:) - LocalOrbs(isite)%rho_OrbSpin(iorb,iorb,ispin) * LocalOrbs(isite)%rho_OrbSpin(jorb,jorb,jspin) )
                        endif
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
            deallocate(n0)
            !
            !
            !Double occupations-------------------------------------------------
            if(Nspin.gt.1)then
               do iorb=1,LocalOrbs(isite)%Norb
                  do jorb=1,LocalOrbs(isite)%Norb
                     LocalOrbs(isite)%Docc(iorb,jorb) = ( NNitau(iorb,jorb,1,2,1)+ NNitau(iorb,jorb,2,1,1) ) / 2d0
                  enddo
               enddo
               call dump_Matrix(LocalOrbs(isite)%Docc,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","Dimp_"//reg(LocalOrbs(isite)%Name)//".DAT")
            endif
            deallocate(NNitau)
            !
            !TEST>>>
            if(Test_flag_1)then
               !
               do ib1=1,ChiCitau%Nbp
                  do ib2=ib1,ChiCitau%Nbp
                     !
                     if(ib1.eq.ib2)then
                        taup = minval(dreal(ChiCitau%screened_local(ib1,ib2,:)))
                        if(taup.lt.0d0) ChiCitau%screened_local(ib1,ib2,:) = ChiCitau%screened_local(ib1,ib2,:) - taup
                     else
                        taup = maxval(dreal(ChiCitau%screened_local(ib1,ib2,:)))
                        if(taup.gt.0d0) ChiCitau%screened_local(ib1,ib2,:) = ChiCitau%screened_local(ib1,ib2,:) - taup
                     endif
                     !
                     ChiCitau%screened_local(ib2,ib1,:) = ChiCitau%screened_local(ib1,ib2,:)
                     !
                  enddo
               enddo
               !
            endif
            !
            allocate(Integral(ChiCitau%Nbp,ChiCitau%Nbp));Integral=0d0
            do itau=1,ChiCitau%Npoints-1
               Integral = Integral + ( ChiCitau%screened_local(:,:,itau) + ChiCitau%screened_local(:,:,itau+1) ) * 0.5d0 * abs(tauB(itau+1)-tauB(itau))
            enddo
            !>>>TEST
            !
            !User-defined modification of local charge susceptibility
            if(alphaChi.ne.1d0)then
               write(*,"(A)") new_line("A")//"     Rescaling ChiC(tau) with factor: "//str(alphaChi,3)
               ChiCitau%screened_local = ChiCitau%screened_local * alphaChi
            endif
            if(ChiDiag)then
               write(*,"(A)") new_line("A")//"     Removing all off-diagonal elements from ChiC(tau)."
               do itau=1,ChiCitau%Npoints
                  ChiCitau%screened_local(:,:,itau) = diag(diagonal(ChiCitau%screened_local(:,:,itau)))
               enddo
            endif
            !
            call dump_BosonicField(ChiCitau,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","ChiC_"//reg(LocalOrbs(isite)%Name)//"_t.DAT",axis=tauB)
            deallocate(tauB)
            !
            !FT to get ChiC(iw)
            call AllocateBosonicField(ChiCmats,LocalOrbs(isite)%Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            call Bitau2mats(Beta,ChiCitau%screened_local,ChiCmats%screened_local,tau_uniform=.true.)
            call DeallocateBosonicField(ChiCitau)
            call isReal(ChiCmats)
            !
            !TEST>>>
            call dump_Matrix(dreal(Integral-ChiCmats%screened_local(:,:,1)),reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","IntDiff_"//reg(LocalOrbs(isite)%Name)//".DAT")
            if(Test_flag_2)then
               ChiCmats%screened_local(:,:,1) = Integral
            endif
            deallocate(Integral)
            !>>>TEST
            !
            !Remove the iw=0 divergency of local charge susceptibility
            if(removeCDW_C)then
               write(*,"(A)") new_line("A")//"     Divergency removal in ChiC(iw=0) of site "//reg(LocalOrbs(isite)%Name)
               call dump_BosonicField(ChiCmats,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","ChiC_CDW_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
               allocate(CDW(ChiCmats%Nbp,ChiCmats%Nbp));CDW=0d0
               CDW = real(ChiCmats%screened_local(:,:,1))
               call remove_CDW(ChiCmats,"imp")
               CDW = CDW - real(ChiCmats%screened_local(:,:,1))
               call dump_Matrix(CDW,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","CDW_"//reg(LocalOrbs(isite)%Name)//".DAT")
               deallocate(CDW)
            endif
            !
            !Print ChiC
            call dump_BosonicField(ChiCmats,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","ChiC_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
            !
            !
            !BOSONIC DYSON -----------------------------------------------------
            !
            !recollect curlyU
            call AllocateBosonicField(curlyU,LocalOrbs(isite)%Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
            call read_BosonicField(curlyU,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/","curlyU_"//reg(LocalOrbs(isite)%Name)//"_w.DAT")
            !
            !Bosonic Dyson equation in the solver basis
            write(*,"(A)") new_line("A")//"     Solving bosonic Dyson of site "//reg(LocalOrbs(isite)%Name)
            call AllocateBosonicField(Pimp,LocalOrbs(isite)%Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            call calc_Pimp(Pimp,curlyU,ChiCmats) !,NaNb=Test_flag_3)
            !
            !Compute convergence benchmark for the interaction
            call AllocateBosonicField(Wimp,LocalOrbs(isite)%Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
            call calc_Wimp(Wimp,curlyU,ChiCmats)
            !
            !Expand to the Lattice basis
            call imp2loc(C_EDMFT,ChiCmats,isite,LocalOrbs,ExpandImpurity,AFMselfcons,RotateHloc,name="Cimp")
            call imp2loc(P_EDMFT,Pimp,isite,LocalOrbs,ExpandImpurity,AFMselfcons,RotateHloc,name="Pimp")
            call imp2loc(W_EDMFT,Wimp,isite,LocalOrbs,ExpandImpurity,AFMselfcons,RotateHloc,name="Wimp")
            call imp2loc(curlyU_EDMFT,curlyU,isite,LocalOrbs,ExpandImpurity,AFMselfcons,RotateUloc,name="curlyU")
            !
            call DeallocateBosonicField(ChiCmats)
            call DeallocateBosonicField(curlyU)
            call DeallocateBosonicField(Pimp)
            call DeallocateBosonicField(Wimp)
            !
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
         !symmetrize and print
         if(verbose)call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_noSym_w.DAT")
         call symmetrize_GW(P_EDMFT,EqvGWndx)
         !
         !Remove the iw=0 divergency of local polarization
         if(removeCDW_P)then
            write(*,"(A)") new_line("A")//"     Divergency removal in Pimp(iw=0)."
            call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_CDW_w.DAT")
            if(RotateUloc)then
               call remove_CDW(P_EDMFT,"lat")
            else
               call remove_CDW(P_EDMFT,"imp_exp")
            endif
         endif
         !
         !Print
         call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_w.DAT")
         call dump_BosonicField(W_EDMFT,reg(PrevItFolder),"Wimp_w.DAT")
         call dump_BosonicField(C_EDMFT,reg(PrevItFolder),"Cimp_w.DAT")
         call dump_BosonicField(curlyU_EDMFT,reg(PrevItFolder),"curlyUimp_w.DAT")
         !
         call dump_MaxEnt(P_EDMFT,"mats",reg(PrevItFolder)//"Convergence/","Pimp",EqvGWndx%SetOrbs)
         call dump_MaxEnt(W_EDMFT,"mats",reg(PrevItFolder)//"Convergence/","Wimp",EqvGWndx%SetOrbs)
         call dump_MaxEnt(curlyU_EDMFT,"mats",reg(PrevItFolder)//"Convergence/","curlyUimp",EqvGWndx%SetOrbs)
         call dump_MaxEnt(C_EDMFT,"mats",reg(PrevItFolder)//"Convergence/","Cimp",EqvGWndx%SetOrbs)
         !
         call DeallocateBosonicField(P_EDMFT)
         call DeallocateBosonicField(W_EDMFT)
         call DeallocateBosonicField(curlyU_EDMFT)
         call DeallocateBosonicField(C_EDMFT)
         !
      endif
      !
      deallocate(Gimp,Simp)
      if(Dyson_Imprvd_F)deallocate(Fimp)
      do isite=1,Nsite
         deallocate(LocalOrbs(isite)%rho_Flav)
         deallocate(LocalOrbs(isite)%rho_OrbSpin)
         deallocate(LocalOrbs(isite)%Docc)
         if(ExpandImpurity.or.AFMselfcons)exit
      enddo
      !
   end subroutine collect_QMC_results


   !---------------------------------------------------------------------------!
   !PURPOSE: Print heterostructure embedding potentials
   !---------------------------------------------------------------------------!
   subroutine print_potentials(printpath)
      implicit none
      character(len=*),intent(in),optional  :: printpath
      type(FermionicField)                  :: Pot
      character(len=255)                    :: printpath_
      !
      printpath_=reg(ItFolder)
      if(present(printpath))printpath_=printpath
      !
      call AllocateFermionicField(Pot,Hetero%Norb,Nmats,Beta=Beta)
      !
      call clear_attributes(Pot)
      if(allocated(Hetero%P_L))then
         Pot%ws = Hetero%P_L
         call dump_FermionicField(Pot,reg(printpath_),"Pot_L_w",paramagnet)
      endif
      !
      call clear_attributes(Pot)
      if(allocated(Hetero%P_R))then
         Pot%ws = Hetero%P_R
         call dump_FermionicField(Pot,reg(printpath_),"Pot_R_w",paramagnet)
      endif
      !
      call DeallocateField(Pot)
      !
   end subroutine print_potentials


   !---------------------------------------------------------------------------!
   !PURPOSE: Deallocate all fields
   !---------------------------------------------------------------------------!
   subroutine DeallocateAllFields()
      implicit none
      if(Glat%status)then
         write(*,"(A)")"     Deallocating Glat"
         call DeallocateFermionicField(Glat)
      endif
      if(S_Full%status)then
         write(*,"(A)") "     Deallocating S_Full"
         call DeallocateFermionicField(S_Full)
      endif
      if(S_G0W0%status)then
         write(*,"(A)") "     Deallocating S_G0W0"
         call DeallocateFermionicField(S_G0W0)
      endif
      if(S_G0W0dc%status)then
         write(*,"(A)") "     Deallocating S_G0W0dc"
         call DeallocateFermionicField(S_G0W0dc)
      endif
      if(S_GW%status)then
         write(*,"(A)") "     Deallocating S_GW"
         call DeallocateFermionicField(S_GW)
      endif
      if(S_GW_C%status)then
         write(*,"(A)") "     Deallocating S_GW_C"
         call DeallocateFermionicField(S_GW_C)
      endif
      if(S_GW_X%status)then
         write(*,"(A)") "     Deallocating S_GW_X"
         call DeallocateFermionicField(S_GW_X)
      endif
      if(S_GWdc%status)then
         write(*,"(A)") "     Deallocating S_GWdc"
         call DeallocateFermionicField(S_GWdc)
      endif
      if(S_GW_Cdc%status)then
         write(*,"(A)") "     Deallocating S_GW_Cdc"
         call DeallocateFermionicField(S_GW_Cdc)
      endif
      if(S_GW_Xdc%status)then
         write(*,"(A)") "     Deallocating S_GW_Xdc"
         call DeallocateFermionicField(S_GW_Xdc)
      endif
      if(S_DMFT%status)then
         write(*,"(A)") "     Deallocating S_DMFT"
         call DeallocateFermionicField(S_DMFT)
      endif
      if(Wlat%status)then
         write(*,"(A)") "     Deallocating Wlat"
         call DeallocateBosonicField(Wlat)
      endif
      if(Ulat%status)then
         write(*,"(A)") "     Deallocating Ulat"
         call DeallocateBosonicField(Ulat)
      endif
      if(Plat%status)then
         write(*,"(A)") "     Deallocating Plat"
         call DeallocateBosonicField(Plat)
      endif
      if(P_EDMFT%status)then
         write(*,"(A)") "     Deallocating P_EDMFT"
         call DeallocateBosonicField(P_EDMFT)
      endif
      if(C_EDMFT%status)then
         write(*,"(A)") "     Deallocating C_EDMFT"
         call DeallocateBosonicField(C_EDMFT)
      endif
      if(Delta_correction%status)then
         write(*,"(A)") "     Deallocating Delta_correction"
         call DeallocateFermionicField(Delta_correction)
      endif
      if(curlyU_correction%status)then
         write(*,"(A)") "     Deallocating curlyU_correction"
         call DeallocateBosonicField(curlyU_correction)
      endif
   end subroutine DeallocateAllFields


   !---------------------------------------------------------------------------!
   !PURPOSE: Print the different density matrices
   !---------------------------------------------------------------------------!
   subroutine show_Densities(Iteration,enlrg)
      use utils_misc
      implicit none
      integer,intent(in)                    :: Iteration
      integer,intent(in),optional           :: enlrg
      integer                               :: Norb,Norb_imp
      integer                               :: iorb,jorb,isite,ispin
      integer                               :: wn,ws,wsi,wnmin
      integer                               :: l1,l2,l3,l4,l5,l6,l7
      integer                               :: enlrg_
      character(len=255)                    :: header1="Lattice density"
      character(len=255)                    :: header2="Impurity density"
      character(len=255)                    :: header3="Solver density"
      character(len=255)                    :: header4="Impurity magnetization"
      character(len=255)                    :: header5="Solver magnetization"
      character(len=255)                    :: header6="Lattice magnetization"
      character(len=255)                    :: header7="Density error"
      !
      Norb = Crystal%Norb
      !
      l1=len(trim(header1)//" up")
      l2=len(trim(header2)//" up")
      l3=len(trim(header3)//" up")
      l4=len(trim(header4))
      l5=len(trim(header5))
      l6=len(trim(header6))
      l7=len(trim(header7)//" up")
      !
      wnmin=max(maxval([l1,l2,l3,l4,l5,l6,l7]),Norb*6) !6 because I have 4 precision
      enlrg_=3
      if(present(enlrg))enlrg_=enlrg
      !
      wn=int(wnmin/Norb)+enlrg_
      ws=2
      !
      select case(reg(CalculationType))
         case default
            !
            call printHeader(Iteration)
            !
            write(*,*)
            if(Nspin.eq.1)write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header1)//" up",wn*Norb)
            if(Nspin.gt.1)write(*,"(2(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header1)//" up",wn*Norb),banner(trim(header1)//" dw",wn*Norb)
            do iorb=1,Norb
               write(*,"("//str(Nspin)//"("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X))") ((dreal(densityGW(iorb,jorb,ispin)),jorb=1,Norb),ispin=1,Nspin)!(dreal(densityGW(iorb,jorb,2)),jorb=1,Norb)
            enddo
            !
            write(*,*)
            if(Nspin.eq.1)write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header2)//" up",wn*Norb)
            if(Nspin.gt.1)write(*,"(2(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header2)//" up",wn*Norb),banner(trim(header2)//" dw",wn*Norb)
            do iorb=1,Norb
               write(*,"("//str(Nspin)//"("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X))") ((dreal(densityDMFT(iorb,jorb,ispin)),jorb=1,Norb),ispin=1,Nspin)!,(dreal(densityDMFT(iorb,jorb,2)),jorb=1,Norb)
            enddo
            !
            write(*,*)
            write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header7)//" up",wn*Norb)
            write(*,"("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X)") (dreal(densityGW(iorb,iorb,1)-densityDMFT(iorb,iorb,1)),iorb=1,Norb)
            !
            if((EqvGWndx%para.eq.0).and.(Nspin.ne.1))then
               write(*,*)
               write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header7)//" dw",wn*Norb)
               write(*,"("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X)") (dreal(densityGW(iorb,iorb,2)-densityDMFT(iorb,iorb,2)),iorb=1,Norb)
               write(*,*)
               write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header6),wn*Norb)
               write(*,"("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X)") (dreal(densityGW(iorb,iorb,1)-densityGW(iorb,iorb,2)),iorb=1,Norb)
               write(*,*)
               write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header4),wn*Norb)
               write(*,"("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X)") (dreal(densityDMFT(iorb,iorb,1)-densityDMFT(iorb,iorb,2)),iorb=1,Norb)
            endif
            !
            if(Nsite.gt.1)then
               do isite=1,Nsite
                  !
                  Norb_imp = LocalOrbs(isite)%Norb
                  wsi = wn*Norb - wn*Norb_imp
                  !
                  write(*,*)
                  if(Nspin.eq.1)write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header3)//" up",wn*Norb)
                  if(Nspin.gt.1)write(*,"("//str(Nspin)//"(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header3)//" up",wn*Norb),banner(trim(header3)//" dw",wn*Norb)
                  do iorb=1,Norb_imp
                     write(*,"("//str(Nspin)//"("//str(wsi)//"X,"//str(Norb_imp)//"F"//str(wn)//".4,"//str(ws)//"X))") ((LocalOrbs(isite)%rho_OrbSpin(iorb,jorb,ispin),jorb=1,Norb_imp),ispin=1,Nspin)
                  enddo
                  !
                  if(ExpandImpurity.or.AFMselfcons)exit
                  !
               enddo
            endif
            !
         case("G0W0","scGW")
            !
            call printHeader(Iteration)
            !
            write(*,*)
            if(Nspin.eq.1)write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header1)//" up",wn*Norb)
            if(Nspin.gt.1)write(*,"(2(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header1)//" up",wn*Norb),banner(trim(header1)//" dw",wn*Norb)
            do iorb=1,Norb
               write(*,"("//str(Nspin)//"("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X))") ((dreal(densityGW(iorb,jorb,1)),jorb=1,Norb),ispin=1,Nspin)!(dreal(densityGW(iorb,jorb,2)),jorb=1,Norb)
            enddo
            !
      end select
      !
      contains
         !
         function banner(txt,width) result(string)
           character(len=*)             :: txt
           integer                      :: width
           character(len=width)         :: string
           if(width.lt.len(txt))stop "banner will be truncated"
           string = txt
           string = adjustr(string)
         end function banner
         !
   end subroutine show_Densities


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate Glat, S_DMFT and P_DMFT from an old iteration
   !---------------------------------------------------------------------------!
   subroutine  interpolate_from_oldBeta()
      !
      implicit none
      !
      complex(8),allocatable                :: HartreeU(:,:,:),Nimp(:,:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- interpolate_from_oldBeta"
      !
      !
      select case(reg(CalculationType))
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0","scGW")
            !
            !Lattice Gf
            write(*,"(A)") "     Interpolating lattice Green function from Beta="//str(Beta_Match%Beta_old,2)//" to Beta="//str(Beta_Match%Beta_new,2)
            call AllocateFermionicField(Glat,Crystal%Norb,Beta_Match%Nmats_old,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta_Match%Beta_old)
            call read_FermionicField(Glat,reg(Beta_Match%Path),"Glat_w",Crystal%kpt)
            call interpolate2Beta(Glat,Beta_Match,"lat",.true.)
            call dump_FermionicField(Glat,reg(PrevItFolder),"Glat_w",.true.,Crystal%kpt,paramagnet)
            call DeallocateFermionicField(Glat)
            !
         case("DMFT+statU","DMFT+dynU")
            !
            !Impurity Self-energy
            write(*,"(A)") "     Interpolating impurity self-energy from Beta="//str(Beta_Match%Beta_old,2)//" to Beta="//str(Beta_Match%Beta_new,2)
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Beta_Match%Nmats_old,Nsite=Nsite,Beta=Beta_Match%Beta_old)
            call read_FermionicField(S_DMFT,reg(Beta_Match%Path),"Simp_w")
            call interpolate2Beta(S_DMFT,Beta_Match,"imp",ExpandImpurity)
            call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w",paramagnet)
            call DeallocateFermionicField(S_DMFT)
            !
         case("EDMFT")
            !
            !Polarization
            write(*,"(A)") "     Interpolating impurity polarization from Beta="//str(Beta_Match%Beta_old,2)//" to Beta="//str(Beta_Match%Beta_new,2)
            call AllocateBosonicField(P_EDMFT,Crystal%Norb,Beta_Match%Nmats_old,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta_Match%Beta_old)
            call read_BosonicField(P_EDMFT,reg(Beta_Match%Path),"Pimp_w.DAT")
            call interpolate2Beta(P_EDMFT,Beta_Match,"imp",.false.)
            !Remove the iw=0 divergency of local charge susceptibility
            if(removeCDW_P)then
               write(*,"(A)") new_line("A")//"     Divergency removal in Pimp(iw=0)."
               call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_CDW_w.DAT")
               if(RotateUloc)then
                  call remove_CDW(P_EDMFT,"lat")
               else
                  call remove_CDW(P_EDMFT,"imp_exp")
               endif
            endif
            call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_w.DAT")
            call DeallocateBosonicField(P_EDMFT)
            !
            !Impurity Self-energy
            write(*,"(A)") "     Interpolating impurity self-energy from Beta="//str(Beta_Match%Beta_old,2)//" to Beta="//str(Beta_Match%Beta_new,2)
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Beta_Match%Nmats_old,Nsite=Nsite,Beta=Beta_Match%Beta_old)
            call read_FermionicField(S_DMFT,reg(Beta_Match%Path),"Simp_w")
            call interpolate2Beta(S_DMFT,Beta_Match,"imp",ExpandImpurity)
            call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w",paramagnet)
            call DeallocateFermionicField(S_DMFT)
            !
            !Read&write instead of execute_command
            allocate(HartreeU(Crystal%Norb,Crystal%Norb,Nspin));HartreeU=czero
            call read_Matrix(HartreeU,reg(Beta_Match%Path)//"Hartree_UNimp",paramagnet)
            call dump_Matrix(HartreeU,reg(PrevItFolder),"Hartree_UNimp",paramagnet)
            deallocate(HartreeU)
            allocate(Nimp(Crystal%Norb,Crystal%Norb,Nspin));Nimp=czero
            call read_Matrix(Nimp,reg(Beta_Match%Path)//"Nimp",paramagnet)
            call dump_Matrix(Nimp,reg(PrevItFolder),"Nimp",paramagnet)
            deallocate(Nimp)
            !
         case("GW+EDMFT")
            !
            !Polarization
            write(*,"(A)") "     Interpolating impurity polarization from Beta="//str(Beta_Match%Beta_old,2)//" to Beta="//str(Beta_Match%Beta_new,2)
            call AllocateBosonicField(P_EDMFT,Crystal%Norb,Beta_Match%Nmats_old,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta_Match%Beta_old)
            call read_BosonicField(P_EDMFT,reg(Beta_Match%Path),"Pimp_w.DAT")
            call interpolate2Beta(P_EDMFT,Beta_Match,"imp",.false.)
            !Remove the iw=0 divergency of local charge susceptibility
            if(removeCDW_P)then
               write(*,"(A)") new_line("A")//"     Divergency removal in Pimp(iw=0)."
               call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_CDW_w.DAT")
               if(RotateUloc)then
                  call remove_CDW(P_EDMFT,"lat")
               else
                  call remove_CDW(P_EDMFT,"imp_exp")
               endif
            endif
            call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_w.DAT")
            call DeallocateBosonicField(P_EDMFT)
            !
            !Impurity Self-energy
            write(*,"(A)") "     Interpolating impurity self-energy from Beta="//str(Beta_Match%Beta_old,2)//" to Beta="//str(Beta_Match%Beta_new,2)
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Beta_Match%Nmats_old,Nsite=Nsite,Beta=Beta_Match%Beta_old)
            call read_FermionicField(S_DMFT,reg(Beta_Match%Path),"Simp_w")
            call interpolate2Beta(S_DMFT,Beta_Match,"imp",ExpandImpurity)
            call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w",paramagnet)
            call DeallocateFermionicField(S_DMFT)
            !
            !Lattice Gf
            write(*,"(A)") "     Interpolating lattice Green function from Beta="//str(Beta_Match%Beta_old,2)//" to Beta="//str(Beta_Match%Beta_new,2)
            call AllocateFermionicField(Glat,Crystal%Norb,Beta_Match%Nmats_old,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta_Match%Beta_old)
            call read_FermionicField(Glat,reg(Beta_Match%Path),"Glat_w",Crystal%kpt)
            call interpolate2Beta(Glat,Beta_Match,"lat",.true.)
            call dump_FermionicField(Glat,reg(PrevItFolder),"Glat_w",.true.,Crystal%kpt,paramagnet)
            call DeallocateFermionicField(Glat)
            !
            !Read&write instead of execute_command
            allocate(HartreeU(Crystal%Norb,Crystal%Norb,Nspin));HartreeU=czero
            call read_Matrix(HartreeU,reg(Beta_Match%Path)//"Hartree_UNimp",paramagnet)
            call dump_Matrix(HartreeU,reg(PrevItFolder),"Hartree_UNimp",paramagnet)
            deallocate(HartreeU)
            allocate(Nimp(Crystal%Norb,Crystal%Norb,Nspin));Nimp=czero
            call read_Matrix(Nimp,reg(Beta_Match%Path)//"Nimp",paramagnet)
            call dump_Matrix(Nimp,reg(PrevItFolder),"Nimp",paramagnet)
            deallocate(Nimp)
            !
      end select
      !
   end subroutine interpolate_from_oldBeta


   !---------------------------------------------------------------------------!
   !PURPOSE: Extract the occupation from a subset of orbtials
   !---------------------------------------------------------------------------!
   function get_Tier_occupation(rho,LocalOrbs) result(occupation)
      implicit none
      complex(8),intent(in)                 :: rho(:,:,:)
      type(LocalOrbitals),allocatable,intent(in):: LocalOrbs(:)
      real(8)                               :: occupation
      integer                               :: Norb,Nspin,Nsite
      integer                               :: iorb,ispin,isite,orbndx
      if(.not.allocated(LocalOrbs)) stop "get_Tier_occupation: LocalOrbs not properly initialized."
      Norb=size(rho,dim=1)
      Nspin=size(rho,dim=3)
      Nsite=size(LocalOrbs)
      if(size(rho,dim=2).ne.Norb)stop"get_Tier_occupation: density matrix not square."
      occupation=0d0
      do isite=1,Nsite
         if(LocalOrbs(isite)%Norb.gt.Norb)stop"get_Tier_occupation: orbital list bigger than density matrix."
         do iorb=1,LocalOrbs(isite)%Norb
            orbndx = LocalOrbs(isite)%Orbs(iorb)
            if(orbndx.gt.Norb)stop"get_Tier_occupation: orbital outside the density matrix space."
            do ispin=1,Nspin
               occupation = occupation + real(rho(orbndx,orbndx,ispin))
            enddo
         enddo
      enddo
   end function get_Tier_occupation


end module utils_main
