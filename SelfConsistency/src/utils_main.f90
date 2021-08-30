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
   !moved to input_vars
   !complex(8),allocatable                   :: OlocSite(:,:,:)
   !complex(8),allocatable                   :: OlocRot(:,:,:)
   !complex(8),allocatable                   :: OlocRotDag(:,:,:)
   !real(8),allocatable                      :: OlocEig(:,:)
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
   type(BosonicField)                       :: Wlat
   type(BosonicField)                       :: W_EDMFT
   !
   type(BosonicField)                       :: Ulat
   type(BosonicField)                       :: Plat
   type(BosonicField)                       :: P_EDMFT
   type(BosonicField)                       :: C_EDMFT
   type(BosonicField)                       :: curlyU_EDMFT
   !
   type(FermionicField)                     :: D_correction
   type(BosonicField)                       :: curlyU_correction
   !
   real(8)                                  :: Ek,Ep
   !
   real(8)                                  :: density2set
   complex(8),allocatable                   :: densityLDA(:,:,:)
   complex(8),allocatable                   :: densityGW(:,:,:)
   complex(8),allocatable                   :: densityDMFT(:,:,:)
   real(8),allocatable                      :: densityQMC(:,:,:,:)
   !
   complex(8),allocatable                   :: Vxc(:,:,:,:)
   complex(8),allocatable                   :: VH(:,:)
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
   real(8)                                  :: HartreeFact=1d0
   logical                                  :: update_curlyG=.true.
   integer                                  :: SigmaMaxMom=3
   !
   logical                                  :: sym_constrained=.false.
   logical                                  :: MultiTier=.false.
   logical                                  :: print_path=.false.
   !
   character(len=255)                       :: ItFolder,PrevItFolder
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
      !
      !Few general checks
      if(ExpandImpurity.and.AFMselfcons) stop "AFM self-consistency and expansion to real space not yet implemented."
      if(ExpandImpurity.and.(Nsite.eq.1)) stop "Cannot expand a single site."
      if(AFMselfcons.and.(Nsite.ne.2)) stop "AFM self-consistency is implemented only for lattices with 2 sites."
      if(RotateUloc.and.(.not.RotateHloc)) stop "Rotate the Bosonic impurity problem without rotating the Ferminic one is not allowed."
      !
      if(Hmodel)HartreeFact=0d0
      !
      causal_D = causal_D .and. (FirstIteration.ne.0)
      causal_U = causal_U .and. (FirstIteration.ne.0)
      !
      calc_Sguess = calc_Sguess .and. (FirstIteration.eq.0) .and. solve_DMFT
      calc_Pguess = calc_Pguess .and. (FirstIteration.eq.0) .and. solve_DMFT
      !
      sym_constrained = ExpandImpurity .and. (reg(VH_type).ne."Ustatic")
      sym_constrained = ExpandImpurity .and. (reg(VH_type).ne."Ubare")
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
               call createDir(reg(Itpath)//"Solver_"//reg(SiteName(isite))//"/fits",verb=verbose)
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
               if(.not.ZeroItexist)stop "G0W0 0th iteration not found, exiting."
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
               if((reg(CalculationType).eq."GW+EDMFT").and.(.not.ZeroItexist)) stop "G0W0 0th iteration not found, exiting."
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
                  call createDir(reg(Itpath)//"/Solver_"//reg(SiteName(isite))//"/fits",verb=verbose)
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
      !Just to be sure
      Ustart = Ustart .and. (ItStart.eq.0)
      calc_Pguess = calc_Pguess .and. (.not.Ustart)
      !
      if((reg(path_funct).ne."G") .and. (reg(path_funct).ne."S") .and. (reg(path_funct).ne."GS") .and. (reg(path_funct).ne."None"))then
         write(*,"(A)")"     Invalid content of PATH_FUNCT. The calculation on the path will be avoided."
         print_path=.false.
      elseif(reg(path_funct).eq."None")then
         print_path=.false.
      endif
      !
      write(*,"(A)") "     HartreeFact: "//str(HartreeFact)
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
      integer                               :: isite,iadd,iset,iorb
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
      Lttc%Nkpt3 = Nkpt3
      !
      !
      select case(reg(CalculationType))
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0","scGW")
            !
            if(Hmodel)then
               !
               call build_Hk(LatticeVec,Norb_model,hopping,Nkpt3,alphaHk,Hetero,&
                             Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc,        &
                             iq_gamma=Lttc%iq_gamma,pathOUTPUT=reg(pathINPUT))
               Lttc%Norb = size(Lttc%Hk,dim=1)
               Lttc%Nkpt = size(Lttc%Hk,dim=3)
               Lttc%Nkpt_irred = Lttc%Nkpt
               !
            else
               !
               call read_Hk(pathINPUT,alphaHk,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc,iq_gamma_Hk)
               Lttc%Norb = size(Lttc%Hk,dim=1)
               Lttc%Nkpt = size(Lttc%Hk,dim=3)
               !
               allocate(Lttc%kptPos(Lttc%Nkpt));Lttc%kptPos=0
               call read_xeps(pathINPUT,Lttc%kpt,Nkpt3,UseXepsKorder, &
               Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,iq_gamma_XEPS,paramagneticSPEX)
               XEPSisread=.true.
               if(iq_gamma_Hk.eq.iq_gamma_XEPS)then
                  Lttc%iq_gamma = iq_gamma_Hk
               else
                  stop "Index of the Gamma point is not the same in H(k) and XEPS."
               endif
               !
            endif
            !
            allocate(Lttc%kptsum(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptsum=0
            allocate(Lttc%kptdif(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptdif=0
            call fill_ksumkdiff(Lttc%kpt,Lttc%kptsum,Lttc%kptdif,Nkpt3)
            !
            allocate(Lttc%small_ik(12,2));Lttc%small_ik=0
            call fill_smallk(Lttc%kpt,Lttc%small_ik)
            !
            Lttc%status=.true.
            !
         case("DMFT+statU","DMFT+dynU")
            !
            if(Hmodel)then
               !
               call build_Hk(LatticeVec,Norb_model,hopping,Nkpt3,alphaHk,Hetero,&
                             Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc,        &
                             pathOUTPUT=reg(pathINPUT))
               Lttc%Norb = size(Lttc%Hk,dim=1)
               Lttc%Nkpt = size(Lttc%Hk,dim=3)
               Lttc%Nkpt_irred = Lttc%Nkpt
               !
            else
               !
               call read_Hk(pathINPUT,alphaHk,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc)
               Lttc%Norb = size(Lttc%Hk,dim=1)
               Lttc%Nkpt = size(Lttc%Hk,dim=3)
               !
            endif
            !
            Lttc%status=.true.
            !
         case("EDMFT","GW+EDMFT")
            !
            if(Hmodel)then
               !
               call build_Hk(LatticeVec,Norb_model,hopping,Nkpt3,alphaHk,Hetero,&
                             Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc,        &
                             iq_gamma=Lttc%iq_gamma,pathOUTPUT=reg(pathINPUT))
               Lttc%Norb = size(Lttc%Hk,dim=1)
               Lttc%Nkpt = size(Lttc%Hk,dim=3)
               Lttc%Nkpt_irred = Lttc%Nkpt
               !
            else
               !
               call read_Hk(pathINPUT,alphaHk,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc,iq_gamma_Hk)
               Lttc%Norb = size(Lttc%Hk,dim=1)
               Lttc%Nkpt = size(Lttc%Hk,dim=3)
               !
               allocate(Lttc%kptPos(Lttc%Nkpt));Lttc%kptPos=0
               call read_xeps(pathINPUT,Lttc%kpt,Nkpt3,UseXepsKorder, &
               Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,iq_gamma_XEPS,paramagneticSPEX)
               XEPSisread=.true.
               if(iq_gamma_Hk.eq.iq_gamma_XEPS)then
                  Lttc%iq_gamma = iq_gamma_Hk
               else
                  stop "Index of the Gamma point is not the same in H(k) and XEPS."
               endif
               !
            endif
            !
            allocate(Lttc%kptsum(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptsum=0
            allocate(Lttc%kptdif(Lttc%Nkpt,Lttc%Nkpt));Lttc%kptdif=0
            call fill_ksumkdiff(Lttc%kpt,Lttc%kptsum,Lttc%kptdif,Nkpt3)
            !
            allocate(Lttc%small_ik(12,2));Lttc%small_ik=0
            call fill_smallk(Lttc%kpt,Lttc%small_ik)
            !
            Lttc%status=.true.
            !
      end select
      if(Lttc%Nkpt.ne.(Nkpt3(1)*Nkpt3(2)*Nkpt3(3)))stop "Total number of K-points does not match with number of K-points per dimension."
      !
      !
      wrealMax = maxval(abs(Lttc%Ek)) + 0.1*maxval(abs(Lttc%Ek))
      !
      !
      !Store the local Hamiltonian
      call dump_Matrix(Lttc%Hloc,reg(pathINPUT),"Hloc.DAT")
      !
      !
      !print some info
      write(*,"(A)") new_line("A")//new_line("A")//"---- site and orbital structure"
      write(*,"(A,1I3)") "     Number of inequivalent sites: ",Nsite
      write(*,"(A,1I3)") "     Number of solved impurities : ",Solver%Nimp
      do isite=1,Nsite
         write(*,"(2(A,I3),A,10I3)")"     site: ",isite,", orbital space: ",SiteNorb(isite),", orbital indexes: ",SiteOrbs(isite,1:SiteNorb(isite))
      enddo
      if(sum(SiteNorb).ne.Lttc%Norb)then
         !
         MultiTier = .true.
         if(Hetero%status) stop "MultiTier construction and Heterostructured setup are not allowed together."
         if(RotateHloc.or.RotateUloc) stop "MultiTier construction and Rotations of the local space are not allowed together."
         if(reg(DC_type).eq."GlocWloc")then
            DC_type="Sloc"
            write(*,"(A,1I3)") "     DC_TYPE updated from GlocWloc to "//reg(DC_type)
         endif
         !
      endif
      !
      !
      !Store the local rotation of each site and add to the input list the Impurity equivalent orbitals
      if(RotateHloc)then
         !
         write(*,"(A)") new_line("A")//new_line("A")//"---- Rotations of the local LDA Hamiltonian (used)"
         call build_rotations("Hloc",OlocSite,OlocEig,OlocRot,OlocRotDag,LatticeOp=Lttc%Hloc)
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
      if(EqvGWndx%Nset.eq.0)then
         !
         EqvGWndx%O=.false.
         EqvGWndx%Ntotset=EqvGWndx%Nset
         EqvGWndx%Gfoffdiag=.false.
         !
      elseif((EqvGWndx%Nset.ne.0).and.(.not.ExpandImpurity))then
         !
         EqvGWndx%O = sym_mode.le.2
         EqvGWndx%Ntotset=EqvGWndx%Nset
         EqvGWndx%Gfoffdiag = sym_mode.le.2
         !
         if(sum(EqvGWndx%SetNorb).lt.Lttc%Norb)then !sum(SiteNorb)
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
      write(*,"(A,L)") "     Syimmetrizing off-diagonal: ",EqvGWndx%Gfoffdiag
      if(sym_mode.eq.1)write(*,"(A)") "     Only lattice quantities will be symmetrized."
      if(sym_mode.eq.2)write(*,"(A)") "     Both lattice and impurity quantities will be symmetrized."
      if(sym_mode.eq.3)write(*,"(A)") "     Only impurity quantities will be symmetrized."
      !
      !
      !Dump some LDA results
      if(ItStart.eq.0)then
         !
         call calc_Glda(0d0,Beta,Lttc)
         !
         allocate(Egrid(Nreal));Egrid=0d0
         Egrid = linspace(-wrealMax,+wrealMax,Nreal)
         call tetrahedron_integration(reg(pathINPUT),Lttc%Ek,Lttc%Nkpt3,Lttc%kpt,Egrid,fact_intp=2,pathOUTPUT=reg(pathINPUT))
         deallocate(Egrid)
         !
         if(reg(structure).ne."None")call interpolateHk2Path(Lttc,reg(structure),Nkpt_path,reg(pathINPUT),doplane=.true.)
         !
      endif
      !
      !
   end subroutine initialize_Lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute or read the site-dependent rotations with respect of the
   !         LatticeOp. Fixed rotations are stored in pathINPUT, iteration
   !         dependent rotations are stored in the iteration folder.
   !---------------------------------------------------------------------------!
   subroutine build_rotations(Opname,SiteOp,SiteEig,SiteRot,SiteRotdag,LatticeOp,read)
      !
      implicit none
      !
      character(len=*),intent(in)           :: Opname
      complex(8),allocatable,intent(out)    :: SiteOp(:,:,:)
      real(8),allocatable,intent(out)       :: SiteEig(:,:)
      complex(8),allocatable,intent(out)    :: SiteRot(:,:,:)
      complex(8),allocatable,intent(out)    :: SiteRotdag(:,:,:)
      complex(8),intent(in),optional        :: LatticeOp(:,:)
      logical,intent(in),optional           :: read
      !
      character(len=255)                    :: Folder
      integer                               :: isite,Norb
      integer,allocatable                   :: Orbs(:)
      complex(8),allocatable                :: Oloc(:,:),Rot(:,:)
      real(8),allocatable                   :: Eig(:),EigR(:,:)
      logical                               :: storeIt,read_,build_
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
      if(allocated(SiteOp))deallocate(SiteOp)
      if(allocated(SiteEig))deallocate(SiteEig)
      if(allocated(SiteRot))deallocate(SiteRot)
      if(allocated(SiteRotdag))deallocate(SiteRotdag)
      !
      if(build_)then
         !
         if(ExpandImpurity)then !Only one set of orbital provided - all matrices with same "Norb" dimension
            !
            allocate(SiteOp(SiteNorb(1),SiteNorb(1),Nsite));SiteOp=czero
            allocate(SiteEig(SiteNorb(1),Nsite));SiteEig=0d0
            allocate(SiteRot(SiteNorb(1),SiteNorb(1),Nsite));SiteRot=czero
            allocate(SiteRotdag(SiteNorb(1),SiteNorb(1),Nsite));SiteRotdag=czero
            !
            do isite=1,Nsite
               !
               if(storeIt)Folder=reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/"
               !
               Norb = SiteNorb(1)
               if(Norb.ne.SiteNorb(isite)) stop "The orbital space is not the same among the different impurities."
               !
               allocate(Oloc(Norb,Norb));Oloc=czero
               allocate(Eig(Norb));Eig=0d0
               allocate(Rot(Norb,Norb));Rot=czero
               allocate(Orbs(Norb));Orbs=0
               !
               Orbs=SiteOrbs(isite,:)
               !
               !Extract then local Hamiltonian for each site
               call loc2imp(Oloc,LatticeOp,Orbs)
               Rot = dreal(Oloc)
               !
               !Rotate
               call eigh(Rot,Eig)
               !
               !Save
               call dump_Matrix(Oloc,reg(Folder),reg(Opname)//"Site_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               call dump_Matrix(diag(Eig),reg(Folder),reg(Opname)//"Eig_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               call dump_Matrix(Rot,reg(Folder),reg(Opname)//"Rot_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               !
               !SO(3)check
               write(*,"(A,"//str(Norb)//"I3)") "     Orbs: ",Orbs
               write(*,"(A,2F)") "     det(Rot) of "//reg(SiteName(1))//"_"//str(isite)//" :",det(Rot)
               !
               SiteOp(1:Norb,1:Norb,isite) = Oloc
               SiteEig(1:Norb,isite) = Eig
               SiteRot(1:Norb,1:Norb,isite) = Rot
               SiteRotdag(1:Norb,1:Norb,isite) = conjg(transpose(Rot))
               !
               deallocate(Oloc,Eig,Rot,Orbs)
               !
            enddo
            !
         else !Potentially different set of orbitals provided - all matrices with different dimension (SiteNorb(isite))
            !
            allocate(SiteOp(maxval(SiteNorb),maxval(SiteNorb),Nsite));SiteOp=czero
            allocate(SiteEig(maxval(SiteNorb),Nsite));SiteEig=0d0
            allocate(SiteRot(maxval(SiteNorb),maxval(SiteNorb),Nsite));SiteRot=czero
            allocate(SiteRotdag(maxval(SiteNorb),maxval(SiteNorb),Nsite));SiteRotdag=czero
            !
            do isite=1,Nsite
               !
               if(storeIt)Folder=reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/"
               !
               Norb = SiteNorb(isite)
               !
               allocate(Oloc(Norb,Norb));Oloc=czero
               allocate(Eig(Norb));Eig=0d0
               allocate(Rot(Norb,Norb));Rot=czero
               allocate(Orbs(Norb));Orbs=0
               !
               Orbs=SiteOrbs(isite,:)
               !
               !Extract then local Hamiltonian for each site
               call loc2imp(Oloc,LatticeOp,Orbs)
               Rot = dreal(Oloc)
               !
               !Rotate
               call eigh(Rot,Eig)
               !
               !Save
               call dump_Matrix(Oloc,reg(Folder),reg(Opname)//"Site_"//reg(SiteName(isite))//".DAT")
               call dump_Matrix(diag(Eig),reg(Folder),reg(Opname)//"Eig_"//reg(SiteName(isite))//".DAT")
               call dump_Matrix(Rot,reg(Folder),reg(Opname)//"Rot_"//reg(SiteName(isite))//".DAT")
               !
               !SO(3)check
               write(*,"(A,"//str(Norb)//"I3)") "     Orbs: ",Orbs
               write(*,"(A,F)") "     det(Rot) of "//reg(SiteName(isite))//" :",det(SiteRot(:,:,isite))
               !
               SiteOp(1:Norb,1:Norb,isite) = Oloc
               SiteEig(1:Norb,isite) = Eig
               SiteRot(1:Norb,1:Norb,isite) = Rot
               SiteRotdag(1:Norb,1:Norb,isite) = conjg(transpose(Rot))
               !
               deallocate(Oloc,Eig,Rot,Orbs)
               !
            enddo
            !
         endif
         !
      elseif(read_)then
         !
         if(ExpandImpurity)then !Only one set of orbital provided - all matrices with same "Norb" dimension
            !
            allocate(SiteOp(SiteNorb(1),SiteNorb(1),Nsite));SiteOp=czero
            allocate(SiteEig(SiteNorb(1),Nsite));SiteEig=0d0
            allocate(SiteRot(SiteNorb(1),SiteNorb(1),Nsite));SiteRot=czero
            allocate(SiteRotdag(SiteNorb(1),SiteNorb(1),Nsite));SiteRotdag=czero
            !
            do isite=1,Nsite
               !
               if(storeIt)Folder=reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/"
               !
               Norb = SiteNorb(1)
               if(Norb.ne.SiteNorb(isite)) stop "The orbital space is not the same among the different impurities."
               !
               allocate(Oloc(Norb,Norb));Oloc=czero
               allocate(EigR(Norb,Norb));EigR=0d0
               allocate(Rot(Norb,Norb));Rot=czero
               allocate(Orbs(Norb));Orbs=0
               !
               Orbs=SiteOrbs(isite,:)
               !
               !Read
               call read_Matrix(Oloc,reg(Folder)//reg(Opname)//"Site_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               call read_Matrix(EigR,reg(Folder)//reg(Opname)//"Eig_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               call read_Matrix(Rot,reg(Folder)//reg(Opname)//"Rot_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               !
               !SO(3)check
               write(*,"(A,"//str(Norb)//"I3)") "     Orbs: ",Orbs
               write(*,"(A,2F)") "     det(Rot) of "//reg(SiteName(1))//"_"//str(isite)//" :",det(Rot)
               !
               SiteOp(1:Norb,1:Norb,isite) = Oloc
               SiteEig(1:Norb,isite) = diagonal(EigR)
               SiteRot(1:Norb,1:Norb,isite) = Rot
               SiteRotdag(1:Norb,1:Norb,isite) = conjg(transpose(Rot))
               !
               deallocate(Oloc,EigR,Rot,Orbs)
               !
            enddo
            !
         else !Potentially different set of orbitals provided - all matrices with different dimension (SiteNorb(isite))
            !
            allocate(SiteOp(maxval(SiteNorb),maxval(SiteNorb),Nsite));SiteOp=czero
            allocate(SiteEig(maxval(SiteNorb),Nsite));SiteEig=0d0
            allocate(SiteRot(maxval(SiteNorb),maxval(SiteNorb),Nsite));SiteRot=czero
            allocate(SiteRotdag(maxval(SiteNorb),maxval(SiteNorb),Nsite));SiteRotdag=czero
            !
            do isite=1,Nsite
               !
               if(storeIt)Folder=reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/"
               !
               Norb = SiteNorb(isite)
               !
               allocate(Oloc(Norb,Norb));Oloc=czero
               allocate(EigR(Norb,Norb));EigR=0d0
               allocate(Rot(Norb,Norb));Rot=czero
               allocate(Orbs(Norb));Orbs=0
               !
               Orbs=SiteOrbs(isite,:)
               !
               !Read
               call read_Matrix(Oloc,reg(Folder)//reg(Opname)//"Site_"//reg(SiteName(isite))//".DAT")
               call read_Matrix(EigR,reg(Folder)//reg(Opname)//"Eig_"//reg(SiteName(isite))//".DAT")
               call read_Matrix(Rot,reg(Folder)//reg(Opname)//"Rot_"//reg(SiteName(isite))//".DAT")
               !
               !SO(3)check
               write(*,"(A,"//str(Norb)//"I3)") "     Orbs: ",Orbs
               write(*,"(A,F)") "     det(Rot) of "//reg(SiteName(isite))//" :",det(SiteRot(:,:,isite))
               !
               SiteOp(1:Norb,1:Norb,isite) = Oloc
               SiteEig(1:Norb,isite) = diagonal(EigR)
               SiteRot(1:Norb,1:Norb,isite) = Rot
               SiteRotdag(1:Norb,1:Norb,isite) = conjg(transpose(Rot))
               !
               deallocate(Oloc,EigR,Rot,Orbs)
               !
            enddo
            !
         endif
         !
      endif
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
            call get_pattern(EqvImpndx(isite)%SetOrbs,OlocEig(:,isite),1e4*eps)
            !
            if(allocated(EqvImpndx(isite)%SetOrbs))then
               !
               EqvImpndx(isite)%Nset = size(EqvImpndx(isite)%SetOrbs,dim=1)
               allocate(EqvImpndx(isite)%SetNorb(EqvImpndx(isite)%Nset))
               do iset=1,EqvImpndx(isite)%Nset
                  EqvImpndx(isite)%SetNorb(iset) = size( pack( EqvImpndx(isite)%SetOrbs(iset,:), EqvImpndx(isite)%SetOrbs(iset,:).gt.0 ) )
               enddo
               !
            else
               !
               EqvImpndx(isite)%Nset = 0
               !
            endif
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
            !
            !loop over lattice sets
            do iset=1,EqvGWndx%Nset
               !
               !scroll elements in the set
               contained=.true.
               do set_ndx=1,EqvGWndx%SetNorb(iset)
                  set_orb = EqvGWndx%SetOrbs(iset,set_ndx)
                  contained = contained.and.any( SiteOrbs(isite,1:SiteNorb(isite)).eq. set_orb )
               enddo
               !
               !the site contains all the indexes within the set
               if(contained.and.(EqvGWndx%SetNorb(iset).gt.1))then
                  EqvImpndx(isite)%Nset = EqvImpndx(isite)%Nset+1
                  EqvImpndx(isite)%SetNorb(EqvImpndx(isite)%Nset) = EqvGWndx%SetNorb(iset)
                  EqvImpndx(isite)%SetOrbs(EqvImpndx(isite)%Nset,:) = EqvGWndx%SetOrbs(iset,:) - (SiteOrbs(isite,1)-1)
               endif
               !
            enddo
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
                  file = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Eqv.DAT"
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
               call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.false.,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal)
               elseif(N_Vnn.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph,LocalOnly=.false.)
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
               endif
            endif
            !
            !Fully screened interaction
            call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            !Polarization
            call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            !
            !Lattice Gf
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
               call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.true.,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal,LocalOnly=.true.)
               elseif(N_Vnn.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
                  call inquireFile(reg(pathINPUTtr)//"Uloc_mats.DAT",filexists,hardstop=.false.,verb=verbose)
                  if(filexists)then
                     call read_BosonicField(Ulat,reg(pathINPUTtr),"Uloc_mats.DAT")
                  else
                     call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph)
                  endif
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
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
               call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.false.,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal)
               elseif(N_Vnn.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph,LocalOnly=.false.)
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
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
            !Polarization
            call AllocateBosonicField(P_EDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
            if(ItStart.ne.0)call read_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_w.DAT")
            !
            !Impurity Self-energy
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w")
               call read_Matrix(S_DMFT%N_s,reg(PrevItFolder)//"HartreeU",paramagnet)
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
               call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.false.,doAC=U_AC,pathOUTPUT=reg(pathINPUTtr))
               call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
            elseif(Umodel)then
               if(Nphonons.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with long-range couplings is frequency-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,1,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,Vnn,Crystal)
               elseif(N_Vnn.eq.0)then
                  write(*,"(A)")"     Warning: the model interaction built with phononic modes is K-independent."
                  call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph,LocalOnly=.false.)
                  call dump_MaxEnt(Ulat,"mats",reg(ItFolder)//"Convergence/","Ulat",EqvGWndx%SetOrbs)
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
            !Impurity Self-energy and Hartree contribution
            call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            if(ItStart.ne.0)then
               call read_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w")
               call read_Matrix(S_DMFT%N_s,reg(PrevItFolder)//"HartreeU",paramagnet)
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
      calc_W = calc_Wedmft .or. calc_Wfull
      !
      !
      !Allocate and initialize different density matrices
      allocate(densityLDA(Crystal%Norb,Crystal%Norb,Nspin));densityLDA=czero
      allocate(densityGW(Crystal%Norb,Crystal%Norb,Nspin));densityGW=czero
      allocate(densityDMFT(Crystal%Norb,Crystal%Norb,Nspin));densityDMFT=czero
      allocate(densityQMC(maxval(SiteNorb),maxval(SiteNorb),Nspin,Nsite));densityQMC=0d0
      if(ItStart.eq.0)then
         !
         densityLDA = Glat%N_s
         densityGW=czero
         densityDMFT=czero
         densityQMC=0d0
         !
         call dump_Matrix(densityLDA,reg(pathINPUT),"Nlda",paramagnet)
         !
      else
         !
         call read_Matrix(densityLDA,reg(pathINPUT)//"Nlda",paramagnet)
         call read_Matrix(densityDMFT,reg(PrevItFolder)//"Nimp",paramagnet)
         densityGW=Glat%N_s
         !
         do isite=1,Nsite
            file = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/Nqmc.DAT"
            call inquireFile(reg(file),filexists,hardstop=.false.,verb=verbose)
            if(filexists)then
               unit = free_unit()
               open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
               read(unit,*) muQMC
               do ib1=1,Nspin*SiteNorb(isite)
                  iorb = (ib1+mod(ib1,2))/2
                  ispin = abs(mod(ib1,2)-2)
                  read(unit,*) idum,densityQMC(iorb,iorb,ispin,isite)
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
      if((reg(CalculationType).eq."GW+EDMFT").and.(.not.S_G0W0%status))then
         call AllocateFermionicField(S_G0W0,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_Sigma_spex(SpexVersion,S_G0W0,Crystal,verbose,recompute=RecomputeG0W0,pathOUTPUT=reg(pathINPUTtr))
      endif
      !
      !
      allocate(Uinst_0th(Crystal%Norb**2,Crystal%Norb**2));Uinst_0th=0d0
      !
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
                     ib1 = iorb + S_DMFT%Norb*(jorb-1)
                     ib2 = korb + S_DMFT%Norb*(lorb-1)
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
   !PURPOSE: Check that the double counting between G0W0 and scGW at the
   !         0th iteration yeld a causal local self-energy. Usually that's the case,
   !         but, given that SPEX is a zero T calculation, some very small numeircal
   !         errors might be present at low freq.
   !         This subroutine, if needed, correct the SPEX self-energy by a
   !         small rescaling factor computed just from the first matsubara
   !         frequency. If a non-causal difference between the local SPEX G0W0
   !         and local scGW is still present it will be removed in join_SigmaFull
   !         This information is present already at the 0th iteration but the
   !         input value of GoWoDC_loc cannot be changed between iterations
   !         therefore the check in join_SigmaFull is required each time the
   !         total self-energy is computed.
   !---------------------------------------------------------------------------!
   subroutine check_S_G0W0()
      !
      implicit none
      integer                               :: iorb,ispin
      real(8)                               :: y1,y2,x1,x2,m,q
      real(8),allocatable                   :: wmats_orig(:)
      type(OldBeta)                         :: Beta_Match
      logical                               :: shift
      !
      if(.not.S_G0W0dc%status) stop "check_S_G0W0: S_G0W0dc not properly initialized."
      if(.not.S_G0W0%status) stop "check_S_G0W0: S_G0W0 not properly initialized."
      !
      allocate(wmats_orig(S_G0W0%Npoints));wmats_orig=0d0
      wmats_orig = FermionicFreqMesh(S_G0W0%Beta,S_G0W0%Npoints)
      !
      !check for any component that does not extrapolate to iw=0 in the imagnary part
      shift=.false.
      checkloop: do ispin=1,Nspin
         do iorb=1,S_G0W0%Norb
            !
            y1 = dimag(S_G0W0%ws(iorb,iorb,1,ispin)) ; x1 = wmats_orig(1)
            y2 = dimag(S_G0W0%ws(iorb,iorb,2,ispin)) ; x2 = wmats_orig(2)
            m = (y2 - y1)/(x2 - x1)
            q = y1 - m*x1
            !
            !shift the S_G0W0 axis in order to have the imaginary part extrapolating to zero
            if(dimag(S_G0W0%ws(iorb,iorb,1,ispin)).ge.dimag(S_G0W0dc%ws(iorb,iorb,1,ispin)))then
               shift=.true.
               wmats_orig = wmats_orig - q/abs(m)
               exit checkloop
            endif
            !
         enddo
      enddo checkloop
      !
      if(shift)then
         !
         write(*,"(A,F)")"     Correcting G0W0 self-energy. Axis shift: ",q/m
         Beta_Match%Nmats_old = S_G0W0%Npoints
         Beta_Match%Nmats_new = S_G0W0%Npoints
         Beta_Match%Beta_old = S_G0W0%Beta
         Beta_Match%Beta_new = S_G0W0%Beta
         call interpolate2Beta(S_G0W0,Beta_Match,"lat",.true.,wmats_in=wmats_orig)
         deallocate(wmats_orig)
         !
         call dump_FermionicField(S_G0W0,reg(pathINPUTtr),"SGoWo_w",.true.,Crystal%kpt,.true.)
         !
      endif
      !
   end subroutine check_S_G0W0


   !---------------------------------------------------------------------------!
   !PURPOSE: Join the all the component of the self-energy.
   !---------------------------------------------------------------------------!
   subroutine join_SigmaFull(Iteration)
      !
      implicit none
      integer,intent(in)                    :: Iteration
      integer                               :: iorb,ik,iw,ispin
      integer                               :: isite,Norb
      integer,allocatable                   :: Orbs(:)
      complex(8),allocatable                :: Vxc_loc(:,:,:)
      real(8)                               :: ImS_1,ImS_2
      logical                               :: causal_G0W0_loc
      type(FermionicField)                  :: S_G0W0_imp
      type(FermionicField)                  :: S_G0W0_DMFT
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
            if(.not.allocated(VH)) stop "join_SigmaFull: VH not allocated."
            if(.not.allocated(Vxc)) stop "join_SigmaFull: Vxc not allocated."
            allocate(Vxc_loc(S_Full%Norb,S_Full%Norb,Nspin));Vxc_loc=czero
            do ik=1,S_Full%Nkpt
               Vxc_loc = Vxc_loc + Vxc(:,:,ik,:)/S_Full%Nkpt
            enddo
            !
            if(.not.S_G0W0%status) stop "join_SigmaFull: S_G0W0 not properly initialized."
            !
            if(Iteration.eq.0)then
               !
               !Keep only single-shot GW
               do ispin=1,Nspin
                  !$OMP PARALLEL DEFAULT(NONE),&
                  !$OMP SHARED(ispin,S_Full,S_G0W0,VH,Vxc),&
                  !$OMP PRIVATE(ik,iw)
                  !$OMP DO
                  do ik=1,S_Full%Nkpt
                     do iw=1,S_Full%Npoints
                        S_Full%wks(:,:,iw,ik,ispin) = + S_G0W0%wks(:,:,iw,ik,ispin) - Vxc(:,:,ik,ispin) + VH(:,:)
                     enddo
                  enddo
                  !$OMP END DO
                  !$OMP END PARALLEL
                  if(paramagnet)then
                     S_Full%wks(:,:,:,:,Nspin) = S_Full%wks(:,:,:,:,1)
                     exit
                  endif
               enddo
               !
            elseif(Iteration.gt.0)then
               !
               if(.not.S_G0W0dc%status) stop "join_SigmaFull: S_G0W0dc not properly initialized."
               if(.not.S_GW%status) stop "join_SigmaFull: S_GW not properly initialized."
               !
               !Remove Dc between G0W0 and scGW self-energies
               !S_G0W0dc contains only diagonal elements in the LDA basis as the S_G0W0 computed from SPEX
               S_G0W0%wks = S_G0W0%wks - S_G0W0dc%wks
               call FermionicKsum(S_G0W0)
               !
               !Compute the G0W0 contribution to the local self-energy
               call AllocateFermionicField(S_G0W0_DMFT,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta,mu=Glat%mu)
               do isite=1,Nsite
                  !
                  Norb = SiteNorb(isite)
                  allocate(Orbs(Norb))
                  Orbs = SiteOrbs(isite,1:Norb)
                  !
                  !Extract the local G0W0 self-energy for each site
                  call AllocateFermionicField(S_G0W0_imp,Norb,Nmats,Beta=Beta)
                  call loc2imp(S_G0W0_imp,S_G0W0,Orbs)
                  !
                  !Put it into an object that contains only the site indexes
                  call imp2loc(S_G0W0_DMFT,S_G0W0_imp,isite,Orbs,.false.,.false.)
                  !
                  call DeallocateField(S_G0W0_imp)
                  deallocate(Orbs)
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
               !Enclose in the EDMFT *ALL* the local contributions to the self-energy
               !From the S_G0W0^{SPEX}_{ij} + S_G0W0^{SPEX}_{i} - S_G0W0^{DC}_{ij} - S_G0W0^{DC}_{i}
               !we remove [ S_G0W0^{SPEX}_{i} - S_G0W0^{DC}_{i} ]
               !here, if Vxc is inside S_G0W0, also the local contribution from Vxc is removed
               if((.not.GoWoDC_loc).or.(.not.causal_G0W0_loc))then
                  write(*,"(A)")"     Local G0W0 self-energy removed."
                  do ik=1,S_G0W0%Nkpt
                     S_G0W0%wks(:,:,:,ik,:) = S_G0W0%wks(:,:,:,ik,:) - S_G0W0_DMFT%ws
                  enddo
               endif
               call DeallocateField(S_G0W0_DMFT)
               !
               !Add S_G0W0 to S_GW (already containing Simp)
               S_GW%wks = S_GW%wks + S_G0W0%wks
               !
               !Put together all the terms
               do ispin=1,Nspin
                  !$OMP PARALLEL DEFAULT(NONE),&
                  !$OMP SHARED(S_Full,S_GW,Glat,VH,Vxc,ispin),&
                  !$OMP PRIVATE(ik,iw)
                  !$OMP DO
                  do ik=1,S_Full%Nkpt
                     do iw=1,S_Full%Npoints
                        S_Full%wks(:,:,iw,ik,ispin) = + S_GW%wks(:,:,iw,ik,ispin) - Vxc(:,:,ik,ispin) + VH(:,:)! + 12d0*(Glat%wks(:,:,iw,ik,ispin)+1d0*deye(S_Full%Norb))
                     enddo
                  enddo
                  !$OMP END DO
                  !$OMP END PARALLEL
                  if(paramagnet)then
                     S_Full%wks(:,:,:,:,Nspin) = S_Full%wks(:,:,:,:,1)
                     exit
                  endif
               enddo
               !
            endif
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
               !$OMP PARALLEL DEFAULT(NONE),&
               !$OMP SHARED(HartreeFact,S_Full,S_DMFT,ispin),&
               !$OMP PRIVATE(ik,iw)
               !$OMP DO
               do ik=1,S_Full%Nkpt
                  do iw=1,S_Full%Npoints
                     S_Full%wks(:,:,iw,ik,ispin) = S_DMFT%ws(:,:,iw,ispin) - HartreeFact*S_DMFT%N_s(:,:,int(Nspin/ispin))
                  enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL
               S_Full%N_s(:,:,ispin) = HartreeFact*S_DMFT%N_s(:,:,int(Nspin/ispin))
               if(paramagnet)then
                  S_Full%wks(:,:,:,:,Nspin) = S_Full%wks(:,:,:,:,1)
                  exit
               endif
            enddo
            call FermionicKsum(S_Full)
            !
      end select
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
      call AllocateFermionicField(D_correction,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
      allocate(GS(S_GW%Norb,S_GW%Norb));GS=czero
      allocate(SG(S_GW%Norb,S_GW%Norb));SG=czero
      allocate(SGS(S_GW%Norb,S_GW%Norb));SGS=czero
      allocate(invG(S_GW%Norb,S_GW%Norb));invG=czero
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(S_GW,Glat,D_correction),&
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
            D_correction%ws(:,:,iw,ispin) = + SGS                             &
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
      call symmetrize_GW(D_correction,EqvGWndx)
      call dump_FermionicField(D_correction,reg(ItFolder),"D_correction_w",paramagnet)
      !
   end subroutine calc_causality_Delta_correction
   !
   subroutine calc_causality_curlyU_correction()
      !
      implicit none
      integer                               :: iq,iw
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
      allocate(WP(Plat%Nbp,Plat%Nbp));WP=czero
      allocate(PW(Plat%Nbp,Plat%Nbp));PW=czero
      allocate(PWP(Plat%Nbp,Plat%Nbp));PWP=czero
      allocate(W_P(Plat%Nbp,Plat%Nbp));W_P=czero
      allocate(P_W(Plat%Nbp,Plat%Nbp));P_W=czero
      allocate(invWa(Plat%Nbp,Plat%Nbp));invWa=czero
      allocate(invWb(Plat%Nbp,Plat%Nbp));invWb=czero
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Plat,Wlat,curlyU_correction),&
      !$OMP PRIVATE(iw,iq,WP,PW,PWP,W_P,P_W,invWa,invWb)
      !$OMP DO
      do iw=1,Plat%Npoints
         !
         WP=czero;PW=czero;PWP=czero
         !
         do iq=1,Plat%Nkpt
            !
            WP = WP + matmul(Wlat%screened(:,:,iw,iq),Plat%screened(:,:,iw,iq))/Plat%Nkpt
            PW = PW + matmul(Plat%screened(:,:,iw,iq),Wlat%screened(:,:,iw,iq))/Plat%Nkpt
            PWP = PWP + matmul(Plat%screened(:,:,iw,iq),matmul(Wlat%screened(:,:,iw,iq),Plat%screened(:,:,iw,iq)))/Plat%Nkpt
            !
         enddo
         !
         W_P=czero;P_W=czero
         W_P = matmul(Wlat%screened_local(:,:,iw),Plat%screened_local(:,:,iw))
         P_W = matmul(Plat%screened_local(:,:,iw),Wlat%screened_local(:,:,iw))
         !
         invWa=czero
         invWa = Plat%screened_local(:,:,iw) + PWP
         call inv(invWa)
         !
         invWb=czero
         invWb = Plat%screened_local(:,:,iw) + matmul(Plat%screened_local(:,:,iw),matmul(Wlat%screened_local(:,:,iw),Plat%screened_local(:,:,iw)))
         call inv(invWb)
         !
         ! Put toghether all the pieces. arxiv:2011.05311 Eq.22
         curlyU_correction%screened_local(:,:,iw) = + matmul(PW,matmul(invWa,WP)) - matmul(P_W,matmul(invWb,W_P))
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(WP,PW,PWP,W_P,P_W,invWa,invWb)
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
      type(FermionicField)                  :: DeltaOld
      integer                               :: Norb,Nflavor,unit
      integer                               :: ispin,iw,iwan,itau,ndx,wndx
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: wmats(:),tau(:),Moments(:,:,:)
      real(8),allocatable                   :: Eloc(:,:),PrintLine(:),coef01(:,:)
      complex(8),allocatable                :: Nloc(:,:,:)
      real(8)                               :: tailShift,CrystalField
      complex(8),allocatable                :: zeta(:,:,:),invG(:,:),Rot(:,:)
      complex(8),allocatable                :: Dfit(:,:,:),Dmats(:,:,:),Ditau(:,:,:)
      complex(8),allocatable                :: invCurlyG(:,:,:)!,DeltaTail(:)
      character(len=255)                    :: file,MomDir
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Delta of "//reg(SiteName(isite))
      !
      !
      if(.not.S_DMFT%status) stop "calc_Delta: S_DMFT not properly initialized."
      if(.not.Glat%status) stop "calc_Delta: Glat not properly initialized."
      if(causal_D.and.(.not.D_correction%status)) stop "calc_Delta: requested causality correction but D_correction not properly initialized."
      !
      Norb = SiteNorb(isite)
      allocate(Orbs(Norb))
      Orbs = SiteOrbs(isite,1:Norb)
      Nflavor = Norb*Nspin
      !
      allocate(Eloc(Norb,Nspin));Eloc=0d0
      allocate(coef01(Norb,Nspin));coef01=0d0
      !
      call AllocateFermionicField(SigmaImp,Norb,Nmats,Beta=Beta)
      call AllocateFermionicField(Gloc,Norb,Nmats,Beta=Beta)
      if(causal_D)call AllocateFermionicField(DeltaCorr,Norb,Nmats,Beta=Beta)
      allocate(invCurlyG(Norb,Nmats,Nspin));invCurlyG=czero
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(Beta,Nmats)
      allocate(zeta(Norb,Nmats,Nspin))
      do iwan=1,Norb
         do iw=1,Nmats
            zeta(iwan,iw,:) = dcmplx( Glat%mu , wmats(iw) )
         enddo
      enddo
      !
      !Extract and rotate from local (non-diagonal) to imp (diagonal) the given sites
      !similar result could be obtained with S_Full%ws(:,:,iw,:) + S_Full%N_s
      !defined in join_SigmaFull but using directly S_DMFT is more precise
      allocate(Nloc(Norb,Norb,Nspin));Nloc=czero
      call clear_attributes(Gloc)
      call clear_attributes(SigmaImp)
      if(causal_D)call clear_attributes(DeltaCorr)
      if(RotateHloc)then
         !
         allocate(Rot(Norb,Norb)); Rot=OlocRot(1:Norb,1:Norb,isite)
         call loc2imp(Gloc,Glat,Orbs,U=Rot)
         call calc_density(Gloc,Nloc)
         call loc2imp(SigmaImp,S_DMFT,Orbs,U=Rot)
         if(causal_D)call loc2imp(DeltaCorr,D_correction,Orbs,U=Rot)
         deallocate(Rot)
         !
         if(sym_mode.gt.1)then
            call symmetrize_imp(Gloc,OlocEig(:,isite))
            call symmetrize_imp(SigmaImp,OlocEig(:,isite))
            if(causal_D)call symmetrize_imp(DeltaCorr,OlocEig(:,isite))
         endif
         !
      else
         !
         call loc2imp(Gloc,Glat,Orbs)
         call calc_density(Gloc,Nloc)
         call loc2imp(SigmaImp,S_DMFT,Orbs)
         if(causal_D)call loc2imp(DeltaCorr,D_correction,Orbs)
         !
      endif
      !
      !Print what's used to compute delta
      call dump_Matrix(Nloc,reg(ItFolder),"Solver_"//reg(SiteName(isite))//"/N_"//reg(SiteName(isite)),paramagnet)
      deallocate(Nloc)
      call dump_FermionicField(Gloc,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","G_"//reg(SiteName(isite))//"_w",paramagnet)
      call dump_FermionicField(SigmaImp,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","S_"//reg(SiteName(isite))//"_w",paramagnet)
      if(causal_D)call dump_FermionicField(DeltaCorr,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","D_correction_"//reg(SiteName(isite))//"_w",paramagnet)
      !
      !Compute the fermionic Weiss field aka the inverse of CurlyG
      allocate(invG(Norb,Norb));invG=czero
      do ispin=1,Nspin
         do iw=1,Nmats
            !
            invG = Gloc%ws(:,:,iw,ispin)
            call inv(invG)
            !
            if(causal_D)then
               do iwan=1,Norb
                  invCurlyG(iwan,iw,ispin) = invG(iwan,iwan) + SigmaImp%ws(iwan,iwan,iw,ispin) - DeltaCorr%ws(iwan,iwan,iw,ispin)
               enddo
            else
               do iwan=1,Norb
                  invCurlyG(iwan,iw,ispin) = invG(iwan,iwan) + SigmaImp%ws(iwan,iwan,iw,ispin)
               enddo
            endif
            !
         enddo
      enddo
      deallocate(invG)
      call DeallocateFermionicField(SigmaImp)
      call DeallocateFermionicField(Gloc)
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
         case("Inf")
            !
            Eloc = real(Dfit(:,Nmats,:))
            !
         case("Analytic")
            !
            file = "DeltaPara_"//reg(SiteName(isite))//".DAT"
            MomDir = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/"
            wndx = minloc(abs(wmats-0.85*wmatsMax),dim=1)
            call fit_Delta(Dfit,Beta,Nfit,reg(MomDir),reg(file),"Shifted",Eloc,filename="DeltaAnd",Wlimit=wndx,coef01=coef01)
            !
         case("Moments")
            !
            file = "DeltaMom_"//reg(SiteName(isite))//".DAT"
            MomDir = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/"
            !
            allocate(Moments(Norb,Nspin,0:4));Moments=0d0
            wndx = minloc(abs(wmats-0.85*wmatsMax),dim=1)
            call fit_moments(Dfit,Beta,reg(MomDir),reg(file),"Sigma",Moments,filename="DeltaMom",Wlimit=wndx)
            !
            Eloc=Moments(:,0,:)
            coef01=Moments(:,1,:)
            !
            deallocate(Moments)
            !
      end select
      deallocate(Dfit)
      !
      !Compute Delta on matsubara
      allocate(Dmats(Norb,Nmats,Nspin));Dmats=czero
      do ispin=1,Nspin
         do iwan=1,Norb
            !
            do iw=1,Nmats
               Dmats(iwan,iw,ispin) = dcmplx( Glat%mu , wmats(iw) ) - Eloc(iwan,ispin) - invCurlyG(iwan,iw,ispin)
            enddo
            !
            !Additional correction
            if(dreal(Dmats(iwan,Nmats,ispin))*real(Dmats(iwan,Nmats-1,ispin)).lt.0d0)then
               tailShift = real(Dmats(iwan,Nmats,ispin)) + (Dmats(iwan,Nmats,ispin)-Dmats(iwan,Nmats-1,ispin))*2d0
               write(*,"(A,E13.3)")"     Additional Eloc shift orb: "//str(iwan)//" spin: "//str(ispin),tailShift
               Eloc(iwan,ispin) = Eloc(iwan,ispin) + tailShift
               Dmats(iwan,:,ispin) = Dmats(iwan,:,ispin) - tailShift
            endif
            !
         enddo
      enddo
      deallocate(invCurlyG)
      !
      !Mixing Delta
      if((Mixing_Delta.gt.0d0).and.(Iteration.gt.0))then
         write(*,"(A)")"     Mixing Delta with "//str(Mixing_Delta,3)//" of old solution."
         call AllocateFermionicField(DeltaOld,Norb,Nmats,Beta=Beta)
         call read_FermionicField(DeltaOld,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","D_"//reg(SiteName(isite))//"_w")
         do ispin=1,Nspin
            do iw=1,Nmats
               do iwan=1,Norb
                  Dmats(iwan,iw,ispin) = (1d0-Mixing_Delta)*Dmats(iwan,iw,ispin) + Mixing_Delta*DeltaOld%ws(iwan,iwan,iw,ispin)
               enddo
            enddo
         enddo
         call DeallocateFermionicField(DeltaOld)
      endif
      !
      !Spin symmetrization for paramagnetic calculations
      if(EqvGWndx%S)then
         Dmats(:,:,1) = (Dmats(:,:,1) + Dmats(:,:,2))/2d0
         Dmats(:,:,2) = Dmats(:,:,1)
         Eloc(:,1) = (Eloc(:,1) + Eloc(:,2))/2d0
         Eloc(:,2) = Eloc(:,1)
      endif
      !
      !Fourier transform to the tau axis
      allocate(Ditau(Norb,Solver%NtauF,Nspin));Ditau=0d0
      do ispin=1,Nspin
         !add correction that allows to use the tail correction in the FT
         do iwan=1,Norb
            if(reg(CalculationType).ne."GW+EDMFT") coef01(iwan,ispin) = abs(dimag(Dmats(iwan,Nmats,ispin)))*wmats(Nmats)
            write(*,"(A,E10.3)")"     First coeff orb: "//str(iwan)//" spin: "//str(ispin),coef01(iwan,ispin)
            Dmats(iwan,:,ispin) = Dmats(iwan,:,ispin)/coef01(iwan,ispin)
         enddo
         !FT
         call Fmats2itau_vec(Beta,Dmats(:,:,ispin),Ditau(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
         !remove correction
         do iwan=1,Norb
            Dmats(iwan,:,ispin) = Dmats(iwan,:,ispin)*coef01(iwan,ispin)
            Ditau(iwan,:,ispin) = Ditau(iwan,:,ispin)*coef01(iwan,ispin)
         enddo
      enddo
      !
      !Check on Delta(tau) causality
      do ispin=1,Nspin
         do iwan=1,Norb
            do itau=1,Solver%NtauF
               if(dreal(Ditau(iwan,itau,ispin)).gt.0d0)then
                  write(*,"(A,E10.3)")"     Warning: Removing non-causality from Delta(tau) at orb: "//str(iwan)//" spin: "//str(ispin)//" itau: "//str(itau)
                  if(dreal(Ditau(iwan,int(Solver%NtauF/2),ispin)).lt.0d0)then
                     Ditau(iwan,itau,ispin) = Ditau(iwan,int(Solver%NtauF/2),ispin)
                  else
                     Ditau(iwan,itau,ispin) = czero
                  endif
               endif
            enddo
         enddo
      enddo
      !
      !Eloc and chemical potential
      if(any(SiteCF(isite,:).ne.0d0))then
         file = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Eloc_noCF.DAT"
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
         write(unit,"(1E20.12)") Glat%mu
         do iwan=1,Norb
            write(unit,"(2E20.12)") Eloc(iwan,1)+EqvGWndx%hseed,Eloc(iwan,2)-EqvGWndx%hseed
         enddo
         close(unit)
      endif
      file = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Eloc.DAT"
      unit = free_unit()
      open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(1E20.12)") Glat%mu
      CrystalField=0d0
      do iwan=1,Norb
         if(iwan.gt.1)CrystalField=SiteCF(isite,iwan-1)
         write(unit,"(2E20.12)") Eloc(iwan,1)+EqvGWndx%hseed+CrystalField,Eloc(iwan,2)-EqvGWndx%hseed+CrystalField
      enddo
      close(unit)
      !
      !Delta(tau)
      allocate(tau(Solver%NtauF));tau=0d0
      tau = linspace(0d0,Beta,Solver%NtauF)
      allocate(PrintLine(Nflavor))
      file = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Delta_t.DAT"
      unit = free_unit()
      open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
      do itau=1,Solver%NtauF
         ndx=1
         PrintLine=0d0
         do iwan=1,Norb
            do ispin=1,Nspin
               PrintLine(ndx) = real(Ditau(iwan,itau,ispin))
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
            do iwan=1,Norb
               FermiPrint%ws(iwan,iwan,iw,ispin) = 1d0/(dcmplx( Glat%mu , wmats(iw) ) - Eloc(iwan,ispin) - Dmats(iwan,iw,ispin))
            enddo
         enddo
      enddo
      FermiPrint%mu=Glat%mu
      call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w",paramagnet)
      !
      !Delta(iw)
      call clear_attributes(FermiPrint)
      do iwan=1,Norb
         FermiPrint%ws(iwan,iwan,:,:) = Dmats(iwan,:,:)
      enddo
      FermiPrint%mu=Glat%mu
      call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","D_"//reg(SiteName(isite))//"_w",paramagnet)
      !
      !Delta(iw) without causality correction
      if(causal_D)then
         call clear_attributes(FermiPrint)
         do ispin=1,Nspin
            do iwan=1,Norb
               FermiPrint%ws(iwan,iwan,:,ispin) = Dmats(iwan,:,ispin) - DeltaCorr%ws(iwan,iwan,:,ispin)
            enddo
         enddo
         call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","D_notCorr_"//reg(SiteName(isite))//"_w",paramagnet)
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
            do iwan=1,Norb
               FermiPrint%ws(iwan,iwan,:,ispin) = Dmats(iwan,:,ispin)
            enddo
         enddo
         call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","testFT_D_"//reg(SiteName(isite))//"_w",paramagnet)
         !
         call clear_attributes(FermiPrint)
         do ispin=1,Nspin
            do iwan=1,Norb
               FermiPrint%ws(iwan,iwan,:,ispin) = Eloc(iwan,ispin) + Dmats(iwan,:,ispin)
            enddo
         enddo
         call dump_FermionicField(FermiPrint,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/fits/","test_Eo+D_"//reg(SiteName(isite))//"_w",paramagnet)
         !
      endif
      call DeallocateFermionicField(FermiPrint)
      !
      deallocate(Orbs,Eloc,zeta,Dmats,Ditau,tau,wmats)
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
      integer                               :: Norb,Nflavor,Nbp
      integer                               :: ib1,ib2,itau,iw
      integer                               :: unit,ndx,isitecheck
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: Uinst(:,:),Ucheck(:,:)
      real(8),allocatable                   :: Kfunct(:,:,:)
      real(8),allocatable                   :: tau(:),PrintLine(:)
      complex(8),allocatable                :: Rot(:,:)
      character(len=255)                    :: file
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Interaction of "//reg(SiteName(isite))
      !
      !
      if(.not.P_EDMFT%status) stop "calc_Interaction: P_EDMFT not properly initialized."
      if(.not.Wlat%status) stop "calc_Interaction: Wlat not properly initialized."
      if(causal_U.and.(.not.curlyU_correction%status)) stop "calc_Interaction: requested causality correction but curlyU_correction not properly initialized."
      !
      Norb = SiteNorb(isite)
      allocate(Orbs(Norb))
      Orbs = SiteOrbs(isite,1:Norb)
      Nbp = Norb**2
      Nflavor = Norb*Nspin
      !
      allocate(tau(Solver%NtauB));tau=0d0
      tau = linspace(0d0,Beta,Solver%NtauB)
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      allocate(Uinst(Nflavor,Nflavor));Uinst=0d0
      call AllocateBosonicField(curlyU,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
      !
      select case(reg(CalculationType))
         case default
            !
            stop "If you got so far somethig is wrong."
            !
         case("DMFT+statU")
            !
            call loc2imp(curlyU,Ulat,Orbs)
            call isReal(curlyU)
            call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
            !
            if(RotateUloc)then
               allocate(Rot(Norb,Norb)); Rot=OlocRot(1:Norb,1:Norb,isite)
               call TransformBosonicField(curlyU,Rot,PhysicalUelements%Full_Map)
               deallocate(Rot)
               call isReal(curlyU)
               call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
            endif
            !
            call calc_QMCinteractions(curlyU,Uinst)
            !
         case("DMFT+dynU")
            !
            call loc2imp(curlyU,Ulat,Orbs)
            call isReal(curlyU)
            call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
            !
            if(RotateUloc)then
               allocate(Rot(Norb,Norb)); Rot=OlocRot(1:Norb,1:Norb,isite)
               call TransformBosonicField(curlyU,Rot,PhysicalUelements%Full_Map)
               deallocate(Rot)
               call isReal(curlyU)
               call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
            endif
            !
            allocate(Kfunct(Nflavor,Nflavor,Solver%NtauB));Kfunct=0d0
            call calc_QMCinteractions(curlyU,Uinst,Kfunct)
            !
         case("EDMFT","GW+EDMFT")
            !
            if(Ustart)then
               !
               write(*,"(A)") "     Using local Ucrpa as effective interaction."
               call loc2imp(curlyU,Ulat,Orbs)
               call isReal(curlyU)
               call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
               !
               if(RotateUloc)then
                  call dump_BosonicField(curlyU,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_notRot_"//reg(SiteName(isite))//"_w.DAT") !remove is the same of curlyU_EDMFT
                  allocate(Rot(Norb,Norb)); Rot=OlocRot(1:Norb,1:Norb,isite)
                  call TransformBosonicField(curlyU,Rot,PhysicalUelements%Full_Map)
                  deallocate(Rot)
                  call isReal(curlyU)
                  call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
               endif
               !
            else
               !
               write(*,"(A)") "     Computing the local effective interaction."
               call AllocateBosonicField(Pimp,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
               call AllocateBosonicField(Wloc,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
               if(causal_U)call AllocateBosonicField(curlyUcorr,Norb,Nmats,Crystal%iq_gamma,Beta=Beta,no_bare=.true.)
               !
               call loc2imp(Pimp,P_EDMFT,Orbs)
               call loc2imp(Wloc,Wlat,Orbs)
               !
               if(causal_U)then
                  call loc2imp(curlyUcorr,curlyU_correction,Orbs)
                  call calc_curlyU(curlyU,Wloc,Pimp,curlyUcorr=curlyUcorr)
               else
                  call calc_curlyU(curlyU,Wloc,Pimp)
               endif
               call isReal(curlyU)
               call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
               if(causal_U)call clear_MatrixElements(curlyUcorr,PhysicalUelements%Full_All)
               !
               if(RotateUloc)then
                  call dump_BosonicField(curlyU,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_notRot_"//reg(SiteName(isite))//"_w.DAT") !remove is the same of curlyU_EDMFT
                  allocate(Rot(Norb,Norb)); Rot=OlocRot(1:Norb,1:Norb,isite)
                  call TransformBosonicField(curlyU,Rot,PhysicalUelements%Full_Map)
                  if(causal_U)call TransformBosonicField(curlyUcorr,Rot,PhysicalUelements%Full_Map)
                  deallocate(Rot)
                  call isReal(curlyU)
                  call clear_MatrixElements(curlyU,PhysicalUelements%Full_All)
               endif
               !
               if(causal_U)then
                  call dump_BosonicField(curlyUcorr,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_correction_"//reg(SiteName(isite))//"_w.DAT")
                  call DeallocateBosonicField(curlyUcorr)
               endif
               !
               call DeallocateBosonicField(Pimp)
               call DeallocateBosonicField(Wloc)
               !
            endif
            !
            !Mixing curlyU
            if((Mixing_curlyU.gt.0d0).and.(Iteration.gt.0))then
               write(*,"(A)")"     Mixing curlyU with "//str(Mixing_curlyU,3)//" of old solution."
               call AllocateBosonicField(curlyUold,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
               call read_BosonicField(curlyUold,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
               curlyU%bare_local = (1d0-Mixing_curlyU)*curlyU%bare_local + Mixing_curlyU*curlyUold%bare_local
               do iw=1,Nmats
                  curlyU%screened_local(:,:,iw) = (1d0-Mixing_curlyU)*curlyU%screened_local(:,:,iw) + Mixing_curlyU*curlyUold%screened_local(:,:,iw)
               enddo
               call DeallocateBosonicField(curlyUold)
            endif
            !
            allocate(Kfunct(Nflavor,Nflavor,Solver%NtauB));Kfunct=0d0
            call calc_QMCinteractions(curlyU,Uinst,Kfunct)
            !
      end select
      deallocate(Orbs)
      !
      !Print curlyU in the solver basis
      call dump_BosonicField(curlyU,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
      !
      !Istantaneous interaction
      call dump_Matrix(Uinst,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","Umat.DAT")
      !
      !Print data for retarded interactions
      if(allocated(Kfunct))then
         !
         !Screening function
         file = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/K_t.DAT"
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="unknown",position="rewind",action="write")
         allocate(PrintLine(Nflavor*(Nflavor+1)/2));PrintLine=0d0
         do itau=1,Solver%NtauB
            ndx=1
            !print diagonal and LT
            do ib1=1,Nflavor
               do ib2=1,ib1
                  !
                  PrintLine(ndx) = Kfunct(ib1,ib2,itau)
                  if(Kdiag) PrintLine(ndx) = Kfunct(ib1,ib1,itau)
                  ndx=ndx+1
                  !
               enddo
            enddo
            write(unit,"(999E20.12)") tau(itau),PrintLine
         enddo
         deallocate(PrintLine)
         close(unit)
         !
         deallocate(Kfunct)
         !
      endif
      deallocate(tau,Uinst)
      !
      !Check if the screened effective local interaction at iw=0 is the same for all the sites. Same Orbital dimension assumed.
      if(checkInvariance)then
         !
         do isitecheck=1,Nsite
            !
            Norb = SiteNorb(isitecheck)
            Nflavor = Nspin*Norb
            allocate(Orbs(Norb))
            Orbs = SiteOrbs(isitecheck,1:Norb)
            !
            write(*,"(A)")"     Checking site "//reg(SiteName(1))//"_"//str(isitecheck)
            write(*,"(A,10I3)")"     Norb: "//str(Norb)//" Orbitals: ",Orbs
            !
            if(isitecheck.eq.1)allocate(Uinst(Nflavor,Nflavor))
            if(isitecheck.ne.1)allocate(Ucheck(Nflavor,Nflavor))
            !
            call clear_attributes(curlyU)
            !
            if(Ustart)then
               call loc2imp(curlyU,Ulat,Orbs)
               call isReal(curlyU)
            else
               call AllocateBosonicField(Pimp,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
               call AllocateBosonicField(Wloc,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
               call loc2imp(Pimp,P_EDMFT,Orbs)
               call loc2imp(Wloc,Wlat,Orbs)
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
               call dump_Matrix(Uinst,reg(ItFolder)//"Solver_"//reg(SiteName(isitecheck))//"/","Umat_noSym_"//str(isitecheck)//".DAT")
            else
               call calc_QMCinteractions(curlyU,Ucheck,sym=.false.)
               call dump_Matrix(Ucheck,reg(ItFolder)//"Solver_"//reg(SiteName(isitecheck))//"/","Umat_noSym_"//str(isitecheck)//".DAT")
               !
               do ib1=1,Nflavor
                  do ib2=1,Nflavor
                     if(abs(Ucheck(ib1,ib2)-Uinst(ib1,ib2)).gt.1e-3)then
                        write(*,"(A,F,A,F)")"     Warning: Element["//str(ib1)//"]["//str(ib2)//"] is different:",Ucheck(ib1,ib2)," instead of: ",Uinst(ib1,ib2)
                     endif
                  enddo
               enddo
               deallocate(Ucheck)
               !
            endif
            !
            deallocate(Orbs)
            !
         enddo !isitecheck
         deallocate(Uinst)
         !
      endif
      !
      call DeallocateBosonicField(curlyU)
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
      integer                               :: Norb,Nflavor,Nbp
      integer                               :: iorb,jorb,ispin,jspin
      integer                               :: ib1,ib2,isite,idum
      integer                               :: unit,ndx,itau,iw,wndx,wndxOpt
      real(8)                               :: taup,muQMC
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: tauF(:),tauB(:),wmats(:)
      real(8),allocatable                   :: ReadLine(:)
      real(8),allocatable                   :: Moments(:,:,:)
      character(len=255)                    :: file,MomDir
      logical                               :: filexists
      type(physicalU)                       :: PhysicalUelements
      !Impurity Green's function
      type(FermionicField)                  :: Gimp
      complex(8),allocatable                :: Gitau(:,:,:)
      complex(8),allocatable                :: Gmats(:,:,:,:)
      !Impurity self-energy and fermionic Dyson equation
      type(FermionicField)                  :: Simp
      type(FermionicField)                  :: G0imp
      type(BosonicField)                    :: curlyU
      complex(8),allocatable                :: Smats(:,:,:,:)
      complex(8),allocatable                :: SmatsTail(:),SmatsNoFit(:,:,:,:)
      !Impurity susceptibilities
      real(8),allocatable                   :: nnt(:,:,:)
      complex(8),allocatable                :: NNitau(:,:,:,:,:)
      complex(8),allocatable                :: ChiMitau(:),ChiMmats(:)
      type(BosonicField)                    :: ChiCitau,ChiCmats
      !Impurity polarization and bosonic Dyson equation
      type(BosonicField)                    :: Pimp
      type(BosonicField)                    :: Wimp
      real(8),allocatable                   :: CDW(:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- collect_QMC_results"
      !
      !
      !
      !
      ! COLLECT IMPURITY OCCUPATION --------------------------------------------
      allocate(densityDMFT(Crystal%Norb,Crystal%Norb,Nspin));densityDMFT=czero
      allocate(densityQMC(maxval(SiteNorb),maxval(SiteNorb),Nspin,Nsite));densityQMC=0d0
      !
      do isite=1,Nsite
         !
         write(*,"(A)") new_line("A")//"     Collecting occupation of site: "//reg(SiteName(isite))
         !
         Norb = SiteNorb(isite)
         allocate(Orbs(Norb))
         Orbs = SiteOrbs(isite,1:Norb)
         Nflavor = Norb*Nspin
         !
         !
         !Read the impurity occupation
         file = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/Nqmc.DAT"
         call inquireFile(reg(file),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
         read(unit,*) muQMC
         do ib1=1,Nflavor
            iorb = (ib1+mod(ib1,2))/2
            ispin = abs(mod(ib1,2)-2)
            read(unit,*) idum,densityQMC(iorb,iorb,ispin,isite)
         enddo
         close(unit)
         !
         !
         !Insert or Expand to the Lattice basis
         if(RotateHloc)then
            if(sym_mode.gt.1)call symmetrize_imp(densityQMC(:,:,:,isite),OlocEig(:,isite))
            call imp2loc(densityDMFT,dcmplx(densityQMC(1:Norb,1:Norb,:,isite),0d0),isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag)
         else
            call imp2loc(densityDMFT,dcmplx(densityQMC(1:Norb,1:Norb,:,isite),0d0),isite,Orbs,ExpandImpurity,AFMselfcons)
         endif
         !
         deallocate(Orbs)
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      !symmetrize GW indexes and print
      if(EqvGWndx%O.or.EqvGWndx%S)then
         !
         if(verbose)call dump_Matrix(densityDMFT,reg(PrevItFolder),"Nimp_noSym",paramagnet)
         !
         call symmetrize_GW(densityDMFT,EqvGWndx)
         !
      endif
      call dump_Matrix(densityDMFT,reg(PrevItFolder),"Nimp",paramagnet)
      deallocate(densityDMFT)
      !
      !
      !
      !
      ! COLLECT IMPURITY GF AND FERMIONIC DYSON --------------------------------
      allocate(Gmats(maxval(SiteNorb),Nmats,Nspin,Nsite));Gmats=czero
      do isite=1,Nsite
         !
         write(*,"(A)") new_line("A")//"     Collecting the impurity Green's function of site: "//reg(SiteName(isite))
         !
         Norb = SiteNorb(isite)
         allocate(Orbs(Norb))
         Orbs = SiteOrbs(isite,1:Norb)
         !
         !Read the impurity Green's function
         allocate(tauF(Solver%NtauF_in));tauF = linspace(0d0,Beta,Solver%NtauF_in)
         allocate(Gitau(Norb,Solver%NtauF_in,Nspin));Gitau=czero
         allocate(ReadLine(Nspin*Norb))
         file = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/Gimp_t.DAT"
         call inquireFile(reg(file),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
         do itau=1,Solver%NtauF_in
            ReadLine=0d0
            read(unit,*) taup,ReadLine
            if(abs(taup-tauF(itau)).gt.eps) stop "Impurity fermionic tau mesh does not coincide."
            !
            ndx=1
            do iorb=1,Norb
               do ispin=1,Nspin
                  Gitau(iorb,itau,ispin) = dcmplx(ReadLine(ndx),0d0)
                  ndx=ndx+1
               enddo
            enddo
            !
         enddo
         deallocate(ReadLine,tauF)
         !
         !FT to the matsubara axis
         do ispin=1,Nspin
            call Fitau2mats_vec(Beta,Gitau(:,:,ispin),Gmats(1:Norb,:,ispin,isite),tau_uniform=.true.)
         enddo
         !
         call dump_MaxEnt(Gitau,"itau",reg(PrevItFolder)//"Convergence/","Gqmc_"//reg(SiteName(isite)))
         deallocate(Gitau)
         !
         deallocate(Orbs)
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      !Save to file in standard format Gimp
      call AllocateFermionicField(G_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
      do isite=1,Nsite
         !
         Norb = SiteNorb(isite)
         allocate(Orbs(Norb))
         Orbs = SiteOrbs(isite,1:Norb)
         call AllocateFermionicField(Gimp,Norb,Nmats,Beta=Beta)
         !
         !Push Gimp to a Local fermionic field
         call clear_attributes(Gimp)
         do iorb=1,Norb
            Gimp%ws(iorb,iorb,:,:) = Gmats(iorb,:,:,isite)
         enddo
         !
         !Expand to the Lattice basis
         if(RotateHloc)then
            if(sym_mode.gt.1)call symmetrize_imp(Gimp,OlocEig(:,isite))
            call imp2loc(G_DMFT,Gimp,isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag)
         else
            call imp2loc(G_DMFT,Gimp,isite,Orbs,ExpandImpurity,AFMselfcons)
         endif
         !
         deallocate(Orbs)
         call DeallocateFermionicField(Gimp)
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      call dump_FermionicField(G_DMFT,reg(PrevItFolder),"Gimp_w",paramagnet)
      call dump_MaxEnt(G_DMFT,"mats",reg(PrevItFolder)//"Convergence/","Gimp",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call dump_MaxEnt(G_DMFT,"mats2itau",reg(PrevItFolder)//"Convergence/","Gimp",EqvGWndx%SetOrbs)
      call DeallocateFermionicField(G_DMFT)
      !
      !Fermionic Dyson equation
      allocate(Smats(maxval(SiteNorb),Nmats,Nspin,Nsite));Smats=czero
      do isite=1,Nsite
         !
         Norb = SiteNorb(isite)
         allocate(Orbs(Norb))
         Orbs = SiteOrbs(isite,1:Norb)
         !
         write(*,"(A)") new_line("A")//"     Collecting curlyG of site: "//reg(SiteName(isite))
         !
         !Read curlyG
         call AllocateFermionicField(G0imp,Norb,Nmats,Beta=Beta)
         call read_FermionicField(G0imp,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w")
         !
         !Adjust with the chemical potential if the solver has changed it
         if((abs(G0imp%mu-muQMC).gt.1e-5).and.update_curlyG)then
            write(*,"(A)") new_line("A")//"     Updating the chemical potential of curlyG from "//str(G0imp%mu)//" to "//str(muQMC)
            do ispin=1,Nspin
               do iw=1,Nmats
                  do iorb=1,Norb
                     G0imp%ws(iorb,iorb,iw,ispin) = 1d0/(1d0/G0imp%ws(iorb,iorb,iw,ispin) - G0imp%mu + muQMC)
                  enddo
               enddo
            enddo
            G0imp%mu=muQMC
            call dump_FermionicField(G0imp,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w",paramagnet)
         endif
         !
         !Fermionic Dyson equation in the solver basis (always diagonal)
         write(*,"(A)") new_line("A")//"     Solving fermionic Dyson of site: "//reg(SiteName(isite))
         do ispin=1,Nspin
            do iorb=1,Norb
               Smats(iorb,:,ispin,isite) = 1d0/G0imp%ws(iorb,iorb,:,ispin) - 1d0/Gmats(iorb,:,ispin,isite)
            enddo
         enddo
         call DeallocateFermionicField(G0imp)
         !
         deallocate(Orbs)
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      deallocate(Gmats)
      !
      !Replace the tail with the fitted self-energy - I'm doing another loop over sites because I want to store all the SmatsNoFit
      if((ReplaceTail_Simp.ne.0d0).and.(ReplaceTail_Simp.lt.wmatsMax))then
         !
         allocate(SmatsNoFit(maxval(SiteNorb),Nmats,Nspin,Nsite))
         SmatsNoFit=Smats
         !
         do isite=1,Nsite
            !
            write(*,"(A)") new_line("A")//"     Fitting self-energy moments of site: "//reg(SiteName(isite))
            !
            Norb = SiteNorb(isite)
            allocate(Orbs(Norb))
            Orbs = SiteOrbs(isite,1:Norb)
            !
            !perform the fit
            file = "SimpMom_"//reg(SiteName(isite))//".DAT"
            MomDir = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/"
            allocate(wmats(Nmats));wmats=FermionicFreqMesh(Beta,Nmats)
            !
            !define the frequency index from which substitute the tail
            wndx = minloc(abs(wmats-ReplaceTail_Simp),dim=1)
            wndxOpt = wndx
            do ispin=1,Nspin
               do iorb=1,Norb
                  wndxOpt = min(wndxOpt,minloc(abs(SmatsNoFit(iorb,:,ispin,isite)-0.3d0*minval(dimag(SmatsNoFit(iorb,:,ispin,isite)))),dim=1))
               enddo
            enddo
            !
            allocate(Moments(Norb,Nspin,0:min(SigmaMaxMom,Nfit)));Moments=0d0
            call fit_moments(Smats(:,:,:,isite),Beta,reg(MomDir),reg(file),"Sigma",Moments,filename="Simp",Wlimit=wndx)
            !
            allocate(SmatsTail(Nmats));SmatsTail=czero
            write(*,"(A,F)") new_line("A")//"     Replacing Sigma tail starting from iw_["//str(wndx)//"]=",wmats(wndx)
            if(abs(wmats(wndx)-wmats(wndxOpt)).gt.1d0)write(*,"(A,F)") "     Optimal would be iw_["//str(wndxOpt)//"]=",wmats(wndxOpt)
            do ispin=1,Nspin
               do iorb=1,Norb
                  SmatsTail = S_Moments(Moments(iorb,ispin,:),wmats)
                  Smats(iorb,wndx:Nmats,ispin,isite) = SmatsTail(wndx:Nmats)
               enddo
            enddo
            deallocate(wmats,Moments,SmatsTail)
            !
            deallocate(Orbs)
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      endif
      !
      !Save to file in standard format the eventually fitted self-energy
      call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
      do isite=1,Nsite
         !
         Norb = SiteNorb(isite)
         allocate(Orbs(Norb))
         Orbs = SiteOrbs(isite,1:Norb)
         call AllocateFermionicField(Simp,Norb,Nmats,Beta=Beta)
         !
         !Push Simp to a Local fermionic field
         call clear_attributes(Simp)
         do iorb=1,Norb
            Simp%ws(iorb,iorb,:,:) = Smats(iorb,:,:,isite)
         enddo
         !
         !Fill up the N_s attribute that correspond to the Hartree term
         !since the impurity interaction contains only the physical interaction elements
         !(particle and spin conserving) I'm neglecting here elements in principle present
         !Because both Simp and densityQMC are diagonal
         !like Hartree_{ab} = Sum_c curlyU_{ab}{cc} * n_{cc}.
         call AllocateBosonicField(curlyU,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
         call read_BosonicField(curlyU,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
         Simp%N_s=czero
         do ispin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  !
                  ib1 = iorb + Norb*(iorb-1)
                  ib2 = jorb + Norb*(jorb-1)
                  !
                  Simp%N_s(iorb,iorb,ispin) = Simp%N_s(iorb,iorb,ispin) + real(curlyU%screened_local(ib1,ib2,1))*densityQMC(jorb,jorb,ispin,isite)
                  !
               enddo
            enddo
         enddo
         call DeallocateBosonicField(curlyU)
         !
         !The magnetization will be given only by the self-energy
         Simp%N_s(:,:,1) = (Simp%N_s(:,:,1)+Simp%N_s(:,:,2))/2d0
         Simp%N_s(:,:,2) = Simp%N_s(:,:,1)
         !
         call dump_Matrix(Simp%N_s,reg(PrevItFolder),"Solver_"//reg(SiteName(isite))//"/HartreeU",paramagnet)
         !
         !Expand to the Lattice basis
         if(RotateHloc)then
            if(sym_mode.gt.1)call symmetrize_imp(Simp,OlocEig(:,isite))
            call imp2loc(S_DMFT,Simp,isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag)
         else
            call imp2loc(S_DMFT,Simp,isite,Orbs,ExpandImpurity,AFMselfcons)
         endif
         !
         call DeallocateFermionicField(Simp)
         deallocate(Orbs)
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      !symmetrize GW indexes and print
      if(EqvGWndx%O.or.EqvGWndx%S)then
         !
         if(verbose)then
            call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_noSym_w",paramagnet)
            call dump_Matrix(S_DMFT%N_s,reg(PrevItFolder),"HartreeU_noSym",paramagnet)
         endif
         !
         call symmetrize_GW(S_DMFT,EqvGWndx)
         !
      endif
      deallocate(Smats)
      !
      call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_w",paramagnet)
      call dump_MaxEnt(S_DMFT,"mats",reg(PrevItFolder)//"Convergence/","Simp",EqvGWndx%SetOrbs,WmaxPade=PadeWlimit)
      call dump_Matrix(S_DMFT%N_s,reg(PrevItFolder),"HartreeU",paramagnet)
      call DeallocateFermionicField(S_DMFT)
      !
      !Save the non-Fitted non-symmetrized self-energy if present
      if(ReplaceTail_Simp.ne.0d0)then
         !
         if(.not.allocated(SmatsNoFit)) stop "SmatsNoFit is not allocated."
         call AllocateFermionicField(S_DMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
         do isite=1,Nsite
            !
            Norb = SiteNorb(isite)
            allocate(Orbs(Norb))
            Orbs = SiteOrbs(isite,1:Norb)
            call AllocateFermionicField(Simp,Norb,Nmats,Beta=Beta)
            !
            !Push Simp to a Local fermionic field
            call clear_attributes(Simp)
            do iorb=1,Norb
               Simp%ws(iorb,iorb,:,:) = SmatsNoFit(iorb,:,:,isite)
            enddo
            !
            !Expand to the Lattice basis
            if(RotateHloc)then
               call imp2loc(S_DMFT,Simp,isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag)
            else
               call imp2loc(S_DMFT,Simp,isite,Orbs,ExpandImpurity,AFMselfcons)
            endif
            !
            call DeallocateFermionicField(Simp)
            deallocate(Orbs)
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         deallocate(SmatsNoFit)
         call dump_FermionicField(S_DMFT,reg(PrevItFolder),"Simp_noFit_w",paramagnet)
         call DeallocateFermionicField(S_DMFT)
         !
      endif
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
            write(*,"(A)") new_line("A")//"     Collecting the impurity susceptibilities of site: "//reg(SiteName(isite))
            !
            Norb = SiteNorb(isite)
            allocate(Orbs(Norb))
            Orbs = SiteOrbs(isite,1:Norb)
            Nbp = Norb**2
            Nflavor = Norb*Nspin
            allocate(wmats(Nmats));wmats=BosonicFreqMesh(Beta,Nmats)
            !
            !Read the impurity N(tau)N(0)
            allocate(tauB(Solver%NtauB_in));tauB=0d0
            tauB = linspace(0d0,Beta,Solver%NtauB_in)
            allocate(nnt(Nflavor,Nflavor,Solver%NtauB_in));nnt=0d0
            allocate(ReadLine(Nflavor*(Nflavor+1)/2))
            file = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/nn_t.DAT"
            call inquireFile(reg(file),filexists,verb=verbose)
            unit = free_unit()
            open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
            do itau=1,Solver%NtauB_in
               ReadLine=0d0
               read(unit,*) taup,ReadLine
               if(abs(taup-tauB(itau)).gt.eps) stop "Impurity bosonic tau mesh does not coincide."
               !
               !this is cumbersome but better be safe
               ndx=1
               do ib1=1,Nflavor
                  do ib2=1,ib1
                     nnt(ib1,ib2,itau) = ReadLine(ndx)
                     if(ib1.ne.ib2)nnt(ib2,ib1,itau) = ReadLine(ndx)
                     ndx=ndx+1
                  enddo
               enddo
            enddo
            deallocate(ReadLine)
            !
            !Reshape N(tau)N(0)for easier handling
            allocate(NNitau(Norb,Norb,Nspin,Nspin,Solver%NtauB_in));NNitau=czero
            do ib1=1,Nflavor
               do ib2=1,Nflavor
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
                  do iorb=1,Norb
                     do jorb=1,Norb
                        !
                        ChiMitau = ChiMitau + NNitau(iorb,jorb,ispin,jspin,:)*(-1d0)**(ispin-jspin)/4d0
                        !
                     enddo
                  enddo
               enddo
            enddo
            call dump_BosonicField(ChiMitau,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiM_"//reg(SiteName(isite))//"_t.DAT",tauB)
            !
            !FT to get ChiM(iw)
            allocate(ChiMmats(Nmats));ChiMmats=czero
            call Bitau2mats(Beta,ChiMitau,ChiMmats,tau_uniform=.true.)
            call dump_BosonicField(ChiMmats,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiM_"//reg(SiteName(isite))//"_w.DAT",wmats)
            deallocate(ChiMitau,ChiMmats,wmats)
            !
            !
            !Charge susceptibility----------------------------------------------
            !ChiC(tau) = Sum_s1s2 <Na(tau)Nb(0)> - <Na><Nb> in the solver basis
            !here I have to use the non-symmetrized Nqmc because <Na><Nb> is not symmetrized
            call AllocateBosonicField(ChiCitau,Norb,Solver%NtauB_in,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            do iorb=1,Norb
               do jorb=1,Norb
                  !
                  ib1 = iorb + Norb*(iorb-1)
                  ib2 = jorb + Norb*(jorb-1)
                  !
                  do ispin=1,Nspin
                     do jspin=1,Nspin
                        !
                        ChiCitau%screened_local(ib1,ib2,:) = ChiCitau%screened_local(ib1,ib2,:) + (NNitau(iorb,jorb,ispin,jspin,:) &
                                                           - densityQMC(iorb,iorb,ispin,isite)*densityQMC(jorb,jorb,jspin,isite))
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
            deallocate(NNitau)
            call isReal(ChiCitau)
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
            call dump_BosonicField(ChiCitau,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiC_"//reg(SiteName(isite))//"_t.DAT",axis=tauB)
            deallocate(tauB)
            !
            !FT to get ChiC(iw)
            call AllocateBosonicField(ChiCmats,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            call Bitau2mats(Beta,ChiCitau%screened_local,ChiCmats%screened_local,tau_uniform=.true.)
            call DeallocateBosonicField(ChiCitau)
            call isReal(ChiCmats)
            !
            !Remove the iw=0 divergency of local charge susceptibility
            if(removeCDW_C)then
               write(*,"(A)") new_line("A")//"     Divergency removal in ChiC(iw=0) of site: "//reg(SiteName(isite))
               call dump_BosonicField(ChiCmats,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiC_CDW_"//reg(SiteName(isite))//"_w.DAT")
               allocate(CDW(ChiCmats%Nbp,ChiCmats%Nbp));CDW=0d0
               CDW = real(ChiCmats%screened_local(:,:,1))
               call remove_CDW(ChiCmats,"diag")
               CDW = CDW - real(ChiCmats%screened_local(:,:,1))
               call dump_Matrix(CDW,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","CDW_"//reg(SiteName(isite))//".DAT")
               deallocate(CDW)
            endif
            !
            !Print ChiC
            call dump_BosonicField(ChiCmats,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiC_"//reg(SiteName(isite))//"_w.DAT")
            !
            !recollect curlyU
            call AllocateBosonicField(curlyU,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
            call read_BosonicField(curlyU,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
            !
            !Bosonic Dyson equation in the solver basis
            write(*,"(A)") new_line("A")//"     Solving bosonic Dyson of site: "//reg(SiteName(isite))
            call AllocateBosonicField(Pimp,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            call calc_Pimp(Pimp,curlyU,ChiCmats)
            !
            !Compute convergence benchmark for the interaction
            call AllocateBosonicField(Wimp,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
            call calc_Wimp(Wimp,curlyU,ChiCmats)
            !
            !Expand to the Lattice basis
            if(RotateHloc)then
               call init_Uelements(Norb,PhysicalUelements)
               call imp2loc(C_EDMFT,ChiCmats,isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag,Map=PhysicalUelements%Full_Map)
               call imp2loc(P_EDMFT,Pimp,isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag,Map=PhysicalUelements%Full_Map)
               call imp2loc(W_EDMFT,Wimp,isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag,Map=PhysicalUelements%Full_Map)
            else
               call imp2loc(C_EDMFT,ChiCmats,isite,Orbs,ExpandImpurity,AFMselfcons)
               call imp2loc(P_EDMFT,Pimp,isite,Orbs,ExpandImpurity,AFMselfcons)
               call imp2loc(W_EDMFT,Wimp,isite,Orbs,ExpandImpurity,AFMselfcons)
            endif
            if(RotateUloc)then
               call imp2loc(curlyU_EDMFT,curlyU,isite,Orbs,ExpandImpurity,AFMselfcons,U=OlocRotDag,Map=PhysicalUelements%Full_Map)
            else
               call imp2loc(curlyU_EDMFT,curlyU,isite,Orbs,ExpandImpurity,AFMselfcons)
            endif
            !
            call DeallocateBosonicField(ChiCmats)
            call DeallocateBosonicField(curlyU)
            call DeallocateBosonicField(Pimp)
            call DeallocateBosonicField(Wimp)
            !
            deallocate(Orbs)
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         deallocate(densityQMC)
         !
         !symmetrize GW indexes
         if(EqvGWndx%O)then
            if(verbose)call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_noSym_w.DAT")
            call symmetrize_GW(P_EDMFT,EqvGWndx)
         endif
         !
         !Remove the iw=0 divergency of local polarization
         if(removeCDW_P)then
            write(*,"(A)") new_line("A")//"     Divergency removal in Pimp(iw=0)."
            call dump_BosonicField(P_EDMFT,reg(PrevItFolder),"Pimp_CDW_w.DAT")
            call remove_CDW(P_EDMFT,"diag")
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
   end subroutine collect_QMC_results


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
      if(D_correction%status)then
         write(*,"(A)") "     Deallocating D_correction"
         call DeallocateFermionicField(D_correction)
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
      integer                               :: iorb,jorb,isite
      integer                               :: wn,ws,wsi,wnmin
      integer                               :: l1,l2,l3,l4,l5,l6
      integer                               :: enlrg_
      character(len=255)                    :: header1="Lattice density"
      character(len=255)                    :: header2="Impurity density"
      character(len=255)                    :: header3="Solver density"
      character(len=255)                    :: header4="Impurity magnetization"
      character(len=255)                    :: header5="Solver magnetization"
      character(len=255)                    :: header6="Lattice magnetization"
      !
      Norb = Crystal%Norb
      !
      l1=len(trim(header1)//" up")
      l2=len(trim(header2)//" up")
      l3=len(trim(header3)//" up")
      l4=len(trim(header4))
      l5=len(trim(header5))
      l6=len(trim(header6))
      !
      wnmin=max(maxval([l1,l2,l3,l4,l5,l6]),Norb*6) !6 because I have 4 precision
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
            write(*,"(2(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header1)//" up",wn*Norb),banner(trim(header1)//" dw",wn*Norb)
            do iorb=1,Norb
               write(*,"(2("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X))") (dreal(densityGW(iorb,jorb,1)),jorb=1,Norb),(dreal(densityGW(iorb,jorb,2)),jorb=1,Norb)
            enddo
            !
            write(*,*)
            write(*,"(2(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header2)//" up",wn*Norb),banner(trim(header2)//" dw",wn*Norb)
            do iorb=1,Norb
               write(*,"(2("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X))") (dreal(densityDMFT(iorb,jorb,1)),jorb=1,Norb),(dreal(densityDMFT(iorb,jorb,2)),jorb=1,Norb)
            enddo
            !
            if(.not.EqvGWndx%S)then
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
                  Norb_imp = SiteNorb(isite)
                  wsi = wn*Norb - wn*Norb_imp
                  !
                  write(*,*)
                  write(*,"(2(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header3)//" up",wn*Norb),banner(trim(header3)//" dw",wn*Norb)
                  do iorb=1,Norb_imp
                     write(*,"(2("//str(wsi)//"X,"//str(Norb_imp)//"F"//str(wn)//".4,"//str(ws)//"X))") &
                                       (densityQMC(iorb,jorb,1,isite),jorb=1,Norb_imp),(densityQMC(iorb,jorb,2,isite),jorb=1,Norb_imp)
                  enddo
                  !
                  if(.not.EqvGWndx%S)then
                     write(*,*)
                     write(*,"(A"//str(wn*Norb)//","//str(ws)//"X)")banner(trim(header5),wn*Norb)
                     !do iorb=1,Norb_imp
                     write(*,"("//str(wsi)//"X,"//str(Norb_imp)//"F"//str(wn)//".4,"//str(ws)//"X)")((densityQMC(iorb,iorb,1,isite)-densityQMC(iorb,iorb,2,isite)),iorb=1,Norb_imp)
                     !enddo
                  endif
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
            write(*,"(2(A"//str(wn*Norb)//","//str(ws)//"X))")banner(trim(header1)//" up",wn*Norb),banner(trim(header1)//" dw",wn*Norb)
            do iorb=1,Norb
               write(*,"(2("//str(Norb)//"F"//str(wn)//".4,"//str(ws)//"X))") (dreal(densityGW(iorb,jorb,1)),jorb=1,Norb),(dreal(densityGW(iorb,jorb,2)),jorb=1,Norb)
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
               call remove_CDW(P_EDMFT,"imp")
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
            call read_Matrix(HartreeU,reg(Beta_Match%Path)//"HartreeU",paramagnet)
            call dump_Matrix(HartreeU,reg(PrevItFolder),"HartreeU",paramagnet)
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
               call remove_CDW(P_EDMFT,"imp")
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
            call read_Matrix(HartreeU,reg(Beta_Match%Path)//"HartreeU",paramagnet)
            call dump_Matrix(HartreeU,reg(PrevItFolder),"HartreeU",paramagnet)
            deallocate(HartreeU)
            allocate(Nimp(Crystal%Norb,Crystal%Norb,Nspin));Nimp=czero
            call read_Matrix(Nimp,reg(Beta_Match%Path)//"Nimp",paramagnet)
            call dump_Matrix(Nimp,reg(PrevItFolder),"Nimp",paramagnet)
            deallocate(Nimp)
            !
      end select
      !
   end subroutine interpolate_from_oldBeta


end module utils_main
