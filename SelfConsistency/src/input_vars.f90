module input_vars

   use parameters
   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface append_to_input_list
      module procedure i_append_to_input_list
      module procedure d_append_to_input_list
      module procedure l_append_to_input_list
      module procedure iv_append_to_input_list
      module procedure dv_append_to_input_list
      module procedure lv_append_to_input_list
      module procedure ch_append_to_input_list
   end interface append_to_input_list

   interface parse_Cmd_variable
      module procedure :: i_parse_cmd_variable
      module procedure :: d_parse_cmd_variable
      module procedure :: l_parse_cmd_variable
      module procedure :: iv_parse_cmd_variable
      module procedure :: dv_parse_cmd_variable
      module procedure :: lv_parse_cmd_variable
      module procedure :: ch_parse_cmd_variable
   end interface parse_Cmd_variable

   interface parse_Input_variable
      module procedure :: i_parse_input
      module procedure :: d_parse_input
      module procedure :: l_parse_input
      module procedure :: iv_parse_input
      module procedure :: dv_parse_input
      module procedure :: lv_parse_input
      module procedure :: ch_parse_input
   end interface parse_Input_variable

   type input_var
      integer,pointer                       :: i
      real(8),pointer                       :: d
      logical,pointer                       :: l
      character(len=:),pointer              :: ch
   end type input_var

   type input_node
      type(input_var),allocatable           :: var(:)
      character(len=3)                      :: type
      character(len=100)                    :: name
      character(len=512)                    :: comment
      type(input_node),pointer              :: next !link to next box
   end type input_node

   type input_list
      logical                               :: status=.false.
      integer                               :: size
      type(input_node),pointer              :: root
   end type input_list

   type input_variable
      character(len=64)                     :: name
      character(len=256)                    :: value
      character(len=20),allocatable         :: args(:)
   end type input_variable

   type input_comment
      integer                               :: Pos
      character(:),allocatable              :: descriptor
   end type input_comment

   type MpiEnv
      integer                               :: MpiComm
      integer                               :: MpiErr
      integer                               :: MpiRank
      integer                               :: MpiSize
      logical                               :: status=.false.
   end type MpiEnv

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   type(input_list),private                 :: default_list
   type(input_comment),allocatable,private  :: comments(:)
   character(len=255),private               :: p_buffer
   character(len=7),private                 :: file_status
   integer,parameter,private                :: pos_comment=45
   logical,private                          :: IOinput=.true.
   !
   character(len=256),public                :: CalculationType
   !
   !OMP parallelization
   integer,public                           :: Nthread
   !
   !K-points
   integer,public                           :: Nkpt3(3)
   logical,public                           :: UseXepsKorder
   !
   !model H(k)
   logical,public                           :: Hmodel
   logical,public                           :: readHr
   logical,public                           :: readHk
   real(8),public                           :: LatticeVec(3,3)
   real(8),allocatable,public               :: ucVec(:,:)
   integer,public                           :: Norb_model !<== remove from public variables
   real(8),allocatable,public               :: hopping(:)
   type(Heterostructures),public            :: Hetero
   !
   !Site and Orbital space
   integer,public                           :: Nsite
   logical,public                           :: ExpandImpurity
   logical,public                           :: RotateHloc
   logical,public                           :: RotateUloc
   logical,public                           :: AFMselfcons
   type(LocalOrbitals),allocatable,public   :: LocalOrbs(:)
   logical,public                           :: addCF
   !
   !Equivalent lattice indexes
   integer,public                           :: sym_mode
   type(Equivalent),public                  :: EqvGWndx
   type(Equivalent),allocatable,public      :: EqvImpndx(:)
   !
   !Imaginary time and frequency meshes
   real(8),public                           :: Beta
   integer,public                           :: Ntau
   logical,public                           :: tau_uniform
   real(8),public                           :: wmatsMax
   integer,public                           :: Nmats
   integer,public                           :: Nreal
   real(8),public                           :: wrealMax
   real(8),public                           :: eta
   real(8),public                           :: PadeWlimit
   !
   !Density lookup
   type(musearch),public                    :: look4dens
   !
   !Interaction variables
   logical,public                           :: UfullStructure
   logical,public                           :: Umodel
   logical,public                           :: Uspex
   logical,public                           :: Ustart
   logical,public                           :: U_AC
   real(8),public                           :: Uthresh
   real(8),public                           :: Uaa=0d0
   real(8),public                           :: Uab=0d0
   real(8),public                           :: J=0d0
   integer,public                           :: Nphonons
   real(8),allocatable,public               :: g_eph(:)
   real(8),allocatable,public               :: wo_eph(:)
   integer,public                           :: N_Vnn
   real(8),allocatable,public               :: Vnn(:,:,:),Vnn_diag(:,:)
   character(len=256),public                :: long_range
   !
   !Double counting types, divergencies, scaling and self-consistency coefficients
   logical,public                           :: Vxc_in
   character(len=256),public                :: SpexVersion
   character(len=256),public                :: VH_type
   character(len=256),public                :: VN_type
   character(len=256),public                :: DC_type_S
   character(len=256),public                :: DC_type_P
   character(len=256),public                :: Embedding
   logical,public                           :: addTierIII
   logical,public                           :: RecomputeG0W0
   integer,public                           :: HandleGammaPoint
   logical,public                           :: calc_Sguess
   logical,public                           :: calc_Pguess
   logical,public                           :: GoWoDC_loc
   logical,public                           :: RemoveHartree
   logical,public                           :: Dyson_Imprvd_F
   logical,public                           :: Dyson_Imprvd_B
   real(8),public                           :: alphaChi
   real(8),public                           :: alphaPi
   real(8),public                           :: alphaSigma
   real(8),public                           :: alphaHk
   logical,public                           :: Mixing_Delta_tau
   real(8),public                           :: Mixing_Delta
   real(8),public                           :: Mixing_curlyU
   integer,public                           :: Mixing_period
   logical,public                           :: causal_D
   logical,public                           :: causal_U
   logical,public                           :: recalc_Hartree
   character(len=256),public                :: causal_U_type
   !
   !Variables for the fit on Delta, Gimp, Simp
   character(len=256),public                :: DeltaFit
   integer,public                           :: Nfit
   real(8),public                           :: ReplaceTail_Simp
   !
   !Paths (directories must end with "/") and loop variables
   integer,public                           :: FirstIteration
   integer,public                           :: LastIteration
   character(len=256),public                :: pathINPUT
   character(len=256),public                :: pathINPUTtr
   character(len=256),public                :: pathDATA
   integer,public                           :: LOGfile
   logical,public                           :: dump_Gk
   logical,public                           :: dump_Sigmak
   logical,public                           :: dump_Chik
   logical,public                           :: dump_Wk
   !
   !Post-processing variables
   character(len=256),public                :: structure
   character(len=256),public                :: path_funct
   integer,private                          :: Nsym_user
   real(8),allocatable,public               :: UserPath(:,:)
   integer,public                           :: Nkpt_path
   integer,public                           :: Nmats_MaxEnt
   integer,public                           :: Ntau_MaxEnt
   integer,public                           :: Nkpt_Fermi
   logical,public                           :: FermiSurf
   real(8),public                           :: FermiCut
   real(8),public                           :: KKcutoff
   !
   !Variables for the matching beta
   type(OldBeta),public                     :: Beta_Match
   !
   !Variables for the gap equation
   type(SCDFT),public                       :: gap_equation
   !
   !Variables related to the impurity solver
   type(QMC),public                         :: Solver
   !
   !Variables not to be readed but used in utils_main
   logical,public                           :: paramagnet=.true.
   logical,public                           :: paramagneticSPEX=.true.
   logical,public                           :: XEPSisread=.false.
   logical,public                           :: solve_DMFT=.true.
   logical,public                           :: bosonicSC=.false.
   !
   !Testing flags -  default values changed only if Testing=.true. at compile time
   logical,private                          :: TESTING=.false.
   !
   logical,public                           :: Kdiag=.false.
   logical,public                           :: Ktilda=.true.
   logical,public                           :: ChiDiag=.false.
   logical,public                           :: removeCDW_C=.false.
   logical,public                           :: removeCDW_P=.false.
   logical,public                           :: Pimp_NaNb=.true.
   real(8),public                           :: dampStail=-10d0
   !
   logical,public                           :: Generic_test_flag_1=.false.
   logical,public                           :: Generic_test_flag_2=.false.
   logical,public                           :: Generic_test_flag_3=.false.
   logical,public                           :: Generic_test_flag_4=.false.
   !
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: save_InputFile
   public :: print_InputFile
   public :: delete_Input
   public :: parse_Cmd_variable
   public :: parse_Input_variable
   public :: append_to_input_list
   public :: read_InputFile
   public :: add_separator

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the Inputfile
   !---------------------------------------------------------------------------!
   subroutine read_InputFile(InputFile)
      use omp_lib
      use linalg, only : diag
      use utils_misc
      use parameters
      implicit none
      character(len=*)                      :: InputFile
      integer                               :: isite,ilayer,iset,iph,irange
      integer                               :: isym_user
      logical                               :: readVnn
      !
      !OMP parallelization.
      !call execute_command_line(" lscpu | grep 'CPU(s):       ' | awk '{print $2}' > Nthread.used ")
      !call execute_command_line(" grep '#$ -pe' submit          | awk '{print $4}' > Nthread.used ")
      !unit=free_unit()
      !open(unit,file=reg("Nthread.used"),form="formatted",action="read",status="old")
      !read(unit,*)Nthread
      !close(unit)
      !call omp_set_num_threads(Nthread) this was called in the main.
      !
      !Done in the submit file via " export OMP_NUM_THREADS= # "
      Nthread = omp_get_max_threads()
      !
      write(*,"(A,1I4)") new_line("A")//"Setting Nthread:",Nthread
      write(*,"(A)") "Reading InputFile: "//reg(InputFile)//new_line("A")
      !
      call add_separator("Calculation type")
      call parse_input_variable(CalculationType,"CALC_TYPE",InputFile,default="GW+EDMFT",comment="Calculation type. Avalibale: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT.")
      if((reg(CalculationType).eq."G0W0").or.(reg(CalculationType).eq."scGW")) solve_DMFT=.false.
      if((reg(CalculationType).eq."EDMFT").or.(reg(CalculationType).eq."GW+EDMFT")) then
         bosonicSC=.true.
         Solver%retarded=1!.true.
      endif
      if(reg(CalculationType).eq."DMFT+dynU")Solver%retarded=1!.true.
      !
      !K-points
      call add_separator("K-points and hopping")
      call parse_input_variable(Nkpt3,"NKPT3",InputFile,default=[8,8,8],comment="Number of K-points per dimension.")
      !
      !model H(k)
      call parse_input_variable(Hmodel,"H_MODEL",InputFile,default=.false.,comment="Flag to build a model non-interacting Hamiltonian.")
      if(Hmodel)then
         call parse_input_variable(Norb_model,"NORB_MODEL",InputFile,default=1,comment="Orbitals in the model non-interacting Hamiltonian (in Hr/Hk if read).")
         call parse_input_variable(readHr,"READ_HR",InputFile,default=.false.,comment="Read W90 real space Hamiltonian (Hr.DAT) from PATH_INPUT.")
         call parse_input_variable(readHk,"READ_HK",InputFile,default=.false.,comment="Read W90 K-space Hamiltonian (Hk.DAT) from PATH_INPUT.")
         call parse_input_variable(Hetero%status,"HETERO",InputFile,default=.false.,comment="Flag to build an heterostructured setup from model the non-interacting Hamiltonian.")
         if(readHr.and.readHk) stop "read_InputFile: Make up your mind, READ_HR or READ_HK?"
         call parse_input_variable(LatticeVec(:,1),"LAT_VEC_1",InputFile,default=[1d0,0d0,0d0],comment="Unit cell vector #1 of the model lattice.")
         call parse_input_variable(LatticeVec(:,2),"LAT_VEC_2",InputFile,default=[0d0,1d0,0d0],comment="Unit cell vector #2 of the model lattice.")
         call parse_input_variable(LatticeVec(:,3),"LAT_VEC_3",InputFile,default=[0d0,0d0,1d0],comment="Unit cell vector #3 of the model lattice.")
         allocate(hopping(Norb_model));hopping=0d0
         if((.not.readHr).and.(.not.readHk))call parse_input_variable(hopping,"HOPPING",InputFile,comment="NN hopping for each orbital of the non-interacting Hamiltonian.")
         !Heterostructured setup
         if(Hetero%status)then
            if(Nkpt3(3).ne.1) stop "read_InputFile: requested Heterostructured non-interacting Hamiltonian but Nk_z is not 1."
            Hetero%Norb = Norb_model
            call parse_input_variable(Hetero%Nslab,"NSLAB",InputFile,default=20,comment="Global dimension fo the slab.")
            call parse_input_variable(Hetero%Explicit,"EXPLICIT",InputFile,default=[1,10],comment="Index boundaries of the impurities explicitly solved.")
            call parse_input_variable(Hetero%tzRange,"TZ_RANGE",InputFile,default=1,comment="Range of the longitudinal hopping.")
            Hetero%Nlayer = Hetero%Explicit(2)-Hetero%Explicit(1)+1
            if(Hetero%Nlayer.le.1) stop "read_InputFile: a single layer heterostructure does not make sense."
            !setting up the indexes of the longitudinal hopping: t_1 is the hopping connecting layer #1 and #2
            Hetero%tzIndex(1) = Hetero%Explicit(1)
            Hetero%tzIndex(2) = Hetero%Explicit(2) - 1
            if(Hetero%Explicit(1).ne.1) Hetero%tzIndex(1) = Hetero%tzIndex(1) - 1              ! hopping to the left potential
            if(Hetero%Explicit(2).ne.Hetero%Nslab) Hetero%tzIndex(2) = Hetero%tzIndex(2) + 1   ! hopping to the right potential
            !allocate(Hetero%tz(Norb_model,Hetero%tzIndex(1):Hetero%tzIndex(2),Hetero%tzRange));Hetero%tz=0d0
            !do ilayer=Hetero%tzIndex(1),Hetero%tzIndex(2)
            allocate(Hetero%tz(Norb_model,Hetero%Nlayer,Hetero%tzRange));Hetero%tz=0d0
            do irange=1,Hetero%tzRange
               do ilayer=1,Hetero%Nlayer
                  !call parse_input_variable(Hetero%tz(:,ilayer,irange),"TZ"//str(irange)//"_"//str(ilayer),InputFile,comment="Longitudinal hopping for each orbital at distance "//str(irange)//" between layer #"//str(ilayer)//" and layer #"//str(ilayer+1))
                  call parse_input_variable(Hetero%tz(:,ilayer,irange),"TZ"//str(irange)//"_"//str(ilayer),InputFile,comment="Longitudinal hopping for each orbital at distance "//str(irange)//" starting from layer #"//str(ilayer))
               enddo
            enddo
         endif
      else
         call parse_input_variable(UseXepsKorder,"XEPS_KORDER",InputFile,default=.true.,comment="Flag to use the K-point ordering of XEPS.DAT if present.")
         readHk=.true.
      endif
      !
      !Site and Orbital space
      call add_separator("Sites and orbitals")
      call append_to_input_list(Nspin,"NSPIN","Number of spins. User cannot set this as its fixed to 2.")
      call parse_input_variable(Nsite,"NSITE",InputFile,default=1,comment="Number of inequivalent sites in the lattice.")
      if(Hetero%status)then
         if(Nsite.ne.(Hetero%Explicit(2)-Hetero%Explicit(1)+1)) stop "read_InputFile: The number of explicit slabs does not match with NSITE."
      endif
      call parse_input_variable(ExpandImpurity,"EXPAND",InputFile,default=.false.,comment="Flag to use a single impurity solution for all the sites of the lattice. Only indexes for site 1 readed.")
      call parse_input_variable(RotateHloc,"ROTATE_F",InputFile,default=.false.,comment="Solve the Fermionic impurity problem in the basis where H(R=0) is diagonal.")
      call parse_input_variable(RotateUloc,"ROTATE_B",InputFile,default=RotateHloc,comment="Solve the Bosonic impurity problem in the basis where H(R=0) is diagonal.")
      call parse_input_variable(AFMselfcons,"AFM",InputFile,default=.false.,comment="Flag to use the AFM self-consistency by flipping the spin. Requires input with doubled unit cell.")
      allocate(LocalOrbs(Nsite))
      do isite=1,Nsite
         if( .not.(ExpandImpurity) .or. (ExpandImpurity.and.(isite.eq.1)) )then
            call parse_input_variable(LocalOrbs(isite)%Name,"NAME_"//str(isite),InputFile,default="El",comment="Chemical species (or 2-character tag) of the site number "//str(isite))
            call parse_input_variable(LocalOrbs(isite)%Norb,"NORB_"//str(isite),InputFile,default=1,comment="Number of orbitals in site number "//str(isite))
            LocalOrbs(isite)%Nflavor = LocalOrbs(isite)%Norb*Nspin
            allocate(LocalOrbs(isite)%Orbs(LocalOrbs(isite)%Norb));LocalOrbs(isite)%Orbs=0
            call parse_input_variable(LocalOrbs(isite)%Orbs,"ORBS_"//str(isite),InputFile,default=LocalOrbs(isite)%Orbs,comment="Lattice orbital indexes of site number "//str(isite))
         elseif(ExpandImpurity.and.(isite.gt.1))then
            LocalOrbs(isite)%Name = LocalOrbs(1)%Name
            LocalOrbs(isite)%Norb = LocalOrbs(1)%Norb
            LocalOrbs(isite)%Nflavor = LocalOrbs(1)%Norb*Nspin
            if(LocalOrbs(1)%Norb.gt.1)then
               if(abs(LocalOrbs(1)%Orbs(2)-LocalOrbs(1)%Orbs(1)).eq.1)then
                  LocalOrbs(isite)%Orbs = LocalOrbs(1)%Orbs + LocalOrbs(1)%Norb*(isite-1)
               elseif(abs(LocalOrbs(1)%Orbs(2)-LocalOrbs(1)%Orbs(1)).eq.Nsite)then
                  LocalOrbs(isite)%Orbs = LocalOrbs(1)%Orbs + isite-1
               endif
            else
               LocalOrbs(isite)%Orbs = LocalOrbs(1)%Orbs + 1*(isite-1)
            endif
         endif
      enddo
      if(Hmodel)then
         allocate(ucVec(3,Nsite))
         ucVec=czero
         do isite=1,Nsite
            call parse_input_variable(ucVec(:,isite),"UC_VEC_"//str(isite),InputFile,default=[0d0,0d0,0d0],comment="Position of site #"//str(isite)//" inside the unit cell.")
         enddo
      endif
      call parse_input_variable(addCF,"ADD_CF",InputFile,default=.false.,comment="Flag to include additional crystal-fields.")
      if(addCF)then
         do isite=1,Nsite
            allocate(LocalOrbs(isite)%CrystalField(LocalOrbs(isite)%Norb));LocalOrbs(isite)%CrystalField=0d0
            call parse_input_variable(LocalOrbs(isite)%CrystalField,"CF_"//str(isite),InputFile,default=LocalOrbs(isite)%CrystalField,comment="Additional crystal-fields on diagonal orbitals of site number "//str(isite))
            if(ExpandImpurity)exit
         enddo
      endif
      !
      !Equivalent lattice indexes
      call add_separator("Symmetrizations")
      call parse_input_variable(EqvGWndx%para,"PARAMAGNET",InputFile,default=1,comment="If =1 spin symmetry is enforced if =0 spin is left free.")
      call parse_input_variable(EqvGWndx%hseed,"H_SEED",InputFile,default=0d0,comment="Seed to break spin symmetry (persistent if non zero).")
      call parse_input_variable(EqvGWndx%Nset,"EQV_SETS",InputFile,default=1,comment="Number of sets of locally equivalent lattice orbitals.")
      paramagnet = EqvGWndx%para.gt.0
      if(EqvGWndx%Nset.gt.0)then
         allocate(EqvGWndx%SetNorb(EqvGWndx%Nset))
         do iset=1,EqvGWndx%Nset
            call parse_input_variable(EqvGWndx%SetNorb(iset),"EQV_NORB_"//str(iset),InputFile,default=1,comment="Number of equivalent lattice orbitals in the set number "//str(iset))
         enddo
         allocate(EqvGWndx%SetOrbs(EqvGWndx%Nset,maxval(EqvGWndx%SetNorb)));EqvGWndx%SetOrbs=0
         do iset=1,EqvGWndx%Nset
            call parse_input_variable(EqvGWndx%SetOrbs(iset,1:EqvGWndx%SetNorb(iset)),"EQV_ORBS_"//str(iset),InputFile,comment="Lattice orbital indexes of equivalent set number "//str(iset))
         enddo
      endif
      if(EqvGWndx%Nset.gt.0)then
         call parse_input_variable(sym_mode,"SYM_MODE",InputFile,default=3,comment="If =1 only the lattice orbitals will be symmetrized, if =2 also the corresponding n(tau) inside the solver, if =3 (PREFERRED) only n(tau). 0 to avoid.")
      else
         sym_mode=0
         call append_to_input_list(sym_mode,"SYM_MODE","No orbital symmetrizations.")
      endif
      !
      !Imaginary time and frequency meshes
      call add_separator("Temperature and frequency meshes")
      call parse_input_variable(Beta,"BETA",InputFile,default=10.d0,comment="Inverse temperature in 1/eV.")
      call parse_input_variable(wmatsMax,"MAX_WMATS",InputFile,default=100.d0,comment="Maximum value of the Matsubara frequency mesh.")
      Nmats = int(Beta*wmatsMax/(2d0*pi))
      call append_to_input_list(Nmats,"NMATS","Number of points on the imaginary frequency axis. User cannot set this as its computed from MAX_WMATS and BETA.")
      call parse_input_variable(Ntau,"NTAU_LAT",InputFile,default=int(2d0*pi*Nmats),comment="Number of points on the imaginary time axis for Fermionic and Bosonic lattice fields. Its gonna be made odd.")
      if(mod(Ntau,2).eq.0)Ntau=Ntau+1
      if(mod(Ntau-1,4).ne.0)Ntau=Ntau+mod(Ntau-1,4)
      call parse_input_variable(tau_uniform,"TAU_UNIF",InputFile,default=.false.,comment="Flag to use a uniform mesh on the imaginary time axis. Only internal for GW.")
      call parse_input_variable(Nreal,"NREAL",InputFile,default=2000,comment="Number of points on the real frequency axis.")
      call parse_input_variable(wrealMax,"MAX_WREAL",InputFile,default=0d0,comment="Maximum absolute value of the real frequency axis.")
      call parse_input_variable(eta,"ETA",InputFile,default=0.04d0,comment="Real frequency broadening.")
      call parse_input_variable(PadeWlimit,"WPADE",InputFile,default=0d0,comment="Higest Matsubara frequency used in pade' analytic continuation. If its =0d0 Pade will not be performed.")
      !
      !Density lookup
      call add_separator("Density lookup")
      call parse_input_variable(look4dens%mu,"MU",InputFile,default=0d0,comment="Absolute chemical potential or shift with respect to the half-filling mu depending on REMOVE_HARTREE.")
      call parse_input_variable(look4dens%TargetDensity,"N_READ_LAT",InputFile,default=0d0,comment="Target density.")
      if(ExpandImpurity.or.AFMselfcons)then
         call parse_input_variable(look4dens%local,"N_READ_LAT_LOC",InputFile,default=.false.,comment="Flag to restrict the lattice density lookup to the ORBS_1 indexes corresponding to the solved impurity.")
         if(look4dens%local)then
            look4dens%orbs = LocalOrbs(1)%Orbs
            look4dens%TargetDensity = look4dens%TargetDensity/Nsite
         endif
      endif
      call parse_input_variable(Solver%mu_scan,"NSCAN_IMP",InputFile,default=1,comment="Integer flag to switch on density lookup within the solver.")
      call parse_input_variable(look4dens%mu_scan,"NSCAN_LAT",InputFile,default=1,comment="Integer flag to switch on density lookup within the lattice problem (self-energy not recomputed).")
      call parse_input_variable(look4dens%densityRelErr,"N_ERR",InputFile,default=0.01d0,comment="Relative error on the target density. Better if not lower than 1e-3.")
      call parse_input_variable(look4dens%muStep,"MU_STEP",InputFile,default=0.2d0,comment="Initial chemical potential step in the density lookup.")
      call parse_input_variable(look4dens%muIter,"MU_ITER",InputFile,default=50,comment="Maximum number of iterations in the density lookup.")
      call parse_input_variable(look4dens%muTime,"MU_TIME",InputFile,default=0.5d0,comment="Minutes of solver runtime in the density lookup.")
      !
      !Interaction variables
      call add_separator("Interaction")
      call parse_input_variable(UfullStructure,"U_FULL",InputFile,default=.true.,comment="Flag to use all the off-diagonal components of SPEX Ucrpa or only the physical ones.")
      call parse_input_variable(Umodel,"U_MODEL",InputFile,default=.false.,comment="Flag to build a model interaction.")
      call parse_input_variable(Uspex,"U_SPEX",InputFile,default=.true.,comment="Flag to read SPEX Ucrpa.")
      call parse_input_variable(Ustart,"U_START",InputFile,default=.true.,comment="Flag to use the local Ucrpa interaction as the effetive interaction in the 0th iteration.")
      call parse_input_variable(U_AC,"U_AC",InputFile,default=.false.,comment="Flag to force the analytic continuation on the SPEX interaction.")
      call parse_input_variable(Uthresh,"U_THRES",InputFile,default=0.001d0,comment="Lowest magnitude considered in SPEX Ucrpa bare interaction (only for local interactions).")
      if((Umodel.and.Uspex).or.((.not.Umodel).and.(.not.Uspex))) stop "read_InputFile: Make up your mind, U_MODEL or U_SPEX?"
      if(Umodel)then
         call parse_input_variable(Uaa,"UAA",InputFile,default=0d0,comment="Interaction between same orbital and opposite spin electrons (orbital independent).")
         if(Norb_model.gt.1)then
            !cahnge this if you want to have them orbital dependent
            call parse_input_variable(Uab,"UAB",InputFile,default=0d0,comment="Interaction between different orbital and opposite spin electrons (orbital independent).")
            call parse_input_variable(J,"JH",InputFile,default=0d0,comment="Hund's coupling (orbital independent).")
         endif
         !phononic model U
         call parse_input_variable(Nphonons,"N_PH",InputFile,default=0,comment="Number of custom phononic modes.")
         if(Nphonons.gt.0)then
            allocate(g_eph(Nphonons));g_eph=0d0
            allocate(wo_eph(Nphonons));wo_eph=0d0
            do iph=1,Nphonons
               call parse_input_variable(g_eph,"G_PH",InputFile,comment="Custom phononic couplings.")
               call parse_input_variable(wo_eph,"WO_PH",InputFile,comment="Custom phononic energies.")
            enddo
         endif
         !long-range model U
         call parse_input_variable(N_Vnn,"N_V",InputFile,default=0,comment="Range of the non-local interaction in real space (orbital independent).")
         if(N_Vnn.gt.0)then
            !long-range coulombian
            call parse_input_variable(long_range,"LONG_RANGE",InputFile,default="Explicit",comment="Avalibale long range interaction: Explicit(reads VNN for each N_V), Coulomb(reads first VNN, max neighbor N_V), Ewald(reads first VNN, unrestricted range).")
            allocate(Vnn(Norb_model,Norb_model,N_Vnn));Vnn=0d0
            call parse_input_variable(readVnn,"READ_LONG_RANGE",InputFile,default=.false.,comment="Flag to read the long-range interaction matrix [NORB,NORB] from file PATH_INPUT/Vnn.DAT. If False the diagonal entries are provided by the user.")
            if(readVnn)then
               call read_Vnn()
            else
               allocate(Vnn_diag(Norb_model,N_Vnn));Vnn_diag=0d0
               if(reg(long_range).eq."Explicit")then
                  do irange=1,N_Vnn
                     call parse_input_variable(Vnn_diag(1:Norb_model,irange),"VNN_"//str(irange),InputFile,comment="Magnitude of the long-range interactions for each orbital at range "//str(irange))
                     Vnn(:,:,irange) = diag(Vnn_diag(1:Norb_model,irange))
                  enddo
               else
                  call parse_input_variable(Vnn_diag(1:Norb_model,1),"VNN_1",InputFile,comment="Magnitude of the long-range interactions for each orbital between nearest neighbor  sites.")
                  Vnn(:,:,1) = diag(Vnn_diag(1:Norb_model,1))
               endif
            endif
         else
            long_range="None"
         endif
         if((Nphonons.gt.0).and.(N_Vnn.gt.0)) stop "read_InputFile: Model interaction with both phonons and non-local couplings not implemented."
         !if((Nphonons.eq.0).and.(N_Vnn.eq.0)) stop "read_InputFile: Model interaction requested buth neither phonons nor long-range couplings provided."
      endif
      !
      !Double counting types, divergencies, scaling coefficients
      call add_separator("Double counting and rescaling coeff")
      call parse_input_variable(addTierIII,"TIER_III",InputFile,default=.true.,comment="Flag to include the Tier-III contribution for ab-initio calculations.")
      if(addTierIII)then
         call parse_input_variable(SpexVersion,"SPEX_VERSION",InputFile,default="Julich",comment="Version of SPEX with which the G0W0 self-energy is computed. Available: Julich, Lund.")
         call parse_input_variable(VH_type,"VH_TYPE",InputFile,default="Ustatic_SPEX",comment="Hartree term mismatch between GoWo and scGW. Available: Ubare, Ustatic, Ubare_SPEX(V_nodiv.DAT required), Ustatic_SPEX(V_nodiv.DAT required).")
         if(Umodel)VH_type="Ustatic"
         call parse_input_variable(VN_type,"VN_TYPE",InputFile,default="Nlat",comment="Density matrix used to compute the Hartree term mismatch between GoWo and scGW. Available: Nlat, Nimp, None to set VH to zero.")
         call parse_input_variable(Vxc_in,"VXC_IN",InputFile,default=.true.,comment="Flag to include the Vxc potential inside the SigmaG0W0.")
         call parse_input_variable(RecomputeG0W0,"RECOMP_G0W0",InputFile,default=.false.,comment="Flag to recompute the G0W0 self-energy from the SPEX input.")
         if(Hmodel)RecomputeG0W0=.false.
         call parse_input_variable(GoWoDC_loc,"G0W0DC_LOC",InputFile,default=.true.,comment="Keep the local contribution of Tier-III. Automatically removed if non-causal.")
      endif
      call parse_input_variable(RemoveHartree,"REMOVE_HARTREE",InputFile,default=(.not.Hmodel),comment="Remove the Hartree term (curlyU(0)*Nimp/2) from the Impurity self-energy and perform the self-consistency only with the remaining part.")
      call parse_input_variable(Dyson_Imprvd_F,"DYSON_F_IMPRVD",InputFile,default=.false.,comment="Perform the fermionic Dyson equation using the improved estimators.") !See PRB,85,205106
      Dyson_Imprvd_B=.false.
      call append_to_input_list(Dyson_Imprvd_B,"DYSON_B_IMPRVD","Perform the bosonic Dyson equation using the improved estimators (NOT IMPLEMENTED).")
      call parse_input_variable(DC_type_S,"DC_TYPE_S",InputFile,default="GlocWloc",comment="Local GW self-energy which is replaced by the DMFT one. Avalibale: GlocWloc, Sloc.")
      call parse_input_variable(DC_type_P,"DC_TYPE_P",InputFile,default="GlocGloc",comment="Local GG polarization which is replaced by the DMFT one. Avalibale: GlocGloc, Ploc.")
      call parse_input_variable(Embedding,"ADD_EMBEDDING",InputFile,default="None",comment="Constant embedding self-energy stored in PATH_INPUT. Avalibale: loc (filename: Semb_w_s[1,2].DAT), nonloc (filename: Semb_w_k_s[1,2].DAT), None to avoid.")
      if(Hmodel.or.Umodel)addTierIII=.false.
      call parse_input_variable(HandleGammaPoint,"SMEAR_GAMMA",InputFile,default=1,comment="If >0 the dielectric function will be averaged on the SMEAR_GAMMA nearest K-points close to Gamma. If <0 the UcRPA will be rescaled like a Lorentzian in the SMEAR_GAMMA nearest K-points close to Gamma. Inactive if =0.")
      if(Umodel)HandleGammaPoint=0
      call parse_input_variable(calc_Sguess,"S_GUESS",InputFile,default=.true.,comment="Use G0W0_loc as a first guess for the DMFT self-energy.")
      call parse_input_variable(calc_Pguess,"P_GUESS",InputFile,default=.false.,comment="Use GG_loc as a first guess for the DMFT polarization. If =T it sets U_START=T.")
      if(.not.solve_DMFT)then
         calc_Sguess=.false.
         calc_Pguess=.false.
      endif
      call parse_input_variable(alphaChi,"ALPHA_CHI",InputFile,default=1d0,comment="Rescaling factor for the local charge susceptibility.")
      call parse_input_variable(alphaPi,"ALPHA_PI",InputFile,default=1d0,comment="Fraction of the EDMFT polarization substituted within the lattice one.")
      call parse_input_variable(alphaSigma,"ALPHA_SIGMA",InputFile,default=1d0,comment="Fraction of the EDMFT self-energy substituted within the lattice one.")
      call parse_input_variable(alphaHk,"ALPHA_HK",InputFile,default=1d0,comment="Rescaling of the non-interacting Hamiltonian.")
      call parse_input_variable(Mixing_Delta_tau,"MIX_D_TAU",InputFile,default=.false.,comment="Flag to mix Delta(tau). If false the mix is done with Delta(iw).")
      call parse_input_variable(Mixing_Delta,"MIX_D",InputFile,default=0.5d0,comment="Fraction of the old iteration Delta.")
      call parse_input_variable(Mixing_curlyU,"MIX_U",InputFile,default=0.5d0,comment="Fraction of the old iteration curlyU.")
      call parse_input_variable(Mixing_period,"MIX_PERIOD",InputFile,default=1,comment="Backward distance with mixing iteration. if=1 mixing with the previous iteration.")
      call parse_input_variable(causal_D,"CAUSAL_D",InputFile,default=.false.,comment="Flag to employ generalized fermionic cavity construction. Active only for GW+EDMFT calculation.")
      if(reg(CalculationType).ne."GW+EDMFT")causal_D=.false.
      call parse_input_variable(causal_U,"CAUSAL_U",InputFile,default=.false.,comment="Flag to employ generalized bosonic cavity construction. Active only for GW+EDMFT calculation.")
      if(reg(CalculationType).ne."GW+EDMFT")causal_U=.false.
      if(causal_U)then
         call parse_input_variable(causal_U_type,"CAUSAL_U_TYPE",InputFile,default="curlyU",comment="Correction mode for generalized bosonic cavity construction. Available: curlyU, Ploc.")
         if((reg(causal_U_type).eq."Ploc").and.((Nsite.gt.1).or.(maxval(LocalOrbs(:)%Norb).gt.1)))causal_U_type="curlyU"
      endif
      !
      !Variables for the fit
      call parse_input_variable(DeltaFit,"DELTA_FIT",InputFile,default="Analytic",comment="Fit to extract the local energy in GW+EDMFT calculations. Available: Analytic(best), Anaderson, Inf, Moments.")
      call parse_input_variable(Nfit,"NFIT",InputFile,default=8,comment="Number of bath levels (Analytic) or coefficient (automatic limit to NFIT=4).")
      call parse_input_variable(ReplaceTail_Simp,"WTAIL_SIMP",InputFile,default=80d0,comment="Frequency value above which the tail of Simp is replaced. If =0d0 the tail is not replaced. Only via moments (automatic limit to NFIT=4).")
      call parse_input_variable(recalc_Hartree,"RECALC_HARTREE",InputFile,default=.false.,comment="Use the lattice density to compute the Hartree term of the impurity self-energy.")
      !
      !Paths and loop variables
      call add_separator("Paths")
      call parse_input_variable(pathINPUT,"PATH_INPUT",InputFile,default="InputFiles",comment="Folder within cwd where to look for input files.")
      call parse_input_variable(pathINPUTtr,"PATH_INPUT_TR",InputFile,default="Iterations",comment="Folder within cwd where to save input files after analytic continuation to the Matsubara axis.")
      call parse_input_variable(pathDATA,"PATH_DATA",InputFile,default="Iterations",comment="Folder within cwd where to store calculation data.")
      call parse_input_variable(FirstIteration,"START_IT",InputFile,default=0,comment="First iteration. If its non zero the code will look for the last item in PATH_DATA/item and start there.")
      call parse_input_variable(LastIteration,"LAST_IT",InputFile,default=100,comment="Last iteration.")
      call parse_input_variable(LOGfile,"LOGFILE",InputFile,default=6,comment="Standard output redirection unit. Use 6 to print to terminal. Not used yet.")
      call parse_input_variable(dump_Gk,"PRINT_GK",InputFile,default=.false.,comment="Print the full k-dependent Green's function (binfmt) at each iteration (mandatory for CALC_TYPE=G0W0,scGW,GW+EDMFT).")
      if(reg(CalculationType).eq."G0W0")dump_Gk=.true.
      if(reg(CalculationType).eq."scGW")dump_Gk=.true.
      if(reg(CalculationType).eq."GW+EDMFT")dump_Gk=.true.
      call parse_input_variable(dump_Sigmak,"PRINT_SIGMAK",InputFile,default=.false.,comment="Print the full k-dependent self-energy (binfmt) at each iteration (always optional).")
      !
      !Post-processing variables
      call add_separator("Post processing")
      call parse_input_variable(structure,"STRUCTURE",InputFile,default="cubic",comment="Available structures: triangular, cubic_[2,3], fcc, bcc, hex, tetragonal, orthorhombic_[1,2], User, None to avoid.")
      if(reg(structure).eq."User")then
         call parse_input_variable(Nsym_user,"NSYM_USER",InputFile,default=4,comment="Number of high-symmetry points provided by the user.")
         allocate(UserPath(3,Nsym_user));UserPath=0d0
         do isym_user=1,Nsym_user
            call parse_input_variable(UserPath(:,isym_user),"KP_"//str(isym_user),InputFile,default=[0d0,0d0,0d0],comment="User provided high symmetry K-point #"//str(isym_user))
         enddo
      endif
      call parse_input_variable(path_funct,"PATH_FUNCT",InputFile,default="None",comment="Print interacting fields on high-symmetry points. Available fields: G=Green's function, S=self-energy, GS=both. None to avoid.")
      call parse_input_variable(Nkpt_path,"NK_PATH",InputFile,default=50,comment="Number of K-points between two hig-symmetry Kpoints.")
      call parse_input_variable(Nmats_MaxEnt,"NMATS_MAXENT",InputFile,default=Nmats,comment="If >NMATS a tail will be attached before the FT to get G(tau) when computing K-resolved spectra.")
      call parse_input_variable(Ntau_MaxEnt,"NTAU_MAXENT",InputFile,default=Ntau,comment="Number of tau-points of the interpolated Field. Ignored otherwise.")
      call parse_input_variable(FermiSurf,"FERMI_SURF",InputFile,default=.false.,comment="Flag to compute the Green's function on the planar {kx,ky} sheet. The mesh is set by NK_PATH/2. Ignored if PATH_FUNCT=None.")
      call parse_input_variable(Nkpt_Fermi,"NK_FERMI",InputFile,default=50,comment="Number of K-points in the side of the planar {kx,ky} sheet.")
      call parse_input_variable(FermiCut,"FERMI_CUT",InputFile,default=0d0,comment="Energy level at which the Fermi surface is computed. Used only in Akf_builder.")
      call parse_input_variable(KKcutoff,"KK_CUTOFF",InputFile,default=50d0,comment="Real frequency cutoff for Kramers Kronig integrals, should be twice the region of interest.")
      KKcutoff=abs(KKcutoff)
      call parse_input_variable(dump_Chik,"PRINT_CHIK",InputFile,default=.false.,comment="Print the k-dependent charge susceptibility along the K-path.")
      call parse_input_variable(dump_Wk,"PRINT_WK",InputFile,default=.false.,comment="Print the k-dependent screened interaction along the K-path.")
      !
      !Variables for the matching beta
      call add_separator("Matching beta")
      call parse_input_variable(Beta_Match%status,"MATCH_BETA",InputFile,default=.false.,comment="Interpolate to new Beta.")
      if(Beta_Match%status)then
         !
         call parse_input_variable(Beta_Match%Beta_old,"OLD_BETA",InputFile,default=Beta,comment="Beta of the calculation from which Pimp and Sigma will be fitted from.")
         call parse_input_variable(Beta_Match%Path,"OLD_BETA_PATH",InputFile,default="None",comment="Folder (within cwd) from which the old Beta Pimp and Sigma will be read from.")
         call parse_input_variable(Beta_Match%wmatsMax,"OLD_BETA_MAX_WMATS",InputFile,default=100.d0,comment="Maximum value of the Matsubara frequency mesh of the old Beta calculation.")
         call append_to_input_list(Beta_Match%Nmats_old,"OLD_BETA_NMATS","Number of points on the imaginary frequency axis of the old Beta calculation. User cannot set this as its computed from OLD_BETA_MAX_WMATS and OLD_BETA.")
         !
         Beta_Match%Beta_new = Beta
         Beta_Match%Nmats_old = int(Beta_Match%Beta_old*Beta_Match%wmatsMax/(2d0*pi))
         Beta_Match%Nmats_new = Nmats
         !
         Mixing_Delta=0d0
         Mixing_curlyU=0d0
         !
      endif
      !
      !Variables for the gap equation
      call add_separator("Gap equation")
      call parse_input_variable(gap_equation%status,"CALC_TC",InputFile,default=.false.,comment="Solve the gap equation to compute the critical Temperature.")
      if(gap_equation%status)then
         call parse_input_variable(gap_equation%Tbounds,"T_BOUNDS",InputFile,default=[0.1d0,10d0],comment="Lower and upper boundaries (Kelvin) of the temperature scan.")
         call parse_input_variable(gap_equation%Tsteps,"T_STEPS",InputFile,default=10,comment="Number of points in the temperature scan.")
         call parse_input_variable(gap_equation%loops,"LOOPS",InputFile,default=100,comment="Maximum number of iteration per each Temperature point.")
         call parse_input_variable(gap_equation%DeltaErr,"DELTA_ERR",InputFile,default=1d-5,comment="Convergence threshold on Delta.")
         call parse_input_variable(gap_equation%DeltaInit,"DELTA_INIT",InputFile,default=0.1d0,comment="Initial guess for Delta[eV].")
         call parse_input_variable(gap_equation%DeltaMix,"DELTA_MIX",InputFile,default=0.5d0,comment="Fraction of the old iteration Delta.")
         call parse_input_variable(gap_equation%HkRenorm,"HK_RENORM",InputFile,default=.true.,comment="Correct the LDA DoS with the self-energy matrix at zero frequency.")
         call parse_input_variable(gap_equation%mode_ph,"MODE_PH",InputFile,default="None",comment="Whether to include phononic Kernel. Available modes: Elk, QEspresso. None to avoid.")
         call parse_input_variable(gap_equation%mode_Zph,"MODE_ZPH",InputFile,default="symrenorm",comment="Low energy limit of Zph. Available modes: symrenorm, sym, asym.")
         call parse_input_variable(gap_equation%mode_el,"MODE_EL",InputFile,default="None",comment="Whether to include electronic Kernel. Available modes: static, static+dynamic. None to avoid.")
         if((reg(gap_equation%mode_ph).eq."None").and.(reg(gap_equation%mode_el).eq."None"))stop "read_InputFile: Tc requested but no mode is choosen."
         call parse_input_variable(gap_equation%Nkpt3_intp_Hk,"NKPT3_HK",InputFile,default=Nkpt3,comment="New interpolation grid for Hk.")
         call parse_input_variable(gap_equation%Nkpt3_intp_Wk,"NKPT3_WK",InputFile,default=Nkpt3,comment="New interpolation grid for Wk.")
         call parse_input_variable(gap_equation%Wk_cutoff,"WK_CUTOFF",InputFile,default=0.98*wmatsMax,comment="Bosonic frequency cutoff for MODE_EL=static+dynamic calculations.")
         call parse_input_variable(gap_equation%printmode_ph,"PRINT_KPH",InputFile,default="None",comment="Printing mode of the phononic Kernel. Available modes: E0, diag, surf, all. None to avoid.")
         call parse_input_variable(gap_equation%printmode_el,"PRINT_KEL",InputFile,default="None",comment="Printing mode of the electronic Kernel. Available modes: E0, 0E, diag, surf, all. None to avoid.")
      endif
      !
      !Variables related to the impurity solver
      call add_separator("Impurity solver")
      Solver%Nimp = Nsite
      if(ExpandImpurity.or.AFMselfcons) Solver%Nimp = 1
      call append_to_input_list(Solver%Nimp,"NIMP","Number of impurities solved. User cannot set this as its deduced from NSITE and EXPAND.")
      call parse_input_variable(Solver%NtauF,"NTAU_F_IMP",InputFile,default=int(2d0*pi*Nmats),comment="Number of points on the imaginary time axis for Fermionic impurity fields. Its gonna be made odd.")
      if(mod(Solver%NtauF,2).eq.0)Solver%NtauF=Solver%NtauF+1
      call parse_input_variable(Solver%NtauB,"NTAU_B_IMP",InputFile,default=int(2d0*pi*Nmats),comment="Number of points on the imaginary time axis for Bosonic impurity fields. Its gonna be made odd.")
      if(mod(Solver%NtauB,2).eq.0)Solver%NtauB=Solver%NtauB+1
      call parse_input_variable(Solver%NtauF_in,"NTAU_F_IMP_IN",InputFile,default=Solver%NtauF,comment="Number of points on the fermionic imaginary time axis used in the previous iteration.")
      call parse_input_variable(Solver%NtauB_in,"NTAU_B_IMP_IN",InputFile,default=Solver%NtauB,comment="Number of points on the bosonic imaginary time axis used in the previous iteration.")
      call parse_input_variable(Solver%NtauF_D,"NTAU_F_IMP_D",InputFile,default=int(2d0*pi*Nmats),comment="Number of points on the imaginary time axis for the hybridization function.")
      call parse_input_variable(Solver%NtauB_K,"NTAU_B_IMP_K",InputFile,default=int(2d0*pi*Nmats),comment="Number of points on the imaginary time axis for the screening function.")
      Solver%TargetDensity = look4dens%TargetDensity
      if((ExpandImpurity.or.AFMselfcons).and.(.not.look4dens%local))Solver%TargetDensity = look4dens%TargetDensity/Nsite
      call append_to_input_list(Solver%TargetDensity,"N_READ_IMP","Target density in the impurity list. User cannot set this as its the the same density on within the impurity orbitals if EXPAND=F otherwise its N_READ_LAT/NSITE.")
      call parse_input_variable(Solver%Norder,"NORDER",InputFile,default=10,comment="Maximum perturbation order measured.")
      call parse_input_variable(Solver%Gexp,"GEXPENSIVE",InputFile,default=0,comment="If =1 the impurity Green's function measurment is considered expensive (Needed at high Beta*Bandwidth).")
      call parse_input_variable(Solver%Nmeas,"NMEAS",InputFile,default=1000,comment="Sweeps where expensive measurments are not performed.")
      call parse_input_variable(Solver%Ntherm,"NTHERM",InputFile,default=100,comment="Thermalization cycles. Each cycle performs NMEAS sweeps.")
      call parse_input_variable(Solver%Nshift,"NSHIFT",InputFile,default=1,comment="Proposed segment shifts at each sweep.")
      call parse_input_variable(Solver%Nswap,"NSWAP",InputFile,default=1,comment="Proposed global spin swaps at each sweep.")
      call parse_input_variable(Solver%Imprvd_F,"IMPRVD_F",InputFile,default=0,comment="If =1 the improved estimator for the self-energy will be computed.")
      if(Dyson_Imprvd_F)Solver%Imprvd_F=1
      call parse_input_variable(Solver%Imprvd_B,"IMPRVD_B",InputFile,default=0,comment="If =1 the improved estimator for the polarization will be computed (NOT IMPLEMENTED).")
      if(Dyson_Imprvd_B)Solver%Imprvd_B=1
      call parse_input_variable(Solver%N_nnt,"N_NNT",InputFile,default=1,comment="Measurment for <n_a(tau)n_b(0)> evaluation. Updated according to CALC_TYPE. Should be either =1 or 2*NTAU_B_IMP if =0 measurment avoided.")
      if(.not.bosonicSC)Solver%N_nnt=0
      call parse_input_variable(Solver%full_ntOrbSym,"NT_FULLSYM",InputFile,default=0,comment="If =1 and orbital symmetrization inside the solver is requested it averages <n_a(tau)> for all the orbitals.")
      if(sym_mode.le.1)Solver%full_ntOrbSym=0
      call parse_input_variable(Solver%PrintTime,"PRINT_TIME",InputFile,default=10,comment="Minutes that have to pass before observables are updated and stored.")
      call parse_input_variable(Solver%binlength,"BINLENGTH",InputFile,default=4,comment="If >0 the Green's function at itau will be the average within +/-binlength.")
      call parse_input_variable(Solver%binstart,"BINSTART",InputFile,default=100,comment="Tau points skipped at the beginning and end of the Green's function average.")
      call append_to_input_list(Solver%retarded,"RETARDED","Integer flag to include the frequency dependent part of the interaction. User cannot set this as its deduced from CALC_TYPE.")
      call parse_input_variable(Solver%removeUhalf,"REMOVE_UHALF",InputFile,default=0,comment="Integer flag to use the particle-hole symmetric interaction by removing the half-filling chemical potential inside the solver.")
      !if(.not.Hmodel)Solver%removeUhalf=0
      if(ExpandImpurity)then
         allocate(Solver%Time(1));Solver%Time=0
         call parse_input_variable(Solver%Time(1),"TIME_1",InputFile,default=15,comment="Minutes of solver runtime for site number 1")
      else
         allocate(Solver%Time(Nsite));Solver%Time=0
         do isite=1,Nsite
            call parse_input_variable(Solver%Time(isite),"TIME_"//str(isite),InputFile,default=15,comment="Minutes of solver runtime for site number "//str(isite))
         enddo
      endif
      !
      if(TESTING)then
         call add_separator("Testing flags")
         call parse_input_variable(Kdiag,"K_DIAG",InputFile,default=.false.,comment="Flag to use only one J-independent screening function.")
         call parse_input_variable(Ktilda,"K_TILDA",InputFile,default=.true.,comment="Flag to U(iw)-U(0) to build the screening function.")
         call parse_input_variable(ChiDiag,"CHI_DIAG",InputFile,default=.false.,comment="Flag to remove the off-diagonal components of the local charge susceptibility.")
         call parse_input_variable(removeCDW_C,"CDW_CHI",InputFile,default=.false.,comment="Flag to remove the iw=0 divergence in the local charge susceptibility.")
         call parse_input_variable(removeCDW_P,"CDW_PI",InputFile,default=.false.,comment="Flag to remove the iw=0 divergence in the local polarization.")
         call parse_input_variable(Pimp_NaNb,"PI_NANB",InputFile,default=.true.,comment="Flag to keep the off-diagonal components o fthe local impurity polarization.")
         call parse_input_variable(dampStail,"DAMP_STAIL",InputFile,default=wmatsMax/2d0,comment="Matsubara Frequency from which powers beyond -1 are removed from ImSigma. Used only when computing K-resolved spectra, ignored if <=0d0 or >MAX_WMATS.")
         call parse_input_variable(Generic_test_flag_1,"GENERIC_TEST_FLAG_1",InputFile,default=.false.,comment="Flag to ___.")
         call parse_input_variable(Generic_test_flag_2,"GENERIC_TEST_FLAG_2",InputFile,default=.false.,comment="Flag to ___.")
         call parse_input_variable(Generic_test_flag_3,"GENERIC_TEST_FLAG_3",InputFile,default=.false.,comment="Flag to ___.")
         call parse_input_variable(Generic_test_flag_4,"GENERIC_TEST_FLAG_4",InputFile,default=.false.,comment="Flag to ___.")
      endif
      !
      !
      !
      !------------------------------------------------------------------------!
      call code_version()
      call save_InputFile(reg(InputFile))
      !
      !The last slash cannot be read for some reason
      pathDATA=trim(pathDATA)//"/"
      if(reg(pathINPUT).eq.reg(pathINPUTtr))then
         pathINPUT="../"//trim(pathINPUT)//"/"
         pathINPUTtr=reg(pathINPUT)
      else
         pathINPUT="../"//trim(pathINPUT)//"/"
         pathINPUTtr=trim(pathINPUTtr)//"/"
      endif
      Beta_Match%Path=trim(Beta_Match%Path)//"/"
      !
   end subroutine read_InputFile


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the non-local interaction matrix from file
   !---------------------------------------------------------------------------!
   subroutine read_Vnn()
      use utils_misc
      implicit none
      integer                               :: unit,idum,idist,iorb
      unit=free_unit()
      open(unit,file=pathINPUT//"Vnn.DAT",form="formatted",status="old",position="rewind",action="read")
      read(unit,*)
      read(unit,*)idum !Number of orbitals
      if(idum.ne.Norb_model) stop "read_InputFile/read_Vnn: unexpected orbital dimension from file pathINPUT/Vnn.DAT"
      do idist=1,N_Vnn
         read(unit,*)
         read(unit,*)idum  !distance index
         if(idum.ne.idist) stop "read_InputFile/read_Vnn: unexpected distance index from file pathINPUT/Vnn.DAT"
         do iorb=1,Norb_model
            read(unit,*) Vnn(iorb,1:Norb_model,idist)
         enddo
      enddo
      close(unit)
   end subroutine read_Vnn


   !---------------------------------------------------------------------------!
   !PURPOSE: Somethig that could save my ass
   !---------------------------------------------------------------------------!
   subroutine code_version()
      implicit none
      include "revision_SelfCons.inc"
      integer(4),dimension(8)               :: dummy
      integer(4)                            :: year
      integer(4)                            :: mese
      integer(4)                            :: day
      integer(4)                            :: h
      integer(4)                            :: m
      integer(4)                            :: s
      integer(4)                            :: ms
      character(len=9),parameter,dimension(12) :: month = (/   &
           'January  ', 'February ', 'March    ', 'April    ', &
           'May      ', 'June     ', 'July     ', 'August   ', &
           'September', 'October  ', 'November ', 'December ' /)
      !
      write(*,"(A)")("CODE VERSION: "//trim(adjustl(trim(revision))))
      call date_and_time(values=dummy)
      year = dummy(1)
      mese = dummy(2)
      day  = dummy(3)
      h    = dummy(5)
      m    = dummy(6)
      s    = dummy(7)
      ms   = dummy(8)
      write(*,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)") "Timestamp: ",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
      write(*,*)""
      if(verbose)then
         open(10,file="code_version.inc")
         write(10,"(A)")"CODE VERSION: "//trim(adjustl(trim(revision)))
         write(10,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")"Timestamp: ",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
         write(10,*)""
         close(10)
      endif
   end subroutine code_version


   !---------------------------------------------------------------------------!
   !PURPOSE: Field Separator
   !---------------------------------------------------------------------------!
   subroutine add_separator(descriptor)
      implicit none
      character(len=*)                      :: descriptor
      integer                               :: sizeOld,sizeNew
      type(input_comment),allocatable       :: prv(:)
      !
      if(allocated(comments))then
         !
         sizeOld=size(comments)
         sizeNew=sizeOld+1
         allocate(prv(sizeOld))
         prv=comments
         deallocate(comments)
         allocate(comments(sizeNew))
         !
         comments(1:sizeOld) = prv
         deallocate(prv)
         !
         comments(sizeNew)%Pos = default_list%size
         comments(sizeNew)%descriptor = trim(descriptor)
         !
      else
         !
         allocate(comments(1))
         comments(1)%Pos = 1
         comments(1)%descriptor = trim(descriptor)
         !
      endif
   end subroutine add_separator


   !---------------------------------------------------------------------------!
   !PURPOSE: Save to file the used inputfile
   !---------------------------------------------------------------------------!
   subroutine save_InputFile(file)
      implicit none
      character(len=*) :: file
      if(IOinput)then
         call print_InputFile(trim(file))
      else
         call print_InputFile(trim(file))
         write(*,*)"input file can not be found. Dumped a *default* version in: used."//trim(file)
         stop
      endif
   end subroutine save_InputFile


   !---------------------------------------------------------------------------!
   !PURPOSE: Print the used inputfile
   !---------------------------------------------------------------------------!
   subroutine print_InputFile(file,list)
      implicit none
      character(len=*),optional             :: file
      type(input_list),optional             :: list
      integer                               :: counter,counterCom,size
      type(input_node),pointer              :: c
      if(present(list))then
         c => list%root%next
      else
         c => default_list%root%next
      endif
      counter = 0
      file_status='replace'
      size=default_list%size
      if(present(list))size=list%size
      !
      counterCom = 1
      if(present(file))then
         call print_separator(comments(1)%descriptor,file)
      else
         call print_separator(comments(counterCom)%descriptor)
      endif
      !
      if(size>0)then
         do
            if(.not.associated(c))exit
            counter=counter+1
            if(present(file))then
               call print_input_node(c,file)
               if(any(comments(:)%Pos.eq.counter))then
                  counterCom=counterCom+1
                  call print_separator(comments(counterCom)%descriptor,file)
               endif
            else
               call print_input_node(c)
               if(any(comments(:)%Pos.eq.counter))then
                  counterCom=counterCom+1
                  call print_separator(comments(counterCom)%descriptor)
               endif
            endif
            c => c%next
         enddo
      else
         write(*,*)"input list: empty"
         return
      endif
      file_status='replace'
      c => null()
   end subroutine print_InputFile


   !---------------------------------------------------------------------------!
   !PURPOSE: delete the list
   !---------------------------------------------------------------------------!
   subroutine delete_Input(list)
      implicit none
      type(input_list),optional             :: list
      type(input_node),pointer              :: p,c
      integer :: i
      if(present(list))then
         do
            p => list%root
            c => p%next
            if(.not.associated(c))exit  !empty list
            p%next => c%next !
            c%next=>null()
            do i=1,size(c%var)
               nullify(c%var(i)%i)
               nullify(c%var(i)%d)
               nullify(c%var(i)%l)
               nullify(c%var(i)%ch)
            enddo
            deallocate(c%var)
            deallocate(c)
         end do
         list%status=.false.
         deallocate(list%root)
      else
         do
            p => default_list%root
            c => p%next
            if(.not.associated(c))exit  !empty list
            p%next => c%next !
            c%next=>null()
            do i=1,size(c%var)
               nullify(c%var(i)%i)
               nullify(c%var(i)%d)
               nullify(c%var(i)%l)
               nullify(c%var(i)%ch)
            enddo
            deallocate(c%var)
            deallocate(c)
         enddo
         default_list%status=.false.
         deallocate(default_list%root)
      endif
      p => null()
      c => null()
   end subroutine delete_Input


   !---------------------------------------------------------------------------!
   !PURPOSE: parse input varibales from file/CMD Line
   !---------------------------------------------------------------------------!
   !
   !SCALAR
   subroutine i_parse_input(variable,name,file,default,comment)
     use utils_misc
     implicit none
     integer                                :: variable
     integer,optional                       :: default
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     character(len=*)                       :: file
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: unit,pos
     integer                                :: status
     logical                                :: bool
     character(len=255)                     :: buffer
     if(present(default))variable=default
     name_=name;call upper_case(name_)
     inquire(file=file,exist=bool)
     if(.not.bool)IOinput=.false.
     if(bool)then
        unit=free_unit()
        open(unit,file=file)
        status=0
        var_search: do while(status>=0)
           read(unit,"(A255)",iostat=status)buffer
           pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
           var = scan_input_variable(trim(buffer))
           if(var%name==name_)then
              read(var%value,*)variable
              exit var_search
           endif
        enddo var_search
        close(unit)
     endif
     call parse_cmd_variable(variable,name_)
     if(present(comment))then
        call append_to_input_list(variable,name_,comment)
     else
        call append_to_input_list(variable,name_)
     endif
     return
   end subroutine i_parse_input
   !
   subroutine d_parse_input(variable,name,file,default,comment)
     use utils_misc
     implicit none
     real(8)                                :: variable
     real(8),optional                       :: default
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     character(len=*)                       :: file
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: unit,pos
     integer                                :: status
     logical                                :: bool
     character(len=255)                     :: buffer
     if(present(default))variable=default
     name_=name;call upper_case(name_)
     inquire(file=file,exist=bool)
     if(.not.bool)IOinput=.false.
     if(bool)then
        unit=free_unit()
        open(unit,file=file)
        status=0
        var_search: do while(status>=0)
           read(unit,"(A255)",iostat=status)buffer
           pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
           var = scan_input_variable(trim(buffer))
           if(var%name==name_)then
              read(var%value,*)variable
              exit var_search
           endif
        enddo var_search
        close(unit)
     endif
     call parse_cmd_variable(variable,name_)
     if(present(comment))then
        call append_to_input_list(variable,name_,comment)
     else
        call append_to_input_list(variable,name_)
     endif
     return
   end subroutine d_parse_input
   !
   subroutine l_parse_input(variable,name,file,default,comment)
     use utils_misc
     implicit none
     logical                                :: variable
     logical,optional                       :: default
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     character(len=*)                       :: file
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: unit,pos
     integer                                :: status
     logical                                :: bool
     character(len=255)                     :: buffer
     if(present(default))variable=default
     name_=name;call upper_case(name_)
     inquire(file=file,exist=bool)
     if(.not.bool)IOinput=.false.
     if(bool)then
        unit=free_unit()
        open(unit,file=file)
        status=0
        var_search: do while(status>=0)
           read(unit,"(A255)",iostat=status)buffer
           pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
           var = scan_input_variable(trim(buffer))
           if(var%name==name_)then
              read(var%value,*)variable
              exit var_search
           endif
        enddo var_search
        close(unit)
     endif
     call parse_cmd_variable(variable,name_)
     if(present(comment))then
        call append_to_input_list(variable,name_,comment)
     else
        call append_to_input_list(variable,name_)
     endif
     return
   end subroutine l_parse_input
   !
   !
   !VECTOR
   subroutine iv_parse_input(variable,name,file,default,comment)
     use utils_misc
     implicit none
     integer,dimension(:)                   :: variable
     integer,dimension(size(variable)),optional :: default
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     character(len=*)                       :: file
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: unit,pos,j,ndim,nargs,pos0,iarg
     integer                                :: status
     logical                                :: bool
     character(len=255)                     :: buffer
     If(present(default))variable=default
     ndim=size(variable)
     name_=name;call upper_case(name_)
     inquire(file=file,exist=bool)
     if(.not.bool)IOinput=.false.
     if(bool)then
        unit=free_unit()
        open(unit,file=file)
        status=0
        var_search: do while(status>=0)
           read(unit,"(A255)",iostat=status)buffer
           pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
           var = scan_input_variable(trim(buffer))
           if(var%name==name_)then
              nargs=check_cmd_vector_size(ndim,var)
              allocate(var%args(nargs))
              iarg=0
              pos0=0
              do j=1,len(var%value)
                 if(var%value(j:j)==",")then
                    iarg=iarg+1
                    var%args(iarg)=var%value(pos0+1:j-1)
                    pos0=j
                 endif
              enddo
              var%args(nargs)=var%value(pos0+1:)
              do iarg=1,nargs
                 read(var%args(iarg),*)variable(iarg)
              enddo
              exit var_search
           endif
        enddo var_search
        close(unit)
     endif
     call parse_cmd_variable(variable,name_)
     if(present(comment))then
        call append_to_input_list(variable,name_,comment)
     else
        call append_to_input_list(variable,name_)
     endif
     return
   end subroutine iv_parse_input
   !
   subroutine dv_parse_input(variable,name,file,default,comment)
     use utils_misc
     implicit none
     real(8),dimension(:)                   :: variable
     real(8),dimension(size(variable)),optional :: default
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     character(len=*)                       :: file
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: unit,pos,j,ndim,nargs,pos0,iarg
     integer                                :: status
     logical                                :: bool
     character(len=255)                     :: buffer
     If(present(default))variable=default
     ndim=size(variable)
     name_=name;call upper_case(name_)
     inquire(file=file,exist=bool)
     if(.not.bool)IOinput=.false.
     if(bool)then
        unit=free_unit()
        open(unit,file=file)
        status=0
        var_search: do while(status>=0)
           read(unit,"(A255)",iostat=status)buffer
           pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
           var = scan_input_variable(trim(buffer))
           if(var%name==name_)then
              nargs=check_cmd_vector_size(ndim,var)
              allocate(var%args(nargs))
              iarg=0
              pos0=0
              do j=1,len(var%value)
                 if(var%value(j:j)==",")then
                    iarg=iarg+1
                    var%args(iarg)=var%value(pos0+1:j-1)
                    pos0=j
                 endif
              enddo
              var%args(nargs)=var%value(pos0+1:)
              do iarg=1,nargs
                 read(var%args(iarg),*)variable(iarg)
              enddo
              exit var_search
           endif
        enddo var_search
        close(unit)
     endif
     call parse_cmd_variable(variable,name_)
     if(present(comment))then
        call append_to_input_list(variable,name_,comment)
     else
        call append_to_input_list(variable,name_)
     endif
     return
   end subroutine dv_parse_input
   !
   subroutine lv_parse_input(variable,name,file,default,comment)
     use utils_misc
     implicit none
     logical,dimension(:)                   :: variable
     logical,dimension(size(variable)),optional :: default
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     character(len=*)                       :: file
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: unit,pos,j,ndim,nargs,pos0,iarg
     integer                                :: status
     logical                                :: bool
     character(len=255)                     :: buffer
     If(present(default))variable=default
     ndim=size(variable)
     name_=name;call upper_case(name_)
     inquire(file=file,exist=bool)
     if(.not.bool)IOinput=.false.
     if(bool)then
        unit=free_unit()
        open(unit,file=file)
        status=0
        var_search: do while(status>=0)
           read(unit,"(A255)",iostat=status)buffer
           pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
           var = scan_input_variable(trim(buffer))
           if(var%name==name_)then
              nargs=check_cmd_vector_size(ndim,var)
              allocate(var%args(nargs))
              iarg=0
              pos0=0
              do j=1,len(var%value)
                 if(var%value(j:j)==",")then
                    iarg=iarg+1
                    var%args(iarg)=var%value(pos0+1:j-1)
                    pos0=j
                 endif
              enddo
              var%args(nargs)=var%value(pos0+1:)
              do iarg=1,nargs
                 read(var%args(iarg),*)variable(iarg)
              enddo
              exit var_search
           endif
        enddo var_search
        close(unit)
     endif
     call parse_cmd_variable(variable,name_)
     if(present(comment))then
        call append_to_input_list(variable,name_,comment)
     else
        call append_to_input_list(variable,name_)
     endif
     return
   end subroutine lv_parse_input
   !
   !
   !STRING
   subroutine ch_parse_input(variable,name,file,default,comment)
     use utils_misc
     implicit none
     character(len=*)                       :: variable
     character(len=*),optional              :: default
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     character(len=*)                       :: file
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: unit,pos
     integer                                :: status
     logical                                :: bool
     character(len=255)                     :: buffer
     if(present(default))variable=default
     name_=name;call upper_case(name_)
     inquire(file=file,exist=bool)
     if(.not.bool)IOinput=.false.
     if(bool)then
        unit=free_unit()
        open(unit,file=file)
        status=0
        var_search: do while(status>=0)
           read(unit,"(A255)",iostat=status)buffer
           pos=scan_comment(buffer);if(pos/=0)buffer=buffer(1:pos-1)
           var = scan_input_variable(trim(buffer))
           if(var%name==name_)then
              read(var%value,*)variable
              exit var_search
           endif
        enddo var_search
        close(unit)
     endif
     call parse_cmd_variable(variable,name_)
     if(present(comment))then
        call append_to_input_list(variable,name_,comment)
     else
        call append_to_input_list(variable,name_)
     endif
     return
   end subroutine ch_parse_input


   !---------------------------------------------------------------------------!
   !PURPOSE: Add variables to the inputFile or parse from CMD line
   !---------------------------------------------------------------------------!
   !
   !SCALAR
   subroutine i_parse_cmd_variable(variable,name,default)
     implicit none
     integer                                :: variable
     integer,optional                       :: default
     character(len=*)                       :: name
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: i
     If(present(default))variable=default
     name_=name;call upper_case(name_)
     do i=1,command_argument_count()
        var = scan_cmd_variable(i)
        if(var%name==name_)then
           read(var%value,*)variable
           write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
        endif
     enddo
   end subroutine i_parse_cmd_variable
   !
   subroutine d_parse_cmd_variable(variable,name,default)
     implicit none
     real(8)                                :: variable
     real(8),optional                       :: default
     character(len=*)                       :: name
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: i
     if(present(default))variable=default
     name_=name;call upper_case(name_)
     do i=1,command_argument_count()
        var = scan_cmd_variable(i)
        if(var%name==name_)then
           read(var%value,*)variable
           write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
        endif
     enddo
   end subroutine d_parse_cmd_variable
   !
   subroutine l_parse_cmd_variable(variable,name,default)
     implicit none
     logical                                :: variable
     logical,optional                       :: default
     character(len=*)                       :: name
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: i
     if(present(default))variable=default
     name_=name;call upper_case(name_)
     do i=1,command_argument_count()
        var = scan_cmd_variable(i)
        if(var%name==name_)then
           read(var%value,*)variable
           write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
        endif
     enddo
   end subroutine l_parse_cmd_variable
   !
   !
   !VECTOR
   subroutine iv_parse_cmd_variable(variable,name,default)
     implicit none
     integer,dimension(:)                   :: variable
     integer,dimension(size(variable)),optional :: default
     character(len=*)                       :: name
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: i,j,ndim,nargs,pos0,iarg
     If(present(default))variable=default
     ndim=size(variable)
     name_=name;call upper_case(name_)
     do i=1,command_argument_count()
        var = scan_cmd_variable(i)
        if(var%name==name_)then
           nargs=check_cmd_vector_size(ndim,var)
           allocate(var%args(nargs))
           iarg=0
           pos0=0
           do j=1,len(var%value)
              if(var%value(j:j)==",")then
                 iarg=iarg+1
                 var%args(iarg)=var%value(pos0+1:j-1)
                 pos0=j
              endif
           enddo
           var%args(nargs)=var%value(pos0+1:)
           do iarg=1,nargs
              read(var%args(iarg),*)variable(iarg)
           enddo
           write(0,"(A,100I6)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
        endif
     enddo
   end subroutine iv_parse_cmd_variable
   !
   subroutine dv_parse_cmd_variable(variable,name,default)
     implicit none
     real(8),dimension(:)                   :: variable
     real(8),dimension(size(variable)),optional :: default
     character(len=*)                       :: name
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: i,j,ndim,nargs,pos0,iarg
     If(present(default))variable=default
     ndim=size(variable)
     name_=name;call upper_case(name_)
     do i=1,command_argument_count()
        var = scan_cmd_variable(i)
        if(var%name==name_)then
           nargs=check_cmd_vector_size(ndim,var)
           allocate(var%args(nargs))
           iarg=0
           pos0=0
           do j=1,len(var%value)
              if(var%value(j:j)==",")then
                 iarg=iarg+1
                 var%args(iarg)=var%value(pos0+1:j-1)
                 pos0=j
              endif
           enddo
           var%args(nargs)=var%value(pos0+1:)
           do iarg=1,nargs
              read(var%args(iarg),*)variable(iarg)
           enddo
           write(0,"(A,100F18.9)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
        endif
     enddo
   end subroutine dv_parse_cmd_variable
   !
   subroutine lv_parse_cmd_variable(variable,name,default)
     implicit none
     logical,dimension(:)                   :: variable
     logical,dimension(size(variable)),optional :: default
     character(len=*)                       :: name
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: i,j,ndim,nargs,pos0,iarg
     If(present(default))variable=default
     ndim=size(variable)
     name_=name;call upper_case(name_)
     do i=1,command_argument_count()
        var = scan_cmd_variable(i)
        if(var%name==name_)then
           nargs=check_cmd_vector_size(ndim,var)
           allocate(var%args(nargs))
           iarg=0
           pos0=0
           do j=1,len(var%value)
              if(var%value(j:j)==",")then
                 iarg=iarg+1
                 var%args(iarg)=var%value(pos0+1:j-1)
                 pos0=j
              endif
           enddo
           var%args(nargs)=var%value(pos0+1:)
           do iarg=1,nargs
              read(var%args(iarg),*)variable(iarg)
           enddo
           write(0,"(A,100L3)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
        endif
     enddo
   end subroutine lv_parse_cmd_variable
   !
   !
   !STRING
   subroutine ch_parse_cmd_variable(variable,name,default)
     implicit none
     character(len=*)                       :: variable
     character(len=*),optional              :: default
     character(len=*)                       :: name
     character(len=len(name))               :: name_
     type(input_variable)                   :: var
     integer                                :: i
     if(present(default))variable=default
     name_=name;call upper_case(name_)
     do i=1,command_argument_count()
        var = scan_cmd_variable(i)
        if(var%name==name_)then
           read(var%value,*)variable
           write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
        endif
     enddo
   end subroutine ch_parse_cmd_variable


   !---------------------------------------------------------------------------!
   !PURPOSE: init the input list
   !---------------------------------------------------------------------------!
   subroutine init_input_list(list)
     implicit none
     type(input_list),optional              :: list
     if(present(list))then
        allocate(list%root)
        list%size=0
        list%status=.true.
        list%root%next=>null()
     else
        allocate(default_list%root)
        default_list%size=0
        default_list%status=.true.
        default_list%root%next=>null()
     endif
   end subroutine init_input_list


   !---------------------------------------------------------------------------!
   !PURPOSE: get list size
   !---------------------------------------------------------------------------!
   function size_input_list(list) result(size)
     implicit none
     type(input_list),optional              :: list
     integer                                :: size
     size=default_list%size
     if(present(list))size=list%size
   end function size_input_list


   !---------------------------------------------------------------------------!
   !PURPOSE: !Append input data to the list:
   !---------------------------------------------------------------------------!
   !
   !SCALAR
   subroutine i_append_to_input_list(variable,name,comment)
     implicit none
     integer,target                         :: variable
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     type(input_node),pointer               :: p,c
     if(.not.default_list%status)call init_input_list()
     p => default_list%root
     c => p%next
     do                            !traverse the list until obj < value (ordered list)
        if(.not.associated(c))exit !empty list or beginning of the list
        p => c
        c => c%next
     end do
     allocate(p%next)                !Create a new element in the list
     !
     allocate(p%next%var(1))
     p%next%var(1)%i=>variable
     p%next%name= name
     p%next%type='i'
     p%next%comment=""
     if(present(comment))p%next%comment=trim(comment)
     !
     default_list%size=default_list%size+1
     if(.not.associated(c))then !end of the list special case (current=>current%next)
        p%next%next  => null()
     else
        p%next%next  => c      !the %next of the new node come to current
     end if
     p=>null()
     c=>null()
   end subroutine i_append_to_input_list
   !
   subroutine d_append_to_input_list(variable,name,comment)
     implicit none
     real(8),target                         :: variable
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     type(input_node),pointer               :: p,c
     if(.not.default_list%status)call init_input_list()
     p => default_list%root
     c => p%next
     do                            !traverse the list until obj < value (ordered list)
        if(.not.associated(c))exit !empty list or beginning of the list
        p => c
        c => c%next
     end do
     allocate(p%next)                !Create a new element in the list
     !
     allocate(p%next%var(1))
     p%next%var(1)%d=>variable
     !
     p%next%name= name
     p%next%type='d'
     p%next%comment=""
     if(present(comment))p%next%comment=trim(comment)
     !
     default_list%size=default_list%size+1
     if(.not.associated(c))then !end of the list special case (current=>current%next)
        p%next%next  => null()
     else
        p%next%next  => c      !the %next of the new node come to current
     end if
     p=>null()
     c=>null()
   end subroutine d_append_to_input_list
   !
   subroutine l_append_to_input_list(variable,name,comment)
     implicit none
     logical,target                         :: variable
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     type(input_node),pointer               :: p,c
     if(.not.default_list%status)call init_input_list()
     p => default_list%root
     c => p%next
     do                            !traverse the list until obj < value (ordered list)
        if(.not.associated(c))exit !empty list or beginning of the list
        p => c
        c => c%next
     end do
     allocate(p%next)                !Create a new element in the list
     !
     ! allocate(p%next%l(1))
     ! p%next%l(1) = variable
     !>NEW
     allocate(p%next%var(1))
     p%next%var(1)%l=>variable
     !<
     p%next%name= name
     p%next%type='l'
     p%next%comment=""
     if(present(comment))p%next%comment=trim(comment)
     !
     default_list%size=default_list%size+1
     if(.not.associated(c))then !end of the list special case (current=>current%next)
        p%next%next  => null()
     else
        p%next%next  => c      !the %next of the new node come to current
     end if
     p=>null()
     c=>null()
   end subroutine l_append_to_input_list
   !
   !
   !VECTOR
   subroutine iv_append_to_input_list(variable,name,comment)
     implicit none
     integer,dimension(:),target            :: variable
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     type(input_node),pointer               :: p,c
     integer :: i
     if(.not.default_list%status)call init_input_list()
     p => default_list%root
     c => p%next
     do                            !traverse the list until obj < value (ordered list)
        if(.not.associated(c))exit !empty list or beginning of the list
        p => c
        c => c%next
     end do
     allocate(p%next)                !Create a new element in the list
     !
     ! allocate(p%next%i(size(variable)))
     ! p%next%i   = variable
     !>NEW
     allocate(p%next%var(size(variable)))
     do i=1,size(variable)
        p%next%var(i)%i=>variable(i)
     enddo
     !<
     p%next%name= name
     p%next%type='i'
     p%next%comment=""
     if(present(comment))p%next%comment=trim(comment)
     !
     default_list%size=default_list%size+1
     if(.not.associated(c))then !end of the list special case (current=>current%next)
        p%next%next  => null()
     else
        p%next%next  => c      !the %next of the new node come to current
     end if
     p=>null()
     c=>null()
   end subroutine iv_append_to_input_list
   !
   subroutine dv_append_to_input_list(variable,name,comment)
     implicit none
     real(8),dimension(:),target            :: variable
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     type(input_node),pointer               :: p,c
     integer :: i
     if(.not.default_list%status)call init_input_list()
     p => default_list%root
     c => p%next
     do                            !traverse the list until obj < value (ordered list)
        if(.not.associated(c))exit !empty list or beginning of the list
        p => c
        c => c%next
     end do
     allocate(p%next)                !Create a new element in the list
     !
     ! allocate(p%next%d(size(variable)))
     ! p%next%d   = variable
     !>NEW
     allocate(p%next%var(size(variable)))
     do i=1,size(variable)
        p%next%var(i)%d=>variable(i)
     enddo
     !<
     p%next%name= name
     p%next%type='d'
     p%next%comment=""
     if(present(comment))p%next%comment=trim(comment)
     !
     default_list%size=default_list%size+1
     if(.not.associated(c))then !end of the list special case (current=>current%next)
        p%next%next  => null()
     else
        p%next%next  => c      !the %next of the new node come to current
     end if
     p=>null()
     c=>null()
   end subroutine dv_append_to_input_list
   !
   subroutine lv_append_to_input_list(variable,name,comment)
     implicit none
     logical,dimension(:),target            :: variable
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     type(input_node),pointer               :: p,c
     integer :: i
     if(.not.default_list%status)call init_input_list()
     p => default_list%root
     c => p%next
     do                            !traverse the list until obj < value (ordered list)
        if(.not.associated(c))exit !empty list or beginning of the list
        p => c
        c => c%next
     end do
     allocate(p%next)                !Create a new element in the list
     !
     ! allocate(p%next%l(size(variable)))
     ! p%next%l   = variable
     !>NEW
     allocate(p%next%var(size(variable)))
     do i=1,size(variable)
        p%next%var(i)%l=>variable(i)
     enddo
     !<
     p%next%name= name
     p%next%type='l'
     p%next%comment=""
     if(present(comment))p%next%comment=trim(comment)
     !
     default_list%size=default_list%size+1
     if(.not.associated(c))then !end of the list special case (current=>current%next)
        p%next%next  => null()
     else
        p%next%next  => c      !the %next of the new node come to current
     end if
     p=>null()
     c=>null()
   end subroutine lv_append_to_input_list
   !
   !
   !STRING
   subroutine ch_append_to_input_list(variable,name,comment)
     implicit none
     character(len=*),target                :: variable
     character(len=*)                       :: name
     character(len=*),optional              :: comment
     type(input_node),pointer               :: p,c
     if(.not.default_list%status)call init_input_list()
     p => default_list%root
     c => p%next
     do                            !traverse the list until obj < value (ordered list)
        if(.not.associated(c))exit !empty list or beginning of the list
        p => c
        c => c%next
     end do
     allocate(p%next)                !Create a new element in the list
     !
     ! allocate(p%next%ch(1))
     ! p%next%ch(1) = variable
     !>NEW
     allocate(p%next%var(1))
     nullify(p%next%var(1)%ch)
     p%next%var(1)%ch=> variable
     !<
     p%next%name= name
     p%next%type='ch'
     p%next%comment=""
     if(present(comment))p%next%comment=trim(comment)
     !
     default_list%size=default_list%size+1
     if(.not.associated(c))then !end of the list special case (current=>current%next)
        p%next%next  => null()
     else
        p%next%next  => c      !the %next of the new node come to current
     end if
     p=>null()
     c=>null()
   end subroutine ch_append_to_input_list


   !---------------------------------------------------------------------------!
   !PURPOSE: auxiliary to print_InputFile
   !---------------------------------------------------------------------------!
   subroutine print_input_node(c,file)
     use utils_misc
     implicit none
     type(input_node)                       :: c
     character(len=*),optional              :: file
     character(len=255)                     :: blank=""
     integer                                :: clen
     integer                                :: unit,i
     !
     call s_blank_delete(c%name)
     select case(c%type)
     case('ch')
        p_buffer=trim(c%name)//"="//trim(adjustl(trim(c%var(1)%ch)))

     case('i')
        if(size(c%var)==1)then   !scalar
           p_buffer=trim(c%name)//"="//str(c%var(1)%i)
        else                     !vector
           p_buffer=trim(c%name)//"="
           do i=1,size(c%var)-1
              p_buffer=trim(p_buffer)//trim(str(c%var(i)%i))//","
           end do
           p_buffer=trim(p_buffer)//trim(str(c%var(size(c%var))%i))
        endif

     case('d')
        if(size(c%var)==1)then   !scalar
           p_buffer=trim(c%name)//"="//str(c%var(1)%d)
        else                     !vector
           p_buffer=trim(c%name)//"="
           do i=1,size(c%var)-1
              p_buffer=trim(p_buffer)//trim(str(c%var(i)%d))//","
           end do
           p_buffer=trim(p_buffer)//trim(str(c%var(size(c%var))%d))
        endif

     case('l')
        if(size(c%var)==1)then   !scalar
           p_buffer=trim(c%name)//"="//str(c%var(1)%l)
        else                     !vector
           p_buffer=trim(c%name)//"="
           do i=1,size(c%var)-1
              p_buffer=trim(p_buffer)//trim(str(c%var(i)%l))//","
           end do
           p_buffer=trim(p_buffer)//trim(str(c%var(size(c%var))%l))
        endif
     end select
     !
     call s_blank_delete(p_buffer)
     clen=pos_comment-len(trim(p_buffer))
     if(clen<=0)clen=1
     p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(c%comment)
     !
     ! write(*,"(1x,A)")trim(p_buffer)
     if(present(file))then
        unit=free_unit()
        open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
        write(unit,"(1x,A)")trim(p_buffer)
        close(unit)
     else
        write(unit,"(1x,A)")trim(p_buffer)
     endif
     p_buffer=""
   end subroutine print_input_node
   !
   subroutine print_separator(descriptor,file)
     use utils_misc
     implicit none
     character(len=*)                       :: descriptor
     character(len=*),optional              :: file
     integer                                :: unit,i
     character(len=pos_comment+2)           :: line
     line="# "//descriptor//" -"
     do i=1,pos_comment-1-len(descriptor)-2-1
        line = trim(line)//"-"
     enddo
     line = trim(line)//"#"
     if(present(file))then
        unit=free_unit()
        open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
     else
        unit=6
     endif
     write(unit,"(1x,A)")
     write(unit,"(1x,A)")trim(line)
     write(unit,"(1x,A)")
     if(present(file))close(unit)
  end subroutine print_separator


   !---------------------------------------------------------------------------!
   !PURPOSE: AUXILIARY ROUTINES
   !---------------------------------------------------------------------------!
   function scan_comment(buffer) result(pos)
     implicit none
     character(len=255),intent(in)           :: buffer
     integer                                 :: pos
     character(len=1),dimension(4),parameter :: comments=["!","c","#","%"]
     integer :: i
     pos=0
     do i=1,4
        pos=scan(buffer,comments(i))
        if(pos/=0)return
     end do
   end function scan_comment

   function scan_cmd_variable(i)  result(var)
     implicit none
     integer                                :: i,ncount,pos
     type(input_variable)                   :: var
     character(len=512)                     :: buffer
     ncount=command_argument_count()
     if(i>ncount)then
        write(*,*)"get_cmd_variable: i > cmdline_argument_count!"
        return
     endif
     call get_command_argument(i,buffer)
     pos      = scan(buffer,"=")
     var%name = buffer(1:pos-1);call upper_case(var%name)
     var%value= buffer(pos+1:)
   end function scan_cmd_variable

   function scan_input_variable(buffer)  result(var)
     implicit none
     character(len=*)                       :: buffer
     integer                                :: pos,len
     type(input_variable)                   :: var
     call s_blank_delete(buffer)
     pos      = scan(buffer,"=")
     var%name = buffer(1:pos-1)
     call upper_case(var%name)
     len=len_trim(buffer)
     if(len.gt.0)then
        if(buffer(len:len)==',')then
          var%value= buffer(pos+1:len-1)
       else
          var%value= buffer(pos+1:)
       endif
     endif
   end function scan_input_variable

   function check_cmd_vector_size(ndim,var) result(nargs)
     implicit none
     integer                                :: ndim,ncount,j,nargs
     type(input_variable)                   :: var
     logical                                :: iscalar
     if(ndim==1)then
        nargs=ndim
        return
     endif
     iscalar=(scan(var%value,",")==0)
     if(iscalar)then
        print*,"warning scalar in parse_cmd array:   ",trim(var%name)
        print*,"expecting a comma separated list of: ",ndim
     endif
     ncount=0
     do j=1,len(var%value)
        if(var%value(j:j)==",")ncount=ncount+1
     enddo
     nargs=ncount+1
     if(nargs/=ndim)then
        print*,"wrong dimensions parsing variable:   ",trim(var%name)
        print*,"expecting a comma separated list of: ",ndim
     endif
   end function check_cmd_vector_size

   subroutine upper_case(s)
     implicit none
     character              ch
     integer   ( kind = 4 ) i
     character ( len = * )  s
     integer   ( kind = 4 ) s_length
     s_length = len_trim ( s )
     do i = 1, s_length
        ch = s(i:i)
        call ch_cap ( ch )
        s(i:i) = ch
     end do
   end subroutine upper_case

   subroutine lower_case(s)
     implicit none
     integer   ( kind = 4 ) i
     character ( len = * )  s
     integer   ( kind = 4 ) s_length
     s_length = len_trim ( s )
     do i = 1, s_length
        call ch_low ( s(i:i) )
     end do
   end subroutine lower_case

   subroutine ch_cap(ch)
     implicit none
     character              ch
     integer   ( kind = 4 ) itemp
     itemp = iachar ( ch )
     if ( 97 <= itemp .and. itemp <= 122 ) then
        ch = achar ( itemp - 32 )
     end if
   end subroutine ch_cap

   subroutine ch_low ( ch )
     implicit none
     character ch
     integer ( kind = 4 ) i
     i = iachar ( ch )
     if ( 65 <= i .and. i <= 90 ) then
        ch = achar ( i + 32 )
     end if
   end subroutine ch_low

   subroutine s_blank_delete ( s )
     !! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
     !    All TAB characters are also removed.
     !    Input/output, character ( len = * ) S, the string to be transformed.
     implicit none
     character              ch
     integer   ( kind = 4 ) get
     integer   ( kind = 4 ) put
     character ( len = * )  s
     integer   ( kind = 4 ) s_length
     character, parameter :: tab = achar ( 9 )
     put = 0
     s_length = len_trim ( s )
     do get = 1, s_length
        ch = s(get:get)
        if ( ch /= ' ' .and. ch /= tab ) then
           put = put + 1
           s(put:put) = ch
        end if
     end do
     s(put+1:s_length) = ' '
     return
   end subroutine s_blank_delete


end module input_vars
