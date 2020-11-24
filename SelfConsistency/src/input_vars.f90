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
   character(len=255),private               :: p_buffer
   character(len=7),private                 :: file_status
   integer,parameter,private                :: pos_comment=46
   logical,private                          :: IOinput=.true.
   integer,allocatable,private              :: commentPos(:)
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
   !Site and Orbital space
   integer,public                           :: Nsite
   logical,public                           :: ExpandImpurity
   logical,public                           :: RotateHloc
   logical,public                           :: AFMselfcons
   integer,public,allocatable               :: SiteNorb(:)
   character(len=2),public,allocatable      :: SiteName(:)
   integer,public,allocatable               :: SiteOrbs(:,:)
   !
   !Equivalent lattice indexes
   type(Equivalent),public                  :: EqvGWndx
   !
   !Imaginary time and frequency meshes
   real(8),public                           :: Beta
   integer,public                           :: NtauF
   integer,public                           :: NtauB
   logical,public                           :: tau_uniform
   real(8),public                           :: wmatsMax
   integer,public                           :: Nmats
   integer,public                           :: Nreal
   real(8),public                           :: wrealMax
   real(8),public                           :: eta
   !
   !Density lookup
   type(musearch),public                    :: look4dens
   !
   !Interaction variables
   logical,public                           :: UfullStructure
   logical,public                           :: Umodel
   logical,public                           :: Uspex
   integer,public                           :: Nphonons
   real(8),public,allocatable               :: g_eph(:)
   real(8),public,allocatable               :: wo_eph(:)
   real(8),public                           :: Uaa
   real(8),public                           :: Uab
   real(8),public                           :: J
   logical,public                           :: Kdiag
   !
   !Double counting types, divergencies, scaling coefficients
   character(len=256),public                :: VH_type
   character(len=256),public                :: DC_type
   logical,public                           :: HandleGammaPoint
   real(8),public                           :: alphaPi
   real(8),public                           :: alphaSigma
   real(8),public                           :: alphaHk
   !
   !Variables for the fit on Delta, Gimp, Simp
   character(len=256),public                :: DeltaFit
   integer,public                           :: Nfit
   real(8),public                           :: ReplaceTail_Gimp
   real(8),public                           :: ReplaceTail_Simp
   !
   !Paths (directories must end with "/") and loop variables
   integer,public                           :: FirstIteration
   integer,public                           :: LastIteration
   character(len=256),public                :: pathINPUT!="InputFiles/"
   character(len=256),public                :: pathDATA!="Iterations/"
   integer,public                           :: LOGfile
   real(8),public                           :: Mixing_curlyG
   real(8),public                           :: Mixing_Delta
   real(8),public                           :: Mixing_curlyU
   logical,public                           :: skipLattice
   logical,public                           :: printSfull
   !
   !Variables related to the impurity solver
   type(QMC),public                         :: Solver
   !
   !Variables not to be readed
   logical,public                           :: paramagneticSPEX=.true.
   logical,public                           :: XEPSisread=.false.
   logical,public                           :: solve_DMFT=.true.
   logical,public                           :: bosonicSC=.false.
   !Once the it=0 problem is solved this can be chosen by the user
   logical,public                           :: Ustart=.true.

   !---------------------------------------------------------------------------!
   !PURPOSE: Internal Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: save_InputFile
   public :: print_InputFile
   public :: delete_Input
   public :: parse_Cmd_variable
   public :: parse_Input_variable
   public :: read_InputFile

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the Inputfile
   !---------------------------------------------------------------------------!
   subroutine read_InputFile(InputFile)
      use utils_misc
      use parameters
      implicit none
      character(len=*)                      :: InputFile
      integer                               :: unit
      integer                               :: isite,iset,iph
      integer,allocatable                   :: tmpOrbs(:)
      !
      write(LOGfile,"(A)") new_line("A")//"Reading InputFile"//new_line("A")
      !
      call parse_input_variable(CalculationType,"CALC_TYPE",InputFile,default="GW+EDMFT",comment="Calculation type. Avalibale: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT.")
      if((reg(CalculationType).eq."G0W0").or.(reg(CalculationType).eq."scGW")) solve_DMFT=.false.
      if((reg(CalculationType).eq."EDMFT").or.(reg(CalculationType).eq."GW+EDMFT")) then
         bosonicSC=.true.
         Solver%retarded=1!.true.
         Solver%nnt_meas=1!.true.
      endif
      if(reg(CalculationType).eq."DMFT+dynU")Solver%retarded=1!.true.
      !
      !OMP parallelization
      call execute_command_line(" lscpu | grep 'CPU(s):       ' | awk '{print $2}' > Nthread.used ")
      unit=free_unit()
      open(unit,file=reg("Nthread.used"),form="formatted",action="read",status="old")
      read(unit,*)Nthread
      close(unit)
      !
      !K-points
      call add_separator()
      call parse_input_variable(Nkpt3,"NKPT3",InputFile,default=[8,8,8],comment="Number of K-points per dimension.")
      call parse_input_variable(UseXepsKorder,"XEPS_KORDER",InputFile,default=.true.,comment="Flag to use the K-point ordering of XEPS.DAT if present.")
      !
      !Site and Orbital space
      call add_separator()
      call append_to_input_list(Nspin,"NSPIN","Number of spins (fixed to 2).")
      call parse_input_variable(Nsite,"NSITE",InputFile,default=1,comment="Number of impurities in the lattice.")
      call parse_input_variable(ExpandImpurity,"EXPAND",InputFile,default=.false.,comment="Flag to use a single impurity solution for all the sites of the lattice. Only indexes for site 1 readed.")
      call parse_input_variable(RotateHloc,"ROTATE",InputFile,default=.false.,comment="Solve the impurity problem in the basis where H(R=0) is diagonal.")
      call parse_input_variable(AFMselfcons,"AFM",InputFile,default=.false.,comment="Flag to use  the AFM self-consistency by flipping the spin. Requires input with doubled unit cell.")
      if(ExpandImpurity)then
         allocate(SiteNorb(1));SiteNorb=0
         allocate(SiteName(1))
         call parse_input_variable(SiteNorb(1),"NORB_1",InputFile,default=1,comment="Number of orbitals in site number 1")
         call parse_input_variable(SiteName(1),"NAME_1",InputFile,default="El",comment="Chemical species of the site number 1")
         allocate(SiteOrbs(1,SiteNorb(1)));SiteOrbs=0
         allocate(tmpOrbs(1:SiteNorb(1)));tmpOrbs=0
         call parse_input_variable(SiteOrbs(1,1:SiteNorb(1)),"ORBS_1",InputFile,default=tmpOrbs,comment="Lattice orbital indexes of site number 1")
         deallocate(tmpOrbs)
      else
         allocate(SiteNorb(Nsite));SiteNorb=0
         allocate(SiteName(Nsite))
         do isite=1,Nsite
            call parse_input_variable(SiteNorb(isite),"NORB_"//str(isite),InputFile,default=1,comment="Number of orbitals in site number "//str(isite))
            call parse_input_variable(SiteName(isite),"NAME_"//str(isite),InputFile,default="El",comment="Chemical species of the site number "//str(isite))
         enddo
         allocate(SiteOrbs(Nsite,maxval(SiteNorb)));SiteOrbs=0
         do isite=1,Nsite
            allocate(tmpOrbs(1:SiteNorb(isite)));tmpOrbs=0
            call parse_input_variable(SiteOrbs(isite,1:SiteNorb(isite)),"ORBS_"//str(isite),InputFile,default=tmpOrbs,comment="Lattice orbital indexes of site number "//str(isite))
            deallocate(tmpOrbs)
         enddo
      endif
      !
      !Equivalent lattice indexes
      call add_separator()
      call parse_input_variable(EqvGWndx%para,"PARAMAGNET",InputFile,default=1,comment="Integer flag to impose spin symmetry.")
      call parse_input_variable(EqvGWndx%hseed,"H_SEED",InputFile,default=0d0,comment="Seed to break spin symmetry (persistent if non zero).")
      call parse_input_variable(EqvGWndx%Nset,"EQV_SETS",InputFile,default=1,comment="Number of sets of locally equivalent lattice orbitals.")
      if(EqvGWndx%Nset.gt.0)then
         allocate(EqvGWndx%SetNorb(EqvGWndx%Nset))
         do iset=1,EqvGWndx%Nset
            call parse_input_variable(EqvGWndx%SetNorb(iset),"EQV_NORB_"//str(iset),InputFile,default=1,comment="Number of equivalent lattice orbitals in the set number "//str(iset))
         enddo
         allocate(EqvGWndx%SetOrbs(EqvGWndx%Nset,maxval(EqvGWndx%SetNorb)));EqvGWndx%SetOrbs=0
         do iset=1,EqvGWndx%Nset
            allocate(tmpOrbs(1:EqvGWndx%SetNorb(iset)));tmpOrbs=0
            call parse_input_variable(EqvGWndx%SetOrbs(iset,1:EqvGWndx%SetNorb(iset)),"EQV_ORBS_"//str(iset),InputFile,default=tmpOrbs,comment="Lattice orbital indexes of equivalent set number "//str(iset))
            deallocate(tmpOrbs)
         enddo
      endif
      if(EqvGWndx%para.eq.1)EqvGWndx%S=.true.
      EqvGWndx%Gfoffdiag=.not.ExpandImpurity
      !
      !Imaginary time and frequency meshes
      call add_separator()
      call parse_input_variable(Beta,"BETA",InputFile,default=10.d0,comment="Inverse temperature in 1/eV.")
      call parse_input_variable(wmatsMax,"MAX_WMATS",InputFile,default=100.d0,comment="Maximum value of the Matsubara frequency mesh.")
      Nmats = int(Beta*wmatsMax/(2d0*pi))
      call append_to_input_list(Nmats,"NMATS","Number of points on the imaginary frequency axis. User cannot set this as its computed from MAX_WMATS and BETA.")
      call parse_input_variable(NtauF,"NTAU_F",InputFile,default=int(100+Nmats),comment="Number of points on the imaginary time axis for Fermionic quantities. Its gonna be made odd. Optimal if >~ NMATS.")
      if(mod(NtauF,2).eq.0)NtauF=NtauF+1
      if(mod(NtauF-1,4).ne.0)NtauF=NtauF+mod(NtauF-1,4)
      call parse_input_variable(NtauB,"NTAU_B",InputFile,default=int(2d0*pi*Nmats),comment="Number of points on the imaginary time axis for Bosonic quantities. Its gonna be made odd. Optimal if >= to 2*PI*NMATS.")
      if(mod(NtauB,2).eq.0)NtauB=NtauB+1
      if(mod(NtauB-1,4).ne.0)NtauB=NtauB+mod(NtauB-1,4)
      call parse_input_variable(tau_uniform,"TAU_UNIF",InputFile,default=.false.,comment="Flag to use a non-tau_uniform mesh on the imaginary time axis.")
      call parse_input_variable(Nreal,"NREAL",InputFile,default=2000,comment="Number of points on the real frequency axis.")
      call parse_input_variable(wrealMax,"MAX_WREAL",InputFile,default=10.d0,comment="Maximum absolute value of the real frequency mesh.")
      call parse_input_variable(eta,"ETA",InputFile,default=0.04d0,comment="Real frequency broadening.")
      !
      !Density lookup
      call add_separator()
      call parse_input_variable(look4dens%mu,"MU",InputFile,default=0d0,comment="Chemical potential.")
      call parse_input_variable(look4dens%TargetDensity,"N_READ",InputFile,default=0d0,comment="Target density per site lookup is switched on to this value if its >0d0. Otherwise mu will be kept fixed.")
      call parse_input_variable(look4dens%quickloops,"N_QUICK",InputFile,default=1,comment="Integer flag to switch on the quick density lookup within the solver.")
      call parse_input_variable(look4dens%densityRelErr,"N_ERR",InputFile,default=0.01d0,comment="Relative error on the target density. Better if not lower than 1e-3.")
      call parse_input_variable(look4dens%muStep,"MU_STEP",InputFile,default=0.2d0,comment="Initial chemical potential step in the density lookup.")
      call parse_input_variable(look4dens%muIter,"MU_ITER",InputFile,default=50,comment="Maximum number of iterations in the density lookup.")
      call parse_input_variable(look4dens%muTime,"MU_TIME",InputFile,default=0.5d0,comment="Minutes of solver runtime in the density lookup.")
      if(ExpandImpurity)look4dens%TargetDensity = look4dens%TargetDensity*Nsite
      !
      !Interaction variables
      call add_separator()
      !call parse_input_variable(Ustart,"USTART",InputFile,default=.false.,comment="Flag to use the local Ucrpa interaction as effective interaction in the 0th iteration.")
      call parse_input_variable(UfullStructure,"U_FULL",InputFile,default=.true.,comment="Flag to check for inverted Re/Im parity in SPEX Ucrpa.")
      call parse_input_variable(Umodel,"U_MODEL",InputFile,default=.false.,comment="Flag to build the screening from user chosen phononic modes.")
      call parse_input_variable(Uspex,"U_SPEX",InputFile,default=.true.,comment="Flag to read SPEX Ucrpa.")
      call parse_input_variable(Kdiag,"K_DIAG",InputFile,default=.false.,comment="Flag to screen only orbital-diagonal channels.")
      if(Umodel.and.Uspex) stop "Make up your mind, U_MODEL or U_SPEX ?"
      if(Umodel)then
         call parse_input_variable(Uaa,"UAA",InputFile,default=5d0,comment="Interaction between same orbital and opposite spin electrons.")
         call parse_input_variable(Uab,"UAB",InputFile,default=4d0,comment="Interaction between different orbital and opposite spin electrons.")
         call parse_input_variable(J,"JH",InputFile,default=0.5d0,comment="Hund's coupling.")
         call parse_input_variable(Nphonons,"N_PH",InputFile,default=3,comment="Number of custom phononic modes.")
         allocate(g_eph(Nphonons));g_eph=0d0
         allocate(wo_eph(Nphonons));wo_eph=0d0
         do iph=1,Nphonons
            call parse_input_variable(g_eph,"G_PH",InputFile,comment="Custom phononic couplings.")
            call parse_input_variable(wo_eph,"WO_PH",InputFile,comment="Custom phononic energies.")
         enddo
      endif
      !
      !Double counting types, divergencies, scaling coefficients
      call add_separator()
      call parse_input_variable(VH_type,"VH_TYPE",InputFile,default="Ustatic",comment="Hartree term mismatch between GoWo and scGW. Available: Ubare, Ustatic, Ubare_SPEX, Ustatic_SPEX.")
      call parse_input_variable(DC_type,"DC_TYPE",InputFile,default="GlocWloc",comment="Local GW self-energy which is replaced by DMFT self-energy. Avalibale: GlocWloc, Sloc.")
      call parse_input_variable(HandleGammaPoint,"SMEAR_GAMMA",InputFile,default=.true.,comment="Remove the interaction divergence at the Gamma point.")
      call parse_input_variable(alphaPi,"ALPHA_PI",InputFile,default=1d0,comment="Fraction of the EDMFT polarization substituted within the lattice one.")
      call parse_input_variable(alphaSigma,"ALPHA_SIGMA",InputFile,default=1d0,comment="Fraction of the EDMFT self-energy substituted within the lattice one.")
      call parse_input_variable(alphaHk,"ALPHA_HK",InputFile,default=1d0,comment="Rescaling of the non-interacting Hamiltonian.")
      !
      !Variables for the fit
      call parse_input_variable(DeltaFit,"DELTA_FIT",InputFile,default="Analytic",comment="Fit to extract the local energy in GW+EDMFT calculations. Available: Analytic, Moments.")
      call parse_input_variable(Nfit,"NFIT",InputFile,default=15,comment="Number of bath levels (Analytic) or coefficient (Moments).")
      call parse_input_variable(ReplaceTail_Gimp,"WTAIL_GIMP",InputFile,default=100d0,comment="Frequency value above which the tail of Gimp is replaced. If =0d0 the tail is not replaced. A good number is around 100.")
      call parse_input_variable(ReplaceTail_Simp,"WTAIL_SIMP",InputFile,default=20d0,comment="Frequency value above which the tail of Simp is replaced. If =0d0 the tail is not replaced. A good number is around 10.")
      !
      !Paths and loop variables
      call add_separator()
      call parse_input_variable(pathINPUT,"PATH_INPUT",InputFile,default="InputFiles",comment="Folder within cwd where to look for input files.")
      call parse_input_variable(pathDATA,"PATH_DATA",InputFile,default="Iterations",comment="Folder within cwd where to store data.")
      call parse_input_variable(FirstIteration,"START_IT",InputFile,default=0,comment="First iteration. If its non zero the code will look for the last item in PATH_DATA/item and start there.")
      call parse_input_variable(LastIteration,"LAST_IT",InputFile,default=100,comment="Last iteration.")
      call parse_input_variable(LOGfile,"LOGFILE",InputFile,default=6,comment="Standard output redirection unit. Use 6 to print to terminal. Not used yet.")
      call parse_input_variable(Mixing_curlyG,"MIX_G",InputFile,default=0d0,comment="Fraction of the old iteration curlyG.")
      call parse_input_variable(Mixing_Delta,"MIX_D",InputFile,default=0.5d0,comment="Fraction of the old iteration Delta.")
      call parse_input_variable(Mixing_curlyU,"MIX_U",InputFile,default=0.5d0,comment="Fraction of the old iteration curlyU.")
      call parse_input_variable(skipLattice,"SKIP_LATT",InputFile,default=.false.,comment="Skip the lattice summation and assuming good the existing Gloc and Wloc. Not used yet.")
      call parse_input_variable(printSfull,"PRINT_SFULL",InputFile,default=.false.,comment="Print the full k-dependent self-energy (binfmt) at each iteration.")
      !
      !Variables related to the impurity solver
      call add_separator()
      call parse_input_variable(Solver%Norder,"NORDER",InputFile,default=10,comment="Maximum perturbation order measured.")
      call parse_input_variable(Solver%Nmeas,"NMEAS",InputFile,default=1000,comment="Sweeps where expensive measurments are not performed.")
      call parse_input_variable(Solver%Ntherm,"NTHERM",InputFile,default=100,comment="Thermalization cycles. Each cycle performs NMEAS sweeps.")
      call parse_input_variable(Solver%Nshift,"NSHIFT",InputFile,default=2,comment="Attempted shift moves at each sweep.")
      call parse_input_variable(Solver%PrintTime,"PRINT_TIME",InputFile,default=10,comment="Minutes that have to pass before observables are updated and stored.")
      call parse_input_variable(Solver%binlength,"BINLENGTH",InputFile,default=4,comment="If >0 the Green's function at itau will be the average within +/-binlength.")
      call parse_input_variable(Solver%binstart,"BINSTART",InputFile,default=100,comment="Tau points skipped at the beginning and end of the Green's function.")
      call append_to_input_list(Solver%retarded,"RETARDED","Integer flag to include the frequency dependent part of the interaction. User cannot set this as its deduced from CALC_TYPE.")
      call append_to_input_list(Solver%nnt_meas,"NNT_MEAS","Integer flag to switch on the measurement of the susceptibility. User cannot set this as its deduced from CALC_TYPE.")
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
      call code_version()
      call save_InputFile(InputFile)
      !
      !The last slash is uneffective for some reason
      pathINPUT=trim(pathINPUT)//"/"
      pathDATA=trim(pathDATA)//"/"
      !
   end subroutine read_InputFile


   !---------------------------------------------------------------------------!
   !PURPOSE: Somethig that could save my ass
   !---------------------------------------------------------------------------!
   subroutine code_version()
      implicit none
      include "revision_SelfCons.inc"
      integer(4),dimension(8)                  :: dummy
      integer(4)                               :: year
      integer(4)                               :: mese
      integer(4)                               :: day
      integer(4)                               :: h
      integer(4)                               :: m
      integer(4)                               :: s
      integer(4)                               :: ms
      character(len=9),parameter,dimension(12) :: month = (/ &
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
      write(*,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
           "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
      write(*,*)""
      open(10,file="code_version.inc")
      write(10,"(A)")"CODE VERSION: "//trim(adjustl(trim(revision)))
      write(10,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
           "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
      write(10,*)""
      close(10)
   end subroutine code_version


   !---------------------------------------------------------------------------!
   !PURPOSE: Field Separator
   !---------------------------------------------------------------------------!
   subroutine add_separator()
      implicit none
      integer             :: sizeOld,sizeNew
      integer,allocatable :: prvPos(:)
      !
      if(allocated(commentPos))then
         sizeOld=size(commentPos)
         sizeNew=sizeOld+1
         allocate(prvPos(sizeOld))
         prvPos=commentPos
         deallocate(commentPos)
         allocate(commentPos(sizeNew))
         commentPos(1:sizeOld) = prvPos
         commentPos(sizeNew) = default_list%size
      else
         allocate(commentPos(1))
         commentPos(1) = default_list%size
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
      character(len=*),optional              :: file
      type(input_list),optional              :: list
      integer                                :: counter,size
      type(input_node),pointer               :: c
      if(present(list))then
         c => list%root%next
      else
         c => default_list%root%next
      endif
      counter = 0
      file_status='replace'
      size=default_list%size
      if(present(list))size=list%size
      if(size>0)then
         do
            if(.not.associated(c))exit
            counter=counter+1
            if(present(file))then
               call print_input_node(c,file)
               if(any(commentPos.eq.counter))call print_separator(file)
            else
               call print_input_node(c)
               if(any(commentPos.eq.counter))call print_separator()
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
      type(input_list),optional              :: list
      type(input_node),pointer               :: p,c
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
     integer,dimension(:)                        :: variable
     integer,dimension(size(variable)),optional  :: default
     character(len=*)                            :: name
     character(len=*),optional                   :: comment
     character(len=*)                            :: file
     character(len=len(name))                    :: name_
     type(input_variable)                        :: var
     integer                                     :: unit,pos,j,ndim,nargs,pos0,iarg
     integer                                     :: status
     logical                                     :: bool
     character(len=255)                          :: buffer
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
     real(8),dimension(:)                        :: variable
     real(8),dimension(size(variable)),optional  :: default
     character(len=*)                            :: name
     character(len=*),optional                   :: comment
     character(len=*)                            :: file
     character(len=len(name))                    :: name_
     type(input_variable)                        :: var
     integer                                     :: unit,pos,j,ndim,nargs,pos0,iarg
     integer                                     :: status
     logical                                     :: bool
     character(len=255)                          :: buffer
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
     logical,dimension(:)                        :: variable
     logical,dimension(size(variable)),optional  :: default
     character(len=*)                            :: name
     character(len=*),optional                   :: comment
     character(len=*)                            :: file
     character(len=len(name))                    :: name_
     type(input_variable)                        :: var
     integer                                     :: unit,pos,j,ndim,nargs,pos0,iarg
     integer                                     :: status
     logical                                     :: bool
     character(len=255)                          :: buffer
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
     integer,dimension(:)                        :: variable
     integer,dimension(size(variable)),optional  :: default
     character(len=*)                            :: name
     character(len=len(name))                    :: name_
     type(input_variable)                        :: var
     integer                                     :: i,j,ndim,nargs,pos0,iarg
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
     real(8),dimension(:)                        :: variable
     real(8),dimension(size(variable)),optional  :: default
     character(len=*)                            :: name
     character(len=len(name))                    :: name_
     type(input_variable)                        :: var
     integer                                     :: i,j,ndim,nargs,pos0,iarg
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
     logical,dimension(:)                        :: variable
     logical,dimension(size(variable)),optional  :: default
     character(len=*)                            :: name
     character(len=len(name))                    :: name_
     type(input_variable)                        :: var
     integer                                     :: i,j,ndim,nargs,pos0,iarg
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
   subroutine print_separator(file)
     use utils_misc
     implicit none
     character(len=*),optional              :: file
     integer                                :: unit,i
     character(len=pos_comment+1)           :: line
     line="#"
     do i=1,pos_comment-1
        line = trim(line)//"-"
     enddo
     line = trim(line)//"#"
     if(present(file))then
        unit=free_unit()
        open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
        write(unit,"(1x,A)")
        write(unit,"(1x,A)")trim(line)
        write(unit,"(1x,A)")
        close(unit)
     else
        write(unit,"(1x,A)")
        write(unit,"(1x,A)")trim(line)
        write(unit,"(1x,A)")
     endif
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
