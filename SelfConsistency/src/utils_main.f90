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
   type(FermionicField)                     :: Gimp
   type(FermionicField)                     :: SigmaFull
   type(FermionicField)                     :: SigmaG0W0,SigmaG0W0dc
   type(FermionicField)                     :: SigmaGW,SigmaGW_C,SigmaGW_X
   type(FermionicField)                     :: SigmaGWdc,SigmaGW_Cdc,SigmaGW_Xdc
   type(FermionicField)                     :: SigmaDMFT
   !
   type(BosonicField)                       :: Wlat
   type(BosonicField)                       :: Ulat
   type(BosonicField)                       :: PiGG
   type(BosonicField)                       :: PiEDMFT
   type(BosonicField)                       :: curlyU
   !
   real(8)                                  :: density2set
   complex(8),allocatable                   :: densityLDA(:,:)
   complex(8),allocatable                   :: densityGW(:,:,:)
   real(8),allocatable                      :: densityDMFT(:,:,:)
   !
   complex(8),allocatable                   :: Vxc(:,:,:,:)
   complex(8),allocatable                   :: VH(:,:)
   !
   real(8),allocatable                      :: Umat(:,:)
   real(8),allocatable                      :: Kfunct(:,:)
   !
   logical                                  :: solve_DMFT=.false.
   !
   logical                                  :: calc_PiGG=.false.
   logical                                  :: merge_Pi=.false.
   logical                                  :: calc_W=.false.
   logical                                  :: calc_Wfull=.false.
   logical                                  :: calc_Wedmft=.false.
   logical                                  :: calc_Sigma=.false.
   logical                                  :: merge_Sigma=.false.

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
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine printHeader(Iteration)
      !
      implicit none
      integer,intent(in),optional           :: Iteration
      character(len=80)                     :: line
      character(:),allocatable              :: header
      integer                               :: i,Lsx,Ldx,Lcalc
      !
      if(present(Iteration))then
         header="<"
         Lcalc=len(" Iteration #: "//str(Iteration,3)//" ")
         Lsx=int((78-Lcalc)/2)-1
         Ldx=78-Lcalc-Lsx
         do i=1,Lsx
            header = header//"="
         enddo
         header = header//" Iteration #: "//str(Iteration,3)//" "
         do i=1,Ldx
            header = header//"="
         enddo
         header = header//">"
         write(LOGfile,"(A)") new_line("A")//header//new_line("A")
      else
         header="*"
         line="*"
         do i=1,79
            line = trim(line)//"*"
         enddo
         !
         Lcalc=len("Calculation type: "//trim(CalculationType))
         !
         Lsx=int((78-Lcalc)/2)-1
         Ldx=78-Lcalc-Lsx
         !
         do i=1,Lsx
            header = header//" "
         enddo
         header = header//"Calculation type: "//trim(CalculationType)
         do i=1,Ldx
            header = header//" "
         enddo
         header = header//"*"
         write(LOGfile,"(A)") line
         write(LOGfile,"(A)") header
         write(LOGfile,"(A)") line//new_line("A")//new_line("A")
      endif
      !
   end subroutine printHeader


   !---------------------------------------------------------------------------!
   !PURPOSE: looks for the current iteration number
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine initialize_DataStructure(ItStart,Itend)
      !
      implicit none
      integer,intent(out)                   :: ItStart
      integer,intent(out)                   :: Itend
      character(len=256)                    :: Itpath
      integer                               :: iter,Itfirst
      logical                               :: Itexist
      !
      do iter=0,1000
         Itpath = reg(pathDATA)//str(iter)
         call inquireDir(reg(Itpath),Itexist,hardstop=.false.,verb=.false.)
         if(.not.Itexist)then
            Itfirst = iter
            Exit
         endif
      enddo
      !
      if(FirstIteration.le.Itfirst)Itfirst=FirstIteration
      !
      Itpath = reg(pathDATA)//str(Itfirst)
      if(Itfirst.eq.0)then
         write(LOGfile,"(A)") "Brand new calculation. Initializing "//reg(Itpath)
      else
         write(LOGfile,"(A)") "Last iteration: "//str(Itfirst-1)//". Initializing "//reg(Itpath)
      endif
      call createDir(reg(Itpath))
      !
      ItStart = Itfirst
      !
      select case(reg(CalculationType))
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0")
            !
            if(ItStart.ne.0) stop "CalculationType is G0W0 but the starting iteration is not 0."
            Itend = 1
            !
         case("scGW")
            !
            Itend = LastIteration
            !
         case("DMFT+statU","DMFT+dynU","EDMFT","GW+EDMFT")
            !
            Itend = ItStart + 1
            !
      end select
      !
   end subroutine initialize_DataStructure


   !---------------------------------------------------------------------------!
   !PURPOSE: Initialize Lattice. I could have used the AllocateLattice in
   !         utils_fields but then als useless attributes would have been allocated
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine initialize_Lattice(Lttc,ItStart)
      !
      implicit none
      type(Lattice),intent(out)             :: Lttc
      integer,intent(in)                    :: ItStart
      integer                               :: m,n,mp,np,ib1,ib2
      integer                               :: iq_gamma_Hk,iq_gamma_XEPS
      !
      !
      write(LOGfile,"(A)") "---- initialize Lattice"
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
            call read_Hk(pathINPUT,alphaHk,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc,iq_gamma_Hk)
            !
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
            call read_Hk(pathINPUT,alphaHk,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc)
            !
            Lttc%Norb = size(Lttc%Hk,dim=1)
            Lttc%Nkpt = size(Lttc%Hk,dim=3)
            !
            Lttc%status=.true.
            !
         case("EDMFT","GW+EDMFT")
            !
            call read_Hk(pathINPUT,alphaHk,Lttc%Hk,Lttc%kpt,Lttc%Ek,Lttc%Zk,Lttc%Hloc,iq_gamma_Hk)
            !
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
      !
      allocate(PhysicalUelement(Lttc%Norb**2,Lttc%Norb**2));PhysicalUelement=.false.
      do m=1,Lttc%Norb
         do n=1,Lttc%Norb
            do mp=1,Lttc%Norb
               do np=1,Lttc%Norb
                  !
                  ib1 = mp + Lttc%Norb*(m-1)
                  ib2 = np + Lttc%Norb*(n-1)
                  !
                  if((mp.eq.m).and.(np.eq.n)) PhysicalUelement(ib1,ib2)=.true.
                  if((mp.eq.np).and.(m.eq.n)) PhysicalUelement(ib1,ib2)=.true.
                  if((mp.eq.n).and.(m.eq.np)) PhysicalUelement(ib1,ib2)=.true.
                  !
               enddo
            enddo
         enddo
      enddo
      !
      if(ItStart.eq.0)call calc_Glda(0d0,Beta,Lttc)
      !
      allocate(densityLDA(Lttc%Norb,Lttc%Norb));densityLDA=czero
      allocate(densityGW(Lttc%Norb,Lttc%Norb,Nspin));densityGW=czero
      allocate(densityDMFT(Lttc%Norb,Lttc%Norb,Nspin));densityDMFT=0d0
      !
   end subroutine initialize_Lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Initialize the Fields depending on the starting iteration
   !TEST ON: 16-10-2020
   !---------------------------------------------------------------------------!
   subroutine initialize_Fields(ItStart)
      !
      implicit none
      integer,intent(in)                    :: ItStart
      logical                               :: filexists
      !
      !
      write(LOGfile,"(A)") "---- initialize Fields"
      !
      !
      select case(reg(CalculationType))
         case default
            !
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
            !
         case("G0W0","scGW")
            !
            !Unscreened interaction
            call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(Umodel)stop "U model is implemented only for non-GW (fully local) screened calculations."
            if(Uspex) call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.false.)
            !
            !Fully screened interaction
            call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            !Polarization
            call AllocateBosonicField(PiGG,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.false.,Beta=Beta)
            !
            !Lattice Gf
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(ItStart.eq.0) call calc_Gmats(Glat,Crystal)
            if(ItStart.ne.0) call read_FermionicField(Glat,reg(pathDATA)//str(ItStart-1)//"/","Glat",kpt=Crystal%kpt)
            !
            !Logical Flags
            calc_PiGG = .true.
            calc_W = .true.
            calc_Wfull = .true.
            calc_Sigma = .true.
            !
         case("DMFT+statU")
            !
            !Hubbard interaction
            allocate(Umat(Crystal%Norb**2,Crystal%Norb**2));Umat=0d0
            if(Uspex) call read_U_spex(Umat)
            if(Umodel)then
               call inquireFile(reg(pathINPUT)//"Umat_model.DAT",filexists,hardstop=.false.)
               if(filexists)then
                  call read_Matrix(Umat,reg(pathINPUT)//"Umat_model.DAT")
               else
                  call build_Uscr(Umat,Uaa,Uab,J)
                  call dump_Matrix(Umat,reg(pathINPUT)//"Umat_model.DAT")
               endif
            endif
            !
            if(ItStart.ne.0)then
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_dw.DAT")
               !
            else
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call calc_Gmats(Glat,Crystal)
               !
            endif
            !
            !Logical Flags
            solve_DMFT = .true.
            !
         case("DMFT+dynU")
            !
            !Unscreened interaction
            call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
            if(Uspex) call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.true.)
            if(Umodel)then
               call inquireFile(reg(pathINPUT)//"Uloc_mats_model.DAT",filexists,hardstop=.false.)
               if(filexists)then
                  call read_BosonicField(Ulat,reg(pathINPUT),"Uloc_mats_model.DAT")
               else
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph)
                  call dump_BosonicField(Ulat,reg(pathINPUT),"Uloc_mats_model.DAT")
               endif
            endif
            !
            if(ItStart.ne.0)then
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_dw.DAT")
               !
            else
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call calc_Gmats(Glat,Crystal)
               !
            endif
            !
            !Logical Flags
            solve_DMFT = .true.
            !
         case("EDMFT")
            !
            !Unscreened interaction
            call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,Beta=Beta)
            if(Uspex) call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.true.)
            if(Umodel)then
               call inquireFile(reg(pathINPUT)//"Uloc_mats_model.DAT",filexists,hardstop=.false.)
               if(filexists)then
                  call read_BosonicField(Ulat,reg(pathINPUT),"Uloc_mats_model.DAT")
               else
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph)
                  call dump_BosonicField(Ulat,reg(pathINPUT),"Uloc_mats_model.DAT")
               endif
            endif
            !
            !Fully screened interaction
            call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            if(ItStart.ne.0)then
               !
               !Polarization
               call AllocateBosonicField(PiEDMFT,Crystal%Norb,Crystal%iq_gamma,Nmats,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               call read_BosonicField(PiEDMFT,reg(pathDATA)//str(ItStart-1)//"/","PiEDMFT")
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_dw.DAT")
               !
            else
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call calc_Gmats(Glat,Crystal)
               !
            endif
            !
            !Logical Flags
            calc_Wedmft = .true.
            solve_DMFT = .true.
            !
         case("GW+EDMFT")
            !
            !Unscreened interaction
            call AllocateBosonicField(Ulat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(Umodel)stop "U model is implemented only for non-GW (fully local) screened calculations."
            if(Uspex) call read_U_spex(Ulat,save2readable=verbose,LocalOnly=.false.)
            !
            !Fully screened interaction
            call AllocateBosonicField(Wlat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            if(ItStart.ne.0)then
               !
               !Polarization
               call AllocateBosonicField(PiGG,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               call AllocateBosonicField(PiEDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               call read_BosonicField(PiEDMFT,reg(pathDATA)//str(ItStart-1)//"/","PiEDMFT.DAT")
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT_dw.DAT")
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(Glat,reg(pathDATA)//str(ItStart-1)//"/","Glat",kpt=Crystal%kpt)
               !
            else
               !
               !Polarization
               call AllocateBosonicField(PiGG,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call calc_Gmats(Glat,Crystal)
               !
            endif
            !
            !Logical Flags
            calc_PiGG = .true.
            calc_Wfull = .true.
            calc_Sigma = .true.
            solve_DMFT = .true.
            if(ItStart.gt.0)then
               merge_Pi = .true.
               merge_Sigma = .true.
            endif
            !
      end select
      !
      calc_W = calc_Wedmft .or. calc_Wfull
      !
      if(ItStart.eq.0)then
         densityLDA = Glat%N_s(:,:,1) + Glat%N_s(:,:,2)
         call dump_Matrix(densityLDA,reg(pathINPUT)//"n_LDA.DAT")
         densityGW=czero
      else
         call read_Matrix(densityLDA,reg(pathINPUT)//"n_LDA.DAT")
         densityGW=Glat%N_s
      endif
      !
      if(look4dens%TargetDensity.eq.0d0)then
         look4dens%TargetDensity = trace(densityLDA)
      else
         if(ItStart.eq.0)call set_density(Glat%mu,Beta,Crystal,look4dens)
      endif
      !
   end subroutine initialize_Fields


   !---------------------------------------------------------------------------!
   !PURPOSE: Join the all the component of the self-energy
   !TEST ON: 27-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_SigmaFull(Iteration)
      !
      implicit none
      integer,intent(in)                    :: Iteration
      integer                               :: iorb,jorb,ik,iw,ispin
      !
      select case(reg(CalculationType))
         case default
            stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
         case("G0W0","scGW","GW+EDMFT")
            !
            call AllocateFermionicField(SigmaFull,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            !
            if(Iteration.eq.0)then
               !
               if(.not.SigmaG0W0%status) stop "calc_SigmaFull: SigmaG0W0 not properly initialized."
               !
               !$OMP PARALLEL DEFAULT(NONE),&
               !$OMP SHARED(SigmaFull,SigmaG0W0,VH,Vxc),&
               !$OMP PRIVATE(iorb,jorb,ik,iw,ispin)
               !$OMP DO
               do ispin=1,Nspin
                  do ik=1,SigmaFull%Nkpt
                     do iw=1,SigmaFull%Npoints
                        do iorb=1,SigmaFull%Norb
                           do jorb=1,SigmaFull%Norb
                              SigmaFull%wks(iorb,jorb,iw,ik,ispin) =  + SigmaG0W0%wks(iorb,jorb,iw,ik,ispin)   &
                                                                      - Vxc(iorb,jorb,ik,ispin)                &
                                                                      + VH(iorb,jorb)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL
               !
            elseif(Iteration.gt.0)then
               !
               if(.not.SigmaGW%status) stop "calc_SigmaFull: SigmaGW not properly initialized."
               if(.not.SigmaG0W0dc%status) stop "calc_SigmaFull: SigmaG0W0dc not properly initialized."
               if(.not.SigmaG0W0%status) stop "calc_SigmaFull: SigmaG0W0 not properly initialized."
               !
               !$OMP PARALLEL DEFAULT(NONE),&
               !$OMP SHARED(SigmaFull,SigmaGW,SigmaG0W0dc,SigmaG0W0,VH,Vxc),&
               !$OMP PRIVATE(iorb,jorb,ik,iw,ispin)
               !$OMP DO
               do ispin=1,Nspin
                  do ik=1,SigmaFull%Nkpt
                     do iw=1,SigmaFull%Npoints
                        do iorb=1,SigmaFull%Norb
                           do jorb=1,SigmaFull%Norb
                              SigmaFull%wks(iorb,jorb,iw,ik,ispin) =  + SigmaGW%wks(iorb,jorb,iw,ik,ispin)     &
                                                                      - SigmaG0W0dc%wks(iorb,jorb,iw,ik,ispin) &
                                                                      + SigmaG0W0%wks(iorb,jorb,iw,ik,ispin)   &
                                                                      - Vxc(iorb,jorb,ik,ispin)                &
                                                                      + VH(iorb,jorb)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL
               !
            endif
            !
            call FermionicKsum(SigmaFull)
            !
            !Dump the local projection of SigmaFull
            call dump_FermionicField(SigmaFull,1,reg(pathDATA)//str(Iteration)//"/","SigmaLoc_up.DAT")
            call dump_FermionicField(SigmaFull,2,reg(pathDATA)//str(Iteration)//"/","SigmaLoc_dw.DAT")
            !
            !Dump the whole SigmaFull
            call dump_FermionicField(SigmaFull,reg(pathDATA)//str(Iteration)//"/","Sigma",.true.,Crystal%kpt)
            !
         case("DMFT+statU","DMFT+dynU","EDMFT")
            !
            call AllocateFermionicField(SigmaFull,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
            !
            if(Iteration.eq.0)then
               !
               call clear_attributes(SigmaFull)
               !
            elseif(Iteration.gt.0)then
               !
               if(.not.SigmaDMFT%status) stop "calc_SigmaFull: SigmaDMFT not properly initialized."
               !
               !$OMP PARALLEL DEFAULT(NONE),&
               !$OMP SHARED(SigmaFull,SigmaDMFT),&
               !$OMP PRIVATE(iorb,jorb,iw,ispin)
               !$OMP DO
               do ispin=1,Nspin
                  do iw=1,SigmaFull%Npoints
                     do iorb=1,SigmaFull%Norb
                        do jorb=1,SigmaFull%Norb
                           SigmaFull%ws(iorb,jorb,iw,ispin) = SigmaDMFT%ws(iorb,jorb,iw,ispin)
                        enddo
                     enddo
                  enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL
               !
            endif
            call DeallocateFermionicField(SigmaDMFT)
            !
            !Not dumping anything since SigmaDMFT is already present
            !
      end select
      !
   end subroutine calc_SigmaFull




end module utils_main











!call dump_FermionicField(Glat,1,reg(pathDATA)//str(ItStart)//"/","Glat_loc_up.DAT")
!call dump_FermionicField(Glat,2,reg(pathDATA)//str(ItStart)//"/","Glat_loc_dn.DAT")
!call dump_FermionicField(Glat,reg(pathDATA)//str(ItStart)//"/","Glat",binfmt=.true.)
!call dump_FermionicField(Glat,reg(pathDATA)//str(ItStart)//"/","Glat",binfmt=.false.)
!call read_FermionicField(Glat,1,reg(pathDATA)//str(ItStart)//"/","Glat_loc_up.DAT")
!call read_FermionicField(Glat,2,reg(pathDATA)//str(ItStart)//"/","Glat_loc_dn.DAT")
!call dump_FermionicField(Glat,1,reg(pathDATA)//str(ItStart)//"/","Glat_loc_up_testread_LOC.DAT")
!call dump_FermionicField(Glat,2,reg(pathDATA)//str(ItStart)//"/","Glat_loc_dn_testread_LOC.DAT")
!call read_FermionicField(Glat,reg(pathDATA)//str(ItStart)//"/","Glat")
!call dump_FermionicField(Glat,reg(pathDATA)//str(ItStart)//"/","Glat_testread_Kdep",binfmt=.true.)
