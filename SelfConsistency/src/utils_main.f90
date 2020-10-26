module utils_main

   use module_container
   implicit none
   public

   !===========================================================================!

   ! COMMENTS:
   !
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
   type(FermionicField)                     :: Glat
   type(FermionicField)                     :: Gimp
   type(FermionicField)                     :: SigmaFull
   type(FermionicField)                     :: SigmaG0W0
   type(FermionicField)                     :: SigmaG0W0dc
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
   subroutine printHeader()
      !
      implicit none
      character(len=80)                     :: line
      character(:),allocatable              :: header
      integer                               :: i,Lsx,Ldx,Lcalc
      !
      header="*"
      line="*"
      do i=1,79
         line = trim(line)//"*"
      enddo
      Lcalc=len("Calculation type: "//trim(CalculationType))
      Lsx=int((78-Lcalc)/2)-1
      Ldx=78-Lcalc-Lsx
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
      !
   end subroutine printHeader


   !---------------------------------------------------------------------------!
   !PURPOSE: looks for the current iteration number
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine initialize_DataStructure(ItStart)
      !
      implicit none
      integer,intent(out)                   :: ItStart
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
      select case(CalculationType)
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
      select case(CalculationType)
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
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
               !
            else
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call calc_Gmats(Glat,Crystal)
               !
            endif
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
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
               !
            else
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call calc_Gmats(Glat,Crystal)
               !
            endif
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
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
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
               call read_FermionicField(SigmaDMFT,1,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
               call read_FermionicField(SigmaDMFT,2,reg(pathDATA)//str(ItStart-1)//"/","SigmaDMFT")
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
            if(ItStart.ne.0) merge_Pi = .true.
            calc_Wfull = .true.
            calc_Sigma = .true.
            merge_Sigma = .true.
            !
      end select
      !
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
      calc_W = calc_Wedmft .or. calc_Wfull
      !
   end subroutine initialize_Fields


   !---------------------------------------------------------------------------!
   !PURPOSE: Join the C and X component of the self-energy
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine join_SigmaCX(SigmaFull,Sigma_C,Sigma_X)
      !
      implicit none
      type(FermionicField),intent(inout)    :: SigmaFull
      type(FermionicField),intent(in)       :: Sigma_C
      type(FermionicField),intent(in)       :: Sigma_X
      real(8)                               :: Beta
      integer                               :: Nkpt,Norb,Nmats
      integer                               :: iw,ik,ispin
      !
      !
      if(verbose)write(LOGfile,"(A)") "---- join_SigmaCX"
      !
      !
      !
      ! Check on the input Fields
      if(.not.SigmaFull%status) stop "SigmaFull not properly initialized."
      if(.not.Sigma_C%status) stop "Sigma_C not properly initialized."
      if(.not.Sigma_X%status) stop "Sigma_X not properly initialized."
      if(SigmaFull%Nkpt.eq.0) stop "SigmaFull k dependent attributes not properly initialized."
      if(Sigma_C%Nkpt.eq.0) stop "Sigma_C k dependent attributes not properly initialized."
      if(Sigma_X%Nkpt.eq.0) stop "Sigma_X k dependent attributes not properly initialized."
      if(Sigma_X%Npoints.ne.0) stop "Sigma_X frequency dependent attributes are supposed to be unallocated."
      !
      Norb = SigmaFull%Norb
      Nkpt = SigmaFull%Nkpt
      Beta = SigmaFull%Beta
      Nmats = SigmaFull%Npoints
      !
      if(all([Sigma_C%Nkpt-Nkpt,Sigma_X%Nkpt-Nkpt].ne.[0,0])) stop "Either Sigma_C or Sigma_X have different number of k-points with respect to SigmaFull."
      if(all([Sigma_C%Beta-Beta,Sigma_X%Beta-Beta].ne.[0d0,0d0])) stop "Either Sigma_C or Sigma_X have different Beta with respect to SigmaFull."
      if(Nmats.ne.Sigma_C%Npoints) stop "Sigma_C has different number of Matsubara points with respect to SigmaFull."
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nmats,Nkpt,SigmaFull,Sigma_C,Sigma_X),&
      !$OMP PRIVATE(iw,ik,ispin)
      !$OMP DO
      do iw=1,Nmats
         do ik=1,Nkpt
            do ispin=1,Nspin
               !
               SigmaFull%wks(:,:,iw,ik,ispin) = Sigma_C%wks(:,:,iw,ik,ispin) + Sigma_X%N_ks(:,:,ik,ispin)
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      call FermionicKsum(SigmaFull)
      !
   end subroutine join_SigmaCX
   !
   subroutine join_SigmaCXdc(SigmaFulldc,Sigma_Cdc,Sigma_Xdc)
      !
      implicit none
      type(FermionicField),intent(inout)    :: SigmaFulldc
      type(FermionicField),intent(in)       :: Sigma_Cdc
      type(FermionicField),intent(in)       :: Sigma_Xdc
      real(8)                               :: Beta
      integer                               :: Norb,Nmats
      integer                               :: iw,ispin
      !
      !
      if(verbose)write(LOGfile,"(A)") "---- join_SigmaCXdc"
      !
      !
      !
      ! Check on the input Fields
      if(.not.SigmaFulldc%status) stop "SigmaFulldc not properly initialized."
      if(.not.Sigma_Cdc%status) stop "Sigma_Cdc not properly initialized."
      if(.not.Sigma_Xdc%status) stop "Sigma_Xdc not properly initialized."
      if(SigmaFulldc%Nkpt.eq.0) stop "SigmaFulldc k dependent attributes not properly initialized."
      if(Sigma_Cdc%Nkpt.ne.0) stop "Sigma_Cdc k dependent attributes are supposed to be unallocated."
      if(Sigma_Xdc%Nkpt.ne.0) stop "Sigma_Xdc k dependent attributes are supposed to be unallocated."
      if(Sigma_Xdc%Npoints.ne.0) stop "Sigma_Xdc frequency dependent attributes are supposed to be unallocated."
      !
      Norb = SigmaFulldc%Norb
      Beta = SigmaFulldc%Beta
      Nmats = SigmaFulldc%Npoints
      !
      if(all([Sigma_Cdc%Beta-Beta,Sigma_Xdc%Beta-Beta].ne.[0d0,0d0])) stop "Either Sigma_Cdc or Sigma_Xdc have different Beta with respect to SigmaFulldc."
      if(Nmats.ne.Sigma_Cdc%Npoints) stop "Sigma_Cdc has different number of Matsubara points with respect to SigmaFulldc."
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nmats,SigmaFulldc,Sigma_Cdc,Sigma_Xdc),&
      !$OMP PRIVATE(iw,ispin)
      !$OMP DO
      do iw=1,Nmats
         do ispin=1,Nspin
            !
            SigmaFulldc%ws(:,:,iw,ispin) = Sigma_Cdc%ws(:,:,iw,ispin) + Sigma_Xdc%N_s(:,:,ispin)
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine join_SigmaCXdc




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
