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
   complex(8),allocatable                   :: HlocSite(:,:,:)
   complex(8),allocatable                   :: HlocRot(:,:,:)
   real(8),allocatable                      :: HlocEig(:,:)
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
   !
   character(len=255)                       :: ItFolder,PrevItFolder

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
         !
         ItFolder = reg(pathDATA)//str(Iteration)//"/"
         PrevItFolder = reg(pathDATA)//str(Iteration-1)//"/"
         !
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
      integer                               :: m,n,mp,np,ib1,ib2,isite,Norb
      integer                               :: iq_gamma_Hk,iq_gamma_XEPS
      complex(8),allocatable                :: Hloc(:,:),Rot(:,:)
      real(8),allocatable                   :: Eig(:)
      integer,allocatable                   :: Orbs(:)
      !
      !
      write(LOGfile,"(A)") "---- initialize_Lattice"
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
      !Store the local rotation of each site
      if(RotateHloc)then
         !
         if(ExpandImpurity)then !Only one set of orbital provided - all matrices with same dimension
            !
            allocate(HlocEig(SiteNorb(1),Nsite));HlocEig=0d0
            allocate(HlocSite(SiteNorb(1),SiteNorb(1),Nsite));HlocSite=czero
            allocate(HlocRot(SiteNorb(1),SiteNorb(1),Nsite));HlocRot=czero
            !
            do isite=1,Nsite
               !
               Norb = SiteNorb(1)
               !
               allocate(Eig(Norb));Eig=0d0
               allocate(Hloc(Norb,Norb));Hloc=czero
               allocate(Rot(Norb,Norb));Rot=czero
               allocate(Orbs(Norb));Orbs=0
               !
               ! only two possible arrangements
               if(abs(SiteOrbs(1,2)-SiteOrbs(1,1)).eq.1)then
                  Orbs = SiteOrbs(1,:) + Norb*(isite-1)
               elseif(abs(SiteOrbs(1,2)-SiteOrbs(1,1)).eq.Nsite)then
                  Orbs = SiteOrbs(1,:) + isite-1
               endif
               !
               !Extract then local Hamiltonian for each site
               call loc2imp(Hloc,Lttc%Hloc,Orbs)
               Rot = Hloc
               !
               !Rotate
               call eigh(Rot,Eig)
               !
               !Save
               call dump_Matrix(Hloc,reg(pathINPUT)//"HlocSite_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               call dump_Matrix(Rot,reg(pathINPUT)//"HlocRot_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               call dump_Matrix(diag(Eig),reg(pathINPUT)//"HlocEig_"//reg(SiteName(1))//"_"//str(isite)//".DAT")
               !
               !SO(3)check
               write(LOGfile,"(A,F)") "     det(Rot) of "//reg(SiteName(isite))//" :",det(HlocRot(:,:,isite))
               !
               HlocSite(1:Norb,1:Norb,isite) = Hloc
               HlocRot(1:Norb,1:Norb,isite) = Rot
               HlocEig(1:Norb,isite) = Eig
               !
               deallocate(Eig,Hloc,Rot,Orbs)
               !
            enddo
            !
         else !Potentially different set of orbitals provided - all matrices with different dimension (SiteNorb(isite))
            !
            allocate(HlocEig(maxval(SiteNorb),Nsite));HlocEig=0d0
            allocate(HlocSite(maxval(SiteNorb),maxval(SiteNorb),Nsite));HlocSite=czero
            allocate(HlocRot(maxval(SiteNorb),maxval(SiteNorb),Nsite));HlocRot=czero
            !
            do isite=1,Nsite
               !
               Norb = SiteNorb(isite)
               !
               allocate(Eig(Norb));Eig=0d0
               allocate(Hloc(Norb,Norb));Hloc=czero
               allocate(Rot(Norb,Norb));Rot=czero
               allocate(Orbs(Norb));Orbs=0
               !
               Orbs=SiteOrbs(isite,:)
               !
               !Extract then local Hamiltonian for each site
               call loc2imp(Hloc,Lttc%Hloc,Orbs)
               Rot = Hloc
               !
               !Rotate
               call eigh(Rot,Eig)
               !
               !Save
               call dump_Matrix(Hloc,reg(pathINPUT)//"HlocSite_"//reg(SiteName(isite))//".DAT")
               call dump_Matrix(Rot,reg(pathINPUT)//"HlocRot_"//reg(SiteName(isite))//".DAT")
               call dump_Matrix(diag(Eig),reg(pathINPUT)//"HlocEig_"//reg(SiteName(isite))//".DAT")
               !
               !SO(3)check
               write(LOGfile,"(A,F)") "     det(Rot) of "//reg(SiteName(isite))//" :",det(HlocRot(:,:,isite))
               !
               HlocSite(1:Norb,1:Norb,isite) = Hloc
               HlocRot(1:Norb,1:Norb,isite) = Rot
               HlocEig(1:Norb,isite) = Eig
               !
               deallocate(Eig,Hloc,Rot,Orbs)
               !
            enddo
            !
         endif
         !
      endif
      !
      !Store the Physical (number and spin conserving) elements of the interaction tensor
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
      !Dump some LDA results
      if(ItStart.eq.0)call calc_Glda(0d0,Beta,Lttc)
      !
      !Allocate different density matrices
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
      write(LOGfile,"(A)") "---- initialize_Fields"
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
            if(ItStart.ne.0) call read_FermionicField(Glat,reg(PrevItFolder),"Glat",kpt=Crystal%kpt)
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
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"SigmaDMFT_dw.DAT")
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
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"SigmaDMFT_dw.DAT")
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
               call read_BosonicField(PiEDMFT,reg(PrevItFolder),"PiEDMFT")
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"SigmaDMFT_dw.DAT")
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
               call read_BosonicField(PiEDMFT,reg(PrevItFolder),"PiEDMFT.DAT")
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"SigmaDMFT_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"SigmaDMFT_dw.DAT")
               !
               !Hartree contribution to the Impurity self-energy
               call read_Matrix(SigmaDMFT%N_s(:,:,1),reg(PrevItFolder)//"HartreeU_up.DAT")
               call read_Matrix(SigmaDMFT%N_s(:,:,2),reg(PrevItFolder)//"HartreeU_dw.DAT")
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(Glat,reg(PrevItFolder),"Glat",kpt=Crystal%kpt)
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
   subroutine join_SigmaFull(Iteration)
      !
      implicit none
      integer,intent(in)                    :: Iteration
      integer                               :: iorb,jorb,ik,iw,ispin
      !
      !
      write(LOGfile,"(A)") "---- join_SigmaFull"
      !
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
               if(.not.SigmaG0W0%status) stop "join_SigmaFull: SigmaG0W0 not properly initialized."
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
               if(.not.SigmaGW%status) stop "join_SigmaFull: SigmaGW not properly initialized."
               if(.not.SigmaG0W0dc%status) stop "join_SigmaFull: SigmaG0W0dc not properly initialized."
               if(.not.SigmaG0W0%status) stop "join_SigmaFull: SigmaG0W0 not properly initialized."
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
            call dump_FermionicField(SigmaFull,1,reg(ItFolder),"SigmaLoc_up.DAT")
            call dump_FermionicField(SigmaFull,2,reg(ItFolder),"SigmaLoc_dw.DAT")
            !
            !Dump the whole SigmaFull
            call dump_FermionicField(SigmaFull,reg(ItFolder),"Sigma",.true.,Crystal%kpt)
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
               if(.not.SigmaDMFT%status) stop "join_SigmaFull: SigmaDMFT not properly initialized."
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
   end subroutine join_SigmaFull


   !---------------------------------------------------------------------------!
   !PURPOSE: Join the all the component of the self-energy
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine calc_Delta(isite)
      !
      implicit none
      integer,intent(in)                    :: isite
      !
      type(FermionicField)                  :: Gimp
      type(FermionicField)                  :: SigmaImp
      integer                               :: Norb,unit
      integer                               :: ispin,iw,iwan,itau,ndx
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: Eloc(:,:),PrintLine(:)
      complex(8),allocatable                :: zeta(:,:,:),invGf(:,:),Rot(:,:)
      complex(8),allocatable                :: Dmats(:,:,:),Ditau(:,:,:)
      complex(8),allocatable                :: invCurlyG(:,:,:)
      character(len=255)                    :: printpath,oldparafile,newparafile
      logical                               :: oldparexist
      !
      !
      write(LOGfile,"(A)") "---- calc_Delta of site "//str(isite)
      !
      !
      Norb = SiteNorb(isite)
      allocate(Orbs(Norb))
      Orbs = SiteOrbs(isite,1:Norb)
      !
      allocate(Eloc(Norb,Nspin));Eloc=0d0
      !
      call AllocateFermionicField(SigmaImp,Norb,Nmats,Beta=Beta)
      call AllocateFermionicField(Gimp,Norb,Nmats,Beta=Beta)
      allocate(invCurlyG(Norb,Nmats,Nspin));invCurlyG=czero
      allocate(Dmats(Norb,Nmats,Nspin));Dmats=czero
      allocate(Ditau(Norb,NtauF,Nspin));Ditau=0d0
      !
      allocate(tau(NtauF));tau=0d0
      tau = linspace(0d0,Beta,NtauF)
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(Beta,Nmats)
      allocate(zeta(Norb,Nmats,Nspin))
      do iwan=1,Norb
         do iw=1,Nmats
            zeta(iwan,iw,:) = dcmplx( Glat%mu , wmats(iw) )
         enddo
      enddo
      !
      ! Extract from local to imp the given sites
      if(RotateHloc)then
         allocate(Rot(Norb,Norb)); Rot=HlocRot(1:Norb,1:Norb,isite)
         call loc2imp(Gimp,Glat,Orbs,U=Rot)
         call loc2imp(SigmaImp,SigmaFull,Orbs,U=Rot)
         deallocate(Rot)
      else
         call loc2imp(Gimp,Glat,Orbs)
         call loc2imp(SigmaImp,SigmaFull,Orbs)
      endif
      !
      !Print if needed
      if(Nsite.gt.1)then
         if(RotateHloc)then
            call dump_FermionicField(Gimp,1,reg(ItFolder),"Gloc_up_diag_Site"//str(isite)//".DAT")
            call dump_FermionicField(Gimp,2,reg(ItFolder),"Gloc_dn_diag_Site"//str(isite)//".DAT")
            call dump_FermionicField(SigmaImp,1,reg(ItFolder),"SigmaLoc_up_diag_Site"//str(isite)//".DAT")
            call dump_FermionicField(SigmaImp,2,reg(ItFolder),"SigmaLoc_dw_diag_Site"//str(isite)//".DAT")
         else
            !THIS IS JUST A PEDANTIC TEST - IF IT WORKS IT CAN BE REMOVED
            call dump_FermionicField(Gimp,1,reg(ItFolder),"Gloc_up_Site"//str(isite)//".DAT")
            call dump_FermionicField(Gimp,2,reg(ItFolder),"Gloc_dn_Site"//str(isite)//".DAT")
            call dump_FermionicField(SigmaImp,1,reg(ItFolder),"SigmaLoc_up_Site"//str(isite)//".DAT")
            call dump_FermionicField(SigmaImp,2,reg(ItFolder),"SigmaLoc_dw_Site"//str(isite)//".DAT")
         endif
      elseif((Nsite.eq.1).and.RotateHloc)then
         call dump_FermionicField(Gimp,1,reg(ItFolder),"Gloc_up_diag.DAT")
         call dump_FermionicField(Gimp,2,reg(ItFolder),"Gloc_dn_diag.DAT")
         call dump_FermionicField(SigmaImp,1,reg(ItFolder),"SigmaLoc_up_diag.DAT")
         call dump_FermionicField(SigmaImp,2,reg(ItFolder),"SigmaLoc_dw_diag.DAT")
      endif
      !
      !Compute the fermionic Weiss field aka the inverse of CurlyG
      allocate(invGf(Norb,Norb));invGf=czero
      do ispin=1,Nspin
         do iw=1,Nmats
            !
            invGf = Gimp%ws(:,:,iw,ispin)
            call inv(invGf)
            !
            do iwan=1,Norb
               invCurlyG(iwan,iw,ispin) = invGf(iwan,iwan) + SigmaImp%ws(iwan,iwan,iw,ispin)
            enddo
            !
         enddo
      enddo
      deallocate(invGf)
      call DeallocateFermionicField(SigmaImp)
      call DeallocateFermionicField(Gimp)
      !
      !Extract the local energy
      select case(reg(CalculationType))
         case default
            !
            stop "If you got so fare somethig is wrong."
            !
         case("GW+EDMFT")
            !
            oldparafile = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/AndPara.DAT"
            newparafile = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/AndPara.DAT"
            call inquireFile(reg(oldparafile),oldparexist,hardstop=.false.)
            if(oldparexist) call execute_command_line(" cp "//reg(oldparafile)//" "//reg(newparafile))
            call fit_Delta(zeta-invCurlyG,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","AndPara.DAT","Shifted",Eloc)
            !
         case("DMFT+statU","DMFT+dynU","EDMFT")
            !
            if(RotateHloc)then
               Eloc(:,1) = HlocEig(1:Norb,isite)
               Eloc(:,2) = HlocEig(1:Norb,isite)
            else
               Eloc(:,1) = diagonal(HlocSite(1:Norb,1:Norb,isite))
               Eloc(:,2) = diagonal(HlocSite(1:Norb,1:Norb,isite))
            endif
            !
      end select
      !
      !Compute Delta on matsubara
      do ispin=1,Nspin
         do iw=1,Nmats
            do iwan=1,Norb
               !
               Dmats(iwan,iw,ispin) = dcmplx( Glat%mu , wmats(iw) ) - Eloc(iwan,ispin) - invCurlyG(iwan,iw,ispin)
               !
            enddo
         enddo
      enddo
      !
      !Fourier transform to the tau axis
      do ispin=1,Nspin
         call Fmats2itau_vec(Beta,Dmats(:,:,ispin),Ditau(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
      enddo
      !
      !Dump Files in the proper directory - The output is printed here as it has to be customized for the solver
      allocate(PrintLine(1+Norb*Nspin))
      call createDir(reg(ItFolder)//"Solver_"//reg(SiteName(isite)))
      !
      !Eloc and chemical potential
      printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Eloc.DAT"
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(1E20.12)") Glat%mu
      do iwan=1,Norb
         write(unit,"(2E20.12)") Eloc(iwan,1),Eloc(iwan,2)
      enddo
      close(unit)
      !
      !Delta(tau)
      printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Delta.DAT"
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do itau=1,NtauF
         ndx=1
         PrintLine=0d0
         PrintLine(ndx) = tau(itau)
         do iwan=1,Norb
            do ispin=1,Nspin
               ndx=ndx+1
               PrintLine(ndx) = real(Ditau(iwan,itau,ispin))
            enddo
         enddo
         write(unit,"(2000E20.12)") PrintLine
      enddo
      close(unit)
      !
      !Delta(iw)
      printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Delta_iw.DAT"
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do iw=1,Nmats
         ndx=1
         PrintLine=0d0
         PrintLine(ndx) = wmats(iw)
         do iwan=1,Norb
            do ispin=1,Nspin
               ndx=ndx+1
               PrintLine(ndx) = real(Dmats(iwan,iw,ispin))
            enddo
         enddo
         write(unit,"(2000E20.12)") PrintLine
      enddo
      close(unit)
      !
      !CurlyG(iw)
      printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/CurlyG_iw.DAT"
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do iw=1,Nmats
         ndx=1
         PrintLine=0d0
         PrintLine(ndx) = wmats(iw)
         do iwan=1,Norb
            do ispin=1,Nspin
               ndx=ndx+1
               PrintLine(ndx) = 1d0/real(invCurlyG(iwan,iw,ispin))
            enddo
         enddo
         write(unit,"(2000E20.12)") PrintLine
      enddo
      close(unit)
      !
      !Eloc _ Delta(iw)
      if(verbose)then
         printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Eloc_Delta_iw.DAT"
         unit = free_unit()
         open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nmats
            ndx=1
            PrintLine=0d0
            PrintLine(ndx) = wmats(iw)
            do iwan=1,Norb
               do ispin=1,Nspin
                  ndx=ndx+1
                  PrintLine(ndx) = dcmplx( Glat%mu , wmats(iw) ) - invCurlyG(iwan,iw,ispin)
               enddo
            enddo
            write(unit,"(2000E20.12)") PrintLine
         enddo
         close(unit)
      endif
      !
      deallocate(Orbs,Eloc,zeta,invCurlyG,Dmats,Ditau,tau,wmats,PrintLine)
      !
   end subroutine calc_Delta












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
