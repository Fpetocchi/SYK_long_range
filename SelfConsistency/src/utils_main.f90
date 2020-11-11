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
   complex(8),allocatable                   :: HlocRotDag(:,:,:)
   real(8),allocatable                      :: HlocEig(:,:)
   !
   type(FermionicField)                     :: Glat
   type(FermionicField)                     :: SigmaFull
   type(FermionicField)                     :: SigmaG0W0,SigmaG0W0dc
   type(FermionicField)                     :: SigmaGW,SigmaGW_C,SigmaGW_X
   type(FermionicField)                     :: SigmaGWdc,SigmaGW_Cdc,SigmaGW_Xdc
   type(FermionicField)                     :: SigmaDMFT
   !
   type(BosonicField)                       :: Wlat
   type(BosonicField)                       :: Ulat
   type(BosonicField)                       :: Plat
   type(BosonicField)                       :: PiEDMFT
   !
   real(8)                                  :: density2set
   complex(8),allocatable                   :: densityLDA(:,:)
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
   logical                                  :: calc_Plat=.false.
   logical                                  :: merge_Pi=.false.
   logical                                  :: calc_W=.false.
   logical                                  :: calc_Wfull=.false.
   logical                                  :: calc_Wedmft=.false.
   logical                                  :: calc_Sigmak=.false.
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
         !
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
         write(LOGfile,"(A)") new_line("A")//new_line("A")//header//new_line("A")
         !
      else
         !
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
         !
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
      integer                               :: iter,Itfirst,isite
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
      call createDir(reg(Itpath),verb=verbose)
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
            Itend = ItStart
            !
      end select
      !
      ItFolder = reg(pathDATA)//str(ItStart)//"/"
      PrevItFolder = reg(pathDATA)//str(ItStart-1)//"/"
      Ustart = Ustart .and. (ItStart.eq.0)
      !
      do isite=1,Nsite
         call createDir(reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/fits",verb=verbose)
         if(ExpandImpurity)exit
      enddo
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
      integer                               :: isite,iorb,Norb,iadd,iset
      integer                               :: iq_gamma_Hk,iq_gamma_XEPS
      complex(8),allocatable                :: Hloc(:,:),Rot(:,:)
      real(8),allocatable                   :: Eig(:)
      integer,allocatable                   :: Orbs(:)
      integer,allocatable                   :: oldSetNorb(:),oldSetOrbs(:,:)
      !
      !
      write(LOGfile,"(A)") new_line("A")//new_line("A")//"---- initialize_Lattice"
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
            allocate(HlocRotDag(SiteNorb(1),SiteNorb(1),Nsite));HlocRotDag=czero
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
               HlocRotDag(1:Norb,1:Norb,isite) = conjg(transpose(Rot))
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
            allocate(HlocRotDag(maxval(SiteNorb),maxval(SiteNorb),Nsite));HlocRotDag=czero
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
               HlocRotDag(1:Norb,1:Norb,isite) = conjg(transpose(Rot))
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
      !Dump some LDA results
      if(ItStart.eq.0)call calc_Glda(0d0,Beta,Lttc)
      !
      !Add the reamining orbitals in the symmetrization list
      if(sum(EqvGWndx%SetNorb).lt.Lttc%Norb)then
         !
         EqvGWndx%Ntotset = EqvGWndx%Nset + (Lttc%Norb-sum(EqvGWndx%SetNorb))
         !
         if(EqvGWndx%Nset.gt.0)then
            !reshape of SetNorb
            allocate(oldSetNorb(size(EqvGWndx%SetNorb)))
            oldSetNorb=EqvGWndx%SetNorb
            deallocate(EqvGWndx%SetNorb)
            allocate(EqvGWndx%SetNorb(EqvGWndx%Ntotset))
            EqvGWndx%SetNorb(1:EqvGWndx%Nset) = oldSetNorb
            EqvGWndx%SetNorb(1+EqvGWndx%Nset:EqvGWndx%Ntotset) = 1
            deallocate(oldSetNorb)
            !
            !reshape of SetOrbs
            allocate(oldSetOrbs(EqvGWndx%Nset,size(EqvGWndx%SetOrbs,dim=2)))
            oldSetOrbs=EqvGWndx%SetOrbs
            deallocate(EqvGWndx%SetOrbs)
            allocate(EqvGWndx%SetOrbs(EqvGWndx%Ntotset,1))
         endif
         !
         iadd=0
         do iorb=1,Lttc%Norb
            do iset=1,EqvGWndx%Nset
               if(allocated(oldSetOrbs).and.any(oldSetOrbs(iset,:).eq.iorb))then
                  cycle
               else
                  iadd=iadd+1
                  EqvGWndx%SetOrbs(EqvGWndx%Nset+iadd,1)=iorb
               endif
            enddo
         enddo
         if(allocated(oldSetOrbs))deallocate(oldSetOrbs)
         !
      else
         EqvGWndx%Ntotset = EqvGWndx%Nset
      endif
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
      integer                               :: unit,isite,iorb
      character(len=255)                    :: filepath
      real(8)                               :: muQMC
      !
      !
      write(LOGfile,"(A)") new_line("A")//new_line("A")//"---- initialize_Fields"
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
            call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.false.,Beta=Beta)
            !
            !Lattice Gf
            call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
            if(ItStart.eq.0) call calc_Gmats(Glat,Crystal)
            if(ItStart.ne.0) call read_FermionicField(Glat,reg(PrevItFolder),"Glat_w",kpt=Crystal%kpt)
            !
            !Logical Flags
            calc_Plat = .true.
            calc_W = .true.
            calc_Wfull = .true.
            calc_Sigmak = .true.
            !
         case("DMFT+statU")
            !
            !Hubbard interaction
            allocate(Umat(Crystal%Norb**2,Crystal%Norb**2));Umat=0d0
            if(Uspex) call read_U_spex(Umat)
            if(Umodel)then
               call inquireFile(reg(pathINPUT)//"Umat_model.DAT",filexists,hardstop=.false.,verb=verbose)
               if(filexists)then
                  call read_Matrix(Umat,reg(pathINPUT)//"Umat_model.DAT")
               else
                  call build_Umat(Umat,Uaa,Uab,J)
                  call dump_Matrix(Umat,reg(pathINPUT)//"Umat_model.DAT")
               endif
            endif
            !
            if(ItStart.ne.0)then
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"Simp_w_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"Simp_w_dw.DAT")
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
               call inquireFile(reg(pathINPUT)//"Uw_model.DAT",filexists,hardstop=.false.,verb=verbose)
               if(filexists)then
                  call read_BosonicField(Ulat,reg(pathINPUT),"Uw_model.DAT")
               else
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph)
                  call dump_BosonicField(Ulat,reg(pathINPUT),"Uw_model.DAT")
               endif
            endif
            !
            if(ItStart.ne.0)then
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"Simp_w_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"Simp_w_dw.DAT")
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
               call inquireFile(reg(pathINPUT)//"Uw_model.DAT",filexists,hardstop=.false.,verb=verbose)
               if(filexists)then
                  call read_BosonicField(Ulat,reg(pathINPUT),"Uw_model.DAT")
               else
                  call build_Uret(Ulat,Uaa,Uab,J,g_eph,wo_eph)
                  call dump_BosonicField(Ulat,reg(pathINPUT),"Uw_model.DAT")
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
               call read_BosonicField(PiEDMFT,reg(PrevItFolder),"Pimp_w")
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"Simp_w_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"Simp_w_dw.DAT")
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
               call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               call AllocateBosonicField(PiEDMFT,Crystal%Norb,Nmats,Crystal%iq_gamma,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               call read_BosonicField(PiEDMFT,reg(PrevItFolder),"Pimp_w.DAT")
               !
               !Impurity Self-energy
               call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"Simp_w_up.DAT")
               call read_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"Simp_w_dw.DAT")
               !
               !Hartree contribution to the Impurity self-energy
               call read_Matrix(SigmaDMFT%N_s(:,:,1),reg(PrevItFolder)//"HartreeU_up.DAT")
               call read_Matrix(SigmaDMFT%N_s(:,:,2),reg(PrevItFolder)//"HartreeU_dw.DAT")
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call read_FermionicField(Glat,reg(PrevItFolder),"Glat_w",kpt=Crystal%kpt)
               !
            else
               !
               !Polarization
               call AllocateBosonicField(Plat,Crystal%Norb,Nmats,Crystal%iq_gamma,Nkpt=Crystal%Nkpt,Nsite=Nsite,no_bare=.true.,Beta=Beta)
               !
               !Lattice Gf
               call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
               call calc_Gmats(Glat,Crystal)
               !
            endif
            !
            !Logical Flags
            calc_Plat = .true.
            calc_Wfull = .true.
            calc_Sigmak = .true.
            if(ItStart.gt.0)then
               merge_Pi = .true.
               merge_Sigma = .true.
            endif
            !
      end select
      !
      calc_W = calc_Wedmft .or. calc_Wfull
      !
      !
      !Allocate and initialize different density matrices
      allocate(densityLDA(Crystal%Norb,Crystal%Norb));densityLDA=czero
      allocate(densityGW(Crystal%Norb,Crystal%Norb,Nspin));densityGW=czero
      allocate(densityDMFT(Crystal%Norb,Crystal%Norb,Nspin));densityDMFT=czero
      allocate(densityQMC(maxval(SiteOrbs),maxval(SiteOrbs),Nspin,Nsite));densityQMC=0d0
      !
      if(ItStart.eq.0)then
         densityLDA = Glat%N_s(:,:,1) + Glat%N_s(:,:,2)
         call dump_Matrix(densityLDA,reg(pathINPUT)//"Nlda.DAT")
         densityGW=czero
         densityDMFT=czero
         densityQMC=0d0
      else
         call read_Matrix(densityLDA,reg(pathINPUT)//"Nlda.DAT")
         densityGW=Glat%N_s
         call read_Matrix(densityDMFT(:,:,1),reg(PrevItFolder)//"Nimp_up.DAT")
         call read_Matrix(densityDMFT(:,:,2),reg(PrevItFolder)//"Nimp_dw.DAT")
         !
         do isite=1,Nsite
            filepath = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/Nqmc.DAT"
            call inquireFile(reg(filepath),filexists,verb=verbose)
            unit = free_unit()
            open(unit,file=reg(filepath),form="formatted",status="old",position="rewind",action="read")
            read(unit,*) muQMC
            do iorb=1,SiteNorb(isite)
               read(unit,*) densityQMC(iorb,iorb,1,isite),densityQMC(iorb,iorb,2,isite)
            enddo
            close(unit)
            if(ExpandImpurity)exit
         enddo
         !
      endif
      !
      if(look4dens%TargetDensity.eq.0d0)then
         Glat%mu = look4dens%mu
      else
         if(ItStart.eq.0)call set_density(Glat%mu,Beta,Crystal,look4dens)
      endif
      !
      write(LOGfile,"(A,F)") new_line("A")//"Lattice chemical potential: ",Glat%mu
      if(ItStart.gt.0)write(LOGfile,"(A,F)") new_line("A")//"Impurity chemical potential: ",muQMC
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
      write(LOGfile,"(A)") new_line("A")//new_line("A")//"---- join_SigmaFull"
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
            !Print full k-dep self-energy: binfmt
            if(printSfull)call dump_FermionicField(SigmaFull,reg(ItFolder),"Sfull_w",.true.,Crystal%kpt)
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
   !PURPOSE: compute, dump and fit the diagonal hybridization function for each
   !         spin-orbital flavor
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine calc_Delta(isite,Iteration)
      !
      implicit none
      integer,intent(in)                    :: isite
      integer,intent(in)                    :: Iteration
      !
      type(FermionicField)                  :: Gloc
      type(FermionicField)                  :: SigmaImp
      type(FermionicField)                  :: FermiPrint
      type(FermionicField)                  :: curlyGold
      integer                               :: Norb,Nflavor,unit
      integer                               :: ispin,iw,iwan,itau,ndx
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: wmats(:),tau(:),Moments(:,:,:)
      real(8),allocatable                   :: Eloc(:,:),PrintLine(:)
      complex(8),allocatable                :: zeta(:,:,:),invGf(:,:),Rot(:,:)
      complex(8),allocatable                :: Dmats(:,:,:),Ditau(:,:,:)
      complex(8),allocatable                :: invCurlyG(:,:,:)
      character(len=255)                    :: printpath,oldparafile,newparafile
      logical                               :: oldparexist
      !
      !
      write(LOGfile,"(A)") new_line("A")//new_line("A")//"---- calc_Delta of "//reg(SiteName(isite))
      !
      !
      Norb = SiteNorb(isite)
      allocate(Orbs(Norb))
      Orbs = SiteOrbs(isite,1:Norb)
      Nflavor = Norb*Nspin
      !
      allocate(Eloc(Norb,Nspin));Eloc=0d0
      !
      call AllocateFermionicField(SigmaImp,Norb,Nmats,Beta=Beta)
      call AllocateFermionicField(Gloc,Norb,Nmats,Beta=Beta)
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
         call loc2imp(Gloc,Glat,Orbs,U=Rot)
         call loc2imp(SigmaImp,SigmaFull,Orbs,U=Rot)
         deallocate(Rot)
      else
         call loc2imp(Gloc,Glat,Orbs)
         call loc2imp(SigmaImp,SigmaFull,Orbs)
      endif
      !
      !Print what's used to compute delta
      call dump_FermionicField(Gloc,1,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","G_"//reg(SiteName(isite))//"_w_up.DAT")
      call dump_FermionicField(Gloc,2,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","G_"//reg(SiteName(isite))//"_w_dw.DAT")
      call dump_FermionicField(SigmaImp,1,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","S_"//reg(SiteName(isite))//"_w_up.DAT")
      call dump_FermionicField(SigmaImp,2,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","S_"//reg(SiteName(isite))//"_w_dw.DAT")
      !
      !Compute the fermionic Weiss field aka the inverse of CurlyG
      allocate(invGf(Norb,Norb));invGf=czero
      do ispin=1,Nspin
         do iw=1,Nmats
            !
            invGf = Gloc%ws(:,:,iw,ispin)
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
      call DeallocateFermionicField(Gloc)
      !
      !Mixing curlyG
      if((Mixing_curlyG.gt.0d0).and.(Iteration.gt.0))then
         call AllocateFermionicField(curlyGold,Norb,Nmats,Beta=Beta)
         call read_FermionicField(curlyGold,1,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_up.DAT")
         call read_FermionicField(curlyGold,2,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_dw.DAT")
         do ispin=1,Nspin
            do iw=1,Nmats
               do iwan=1,Norb
                  invCurlyG(iwan,iw,ispin) = 1d0/((1d0-Mixing_curlyG)*(1d0/invCurlyG(iwan,iw,ispin)) + Mixing_curlyG*curlyGold%ws(iwan,iwan,iw,ispin))
               enddo
            enddo
         enddo
         call DeallocateFermionicField(curlyGold)
      endif
      !
      !Extract the local energy
      select case(reg(CalculationType))
         case default
            !
            stop "If you got so far somethig is wrong."
            !
         case("GW+EDMFT")
            !
            select case(reg(DeltaFit))
               case default
                  !
                  stop "Available modes for Delta fitting: Analytic, Moments."
                  !
               case("Analytic")
                  !
                  oldparafile = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/AndPara_"//reg(SiteName(isite))//".DAT"
                  newparafile = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/AndPara_"//reg(SiteName(isite))//".DAT"
                  call inquireFile(reg(oldparafile),oldparexist,hardstop=.false.,verb=verbose)
                  if(oldparexist) call execute_command_line(" cp "//reg(oldparafile)//" "//reg(newparafile))
                  !
                  call fit_Delta(zeta-invCurlyG,Nbath,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","AndPara_"//reg(SiteName(isite))//".DAT","Shifted",Eloc)
                  !
               case("Moments")
                  !
                  oldparafile = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/Moments_"//reg(SiteName(isite))//".DAT"
                  newparafile = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Moments_"//reg(SiteName(isite))//".DAT"
                  call inquireFile(reg(oldparafile),oldparexist,hardstop=.false.,verb=verbose)
                  if(oldparexist) call execute_command_line(" cp "//reg(oldparafile)//" "//reg(newparafile))
                  !
                  allocate(Moments(Norb,Nbath,Nspin));Moments=0d0
                  call fit_moments(zeta-invCurlyG,Nbath,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","Moments_"//reg(SiteName(isite))//".DAT","Generic",Moments)
                  Eloc=Moments(:,1,:)
                  deallocate(Moments)
                  !
            end select
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
               Dmats(iwan,iw,ispin) = dcmplx( Glat%mu , wmats(iw) ) - Eloc(iwan,ispin) - invCurlyG(iwan,iw,ispin) + EqvGWndx%hseed*(-1d0)**ispin
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
      call createDir(reg(ItFolder)//"Solver_"//reg(SiteName(isite)),verb=verbose)
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
      allocate(PrintLine(Nflavor))
      printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Delta_t.DAT"
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do itau=1,NtauF
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
      !fields that is better to save
      call AllocateFermionicField(FermiPrint,Norb,Nmats,Beta=Beta)
      !
      !Delta(iw)
      do iwan=1,Norb
         FermiPrint%ws(iwan,iwan,:,:) = Dmats(iwan,:,:)
      enddo
      call dump_FermionicField(FermiPrint,1,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","D_"//reg(SiteName(isite))//"_w_up.DAT")
      call dump_FermionicField(FermiPrint,2,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","D_"//reg(SiteName(isite))//"_w_dw.DAT")
      call clear_attributes(FermiPrint)
      !
      !CurlyG(iw)
      do iwan=1,Norb
         FermiPrint%ws(iwan,iwan,:,:) = 1d0/invCurlyG(iwan,:,:)
      enddo
      call dump_FermionicField(FermiPrint,1,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_up.DAT")
      call dump_FermionicField(FermiPrint,2,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_dw.DAT")
      call clear_attributes(FermiPrint)
      !
      !Eo+Delta(iw)
      if((reg(CalculationType).eq."GW+EDMFT").and.verbose)then
         do ispin=1,Nspin
            do iwan=1,Norb
               FermiPrint%ws(iwan,iwan,:,ispin) = zeta(iwan,:,ispin) - invCurlyG(iwan,:,ispin)
            enddo
         enddo
         call dump_FermionicField(FermiPrint,1,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","Eo+D_"//reg(SiteName(isite))//"_w_up.DAT")
         call dump_FermionicField(FermiPrint,2,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","Eo+D_"//reg(SiteName(isite))//"_w_dw.DAT")
      endif
      call DeallocateFermionicField(FermiPrint)
      !
      !test of the FT procedure
      if(verbose)then
         !
         Dmats=czero
         do ispin=1,Nspin
            call Fitau2mats_vec(Beta,Ditau(:,:,ispin),Dmats(:,:,ispin),tau_uniform=.true.)
         enddo
         !
         allocate(PrintLine(Nflavor))
         printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/ReDelta_w_test.DAT"
         unit = free_unit()
         open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nmats
            ndx=1
            PrintLine=0d0
            do iwan=1,Norb
               do ispin=1,Nspin
                  PrintLine(ndx) = real(Dmats(iwan,iw,ispin))
                  ndx=ndx+1
               enddo
            enddo
            write(unit,"(2000E20.12)") wmats(iw),PrintLine
         enddo
         close(unit)
         deallocate(PrintLine)
         !
         allocate(PrintLine(Nflavor))
         printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/ImDelta_w_test.DAT"
         unit = free_unit()
         open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nmats
            ndx=1
            PrintLine=0d0
            do iwan=1,Norb
               do ispin=1,Nspin
                  PrintLine(ndx) = aimag(Dmats(iwan,iw,ispin))
                  ndx=ndx+1
               enddo
            enddo
            write(unit,"(2000E20.12)") wmats(iw),PrintLine
         enddo
         close(unit)
         deallocate(PrintLine)
         !
      endif
      !
      deallocate(Orbs,Eloc,zeta,invCurlyG,Dmats,Ditau,tau,wmats)
      !
   end subroutine calc_Delta


   !---------------------------------------------------------------------------!
   !PURPOSE: compute and dump the istantanous and retarded interactions. Since
   !         if(ExpandImpurity) I have to be sure that all the local interactions
   !         are the same, here I'm extracting the tensor for all the sites and
   !         using the isite input just to print out the wanted one.
   !TEST ON:
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
      integer                               :: Norb,Nflavor,Nbp
      integer                               :: ib1,ib2,itau,iw
      integer                               :: unit,ndx,isitecheck
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: Uinst(:,:),Ucheck(:,:)
      real(8),allocatable                   :: Kfunct(:,:,:)
      real(8),allocatable                   :: tau(:),PrintLine(:)
      character(len=255)                    :: printpath
      !
      !
      write(LOGfile,"(A)") new_line("A")//new_line("A")//"---- calc_Interaction of "//reg(SiteName(isite))
      !
      !
      Norb = SiteNorb(isite)
      allocate(Orbs(Norb))
      Orbs = SiteOrbs(isite,1:Norb)
      Nbp = Norb**2
      Nflavor = Norb*Nspin
      !
      allocate(tau(NtauB));tau=0d0
      tau = linspace(0d0,Beta,NtauB)
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
            call calc_QMCinteractions(curlyU,Uinst)
            !
         case("DMFT+dynU")
            !
            call loc2imp(curlyU,Ulat,Orbs)
            !
            allocate(Kfunct(Nflavor,Nflavor,NtauB));Kfunct=0d0
            call calc_QMCinteractions(curlyU,Uinst,Kfunct)
            !
         case("EDMFT","GW+EDMFT")
            !
            if(Ustart)then
               !
               write(LOGfile,"(A)") "     Using local Ucrpa as effective interaction."
               call loc2imp(curlyU,Ulat,Orbs)
               call DeallocateBosonicField(Ulat)
               !
            else
               !
               write(LOGfile,"(A)") "     Computing the local effective interaction."
               call AllocateBosonicField(Pimp,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
               call AllocateBosonicField(Wloc,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
               call loc2imp(Pimp,Plat,Orbs)
               call loc2imp(Wloc,Wlat,Orbs)
               !
               call calc_curlyU(curlyU,Wloc,Pimp)
               !
               call DeallocateBosonicField(Pimp)
               if(.not.checkInvariance)call DeallocateBosonicField(Wloc)
               !
            endif
            !
            !Mixing curlyU
            if((Mixing_curlyU.gt.0d0).and.(Iteration.gt.0))then
               call AllocateBosonicField(curlyUold,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
               call read_BosonicField(curlyUold,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
               curlyU%bare_local = (1d0-Mixing_curlyU)*curlyU%bare_local + Mixing_curlyU*curlyUold%bare_local
               do iw=1,Nmats
                  curlyU%screened_local(:,:,iw) = (1d0-Mixing_curlyU)*curlyU%screened_local(:,:,iw) + Mixing_curlyU*curlyUold%screened_local(:,:,iw)
               enddo
               call DeallocateBosonicField(curlyU)
            endif
            !
            allocate(Kfunct(Nflavor,Nflavor,NtauB));Kfunct=0d0
            call calc_QMCinteractions(curlyU,Uinst,Kfunct)
            !
      end select
      deallocate(Orbs)
      !
      !Dump Files in the proper directory - The output is printed here as it has to be customized for the solver
      call createDir(reg(ItFolder)//"Solver_"//reg(SiteName(isite)),verb=verbose)
      !
      !Istantaneous interaction
      printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/Umat.DAT"
      call dump_Matrix(Uinst,reg(printpath))
      !
      !Screening function and effective local interaction
      if(allocated(Kfunct))then
         !
         printpath = reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/K_t.DAT"
         unit = free_unit()
         open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
         !
         allocate(PrintLine(Nflavor*(Nflavor+1)/2));PrintLine=0d0
         do itau=1,NtauB
            ndx=1
            !print diagonal and LT
            do ib1=1,Nflavor
               do ib2=1,ib1
                  !
                  PrintLine(ndx) = Kfunct(ib1,ib2,itau)
                  if(Kdiag)then
                     if((mod(ib1,2).eq.0).and.(abs(ib1-ib2).gt.1))PrintLine(ndx) = 0d0
                     if((mod(ib1,2).eq.1).and.(abs(ib1-ib2).ne.0))PrintLine(ndx) = 0d0
                  endif
                  ndx=ndx+1
                  !
               enddo
            enddo
            write(unit,"(999E20.12)") tau(itau),PrintLine
         enddo
         deallocate(PrintLine)
         close(unit)
         !
         call dump_BosonicField(curlyU,reg(ItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
         !
      endif
      !
      deallocate(tau)
      !
      !Check if the fully screened local interaction at iw=0 is the same for all the sites
      if(checkInvariance)then
         !
         do isitecheck=2,Nsite
            !
            Norb = SiteNorb(isitecheck)
            allocate(Orbs(Norb))
            Orbs = SiteOrbs(isitecheck,1:Norb)
            allocate(Ucheck(Nflavor,Nflavor));Ucheck=0d0
            !
            call loc2imp(Wloc,Wlat,Orbs)
            call calc_QMCinteractions(Wloc,Ucheck)
            !
            write(LOGfile,"(A)")"Checking site "//reg(SiteName(1))//"_"//str(isitecheck)
            do ib1=1,Nflavor
               do ib2=1,Nflavor
                  if(abs(Ucheck(ib1,ib2)-Uinst(ib1,ib2)).gt.1e-6)then
                     write(LOGfile,"(A,F,A,F)")"Warning: Element["//str(ib1)//"]["//str(ib2)//"] is different:",Ucheck(ib1,ib2)," instead of: ",Uinst(ib1,ib2)
                  endif
               enddo
            enddo
            !
            deallocate(Orbs,Ucheck)
            !
         enddo
         call DeallocateBosonicField(Wloc)
         !
      endif
      !
      deallocate(Uinst,Kfunct)
      call DeallocateBosonicField(curlyU)
      !
   end subroutine calc_Interaction




   !---------------------------------------------------------------------------!
   !PURPOSE: fetch the impurity fields for each site solved by the QMC solver
   !         and plugs into a container with the same dimension of the lattice
   !         to facilitate the merge. If(ExpandImpurity) the impurity Gf is
   !         rotated back from diagonal to Wannier basis.
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine collect_QMC_results()
      !
      implicit none
      !
      type(FermionicField)                  :: SigmaDMFT
      type(FermionicField)                  :: SigmaImp
      type(FermionicField)                  :: G0imp
      type(BosonicField)                    :: curlyU
      type(BosonicField)                    :: PiEDMFT
      type(BosonicField)                    :: Pimp
      type(BosonicField)                    :: Wimp
      type(BosonicField)                    :: ChiCitau,ChiCmats
      integer                               :: Norb,Nflavor,Nbp
      integer                               :: iorb,jorb,ispin,jspin
      integer                               :: ib1,ib2,isite
      integer                               :: unit,ndx,itau,iw
      real(8)                               :: taup,muQMC
      integer,allocatable                   :: Orbs(:)
      real(8),allocatable                   :: ReadLine(:)
      real(8),allocatable                   :: tauF(:),tauB(:),wmats(:)
      complex(8),allocatable                :: Gitau(:,:,:),Gmats(:,:,:)
      real(8),allocatable                   :: nnt(:,:,:)
      complex(8),allocatable                :: NNitau(:,:,:,:,:),invP(:,:)
      complex(8),allocatable                :: ChiMitau(:),ChiMmats(:)
      logical                               :: filexists
      character(len=255)                    :: filepath
      !
      !
      write(LOGfile,"(A)") "---- collect_QMC_results"
      !
      !
      !
      ! COLLECT IMPURITY OCCUPATION
      allocate(densityDMFT(Crystal%Norb,Crystal%Norb,Nspin));densityDMFT=czero
      allocate(densityQMC(maxval(SiteOrbs),maxval(SiteOrbs),Nspin,Nsite));densityQMC=0d0
      !
      do isite=1,Nsite
         !
         write(LOGfile,"(A)") "     Collecting occupation of site: "//reg(SiteName(isite))
         !
         Norb = SiteNorb(isite)
         allocate(Orbs(Norb))
         Orbs = SiteOrbs(isite,1:Norb)
         Nflavor = Norb*Nspin
         !
         !Read the impurity occupation
         filepath = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/Nqmc.DAT"
         call inquireFile(reg(filepath),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(filepath),form="formatted",status="old",position="rewind",action="read")
         read(unit,*) muQMC
         do iorb=1,Norb
            read(unit,*) densityQMC(iorb,iorb,1,isite),densityQMC(iorb,iorb,2,isite)
         enddo
         close(unit)
         !
         !Insert or Expand to the Lattice basis
         do ispin=1,Nspin
            if(RotateHloc)then
               call imp2loc(densityDMFT(:,:,ispin),dcmplx(densityQMC(:,:,ispin,isite),0d0),Orbs,U=HlocRotDag,expand=ExpandImpurity)
            else
               call imp2loc(densityDMFT(:,:,ispin),dcmplx(densityQMC(:,:,ispin,isite),0d0),Orbs,expand=ExpandImpurity)
            endif
         enddo
         !
         if(ExpandImpurity)exit
         !
      enddo
      !
      !Symmetrize
      if(verbose)call dump_Matrix(densityDMFT(:,:,1),reg(PrevItFolder)//"Nimp_up_notsymm.DAT")
      call symmetrize(densityDMFT,EqvGWndx)
      !
      !Save to the proper iteration folder
      call dump_Matrix(densityDMFT(:,:,1),reg(PrevItFolder)//"Nimp_up.DAT")
      call dump_Matrix(densityDMFT(:,:,2),reg(PrevItFolder)//"Nimp_dw.DAT")
      deallocate(densityDMFT)
      !
      !
      !
      !
      ! COLLECT IMPURITY GF AND FERMIONIC DYSON
      call AllocateFermionicField(SigmaDMFT,Crystal%Norb,Nmats,Nsite=Nsite,Beta=Beta)
      do isite=1,Nsite
         !
         write(LOGfile,"(A)") "     Solving fermionic Dyson of site: "//reg(SiteName(isite))
         !
         Norb = SiteNorb(isite)
         allocate(Orbs(Norb))
         Orbs = SiteOrbs(isite,1:Norb)
         !
         !Read the impurity Green's function
         allocate(tauF(NtauF));tauF=0d0
         tauF = linspace(0d0,Beta,NtauF)
         allocate(Gitau(Norb,NtauF,Nspin));Gitau=czero
         allocate(ReadLine(Nspin*Norb))
         filepath = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/Gimp_t.DAT"
         call inquireFile(reg(filepath),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(filepath),form="formatted",status="old",position="rewind",action="read")
         do itau=1,NtauF
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
         allocate(Gmats(Norb,Nmats,Nspin));Gmats=czero
         do ispin=1,Nspin
            call Fitau2mats_vec(Beta,Gitau(:,:,ispin),Gmats(:,:,ispin),tau_uniform=.true.)
         enddo
         deallocate(Gitau)
         !
         !Read curlyG
         call AllocateFermionicField(G0imp,Norb,Nmats,Beta=Beta)
         call read_FermionicField(G0imp,1,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_up.DAT")
         call read_FermionicField(G0imp,2,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_dw.DAT")
         !
         !Adjust with the chemical potential used by the solver
         if(G0imp%mu.ne.muQMC)then
            write(LOGfile,"(A)") "     Updating the chemical potential of curlyG."
            do ispin=1,Nspin
               do iw=1,Nmats
                  do iorb=1,Norb
                     G0imp%ws(iorb,iorb,iw,ispin) = 1d0/(1d0/G0imp%ws(iorb,iorb,iw,ispin) - G0imp%mu + muQMC)
                  enddo
               enddo
            enddo
         endif
         call dump_FermionicField(G0imp,1,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_up.DAT")
         call dump_FermionicField(G0imp,2,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","G0_"//reg(SiteName(isite))//"_w_dw.DAT")
         !
         !Fermionic Dyson equation in the solver basis (always diagonal)
         call AllocateFermionicField(SigmaImp,Norb,Nmats,Beta=Beta)
         do ispin=1,Nspin
            do iorb=1,Norb
               SigmaImp%ws(iorb,iorb,:,ispin) = 1d0/G0imp%ws(iorb,iorb,:,ispin) - 1d0/Gmats(iorb,:,ispin)
            enddo
         enddo
         deallocate(Gmats)
         call DeallocateFermionicField(G0imp)
         !
         !Fill up the N_s attribute that correspond to the Hartree term
         !since the impurity interaction contains only the physical interaction elements
         !(particle and spin conserving) I'm neglecting here elements in principle present
         !like Hartree_{ab} = Sum_c curlyU_{ab}{cc} * n_{cc}
         call AllocateBosonicField(curlyU,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
         call read_BosonicField(curlyU,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
         do ispin=1,Nspin
            do iorb=1,Norb
               do jorb=1,Norb
                  !
                  ib1 = iorb + Norb*(iorb-1)
                  ib2 = jorb + Norb*(jorb-1)
                  !
                  SigmaImp%N_s(iorb,iorb,ispin) = SigmaImp%N_s(iorb,iorb,ispin) + curlyU%screened_local(ib1,ib2,1)*densityQMC(jorb,jorb,ispin,isite)
                  !
               enddo
            enddo
         enddo
         call DeallocateBosonicField(curlyU)
         !
         !Expand to the Lattice basis
         if(RotateHloc)then
            call imp2loc(SigmaDMFT,SigmaImp,Orbs,U=HlocRotDag,expand=ExpandImpurity,AFM=AFMselfcons)
         else
            call imp2loc(SigmaDMFT,SigmaImp,Orbs,expand=ExpandImpurity,AFM=AFMselfcons)
         endif
         call DeallocateFermionicField(SigmaImp)
         deallocate(Orbs)
         !
         if(ExpandImpurity)exit
         !
      enddo
      !
      !Symmetrize
      if(verbose)call dump_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"Simp_w_up_notsymm.DAT")
      call symmetrize(SigmaDMFT,EqvGWndx)
      !
      !Save to the proper iteration folder
      call dump_FermionicField(SigmaDMFT,1,reg(PrevItFolder),"Simp_w_up.DAT")
      call dump_FermionicField(SigmaDMFT,2,reg(PrevItFolder),"Simp_w_dw.DAT")
      call dump_Matrix(SigmaDMFT%N_s(:,:,1),reg(PrevItFolder)//"HartreeU_up.DAT")
      call dump_Matrix(SigmaDMFT%N_s(:,:,2),reg(PrevItFolder)//"HartreeU_dw.DAT")
      call DeallocateFermionicField(SigmaDMFT)
      !
      !
      !
      !
      ! COLLECT IMPURITY CHARGE SUSCEPTIBILITY AND BOSONIC DYSON
      if(bosonicSC)then
         !
         call AllocateBosonicField(PiEDMFT,Crystal%Norb,Crystal%iq_gamma,Nmats,Nsite=Nsite,no_bare=.true.,Beta=Beta)
         !
         do isite=1,Nsite
            !
            write(LOGfile,"(A)") "     Solving bosonic Dyson of site: "//reg(SiteName(isite))
            !
            Norb = SiteNorb(isite)
            allocate(Orbs(Norb))
            Orbs = SiteOrbs(isite,1:Norb)
            Nbp = Norb**2
            !
            !Read the impurity N(tau)N(0)
            allocate(tauB(NtauB));tauB=0d0
            tauB = linspace(0d0,Beta,NtauB)
            allocate(nnt(Nflavor,Nflavor,NtauB));nnt=0d0
            allocate(ReadLine(Nflavor*(Nflavor+1)/2))
            filepath = reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/resultsQMC/nn_t.DAT"
            call inquireFile(reg(filepath),filexists,verb=verbose)
            unit = free_unit()
            open(unit,file=reg(filepath),form="formatted",status="old",position="rewind",action="read")
            do itau=1,NtauB
               ReadLine=0d0
               read(unit,*) taup,ReadLine
               if(abs(taup-tauB(itau)).gt.eps) stop "Impurity bosonic tau mesh does not coincide."
               !
               !this is cumbersome but better be safe
               ndx=1
               do ib1=1,Nflavor
                  do ib2=1,ib1
                     nnt(ib1,ib2,itau) = ReadLine(ndx)
                     nnt(ib2,ib1,itau) = ReadLine(ndx)
                     ndx=ndx+1
                  enddo
               enddo
            enddo
            deallocate(ReadLine)
            !
            !Reshape N(tau)N(0)for easier handling
            allocate(NNitau(Norb,Norb,Nspin,Nspin,NtauB));NNitau=czero
            do ib1=1,Nspin*Norb
               do ib2=1,Nspin*Norb
                  !
                  iorb = (ib1+mod(ib1,2))/2
                  jorb = (ib2+mod(ib2,2))/2
                  ispin = abs(mod(ib1,2)-2)
                  jspin = abs(mod(ib2,2)-2)
                  !
                  NNitau(iorb,jorb,ispin,jspin,:) = dcmplx(nnt(ib1,ib2,:),0d0)
                  !
               enddo
            enddo
            deallocate(nnt)
            !
            !Compute the spin susceptibility Sum_ab <S(tau)S(0)>
            allocate(ChiMitau(NtauB));ChiMitau=czero
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        !
                        ib1 = iorb + Norb*(iorb-1)
                        ib2 = jorb + Norb*(jorb-1)
                        !
                        ChiMitau = ChiMitau + NNitau(iorb,jorb,ispin,jspin,:)*(-1d0)**(ispin-jspin)
                        !
                     enddo
                  enddo
               enddo
            enddo
            allocate(ChiMmats(Nmats));ChiMmats=czero
            call Bitau2mats(Beta,ChiMitau,ChiMmats,tau_uniform=.true.)
            call dump_BosonicField(ChiMitau,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiM_"//reg(SiteName(isite))//"_t.DAT",tauB)
            allocate(wmats(Nmats));wmats=BosonicFreqMesh(Beta,Nmats)
            call dump_BosonicField(ChiMmats,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiM_"//reg(SiteName(isite))//"_w.DAT",wmats)
            deallocate(ChiMitau,ChiMmats,wmats)
            !
            !Compute the charge susceptibility Sum_s1s2 <Na(tau)Nb(0)> - <Na><Nb>
            call AllocateBosonicField(ChiCitau,Norb,NtauB,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        !
                        ib1 = iorb + Norb*(iorb-1)
                        ib2 = jorb + Norb*(jorb-1)
                        !
                        ChiCitau%screened_local(ib1,ib2,:) = ChiCitau%screened_local(ib1,ib2,:) + &
                        NNitau(iorb,jorb,ispin,jspin,:) - densityQMC(iorb,iorb,ispin,isite)*densityQMC(jorb,jorb,jspin,isite)
                        !
                     enddo
                  enddo
               enddo
            enddo
            call AllocateBosonicField(ChiCmats,Norb,Nmats,Crystal%iq_gamma,no_bare=.true.,Beta=Beta)
            call Bitau2mats(Beta,ChiCitau%screened_local,ChiCmats%screened_local,tau_uniform=.true.)
            call dump_BosonicField(ChiCitau,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiC_"//reg(SiteName(isite))//"_t.DAT",axis=tauB)
            call dump_BosonicField(ChiCmats,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","ChiC_"//reg(SiteName(isite))//"_w.DAT")
            call DeallocateBosonicField(ChiCitau)
            deallocate(tauB,densityQMC)
            !
            !Bosonic Dyson equation in the solver basis
            call AllocateBosonicField(curlyU,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
            call read_BosonicField(curlyU,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","curlyU_"//reg(SiteName(isite))//"_w.DAT")
            call AllocateBosonicField(Pimp,Norb,Crystal%iq_gamma,Nmats,no_bare=.true.,Beta=Beta)
            allocate(invP(Nbp,Nbp));invP=czero
            do iw=1,Nmats
               !
               invP = matmul(curlyU%screened_local(:,:,iw),ChiCmats%screened_local(:,:,iw)) - zeye(Nbp)
               !
               call inv(invP)
               !
               Pimp%screened_local(:,:,iw) = matmul(ChiCmats%screened_local(:,:,iw),invP)
               !
            enddo
            deallocate(invP)
            !
            !Compute Wimp to check convergence on the interaction
            call AllocateBosonicField(Wimp,Norb,Nmats,Crystal%iq_gamma,Beta=Beta)
            call calc_Wimp(Wimp,curlyU,ChiCmats)
            call dump_BosonicField(Wimp,reg(PrevItFolder)//"Solver_"//reg(SiteName(isite))//"/","Wimp_"//reg(SiteName(isite))//"_w.DAT")
            !
            call DeallocateBosonicField(Wimp)
            call DeallocateBosonicField(curlyU)
            call DeallocateBosonicField(ChiCmats)
            !
            !Expand to the Lattice basis
            call imp2loc(PiEDMFT,Pimp,Orbs,expand=ExpandImpurity)
            call DeallocateBosonicField(Pimp)
            deallocate(Orbs)
            !
            if(ExpandImpurity)exit
            !
         enddo
         !
         !Symmetrize
         if(verbose)call dump_BosonicField(PiEDMFT,reg(PrevItFolder),"Pimp_w_notsymm")
         call symmetrize(PiEDMFT,EqvGWndx)
         !
         !Save to the proper iteration folder
         call dump_BosonicField(PiEDMFT,reg(PrevItFolder),"Pimp_w")
         call DeallocateBosonicField(PiEDMFT)
         !
      endif
      !
   end subroutine collect_QMC_results


   !---------------------------------------------------------------------------!
   !PURPOSE: Deallocate all fields
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine DeallocateAllFields()
      implicit none
      if(Glat%status) call DeallocateFermionicField(Glat)
      if(SigmaFull%status) call DeallocateFermionicField(SigmaFull)
      if(SigmaG0W0%status) call DeallocateFermionicField(SigmaG0W0)
      if(SigmaG0W0dc%status) call DeallocateFermionicField(SigmaG0W0dc)
      if(SigmaGW%status) call DeallocateFermionicField(SigmaGW)
      if(SigmaGW_C%status) call DeallocateFermionicField(SigmaGW_C)
      if(SigmaGW_X%status) call DeallocateFermionicField(SigmaGW_X)
      if(SigmaGWdc%status) call DeallocateFermionicField(SigmaGWdc)
      if(SigmaGW_Cdc%status) call DeallocateFermionicField(SigmaGW_Cdc)
      if(SigmaGW_Xdc%status) call DeallocateFermionicField(SigmaGW_Xdc)
      if(SigmaDMFT%status) call DeallocateFermionicField(SigmaDMFT)
      if(Wlat%status) call DeallocateBosonicField(Wlat)
      if(Ulat%status) call DeallocateBosonicField(Ulat)
      if(Plat%status) call DeallocateBosonicField(Plat)
      if(PiEDMFT%status) call DeallocateBosonicField(PiEDMFT)
   end subroutine DeallocateAllFields


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
