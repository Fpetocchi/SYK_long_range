module crystal

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !metto i vettori in spazio reale qui perche sevono anche a routine interne e
   !in genere mai al codice che lavora in spazio K.
   !cmq possono essere accessibili e sbattuti nell"attributo lattice con suroutine dedicate
   !il calc wigsiez e smallk hanno bisogno che LATTC sia letto
   ! qui mantengo il path allinputfile dato dall'utente perche cmq se leggo la Hk la devo leggere da qulahce parte
   ! per funzionare wannier intp devo leggere Hk perche' e' l'unico modo per dare il path

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface wannierinterpolation
      module procedure wannierinterpolation_mat                                 !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(Norb,Norb,Ndat,Nkpt_orig),mat_intp(Norb,Norb,Ndat,Nkpt_intp))
      module procedure wannierinterpolation_vec                                 !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(Norb,Ndat,Nkpt_orig),mat_intp(Norb,Ndat,Nkpt_intp))
   end interface wannierinterpolation

   interface wannier_K2R
      module procedure wannier_K2R_mat                                          !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_K(Norb,Norb,Ndat,Nkpt_orig),mat_R(internally allocated))
      module procedure wannier_K2R_vec                                          !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_K(Norb,Ndat,Nkpt_orig),mat_R(internally allocated))
   end interface wannier_K2R

   interface wannier_R2K
      module procedure wannier_R2K_mat                                          !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_R(Norb,Norb,Ndat,Nwig),mat_K(Norb,Norb,Ndat,Nkpt_orig))
      module procedure wannier_R2K_vec                                          !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_R(Norb,Ndat,Nwig),mat_K(Norb,Ndat,Nkpt_orig))
   end interface wannier_R2K

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),parameter,private                :: H2eV=27.2113831243217d0
   real(8),parameter,private                :: pi=3.14159265358979323846d0
   complex(8),parameter,private             :: czero=dcmplx(0.d0,0.d0)
   !
   real(8),parameter,private                :: eps=1e-9
   real(8),parameter,private                :: epsWig=1e-5
   !
   real(8),allocatable                      :: rsite(:,:)
   real(8)                                  :: rlat(3,3)
   integer,allocatable                      :: Kprint(:)
   !
   real(8),private                          :: lat(3,3)
   real(8),private                          :: vol
   real(8),private                          :: rvol
   !
   integer,private                          :: Nwig
   integer,allocatable,private              :: rvecwig(:,:)
   integer,allocatable,private              :: nrdegwig(:)
   !
   logical                                  :: Hk_stored=.false.
   logical                                  :: Ruc_stored=.false.               !Global flag for routines that need positions within the u.c.
   logical,private                          :: Lat_stored=.false.               !Internal flag for routines that need rlat
   logical,private                          :: Wig_stored=.false.               !Internal flag for routines performing Wannier interpolation

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !variables
   public :: Hk_stored
   public :: Ruc_stored
   public :: rlat
   !subroutines
   public :: read_lattice
   public :: read_xeps
   public :: read_Hk
   public :: fill_ksumkdiff
   public :: fill_smallk
   public :: set_siteposition
   public :: wannierinterpolation                                               ![Nkpt3_orig(3),Kpt_orig(3,Nkpt_orig),Kpt_intp(3,Nkpt_intp),Mat_orig(n,n,Npoins,Nkpt_orig),Mat_intp(n,n,Npoins,Nkpt_intp)]
   public :: wannier_K2R_NN                                                     ![Nkpt3_orig(3),Kpt_orig(3,Nkpt_orig),Mat_orig(n,n,Npoins,Nkpt_orig),mat_R_nn(n,n,Npoins,3)]
   public :: wannier_K2R
   public :: wannier_R2K
   !public :: add_crystalfields

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the Lattice vectors
   !---------------------------------------------------------------------------!
   subroutine read_lattice(pathINPUT)
      !
      use utils_misc
      use linalg
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=256)                    :: path
      integer                               :: unit
      logical                               :: filexists
      !
      !
      write(*,"(A)") "--- read_lattice ---"
      !
      !
      ! Look for LATTC
      path=reg(pathINPUT)//"LATTC.DAT"
      call inquireFile(reg(path),filexists)
      !
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="old",position="rewind",action="read")
      read(unit,*)
      read(unit,*) lat(1:3,1)
      read(unit,*) lat(1:3,2)
      read(unit,*) lat(1:3,3)
      close(unit)
      !
      rlat = lat
      call inv_sym(rlat)
      !
      rlat = 2*pi*transpose(rlat)
      !
      vol = det(lat)
      rvol = 8*pi**3 / vol
      write(*,"(A,1F10.6)")"Unit cell volume: ",vol
      !
      Lat_stored=.true.
      !
   end subroutine read_lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Read XEPS.DAT file
   !---------------------------------------------------------------------------!
   subroutine read_xeps(pathINPUT,kpt,Nkpt3,UseXepsKorder,kptPos,Nkpt_irred,UseDisentangledBS,spex_para)
      !
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      real(8),allocatable,intent(in)        :: kpt(:,:)
      integer,intent(in)                    :: Nkpt3(3)
      logical,intent(in)                    :: UseXepsKorder
      integer,allocatable,intent(inout)     :: kptPos(:)
      integer,intent(out)                   :: Nkpt_irred
      logical,intent(out)                   :: UseDisentangledBS
      logical,intent(in),optional           :: spex_para
      !
      character(len=256)                    :: path
      real(8),allocatable                   :: Ene_xeps(:,:,:)
      real(8),allocatable                   :: kpt_xeps(:,:)
      integer,allocatable                   :: kptPos_xeps(:)
      integer                               :: Nspin_xeps
      integer                               :: Nkpt3_xeps(3)
      integer                               :: Nkpt,Nkpt_xeps
      integer                               :: Nkpt_xeps_irred
      integer                               :: Nband_xeps
      real(8)                               :: Efermi_xeps
      logical                               :: spex_para_
      integer                               :: ik,unit
      logical                               :: dumlogical,filexists
      !
      !
      write(*,"(A)") "--- read_xeps ---"
      !
      !
      path = reg(pathINPUT)//"XEPS.DAT"
      call inquireFile(reg(path),filexists)
      !
      unit = free_unit()
      open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
      read(unit) Nspin_xeps,Nkpt3_xeps(:),Nkpt_xeps,Nkpt_xeps_irred, &
                 Nband_xeps,Efermi_xeps,dumlogical,UseDisentangledBS
      !
      write(*,"(A,1F10.5)") "Fermi energy in XEPS: ",Efermi_xeps
      !
      allocate(kpt_xeps(3,Nkpt_xeps));kpt_xeps=0d0
      allocate(kptPos_xeps(Nkpt_xeps));kptPos_xeps=0
      allocate(Ene_xeps(Nband_xeps,Nkpt_xeps_irred,Nspin_xeps));Ene_xeps=0d0
      !
      read(unit) kpt_xeps,kptPos_xeps
      read(unit) Ene_xeps
      !
      close(unit)
      !
      Nkpt = size(kpt,dim=2)
      Nkpt_irred = Nkpt
      if(UseXepsKorder) Nkpt_irred = Nkpt_xeps_irred
      !
      spex_para_=.true.
      if(present(spex_para))spex_para_=spex_para
      !
      ! Global checks
      if(spex_para_.and.(Nspin_xeps.ne.1)) stop "Nspin_xeps.ne.1 in XEPS.DAT"
      if(Nkpt_xeps.ne.Nkpt)         stop "Nkpt_xeps.ne.Nkpt in XEPS.DAT"
      if(Nkpt3_xeps(1).ne.Nkpt3(1)) stop "Nkpt(1)_xeps.ne.Nkpt(1) in XEPS.DAT"
      if(Nkpt3_xeps(2).ne.Nkpt3(2)) stop "Nkpt(2)_xeps.ne.Nkpt(2) in XEPS.DAT"
      if(Nkpt3_xeps(3).ne.Nkpt3(3)) stop "Nkpt(3)_xeps.ne.Nkpt(3) in XEPS.DAT"
      !
      ! Check of the K-point ordering
      do ik=1,Nkpt
         if (.not.keq(kpt_xeps(:,ik),kpt(:,ik))) then
            write(*,"(A)")"ik=",ik,"kpt(:,ik)=",kpt(:,ik),"kpt_loc(:,ik=)",kpt_xeps(:,ik)
            !write(*,"(A)") "kptp(ik)=",kptPos(ik),"kptp_loc(ik)=",kptPos_xeps(ik)
            stop "K-points grid does not match"
         endif
      enddo
      !
      call assert_shape(kptPos,[Nkpt],"read_xeps","kptPos")
      kptPos=0
      if(UseXepsKorder)then
         kptPos=kptPos_xeps
      else
         do ik=1,Nkpt
            kptPos(ik)=ik
         enddo
      endif
      !
   end subroutine read_xeps


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the Hamiltonian and kpoints providing Eigen-values/vectors
   !by now only for paramagnetic Hk
   !---------------------------------------------------------------------------!
   subroutine read_Hk(pathINPUT,Hk,kpt,Ek,Zk,Hloc)
      !
      use utils_misc
      use linalg, only :  eigh
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      complex(8),allocatable,intent(out)    :: Hk(:,:,:)
      real(8),allocatable,intent(out)       :: kpt(:,:)
      real(8),allocatable,intent(out)       :: Ek(:,:)
      complex(8),allocatable,intent(out)    :: Zk(:,:,:)
      complex(8),allocatable,intent(out)    :: Hloc(:,:)
      !
      character(len=256)                    :: path
      integer                               :: unit,Nkpt,Norb
      integer                               :: iwan1,iwan2,ik
      integer                               :: idum1,idum2
      real(8)                               :: ReHk,ImHk
      logical                               :: filexists
      !
      !
      write(*,"(A)") "--- read_Hk ---"
      !
      !
      ! Look for Hk.DAT
      path=reg(pathINPUT)//"Hk.DAT"
      call inquireFile(reg(path),filexists)
      !
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) idum1,Nkpt,Norb
      !
      allocate(Hk(Norb,Norb,Nkpt));Hk=czero
      allocate(kpt(3,Nkpt));kpt=0d0
      allocate(Ek(Norb,Nkpt));Ek=0d0
      allocate(Zk(Norb,Norb,Nkpt));Zk=czero
      allocate(Hloc(Norb,Norb));Hloc=czero
      !
      Hk=czero
      Zk=czero
      Hloc=czero
      Ek=0d0
      do ik=1,nkpt
         read(unit,*) idum1,idum2,kpt(:,ik)
         if (idum2.ne.ik) stop "ik"
         do iwan1=1,Norb
            do iwan2=1,Norb
               read(unit,*) idum1,idum2,ReHk,ImHk
               if (idum1.ne.iwan1) stop "iwan1"
               if (idum2.ne.iwan2) stop "iwan2"
               Hk(iwan1,iwan2,ik) = dcmplx(ReHk,ImHk)*H2eV
            enddo
         enddo
         Hloc = Hloc + Hk(:,:,ik)/nkpt
         !
         call check_Hermiticity(Hk(:,:,ik),eps)
         !
         Ek(:,ik) = 0d0
         Zk(:,:,ik) = Hk(:,:,ik)
         call eigh(Zk(:,:,ik),Ek(:,ik))
         !
      enddo
      !
      Hk_stored=.true.
      call read_lattice(reg(pathINPUT))
      !
   end subroutine read_Hk


   !---------------------------------------------------------------------------!
   !PURPOSE: Fill up the list likning indexes of the sum and diff of K-points
   !---------------------------------------------------------------------------!
   subroutine fill_ksumkdiff(kpt,kptsum,kptdif,nkpt3,pkpt)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: kpt(:,:)
      integer,intent(inout)                 :: kptsum(:,:)
      integer,intent(inout)                 :: kptdif(:,:)
      integer,intent(in),optional           :: nkpt3(:)
      integer,intent(inout),optional        :: pkpt(:,:,:)
      !
      integer                               :: Nkpt
      integer                               :: iq,ik1,ik2
      integer                               :: i1,i2,i3
      integer                               :: k(3)
      real(8)                               :: dk(3)
      integer,allocatable                   :: pkpt_(:,:,:)
      !
      !
      write(*,"(A)") "--- fill_ksumkdiff ---"
      !
      !
      Nkpt = size(kpt,dim=2)
      call assert_shape(kptsum,[Nkpt,Nkpt],"fill_ksumkdiff","kptsum")
      call assert_shape(kptdif,[Nkpt,Nkpt],"fill_ksumkdiff","kptdif")
      !
      ! k1+q=k2
      kptsum=0
      do iq=1,Nkpt
         do ik1=1,iq
            do ik2=1,Nkpt
               dk(:)=kpt(:,ik1)+kpt(:,iq)-kpt(:,ik2)
               dk(:)=dk(:)-nint(dk(:))
               if (all(abs(dk(:)).lt.eps)) then
                  kptsum(ik1,iq)=ik2
                  if (ik1.ne.iq) kptsum(iq,ik1)=kptsum(ik1,iq)
                  exit
               endif
            enddo ! ik2
            if (kptsum(ik1,iq).eq.0) stop "kptsum"
         enddo ! ik1
      enddo ! iq
      !
      ! k1-q=k2
      kptdif=0
      do iq=1,nkpt
         do ik1=1,nkpt
            do ik2=1,nkpt
               dk(:)=kpt(:,ik1)-kpt(:,iq)-kpt(:,ik2)
               dk(:)=dk(:)-nint(dk(:))
               if (all(abs(dk(:)).lt.eps)) then
                  kptdif(ik1,iq)=ik2
                  exit
               endif
            enddo ! ik2
            if (kptdif(ik1,iq).eq.0) stop "kptdif"
            if (kptdif(ik1,iq).gt.nkpt) stop "kptdif2"
         enddo ! ik1
      enddo ! iq
      !
      ! Not sure if the following is needed
      if(present(nkpt3))then
         call assert_shape(nkpt3,[3],"fill_ksumkdiff","nkpt3")
         allocate(pkpt_(nkpt3(1)+1,nkpt3(2)+1,nkpt3(3)+1));pkpt_=0
         do i1=1,nkpt3(1)+1
            do i2=1,nkpt3(2)+1
               do i3=1,nkpt3(3)+1
                  k(1)=mod(i1-1,nkpt3(1))
                  k(2)=mod(i2-1,nkpt3(2))
                  k(3)=mod(i3-1,nkpt3(3))
                  k(:)=k(:)/nkpt3(:)
                  do ik1=1,nkpt
                     dk(:)=kpt(:,ik1)-k(:)
                     dk(:)=dk(:)-nint(dk(:))
                     if (all(abs(dk(:)).lt.eps)) then
                        pkpt_(i1,i2,i3)=ik1
                        exit
                     endif
                  enddo
                  if (pkpt_(i1,i2,i3).eq.0) stop "pkpt"
               enddo
            enddo
         enddo
         if(present(pkpt))then
            call assert_shape(pkpt,[nkpt3(1)+1,nkpt3(2)+1,nkpt3(3)+1],"fill_ksumkdiff","pkpt")
            pkpt = pkpt_
         endif
      endif
      !
   end subroutine fill_ksumkdiff



   !---------------------------------------------------------------------------!
   !PURPOSE: find the smallest k vectors for interaction interpolation
   !---------------------------------------------------------------------------!
   subroutine fill_smallk(kpt,small_ik)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: kpt(:,:)
      integer,intent(inout)                 :: small_ik(:,:)
      !
      integer                               :: Nkpt,ik,i,j
      real(8)                               :: kreal1(3),kreal(3,12)
      !
      !
      write(*,"(A)") "--- fill_smallk ---"
      if(.not.Lat_stored)stop "Lattice positions not stored. Either call read_lattice(path) or read_Hk(path,Hk,kpt)"
      !
      !
      Nkpt = size(kpt,dim=2)
      call assert_shape(small_ik,[12,2],"fill_smallk","small_ik")
      !
      small_ik(:,1)=nkpt
      do ik=1,nkpt
         !
         !avoid the gamma point
         if(all(kpt(:,ik).eq.[0d0,0d0,0d0]))cycle
         !
         kreal1 = kpt(1,ik)*rlat(:,1) + kpt(2,ik)*rlat(:,2) + kpt(3,ik)*rlat(:,3)
         do i=1,12
            kreal(:,i) = kpt(1,small_ik(i,1))*rlat(:,1) + kpt(2,small_ik(i,1))*rlat(:,2) + kpt(3,small_ik(i,1))*rlat(:,3)
            if ((kreal1(1)**2+kreal1(2)**2+kreal1(3)**2).lt.(kreal(1,i)**2+kreal(2,i)**2+kreal(3,i)**2)) then
               do j=12,i+1,-1
                  small_ik(j,1)=small_ik(j-1,1)
               enddo
               small_ik(i,1)=ik
               exit
            endif
         enddo
         !
      enddo
      !
      !find kpoints that are equally far from the origin
      small_ik(1,2)=1
      do i=2,12
         if (abs((kreal(1,i-1)**2+kreal(2,i-1)**2+kreal(3,i-1)**2)-(kreal(1,i)**2+kreal(2,i)**2+kreal(3,i)**2)).lt.eps) then
            small_ik(i,2)=small_ik(i-1,2)
         else
            small_ik(i,2)=small_ik(i-1,2)+1
         endif
      enddo
      write(*,"(A)")"12 smallest k vectors are:"
      do i=1,12
         write(*,"(12(3F5.2,1X))")kpt(:,small_ik(i,1))
         write(*,"(3F5.2,1X)")kreal(:,i)
         write(*,"(1I5)")small_ik(i,2)
      enddo
      !
   end subroutine fill_smallk


   !---------------------------------------------------------------------------!
   !PURPOSE: calculate lattice points inside Wigner-Seitz cell of given supercell
   !NEEDED:
   ! lat(3,3) : premitive lattice vectors
   ! nkpt : number of k-points
   ! nkpt3(3) : number of k-points along each direction
   !OUTPUT:
   ! nwig : number of points inside wigner-seitz supercell
   ! rvec(3,2*nkpt) : (INTEGER) lattice vectors
   ! nrdeg(2*nkpt) : degeneracy
   !---------------------------------------------------------------------------!
   subroutine calc_wignerseiz(nkpt,nkpt3)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)           :: nkpt,nkpt3(3)
      integer                      :: ir1,ir2,ir3,irsc1,irsc2,irsc3
      integer,parameter            :: nshell=2
      double precision             :: rtmp(3),rtmpsc(3),dr(3)
      integer                      :: i,i0
      double precision,allocatable :: dist(:)
      double precision             :: distmin
      !
      !
      write(*,"(A)") "--- calc_wignerseiz ---"
      if(.not.Lat_stored)stop "Lattice positions not stored. Either call read_lattice(path) or read_Hk(path,Hk,kpt)."
      !
      !
      allocate(rvecwig(3,10*nkpt));rvecwig=0
      allocate(nrdegwig(10*nkpt));nrdegwig=0
      allocate(dist((2*nshell+1)**3))
      ! this i0 corresponds to irsc1=irsc2=irsc3=0
      i0=nshell*(1+(2*nshell+1)*(1+(2*nshell+1)))+1
      !
      nwig=0
      do ir1=-nkpt3(1),+nkpt3(1)
         do ir2=-nkpt3(2),+nkpt3(2)
            do ir3=-nkpt3(3),+nkpt3(3)
               rtmp(:)=matmul(lat,(/ir1,ir2,ir3/))
               i=0
               !
               do irsc1=-nshell,+nshell
                  do irsc2=-nshell,+nshell
                     do irsc3=-nshell,+nshell
                        i=i+1
                        rtmpsc(:)=matmul(lat,(/nkpt3(1)*irsc1,nkpt3(2)*irsc2,nkpt3(3)*irsc3/))
                        dr(:)=rtmp(:)-rtmpsc(:)
                        dist(i)=sum(dr(:)**2)
                     enddo ! irsc3
                  enddo ! irsc2
               enddo ! irsc1
               !
               distmin=minval(dist(:))
               if (abs(distmin-dist(i0)).le.epsWig) then
                  nwig=nwig+1
                  if (nwig.gt.10*nkpt) stop "nwig>10*nkpt"
                  rvecwig(:,nwig)=(/ir1,ir2,ir3/)
                  nrdegwig(nwig)=count(abs(distmin-dist(:)).le.epsWig)
                  !write(*,"(A)") nwig,rvecwig(:,nwig),nrdegwig(nwig)
               endif
               !
            enddo
         enddo
      enddo
      deallocate(dist)
      !
      if (abs(sum(1d0/nrdegwig(1:nwig))-nkpt).gt.epsWig) then
         write(*,"(A)") "Error: sum(1/nrdeg(:))=",sum(1d0/nrdegwig(1:nwig))
         stop "nrdeg"
      endif
      !
      Wig_stored=.true.
      !
   end subroutine calc_wignerseiz


   !---------------------------------------------------------------------------!
   !PURPOSE: Set the positions of the non-equivalent sites within the unit cell
   !---------------------------------------------------------------------------!
   subroutine set_siteposition(Rsites)
      use utils_misc
      implicit none
      real(8),intent(in)           :: Rsites(:,:)
      integer                      :: isite,Nsite
      write(*,"(A)") "--- calc_wignerseiz ---"
      Nsite = size(Rsites,dim=1)
      if(allocated(rsite))deallocate(rsite)
      allocate(rsite(Nsite,3));rsite=0d0
      do isite=1,Nsite
         rsite(isite,:) = Rsites(isite,:)
      enddo
      Ruc_stored=.true.
   end subroutine set_siteposition


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolates a K-dependent matrix between two different K meshes
   !---------------------------------------------------------------------------!
   subroutine wannierinterpolation_mat(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_orig(:,:,:,:)
      complex(8),intent(inout)              :: mat_intp(:,:,:,:)
      !
      integer                               :: Nkpt_orig,Nkpt_intp
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id,i1,i2
      real(8)                               :: kR
      complex(8)                            :: cfac
      complex(8),allocatable                :: mat_R(:,:,:,:)
      !
      !
      write(*,"(A)") "--- wannierinterpolation_mat ---"
      if(.not.Wig_stored)then!stop "Wigner Seiz cell not initialized. Call calc_wignerseiz."
         write(*,"(A)") "Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_mat","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_orig,dim=3)
      if(size(mat_orig,dim=1).ne.size(mat_orig,dim=2)) stop "mat_orig not square."
      Nsize = size(mat_orig,dim=1)
      call assert_shape(mat_intp,[Nsize,Nsize,Npoints,Nkpt_intp],"wannierinterpolation_mat","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize,Nsize,Npoints,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Npoints,Nsize,kpt_orig,rvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,id,i1,i2,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do ik=1,Nkpt_orig
            !
            kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),-sin(kR))
            !
            do id=1,Npoints
               do i1=1,Nsize
                  do i2=1,Nsize
                     !
                     mat_R(i1,i2,id,ir) = mat_R(i1,i2,id,ir) + mat_orig(i1,i2,id,ik)*cfac
                     !
                  enddo
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp(:,:,:,:)=0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Npoints,Nsize,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,id,i1,i2,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do ir=1,Nwig
            !
            kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
            !
            do id=1,Npoints
               do i1=1,Nsize
                  do i2=1,Nsize
                     !
                     mat_intp(i1,i2,id,ik) = mat_intp(i1,i2,id,ik) + mat_R(i1,i2,id,ir)*cfac
                     !
                  enddo
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_mat
   !
   subroutine wannierinterpolation_vec(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_orig(:,:,:)
      complex(8),intent(inout)              :: mat_intp(:,:,:)
      !
      integer                               :: Nkpt_orig,Nkpt_intp
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id,i1
      real(8)                               :: kR
      complex(8)                            :: cfac
      complex(8),allocatable                :: mat_R(:,:,:)
      !
      !
      write(*,"(A)") "--- wannierinterpolation_vec ---"
      if(.not.Wig_stored)then!stop "Wigner Seiz cell not initialized. Call calc_wignerseiz."
         write(*,"(A)") "Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_vec","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_orig,dim=2)
      Nsize = size(mat_orig,dim=1)
      call assert_shape(mat_intp,[Nsize,Npoints,Nkpt_intp],"wannierinterpolation_vec","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize,Npoints,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Npoints,Nsize,kpt_orig,rvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,id,i1,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do ik=1,Nkpt_orig
            !
            kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),-sin(kR))
            !
            do id=1,Npoints
               do i1=1,Nsize
                  !
                  mat_R(i1,id,ir) = mat_R(i1,id,ir) + mat_orig(i1,id,ik)*cfac
                  !
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Npoints,Nsize,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,id,i1,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do ir=1,Nwig
            !
            kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
            !
            do id=1,Npoints
               do i1=1,Nsize
                  !
                  mat_intp(i1,id,ik) = mat_intp(i1,id,ik) + mat_R(i1,id,ir)*cfac
                  !
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_vec


   !---------------------------------------------------------------------------!
   !PURPOSE: Extract from a K-dependent matrix the next-neighbor components
   !---------------------------------------------------------------------------!
   subroutine wannier_K2R_NN(nkpt3_orig,kpt_orig,mat_K,mat_R_nn)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      complex(8),intent(in)                 :: mat_K(:,:,:,:)
      complex(8),intent(inout)              :: mat_R_nn(:,:,:,:)
      !
      integer                               :: Nkpt_orig
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id,i1,i2,ir2
      real(8)                               :: kR
      logical                               :: Rx,Ry,Rz
      complex(8)                            :: cfac
      !
      !
      write(*,"(A)") "--- wannier_K2R_NN ---"
      if(.not.Wig_stored)then!stop "Wigner Seiz cell not initialized. Call calc_wignerseiz."
         write(*,"(A)") "Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannier_K2R_NN","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_K,dim=3)
      if(size(mat_K,dim=1).ne.size(mat_K,dim=2)) stop "mat_K not square."
      Nsize = size(mat_K,dim=1)
      call assert_shape(mat_R_nn,[Nsize,Nsize,Npoints,3],"wannierinterpolation","mat_R_nn")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Npoints,Nsize,kpt_orig,rvecwig,mat_K,mat_R_nn),&
      !$OMP PRIVATE(Rx,Ry,Rz,ir2,ir,ik,id,i1,i2,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         !
         Rx = all(rvecwig(:,ir).eq.[1,0,0])
         Ry = all(rvecwig(:,ir).eq.[0,1,0])
         Rz = all(rvecwig(:,ir).eq.[0,0,1])
         !
         if(Rx)then
            ir2 = 1
         elseif(Ry)then
            ir2 = 2
         elseif(Rz)then
            ir2 = 3
         else
            cycle
         endif
         !
         do ik=1,Nkpt_orig
            !
            kR=2*pi*dot_product(kpt_orig(:,ik),rvecwig(:,ir))
            cfac=dcmplx(cos(kR),-sin(kR))
            !
            do id=1,Npoints
               do i1=1,Nsize
                  do i2=1,Nsize
                     !
                     mat_R_nn(i1,i2,id,ir2) = mat_R_nn(i1,i2,id,ir2) + mat_K(i1,i2,id,ik)*cfac
                     !
                  enddo
               enddo
            enddo
            !
         enddo ! ik
         !
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_K2R_NN


   !---------------------------------------------------------------------------!
   !PURPOSE: Transforms a K-dependent matrix into Wannier basis
   !Nwig might not be already available in the main program
   !so I'm allocating here the mat_R which does not need to be allocated in
   !the calling routine.
   !---------------------------------------------------------------------------!
   subroutine wannier_K2R_mat(nkpt3_orig,kpt_orig,mat_K,mat_R)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      complex(8),intent(in)                 :: mat_K(:,:,:,:)
      complex(8),allocatable,intent(out)    :: mat_R(:,:,:,:)
      !
      integer                               :: Nkpt_orig
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id,i1,i2
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      write(*,"(A)") "--- wannier_K2R_mat ---"
      if(.not.Wig_stored)then!stop "Wigner Seiz cell not initialized. Call calc_wignerseiz."
         write(*,"(A)") "Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_K,dim=3)
      if(size(mat_K,dim=1).ne.size(mat_K,dim=2)) stop "mat_K not square."
      Nsize = size(mat_K,dim=1)
      if(.not.allocated(mat_R))then
         allocate(mat_R(Nsize,Nsize,Npoints,Nwig))
         mat_R=czero
      endif
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Npoints,Nsize,kpt_orig,rvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,id,i1,i2,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do ik=1,Nkpt_orig
            !
            kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),-sin(kR))
            !
            do id=1,Npoints
               do i1=1,Nsize
                  do i2=1,Nsize
                     !
                     mat_R(i1,i2,id,ir) = mat_R(i1,i2,id,ir) + mat_K(i1,i2,id,ik)*cfac
                     !
                  enddo
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
   end subroutine wannier_K2R_mat
   !
   subroutine wannier_K2R_vec(nkpt3_orig,kpt_orig,mat_K,mat_R)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      complex(8),intent(in)                 :: mat_K(:,:,:)
      complex(8),allocatable,intent(out)    :: mat_R(:,:,:)
      !
      integer                               :: Nkpt_orig
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id,i1
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      write(*,"(A)") "--- wannier_K2R_vec ---"
      if(.not.Wig_stored)then!stop "Wigner Seiz cell not initialized. Call calc_wignerseiz."
         write(*,"(A)") "Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_K,dim=2)
      Nsize = size(mat_K,dim=1)
      if(.not.allocated(mat_R))then
         allocate(mat_R(Nsize,Npoints,Nwig))
         mat_R=czero
      endif
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Npoints,Nsize,kpt_orig,rvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,id,i1,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do ik=1,Nkpt_orig
            !
            kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),-sin(kR))
            !
            do id=1,Npoints
               do i1=1,Nsize
                  !
                  mat_R(i1,id,ir) = mat_R(i1,id,ir) + mat_K(i1,id,ik)*cfac
                  !
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
   end subroutine wannier_K2R_vec


   !---------------------------------------------------------------------------!
   !PURPOSE: Transforms matrix in Wannier basis into K-space
   !---------------------------------------------------------------------------!
   subroutine wannier_R2K_mat(nkpt3_orig,kpt_intp,mat_R,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_R(:,:,:,:)
      complex(8),intent(inout)              :: mat_intp(:,:,:,:)
      !
      integer                               :: Nkpt_intp
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id,i1,i2
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      write(*,"(A)") "--- wannier_R2K_mat ---"
      if(.not.Wig_stored)then!stop "Wigner Seiz cell not initialized. Call calc_wignerseiz."
         write(*,"(A)") "Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_intp,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_R,dim=3)
      if(size(mat_R,dim=1).ne.size(mat_R,dim=2)) stop "mat_K not square."
      Nsize = size(mat_R,dim=1)
      call assert_shape(mat_intp,[Nsize,Nsize,Npoints,Nkpt_intp],"wannierinterpolation","mat_intp")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Npoints,Nsize,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,id,i1,i2,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do ir=1,Nwig
            !
            kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
            !
            do id=1,Npoints
               do i1=1,Nsize
                  do i2=1,Nsize
                     !
                     mat_intp(i1,i2,id,ik) = mat_intp(i1,i2,id,ik) + mat_R(i1,i2,id,ir)*cfac
                     !
                  enddo
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_R2K_mat
   !
   subroutine wannier_R2K_vec(nkpt3_orig,kpt_intp,mat_R,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_R(:,:,:)
      complex(8),intent(inout)              :: mat_intp(:,:,:)
      !
      integer                               :: Nkpt_intp
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id,i1
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      write(*,"(A)") "--- wannier_R2K_vec ---"
      if(.not.Wig_stored)then!stop "Wigner Seiz cell not initialized. Call calc_wignerseiz."
         write(*,"(A)") "Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_intp,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_R,dim=2)
      Nsize = size(mat_R,dim=1)
      call assert_shape(mat_intp,[Nsize,Npoints,Nkpt_intp],"wannierinterpolation","mat_intp")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Npoints,Nsize,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,id,i1,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do ir=1,Nwig
            !
            kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
            cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
            !
            do id=1,Npoints
               do i1=1,Nsize
                  !
                  mat_intp(i1,id,ik) = mat_intp(i1,id,ik) + mat_R(i1,id,ir)*cfac
                  !
               enddo
            enddo
            !
         enddo ! ik
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_R2K_vec


end module crystal
