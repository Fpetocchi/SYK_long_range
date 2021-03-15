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
      module procedure wannierinterpolation_d1                                  !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(d1,Nkpt_orig),mat_intp(d1,Nkpt_intp))
      module procedure wannierinterpolation_d2                                  !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(d1,d2,Nkpt_orig),mat_intp(d1,d2,Nkpt_intp))
      module procedure wannierinterpolation_d3                                  !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(d1,d2,d3,Nkpt_orig),mat_intp(d1,d2,d3,Nkpt_intp))
   end interface wannierinterpolation

   interface wannier_K2R
      module procedure wannier_K2R_d1                                           !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_K(d1,Nkpt_orig),mat_R(internally allocated))
      module procedure wannier_K2R_d2                                           !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_K(d1,d2,Nkpt_orig),mat_R(internally allocated))
      module procedure wannier_K2R_d3                                           !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_K(d1,d2,d3,Nkpt_orig),mat_R(internally allocated))
   end interface wannier_K2R

   interface wannier_R2K
      module procedure wannier_R2K_d1                                           !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_R(d1,Nwig),mat_K(d1,Nkpt_orig))
      module procedure wannier_R2K_d2                                           !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_R(d1,d2,Nwig),mat_K(d1,d2,Nkpt_orig))
      module procedure wannier_R2K_d3                                           !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),mat_R(d1,d2,d3,Nwig),mat_K(d1,d2,d3,Nkpt_orig))
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
   integer,public,protected                 :: Nwig=0
   integer,public,protected                 :: wig0=0
   integer,allocatable,public,protected     :: rvecwig(:,:)
   integer,allocatable,public,protected     :: nrdegwig(:)
   real(8),allocatable,public,protected     :: radiuswig(:)
   !
   logical,public,protected                 :: Hk_stored=.false.
   logical,public,protected                 :: Ruc_stored=.false.               !Global flag for routines that need positions within the u.c.
   logical,public,protected                 :: Lat_stored=.false.               !Internal flag for routines that need rlat
   logical,public,protected                 :: Wig_stored=.false.               !Internal flag for routines performing Wannier interpolation
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
   public :: read_lattice
   public :: read_xeps
   public :: read_Hk
   public :: build_Hk
   public :: fill_ksumkdiff
   public :: fill_smallk
   public :: set_siteposition
   public :: calc_wignerseiz
   public :: wannierinterpolation                                               ![Nkpt3_orig(3),Kpt_orig(3,Nkpt_orig),Kpt_intp(3,Nkpt_intp),Mat_orig(n,n,Npoins,Nkpt_orig),Mat_intp(n,n,Npoins,Nkpt_intp)]
   public :: wannier_K2R_NN                                                     ![Nkpt3_orig(3),Kpt_orig(3,Nkpt_orig),Mat_orig(n,n,Npoins,Nkpt_orig),mat_R_nn(n,n,Npoins,3)]
   public :: wannier_K2R
   public :: wannier_R2K
   public :: interpolateHk2Path
   public :: calc_path
   !public :: add_crystalfields

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the Lattice vectors
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine read_lattice(pathINPUT)
      !
      use utils_misc
      use linalg, only : det, inv_sym
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=256)                    :: path
      integer                               :: unit
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_lattice"
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
      if(verbose)write(*,"(A,F)")"     Unit cell volume: ",vol
      !
      Lat_stored=.true.
      !
   end subroutine read_lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Build the lattice vectors
   !---------------------------------------------------------------------------!
   subroutine set_lattice(Rinput)
      !
      use utils_misc
      use linalg, only : det, inv_sym
      implicit none
      !
      real(8),intent(in)                    :: Rinput(3,3)
      !
      !
      if(verbose)write(*,"(A)") "---- set_lattice"
      !
      !
      lat(:,1) = Rinput(:,1)
      lat(:,2) = Rinput(:,2)
      lat(:,3) = Rinput(:,3)
      !
      Lat_stored=.true.
      !
   end subroutine set_lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Build the k-point mesh
   !---------------------------------------------------------------------------!
   subroutine build_kpt(Nkpt3,kpt)
      !
      use utils_misc
      implicit none
      !
      integer,dimension(3),intent(in)       :: Nkpt3
      real(8),allocatable,intent(inout)     :: kpt(:,:)
      !
      integer                               :: ikx,iky,ikz,ik,Nkpt
      real(8)                               :: Dkx,Dky,Dkz
      real(8)                               :: kx,ky,kz
      !
      !
      if(verbose)write(*,"(A)") "---- build_kpt"
      !
      !
      Nkpt = Nkpt3(1)*Nkpt3(2)*Nkpt3(3)
      !
      if(allocated(kpt))deallocate(kpt)
      allocate(kpt(3,Nkpt));kpt=0d0
      !
      Dkx=0d0;Dky=0d0;Dkz=0d0
      if(Nkpt3(1).ne.1) Dkx = 1d0/dble(Nkpt3(1)-1)
      if(Nkpt3(2).ne.1) Dky = 1d0/dble(Nkpt3(2)-1)
      if(Nkpt3(3).ne.1) Dkz = 1d0/dble(Nkpt3(3)-1)
      !
      ik=0
      do ikx=1,Nkpt3(1)
         do iky=1,Nkpt3(2)
            do ikz=1,Nkpt3(3)
               !
               kx = (ikx-1)*Dkx
               ky = (iky-1)*Dky
               kz = (ikz-1)*Dkz
               !
               ik=ik+1
               kpt(:,ik) = [kx,ky,kz]
               !
            enddo
         enddo
      enddo
      !
   end subroutine build_kpt


   !---------------------------------------------------------------------------!
   !PURPOSE: Read XEPS.DAT file
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine read_xeps(pathINPUT,kpt,Nkpt3,UseXepsKorder,kptPos,Nkpt_irred,UseDisentangledBS,iq_gamma,spex_para)
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
      integer,intent(out)                   :: iq_gamma
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
      if(verbose)write(*,"(A)") "---- read_xeps"
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
      if(verbose)write(*,"(A,F)") "     Fermi energy in XEPS: ",Efermi_xeps
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
      iq_gamma = find_kpt([0d0,0d0,0d0],kpt_xeps,eps)
      !
      spex_para_=.true.
      if(present(spex_para))spex_para_=spex_para
      !
      ! Global checks
      if(spex_para_.and.(Nspin_xeps.ne.1)) stop "read_xeps: Nspin_xeps.ne.1 in XEPS.DAT"
      if(Nkpt_xeps.ne.Nkpt)         stop "read_xeps: Nkpt_xeps.ne.Nkpt in XEPS.DAT"
      if(Nkpt3_xeps(1).ne.Nkpt3(1)) stop "read_xeps: Nkpt(1)_xeps.ne.Nkpt(1) in XEPS.DAT"
      if(Nkpt3_xeps(2).ne.Nkpt3(2)) stop "read_xeps: Nkpt(2)_xeps.ne.Nkpt(2) in XEPS.DAT"
      if(Nkpt3_xeps(3).ne.Nkpt3(3)) stop "read_xeps: Nkpt(3)_xeps.ne.Nkpt(3) in XEPS.DAT"
      !
      ! Check of the K-point ordering
      do ik=1,Nkpt
         if (.not.keq(kpt_xeps(:,ik),kpt(:,ik))) then
            write(*,"(A)")"ik=",ik,"kpt(:,ik)=",kpt(:,ik),"kpt_loc(:,ik=)",kpt_xeps(:,ik)
            !write(*,"(A)") "kptp(ik)=",kptPos(ik),"kptp_loc(ik)=",kptPos_xeps(ik)
            stop "read_xeps: K-points grid does not match"
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
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine read_Hk(pathINPUT,alphaHk,Hk,kpt,Ek,Zk,Hloc,iq_gamma)
      !
      use utils_misc
      use linalg, only :  eigh
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      real(8),intent(in)                    :: alphaHk
      complex(8),allocatable,intent(out)    :: Hk(:,:,:)
      real(8),allocatable,intent(out)       :: kpt(:,:)
      real(8),allocatable,intent(out)       :: Ek(:,:)
      complex(8),allocatable,intent(out)    :: Zk(:,:,:)
      complex(8),allocatable,intent(out)    :: Hloc(:,:)
      integer,intent(out),optional          :: iq_gamma
      !
      character(len=256)                    :: path
      integer                               :: unit,Nkpt,Norb
      integer                               :: iwan1,iwan2,ik
      integer                               :: idum1,idum2
      real(8)                               :: ReHk,ImHk
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_Hk"
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
      if(allocated(Hk))deallocate(Hk)
      if(allocated(kpt))deallocate(kpt)
      if(allocated(Ek))deallocate(Ek)
      if(allocated(Zk))deallocate(Zk)
      if(allocated(Hloc))deallocate(Hloc)
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
               if (idum1.ne.iwan1) stop "read_Hk: wrong index iwan1."
               if (idum2.ne.iwan2) stop "read_Hk: wrong index iwan2."
               Hk(iwan1,iwan2,ik) = dcmplx(ReHk,ImHk)*H2eV*alphaHk
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
      if(present(iq_gamma))iq_gamma = find_kpt([0d0,0d0,0d0],kpt,eps)
      Hk_stored=.true.
      call read_lattice(reg(pathINPUT))
      !
   end subroutine read_Hk


   !---------------------------------------------------------------------------!
   !PURPOSE: Build the Hamiltonian and kpoints from user-given parameters
   !---------------------------------------------------------------------------!
   subroutine build_Hk(Rinput,Norb,hopping,Nkpt3,alphaHk,Hk,kpt,Ek,Zk,Hloc,iq_gamma,pathOUTPUT)
      !
      use utils_misc
      use linalg, only : zeye, diagonal
      implicit none
      !
      real(8),intent(in)                    :: Rinput(3,3)
      integer,intent(in)                    :: Norb
      real(8),intent(in)                    :: hopping(:)
      integer,intent(in)                    :: Nkpt3(3)
      real(8),intent(in)                    :: alphaHk
      complex(8),allocatable,intent(out)    :: Hk(:,:,:)
      real(8),allocatable,intent(out)       :: kpt(:,:)
      real(8),allocatable,intent(out)       :: Ek(:,:)
      complex(8),allocatable,intent(out)    :: Zk(:,:,:)
      complex(8),allocatable,intent(out)    :: Hloc(:,:)
      integer,intent(out),optional          :: iq_gamma
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      integer                               :: unit,Nkpt
      integer                               :: iwan1,iwan2,ik
      integer                               :: Trange,idist,iwig
      real(8),allocatable                   :: Rsorted(:)
      integer,allocatable                   :: Rorder(:)
      complex(8),allocatable                :: Hr(:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- build_Hk"
      !
      !
      Nkpt = Nkpt3(1)*Nkpt3(2)*Nkpt3(3)
      call assert_shape(hopping,[Norb],"build_Hk","hopping")
      !
      if(allocated(Hk))deallocate(Hk)
      if(allocated(kpt))deallocate(kpt)
      if(allocated(Ek))deallocate(Ek)
      if(allocated(Zk))deallocate(Zk)
      if(allocated(Hloc))deallocate(Hloc)
      !
      allocate(Hk(Norb,Norb,Nkpt));Hk=czero
      allocate(kpt(3,Nkpt));kpt=0d0
      allocate(Ek(Norb,Nkpt));Ek=0d0
      allocate(Zk(Norb,Norb,Nkpt));Zk=czero
      allocate(Hloc(Norb,Norb));Hloc=czero
      !
      call set_lattice(Rinput)
      call build_kpt(Nkpt3,kpt)
      !
      !recover the vectors in real space and allocate hopping in real space
      if(.not.Wig_stored)call calc_wignerseiz(Nkpt,Nkpt3)
      allocate(Rsorted(Nwig));Rsorted = radiuswig
      allocate(Rorder(Nwig))
      call sort_array(Rsorted,Rorder)
      allocate(Hr(Norb,Norb,Nwig));Hr=czero
      !
      !Hopping is only nearest neighbor by now
      Trange=1
      !
      !loop over the sorted Wigner Seiz positions
      idist=1
      loopwig:do iwig=1,Nwig
         !
         !setting the local energy
         if(Rsorted(Rorder(iwig)).eq.0d0)then
            if(Rorder(iwig).ne.wig0)stop "build_Hk: wrong index of R=0 vector."
            cycle
         endif
         !
         !increasing range
         if(iwig.gt.2)then
            if(Rsorted(Rorder(iwig)).gt.Rsorted(Rorder(iwig-1))) idist=idist+1
            if(idist.gt.Trange) exit loopwig
         endif
         !
         !setting matrix element
         do iwan1=1,Norb
            Hr(iwan1,iwan1,Rorder(iwig)) = -dcmplx(hopping(iwan1),0d0)
         enddo
         !
      enddo loopwig
      !
      if(verbose)then
         write(*,*)"     Real-space hopping elements:"
         do iwig=1,Nwig
            write(*,"(5X,2F12.6)")Rsorted(Rorder(iwig)),real(Hr(1,1,Rorder(iwig)))
         enddo
      endif
      !
      call wannier_R2K(Nkpt3,kpt,Hr,Hk)
      deallocate(Hr,Rorder,Rsorted)
      !
      Hk = Hk*alphaHk
      !
      do ik=1,Nkpt
         Ek(:,ik) = diagonal(Hk(:,:,ik))
         Zk(:,:,ik) = zeye(Norb)
      enddo
      Hloc = sum(Hk,dim=3)/Nkpt
      !
      if(present(iq_gamma))iq_gamma = 1
      Hk_stored=.true.
      !
      if(present(pathOUTPUT))then
         !
         unit = free_unit()
         open(unit,file=reg(pathOUTPUT)//"Hk.DAT",form="formatted",status="unknown",position="rewind",action="write")
         write(unit,("(3I10)")) 1,Nkpt,Norb
         do ik=1,Nkpt
            write(unit,("(3F14.8)")) kpt(:,ik)
            do iwan1=1,Norb
               do iwan2=1,Norb
                  write(unit,("(2I4,2E20.12)")) iwan1,iwan2,dreal(Hk(iwan1,iwan2,ik)),dimag(Hk(iwan1,iwan2,ik))
               enddo
            enddo
         enddo
         !
      endif
      !
   end subroutine build_Hk


   !---------------------------------------------------------------------------!
   !PURPOSE: Fill up the list likning indexes of the sum and diff of K-points
   !TEST ON: 14-10-2020
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
      if(verbose)write(*,"(A)") "---- fill_ksumkdiff"
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
            !
            !first attempt
            do ik2=1,Nkpt
               dk(:)=kpt(:,ik1)+kpt(:,iq)-kpt(:,ik2)
               dk(:)=dk(:)-nint(dk(:))
               if (all(abs(dk(:)).lt.eps)) then
                  kptsum(ik1,iq)=ik2
                  if (ik1.ne.iq) kptsum(iq,ik1)=kptsum(ik1,iq)
                  exit
               endif
            enddo ! ik2
            !
            !second attempt with larger tolerance
            if (kptsum(ik1,iq).eq.0)then
               do ik2=1,Nkpt
                  dk(:)=kpt(:,ik1)+kpt(:,iq)-kpt(:,ik2)
                  dk(:)=dk(:)-nint(dk(:))
                  if (all(abs(dk(:)).lt.1e-6)) then
                     kptsum(ik1,iq)=ik2
                     if (ik1.ne.iq) kptsum(iq,ik1)=kptsum(ik1,iq)
                     exit
                  endif
               enddo ! ik2
            endif
            !
            !missing sum
            if (kptsum(ik1,iq).eq.0) stop "fill_ksumkdiff: kptsum failed."
            !
         enddo ! ik1
      enddo ! iq
      !
      ! k1-q=k2
      kptdif=0
      do iq=1,nkpt
         do ik1=1,nkpt
            !
            !first attempt
            do ik2=1,nkpt
               dk(:)=kpt(:,ik1)-kpt(:,iq)-kpt(:,ik2)
               dk(:)=dk(:)-nint(dk(:))
               if (all(abs(dk(:)).lt.eps)) then
                  kptdif(ik1,iq)=ik2
                  exit
               endif
            enddo ! ik2
            !
            !second attempt with larger tolerance
            if (kptdif(ik1,iq).eq.0)then
               do ik2=1,nkpt
                  dk(:)=kpt(:,ik1)-kpt(:,iq)-kpt(:,ik2)
                  dk(:)=dk(:)-nint(dk(:))
                  if (all(abs(dk(:)).lt.1e-6)) then
                     kptdif(ik1,iq)=ik2
                     exit
                  endif
               enddo ! ik2
            endif
            !
            !missing difference
            if (kptdif(ik1,iq).eq.0) stop "fill_ksumkdiff: kptdif failed."
            if (kptdif(ik1,iq).gt.nkpt) stop "fill_ksumkdiff: kptdif2 failed."
            !
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
                  !
                  !first attempt
                  do ik1=1,nkpt
                     dk(:)=kpt(:,ik1)-k(:)
                     dk(:)=dk(:)-nint(dk(:))
                     if (all(abs(dk(:)).lt.eps)) then
                        pkpt_(i1,i2,i3)=ik1
                        exit
                     endif
                  enddo
                  !
                  !second attempt with larger tolerance
                  if (pkpt_(i1,i2,i3).eq.0)then
                     do ik1=1,nkpt
                        dk(:)=kpt(:,ik1)-k(:)
                        dk(:)=dk(:)-nint(dk(:))
                        if (all(abs(dk(:)).lt.1e-6)) then
                           pkpt_(i1,i2,i3)=ik1
                           exit
                        endif
                     enddo
                  endif
                  !
                  !missing positon
                  if (pkpt_(i1,i2,i3).eq.0) stop "fill_ksumkdiff: pkpt failed."
                  !
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
   !TEST ON: 14-10-2020
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
      if(verbose)write(*,"(A)") "---- fill_smallk"
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
      if(verbose)then
         write(*,"(A)")"     12 smallest k vectors are:"
         do i=1,12
            write(*,"(5X,12(3F6.3,1X))")kpt(:,small_ik(i,1))
            write(*,"(5X,3F5.2,1X)")kreal(:,i)
            write(*,"(5X,1I5)")small_ik(i,2)
         enddo
      endif
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
   !TEST ON: 16-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_wignerseiz(nkpt,nkpt3)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt
      integer,intent(in)                    :: nkpt3(3)
      !
      integer                               :: ir1,ir2,ir3,irsc1,irsc2,irsc3
      integer                               :: i,i0
      integer                               :: rscan(3),rshel(3),dims(3)
      integer,parameter                     :: nshell=2
      real(8)                               :: rtmp(3),rtmpsc(3),dr(3)
      real(8),allocatable                   :: dist(:),radiuswig_tmp(:)
      integer,allocatable                   :: rvecwig_tmp(:,:),nrdegwig_tmp(:)
      real(8)                               :: distmin
      !
      !
      if(verbose)write(*,"(A)") "---- calc_wignerseiz"
      if(.not.Lat_stored)stop "calc_wignerseiz: Lattice positions not stored. Either call read_lattice(path) or read_Hk(path,Hk,kpt)."
      !
      !public,protected
      allocate(rvecwig_tmp(3,10*nkpt));rvecwig_tmp=0
      allocate(nrdegwig_tmp(10*nkpt));nrdegwig_tmp=0
      allocate(radiuswig_tmp(10*nkpt));radiuswig_tmp=0d0
      !
      rscan=nkpt3
      if(nkpt3(1).eq.1)rscan(1)=0
      if(nkpt3(2).eq.1)rscan(2)=0
      if(nkpt3(3).eq.1)rscan(3)=0
      rshel=nshell
      if(nkpt3(1).eq.1)rshel(1)=0
      if(nkpt3(2).eq.1)rshel(2)=0
      if(nkpt3(3).eq.1)rshel(3)=0
      dims=(2*nshell+1)
      if(nkpt3(1).eq.1)dims(1)=1
      if(nkpt3(2).eq.1)dims(2)=1
      if(nkpt3(3).eq.1)dims(3)=1
      !
      !OLD
      !this i0 corresponds to irsc1=irsc2=irsc3=0
      !allocate(dist((2*nshell+1)**3))
      !i0=nshell*(1+(2*nshell+1)*(1+(2*nshell+1)))+1
      !
      !NEW
      allocate(dist(dims(1)*dims(2)*dims(3)))
      i0=int((dims(1)*dims(2)*dims(3))/2)+1
      !
      !
      nwig=0
      do ir1=-rscan(1),+rscan(1)
         do ir2=-rscan(2),+rscan(2)
            do ir3=-rscan(3),+rscan(3)
               rtmp(:)=matmul(lat,(/ir1,ir2,ir3/))
               i=0
               !
               do irsc1=-rshel(1),+rshel(1)
                  do irsc2=-rshel(2),+rshel(2)
                     do irsc3=-rshel(3),+rshel(3)
                        i=i+1
                        rtmpsc(:)=matmul(lat,(/rscan(1)*irsc1,rscan(2)*irsc2,rscan(3)*irsc3/))
                        dr(:)=rtmp(:)-rtmpsc(:)
                        dist(i)=sum(dr(:)**2)
                        if((i.eq.i0).and.(.not.all([irsc1,irsc2,irsc3].eq.[0,0,0])))stop "calc_wignerseiz: wrong index of R=0 vector."
                     enddo ! irsc3
                  enddo ! irsc2
               enddo ! irsc1
               !
               distmin=minval(dist(:))
               if (abs(distmin-dist(i0)).le.epsWig) then
                  nwig=nwig+1
                  if (nwig.gt.10*nkpt) stop "calc_wignerseiz: nwig>10*nkpt."
                  rvecwig_tmp(:,nwig)=(/ir1,ir2,ir3/)
                  nrdegwig_tmp(nwig)=count(abs(distmin-dist(:)).le.epsWig)
                  radiuswig_tmp(nwig)=sqrt(dble(dot_product([ir1,ir2,ir3],[ir1,ir2,ir3])))
                  if(all([ir1,ir2,ir3].eq.[0,0,0]))wig0=nwig
                  !if(verbose)write(*,*) nwig,rvecwig_tmp(:,nwig),nrdegwig(nwig)
               endif
               !
            enddo
         enddo
      enddo
      deallocate(dist)
      !
      if (abs(sum(1d0/nrdegwig_tmp(1:nwig))-nkpt).gt.epsWig) then
         write(*,"(A,F)") "Error: sum(1/nrdeg(:))=",sum(1d0/nrdegwig_tmp(1:nwig))
         stop "calc_wignerseiz: nrdeg failed."
      endif
      !
      !public,protected
      allocate(rvecwig(3,nwig)) ; rvecwig=rvecwig_tmp(:,1:nwig)
      allocate(nrdegwig(nwig))  ; nrdegwig=nrdegwig_tmp(1:nwig)
      allocate(radiuswig(nwig)) ; radiuswig=radiuswig_tmp(1:nwig)
      Wig_stored=.true.
      !
      deallocate(rvecwig_tmp,nrdegwig_tmp,radiuswig_tmp)
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
      if(verbose)write(*,"(A)") "---- set_siteposition"
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
   subroutine wannierinterpolation_d1(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_orig(:,:)
      complex(8),intent(inout)              :: mat_intp(:,:)
      !
      integer                               :: Nkpt_orig,Nkpt_intp
      integer                               :: Nsize1
      integer                               :: i1
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      complex(8),allocatable                :: mat_R(:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- wannierinterpolation_d1"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_d1","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_d1: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_d1: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      call assert_shape(mat_intp,[Nsize1,Nkpt_intp],"wannierinterpolation_d1","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,kpt_orig,rvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         !
         do ir=1,Nwig
            do ik=1,Nkpt_orig
               !
               kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
               cfac = dcmplx(cos(kR),-sin(kR))
               !
               mat_R(i1,ir) = mat_R(i1,ir) + mat_orig(i1,ik)*cfac
               !
            enddo ! ik
         enddo ! ir
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         !
         do ik=1,Nkpt_intp
            do ir=1,Nwig
               !
               kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
               cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
               !
               mat_intp(i1,ik) = mat_intp(i1,ik) + mat_R(i1,ir)*cfac
               !
            enddo ! ir
         enddo ! ik
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_d1
   !
   subroutine wannierinterpolation_d2(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
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
      integer                               :: Nsize1,Nsize2
      integer                               :: i1,i2
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      complex(8),allocatable                :: mat_R(:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- wannierinterpolation_d2"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_d2","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_d2: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_d2: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      Nsize2 = size(mat_orig,dim=2)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nkpt_intp],"wannierinterpolation_d2","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nsize2,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,kpt_orig,rvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            !
            do ir=1,Nwig
               do ik=1,Nkpt_orig
                  !
                  kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
                  cfac = dcmplx(cos(kR),-sin(kR))
                  !
                  mat_R(i1,i2,ir) = mat_R(i1,i2,ir) + mat_orig(i1,i2,ik)*cfac
                  !
               enddo ! ik
            enddo ! ir
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            !
            do ik=1,Nkpt_intp
               do ir=1,Nwig
                  !
                  kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
                  cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                  !
                  mat_intp(i1,i2,ik) = mat_intp(i1,i2,ik) + mat_R(i1,i2,ir)*cfac
                  !
               enddo ! ir
            enddo ! ik
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_d2
   !
   subroutine wannierinterpolation_d3(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
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
      integer                               :: Nsize1,Nsize2,Nsize3
      integer                               :: i1,i2,i3
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      complex(8),allocatable                :: mat_R(:,:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- wannierinterpolation_d3"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_d3","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_d3: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_d3: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      Nsize2 = size(mat_orig,dim=2)
      Nsize3 = size(mat_orig,dim=3)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nsize3,Nkpt_intp],"wannierinterpolation_d3","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nsize2,Nsize3,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,Nsize3,kpt_orig,rvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            do i3=1,Nsize3
               !
               do ir=1,Nwig
                  do ik=1,Nkpt_orig
                     !
                     kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
                     cfac = dcmplx(cos(kR),-sin(kR))
                     !
                     mat_R(i1,i2,i3,ir) = mat_R(i1,i2,i3,ir) + mat_orig(i1,i2,i3,ik)*cfac
                     !
                  enddo ! ik
               enddo ! ir
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp(:,:,:,:)=0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,Nsize3,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            do i3=1,Nsize3
               !
               do ik=1,Nkpt_intp
                  do ir=1,Nwig
                     !
                     kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
                     cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                     !
                     mat_intp(i1,i2,i3,ik) = mat_intp(i1,i2,i3,ik) + mat_R(i1,i2,i3,ir)*cfac
                     !
                  enddo ! ir
               enddo ! ik
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_d3


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
      if(verbose)write(*,"(A)") "---- wannier_K2R_NN"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannier_K2R_NN","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_NN: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Npoints = size(mat_K,dim=3)
      if(size(mat_K,dim=1).ne.size(mat_K,dim=2)) stop "wannier_K2R_NN: mat_K not square."
      Nsize = size(mat_K,dim=1)
      call assert_shape(mat_R_nn,[Nsize,Nsize,Npoints,3],"wannierinterpolation","mat_R_nn")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Npoints,Nsize,kpt_orig,rvecwig,mat_K,mat_R_nn),&
      !$OMP PRIVATE(Rx,Ry,Rz,ir2,ir,ik,id,i1,i2,kR,cfac)
      !$OMP DO
      do i1=1,Nsize
         do i2=1,Nsize
            do id=1,Npoints
               !
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
                     mat_R_nn(i1,i2,id,ir2) = mat_R_nn(i1,i2,id,ir2) + mat_K(i1,i2,id,ik)*cfac
                     !
                  enddo ! ik
               enddo ! ir
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_K2R_NN


   !---------------------------------------------------------------------------!
   !PURPOSE: Transforms a K-dependent matrix into Wannier basis
   !Nwig might not be already available in the main program
   !so I'm allocating here the mat_R which does not need to be allocated in
   !the calling routine.
   !TEST ON: 16-10-2020(mat)
   !---------------------------------------------------------------------------!
   subroutine wannier_K2R_d1(nkpt3_orig,kpt_orig,mat_K,mat_R)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      complex(8),intent(in)                 :: mat_K(:,:)
      complex(8),allocatable,intent(out)    :: mat_R(:,:)
      !
      integer                               :: Nkpt_orig
      integer                               :: Nsize1
      integer                               :: i1
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      if(verbose)write(*,"(A)") "---- wannier_K2R_d1"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_d1: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_K,dim=1)
      !
      if(allocated(mat_R))deallocate(mat_R)
      allocate(mat_R(Nsize1,Nwig));mat_R=czero
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,kpt_orig,rvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         !
         do ir=1,Nwig
            do ik=1,Nkpt_orig
               !
               kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
               cfac = dcmplx(cos(kR),-sin(kR))
               !
               mat_R(i1,ir) = mat_R(i1,ir) + mat_K(i1,ik)*cfac
               !
            enddo ! ik
         enddo ! ir
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
   end subroutine wannier_K2R_d1
   !
   subroutine wannier_K2R_d2(nkpt3_orig,kpt_orig,mat_K,mat_R)
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
      integer                               :: Nsize1,Nsize2
      integer                               :: i1,i2
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      if(verbose)write(*,"(A)") "---- wannier_K2R_d2"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_d2: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_K,dim=1)
      Nsize2 = size(mat_K,dim=2)
      !
      if(allocated(mat_R))deallocate(mat_R)
      allocate(mat_R(Nsize1,Nsize2,Nwig));mat_R=czero
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,kpt_orig,rvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            !
            do ir=1,Nwig
               do ik=1,Nkpt_orig
                  !
                  kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
                  cfac = dcmplx(cos(kR),-sin(kR))
                  !
                  mat_R(i1,i2,ir) = mat_R(i1,i2,ir) + mat_K(i1,i2,ik)*cfac
                  !
               enddo ! ik
            enddo ! ir
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
   end subroutine wannier_K2R_d2
   !
   subroutine wannier_K2R_d3(nkpt3_orig,kpt_orig,mat_K,mat_R)
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
      integer                               :: Nsize1,Nsize2,Nsize3
      integer                               :: i1,i2,i3
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      if(verbose)write(*,"(A)") "---- wannier_K2R_d3"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_orig,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_d3: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_K,dim=1)
      Nsize2 = size(mat_K,dim=2)
      Nsize3 = size(mat_K,dim=3)
      !
      if(allocated(mat_R))deallocate(mat_R)
      allocate(mat_R(Nsize1,Nsize2,Nsize3,Nwig));mat_R=czero
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,Nsize3,kpt_orig,rvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            do i3=1,Nsize3
               !
               do ir=1,Nwig
                  do ik=1,Nkpt_orig
                     !
                     kR = 2*pi * dot_product(kpt_orig(:,ik),rvecwig(:,ir))
                     cfac = dcmplx(cos(kR),-sin(kR))
                     !
                     mat_R(i1,i2,i3,ir) = mat_R(i1,i2,i3,ir) + mat_K(i1,i2,i3,ik)*cfac
                     !
                  enddo ! ik
               enddo ! ir
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
   end subroutine wannier_K2R_d3


   !---------------------------------------------------------------------------!
   !PURPOSE: Transforms matrix in Wannier basis into K-space
   !TEST ON: 16-10-2020(mat)
   !---------------------------------------------------------------------------!
   subroutine wannier_R2K_d1(nkpt3_orig,kpt_intp,mat_R,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_R(:,:)
      complex(8),intent(inout)              :: mat_intp(:,:)
      !
      integer                               :: Nkpt_intp
      integer                               :: Nsize1
      integer                               :: i1
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      if(verbose)write(*,"(A)") "---- wannier_R2K_d1"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_intp,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "wannier_R2K_d1: size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_R,dim=1)
      call assert_shape(mat_intp,[Nsize1,Nkpt_intp],"wannierinterpolation","mat_intp")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         !
         do ik=1,Nkpt_intp
            do ir=1,Nwig
               !
               kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
               cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
               !
               mat_intp(i1,ik) = mat_intp(i1,ik) + mat_R(i1,ir)*cfac
               !
            enddo ! ir
         enddo ! ik
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_R2K_d1
   !
   subroutine wannier_R2K_d2(nkpt3_orig,kpt_intp,mat_R,mat_intp)
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
      integer                               :: Nsize1,Nsize2
      integer                               :: i1,i2
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      if(verbose)write(*,"(A)") "---- wannier_R2K_d2"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_intp,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "wannier_R2K_d2: size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_R,dim=1)
      Nsize2 = size(mat_R,dim=2)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nkpt_intp],"wannierinterpolation","mat_intp")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            !
            do ik=1,Nkpt_intp
               do ir=1,Nwig
                  !
                  kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
                  cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                  !
                  mat_intp(i1,i2,ik) = mat_intp(i1,i2,ik) + mat_R(i1,i2,ir)*cfac
                  !
               enddo ! ir
            enddo ! ik
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_R2K_d2
   !
   subroutine wannier_R2K_d3(nkpt3_orig,kpt_intp,mat_R,mat_intp)
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
      integer                               :: Nsize1,Nsize2,Nsize3
      integer                               :: i1,i2,i3
      integer                               :: ik,ir
      real(8)                               :: kR
      complex(8)                            :: cfac
      !
      !
      if(verbose)write(*,"(A)") "---- wannier_R2K_d3"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation","nkpt3_orig")
         call calc_wignerseiz(size(kpt_intp,dim=2),nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "wannier_R2K_d3: size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop 'nkpt'
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_R,dim=1)
      Nsize2 = size(mat_R,dim=2)
      Nsize3 = size(mat_R,dim=3)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nsize3,Nkpt_intp],"wannierinterpolation","mat_intp")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,Nsize3,kpt_intp,rvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do i1=1,Nsize1
         do i2=1,Nsize2
            do i3=1,Nsize3
               !
               do ik=1,Nkpt_intp
                  do ir=1,Nwig
                     !
                     kR = 2*pi * dot_product(kpt_intp(:,ik),rvecwig(:,ir))
                     cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                     !
                     mat_intp(i1,i2,i3,ik) = mat_intp(i1,i2,i3,ik) + mat_R(i1,i2,i3,ir)*cfac
                     !
                  enddo ! ir
               enddo ! ik
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_R2K_d3


   !---------------------------------------------------------------------------!
   !PURPOSE: generates thew K-points along some pre-stored high-symmetry points.
   !         taken from https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#high-symmetry-pathsv
   !---------------------------------------------------------------------------!
   subroutine calc_path(kpt_path,structure,Nkpt_path,Kaxis,KaxisPoints)
      !
      use utils_misc
      implicit none
      !
      real(8),allocatable,intent(out)       :: kpt_path(:,:)
      character(len=*),intent(in)           :: structure
      integer,intent(in)                    :: Nkpt_path
      real(8),allocatable,intent(out),optional :: Kaxis(:)
      real(8),allocatable,intent(out),optional :: KaxisPoints(:)
      !
      real(8),dimension(3)                  :: Gamma,M,R,X,K,L,U,W,H,N,P,A,Z,S,T,Y
      real(8),dimension(3)                  :: Kdiff
      real(8),allocatable                   :: Kpoints(:,:),Kdist(:),Kturn(:)
      integer                               :: idir,Ndir,idk,ik,lastK
      real(8)                               :: dKtot,theta,phi,dk,kx,ky,kz
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_path"
      !
      !
      !
      if(allocated(kpt_path))deallocate(kpt_path)
      select case(reg(structure))
         case default
            !
            stop "Available structures: cubic, fcc, bcc, hex, tetragonal, orthorhombic."
            !
         case("cubic")
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            M     = [ 1d0/2d0, 1d0/2d0,     0d0 ]
            R     = [ 1d0/2d0, 1d0/2d0, 1d0/2d0 ]
            X     = [ 1d0/2d0,     0d0,     0d0 ]
            !
            allocate(Kpoints(3,8));Kpoints=0d0
            Kpoints(:,1) = Gamma
            Kpoints(:,2) = X
            Kpoints(:,3) = M
            Kpoints(:,4) = Gamma
            Kpoints(:,5) = R
            Kpoints(:,6) = X
            Kpoints(:,7) = M
            Kpoints(:,8) = R
            write(*,"(A)") "     structure: simple-cubic."
            write(*,"(A)") "     path: GXMGRX,MR"
            !
         case("fcc")
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            K     = [ 3d0/8d0, 3d0/8d0, 3d0/4d0 ]
            L     = [ 1d0/2d0, 1d0/2d0, 1d0/2d0 ]
            U     = [ 1d0/4d0, 5d0/8d0, 5d0/8d0 ]
            W     = [ 1d0/4d0, 1d0/2d0, 3d0/4d0 ]
            X     = [     0d0, 1d0/2d0, 1d0/2d0 ]
            !
            allocate(Kpoints(3,12));Kpoints=0d0
            Kpoints(:,1)  = Gamma
            Kpoints(:,2)  = X
            Kpoints(:,3)  = W
            Kpoints(:,4)  = K
            Kpoints(:,5)  = Gamma
            Kpoints(:,6)  = L
            Kpoints(:,7)  = U
            Kpoints(:,8)  = W
            Kpoints(:,9)  = L
            Kpoints(:,10) = K
            Kpoints(:,11) = U
            Kpoints(:,12) = X
            write(*,"(A)") "     structure: face-centred cubic."
            write(*,"(A)") "     path: GXWKGLUWLK,UX"
            !
         case("bcc")
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            H     = [-1d0/2d0, 1d0/2d0, 1d0/2d0 ]
            N     = [     0d0,     0d0, 1d0/2d0 ]
            P     = [ 1d0/4d0, 1d0/4d0, 1d0/4d0 ]
            !
            allocate(Kpoints(3,8));Kpoints=0d0
            Kpoints(:,1) = Gamma
            Kpoints(:,2) = H
            Kpoints(:,3) = N
            Kpoints(:,4) = Gamma
            Kpoints(:,5) = P
            Kpoints(:,6) = H
            Kpoints(:,7) = P
            Kpoints(:,8) = N
            write(*,"(A)") "     structure: body-centred cubic."
            write(*,"(A)") "     path: GHNGPH,PN"
            !
         case("hex")
            !
            A     = [     0d0,     0d0, 1d0/2d0 ]
            Gamma = [     0d0,     0d0,     0d0 ]
            H     = [ 1d0/3d0, 1d0/3d0, 1d0/2d0 ]
            K     = [ 1d0/3d0, 1d0/3d0,     0d0 ]
            L     = [     0d0, 1d0/2d0, 1d0/2d0 ]
            M     = [     0d0, 1d0/2d0,     0d0 ]
            !
            allocate(Kpoints(3,12));Kpoints=0d0
            Kpoints(:,1)  = Gamma
            Kpoints(:,2)  = M
            Kpoints(:,3)  = K
            Kpoints(:,4)  = Gamma
            Kpoints(:,5)  = A
            Kpoints(:,6)  = L
            Kpoints(:,7)  = H
            Kpoints(:,8)  = A
            Kpoints(:,9)  = L
            Kpoints(:,10) = M
            Kpoints(:,11) = K
            Kpoints(:,12) = H
            write(*,"(A)") "     structure: hexagonal."
            write(*,"(A)") "     path: GMKGALHA,LM,KH"
            !
         case("tetragonal")
            !
            A     = [ 1d0/2d0, 1d0/2d0, 1d0/2d0 ]
            Gamma = [     0d0,     0d0,     0d0 ]
            M     = [ 1d0/2d0, 1d0/2d0,     0d0 ]
            R     = [ 1d0/2d0,     0d0, 1d0/2d0 ]
            X     = [ 1d0/2d0,     0d0,     0d0 ]
            Z     = [     0d0,     0d0, 1d0/2d0 ]
            !
            allocate(Kpoints(3,12));Kpoints=0d0
            Kpoints(:,1)  = Gamma
            Kpoints(:,2)  = X
            Kpoints(:,3)  = M
            Kpoints(:,4)  = Gamma
            Kpoints(:,5)  = Z
            Kpoints(:,6)  = R
            Kpoints(:,7)  = A
            Kpoints(:,8)  = Z
            Kpoints(:,9)  = X
            Kpoints(:,10) = R
            Kpoints(:,11) = M
            Kpoints(:,12) = A
            write(*,"(A)") "     structure: tetragonal."
            write(*,"(A)") "     path: GXMGZRAZ,XR,MA"
            !
         case("orthorhombic_1")
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            R     = [ 1d0/2d0, 1d0/2d0, 1d0/2d0 ]
            S     = [ 1d0/2d0, 1d0/2d0,     0d0 ]
            T     = [     0d0, 1d0/2d0, 1d0/2d0 ]
            U     = [ 1d0/2d0,     0d0, 1d0/2d0 ]
            X     = [ 1d0/2d0,     0d0,     0d0 ]
            Y     = [     0d0, 1d0/2d0,     0d0 ]
            Z     = [     0d0,     0d0, 1d0/2d0 ]
            !
            allocate(Kpoints(3,16));Kpoints=0d0
            Kpoints(:,1)  = Gamma
            Kpoints(:,2)  = X
            Kpoints(:,3)  = S
            Kpoints(:,4)  = Y
            Kpoints(:,5)  = Gamma
            Kpoints(:,6)  = Z
            Kpoints(:,7)  = U
            Kpoints(:,8)  = R
            Kpoints(:,9)  = T
            Kpoints(:,10) = Z
            Kpoints(:,11) = Y
            Kpoints(:,12) = T
            Kpoints(:,13) = U
            Kpoints(:,14) = X
            Kpoints(:,15) = S
            Kpoints(:,16) = R
            write(*,"(A)") "     structure: orthorhombic - version 1."
            write(*,"(A)") "     path: GXSYGZURTZ,YT,UX,SR"
            !
         case("orthorhombic_2")
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            R     = [ 1d0/2d0, 1d0/2d0, 1d0/2d0 ]
            S     = [ 1d0/2d0, 1d0/2d0,     0d0 ]
            T     = [     0d0, 1d0/2d0, 1d0/2d0 ]
            U     = [ 1d0/2d0,     0d0, 1d0/2d0 ]
            X     = [ 1d0/2d0,     0d0,     0d0 ]
            Y     = [     0d0, 1d0/2d0,     0d0 ]
            Z     = [     0d0,     0d0, 1d0/2d0 ]
            !
            allocate(Kpoints(3,8));Kpoints=0d0
            Kpoints(:,1) = Gamma
            Kpoints(:,2) = X
            Kpoints(:,3) = S
            Kpoints(:,4) = Gamma
            Kpoints(:,5) = Y
            Kpoints(:,6) = T
            Kpoints(:,7) = Gamma
            Kpoints(:,8) = Z
            write(*,"(A)") "     structure: orthorhombic - version 2."
            write(*,"(A)") "     path: GXSGYTGZ"
            !
      end select
      !
      !
      Ndir = size(Kpoints,dim=2)
      allocate(Kdist((Ndir-1)*Nkpt_path+1));Kdist=0d0
      allocate(kpt_path(3,(Ndir-1)*Nkpt_path+1));kpt_path=0d0
      allocate(Kturn(Ndir));Kturn=0d0
      !
      ik=0
      do idir=2,Ndir
         !
         Kdiff = Kpoints(:,idir) - Kpoints(:,idir-1)
         !
         dKtot = sqrt(dot_product(Kdiff,Kdiff))
         theta = acos(Kdiff(3)/dKtot)
         phi = atan2(Kdiff(2),Kdiff(1))
         !
         dk = dKtot/Nkpt_path
         !
         lastK=0
         if(idir.eq.Ndir)lastK=1
         do idk=1,Nkpt_path+lastK
            !
            kx = Kpoints(1,idir-1) + (idk-1)*dk*sin(theta)*cos(phi)
            ky = Kpoints(2,idir-1) + (idk-1)*dk*sin(theta)*sin(phi)
            kz = Kpoints(3,idir-1) + (idk-1)*dk*cos(theta)
            !
            ik=ik+1
            !
            kpt_path(:,ik) = [kx,ky,kz]
            if(ik.gt.1)Kdist(ik) = Kdist(ik-1) + dk
            if(idk.eq.(Nkpt_path+lastK)) Kturn(idir) = Kdist(ik)
            !
         enddo
         !
      enddo
      Kdist = Kdist/Kdist((Ndir-1)*Nkpt_path+1)
      Kturn = Kturn/Kturn(Ndir)
      if(present(Kaxis))Kaxis=Kdist
      if(present(KaxisPoints))KaxisPoints=Kturn
      !
   end subroutine calc_path


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate to a user provided K-point path the Hamiltonian
   !---------------------------------------------------------------------------!
   subroutine interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT)
      !
      use parameters
      use utils_misc
      use linalg, only : eigh, inv
      implicit none
      !
      type(Lattice),intent(inout)           :: Lttc
      character(len=*),intent(in)           :: structure
      integer,intent(in)                    :: Nkpt_path
      character(len=*),intent(in)           :: pathOUTPUT
      !
      character(len=256)                    :: path
      integer                               :: ik,iorb,unit
      integer                               :: Norb
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- interpolateHk2Path"
      !
      !
      ! Check on the input Fields
      if(.not.Lttc%status) stop "interpolateHk2Path: Lttc not properly initialized."
      Norb = Lttc%Norb
      !
      !
      !Create K-points along high-symmetry points
      if(allocated(Lttc%kptpath))deallocate(Lttc%kptpath)
      if(allocated(Lttc%Kpathaxis))deallocate(Lttc%Kpathaxis)
      call calc_path(Lttc%kptpath,reg(structure),Nkpt_path,Kaxis=Lttc%Kpathaxis,KaxisPoints=Lttc%KpathaxisPoints)
      Lttc%Nkpt_path = size(Lttc%kptpath,dim=2)
      !
      !Fill in Hk along points
      if(allocated(Lttc%Hk_path))deallocate(Lttc%Hk_path)
      allocate(Lttc%Hk_path(Norb,Norb,Lttc%Nkpt_path));Lttc%Hk_path=czero
      call cpu_time(start)
      call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath,Lttc%Hk,Lttc%Hk_path)
      call cpu_time(finish)
      write(*,"(A,F)") "     H(fullBZ) --> H(Kpath) cpu timing:", finish-start
      !
      !Fill in Ek along points
      if(allocated(Lttc%Ek_path))deallocate(Lttc%Ek_path)
      allocate(Lttc%Ek_path(Norb,Lttc%Nkpt_path));Lttc%Ek_path=0d0
      if(allocated(Lttc%Zk_path))deallocate(Lttc%Zk_path)
      allocate(Lttc%Zk_path(Norb,Norb,Lttc%Nkpt_path));Lttc%Zk_path=czero
      do ik=1,Lttc%Nkpt_path
         Lttc%Zk_path(:,:,ik) = Lttc%Hk_path(:,:,ik)
         call eigh(Lttc%Zk_path(:,:,ik),Lttc%Ek_path(:,ik))
      enddo
      !
      !Print bands
      path = reg(pathOUTPUT)//"Bands.DAT"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      do ik=1,Lttc%Nkpt_path
         write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),(Lttc%Ek_path(:,ik),iorb=1,Norb)
      enddo
      close(unit)
      write(*,"(A,I)") "     Total number of K-points along path:",Lttc%Nkpt_path
      !
      !Print position of High-symmetry points in the same folder where the function is
      path = reg(pathOUTPUT)//"Kpoints.DAT"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      do ik=1,size(Lttc%KpathaxisPoints,dim=1)
         write(unit,"(1I5,200E20.12)") ik,Lttc%KpathaxisPoints(ik)
      enddo
      close(unit)
      write(*,"(A,I)") "     Total number of High symmetry points:",size(Lttc%KpathaxisPoints,dim=1)
      !
      Lttc%pathStored=.true.
      !
   end subroutine interpolateHk2Path


end module crystal
