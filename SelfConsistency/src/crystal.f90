module crystal

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !metto i vettori in spazio reale qui perche sevono anche a routine interne e
   !in genere mai al codice che lavora in spazio K.
   !cmq possono essere accessibili e sbattuti nell"attributo lattice con suroutine dedicate
   !il calc wigsiez e smallk hanno bisogno che LATTC sia letto
   ! qui mantengo il path allinputfile dato dall"utente perche cmq se leggo la Hk la devo leggere da qulahce parte
   ! per funzionare wannier intp devo leggere Hk perche" e" l"unico modo per dare il path

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface wannierinterpolation
      module procedure wannierinterpolation_D1_d
      module procedure wannierinterpolation_D2_d
      module procedure wannierinterpolation_D3_d
      module procedure wannierinterpolation_D1_z                                !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(d1,Nkpt_orig),mat_intp(d1,Nkpt_intp))
      module procedure wannierinterpolation_D2_z                                !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(d1,d2,Nkpt_orig),mat_intp(d1,d2,Nkpt_intp))
      module procedure wannierinterpolation_D3_z                                !(nkpt3_orig(3),kpt_orig(3,Nkpt_orig),kpt_intp(3,Nkpt_intp),mat_K(d1,d2,d3,Nkpt_orig),mat_intp(d1,d2,d3,Nkpt_intp))
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
   !PURPOSE: Module custom types
   !---------------------------------------------------------------------------!
   type symtype
      integer                               :: rot(3,3)                         !rotation in real space(lattice coords.)
      integer                               :: rrot(3,3)                        !rotation in recip. space(inverse of transpose(rot))
      real(8)                               :: transl(3)
      integer                               :: inv
      logical                               :: symmor
   end type symtype

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
   integer,allocatable                      :: Kprint(:)
   !
   real(8),private                          :: Rlat(3,3)
   real(8),private                          :: Blat(3,3)
   real(8),private                          :: vol
   !
   integer,public,protected                 :: Nwig=0
   integer,public,protected                 :: wig0=0
   real(8),allocatable,public,protected     :: radiuswig(:)
   integer,allocatable,public,protected     :: Nvecwig(:,:)
   real(8),allocatable,public,protected     :: Rvecwig(:,:)
   integer,allocatable,public,protected     :: nrdegwig(:)
   !
   real(8),allocatable,public,protected     :: UserPath(:,:)
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
#ifdef _Klib
   logical,private                          :: Klib=.true.
#else
   logical,private                          :: Klib=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: read_lattice
   public :: build_kpt
   public :: read_xeps
   public :: read_Hk
   public :: build_Hk
   public :: calc_irredBZ
   public :: fill_ksumkdiff
   public :: fill_smallk
   public :: set_siteposition
   public :: calc_wignerseiz
   public :: clear_wignerseiz
   public :: wannierinterpolation                                               ![Nkpt3_orig(3),Kpt_orig(3,Nkpt_orig),Kpt_intp(3,Nkpt_intp),Mat_orig(n,n,Npoins,Nkpt_orig),Mat_intp(n,n,Npoins,Nkpt_intp)]
   public :: wannier_K2R_NN                                                     ![Nkpt3_orig(3),Kpt_orig(3,Nkpt_orig),Mat_orig(n,n,Npoins,Nkpt_orig),mat_R_nn(n,n,Npoins,3)]
   public :: wannier_K2R
   public :: wannier_R2K
   public :: set_UserPath
   public :: get_Rlat,get_Blat
   public :: interpolateHk2Path
   public :: calc_Kpath
   public :: calc_Kplane
   public :: calc_Ewald
   public :: tetrahedron_integration

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the Lattice vectors
   !---------------------------------------------------------------------------!
   subroutine read_lattice(pathINPUT)
      !
      use utils_misc
      use linalg, only : det, inv_sym, cross_product
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=256)                    :: path
      integer                               :: unit,ir
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
      read(unit,*) Rlat(1:3,1)
      read(unit,*) Rlat(1:3,2)
      read(unit,*) Rlat(1:3,3)
      close(unit)
      vol = dot_product(cross_product(Rlat(:,1),Rlat(:,2)),Rlat(:,3))
      if(verbose)write(*,"(A,F)")"     Unit cell volume: ",vol
      Blat(:,1) = cross_product(Rlat(:,2),Rlat(:,3))/vol
      Blat(:,2) = cross_product(Rlat(:,3),Rlat(:,1))/vol
      Blat(:,3) = cross_product(Rlat(:,1),Rlat(:,2))/vol
      !
      write(*,"(A)")new_line("A")//"     Unit cell vectors: "
      do ir=1,3
         write(*,"(A)")"     R_"//str(ir)//": [ "//str(Rlat(1,ir),3)//" , "//str(Rlat(2,ir),3)//" , "//str(Rlat(3,ir),3)//" ]"
      enddo
      write(*,"(A)")new_line("A")//"     Reciprocal lattice vectors: "
      do ir=1,3
         write(*,"(A)")"     B_"//str(ir)//": [ "//str(Blat(1,ir),3)//" , "//str(Blat(2,ir),3)//" , "//str(Blat(3,ir),3)//" ]*2pi"
      enddo
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
      use linalg, only : det, inv_sym, cross_product
      implicit none
      !
      real(8),intent(in)                    :: Rinput(3,3)
      integer                               :: ir
      !
      !
      if(verbose)write(*,"(A)") "---- set_lattice"
      !
      !
      Rlat(:,1) = Rinput(:,1)
      Rlat(:,2) = Rinput(:,2)
      Rlat(:,3) = Rinput(:,3)
      vol = dot_product(cross_product(Rlat(:,1),Rlat(:,2)),Rlat(:,3))
      if(verbose)write(*,"(A,F)")"     Unit cell volume: ",vol
      Blat(:,1) = cross_product(Rlat(:,2),Rlat(:,3))/vol
      Blat(:,2) = cross_product(Rlat(:,3),Rlat(:,1))/vol
      Blat(:,3) = cross_product(Rlat(:,1),Rlat(:,2))/vol
      !
      write(*,"(A)")new_line("A")//"     Unit cell vectors: "
      do ir=1,3
         write(*,"(A)")"     R_"//str(ir)//": [ "//str(Rlat(1,ir),3)//" , "//str(Rlat(2,ir),3)//" , "//str(Rlat(3,ir),3)//" ]"
      enddo
      write(*,"(A)")new_line("A")//"     Reciprocal lattice vectors: "
      do ir=1,3
         write(*,"(A)")"     B_"//str(ir)//": [ "//str(Blat(1,ir),3)//" , "//str(Blat(2,ir),3)//" , "//str(Blat(3,ir),3)//" ]*2pi"
      enddo
      !
      Lat_stored=.true.
      !
   end subroutine set_lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Access private-protected variables
   !---------------------------------------------------------------------------!
   subroutine get_Rlat(Rlat_out)
      implicit none
      real(8),intent(out)                   :: Rlat_out(3,3)
      if(Lat_stored)then
         Rlat_out = Rlat
      else
         write(*,"(A)")"     Warning: requested lattice vetors but lattice is not stored."
      endif
   end subroutine get_Rlat
   !
   subroutine get_Blat(Blat_out)
      implicit none
      real(8),intent(out)                   :: Blat_out(3,3)
      if(Lat_stored)then
         Blat_out = Blat
      else
         write(*,"(A)")"     Warning: requested reciprocal lattice vetors but lattice is not stored."
      endif
   end subroutine get_Blat


   !---------------------------------------------------------------------------!
   !PURPOSE: Build the k-point mesh
   !---------------------------------------------------------------------------!
   subroutine build_kpt(Nkpt3,kpt,pathOUTPUT)
      !
      use utils_misc
      implicit none
      !
      integer,dimension(3),intent(in)       :: Nkpt3
      real(8),allocatable,intent(inout)     :: kpt(:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      integer                               :: ik,Nkpt,unit
      integer                               :: k1,k2,k3
      !
      !
      if(verbose)write(*,"(A)") "---- build_kpt"
      !
      !
      !if(.not.Klib)stop "build_kpt: lthis routine needs an outer library. Recompile the code with the option GT=T."
      if(.not.Lat_stored)stop "build_kpt: lattice vectors are not stored."
      !
      Nkpt = product(Nkpt3)
      if(allocated(kpt))deallocate(kpt)
      allocate(kpt(3,Nkpt));kpt=0d0
      !
      ik=0
      do k1=1,Nkpt3(1)
         do k2=1,Nkpt3(2)
            do k3=1,Nkpt3(3)
               !
               ik=ik+1
               kpt(:,ik) = [ dble(k1-1)/nkpt3(1), dble(k2-1)/nkpt3(2), dble(k3-1)/nkpt3(3) ]
               !
            enddo
         enddo
      enddo
      !
#ifdef _Klib
      kpt=0d0
      call build_kptGT(kpt,Nkpt3,Blat,[0d0,0d0,0d0])
#endif
      !
      if(present(pathOUTPUT))then
         unit = free_unit()
         open(unit,file=reg(pathOUTPUT)//"Kpoints_BZ.DAT",form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Nkpt
            write(unit,"(1I8,6F20.10)")ik,kpt(:,ik),kpt(1,ik)*Blat(:,1)+kpt(2,ik)*Blat(:,2)+kpt(3,ik)*Blat(:,3)
         enddo
         close(unit)
      endif
      !
   end subroutine build_kpt
   !
#ifdef _Klib
   !
   subroutine build_kptGT(kpt,Nkpt3,Blat,shift)
      use kpointGeneration
      use num_types
      use vector_matrix_utilities
      use utils_misc
      implicit none
      real(8),intent(inout)           :: kpt(:,:)
      integer,intent(in)              :: Nkpt3(3)
      real(8),intent(in)              :: Blat(:,:)
      real(8),intent(in)              :: shift(3)
      !
      real(dp)                        :: B(3,3),H(3,3)
      real(dp)                        :: K(3,3),Hinv(3,3)
      real(dp)                        :: reps,aeps
      real(dp),pointer                :: klist(:,:)
      integer                         :: ik
      !
      !
      if(verbose)write(*,"(A)") "---- build_kptGT"
      !
      !
      call assert_shape(kpt,[3,Nkpt3(1)*Nkpt3(2)*Nkpt3(3)],"build_kptGT","kpt")
      !
      ! Finite precision tolerance (same as default value)
      reps = 1e-8_dp
      aeps = 1e-10_dp
      !
      ! Reciprocal lattice vectors
      B(:,1) = Blat(:,1)
      B(:,2) = Blat(:,2)
      B(:,3) = Blat(:,3)
      !
      ! HNF Matrix. See arXiv:1009.4826 and arXiv:0804.3544v1
      H(:,1) = [ Nkpt3(1), 0, 0 ]
      H(:,2) = [ 0, Nkpt3(2), 0 ]
      H(:,3) = [ 0, 0, Nkpt3(3) ]
      !
      ! Columns of K are the grid generating vectors.
      call matrix_inverse(real(H,dp), Hinv, eps_=aeps)
      K = matmul(Blat,Hinv)
      !
      ! create the kpt list
      call generateFullKpointList(K, B, real(shift,dp), klist, reps_=reps,aeps_=aeps)
      do ik = 1,determinant(H)
         kpt(:,ik) = klist(ik,:)
         if(verbose) write(*,"(5X,1I5,3F15.7)")ik,kpt(:,ik)
      enddo
      !
   end subroutine build_kptGT
   !
#endif


   !---------------------------------------------------------------------------!
   !PURPOSE: Read XEPS.DAT file
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
      iq_gamma = find_vec([0d0,0d0,0d0],kpt_xeps,eps)
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
            !Hk(iwan1,iwan1,ik) = dcmplx(dreal(Hk(iwan1,iwan1,ik)),0d0)
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
      if(present(iq_gamma))iq_gamma = find_vec([0d0,0d0,0d0],kpt,eps)
      write(*,"(A,I4)")"     Gamma point index: ",iq_gamma
      Hk_stored=.true.
      call read_lattice(reg(pathINPUT))
      !
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"Kpoints_BZ.DAT",form="formatted",status="unknown",position="rewind",action="write")
      do ik=1,Nkpt
         write(unit,"(1I8,6F20.10)")ik,kpt(:,ik),kpt(1,ik)*Blat(:,1)+kpt(2,ik)*Blat(:,2)+kpt(3,ik)*Blat(:,3)
      enddo
      close(unit)
      !
   end subroutine read_Hk


   !---------------------------------------------------------------------------!
   !PURPOSE: Build the Hamiltonian and kpoints from user-given parameters
   !---------------------------------------------------------------------------!
   subroutine build_Hk(Rinput,Norb,hopping,Nkpt3,alphaHk,readHr,Hetero,Hk,kpt,Ek,Zk,Hloc,iq_gamma,pathOUTPUT)
      !
      use utils_misc
      use parameters, only : Heterostructures !WHY IS THIS WORKING?
      use linalg, only : zeye, diagonal, diag, eigh, dag
      implicit none
      !
      real(8),intent(in)                    :: Rinput(3,3)
      integer,intent(in)                    :: Norb
      real(8),intent(in)                    :: hopping(:)
      integer,intent(in)                    :: Nkpt3(3)
      real(8),intent(in)                    :: alphaHk
      logical,intent(in)                    :: readHr
      type(Heterostructures),intent(inout)  :: Hetero
      complex(8),allocatable,intent(out)    :: Hk(:,:,:)
      real(8),allocatable,intent(out)       :: kpt(:,:)
      real(8),allocatable,intent(out)       :: Ek(:,:)
      complex(8),allocatable,intent(out)    :: Zk(:,:,:)
      complex(8),allocatable,intent(out)    :: Hloc(:,:)
      integer,intent(out),optional          :: iq_gamma
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      !User
      integer                               :: unit,Nkpt
      integer                               :: iwan1,iwan2,ik
      integer                               :: Trange,idist,iwig
      !W90
      integer,parameter                     :: W90NumCol=15
      integer                               :: Num_wann,Nrpts
      integer                               :: Qst,Rst,i,j,ir
      integer                               :: nx,ny,nz
      integer,allocatable                   :: Ndegen(:)
      real(8)                               :: ReHr,ImHr
      character(len=256)                    :: path
      logical                               :: filexists,Tcond
      !Hetero
      integer                               :: isite,Nsite,na,nb
      integer                               :: tzl,tzr,ilayer
      logical,allocatable                   :: inHomo(:)
      real(8)                               :: tzRatio,angle,Rvec(3)
      real(8),allocatable                   :: Rsorted(:)
      integer,allocatable                   :: Rorder(:),itz(:)
      complex(8),allocatable                :: Hr(:,:,:),Hk_single(:,:,:),Hk_single_offdiag(:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- build_Hk"
      !
      !
      if(readHr.and.(.not.present(pathOUTPUT))) stop "build_Hk: reading of Hr.DAT requested but missing path."
      !
      Nkpt = Nkpt3(1)*Nkpt3(2)*Nkpt3(3)
      call assert_shape(hopping,[Norb],"build_Hk","hopping")
      !
      if(allocated(Hk))deallocate(Hk)
      allocate(Hk(Norb,Norb,Nkpt));Hk=czero
      !
      if(allocated(kpt))deallocate(kpt)
      allocate(kpt(3,Nkpt));kpt=0d0
      !
      call set_lattice(Rinput)
      call build_kpt(Nkpt3,kpt,pathOUTPUT=reg(pathOUTPUT))
      !
      !recover the vectors in real space and allocate hopping in real space
      if(.not.Wig_stored)call calc_wignerseiz(Nkpt3)
      allocate(Rsorted(Nwig));Rsorted = radiuswig
      allocate(Rorder(Nwig))
      call sort_array(Rsorted,Rorder)
      allocate(Hr(Norb,Norb,Nwig));Hr=czero
      !
      if(readHr)then
         !
         ! Look for Hk.DAT
         path=reg(pathOUTPUT)//"Hr.DAT"
         call inquireFile(reg(path),filexists)
         !
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
         read(unit,*)                      !skip first line
         read(unit,*) Num_wann !Number of Wannier orbitals
         read(unit,*) Nrpts    !Number of Wigner-Seitz vectors
         !
         if(Num_wann.ne.Norb) stop "build_Hk: number of Wannier orbital in Hr.DAT and model orbital space does not coincide."
         !
         Qst = int(Nrpts/W90NumCol)
         Rst = mod(Nrpts,W90NumCol)
         !
         allocate(Ndegen(Nrpts));Ndegen=0
         do i=1,Qst
            read(unit,*)(Ndegen(j+(i-1)*W90NumCol),j=1,W90NumCol)
         enddo
         if(Rst.ne.0)read(unit,*)(Ndegen(j+Qst*W90NumCol),j=1,Rst)
         !
         !Read W90 TB hoppings in real space. Assumed paramagnetic
         do ir=1,Nrpts
            do i=1,Num_wann
               do j=1,Num_wann
                  !
                  read(unit,*) nx, ny, nz, iwan1, iwan2, ReHr, ImHr
                  !
                  iwig = find_vec([nx,ny,nz],Nvecwig)
                  !
                  Hr(iwan1,iwan2,iwig) = dcmplx(ReHr,ImHr)/Ndegen(ir)
                  !nrdegwig(iwig) = Ndegen(ir) <-- this would mess-up things in the FT
                  !
               enddo
            enddo
         enddo
         close(unit)
         deallocate(Ndegen)
         !
      else
         !
         !User-provided hopping is only nearest neighbor by now
         Trange=1
         !
         !loop over the sorted Wigner Seiz positions
         idist=1
         loopwigD:do iwig=1,Nwig
            !
            !setting the local energy
            if(Rsorted(Rorder(iwig)).eq.0d0)then
               if(Rorder(iwig).ne.wig0)stop "build_Hk: wrong index of R=0 vector."
               cycle
            endif
            !
            !increasing range
            if(iwig.gt.2)then
               if((Rsorted(Rorder(iwig))-Rsorted(Rorder(iwig-1))).gt.1e-5) idist=idist+1  !if(Rsorted(Rorder(iwig)).gt.Rsorted(Rorder(iwig-1))) idist=idist+1
               if(idist.gt.Trange) exit loopwigD
            endif
            !
            !setting matrix element
            do iwan1=1,Norb
               Hr(iwan1,iwan1,Rorder(iwig)) = -dcmplx(hopping(iwan1),0d0)
            enddo
            !
         enddo loopwigD
         !
      endif
      !
      if(verbose)then
         if(readHr)then
            write(*,'(1A)')        "     H_W90:"
            write(*,'(A,I6)')      "     Number of Wannier functions:   ",Num_wann
            write(*,'(A,I6)')      "     Number of Wigner-Seitz vectors:",Nrpts
            write(*,'(A,I6,A,I6)') "     Deg rows:",Qst," N last row   :",Rst
         endif
         write(*,'(1A)')"     Real-space hopping elements:"
         write(*,"(A6,3A12,1A4)") "  i  ","  Ri  ","  H(Ri)  "," [n1,n2,n3] "," Ndeg "
         do iwig=1,Nwig
            write(*,"(1I6,2F12.4,5I4)")Rorder(iwig),Rsorted(Rorder(iwig)),real(Hr(1,1,Rorder(iwig))),Nvecwig(:,Rorder(iwig)),nrdegwig(Rorder(iwig))
         enddo
      endif
      !
      !FT Hr-->Hk
      call wannier_R2K(Nkpt3,kpt,Hr,Hk)
      deallocate(Hr)
      !
      do ik=1,nkpt
         do iwan1=1,Norb
            Hk(iwan1,iwan1,ik) = dcmplx(dreal(Hk(iwan1,iwan1,ik)),0d0)
         enddo
         if(Norb.gt.1)call check_Hermiticity(Hk(:,:,ik),eps)
      enddo
      !
      !Build up the Heterostructure Hamiltonian
      Nsite = 1
      if(Hetero%status)then
         !
         !this should be already been checked in input_vars
         Nsite = Hetero%Explicit(2)-Hetero%Explicit(1)+1
         !
         !Setting up off-diagonal dispersion if requested
         if(Hetero%offDiagEk)then
            !
            allocate(Hk_single_offdiag(Norb,Norb,Nkpt));Hk_single_offdiag=czero
            allocate(Hr(Norb,Norb,Nwig));Hr=czero
            !
            !User-provided hopping is only nearest neighbor by now
            Trange=1
            !
            !loop over the sorted Wigner Seiz positions
            idist=1
            loopwigOD:do iwig=1,Nwig
               !
               !setting the local energy
               if(Rsorted(Rorder(iwig)).eq.0d0)then
                  if(Rorder(iwig).ne.wig0)stop "build_Hk: wrong index of R=0 vector."
                  cycle
               endif
               !
               !increasing range
               if(iwig.gt.2)then
                  if((Rsorted(Rorder(iwig))-Rsorted(Rorder(iwig-1))).gt.1e-5) idist=idist+1
                  if(idist.gt.Trange) exit loopwigOD
               endif
               !
               !setting matrix element
               !PROJECT SPECIFIC (TaS2)>>>
               !do iwan1=1,Norb
               !   Hr(iwan1,iwan1,Rorder(iwig)) = dcmplx(1d0,0d0)
               !enddo
               Rvec = Nvecwig(1,Rorder(iwig))*Rlat(:,1) + Nvecwig(2,Rorder(iwig))*Rlat(:,2) + Nvecwig(3,Rorder(iwig))*Rlat(:,3)
               angle = atan2(Rvec(2),Rvec(1))
               if(angle.lt.0d0) angle = angle + 2d0*pi
               Tcond = (mod(nint(angle*180/pi)/60,2)-1) .eq. 0
               if(Tcond)then
                  !write(*,*)angle,angle*180/pi,nint(angle*180/pi),mod(nint(angle*180/pi)/60,2),(mod(nint(angle*180/pi)/60,2)-1)
                  !write(*,*)Nvecwig(:,Rorder(iwig))
                  do iwan1=1,Norb
                     Hr(iwan1,iwan1,Rorder(iwig)) = dcmplx(1d0,0d0)
                  enddo
               endif
               !>>>PROJECT SPECIFIC (TaS2)
               !
            enddo loopwigOD
            !
            call wannier_R2K(Nkpt3,kpt,Hr,Hk_single_offdiag)
            deallocate(Hr)
            !
         endif
         !
         !Setting up the out-of-plane hopping array
         tzl = 0 ; tzr = 0
         if(Hetero%Explicit(1).ne.1) tzl = 1              ! hopping to the left potential
         if(Hetero%Explicit(2).ne.Hetero%Nslab) tzr = 1   ! hopping to the right potential
         allocate(Hetero%tz(Norb,Norb,(Hetero%Explicit(1)-tzl):(Hetero%Explicit(2)-1+tzr)));Hetero%tz=czero
         allocate(inHomo((Hetero%Explicit(1)-tzl):(Hetero%Explicit(2)-1+tzr)));inHomo=.false.
         write(*,"(A)")new_line("A")//"     Hetero:"
         do ilayer = Hetero%Explicit(1)-tzl,Hetero%Explicit(2)-1+tzr
            !
            inHomo(ilayer) = (Hetero%NtzExplicit.gt.0) !.and. any(Hetero%ExplicitTzPos.eq.ilayer)
            if(inHomo(ilayer)) inHomo(ilayer) = inHomo(ilayer) .and. any(Hetero%ExplicitTzPos.eq.ilayer)
            !
            tzRatio = 1d0
            if(inHomo(ilayer))then
               allocate(itz(Hetero%NtzExplicit));itz=0
               itz = findloc(Hetero%ExplicitTzPos,value=ilayer)
               if(itz(1).eq.0) stop "build_Hk: something wrong with the Hetero%ExplicitTzPos"
               tzRatio = Hetero%ExplicitTzRatios(itz(1))
               deallocate(itz)
            else
               tzRatio = Hetero%GlobalTzRatio
            endif
            !
            Hetero%tz(:,:,ilayer) = diag(hopping)*tzRatio
            write(*,"(A,F)")"     tz/tplane ["//str(ilayer)//"-"//str(ilayer+1)//"]:",tzRatio
            !
         enddo
         !
         !Setting up multi-site H(k)
         allocate(Hk_single(Norb,Norb,Nkpt));Hk_single=czero
         Hk_single = Hk
         deallocate(Hk)
         allocate(Hk(Norb*Nsite,Norb*Nsite,Nkpt));Hk=czero
         !
         !adding non-diaongonal part
         do isite=1,Nsite
            !
            !Index of the layer inside the slab - Needed because Hetero%tz has a different indexing
            ilayer = Hetero%Explicit(1) + (isite-1)
            !
            !In-plane orbital block
            na = 1+(isite-1)*Norb
            nb = isite*Norb
            !
            !In-plane Hk
            Hk(na:nb,na:nb,:) = Hk_single
            !
            !Out-of-plane hopping
            if(isite.ne.Nsite)then
               !
               !PROJECT SPECIFIC (TaS2)>>>
               !if(Hetero%offDiagEk)then
               if(Hetero%offDiagEk.and.(.not.inHomo(ilayer)))then
               !>>>PROJECT SPECIFIC (TaS2)
                  do ik=1,Nkpt
                     Hk(na:nb,na+Norb:nb+Norb,ik) = matmul(Hetero%tz(:,:,ilayer),Hk_single_offdiag(:,:,ik))
                     Hk(na+Norb:nb+Norb,na:nb,ik) = dag(Hk(na:nb,na+Norb:nb+Norb,ik))
                  enddo
               else
                  do ik=1,Nkpt
                     Hk(na:nb,na+Norb:nb+Norb,ik) = Hetero%tz(:,:,ilayer)
                     Hk(na+Norb:nb+Norb,na:nb,ik) = dag(Hk(na:nb,na+Norb:nb+Norb,ik))
                  enddo
               endif
               !
            endif
            !
         enddo
         deallocate(Hk_single,inHomo)
         !
      endif
      deallocate(Rorder,Rsorted)
      if(Hetero%offDiagEk)deallocate(Hk_single_offdiag)
      !
      Hk = Hk*alphaHk
      !
      if(allocated(Ek))deallocate(Ek)
      if(allocated(Zk))deallocate(Zk)
      if(allocated(Hloc))deallocate(Hloc)
      allocate(Ek(Norb*Nsite,Nkpt));Ek=0d0
      allocate(Zk(Norb*Nsite,Norb*Nsite,Nkpt));Zk=czero
      allocate(Hloc(Norb*Nsite,Norb*Nsite));Hloc=czero
      !
      do ik=1,Nkpt
         !
         call check_Hermiticity(Hk(:,:,ik),eps)
         !
         Ek(:,ik) = 0d0
         Zk(:,:,ik) = Hk(:,:,ik)
         call eigh(Zk(:,:,ik),Ek(:,ik))
         !
      enddo
      Hloc = sum(Hk,dim=3)/Nkpt
      !
      if(present(iq_gamma))iq_gamma = find_vec([0d0,0d0,0d0],kpt,eps)
      write(*,"(A,I4)")"     Gamma point index: ",iq_gamma
      Hk_stored=.true.
      !
      if(present(pathOUTPUT))then
         !
         unit = free_unit()
         open(unit,file=reg(pathOUTPUT)//"Hk.DAT",form="formatted",status="unknown",position="rewind",action="write")
         write(unit,("(3I10)")) 1,Nkpt,Norb
         do ik=1,Nkpt
            write(unit,("(3F14.8)")) kpt(:,ik)
            do iwan1=1,Norb*Nsite
               do iwan2=1,Norb*Nsite
                  write(unit,("(2I4,2E20.12)")) iwan1,iwan2,dreal(Hk(iwan1,iwan2,ik)),dimag(Hk(iwan1,iwan2,ik))
               enddo
            enddo
         enddo
         !
      endif
      !
   end subroutine build_Hk


   !---------------------------------------------------------------------------!
   !PURPOSE: Generate K-points in the irreducible BZ
   !         symkpt(1:nkpt) contains symmetry operation indices that transforms
   !         a k-point in IBZ to the current k-point: R(symkpt(ik)) * kptp(ik) = kpt(ik)
   !         kptsym(1:nkpt,nsym): k-point indices after performing a symmetry operation
   !         to each k: R(isym)*kpt(ik) = kpt(kptsym(ik,isym)) + gkptsym(1:3,ik,isym)
   !---------------------------------------------------------------------------!
   subroutine calc_irredBZ(pathINPUT,nkpt3,nkpti,kptp,pkpt,nkstar,kpt_out,store)
      !
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      integer,intent(in)                    :: nkpt3(3)
      integer,intent(out)                   :: nkpti
      integer,allocatable,intent(out)       :: kptp(:)
      integer,allocatable,intent(out)       :: pkpt(:,:,:)
      real(8),allocatable,intent(out)       :: nkstar(:)
      real(8),allocatable,intent(out),optional :: kpt_out(:,:)
      logical,intent(in),optional           :: store
      !
      type(symtype),allocatable             :: sym(:)
      integer                               :: Nsym,Nkpt,ik,k1,k2,k3
      integer                               :: i,j,k,l,idum
      integer                               :: iarr2(3)
      real(8)                               :: scale,rarr(3),kvec(3)
      real(8),allocatable                   :: kpt(:,:),kptw(:)
      integer,allocatable                   :: symkpt(:),iarr(:),kptsym(:,:)!,gkptsym(:,:,:)
      integer                               :: unit
      logical                               :: store_,filexists,cond1,cond2,err
      !
      !
      if(verbose)write(*,"(A)") "---- calc_irredBZ"
      !
      !
      store_=.false.
      if(present(store))store_=store
      !
      !reading symmetry file from SPEX
      call inquireFile(reg(pathINPUT)//"sym.DAT",filexists)
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"sym.DAT",form="formatted",action="read",position="rewind")
      read(unit,*) Nsym
      allocate(sym(Nsym))
      do i=1,Nsym
         read(unit,*)
         read(unit,*) ( sym(i)%rot(j,:),sym(i)%transl(j), j=1,3 )
      enddo
      close(unit)
      !
      do i=1,Nsym
         !
         if(any(sym(i)%transl.ne.0)) then
            sym(i)%symmor = .false.
         else
            sym(i)%symmor = .true.
         endif
         !
         sym(i)%inv = 0
         do j=1,Nsym
            !
            cond1 = all(matmul(sym(i)%rot,sym(j)%rot) .eq. reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/)))
            cond2 = all(modulo1r(matmul(sym(i)%rot,sym(j)%transl)+sym(i)%transl).lt.1d-10)
            !
            if(cond1.and.cond2)then
               if(sym(i)%inv.ne.0) stop "calc_irredBZ: inverse operation already defined."
               sym(i)%inv = j
               !sym(i)%rrot = transpose_int ( sym(j)%rot ) ! temporary fix for ifc
               sym(i)%rrot = transpose( sym(j)%rot )
            endif
            !
         enddo
         if(sym(i)%inv.eq.0) stop "calc_irredBZ: inverse operation not found."
         !
      enddo
      !
      Nkpt = product(Nkpt3)
      !
      !internal K-point arrays
      allocate(kpt(3,Nkpt));kpt=0d0
      allocate(symkpt(Nkpt));symkpt=0
      !
      !output K-point arrays
      if(allocated(kptp))deallocate(kptp)
      if(allocated(pkpt))deallocate(pkpt)
      allocate(kptp(Nkpt));kptp=0
      allocate(pkpt(Nkpt3(1)+1,Nkpt3(2)+1,Nkpt3(3)+1));pkpt=0
      !
      !fill in kpt exactly as in build_kpt
      ik=0
      do k1 = 1,nkpt3(1)
         do k2 = 1,nkpt3(2)
            do k3 = 1,nkpt3(3)
               ik=ik+1
               kpt(:,ik) = [ dble(k1-1)/nkpt3(1), dble(k2-1)/nkpt3(2), dble(k3-1)/nkpt3(3) ]
            enddo
         enddo
      enddo
      !
      !fill in the kpt pointer
      do ik=1,Nkpt
        iarr2 = nint( kpt(:,ik) * nkpt3 ) + 1
        pkpt(iarr2(1),iarr2(2),iarr2(3)) = ik
      enddo
      pkpt(nkpt3(1)+1,    :     ,    :     ) = pkpt(1,:,:)
      pkpt(    :     ,nkpt3(2)+1,    :     ) = pkpt(:,1,:)
      pkpt(    :     ,    :     ,nkpt3(3)+1) = pkpt(:,:,1)
      if(any(pkpt.eq.0)) stop "calc_irredBZ: Definition of pkpt-pointer failed."
      !
      !Very obscure part
      allocate(iarr(Nkpt));iarr=1
      err = .false.
      do i=1,Nkpt
         !
         if(iarr(i).eq.0) cycle
         kptp(i) = i
         symkpt(i) = 1
         !
         do k=2,Nsym
            !
            l = k
            rarr  = matmul(sym(l)%rrot,kpt(:,i)) * nkpt3
            iarr2 = nint(rarr)
            if(any(abs(iarr2-rarr).gt.1d-10)) then
               write(0,"(A)") "calc_irredBZ: Symmetry operation "//str(l)//" incompatible with k-point set."
               err = .true.
            endif
            iarr2 = modulo(iarr2,nkpt3) + 1
            !
            if(any(iarr2.gt.nkpt3)) stop "calc_irredBZ: pointer indices exceed pointer dimensions. (bug?)"
            j = pkpt(iarr2(1),iarr2(2),iarr2(3))
            !
            if(j.eq.0) stop "calc_irredBZ: k-point index is zero (bug?)"
            if((iarr(j).eq.0).or.(j.eq.i)) cycle
            iarr(j) = 0
            kptp(j) = i
            symkpt(j) = l
         enddo
      enddo
      !
      if(err) then
        stop "calc_irredBZ: Some symmetry operations are incompatible with k-point set."
      endif
      !
      i=0
      do ik=1,nkpt
         if(iarr(ik).eq.1) then
            i = i + 1
            iarr(ik) = i
         endif
      enddo
      nkpti = i
      !
      do ik=1,nkpt
         if(iarr(ik).eq.0) then
            i = i + 1
            iarr(ik) = i
         endif
      enddo
      kpt(:,iarr)  = kpt
      kptp         = iarr(kptp)
      kptp(iarr)   = kptp
      symkpt(iarr) = symkpt
      !
      do i=1,nkpt3(1)+1
         do j=1,nkpt3(2)+1
            do k=1,nkpt3(3)+1
               pkpt(i,j,k) = iarr(pkpt(i,j,k))
            enddo
         enddo
      enddo
      !
      !Define k-point mapping wrt symmetry operations (->kptsym,gkptsym)
      !allocate(gkptsym(3,Nkpt,Nsym));gkptsym=0  !not used
      allocate(kptsym(Nkpt,Nsym));kptsym=0
      do ik=1,Nkpt
         do k=1,Nsym
            !
            if(ik.le.nkpt) then
               rarr = matmul(sym(k)%rrot,kpt(:,ik))
               iarr(1:3) = nint ( modulo1(rarr) * nkpt3 ) + 1
               idum = pkpt(iarr(1),iarr(2),iarr(3))
               kvec = kpt(:,idum)
            else
               stop "calc_irredBZ: Bug?"
            endif
            !
            if(idum.eq.0) then
              write(*,"(A,I3)") "calc_irredBZ: K-point mapping failed for symmetry operation",k
              write(*,"(A)")    "              (K-point mesh and symmetry consistent?)"
              stop              "calc_irredBZ: K-point mapping failed."
            endif
            !
            iarr(1:3) = nint ( rarr - kvec )
            kptsym(ik,k) = idum
            !gkptsym(:,ik,k) = iarr(1:3) !not used
            !
            if(any(abs(rarr-iarr(1:3)-kvec).gt.1d-8)) then
               write(*,"(A,I3)") "calc_irredBZ: K-point mapping failed for symmetry operation",k
               write(*,"(A)")    "              (K-point mesh and symmetry consistent?)"
               stop              "calc_irredBZ: K-point mapping failed."
            endif
            !
         enddo ! k
      enddo ! ik
      !
      scale = kgv(nkpt3,3)
      !
      if(allocated(nkstar))deallocate(nkstar)
      allocate(nkstar(Nkpti));nkstar=0d0
      allocate(kptw(Nkpti));kptw=0d0
      do ik=1,Nkpti
        iarr=0
        do j=1,Nsym
           iarr(kptsym(ik,j)) = 1
        enddo
        kptw(ik) = dble(sum(iarr))/dble(Nkpt)
        nkstar(ik) = kptw(ik)*Nkpt
      enddo
      !
      write(*,"(A,I)") "     Total number of k-points: ",Nkpt
      write(*,"(A,I)") "     Number of k-points in IBZ: ",Nkpti
      if(store_)then
         unit = free_unit()
         open(unit,file=reg(pathINPUT)//"Kpoints_BZirred.DAT",form="formatted",status="unknown",position="rewind",action="write")
         write(unit,"(I5,F20.10)") Nkpti,scale
         do ik=1,Nkpti
            write(unit,"(1I8,6F20.10)")ik,kpt(:,ik)*scale,count(iarr.eq.1)*1d0,kptw(ik),nkstar(ik)
         enddo
         close(unit)
         write(*,"(A)") "     Written to file: "//reg(pathINPUT)//"Kpoints_BZirred.DAT"
      endif
      !
      if(present(kpt_out))kpt_out=kpt
      deallocate(sym,kpt,symkpt,iarr,kptw)
      !
      !
   contains
      !
      !
      !
      integer function kgv(iarr,n)
         implicit none
         integer,intent(in)                 :: n,iarr(n)
         logical                            :: lprim(2:maxval(iarr))
         integer,allocatable                :: prim(:),expo(:)
         integer                            :: nprim,marr
         integer                            :: i,j,ia,k
         !
         !Determine prime numbers
         marr=maxval(iarr)
         lprim=.true.
         do i=2,marr
            j = 2
            do while (i*j.le.marr)
               lprim(i*j) = .false.
               j = j + 1
            enddo
         enddo
         nprim=count(lprim)
         allocate(prim(nprim),expo(nprim))
         j=0
         do i=2,marr
            if(lprim(i))then
               j = j + 1
               prim(j) = i
            endif
         enddo
         !
         !Determine least common multiple
         expo=0
         do i=1,n
            ia = iarr(i)
            if(ia.eq.0) cycle
            do j=1,nprim
               k=0
               do while(ia/prim(j)*prim(j).eq.ia)
                  k  = k + 1
                  ia = ia / prim(j)
               enddo
               expo(j) = max(expo(j),k)
            enddo
         enddo
         kgv=1
         do j=1,nprim
            kgv = kgv * prim(j)**expo(j)
         enddo
         deallocate(prim,expo)
      end function kgv
      !
      !Replaces the function modulo(kpoint,1d0) for a kpoint in kpt(:,ikpt).
      function modulo1(kpoint)
         implicit none
         real(8),intent(in)                 :: kpoint(3)
         real(8)                            :: modulo1(3)
         integer                            :: help(3)
         modulo1 = kpoint*nkpt3
         help = nint(modulo1)
         if(any(abs(help-modulo1).gt.1d-10)) then
            write(*,"(A)") "     modulo1: argument ("//str(kpoint(1),5)//","//str(kpoint(2),5)//","//str(kpoint(3),5)//") is not an element of the k-point set."
            stop "modulo1: argument not an element of k-point set."
         endif
         modulo1 = modulo(help,nkpt3)*1d0/nkpt3
      end function modulo1
      !
      !Same for shifted k-point set
      function modulo1r(kpoint)
         implicit none
         real(8),intent(in)                 :: kpoint(3)
         real(8)                            :: modulo1r(3)
         integer                            :: i
         modulo1r = modulo( kpoint , 1d0 )
         do i=1,3
          if(abs(1-abs(modulo1r(i))).lt.1d-13) modulo1r(i) = 0d0
        enddo
      end function modulo1r
      !
      !
   end subroutine calc_irredBZ


   !---------------------------------------------------------------------------!
   !PURPOSE: Fill up the list likning indexes of the sum and diff of K-points
   !         that was originally meant for cubic systems and the points external
   !         to the 1st BZ recovered by just removing nint(dk) to dk.
   !         For generic lattice the BZ vector is not simply the versor as implied
   !         by nint.
   !---------------------------------------------------------------------------!
   subroutine fill_ksumkdiff(kpt,kptsum,kptdif,Nkpt3,pkpt)
      !
      use utils_misc
      use linalg, only : det3
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
      !k1 + q = k2
      kptsum=0
      do iq=1,Nkpt
         do ik1=1,iq
            !
            !first attempt
            do ik2=1,Nkpt
               dk = kpt(:,ik1) + kpt(:,iq) - kpt(:,ik2)
               dk = dk - nint(dk)
               if(all(abs(dk).lt.eps))then
                  kptsum(ik1,iq)=ik2
                  if (ik1.ne.iq) kptsum(iq,ik1)=kptsum(ik1,iq)
                  exit
               endif
            enddo
            !
            !second attempt with larger tolerance
            if (kptdif(ik1,iq).eq.0)then
               do ik2=1,Nkpt
                  dk = kpt(:,ik1) + kpt(:,iq) - kpt(:,ik2)
                  dk = dk - nint(dk)
                  if(all(abs(dk).lt.1e3*eps))then
                     kptsum(ik1,iq)=ik2
                     if (ik1.ne.iq) kptsum(iq,ik1)=kptsum(ik1,iq)
                     exit
                  endif
               enddo
            endif
            !
            !missing sum
            if (kptsum(ik1,iq).eq.0) stop "fill_ksumkdiff: kptsum failed."
            !
         enddo ! ik1
      enddo ! iq
      !
      !k1 - q = k2
      kptdif=0
      do iq=1,nkpt
         do ik1=1,nkpt
            !
            !first attempt
            do ik2=1,Nkpt
               dk = kpt(:,ik1) - kpt(:,iq) - kpt(:,ik2)
               dk = dk - nint(dk)
               if(all(abs(dk).lt.eps))then
                  kptdif(ik1,iq)=ik2
                  exit
               endif
            enddo
            !
            !second attempt with larger tolerance
            if (kptdif(ik1,iq).eq.0)then
               do ik2=1,Nkpt
                  dk = kpt(:,ik1) - kpt(:,iq) - kpt(:,ik2)
                  dk = dk - nint(dk)
                  if(all(abs(dk).lt.1e3*eps))then
                     kptdif(ik1,iq)=ik2
                     exit
                  endif
               enddo
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
                     dk = kpt(:,ik1) - k
                     dk = dk - nint(dk)
                     if(all(abs(dk).lt.eps))then
                        pkpt_(i1,i2,i3)=ik1
                        exit
                     endif
                  enddo
                  !
                  !second attempt with larger tolerance
                  if (pkpt_(i1,i2,i3).eq.0)then
                     do ik1=1,nkpt
                        dk = kpt(:,ik1) - k
                        dk = dk - nint(dk)
                        if(all(abs(dk).lt.1e3*eps))then
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
      !
   contains
      !
      !
      !
      !Not used
      logical function moduloG(Dk,tol)
         implicit none
         real(8),intent(in)                 :: Dk(3)
         real(8),intent(in)                 :: tol
         real(8)                            :: M1(3,3),M2(3,3),M3(3,3)
         real(8)                            :: n1,n2,n3
         !
         !that"s simply the Cramer rule
         M1 = Blat
         M2 = Blat
         M3 = Blat
         !
         M1(:,1) = Dk
         M2(:,2) = Dk
         M3(:,3) = Dk
         !
         n1 = det3(M1)/det3(Blat)
         n2 = det3(M2)/det3(Blat)
         n3 = det3(M3)/det3(Blat)
         !
         !true if all are equal to 0,+1,-1
         moduloG = (abs(n1*(n1+1)*(n1-1)).lt.tol) .and. &
                   (abs(n2*(n2+1)*(n2-1)).lt.tol) .and. &
                   (abs(n3*(n3+1)*(n3-1)).lt.tol)
         !
      end function moduloG
      !
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
   ! Rlat(3,3) : premitive lattice vectors
   ! nkpt : number of k-points
   ! nkpt3(3) : number of k-points along each direction
   !OUTPUT:
   ! nwig : number of points inside wigner-seitz supercell
   ! rvec(3,2*nkpt) : (INTEGER) lattice vectors
   ! nrdeg(2*nkpt) : degeneracy
   !---------------------------------------------------------------------------!
   subroutine calc_wignerseiz(nkpt3)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3(3)
      !
      integer                               :: Nkpt
      integer                               :: ir1,ir2,ir3,irsc1,irsc2,irsc3
      integer                               :: i,i0,iwig
      integer,parameter                     :: nshell=2
      real(8)                               :: rtmp(3),rtmpsc(3),dr(3)
      real(8),allocatable                   :: dist(:),radiuswig_tmp(:)
      real(8),allocatable                   :: Rvecwig_tmp(:,:)
      integer,allocatable                   :: Nvecwig_tmp(:,:),nrdegwig_tmp(:)
      real(8)                               :: distmin
      !
      !
      if(verbose)write(*,"(A)") "---- calc_wignerseiz"
      if(.not.Lat_stored)stop "calc_wignerseiz: Lattice positions not stored. Either call read_lattice(path) or read_Hk(path,Hk,kpt)."
      !
      Nkpt = nkpt3(1)*nkpt3(2)*nkpt3(3)
      !
      allocate(Nvecwig_tmp(3,100*nkpt));Nvecwig_tmp=0
      allocate(Rvecwig_tmp(3,100*nkpt));Rvecwig_tmp=0d0
      allocate(nrdegwig_tmp(100*nkpt));nrdegwig_tmp=0
      allocate(radiuswig_tmp(100*nkpt));radiuswig_tmp=0d0
      !
      !this i0 corresponds to irsc1=irsc2=irsc3=0
      allocate(dist((2*nshell+1)**3))
      i0=nshell*(1+(2*nshell+1)*(1+(2*nshell+1)))+1
      !
      iwig=0
      do ir1=-nkpt3(1),+nkpt3(1)
         do ir2=-nkpt3(2),+nkpt3(2)
            do ir3=-nkpt3(3),+nkpt3(3)
               !
               rtmp(:)=matmul(Rlat,(/ir1,ir2,ir3/))
               i=0
               !
               do irsc1=-nshell,+nshell
                  do irsc2=-nshell,+nshell
                     do irsc3=-nshell,+nshell
                        !
                        i=i+1
                        rtmpsc(:)=matmul(Rlat,(/nkpt3(1)*irsc1,nkpt3(2)*irsc2,nkpt3(3)*irsc3/))
                        dr(:)=rtmp(:)-rtmpsc(:)
                        dist(i)=sum(dr(:)**2)
                        if((i.eq.i0).and.(.not.all([irsc1,irsc2,irsc3].eq.[0,0,0])))stop "calc_wignerseiz: wrong index of R=0 vector."
                        !
                     enddo ! irsc3
                  enddo ! irsc2
               enddo ! irsc1
               !
               distmin=minval(dist(:))
               if (abs(distmin-dist(i0)).le.epsWig) then
                  !
                  iwig=iwig+1
                  if (iwig.gt.100*nkpt) stop "calc_wignerseiz: iwig>100*nkpt."
                  !
                  Nvecwig_tmp(:,iwig)=(/ir1,ir2,ir3/)
                  Rvecwig_tmp(:,iwig)=matmul(Rlat,(/ir1,ir2,ir3/))
                  nrdegwig_tmp(iwig)=count(abs(distmin-dist(:)).le.epsWig)
                  radiuswig_tmp(iwig)=sqrt(dble(dot_product(rtmp,rtmp)))
                  if(all([ir1,ir2,ir3].eq.[0,0,0]))wig0=iwig
                  !
                  !if(verbose)then
                  !   write(*,"(A12,I10,A5,I4)") "     iwig:",iwig,"deg:",nrdegwig_tmp(iwig)
                  !   write(*,"(A,3I8)") "     Nvecwig:",Nvecwig_tmp(:,iwig)
                  !   write(*,"(A,3F8.2)") "     Rvecwig:",Rvecwig_tmp(:,iwig)
                  !endif
                  !
               endif
               !
            enddo
         enddo
      enddo
      deallocate(dist)
      !
      if (abs(sum(1d0/nrdegwig_tmp(1:iwig))-nkpt).gt.epsWig) then
         write(*,"(A,F)") "Error: sum(1/nrdeg(:))=",sum(1d0/nrdegwig_tmp(1:iwig))
         stop "calc_wignerseiz: nrdeg failed."
      endif
      !
      if(allocated(radiuswig))deallocate(radiuswig)
      if(allocated(Nvecwig))  deallocate(Nvecwig)
      if(allocated(Rvecwig))  deallocate(Rvecwig)
      if(allocated(nrdegwig)) deallocate(nrdegwig)
      !
      !public,protected
      Nwig = iwig
      allocate(radiuswig(Nwig)) ; radiuswig=radiuswig_tmp(1:Nwig)
      allocate(Nvecwig(3,Nwig)) ; Nvecwig=Nvecwig_tmp(:,1:Nwig)
      allocate(Rvecwig(3,Nwig)) ; Rvecwig=Rvecwig_tmp(:,1:Nwig)
      allocate(nrdegwig(Nwig))  ; nrdegwig=nrdegwig_tmp(1:Nwig)
      !
      Wig_stored=.true.
      !
      deallocate(Nvecwig_tmp,Rvecwig_tmp,nrdegwig_tmp,radiuswig_tmp)
      !
   end subroutine calc_wignerseiz
   !
   subroutine clear_wignerseiz
      implicit none
      !
      Nwig = 0
      if(allocated(radiuswig))deallocate(radiuswig)
      if(allocated(Nvecwig))  deallocate(Nvecwig)
      if(allocated(Rvecwig))  deallocate(Rvecwig)
      if(allocated(nrdegwig)) deallocate(nrdegwig)
      Wig_stored=.false.
      !
   end subroutine clear_wignerseiz



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
   subroutine wannierinterpolation_D1_d(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      real(8),intent(in)                    :: mat_orig(:,:)
      real(8),intent(inout)                 :: mat_intp(:,:)
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
      if(verbose)write(*,"(A)") "---- wannierinterpolation_D1_d"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_D1_d","nkpt3_orig")
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_D1_d: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_D1_d: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      call assert_shape(mat_intp,[Nsize1,Nkpt_intp],"wannierinterpolation_D1_d","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nwig));mat_R=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,kpt_orig,Nvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            !
            do ik=1,Nkpt_orig
               !
               kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
               cfac = dcmplx(cos(kR),-sin(kR))
               !
               mat_R(i1,ir) = mat_R(i1,ir) + mat_orig(i1,ik)*cfac
               !
            enddo ! ik
            !
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            !
            do ir=1,Nwig
               !
               kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
               cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
               !
               mat_intp(i1,ik) = mat_intp(i1,ik) + mat_R(i1,ir)*cfac
               !
            enddo ! ir
            !
         enddo ! i1
      enddo ! ik
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_D1_d
   !
   subroutine wannierinterpolation_D1_z(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
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
      if(verbose)write(*,"(A)") "---- wannierinterpolation_D1_z"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_D1_z","nkpt3_orig")
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_D1_z: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_D1_z: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      call assert_shape(mat_intp,[Nsize1,Nkpt_intp],"wannierinterpolation_D1_z","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,kpt_orig,Nvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            !
            do ik=1,Nkpt_orig
               !
               kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
               cfac = dcmplx(cos(kR),-sin(kR))
               !
               mat_R(i1,ir) = mat_R(i1,ir) + mat_orig(i1,ik)*cfac
               !
            enddo ! ik
            !
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            !
            do ir=1,Nwig
               !
               kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
               cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
               !
               mat_intp(i1,ik) = mat_intp(i1,ik) + mat_R(i1,ir)*cfac
               !
            enddo ! ir
            !
         enddo ! i1
      enddo ! ik
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_D1_z
   !
   !
   subroutine wannierinterpolation_D2_d(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      real(8),intent(in)                    :: mat_orig(:,:,:)
      real(8),intent(inout)                 :: mat_intp(:,:,:)
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
      if(verbose)write(*,"(A)") "---- wannierinterpolation_D2_d"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_D2_d","nkpt3_orig")
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_D2_d: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_D2_d: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      Nsize2 = size(mat_orig,dim=2)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nkpt_intp],"wannierinterpolation_D2_d","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nsize2,Nwig));mat_R=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,kpt_orig,Nvecwig,mat_orig,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            do i2=1,Nsize2
               !
               do ik=1,Nkpt_orig
                  !
                  kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
                  cfac = dcmplx(cos(kR),-sin(kR))
                  !
                  mat_R(i1,i2,ir) = mat_R(i1,i2,ir) + mat_orig(i1,i2,ik)*cfac
                  !
               enddo ! ik
               !
            enddo ! i2
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            do i2=1,Nsize2
               !
               do ir=1,Nwig
                  !
                  kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
                  cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                  !
                  mat_intp(i1,i2,ik) = mat_intp(i1,i2,ik) + mat_R(i1,i2,ir)*cfac
                  !
               enddo ! ir
               !
            enddo ! i2
         enddo ! i1
      enddo ! ik
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_D2_d
   !
   subroutine wannierinterpolation_D2_z(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
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
      if(verbose)write(*,"(A)") "---- wannierinterpolation_D2_z"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_D2_z","nkpt3_orig")
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_D2_z: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_D2_z: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      Nsize2 = size(mat_orig,dim=2)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nkpt_intp],"wannierinterpolation_D2_z","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nsize2,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,kpt_orig,Nvecwig,mat_orig,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            do i2=1,Nsize2
               !
               do ik=1,Nkpt_orig
                  !
                  kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
                  cfac = dcmplx(cos(kR),-sin(kR))
                  !
                  mat_R(i1,i2,ir) = mat_R(i1,i2,ir) + mat_orig(i1,i2,ik)*cfac
                  !
               enddo ! ik
               !
            enddo ! i2
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            do i2=1,Nsize2
               !
               do ir=1,Nwig
                  !
                  kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
                  cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                  !
                  mat_intp(i1,i2,ik) = mat_intp(i1,i2,ik) + mat_R(i1,i2,ir)*cfac
                  !
               enddo ! ir
               !
            enddo ! i2
         enddo ! i1
      enddo ! ik
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_D2_z
   !
   !
   subroutine wannierinterpolation_D3_d(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_orig(:,:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      real(8),intent(in)                    :: mat_orig(:,:,:,:)
      real(8),intent(inout)                 :: mat_intp(:,:,:,:)
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
      if(verbose)write(*,"(A)") "---- wannierinterpolation_D3_d"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_D3_d","nkpt3_orig")
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_D3_d: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_D3_d: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      Nsize2 = size(mat_orig,dim=2)
      Nsize3 = size(mat_orig,dim=3)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nsize3,Nkpt_intp],"wannierinterpolation_D3_d","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nsize2,Nsize3,Nwig));mat_R=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,Nsize3,kpt_orig,Nvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            do i2=1,Nsize2
               do i3=1,Nsize3
                  !
                  do ik=1,Nkpt_orig
                     !
                     kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
                     cfac = dcmplx(cos(kR),-sin(kR))
                     !
                     mat_R(i1,i2,i3,ir) = mat_R(i1,i2,i3,ir) + mat_orig(i1,i2,i3,ik)*cfac
                     !
                  enddo ! ik
                  !
               enddo ! i3
            enddo ! i2
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp(:,:,:,:)=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,Nsize3,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            do i2=1,Nsize2
               do i3=1,Nsize3
                  !
                  do ir=1,Nwig
                     !
                     kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
                     cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                     !
                     mat_intp(i1,i2,i3,ik) = mat_intp(i1,i2,i3,ik) + mat_R(i1,i2,i3,ir)*cfac
                     !
                  enddo ! ir
                  !
               enddo ! i3
            enddo ! i2
         enddo ! i1
      enddo ! ik
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_D3_d
   !
   subroutine wannierinterpolation_D3_z(nkpt3_orig,kpt_orig,kpt_intp,mat_orig,mat_intp)
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
      if(verbose)write(*,"(A)") "---- wannierinterpolation_D3_z"
      if(.not.Wig_stored)then
         if(verbose)write(*,"(A)") "     Calculating Wigner Seiz."
         call assert_shape(nkpt3_orig,[3],"wannierinterpolation_D3_z","nkpt3_orig")
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannierinterpolation_D3_z: size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "wannierinterpolation_D3_z: size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_orig,dim=1)
      Nsize2 = size(mat_orig,dim=2)
      Nsize3 = size(mat_orig,dim=3)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nsize3,Nkpt_intp],"wannierinterpolation_D3_z","mat_intp")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      allocate(mat_R(Nsize1,Nsize2,Nsize3,Nwig));mat_R=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,Nsize3,kpt_orig,Nvecwig,mat_orig,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            do i2=1,Nsize2
               do i3=1,Nsize3
                  !
                  do ik=1,Nkpt_orig
                     !
                     kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
                     cfac = dcmplx(cos(kR),-sin(kR))
                     !
                     mat_R(i1,i2,i3,ir) = mat_R(i1,i2,i3,ir) + mat_orig(i1,i2,i3,ik)*cfac
                     !
                  enddo ! ik
                  !
               enddo ! i3
            enddo ! i2
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp(:,:,:,:)=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,Nsize3,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            do i2=1,Nsize2
               do i3=1,Nsize3
                  !
                  do ir=1,Nwig
                     !
                     kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
                     cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                     !
                     mat_intp(i1,i2,i3,ik) = mat_intp(i1,i2,i3,ik) + mat_R(i1,i2,i3,ir)*cfac
                     !
                  enddo ! ir
                  !
               enddo ! i3
            enddo ! i2
         enddo ! i1
      enddo ! ik
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(mat_R)
      !
   end subroutine wannierinterpolation_D3_z


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
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_NN: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Npoints = size(mat_K,dim=3)
      if(size(mat_K,dim=1).ne.size(mat_K,dim=2)) stop "wannier_K2R_NN: mat_K not square."
      Nsize = size(mat_K,dim=1)
      call assert_shape(mat_R_nn,[Nsize,Nsize,Npoints,3],"wannierinterpolation","mat_R_nn")
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Npoints,Nsize,kpt_orig,Nvecwig,mat_K,mat_R_nn),&
      !$OMP PRIVATE(Rx,Ry,Rz,ir2,ir,ik,id,i1,i2,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize
            do i2=1,Nsize
               do id=1,Npoints
                  !
                  Rx = all(Nvecwig(:,ir).eq.[1,0,0])
                  Ry = all(Nvecwig(:,ir).eq.[0,1,0])
                  Rz = all(Nvecwig(:,ir).eq.[0,0,1])
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
                     kR=2*pi*dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
                     cfac=dcmplx(cos(kR),-sin(kR))
                     !
                     mat_R_nn(i1,i2,id,ir2) = mat_R_nn(i1,i2,id,ir2) + mat_K(i1,i2,id,ik)*cfac
                     !
                  enddo ! ik
                  !
               enddo ! id
            enddo ! i2
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_K2R_NN


   !---------------------------------------------------------------------------!
   !PURPOSE: Transforms a K-dependent matrix into Wannier basis
   !Nwig might not be already available in the main program
   !so I"m allocating here the mat_R which does not need to be allocated in
   !the calling routine.
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
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_d1: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_K,dim=1)
      !
      if(allocated(mat_R))deallocate(mat_R)
      allocate(mat_R(Nsize1,Nwig));mat_R=czero
      !
      ! M(R)=\sum_{k} M(k)*exp[-ik*R]
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,kpt_orig,Nvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            !
            do ik=1,Nkpt_orig
               !
               kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
               cfac = dcmplx(cos(kR),-sin(kR))
               !
               mat_R(i1,ir) = mat_R(i1,ir) + mat_K(i1,ik)*cfac
               !
            enddo ! ik
            !
         enddo ! i1
      enddo! ir
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
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_d2: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
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
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,kpt_orig,Nvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            do i2=1,Nsize2
               !
               do ik=1,Nkpt_orig
                  !
                  kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
                  cfac = dcmplx(cos(kR),-sin(kR))
                  !
                  mat_R(i1,i2,ir) = mat_R(i1,i2,ir) + mat_K(i1,i2,ik)*cfac
                  !
               enddo ! ik
               !
            enddo ! i2
         enddo ! i1
      enddo ! ir
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
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "wannier_K2R_d3: size(kpt_orig,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
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
      !$OMP SHARED(Nwig,Nkpt_orig,Nsize1,Nsize2,Nsize3,kpt_orig,Nvecwig,mat_K,mat_R),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do ir=1,Nwig
         do i1=1,Nsize1
            do i2=1,Nsize2
               do i3=1,Nsize3
                  !
                  do ik=1,Nkpt_orig
                     !
                     kR = 2*pi * dot_product(kpt_orig(:,ik),Nvecwig(:,ir))
                     cfac = dcmplx(cos(kR),-sin(kR))
                     !
                     mat_R(i1,i2,i3,ir) = mat_R(i1,i2,i3,ir) + mat_K(i1,i2,i3,ik)*cfac
                     !
                  enddo ! ik
                  !
               enddo ! i3
            enddo ! i2
         enddo ! i1
      enddo ! ir
      !$OMP END DO
      !$OMP END PARALLEL
      mat_R = mat_R/Nkpt_orig
      !
   end subroutine wannier_K2R_d3


   !---------------------------------------------------------------------------!
   !PURPOSE: Transforms matrix in Wannier basis into K-space
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
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "wannier_R2K_d1: size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_R,dim=1)
      call assert_shape(mat_intp,[Nsize1,Nkpt_intp],"wannierinterpolation","mat_intp")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            !
            do ir=1,Nwig
               !
               kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
               cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
               !
               mat_intp(i1,ik) = mat_intp(i1,ik) + mat_R(i1,ir)*cfac
               !
            enddo ! ir
            !
         enddo ! i1
      enddo ! ik
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
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "wannier_R2K_d2: size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
      !
      ! Size checks on Matrices
      Nsize1 = size(mat_R,dim=1)
      Nsize2 = size(mat_R,dim=2)
      call assert_shape(mat_intp,[Nsize1,Nsize2,Nkpt_intp],"wannierinterpolation","mat_intp")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_intp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            do i2=1,Nsize2
               !
               do ir=1,Nwig
                  !
                  kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
                  cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                  !
                  mat_intp(i1,i2,ik) = mat_intp(i1,i2,ik) + mat_R(i1,i2,ir)*cfac
                  !
               enddo ! ir
               !
            enddo ! i2
         enddo ! i1
      enddo ! ik
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
         call calc_wignerseiz(nkpt3_orig)
      endif
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_intp,dim=1).ne.3) stop "wannier_R2K_d3: size(kpt_intp,dim=1).ne.3"
      Nkpt_intp = size(kpt_intp,dim=2)
      !if (Nkpt_orig.ne.size(nrdegwig)/2) stop "nkpt"
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
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,Nsize3,kpt_intp,Nvecwig,mat_intp,mat_R,nrdegwig),&
      !$OMP PRIVATE(ir,ik,i1,i2,i3,kR,cfac)
      !$OMP DO
      do ik=1,Nkpt_intp
         do i1=1,Nsize1
            do i2=1,Nsize2
               do i3=1,Nsize3
                  !
                  do ir=1,Nwig
                     !
                     kR = 2*pi * dot_product(kpt_intp(:,ik),Nvecwig(:,ir))
                     cfac = dcmplx(cos(kR),+sin(kR))/nrdegwig(ir)
                     !
                     mat_intp(i1,i2,i3,ik) = mat_intp(i1,i2,i3,ik) + mat_R(i1,i2,i3,ir)*cfac
                     !
                  enddo ! ir
                  !
               enddo ! i3
            enddo ! i2
         enddo ! i1
      enddo ! ik
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine wannier_R2K_d3


   !---------------------------------------------------------------------------!
   !PURPOSE: Save in the module variables the user-defined high symmetry points
   !---------------------------------------------------------------------------!
   subroutine set_UserPath(UserPath_input)
      implicit none
      real(8),allocatable,intent(in)        :: UserPath_input(:,:)
      if(allocated(UserPath))deallocate(UserPath)
      UserPath = UserPath_input
   end subroutine set_UserPath


   !---------------------------------------------------------------------------!
   !PURPOSE: generates thew K-points along some pre-stored high-symmetry points.
   !         taken from http://lamp.tu-graz.ac.at/~hadley/ss1/bzones/
   !         and https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#high-symmetry-pathsv
   !---------------------------------------------------------------------------!
   subroutine calc_Kpath(kpt_path,structure,Nkpt_path,Kaxis,KaxisPoints,hetero)
      !
      use utils_misc
      implicit none
      !
      real(8),allocatable,intent(out)       :: kpt_path(:,:)
      character(len=*),intent(in)           :: structure
      integer,intent(in)                    :: Nkpt_path
      logical,intent(in),optional           :: hetero
      real(8),allocatable,intent(out)       :: Kaxis(:)
      real(8),allocatable,intent(out)       :: KaxisPoints(:)
      !
      real(8),dimension(3)                  :: Gamma,M,R,X,K,L,U,W,H,N,P,A,Z,S,T,Y
      real(8),dimension(3)                  :: Kdiff
      real(8),allocatable                   :: Kpoints(:,:),Kdist(:),Kturn(:)
      integer                               :: idir,Ndir,idk,ik,lastK
      real(8)                               :: dKtot,theta,phi,dk,kx,ky,kz
      logical                               :: hetero_
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Kpath"
      !
      !
      hetero_=.false.
      if(present(hetero))hetero_=hetero
      !
      if(allocated(kpt_path))deallocate(kpt_path)
      select case(reg(structure))
         case default
            !
            stop "calc_Kpath: Available structures: triangular, cubic_[2,3], fcc, bcc, hex, tetragonal, orthorhombic_[1,2], User."
            !
         case("triangular")
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            M     = [ 1d0/2d0,     0d0,     0d0 ]
            K     = [ 2d0/3d0, 1d0/3d0,     0d0 ]
            A     = [     0d0,     0d0, 1d0/2d0 ]
            !
            if(hetero_)then
               allocate(Kpoints(3,5));Kpoints=0d0
            else
               allocate(Kpoints(3,4));Kpoints=0d0
            endif
            Kpoints(:,1) = Gamma
            Kpoints(:,2) = M
            Kpoints(:,3) = K
            Kpoints(:,4) = Gamma
            if(hetero_)then
               Kpoints(:,5) = A
               write(*,"(A)") "     structure: simple-triangular (2D+A)."
               write(*,"(A)") "     path: GMKG,A"
            else
               write(*,"(A)") "     structure: simple-triangular (2D)."
               write(*,"(A)") "     path: GMKG"
            endif
            !
         case("cubic_2")
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            X     = [     0d0, 1d0/2d0,     0d0 ]
            M     = [ 1d0/2d0, 1d0/2d0,     0d0 ]
            A     = [     0d0,     0d0, 1d0/2d0 ]
            !
            if(hetero_)then
               allocate(Kpoints(3,5));Kpoints=0d0
            else
               allocate(Kpoints(3,4));Kpoints=0d0
            endif
            Kpoints(:,1) = Gamma
            Kpoints(:,2) = X
            Kpoints(:,3) = M
            Kpoints(:,4) = Gamma
            if(hetero_)then
               Kpoints(:,5) = A
               write(*,"(A)") "     structure: simple-cubic (2D+a)."
               write(*,"(A)") "     path: GXMG,A"
            else
               write(*,"(A)") "     structure: simple-cubic (2D)."
               write(*,"(A)") "     path: GXMG"
            endif
            !
         case("cubic_3")
            !
            if(hetero_) stop "calc_Kpath: [cubic_3] K-path does not allow for heterostructure."
            !
            Gamma = [     0d0,     0d0,     0d0 ]
            M     = [ 1d0/2d0, 1d0/2d0,     0d0 ]
            R     = [ 1d0/2d0, 1d0/2d0, 1d0/2d0 ]
            X     = [     0d0, 1d0/2d0,     0d0 ]
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
            write(*,"(A)") "     structure: simple-cubic (3D)."
            write(*,"(A)") "     path: GXMGRX,MR"
            !
         case("fcc")
            !
            if(hetero_) stop "calc_Kpath: [fcc] K-path does not allow for heterostructure."
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
            if(hetero_) stop "calc_Kpath: [bcc] K-path does not allow for heterostructure."
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
            if(hetero_) stop "calc_Kpath: [hex] K-path does not allow for heterostructure."
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
            if(hetero_) stop "calc_Kpath: [tetragonal] K-path does not allow for heterostructure."
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
            if(hetero_) stop "calc_Kpath: [orthorhombic_1] K-path does not allow for heterostructure."
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
            if(hetero_) stop "calc_Kpath: [orthorhombic_2] K-path does not allow for heterostructure."
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
         case("User")
            !
            if(.not.allocated(UserPath)) stop "calc_Kpath: UserPath not defined."
            Kpoints = UserPath
            write(*,"(A)") "     structure: User provided."
            write(*,"(A)") "     path: "//str(size(Kpoints,dim=2))//" high-symmetry points."
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
      !
      if(allocated(Kaxis))deallocate(Kaxis)
      Kaxis=Kdist
      if(allocated(KaxisPoints))deallocate(KaxisPoints)
      KaxisPoints=Kturn
      !
   end subroutine calc_Kpath


   !---------------------------------------------------------------------------!
   !PURPOSE: generates a uniform K mesh in the kx,ky plane
   !---------------------------------------------------------------------------!
   subroutine calc_Kplane(kpt_plane,Nkpt_Kside)
      !
      use utils_misc
      implicit none
      !
      real(8),allocatable,intent(out)       :: kpt_plane(:,:)
      integer,intent(in)                    :: Nkpt_Kside
      !
      integer                               :: ikx,iky,ik
      real(8)                               :: Kmax=1d0
      real(8)                               :: kx,ky,dK
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Kplane"
      !
      !
      if(allocated(kpt_plane))deallocate(kpt_plane)
      allocate(kpt_plane(3,Nkpt_Kside**2))
      !
      dk=Kmax/(Nkpt_Kside-1)
      !
      ik=0
      do ikx=1,Nkpt_Kside
         do iky=1,Nkpt_Kside
            !
            ik=ik+1
            !
            kx = (ikx-1)*dk - Kmax/2d0
            ky = (iky-1)*dk - Kmax/2d0
            !
            kpt_plane(:,ik) = [kx,ky,0d0]
            !
         enddo
      enddo
      !
   end subroutine calc_Kplane


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate to a user provided K-point path the Hamiltonian
   !---------------------------------------------------------------------------!
   subroutine interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT            &
                                                         ,filename,data_in      &
                                                         ,corrname,correction   &
                                                         ,doplane,Nkpt_Kside    &
                                                         ,hetero)
      !
      use parameters !WHY IS THIS WORKING?
      use utils_misc
      use linalg, only : eigh, inv, zeye
      implicit none
      !
      type(Lattice),intent(inout),target    :: Lttc
      character(len=*),intent(in)           :: structure
      integer,intent(in)                    :: Nkpt_path
      character(len=*),intent(in),optional  :: pathOUTPUT
      character(len=*),intent(in),optional  :: filename
      complex(8),intent(in),target,optional :: data_in(:,:,:)
      character(len=*),intent(in),optional  :: corrname
      complex(8),intent(in),optional        :: correction(:,:,:)
      logical,intent(in),optional           :: doplane
      integer,intent(in),optional           :: Nkpt_Kside
      logical,intent(in),optional           :: hetero
      !
      character(len=256)                    :: path,label,filename_,corrname_
      integer                               :: ik,ikz,iorb,unit
      integer                               :: Norb,Nkpt_Kside_,ikx,iky
      integer                               :: Norb_layer,Nreal,iw,ndx
      real(8)                               :: kp,kx,ky,Bvec(3),wrealMax,eta
      complex(8),pointer                    :: data(:,:,:)
      complex(8),allocatable                :: invGf(:,:)
      logical                               :: Hamiltonian,doplane_,hetero_,printout
      real                                  :: start,finish
      !Hetero
      complex(8),allocatable                :: data_intp_kpkz(:,:,:,:),dataZk_kpkz(:,:,:,:)
      real(8),allocatable                   :: dataEk_kpkz(:,:,:)
      !Interp
      complex(8),allocatable                :: data_intp(:,:,:),dataZk(:,:,:)
      real(8),allocatable                   :: dataEk(:,:)
      !Plots
      real(8),allocatable                   :: Fk(:,:),Akw(:,:,:),wreal(:)
      complex(8),allocatable                :: zeta(:,:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- interpolateHk2Path"
      !
      !
      if(.not.Lttc%status) stop "interpolateHk2Path: Lttc not properly initialized."
      !
      printout=.false.
      if(present(pathOUTPUT))printout=.true.
      !
      filename_="Bands"
      if(present(filename))filename_=reg(filename)
      !
      if(present(data_in))then
         label = reg(filename_)
         data => data_in
         Hamiltonian = .false.
      else
         label = "Hk"
         data => Lttc%Hk
         Hamiltonian = .true.
      endif
      !
      doplane_=.false.
      if(present(doplane))doplane_=doplane
      !
      hetero_=.false.
      if(present(hetero))hetero_=hetero
      !
      !Create path along high-symmetry points-------------------------------
      if(allocated(Lttc%kptpath))deallocate(Lttc%kptpath)
      if(allocated(Lttc%Kpathaxis))deallocate(Lttc%Kpathaxis)
      call calc_Kpath(Lttc%kptpath,reg(structure),Nkpt_path,Lttc%Kpathaxis,Lttc%KpathaxisPoints,hetero=hetero_)
      !
      !path in the bulk
      Lttc%Nkpt_path = size(Lttc%kptpath,dim=2)
      !first Nkpt_path*SymmetryPoints then other Nkpt_path along Gamma-A
      if(hetero_)then
         Norb_layer = int(Lttc%Norb/Lttc%Nsite)
         if((Norb_layer*Lttc%Nsite).ne.Lttc%Norb) stop "interpolateHk2Path: Orbital dimension of Hk is not a multiple of the number of sites."
         Lttc%Nkpt_path = Lttc%Nkpt_path - Nkpt_path
      endif
      !
      Norb = size(data,dim=1)
      call assert_shape(data,[Norb,Norb,Lttc%Nkpt],"interpolateHk2Path",reg(label))
      !
      corrname_="nonInt"
      if(present(correction))then
         call assert_shape(correction,[Norb,Norb,Lttc%Nkpt],"interpolateHk2Path","correction")
         data = data + correction
         corrname_="Corrected"
         if(present(corrname))corrname_=reg(corrname)
         write(*,"(A)")"     Correction: "//reg(corrname_)
      endif
      !
      !Interpolate along path---------------------------------------------------
      allocate(data_intp(Norb,Norb,Lttc%Nkpt_path));data_intp=czero
      call cpu_time(start)
      call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),data,data_intp)
      call cpu_time(finish)
      write(*,"(A,F)") "     "//reg(label)//"(fullBZ) --> "//reg(label)//"(Kpath) cpu timing:", finish-start
      if(Hamiltonian)then
         if(allocated(Lttc%Hk_path))deallocate(Lttc%Hk_path)
         Lttc%Hk_path = data_intp
      endif
      !
      !Compute eigenvalues along path-------------------------------------------
      allocate(dataEk(Norb,Lttc%Nkpt_path));dataEk=0d0
      allocate(dataZk(Norb,Norb,Lttc%Nkpt_path));dataZk=czero
      dataZk = data_intp
      do ik=1,Lttc%Nkpt_path
         call eigh(dataZk(:,:,ik),dataEk(:,ik))
      enddo
      if(Hamiltonian)then
         if(allocated(Lttc%Zk_path))deallocate(Lttc%Zk_path)
         if(allocated(Lttc%Ek_path))deallocate(Lttc%Ek_path)
         Lttc%Zk_path = dataZk
         Lttc%Ek_path = dataEk
         Lttc%pathStored = .true.
      endif
      deallocate(dataZk)
      !
      !Default parameters on the real frequency axis
      wrealMax = 1.2*maxval(abs(dataEk))
      Nreal = 2000
      eta = wrealMax/100
      allocate(wreal(Nreal));wreal=0d0
      wreal = linspace(-wrealMax,+wrealMax,Nreal)
      !
      !Print eigenvalues along path
      if(printout)then
         path = reg(pathOUTPUT)//reg(filename_)//"_"//reg(corrname_)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Lttc%Nkpt_path
            write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik)/Lttc%Kpathaxis(Lttc%Nkpt_path),(dataEk(:,ik),iorb=1,Norb)
         enddo
         close(unit)
      endif
      write(*,"(A,I)") "     Total number of K-points along path:",Lttc%Nkpt_path
      deallocate(dataEk)
      !
      !Compute non-interacting spectral function along path---------------------
      allocate(zeta(Norb,Norb,Nreal));zeta=czero
      do iorb=1,Norb
         do iw=1,Nreal
            zeta(iorb,iorb,iw) = dcmplx(wreal(iw),eta)
         enddo
      enddo
      !
      allocate(Akw(Norb,Nreal,Lttc%Nkpt_path));Akw=0d0
      allocate(invGf(Norb,Norb));invGf=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nreal,wreal,Norb,zeta,Lttc,data_intp,Akw),&
      !$OMP PRIVATE(ik,iw,iorb,invGf)
      !$OMP DO
      do ik=1,Lttc%Nkpt_path
         do iw=1,Nreal
            invGf = zeta(:,:,iw) - data_intp(:,:,ik)
            call inv(invGf)
            do iorb=1,Norb
               Akw(iorb,iw,ik) = dimag(invGf(iorb,iorb))
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      do ik=1,Lttc%Nkpt_path
         do iorb=1,Norb
            Akw(iorb,:,ik) = Akw(iorb,:,ik)/(sum(Akw(iorb,:,ik))*abs(wreal(2)-wreal(1)))
         enddo
      enddo
      deallocate(zeta,invGf)
      !
      !Print k-resolved spectral function
      if(printout)then
         path = reg(pathOUTPUT)//"Akw_"//reg(label)//"_"//reg(corrname_)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Lttc%Nkpt_path
            do iw=1,Nreal
                write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik)/Lttc%Kpathaxis(Lttc%Nkpt_path),wreal(iw),(Akw(iorb,iw,ik),iorb=1,Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
      endif
      deallocate(Akw)
      !
      !Compute dispersion in the Gamma-A direction------------------------------
      if(hetero_)then
         !
         !data for the first Nkpt_path*SymmetryPoints at each kz
         allocate(data_intp_kpkz(Norb_layer,Norb_layer,Lttc%Nkpt_path,0:Nkpt_path));data_intp_kpkz=czero
         call fill_Gamma_A(data_intp,data_intp_kpkz)
         !
         !Compute eigenvalues along path for each kz
         allocate(dataEk_kpkz(Norb_layer,Lttc%Nkpt_path,0:Nkpt_path));dataEk_kpkz=0d0
         allocate(dataZk_kpkz(Norb_layer,Norb_layer,Lttc%Nkpt_path,0:Nkpt_path));dataZk_kpkz=czero
         dataZk_kpkz = data_intp_kpkz
         do ik=1,Lttc%Nkpt_path
            do ikz=0,Nkpt_path
               call eigh(dataZk_kpkz(:,:,ik,ikz),dataEk_kpkz(:,ik,ikz))
            enddo
         enddo
         deallocate(dataZk_kpkz)
         !
         !Print Data along the path at kz=0 plus the Gamma-A direction
         if(printout)then
            path = reg(pathOUTPUT)//reg(filename_)//"_Hetero_"//reg(corrname_)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,Lttc%Nkpt_path
               write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),(dataEk_kpkz(:,ik,0),iorb=1,Norb_layer)
            enddo
            do ikz=1,Nkpt_path
               write(unit,"(1I5,200E20.12)") ik+ikz,Lttc%Kpathaxis(Lttc%Nkpt_path+ikz),(dataEk_kpkz(:,Lttc%iq_gamma,ikz),iorb=1,Norb_layer)
            enddo
            close(unit)
         endif
         write(*,"(A,I)") "     Total number of K-points along path (hetero):",Lttc%Nkpt_path+Nkpt_path
         deallocate(dataEk_kpkz)
         !
         !Compute non-interacting spectral function along path at kz=0 plus the Gamma-A direction
         allocate(zeta(Norb_layer,Norb_layer,Nreal));zeta=czero
         do iorb=1,Norb_layer
            do iw=1,Nreal
               zeta(iorb,iorb,iw) = dcmplx(wreal(iw),eta)
            enddo
         enddo
         !
         allocate(Akw(Norb_layer,Nreal,Lttc%Nkpt_path+Nkpt_path));Akw=0d0
         allocate(invGf(Norb_layer,Norb_layer));invGf=czero
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nreal,wreal,Norb_layer,Nkpt_path,zeta,Lttc,data_intp_kpkz,Akw),&
         !$OMP PRIVATE(ik,iw,iorb,invGf)
         !
         !path in the layer
         !$OMP DO
         do ik=1,Lttc%Nkpt_path
            do iw=1,Nreal
               invGf = zeta(:,:,iw) - data_intp_kpkz(:,:,ik,0)
               call inv(invGf)
               do iorb=1,Norb_layer
                  Akw(iorb,iw,ik) = dimag(invGf(iorb,iorb))
               enddo
            enddo
         enddo
         !$OMP END DO
         !
         !Gamma-A direction
         !$OMP DO
         do ik=1,Nkpt_path
            do iw=1,Nreal
               invGf = zeta(:,:,iw) - data_intp_kpkz(:,:,Lttc%iq_gamma,ik)
               call inv(invGf)
               do iorb=1,Norb_layer
                  Akw(iorb,iw,Lttc%Nkpt_path+ik) = dimag(invGf(iorb,iorb))
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         do ik=1,Lttc%Nkpt_path+Nkpt_path
            do iorb=1,Norb_layer
               Akw(iorb,:,ik) = Akw(iorb,:,ik)/(sum(Akw(iorb,:,ik))*abs(wreal(2)-wreal(1)))
            enddo
         enddo
         deallocate(zeta,invGf,data_intp_kpkz)
         !
         !Print k-resolved spectral function
         if(printout)then
            path = reg(pathOUTPUT)//"Akw_"//reg(label)//"_Hetero_"//reg(corrname_)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,(Lttc%Nkpt_path+Nkpt_path)
               do iw=1,Nreal
                   write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),wreal(iw),(Akw(iorb,iw,ik),iorb=1,Norb_layer)
               enddo
               write(unit,*)
            enddo
            close(unit)
         endif
         deallocate(Akw)
         !
      endif
      !
      deallocate(data_intp,wreal)
      !
      if(doplane_)then
         !
         Nkpt_Kside_ = 201 !fixed
         if(present(Nkpt_Kside)) Nkpt_Kside_ = Nkpt_Kside
         !
         !Create K-points inside the kx,ky plane
         if(allocated(Lttc%kptPlane))deallocate(Lttc%kptPlane)
         call calc_Kplane(Lttc%kptPlane,Nkpt_Kside_)
         Lttc%Nkpt_Plane = size(Lttc%kptPlane,dim=2)
         !
         !Interpolate inside the kx,ky plane-------------------------------------
         allocate(data_intp(Norb,Norb,Lttc%Nkpt_Plane));data_intp=czero
         call cpu_time(start)
         call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,data,data_intp)
         call cpu_time(finish)
         write(*,"(A,F)") "     "//reg(label)//"(fullBZ) --> "//reg(label)//"(kx,ky) cpu timing:", finish-start
         if(Hamiltonian)then
            if(allocated(Lttc%Hk_Plane))deallocate(Lttc%Hk_Plane)
            Lttc%Hk_Plane = data_intp
            Lttc%planeStored = .true.
         endif
         !
         !Compute non-interacting Fermi surface---------------------------------
         allocate(Fk(Norb,Lttc%Nkpt_Plane));Fk=0d0
         allocate(invGf(Norb,Norb));invGf=czero
         do ik=1,Lttc%Nkpt_Plane
            invGf = zeye(Norb)*dcmplx(EcutSheet,eta) - data_intp(:,:,ik)
            call inv(invGf)
            do iorb=1,Norb
               Fk(iorb,ik) = -dimag(invGf(iorb,iorb))
            enddo
         enddo
         deallocate(invGf)
         !
         !Print non-interacting Fermi surface
         if(printout)then
            path = reg(pathOUTPUT)//"Fk_"//reg(label)//"_"//reg(corrname_)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,Lttc%Nkpt_Plane
               ikx = int(ik/(Nkpt_Kside_+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside_-1) - 0.5d0
               iky = ik - (ikx-1)*Nkpt_Kside_      ; ky = (iky-1)/dble(Nkpt_Kside_-1) - 0.5d0
               Bvec = kx*Blat(:,1) + ky*Blat(:,2)
               write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Fk(iorb,ik),iorb=1,Norb)
               if(iky.eq.Nkpt_Kside_)write(unit,*)
            enddo
            close(unit)
         endif
         deallocate(Fk)
         !
         !Compute Fermi surface at given points in the Gamma-A direction--------
         if(hetero_)then
            !
            !data inside the kx,ky plane at each kz
            allocate(data_intp_kpkz(Norb_layer,Norb_layer,Lttc%Nkpt_Plane,0:Nkpt_path));data_intp_kpkz=czero
            call fill_Gamma_A(data_intp,data_intp_kpkz)
            !
            !Compute Fermi surface at kz=0
            allocate(Fk(Norb_layer,Lttc%Nkpt_Plane));Fk=0d0
            allocate(invGf(Norb_layer,Norb_layer));invGf=czero
            do ik=1,Lttc%Nkpt_Plane
               invGf = zeye(Norb_layer)*dcmplx(EcutSheet,eta) - data_intp_kpkz(:,:,ik,0)
               call inv(invGf)
               do iorb=1,Norb_layer
                  Fk(iorb,ik) = -dimag(invGf(iorb,iorb))
               enddo
            enddo
            deallocate(invGf)
            !
            !Print non-interacting Fermi surface at kz=0
            if(printout)then
               path = reg(pathOUTPUT)//"Fk_"//reg(label)//"_Hetero_"//reg(corrname_)//"_kz0.000.DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_Plane
                  ikx = int(ik/(Nkpt_Kside_+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  iky = ik - (ikx-1)*Nkpt_Kside_      ; ky = (iky-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  Bvec = kx*Blat(:,1) + ky*Blat(:,2)
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Fk(iorb,ik),iorb=1,Norb_layer)
                  if(iky.eq.Nkpt_Kside_)write(unit,*)
               enddo
               close(unit)
            endif
            deallocate(Fk)
            !
            !find the kz where to compute the Fermi-surface
            ikz = minloc(abs(Lttc%kptpath(3,1+Lttc%Nkpt_path:Lttc%Nkpt_path+Nkpt_path)-kz_cut),dim=1)
            !
            !Compute Fermi surface at kz=0
            allocate(Fk(Norb_layer,Lttc%Nkpt_Plane));Fk=0d0
            allocate(invGf(Norb_layer,Norb_layer));invGf=czero
            do ik=1,Lttc%Nkpt_Plane
               invGf = zeye(Norb_layer)*dcmplx(EcutSheet,eta) - data_intp_kpkz(:,:,ik,ikz)
               call inv(invGf)
               do iorb=1,Norb_layer
                  Fk(iorb,ik) = -dimag(invGf(iorb,iorb))
               enddo
            enddo
            deallocate(invGf)
            !
            !Print non-interacting Fermi surface at kz=0
            if(printout)then
               path = reg(pathOUTPUT)//"Fk_"//reg(label)//"_Hetero_"//reg(corrname_)//"_kz"//str(kz_cut,3)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_Plane
                  ikx = int(ik/(Nkpt_Kside_+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  iky = ik - (ikx-1)*Nkpt_Kside_      ; ky = (iky-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  Bvec = kx*Blat(:,1) + ky*Blat(:,2)
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Fk(iorb,ik),iorb=1,Norb_layer)
                  if(iky.eq.Nkpt_Kside_)write(unit,*)
               enddo
               close(unit)
            endif
            deallocate(Fk)
            !
         endif
         !
      endif
      if(associated(data))nullify(data)
      !
      !Print position of High-symmetry points in the same folder where the function is
      if(printout)then
         path = reg(pathOUTPUT)//"Kpoints_labels.DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,size(Lttc%KpathaxisPoints,dim=1)
            write(unit,"(1I5,200E20.12)") ik,Lttc%KpathaxisPoints(ik)
         enddo
         close(unit)
         !
         path = reg(pathOUTPUT)//"Kpoints_path.DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,size(Lttc%kptpath,dim=2)
            kp = Lttc%Kpathaxis(ik)
            if(any(abs(Lttc%KpathaxisPoints-kp).lt.eps))then
               ndx = minloc(abs(Lttc%KpathaxisPoints-kp),dim=1)
               write(unit,"(1I5,200E20.12)") ik,kp,Lttc%kptpath(:,ik),Lttc%KpathaxisPoints(ndx)
            else
               write(unit,"(1I5,200E20.12)") ik,kp,Lttc%kptpath(:,ik)
            endif
         enddo
         close(unit)
      endif
      !
      write(*,"(A,I)") "     Total number of High symmetry points:",size(Lttc%KpathaxisPoints,dim=1)
      if(doplane_)write(*,"(A,I)") "     Total number of K-points along {kx,ky} plane:",Lttc%Nkpt_Plane
      !
      !
      !
   contains
      !
      !
      subroutine fill_Gamma_A(data_in,data_out)
         !
         implicit none
         !
         complex(8),intent(in)              :: data_in(:,:,:)
         complex(8),intent(out)             :: data_out(:,:,:,0:)
         !
         integer                            :: ra,rb,ca,cb
         integer                            :: isite,jsite
         integer                            :: Nkpt_layer
         real(8)                            :: kR
         complex(8)                         :: cfac
         !
         Nkpt_layer = size(data_in,dim=3)
         !
         !$OMP PARALLEL DEFAULT(PRIVATE),&
         !$OMP SHARED(Lttc,Nkpt_layer,Nkpt_path,Norb_layer,data_out,data_in)
         !$OMP DO
         do ik=1,Nkpt_layer
            do ikz=0,Nkpt_path
               do isite=1,Lttc%Nsite
                  do jsite=1,Lttc%Nsite
                     !
                     ra = 1+(isite-1)*Norb_layer ; rb = ra + Norb_layer-1
                     ca = 1+(jsite-1)*Norb_layer ; cb = ca + Norb_layer-1
                     !
                     kR = 2*pi * Lttc%kptpath(3,Lttc%Nkpt_path+ikz) * (isite-jsite)
                     cfac = dcmplx(cos(kR),+sin(kR))
                     !
                     data_out(:,:,ik,ikz) = data_out(:,:,ik,ikz) + data_in(ra:rb,ca:cb,ik)*cfac / Lttc%Nsite
                     !
                  enddo
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
      end subroutine fill_Gamma_A
      !
      !
   end subroutine interpolateHk2Path


   !---------------------------------------------------------------------------!
   !PURPOSE: computes the Ewald corrections to the infinite range coulombian
   !---------------------------------------------------------------------------!
   subroutine calc_Ewald(EwaldShift,kpt,eta,mode)
      !
      use parameters !WHY IS THIS WORKING?
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)              :: EwaldShift(:)
      real(8),intent(in)                    :: kpt(:,:)
      real(8),intent(in)                    :: eta
      character(len=*),intent(in)           :: mode
      !
      integer                               :: ik,ir,Nkpt
      real(8)                               :: kR,Lk,funct
      complex(8)                            :: cfac
      !
      !
      if(verbose)write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Ewald"
      if(verbose)write(*,"(A)") "     mode: "//str(reg(mode))
      !
      !
      if(.not.Wig_stored) stop "calc_Ewald: the Wigner-Seiz positions are not stored."
      if((reg(mode).ne."2D").and.(reg(mode).ne."3D")) stop "calc_Ewald: Available modes: 2D, 3D."
      call assert_shape(EwaldShift,[Nwig],"calc_Ewald","EwaldShift")
      EwaldShift=czero
      !
      Nkpt = size(kpt,dim=2)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt,kpt,Nvecwig,EwaldShift,eta,mode),&
      !$OMP PRIVATE(ir,ik,kR,Lk,cfac,funct)
      !$OMP DO
      do ir=1,Nwig
         do ik=1,Nkpt
            !
            kR = 2*pi * dot_product(kpt(:,ik),Nvecwig(:,ir))
            cfac = dcmplx(cos(kR),sin(kR))
            !
            Lk = 2*pi * sqrt(dot_product(kpt(:,ik),kpt(:,ik)))
            if(Lk.lt.eps)cycle
            !
            if(reg(mode).eq."2D") funct = 2d0*pi*erfc(Lk*eta)/Lk
            if(reg(mode).eq."3D") funct = 4d0*pi*exp(-eta*Lk**2)/(Lk**2)
            !
            EwaldShift(ir) = EwaldShift(ir) + funct*cfac/Nkpt
            !
         enddo ! ik
      enddo! ir
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine calc_Ewald


   !---------------------------------------------------------------------------!
   !PURPOSE: computes weights for integrations over the reciprocal space
   !         using the tetrahedron method. In this method "reciprocal cubes"
   !         defined by 8 k points as cube edges are decomposed into six
   !         tetrahedra inside which the electron bands are linearized.
   !
   ! Symmetry considerations:
   ! Unfortunately, the definition of the tetrahedra
   ! is not unique. Any group of six tetrahedra out of the 24 tetrahedra
   ! ( 1,2,3,6, 5,7,3,6, 1,5,3,6, 2,4,3,6, 4,8,3,6, 7,8,3,6,
   !   5,6,2,7, 1,5,2,7, 1,3,2,7, 8,6,2,7, 4,3,2,7, 8,4,2,7,
   !   2,6,1,8, 2,4,1,8, 3,4,1,8, 3,7,1,8, 5,7,1,8, 5,6,1,8,
   !   2,6,4,5, 1,2,4,5, 1,3,4,5, 3,7,4,5, 7,8,4,5, 6,8,4,5 )
   ! that fills up a given cube completely can be chosen. This non-uniqueness
   ! breaks the symmetry of the system. As a result, quantities of symmetry-equivalent
   ! k points might not be identical (but similar). However, even if we sum the
   ! contribution of all 24 tetrahedra and divide by 4, there can still be
   ! symmetry-equivalent k points (e.g. (0,0,1/4) and (1/4,1/4,1/4) in fcc)
   ! where quantities are not identical. This comes from the definition of the
   ! cubes which also break the symmetry.
   ! (For testing, the Gauss method should be used. It is not as accurate but does
   ! not suffer from these symmetry-related problems.)
   !
   ! --------------------------------------------------------------------------!
   !
   !     Initializes the tetrahedron integration on the Fermi surface,
   !     i.e. calculates (for any f) the weights wintgr4 in the sum
   !
   !                                   occ                     2
   !     SUM  wintgr4(n,k) * f(n,k) = SUM     INT     f(n,k) d k
   !      k,n                           n   e (k) = E
   !                                         n       F
   !---------------------------------------------------------------------------!
   subroutine tetrahedron_integration(pathINPUT,Ek_orig,nkpt3,kpt,Egrid,weights_out,DoS_out,fact_intp,pathOUTPUT,store_weights)
      !
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      real(8),intent(in)                    :: Ek_orig(:,:)
      integer,intent(in)                    :: nkpt3(3)
      real(8),intent(in)                    :: kpt(:,:)
      real(8),intent(in)                    :: Egrid(:)
      integer,intent(in),optional           :: fact_intp
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: store_weights
      real(8),intent(out),optional          :: weights_out(:,:,:)
      real(8),intent(out),optional          :: DoS_out(:)
      !
      integer                               :: igrid,Ngrid
      integer                               :: iorb,Norb
      integer                               :: Nkpt,Nkpti,ik,ik0
      integer,allocatable                   :: kptp(:),pkpt(:,:,:)
      real(8),allocatable                   :: Ek(:,:),Ek_intp(:,:)
      real(8),allocatable                   :: nkstar(:),DoS(:,:)
      real(8),allocatable                   :: kpt_intp(:,:)
      real(8),allocatable                   :: weights(:,:,:)
      integer                               :: k1,k2,k3,i,kk
      integer                               :: itria,ntria,itetra,unit
      integer                               :: pnt(4),kindx(8),Nkpt3_used(3)
      real(8)                               :: rdum,dE,DoSnorm
      real(8)                               :: Ecube(8),Etetra(4)
      real(8)                               :: f(4),atria(2),wtria(4,2)
      logical                               :: filexists,store_weights_
      integer,parameter                     :: tetra(4,6)=reshape( (/ 1,2,3,6, 5,3,6,7, 1,5,3,6, 8,6,7,2, 4,7,2,3, 8,4,7,2 /),(/ 4,6 /) )
      !
      !
      if(verbose)write(*,"(A)") "---- tetrahedron_integration"
      !
      !
      call inquireFile(reg(pathINPUT)//"sym.DAT",filexists,hardstop=.false.)
      if(.not.filexists)then
         write(*,"(A)") "     Tetrahedron integration skipped."
         return
      endif
      !
      Norb = size(Ek_orig,dim=1)
      Ngrid = size(Egrid)
      dE = abs(Egrid(10)-Egrid(9))
      !
      store_weights_=.false.
      if(present(store_weights))store_weights_=store_weights
      !
      if(present(fact_intp))then
         !
         Nkpt3_used = Nkpt3*fact_intp
         Nkpt = product(Nkpt3_used)
         !
         !Generate K-points in the irreducible BZ (finer K-mesh)
         call calc_irredBZ(reg(pathINPUT),Nkpt3_used,Nkpti,kptp,pkpt,nkstar,kpt_out=kpt_intp)
         !
         !Interpolate Ek to the new K mesh
         allocate(Ek_intp(Norb,Nkpt));Ek_intp=0d0
         call wannierinterpolation(Nkpt3,kpt,kpt_intp,Ek_orig,Ek_intp)
         !
         !The irreducible K-points are always the first Nkpti found by kptp
         allocate(Ek(Norb,Nkpt));Ek=0d0
         do ik=1,Nkpt
            do iorb=1,Norb
               Ek(iorb,ik) = Ek_intp(iorb,kptp(ik))
            enddo
         enddo
         deallocate(Ek_intp)
         !
      else
         !
         !
         Nkpt3_used = Nkpt3
         Nkpt = product(Nkpt3_used)
         !
         !Generate K-points in the irreducible BZ
         call calc_irredBZ(reg(pathINPUT),Nkpt3_used,Nkpti,kptp,pkpt,nkstar,store=store_weights_)
         !
         !allocate dispersion
         allocate(Ek(Norb,Nkpt));Ek=0d0
         do ik=1,Nkpt
            do iorb=1,Norb
               Ek(iorb,ik) = Ek_orig(iorb,kptp(ik))
            enddo
         enddo
         !
      endif
      !
      allocate(weights(Ngrid,Norb,Nkpt));weights=0d0
      !$OMP PARALLEL DEFAULT(PRIVATE),&
      !$OMP SHARED(weights,Ek,Egrid,Nkpt3_used,Nkpti,kptp,pkpt,nkstar,Ngrid,Norb,tetra)
      !$OMP DO
      do igrid=1,Ngrid
         do iorb=1,Norb
            !
            !boundaries of the energy grid
            if(minval(Ek(iorb,:)).ge.Egrid(igrid)) cycle
            if(maxval(Ek(iorb,:)).le.Egrid(igrid)) cycle
            !
            !Loop over cubes in BZ
            do k3 = 1,Nkpt3_used(3)
               do k2 = 1,Nkpt3_used(2)
                  do k1 = 1,Nkpt3_used(1)
                     !
                     !Definition of the cube corners
                     kindx(1) = pkpt(k1  ,k2  ,k3  )
                     kindx(2) = pkpt(k1+1,k2  ,k3  )
                     kindx(3) = pkpt(k1  ,k2+1,k3  )
                     kindx(4) = pkpt(k1+1,k2+1,k3  )
                     kindx(5) = pkpt(k1  ,k2  ,k3+1)
                     kindx(6) = pkpt(k1+1,k2  ,k3+1)
                     kindx(7) = pkpt(k1  ,k2+1,k3+1)
                     kindx(8) = pkpt(k1+1,k2+1,k3+1)
                     !
                     !Get energies at cube corners
                     Ecube = Ek(iorb,kindx)
                     !
                     do itetra=1,6
                        !
                        !Get energies at tetrahedron corners
                        Etetra = Ecube(tetra(:,itetra))
                        !
                        !Nothing to be done if tetrahedron-corner energies larger/smaller than Egrid(igrid)
                        if(all(Etetra.ge.Egrid(igrid))) cycle
                        if(all(Etetra.le.Egrid(igrid))) cycle
                        !
                        !size-ordered Etetra
                        call rorderp(pnt,Etetra,4)
                        Etetra = Etetra(pnt)
                        !
                        !Intersecting plane is a triangle (ntria=1) or a quadrangle (ntria=2).
                        !tria contains the triangle points and atria the triangle area(s).
                        if(Egrid(igrid).gt.Etetra(3))then
                           ntria       = 1
                           f(1:3)      = (Egrid(igrid)-Etetra(1:3))/(Etetra(4)-Etetra(1:3))
                           atria(1)    = 3 * (Egrid(igrid)-Etetra(4))**2 / ( (Etetra(4)-Etetra(1)) * (Etetra(4)-Etetra(2)) * (Etetra(4)-Etetra(3)) )
                           wtria(:3,1) = 1-f(:3)
                           wtria(4,1)  = f(1)+f(2)+f(3)
                        elseif(Egrid(igrid).gt.Etetra(2))then
                           ntria       = 2
                           f(1:2)      = (Egrid(igrid)-Etetra(1))/(Etetra(3:4)-Etetra(1))
                           f(3:4)      = (Egrid(igrid)-Etetra(2))/(Etetra(3:4)-Etetra(2))
                           rdum        = (Etetra(4)-Etetra(2)) * (Etetra(3)-Etetra(1))
                           atria(1)    = 3*f(3) * (Etetra(3)-Egrid(igrid)) / rdum
                           atria(2)    = 3*f(2) * (Etetra(4)-Egrid(igrid)) / rdum
                           wtria(1,1)  = 1-f(1)      ; wtria(1,2) = 2-f(1)-f(2)
                           wtria(2,1)  = 2-f(3)-f(4) ; wtria(2,2) = 1-f(4)
                           wtria(3,1)  = f(1)+f(3)   ; wtria(3,2) = f(1)
                           wtria(4,1)  = f(4)        ; wtria(4,2) = f(2)+f(4)
                        else ! if(h0.gt.htetra(1)) then
                           ntria       = 1
                           f(1:3)      = (Egrid(igrid)-Etetra(1))/(Etetra(2:4)-Etetra(1))
                           atria(1)    = 3*f(1)*f(2) / (Etetra(4)-Etetra(1))
                           wtria(1,1)  = 3-f(1)-f(2)-f(3)
                           wtria(2:,1) = f(:3)
                        endif
                        wtria(pnt,:) = wtria
                        !
                        do itria=1,ntria
                           do i=1,4
                              kk = kindx(tetra(i,itetra))
                              weights(igrid,iorb,kk) = weights(igrid,iorb,kk) + atria(itria) * wtria(i,itria)
                           enddo
                        enddo
                        !
                     enddo ! itetra
                     !
                  enddo ! k1
               enddo ! k2
            enddo ! k3
         enddo ! iorb
      enddo ! igrid
      !$OMP END DO
      !$OMP END PARALLEL
      !
      !take average over three triangle edges, six tetrahedra per cube, nkpt cubes
      do igrid=1,Ngrid
         !
         weights(igrid,:,:) = weights(igrid,:,:)/3d0/6d0/dble(nkpt)
         !
         !Average over symmetry-equivalent k-points
         do iorb=1,Norb
            do ik=nkpti+1,nkpt
               ik0=kptp(ik)
               weights(igrid,iorb,ik0)=weights(igrid,iorb,ik0)+weights(igrid,iorb,ik)
            enddo
            do ik0=1,nkpti
               weights(igrid,iorb,ik0)=weights(igrid,iorb,ik0)/nkstar(ik0)
            enddo
            do ik=nkpti+1,nkpt
               ik0=kptp(ik)
               weights(igrid,iorb,ik)=weights(igrid,iorb,ik0)
            enddo
         enddo
         !
      enddo
      deallocate(kptp,pkpt,nkstar,Ek)
      !
      if(present(pathOUTPUT).and.store_weights_)then
         unit = free_unit()
         open(unit,file=reg(pathOUTPUT)//"Weights_BZ.DAT",form="formatted",status="unknown",position="rewind",action="write")
         do igrid=1,Ngrid
            do ik=1,Nkpt
               write(unit,"(1F20.10,1I8,200F20.10)")Egrid(igrid),ik,(weights(igrid,iorb,ik),iorb=1,Norb)
            enddo
         enddo
         close(unit)
      endif
      !
      !Compute the DoS - the grid might not be uniform
      allocate(DoS(Ngrid,Norb));DoS=0d0
      do iorb=1,Norb
         DoSnorm=0d0
         do igrid=1,Ngrid
            DoS(igrid,iorb) = sum(weights(igrid,iorb,:))
            if(igrid.ge.2)then
               dE = abs(Egrid(igrid)-Egrid(igrid-1))
               DoSnorm = DoSnorm + ( DoS(igrid-1,iorb)+DoS(igrid,iorb) ) * (dE/2d0)
            endif
         enddo
         write(*,"(A,F)") "     Normalization DoS_"//str(iorb)//":", DoSnorm
      enddo
      if(present(pathOUTPUT))then
         unit = free_unit()
         open(unit,file=reg(pathOUTPUT)//"DoS.DAT",form="formatted",status="unknown",position="rewind",action="write")
         do igrid=1,Ngrid
            write(unit,"(200F20.10)")Egrid(igrid),(DoS(igrid,iorb),iorb=1,Norb)
         enddo
         close(unit)
      endif
      !
      if(present(weights_out))then
         call assert_shape(weights_out,[Ngrid,Norb,Nkpt],"tetrahedron_integration","weights_out")
         weights_out = weights
      endif
      deallocate(weights)
      !
      if(present(DoS_out))then
         call assert_shape(DoS_out,[Ngrid],"tetrahedron_integration","DoS_out")
         DoS_out = sum(DoS,dim=2)
      endif
      deallocate(DoS)
      !
      !
   contains
      !
      !
      !
      !Orders rarr(1:n) according to size and returns a correspondingly defined pointer in pnt.
      subroutine rorderp(pnt,rarr,n)
          implicit none
          integer,intent(in)                :: n
          integer,intent(out)               :: pnt(n)
          real(8),intent(in)                :: rarr(n)
          integer                           :: i,j,k
          !
          do i=1,n
             pnt(i) = i
             do j=1,i-1
                if(rarr(pnt(j)).gt.rarr(i)) then
                   do k=i,j+1,-1
                      pnt(k) = pnt(k-1)
                   enddo
                   pnt(j) = i
                   exit
                endif
             enddo
          enddo
          !
      end subroutine rorderp
      !
      !
      !
   end subroutine tetrahedron_integration

end module crystal
