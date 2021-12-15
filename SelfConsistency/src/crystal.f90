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
   real(8),allocatable,private,protected    :: Ruc(:,:)
   real(8),private,protected                :: Rlat(3,3)
   real(8),private,protected                :: Blat(3,3)
   real(8),private,protected                :: vol
   !
   !Available outside the module
   !
   integer,public,protected                 :: Nwig=0
   integer,public,protected                 :: wig0=0
   real(8),allocatable,public,protected     :: radiuswig(:)
   integer,allocatable,public,protected     :: Nvecwig(:,:)
   real(8),allocatable,public,protected     :: Rvecwig(:,:)
   integer,allocatable,public,protected     :: nrdegwig(:)
   !
   real(8),allocatable,public,protected     :: KvecPolar(:,:)
   integer,allocatable,public,protected     :: small_ik(:,:)
   !
   real(8),allocatable,public,protected     :: UserPath(:,:)
   !
   logical,public,protected                 :: Ruc_stored=.false.               !Global flag for routines that need positions within the u.c.
   logical,public,protected                 :: Lat_stored=.false.               !Global flag for routines that need rlat
   logical,public,protected                 :: Wig_stored=.false.               !Global flag for routines performing Wannier interpolation
   logical,public,protected                 :: small_ik_stored=.false.          !Global flag for routines that need sorted k-points
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
   public :: set_lattice
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
   public :: get_Ruc,get_Rlat,get_Blat
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
   subroutine set_lattice(Rlat_input,Ruc_input)
      !
      use utils_misc
      use linalg, only : det, inv_sym, cross_product
      implicit none
      !
      real(8),intent(in)                    :: Rlat_input(3,3)
      real(8),intent(in)                    :: Ruc_input(:,:)
      integer                               :: ir
      !
      !
      if(verbose)write(*,"(A)") "---- set_lattice"
      !
      !
      Rlat(:,1) = Rlat_input(:,1)
      Rlat(:,2) = Rlat_input(:,2)
      Rlat(:,3) = Rlat_input(:,3)
      vol = dot_product(cross_product(Rlat(:,1),Rlat(:,2)),Rlat(:,3))
      if(verbose)write(*,"(A,F)")"     Unit cell volume: ",vol
      Blat(:,1) = cross_product(Rlat(:,2),Rlat(:,3))/vol
      Blat(:,2) = cross_product(Rlat(:,3),Rlat(:,1))/vol
      Blat(:,3) = cross_product(Rlat(:,1),Rlat(:,2))/vol
      !
      if(allocated(Ruc))deallocate(Ruc)
      allocate(Ruc(3,size(Ruc_input,dim=2)));Ruc=0d0
      Ruc = Ruc_input
      !
      write(*,"(A)")new_line("A")//"     Unit cell vectors: "
      do ir=1,3
         write(*,"(A)")"     R_"//str(ir)//": [ "//str(Rlat(1,ir),3)//" , "//str(Rlat(2,ir),3)//" , "//str(Rlat(3,ir),3)//" ]"
      enddo
      write(*,"(A)")new_line("A")//"     Positions inside the unit cell: "
      do ir=1,size(Ruc,dim=2)
         write(*,"(A)")"     r_"//str(ir)//": [ "//str(Ruc(1,ir),3)//" , "//str(Ruc(2,ir),3)//" , "//str(Ruc(3,ir),3)//" ]"
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
   subroutine get_Ruc(Ruc_out)
      implicit none
      real(8),allocatable,intent(inout)     :: Ruc_out(:,:)
      if(Lat_stored)then
         if(allocated(Ruc_out))deallocate(Ruc_out)
         allocate(Ruc_out(size(Ruc,dim=1),size(Ruc,dim=2)));Ruc_out=0d0
         Ruc_out = Ruc
      else
         write(*,"(A)")"     Warning: requested unit cell vetors but lattice is not stored."
      endif
   end subroutine get_Ruc
   !
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
      ! create the kpt List
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
            write(*,"(A,1I4,2(A,3I3))")"     ik= ",ik," kpt(:,ik)= ",kpt(:,ik)," kpt_xeps(:,ik)= ",kpt_xeps(:,ik)
            !write(*,"(A)") "kptp(ik)=",kptPos(ik),"kptp_loc(ik)=",kptPos_xeps(ik)
            stop "read_xeps: K-points grid does not match."
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
   subroutine read_Hk(Hk,kpt,pathINPUT,filename)
      !
      use utils_misc
      implicit none
      !
      complex(8),allocatable,intent(out)    :: Hk(:,:,:)
      real(8),allocatable,intent(out)       :: kpt(:,:)
      character(len=*),intent(in)           :: pathINPUT
      character(len=*),intent(in),optional  :: filename
      !
      character(len=256)                    :: path
      integer                               :: unit,Nkpt,Norb
      integer                               :: iorb,jorb,ik
      integer                               :: idum1,idum2
      real(8)                               :: ReHk,ImHk
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_Hk"
      !
      !
      if(.not.Lat_stored) stop "read_Hk: Lattice vectors not stored."
      !
      ! Look for Hk.DAT
      path=reg(pathINPUT)//"Hk.DAT"
      if(present(filename))path=reg(pathINPUT)//reg(filename)
      write(*,"(A)")"     Reading "//reg(path)
      call inquireFile(reg(path),filexists)
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) idum1,Nkpt,Norb
      !
      write(*,"(A)") "     Reading "//reg(path)
      write(*,"(A)") "     Number of K-points from file: "//str(Nkpt)
      write(*,"(A)") "     Orbital space from file: "//str(Norb)
      !
      if(allocated(Hk))deallocate(Hk)
      if(allocated(kpt))deallocate(kpt)
      !
      allocate(Hk(Norb,Norb,Nkpt));Hk=czero
      allocate(kpt(3,Nkpt));kpt=0d0
      !
      Hk=czero
      do ik=1,Nkpt
         read(unit,*) idum1,idum2,kpt(:,ik)
         if (idum2.ne.ik) stop "read_Hk: wrong index ik"
         do iorb=1,Norb
            do jorb=1,Norb
               !
               read(unit,*) idum1,idum2,ReHk,ImHk
               if (idum1.ne.iorb) stop "read_Hk: wrong index iwan1."
               if (idum2.ne.jorb) stop "read_Hk: wrong index iwan2."
               !
               Hk(iorb,jorb,ik) = dcmplx(ReHk,ImHk)
               !
            enddo
         enddo
         !
         call check_Hermiticity(Hk(:,:,ik),eps)
         !
      enddo
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
   subroutine build_Hk(Hk,kpt,hopping,Nkpt3,readHr,Hetero,pathOUTPUT)
      !
      use utils_misc
      use parameters, only : Heterostructures !WHY IS THIS WORKING?
      use linalg, only : zeye, diagonal, diag, eigh, dag
      implicit none
      !
      complex(8),allocatable,intent(out)    :: Hk(:,:,:)
      real(8),allocatable,intent(out)       :: kpt(:,:)
      real(8),intent(in)                    :: hopping(:)
      integer,intent(in)                    :: Nkpt3(3)
      logical,intent(in)                    :: readHr
      type(Heterostructures),intent(inout)  :: Hetero
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      integer                               :: Nkpt,Norb
      integer                               :: Nsite,Nsite_bulk
      integer                               :: iorb,jorb,io,jo,ik
      integer                               :: Trange,iwig,iD
      integer                               :: unit
      complex(8),allocatable                :: Hr(:,:,:),Hr_bulk(:,:,:)
      !W90
      integer,parameter                     :: W90NumCol=15
      integer                               :: Num_wann,Nrpts
      integer                               :: Qst,Rst,i,j,ir
      integer                               :: nx,ny,nz
      integer,allocatable                   :: Ndegen(:)
      real(8)                               :: ReHr,ImHr
      character(len=256)                    :: path
      logical                               :: filexists,rebuild
      !multi-site
      integer                               :: isite,jsite
      integer                               :: ilayer,jlayer,islab
      integer                               :: na,nb,site_i,site_j,io_l,jo_l
      real(8)                               :: Rvec(3),Rdist!,fact
      real(8),allocatable                   :: Rsorted(:,:),Rsorted_bkp(:,:)
      integer,allocatable                   :: Rorder(:),Dist(:,:),DistList(:)
      logical,allocatable                   :: layerDone(:)
      !
      !
      if(verbose)write(*,"(A)") "---- build_Hk"
      !
      !
      path=reg(pathOUTPUT)//"Hk_built.DAT"
      call inquireFile(reg(path),filexists,hardstop=.false.)
      !
      Nkpt = product(Nkpt3)
      Norb = size(hopping)
      Nsite = size(Ruc,dim=2)
      Nsite_bulk = Nsite
      if(Hetero%status) Nsite_bulk = int(Nsite/Hetero%Nlayer)
      !
      Trange=1
      !
      if(filexists)then
         !
         rebuild=.false.
         write(*,"(A)")"     Reading Hk_built.DAT from "//reg(pathOUTPUT)
         call read_Hk(Hk,kpt,reg(pathOUTPUT),filename="Hk_built.DAT")
         if(size(Hk,dim=3).ne.Nkpt)then
            write(*,"(A)")"     Hk_built.DAT has the wrong number of K-points. Rebuiuding."
            rebuild=.true.
         elseif(size(Hk,dim=1).ne.(Nsite*Norb))then
             write(*,"(A)")"     Hk_built.DAT has the wrong orbital dimension. Rebuiuding."
             rebuild=.true.
          endif
         !
      endif
      !
      if((.not.filexists).or.rebuild)then
         !
         write(*,"(A)")"     Building Hk.DAT from input parameters."
         if(.not.Lat_stored) stop "build_Hk: Lattice vectors not stored."
         if(readHr.and.(.not.present(pathOUTPUT))) stop "build_Hk: reading of Hr.DAT requested but missing path."
         !
         if(allocated(Hk))deallocate(Hk)
         allocate(Hk(Norb*Nsite,Norb*Nsite,Nkpt));Hk=czero
         !
         if(allocated(kpt))deallocate(kpt)
         allocate(kpt(3,Nkpt));kpt=0d0
         !
         call build_kpt(Nkpt3,kpt,pathOUTPUT=reg(pathOUTPUT))
         !
         !recover the vectors in real space
         if(.not.Wig_stored)call calc_wignerseiz(Nkpt3)
         !
         if(readHr)then
            !
            ! Look for Hr.DAT
            path=reg(pathOUTPUT)//"Hr.DAT"
            call inquireFile(reg(path),filexists)
            !
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            read(unit,*)                      !skip first line
            read(unit,*) Num_wann !Number of Wannier orbitals
            read(unit,*) Nrpts    !Number of Wigner-Seitz vectors
            !
            if(Norb*Nsite_bulk.ne.Num_wann) stop "build_Hk: Hr.DAT does not have the correct dimension."
            allocate(Hr_bulk(Norb*Nsite_bulk,Norb*Nsite_bulk,Nwig));Hr_bulk=czero
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
                     read(unit,*) nx, ny, nz, iorb, jorb, ReHr, ImHr
                     !
                     iwig = find_vec([nx,ny,nz],Nvecwig)
                     !
                     Hr_bulk(iorb,jorb,iwig) = dcmplx(ReHr,ImHr)/Ndegen(ir)
                     !
                  enddo
               enddo
            enddo
            close(unit)
            deallocate(Ndegen)
            !
         else
            !
            !Get all the possible positions
            allocate(Rsorted(Nwig*Nsite_bulk*Nsite_bulk,4));Rsorted=0d0
            iR=0
            do iwig=1,Nwig
               do jsite=1,Nsite_bulk
                  do isite=1,Nsite_bulk
                     !
                     Rvec = Rvecwig(:,iwig) + Ruc(:,jsite) - Ruc(:,isite)
                     Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                     !
                     iR = iR +1
                     !
                     Rsorted(iR,1) = Rdist
                     Rsorted(iR,2) = iwig
                     Rsorted(iR,3) = jsite
                     Rsorted(iR,4) = isite
                     !
                  enddo
               enddo
            enddo
            !
            !Sorting the positions according to distance
            allocate(Rorder(Nwig*Nsite_bulk*Nsite_bulk));Rorder=0
            allocate(Rsorted_bkp(Nwig*Nsite_bulk*Nsite_bulk,4));Rsorted_bkp=0d0
            Rsorted_bkp = Rsorted
            call sort_array(Rsorted(:,1),Rorder)
            Rsorted=0d0
            do iR=1,Nwig*Nsite_bulk*Nsite_bulk
               Rsorted(iR,:) = Rsorted_bkp(Rorder(iR),:)
            enddo
            deallocate(Rsorted_bkp,Rorder)
            !
            !Regroup according to distance. The list contains the indexes of all the positions with a given distance
            call get_pattern(Dist,Rsorted(:,1),1e4*eps,listDim=DistList,IncludeSingle=.true.)
            !
            !User-provided hopping is only nearest neighbor and its the same for all the sites
            allocate(Hr_bulk(Norb*Nsite_bulk,Norb*Nsite_bulk,Nwig));Hr_bulk=czero
            !all the possible ranges, iD=1: local energy iD=2 nearest neighbor hopping
            do iD=1,Trange+1
               !all the indexes within that range
               do iR=1,DistList(iD)
                  !
                  !retrieve indexes from sorted list
                  iwig = Rsorted(Dist(iD,iR),2)
                  jsite = Rsorted(Dist(iD,iR),3)
                  isite = Rsorted(Dist(iD,iR),4)
                  !
                  if(iD.eq.1)then
                     !
                     !Local energy
                     if(Rsorted(Dist(iD,iR),2).ne.wig0) stop "build_Hk: wrong R=0 local index in position list."
                     if(Rsorted(Dist(iD,iR),1).ne.0d0) stop "build_Hk: wrong R=0 distance in position list."
                     !
                  else
                     !
                     !nearest neighbor hopping (in the plane)
                     do iorb=1,Norb
                        !
                        !site-orbital indexes of Hr
                        io = iorb + Norb*(isite-1)
                        jo = iorb + Norb*(jsite-1)
                        !
                        Hr_bulk(io,jo,iwig) = dcmplx(hopping(iorb),0d0)
                        !
                     enddo
                     !
                  endif
                  !
               enddo
            enddo
            !
         endif
         !
         if(verbose.and.(.not.readHr))then
            path="Hr_report.DAT"
            if(present(pathOUTPUT))then
               path=reg(pathOUTPUT)//"Hr_report.DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            else
               unit=6
            endif
            do iD=1,size(Dist,dim=1)
               write(unit,"(A)") "     Dist: "//str(iD)
               write(unit,"(A8,6A6,20A12)") "ndx" , "n1" , "n2" , "n3" , "iwig" , "jsite" , "isite" , "R" , "H(Ri)"
               do iR=1,DistList(iD)
                  iwig = Rsorted(Dist(iD,iR),2)
                  jsite = Rsorted(Dist(iD,iR),3)
                  isite = Rsorted(Dist(iD,iR),4)
                  Rvec = Rvecwig(:,iwig) + Ruc(:,jsite) - Ruc(:,isite)
                  Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                  write(unit,"(I8,6I6,20F12.4)") iR,Nvecwig(:,iwig),iwig,jsite,isite,Rsorted(Dist(iD,iR),1),&
                  (dreal(Hr_bulk(iorb + Norb*(isite-1),iorb + Norb*(jsite-1),iwig)),iorb=1,Norb)
                  if(Rsorted(Dist(iD,iR),1).ne.Rdist)then
                     write(unit,"(A,2F12.4)") "ERROR: Rsorted(Dist(iD,iR),1).ne.Rdist",Rsorted(Dist(iD,iR),1),Rdist
                     stop "build_Hk:  check Hr_report.DAT"
                  endif
               enddo
            enddo
            deallocate(Rsorted,Dist,DistList)
            if(present(pathOUTPUT))close(unit)
         endif
         !
         !Set up the heterostructure
         if(Hetero%status)then
            !
            allocate(Hr(Norb*Nsite,Norb*Nsite,Nwig));Hr=czero
            !
            !reshuffle the single-layer Hr into the heterostructured one
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  do jsite=1,Nsite_bulk
                     do iorb=1,Norb
                        !
                        io = iorb + Norb*(isite-1)
                        jo = iorb + Norb*(jsite-1)
                        !
                        io_l = io + Norb*Nsite_bulk*(ilayer-1)
                        jo_l = jo + Norb*Nsite_bulk*(ilayer-1)
                        !
                        Hr(io_l,jo_l,:) = Hr_bulk(io,jo,:)
                        !
                     enddo
                  enddo
               enddo
            enddo
            !
            !adding the hopping between the layers. This implies that the Ruc are ordered correctly!
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  !
                  allocate(Rsorted(Nsite*Nwig,4));Rsorted=0d0
                  allocate(Rsorted_bkp(Nsite*Nwig,4));Rsorted_bkp=0d0
                  allocate(Rorder(Nsite*Nwig));Rorder=0
                  !
                  site_i = isite + (ilayer-1)*Nsite_bulk
                  !
                  iR=0
                  do jlayer=1,Hetero%Nlayer
                     do jsite=1,Nsite_bulk
                        !
                        site_j = jsite + (jlayer-1)*Nsite_bulk
                        !
                        do iwig=1,Nwig
                           !
                           Rvec = Rvecwig(:,iwig) + Ruc(:,site_j) - Ruc(:,site_i)
                           Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                           !
                           iR = iR +1
                           !
                           Rsorted(iR,1) = Rdist
                           if(ilayer.eq.jlayer) Rsorted(iR,1)=1e6               !push away the hopping on the same layer
                           Rsorted(iR,2) = iwig
                           Rsorted(iR,3) = jsite
                           Rsorted(iR,4) = jlayer
                           !
                        enddo
                     enddo
                  enddo
                  !
                  !Re-ordering distances for each site in the system
                  Rorder=0
                  Rsorted_bkp = Rsorted
                  call sort_array(Rsorted(:,1),Rorder)
                  do iR=1,Nwig*Nsite
                     Rsorted(iR,:) = Rsorted_bkp(Rorder(iR),:)
                  enddo
                  !
                  !Regrouping according to distance. The list contains the indexes of all the positions with a given distance
                  call get_pattern(Dist,Rsorted(:,1),1e4*eps,listDim=DistList,IncludeSingle=.true.)
                  !
                  !add the inter-layer hopping up to the first neighbors
                  !QUI ID=1 È IL PRIMO VICINO MENTRE ID=2 È IL SECONDO VICINO DATO TUTTO TRA LAYER DIVERSI DATO CHE HO TOLTO
                  !L INDICE DENTRO LO STESSO LAYER ORA DEVO AGGIUNGERE UN QUALCHE SCALING A SECONDA DELLA DISTANZA
                  !
                  !QUI LO DEVO TOGLIERE PERCHE VOGLIO I PIRMI DUE TRA LAYER DIVERSI
                  !MENTRE PER L INTERAZIONE NON POSSO TENERE TUTTE LE DISTANZE TANTO IL CUTOFF È DATO DAL VRANGE CHE
                  !MI DICE SE HO, O MENO, INTERAZIONE TRA I LAYER.
                  !
                  allocate(layerDone(Hetero%Nlayer));layerDone=.false.
                  do iD=1,Hetero%tzRange
                     !
                     !Rescaling factor proportional to the relative distance
                     !fact = Rsorted(Dist(1,1),1) / Rsorted(Dist(iD,1),1)
                     !write(*,*)"fact",fact,iD
                     !
                     !all the indexes within that range
                     do iR=1,DistList(iD)
                        !
                        !retrieve indexes from sorted list
                        iwig = Rsorted(Dist(iD,iR),2)
                        jsite = Rsorted(Dist(iD,iR),3)
                        jlayer = Rsorted(Dist(iD,iR),4)
                        !
                        do iorb=1,Norb
                           !
                           !site-orbital indexes of Hr
                           io = iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1)
                           jo = iorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1)
                           !
                           islab = min(ilayer,jlayer) + Hetero%Explicit(1) - 1
                           !PROJECT SPECIFIC (TaS2) >>>
                           !if(abs(ilayer-jlayer).eq.1) Hr(io,jo,iwig) = dcmplx(Hetero%tz(iorb,islab),0d0)*fact
                           if((abs(ilayer-jlayer).eq.1).and.(.not.layerDone(jlayer))) Hr(io,jo,iwig) = dcmplx(Hetero%tz(iorb,islab),0d0)
                           !>>> PROJECT SPECIFIC (TaS2)
                           !
                        enddo
                        !
                     enddo
                     layerDone(jlayer)=.true.
                  enddo
                  deallocate(layerDone)
                  !
                  if(verbose)then
                     path="Hr_report_Hetero_"//str(site_i)//".DAT"
                     if(present(pathOUTPUT))then
                        path=reg(pathOUTPUT)//"Hr_report_Hetero_"//str(site_i)//".DAT"
                        unit = free_unit()
                        open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
                     else
                        unit=6
                     endif
                     write(unit,"(A)") "     Global site: "//str(site_i)
                     do iD=1,size(Dist,dim=1)
                        write(unit,"(A)") "     Dist: "//str(iD)
                        write(unit,"(A8,8A6,20A12)") "ndx" , "n1" , "n2" , "n3" , "iwig" , "ilay" , "isite" , "jlay" , "jsite" , "R", "H(Ri)"
                        do iR=1,DistList(iD)
                           iwig = Rsorted(Dist(iD,iR),2)
                           jsite = Rsorted(Dist(iD,iR),3)
                           jlayer = Rsorted(Dist(iD,iR),4)
                           Rvec = Rvecwig(:,iwig) + Ruc(:,jsite+Nsite_bulk*(jlayer-1)) - Ruc(:,isite+Nsite_bulk*(ilayer-1))
                           Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                           write(unit,"(I8,8I6,20F12.4)") iR,Nvecwig(:,iwig),iwig,ilayer,isite,jlayer,jsite,Rsorted(Dist(iD,iR),1),&
                           (dreal(Hr(iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1),iorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1),iwig)),iorb=1,Norb)
                           if(Rsorted(Dist(iD,iR),1).ne.Rdist)then
                              write(unit,"(A,2F20.4)") "Warning: Rsorted(Dist(iD,iR),1).ne.Rdist",Rsorted(Dist(iD,iR),1),Rdist
                           endif
                        enddo
                     enddo
                     if(present(pathOUTPUT))close(unit)
                  endif
                  !
                  deallocate(Rsorted,Rsorted_bkp,Rorder,Dist,DistList)
                  !
               enddo !isite
            enddo !ilayer
            !
         else
            !
            Hr = Hr_bulk
            deallocate(Hr_bulk)
            !
         endif
         !
         !Store Hr in w90 format
         if(present(pathOUTPUT))then
            !
            unit = free_unit()
            open(unit,file=reg(pathOUTPUT)//"Hr_built.DAT",form="formatted",status="unknown",position="rewind",action="write")
            write(unit,*)                  !skip first line
            write(unit,"(1I5)") Norb*Nsite !Number of Wannier orbitals
            write(unit,"(1I5)") Nwig       !Number of Wigner-Seitz vectors
            Qst = int(Nwig/W90NumCol)
            Rst = mod(Nwig,W90NumCol)
            !I'm going to find later the proper way to print the degeneracy only for the no-vanishing H(r)
            !do i=1,Qst
            !   write(unit,"("//str(W90NumCol)//"I5)")(nrdegwig(j+(i-1)*W90NumCol),j=1,W90NumCol)
            !enddo
            !if(Rst.ne.0)write(unit,"("//str(Rst)//"I5)")(nrdegwig(j+Qst*W90NumCol),j=1,Rst)
            do iwig=1,Nwig
               do iorb=1,Norb*Nsite
                  do jorb=1,Norb*Nsite
                     if(abs(Hr(iorb,jorb,iwig)).gt.0d0)write(unit,"(5I5,2F12.6)") Nvecwig(:,iwig), iorb, jorb, dreal(Hr(iorb,jorb,iwig)), dimag(Hr(iorb,jorb,iwig))
                  enddo
               enddo
            enddo
            close(unit)
            !
         endif
         !
         !FT Hr-->Hk
         call wannier_R2K(Nkpt3,kpt,Hr,Hk)
         deallocate(Hr)
         do ik=1,nkpt
            if(Norb.gt.1)call check_Hermiticity(Hk(:,:,ik),eps)
         enddo
         where(abs((Hk)).lt.eps) Hk=czero
         !
      endif
      !
      !Extract the off-diagonal dispersion
      if(Hetero%status)then
         !
         allocate(Hetero%tkz(Norb*Nsite_bulk,Norb*Nsite_bulk,Nkpt,Hetero%tzIndex(1):Hetero%tzIndex(2)));Hetero%tkz=czero
         do ilayer=1,Hetero%Nlayer-1
            !
            islab = ilayer + Hetero%Explicit(1) - 1
            na = 1 + (ilayer-1)*Norb*Nsite_bulk
            nb = ilayer*Norb*Nsite_bulk
            !
            do ik=1,Nkpt
               Hetero%tkz(:,:,ik,islab) = Hk(na:nb,na+Norb*Nsite_bulk:nb+Norb*Nsite_bulk,ik)
            enddo
            !
         enddo
         !
         if(Hetero%Explicit(1).ne.1)then
            !
            ilayer = Hetero%Explicit(1) + 1
            na = 1 + (ilayer-1)*Norb*Nsite_bulk
            nb = ilayer*Norb*Nsite_bulk
            !
            do ik=1,Nkpt
               Hetero%tkz(:,:,ik,Hetero%tzIndex(1)) = Hk(na:nb,na+Norb*Nsite_bulk:nb+Norb*Nsite_bulk,ik)
            enddo
            !
         endif
         if(Hetero%Explicit(2).ne.Hetero%Nslab)then
            !
            ilayer = Hetero%Explicit(2) - 2
            na = 1 + (ilayer-1)*Norb*Nsite_bulk
            nb = ilayer*Norb*Nsite_bulk
            !
            do ik=1,Nkpt
               Hetero%tkz(:,:,ik,Hetero%tzIndex(2)) = Hk(na:nb,na+Norb*Nsite_bulk:nb+Norb*Nsite_bulk,ik)
            enddo
            !
         endif
         !
      endif
      !
      if(present(pathOUTPUT))then
         !
         unit = free_unit()
         open(unit,file=reg(pathOUTPUT)//"Hk_built.DAT",form="formatted",status="unknown",position="rewind",action="write")
         write(unit,("(3I10)")) 1,Nkpt,Norb*Nsite
         do ik=1,Nkpt
            write(unit,("(2I6,3F14.8)")) 1,ik,kpt(:,ik)
            do iorb=1,Norb*Nsite
               do jorb=1,Norb*Nsite
                  write(unit,("(2I4,2E20.12)")) iorb,jorb,dreal(Hk(iorb,jorb,ik)),dimag(Hk(iorb,jorb,ik))
               enddo
            enddo
         enddo
         close(unit)
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
   !PURPOSE: Fill up the List likning indexes of the sum and diff of K-points
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
      integer,allocatable,intent(inout)     :: kptsum(:,:)
      integer,allocatable,intent(inout)     :: kptdif(:,:)
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
      !
      if(allocated(kptsum))deallocate(kptsum)
      if(allocated(kptdif))deallocate(kptdif)
      allocate(kptsum(Nkpt,Nkpt));kptsum=0
      allocate(kptdif(Nkpt,Nkpt));kptdif=0
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
   subroutine fill_smallk(kpt)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: kpt(:,:)
      !
      integer                               :: Nkpt,ik,idist,ismall
      real(8)                               :: Bvec(3),kvec(3)
      real(8),allocatable                   :: Kradius(:)
      integer,allocatable                   :: Korder(:)
      !
      !
      if(verbose)write(*,"(A)") "---- fill_smallk"
      if(.not.Lat_stored)stop "fill_smallk: Lattice positions not stored."
      !
      Nkpt = size(kpt,dim=2)
      !
      if(allocated(small_ik))deallocate(small_ik)
      if(allocated(KvecPolar))deallocate(KvecPolar)
      allocate(small_ik(Nkpt-1,2));small_ik=0
      allocate(KvecPolar(Nkpt-1,3));KvecPolar=0d0
      !
      allocate(Korder(Nkpt));Korder=0
      allocate(Kradius(Nkpt));Kradius=0d0
      do ik=1,Nkpt
         !
         kvec = kpt(:,ik) - nint(kpt(:,ik))
         Bvec = kvec(1)*Blat(:,1) + kvec(2)*Blat(:,2) + kvec(3)*Blat(:,3)
         !
         Kradius(ik) = sqrt(dble(dot_product(Bvec,Bvec)))
         !
      enddo
      !
      call sort_array(Kradius,Korder)
      !
      idist=0
      ismall=0
      if(verbose)write(*,"(A)")"     Smallest k vectors:"
      do ik=1,Nkpt
         !
         !skip the gamma point
         if(Kradius(Korder(ik)).ne.0d0)then
            !
            ismall = ismall + 1
            if((Kradius(Korder(ik))-Kradius(Korder(ik-1))).gt.eps) idist=idist+1
            !
            kvec = kpt(:,Korder(ik)) - nint(kpt(:,Korder(ik)))
            Bvec = kvec(1)*Blat(:,1) + kvec(2)*Blat(:,2) + kvec(3)*Blat(:,3)
            !
            small_ik(ismall,1) = Korder(ik)
            small_ik(ismall,2) = idist
            !
            !Polar coordinates sorted according to |k|
            KvecPolar(ismall,1) = sqrt(dble(dot_product(Bvec,Bvec)))                !r = sqrt( x^2 + y^2 + z^2 )
            KvecPolar(ismall,2) = atan2(Bvec(2),Bvec(1))                            !theta = atan( y/x )
            KvecPolar(ismall,3) = acos(Bvec(3)/sqrt(dble(dot_product(Bvec,Bvec))))  !phi = acos( z/r )
            !
            if(verbose.and.(small_ik(ismall,2).le.5)) write(*,"(A,I5,A,3F12.6,A,I3,A,3F12.6,3(A,1F20.10))") "     ikpt: ",small_ik(ismall,1)       &
                                                                         ,"   kpt: ",kpt(:,Korder(ik))         &
                                                                         ,"   dist: ",small_ik(ismall,2)       &
                                                                         ,"   Kvec: ",Bvec                     &
                                                                         ,"   r: ",KvecPolar(ismall,1)         &
                                                                         ,"   theta: ",KvecPolar(ismall,2)     &
                                                                         ,"   phi: ",KvecPolar(ismall,3)
            !
         endif
         !
      enddo
      deallocate(Kradius,Korder)
      !
      small_ik_stored=.true.
      !
     !if(verbose)then
     !   !
     !   !sorting dist
     !   call get_pattern(List,small_ik(:,2),listDim=DimList)
     !   write(*,*)new_line("A")//new_line("A")//"Lists: ",size(List,dim=1),size(DimList,dim=1)
     !   do ilist=1,size(List,dim=1)
     !      write(*,"(A,1F12.5,A,1I3,A,100I4)")"D:",dble(small_ik(List(ilist,1),2)),"Dim:", DimList(ilist)," List["//str(ilist)//",:] ",List(ilist,1:DimList(ilist))
     !   enddo
     !   deallocate(List,DimList)
     !   !
     !   !sorting R
     !   call get_pattern(List,Polar(:,1),eps,listDim=DimList)
     !   write(*,*)new_line("A")//new_line("A")//"Lists: ",size(List,dim=1),size(DimList,dim=1)
     !   do ilist=1,size(List,dim=1)
     !      write(*,"(A,1F12.5,A,1I3,A,100I4)")"R:",Polar(List(ilist,1),1),"Dim:", DimList(ilist)," List["//str(ilist)//",:] ",List(ilist,1:DimList(ilist))
     !   enddo
     !   deallocate(List,DimList)
     !   !
     !   !sorting theta
     !   call get_pattern(List,Polar(:,2),eps,listDim=DimList)
     !   write(*,*)new_line("A")//new_line("A")//"Lists: ",size(List,dim=1),size(DimList,dim=1)
     !   do ilist=1,size(List,dim=1)
     !      write(*,"(A,1F12.5,A,1I3,A,100I4)")"T:",Polar(List(ilist,1),2),"Dim:", DimList(ilist)," List["//str(ilist)//",:] ",List(ilist,1:DimList(ilist))
     !   enddo
     !   deallocate(List,DimList)
     !   !
     !   !sorting phi
     !   call get_pattern(List,Polar(:,3),eps,listDim=DimList)
     !   write(*,*)new_line("A")//new_line("A")//"Lists: ",size(List,dim=1),size(DimList,dim=1)
     !   do ilist=1,size(List,dim=1)
     !      write(*,"(A,1F12.5,A,1I3,A,100I4)")"P:",Polar(List(ilist,1),3),"Dim:", DimList(ilist)," List["//str(ilist)//",:] ",List(ilist,1:DimList(ilist))
     !   enddo
     !   deallocate(List,DimList)
     !   !
     !endif
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
                  Nvecwig_tmp(:,iwig) = (/ir1,ir2,ir3/)
                  Rvecwig_tmp(:,iwig) = matmul(Rlat,(/ir1,ir2,ir3/))
                  nrdegwig_tmp(iwig) = count(abs(distmin-dist(:)).le.epsWig)
                  radiuswig_tmp(iwig) = sqrt(dble(dot_product(rtmp,rtmp)))
                  if(all([ir1,ir2,ir3].eq.[0,0,0]))wig0=iwig
                  !
                  if(verbose)then
                     write(*,"(A,I10,A,I4)") "     iwig: ",iwig,"   deg:",nrdegwig_tmp(iwig)
                     write(*,"(A,3I8)") "     Nvecwig: ",Nvecwig_tmp(:,iwig)
                     write(*,"(A,3F)") "     Rvecwig: ",Rvecwig_tmp(:,iwig)
                     write(*,"(A,F)") "     R: ",radiuswig_tmp(iwig)
                  endif
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
         stop "calc_wignerseiz: nrdeg failed one of the lattice vectors might be wrong."
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
      complex(8),intent(inout)              :: mat_R(:,:)
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
      call assert_shape(mat_R,[Nsize1,Nwig],"wannier_K2R_d1","mat_R")
      !
      mat_R=czero
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
      complex(8),intent(inout)              :: mat_R(:,:,:)
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
      call assert_shape(mat_R,[Nsize1,Nsize2,Nwig],"wannier_K2R_d2","mat_R")
      !
      mat_R=czero
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
      complex(8),intent(inout)              :: mat_R(:,:,:,:)
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
      call assert_shape(mat_R,[Nsize1,Nsize2,Nsize3,Nwig],"wannier_K2R_d3","mat_R")
      !
      mat_R=czero
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
   subroutine wannier_R2K_d1(nkpt3_orig,kpt_intp,mat_R,mat_K)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_R(:,:)
      complex(8),intent(inout)              :: mat_K(:,:)
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
      call assert_shape(mat_K,[Nsize1,Nkpt_intp],"wannierinterpolation","mat_K")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_K=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,kpt_intp,Nvecwig,mat_K,mat_R,nrdegwig),&
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
               mat_K(i1,ik) = mat_K(i1,ik) + mat_R(i1,ir)*cfac
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
   subroutine wannier_R2K_d2(nkpt3_orig,kpt_intp,mat_R,mat_K)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_R(:,:,:)
      complex(8),intent(inout)              :: mat_K(:,:,:)
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
      call assert_shape(mat_K,[Nsize1,Nsize2,Nkpt_intp],"wannierinterpolation","mat_K")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_K=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,kpt_intp,Nvecwig,mat_K,mat_R,nrdegwig),&
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
                  mat_K(i1,i2,ik) = mat_K(i1,i2,ik) + mat_R(i1,i2,ir)*cfac
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
   subroutine wannier_R2K_d3(nkpt3_orig,kpt_intp,mat_R,mat_K)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)                    :: nkpt3_orig(:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_R(:,:,:,:)
      complex(8),intent(inout)              :: mat_K(:,:,:,:)
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
      call assert_shape(mat_K,[Nsize1,Nsize2,Nsize3,Nkpt_intp],"wannierinterpolation","mat_K")
      !
      ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
      mat_K=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nwig,Nkpt_intp,Nsize1,Nsize2,Nsize3,kpt_intp,Nvecwig,mat_K,mat_R,nrdegwig),&
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
                     mat_K(i1,i2,i3,ik) = mat_K(i1,i2,i3,ik) + mat_R(i1,i2,i3,ir)*cfac
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
      if(verbose)write(*,"(A)")"---- calc_Kpath"
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
      if(verbose)write(*,"(A)") "---- calc_Kplane"
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
   subroutine interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT                   &
                                                         ,filename,data_in,data_out    &
                                                         ,corrname,correction          &
                                                         ,doplane,Nkpt_Kside           &
                                                         ,hetero,store)
      !
      use parameters !WHY IS THIS WORKING?
      use utils_misc
      use linalg, only : eigh, inv, zeye
      implicit none
      !
      type(Lattice),intent(inout)           :: Lttc
      character(len=*),intent(in)           :: structure
      integer,intent(in)                    :: Nkpt_path
      character(len=*),intent(in),optional  :: pathOUTPUT
      character(len=*),intent(in),optional  :: filename
      complex(8),intent(in),optional        :: data_in(:,:,:)
      complex(8),intent(out),allocatable,optional :: data_out(:,:,:)
      character(len=*),intent(in),optional  :: corrname
      complex(8),intent(in),optional        :: correction(:,:,:)
      logical,intent(in),optional           :: doplane
      integer,intent(in),optional           :: Nkpt_Kside
      logical,intent(in),optional           :: store
      type(Heterostructures),intent(inout),optional :: hetero
      !
      character(len=256)                    :: path,label,filename_,corrname_
      integer                               :: ik,ikz,iorb,unit,ilayer
      integer                               :: Norb,Nkpt_Kside_,ikx,iky
      integer                               :: Nreal,iw,ndx
      real(8)                               :: kp,kx,ky,Bvec(3),wrealMax,eta,kz_cut
      complex(8),pointer                    :: data(:,:,:)
      complex(8),allocatable,target         :: data_tmp(:,:,:)
      complex(8),allocatable                :: invGf(:,:)
      logical                               :: Hamiltonian,doplane_,hetero_,printout,store_
      real                                  :: start,finish
      !Interp
      complex(8),allocatable                :: data_intp(:,:,:),dataZk(:,:,:)
      real(8),allocatable                   :: dataEk(:,:)
      !Plots
      real(8),allocatable                   :: Fk(:,:,:),Akw(:,:,:,:),wreal(:)
      complex(8),allocatable                :: zeta(:,:,:)
      !Hetero
      integer                               :: NbulkL,NbulkR
      integer                               :: Ln(2),Rn(2)
      complex(8),allocatable                :: Fk_kz(:,:,:,:),Akw_kz(:,:,:,:,:)
      complex(8),allocatable                :: Potential_L(:,:,:,:,:),Potential_R(:,:,:,:,:)
      real(8),allocatable                   :: Akw_print(:,:,:),Akw_kz_print(:,:,:,:)
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
         data_tmp = data_in
         Hamiltonian = .false.
      else
         label = "Hk"
         data_tmp = Lttc%Hk
         Hamiltonian = .true.
      endif
      data => data_tmp
      !
      doplane_=.false.
      if(present(doplane))doplane_=doplane
      !
      hetero_=.false.
      if(present(hetero))hetero_=Hetero%status
      !
      Norb = size(data,dim=1)
      call assert_shape(data,[Norb,Norb,Lttc%Nkpt],"interpolateHk2Path",reg(label))
      !
      store_=.true.
      if(present(store))store_=store
      !
      !static correction to the input data--------------------------------------
      corrname_="nonInt"
      if(present(correction))then
         call assert_shape(correction,[Norb,Norb,Lttc%Nkpt],"interpolateHk2Path","correction")
         data = data + correction
         corrname_="Corrected"
         if(present(corrname))corrname_=reg(corrname)
         write(*,"(A)")"     Correction: "//reg(corrname_)
         store_=.false.
      endif
      if(store_) write(*,"(A)")"     Storing Lttc interpolated attributes."
      !
      !
      !Create path along high-symmetry points-----------------------------------
      if(allocated(Lttc%kptpath))deallocate(Lttc%kptpath)
      if(allocated(Lttc%Kpathaxis))deallocate(Lttc%Kpathaxis)
      call calc_Kpath(Lttc%kptpath,reg(structure),Nkpt_path,Lttc%Kpathaxis,Lttc%KpathaxisPoints,hetero=hetero_)
      !
      !path in the bulk
      Lttc%Nkpt_path = size(Lttc%kptpath,dim=2)
      !
      !path for the Heterostructure: first Nkpt_path*SymmetryPoints then other Nkpt_path along Gamma-A
      if(hetero_)then
         if((Hetero%Norb*Lttc%Nsite).ne.Lttc%Norb) stop "interpolateHk2Path: Orbital dimension of Hk is not a multiple of the number of sites."
         Lttc%Nkpt_path = Lttc%Nkpt_path - Nkpt_path
      endif
      !
      !
      !Interpolate input data along path---------------------------------------------------
      allocate(data_intp(Norb,Norb,Lttc%Nkpt_path));data_intp=czero
      call cpu_time(start)
      call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),data,data_intp)
      call cpu_time(finish)
      write(*,"(A,F)") "     "//reg(label)//"(fullBZ) --> "//reg(label)//"(Kpath) cpu timing:", finish-start
      if(Hamiltonian.and.store_)then
         if(allocated(Lttc%Hk_path))deallocate(Lttc%Hk_path)
         Lttc%Hk_path = data_intp
      endif
      if(present(data_out))then
         if(allocated(data_out))deallocate(data_out)
         data_out = data_intp
      endif
      !
      !
      !Compute eigenvalues along path-------------------------------------------
      allocate(dataEk(Norb,Lttc%Nkpt_path));dataEk=0d0
      allocate(dataZk(Norb,Norb,Lttc%Nkpt_path));dataZk=czero
      dataZk = data_intp
      do ik=1,Lttc%Nkpt_path
         call eigh(dataZk(:,:,ik),dataEk(:,ik))
      enddo
      if(Hamiltonian.and.store_)then
         if(allocated(Lttc%Zk_path))deallocate(Lttc%Zk_path)
         if(allocated(Lttc%Ek_path))deallocate(Lttc%Ek_path)
         Lttc%Zk_path = dataZk
         Lttc%Ek_path = dataEk
         Lttc%pathStored = .true.
      endif
      deallocate(dataZk)
      !
      !Print eigenvalues along path
      if(printout)then
         path = reg(pathOUTPUT)//reg(filename_)//"_"//reg(corrname_)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Lttc%Nkpt_path
            write(unit,"(1I5,10000E20.12)") ik,Lttc%Kpathaxis(ik)/Lttc%Kpathaxis(Lttc%Nkpt_path),(dataEk(:,ik),iorb=1,Norb)
         enddo
         close(unit)
      endif
      write(*,"(A,I)") "     Total number of K-points along path:",Lttc%Nkpt_path
      wrealMax = 2.0*maxval(abs(dataEk))
      deallocate(dataEk)
      !
      !
      !Non-interacting spectral function along path-----------------------------
      if(hamiltonian)then
         !
         !Default parameters on the real frequency axis
         Nreal = 2000
         eta = wrealMax/200
         allocate(wreal(Nreal));wreal=0d0
         wreal = linspace(-wrealMax,+wrealMax,Nreal)
         !
         allocate(zeta(Norb,Norb,Nreal));zeta=czero
         do iorb=1,Norb
            do iw=1,Nreal
               zeta(iorb,iorb,iw) = dcmplx(wreal(iw),eta)
            enddo
         enddo
         !
         !Interpolate longitudinal tz along the path and compute potentials
         Ln=0;Rn=0
         if(hetero_)then
            !
            if(allocated(Hetero%tkz_path))deallocate(Hetero%tkz_path)
            allocate(Hetero%tkz_path(Hetero%Norb,Hetero%Norb,Lttc%Nkpt_path,Hetero%tzIndex(1):Hetero%tzIndex(2)));Hetero%tkz_path=czero
            do ilayer = Hetero%tzIndex(1),Hetero%tzIndex(2)
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),Hetero%tkz(:,:,:,ilayer),Hetero%tkz_path(:,:,:,ilayer))
            enddo
            !
            !Non-interacting potential to the left/upper side of the Heterostructure
            if(Hetero%Explicit(1).ne.1)then
               allocate(Potential_L(Hetero%Norb,Hetero%Norb,Nreal,Lttc%Nkpt_path,Nspin));Potential_L=czero
               call build_Potential(Potential_L,Hetero,Ln,NbulkL,zeta,Lttc%Hk_path,Hetero%tkz_path,"left",.true.)
               write(*,"(2(A,2I4))") "     Left potential (path) orbital lattice indexes: ",Ln(1),Ln(2)," thickness: ",NbulkL
            endif
            !
            !Non-interacting potential to the right/lower side of the Heterostructure
            if(Hetero%Explicit(2).ne.Hetero%Nslab)then
               allocate(Potential_R(Hetero%Norb,Hetero%Norb,Nreal,Lttc%Nkpt_path,Nspin));Potential_R=czero
               call build_Potential(Potential_R,Hetero,Rn,NbulkR,zeta,Lttc%Hk_path,Hetero%tkz_path,"right",.true.)
               write(*,"(2(A,2I4))") "     Right potential (path) orbital lattice indexes: ",Rn(1),Rn(2)," thickness: ",NbulkR
            endif
            !
         endif
         !
         !Compute non-interacting spectral function along path
         allocate(Akw(Norb,Norb,Nreal,Lttc%Nkpt_path));Akw=0d0
         allocate(invGf(Norb,Norb));invGf=czero
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(ik,iw,invGf)
         !$OMP DO
         do ik=1,Lttc%Nkpt_path
            do iw=1,Nreal
               !
               invGf = zeta(:,:,iw) - data_intp(:,:,ik)
               !
               if(allocated(Potential_L)) invGf(Ln(1):Ln(2),Ln(1):Ln(2)) = invGf(Ln(1):Ln(2),Ln(1):Ln(2)) - Potential_L(:,:,iw,ik,1)
               if(allocated(Potential_R)) invGf(Rn(1):Rn(2),Rn(1):Rn(2)) = invGf(Rn(1):Rn(2),Rn(1):Rn(2)) - Potential_R(:,:,iw,ik,1)
               !
               call inv(invGf)
               Akw(:,:,iw,ik) = dimag(invGf)
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         if(allocated(Potential_L))deallocate(Potential_L)
         if(allocated(Potential_R))deallocate(Potential_R)
         deallocate(zeta,invGf)
         !
         !Print non-interacting spectral function along path
         if(printout)then
            !
            !Normalization
            allocate(Akw_print(Norb,Nreal,Lttc%Nkpt_path));Akw_print=0d0
            do ik=1,Lttc%Nkpt_path
               do iorb=1,Norb
                  Akw_print(iorb,:,ik) = Akw(iorb,iorb,:,ik)/(sum(Akw(iorb,iorb,:,ik))*abs(wreal(2)-wreal(1)))
               enddo
            enddo
            !
            !print
            path = reg(pathOUTPUT)//"Akw_"//reg(label)//"_"//reg(corrname_)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,Lttc%Nkpt_path
               do iw=1,Nreal
                   write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik)/Lttc%Kpathaxis(Lttc%Nkpt_path),wreal(iw),(Akw_print(iorb,iw,ik),iorb=1,Norb)
               enddo
               write(unit,*)
            enddo
            close(unit)
            deallocate(Akw_print)
            !
         endif
         !
         !Compute non-interacting spectral function along the Gamma-A direction
         if(hetero_)then
            !
            allocate(Akw_kz(Hetero%Norb,Hetero%Norb,Nreal,Lttc%Nkpt_path,0:Nkpt_path));Akw_kz=czero
            do iw=1,Nreal
               call fill_Gamma_A(Akw(:,:,iw,:),Akw_kz(:,:,iw,:,:))
            enddo
            !
            !Normalization
            allocate(Akw_kz_print(Hetero%Norb,Nreal,Lttc%Nkpt_path,0:Nkpt_path));Akw_kz_print=0d0
            do ik=1,Lttc%Nkpt_path
               do ikz=0,Nkpt_path
                  do iorb=1,Hetero%Norb
                     Akw_kz_print(iorb,:,ik,ikz) = Akw_kz(iorb,iorb,:,ik,ikz)/(sum(Akw_kz(iorb,iorb,:,ik,ikz))*abs(wreal(2)-wreal(1)))
                  enddo
               enddo
            enddo
            deallocate(Akw_kz)
            !
            !Print non-interacting spectral function along path with the Gamma-A direction
            if(printout)then
               path = reg(pathOUTPUT)//"Akw_"//reg(label)//"_Hetero_"//reg(corrname_)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_path
                  do iw=1,Nreal
                      write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),wreal(iw),(Akw_kz_print(iorb,iw,ik,0),iorb=1,Hetero%Norb)
                  enddo
                  write(unit,*)
               enddo
               do ikz=1,Nkpt_path
                  do iw=1,Nreal
                      write(unit,"(1I5,200E20.12)") ik+ikz,Lttc%Kpathaxis(Lttc%Nkpt_path+ikz),wreal(iw),(Akw_kz_print(iorb,iw,Lttc%iq_gamma,ikz),iorb=1,Hetero%Norb)
                  enddo
                  write(unit,*)
               enddo
               close(unit)
            endif
            deallocate(Akw_kz_print)
            !
         endif
         deallocate(wreal,Akw)
         !
      endif
      deallocate(data_intp)
      !
      !
      !Non-interacting Fermi surface--------------------------------------------
      if(doplane_.and.hamiltonian)then
         !
         !
         !Create K-points inside the kx,ky plane
         Nkpt_Kside_ = 201
         if(present(Nkpt_Kside)) Nkpt_Kside_ = Nkpt_Kside
         !
         if(allocated(Lttc%kptPlane))deallocate(Lttc%kptPlane)
         call calc_Kplane(Lttc%kptPlane,Nkpt_Kside_)
         Lttc%Nkpt_Plane = size(Lttc%kptPlane,dim=2)
         !
         !
         !Interpolate hamiltonian inside the kx,ky plane
         allocate(data_intp(Norb,Norb,Lttc%Nkpt_Plane));data_intp=czero
         call cpu_time(start)
         call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,data,data_intp)
         call cpu_time(finish)
         write(*,"(A,F)") "     "//reg(label)//"(fullBZ) --> "//reg(label)//"(kx,ky) cpu timing:", finish-start
         if(Hamiltonian.and.store_)then
            if(allocated(Lttc%Hk_Plane))deallocate(Lttc%Hk_Plane)
            Lttc%Hk_Plane = data_intp
            Lttc%planeStored = .true.
         endif
         !
         !Create zeta array for compatibility
         eta = wrealMax/100 !same as before for Akw
         allocate(zeta(Norb,Norb,1));zeta=czero
         do iorb=1,Norb
            zeta(iorb,iorb,1) = dcmplx(EcutSheet,eta)
         enddo
         !
         !Interpolate longitudinal tz inside the kx,ky plane and compute potentials
         Ln=0;Rn=0
         if(hetero_)then
            !
            if(allocated(Hetero%tkz_Plane))deallocate(Hetero%tkz_Plane)
            allocate(Hetero%tkz_Plane(Hetero%Norb,Hetero%Norb,Lttc%Nkpt_Plane,Hetero%tzIndex(1):Hetero%tzIndex(2)));Hetero%tkz_Plane=czero
            do ilayer = Hetero%tzIndex(1),Hetero%tzIndex(2)
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,Hetero%tkz(:,:,:,ilayer),Hetero%tkz_Plane(:,:,:,ilayer))
            enddo
            !
            !Non-interacting potential to the left/upper side of the Heterostructure
            if(Hetero%Explicit(1).ne.1)then
               allocate(Potential_L(Hetero%Norb,Hetero%Norb,1,Lttc%Nkpt_Plane,Nspin));Potential_L=czero
               call build_Potential(Potential_L,Hetero,Ln,NbulkL,zeta,Lttc%Hk_Plane,Hetero%tkz_Plane,"left",.true.)
               write(*,"(2(A,2I4))") "     Left potential (plane) orbital lattice indexes: ",Ln(1),Ln(2)," thickness: ",NbulkL
            endif
            !
            !Non-interacting potential to the right/lower side of the Heterostructure
            if(Hetero%Explicit(2).ne.Hetero%Nslab)then
               allocate(Potential_R(Hetero%Norb,Hetero%Norb,1,Lttc%Nkpt_Plane,Nspin));Potential_R=czero
               call build_Potential(Potential_R,Hetero,Rn,NbulkR,zeta,Lttc%Hk_Plane,Hetero%tkz_Plane,"right",.true.)
               write(*,"(2(A,2I4))") "     Right potential (plane) orbital lattice indexes: ",Rn(1),Rn(2)," thickness: ",NbulkR
            endif
            !
         endif
         !
         !Compute non-interacting Fermi surface
         allocate(Fk(Norb,Norb,Lttc%Nkpt_Plane));Fk=0d0
         allocate(invGf(Norb,Norb));invGf=czero
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(ik,invGf)
         !$OMP DO
         do ik=1,Lttc%Nkpt_Plane
            !
            invGf = zeta(:,:,1) - data_intp(:,:,ik)
            !
            if(allocated(Potential_L)) invGf(Ln(1):Ln(2),Ln(1):Ln(2)) = invGf(Ln(1):Ln(2),Ln(1):Ln(2)) - Potential_L(:,:,1,ik,1)
            if(allocated(Potential_R)) invGf(Rn(1):Rn(2),Rn(1):Rn(2)) = invGf(Rn(1):Rn(2),Rn(1):Rn(2)) - Potential_R(:,:,1,ik,1)
            !
            call inv(invGf)
            Fk(:,:,ik) = -dimag(invGf)
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         if(allocated(Potential_L))deallocate(Potential_L)
         if(allocated(Potential_R))deallocate(Potential_R)
         deallocate(zeta,invGf,data_intp)
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
               write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Fk(iorb,iorb,ik),iorb=1,Norb)
               if(iky.eq.Nkpt_Kside_)write(unit,*)
            enddo
            close(unit)
         endif
         !
         !Compute non-interacting Fermi surface along the Gamma-A direction
         if(hetero_)then
            !
            allocate(Fk_kz(Hetero%Norb,Hetero%Norb,Lttc%Nkpt_Plane,0:Nkpt_path));Fk_kz=czero
            call fill_Gamma_A(Fk,Fk_kz)
            !
            !find the kz where to compute the Fermi-surface and print
            if(printout)then
               !
               !Gamma
               kz_cut = 0d0
               ikz = minloc(abs(Lttc%kptpath(3,1+Lttc%Nkpt_path:Lttc%Nkpt_path+Nkpt_path)-kz_cut),dim=1)
               !
               path = reg(pathOUTPUT)//"Fk_"//reg(label)//"_Hetero_"//reg(corrname_)//"_kz"//str(kz_cut,3)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_Plane
                  ikx = int(ik/(Nkpt_Kside_+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  iky = ik - (ikx-1)*Nkpt_Kside_      ; ky = (iky-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  Bvec = kx*Blat(:,1) + ky*Blat(:,2)
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(dreal(Fk_kz(iorb,iorb,ik,ikz)),iorb=1,Hetero%Norb)
                  if(iky.eq.Nkpt_Kside_)write(unit,*)
               enddo
               close(unit)
               !
               !A
               kz_cut = 0.5d0
               ikz = minloc(abs(Lttc%kptpath(3,1+Lttc%Nkpt_path:Lttc%Nkpt_path+Nkpt_path)-kz_cut),dim=1)
               !
               path = reg(pathOUTPUT)//"Fk_"//reg(label)//"_Hetero_"//reg(corrname_)//"_kz"//str(kz_cut,3)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_Plane
                  ikx = int(ik/(Nkpt_Kside_+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  iky = ik - (ikx-1)*Nkpt_Kside_      ; ky = (iky-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  Bvec = kx*Blat(:,1) + ky*Blat(:,2)
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(dreal(Fk_kz(iorb,iorb,ik,ikz)),iorb=1,Hetero%Norb)
                  if(iky.eq.Nkpt_Kside_)write(unit,*)
               enddo
               close(unit)
               !
               !half-way between Gamma-A
               kz_cut = 0.25d0
               ikz = minloc(abs(Lttc%kptpath(3,1+Lttc%Nkpt_path:Lttc%Nkpt_path+Nkpt_path)-kz_cut),dim=1)
               !
               path = reg(pathOUTPUT)//"Fk_"//reg(label)//"_Hetero_"//reg(corrname_)//"_kz"//str(kz_cut,3)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_Plane
                  ikx = int(ik/(Nkpt_Kside_+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  iky = ik - (ikx-1)*Nkpt_Kside_      ; ky = (iky-1)/dble(Nkpt_Kside_-1) - 0.5d0
                  Bvec = kx*Blat(:,1) + ky*Blat(:,2)
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(dreal(Fk_kz(iorb,iorb,ik,ikz)),iorb=1,Hetero%Norb)
                  if(iky.eq.Nkpt_Kside_)write(unit,*)
               enddo
               close(unit)
               !
            endif
            deallocate(Fk_kz)
            !
         endif
         deallocate(Fk)
         !
      endif
      if(associated(data))nullify(data)
      if(allocated(data_tmp))deallocate(data_tmp)
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
            Bvec = Lttc%kptpath(1,ik)*Blat(:,1) + Lttc%kptpath(2,ik)*Blat(:,2) + Lttc%kptpath(3,ik)*Blat(:,3)
            if(any(abs(Lttc%KpathaxisPoints-kp).lt.eps))then
               ndx = minloc(abs(Lttc%KpathaxisPoints-kp),dim=1)
               write(unit,"(1I5,200E20.12)") ik,kp,Lttc%kptpath(:,ik),Bvec,Lttc%KpathaxisPoints(ndx)
            else
               write(unit,"(1I5,200E20.12)") ik,kp,Lttc%kptpath(:,ik),Bvec
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
      !
      subroutine fill_Gamma_A(data_in,data_out)
         !
         implicit none
         !
         real(8),intent(in)                 :: data_in(:,:,:)
         complex(8),intent(inout)           :: data_out(:,:,:,0:)
         !
         integer                            :: ra,rb,ca,cb
         integer                            :: isite,jsite
         integer                            :: Nkpt_layer
         real(8)                            :: kR
         complex(8)                         :: cfac
         !
         Nkpt_layer = size(data_in,dim=3)
         if(Nkpt_layer.ne.size(data_out,dim=3)) stop "fill_Gamma_A: planar K-mesh does not coincide between layer-resolved and kz integrated."
         !
         data_out=czero
         !$OMP PARALLEL DEFAULT(PRIVATE),&
         !$OMP SHARED(Lttc,Nkpt_layer,Nkpt_path,Hetero,data_out,data_in)
         !$OMP DO
         do ik=1,Nkpt_layer
            do ikz=0,Nkpt_path
               do isite=1,Lttc%Nsite
                  do jsite=1,Lttc%Nsite
                     !
                     ra = 1+(isite-1)*Hetero%Norb ; rb = ra + Hetero%Norb-1
                     ca = 1+(jsite-1)*Hetero%Norb ; cb = ca + Hetero%Norb-1
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
      subroutine fill_Gamma_A_noPot(data_in,data_out)
         !
         implicit none
         !
         real(8),intent(in)                 :: data_in(:,:,:)
         complex(8),intent(inout)           :: data_out(:,:,:,0:)
         !
         integer                            :: ra,rb,ca,cb
         integer                            :: isite,jsite
         integer                            :: idist,jdist
         integer                            :: islab,jslab
         integer                            :: Nkpt_layer
         logical                            :: Explicit,BulkL,BulkR
         real(8)                            :: kR
         complex(8)                         :: cfac
         !
         Nkpt_layer = size(data_in,dim=3)
         if(Nkpt_layer.ne.size(data_out,dim=3)) stop "fill_Gamma_A_noPot: planar K-mesh does not coincide between layer-resolved and kz integrated."
         !
         !$OMP PARALLEL DEFAULT(PRIVATE),&
         !$OMP SHARED(Lttc,Nkpt_layer,Nkpt_path,Hetero,data_out,data_in)
         !$OMP DO
         do ik=1,Nkpt_layer
            do ikz=0,Nkpt_path
               !
               do islab=1,Hetero%Nslab
                  do jslab=1,Hetero%Nslab
                     !
                     Explicit = (islab.ge.Hetero%Explicit(1)) .and. (jslab.ge.Hetero%Explicit(1)) .and. &
                                (islab.le.Hetero%Explicit(2)) .and. (jslab.le.Hetero%Explicit(2))
                     if((.not.Explicit).and.abs(islab-jslab).gt.1)cycle
                     BulkL = (.not.Explicit) .and. (islab.lt.Hetero%Explicit(1)) .and. (jslab.lt.Hetero%Explicit(1))
                     BulkR = (.not.Explicit) .and. (islab.gt.Hetero%Explicit(2)) .and. (jslab.gt.Hetero%Explicit(2))
                     !
                     if(Explicit)then
                        !
                        isite = islab - (Hetero%Explicit(1)-1)
                        jsite = jslab - (Hetero%Explicit(1)-1)
                        !
                     elseif(BulkL)then
                        !
                        idist = abs(islab-Hetero%Explicit(1))
                        jdist = abs(jslab-Hetero%Explicit(1))
                        !
                        isite = islab + int((max(idist,jdist)+1)/2)*2
                        jsite = jslab + int((max(idist,jdist)+1)/2)*2
                        !
                     elseif(BulkR)then
                        !
                        idist = abs(islab-Hetero%Explicit(2))
                        jdist = abs(jslab-Hetero%Explicit(2))
                        !
                        isite = islab - int((max(idist,jdist)+1)/2)*2
                        jsite = jslab - int((max(idist,jdist)+1)/2)*2
                        !
                     endif
                     !
                     ra = 1+(isite-1)*Hetero%Norb ; rb = ra + Hetero%Norb-1
                     ca = 1+(jsite-1)*Hetero%Norb ; cb = ca + Hetero%Norb-1
                     !
                     kR = 2*pi * Lttc%kptpath(3,Lttc%Nkpt_path+ikz) * (islab-jslab)
                     cfac = dcmplx(cos(kR),+sin(kR))
                     !
                     data_out(:,:,ik,ikz) = data_out(:,:,ik,ikz) + data_in(ra:rb,ca:cb,ik)*cfac / Hetero%Nslab
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
      end subroutine fill_Gamma_A_noPot
      !
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


   !---------------------------------------------------------------------------!
   !PURPOSE: Interface with the potential subroutine in utils_mist.
   !         This subroutine takes the frequency mesh, Hk, vertical hoppig array
   !         and optionally the self-energy of the whole heterostructure and
   !         extract the slices needed to construct the embedding potentials
   !         via a continued fraction representation with the periodicity of
   !         the two extremals layers.
   !---------------------------------------------------------------------------!
   subroutine build_Potential(Potential,Hetero,ndx,Nbulk,zeta,Hk,tz,mode,paramagnet,Smats)
      !
      use parameters
      use linalg, only : inv, rotate
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)              :: Potential(:,:,:,:,:)
      integer,intent(out)                   :: Nbulk
      integer,intent(out)                   :: ndx(2)
      type(Heterostructures),intent(inout)  :: Hetero
      complex(8),intent(in)                 :: zeta(:,:,:)
      complex(8),intent(in)                 :: Hk(:,:,:)
      complex(8),intent(in)                 :: tz(:,:,:,:)
      character(len=*),intent(in)           :: mode
      logical,intent(in)                    :: paramagnet
      complex(8),intent(in),optional        :: Smats(:,:,:,:,:)
      !
      integer                               :: Norb,Nmats,Nkpt,ik
      complex(8),allocatable                :: zeta_(:,:,:),Hk_(:,:,:)
      complex(8),allocatable                :: ta(:,:,:),tb(:,:,:)
      complex(8),allocatable                :: Sa(:,:,:,:,:),Sb(:,:,:,:,:)
      complex(8),allocatable                :: Potential_loc(:,:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- build_Potential"
      !
      !
      if(.not.Hetero%status) stop "build_Potential: Heterostructure not properly initialized."
      if(size(Potential,dim=1).ne.Hetero%Norb) stop "build_Potential: Potential and Heterostructure have different orbital dimension."
      !
      Norb = size(Hk,dim=1)
      Nmats = size(Potential,dim=3)
      Nkpt = size(Potential,dim=4)
      !
      call assert_shape(Potential,[Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin],"build_Potential","Potential")
      call assert_shape(zeta,[Norb,Norb,Nmats],"build_Potential","zeta")
      call assert_shape(Hk,[Norb,Norb,Nkpt],"build_Potential","Hk")
      call assert_shape(tz,[Hetero%Norb,Hetero%Norb,Nkpt,(Hetero%tzIndex(2)-Hetero%tzIndex(1)+1)],"build_Potential","tz")
      !
      if(present(Smats))then
         if(size(Smats,dim=1).ne.Norb) stop "build_Potential: Smats has wrong orbital 1st dimension."
         if(size(Smats,dim=2).ne.Norb) stop "build_Potential: Smats has wrong orbital 2nd dimension."
         if(size(Smats,dim=3).ne.Nmats) stop "build_Potential: Smats has wrong Matsubara frequency mesh."
         if(size(Smats,dim=4).ne.Nkpt) stop "build_Potential: Smats has wrong k-point mesh."
      endif
      !
      Potential=czero
      allocate(ta(Hetero%Norb,Hetero%Norb,Nkpt));ta=czero
      allocate(tb(Hetero%Norb,Hetero%Norb,Nkpt));tb=czero
      allocate(Sa(Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin));Sa=czero
      allocate(Sb(Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin));Sb=czero
      !
      select case(reg(mode))
         case default
            !
            stop "build_Potential: Available modes: left, right."
            !
         case("left")
            !
            if(Hetero%Explicit(1).eq.1) stop "build_Potential: requested left potential but no bulk present."
            ndx(1) = 1
            ndx(2) = Hetero%Norb
            Nbulk = Hetero%Explicit(1)-1
            !
            !Connection-to and self-energy-of the first layer explicitly solved
            ta = tz(:,:,:,Hetero%Explicit(1)-1)
            if(present(Smats)) Sa = Smats(ndx(1):ndx(2),ndx(1):ndx(2),:,:,:)
            !
            !Connection-to and self-energy-of the second layer explicitly solved
            tb = tz(:,:,:,Hetero%Explicit(1))
            if(present(Smats)) Sb = Smats(ndx(1)+Hetero%Norb:ndx(2)+Hetero%Norb,ndx(1)+Hetero%Norb:ndx(2)+Hetero%Norb,:,:,:)
            !
         case("right")
            !
            if(Hetero%Explicit(2).eq.Hetero%Nslab) stop "build_Potential: requested right potential but no bulk present."
            ndx(1) = 1+ Norb - Hetero%Norb
            ndx(2) = Norb
            Nbulk = Hetero%Nslab-Hetero%Explicit(2)
            !
            !Connection-to and self-energy-of the last layer explicitly solved
            ta = tz(:,:,:,Hetero%Explicit(2))
            if(present(Smats)) Sa = Smats(ndx(1):ndx(2),ndx(1):ndx(2),:,:,:)
            !
            !Connection-to and self-energy-of the semi-last layer explicitly solved
            tb = tz(:,:,:,Hetero%Explicit(2)-1)
            if(present(Smats)) Sb = Smats(ndx(1)-Hetero%Norb:ndx(2)-Hetero%Norb,ndx(1)-Hetero%Norb:ndx(2)-Hetero%Norb,:,:,:)
            !
      end select
      !
      if(mod(Nbulk,2).ne.0)Nbulk=Nbulk+1
      !
      !Diagonal arrays of the first/last layer explicitly solved
      allocate(zeta_(Hetero%Norb,Hetero%Norb,Nmats));  zeta_ = zeta(ndx(1):ndx(2),ndx(1):ndx(2),:)
      allocate(Hk_(Hetero%Norb,Hetero%Norb,Nkpt))   ;  Hk_ = Hk(ndx(1):ndx(2),ndx(1):ndx(2),:)
      !
      call Embedding_ContinuedFraction(Potential,Nbulk,zeta_,Hk_,ta,tb,Sa,Sb,paramagnet)
      !
      allocate(Potential_loc(Hetero%Norb,Hetero%Norb,Nmats,Nspin));Potential_loc=czero
      do ik=1,Nkpt
         Potential_loc = Potential_loc + Potential(:,:,:,ik,:)/Nkpt
      enddo
      !
      !this is to be able to print it form main
      if(reg(mode).eq."left")then
         if(allocated(Hetero%P_L))deallocate(Hetero%P_L)
         Hetero%P_L = Potential_loc
      elseif(reg(mode).eq."right")then
         if(allocated(Hetero%P_R))deallocate(Hetero%P_R)
         Hetero%P_R = Potential_loc
      endif
      deallocate(Potential_loc)
      !
   end subroutine build_Potential


   !---------------------------------------------------------------------------!
   !PURPOSE: This subroutine takes the frequency mesh, Hk, vertical hoppig array
   !         and the self-energy of the first/last (depending on "mode") layers
   !         of the heterostructure and build the embedding potentials with a
   !         continued fraction representation.
   !---------------------------------------------------------------------------!
   subroutine Embedding_ContinuedFraction(Potential,Npot,zeta,Hk,tz_a,tz_b,Smats_a,Smats_b,paramagnet)
      !
      use parameters
      use utils_misc
      use linalg, only : inv, rotate
      implicit none
      !
      complex(8),intent(inout)              :: Potential(:,:,:,:,:)
      integer,intent(in)                    :: Npot
      complex(8),intent(in)                 :: zeta(:,:,:)
      complex(8),intent(in)                 :: Hk(:,:,:)
      complex(8),intent(in)                 :: tz_a(:,:,:)
      complex(8),intent(in)                 :: tz_b(:,:,:)
      complex(8),intent(in)                 :: Smats_a(:,:,:,:,:)
      complex(8),intent(in)                 :: Smats_b(:,:,:,:,:)
      logical,intent(in)                    :: paramagnet
      !
      complex(8),allocatable                :: invGbulk(:,:)
      complex(8),allocatable                :: Gbulk(:,:)
      complex(8),allocatable                :: Ptmp(:,:)
      complex(8),allocatable                :: tkz(:,:,:)
      complex(8),allocatable                :: Swks(:,:,:)
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iw,ik,ispin,ibulk,Pndx
      !
      !
      if(verbose)write(*,"(A)") "---- Embedding_ContinuedFraction"
      !
      !
      Norb = size(Potential,dim=1)
      Nmats = size(Potential,dim=3)
      Nkpt = size(Potential,dim=4)
      !
      call assert_shape(Potential,[Norb,Norb,Nmats,Nkpt,Nspin],"Embedding_ContinuedFraction","Potential")
      call assert_shape(zeta,[Norb,Norb,Nmats],"Embedding_ContinuedFraction","zeta")
      call assert_shape(Hk,[Norb,Norb,Nkpt],"Embedding_ContinuedFraction","Hk")
      call assert_shape(tz_a,[Norb,Norb,Nkpt],"Embedding_ContinuedFraction","tz_a")
      call assert_shape(tz_b,[Norb,Norb,Nkpt],"Embedding_ContinuedFraction","tz_b")
      call assert_shape(Smats_a,[Norb,Norb,Nmats,Nkpt,Nspin],"Embedding_ContinuedFraction","Smats_a")
      call assert_shape(Smats_b,[Norb,Norb,Nmats,Nkpt,Nspin],"Embedding_ContinuedFraction","Smats_b")
      !
      Potential = czero
      !
      !G and invG of each layer are constant
      allocate(invGbulk(Norb,Norb));invGbulk=czero
      allocate(Gbulk(Norb,Norb));Gbulk=czero
      allocate(Ptmp(Norb,Norb));Ptmp=czero
      allocate(tkz(Norb,Norb,0:1));tkz=czero
      allocate(Swks(Norb,Norb,0:1));Swks=czero
      spinPloop: do ispin=1,Nspin
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(ispin,Nmats,Nkpt,zeta,Hk,Smats_b,Smats_a,Npot,tz_a,tz_b,Potential),&
         !$OMP PRIVATE(ik,iw,ibulk,invGbulk,Gbulk,Ptmp,tkz,Swks,Pndx)
         !$OMP DO
         do iw=1,Nmats
            do ik=1,Nkpt
               !
               !connecting hopping: even Nbulk--> ta...tb odd Nbulk--> ta...ta
               tkz(:,:,1) = tz_a(:,:,ik)
               tkz(:,:,0) = tz_b(:,:,ik)
               !
               !self-energy repetition: even Nbulk--> Sb...Sa odd Nbulk--> Sb...Sb
               Swks(:,:,1) = Smats_b(:,:,iw,ik,ispin)
               Swks(:,:,0) = Smats_a(:,:,iw,ik,ispin)
               !
               !first t*G*t is the farthest layer
               Pndx = mod(Npot,2)
               invGbulk = zeta(:,:,iw) - Hk(:,:,ik) - Swks(:,:,Pndx)
               Gbulk = invGbulk
               call inv(Gbulk)
               Potential(:,:,iw,ik,ispin) = rotate(Gbulk,tkz(:,:,Pndx))
               !
               !all the other
               do ibulk=2,Npot
                  !
                  Pndx = mod(Npot-ibulk+1,2)
                  Ptmp = zeta(:,:,iw) - Hk(:,:,ik) - Swks(:,:,Pndx) - Potential(:,:,iw,ik,ispin)
                  call inv(Ptmp)
                  Potential(:,:,iw,ik,ispin) = rotate(Ptmp,tkz(:,:,Pndx))
                  !
               enddo
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         if(paramagnet)then
            Potential(:,:,:,:,Nspin) = Potential(:,:,:,:,1)
            exit spinPloop
         endif
      enddo spinPloop
      deallocate(invGbulk,Gbulk,Ptmp,tkz,Swks)
      !
   end subroutine Embedding_ContinuedFraction


end module crystal
