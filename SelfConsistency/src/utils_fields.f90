module utils_fields

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface clear_attributes
      module procedure clear_attributes_Fermion
      module procedure clear_attributes_Boson
   end interface clear_attributes

   interface loc2imp
      module procedure loc2imp_Fermionic
      module procedure loc2imp_Matrix
      module procedure loc2imp_Bosonic
   end interface loc2imp

   interface imp2loc
      module procedure imp2loc_Fermionic
      module procedure imp2loc_Matrix
      module procedure imp2loc_Bosonic
   end interface imp2loc

   interface MergeFields
      module procedure MergeSelfEnergy
      module procedure MergePolarization
   end interface MergeFields

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: FermionicKsum
   public :: BosonicKsum
   public :: AllocateLattice
   public :: DeallocateLattice
   public :: AllocateFermionicField
   public :: DeallocateFermionicField
   public :: AllocateBosonicField
   public :: DeallocateBosonicField
   public :: loc2imp
   public :: imp2loc
   public :: MergeFields
   public :: clear_attributes

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Fill the local attributes of a Fermionic Field
   !---------------------------------------------------------------------------!
   subroutine FermionicKsum(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      integer                               :: ik,in,ispin
      !
      if(.not.G%status)stop "FermionicKsum. Field not properly initialized."
      if(G%Nkpt.eq.0)stop "FermionicKsum. Field k dependent attributes not properly initialized."
      !
      G%N_s=czero
      do ispin=1,Nspin
         do ik=1,G%Nkpt
            G%N_s(:,:,ispin) = G%N_s(:,:,ispin) + G%N_ks(:,:,ik,ispin)/G%Nkpt
         enddo
      enddo
      !
      if(G%Npoints.ne.0)then
         G%ws=czero
         do ispin=1,Nspin
            do ik=1,G%Nkpt
               do in=1,G%Npoints
                  G%ws(:,:,in,ispin) = G%ws(:,:,in,ispin) + G%wks(:,:,in,ik,ispin)/G%Nkpt
               enddo
            enddo
         enddo
      endif
      !
      G%local_filled=.true.
      !
   end subroutine FermionicKsum


   !---------------------------------------------------------------------------!
   !PURPOSE: Fill the local attributes of a Bosonic Field
   !---------------------------------------------------------------------------!
   subroutine BosonicKsum(W)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      integer                               :: ik,in
      !
      if(.not.W%status)stop "BosonicKsum. Field not properly initialized."
      if(W%Nkpt.eq.0)stop "BosonicKsum. Field k dependent attributes not properly initialized."
      !
      W%screened_local=czero
      do ik=1,W%Nkpt
         do in=1,W%Npoints
            W%screened_local(:,:,in) = W%screened_local(:,:,in) + W%screened(:,:,in,ik)/W%Nkpt
         enddo
      enddo
      !
      if(allocated(W%bare_local).and.allocated(W%bare))then
         do ik=1,W%Nkpt
            W%bare_local = W%bare_local + W%bare(:,:,ik)/W%Nkpt
         enddo
      endif
      !
      W%local_filled=.true.
      !
   end subroutine BosonicKsum


   !---------------------------------------------------------------------------!
   !PURPOSE: Allocate/deallocate Lattice attributes in a consistent way
   !---------------------------------------------------------------------------!
   subroutine AllocateLattice(lttc,Norb,Nkpt,name)
      use parameters
      implicit none
      type(Lattice),intent(inout)           :: lttc
      integer,intent(in)                    :: Norb,Nkpt
      character(len=*),intent(in),optional  :: name
      !
      if(present(name)) write(*,"(A)") "Allocation of "//trim(name)
      if(lttc%status) stop "AllocateLattice: container already allocated."
      !
      if(allocated(lttc%kpt))deallocate(lttc%kpt)
      allocate(lttc%kpt(2,Nkpt));lttc%kpt=0d0
      !
      if(allocated(lttc%Hk))deallocate(lttc%Hk)
      allocate(lttc%Hk(Norb,Norb,Nkpt));lttc%Hk=czero
      !
      if(allocated(lttc%Zk))deallocate(lttc%Zk)
      allocate(lttc%Zk(Norb,Norb,Nkpt));lttc%Zk=czero
      !
      if(allocated(lttc%Ek))deallocate(lttc%Ek)
      allocate(lttc%Ek(Norb,Nkpt));lttc%Ek=0d0
      !
      if(allocated(lttc%kptPos))deallocate(lttc%kptPos)
      allocate(lttc%kptPos(Nkpt));lttc%kptPos=0
      !
      if(allocated(lttc%kptsum))deallocate(lttc%kptsum)
      allocate(lttc%kptsum(Nkpt,Nkpt));lttc%kptsum=0
      !
      if(allocated(lttc%kptdif))deallocate(lttc%kptdif)
      allocate(lttc%kptdif(Nkpt,Nkpt));lttc%kptdif=0
      !
      if(allocated(lttc%small_ik))deallocate(lttc%small_ik)
      allocate(lttc%small_ik(12,2));lttc%small_ik=0
      !
      lttc%status=.true.
      !
   end subroutine AllocateLattice
   !
   subroutine DeallocateLattice(lttc,name)
      use parameters
      implicit none
      type(Lattice),intent(inout)           :: lttc
      character(len=*),intent(in),optional  :: name
      !
      if(present(name)) write(*,"(A)") "Deallocation of "//trim(name)
      if(.not.lttc%status) stop "AllocateLattice: container is unallocated."
      !
      if(allocated(lttc%kpt))deallocate(lttc%kpt)
      if(allocated(lttc%Hk))deallocate(lttc%Hk)
      if(allocated(lttc%Zk))deallocate(lttc%Zk)
      if(allocated(lttc%Ek))deallocate(lttc%Ek)
      if(allocated(lttc%kptPos))deallocate(lttc%kptPos)
      if(allocated(lttc%kptsum))deallocate(lttc%kptsum)
      if(allocated(lttc%kptdif))deallocate(lttc%kptdif)
      if(allocated(lttc%kprint))deallocate(lttc%kprint)
      if(allocated(lttc%small_ik))deallocate(lttc%small_ik)
      lttc%Nkpt3=0
      lttc%Nkpt=0
      lttc%Nkpt_irred=0
      lttc%Norb=0
      lttc%mu=0d0
      lttc%density=0d0
      lttc%UseDisentangledBS=.false.
      lttc%status=.false.
      !
   end subroutine DeallocateLattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Allocate/deallocate Fermionic attributes in a consistent way
   !---------------------------------------------------------------------------!
   subroutine AllocateFermionicField(G,Norb,Npoints,Nkpt,Nsite,name)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      integer,intent(in)                    :: Norb,Npoints
      integer,intent(in),optional           :: Nkpt
      integer,intent(in),optional           :: Nsite
      character(len=*),intent(in),optional  :: name
      integer                               :: Nkpt_
      !
      Nkpt_=0
      if(present(Nkpt))Nkpt_=Nkpt
      if(present(name)) write(*,"(A)") "Allocation of "//trim(name)
      if(G%status) stop "AllocateFermionicField: container already allocated."
      if(Norb.eq.0) stop "AllocateFermionicField: Norb not defined."
      if(Npoints.eq.0) write(*,"(A)") "AllocateFermionicField: frequency dependent attributes are not going to be allocated."
      if(Nkpt_.eq.0) write(*,"(A)") "AllocateFermionicField: K-dependent attributes are not going to be allocated."
      !
      if(allocated(G%N_s))deallocate(G%N_s)
      if(allocated(G%N_ks))deallocate(G%N_ks)
      if(allocated(G%ws))deallocate(G%ws)
      if(allocated(G%wks))deallocate(G%wks)
      !
      allocate(G%N_s(Norb,Norb,Nspin));G%N_s=czero
      !
      if((Npoints.eq.0).and.(Nkpt_.ne.0))then
         allocate(G%N_ks(Norb,Norb,Nkpt_,Nspin));G%N_ks=czero
      elseif((Npoints.ne.0).and.(Nkpt_.eq.0))then
         allocate(G%ws(Norb,Norb,Npoints,Nspin));G%ws=czero
      else
         allocate(G%N_ks(Norb,Norb,Nkpt_,Nspin));G%N_ks=czero
         allocate(G%ws(Norb,Norb,Npoints,Nspin));G%ws=czero
         allocate(G%wks(Norb,Norb,Npoints,Nkpt_,Nspin));G%wks=czero
      endif
      !
      if(present(Nsite))G%Nsite=Nsite
      G%status=.true.
      !
   end subroutine AllocateFermionicField
   !
   subroutine DeallocateFermionicField(G,name)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      character(len=*),intent(in),optional  :: name
      !
      if(present(name)) write(*,"(A)") "Deallocation of "//trim(name)
      if(.not.G%status) stop "DeallocateFermionicField: container is unallocated."
      !
      if(allocated(G%N_s))deallocate(G%N_s)
      if(allocated(G%N_ks))deallocate(G%N_ks)
      if(allocated(G%ws))deallocate(G%ws)
      if(allocated(G%wks))deallocate(G%wks)
      !
      G%Norb=0
      G%Npoints=0
      G%Nkpt=0
      G%Nsite=0
      G%Beta=0d0
      G%mu=0d0
      G%status=.false.
      !
   end subroutine DeallocateFermionicField


   !---------------------------------------------------------------------------!
   !PURPOSE: Allocate/deallocate Bosonic attributes in a consistent way
   !---------------------------------------------------------------------------!
   subroutine AllocateBosonicField(W,Norb,Npoints,Nkpt,Nsite,name,no_bare)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      integer,intent(in)                    :: Norb,Npoints
      integer,intent(in),optional           :: Nkpt
      integer,intent(in),optional           :: Nsite
      character(len=*),intent(in),optional  :: name
      logical,intent(in),optional           :: no_bare
      integer                               :: Nbp,Nkpt_
      logical                               :: no_bare_
      !
      Nbp=Norb**2
      Nkpt_=0
      if(present(Nkpt))Nkpt_=Nkpt
      if(present(name)) write(*,"(A)") "Allocation of "//trim(name)
      if(W%status) stop "AllocateBosonicField: container already allocated."
      if(Nbp.eq.0) stop "AllocateBosonicField: Nbp not defined."
      if(Npoints.eq.0) stop "AllocateBosonicField: Npoints not defined."
      if(Nkpt_.eq.0) write(*,"(A)") "AllocateBosonicField: K-dependent attributes are not going to be allocated."
      no_bare_=.false.
      if(present(no_bare))no_bare_=no_bare
      if(no_bare_) write(*,"(A)") "AllocateBosonicField: the bare attributes are not going to be allocated."
      !
      if(allocated(W%bare_local))deallocate(W%bare_local)
      if(allocated(W%screened_local))deallocate(W%screened_local)
      if(allocated(W%bare))deallocate(W%bare)
      if(allocated(W%screened))deallocate(W%screened)
      !
      allocate(W%screened_local(Nbp,Nbp,Npoints));W%screened_local=czero
      if(Nkpt_.ne.0)then
         allocate(W%screened(Nbp,Nbp,Npoints,Nkpt_));W%screened=czero
      endif
      !
      if(.not.no_bare_)then
         allocate(W%bare_local(Nbp,Nbp));W%bare_local=czero
         if(Nkpt_.ne.0)then
            allocate(W%bare(Nbp,Nbp,Nkpt_));W%bare=czero
         endif
      endif
      !
      if(present(Nsite))W%Nsite=Nsite
      W%status=.true.
      !
   end subroutine AllocateBosonicField
   !
   subroutine DeallocateBosonicField(W,name)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      character(len=*),intent(in),optional  :: name
      !
      if(present(name)) write(*,"(A)") "Deallocation of "//trim(name)
      if(.not.W%status) stop "DeallocateBosonicField: container is unallocated."
      !
      if(allocated(W%bare_local))deallocate(W%bare_local)
      if(allocated(W%screened_local))deallocate(W%screened_local)
      if(allocated(W%bare))deallocate(W%bare)
      if(allocated(W%screened))deallocate(W%screened)
      !
      W%Nbp=0
      W%Npoints=0
      W%Nkpt=0
      W%Nsite=0
      W%iq_gamma=-1
      W%Beta=0d0
      W%status=.false.
      !
   end subroutine DeallocateBosonicField


   !---------------------------------------------------------------------------!
   !PURPOSE: Extract from a lattice local projection the subset of orbital
   !         indexes of corresponging to a given site
   !---------------------------------------------------------------------------!
   subroutine loc2imp_Fermionic(Gimp,Gloc,orbs,sitename,U)
      !
      use parameters
      use linalg, only : rotate
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gimp
      type(FermionicField),intent(in)       :: Gloc
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:)
      !
      integer                               :: ip,ispin
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      write(*,"(A)") "--- loc2imp(F) ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Extraction of Gloc for site: "//trim(sitename)
      if(.not.Gloc%status) stop "loc2imp(F): Gloc not properly initialized."
      if(.not.Gimp%status) stop "loc2imp(F): Gimp not properly initialized."
      if(Gloc%Norb.eq.0) stop "loc2imp(F): Norb of Gloc not defined."
      if(Gimp%Norb.eq.0) stop "loc2imp(F): Norb of Gimp not defined."
      if(Gloc%Npoints.eq.0) stop "loc2imp(F): Npoints of Gloc not defined."
      if(Gimp%Npoints.eq.0) stop "loc2imp(F): Npoints of Gimp not defined."
      if(Gimp%Beta.ne.Gloc%Beta) stop "Gimp2Gimp: Gimp and Gloc have different beta."
      if(Gimp%Npoints.ne.Gloc%Npoints) stop "Gimp2Gimp: Gimp and Gloc have different number of Matsubara points."
      if(Gimp%Nkpt.ne.0) stop "loc2imp(F): Gimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Gloc%ws)) stop "loc2imp(F): Gloc local projection not allocated."
      if(.not.allocated(Gimp%ws)) stop "loc2imp(F): Gimp local projection not allocated."
      if(size(orbs).ne.Gimp%Norb) stop "loc2imp(F): can't fit the requested orbitals inside Gimp."
      if(size(orbs).gt.Gloc%Norb) stop "loc2imp(F): number of requested orbitals greater than Gloc size."
      if(present(U))then
         write(*,"(A)") "The local orbital space will be rotated during extraction."
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if(size(U,dim=1).ne.size(orbs)) stop "Rotation matrix has the wrong dimension."
         if(size(U,dim=1).ne.3) write(*,"(A)") "Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
      endif
      !
      do i_imp=1,size(orbs)
         do j_imp=1,size(orbs)
            !
            i_loc = orbs(i_imp)
            j_loc = orbs(j_imp)
            !
            do ispin=1,Nspin
               do ip=1,Gimp%Npoints
                  Gimp%ws(i_imp,j_imp,ip,ispin) = Gloc%ws(i_loc,j_loc,ip,ispin)
               enddo
            enddo
            !
         enddo
      enddo
      !
      if(present(U))then
         do ispin=1,Nspin
            do ip=1,Gimp%Npoints
               Gimp%ws(:,:,ip,ispin) = rotate(Gimp%ws(:,:,ip,ispin),U)
            enddo
         enddo
      endif
      !
   end subroutine loc2imp_Fermionic
   !
   subroutine loc2imp_Matrix(Oimp,Oloc,orbs,sitename,U)
      !
      use parameters
      use linalg, only : rotate
      implicit none
      !
      complex(8),allocatable,intent(inout)  :: Oimp(:,:)
      complex(8),allocatable,intent(in)     :: Oloc(:,:)
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:)
      !
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      write(*,"(A)") "--- loc2imp(O) ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Extraction of Operator for site: "//trim(sitename)
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "loc2imp(O): Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "loc2imp(O): Oimp not square."
      if(size(orbs).ne.size(Oimp,dim=1)) stop "loc2imp(O): can't fit the requested orbitals inside Oimp."
      if(size(orbs).gt.size(Oloc,dim=1)) stop "loc2imp(O): number of requested orbitals greater than Oloc size."
      if(present(U))then
         write(*,"(A)") "The local orbital space will be rotated during extraction."
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if(size(U,dim=1).ne.size(orbs)) stop "Rotation matrix has the wrong dimension."
         if(size(U,dim=1).ne.3) write(*,"(A)") "Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
      endif
      !
      do i_imp=1,size(orbs)
         do j_imp=1,size(orbs)
            !
            i_loc = orbs(i_imp)
            j_loc = orbs(j_imp)
            !
            Oimp(i_imp,j_imp) = Oloc(i_loc,j_loc)
            !
         enddo
      enddo
      !
      if(present(U)) Oimp = rotate(Oimp,U)
      !
   end subroutine loc2imp_Matrix
   !
   subroutine loc2imp_Bosonic(Wimp,Wloc,orbs,sitename)
      !
      use parameters
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wimp
      type(BosonicField),intent(in)         :: Wloc
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      !
      integer                               :: ip
      integer                               :: Norb_imp,Norb_loc
      integer                               :: ib_imp,jb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: i_imp,j_imp,k_imp,l_imp
      !
      !
      write(*,"(A)") "--- loc2imp(B) ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Extraction of Wloc for site: "//trim(sitename)
      if(.not.Wloc%status) stop "loc2imp(B): Wloc not properly initialized."
      if(.not.Wimp%status) stop "loc2imp(B): Wimp not properly initialized."
      if(Wloc%Nbp.eq.0) stop "loc2imp(B): Norb of Wloc not defined."
      if(Wimp%Nbp.eq.0) stop "loc2imp(B): Norb of Wimp not defined."
      if(Wloc%Npoints.eq.0) stop "loc2imp(B): Npoints of Wloc not defined."
      if(Wimp%Npoints.eq.0) stop "loc2imp(B): Npoints of Wimp not defined."
      if(Wimp%Beta.ne.Wloc%Beta) stop "Wimp2Wimp: Wimp and Wloc have different beta."
      if(Wimp%Npoints.ne.Wloc%Npoints) stop "Wimp2Wimp: Wimp and Wloc have different number of Matsubara points."
      if(Wimp%Nkpt.ne.0) stop "loc2imp(B): Wimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Wloc%screened_local)) stop "loc2imp(B): Wloc screened_local attribute not allocated."
      if(.not.allocated(Wloc%bare_local)) stop "loc2imp(B): Wloc bare_local attribute not allocated."
      if(.not.allocated(Wimp%screened_local)) stop "loc2imp(B): Wimp screened_local attribute not allocated."
      if(.not.allocated(Wimp%bare_local)) stop "loc2imp(B): Wimp bare_local attribute not allocated."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      !
      if(size(orbs).ne.Norb_imp) stop "loc2imp(B): can't fit the requested orbitals inside Wimp."
      if(size(orbs).gt.Norb_loc) stop "loc2imp(B): number of requested orbitals greater than Wloc size."
      !
      do i_imp=1,Norb_imp
         do j_imp=1,Norb_imp
            do k_imp=1,Norb_imp
               do l_imp=1,Norb_imp
                  !
                  ! bosonic indexes of the impurity
                  ib_imp = k_imp + Norb_imp*(i_imp-1)
                  jb_imp = l_imp + Norb_imp*(j_imp-1)
                  !
                  ! mapping
                  i_loc = orbs(i_imp)
                  j_loc = orbs(j_imp)
                  k_loc = orbs(k_imp)
                  l_loc = orbs(l_imp)
                  !
                  ! corresponding bosonic indexes on the lattice
                  ib_loc = k_loc + Norb_loc*(i_loc-1)
                  jb_loc = l_loc + Norb_loc*(j_loc-1)
                  !
                  Wimp%bare_local(ib_imp,jb_imp) =  Wloc%bare_local(ib_loc,jb_loc)
                  do ip=1,Wimp%Npoints
                     Wimp%screened_local(ib_imp,jb_imp,ip) = Wloc%screened_local(ib_loc,jb_loc,ip)
                  enddo
                  !
               enddo
            enddo
         enddo
      enddo
      !
   end subroutine loc2imp_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE: Expand the the subset of orbital indexes of corresponging to a
   !         given site to the full lattice quantity
   !---------------------------------------------------------------------------!
   subroutine imp2loc_Fermionic(Gloc,Gimp,orbs,sitename,U,expand,AFM)
      use parameters
      use linalg, only : rotate
      implicit none
      type(FermionicField),intent(inout)    :: Gloc
      type(FermionicField),intent(in)       :: Gimp
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:,:)
      logical,intent(in),optional           :: expand
      logical,intent(in),optional           :: AFM
      !
      complex(8),allocatable                :: Gtmp(:,:,:,:,:)
      integer                               :: ip,ispin,ispin_DMFT,isite,Nsite,shift
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      logical                               :: expand_,AFM_
      !
      !
      if(verbose)write(*,"(A)") "--- imp2loc(F) ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Insertion of Gimp of site: "//trim(sitename)
      if(.not.Gloc%status) stop "imp2loc(F): Gloc not properly initialized."
      if(.not.Gimp%status) stop "imp2loc(F): Gimp not properly initialized."
      if(Gloc%Norb.eq.0) stop "imp2loc(F): Norb of Gloc not defined."
      if(Gimp%Norb.eq.0) stop "imp2loc(F): Norb of Gimp not defined."
      if(Gloc%Npoints.eq.0) stop "imp2loc(F): Npoints of Gloc not defined."
      if(Gimp%Npoints.eq.0) stop "imp2loc(F): Npoints of Gimp not defined."
      if(Gloc%Beta.ne.Gimp%Beta) stop "imp2loc(F): Gimp and Gloc have different beta."
      if(Gloc%Npoints.ne.Gimp%Npoints) stop "imp2loc(F): Gimp and Gloc have different number of Matsubara points."
      if(Gimp%Nkpt.ne.0) stop "imp2loc(F): Gimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Gloc%ws)) stop "imp2loc(F): Gloc local projection not allocated."
      if(.not.allocated(Gimp%ws)) stop "imp2loc(F): Gimp local projection not allocated."
      if(size(orbs).ne.Gimp%Norb) stop "imp2loc(F): can't fit the requested orbitals from Gimp."
      if(size(orbs).gt.Gloc%Norb) stop "imp2loc(F): number of requested orbitals greater than Gloc size."
      !
      Norb_loc = Gloc%Norb
      Norb_imp = Gimp%Norb
      expand_=.false.
      if(present(expand))expand_=expand
      AFM_=.false.
      if(present(AFM))AFM_=AFM
      Nsite=1
      !
      if(present(U))then
         write(*,"(A)") "The local orbital space will be rotated during extraction."
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if(size(U,dim=1).ne.size(orbs)) stop "Rotation matrix has the wrong dimension."
         if(size(U,dim=3).ne.Gloc%Nsite) stop "Number of rotation matrices and number of sites does not match."
         if(size(U,dim=1).ne.3) write(*,"(A)") "Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
      endif
      !
      if(AFM)then
         if(Gloc%Nsite.ne.2) stop "AFM is implemented only for a two site lattice."
         if(Norb/Gimp%Norb.ne.2) stop "Lattice indexes are not twice the impurity ones."
         if(expand) stop "AFM condition and expansion to real space not yet implemented."
         Nsite = 2
      endif
      !
      if(expand)then
         if(mod(Norb_loc,size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of the lattice field."
         if(AFM) stop "Expansion to real space and AFM condition not yet implemented."
         write(*,"(A)") "The impurity field will be expanded to match the lattice orbital space."
         Nsite = Gloc%Nsite
      else
         write(*,"(A,15I3)") "Impurity field will be inserted into the lattice orbital indexes: ",orbs
         Nsite = 1
      endif
      !
      allocate(Gtmp(Gimp%Norb,Gimp%Norb,Gimp%Npoints,Nspin,Nsite));Gtmp=czero
      !
      ! Rotating either one site or all of them depending on expand_
      do isite=1,Nsite
         !
         if(present(U))then
            do ispin=1,Nspin
               do ip=1,Gimp%Npoints
                  Gtmp(:,:,ip,ispin,isite) = rotate(Gimp%ws(:,:,ip,ispin),U(:,:,isite))
               enddo
            enddo
         else
            Gtmp(:,:,:,:,isite) = Gimp%ws
         endif
         !
      enddo
      !
      do isite=1,Nsite
         !
         ! only two possible arrangements
         if(abs(orbs(2)-orbs(1)).eq.1)then
            shift = size(orbs)*(isite-1)
         elseif(abs(orbs(2)-orbs(1)).eq.Nsite)then
            shift = isite-1
         endif
         !
         do i_imp=1,size(orbs)
            do j_imp=1,size(orbs)
               !
               i_loc = orbs(i_imp) + shift
               j_loc = orbs(j_imp) + shift
               !
               do ispin=1,Nspin
                  !
                  ispin_DMFT=ispin
                  if(isite.eq.2) ispin_DMFT=int(Nspin/ispin)
                  !
                  do ip=1,Gimp%Npoints
                     Gloc%ws(i_loc,j_loc,ip,ispin) = Gtmp(i_imp,j_imp,ip,ispin,isite)
                  enddo
               enddo
               !
            enddo
         enddo
         !
      enddo
      deallocate(Gtmp)
      !
   end subroutine imp2loc_Fermionic
   !
   subroutine imp2loc_Matrix(Oloc,Oimp,orbs,sitename,U,expand)
      use parameters
      use linalg, only : rotate
      implicit none
      complex(8),allocatable,intent(inout)  :: Oloc(:,:)
      complex(8),allocatable,intent(in)     :: Oimp(:,:)
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:,:)
      logical,intent(in),optional           :: expand
      complex(8),allocatable                :: Otmp(:,:,:)
      integer                               :: isite,Nsite,shift
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      logical                               :: expand_
      !
      !
      write(*,"(A)") "--- imp2loc(O) ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Insertion of Operator for site: "//trim(sitename)
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "imp2loc(O): Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "imp2loc(O): Oimp not square."
      if(size(orbs).ne.size(Oimp,dim=1)) stop "imp2loc(O): can't fit the requested orbitals from Gimp."
      if(size(orbs).gt.size(Oloc,dim=1)) stop "imp2loc(O): number of requested orbitals greater than Gloc size."
      if(present(U))then
         write(*,"(A)") "The local orbital space will be rotated during extraction."
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if(size(U,dim=1).ne.size(orbs)) stop "Rotation matrix has the wrong dimension."
         if(size(U,dim=1).ne.3) write(*,"(A)") "Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
      endif
      expand_=.false.
      Nsite=1
      if(present(expand))then
         expand_ = expand
         Nsite = size(Oloc,dim=1)/size(orbs)
         if(present(U).and.(size(U,dim=3).ne.Nsite))stop "Number of rotation matrices and number of sites does not match."
         if(mod(size(Oloc,dim=1),size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of Gloc."
      endif
      !
      allocate(Otmp(size(Oloc,dim=1),size(Oloc,dim=1),Nsite));Otmp=czero
      !
      ! Rotating either one site or all of them depending on expand_
      do isite=1,Nsite
         !
         if(present(U))then
            Otmp(:,:,isite) = rotate(Oimp(:,:),U(:,:,isite))
         else
            Otmp(:,:,isite) = Oimp
         endif
         !
      enddo
      !
      do isite=1,Nsite
         !
         ! only two possible arrangements
         if(abs(orbs(2)-orbs(1)).eq.1)then
            shift = size(orbs)*(isite-1)
         elseif(abs(orbs(2)-orbs(1)).eq.Nsite)then
            shift = isite-1
         endif
         !
         do i_imp=1,size(orbs)
            do j_imp=1,size(orbs)
               !
               i_loc = orbs(i_imp) + shift
               j_loc = orbs(j_imp) + shift
               !
               Oloc(i_loc,j_loc) = Otmp(i_imp,j_imp,isite)
               !
            enddo
         enddo
         !
      enddo
      deallocate(Otmp)
      !
   end subroutine imp2loc_Matrix
   !
   subroutine imp2loc_Bosonic(Wloc,Wimp,orbs,sitename,expand)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: Wloc
      type(BosonicField),intent(in)         :: Wimp
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      logical,intent(in),optional           :: expand
      integer                               :: ip,isite,Nsite,shift
      integer                               :: Norb_imp,Norb_loc
      integer                               :: ib_imp,jb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: i_imp,j_imp,k_imp,l_imp
      logical                               :: expand_
      !
      !
      write(*,"(A)") "--- imp2loc(B) ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Insertion of Wimp of site: "//trim(sitename)
      if(.not.Wloc%status) stop "imp2loc(B): Wloc not properly initialized."
      if(.not.Wimp%status) stop "imp2loc(B): Wimp not properly initialized."
      if(Wloc%Nbp.eq.0) stop "imp2loc(B): Norb of Wloc not defined."
      if(Wimp%Nbp.eq.0) stop "imp2loc(B): Norb of Wimp not defined."
      if(Wloc%Npoints.eq.0) stop "imp2loc(B): Npoints of Wloc not defined."
      if(Wimp%Npoints.eq.0) stop "imp2loc(B): Npoints of Wimp not defined."
      if(Wloc%Beta.ne.Wimp%Beta) stop "imp2loc(B): Wimp and Wloc have different beta."
      if(Wloc%Npoints.ne.Wimp%Npoints) stop "imp2loc(B): Wimp and Wloc have different number of Matsubara points."
      if(Wimp%Nkpt.ne.0) stop "imp2loc(B): Wimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Wloc%screened_local)) stop "imp2loc(B): Wloc screened_local attribute not allocated."
      if(.not.allocated(Wloc%bare_local)) stop "imp2loc(B): Wloc bare_local attribute not allocated."
      if(.not.allocated(Wimp%screened_local)) stop "imp2loc(B): Wimp screened_local attribute not allocated."
      if(.not.allocated(Wimp%bare_local)) stop "imp2loc(B): Wimp bare_local attribute not allocated."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      !
      if(size(orbs).ne.Norb_imp) stop "imp2loc(B): can't fit the requested orbitals from Wimp."
      if(size(orbs).gt.Norb_loc) stop "imp2loc(B): number of requested orbitals greater than Wloc size."
      expand_=.false.
      Nsite=1
      if(present(expand))then
         expand_ = expand
         Nsite = Wloc%Nsite
         if(mod(Norb_loc,size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of Wloc."
      endif
      !
      do isite=1,Nsite
         !
         ! only two possible arrangements
         if(abs(orbs(2)-orbs(1)).eq.1)then
            shift = size(orbs)*(isite-1)
         elseif(abs(orbs(2)-orbs(1)).eq.Nsite)then
            shift = isite-1
         endif
         !
         do i_imp=1,Norb_imp
            do j_imp=1,Norb_imp
               do k_imp=1,Norb_imp
                  do l_imp=1,Norb_imp
                     !
                     ! bosonic indexes of the impurity
                     ib_imp = k_imp + Norb_imp*(i_imp-1)
                     jb_imp = l_imp + Norb_imp*(j_imp-1)
                     !
                     ! mapping
                     i_loc = orbs(i_imp) + shift
                     j_loc = orbs(j_imp) + shift
                     k_loc = orbs(k_imp) + shift
                     l_loc = orbs(l_imp) + shift
                     !
                     ! corresponding bosonic indexes on the lattice
                     ib_loc = k_loc + Norb_loc*(i_loc-1)
                     jb_loc = l_loc + Norb_loc*(j_loc-1)
                     !
                     Wloc%bare_local(ib_loc,jb_loc) = Wimp%bare_local(ib_imp,jb_imp)
                     do ip=1,Wimp%Npoints
                        Wloc%screened_local(ib_loc,jb_loc,ip) = Wimp%screened_local(ib_imp,jb_imp,ip)
                     enddo
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      enddo
      !
   end subroutine imp2loc_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !---------------------------------------------------------------------------!
   subroutine MergeSelfEnergy(SigmaGW,SigmaImp,coeff,orbs,SigmaGW_DC)
      use parameters
      use linalg, only : rotate
      implicit none
      type(FermionicField),intent(inout)    :: SigmaGW
      type(FermionicField),intent(in),optional :: SigmaGW_DC
      type(FermionicField),intent(in),target :: SigmaImp
      real(8),intent(in)                    :: coeff
      integer,allocatable,intent(in)        :: orbs(:,:)
      !
      real(8)                               :: Beta
      integer                               :: iw,ik,isite,iorb,jorb
      integer                               :: ispin,ispin_DMFT,shift
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      integer                               :: Nkpt,Norb,Nmats,Nsite
      logical                               :: localDC
      type(FermionicField),target           :: SigmaImpFull
      type(FermionicField),pointer          :: SigmaImpUsed
      !
      !
      write(*,"(A)") "--- MergeSelfEnergy ---"
      !
      !
      ! Check on the input Fields
      if(.not.SigmaGW%status) stop "SigmaGW not properly initialized."
      if(.not.SigmaGW_DC%status) stop "SigmaGW_DC not properly initialized."
      if(.not.SigmaImp%status) stop "SigmaImp not properly initialized."
      if(SigmaGW%Nkpt.eq.0) stop "SigmaGW k dependent attributes not properly initialized."
      if(present(SigmaGW_DC).and.(SigmaGW_DC%Nkpt.ne.0)) stop "SigmaGW_DC k dependent attributes are supposed to be unallocated."
      if(SigmaImp%Nkpt.ne.0) stop "SigmaImp k dependent attributes are supposed to be unallocated."
      !
      Norb = SigmaGW%Norb
      Nkpt = SigmaGW%Nkpt
      Beta = SigmaGW%Beta
      Nmats = SigmaGW%Npoints
      localDC = .true.
      !
      if(SigmaImp%Beta.ne.Beta) stop "SigmaImp has different Beta with respect to SigmaGW."
      if(SigmaImp%Npoints.ne.Nmats) stop "SigmaImp has different number of Matsubara points with respect to SigmaGW."
      if(present(SigmaGW_DC))then
         if(SigmaGW_DC%Norb.ne.Norb) stop "SigmaGW_DC has different orbital dimension with respect to SigmaGW."
         if(SigmaGW_DC%Beta.ne.Beta) stop "SigmaGW_DC has different Beta with respect to SigmaGW."
         if(SigmaGW_DC%Npoints.ne.Nmats) stop "SigmaGW_DC has different number of Matsubara points with respect to SigmaGW."
         localDC = .false.
      endif
      !
      if(AFM)then
         if(SigmaGW%Nsite.ne.2) stop "AFM is implemented only for a two site lattice."
         if(Norb/SigmaImp%Norb.ne.2) stop "SigmaGW and SigmaImp does not match the AFM condition."
         if(expand) stop "AFM condition and expansion to real space not yet implemented."
         Nsite = 2
      endif
      !
      if(expand)then
         if(size(U,dim=3).ne.SigmaGW%Nsite)stop "Number of rotation matrices and number of sites does not match."
         if(mod(Norb,size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of SigmaGW."
         if(AFM) stop "Expansion to real space and AFM condition not yet implemented."
         write(*,"(A)") "SigmaImp will be expanded to match the SigmaGW orbital space."
         call AllocateFermionicField(SigmaImpFull,Norb,Nmats,0)
         call imp2loc_Fermionic(SigmaImpFull,SigmaImp,orbs,U=U,expand=.true.)
         SigmaImpUsed => SigmaImpFull
         Nsite = SigmaGW%Nsite
      else
         if(size(orbs).ne.SigmaImp%Norb) stop "Can't fit the requested orbitals from SigmaImp."
         write(*,"(A,15I3)") "SigmaImp will be inserted into the SigmaGW orbital indexes: ",orbs
         SigmaImpUsed => SigmaImp
         Nsite = 1
      endif
      !
      !all sites if(expand.or.AFM) otherwise only one site and the orbitals within orbs
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nsite,Nmats,Nkpt,orbs,expand,AFM,coeff,SigmaGW,SigmaGW_DC,SigmaImpUsed,localDC),&
      !$OMP PRIVATE(isite,shift,ispin_DMFT,ispin,iorb,jorb,i_loc,j_loc,i_imp,j_imp,iw,ik)
      !$OMP DO
      do isite=1,Nsite
         !
         ! only two possible arrangements
         if(abs(orbs(2)-orbs(1)).eq.1)then
            shift = size(orbs)*(isite-1)
         elseif(abs(orbs(2)-orbs(1)).eq.Nsite)then
            shift = isite-1
         endif
         !
         do iorb=1,size(orbs)
            do jorb=1,size(orbs)
               !
               do ispin=1,Nspin
                  !
                  ispin_DMFT=ispin
                  !
                  if(expand)then
                     !SigmaImp has the same dimension of SigmaGW becasue it has been expanded
                     i_loc = iorb + shift
                     j_loc = jorb + shift
                     i_imp = i_loc
                     j_imp = j_loc
                  elseif(AFM)then
                     !SigmaImp has the dimension of orbs and I'm just flipping the spin
                     i_loc = iorb + shift
                     j_loc = jorb + shift
                     i_imp = iorb
                     j_imp = jorb
                     if(isite.eq.2) ispin_DMFT=int(Nspin/ispin)
                  else
                     !SigmaImp has the dimension of orbs
                     i_loc = orbs(iorb)
                     j_loc = orbs(jorb)
                     i_imp = iorb
                     j_imp = jorb
                  endif
                  !
                  do iw=1,Nmats
                     do ik=1,Nkpt
                        !
                        if(localDC)then
                           SigmaGW%wks(i_loc,j_loc,iw,ik,ispin) = SigmaGW%wks(i_loc,j_loc,iw,ik,ispin)              &
                                                                - coeff*SigmaGW%ws(i_loc,j_loc,iw,ispin)            &
                                                                + coeff*SigmaImpUsed%ws(i_imp,j_imp,iw,ispin_DMFT)
                        else
                           SigmaGW%wks(i_loc,j_loc,iw,ik,ispin) = SigmaGW%wks(i_loc,j_loc,iw,ik,ispin)              &
                                                                - coeff*SigmaGW_DC%ws(i_loc,j_loc,iw,ispin)         &
                                                                + coeff*SigmaImpUsed%ws(i_imp,j_imp,iw,ispin_DMFT)
                        endif
                        !
                     enddo
                  enddo
                  !
               enddo !ispin
               !
            enddo !jorb
         enddo !iorb
      enddo !isite
      !$OMP END DO
      !$OMP END PARALLEL
      !
      call FermionicKsum(SigmaGW)
      !
      if(expand)call DeallocateFermionicField(SigmaImpFull)
      if(associated(SigmaImpUsed))nullify(SigmaImpUsed)
      !
   end subroutine MergeSelfEnergy


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !---------------------------------------------------------------------------!
   subroutine MergePolarization(PiGW,PiImp,coeff,orbs,expand)
      use parameters
      use linalg, only : rotate
      implicit none
      type(BosonicField),intent(inout)      :: PiGW
      type(BosonicField),intent(in),target  :: PiImp
      real(8),intent(in)                    :: coeff
      integer,allocatable,intent(in)        :: orbs(:)
      logical,intent(in)                    :: expand
      !
      real(8)                               :: Beta
      integer                               :: iw,ik,isite,shift
      integer                               :: iorb,jorb,korb,lorb
      integer                               :: ib_imp,jb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: Nkpt,Nmats,Nsite
      integer                               :: Norb_imp,Norb_loc
      type(BosonicField),target             :: PiImpFull
      type(BosonicField),pointer            :: PiImpUsed
      !
      !
      write(*,"(A)") "--- MergePolarization ---"
      !
      !
      ! Check on the input Fields
      if(.not.PiGW%status) stop "PiGW not properly initialized."
      if(.not.PiImp%status) stop "PiImp not properly initialized."
      if(PiGW%Nkpt.eq.0) stop "PiGW k dependent attributes not properly initialized."
      if(PiImp%Nkpt.ne.0) stop "PiImp k dependent attributes are supposed to be unallocated."
      !
      Norb_loc = int(sqrt(dble(PiGW%Nbp)))
      Norb_imp = int(sqrt(dble(PiImp%Nbp)))
      Nkpt = PiGW%Nkpt
      Beta = PiGW%Beta
      Nmats = PiGW%Npoints
      !
      if(PiImp%Beta.ne.Beta) stop "PiImp has different Beta with respect to PiGW."
      if(PiImp%Npoints.ne.Nmats) stop "PiImp has different number of Matsubara points with respect to PiGW."
      !
      if(expand)then
         if(mod(Norb_loc,size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of PiGW."
         write(*,"(A)") "PiImp will be expanded to match the PiGW orbital space."
         call AllocateBosonicField(PiImpFull,Norb_loc,Nmats,0)
         call imp2loc_Bosonic(PiImpFull,PiImp,orbs,expand=.true.)
         PiImpUsed => PiImpFull
         Nsite = PiGW%Nsite
      else
         if(size(orbs).ne.Norb_imp) stop "Can't fit the requested orbitals from PiImp."
         write(*,"(A,15I3)") "PiImp will be inserted into the PiGW orbital indexes: ",orbs
         PiImpUsed => PiImp
         Nsite = 1
      endif
      !
      !all sites if(expand) otherwise only one site and the orbitals within orbs
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nsite,Nmats,Nkpt,Norb_loc,Norb_imp,orbs,expand,coeff,PiGW,PiImpUsed),&
      !$OMP PRIVATE(isite,shift,iorb,jorb,korb,lorb,iw,ik),&
      !$OMP PRIVATE(i_loc,j_loc,k_loc,l_loc,ib_imp,jb_imp,ib_loc,jb_loc)
      !$OMP DO
      do isite=1,Nsite
         !
         ! only two possible arrangements
         if(abs(orbs(2)-orbs(1)).eq.1)then
            shift = size(orbs)*(isite-1)
         elseif(abs(orbs(2)-orbs(1)).eq.Nsite)then
            shift = isite-1
         endif
         !
         do iorb=1,size(orbs)
            do jorb=1,size(orbs)
               do korb=1,size(orbs)
                  do lorb=1,size(orbs)
                     !
                     if(expand)then
                        !mapping
                        i_loc = iorb + shift
                        j_loc = jorb + shift
                        k_loc = korb + shift
                        l_loc = lorb + shift
                        !bosonic indexes on the lattice
                        ib_loc = k_loc + Norb_loc*(i_loc-1)
                        jb_loc = l_loc + Norb_loc*(j_loc-1)
                        !PiImp has the same dimension of PiGW becasue it has been expanded
                        ib_imp = ib_loc
                        jb_imp = jb_loc
                     else
                        !mapping
                        i_loc = orbs(iorb)
                        j_loc = orbs(jorb)
                        k_loc = orbs(korb)
                        l_loc = orbs(lorb)
                        !bosonic indexes on the lattice
                        ib_loc = k_loc + Norb_loc*(i_loc-1)
                        jb_loc = l_loc + Norb_loc*(j_loc-1)
                        !PiImp has the dimension of orbs
                        ib_imp = korb + Norb_imp*(iorb-1)
                        jb_imp = lorb + Norb_imp*(jorb-1)
                     endif
                     !
                     if((i_loc.eq.k_loc).and.(l_loc.eq.j_loc))then
                        !
                        do iw=1,Nmats
                           do ik=1,Nkpt
                              !
                              PiGW%screened(ib_loc,jb_loc,iw,ik) = PiGW%screened(ib_loc,jb_loc,iw,ik)                &
                                                                 - coeff*PiGW%screened_local(ib_loc,jb_loc,iw)       &
                                                                 + coeff*PiImpUsed%screened_local(ib_imp,jb_imp,iw)
                              !
                           enddo
                        enddo
                        !
                     endif
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      enddo !isite
      !$OMP END DO
      !$OMP END PARALLEL
      !
      call BosonicKsum(PiGW)
      !
      if(expand)call DeallocateBosonicField(PiImpFull)
      if(associated(PiImpUsed))nullify(PiImpUsed)
      !
   end subroutine MergePolarization


   !---------------------------------------------------------------------------!
   !PURPOSE: Clear the internal attributes of a Fermionic/Bosonic field
   !---------------------------------------------------------------------------!
   subroutine clear_attributes_Fermion(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      if(allocated(G%N_s))G%N_s=czero
      if(allocated(G%ws))G%ws=czero
      if(allocated(G%N_ks))G%N_ks=czero
      if(allocated(G%wks))G%wks=czero
   end subroutine clear_attributes_Fermion
   !
   subroutine clear_attributes_Boson(W)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      if(allocated(W%bare_local))W%bare_local=czero
      if(allocated(W%screened_local))W%screened_local=czero
      if(allocated(W%bare))W%bare=czero
      if(allocated(W%screened))W%screened=czero
   end subroutine clear_attributes_Boson



end module utils_fields
