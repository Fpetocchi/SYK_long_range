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
      module procedure Gloc2Gimp
      module procedure Oloc2Oimp
      module procedure Wloc2Wimp
   end interface loc2imp

   interface imp2loc
      module procedure Gimp2Gloc
      module procedure Oimp2Oloc
      module procedure Wimp2Wloc
   end interface imp2loc

   interface MergeField
      module procedure MergeSelfEnergy
      !module procedure MergePolarization
   end interface MergeField

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
   public :: MergeField
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
      integer,intent(in)                    :: Norb,Npoints,Nkpt
      integer,intent(in),optional           :: Nsite
      character(len=*),intent(in),optional  :: name
      !
      if(present(name)) write(*,"(A)") "Allocation of "//trim(name)
      if(G%status) stop "AllocateFermionicField: container already allocated."
      if(Norb.eq.0) stop "AllocateFermionicField: Norb not defined."
      if(Npoints.eq.0) write(*,"(A)") "AllocateFermionicField: frequency dependent attributes are not going to be allocated."
      if(Nkpt.eq.0) write(*,"(A)") "AllocateFermionicField: K-dependent attributes are not going to be allocated."
      !
      if(allocated(G%N_s))deallocate(G%N_s)
      if(allocated(G%N_ks))deallocate(G%N_ks)
      if(allocated(G%ws))deallocate(G%ws)
      if(allocated(G%wks))deallocate(G%wks)
      !
      allocate(G%N_s(Norb,Norb,Nspin));G%N_s=czero
      !
      if((Npoints.eq.0).and.(Nkpt.ne.0))then
         allocate(G%N_ks(Norb,Norb,Nkpt,Nspin));G%N_ks=czero
      elseif((Npoints.ne.0).and.(Nkpt.eq.0))then
         allocate(G%ws(Norb,Norb,Npoints,Nspin));G%ws=czero
      else
         allocate(G%N_ks(Norb,Norb,Nkpt,Nspin));G%N_ks=czero
         allocate(G%ws(Norb,Norb,Npoints,Nspin));G%ws=czero
         allocate(G%wks(Norb,Norb,Npoints,Nkpt,Nspin));G%wks=czero
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
   subroutine AllocateBosonicField(W,Nbp,Npoints,Nkpt,Nsite,name,no_bare)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      integer,intent(in)                    :: Nbp,Npoints,Nkpt
      integer,intent(in),optional           :: Nsite
      character(len=*),intent(in),optional  :: name
      logical,intent(in),optional           :: no_bare
      logical                               :: no_bare_
      !
      if(present(name)) write(*,"(A)") "Allocation of "//trim(name)
      if(W%status) stop "AllocateBosonicField: container already allocated."
      if(Nbp.eq.0) stop "AllocateBosonicField: Nbp not defined."
      if(Npoints.eq.0) stop "AllocateBosonicField: Npoints not defined."
      if(Nkpt.eq.0) write(*,"(A)") "AllocateBosonicField: K-dependent attributes are not going to be allocated."
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
      if(Nkpt.ne.0)then
         allocate(W%screened(Nbp,Nbp,Npoints,Nkpt));W%screened=czero
      endif
      !
      if(.not.no_bare_)then
         allocate(W%bare_local(Nbp,Nbp));W%bare_local=czero
         if(Nkpt.ne.0)then
            allocate(W%bare(Nbp,Nbp,Nkpt));W%bare=czero
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
   subroutine Gloc2Gimp(Gimp,Gloc,orbs,sitename,U)
      use parameters
      use linalg, only : rotate
      implicit none
      type(FermionicField),intent(inout)    :: Gimp
      type(FermionicField),intent(in)       :: Gloc
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:)
      integer                               :: ip,ispin
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      write(*,"(A)") "--- Gloc2Gimp ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Extraction of Gloc for site: "//trim(sitename)
      if(.not.Gloc%status) stop "Gloc2Gimp: Gloc not properly initialized."
      if(.not.Gimp%status) stop "Gloc2Gimp: Gimp not properly initialized."
      if(Gloc%Norb.eq.0) stop "Gloc2Gimp: Norb of Gloc not defined."
      if(Gimp%Norb.eq.0) stop "Gloc2Gimp: Norb of Gimp not defined."
      if(Gloc%Npoints.eq.0) stop "Gloc2Gimp: Npoints of Gloc not defined."
      if(Gimp%Npoints.eq.0) stop "Gloc2Gimp: Npoints of Gimp not defined."
      if(Gimp%Beta.ne.Gloc%Beta) stop "Gimp2Gimp: Gimp and Gloc have different beta."
      if(Gimp%Npoints.ne.Gloc%Npoints) stop "Gimp2Gimp: Gimp and Gloc have different number of Matsubara points."
      if(Gimp%Nkpt.ne.0) stop "Gloc2Gimp: Gimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Gloc%ws)) stop "Gloc2Gimp: Gloc local projection not allocated."
      if(.not.allocated(Gimp%ws)) stop "Gloc2Gimp: Gimp local projection not allocated."
      if(size(orbs).ne.Gimp%Norb) stop "Gloc2Gimp: can't fit the requested orbitals inside Gimp."
      if(size(orbs).gt.Gloc%Norb) stop "Gloc2Gimp: number of requested orbitals greater than Gloc size."
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
   end subroutine Gloc2Gimp
   !
   subroutine Oloc2Oimp(Oimp,Oloc,orbs,sitename,U)
      use parameters
      use linalg, only : rotate
      implicit none
      complex(8),allocatable,intent(inout)  :: Oimp(:,:)
      complex(8),allocatable,intent(in)     :: Oloc(:,:)
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:)
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      write(*,"(A)") "--- Oloc2Oimp ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Extraction of Operator for site: "//trim(sitename)
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "Oloc2Oimp: Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "Oloc2Oimp: Oimp not square."
      if(size(orbs).ne.size(Oimp,dim=1)) stop "Oloc2Oimp: can't fit the requested orbitals inside Oimp."
      if(size(orbs).gt.size(Oloc,dim=1)) stop "Oloc2Oimp: number of requested orbitals greater than Oloc size."
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
   end subroutine Oloc2Oimp
   !
   subroutine Wloc2Wimp(Wimp,Wloc,orbs,sitename)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: Wimp
      type(BosonicField),intent(in)         :: Wloc
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      integer                               :: ip
      integer                               :: Norb_imp,Norb_loc
      integer                               :: ib_imp,jb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: i_imp,j_imp,k_imp,l_imp
      !
      !
      write(*,"(A)") "--- Wloc2Wimp ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Extraction of Wloc for site: "//trim(sitename)
      if(.not.Wloc%status) stop "Wloc2Wimp: Wloc not properly initialized."
      if(.not.Wimp%status) stop "Wloc2Wimp: Wimp not properly initialized."
      if(Wloc%Nbp.eq.0) stop "Wloc2Wimp: Norb of Wloc not defined."
      if(Wimp%Nbp.eq.0) stop "Wloc2Wimp: Norb of Wimp not defined."
      if(Wloc%Npoints.eq.0) stop "Wloc2Wimp: Npoints of Wloc not defined."
      if(Wimp%Npoints.eq.0) stop "Wloc2Wimp: Npoints of Wimp not defined."
      if(Wimp%Beta.ne.Wloc%Beta) stop "Wimp2Wimp: Wimp and Wloc have different beta."
      if(Wimp%Npoints.ne.Wloc%Npoints) stop "Wimp2Wimp: Wimp and Wloc have different number of Matsubara points."
      if(Wimp%Nkpt.ne.0) stop "Wloc2Wimp: Wimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Wloc%screened_local)) stop "Wloc2Wimp: Wloc screened_local attribute not allocated."
      if(.not.allocated(Wloc%bare_local)) stop "Wloc2Wimp: Wloc bare_local attribute not allocated."
      if(.not.allocated(Wimp%screened_local)) stop "Wloc2Wimp: Wimp screened_local attribute not allocated."
      if(.not.allocated(Wimp%bare_local)) stop "Wloc2Wimp: Wimp bare_local attribute not allocated."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      !
      if(size(orbs).ne.Norb_imp) stop "Wloc2Wimp: can't fit the requested orbitals inside Wimp."
      if(size(orbs).gt.Norb_loc) stop "Wloc2Wimp: number of requested orbitals greater than Wloc size."
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
   end subroutine Wloc2Wimp


   !---------------------------------------------------------------------------!
   !PURPOSE: Expand the the subset of orbital indexes of corresponging to a
   !         given site to the full lattice quantity
   !---------------------------------------------------------------------------!
   subroutine Gimp2Gloc(Gloc,Gimp,orbs,sitename,U,expand)
      use parameters
      use linalg, only : rotate
      implicit none
      type(FermionicField),intent(inout)    :: Gloc
      type(FermionicField),intent(in)       :: Gimp
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:,:)
      logical,intent(in),optional           :: expand
      complex(8),allocatable                :: Gtmp(:,:,:,:,:)
      integer                               :: ip,ispin,isite,Nsite,shift
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      logical                               :: expand_
      !
      !
      write(*,"(A)") "--- Gimp2Gloc ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Insertion of Gimp of site: "//trim(sitename)
      if(.not.Gloc%status) stop "Gimp2Gloc: Gloc not properly initialized."
      if(.not.Gimp%status) stop "Gimp2Gloc: Gimp not properly initialized."
      if(Gloc%Norb.eq.0) stop "Gimp2Gloc: Norb of Gloc not defined."
      if(Gimp%Norb.eq.0) stop "Gimp2Gloc: Norb of Gimp not defined."
      if(Gloc%Npoints.eq.0) stop "Gimp2Gloc: Npoints of Gloc not defined."
      if(Gimp%Npoints.eq.0) stop "Gimp2Gloc: Npoints of Gimp not defined."
      if(Gloc%Beta.ne.Gimp%Beta) stop "Gloc2Gimp: Gimp and Gloc have different beta."
      if(Gloc%Npoints.ne.Gimp%Npoints) stop "Gloc2Gimp: Gimp and Gloc have different number of Matsubara points."
      if(Gimp%Nkpt.ne.0) stop "Gimp2Gloc: Gimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Gloc%ws)) stop "Gimp2Gloc: Gloc local projection not allocated."
      if(.not.allocated(Gimp%ws)) stop "Gimp2Gloc: Gimp local projection not allocated."
      if(size(orbs).ne.Gimp%Norb) stop "Gimp2Gloc: can't fit the requested orbitals from Gimp."
      if(size(orbs).gt.Gloc%Norb) stop "Gimp2Gloc: number of requested orbitals greater than Gloc size."
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
         Nsite = Gloc%Nsite
         if(present(U).and.(size(U,dim=3).ne.Nsite))stop "Number of rotation matrices and number of sites does not match."
         if(size(orbs).ne.(Gloc%Norb/Nsite))stop "Number of requested orbitals is not a commensurate subset of Gloc."
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
   end subroutine Gimp2Gloc
   !
   subroutine Oimp2Oloc(Oloc,Oimp,orbs,sitename,U,expand)
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
      write(*,"(A)") "--- Oimp2Oloc ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Insertion of Operator for site: "//trim(sitename)
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "Oimp2Oloc: Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "Oimp2Oloc: Oimp not square."
      if(size(orbs).ne.size(Oimp,dim=1)) stop "Oimp2Oloc: can't fit the requested orbitals from Gimp."
      if(size(orbs).gt.size(Oloc,dim=1)) stop "Oimp2Oloc: number of requested orbitals greater than Gloc size."
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
         if(size(orbs).ne.(size(Oloc,dim=1)/Nsite))stop "Number of requested orbitals is not a commensurate subset of Gloc."
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
   end subroutine Oimp2Oloc
   !
   subroutine Wimp2Wloc(Wloc,Wimp,orbs,sitename,expand)
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
      write(*,"(A)") "--- Wimp2Wloc ---"
      !
      !
      if(present(sitename)) write(*,"(A)") "Insertion of Wimp of site: "//trim(sitename)
      if(.not.Wloc%status) stop "Wloc2Wimp: Wloc not properly initialized."
      if(.not.Wimp%status) stop "Wloc2Wimp: Wimp not properly initialized."
      if(Wloc%Nbp.eq.0) stop "Wloc2Wimp: Norb of Wloc not defined."
      if(Wimp%Nbp.eq.0) stop "Wloc2Wimp: Norb of Wimp not defined."
      if(Wloc%Npoints.eq.0) stop "Wloc2Wimp: Npoints of Wloc not defined."
      if(Wimp%Npoints.eq.0) stop "Wloc2Wimp: Npoints of Wimp not defined."
      if(Wloc%Beta.ne.Wimp%Beta) stop "Wloc2Wimp: Wimp and Wloc have different beta."
      if(Wloc%Npoints.ne.Wimp%Npoints) stop "Wloc2Wimp: Wimp and Wloc have different number of Matsubara points."
      if(Wimp%Nkpt.ne.0) stop "Wloc2Wimp: Wimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Wloc%screened_local)) stop "Wloc2Wimp: Wloc screened_local attribute not allocated."
      if(.not.allocated(Wloc%bare_local)) stop "Wloc2Wimp: Wloc bare_local attribute not allocated."
      if(.not.allocated(Wimp%screened_local)) stop "Wloc2Wimp: Wimp screened_local attribute not allocated."
      if(.not.allocated(Wimp%bare_local)) stop "Wloc2Wimp: Wimp bare_local attribute not allocated."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      !
      if(size(orbs).ne.Norb_imp) stop "Wloc2Wimp: can't fit the requested orbitals from Wimp."
      if(size(orbs).gt.Norb_loc) stop "Wloc2Wimp: number of requested orbitals greater than Wloc size."
      expand_=.false.
      Nsite=1
      if(present(expand))then
         expand_ = expand
         Nsite = Wloc%Nsite
         if(size(orbs).ne.Norb_loc/Nsite)stop "Number of requested orbitals is not a commensurate subset of Wloc."
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
   end subroutine Wimp2Wloc


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !---------------------------------------------------------------------------!
   subroutine MergeSelfEnergy(SigmaGW,SigmaGW_DC,SigmaImp,coeff,AFM,orbs,U)
      use parameters
      use linalg, only : rotate
      implicit none
      type(FermionicField),intent(inout)    :: SigmaGW
      type(FermionicField),intent(in)       :: SigmaGW_DC
      type(FermionicField),intent(in),target :: SigmaImp
      real(8),intent(in)                    :: coeff
      logical,intent(in)                    :: AFM
      complex(8),allocatable,optional       :: U(:,:,:)
      integer,allocatable,optional          :: orbs(:)
      !
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Norb,Nmats
      type(FermionicField),target           :: SigmaImpFull
      type(FermionicField),pointer          :: SigmaImpUsed
      !
      !
      write(*,"(A)") "--- MergeSelfEnergy ---"
      !
      !
      ! Check on the input Fields
      if(.not.SigmaGW%status) stop "Smats_C not properly initialized."
      if(.not.SigmaGW_DC%status) stop "Smats_X not properly initialized."
      if(.not.SigmaImp%status) stop "Gmats not properly initialized."
      if(SigmaGW%Nkpt.eq.0) stop "SigmaGW_DC k dependent attributes not properly initialized."
      if(SigmaGW_DC%Nkpt.ne.0) stop "SigmaGW_DC k dependent attributes are supposed to be unallocated."
      if(SigmaImp%Nkpt.ne.0) stop "SigmaImp k dependent attributes are supposed to be unallocated."
      !
      Norb = SigmaGW%Norb
      Nkpt = SigmaGW%Nkpt
      Beta = SigmaGW%Beta
      Nmats = SigmaGW%Npoints
      !
      if(SigmaGW_DC%Norb.ne.Norb) stop "SigmaGW_DC has different orbital dimension with respect to SigmaGW."
      if(all([SigmaGW_DC%Beta-Beta,SigmaImp%Beta-Beta].ne.[0d0,0d0])) stop "Either Smats_X, Gmats or Wmats have different Beta with respect to SigmaGW."
      if(all([SigmaGW_DC%Npoints-Nmats,SigmaImp%Npoints-Nmats].ne.[0,0])) stop "Either Smats_C, Gmats or Wmats have different number of Matsubara points with respect to SigmaGW."
      if(present(orbs))then
         write(*,"(A)") "Warning: SigmaImp has different orbital dimension with respect to SigmaGW. Trying to expand it."
         call AllocateFermionicField(SigmaImpFull,Norb,Nmats,0)
         if(present(U))then
            call Gimp2Gloc(SigmaImpFull,SigmaImp,orbs,U=U,expand=.true.)
         else
            call Gimp2Gloc(SigmaImpFull,SigmaImp,orbs,expand=.true.)
         endif
         SigmaImpUsed => SigmaImpFull
      else
         if(SigmaImp%Norb.ne.Norb) stop "SigmaImp has different orbital dimension with respect to SigmaGW."
         SigmaImpUsed => SigmaImp
      endif
      !



      !
      if(present(orbs))call DeallocateFermionicField(SigmaImpFull)
      if(associated(SigmaImpUsed))nullify(SigmaImpUsed)
      !
   end subroutine MergeSelfEnergy


























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
