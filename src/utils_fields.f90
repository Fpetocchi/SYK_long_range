module utils_fields

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   interface clear_attributes
      module procedure clear_attributes_Fermion
      module procedure clear_attributes_Boson
   end interface clear_attributes

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: FermionicKsum
   public :: BosonicKsum
   public :: selfAllocateLattice
   public :: selfDeallocateLattice
   public :: selfAllocateFermionicField
   public :: selfDeallocateFermionicField
   public :: selfAllocateBosonicField
   public :: selfDeallocateBosonicField
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
      if(.not.G%status)stop "FermionicKsum. Field not allocated."
      !
      G%w=czero
      do ispin=1,Nspin
         do ik=1,G%Nkpt
            do in=1,G%Npoints
               G%w(:,:,in,ispin) = G%w(:,:,in,ispin) + G%wk(:,:,in,ik,ispin)/G%Nkpt
            enddo
         enddo
      enddo
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
      if(.not.W%status)stop "BosonicKsum. Field not allocated."
      !
      if(allocated(W%bare_local))W%bare_local=czero
      W%screened_local=czero
      do ik=1,W%Nkpt
         if(allocated(W%bare_local))W%bare_local = W%bare_local + W%bare(:,:,ik)/W%Nkpt
         do in=1,W%Npoints
            W%screened_local(:,:,in) = W%screened_local(:,:,in) + W%screened(:,:,in,ik)/W%Nkpt
         enddo
      enddo
      !
   end subroutine BosonicKsum


   !---------------------------------------------------------------------------!
   !PURPOSE: Allocate/deallocate Lattice attributes in a consistent way
   !---------------------------------------------------------------------------!
   subroutine selfAllocateLattice(lttc)
      use parameters
      implicit none
      type(Lattice),intent(inout)             :: lttc
      !
      if(.not.lttc%status) stop "selfAllocateLattice: Field not properly initialized."
      if(lttc%Norb.eq.0) stop "selfAllocateLattice: Norb not defined."
      if(lttc%Nkpt.eq.0) stop "selfAllocateLattice: Nkpt not defined."
      !
      if(allocated(lttc%kpt))deallocate(lttc%kpt)
      allocate(lttc%kpt(2,lttc%Nkpt));lttc%kpt=0d0
      !
      if(allocated(lttc%Hk))deallocate(lttc%Hk)
      allocate(lttc%Hk(lttc%Norb,lttc%Norb,lttc%Nkpt));lttc%Hk=czero
      !
      if(allocated(lttc%Zk))deallocate(lttc%Zk)
      allocate(lttc%Zk(lttc%Norb,lttc%Norb,lttc%Nkpt));lttc%Zk=czero
      !
      if(allocated(lttc%Ek))deallocate(lttc%Ek)
      allocate(lttc%Ek(lttc%Norb,lttc%Nkpt));lttc%Ek=0d0
      !
      if(allocated(lttc%kptPos))deallocate(lttc%kptPos)
      allocate(lttc%kptPos(lttc%Nkpt));lttc%kptPos=0
      !
      if(allocated(lttc%kptsum))deallocate(lttc%kptsum)
      allocate(lttc%kptsum(lttc%Nkpt,lttc%Nkpt));lttc%kptsum=0
      !
      if(allocated(lttc%kptdif))deallocate(lttc%kptdif)
      allocate(lttc%kptdif(lttc%Nkpt,lttc%Nkpt));lttc%kptdif=0
      !
      if(allocated(lttc%small_ik))deallocate(lttc%small_ik)
      allocate(lttc%small_ik(12,2));lttc%small_ik=0
      !
   end subroutine selfAllocateLattice
   !
   subroutine selfDeallocateLattice(lttc)
      use parameters
      implicit none
      type(Lattice),intent(inout)           :: lttc
      !
      if(.not.lttc%status) stop "selfDeallocateLattice: Field not properly initialized."
      if(lttc%Norb.eq.0) stop "selfDeallocateLattice: Norb not defined."
      if(lttc%Nkpt.eq.0) stop "selfDeallocateLattice: Nkpt not defined."
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
   end subroutine selfDeallocateLattice


   !---------------------------------------------------------------------------!
   !PURPOSE: Allocate/deallocate Fermionic attributes in a consistent way
   !---------------------------------------------------------------------------!
   subroutine selfAllocateFermionicField(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      !
      if(.not.G%status) stop "selfAllocateFermionicField: Field not properly initialized."
      if(G%Norb.eq.0) stop "selfAllocateFermionicField: Norb not defined."
      if(G%Npoints.eq.0) stop "selfAllocateFermionicField: Npoints not defined."
      !
      if(allocated(G%w))deallocate(G%w)
      allocate(G%w(G%Norb,G%Norb,G%Npoints,Nspin));G%w=czero
      !
      if(G%Nkpt.ne.0)then
         if(allocated(G%wk))deallocate(G%wk)
         allocate(G%wk(G%Norb,G%Norb,G%Npoints,G%Nkpt,Nspin));G%wk=czero
      endif
      !
   end subroutine selfAllocateFermionicField
   !
   subroutine selfDeallocateFermionicField(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      !
      if(.not.G%status) stop "selfDeallocateFermionicField: Field not properly initialized."
      if(G%Norb.eq.0) stop "selfDeallocateFermionicField: Norb not defined."
      if(G%Npoints.eq.0) stop "selfDeallocateFermionicField: Npoints not defined."
      !
      if(allocated(G%w))deallocate(G%w)
      if(allocated(G%wk))deallocate(G%wk)
      G%Norb=0
      G%Npoints=0
      G%Nkpt=0
      G%Nsite=1
      G%Beta=0d0
      G%mu=0d0
      G%status=.false.
      !
   end subroutine selfDeallocateFermionicField


   !---------------------------------------------------------------------------!
   !PURPOSE: Allocate/deallocate Bosonic attributes in a consistent way
   !---------------------------------------------------------------------------!
   subroutine selfAllocateBosonicField(W,zerobare)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      logical,intent(in),optional           :: zerobare
      logical                               :: zerobare_
      !
      if(.not.W%status) stop "selfAllocateBosonicField: Field not properly initialized."
      if(W%Nbp.eq.0) stop "selfAllocateBosonicField: Nbp not defined."
      if(W%Npoints.eq.0) stop "selfAllocateBosonicField: Npoints not defined."
      zerobare_=.false.
      if(present(zerobare))zerobare_=zerobare
      !
      if(allocated(W%bare_local))deallocate(W%bare_local)
      if(.not.zerobare_)then
         allocate(W%bare_local(W%Nbp,W%Nbp))
         W%bare_local=czero
      endif
      !
      if(allocated(W%screened_local))deallocate(W%screened_local)
      allocate(W%screened_local(W%Nbp,W%Nbp,W%Npoints));W%screened_local=czero
      !
      if(W%Nkpt.ne.0)then
         if(allocated(W%bare))deallocate(W%bare)
         if(.not.zerobare_)then
            allocate(W%bare(W%Nbp,W%Nbp,W%Nkpt))
            W%bare=czero
         endif
         if(allocated(W%screened))deallocate(W%screened)
         allocate(W%screened(W%Nbp,W%Nbp,W%Npoints,W%Nkpt));W%screened=czero
      endif
      !
   end subroutine selfAllocateBosonicField
   !
   subroutine selfDeallocateBosonicField(W)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      !
      if(.not.W%status) stop "selfDeallocateBosonicField: Field not properly initialized."
      if(W%Nbp.eq.0) stop "selfDeallocateBosonicField: Nbp not defined."
      if(W%Npoints.eq.0) stop "selfDeallocateBosonicField: Npoints not defined."
      !
      if(allocated(W%bare_local))deallocate(W%bare_local)
      if(allocated(W%screened_local))deallocate(W%screened_local)
      if(allocated(W%bare))deallocate(W%bare)
      if(allocated(W%screened))deallocate(W%screened)
      W%Nbp=0
      W%Npoints=0
      W%Nkpt=0
      W%Nsite=1
      W%iq_gamma=-1
      W%Beta=0d0
      W%status=.false.
      !
   end subroutine selfDeallocateBosonicField


   !---------------------------------------------------------------------------!
   !PURPOSE: Clear the internal attributes of a Fermionic/Bosonic field
   !---------------------------------------------------------------------------!
   subroutine clear_attributes_Fermion(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      if(allocated(G%w))G%w=czero
      if(allocated(G%wk))G%wk=czero
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
