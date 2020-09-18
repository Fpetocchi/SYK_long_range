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
      G%s=czero
      do ispin=1,Nspin
         do ik=1,G%Nkpt
            G%s(:,:,ispin) = G%s(:,:,ispin) + G%ks(:,:,ik,ispin)/G%Nkpt
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
   subroutine AllocateFermionicField(G,Norb,Npoints,Nkpt,name)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      integer,intent(in)                    :: Norb,Npoints,Nkpt
      character(len=*),intent(in),optional  :: name
      !
      if(present(name)) write(*,"(A)") "Allocation of "//trim(name)
      if(G%status) stop "AllocateFermionicField: container already allocated."
      if(Norb.eq.0) stop "AllocateFermionicField: Norb not defined."
      if(Npoints.eq.0) write(*,"(A)") "AllocateFermionicField: frequency dependent attributes are not going to be allocated."
      if(Nkpt.eq.0) write(*,"(A)") "AllocateFermionicField: K-dependent attributes are not going to be allocated."
      !
      if(allocated(G%s))deallocate(G%s)
      if(allocated(G%ks))deallocate(G%ks)
      if(allocated(G%ws))deallocate(G%ws)
      if(allocated(G%wks))deallocate(G%wks)
      !
      allocate(G%s(Norb,Norb,Nspin));G%s=czero
      !
      if((Npoints.eq.0).and.(Nkpt.ne.0))then
         allocate(G%ks(Norb,Norb,Nkpt,Nspin));G%ks=czero
      elseif((Npoints.ne.0).and.(Nkpt.eq.0))then
         allocate(G%ws(Norb,Norb,Npoints,Nspin));G%ws=czero
      else
         allocate(G%ks(Norb,Norb,Nkpt,Nspin));G%ks=czero
         allocate(G%ws(Norb,Norb,Npoints,Nspin));G%ws=czero
         allocate(G%wks(Norb,Norb,Npoints,Nkpt,Nspin));G%wks=czero
      endif
      !
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
      if(allocated(G%s))deallocate(G%s)
      if(allocated(G%ks))deallocate(G%ks)
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
   subroutine AllocateBosonicField(W,Nbp,Npoints,Nkpt,name,no_bare)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      integer,intent(in)                    :: Nbp,Npoints,Nkpt
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
      W%Nsite=1
      W%iq_gamma=-1
      W%Beta=0d0
      W%status=.false.
      !
   end subroutine DeallocateBosonicField


   !---------------------------------------------------------------------------!
   !PURPOSE: Clear the internal attributes of a Fermionic/Bosonic field
   !---------------------------------------------------------------------------!
   subroutine clear_attributes_Fermion(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      if(allocated(G%s))G%s=czero
      if(allocated(G%ws))G%ws=czero
      if(allocated(G%ks))G%ks=czero
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
