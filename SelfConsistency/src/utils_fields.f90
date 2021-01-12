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
      module procedure loc2imp_Matrix
      module procedure loc2imp_Fermionic
      module procedure loc2imp_Bosonic
   end interface loc2imp

   interface imp2loc
      module procedure imp2loc_Matrix
      module procedure imp2loc_Fermionic
      module procedure imp2loc_Bosonic
   end interface imp2loc

   interface symmetrize
      module procedure symmetrize_Matrix
      module procedure symmetrize_Fermionic
      module procedure symmetrize_Bosonic
   end interface symmetrize

   interface MergeFields
      module procedure MergeSelfEnergy
      module procedure MergePolarization
   end interface MergeFields

   interface isReal
      module procedure isReal_Matrix
      module procedure isReal_Fermionic
      module procedure isReal_Bosonic
   end interface isReal

   interface duplicate
      module procedure duplicate_Fermionic
      module procedure duplicate_Bosonic
   end interface duplicate

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
   public :: clear_attributes
   public :: duplicate
   public :: isReal
   public :: loc2imp
   public :: imp2loc
   public :: symmetrize
   public :: join_SigmaCX
   public :: MergeFields

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Fill the local attributes of a Fermionic Field
   !TEST ON: 16-10-2020
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
      G%N_s = czero
      do ispin=1,Nspin
         do ik=1,G%Nkpt
            G%N_s(:,:,ispin) = G%N_s(:,:,ispin) + G%N_ks(:,:,ik,ispin)/G%Nkpt
         enddo
      enddo
      !
      if(G%Npoints.ne.0)then
         G%ws = czero
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
   !TEST ON: 21-10-2020
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
      W%screened_local = czero
      do ik=1,W%Nkpt
         do in=1,W%Npoints
            W%screened_local(:,:,in) = W%screened_local(:,:,in) + W%screened(:,:,in,ik)/W%Nkpt
         enddo
      enddo
      !
      if(allocated(W%bare_local).and.allocated(W%bare))then
         W%bare_local = czero
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
   !TEST ON: 14-10-2020(both)
   !---------------------------------------------------------------------------!
   subroutine AllocateLattice(lttc,Norb,Nkpt,name)
      use parameters
      implicit none
      type(Lattice),intent(inout)           :: lttc
      integer,intent(in)                    :: Norb,Nkpt
      character(len=*),intent(in),optional  :: name
      !
      if(present(name)) write(*,"(A)") "     Allocation of "//trim(name)
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
      if(present(name)) write(*,"(A)") "     Deallocation of "//trim(name)
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
   !TEST ON: 16-10-2020(both)
   !---------------------------------------------------------------------------!
   subroutine AllocateFermionicField(G,Norb,Npoints,Nkpt,Nsite,name,Beta,mu)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      integer,intent(in)                    :: Norb,Npoints
      integer,intent(in),optional           :: Nkpt
      integer,intent(in),optional           :: Nsite
      character(len=*),intent(in),optional  :: name
      real(8),intent(in),optional           :: Beta
      real(8),intent(in),optional           :: mu
      integer                               :: Nkpt_
      !
      Nkpt_=0
      if(present(Nkpt))Nkpt_=Nkpt
      if(present(name).and.verbose) write(*,"(A)") "     Allocation of "//trim(name)
      if(G%status) stop "AllocateFermionicField: container already allocated."
      if(Norb.eq.0) stop "AllocateFermionicField: Norb not defined."
      if(Npoints.eq.0.and.verbose) write(*,"(A)") "     AllocateFermionicField: frequency dependent attributes are not going to be allocated."
      if(Nkpt_.eq.0.and.verbose) write(*,"(A)") "     AllocateFermionicField: K-dependent attributes are not going to be allocated."
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
      G%Norb=Norb
      G%Npoints=Npoints
      G%Nkpt=Nkpt_
      if(present(mu))G%mu=mu
      if(present(Beta))G%Beta=Beta
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
      if(present(name).and.verbose) write(*,"(A)") "     Deallocation of "//trim(name)
      !if(.not.G%status) stop "DeallocateFermionicField: container is unallocated."
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
   !TEST ON: 21-10-2020(both)
   !---------------------------------------------------------------------------!
   subroutine AllocateBosonicField(W,Norb,Npoints,iq_gamma,Nkpt,Nsite,name,no_bare,Beta)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      integer,intent(in)                    :: Norb,Npoints,iq_gamma
      integer,intent(in),optional           :: Nkpt
      integer,intent(in),optional           :: Nsite
      character(len=*),intent(in),optional  :: name
      logical,intent(in),optional           :: no_bare
      real(8),intent(in),optional           :: Beta
      integer                               :: Nbp,Nkpt_
      logical                               :: no_bare_
      !
      Nbp=Norb**2
      Nkpt_=0
      if(present(Nkpt))Nkpt_=Nkpt
      if(present(name).and.verbose) write(*,"(A)") "     Allocation of "//trim(name)
      if(W%status) stop "AllocateBosonicField: container already allocated."
      if(Nbp.eq.0) stop "AllocateBosonicField: Nbp not defined."
      if(Npoints.eq.0) stop "AllocateBosonicField: Npoints not defined."
      if(Nkpt_.eq.0.and.verbose) write(*,"(A)") "     AllocateBosonicField: K-dependent attributes are not going to be allocated."
      no_bare_=.false.
      if(present(no_bare))no_bare_=no_bare
      if(no_bare_.and.verbose) write(*,"(A)") "     AllocateBosonicField: the bare attributes are not going to be allocated."
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
      W%iq_gamma=iq_gamma
      W%Nbp=Nbp
      W%Npoints=Npoints
      W%Nkpt=Nkpt_
      if(present(Beta))W%Beta=Beta
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
      if(present(name).and.verbose) write(*,"(A)") "     Deallocation of "//trim(name)
      !if(.not.W%status) stop "DeallocateBosonicField: container is unallocated."
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
   !PURPOSE: Clear the internal attributes of a Fermionic/Bosonic field
   !TEST ON: 16-10-2020(both)
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


   !---------------------------------------------------------------------------!
   !PURPOSE: Duplicate a Fermionic/Bosonic field
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine duplicate_Fermionic(Gnew,Gold)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: Gnew
      type(FermionicField),intent(in)       :: Gold
      if(Gnew%status.and.verbose) write(*,"(A)") "Warning duplicate_Fermionic: the new fermionic field will be reinitilized."
      call DeallocateFermionicField(Gnew)
      call AllocateFermionicField(Gnew,Gold%Norb,Gold%Npoints,Nkpt=Gold%Nkpt,Nsite=Gold%Nsite,Beta=Gold%Beta,mu=Gold%mu)
      if(allocated(Gold%N_s)) Gnew%N_s = Gold%N_s
      if(allocated(Gold%ws))  Gnew%ws  = Gold%ws
      if(allocated(Gold%N_ks))Gnew%N_ks= Gold%N_ks
      if(allocated(Gold%wks)) Gnew%wks = Gold%wks
   end subroutine duplicate_Fermionic
   !
   subroutine duplicate_Bosonic(Wnew,Wold)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: Wnew
      type(BosonicField),intent(in)         :: Wold
      if(Wnew%status.and.verbose) write(*,"(A)") "Warning duplicate_Bosonic: the new bosonic field will be reinitilized."
      call DeallocateBosonicField(Wnew)
      call AllocateBosonicField(Wnew,int(sqrt(dble(Wold%Nbp))),Wold%Npoints,Wold%iq_gamma,Nkpt=Wold%Nkpt,Nsite=Wold%Nsite,no_bare=(.not.allocated(Wold%bare)),Beta=Wold%Beta)
      if(allocated(Wold%bare_local))    Wnew%bare_local    = Wold%bare_local
      if(allocated(Wold%screened_local))Wnew%screened_local= Wold%screened_local
      if(allocated(Wold%bare))          Wnew%bare          = Wold%bare
      if(allocated(Wold%screened))      Wnew%screened      = Wold%screened
   end subroutine duplicate_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE: Remove the complex part of the allocated attributes
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine isReal_Matrix(A)
      use parameters
      implicit none
      complex(8),allocatable,intent(inout)  :: A(:,:)
      if(allocated(A))A=dcmplx(real(A),0d0)
   end subroutine isReal_Matrix
   !
   subroutine isReal_Fermionic(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      if(allocated(G%N_s))   G%N_s  = dcmplx(real(G%N_s),0d0)
      if(allocated(G%ws))    G%ws   = dcmplx(real(G%ws),0d0)
      if(allocated(G%N_ks))  G%N_ks = dcmplx(real(G%N_ks),0d0)
      if(allocated(G%wks))   G%wks  = dcmplx(real(G%wks),0d0)
   end subroutine isReal_Fermionic
   !
   subroutine isReal_Bosonic(W)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      if(allocated(W%bare_local))     W%bare_local     = dcmplx(real(W%bare_local),0d0)
      if(allocated(W%screened_local)) W%screened_local = dcmplx(real(W%screened_local),0d0)
      if(allocated(W%bare))           W%bare           = dcmplx(real(W%bare),0d0)
      if(allocated(W%screened))       W%screened       = dcmplx(real(W%screened),0d0)
   end subroutine isReal_Bosonic


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
      if(verbose)write(*,"(A)") "---- loc2imp(F)"
      !
      !
      if(present(sitename)) write(*,"(A)") "     Extraction of Gloc for site: "//trim(sitename)
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
         write(*,"(A)") "     The local orbital space will be rotated during extraction."
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if(size(U,dim=1).ne.size(orbs)) stop "Rotation matrix has the wrong dimension."
         if(size(U,dim=1).ne.3) write(*,"(A)") "     Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
      endif
      !
      call clear_attributes(Gimp)
      !
      do i_imp=1,size(orbs)
         do j_imp=1,size(orbs)
            !
            i_loc = orbs(i_imp)
            j_loc = orbs(j_imp)
            !
            do ispin=1,Nspin
               Gimp%N_s(i_imp,j_imp,ispin) = Gloc%N_s(i_loc,j_loc,ispin)
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
            Gimp%N_s(:,:,ispin) = rotate(Gimp%N_s(:,:,ispin),U)
            do ip=1,Gimp%Npoints
               Gimp%ws(:,:,ip,ispin) = rotate(Gimp%ws(:,:,ip,ispin),U)
            enddo
         enddo
      endif
      !
      Gimp%mu = Gloc%mu
      !
   end subroutine loc2imp_Fermionic
   !
   subroutine loc2imp_Matrix(Oimp,Oloc,orbs,sitename,U)
      !
      use parameters
      use linalg, only : rotate
      implicit none
      !
      complex(8),intent(inout)              :: Oimp(:,:)
      complex(8),intent(in)                 :: Oloc(:,:)
      integer,intent(in)                    :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),intent(in),optional        :: U(:,:)
      !
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      if(verbose)write(*,"(A)") "---- loc2imp(O)"
      !
      !
      if(present(sitename)) write(*,"(A)") "     Extraction of Operator for site: "//trim(sitename)
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "loc2imp(O): Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "loc2imp(O): Oimp not square."
      if(size(orbs).ne.size(Oimp,dim=1)) stop "loc2imp(O): can't fit the requested orbitals inside Oimp."
      if(size(orbs).gt.size(Oloc,dim=1)) stop "loc2imp(O): number of requested orbitals greater than Oloc size."
      if(present(U))then
         write(*,"(A)") "     The local orbital space will be rotated during extraction."
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if(size(U,dim=1).ne.size(orbs)) stop "Rotation matrix has the wrong dimension."
         if(size(U,dim=1).ne.3) write(*,"(A)") "     Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
      endif
      !
      Oimp=czero
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
      logical                               :: doBare
      !
      !
      if(verbose)write(*,"(A)") "---- loc2imp(B)"
      !
      !
      if(present(sitename)) write(*,"(A)") "     Extraction of Wloc for site: "//trim(sitename)
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
      if(.not.allocated(Wimp%screened_local)) stop "loc2imp(B): Wimp screened_local attribute not allocated."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      doBare = allocated(Wimp%bare_local)
      !
      if(size(orbs).ne.Norb_imp) stop "loc2imp(B): can't fit the requested orbitals inside Wimp."
      if(size(orbs).gt.Norb_loc) stop "loc2imp(B): number of requested orbitals greater than Wloc size."
      !
      call clear_attributes(Wimp)
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
                  if(doBare)Wimp%bare_local(ib_imp,jb_imp) =  Wloc%bare_local(ib_loc,jb_loc)
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
   subroutine imp2loc_Fermionic(Gloc,Gimp,orbs,expand,AFM,U)
      !
      use parameters
      use utils_misc
      use linalg, only : rotate
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gloc
      type(FermionicField),intent(in)       :: Gimp
      integer,allocatable,intent(in)        :: orbs(:)
      logical,intent(in)                    :: expand
      logical,intent(in)                    :: AFM
      complex(8),allocatable,optional       :: U(:,:,:)
      !
      complex(8),allocatable                :: Rot(:,:)
      complex(8),allocatable                :: Gtmp(:,:,:,:,:),Ntmp(:,:,:,:)
      integer                               :: ip,ispin,ispin_imp,isite,Nsite,shift
      integer                               :: Norb_loc,Norb_imp
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      if(verbose)write(*,"(A)") "---- imp2loc(F)"
      !
      !
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
      Nsite=1
      !
      if(AFM)then
         if(Gloc%Nsite.ne.2) stop "AFM is implemented only for a two site lattice."
         if(Norb_loc/Norb_imp.ne.2) stop "Lattice indexes are not twice the impurity ones."
         if(expand) stop "AFM condition and expansion to real space not yet implemented."
         Nsite = 2
      endif
      !
      if(expand)then
         Nsite = Gloc%Nsite
         if(mod(Norb_loc,size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of the lattice field."
         if(AFM) stop "Expansion to real space and AFM condition not yet implemented."
         write(*,"(A)") "     The impurity field will be expanded to match the lattice orbital space."
      else
         write(*,"(A,15I3)") "     Impurity field will be inserted into the lattice orbital indexes: ",orbs
         Nsite = 1
      endif
      !
      if(present(U))then
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if(size(U,dim=3).ne.Gloc%Nsite) stop "Number of rotation matrices and number of sites does not match."
         if(size(U,dim=1).ne.3) write(*,"(A)") "     Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
         write(*,"(A)") "     The impurity orbital space will be rotated during insertion in "//str(Nsite)//" sites."
      endif
      !
      allocate(Ntmp(Gimp%Norb,Gimp%Norb,Nspin,Nsite));Ntmp=czero
      allocate(Gtmp(Gimp%Norb,Gimp%Norb,Gimp%Npoints,Nspin,Nsite));Gtmp=czero
      !
      ! Rotating either one site or all of them depending on expand
      do isite=1,Nsite
         !
         if(present(U))then
            allocate(Rot(Norb_imp,Norb_imp)); Rot=U(1:Norb_imp,1:Norb_imp,isite)
            do ispin=1,Nspin
               Ntmp(:,:,ispin,isite) = rotate(Gimp%N_s(:,:,ispin),Rot)
               do ip=1,Gimp%Npoints
                  Gtmp(:,:,ip,ispin,isite) = rotate(Gimp%ws(:,:,ip,ispin),Rot)
               enddo
            enddo
            deallocate(Rot)
         else
            Ntmp(:,:,:,isite) = Gimp%N_s
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
                  ispin_imp=ispin
                  if(AFM.and.(isite.eq.2)) ispin_imp=int(Nspin/ispin)
                  !
                  Gloc%N_s(i_loc,j_loc,ispin) = Ntmp(i_imp,j_imp,ispin_imp,isite)
                  do ip=1,Gimp%Npoints
                     Gloc%ws(i_loc,j_loc,ip,ispin) = Gtmp(i_imp,j_imp,ip,ispin_imp,isite)
                  enddo
               enddo
               !
            enddo
         enddo
         !
      enddo
      deallocate(Gtmp,Ntmp)
      !
   end subroutine imp2loc_Fermionic
   !
   subroutine imp2loc_Matrix(Oloc,Oimp,orbs,expand,AFM,U)
      !
      use parameters
      use utils_misc
      use linalg, only : rotate
      implicit none
      !
      complex(8),intent(inout)              :: Oloc(:,:,:)
      complex(8),intent(in)                 :: Oimp(:,:,:)
      integer,intent(in)                    :: orbs(:)
      logical,intent(in)                    :: expand
      logical,intent(in)                    :: AFM
      complex(8),intent(in),optional        :: U(:,:,:)
      !
      complex(8),allocatable                :: Rot(:,:)
      complex(8),allocatable                :: Otmp(:,:,:,:)
      integer                               :: ispin,ispin_imp,isite,Nsite,shift
      integer                               :: Norb_loc,Norb_imp
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      if(verbose)write(*,"(A)") "---- imp2loc(O)"
      !
      !
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "imp2loc(O): Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "imp2loc(O): Oimp not square."
      if(size(Oimp,dim=3).ne.Nspin) stop "imp2loc(O): Oimp third dimension is not Nspin."
      if(size(Oloc,dim=3).ne.Nspin) stop "imp2loc(O): Oloc third dimension is not Nspin."
      if(size(orbs).ne.size(Oimp,dim=1)) stop "imp2loc(O): can't fit the requested orbitals from Gimp."
      if(size(orbs).gt.size(Oloc,dim=1)) stop "imp2loc(O): number of requested orbitals greater than Gloc size."
      !
      Norb_loc = size(Oloc,dim=1)
      Norb_imp = size(orbs)
      Nsite=1
      !
      if(AFM)then
         if(Norb_loc/Norb_imp.ne.2) stop "Lattice indexes are not twice the impurity ones."
         if(expand) stop "AFM condition and expansion to real space not yet implemented."
         Nsite = 2
      endif
      !
      if(expand)then
         Nsite = size(Oloc,dim=1)/size(orbs)
         if(mod(size(Oloc,dim=1),size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of the lattice observable."
         if(AFM) stop "Expansion to real space and AFM condition not yet implemented."
         write(*,"(A)") "     The impurity field will be expanded to match the lattice orbital space."
      else
         write(*,"(A,15I3)") "     Impurity field will be inserted into the lattice orbital indexes: ",orbs
         Nsite = 1
      endif
      !
      if(present(U))then
         if(size(U,dim=1).ne.size(U,dim=2)) stop "Rotation matrix not square."
         if((size(U,dim=3).ne.Nsite))stop "Number of rotation matrices and number of sites does not match."
         if(size(U,dim=1).ne.3) write(*,"(A)") "     Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
         write(*,"(A)") "     The impurity orbital space will be rotated during insertion in "//str(Nsite)//" sites."
      endif
      !
      allocate(Otmp(size(Oloc,dim=1),size(Oloc,dim=1),Nspin,Nsite));Otmp=czero
      !
      ! Rotating either one site or all of them depending on expand
      do isite=1,Nsite
         !
         if(present(U))then
            allocate(Rot(Norb_imp,Norb_imp)); Rot=U(1:Norb_imp,1:Norb_imp,isite)
            do ispin=1,Nspin
               Otmp(:,:,ispin,isite) = rotate(Oimp(:,:,ispin),Rot)
            enddo
            deallocate(Rot)
         else
            Otmp(:,:,:,isite) = Oimp
         endif
         !
      enddo
      !
      do isite=1,Nsite
         !
         ! only two possible arrangements
         if(abs(orbs(2)-orbs(1)).eq.1)then
            shift = Norb_imp*(isite-1)
         elseif(abs(orbs(2)-orbs(1)).eq.Nsite)then
            shift = isite-1
         endif
         !
         do i_imp=1,Norb_imp
            do j_imp=1,Norb_imp
               !
               i_loc = orbs(i_imp) + shift
               j_loc = orbs(j_imp) + shift
               !
               do ispin=1,Nspin
                  !
                  ispin_imp=ispin
                  if(AFM.and.(isite.eq.2)) ispin_imp=int(Nspin/ispin)
                  !
                  Oloc(i_loc,j_loc,ispin) = Otmp(i_imp,j_imp,ispin_imp,isite)
                  !
               enddo
               !
            enddo
         enddo
         !
      enddo
      deallocate(Otmp)
      !
   end subroutine imp2loc_Matrix
   !
   subroutine imp2loc_Bosonic(Wloc,Wimp,orbs,expand,AFM)
      !
      use parameters
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wloc
      type(BosonicField),intent(in)         :: Wimp
      integer,allocatable,intent(in)        :: orbs(:)
      logical,intent(in)                    :: expand
      logical,intent(in)                    :: AFM
      !
      integer                               :: ip,isite,Nsite,shift
      integer                               :: Norb_imp,Norb_loc
      integer                               :: ib_imp,jb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: i_imp,j_imp,k_imp,l_imp
      logical                               :: doBare
      !
      !
      if(verbose)write(*,"(A)") "---- imp2loc(B)"
      !
      !
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
      if(.not.allocated(Wimp%screened_local)) stop "imp2loc(B): Wimp screened_local attribute not allocated."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      doBare = allocated(Wimp%bare_local)
      Nsite=1
      !
      if(size(orbs).ne.Norb_imp) stop "imp2loc(B): can't fit the requested orbitals from Wimp."
      if(size(orbs).gt.Norb_loc) stop "imp2loc(B): number of requested orbitals greater than Wloc size."
      !
      if(AFM)then
         if(Wloc%Nsite.ne.2) stop "AFM is implemented only for a two site lattice."
         if(Norb_loc/Norb_imp.ne.2) stop "Lattice indexes are not twice the impurity ones."
         if(expand) stop "AFM condition and expansion to real space not yet implemented."
         Nsite = 2
      endif
      !
      if(expand)then
         Nsite = Wloc%Nsite
         if(mod(Norb_loc,size(orbs)).ne.0)stop "Number of requested orbitals is not a commensurate subset of Wloc."
         if(AFM) stop "Expansion to real space and AFM condition not yet implemented."
         write(*,"(A)") "     The impurity field will be expanded to match the lattice orbital space."
      else
         write(*,"(A,15I3)") "     Impurity field will be inserted into the lattice orbital indexes: ",orbs
         Nsite = 1
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
                     if(doBare)Wloc%bare_local(ib_loc,jb_loc) = Wimp%bare_local(ib_imp,jb_imp)
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
   !PURPOSE: Symmetrize the Matrix/Field with respect to some lattice indexes.
   !         The Fermionic simmetrization is done only on the local attributes
   !         The Bosonic simmetrization is done only on the physical elements
   !---------------------------------------------------------------------------!
   subroutine symmetrize_Matrix(Mat,Eqv)
      !
      use parameters
      implicit none
      !
      complex(8),intent(inout)              :: Mat(:,:,:)
      type(Equivalent)                      :: Eqv
      !
      real(8)                               :: dimdiag,dimoffdiag
      complex(8)                            :: Delem,Uelem,Lelem
      integer                               :: iset,jset,iorb,jorb,ispin
      integer                               :: i,j
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_Matrix"
      !
      !
      ! Check on the input Fields
      if(size(Mat,dim=1).ne.size(Mat,dim=2)) stop "symmetrize_Matrix: Matix not square."
      !
      if(Eqv%para.eq.1)then
         Mat(:,:,1) = (Mat(:,:,1) + Mat(:,:,2))/2d0
         Mat(:,:,2) = Mat(:,:,1)
      endif
      !
      if(Eqv%O)then
         !
         do ispin=1,Nspin
            !
            !symmetrization of the diagonal sets
            do iset=1,Eqv%Ntotset
               !
               dimdiag = Eqv%SetNorb(iset)
               dimoffdiag = Eqv%SetNorb(iset)*(Eqv%SetNorb(iset)-1)/2
               !
               if(dimdiag.eq.1)cycle
               if(verbose)write(*,"(3(A,I4))") "     Diagonal set: ",iset," dimD: ",int(dimdiag)," dimOD: ",int(dimoffdiag)
               !
               !Average elements
               Delem=czero;Lelem=czero;Uelem=czero
               do iorb=1,Eqv%SetNorb(iset)
                  do jorb=1,Eqv%SetNorb(iset)
                     !
                     i = Eqv%SetOrbs(iset,iorb)
                     j = Eqv%SetOrbs(iset,jorb)
                     if(verbose)write(*,"(A,2I4)")    "     Component: ",i,j
                     !
                     if(i.eq.j) Delem = Delem + Mat(i,j,ispin)/dimdiag
                     if(i.gt.j) Lelem = Lelem + Mat(i,j,ispin)/dimoffdiag
                     if(i.lt.j) Uelem = Uelem + Mat(i,j,ispin)/dimoffdiag
                  enddo
               enddo
               !Re-insert elements
               do iorb=1,Eqv%SetNorb(iset)
                  do jorb=1,Eqv%SetNorb(iset)
                     !
                     i = Eqv%SetOrbs(iset,iorb)
                     j = Eqv%SetOrbs(iset,jorb)
                     !
                     if(i.eq.j) Mat(i,j,ispin) = Delem
                     if((i.gt.j).and.Eqv%Gfoffdiag) Mat(i,j,ispin) = Lelem
                     if((i.lt.j).and.Eqv%Gfoffdiag) Mat(i,j,ispin) = Uelem
                  enddo
               enddo
               !
            enddo !iset
            !
            !symmetrization of the off-diagonal sets
            if(Eqv%Gfoffdiag)then
               do iset=1,Eqv%Ntotset
                  do jset=1,Eqv%Ntotset
                     !
                     if(iset.eq.jset)cycle
                     dimoffdiag = Eqv%SetNorb(iset)*Eqv%SetNorb(jset)
                     !
                     if(verbose)write(*,"(A,2I4,A,I4)") "     Off-diagonal set: ",iset,jset," dimOD: ",int(dimoffdiag)
                     !
                     !Average elements
                     Lelem=czero;Uelem=czero
                     do iorb=1,Eqv%SetNorb(iset)
                        do jorb=1,Eqv%SetNorb(jset)
                           !
                           i = Eqv%SetOrbs(iset,iorb)
                           j = Eqv%SetOrbs(jset,jorb)
                           if(verbose)write(*,"(A,2I4)")    "     Component: ",i,j
                           !
                           if(iset.gt.jset) Lelem = Lelem + Mat(i,j,ispin)/dimoffdiag
                           if(iset.lt.jset) Uelem = Uelem + Mat(i,j,ispin)/dimoffdiag
                        enddo
                     enddo
                     !Re-insert elements
                     do iorb=1,Eqv%SetNorb(iset)
                        do jorb=1,Eqv%SetNorb(jset)
                           !
                           i = Eqv%SetOrbs(iset,iorb)
                           j = Eqv%SetOrbs(jset,jorb)
                           !
                           if(iset.gt.jset) Mat(i,j,ispin) = Lelem
                           if(iset.lt.jset) Mat(i,j,ispin) = Uelem
                        enddo
                     enddo
                     !
                  enddo
               enddo
            endif
            !
         enddo !ispin
         !
      endif
      !
   end subroutine symmetrize_Matrix
   !
   subroutine symmetrize_Fermionic(G,Eqv)
      !
      use parameters
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      type(Equivalent)                      :: Eqv
      !
      real(8)                               :: dimdiag,dimoffdiag
      complex(8)                            :: Delem,Uelem,Lelem
      integer                               :: iset,jset,iorb,jorb,ispin
      integer                               :: i,j,ip
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_Fermionic"
      !
      !
      ! Check on the input Fields
      if(.not.G%status) stop "symmetrize_Fermionic: field not properly initialized."
      !
      if(Eqv%para.eq.1)then
         !
         G%N_s(:,:,1) = (G%N_s(:,:,1) + G%N_s(:,:,2))/2d0
         G%N_s(:,:,2) = G%N_s(:,:,1)
         !
         do ip=1,G%Npoints
            G%ws(:,:,ip,1) = (G%ws(:,:,ip,1) + G%ws(:,:,ip,2))/2d0
            G%ws(:,:,ip,2) = G%ws(:,:,ip,1)
         enddo
         !
      endif
      !
      if(Eqv%O)then
         !
         do ispin=1,Nspin
            !
            !symmetrization of the diagonal sets
            do iset=1,Eqv%Ntotset
               !
               dimdiag = Eqv%SetNorb(iset)
               dimoffdiag = Eqv%SetNorb(iset)*(Eqv%SetNorb(iset)-1)/2
               !
               if(dimdiag.eq.1)cycle
               if(verbose.and.(ispin.eq.1))write(*,"(3(A,I4))") "     Diagonal set: ",iset," dimD: ",int(dimdiag)," dimOD: ",int(dimoffdiag)
               !
               !Average elements
               Delem=czero;Lelem=czero;Uelem=czero
               do iorb=1,Eqv%SetNorb(iset)
                  do jorb=1,Eqv%SetNorb(iset)
                     !
                     i = Eqv%SetOrbs(iset,iorb)
                     j = Eqv%SetOrbs(iset,jorb)
                     if(verbose.and.(ispin.eq.1))write(*,"(A,2I4)")    "     Component: ",i,j
                     !
                     if(i.eq.j) Delem = Delem + G%N_s(i,j,ispin)/dimdiag
                     if(i.gt.j) Lelem = Lelem + G%N_s(i,j,ispin)/dimoffdiag
                     if(i.lt.j) Uelem = Uelem + G%N_s(i,j,ispin)/dimoffdiag
                     !
                  enddo
               enddo
               !Re-insert elements
               do iorb=1,Eqv%SetNorb(iset)
                  do jorb=1,Eqv%SetNorb(iset)
                     !
                     i = Eqv%SetOrbs(iset,iorb)
                     j = Eqv%SetOrbs(iset,jorb)
                     !
                     if(i.eq.j) G%N_s(i,j,ispin) = Delem
                     if((i.gt.j).and.Eqv%Gfoffdiag) G%N_s(i,j,ispin) = Lelem
                     if((i.lt.j).and.Eqv%Gfoffdiag) G%N_s(i,j,ispin) = Uelem
                     !
                  enddo
               enddo
               !
               do ip=1,G%Npoints
                  !Average elements
                  Delem=czero;Lelem=czero;Uelem=czero
                  do iorb=1,Eqv%SetNorb(iset)
                     do jorb=1,Eqv%SetNorb(iset)
                        !
                        i = Eqv%SetOrbs(iset,iorb)
                        j = Eqv%SetOrbs(iset,jorb)
                        !
                        if(i.eq.j) Delem = Delem + G%ws(i,j,ip,ispin)/dimdiag
                        if(i.gt.j) Lelem = Lelem + G%ws(i,j,ip,ispin)/dimoffdiag
                        if(i.lt.j) Uelem = Uelem + G%ws(i,j,ip,ispin)/dimoffdiag
                        !
                     enddo
                  enddo
                  !Re-insert elements
                  do iorb=1,Eqv%SetNorb(iset)
                     do jorb=1,Eqv%SetNorb(iset)
                        !
                        i = Eqv%SetOrbs(iset,iorb)
                        j = Eqv%SetOrbs(iset,jorb)
                        !
                        if(i.eq.j) G%ws(i,j,ip,ispin) = Delem
                        if((i.gt.j).and.Eqv%Gfoffdiag) G%ws(i,j,ip,ispin) = Lelem
                        if((i.lt.j).and.Eqv%Gfoffdiag) G%ws(i,j,ip,ispin) = Uelem
                        !
                     enddo
                  enddo
               enddo! ip
               !
            enddo !iset
            !
            !symmetrization of the off-diagonal sets
            if(Eqv%Gfoffdiag)then
               do iset=1,Eqv%Ntotset
                  do jset=1,Eqv%Ntotset
                     !
                     if(iset.eq.jset)cycle
                     dimoffdiag = Eqv%SetNorb(iset)*Eqv%SetNorb(jset)
                     !
                     if(verbose.and.(ispin.eq.1))write(*,"(A,2I4,A,I4)") "     Off-diagonal set: ",iset,jset," dimOD: ",int(dimoffdiag)
                     !
                     !Average elements
                     Lelem=czero;Uelem=czero
                     do iorb=1,Eqv%SetNorb(iset)
                        do jorb=1,Eqv%SetNorb(jset)
                           !
                           i = Eqv%SetOrbs(iset,iorb)
                           j = Eqv%SetOrbs(jset,jorb)
                           !
                           if(iset.gt.jset) Lelem = Lelem + G%N_s(i,j,ispin)/dimoffdiag
                           if(iset.lt.jset) Uelem = Uelem + G%N_s(i,j,ispin)/dimoffdiag
                           !
                        enddo
                     enddo
                     !Re-insert elements
                     do iorb=1,Eqv%SetNorb(iset)
                        do jorb=1,Eqv%SetNorb(jset)
                           !
                           i = Eqv%SetOrbs(iset,iorb)
                           j = Eqv%SetOrbs(jset,jorb)
                           if(verbose.and.(ispin.eq.1))write(*,"(A,2I4)")    "     Component: ",i,j
                           !
                           if(iset.gt.jset) G%N_s(i,j,ispin) = Lelem
                           if(iset.lt.jset) G%N_s(i,j,ispin) = Uelem
                           !
                        enddo
                     enddo
                     !
                     do ip=1,G%Npoints
                        !Average elements
                        Lelem=czero;Uelem=czero
                        do iorb=1,Eqv%SetNorb(iset)
                           do jorb=1,Eqv%SetNorb(jset)
                              !
                              i = Eqv%SetOrbs(iset,iorb)
                              j = Eqv%SetOrbs(jset,jorb)
                              !
                              if(iset.gt.jset) Lelem = Lelem + G%ws(i,j,ip,ispin)/dimoffdiag
                              if(iset.lt.jset) Uelem = Uelem + G%ws(i,j,ip,ispin)/dimoffdiag
                              !
                           enddo
                        enddo
                        !Re-insert elements
                        do iorb=1,Eqv%SetNorb(iset)
                           do jorb=1,Eqv%SetNorb(jset)
                              !
                              i = Eqv%SetOrbs(iset,iorb)
                              j = Eqv%SetOrbs(jset,jorb)
                              !
                              if(iset.gt.jset) G%ws(i,j,ip,ispin) = Lelem
                              if(iset.lt.jset) G%ws(i,j,ip,ispin) = Uelem
                              !
                           enddo
                        enddo
                     enddo !ip
                     !
                  enddo
               enddo
            endif
            !
         enddo !ispin
         !
      endif
      !
   end subroutine symmetrize_Fermionic
   !
   !
   subroutine symmetrize_Bosonic(W,Eqv)
      !
      use parameters
      implicit none
      !
      type(BosonicField),intent(inout)      :: W
      type(Equivalent)                      :: Eqv
      !
      real(8)                               :: dimdiag,dimoffdiag
      complex(8)                            :: Waaaa,Waabb
      complex(8)                            :: Wabab,Wabba
      integer                               :: iset,jset,ip,Norb
      integer                               :: iorb,jorb,i,j,ib1_aa
      integer                               :: ib1_ab,ib2_ab,ib1_sf,ib2_sf
      integer                               :: ib1_pa,ib2_pa,ib1_pb,ib2_pb
      logical                               :: doBare
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_Bosonic"
      !
      !
      ! Check on the input Fields
      if(.not.W%status) stop "symmetrize_Bosonic: field not properly initialized."
      Norb = sqrt(dble(W%Nbp))
      doBare = allocated(W%bare_local)
      !
      if(Eqv%O)then
         !
         !symmetrization of the diagonal sets
         do iset=1,Eqv%Ntotset
            !
            dimdiag  = Eqv%SetNorb(iset)
            dimoffdiag = Eqv%SetNorb(iset)*(Eqv%SetNorb(iset)-1)
            !
            if(verbose)write(*,"(3(A,I3))") "     Diagonal set: ",iset," dimD: ",int(dimdiag)," dimOD: ",int(dimoffdiag)
            !
            !Bare attribute
            if(doBare)then
               !
               !Average elements
               Waaaa=czero;Waabb=czero;Wabba=czero;Wabab=czero
               do iorb=1,Eqv%SetNorb(iset)
                  i = Eqv%SetOrbs(iset,iorb)
                  !
                  ib1_aa = i + Norb*(i-1)
                  if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(aa)(aa): (",i,i,"),(",i,i,") = [",ib1_aa,",",ib1_aa,"]"
                  Waaaa = Waaaa + W%bare_local(ib1_aa,ib1_aa) / dimdiag
                  !
                  do jorb=1+iorb,Eqv%SetNorb(iset)
                     j = Eqv%SetOrbs(iset,jorb)
                     !
                     ib1_ab = i + Norb*(i-1)
                     ib2_ab = j + Norb*(j-1)
                     if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(aa)(bb): (",i,i,"),(",j,j,") = [",ib1_ab,",",ib2_ab,"]"
                     !
                     ib1_sf = i + Norb*(j-1)
                     ib2_sf = j + Norb*(i-1)
                     if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ba): (",i,j,"),(",j,i,") = [",ib1_sf,",",ib2_sf,"]"
                     !
                     ib1_pa = i + Norb*(j-1)
                     ib2_pa = i + Norb*(j-1)
                     if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ab): (",i,j,"),(",i,j,") = [",ib1_pa,",",ib2_pa,"]"
                     ib1_pb = j + Norb*(i-1)
                     ib2_pb = j + Norb*(i-1)
                     if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ba)(ba): (",j,i,"),(",j,i,") = [",ib1_pb,",",ib2_pb,"]"
                     !
                     Waabb = Waabb + (W%bare_local(ib1_ab,ib2_ab)+W%bare_local(ib2_ab,ib1_ab)) / dimoffdiag
                     Wabba = Wabba + (W%bare_local(ib1_sf,ib2_sf)+W%bare_local(ib2_sf,ib1_sf)) / dimoffdiag
                     Wabab = Wabab + (W%bare_local(ib1_pa,ib2_pa)+W%bare_local(ib1_pb,ib2_pb)) / dimoffdiag
                     !
                  enddo
               enddo
               !Re-insert elements
               do iorb=1,Eqv%SetNorb(iset)
                  i = Eqv%SetOrbs(iset,iorb)
                  !
                  ib1_aa = i + Norb*(i-1)
                  W%bare_local(ib1_aa,ib1_aa) = Waaaa
                  !
                  do jorb=1+iorb,Eqv%SetNorb(iset)
                     j = Eqv%SetOrbs(iset,jorb)
                     !
                     ib1_ab = i + Norb*(i-1)
                     ib2_ab = j + Norb*(j-1)
                     !
                     ib1_sf = i + Norb*(j-1)
                     ib2_sf = j + Norb*(i-1)
                     !
                     ib1_pa = i + Norb*(j-1)
                     ib2_pa = i + Norb*(j-1)
                     ib1_pb = j + Norb*(i-1)
                     ib2_pb = j + Norb*(i-1)
                     !
                     W%bare_local(ib1_ab,ib2_ab) = Waabb
                     W%bare_local(ib2_ab,ib1_ab) = Waabb
                     W%bare_local(ib1_sf,ib2_sf) = Wabba
                     W%bare_local(ib2_sf,ib1_sf) = Wabba
                     W%bare_local(ib1_pa,ib2_pa) = Wabab
                     W%bare_local(ib1_pb,ib2_pb) = Wabab
                     !
                  enddo
               enddo
               !
            endif
            !
            !Screened attribute
            do ip=1,W%Npoints
               !
               !Average elements
               Waaaa=czero;Waabb=czero;Wabba=czero;Wabab=czero
               do iorb=1,Eqv%SetNorb(iset)
                  i = Eqv%SetOrbs(iset,iorb)
                  !
                  ib1_aa = i + Norb*(i-1)
                  Waaaa = Waaaa + W%screened_local(ib1_aa,ib1_aa,ip) / dimdiag
                  !
                  do jorb=1+iorb,Eqv%SetNorb(iset)
                     j = Eqv%SetOrbs(iset,jorb)
                     !
                     ib1_ab = i + Norb*(i-1)
                     ib2_ab = j + Norb*(j-1)
                     !
                     ib1_sf = i + Norb*(j-1)
                     ib2_sf = j + Norb*(i-1)
                     !
                     ib1_pa = i + Norb*(j-1)
                     ib2_pa = i + Norb*(j-1)
                     ib1_pb = j + Norb*(i-1)
                     ib2_pb = j + Norb*(i-1)
                     !
                     Waabb = Waabb + (W%screened_local(ib1_ab,ib2_ab,ip)+W%screened_local(ib2_ab,ib1_ab,ip)) / dimoffdiag
                     Wabba = Wabba + (W%screened_local(ib1_sf,ib2_sf,ip)+W%screened_local(ib2_sf,ib1_sf,ip)) / dimoffdiag
                     Wabab = Wabab + (W%screened_local(ib1_pa,ib2_pa,ip)+W%screened_local(ib1_pb,ib2_pb,ip)) / dimoffdiag
                     !
                  enddo
               enddo
               !Re-insert elements
               do iorb=1,Eqv%SetNorb(iset)
                  i = Eqv%SetOrbs(iset,iorb)
                  !
                  ib1_aa = i + Norb*(i-1)
                  W%screened_local(ib1_aa,ib1_aa,ip) = Waaaa
                  !
                  do jorb=1+iorb,Eqv%SetNorb(iset)
                     j = Eqv%SetOrbs(iset,jorb)
                     !
                     ib1_ab = i + Norb*(i-1)
                     ib2_ab = j + Norb*(j-1)
                     !
                     ib1_sf = i + Norb*(j-1)
                     ib2_sf = j + Norb*(i-1)
                     !
                     ib1_pa = i + Norb*(j-1)
                     ib2_pa = i + Norb*(j-1)
                     ib1_pb = j + Norb*(i-1)
                     ib2_pb = j + Norb*(i-1)
                     !
                     W%screened_local(ib1_ab,ib2_ab,ip) = Waabb
                     W%screened_local(ib2_ab,ib1_ab,ip) = Waabb
                     W%screened_local(ib1_sf,ib2_sf,ip) = Wabba
                     W%screened_local(ib2_sf,ib1_sf,ip) = Wabba
                     W%screened_local(ib1_pa,ib2_pa,ip) = Wabab
                     W%screened_local(ib1_pb,ib2_pb,ip) = Wabab
                     !
                  enddo
               enddo
               !
            enddo !ip
            !
         enddo !iset - diagonal
         !
         !symmetrization of the off-diagonal sets
         do iset=1,Eqv%Ntotset
            do jset=1+iset,Eqv%Ntotset
               !
               dimoffdiag = Eqv%SetNorb(iset)*Eqv%SetNorb(jset)*2
               !
               if(verbose)write(*,"(A,2I3,A,I4)") "     Off-diagonal set: ",iset,jset," dimOD: ",int(dimoffdiag)
               !
               !Bare attribute
               if(doBare)then
                  !
                  !Average elements
                  Waabb=czero;Wabba=czero;Wabab=czero
                  do iorb=1,Eqv%SetNorb(iset)
                     do jorb=1,Eqv%SetNorb(jset)
                        !
                        i = Eqv%SetOrbs(iset,iorb)
                        j = Eqv%SetOrbs(jset,jorb)
                        if(i.eq.j) stop "symmetrize_Bosonic: different sets cannot contain the same orbital index."
                        !
                        ib1_ab = i + Norb*(i-1)
                        ib2_ab = j + Norb*(j-1)
                        if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(aa)(bb): (",i,i,"),(",j,j,") = [",ib1_ab,",",ib2_ab,"]"
                        !
                        ib1_sf = i + Norb*(j-1)
                        ib2_sf = j + Norb*(i-1)
                        if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ba): (",i,j,"),(",j,i,") = [",ib1_sf,",",ib2_sf,"]"
                        !
                        ib1_pa = i + Norb*(j-1)
                        ib2_pa = i + Norb*(j-1)
                        if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ab): (",i,j,"),(",i,j,") = [",ib1_pa,",",ib2_pa,"]"
                        ib1_pb = j + Norb*(i-1)
                        ib2_pb = j + Norb*(i-1)
                        if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ba)(ba): (",j,i,"),(",j,i,") = [",ib1_pb,",",ib2_pb,"]"
                        !
                        Waabb = Waabb + (W%bare_local(ib1_ab,ib2_ab)+W%bare_local(ib2_ab,ib1_ab)) / dimoffdiag
                        Wabba = Wabba + (W%bare_local(ib1_sf,ib2_sf)+W%bare_local(ib2_sf,ib1_sf)) / dimoffdiag
                        Wabab = Wabab + (W%bare_local(ib1_pa,ib2_pa)+W%bare_local(ib1_pb,ib2_pb)) / dimoffdiag
                        !
                     enddo
                  enddo
                  !Re-insert elements
                  do iorb=1,Eqv%SetNorb(iset)
                     do jorb=1,Eqv%SetNorb(jset)
                        !
                        i = Eqv%SetOrbs(iset,iorb)
                        j = Eqv%SetOrbs(jset,jorb)
                        !
                        ib1_ab = i + Norb*(i-1)
                        ib2_ab = j + Norb*(j-1)
                        !
                        ib1_sf = i + Norb*(j-1)
                        ib2_sf = j + Norb*(i-1)
                        !
                        ib1_pa = i + Norb*(j-1)
                        ib2_pa = i + Norb*(j-1)
                        ib1_pb = j + Norb*(i-1)
                        ib2_pb = j + Norb*(i-1)
                        !
                        W%bare_local(ib1_ab,ib2_ab) = Waabb
                        W%bare_local(ib2_ab,ib1_ab) = Waabb
                        W%bare_local(ib1_sf,ib2_sf) = Wabba
                        W%bare_local(ib2_sf,ib1_sf) = Wabba
                        W%bare_local(ib1_pa,ib2_pa) = Wabab
                        W%bare_local(ib1_pb,ib2_pb) = Wabab
                        !
                     enddo
                  enddo
                  !
               endif
               !
               !Screened attribute
               do ip=1,W%Npoints
                  !
                  !Average elements
                  Waabb=czero;Wabba=czero;Wabab=czero
                  do iorb=1,Eqv%SetNorb(iset)
                     do jorb=1,Eqv%SetNorb(jset)
                        !
                        i = Eqv%SetOrbs(iset,iorb)
                        j = Eqv%SetOrbs(jset,jorb)
                        if(i.eq.j) stop "symmetrize_Bosonic: different sets cannot contain the same orbital index."
                        !
                        ib1_ab = i + Norb*(i-1)
                        ib2_ab = j + Norb*(j-1)
                        !
                        ib1_sf = i + Norb*(j-1)
                        ib2_sf = j + Norb*(i-1)
                        !
                        ib1_pa = i + Norb*(j-1)
                        ib2_pa = i + Norb*(j-1)
                        ib1_pb = j + Norb*(i-1)
                        ib2_pb = j + Norb*(i-1)
                        !
                        Waabb = Waabb + (W%screened_local(ib1_ab,ib2_ab,ip)+W%screened_local(ib2_ab,ib1_ab,ip)) / dimoffdiag
                        Wabba = Wabba + (W%screened_local(ib1_sf,ib2_sf,ip)+W%screened_local(ib2_sf,ib1_sf,ip)) / dimoffdiag
                        Wabab = Wabab + (W%screened_local(ib1_pa,ib2_pa,ip)+W%screened_local(ib1_pb,ib2_pb,ip)) / dimoffdiag
                        !
                     enddo
                  enddo
                  !Re-insert elements
                  do iorb=1,Eqv%SetNorb(iset)
                     do jorb=1,Eqv%SetNorb(jset)
                        !
                        i = Eqv%SetOrbs(iset,iorb)
                        j = Eqv%SetOrbs(jset,jorb)
                        !
                        ib1_ab = i + Norb*(i-1)
                        ib2_ab = j + Norb*(j-1)
                        !
                        ib1_sf = i + Norb*(j-1)
                        ib2_sf = j + Norb*(i-1)
                        !
                        ib1_pa = i + Norb*(j-1)
                        ib2_pa = i + Norb*(j-1)
                        ib1_pb = j + Norb*(i-1)
                        ib2_pb = j + Norb*(i-1)
                        !
                        W%screened_local(ib1_ab,ib2_ab,ip) = Waabb
                        W%screened_local(ib2_ab,ib1_ab,ip) = Waabb
                        W%screened_local(ib1_sf,ib2_sf,ip) = Wabba
                        W%screened_local(ib2_sf,ib1_sf,ip) = Wabba
                        W%screened_local(ib1_pa,ib2_pa,ip) = Wabab
                        W%screened_local(ib1_pb,ib2_pb,ip) = Wabab
                        !
                     enddo
                  enddo
                  !
               enddo !ip
               !
            enddo !jset
         enddo !iset
         !
      endif
      !
   end subroutine symmetrize_Bosonic




   !---------------------------------------------------------------------------!
   !PURPOSE: Join the C and X component of the self-energy
   !TEST ON: 27-10-2020
   !---------------------------------------------------------------------------!
   subroutine join_SigmaCX(SigmaFull,Sigma_C,Sigma_X)
      !
      use parameters
      implicit none
      !
      type(FermionicField),intent(inout)    :: SigmaFull
      type(FermionicField),intent(in)       :: Sigma_C
      type(FermionicField),intent(in)       :: Sigma_X
      !
      real(8)                               :: Beta
      integer                               :: Nkpt,Norb,Nmats
      integer                               :: iw,ik,ispin
      logical                               :: local
      !
      !
      if(verbose)write(*,"(A)") "---- join_SigmaCX"
      !
      !
      ! Check on the input Fields
      if(.not.SigmaFull%status) stop "SigmaFull not properly initialized."
      if(.not.Sigma_C%status) stop "Sigma_C not properly initialized."
      if(.not.Sigma_X%status) stop "Sigma_X not properly initialized."
      !
      if(SigmaFull%Nkpt.eq.0) local=.true.
      if(SigmaFull%Nkpt.ne.0) local=.false.
      !
      if(local)then
         if(Sigma_C%Nkpt.ne.0) stop "Sigma_C k dependent attributes are supposed to be unallocated."
         if(Sigma_X%Nkpt.ne.0) stop "Sigma_X k dependent attributes are supposed to be unallocated."
      else
         if(Sigma_C%Nkpt.eq.0) stop "Sigma_C k dependent attributes not properly initialized."
         if(Sigma_X%Nkpt.eq.0) stop "Sigma_X k dependent attributes not properly initialized."
      endif
      if(Sigma_X%Npoints.ne.0) stop "Sigma_X frequency dependent attributes are supposed to be unallocated."
      !
      Norb = SigmaFull%Norb
      Beta = SigmaFull%Beta
      Nmats = SigmaFull%Npoints
      !
      if(.not.local)then
         Nkpt = SigmaFull%Nkpt
         if(all([Sigma_C%Nkpt-Nkpt,Sigma_X%Nkpt-Nkpt].ne.[0,0])) stop "Either Sigma_C or Sigma_X have different number of k-points with respect to SigmaFull."
      endif
      !
      if(all([Sigma_C%Beta-Beta,Sigma_X%Beta-Beta].ne.[0d0,0d0])) stop "Either Sigma_C or Sigma_X have different Beta with respect to SigmaFull."
      if(Nmats.ne.Sigma_C%Npoints) stop "Sigma_C has different number of Matsubara points with respect to SigmaFull."
      !
      if(local)then
         !
         if(verbose)write(*,"(A)") "     Join of local attributes."
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nmats,SigmaFull,Sigma_C,Sigma_X),&
         !$OMP PRIVATE(iw,ispin)
         !$OMP DO
         do iw=1,Nmats
            do ispin=1,Nspin
               !
               SigmaFull%ws(:,:,iw,ispin) = Sigma_C%ws(:,:,iw,ispin) + Sigma_X%N_s(:,:,ispin)
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
      else
         !
         if(verbose)write(*,"(A)") "     Join of non-local attributes."
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nmats,Nkpt,SigmaFull,Sigma_C,Sigma_X),&
         !$OMP PRIVATE(iw,ik,ispin)
         !$OMP DO
         do iw=1,Nmats
            do ik=1,Nkpt
               do ispin=1,Nspin
                  !
                  SigmaFull%wks(:,:,iw,ik,ispin) = Sigma_C%wks(:,:,iw,ik,ispin) + Sigma_X%N_ks(:,:,ik,ispin)
                  !
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         call FermionicKsum(SigmaFull)
         !
      endif
      !
   end subroutine join_SigmaCX


   !---------------------------------------------------------------------------!
   !PURPOSE: Replace SigmaImp in SigmaGW at the indexes contained in orbs.
   !         The Hartee contribution computed as N*curlyU is stored in the
   !         SigmaImp%N_s attribute and it is removed during the merge.
   !TEST ON: 27-10-2020
   !---------------------------------------------------------------------------!
   subroutine MergeSelfEnergy(SigmaGW,SigmaGW_DC,SigmaImp,coeff,orbs,DC_type,Rot)
      !
      use parameters
      use utils_misc
      use linalg, only : rotate
      implicit none
      !
      type(FermionicField),intent(inout)    :: SigmaGW
      type(FermionicField),intent(in)       :: SigmaGW_DC
      type(FermionicField),intent(in)       :: SigmaImp
      real(8),intent(in)                    :: coeff
      integer,allocatable,intent(in)        :: orbs(:,:)
      character(len=*),intent(in)           :: DC_type
      logical,intent(in)                    :: Rot
      !
      real(8)                               :: Beta
      integer                               :: iw,ik,isite,iorb,jorb
      integer                               :: ispin,Norb_imp
      integer                               :: i_loc,j_loc
      integer                               :: Nkpt,Norb,Nmats,Nsite
      logical                               :: localDC
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- Merge SelfEnergy"
      !
      !
      ! Check on the input Fields
      if(.not.SigmaGW%status) stop "SigmaGW not properly initialized."
      if(.not.SigmaImp%status) stop "SigmaImp not properly initialized."
      if(SigmaGW%Nkpt.eq.0) stop "SigmaGW k dependent attributes not properly initialized."
      if(SigmaImp%Nkpt.ne.0) stop "SigmaImp k dependent attributes are supposed to be unallocated."
      !
      Norb = SigmaGW%Norb
      Nkpt = SigmaGW%Nkpt
      Beta = SigmaGW%Beta
      Nmats = SigmaGW%Npoints
      Nsite = SigmaGW%Nsite
      !
      if(SigmaImp%Norb.ne.Norb) stop "SigmaImp has different orbital dimension with respect to SigmaGW."
      if(SigmaImp%Beta.ne.Beta) stop "SigmaImp has different Beta with respect to SigmaGW."
      if(SigmaImp%Npoints.ne.Nmats) stop "SigmaImp has different number of Matsubara points with respect to SigmaGW."

      if(size(orbs,dim=1).ne.Nsite) stop "Number of orbital lists does not match the number of sites."
      Norb_imp=0
      do isite=1,Nsite
         do iorb=1,size(orbs(isite,:))
            if(orbs(isite,iorb).ne.0) Norb_imp=Norb_imp+1
         enddo
      enddo
      if(Norb_imp.gt.Norb) stop "Number of orbital to be inserted is bigger than the total lattice orbital space."
      !
      select case(DC_type)
         case default
            stop "Available DC types for the self-energy: Sloc or GlocWloc."
         case("Sloc")
            localDC = .true.
         case("GlocWloc")
            if(.not.SigmaGW_DC%status) stop "SigmaGW_DC not properly initialized."
            if(SigmaGW_DC%Nkpt.ne.0) stop "SigmaGW_DC k dependent attributes are supposed to be unallocated."
            if(SigmaGW_DC%Norb.ne.Norb) stop "SigmaGW_DC has different orbital dimension with respect to SigmaGW."
            if(SigmaGW_DC%Beta.ne.Beta) stop "SigmaGW_DC has different Beta with respect to SigmaGW."
            if(SigmaGW_DC%Npoints.ne.Nmats) stop "SigmaGW_DC has different number of Matsubara points with respect to SigmaGW."
            localDC = .false.
      end select
      !
      !Fill the local attributes so as to fully replace the local GW contibution
      call FermionicKsum(SigmaGW)
      !
      !all sites if(expand.or.AFM) otherwise only one site and the orbitals within orbs
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nsite,Nmats,Nkpt,orbs,coeff,SigmaGW,SigmaGW_DC,SigmaImp,localDC,Rot),&
      !$OMP PRIVATE(isite,ispin,iorb,jorb,i_loc,j_loc,iw,ik)
      !$OMP DO
      do isite=1,Nsite
         do iorb=1,size(orbs(isite,:))
            do jorb=1,size(orbs(isite,:))
               !
               if((.not.Rot).and.(iorb.ne.jorb))cycle
               !
               !SigmaImp and SigmaGW have the same arrangements
               i_loc = orbs(isite,iorb)
               j_loc = orbs(isite,jorb)
               !
               do ispin=1,Nspin
                  do iw=1,Nmats
                     do ik=1,Nkpt
                        !
                        if(localDC)then
                           SigmaGW%wks(i_loc,j_loc,iw,ik,ispin) = SigmaGW%wks(i_loc,j_loc,iw,ik,ispin)              &
                                                                - coeff*SigmaGW%ws(i_loc,j_loc,iw,ispin)            &
                                                                + coeff*(SigmaImp%ws(i_loc,j_loc,iw,ispin)-2d0*SigmaImp%N_s(i_loc,j_loc,ispin))
                        else
                           SigmaGW%wks(i_loc,j_loc,iw,ik,ispin) = SigmaGW%wks(i_loc,j_loc,iw,ik,ispin)              &
                                                                - coeff*SigmaGW_DC%ws(i_loc,j_loc,iw,ispin)         &
                                                                + coeff*(SigmaImp%ws(i_loc,j_loc,iw,ispin)-2d0*SigmaImp%N_s(i_loc,j_loc,ispin))
                        endif
                        !
                     enddo
                  enddo
               enddo
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      !Put SigmaImp in the local attribute
      call FermionicKsum(SigmaGW)
      !
   end subroutine MergeSelfEnergy


   !---------------------------------------------------------------------------!
   !PURPOSE: Replace PiImp in PiGW at the indexes contained in orbs
   !TEST ON: 23-10-2020
   !---------------------------------------------------------------------------!
   subroutine MergePolarization(PiGW,PiImp,coeff,orbs)
      !
      use parameters
      use linalg, only : rotate
      implicit none
      !
      type(BosonicField),intent(inout)      :: PiGW
      type(BosonicField),intent(in)         :: PiImp
      real(8),intent(in)                    :: coeff
      integer,allocatable,intent(in)        :: orbs(:,:)
      !
      real(8)                               :: Beta
      integer                               :: iw,ik,isite
      integer                               :: iorb,jorb,korb,lorb
      integer                               :: Norb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: Nkpt,Norb,Nmats,Nsite
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- Merge Polarization"
      !
      !
      ! Check on the input Fields
      if(.not.PiGW%status) stop "PiGW not properly initialized."
      if(.not.PiImp%status) stop "PiImp not properly initialized."
      if(PiGW%Nkpt.eq.0) stop "PiGW k dependent attributes not properly initialized."
      if(PiImp%Nkpt.ne.0) stop "PiImp k dependent attributes are supposed to be unallocated."
      !
      Norb = int(sqrt(dble(PiGW%Nbp)))
      Nkpt = PiGW%Nkpt
      Beta = PiGW%Beta
      Nmats = PiGW%Npoints
      Nsite = PiGW%Nsite
      !
      if(PiImp%Beta.ne.Beta) stop "PiImp has different Beta with respect to PiGW."
      if(PiImp%Npoints.ne.Nmats) stop "PiImp has different number of Matsubara points with respect to PiGW."
      !
      if(size(orbs,dim=1).ne.Nsite) stop "Number of orbital lists does not match the number of sites."
      Norb_imp=0
      do isite=1,Nsite
         do iorb=1,size(orbs(isite,:))
            if(orbs(isite,iorb).ne.0) Norb_imp=Norb_imp+1
         enddo
      enddo
      if(Norb_imp.gt.Norb) stop "Number of orbital to be inserted is bigger than the total lattice orbital space."
      !
      !Fill the local attributes so as to fully replace the local GW contibution
      call BosonicKsum(PiGW)
      !
      !all sites if(expand) otherwise only one site and the orbitals within orbs
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nsite,Nmats,Nkpt,Norb,orbs,coeff,PiGW,PiImp),&
      !$OMP PRIVATE(isite,iorb,jorb,korb,lorb,iw,ik),&
      !$OMP PRIVATE(i_loc,j_loc,k_loc,l_loc,ib_loc,jb_loc)
      !$OMP DO
      do isite=1,Nsite
         !
         do iorb=1,size(orbs(isite,:))
            do jorb=1,size(orbs(isite,:))
               do korb=1,size(orbs(isite,:))
                  do lorb=1,size(orbs(isite,:))
                     !
                     !PiImp and PiGW
                     i_loc = orbs(isite,iorb)
                     j_loc = orbs(isite,jorb)
                     k_loc = orbs(isite,korb)
                     l_loc = orbs(isite,lorb)
                     !
                     ib_loc = k_loc + Norb*(i_loc-1)
                     jb_loc = l_loc + Norb*(j_loc-1)
                     !
                     if((i_loc.eq.k_loc).and.(l_loc.eq.j_loc))then
                        !
                        do iw=1,Nmats
                           do ik=1,Nkpt
                              !
                              PiGW%screened(ib_loc,jb_loc,iw,ik) = PiGW%screened(ib_loc,jb_loc,iw,ik)                &
                                                                 - coeff*PiGW%screened_local(ib_loc,jb_loc,iw)       &
                                                                 + coeff*PiImp%screened_local(ib_loc,jb_loc,iw)
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
      !Put PiImp in the local attribute
      call BosonicKsum(PiGW)
      !
   end subroutine MergePolarization


end module utils_fields
