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
   interface FieldKsum
      module procedure FermionicKsum
      module procedure BosonicKsum
   end interface FieldKsum

   interface clear_attributes
      module procedure clear_attributes_Fermion
      module procedure clear_attributes_Boson
   end interface clear_attributes

   interface clear_MatrixElements
      module procedure clear_MatrixElements_Matrix
      module procedure clear_MatrixElements_Fermion
      module procedure clear_MatrixElements_Boson
   end interface clear_MatrixElements

   interface loc2imp
      module procedure loc2imp_Matrix
      module procedure loc2imp_Fermionic
      module procedure loc2imp_Bosonic
   end interface loc2imp

   interface imp2loc
      module procedure imp2loc_Matrix
      module procedure imp2loc_Matrix_s
      module procedure imp2loc_Fermionic
      module procedure imp2loc_Bosonic
   end interface imp2loc

   interface Expand2Nsite
      module procedure Expand2Nsite_Matrix
      module procedure Expand2Nsite_Matrix_s
      module procedure Expand2Nsite_Fermionic
      module procedure Expand2Nsite_Bosonic
   end interface Expand2Nsite

   interface symmetrize_GW
      module procedure symmetrize_GW_Matrix_d
      module procedure symmetrize_GW_Matrix_z
      module procedure symmetrize_GW_Fermionic
      module procedure symmetrize_GW_Bosonic
   end interface symmetrize_GW

   interface symmetrize_imp
      module procedure symmetrize_imp_Matrix_d
      module procedure symmetrize_imp_Matrix_z
      module procedure symmetrize_imp_Fermionic
   end interface symmetrize_imp

   interface MergeFields
      module procedure Merge_SelfEnergy
      module procedure Merge_Polarization
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

   interface DeallocateField
      module procedure DeallocateFermionicField
      module procedure DeallocateBosonicField
   end interface DeallocateField

   interface build_Potential
      module procedure build_Potential_Field
      module procedure build_Potential_Mat
   end interface build_Potential

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
   public :: FieldKsum
   public :: FermionicKsum
   public :: BosonicKsum
   public :: AllocateLattice
   public :: DeallocateLattice
   public :: AllocateFermionicField
   public :: AllocateBosonicField
   public :: DeallocateFermionicField
   public :: DeallocateBosonicField
   public :: DeallocateField
   public :: TransformBosonicField
   public :: clear_attributes
   public :: clear_MatrixElements
   public :: duplicate
   public :: isReal
   public :: loc2imp
   public :: imp2loc
   public :: Expand2Nsite
   public :: symmetrize_GW
   public :: symmetrize_imp
   public :: join_SigmaCX
   public :: MergeFields
   public :: build_Potential
   public :: product2NN
   public :: NN2product
   public :: check_QP_poles
   !functions
   public :: calc_Ek
   public :: calc_Ep

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Fill the local attributes of a Field
   !---------------------------------------------------------------------------!
   subroutine FermionicKsum(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      integer                               :: ik,in,ispin
      !
      if(.not.G%status)stop "FermionicKsum: Field not properly initialized."
      if(G%Nkpt.eq.0)stop "FermionicKsum: Field k dependent attributes not properly initialized."
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
   !
   subroutine BosonicKsum(W)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      integer                               :: ik,in
      !
      if(.not.W%status)stop "BosonicKsum: Field not properly initialized."
      if(W%Nkpt.eq.0)stop "BosonicKsum: Field k dependent attributes not properly initialized."
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
   !PURPOSE: This is just a wrapper around tensor_transform defined in linalg
   !         in order to deal with the frequency axis
   !---------------------------------------------------------------------------!
   subroutine TransformBosonicField(W,Map,UL,UR,type)
      !
      use parameters
      use utils_misc
      use linalg, only : tensor_transform
      implicit none
      !
      type(BosonicField),intent(inout)      :: W
      integer,intent(in)                    :: Map(:,:,:)
      complex(8),intent(in)                 :: UL(:,:)
      complex(8),intent(in),optional        :: UR(:,:)
      character(len=2),intent(in),optional  :: type
      !
      integer                               :: iw,Norb
      complex(8),allocatable                :: UR_(:,:)
      character(len=2)                      :: type_
      !
      if(.not.W%status) stop "TransformLocalBosonicField: field not properly initialized."
      call assert_shape(Map,[W%Nbp,W%Nbp,4],"TransformLocalBosonicField","Map")
      Norb = int(sqrt(dble(W%Nbp)))
      call assert_shape(UL,[Norb,Norb],"TransformLocalBosonicField","UL")
      if(present(UR))then
         call assert_shape(UR,[Norb,Norb],"TransformLocalBosonicField","UR")
         UR_ = UR
      else
         UR_ = UL
      endif
      !
      type_="NN"
      if(present(type))type_=reg(type)
      !
      if(allocated(W%bare_local))call tensor_transform(reg(type_),W%bare_local,Map,UL,UR_)
      if(allocated(W%screened_local))then
         do iw=1,W%Npoints
            call tensor_transform(reg(type_),W%screened_local(:,:,iw),Map,UL,UR_)
         enddo
      endif
      deallocate(UR_)
      !
   end subroutine TransformBosonicField


   !---------------------------------------------------------------------------!
   !PURPOSE: Clear the internal attributes of a Fermionic/Bosonic field
   !---------------------------------------------------------------------------!
   subroutine clear_attributes_Fermion(G,thresh)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      real(8),intent(in),optional           :: thresh
      if(.not.G%status) stop "clear_attributes_Fermion: field not properly initialized."
      if(present(thresh))then
         if(allocated(G%N_s)) where(abs(G%N_s)<abs(thresh)) G%N_s=czero
         if(allocated(G%ws)) where(abs(G%ws)<abs(thresh)) G%ws=czero
         if(allocated(G%N_ks)) where(abs(G%N_ks)<abs(thresh)) G%N_ks=czero
         if(allocated(G%wks)) where(abs(G%wks)<abs(thresh)) G%wks=czero
      else
         if(allocated(G%N_s)) G%N_s=czero
         if(allocated(G%ws)) G%ws=czero
         if(allocated(G%N_ks)) G%N_ks=czero
         if(allocated(G%wks)) G%wks=czero
      endif
   end subroutine clear_attributes_Fermion
   !
   subroutine clear_attributes_Boson(W,thresh)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      real(8),intent(in),optional           :: thresh
      if(.not.W%status) stop "clear_attributes_Boson: field not properly initialized."
      if(present(thresh))then
         if(allocated(W%bare_local)) where(abs(W%bare_local)<abs(thresh)) W%bare_local=czero
         if(allocated(W%screened_local)) where(abs(W%screened_local)<abs(thresh)) W%screened_local=czero
         if(allocated(W%bare)) where(abs(W%bare)<abs(thresh)) W%bare=czero
         if(allocated(W%screened)) where(abs(W%screened)<abs(thresh)) W%screened=czero
      else
         if(allocated(W%bare_local)) W%bare_local=czero
         if(allocated(W%screened_local)) W%screened_local=czero
         if(allocated(W%bare)) W%bare=czero
         if(allocated(W%screened)) W%screened=czero
      endif
   end subroutine clear_attributes_Boson


   !---------------------------------------------------------------------------!
   !PURPOSE: Clear Matix elements of a Matrix/Fermionic/Bosonic field attributes
   !         from user-provided Map
   !---------------------------------------------------------------------------!
   subroutine clear_MatrixElements_Matrix(O,Map)
      use parameters
      use utils_misc
      implicit none
      complex(8),intent(inout)              :: O(:,:)
      logical,intent(in)                    :: Map(:,:)
      integer                               :: iorb,jorb
      logical                               :: erase
      !
      call assert_shape(Map,[size(O,dim=1),size(O,dim=2)],"clear_MatrixElements_Matrix","Map")
      !
      do iorb=1,size(O,dim=1)
         do jorb=1,size(O,dim=2)
            erase = Map(iorb,jorb)
            if(erase) O(iorb,jorb)=czero
         enddo
      enddo
      !
   end subroutine clear_MatrixElements_Matrix
   !
   subroutine clear_MatrixElements_Fermion(G,Map,LocalOnly)
      use parameters
      use utils_misc
      implicit none
      type(FermionicField),intent(inout)    :: G
      logical,intent(in)                    :: Map(:,:)
      logical,intent(in),optional           :: LocalOnly
      integer                               :: ispin,iorb,jorb
      logical                               :: erase,LocalOnly_
      !
      if(.not.G%status) stop "clear_MatrixElements_Fermion: field not properly initialized."
      call assert_shape(Map,[G%Norb,G%Norb],"clear_MatrixElements_Fermion","Map")
      !
      LocalOnly_=.true.
      if(present(LocalOnly))LocalOnly_=LocalOnly
      !
      do ispin=1,Nspin
         do iorb=1,G%Norb
            do jorb=1,G%Norb
               erase = Map(iorb,jorb)
               if(allocated(G%N_s).and.erase) G%N_s(iorb,jorb,ispin)=czero
               if(allocated(G%ws).and.erase) G%ws(iorb,jorb,:,ispin)=czero
               if(.not.LocalOnly_)then
                  if(allocated(G%N_ks).and.erase) G%N_ks(iorb,jorb,:,ispin)=czero
                  if(allocated(G%wks).and.erase) G%wks(iorb,jorb,:,:,ispin)=czero
               endif
            enddo
         enddo
      enddo
      !
   end subroutine clear_MatrixElements_Fermion
   !
   subroutine clear_MatrixElements_Boson(W,Map,LocalOnly)
      use parameters
      use utils_misc
      implicit none
      type(BosonicField),intent(inout)      :: W
      logical,intent(in)                    :: Map(:,:)
      logical,intent(in),optional           :: LocalOnly
      integer                               :: ib1,ib2
      integer                               :: iorb,jorb,Norb
      logical                               :: prod,erase,LocalOnly_
      !
      if(.not.W%status) stop "clear_MatrixElements_Boson: field not properly initialized."
      Norb = int(sqrt(dble(W%Nbp)))
      if((size(Map,dim=1).eq.W%Nbp) .and. (size(Map,dim=2).eq.W%Nbp))then
         prod = .true.
      elseif((size(Map,dim=1).eq.Norb) .and. (size(Map,dim=2).eq.Norb))then
         prod = .false.
      else
         stop "clear_MatrixElements_Boson: wrong Map dimensions."
      endif
      !
      LocalOnly_=.true.
      if(present(LocalOnly))LocalOnly_=LocalOnly
      !
      if(prod)then
         !this deletes all the possible elements
         do ib1=1,W%Nbp
            do ib2=1,W%Nbp
               erase = Map(ib1,ib2)
               if(allocated(W%bare_local).and.erase) W%bare_local(ib1,ib2)=czero
               if(allocated(W%screened_local).and.erase) W%screened_local(ib1,ib2,:)=czero
               if(.not.LocalOnly_)then
                  if(allocated(W%bare).and.erase) W%bare(ib1,ib2,:)=czero
                  if(allocated(W%screened).and.erase) W%screened(ib1,ib2,:,:)=czero
               endif
            enddo
         enddo
      else
         !this deletes only density density terms
         do iorb=1,Norb
            do jorb=1,Norb
               erase = Map(iorb,jorb)
               call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
               if(allocated(W%bare_local).and.erase) W%bare_local(ib1,ib2)=czero
               if(allocated(W%screened_local).and.erase) W%screened_local(ib1,ib2,:)=czero
               if(.not.LocalOnly_)then
                  if(allocated(W%bare).and.erase) W%bare(ib1,ib2,:)=czero
                  if(allocated(W%screened).and.erase) W%screened(ib1,ib2,:,:)=czero
               endif
            enddo
         enddo
      endif
      !
   end subroutine clear_MatrixElements_Boson


   !---------------------------------------------------------------------------!
   !PURPOSE: Duplicate a Fermionic/Bosonic field
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
      if(allocated(Gold%ws)) Gnew%ws = Gold%ws
      if(allocated(Gold%N_ks)) Gnew%N_ks= Gold%N_ks
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
   !---------------------------------------------------------------------------!
   subroutine isReal_Matrix(O)
      use parameters
      implicit none
      complex(8),intent(inout)              :: O(:,:)
      O=dcmplx(dreal(O),0d0)
   end subroutine isReal_Matrix
   !
   subroutine isReal_Fermionic(G)
      use parameters
      implicit none
      type(FermionicField),intent(inout)    :: G
      if(allocated(G%N_s)) G%N_s = dcmplx(dreal(G%N_s),0d0)
      if(allocated(G%ws)) G%ws = dcmplx(dreal(G%ws),0d0)
      if(allocated(G%N_ks)) G%N_ks = dcmplx(dreal(G%N_ks),0d0)
      if(allocated(G%wks)) G%wks = dcmplx(dreal(G%wks),0d0)
   end subroutine isReal_Fermionic
   !
   subroutine isReal_Bosonic(W,LocalOnly)
      use parameters
      implicit none
      type(BosonicField),intent(inout)      :: W
      integer,intent(in),optional           :: LocalOnly
      integer                               :: LocalOnly_
      LocalOnly_=.true.
      if(present(LocalOnly))LocalOnly_=LocalOnly
      if(allocated(W%bare_local)) W%bare_local = dcmplx(dreal(W%bare_local),0d0)
      if(allocated(W%screened_local)) W%screened_local = dcmplx(dreal(W%screened_local),0d0)
      if(allocated(W%bare).and.(.not.LocalOnly_)) W%bare = dcmplx(dreal(W%bare),0d0)
      if(allocated(W%screened).and.(.not.LocalOnly_)) W%screened = dcmplx(dreal(W%screened),0d0)
   end subroutine isReal_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE: Extract from a lattice local projection the subset of orbital
   !         indexes of corresponging to a given site
   !---------------------------------------------------------------------------!
   subroutine loc2imp_Fermionic(Gimp,Gloc,orbs,sitename,U)
      !
      use parameters
      use utils_misc
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
      if(verbose)write(*,"(A)") "---- loc2imp_Fermionic"
      !
      !
      if(present(sitename)) write(*,"(A)") "     Extraction of Gloc for site: "//trim(sitename)
      if(.not.Gloc%status) stop "loc2imp_Fermionic: Gloc not properly initialized."
      if(.not.Gimp%status) stop "loc2imp_Fermionic: Gimp not properly initialized."
      if(Gloc%Norb.eq.0) stop "loc2imp_Fermionic: Norb of Gloc not defined."
      if(Gimp%Norb.eq.0) stop "loc2imp_Fermionic: Norb of Gimp not defined."
      if(Gloc%Npoints.eq.0) stop "loc2imp_Fermionic: Npoints of Gloc not defined."
      if(Gimp%Npoints.eq.0) stop "loc2imp_Fermionic: Npoints of Gimp not defined."
      if(Gimp%Beta.ne.Gloc%Beta) stop "loc2imp_Fermionic: Gimp and Gloc have different beta."
      if(Gimp%Npoints.ne.Gloc%Npoints) stop "loc2imp_Fermionic: Gimp and Gloc have different number of Matsubara points."
      if(Gimp%Nkpt.ne.0) stop "loc2imp_Fermionic: Gimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Gloc%ws)) stop "loc2imp_Fermionic: Gloc local projection not allocated."
      if(.not.allocated(Gimp%ws)) stop "loc2imp_Fermionic: Gimp local projection not allocated."
      if(size(orbs).ne.Gimp%Norb) stop "loc2imp_Fermionic: can't fit the requested orbitals inside Gimp."
      if(size(orbs).gt.Gloc%Norb) stop "loc2imp_Fermionic: number of requested orbitals greater than Gloc size."
      if(present(U))then
         call assert_shape(U,[size(orbs),size(orbs)],"loc2imp_Fermionic","U")
         !if(size(U,dim=1).ne.3) write(*,"(A)") "     Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
         write(*,"(A)") "     The local orbital space will be rotated during extraction."
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
      use utils_misc
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
      if(verbose)write(*,"(A)") "---- loc2imp_Matrix"
      !
      !
      if(present(sitename)) write(*,"(A)") "     Extraction of Operator for site: "//trim(sitename)
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "loc2imp_Matrix: Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "loc2imp_Matrix: Oimp not square."
      if(size(orbs).ne.size(Oimp,dim=1)) stop "loc2imp_Matrix: can't fit the requested orbitals inside Oimp."
      if(size(orbs).gt.size(Oloc,dim=1)) stop "loc2imp_Matrix: number of requested orbitals greater than Oloc size."
      if(present(U))then
         call assert_shape(U,[size(orbs),size(orbs)],"loc2imp_Matrix","U")
         !if(size(U,dim=1).ne.3) write(*,"(A)") "     Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
         write(*,"(A)") "     The local orbital space will be rotated during extraction."
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
   subroutine loc2imp_Bosonic(Wimp,Wloc,orbs,sitename,U,Map,type)
      !
      use parameters
      use utils_misc
      use linalg, only : tensor_transform
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wimp
      type(BosonicField),intent(in)         :: Wloc
      integer,allocatable,intent(in)        :: orbs(:)
      character(len=*),intent(in),optional  :: sitename
      complex(8),allocatable,optional       :: U(:,:)
      integer,allocatable,optional          :: Map(:,:,:)
      character(len=2),intent(in),optional  :: type
      !
      integer                               :: ip,ip_loc
      integer                               :: Norb_imp,Norb_loc
      integer                               :: ib_imp,jb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: i_imp,j_imp,k_imp,l_imp
      logical                               :: doBare,rot,WlocStatic
      character(len=2)                      :: type_
      !
      !
      if(verbose)write(*,"(A)") "---- loc2imp_Bosonic"
      !
      !
      if(present(sitename)) write(*,"(A)") "     Extraction of Wloc for site: "//trim(sitename)
      if(.not.Wloc%status) stop "loc2imp_Bosonic: Wloc not properly initialized."
      if(.not.Wimp%status) stop "loc2imp_Bosonic: Wimp not properly initialized."
      if(Wloc%Nbp.eq.0) stop "loc2imp_Bosonic: Norb of Wloc not defined."
      if(Wimp%Nbp.eq.0) stop "loc2imp_Bosonic: Norb of Wimp not defined."
      if(Wloc%Npoints.eq.0) stop "loc2imp_Bosonic: Npoints of Wloc not defined."
      if(Wimp%Npoints.eq.0) stop "loc2imp_Bosonic: Npoints of Wimp not defined."
      if(Wimp%Beta.ne.Wloc%Beta) stop "loc2imp_Bosonic: Wimp and Wloc have different beta."
      if(Wimp%Nkpt.ne.0) stop "loc2imp_Bosonic: Wimp k-dependent attributes attributes are supposed to be unallocated."
      if((Wloc%Npoints.ne.1).and.(Wimp%Npoints.ne.Wloc%Npoints))stop "loc2imp_Bosonic: Wimp and Wloc have different number of Matsubara points."
      if(.not.allocated(Wloc%screened_local)) stop "loc2imp_Bosonic: Wloc screened_local attribute not allocated."
      if(.not.allocated(Wimp%screened_local)) stop "loc2imp_Bosonic: Wimp screened_local attribute not allocated."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      doBare = allocated(Wimp%bare_local)
      WlocStatic = Wloc%Npoints .eq. 1
      if(.not.allocated(Wloc%bare_local).and.doBare) stop "loc2imp_Bosonic: Wloc bare_local attribute not allocated."
      !
      if(size(orbs).ne.Norb_imp) stop "loc2imp_Bosonic: can't fit the requested orbitals inside Wimp."
      if(size(orbs).gt.Norb_loc) stop "loc2imp_Bosonic: number of requested orbitals greater than Wloc size."
      if(present(U))then
         type_="NN"
         if(present(type))type_=reg(type)
         if(.not.present(Map)) stop "loc2imp_Bosonic: requested orbital rotation but map not provided."
         call assert_shape(U,[size(orbs),size(orbs)],"loc2imp_Bosonic","U")
         !if(size(U,dim=1).ne.3) write(*,"(A)") "     Warning: The local orbital space rotation is well defined only for a t2g sub-shell."
         write(*,"(A)") "     The local orbital space will be rotated during extraction."
      endif
      if(present(Map))then
         if(.not.present(U)) stop "loc2imp_Bosonic: Also the rotation must be provided."
         call assert_shape(Map,[Wimp%Nbp,Wimp%Nbp,4],"loc2imp_Bosonic","Map")
      endif
      rot=.false.
      if(present(U).and.present(Map))rot=.true.
      !
      call clear_attributes(Wimp)
      !
      do i_imp=1,Norb_imp
         do j_imp=1,Norb_imp
            do k_imp=1,Norb_imp
               do l_imp=1,Norb_imp
                  !
                  ! bosonic indexes of the impurity
                  ib_imp = j_imp + Norb_imp*(i_imp-1)
                  jb_imp = l_imp + Norb_imp*(k_imp-1)
                  !
                  ! mapping
                  i_loc = orbs(i_imp)
                  j_loc = orbs(j_imp)
                  k_loc = orbs(k_imp)
                  l_loc = orbs(l_imp)
                  !
                  ! corresponding bosonic indexes on the lattice
                  call F2Bindex(Norb_loc,[i_loc,j_loc],[k_loc,l_loc],ib_loc,jb_loc)
                  !
                  if(doBare)Wimp%bare_local(ib_imp,jb_imp) =  Wloc%bare_local(ib_loc,jb_loc)
                  do ip=1,Wimp%Npoints
                     ip_loc = ip
                     if(WlocStatic) ip_loc=1
                     Wimp%screened_local(ib_imp,jb_imp,ip) = Wloc%screened_local(ib_loc,jb_loc,ip_loc)
                  enddo
                  !
               enddo
            enddo
         enddo
      enddo
      !
      if(rot)then
         if(doBare) call tensor_transform(reg(type_),Wimp%bare_local,Map,U,U)
         do ip=1,Wimp%Npoints
            call tensor_transform(reg(type_),Wimp%screened_local(:,:,ip),Map,U,U)
         enddo
      endif
      !
   end subroutine loc2imp_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE: Expand the the subset of orbital indexes of corresponging to a
   !         given site to the full lattice quantity
   !---------------------------------------------------------------------------!
   subroutine imp2loc_Fermionic(Gloc,Gimp,impndx,LocalOrbs,expand,AFM,rot,name)
      !
      use parameters
      use utils_misc
      use linalg, only : rotate
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gloc
      type(FermionicField),intent(in)       :: Gimp
      integer,intent(in)                    :: impndx
      type(LocalOrbitals),allocatable,intent(in):: LocalOrbs(:)
      logical,intent(in)                    :: expand
      logical,intent(in)                    :: AFM
      logical,intent(in)                    :: rot
      character(len=*),intent(in),optional  :: name
      !
      complex(8),allocatable                :: Gtmp(:,:,:,:),Ntmp(:,:,:)
      integer                               :: ip,ispin,ispin_imp,isite,Nsite
      integer                               :: Norb_loc,Norb_imp,sitendx
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      if(verbose)write(*,"(A)") "---- imp2loc_Fermionic"
      !
      !
      if(expand.and.AFM) stop "imp2loc_Fermionic: AFM condition and expansion to real space not implemented."
      if(.not.Gloc%status) stop "imp2loc_Fermionic: Gloc not properly initialized."
      if(.not.Gimp%status) stop "imp2loc_Fermionic: Gimp not properly initialized."
      if(Gloc%Norb.eq.0) stop "imp2loc_Fermionic: Norb of Gloc not defined."
      if(Gimp%Norb.eq.0) stop "imp2loc_Fermionic: Norb of Gimp not defined."
      if(Gloc%Npoints.eq.0) stop "imp2loc_Fermionic: Npoints of Gloc not defined."
      if(Gimp%Npoints.eq.0) stop "imp2loc_Fermionic: Npoints of Gimp not defined."
      if(Gloc%Beta.ne.Gimp%Beta) stop "imp2loc_Fermionic: Gimp and Gloc have different beta."
      if(Gloc%Npoints.ne.Gimp%Npoints) stop "imp2loc_Fermionic: Gimp and Gloc have different number of Matsubara points."
      if(Gimp%Nkpt.ne.0) stop "imp2loc_Fermionic: Gimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Gloc%ws)) stop "imp2loc_Fermionic: Gloc local projection not allocated."
      if(.not.allocated(Gimp%ws)) stop "imp2loc_Fermionic: Gimp local projection not allocated."
      if(.not.allocated(LocalOrbs)) stop "imp2loc_Fermionic: LocalOrbs not properly initialized."
      !
      Norb_loc = Gloc%Norb
      Norb_imp = Gimp%Norb
      !
      if(expand)then
         if(mod(Norb_loc,LocalOrbs(1)%Norb).ne.0)stop "imp2loc_Fermionic: Number of requested orbitals is not a commensurate subset of the lattice field."
         if(LocalOrbs(1)%Norb.ne.Norb_imp)stop "imp2loc_Fermionic: orbital dimension of impurity field does not match the LocalOrbs list."
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         else
            write(*,"(A)") "     The impurity field of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         endif
         Nsite = Gloc%Nsite
      elseif(AFM)then
         if(Gloc%Nsite.ne.2) stop "imp2loc_Fermionic: AFM is implemented only for a two site lattice."
         if(Norb_loc/Norb_imp.ne.2) stop "imp2loc_Fermionic: Lattice indexes are not twice the impurity ones."
         if(present(name))then
            write(*,"(2(A,"//str(Norb_imp)//"I3))") "     "//reg(name)//" of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs," and spin-flipped into: ",LocalOrbs(impndx)%Orbs+Norb_imp
         else
            write(*,"(2(A,"//str(Norb_imp)//"I3))") "     The impurity field of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs," and spin-flipped into: ",LocalOrbs(impndx)%Orbs+Norb_imp
         endif
         Nsite = 2
      else
         if(present(name))then
            write(*,"(A,15I3)") "     "//reg(name)//" of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         else
            write(*,"(A,15I3)") "     The impurity field of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         endif
         Nsite = 1
      endif
      !
      if(Nsite.gt.size(LocalOrbs)) stop "imp2loc_Fermionic: number of target sites is larger LocalOrbs list."
      if(impndx.gt.size(LocalOrbs)) stop "imp2loc_Fermionic: target site is outside LocalOrbs list."
      !
      if(rot)then
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         else
            write(*,"(A)") "     The impurity orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         endif
      endif
      !
      ! Rotating either one site or all of them depending on expand
      allocate(Ntmp(Norb_imp,Norb_imp,Nspin))
      allocate(Gtmp(Norb_imp,Norb_imp,Gimp%Npoints,Nspin))
      do isite=1,Nsite
         !
         sitendx = impndx
         if(expand.or.AFM) sitendx = isite
         !
         if(LocalOrbs(sitendx)%Norb.ne.Norb_imp) stop "imp2loc_Fermionic: can't fit the requested orbitals from Gimp."
         if(LocalOrbs(sitendx)%Norb.gt.Norb_loc) stop "imp2loc_Fermionic: number of requested orbitals greater than Gloc size."
         !
         Ntmp=czero;Gtmp=czero
         !
         if(rot)then
            !
            if(.not.allocated(LocalOrbs(sitendx)%RotDag)) stop "imp2loc_Fermionic: One of the rotations is not allocated."
            call assert_shape(LocalOrbs(sitendx)%RotDag,[Norb_imp,Norb_imp],"imp2loc_Fermionic","LocalOrbs("//str(sitendx)//")%RotDag")
            !
            do ispin=1,Nspin
               Ntmp(:,:,ispin) = rotate(Gimp%N_s(:,:,ispin),LocalOrbs(sitendx)%RotDag)
               do ip=1,Gimp%Npoints
                  Gtmp(:,:,ip,ispin) = rotate(Gimp%ws(:,:,ip,ispin),LocalOrbs(sitendx)%RotDag)
               enddo
            enddo
            !
         else
            !
            Ntmp = Gimp%N_s
            Gtmp = Gimp%ws
            !
         endif
         !
         !Insert into local indexes
         do i_imp=1,LocalOrbs(sitendx)%Norb
            do j_imp=1,LocalOrbs(sitendx)%Norb
               !
               i_loc = LocalOrbs(sitendx)%Orbs(i_imp)
               j_loc = LocalOrbs(sitendx)%Orbs(j_imp)
               !
               do ispin=1,Nspin
                  !
                  ispin_imp=ispin
                  if(AFM.and.(isite.eq.2)) ispin_imp=int(Nspin/ispin)
                  !
                  Gloc%N_s(i_loc,j_loc,ispin) = Ntmp(i_imp,j_imp,ispin_imp)
                  do ip=1,Gimp%Npoints
                     Gloc%ws(i_loc,j_loc,ip,ispin) = Gtmp(i_imp,j_imp,ip,ispin_imp)
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
   subroutine imp2loc_Matrix(Oloc,Oimp,impndx,LocalOrbs,expand,rot,name)
      !
      use parameters
      use utils_misc
      use linalg, only : rotate
      implicit none
      !
      complex(8),intent(inout)              :: Oloc(:,:)
      complex(8),intent(in)                 :: Oimp(:,:)
      integer,intent(in)                    :: impndx
      type(LocalOrbitals),allocatable,intent(in):: LocalOrbs(:)
      logical,intent(in)                    :: expand
      logical,intent(in)                    :: rot
      character(len=*),intent(in),optional  :: name
      !
      complex(8),allocatable                :: Otmp(:,:)
      integer                               :: isite,Nsite
      integer                               :: Norb_loc,Norb_imp,sitendx
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      if(verbose)write(*,"(A)") "---- imp2loc_Matrix"
      !
      !
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "imp2loc_Matrix: Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "imp2loc_Matrix: Oimp not square."
      if(.not.allocated(LocalOrbs)) stop "imp2loc_Matrix: LocalOrbs not properly initialized."
      !
      Norb_loc = size(Oloc,dim=1)
      Norb_imp = size(Oimp,dim=1)
      !
      if(expand)then
         if(mod(Norb_loc,LocalOrbs(1)%Norb).ne.0)stop "imp2loc_Matrix: Number of requested orbitals is not a commensurate subset of the lattice observable."
         if(LocalOrbs(1)%Norb.ne.Norb_imp)stop "imp2loc_Matrix: orbital dimension of impurity matrix does not match the LocalOrbs list."
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         else
            write(*,"(A)") "     The matrix of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         endif
         Nsite = Norb_loc/Norb_imp
      else
         if(present(name))then
            write(*,"(A,15I3)") "     "//reg(name)//" of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         else
            write(*,"(A,15I3)") "     The matrix of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         endif
         Nsite = 1
      endif
      !
      if(Nsite.gt.size(LocalOrbs)) stop "imp2loc_Matrix: number of target sites is larger LocalOrbs list."
      if(impndx.gt.size(LocalOrbs)) stop "imp2loc_Matrix: target site is outside LocalOrbs list."
      !
      if(rot)then
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         else
            write(*,"(A)") "     The matrix orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         endif
      endif
      !
      ! Rotating either one site or all of them depending on expand
      allocate(Otmp(Norb_imp,Norb_imp))
      do isite=1,Nsite
         !
         sitendx = impndx
         if(expand) sitendx = isite
         !
         if(LocalOrbs(sitendx)%Norb.ne.Norb_imp) stop "imp2loc_Matrix: can't fit the requested orbitals from Gimp."
         if(LocalOrbs(sitendx)%Norb.gt.Norb_loc) stop "imp2loc_Matrix: number of requested orbitals greater than Gloc size."
         !
         Otmp=czero
         !
         if(rot)then
            !
            if(.not.allocated(LocalOrbs(sitendx)%RotDag)) stop "imp2loc_Matrix: One of the rotations is not allocated."
            call assert_shape(LocalOrbs(sitendx)%RotDag,[Norb_imp,Norb_imp],"imp2loc_Matrix","LocalOrbs("//str(sitendx)//")%RotDag")
            !
            Otmp = rotate(Oimp,LocalOrbs(sitendx)%RotDag)
            !
         else
            !
            Otmp = Oimp
            !
         endif
         !
         !Insert into local indexes
         do i_imp=1,LocalOrbs(sitendx)%Norb
            do j_imp=1,LocalOrbs(sitendx)%Norb
               !
               i_loc = LocalOrbs(sitendx)%Orbs(i_imp)
               j_loc = LocalOrbs(sitendx)%Orbs(j_imp)
               !
               Oloc(i_loc,j_loc) = Otmp(i_imp,j_imp)
               !
            enddo
         enddo
         !
      enddo
      deallocate(Otmp)
      !
   end subroutine imp2loc_Matrix
   !
   subroutine imp2loc_Matrix_s(Oloc,Oimp,impndx,LocalOrbs,expand,AFM,rot,name)
      !
      use parameters
      use utils_misc
      use linalg, only : rotate
      implicit none
      !
      complex(8),intent(inout)              :: Oloc(:,:,:)
      complex(8),intent(in)                 :: Oimp(:,:,:)
      integer,intent(in)                    :: impndx
      type(LocalOrbitals),allocatable,intent(in):: LocalOrbs(:)
      logical,intent(in)                    :: expand
      logical,intent(in)                    :: AFM
      logical,intent(in)                    :: rot
      character(len=*),intent(in),optional  :: name
      !
      complex(8),allocatable                :: Otmp(:,:,:)
      integer                               :: ispin,ispin_imp,isite,Nsite
      integer                               :: Norb_loc,Norb_imp,sitendx
      integer                               :: i_loc,j_loc
      integer                               :: i_imp,j_imp
      !
      !
      if(verbose)write(*,"(A)") "---- imp2loc_Matrix_s"
      !
      !
      if(expand.and.AFM) stop "imp2loc_Matrix_s: AFM condition and expansion to real space not implemented."
      if(size(Oloc,dim=1).ne.size(Oloc,dim=2)) stop "imp2loc_Matrix_s: Oloc not square."
      if(size(Oimp,dim=1).ne.size(Oimp,dim=2)) stop "imp2loc_Matrix_s: Oimp not square."
      if(size(Oimp,dim=3).ne.Nspin) stop "imp2loc_Matrix_s: Oimp third dimension is not Nspin."
      if(size(Oloc,dim=3).ne.Nspin) stop "imp2loc_Matrix_s: Oloc third dimension is not Nspin."
      if(.not.allocated(LocalOrbs)) stop "imp2loc_Matrix_s: LocalOrbs not properly initialized."
      !
      Norb_loc = size(Oloc,dim=1)
      Norb_imp = size(Oimp,dim=1)
      !
      if(expand)then
         if(mod(Norb_loc,LocalOrbs(1)%Norb).ne.0)stop "imp2loc_Matrix_s: Number of requested orbitals is not a commensurate subset of the lattice observable."
         if(LocalOrbs(1)%Norb.ne.Norb_imp)stop "imp2loc_Matrix_s: orbital dimension of impurity matrix does not match the LocalOrbs list."
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         else
            write(*,"(A)") "     The matrix of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         endif
         Nsite = Norb_loc/Norb_imp
      elseif(AFM)then
         if(Norb_loc/Norb_imp.ne.2) stop "imp2loc_Matrix_s: Lattice indexes are not twice the impurity ones."
         if(present(name))then
            write(*,"(2(A,"//str(Norb_imp)//"I3))") "     "//reg(name)//" of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs," and spin-flipped into: ",LocalOrbs(impndx)%Orbs+Norb_imp
         else
            write(*,"(2(A,"//str(Norb_imp)//"I3))") "     The matrix of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs," and spin-flipped into: ",LocalOrbs(impndx)%Orbs+Norb_imp
         endif
         Nsite = 2
      else
         if(present(name))then
            write(*,"(A,15I3)") "     "//reg(name)//" of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         else
            write(*,"(A,15I3)") "     The matrix of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         endif
         Nsite = 1
      endif
      !
      if(Nsite.gt.size(LocalOrbs)) stop "imp2loc_Matrix_s: number of target sites is larger LocalOrbs list."
      if(impndx.gt.size(LocalOrbs)) stop "imp2loc_Matrix_s: target site is outside LocalOrbs list."
      !
      if(rot)then
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         else
            write(*,"(A)") "     The matrix orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         endif
      endif
      !
      !
      ! Rotating either one site or all of them depending on expand
      allocate(Otmp(Norb_imp,Norb_imp,Nspin))
      do isite=1,Nsite
         !
         sitendx = impndx
         if(expand) sitendx = isite
         !
         if(LocalOrbs(sitendx)%Norb.ne.Norb_imp) stop "imp2loc_Matrix_s: can't fit the requested orbitals from Gimp."
         if(LocalOrbs(sitendx)%Norb.gt.Norb_loc) stop "imp2loc_Matrix_s: number of requested orbitals greater than Gloc size."
         !
         Otmp=czero
         !
         if(rot)then
            !
            if(.not.allocated(LocalOrbs(sitendx)%RotDag)) stop "imp2loc_Matrix_s: One of the rotations is not allocated."
            call assert_shape(LocalOrbs(sitendx)%RotDag,[Norb_imp,Norb_imp],"imp2loc_Matrix_s","LocalOrbs("//str(sitendx)//")%RotDag")
            !
            do ispin=1,Nspin
               Otmp(:,:,ispin) = rotate(Oimp(:,:,ispin),LocalOrbs(sitendx)%RotDag)
            enddo
            !
         else
            !
            Otmp = Oimp
            !
         endif
         !
         !Insert into local indexes
         do i_imp=1,LocalOrbs(sitendx)%Norb
            do j_imp=1,LocalOrbs(sitendx)%Norb
               !
               i_loc = LocalOrbs(sitendx)%Orbs(i_imp)
               j_loc = LocalOrbs(sitendx)%Orbs(j_imp)
               !
               do ispin=1,Nspin
                  !
                  ispin_imp=ispin
                  if(AFM.and.(isite.eq.2)) ispin_imp=int(Nspin/ispin)
                  !
                  Oloc(i_loc,j_loc,ispin) = Otmp(i_imp,j_imp,ispin_imp)
                  !
               enddo
               !
            enddo
         enddo
         !
      enddo
      deallocate(Otmp)
      !
   end subroutine imp2loc_Matrix_s
   !
   subroutine imp2loc_Bosonic(Wloc,Wimp,impndx,LocalOrbs,expand,AFM,rot,name,type)
      !
      use parameters
      use utils_misc
      use linalg, only : tensor_transform
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wloc
      type(BosonicField),intent(in)         :: Wimp
      integer,intent(in)                    :: impndx
      type(LocalOrbitals),allocatable,intent(in):: LocalOrbs(:)
      logical,intent(in)                    :: expand
      logical,intent(in)                    :: AFM
      logical,intent(in)                    :: rot
      character(len=*),intent(in),optional  :: name
      character(len=2),intent(in),optional  :: type
      !
      complex(8),allocatable                :: Wbtmp(:,:),Wstmp(:,:,:)
      type(physicalU)                       :: PhysicalUelements
      integer                               :: ip,isite,Nsite
      integer                               :: Norb_imp,Norb_loc,sitendx
      integer                               :: ib_imp,jb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: i_imp,j_imp,k_imp,l_imp
      logical                               :: doBare
      character(len=2)                      :: type_
      !
      !
      if(verbose)write(*,"(A)") "---- imp2loc_Bosonic"
      !
      !
      if(expand.and.AFM) stop "imp2loc_Bosonic: AFM condition and expansion to real space not implemented."
      if(.not.Wloc%status) stop "imp2loc_Bosonic: Wloc not properly initialized."
      if(.not.Wimp%status) stop "imp2loc_Bosonic: Wimp not properly initialized."
      if(Wloc%Nbp.eq.0) stop "imp2loc_Bosonic: Norb of Wloc not defined."
      if(Wimp%Nbp.eq.0) stop "imp2loc_Bosonic: Norb of Wimp not defined."
      if(Wloc%Npoints.eq.0) stop "imp2loc_Bosonic: Npoints of Wloc not defined."
      if(Wimp%Npoints.eq.0) stop "imp2loc_Bosonic: Npoints of Wimp not defined."
      if(Wloc%Beta.ne.Wimp%Beta) stop "imp2loc_Bosonic: Wimp and Wloc have different beta."
      if(Wloc%Npoints.ne.Wimp%Npoints) stop "imp2loc_Bosonic: Wimp and Wloc have different number of Matsubara points."
      if(Wimp%Nkpt.ne.0) stop "imp2loc_Bosonic: Wimp k-dependent attributes attributes are supposed to be unallocated."
      if(.not.allocated(Wloc%screened_local)) stop "imp2loc_Bosonic: Wloc screened_local attribute not allocated."
      if(.not.allocated(Wimp%screened_local)) stop "imp2loc_Bosonic: Wimp screened_local attribute not allocated."
      if(.not.allocated(LocalOrbs)) stop "imp2loc_Bosonic: LocalOrbs not properly initialized."
      !
      Norb_imp = int(sqrt(dble(Wimp%Nbp)))
      Norb_loc = int(sqrt(dble(Wloc%Nbp)))
      doBare = allocated(Wimp%bare_local)
      !
      if(expand)then
         if(mod(Norb_loc,LocalOrbs(1)%Norb).ne.0)stop "imp2loc_Bosonic: Number of requested orbitals is not a commensurate subset of Wloc."
         if(LocalOrbs(1)%Norb.ne.Norb_imp)stop "imp2loc_Bosonic: orbital dimension of impurity field does not match the LocalOrbs list."
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         else
            write(*,"(A)") "     The impurity field of site "//str(impndx)//" will be expanded to match the lattice orbital space."
         endif
         Nsite = Wloc%Nsite
      elseif(AFM)then
         if(Wloc%Nsite.ne.2) stop "imp2loc_Bosonic: AFM is implemented only for a two site lattice."
         if(Norb_loc/Norb_imp.ne.2) stop "imp2loc_Bosonic: Lattice indexes are not twice the impurity ones."
         if(present(name))then
            write(*,"(2(A,"//str(Norb_imp)//"I3))") "     "//reg(name)//" of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs," and into: ",LocalOrbs(impndx)%Orbs+Norb_imp
         else
            write(*,"(2(A,"//str(Norb_imp)//"I3))") "     The impurity field of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs," and into: ",LocalOrbs(impndx)%Orbs+Norb_imp
         endif
         Nsite = 2
      else
         if(present(name))then
            write(*,"(A,15I3)") "     "//reg(name)//" of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         else
            write(*,"(A,15I3)") "     The impurity field of site "//str(impndx)//" will be inserted into the lattice orbital indexes: ",LocalOrbs(impndx)%Orbs
         endif
         Nsite = 1
      endif
      !
      if(Nsite.gt.size(LocalOrbs)) stop "imp2loc_Bosonic: number of target sites is larger LocalOrbs list."
      if(impndx.gt.size(LocalOrbs)) stop "imp2loc_Bosonic: target site is outside LocalOrbs list."
      !
      if(rot)then
         !
         if(present(name))then
            write(*,"(A)") "     "//reg(name)//" orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         else
            write(*,"(A)") "     The impurity orbital space will be rotated during insertion in "//str(Nsite)//" sites."
         endif
         !
         type_="NN"
         if(present(type))type_=reg(type)
         !
      endif
      !
      ! Rotating either one site or all of them depending on expand.
      ! I'm allocating another one because Wimp has intent(in)
      if(doBare) allocate(Wbtmp(Wimp%Nbp,Wimp%Nbp))
      allocate(Wstmp(Wimp%Nbp,Wimp%Nbp,Wimp%Npoints))
      do isite=1,Nsite
         !
         sitendx = impndx
         if(expand.or.AFM) sitendx = isite
         !
         if(LocalOrbs(sitendx)%Norb.ne.Norb_imp) stop "imp2loc_Bosonic: can't fit the requested orbitals from Gimp."
         if(LocalOrbs(sitendx)%Norb.gt.Norb_loc) stop "imp2loc_Bosonic: number of requested orbitals greater than Gloc size."
         !
         if(doBare) Wbtmp = Wimp%bare_local
         Wstmp = Wimp%screened_local
         !
         if(rot)then
            !
            call init_Uelements(LocalOrbs(sitendx)%Norb,PhysicalUelements)
            !
            if(.not.allocated(LocalOrbs(sitendx)%RotDag)) stop "imp2loc_Fermionic: One of the rotations is not allocated."
            call assert_shape(LocalOrbs(sitendx)%RotDag,[Norb_imp,Norb_imp],"imp2loc_Fermionic","LocalOrbs("//str(sitendx)//")%RotDag")
            !
            if(doBare) call tensor_transform(reg(type_),Wbtmp,PhysicalUelements%Full_Map,LocalOrbs(sitendx)%RotDag,LocalOrbs(sitendx)%RotDag)
            do ip=1,Wimp%Npoints
               call tensor_transform(reg(type_),Wstmp(:,:,ip),PhysicalUelements%Full_Map,LocalOrbs(sitendx)%RotDag,LocalOrbs(sitendx)%RotDag)
            enddo
            !
         endif
         !
         !Insert into local indexes
         do i_imp=1,LocalOrbs(sitendx)%Norb
            do j_imp=1,LocalOrbs(sitendx)%Norb
               do k_imp=1,LocalOrbs(sitendx)%Norb
                  do l_imp=1,LocalOrbs(sitendx)%Norb
                     !
                     ! bosonic indexes of the impurity
                     call F2Bindex(Norb_imp,[i_imp,j_imp],[k_imp,l_imp],ib_imp,jb_imp)
                     !
                     ! mapping
                     i_loc = LocalOrbs(sitendx)%Orbs(i_imp)
                     j_loc = LocalOrbs(sitendx)%Orbs(j_imp)
                     k_loc = LocalOrbs(sitendx)%Orbs(k_imp)
                     l_loc = LocalOrbs(sitendx)%Orbs(l_imp)
                     !
                     ! corresponding bosonic indexes on the lattice
                     call F2Bindex(Norb_loc,[i_loc,j_loc],[k_loc,l_loc],ib_loc,jb_loc)
                     !
                     if(doBare)Wloc%bare_local(ib_loc,jb_loc) = Wbtmp(ib_imp,jb_imp)
                     do ip=1,Wimp%Npoints
                        Wloc%screened_local(ib_loc,jb_loc,ip) = Wstmp(ib_imp,jb_imp,ip)
                     enddo
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      enddo
      if(doBare) deallocate(Wbtmp)
      deallocate(Wstmp)
      !
   end subroutine imp2loc_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE: Expand an impurity matrix/field to a given number of sites
   !---------------------------------------------------------------------------!
   subroutine Expand2Nsite_Fermionic(G,Gsite,Nsite)
      !
      use parameters
      use utils_misc
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      type(FermionicField),intent(in)       :: Gsite
      integer,intent(in)                    :: Nsite
      !
      integer                               :: isite,shift
      integer                               :: Norb_site,Norb
      logical                               :: condGsite,condG
      !
      !
      if(verbose)write(*,"(A)") "---- Expand2Nsite_Fermionic"
      !
      !
      if(.not.G%status) stop "Expand2Nsite_Fermionic: G not properly initialized."
      if(.not.Gsite%status) stop "Expand2Nsite_Fermionic: Gsite not properly initialized."
      if(G%Norb.eq.0) stop "Expand2Nsite_Fermionic: Norb of G not defined."
      if(Gsite%Norb.eq.0) stop "Expand2Nsite_Fermionic: Norb of Gsite not defined."
      if(G%Npoints.eq.0) stop "Expand2Nsite_Fermionic: Npoints of G not defined."
      if(Gsite%Npoints.eq.0) stop "Expand2Nsite_Fermionic: Npoints of Gsite not defined."
      if(G%Beta.ne.Gsite%Beta) stop "Expand2Nsite_Fermionic: Gsite and G have different beta."
      if(G%Npoints.ne.Gsite%Npoints) stop "Expand2Nsite_Fermionic: Gsite and G have different number of Matsubara points."
      !
      Norb = G%Norb
      Norb_site = Gsite%Norb
      !
      if(Norb.ne.(Norb_site*Nsite)) stop "Expand2Nsite_Fermionic: output G dimension does not match Nsite times Gsite orbitals."
      !
      call clear_attributes(G)
      do isite=1,Nsite
         !
         shift = (isite-1)*Norb_site
         !
         condG = allocated(G%wks)
         condGsite = allocated(Gsite%wks)
         if(condG.and.condGsite) G%wks(1+shift:Norb_site+shift,1+shift:Norb_site+shift,:,:,:) = Gsite%wks
         !
         condG = allocated(G%ws)
         condGsite = allocated(Gsite%ws)
         if(condG.and.condGsite) G%ws(1+shift:Norb_site+shift,1+shift:Norb_site+shift,:,:) = Gsite%ws
         !
         condG = allocated(G%N_ks)
         condGsite = allocated(Gsite%N_ks)
         if(condG.and.condGsite) G%N_ks(1+shift:Norb_site+shift,1+shift:Norb_site+shift,:,:) = Gsite%N_ks
         !
         condG = allocated(G%N_s)
         condGsite = allocated(Gsite%N_s)
         if(condG.and.condGsite) G%N_s(1+shift:Norb_site+shift,1+shift:Norb_site+shift,:) = Gsite%N_s
         !
      enddo
      !
   end subroutine Expand2Nsite_Fermionic
   !
   subroutine Expand2Nsite_Matrix(O,Osite,Nsite)
      !
      use parameters
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)              :: O(:,:)
      complex(8),intent(in)                 :: Osite(:,:)
      integer,intent(in)                    :: Nsite
      !
      integer                               :: isite,shift
      integer                               :: Norb_site,Norb
      !
      !
      if(verbose)write(*,"(A)") "---- Expand2Nsite_Matrix"
      !
      !
      if(size(O,dim=1).ne.size(O,dim=2)) stop "Expand2Nsite_Matrix: O not square."
      if(size(Osite,dim=1).ne.size(Osite,dim=2)) stop "Expand2Nsite_Matrix: Osite not square."
      !
      Norb = size(O,dim=1)
      Norb_site = size(Osite,dim=1)
      !
      if(Norb.ne.(Norb_site*Nsite)) stop "Expand2Nsite_Matrix: output Matrix dimension does not match Nsite times site Matrix orbitals."
      !
      O=czero
      do isite=1,Nsite
         !
         shift = (isite-1)*Norb_site
         !
         O(1+shift:Norb_site+shift,1+shift:Norb_site+shift) = Osite
         !
      enddo
      !
   end subroutine Expand2Nsite_Matrix
   !
   subroutine Expand2Nsite_Matrix_s(O,Osite,Nsite)
      !
      use parameters
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)              :: O(:,:,:)
      complex(8),intent(in)                 :: Osite(:,:,:)
      integer,intent(in)                    :: Nsite
      !
      integer                               :: isite,shift
      integer                               :: Norb_site,Norb
      !
      !
      if(verbose)write(*,"(A)") "---- Expand2Nsite_Matrix_s"
      !
      !
      if(size(O,dim=1).ne.size(O,dim=2)) stop "Expand2Nsite_Matrix_s: O not square."
      if(size(Osite,dim=1).ne.size(Osite,dim=2)) stop "Expand2Nsite_Matrix_s: Osite not square."
      if(size(O,dim=3).ne.size(Osite,dim=3)) stop "Expand2Nsite_Matrix_s: third dimension is different between O and Osite."
      !
      Norb = size(O,dim=1)
      Norb_site = size(Osite,dim=1)
      !
      if(Norb.ne.(Norb_site*Nsite)) stop "Expand2Nsite_Matrix_s: output Matrix dimension does not match Nsite times site Matrix orbitals."
      !
      O=czero
      do isite=1,Nsite
         !
         shift = (isite-1)*Norb_site
         !
         O(1+shift:Norb_site+shift,1+shift:Norb_site+shift,:) = Osite
         !
      enddo
      !
   end subroutine Expand2Nsite_Matrix_s
   !
   subroutine Expand2Nsite_Bosonic(W,Wsite,Nsite)
      !
      use parameters
      use utils_misc
      implicit none
      !
      type(BosonicField),intent(inout)      :: W
      type(BosonicField),intent(in)         :: Wsite
      integer,intent(in)                    :: Nsite
      !
      integer                               :: isite,shift
      integer                               :: Norb_site,Norb
      integer                               :: i_site,j_site,k_site,l_site
      integer                               :: ib_site,jb_site
      integer                               :: i,j,k,l
      integer                               :: ib,jb
      logical                               :: condWsite,condW
      !
      !
      if(verbose)write(*,"(A)") "---- Expand2Nsite_Bosonic"
      !
      !
      if(.not.W%status) stop "Expand2Nsite_Bosonic: W not properly initialized."
      if(.not.Wsite%status) stop "Expand2Nsite_Bosonic: Wsite not properly initialized."
      if(W%Nbp.eq.0) stop "Expand2Nsite_Bosonic: Norb of W not defined."
      if(Wsite%Nbp.eq.0) stop "Expand2Nsite_Bosonic: Norb of Wsite not defined."
      if(W%Npoints.eq.0) stop "Expand2Nsite_Bosonic: Npoints of W not defined."
      if(Wsite%Npoints.eq.0) stop "Expand2Nsite_Bosonic: Npoints of Wsite not defined."
      if(W%Beta.ne.Wsite%Beta) stop "Expand2Nsite_Bosonic: Wsite and W have different beta."
      if(W%Npoints.ne.Wsite%Npoints) stop "Expand2Nsite_Bosonic: Wsite and W have different number of Matsubara points."
      !
      Norb = int(sqrt(dble(W%Nbp)))
      Norb_site = int(sqrt(dble(Wsite%Nbp)))
      !
      if(Norb.ne.(Norb_site*Nsite)) stop "Expand2Nsite_Bosonic: output W dimension does not match Nsite times Wsite orbitals."
      !
      call clear_attributes(W)
      do isite=1,Nsite
         !
         shift = (isite-1)*Norb_site
         !
         do i_site=1,Norb_site
            do j_site=1,Norb_site
               do k_site=1,Norb_site
                  do l_site=1,Norb_site
                     !
                     ! bosonic indexes of the impurity
                     ib_site = j_site + Norb_site*(i_site-1)
                     jb_site = l_site + Norb_site*(k_site-1)
                     !
                     ! mapping
                     i = i_site + shift
                     j = j_site + shift
                     k = k_site + shift
                     l = l_site + shift
                     !
                     ! corresponding bosonic indexes on the lattice
                     call F2Bindex(Norb,[i,j],[k,l],ib,jb)
                     !
                     condW = allocated(W%screened)
                     condWsite = allocated(Wsite%screened)
                     if(condW.and.condWsite) W%screened(ib,jb,:,:) = Wsite%screened(ib_site,jb_site,:,:)
                     !
                     condW = allocated(W%bare)
                     condWsite = allocated(Wsite%bare)
                     if(condW.and.condWsite) W%bare(ib,jb,:) = Wsite%bare(ib_site,jb_site,:)
                     !
                     condW = allocated(W%screened_local)
                     condWsite = allocated(Wsite%screened_local)
                     if(condW.and.condWsite) W%screened_local(ib,jb,:) = Wsite%screened_local(ib_site,jb_site,:)
                     !
                     condW = allocated(W%bare_local)
                     condWsite = allocated(Wsite%bare_local)
                     if(condW.and.condWsite) W%bare_local(ib,jb) = Wsite%bare_local(ib_site,jb_site)
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      enddo
      !
   end subroutine Expand2Nsite_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE: Symmetrize the Matrix/Field with respect to some lattice indexes.
   !         The Fermionic simmetrization is done only on the local attributes
   !         The Bosonic simmetrization is done only on the physical elements
   !---------------------------------------------------------------------------!
   subroutine symmetrize_GW_Matrix_d(Mat,Eqv,override)
      !
      use parameters
      implicit none
      !
      real(8),intent(inout)                 :: Mat(:,:,:)
      type(Equivalent),intent(in)           :: Eqv
      logical,intent(in),optional           :: override
      !
      real(8)                               :: dimdiag,dimoffdiag
      real(8)                               :: Delem,Uelem,Lelem
      integer                               :: iset,jset,iorb,jorb,ispin
      integer                               :: i,j
      logical                               :: orbsym
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_GW_Matrix_d"
      !
      !
      ! Check on the input Fields
      if(size(Mat,dim=1).ne.size(Mat,dim=2)) stop "symmetrize_GW_Matrix_d: Matix not square."
      orbsym = Eqv%O
      if(present(override))orbsym = override
      !
      if(Eqv%para.eq.1)then
         Mat(:,:,1) = (Mat(:,:,1) + Mat(:,:,Nspin))/2d0
         Mat(:,:,Nspin) = Mat(:,:,1)
      endif
      !
      if(orbsym)then
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
   end subroutine symmetrize_GW_Matrix_d
   !
   subroutine symmetrize_GW_Matrix_z(Mat,Eqv,override)
      !
      use parameters
      implicit none
      !
      complex(8),intent(inout)              :: Mat(:,:,:)
      type(Equivalent),intent(in)           :: Eqv
      logical,intent(in),optional           :: override
      !
      real(8)                               :: dimdiag,dimoffdiag
      complex(8)                            :: Delem,Uelem,Lelem
      integer                               :: iset,jset,iorb,jorb,ispin
      integer                               :: i,j
      logical                               :: orbsym
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_GW_Matrix_z"
      !
      !
      ! Check on the input Fields
      if(size(Mat,dim=1).ne.size(Mat,dim=2)) stop "symmetrize_GW_Matrix_z: Matix not square."
      orbsym = Eqv%O
      if(present(override))orbsym = override
      !
      if(Eqv%para.eq.1)then
         Mat(:,:,1) = (Mat(:,:,1) + Mat(:,:,Nspin))/2d0
         Mat(:,:,Nspin) = Mat(:,:,1)
      endif
      !
      if(orbsym)then
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
   end subroutine symmetrize_GW_Matrix_z
   !
   subroutine symmetrize_GW_Fermionic(G,Eqv,override)
      !
      use parameters
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      type(Equivalent),intent(in)           :: Eqv
      logical,intent(in),optional           :: override
      !
      real(8)                               :: dimdiag,dimoffdiag
      complex(8)                            :: Delem,Uelem,Lelem
      integer                               :: iset,jset,iorb,jorb,ispin
      integer                               :: i,j,ip
      logical                               :: orbsym
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_GW_Fermionic"
      !
      !
      ! Check on the input Fields
      if(.not.G%status) stop "symmetrize_GW_Fermionic: field not properly initialized."
      orbsym = Eqv%O
      if(present(override))orbsym = override
      !
      if(Eqv%para.eq.1)then
         !
         G%N_s(:,:,1) = (G%N_s(:,:,1) + G%N_s(:,:,Nspin))/2d0
         G%N_s(:,:,Nspin) = G%N_s(:,:,1)
         !
         do ip=1,G%Npoints
            G%ws(:,:,ip,1) = (G%ws(:,:,ip,1) + G%ws(:,:,ip,Nspin))/2d0
            G%ws(:,:,ip,Nspin) = G%ws(:,:,ip,1)
         enddo
         !
      endif
      !
      if(orbsym)then
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
   end subroutine symmetrize_GW_Fermionic
   !
   subroutine symmetrize_GW_Bosonic(W,Eqv,override)
      !
      use parameters
      use utils_misc
      implicit none
      !
      type(BosonicField),intent(inout)      :: W
      type(Equivalent),intent(in)           :: Eqv
      logical,intent(in),optional           :: override
      !
      real(8)                               :: dimdiag,dimoffdiag
      complex(8)                            :: Waaaa,Waabb
      complex(8)                            :: Wabab,Wabba
      integer                               :: iset,jset,ip,Norb
      integer                               :: iorb,jorb,i,j,ib1_aa
      integer                               :: ib1_ab,ib2_ab,ib1_sf,ib2_sf
      integer                               :: ib1_pa,ib2_pa,ib1_pb,ib2_pb
      logical                               :: doBare,orbsym
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_GW_Bosonic"
      !
      !
      ! Check on the input Fields
      if(.not.W%status) stop "symmetrize_GW_Bosonic: field not properly initialized."
      Norb = sqrt(dble(W%Nbp))
      doBare = allocated(W%bare_local)
      orbsym = Eqv%O
      if(present(override))orbsym = override
      !
      if(orbsym)then
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
                     call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                     if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(aa)(bb): (",i,i,"),(",j,j,") = [",ib1_ab,",",ib2_ab,"]"
                     !
                     call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                     if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ba): (",i,j,"),(",j,i,") = [",ib1_sf,",",ib2_sf,"]"
                     !
                     call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                     if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ab): (",i,j,"),(",i,j,") = [",ib1_pa,",",ib2_pa,"]"
                     call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
                     call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                     call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                     call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                     call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
                     call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                     call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                     call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                     call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
                     call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                     call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                     call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                     call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
                        if(i.eq.j) stop "symmetrize_GW_Bosonic: different sets cannot contain the same orbital index."
                        !
                        call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                        if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(aa)(bb): (",i,i,"),(",j,j,") = [",ib1_ab,",",ib2_ab,"]"
                        !
                        call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                        if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ba): (",i,j,"),(",j,i,") = [",ib1_sf,",",ib2_sf,"]"
                        !
                        call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                        if(verbose)write(*,"(2(A,2I3),2(A,I4),A)")    "     W(ab)(ab): (",i,j,"),(",i,j,") = [",ib1_pa,",",ib2_pa,"]"
                        call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
                        call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                        call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                        call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                        call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
                        if(i.eq.j) stop "symmetrize_GW_Bosonic: different sets cannot contain the same orbital index."
                        !
                        call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                        call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                        call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                        call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
                        call F2Bindex(Norb,[i,i],[j,j],ib1_ab,ib2_ab)
                        call F2Bindex(Norb,[i,j],[j,i],ib1_sf,ib2_sf)
                        call F2Bindex(Norb,[i,j],[i,j],ib1_pa,ib2_pa)
                        call F2Bindex(Norb,[j,i],[j,i],ib1_pb,ib2_pb)
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
   end subroutine symmetrize_GW_Bosonic


   !---------------------------------------------------------------------------!
   !PURPOSE: Symmetrize the Matrix/Field with respect to some impurity indexes.
   !         The given pattern is recognized and the GW symmetrization is called
   !         only for diagonal indexes.
   !---------------------------------------------------------------------------!
   subroutine symmetrize_imp_Matrix_d(Mat,Pattern)
      !
      use parameters
      use utils_misc
      use linalg, only: diag, diagonal
      implicit none
      !
      real(8),intent(inout)                 :: Mat(:,:,:)
      real(8),intent(in)                    :: Pattern(:)
      !
      type(Equivalent)                      :: Eqv_diag
      integer                               :: ispin,Norb!,iset
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_imp_Matrix_d"
      !
      !
      ! Check on the input Fields
      Norb = size(Mat,dim=1)
      call assert_shape(Mat,[Norb,Norb,Nspin],"symmetrize_imp_Matrix_d","Mat")
      call assert_shape(Pattern,[Norb],"symmetrize_imp_Matrix_d","Mat")
      !
      !remove accidental off-diagonal stuff
      do ispin=1,Nspin
         Mat(:,:,ispin) = diag(diagonal(Mat(:,:,ispin)))
      enddo
      !
      !The spin symmetrization is done in utils_main
      Eqv_diag%para = 0
      !
      !The off-diagonal components are absent by constuction
      Eqv_diag%Gfoffdiag = .false.
      !
      !get the equivalent sets from the input pattern
      call get_pattern(Eqv_diag%SetOrbs,Pattern,1e4*eps,listDim=Eqv_diag%SetNorb)
      !
      !fill in other stuff
      if(allocated(Eqv_diag%SetOrbs))then
         !
         Eqv_diag%Nset = size(Eqv_diag%SetOrbs,dim=1)
         Eqv_diag%Ntotset = Eqv_diag%Nset
         Eqv_diag%O = .true.
         !
         !allocate(Eqv_diag%SetNorb(Eqv_diag%Nset))
         !do iset=1,Eqv_diag%Nset
         !   Eqv_diag%SetNorb(iset) = size( pack( Eqv_diag%SetOrbs(iset,:), Eqv_diag%SetOrbs(iset,:).gt.0 ) )
         !enddo
         !
         !call the GW subroutine
         call symmetrize_GW(Mat,Eqv_diag)
         !
      else
         !
         write(*,"(A)") "     Warning symmetrize_imp_Matrix_d: nothing to symmetrize from the pattern."
         !
      endif
      !
   end subroutine symmetrize_imp_Matrix_d
   !
   subroutine symmetrize_imp_Matrix_z(Mat,Pattern)
      !
      use parameters
      use utils_misc
      use linalg, only: diag, diagonal
      implicit none
      !
      complex(8),intent(inout)              :: Mat(:,:,:)
      real(8),intent(in)                    :: Pattern(:)
      !
      type(Equivalent)                      :: Eqv_diag
      integer                               :: ispin,Norb!,iset
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_imp_Matrix_z"
      !
      !
      ! Check on the input Fields
      Norb = size(Mat,dim=1)
      call assert_shape(Mat,[Norb,Norb,Nspin],"symmetrize_imp_Matrix_z","Mat")
      call assert_shape(Pattern,[Norb],"symmetrize_imp_Matrix_z","Mat")
      !
      !remove accidental off-diagonal stuff
      do ispin=1,Nspin
         Mat(:,:,ispin) = diag(diagonal(Mat(:,:,ispin)))
      enddo
      !
      !The spin symmetrization is done in utils_main
      Eqv_diag%para = 0
      !
      !The off-diagonal components are absent by constuction
      Eqv_diag%Gfoffdiag = .false.
      !
      !get the equivalent sets from the input pattern
      call get_pattern(Eqv_diag%SetOrbs,Pattern,1e4*eps,listDim=Eqv_diag%SetNorb)
      !
      !fill in other stuff
      if(allocated(Eqv_diag%SetOrbs))then
         !
         Eqv_diag%Nset = size(Eqv_diag%SetOrbs,dim=1)
         Eqv_diag%Ntotset = Eqv_diag%Nset
         Eqv_diag%O = .true.
         !
         !allocate(Eqv_diag%SetNorb(Eqv_diag%Nset))
         !do iset=1,Eqv_diag%Nset
         !   Eqv_diag%SetNorb(iset) = size( pack( Eqv_diag%SetOrbs(iset,:), Eqv_diag%SetOrbs(iset,:).gt.0 ) )
         !enddo
         !
         !call the GW subroutine
         call symmetrize_GW(Mat,Eqv_diag)
         !
      else
         !
         write(*,"(A)") "     Warning symmetrize_imp_Matrix_z: nothing to symmetrize from the pattern."
         !
      endif
      !
   end subroutine symmetrize_imp_Matrix_z
   !
   subroutine symmetrize_imp_Fermionic(G,Pattern)
      !
      use parameters
      use utils_misc
      use linalg, only: diag, diagonal
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      real(8),intent(in)                    :: Pattern(:)
      !
      type(Equivalent)                      :: Eqv_diag
      integer                               :: iw,ispin,iset
      !
      !
      if(verbose)write(*,"(A)") "---- symmetrize_imp_Fermionic"
      !
      !
      ! Check on the input Fields
      if(.not.G%status) stop "symmetrize_imp_Fermionic: field not properly initialized."
      if(G%Norb.ne.size(Pattern)) stop "symmetrize_imp_Fermionic: field and pattern dimension mismatch."
      !
      !remove accidental off-diagonal stuff
      do ispin=1,Nspin
         G%N_s(:,:,ispin) = diag(diagonal(G%N_s(:,:,ispin)))
         do iw=1,G%Npoints
            G%ws(:,:,iw,ispin) = diag(diagonal(G%ws(:,:,iw,ispin)))
         enddo
      enddo
      !
      !The spin symmetrization is done in utils_main
      Eqv_diag%para = 0
      !
      !The off-diagonal components are absent by constuction
      Eqv_diag%Gfoffdiag = .false.
      !
      !get the equivalent sets from the input pattern
      call get_pattern(Eqv_diag%SetOrbs,Pattern,1e4*eps)
      !
      !fill in other stuff
      if(allocated(Eqv_diag%SetOrbs))then
         !
         Eqv_diag%Nset = size(Eqv_diag%SetOrbs,dim=1)
         Eqv_diag%Ntotset = Eqv_diag%Nset
         Eqv_diag%O = .true.
         !
         allocate(Eqv_diag%SetNorb(Eqv_diag%Nset))
         do iset=1,Eqv_diag%Nset
            Eqv_diag%SetNorb(iset) = size( pack( Eqv_diag%SetOrbs(iset,:), Eqv_diag%SetOrbs(iset,:).gt.0 ) )
         enddo
         !
         !call the GW subroutine
         call symmetrize_GW(G,Eqv_diag)
         !
      else
         !
         write(*,"(A)") "     Warning symmetrize_imp_Fermionic: nothing to symmetrize from the pattern."
         !
      endif
      !
   end subroutine symmetrize_imp_Fermionic


   !---------------------------------------------------------------------------!
   !PURPOSE: Join the C and X component of the self-energy
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
      if(.not.SigmaFull%status) stop "join_SigmaCX: SigmaFull not properly initialized."
      if(.not.Sigma_C%status) stop "join_SigmaCX: Sigma_C not properly initialized."
      if(.not.Sigma_X%status) stop "join_SigmaCX: Sigma_X not properly initialized."
      !
      if(SigmaFull%Nkpt.eq.0) local=.true.
      if(SigmaFull%Nkpt.ne.0) local=.false.
      !
      if(local)then
         if(Sigma_C%Nkpt.ne.0) stop "join_SigmaCX: Sigma_C k dependent attributes are supposed to be unallocated."
         if(Sigma_X%Nkpt.ne.0) stop "join_SigmaCX: Sigma_X k dependent attributes are supposed to be unallocated."
      else
         if(Sigma_C%Nkpt.eq.0) stop "join_SigmaCX: Sigma_C k dependent attributes not properly initialized."
         if(Sigma_X%Nkpt.eq.0) stop "join_SigmaCX: Sigma_X k dependent attributes not properly initialized."
      endif
      if(Sigma_X%Npoints.ne.0) stop "join_SigmaCX: Sigma_X frequency dependent attributes are supposed to be unallocated."
      !
      Norb = SigmaFull%Norb
      Beta = SigmaFull%Beta
      Nmats = SigmaFull%Npoints
      !
      if(.not.local)then
         Nkpt = SigmaFull%Nkpt
         if(all([Sigma_C%Nkpt-Nkpt,Sigma_X%Nkpt-Nkpt].ne.[0,0])) stop "join_SigmaCX: Either Sigma_C or Sigma_X have different number of k-points with respect to SigmaFull."
      endif
      !
      if(all([Sigma_C%Beta-Beta,Sigma_X%Beta-Beta].ne.[0d0,0d0])) stop "join_SigmaCX: Either Sigma_C or Sigma_X have different Beta with respect to SigmaFull."
      if(Nmats.ne.Sigma_C%Npoints) stop "join_SigmaCX: Sigma_C has different number of Matsubara points with respect to SigmaFull."
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
   !PURPOSE: + Replace SigmaImp in SigmaGW at the indexes contained in orbs.
   !           The Hartee contribution computed as N*curlyU is stored in the
   !           SigmaImp%N_s attribute and it is removed during the merge.
   !         + Replace PiImp in PiGW at the indexes contained in orbs
   !         + The only difference between v1 and v2 is the optional SigmaGW_DC
   !---------------------------------------------------------------------------!
   subroutine Merge_SelfEnergy(SigmaGW,SigmaImp,coeff,LocalOrbs,OffDiag,SigmaGW_DC)
      !
      use parameters
      use utils_misc
      implicit none
      !
      type(FermionicField),intent(inout)    :: SigmaGW
      type(FermionicField),intent(in)       :: SigmaImp
      real(8),intent(in)                    :: coeff(2)
      type(LocalOrbitals),allocatable,intent(in):: LocalOrbs(:)
      logical,intent(in),optional           :: OffDiag
      type(FermionicField),intent(in),optional :: SigmaGW_DC
      !
      real(8)                               :: Beta
      integer                               :: iw,ik,isite,iorb,jorb
      integer                               :: ispin,Norb_imp
      integer                               :: i_loc,j_loc
      integer                               :: Nkpt,Norb,Nmats,Nsite,Nimp
      logical                               :: impDC,OffDiag_
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- Merge_SelfEnergy"
      !
      !
      ! Check on the input Fields
      if(.not.SigmaGW%status) stop "Merge_SelfEnergy: SigmaGW not properly initialized."
      if(.not.SigmaImp%status) stop "Merge_SelfEnergy: SigmaImp not properly initialized."
      if(SigmaGW%Nkpt.eq.0) stop "Merge_SelfEnergy: SigmaGW k dependent attributes not properly initialized."
      if(SigmaImp%Nkpt.ne.0) stop "Merge_SelfEnergy: SigmaImp k dependent attributes are supposed to be unallocated."
      if(.not.allocated(LocalOrbs)) stop "Merge_SelfEnergy: LocalOrbs not properly initialized."
      !
      Norb = SigmaGW%Norb
      Nkpt = SigmaGW%Nkpt
      Beta = SigmaGW%Beta
      Nmats = SigmaGW%Npoints
      Nsite = SigmaGW%Nsite
      Nimp = size(LocalOrbs)
      !
      if(SigmaImp%Norb.ne.Norb) stop "Merge_SelfEnergy: SigmaImp has different orbital dimension with respect to SigmaGW."
      if(SigmaImp%Beta.ne.Beta) stop "Merge_SelfEnergy: SigmaImp has different Beta with respect to SigmaGW."
      if(SigmaImp%Npoints.ne.Nmats) stop "Merge_SelfEnergy: SigmaImp has different number of Matsubara points with respect to SigmaGW."
      !
      Norb_imp=0
      do isite=1,Nimp
         Norb_imp = Norb_imp + LocalOrbs(isite)%Norb
      enddo
      if(Norb_imp.gt.Norb) stop "Merge_SelfEnergy: Number of orbital to be inserted is bigger than the total lattice orbital space."
      !
      impDC = .true.
      if(present(SigmaGW_DC))then
         if(.not.SigmaGW_DC%status) stop "Merge_SelfEnergy: SigmaGW_DC not properly initialized."
         if(SigmaGW_DC%Nkpt.ne.0) stop "Merge_SelfEnergy: SigmaGW_DC k dependent attributes are supposed to be unallocated."
         if(SigmaGW_DC%Norb.ne.Norb) stop "Merge_SelfEnergy: SigmaGW_DC has different orbital dimension with respect to SigmaGW."
         if(SigmaGW_DC%Beta.ne.Beta) stop "Merge_SelfEnergy: SigmaGW_DC has different Beta with respect to SigmaGW."
         if(SigmaGW_DC%Npoints.ne.Nmats) stop "Merge_SelfEnergy: SigmaGW_DC has different number of Matsubara points with respect to SigmaGW."
         impDC = .false.
      endif
      !
      OffDiag_=.true.
      if(present(OffDiag)) OffDiag_=OffDiag
      !
      !Fill the local attributes so as to fully replace the local GW contibution
      call FermionicKsum(SigmaGW)
      !
      !all sites if(expand.or.AFM) otherwise only one site and the orbitals within orbs
      do isite=1,Nimp
         do iorb=1,LocalOrbs(isite)%Norb
            do jorb=1,LocalOrbs(isite)%Norb
               !
               if((.not.OffDiag_).and.(iorb.ne.jorb))cycle
               !
               !SigmaGW, SigmaImp and SigmaGW_DC have the same arrangements
               i_loc = LocalOrbs(isite)%Orbs(iorb)
               j_loc = LocalOrbs(isite)%Orbs(jorb)
               !
               !$OMP PARALLEL DEFAULT(SHARED),&
               !$OMP PRIVATE(ispin,iw,ik)
               !$OMP DO
               do iw=1,Nmats
                  do ik=1,Nkpt
                     do ispin=1,Nspin
                        !
                        if(impDC)then
                           SigmaGW%wks(i_loc,j_loc,iw,ik,ispin) = SigmaGW%wks(i_loc,j_loc,iw,ik,ispin)                 &
                                                                - coeff(1)*SigmaGW%ws(i_loc,j_loc,iw,ispin)            &
                                                                + coeff(1)*(SigmaImp%ws(i_loc,j_loc,iw,ispin)-coeff(2)*SigmaImp%N_s(i_loc,j_loc,ispin))
                        else
                           SigmaGW%wks(i_loc,j_loc,iw,ik,ispin) = SigmaGW%wks(i_loc,j_loc,iw,ik,ispin)                 &
                                                                - coeff(1)*SigmaGW_DC%ws(i_loc,j_loc,iw,ispin)         &
                                                                + coeff(1)*(SigmaImp%ws(i_loc,j_loc,iw,ispin)-coeff(2)*SigmaImp%N_s(i_loc,j_loc,ispin))
                        endif
                        !
                     enddo
                  enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL
               !
            enddo
         enddo
      enddo
      !
      !Put SigmaImp in the local attribute
      call FermionicKsum(SigmaGW)
      !
   end subroutine Merge_SelfEnergy
   !
   subroutine Merge_Polarization(PiGW,PiImp,coeff,LocalOrbs,OffDiag,PiGG_DC)
      !
      use parameters
      use utils_misc
      implicit none
      !
      type(BosonicField),intent(inout)      :: PiGW
      type(BosonicField),intent(in)         :: PiImp
      real(8),intent(in)                    :: coeff
      type(LocalOrbitals),allocatable,intent(in):: LocalOrbs(:)
      logical,intent(in),optional            :: OffDiag
      type(BosonicField),intent(in),optional :: PiGG_DC
      !
      real(8)                               :: Beta
      integer                               :: iw,ik,isite
      integer                               :: iorb,jorb,korb,lorb
      integer                               :: Norb_imp,ib_loc,jb_loc
      integer                               :: i_loc,j_loc,k_loc,l_loc
      integer                               :: Nkpt,Norb,Nmats,Nsite,Nimp
      logical                               :: impDC,OffDiag_
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- Merge_Polarization"
      !
      !
      ! Check on the input Fields
      if(.not.PiGW%status) stop "Merge_Polarization: PiGW not properly initialized."
      if(.not.PiImp%status) stop "Merge_Polarization: PiImp not properly initialized."
      if(PiGW%Nkpt.eq.0) stop "Merge_Polarization: PiGW k dependent attributes not properly initialized."
      if(PiImp%Nkpt.ne.0) stop "Merge_Polarization: PiImp k dependent attributes are supposed to be unallocated."
      if(.not.allocated(LocalOrbs)) stop "Merge_Polarization: LocalOrbs not properly initialized."
      !
      Norb = int(sqrt(dble(PiGW%Nbp)))
      Nkpt = PiGW%Nkpt
      Beta = PiGW%Beta
      Nmats = PiGW%Npoints
      Nsite = PiGW%Nsite
      Nimp = size(LocalOrbs)
      !
      if(PiImp%Beta.ne.Beta) stop "Merge_Polarization: PiImp has different Beta with respect to PiGW."
      if(PiImp%Npoints.ne.Nmats) stop "Merge_Polarization: PiImp has different number of Matsubara points with respect to PiGW."
      !
      Norb_imp=0
      do isite=1,Nimp
         Norb_imp = Norb_imp + LocalOrbs(isite)%Norb
      enddo
      if(Norb_imp.gt.Norb) stop "Merge_Polarization: Number of orbital to be inserted is bigger than the total lattice orbital space."
      !
      impDC = .true.
      if(present(PiGG_DC))then
         if(.not.PiGG_DC%status) stop "Merge_Polarization: PiGG_DC not properly initialized."
         if(PiGG_DC%Nkpt.ne.0) stop "Merge_Polarization: PiGG_DC k dependent attributes are supposed to be unallocated."
         if(PiGG_DC%Nbp.ne.PiGW%Nbp) stop "Merge_Polarization: PiGG_DC has different orbital dimension with respect to PiGG."
         if(PiGG_DC%Beta.ne.Beta) stop "Merge_Polarization: PiGG_DC has different Beta with respect to PiGG."
         if(PiGG_DC%Npoints.ne.Nmats) stop "Merge_Polarization: PiGG_DC has different number of Matsubara points with respect to PiGG."
         impDC = .false.
      endif
      !
      OffDiag_=.false. ! this refers to the components different than (aa)(bb)
      if(present(OffDiag)) OffDiag_=OffDiag
      !
      !Fill the local attributes so as to fully replace the local GW contibution
      call BosonicKsum(PiGW)
      !
      !all sites if(expand) otherwise only one site and the orbitals within orbs
      do isite=1,Nimp
         !
         do iorb=1,LocalOrbs(isite)%Norb
            do jorb=1,LocalOrbs(isite)%Norb
               do korb=1,LocalOrbs(isite)%Norb
                  do lorb=1,LocalOrbs(isite)%Norb
                     !
                     !mapping
                     i_loc = LocalOrbs(isite)%Orbs(iorb)
                     j_loc = LocalOrbs(isite)%Orbs(jorb)
                     k_loc = LocalOrbs(isite)%Orbs(korb)
                     l_loc = LocalOrbs(isite)%Orbs(lorb)
                     !
                     call F2Bindex(Norb,[i_loc,j_loc],[k_loc,l_loc],ib_loc,jb_loc)
                     !
                     if((.not.OffDiag_).and.(.not.((i_loc.eq.j_loc).and.(k_loc.eq.l_loc))))cycle
                     !if(.not.((i_loc.eq.k_loc).and.(l_loc.eq.j_loc)))cycle
                     !
                     !$OMP PARALLEL DEFAULT(SHARED),&
                     !$OMP PRIVATE(iw,ik)
                     !$OMP DO
                     do iw=1,Nmats
                        do ik=1,Nkpt
                           !
                           if(impDC)then
                              PiGW%screened(ib_loc,jb_loc,iw,ik) = PiGW%screened(ib_loc,jb_loc,iw,ik)                &
                                                                 - coeff*PiGW%screened_local(ib_loc,jb_loc,iw)       &
                                                                 + coeff*PiImp%screened_local(ib_loc,jb_loc,iw)
                           else
                              PiGW%screened(ib_loc,jb_loc,iw,ik) = PiGW%screened(ib_loc,jb_loc,iw,ik)                &
                                                                 - coeff*PiGG_DC%screened_local(ib_loc,jb_loc,iw)    &
                                                                 + coeff*PiImp%screened_local(ib_loc,jb_loc,iw)
                           endif
                           !
                        enddo
                     enddo
                     !$OMP END DO
                     !$OMP END PARALLEL
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      enddo !isite
      !
      !Put PiImp in the local attribute
      call BosonicKsum(PiGW)
      !
   end subroutine Merge_Polarization


   !---------------------------------------------------------------------------!
   !PURPOSE: Resize matrix between full product basis and just NaNb components
   !---------------------------------------------------------------------------!
   subroutine product2NN(Mfull,Mnanb)
      use parameters
      use utils_misc
      implicit none
      complex(8),intent(in)                 :: Mfull(:,:)
      complex(8),intent(inout)              :: Mnanb(:,:)
      integer                               :: Nbp,Norb
      integer                               :: i,j,ib,jb
      !
      if(size(Mfull,dim=1).ne.size(Mfull,dim=2)) stop "product2NN: input matrix not square."
      Nbp = size(Mfull,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      call assert_shape(Mnanb,[Norb,Norb],"product2NN","Mnanb")
      !
      do i=1,Norb
         do j=1,Norb
            !
            call F2Bindex(Norb,[i,i],[j,j],ib,jb)
            !
            Mnanb(i,j) = Mfull(ib,jb)
            !
         enddo
      enddo
      !
   end subroutine product2NN
   !
   subroutine NN2product(Mnanb,Mfull)
      use parameters
      use utils_misc
      implicit none
      complex(8),intent(in)                 :: Mnanb(:,:)
      complex(8),intent(inout)              :: Mfull(:,:)
      integer                               :: Nbp,Norb
      integer                               :: i,j,ib,jb
      !
      if(size(Mnanb,dim=1).ne.size(Mnanb,dim=2)) stop "NN2product: input matrix not square."
      Norb = size(Mnanb,dim=1)
      Nbp = Norb**2
      call assert_shape(Mfull,[Nbp,Nbp],"NN2product","Mfull")
      !
      do i=1,Norb
         do j=1,Norb
            !
            call F2Bindex(Norb,[i,i],[j,j],ib,jb)
            !
            Mfull(ib,jb) = Mnanb(i,j)
            !
         enddo
      enddo
      !
   end subroutine NN2product


   !---------------------------------------------------------------------------!
   !PURPOSE: compute the kinetic and potential energy
   !---------------------------------------------------------------------------!
   function calc_Ek(Gmats,Lttc) result(Ek)
      !
      use parameters
      use linalg, only : dag
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      real(8)                               :: Ek(Gmats%Norb,Gmats%Norb)
      real(8)                               :: Ek_T(Gmats%Norb,Gmats%Norb)
      integer                               :: ik,iw,ispin
      complex(8),allocatable                :: E(:,:),G(:,:)
      !
      if(.not.Gmats%status) stop "calc_Ek: Gmats not properly initialized."
      if(.not.Lttc%status) stop "calc_Ek: Lttc not properly initialized."
      if(Lttc%Nkpt.ne.Gmats%Nkpt) stop "calc_Ek: number of K-points does not match between Lttc and Gmats."
      if(Lttc%Norb.ne.Gmats%Norb) stop "calc_Ek: orbital dimension does not match between Lttc and Gmats."
      !
      Ek=0d0
      allocate(E(Gmats%Norb,Gmats%Norb));E=czero
      allocate(G(Gmats%Norb,Gmats%Norb));G=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Gmats,Lttc,Ek),&
      !$OMP PRIVATE(iw,ik,ispin,E,G,Ek_T)
      !$OMP DO
      do iw=1,Gmats%Npoints
         Ek_T=0d0
         do ik=1,Gmats%Nkpt
            do ispin=1,Nspin
               !
               E = Lttc%Hk(:,:,ik)
               G = Gmats%wks(:,:,iw,ik,ispin)
               !
               Ek_T = Ek_T + ( matmul(E,G) + matmul(dag(G),dag(E)) ) / (Gmats%Nkpt*Nspin)
               !
            enddo
         enddo
         Ek = Ek + Ek_T / Gmats%Beta
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(E,G)
      !
   end function calc_Ek
   !
   function calc_Ep(Gmats,Smats) result(Ep)
      !
      use parameters
      use linalg, only : dag, deye
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats
      type(FermionicField),intent(in)       :: Smats
      real(8)                               :: Ep(Gmats%Norb,Gmats%Norb)
      real(8)                               :: Ep_T(Gmats%Norb,Gmats%Norb)
      integer                               :: ik,iw,ispin
      complex(8),allocatable                :: U(:,:),G(:,:)
      !
      if(.not.Gmats%status) stop "calc_Ep: Gmats not properly initialized."
      if(.not.Smats%status) stop "calc_Ep: Smats not properly initialized."
      if(Smats%Nkpt.ne.Gmats%Nkpt) stop "calc_Ep: number of K-points does not match between Smats and Gmats."
      if(Smats%Beta.ne.Gmats%Beta) stop "calc_Ep: Beta does not match between Smats and Gmats."
      if(Smats%Norb.ne.Gmats%Norb) stop "calc_Ep: orbital dimension does not match between Smats and Gmats."
      !
      Ep=0d0
      allocate(U(Gmats%Norb,Gmats%Norb));U=czero
      allocate(G(Gmats%Norb,Gmats%Norb));G=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Gmats,Smats,Ep),&
      !$OMP PRIVATE(iw,ik,ispin,U,G,Ep_T)
      !$OMP DO
      do iw=1,Gmats%Npoints
         Ep_T=0d0
         do ik=1,Gmats%Nkpt
            do ispin=1,Nspin
               !
               U = Smats%wks(:,:,iw,ik,ispin)-deye(Smats%Norb)*Gmats%mu
               G = Gmats%wks(:,:,iw,ik,ispin)
               !
               Ep_T = Ep_T + ( matmul(U,G) + matmul(dag(G),dag(U)) ) / (Gmats%Nkpt*Nspin)
               !
            enddo
         enddo
         Ep = Ep + Ep_T / Gmats%Beta
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(U,G)
      !
   end function calc_Ep


   !---------------------------------------------------------------------------!
   !PURPOSE: check if the quasiparticle equatin has solutions in w=0
   !---------------------------------------------------------------------------!
   subroutine check_QP_poles(Lttc,Smats)
      !
      use parameters
      use utils_misc
      use linalg, only : diagonal, rotate
      implicit none
      !
      type(Lattice),intent(in)              :: Lttc
      type(FermionicField),intent(in)       :: Smats
      integer                               :: ispin,ik,iorb
      real(8),allocatable                   :: mu_ReSigma0(:)
      logical                               :: pole
      !
      if(.not.Lttc%status) stop "check_QP_poles: Lttc not properly initialized."
      if(.not.Smats%status) stop "check_QP_poles: Smats not properly initialized."
      if(Lttc%Nkpt.ne.Smats%Nkpt) stop "check_QP_poles: number of K-points does not match between Lttc and Smats."
      if(Lttc%Norb.ne.Smats%Norb) stop "check_QP_poles: orbital dimension does not match between Lttc and Smats."
      !
      if(Smats%mu.eq.0d0)write(*,"(A)")"     check_QP_poles: warning, the chemical potential is exactly zero."
      !
      allocate(mu_ReSigma0(Lttc%Norb));mu_ReSigma0=0d0
      do ispin=1,Nspin
         kloop: do ik=1,Lttc%Nkpt
            !
            mu_ReSigma0 = Smats%mu - dreal(diagonal(rotate(Smats%wks(:,:,1,ik,ispin),Lttc%Zk(:,:,ik))))
            !
            do iorb=1,Lttc%Norb
               pole = (Lttc%D_lower.lt.mu_ReSigma0(iorb)) .and. (mu_ReSigma0(iorb).lt.Lttc%D_upper)
               if(pole) exit kloop
            enddo
            !
         enddo kloop
         !
         if(pole)then
            write(*,"(A)")"     Spin #"//str(ispin)//" is a metal."
         else
            write(*,"(A)")"     Spin #"//str(ispin)//" is an insulator."
         endif
         !
      enddo
      deallocate(mu_ReSigma0)
      !
   end subroutine check_QP_poles


   !---------------------------------------------------------------------------!
   !PURPOSE: Interface with the potential subroutine in utils_mist.
   !         This subroutine takes the frequency mesh, Hk, vertical hoppig array
   !         and optionally the self-energy of the whole heterostructure and
   !         extract the slices needed to construct the embedding potentials
   !         via a continued fraction representation with the periodicity of
   !         the two extremals layers.
   !---------------------------------------------------------------------------!
   subroutine build_Potential_Field(Potential,Hetero,ndx,Nbulk,zeta,Hk,tz,mode,paramagnet,Smats)
      !
      use parameters
      use linalg, only : inv, rotate
      use utils_misc
      implicit none
      !
      type(FermionicField),intent(inout)    :: Potential
      integer,intent(out)                   :: Nbulk
      integer,intent(out)                   :: ndx(2)
      type(Heterostructures),intent(inout)  :: Hetero
      complex(8),intent(in)                 :: zeta(:,:,:)
      complex(8),intent(in)                 :: Hk(:,:,:)
      complex(8),allocatable,intent(in)     :: tz(:,:,:,:)
      character(len=*),intent(in)           :: mode
      logical,intent(in)                    :: paramagnet
      type(FermionicField),intent(in),optional :: Smats
      !
      integer                               :: Norb,Nmats,Nkpt
      complex(8),allocatable                :: zeta_(:,:,:),Hk_(:,:,:)
      complex(8),allocatable                :: ta(:,:,:),tb(:,:,:)
      complex(8),allocatable                :: Sa(:,:,:,:,:),Sb(:,:,:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- build_Potential_Field"
      !
      !
      if(.not.Potential%status) stop "build_Potential_Field: Potential not properly initialized."
      if(.not.Hetero%status) stop "build_Potential_Field: Heterostructure not properly initialized."
      if(Potential%Norb.ne.Hetero%Norb) stop "build_Potential_Field: Potential and Heterostructure have different orbital dimension."
      !
      Norb = size(Hk,dim=1)
      Nmats = Potential%Npoints
      Nkpt = Potential%Nkpt
      !
      call assert_shape(zeta,[Norb,Norb,Nmats],"build_Potential_Field","zeta")
      call assert_shape(Hk,[Norb,Norb,Nkpt],"build_Potential_Field","Hk")
      call assert_shape(tz,[Hetero%Norb,Hetero%Norb,Nkpt,(Hetero%tzIndex(2)-Hetero%tzIndex(1)+1)],"build_Potential_Field","tz")
      !
      if(present(Smats))then
         if(.not.Smats%status) stop "build_Potential_Field: Smats present but not properly initialized."
         if(Smats%Norb.ne.Norb) stop "build_Potential_Field: Smats has wrong orbital dimension."
         if(Smats%Npoints.ne.Nmats) stop "build_Potential_Field: Smats has wrong Matsubara frequency mesh."
         if(Smats%Nkpt.ne.Nkpt) stop "build_Potential_Field: Smats has wrong k-point mesh."
         if(Smats%Beta.ne.Potential%Beta) stop "build_Potential_Field: Smats has different Beta with respect of the Potential."
      endif
      !
      call clear_attributes(Potential)
      allocate(ta(Hetero%Norb,Hetero%Norb,Nkpt));ta=czero
      allocate(tb(Hetero%Norb,Hetero%Norb,Nkpt));tb=czero
      allocate(Sa(Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin));Sa=czero
      allocate(Sb(Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin));Sb=czero
      !
      select case(reg(mode))
         case default
            !
            stop "build_Potential_Field: Available modes: left, right."
            !
         case("left")
            !
            if(Hetero%Explicit(1).eq.1) stop "build_Potential_Field: requested left potential but no bulk present."
            ndx(1) = 1
            ndx(2) = Hetero%Norb
            Nbulk = Hetero%Explicit(1)-1
            !
            !Connection-to and self-energy-of the first layer explicitly solved
            ta = tz(:,:,:,Hetero%Explicit(1)-1)
            if(present(Smats)) Sa = Smats%wks(ndx(1):ndx(2),ndx(1):ndx(2),:,:,:)
            !
            !Connection-to and self-energy-of the second layer explicitly solved
            tb = tz(:,:,:,Hetero%Explicit(1))
            if(present(Smats)) Sb = Smats%wks(ndx(1)+Hetero%Norb:ndx(2)+Hetero%Norb,ndx(1)+Hetero%Norb:ndx(2)+Hetero%Norb,:,:,:)
            !
         case("right")
            !
            if(Hetero%Explicit(2).eq.Hetero%Nslab) stop "build_Potential_Field: requested right potential but no bulk present."
            ndx(1) = 1+ Norb - Hetero%Norb
            ndx(2) = Norb
            Nbulk = Hetero%Nslab-Hetero%Explicit(2)
            !
            !Connection-to and self-energy-of the last layer explicitly solved
            ta = tz(:,:,:,Hetero%Explicit(2))
            if(present(Smats)) Sa = Smats%wks(ndx(1):ndx(2),ndx(1):ndx(2),:,:,:)
            !
            !Connection-to and self-energy-of the semi-last layer explicitly solved
            tb = tz(:,:,:,Hetero%Explicit(2)-1)
            if(present(Smats)) Sb = Smats%wks(ndx(1)-Hetero%Norb:ndx(2)-Hetero%Norb,ndx(1)-Hetero%Norb:ndx(2)-Hetero%Norb,:,:,:)
            !
      end select
      !
      if(mod(Nbulk,2).ne.0)Nbulk=Nbulk+1
      !
      !Diagonal arrays of the first/last layer explicitly solved
      allocate(zeta_(Hetero%Norb,Hetero%Norb,Nmats));  zeta_ = zeta(ndx(1):ndx(2),ndx(1):ndx(2),:)
      allocate(Hk_(Hetero%Norb,Hetero%Norb,Nkpt))   ;  Hk_ = Hk(ndx(1):ndx(2),ndx(1):ndx(2),:)
      !
      call Embedding_ContinuedFraction(Potential%wks,Nbulk,zeta_,Hk_,ta,tb,Sa,Sb,paramagnet)
      !
      call FermionicKsum(Potential)
      !
      !this is to be able to print it form main
      if(reg(mode).eq."left")then
         if(allocated(Hetero%P_L))deallocate(Hetero%P_L)
         Hetero%P_L = Potential%ws
      elseif(reg(mode).eq."right")then
         if(allocated(Hetero%P_R))deallocate(Hetero%P_R)
         Hetero%P_R = Potential%ws
      endif
      !
   end subroutine build_Potential_Field
   !
   subroutine build_Potential_Mat(Potential,Hetero,ndx,Nbulk,zeta,Hk,tz,mode,paramagnet,Smats)
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
      complex(8),allocatable,intent(in)     :: tz(:,:,:,:)
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
      if(verbose)write(*,"(A)") "---- build_Potential_Mat"
      !
      !
      if(.not.Hetero%status) stop "build_Potential_Mat: Heterostructure not properly initialized."
      if(size(Potential,dim=1).ne.Hetero%Norb) stop "build_Potential_Mat: Potential and Heterostructure have different orbital dimension."
      !
      Norb = size(Hk,dim=1)
      Nmats = size(Potential,dim=3)
      Nkpt = size(Potential,dim=4)
      !
      call assert_shape(Potential,[Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin],"build_Potential_Mat","Potential")
      call assert_shape(zeta,[Norb,Norb,Nmats],"build_Potential_Mat","zeta")
      call assert_shape(Hk,[Norb,Norb,Nkpt],"build_Potential_Mat","Hk")
      call assert_shape(tz,[Hetero%Norb,Hetero%Norb,Nkpt,(Hetero%tzIndex(2)-Hetero%tzIndex(1)+1)],"build_Potential_Mat","tz")
      !
      if(present(Smats))then
         if(size(Smats,dim=1).ne.Norb) stop "build_Potential_Mat: Smats has wrong orbital 1st dimension."
         if(size(Smats,dim=2).ne.Norb) stop "build_Potential_Mat: Smats has wrong orbital 2nd dimension."
         if(size(Smats,dim=3).ne.Nmats) stop "build_Potential_Mat: Smats has wrong Matsubara frequency mesh."
         if(size(Smats,dim=4).ne.Nkpt) stop "build_Potential_Mat: Smats has wrong k-point mesh."
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
            stop "build_Potential_Mat: Available modes: left, right."
            !
         case("left")
            !
            if(Hetero%Explicit(1).eq.1) stop "build_Potential_Mat: requested left potential but no bulk present."
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
            if(Hetero%Explicit(2).eq.Hetero%Nslab) stop "build_Potential_Mat: requested right potential but no bulk present."
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
   end subroutine build_Potential_Mat


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


end module utils_fields
