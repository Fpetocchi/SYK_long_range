!------------------------------------------------------------------------------!
!PURPOSE: Rotate a tensor diven by the kroenecker product of square matrices
!------------------------------------------------------------------------------!
subroutine tensor_transform_NNNN_c(mode,Utensor,rotL,rotR)
   !
   implicit none
   character(len=2),intent(in)              :: mode
   complex(8),intent(inout)                 :: Utensor(:,:,:,:)
   complex(8),intent(in)                    :: rotL(:,:)
   complex(8),intent(in)                    :: rotR(:,:)
   !
   complex(8),allocatable                   :: Umatrix(:,:)
   complex(8),allocatable                   :: X(:,:,:),Y(:,:,:)
   integer                                  :: Norb,Nbp,ib,ib1,ib2
   integer                                  :: i,j,k,l
   integer                                  :: a,b,c,d
   !
   Norb = size(Utensor,dim=1)
   Nbp = Norb*Norb
   call assert_shape(Utensor,[Norb,Norb,Norb,Norb],"tensor_transform_NNNN_c","Utensor")
   call assert_shape(rotL,[Norb,Norb],"tensor_transform_NNNN_c","rotL")
   call assert_shape(rotR,[Norb,Norb],"tensor_transform_NNNN_c","rotR")
   !
   allocate(X(Nbp,Norb,Norb));X=zero
   allocate(Y(Nbp,Norb,Norb));Y=zero
   !
   select case(mode)
      case default
         !
         stop "tensor_transform_NNNN_c: Available modes: GG, NN."
         !
      case("GG")
         !
         ! U_ijkl = sum_abcd [ (L*_ai G_ac L_ck) * (R*_dl G_db R_bj) ]
         ! X_ij(a,b) = L*_ai R_bj
         ! Y_kl(c,d) = L*_ck R_dl
         ! U_ijkl = sum_abcd G_acG_db X_ij(a,b) * Y*_kl(c,d)
         !
         do a=1,Norb
            do b=1,Norb
               !
               do i=1,Norb
                  do j=1,Norb
                     !
                     call F2Bindex(Norb,[i,j],ib)
                     !
                     X(ib,a,b) = dconjg(rotL(a,i))*rotR(b,j)
                     Y(ib,a,b) = dconjg(rotL(a,i))*rotR(b,j)
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !
         !
      case("NN")
         !
         ! U_ijkl = sum_abcd [ (L*_ai n_ab L_bj) * (R*_ck n_cd R_dl) ] U_abcd
         ! X_ij(a,b) = L*_ai n_ab L_bj
         ! Y_kl(c,d) = R_ck n_cd R*_dl
         ! U_ijkl = sum_abcd U_abcd X_ij(a,b) * Y*_kl(c,d)
         !
         do a=1,Norb
            do b=1,Norb
               !
               do i=1,Norb
                  do j=1,Norb
                     !
                     call F2Bindex(Norb,[i,j],ib)
                     !
                     X(ib,a,b) = dconjg(rotL(a,i))*rotL(b,j)
                     Y(ib,a,b) = rotR(a,i)*dconjg(rotR(b,j))
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !
         !
   end select
   !
   !performs the hermitian rank 2 operation: A := alpha*x*y**H + conjg( alpha )*y*x**H + A
   allocate(Umatrix(Nbp,Nbp));Umatrix = zero
   do a=1,Norb
      do b=1,Norb
         do c=1,Norb
            do d=1,Norb
               call ZHER2('U',Nbp,Utensor(a,b,c,d),X(:,a,b),1,Y(:,c,d),1,Umatrix,Nbp)
            enddo
         enddo
      enddo
   enddo
   deallocate(X,Y)
   do ib1=1,Nbp
      do ib2=1+ib1,Nbp
         Umatrix(ib2,ib1) = conjg(Umatrix(ib1,ib2))
      enddo
   enddo
   !
   !provide U_ijkl
   Utensor = zero
   do i=1,Norb
      do j=1,Norb
         do k=1,Norb
            do l=1,Norb
               !
               call F2Bindex(Norb,[i,j],ib1)
               call F2Bindex(Norb,[k,l],ib2)
               !
               Utensor(i,j,k,l) = Umatrix(ib1,ib2)/2d0
               !
            enddo
         enddo
      enddo
   enddo
   deallocate(Umatrix)
   !
   contains
      !
      ! internal product basis ordering
      subroutine F2Bindex(Norb,orbs,ib)
         implicit none
         integer,intent(in)                    :: orbs(2)
         integer,intent(in)                    :: Norb
         integer,intent(out)                   :: ib
         integer                               :: i,j
         !
         i = orbs(1)
         j = orbs(2)
         !
         ib = j + Norb*(i-1)
         !
      end subroutine F2Bindex
      !
end subroutine tensor_transform_NNNN_c
!
subroutine tensor_transform_NNNN_d(mode,Utensor,rotL,rotR)
   !
   implicit none
   character(len=2),intent(in)              :: mode
   real(8),intent(inout)                    :: Utensor(:,:,:,:)
   real(8),intent(in)                       :: rotL(:,:)
   real(8),intent(in)                       :: rotR(:,:)
   !
   real(8),allocatable                      :: Umatrix(:,:)
   real(8),allocatable                      :: X(:,:,:),Y(:,:,:)
   integer                                  :: Norb,Nbp,ib,ib1,ib2
   integer                                  :: i,j,k,l
   integer                                  :: a,b,c,d
   !
   Norb = size(Utensor,dim=1)
   Nbp = Norb*Norb
   call assert_shape(Utensor,[Norb,Norb,Norb,Norb],"tensor_transform_NNNN_d","Utensor")
   call assert_shape(rotL,[Norb,Norb],"tensor_transform_NNNN_d","rotL")
   call assert_shape(rotR,[Norb,Norb],"tensor_transform_NNNN_d","rotR")
   !
   allocate(X(Nbp,Norb,Norb));X=0d0
   allocate(Y(Nbp,Norb,Norb));Y=0d0
   !
   select case(mode)
      case default
         !
         stop "tensor_transform_NNNN_d: Available modes: GG, NN."
         !
      case("GG")
         !
         ! U_ijkl = sum_abcd [ (L*_ai G_ac L_ck) * (R*_dl G_db R_bj) ]
         ! X_ij(a,b) = L*_ai R_bj
         ! Y_kl(c,d) = L*_ck R_dl
         ! U_ijkl = sum_abcd G_acG_db X_ij(a,b) * Y*_kl(c,d)
         !
         do a=1,Norb
            do b=1,Norb
               !
               do i=1,Norb
                  do j=1,Norb
                     !
                     call F2Bindex(Norb,[i,j],ib)
                     !
                     X(ib,a,b) = rotL(a,i)*rotR(b,j)
                     Y(ib,a,b) = rotL(a,i)*rotR(b,j)
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !
         !
      case("NN")
         !
         ! U_ijkl = sum_abcd [ (L*_ai n_ab L_bj) * (R*_ck n_cd R_dl) ] U_abcd
         ! X_ij(a,b) = L*_ai n_ab L_bj
         ! Y_kl(c,d) = R_ck n_cd R*_dl
         ! U_ijkl = sum_abcd U_abcd X_ij(a,b) * Y*_kl(c,d)
         !
         do a=1,Norb
            do b=1,Norb
               !
               do i=1,Norb
                  do j=1,Norb
                     !
                     call F2Bindex(Norb,[i,j],ib)
                     !
                     X(ib,a,b) = rotL(a,i)*rotL(b,j)
                     Y(ib,a,b) = rotR(a,i)*rotR(b,j)
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !
         !
   end select
   !
   !performs the symmetric rank 2 operation: A := alpha*x*y**T + alpha*y*x**T + A
   allocate(Umatrix(Nbp,Nbp));Umatrix = 0d0
   do a=1,Norb
      do b=1,Norb
         do c=1,Norb
            do d=1,Norb
               call DSPR2('U',Nbp,Utensor(a,b,c,d),X(:,a,b),1,Y(:,c,d),1,Umatrix,Nbp)
            enddo
         enddo
      enddo
   enddo
   deallocate(X,Y)
   do ib1=1,Nbp
      do ib2=1+ib1,Nbp
         Umatrix(ib2,ib1) = Umatrix(ib1,ib2)
      enddo
   enddo
   !
   !provide U_ijkl
   Utensor = 0d0
   do i=1,Norb
      do j=1,Norb
         do k=1,Norb
            do l=1,Norb
               !
               call F2Bindex(Norb,[i,j],ib1)
               call F2Bindex(Norb,[k,l],ib2)
               !
               Utensor(i,j,k,l) = Umatrix(ib1,ib2)/2d0
               !
            enddo
         enddo
      enddo
   enddo
   deallocate(Umatrix)
   !
   contains
      !
      ! internal product basis ordering
      subroutine F2Bindex(Norb,orbs,ib)
         implicit none
         integer,intent(in)                    :: orbs(2)
         integer,intent(in)                    :: Norb
         integer,intent(out)                   :: ib
         integer                               :: i,j
         !
         i = orbs(1)
         j = orbs(2)
         !
         ib = j + Norb*(i-1)
         !
      end subroutine F2Bindex
      !
end subroutine tensor_transform_NNNN_d


!------------------------------------------------------------------------------!
!PURPOSE: If the tensor is in some given square-matrix representation the user
!         must provide the representation mapping of the input matrix in the
!         form of an index list pointing at a given element of the tensor.
!         - Map(io,jo,1) --> iorb
!         - Map(io,jo,2) --> jorb
!         - Map(io,jo,3) --> korb
!         - Map(io,jo,4) --> lorb
!------------------------------------------------------------------------------!
subroutine tensor_transform_NN_c(mode,Umatrix,Map,rotL,rotR)
   !
   implicit none
   character(len=2),intent(in)              :: mode
   complex(8),intent(inout)                 :: Umatrix(:,:)
   integer,intent(in)                       :: Map(:,:,:)
   complex(8),intent(in)                    :: rotL(:,:)
   complex(8),intent(in)                    :: rotR(:,:)
   !
   complex(8),allocatable                   :: Utensor(:,:,:,:)
   integer                                  :: Norb,Nbp
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo
   !
   Nbp = size(Umatrix,dim=1)
   call assert_shape(Umatrix,[Nbp,Nbp],"tensor_transform_NN_c","Umatrix")
   call assert_shape(Map,[Nbp,Nbp,4],"tensor_transform_NN_c","Map")
   !
   Norb = int(sqrt(dble(Nbp)))
   call assert_shape(rotL,[Norb,Norb],"tensor_transform_NN_c","rotL")
   call assert_shape(rotR,[Norb,Norb],"tensor_transform_NN_c","rotR")
   !
   allocate(Utensor(Norb,Norb,Norb,Norb));Utensor=zero
   do io=1,Nbp
      do jo=1,Nbp
         !
         iorb = Map(io,jo,1)
         jorb = Map(io,jo,2)
         korb = Map(io,jo,3)
         lorb = Map(io,jo,4)
         !
         Utensor(iorb,jorb,korb,lorb) = Umatrix(io,jo)
         !
      enddo
   enddo
   !
   call tensor_transform_NNNN_c(mode,Utensor,rotL,rotR)
   !
   Umatrix=zero
   do io=1,Nbp
      do jo=1,Nbp
         !
         iorb = Map(io,jo,1)
         jorb = Map(io,jo,2)
         korb = Map(io,jo,3)
         lorb = Map(io,jo,4)
         !
         Umatrix(io,jo) = Utensor(iorb,jorb,korb,lorb)
         !
      enddo
   enddo
   deallocate(Utensor)
   !
end subroutine tensor_transform_NN_c
!
subroutine tensor_transform_NN_d(mode,Umatrix,Map,rotL,rotR)
   !
   implicit none
   character(len=2),intent(in)              :: mode
   real(8),intent(inout)                    :: Umatrix(:,:)
   integer,intent(in)                       :: Map(:,:,:)
   real(8),intent(in)                       :: rotL(:,:)
   real(8),intent(in)                       :: rotR(:,:)
   !
   real(8),allocatable                      :: Utensor(:,:,:,:)
   integer                                  :: Norb,Nbp
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo
   !
   Nbp = size(Umatrix,dim=1)
   call assert_shape(Umatrix,[Nbp,Nbp],"tensor_transform_NN_d","Umatrix")
   call assert_shape(Map,[Nbp,Nbp,4],"tensor_transform_NN_d","Map")
   !
   Norb = int(sqrt(dble(Nbp)))
   call assert_shape(rotL,[Norb,Norb],"tensor_transform_NN_d","rotL")
   call assert_shape(rotR,[Norb,Norb],"tensor_transform_NN_d","rotR")
   !
   allocate(Utensor(Norb,Norb,Norb,Norb));Utensor=0d0
   do io=1,Nbp
      do jo=1,Nbp
         !
         iorb = Map(io,jo,1)
         jorb = Map(io,jo,2)
         korb = Map(io,jo,3)
         lorb = Map(io,jo,4)
         !
         Utensor(iorb,jorb,korb,lorb) = Umatrix(io,jo)
         !
      enddo
   enddo
   !
   call tensor_transform_NNNN_d(mode,Utensor,rotL,rotR)
   !
   Umatrix=0d0
   do io=1,Nbp
      do jo=1,Nbp
         !
         iorb = Map(io,jo,1)
         jorb = Map(io,jo,2)
         korb = Map(io,jo,3)
         lorb = Map(io,jo,4)
         !
         Umatrix(io,jo) = Utensor(iorb,jorb,korb,lorb)
         !
      enddo
   enddo
   deallocate(Utensor)
   !
end subroutine tensor_transform_NN_d
