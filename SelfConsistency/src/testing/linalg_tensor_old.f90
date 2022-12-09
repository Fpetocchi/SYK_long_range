!------------------------------------------------------------------------------!
!PURPOSE: Rotate a tensor diven by the kroenecker product of square matrices
!------------------------------------------------------------------------------!
subroutine tensor_transform_NNNN_d(mode,Utensor,rotL,rotR)
   !
   implicit none
   character(len=2),intent(in)              :: mode
   real(8),intent(inout)                    :: Utensor(:,:,:,:)
   complex(8),intent(in)                    :: rotL(:,:)
   complex(8),intent(in)                    :: rotR(:,:)
   !
   complex(8),allocatable                   :: Umatrix(:,:)
   complex(8),allocatable                   :: Nab(:,:,:),NabL(:,:,:),NabR(:,:,:)
   integer,allocatable                      :: stride(:,:)
   integer                                  :: Norb
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo,count
   !
   Norb = size(Utensor,dim=1)
   call assert_shape(Utensor,[Norb,Norb,Norb,Norb],"tensor_transform_NNNN_d","Utensor")
   call assert_shape(rotL,[Norb,Norb],"tensor_transform_NNNN_d","rotL")
   call assert_shape(rotR,[Norb,Norb],"tensor_transform_NNNN_d","rotR")
   !
   allocate(Nab(Norb,Norb,Norb*Norb));Nab=zero
   allocate(stride(Norb,Norb));stride=0
   count=0
   do iorb=1,Norb
      do jorb=1,Norb
         count=count+1
         Nab(jorb,iorb,count) = one !N=c^+,c = {col,row} --> transpose
         stride(iorb,jorb)=count
      enddo
   enddo
   !
   allocate(NabL(Norb,Norb,Norb*Norb));NabL=zero
   allocate(NabR(Norb,Norb,Norb*Norb));NabR=zero
   do io=1,Norb*Norb
      NabL(:,:,io) = rotate_zz(Nab(:,:,io),rotL)
      NabR(:,:,io) = rotate_zz(Nab(:,:,io),rotR)
   enddo
   deallocate(Nab)
   !
   allocate(Umatrix(Norb*Norb,Norb*Norb));Umatrix=zero
   select case(mode)
      case default
         !
         stop "tensor_transform_NNNN_d: Available modes: ZHER, TT."
         !
      case("GG")
         !
         !
         do io=1,Norb*Norb
            do jo=1,Norb*Norb
               !
               call get_element_from_stride(iorb,jorb,io)
               call get_element_from_stride(korb,lorb,jo)
               !
               ! U_(il)(jk) = G_ij * G_kl * U_ijkl
               Umatrix = Umatrix + kronecker_product(NabL(:,:,io),NabR(:,:,jo))*Utensor(iorb,lorb,jorb,korb)
               !
            enddo
         enddo
         !
         !
      case("NN")
         !
         !
         do io=1,Norb*Norb
            do jo=1,Norb*Norb
               !
               call get_element_from_stride(iorb,jorb,io)
               call get_element_from_stride(korb,lorb,jo)
               !
               ! U_(ij)(kl) = n_ij * n_kl
               Umatrix = Umatrix + kronecker_product(NabL(:,:,io),NabR(:,:,jo))*Utensor(iorb,jorb,korb,lorb)
               !
            enddo
         enddo
         !
         !
   end select
   deallocate(NabL,NabR,stride)
   !
   Utensor = 0d0
   do iorb=1,Norb !row of external submatrix
      do jorb=1,Norb !col of external submatrix
         do korb=1,Norb !row of internal submatrix
            do lorb=1,Norb !col of internal submatrix
               !
               ! this ordering reflects the kronecker_product
               io = korb + (iorb-1)*Norb
               jo = lorb + (jorb-1)*Norb
               !
               Utensor(iorb,jorb,korb,lorb) = dreal(Umatrix(io,jo))
               !
            enddo
         enddo
      enddo
   enddo
   deallocate(Umatrix)
   !
   contains
      !
      subroutine get_element_from_stride(i,j,ndx)
        implicit none
        integer,intent(in)                  :: ndx
        integer,intent(out)                 :: i,j
        !
        strideloop:do i=1,Norb
           do j=1,Norb
              if(stride(i,j)==ndx) exit strideloop
           enddo
        enddo strideloop
        !
      end subroutine get_element_from_stride
      !
end subroutine tensor_transform_NNNN_d
!
subroutine tensor_transform_NNNN_c(mode,Utensor,rotL,rotR)
   !
   implicit none
   character(len=2),intent(in)              :: mode
   complex(8),intent(inout)                 :: Utensor(:,:,:,:)
   complex(8),intent(in)                    :: rotL(:,:)
   complex(8),intent(in)                    :: rotR(:,:)
   !
   complex(8),allocatable                   :: Umatrix(:,:)
   complex(8),allocatable                   :: Nab(:,:,:),NabL(:,:,:),NabR(:,:,:)
   integer,allocatable                      :: stride(:,:)
   integer                                  :: Norb
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo,count
   !
   Norb = size(Utensor,dim=1)
   call assert_shape(Utensor,[Norb,Norb,Norb,Norb],"tensor_transform_NNNN_c","Utensor")
   call assert_shape(rotL,[Norb,Norb],"tensor_transform_NNNN_c","rotL")
   call assert_shape(rotR,[Norb,Norb],"tensor_transform_NNNN_c","rotR")
   !
   allocate(Nab(Norb,Norb,Norb*Norb));Nab=zero
   allocate(stride(Norb,Norb));stride=0
   count=0
   do iorb=1,Norb
      do jorb=1,Norb
         count=count+1
         Nab(jorb,iorb,count) = one !N=c^+,c = {col,row} --> transpose
         stride(iorb,jorb)=count
      enddo
   enddo
   !
   allocate(NabL(Norb,Norb,Norb*Norb));NabL=zero
   allocate(NabR(Norb,Norb,Norb*Norb));NabR=zero
   do io=1,Norb*Norb
      NabL(:,:,io) = rotate_zz(Nab(:,:,io),rotL)
      NabR(:,:,io) = rotate_zz(Nab(:,:,io),rotR)
   enddo
   deallocate(Nab)
   !
   allocate(Umatrix(Norb*Norb,Norb*Norb));Umatrix=zero
   select case(mode)
      case default
         !
         stop "tensor_transform_NNNN_d: Available modes: ZHER, TT."
         !
      case("GG")
         !
         !
         do io=1,Norb*Norb
            do jo=1,Norb*Norb
               !
               call get_element_from_stride(iorb,jorb,io)
               call get_element_from_stride(korb,lorb,jo)
               !
               ! U_(il)(jk) = G_ij * G_kl * U_ijkl
               Umatrix = Umatrix + kronecker_product(NabL(:,:,io),NabR(:,:,jo))*Utensor(iorb,lorb,jorb,korb)
               !
            enddo
         enddo
         !
         !
      case("NN")
         !
         !
         do io=1,Norb*Norb
            do jo=1,Norb*Norb
               !
               call get_element_from_stride(iorb,jorb,io)
               call get_element_from_stride(korb,lorb,jo)
               !
               ! U_(ij)(kl) = n_ij * n_kl
               Umatrix = Umatrix + kronecker_product(NabL(:,:,io),NabR(:,:,jo))*Utensor(iorb,jorb,korb,lorb)
               !
            enddo
         enddo
         !
         !
   end select
   deallocate(NabL,NabR,stride)
   !
   Utensor = zero
   do iorb=1,Norb !row of external submatrix
      do jorb=1,Norb !col of external submatrix
         do korb=1,Norb !row of internal submatrix
            do lorb=1,Norb !col of internal submatrix
               !
               ! this ordering reflects the kronecker_product
               io = korb + (iorb-1)*Norb
               jo = lorb + (jorb-1)*Norb
               !
               Utensor(iorb,jorb,korb,lorb) = Umatrix(io,jo)
               !
            enddo
         enddo
      enddo
   enddo
   deallocate(Umatrix)
   !
   contains
      !
      subroutine get_element_from_stride(i,j,ndx)
        implicit none
        integer,intent(in)                  :: ndx
        integer,intent(out)                 :: i,j
        !
        strideloop:do i=1,Norb
           do j=1,Norb
              if(stride(i,j)==ndx) exit strideloop
           enddo
        enddo strideloop
        !
      end subroutine get_element_from_stride
      !
end subroutine tensor_transform_NNNN_c


!------------------------------------------------------------------------------!
!PURPOSE: If the tensor is in some given square-matrix representation the user
!         must provide the representation mapping of the input matrix in the
!         form of an index list pointing at a given element of the tensor.
!         - Map(io,jo,1) --> iorb
!         - Map(io,jo,2) --> jorb
!         - Map(io,jo,3) --> korb
!         - Map(io,jo,4) --> lorb
!------------------------------------------------------------------------------!
subroutine tensor_transform_NN_d(mode,Umatrix,Map,rotL,rotR)
   !
   implicit none
   character(len=2),intent(in)              :: mode
   real(8),intent(inout)                    :: Umatrix(:,:)
   integer,intent(in)                       :: Map(:,:,:)
   complex(8),intent(in)                    :: rotL(:,:)
   complex(8),intent(in)                    :: rotR(:,:)
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
   call tensor_transform(mode,Utensor,rotL,rotR)
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
!
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
   call tensor_transform(mode,Utensor,rotL,rotR)
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
