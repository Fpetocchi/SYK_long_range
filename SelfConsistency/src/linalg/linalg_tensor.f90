!------------------------------------------------------------------------------!
!PURPOSE: Rotate a tensor diven by the kroenecker product of square matrices
!------------------------------------------------------------------------------!
subroutine tensor_transform_NNNN_d(Utensor,rot,onlyNaNb)
   !
   implicit none
   real(8),intent(inout)                    :: Utensor(:,:,:,:)
   complex(8),intent(in)                    :: rot(:,:)
   logical,intent(in),optional              :: onlyNaNb
   !
   real(8),allocatable                      :: Umatrix(:,:)
   complex(8),allocatable                   :: Nab(:,:,:)
   complex(8),allocatable                   :: rotdag(:,:)
   integer,allocatable                      :: stride(:,:)
   integer                                  :: Norb
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo,count
   logical                                  :: onlyNaNb_
   !
   Norb = size(Utensor,dim=1)
   call assert_shape(Utensor,[Norb,Norb,Norb,Norb],"tensor_transform_NNNN_d","Utensor")
   if(size(rot,1).ne.Norb) stop "tensor_transform_NNNN_d: Rotation has a wrong dimension."
   !
   onlyNaNb_=.false.
   if(present(onlyNaNb))onlyNaNb_=onlyNaNb
   !
   allocate(Umatrix(Norb*Norb,Norb*Norb));Umatrix=0d0
   allocate(Nab(Norb,Norb,Norb*Norb));Nab=zero
   allocate(stride(Norb,Norb));stride=0
   !
   allocate(rotdag(Norb,Norb));rotdag=zero
   rotdag=transpose(conjg(rot))
   !
   count=0
   do iorb=1,Norb
      do jorb=1,Norb
         count=count+1
         Nab(jorb,iorb,count) = one  !N=c^+,c = {col,row} --> transpose
         stride(iorb,jorb)=count
      enddo
   enddo
   !
   do io=1,Norb*Norb
      Nab(:,:,io) = matmul(rotdag,matmul(Nab(:,:,io),rot))
   enddo
   !
   do io=1,Norb*Norb
      do jo=1,Norb*Norb
         !
         call get_element_from_stride(iorb,jorb,io)
         call get_element_from_stride(korb,lorb,jo)
         !
         Umatrix = Umatrix + dreal(kronecker_product(Nab(:,:,io),Nab(:,:,jo)))*Utensor(iorb,jorb,korb,lorb)
         !
      enddo
   enddo
   !
   do iorb=1,Norb !row of external submatrix
      do jorb=1,Norb !col of external submatrix
         if(onlyNaNb_.and.(iorb.ne.iorb))cycle
         do korb=1,Norb !row of internal submatrix
            do lorb=1,Norb !col of internal submatrix
               if(onlyNaNb_.and.(korb.ne.lorb))cycle
               !
               io = korb + (iorb-1)*Norb
               jo = lorb + (jorb-1)*Norb
               !
               Utensor(iorb,jorb,korb,lorb) = Umatrix(io,jo)
               !
            enddo
         enddo
      enddo
   enddo
   deallocate(Umatrix,Nab,stride)
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
subroutine tensor_transform_NNNN_c(Utensor,rot,onlyNaNb)
   !
   implicit none
   complex(8),intent(inout)                 :: Utensor(:,:,:,:)
   complex(8),intent(in)                    :: rot(:,:)
   logical,intent(in),optional              :: onlyNaNb
   !
   real(8),allocatable                      :: Umatrix(:,:)
   complex(8),allocatable                   :: Nab(:,:,:)
   complex(8),allocatable                   :: rotdag(:,:)
   integer,allocatable                      :: stride(:,:)
   integer                                  :: Norb
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo,count
   logical                                  :: onlyNaNb_
   !
   Norb = size(Utensor,dim=1)
   call assert_shape(Utensor,[Norb,Norb,Norb,Norb],"tensor_transform_NNNN_c","Utensor")
   if(size(rot,1).ne.Norb) stop "tensor_transform_NNNN_c: Rotation has a wrong dimension."
   !
   onlyNaNb_=.false.
   if(present(onlyNaNb))onlyNaNb_=onlyNaNb
   !
   allocate(Umatrix(Norb*Norb,Norb*Norb));Umatrix=zero
   allocate(Nab(Norb,Norb,Norb*Norb));Nab=zero
   allocate(stride(Norb,Norb));stride=0
   !
   allocate(rotdag(Norb,Norb));rotdag=zero
   rotdag=transpose(conjg(rot))
   !
   count=0
   do iorb=1,Norb
      do jorb=1,Norb
         count=count+1
         Nab(jorb,iorb,count) = one  !N=c^+,c = {col,row} --> transpose
         stride(iorb,jorb)=count
      enddo
   enddo
   !
   do io=1,Norb*Norb
      Nab(:,:,io) = matmul(rotdag,matmul(Nab(:,:,io),rot))
   enddo
   !
   do io=1,Norb*Norb
      do jo=1,Norb*Norb
         !
         call get_element_from_stride(iorb,jorb,io)
         call get_element_from_stride(korb,lorb,jo)
         !
         Umatrix = Umatrix + kronecker_product(Nab(:,:,io),Nab(:,:,jo))*Utensor(iorb,jorb,korb,lorb)
         !
      enddo
   enddo
   !
   do iorb=1,Norb !row of external submatrix
      do jorb=1,Norb !col of external submatrix
         if(onlyNaNb_.and.(iorb.ne.iorb))cycle
         do korb=1,Norb !row of internal submatrix
            do lorb=1,Norb !col of internal submatrix
               if(onlyNaNb_.and.(korb.ne.lorb))cycle
               !
               io = korb + (iorb-1)*Norb
               jo = lorb + (jorb-1)*Norb
               !
               Utensor(iorb,jorb,korb,lorb) = Umatrix(io,jo)
               !
            enddo
         enddo
      enddo
   enddo
   deallocate(Umatrix,Nab,stride)
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
subroutine tensor_transform_NN_d(Umatrix,Map,rot,onlyNaNb)
   !
   implicit none
   real(8),intent(inout)                    :: Umatrix(:,:)
   integer,intent(in)                       :: Map(:,:,:)
   complex(8),intent(in)                    :: rot(:,:)
   logical,intent(in),optional              :: onlyNaNb
   !
   real(8),allocatable                      :: Utensor(:,:,:,:)
   integer                                  :: Norb,Nbp
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo
   logical                                  :: onlyNaNb_
   !
   Nbp = size(Umatrix,dim=1)
   call assert_shape(Umatrix,[Nbp,Nbp],"tensor_transform_NN_d","Umatrix")
   call assert_shape(Map,[Nbp,Nbp,4],"tensor_transform_NN_d","Map")
   !
   Norb = int(sqrt(dble(Nbp)))
   if(size(rot,1).ne.Norb) stop "tensor_transform_NN_d: Rotation has a wrong dimension."
   !
   onlyNaNb_=.false.
   if(present(onlyNaNb))onlyNaNb_=onlyNaNb
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
   call tensor_transform(Utensor,rot,onlyNaNb=onlyNaNb_)
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
subroutine tensor_transform_NN_c(Umatrix,Map,rot,onlyNaNb)
   !
   implicit none
   complex(8),intent(inout)                 :: Umatrix(:,:)
   integer,intent(in)                       :: Map(:,:,:)
   complex(8),intent(in)                    :: rot(:,:)
   logical,intent(in),optional              :: onlyNaNb
   !
   complex(8),allocatable                   :: Utensor(:,:,:,:)
   integer                                  :: Norb,Nbp
   integer                                  :: iorb,jorb,korb,lorb
   integer                                  :: io,jo
   logical                                  :: onlyNaNb_
   !
   Nbp = size(Umatrix,dim=1)
   call assert_shape(Umatrix,[Nbp,Nbp],"tensor_transform_NN_c","Umatrix")
   call assert_shape(Map,[Nbp,Nbp,4],"tensor_transform_NN_c","Map")
   !
   Norb = int(sqrt(dble(Nbp)))
   if(size(rot,1).ne.Norb) stop "tensor_transform_NN_c: Rotation has a wrong dimension."
   !
   onlyNaNb_=.false.
   if(present(onlyNaNb))onlyNaNb_=onlyNaNb
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
   call tensor_transform(Utensor,rot,onlyNaNb=onlyNaNb_)
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
end subroutine tensor_transform_NN_c
