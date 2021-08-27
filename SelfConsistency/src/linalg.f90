module linalg

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface eig
      module procedure deig
      module procedure zeig
   end interface eig

   interface eigh
      module procedure deigh
      module procedure zeigh
   end interface eigh

   interface inv
      module procedure dinv
      module procedure zinv
   end interface inv

   interface inv_sym
      module procedure dinv_sym
      module procedure zinv_sym
   end interface inv_sym

   interface inv_her
      module procedure zinv_her
   end interface inv_her

   interface det
      module procedure ddet
      module procedure zdet
   end interface det

   interface det3
      module procedure ddet3
      module procedure zdet3
   end interface det3

   interface diag
      module procedure ddiag
      module procedure zdiag
   end interface diag

   interface diagonal
      module procedure d_diagonal
      module procedure z_diagonal
   end interface diagonal

   interface trace
      module procedure dtrace
      module procedure ztrace
   end interface trace

   interface rotate
      module procedure rotate_d
      module procedure rotate_z
   end interface rotate

   interface Cramer
      module procedure Cramer_2_d
      module procedure Cramer_2_z
      module procedure Cramer_3_d
      module procedure Cramer_3_z
   end interface Cramer

   interface zeros
      module procedure zzeros_1
      module procedure zzeros_2
      module procedure zzeros_3
      module procedure zzeros_4
      module procedure zzeros_5
      module procedure zzeros_6
      module procedure zzeros_7
   end interface zeros

   interface ones
      module procedure zones_1
      module procedure zones_2
      module procedure zones_3
      module procedure zones_4
      module procedure zones_5
      module procedure zones_6
      module procedure zones_7
   end interface ones

   interface diag_factor
      module procedure diag_factor_d
      module procedure diag_factor_z
   end interface diag_factor

   interface kronecker_product
      module procedure i_kronecker_product
      module procedure d_kronecker_product
      module procedure c_kronecker_product
   end interface kronecker_product

   interface outerprod
      module procedure outerprod_d,outerprod_c
   end interface outerprod

   interface cross_product
      module procedure cross_3d_d
      module procedure cross_3d_c
   end interface cross_product

   interface s3_product
      module procedure s3_product_d
      module procedure s3_product_c
   end interface s3_product

   interface tensor_transform
      module procedure tensor_transform_NNNN_d
      module procedure tensor_transform_NNNN_c
      module procedure tensor_transform_NN_d
      module procedure tensor_transform_NN_c
   end interface tensor_transform

   !NOT PUBLIC:
   interface assert_shape
      module procedure iassert_shape_2
      module procedure dassert_shape_2
      module procedure zassert_shape_2
      module procedure iassert_shape_3
      module procedure dassert_shape_3
      module procedure zassert_shape_3
      module procedure iassert_shape_4
      module procedure dassert_shape_4
      module procedure zassert_shape_4
   end interface assert_shape
   interface
      integer function ilaenv( ispec, name, opts, n1, n2, n3, n4 )
        character*(*) name, opts
        integer       ispec, n1, n2, n3, n4
      end function ilaenv
   end interface

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   complex(8),parameter,private             :: zero=(0.d0,0.d0)
   complex(8),parameter,private             :: xi=(0.d0,1.d0)
   complex(8),parameter,private             :: one=(1.d0,0.d0)

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: eig                                                                ![GenericMatrix,Evalues,Evectors])
   public :: eigh                                                               ![HermitainMatrix(replaced by Evectors),Evalues])
   public :: inv                                                                ![GenericMatrix(replaced with inverse)])
   public :: inv_sym                                                            ![SymmetricMatrix(replaced with inverse),uplo(optional for upper or lower)])
   public :: inv_her                                                            ![HermitainMatrix(replaced with inverse),uplo(optional for upper or lower)])
   public :: tensor_transform                                                   ![Tensor, either (Norb,Norb,Norb,Norb) or (Norb**2,Norb**2), Rotation(Norb,Norb)]
   !functions
   public :: det, det3, trace                                                   !numpy-like. Takes matrix returns scalar.
   public :: Cramer                                                             !numpy-like. Takes matrix returns solution array.
   public :: deye, zeye, zeros, ones                                            !numpy-like. Takes integer returns square matrix.
   public :: diag                                                               !numpy-like. Takes vector returns diagonal square matrix.
   public :: diagonal                                                           !numpy-like. Takes square matrix returns diagonal vector.
   public :: rotate                                                             !numpy-like. Takes two square matrix returns square matrix.
   public :: diag_factor                                                        !numpy-like. Takes square matrix and real factor returns square matrix.
   public :: kronecker_product                                                  !Takes two square matrix (size Norb) returns square matrix (size Norb**2).
   !
   public :: outerprod                                                          !Form a matrix A(:,:) from the outerproduct of two 1d arrays: A(i,j) = a_i*b_j
   public :: cross_product                                                      !cross or vector product for 2d and 3d vectors.
   public :: s3_product                                                         !evaluate the S3 product A.(BxC) for 3d vectors

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE:  SOLUTION TO EIGEN PROBLEMS
   ! - general matrices (in general non-symmetric or non-complex-hermitian matrices).
   ! - real symmetric/complex hermitian matrices
   !---------------------------------------------------------------------------!
   include "linalg/linalg_eig.f90"
   include "linalg/linalg_eigh.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE: INVERSION OF A MATRIX USING LAPACK LIBRARY
   ! - General M*N (D,C)
   ! - Symmetric N*N (D,C)
   ! - Hermitial (C)
   ! note: M is destroyed and replaces by its inverse M^-1
   !---------------------------------------------------------------------------!
   include "linalg/linalg_inv.f90"
   include "linalg/linalg_inv_sym.f90"
   include "linalg/linalg_inv_her.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE: AUXILIARY AND COMPUTATIONAL ROUTINES
   ! - det: compute determinant of a real/complex matrix
   ! - diag: construct real matrix from diagonal elements
   ! - trace: return trace along the main diagonal
   ! - Xeye: returns the identity matrix of size n x n
   ! - Rotate: returns the rotated matrix of size n x n
   !---------------------------------------------------------------------------!
   include "linalg/linalg_auxiliary.f90"


   !+--------------------------------------------------------------------------!
   !PURPOSE: EXTERNAL PRODUCTS ROUTINES
   ! - kronecker:  compute the tensor product (M1_kp_M2) of
   ! two complex matrices M1 and M2. nr1(nr2) and nc1(nc2) are
   ! the number of rows and columns of the Matrix M1 and M2
   ! - outerprod: Form a matrix A(:,:) from the outerproduct of two 1d arrays:
   ! A(i,j) = a_i*b_j
   ! - cross: cross or vector product for 2d and 3d vectors.
   ! - s3_product: evaluate the S3 product A.(BxC) for 3d vectors
   !+--------------------------------------------------------------------------!
   include "linalg/linalg_external_products.f90"


   !+--------------------------------------------------------------------------!
   !PURPOSE: TENSOR ROUTINES
   ! - tensor_transform: rotate a tensor that's given by the kroenecker product
   ! of two square matrices into a different basis
   !+--------------------------------------------------------------------------!
   include "linalg/linalg_tensor.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE: Routines needed to make the module self-sufficient
   !---------------------------------------------------------------------------!
   subroutine iassert_shape_2(A, shap, routine, matname)
     integer,intent(in)  :: A(:,:)
     integer,intent(in)  :: shap(:)
     character(len=*)    :: routine, matname
     if(any(shape(A) /= shap)) then
        print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
        print*, "Shape should be ", shap
        stop "Aborting due to illegal matrix operation"
     end if
  end subroutine iassert_shape_2
  subroutine dassert_shape_2(A, shap, routine, matname)
    real(8),intent(in)  :: A(:,:)
    integer,intent(in)  :: shap(:)
    character(len=*)    :: routine, matname
    if(any(shape(A) /= shap)) then
      print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
      print*, "Shape should be ", shap
      stop "Aborting due to illegal matrix operation"
    end if
   end subroutine dassert_shape_2
   subroutine zassert_shape_2(A, shap, routine, matname)
     complex(8),intent(in) :: A(:,:)
     integer,intent(in)    :: shap(:)
     character(len=*)      :: routine, matname
     if(any(shape(A) /= shap)) then
        print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
        print*, "Shape should be ", shap
        stop "Aborting due to illegal matrix operation"
     end if
  end subroutine zassert_shape_2
  subroutine iassert_shape_3(A, shap, routine, matname)
    integer,intent(in)  :: A(:,:,:)
    integer,intent(in)  :: shap(:)
    character(len=*)    :: routine, matname
    if(any(shape(A) /= shap)) then
      print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
      print*, "Shape should be ", shap
      stop "Aborting due to illegal matrix operation"
    end if
end subroutine iassert_shape_3
  subroutine dassert_shape_3(A, shap, routine, matname)
    real(8),intent(in)  :: A(:,:,:)
    integer,intent(in)  :: shap(:)
    character(len=*)    :: routine, matname
    if(any(shape(A) /= shap)) then
      print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
      print*, "Shape should be ", shap
      stop "Aborting due to illegal matrix operation"
    end if
  end subroutine dassert_shape_3
  subroutine zassert_shape_3(A, shap, routine, matname)
    complex(8),intent(in) :: A(:,:,:)
    integer,intent(in)    :: shap(:)
    character(len=*)      :: routine, matname
    if(any(shape(A) /= shap)) then
      print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
      print*, "Shape should be ", shap
      stop "Aborting due to illegal matrix operation"
    end if
  end subroutine zassert_shape_3
  subroutine iassert_shape_4(A, shap, routine, matname)
     integer,intent(in)  :: A(:,:,:,:)
     integer,intent(in)  :: shap(:)
     character(len=*)    :: routine, matname
     if(any(shape(A) /= shap)) then
      print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
      print*, "Shape should be ", shap
      stop "Aborting due to illegal matrix operation"
     end if
  end subroutine iassert_shape_4
  subroutine dassert_shape_4(A, shap, routine, matname)
     real(8),intent(in)  :: A(:,:,:,:)
     integer,intent(in)  :: shap(:)
     character(len=*)    :: routine, matname
     if(any(shape(A) /= shap)) then
      print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
      print*, "Shape should be ", shap
      stop "Aborting due to illegal matrix operation"
     end if
  end subroutine dassert_shape_4
  subroutine zassert_shape_4(A, shap, routine, matname)
    complex(8),intent(in) :: A(:,:,:,:)
    integer,intent(in)    :: shap(:)
    character(len=*)      :: routine, matname
    if(any(shape(A) /= shap)) then
      print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
      print*, "Shape should be ", shap
      stop "Aborting due to illegal matrix operation"
    end if
  end subroutine zassert_shape_4




end module linalg
