module linalg

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

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

   !NOT PUBLIC:
   interface assert_shape
      module procedure dassert_shape
      module procedure zassert_shape
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
   !functions
   public :: det,trace                                                          !numpy-like. Takes matrix returns scalar.
   public :: deye, zeye, zeros, ones                                            !numpy-like. Takes integer returns square matrix.
   public :: diag                                                               !numpy-like. Takes vector returns diagonal sqaure matrix.
   public :: diagonal                                                           !numpy-like. Takes sqaure matrix returns diagonal vector.
   public :: rotate                                                             !numpy-like. Takes two sqaure matrix returns sqaure matrix.

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
   !PURPOSE  : AUXILIARY AND COMPUTATIONAL ROUTINES
   ! - det: compute determinant of a real/complex matrix
   ! - diag: construct real matrix from diagonal elements
   ! - trace: return trace along the main diagonal
   ! - Xeye: returns the identity matrix of size n x n
   ! - Rotate: returns the rotated matrix of size n x n
   !---------------------------------------------------------------------------!
   include "linalg/linalg_auxiliary.f90"


   !---------------------------------------------------------------------------!
   !PURPOSE: Routines needed to make the module self-sufficient
   !---------------------------------------------------------------------------!
   subroutine dassert_shape(A, shap, routine, matname)
     real(8),intent(in)  :: A(:,:)
     integer,intent(in)  :: shap(:)
     character(len=*)    :: routine, matname
     if(any(shape(A) /= shap)) then
        print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
        print*, "Shape should be ", shap
        stop "Aborting due to illegal matrix operation"
     end if
   end subroutine dassert_shape
   subroutine zassert_shape(A, shap, routine, matname)
     complex(8),intent(in) :: A(:,:)
     integer,intent(in)    :: shap(:)
     character(len=*)      :: routine, matname
     if(any(shape(A) /= shap)) then
        print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
        print*, "Shape should be ", shap
        stop "Aborting due to illegal matrix operation"
     end if
   end subroutine zassert_shape

end module linalg
