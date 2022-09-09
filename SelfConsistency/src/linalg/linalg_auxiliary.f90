!-------------------------------------------------------------------------------
!PURPOSE: compute the determinant of a real matrix using an LU factorization
!-------------------------------------------------------------------------------
function ddet(A) result(x)
  ! compute the determinant of a real matrix using an LU factorization
  real(8), intent(in)  :: A(:, :)
  real(8)              :: x
  integer              :: i
  integer              :: info, n
  integer, allocatable :: ipiv(:)
  real(8), allocatable :: At(:,:)
  n = size(A(1,:))
  call assert_shape(A, [n, n], "det", "A")
  allocate(At(n,n), ipiv(n))
  At = A
  call dgetrf(n, n, At, n, ipiv, info)
  if(info /= 0) then
     print *, "dgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop 'det: dgetrf error'
  end if
  ! At now contains the LU of the factorization A = PLU
  ! as L has unit diagonal entries, the determinant can be computed
  ! from the product of U's diagonal entries. Additional sign changes
  ! stemming from the permutations P have to be taken into account as well.
  x = 1d0
  do i = 1,n
     if(ipiv(i) /= i) then  ! additional sign change
        x = -x*At(i,i)
     else
        x = x*At(i,i)
     endif
  end do
end function ddet
!
function zdet(A) result(x)
  ! compute the determinant of a real matrix using an LU factorization
  complex(8), intent(in)  :: A(:, :)
  complex(8)              :: x
  integer                 :: i
  integer                 :: info, n
  integer, allocatable    :: ipiv(:)
  complex(8), allocatable :: At(:,:)
  n = size(A(1,:))
  call assert_shape(A, [n, n], "det", "A")
  allocate(At(n,n), ipiv(n))
  At = A
  call zgetrf(n, n, At, n, ipiv, info)
  if(info /= 0) then
     print *, "zgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop 'zdet error: zgetrf '
  end if
  ! for details on the computation, compare the comment in ddet().
  x = one
  do i = 1,n
     if(ipiv(i) /= i) then  ! additional sign change
        x = -x*At(i,i)
     else
        x = x*At(i,i)
     endif
  end do
end function zdet
!
function ddet3(A) result(x)
  real(8), intent(in)  :: A(3,3)
  real(8)              :: x
  x =     A(1,1)*A(2,2)*A(3,3)
  x = x + A(1,2)*A(2,3)*A(3,1)
  x = x + A(1,3)*A(2,1)*A(3,2)
  x = x - A(1,3)*A(2,2)*A(3,1)
  x = x - A(1,2)*A(2,1)*A(3,3)
  x = x - A(1,1)*A(2,3)*A(3,2)
end function ddet3
!
function zdet3(A) result(x)
  complex(8), intent(in)  :: A(3,3)
  complex(8)              :: x
  x =     A(1,1)*A(2,2)*A(3,3)
  x = x + A(1,2)*A(2,3)*A(3,1)
  x = x + A(1,3)*A(2,1)*A(3,2)
  x = x - A(1,3)*A(2,2)*A(3,1)
  x = x - A(1,2)*A(2,1)*A(3,3)
  x = x - A(1,1)*A(2,3)*A(3,2)
end function zdet3





!-------------------------------------------------------------------------------
!PURPOSE:  construct real matrix from diagonal elements
!-------------------------------------------------------------------------------
pure function ddiag(x) result(A)
  real(8), intent(in)  :: x(:)
  real(8), allocatable :: A(:,:)
  integer              :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = 0d0
  forall(i=1:n) A(i,i) = x(i)
end function ddiag
!
pure function zdiag(x) result(A)
  complex(8), intent(in)  :: x(:)
  complex(8), allocatable :: A(:,:)
  integer                 :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = zero
  forall(i=1:n) A(i,i) = x(i)
end function zdiag





!-------------------------------------------------------------------------------
!PURPOSE:  Create the hermitian conjugate of a matrix
!-------------------------------------------------------------------------------
function dag_d(A) result(Adag)
  real(8),intent(in)           :: A(:,:)
  real(8),allocatable          :: Adag(:,:)
  integer                      :: NA
  NA = size(A,dim=1) ; if(size(A,dim=2).ne.NA) stop "dag_d: input matrix not square."
  allocate(Adag(NA,NA));Adag=0d0
  Adag = transpose(A)
end function dag_d
function dag_z(A) result(Adag)
  complex(8),intent(in)        :: A(:,:)
  complex(8),allocatable       :: Adag(:,:)
  integer                      :: NA
  NA = size(A,dim=1) ; if(size(A,dim=2).ne.NA) stop "dag_z: input matrix not square."
  allocate(Adag(NA,NA));Adag=dcmplx(0d0,0d0)
  Adag = transpose(conjg(A))
end function dag_z




!-------------------------------------------------------------------------------
!PURPOSE:  Rotate matrix with a given rotation
!-------------------------------------------------------------------------------
function rotate_dz(A,U) result(Arot)
  real(8),intent(in)           :: A(:,:)
  complex(8),intent(in)        :: U(:,:)
  complex(8),allocatable       :: Arot(:,:)
  integer                      :: NA,NU
  NA = size(A,dim=1) ; if(size(A,dim=2).ne.NA) stop "rotate_dz: input matrix not square."
  NU = size(U,dim=1) ; if(size(U,dim=2).ne.NU) stop "rotate_dz: rotation matrix not square."
  allocate(Arot(NA,NA));Arot=dcmplx(0d0,0d0)
  Arot = matmul(transpose(conjg(U)),matmul(A,U))
end function rotate_dz
function rotate_zz(A,U) result(Arot)
  complex(8),intent(in)        :: A(:,:)
  complex(8),intent(in)        :: U(:,:)
  complex(8),allocatable       :: Arot(:,:)
  integer                      :: NA,NU
  NA = size(A,dim=1) ; if(size(A,dim=2).ne.NA) stop "rotate_zz: input matrix not square."
  NU = size(U,dim=1) ; if(size(U,dim=2).ne.NU) stop "rotate_zz: rotation matrix not square."
  allocate(Arot(NA,NA));Arot=dcmplx(0d0,0d0)
  Arot = matmul(transpose(conjg(U)),matmul(A,U))
end function rotate_zz
function rotate_dd(A,U) result(Arot)
  real(8),intent(in)           :: A(:,:)
  real(8),intent(in)           :: U(:,:)
  real(8),allocatable          :: Arot(:,:)
  integer                      :: NA,NU
  NA = size(A,dim=1) ; if(size(A,dim=2).ne.NA) stop "rotate_dd: input matrix not square."
  NU = size(U,dim=1) ; if(size(U,dim=2).ne.NU) stop "rotate_dd: rotation matrix not square."
  allocate(Arot(NA,NA));Arot=0d0
  Arot = matmul(transpose(U),matmul(A,U))
end function rotate_dd
function rotate_zd(A,U) result(Arot)
  complex(8),intent(in)        :: A(:,:)
  real(8),intent(in)           :: U(:,:)
  complex(8),allocatable       :: Arot(:,:)
  integer                      :: NA,NU
  NA = size(A,dim=1) ; if(size(A,dim=2).ne.NA) stop "rotate_zd: input matrix not square."
  NU = size(U,dim=1) ; if(size(U,dim=2).ne.NU) stop "rotate_zd: rotation matrix not square."
  allocate(Arot(NA,NA));Arot=dcmplx(0d0,0d0)
  Arot = matmul(transpose(U),matmul(A,U))
end function rotate_zd





!-------------------------------------------------------------------------------
!PURPOSE:  Rotate matrix with a given rotation
!-------------------------------------------------------------------------------
function diag_factor_d(A,fact) result(B)
  real(8),intent(in)           :: A(:,:)
  real(8),intent(in)           :: fact
  real(8),allocatable          :: B(:,:)
  integer                      :: i,NA
  NA = size(A,dim=1)
  if(size(A,dim=2).ne.NA) stop "diag_factor_d: input matrix not square."
  B=A
  do i=1,NA
     B(i,i)=fact*A(i,i)
  enddo
end function diag_factor_d
function diag_factor_z(A,fact) result(B)
  complex(8),intent(in)        :: A(:,:)
  real(8),intent(in)           :: fact
  complex(8),allocatable       :: B(:,:)
  integer                      :: i,NA
  NA = size(A,dim=1)
  if(size(A,dim=2).ne.NA) stop "diag_factor_z: input matrix not square."
  B=A
  do i=1,NA
     B(i,i)=fact*A(i,i)
  enddo
end function diag_factor_z






!-------------------------------------------------------------------------------
!PURPOSE:  return the diagonal of a matrix
!-------------------------------------------------------------------------------
pure function d_diagonal(A) result(dd)
  real(8),intent(in)           :: A(:,:)
  real(8),dimension(size(A,1)) :: dd
  integer                      :: i
  do i = 1,size(A,1)
     dd(i) = A(i,i)
  end do
end function d_diagonal

pure function z_diagonal(A) result(dd)
  complex(8),intent(in)           :: A(:,:)
  complex(8),dimension(size(A,1)) :: dd
  integer                         :: i
  do i = 1,size(A,1)
     dd(i) = A(i,i)
  end do
end function z_diagonal





!-------------------------------------------------------------------------------
!PURPOSE:  return trace along the main diagonal
!-------------------------------------------------------------------------------
pure function dtrace(A) result(t)
  real(8), intent(in) :: A(:,:)
  real(8)             :: t
  integer             :: i
  t = 0d0
  do i = 1,minval(shape(A))
     t = t + A(i,i)
  end do
end function dtrace

pure function ztrace(A) result(t)
  complex(8), intent(in) :: A(:,:)
  complex(8)             :: t
  integer                :: i
  t = zero
  do i = 1,minval(shape(A))
     t = t + A(i,i)
  end do
end function ztrace





!-------------------------------------------------------------------------------
!PURPOSE:  Returns the identity matrix of size n x n and type real.
!-------------------------------------------------------------------------------
pure function deye(n) result(A)
  integer, intent(in) :: n
  real(8)             :: A(n, n)
  integer             :: i
  A = 0d0
  do i = 1, n
     A(i,i) = 1d0
  end do
end function deye
!
pure function zeye(n) result(A)
  integer, intent(in) :: n
  complex(8)          :: A(n, n)
  integer             :: i
  A = zero
  do i = 1, n
     A(i,i) = one
  end do
end function zeye





!-------------------------------------------------------------------------------
!PURPOSE: Cramer
!-------------------------------------------------------------------------------
function Cramer_3_d(x,y,z,b) result(M)
   implicit none
   real(8),intent(in)      :: x(3),y(3),z(3),b(3)
   real(8)                 :: M(3)
   integer                 :: i
   real(8)                 :: num(3,3),den(3,3)
   !
   den(:,1) = x
   den(:,2) = y
   den(:,3) = z
   if(ddet3(den).eq.0d0)stop"Cramer_3_d: denominator"
   !
   do i=1,3
      num = den
      num(:,i) = b
      M(i) = ddet3(num) / ddet3(den)
   enddo
   !
end function Cramer_3_d
!
function Cramer_3_z(x,y,z,b) result(M)
   implicit none
   complex(8),intent(in)   :: x(3),y(3),z(3),b(3)
   complex(8)              :: M(3)
   integer                 :: i
   complex(8)              :: num(3,3),den(3,3)
   !
   den(:,1) = x
   den(:,2) = y
   den(:,3) = z
   if(zdet3(den).eq.0d0)stop"Cramer_3_z: denominator"
   !
   do i=1,3
      num = den
      num(:,i) = b
      M(i) = zdet3(num) / zdet3(den)
   enddo
   !
end function Cramer_3_z
!
function Cramer_2_d(x,y,b) result(M)
   implicit none
   real(8),intent(in)      :: x(2),y(2),b(2)
   real(8)                 :: M(2)
   integer                 :: i
   real(8)                 :: num(2,2),den(2,2),detden
   !
   den(:,1) = x
   den(:,2) = y
   detden = den(1,1)*den(2,2) - den(1,2)*den(2,1)
   if(detden.eq.0d0)stop"Cramer_2_d: denominator"
   !
   do i=1,2
      num = den
      num(:,i) = b
      M(i) = (num(1,1)*num(2,2) - num(1,2)*num(2,1)) / detden
   enddo
   !
end function Cramer_2_d
!
function Cramer_2_z(x,y,b) result(M)
   implicit none
   complex(8),intent(in)   :: x(2),y(2),b(2)
   complex(8)              :: M(2)
   integer                 :: i
   complex(8)              :: num(2,2),den(2,2),detden
   !
   den(:,1) = x
   den(:,2) = y
   detden = den(1,1)*den(2,2) - den(1,2)*den(2,1)
   if(detden.eq.0d0)stop"Cramer_2_z: denominator"
   !
   do i=1,2
      num = den
      num(:,i) = b
      M(i) = (num(1,1)*num(2,2) - num(1,2)*num(2,1)) / detden
   enddo
   !
end function Cramer_2_z






!-------------------------------------------------------------------------------
!PURPOSE:  Returns an array of zeros of specified size from 1 to 7 dimension
!-------------------------------------------------------------------------------
pure function zzeros_1(n) result(A)
  integer, intent(in) :: n
  complex(8)          :: A(n)
  A = zero
end function zzeros_1
!
pure function zzeros_2(n1,n2) result(A)
  integer, intent(in) :: n1,n2
  complex(8)          :: A(n1,n2)
  A = zero
end function zzeros_2
!
pure function zzeros_3(n1,n2,n3) result(A)
  integer, intent(in) :: n1,n2,n3
  complex(8)          :: A(n1,n2,n3)
  A = zero
end function zzeros_3
!
pure function zzeros_4(n1,n2,n3,n4) result(A)
  integer, intent(in) :: n1,n2,n3,n4
  complex(8)          :: A(n1,n2,n3,n4)
  A = zero
end function zzeros_4
!
pure function zzeros_5(n1,n2,n3,n4,n5) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5
  complex(8)          :: A(n1,n2,n3,n4,n5)
  A = zero
end function zzeros_5
!
pure function zzeros_6(n1,n2,n3,n4,n5,n6) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6
  complex(8)          :: A(n1,n2,n3,n4,n5,n6)
  A = zero
end function zzeros_6
!
pure function zzeros_7(n1,n2,n3,n4,n5,n6,n7) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6,n7
  complex(8)          :: A(n1,n2,n3,n4,n5,n6,n7)
  A = zero
end function zzeros_7



pure function zones_1(n) result(A)
  integer, intent(in) :: n
  complex(8)          :: A(n)
  A = one
end function zones_1
!
pure function zones_2(n1,n2) result(A)
  integer, intent(in) :: n1,n2
  complex(8)          :: A(n1,n2)
  A = one
end function zones_2
!
pure function zones_3(n1,n2,n3) result(A)
  integer, intent(in) :: n1,n2,n3
  complex(8)          :: A(n1,n2,n3)
  A = one
end function zones_3
!
pure function zones_4(n1,n2,n3,n4) result(A)
  integer, intent(in) :: n1,n2,n3,n4
  complex(8)          :: A(n1,n2,n3,n4)
  A = one
end function zones_4
!
pure function zones_5(n1,n2,n3,n4,n5) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5
  complex(8)          :: A(n1,n2,n3,n4,n5)
  A = one
end function zones_5
!
pure function zones_6(n1,n2,n3,n4,n5,n6) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6
  complex(8)          :: A(n1,n2,n3,n4,n5,n6)
  A = one
end function zones_6
!
pure function zones_7(n1,n2,n3,n4,n5,n6,n7) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6,n7
  complex(8)          :: A(n1,n2,n3,n4,n5,n6,n7)
  A = one
end function zones_7
