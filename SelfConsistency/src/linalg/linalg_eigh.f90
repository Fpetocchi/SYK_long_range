

subroutine deigh(A,E,method,jobz,uplo,vl,vu,il,iu,tol)
  real(8),dimension(:,:),intent(inout)       :: A ! M v = E v/v(i,j) = ith component of jth vec.
  real(8),dimension(size(A,2)),intent(inout) :: E ! eigenvalues
  character(len=*),optional                  :: method
  character(len=1),optional                  :: jobz,uplo
  character(len=1)                           :: jobz_,uplo_,range
  character(len=20)                          :: method_
  real(8),optional                           :: vl,vu,tol
  integer,optional                           :: il,iu
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: info!i,j
  integer                                    :: lwork,liwork,mE
  integer                                    :: guess_liwork(1)
  real(8)                                    :: guess_lwork(1)
  logical                                    :: boolV,boolI
  real(8),dimension(:,:),allocatable         :: Z
  real(8),dimension(:),allocatable           :: work,rwork
  integer,dimension(:),allocatable           :: Isuppz,Ifail
  integer,dimension(:),allocatable           :: Iwork
  !
  method_='dsyevd';if(present(method))method_=trim(method)
  jobz_='V'  ;if(present(jobz))jobz_=jobz
  uplo_='U'  ;if(present(uplo))uplo_=uplo
  vl_  = 1d0 ;if(present(vL))vL_=vL
  vu_  = 1d0 ;if(present(vU))vU_=vU
  iL_  = 1   ;if(present(iL))iL_=iL
  iU_  = 1   ;if(present(iU))iU_=iU
  tol_ = 0d0 ;if(present(tol))tol_=tol
  !
  E=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  !
  if(jobz_/='V'.AND.jobz_/='N')stop "deigh_simple error: jobz has illegal value"
  if(uplo_/='U'.AND.uplo_/='L')stop "deigh_simple error: uplo has illegal value"
  !
  Ns = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "deigh_simple error: A has illegal shape"
  !
  select case(method_)
  case ("dsyevr")
     allocate(Isuppz(2*Ns))
     allocate(Z(Ns,Ns))
     call dsyevr(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,Isuppz,guess_lwork,-1,guess_liwork,-1,info)
     lwork = int(guess_lwork(1))
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call dsyevr(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,Isuppz,work,lwork,iwork,liwork,info)
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)
     !
  case ("dsyev")
     call dsyev(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,info)
     lwork = int(guess_lwork(1))
     allocate(work(lwork))
     call dsyev(jobz_,uplo_,Ns,A,Ns,E,work,lwork,info)
     !
  case default
     ! case("dsyevd")
     call dsyevd(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,guess_liwork,-1,info)
     lwork = int(guess_lwork(1))
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call dsyevd(jobz_,uplo_,Ns,A,Ns,E,work,lwork,iwork,liwork,info)
     !
  case ("dsyevx")
     allocate(Z(Ns,Ns))
     allocate(Ifail(Ns))
     allocate(rwork(7*Ns))
     allocate(iwork(5*Ns))
     call dsyevx(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,guess_lwork,-1,rwork,iwork,ifail,info)
     lwork = int(guess_lwork(1))
     allocate(work(lwork))
     call dsyevx(jobz_,range,uplo_,Ns,A,Ns,vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,work,lwork,rwork,iwork,ifail,info)
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)
  end select
  return
end subroutine deigh


subroutine zeigh(A,E,method,jobz,uplo,vl,vu,il,iu,tol)
  complex(8),dimension(:,:),intent(inout)       :: A ! M v = E v/v(i,j) = ith cmpt of jth vec.
  real(8),dimension(size(A,2)),intent(inout) :: E ! eigenvalues
  character(len=*),optional                  :: method
  character(len=1),optional                  :: jobz,uplo
  character(len=1)                           :: jobz_,uplo_,range
  character(len=20)                          :: method_
  real(8),optional                           :: vl,vu,tol
  integer,optional                           :: il,iu
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: info!i,j,
  integer                                    :: lwork,liwork,lrwork,mE
  complex(8)                                 :: guess_lwork(1)
  real(8)                                    :: guess_lrwork(1)
  integer                                    :: guess_liwork(1)
  logical                                    :: boolV,boolI
  complex(8),dimension(:,:),allocatable      :: Z
  complex(8),dimension(:),allocatable        :: Work
  real(8),dimension(:),allocatable           :: Rwork
  integer,dimension(:),allocatable           :: Iwork
  integer,dimension(:),allocatable           :: Isuppz
  integer,dimension(:),allocatable           :: Ifail
  !
  method_='zheevd';if(present(method))method_=trim(method)
  jobz_='V'  ;if(present(jobz))jobz_=jobz
  uplo_='U'  ;if(present(uplo))uplo_=uplo
  vl_  = 1d0 ;if(present(vL))vL_=vL
  vu_  = 1d0 ;if(present(vU))vU_=vU
  iL_  = 1   ;if(present(iL))iL_=iL
  iU_  = 1   ;if(present(iU))iU_=iU
  tol_ = 0d0 ;if(present(tol))tol_=tol
  !
  E=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  !
  if(jobz_/='V'.AND.jobz_/='N')stop "zeigh_simple error: jobz has illegal value"
  if(uplo_/='U'.AND.uplo_/='L')stop "zeigh_simple error: uplo has illegal value"
  !
  Ns = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "zeigh_simple error: A has illegal shape"
  !
  mE = Ns
  select case(method_)
  case ("zheevr")
     allocate(Isuppz(2*Ns))
     allocate(Z(Ns,Ns))
     call zheevr(jobz_,range,uplo_,&
          Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,&
          ME,E,Z,Ns,&
          Isuppz,guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,info)
     lwork = int(guess_lwork(1))
     lrwork= int(guess_lrwork(1))
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call zheevr(jobz_,range,uplo_,&
          Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,&
          ME,E,Z,Ns,&
          Isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
     !<copy the Evecs from Z to the input matrix A
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)
     !
  case ("zheev")
     allocate(rwork(max(1,3*Ns)))
     call zheev(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,rwork,info)
     lwork = int(guess_lwork(1))
     allocate(work(lwork))
     call zheev(jobz_,uplo_,Ns,A,Ns,E,work,lwork,rwork,info)
     !
  case default
     ! case("zheevd")
     call zheevd(jobz_,uplo_,Ns,A,Ns,E,guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,info)
     lwork = int(guess_lwork(1)) ; lrwork= int(guess_lrwork(1)) ; liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call zheevd(jobz_,uplo_,Ns,A,Ns,E,work,lwork,rwork,lrwork,iwork,liwork,info)
     !
  case ("zheevx")
     allocate(Z(Ns,Ns))
     allocate(Ifail(Ns))
     allocate(rwork(7*Ns))
     allocate(iwork(5*Ns))
     call zheevx(jobz_,range,uplo_,Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,guess_lwork,-1,rwork,iwork,ifail,info)
     lwork = int(guess_lwork(1))
     allocate(work(lwork))
     call zheevx(jobz_,range,uplo_,Ns,A,Ns,&
          vl_,vu_,iL_,iU_,tol_,mE,E,Z,Ns,work,lwork,rwork,iwork,ifail,info)
     if(jobz_=='V') A(:,1:mE) = Z(:,1:mE)
  end select
  return
end subroutine zeigh
