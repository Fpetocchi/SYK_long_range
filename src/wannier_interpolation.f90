module wannier_interpolation

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   ! qui lo devo mettere 'input path? se no sbatto tutto u '
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),parameter,private                :: pi=3.14159265358979323846d0
   !
   real(8),allocatable                      :: rsite(:,:)
   real(8)                                  :: rlat(3,3)
   real(8),private                          :: lat(3,3)
   real(8),private                          :: vol
   real(8),private                          :: rvol
   !
   real(8),parameter,private                :: epsWig=1e-5
   integer,private                          :: Nwig
   integer,allocatable,private              :: rvecwig(:,:)
   integer,allocatable,private              :: nrdegwig(:)
   !
   logical,private                          :: Lat_stored=.false.               !Internal flag for routines that need rlat
   logical,private                          :: Wig_stored=.false.               !Internal flag for routines performing Wannier interpolation

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: wannierinterpolation_matrix                                        ![Kpt_orig(3,Nkpt_orig),Kpt_intp(3,Nkpt_intp),Mat_orig(n,n,Npoins,Nkpt_orig),Mat_intp(n,n,Npoins,Nkpt_intp)]
   public :: wannierinterpolation_WNN

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the Lattice vectors
   !---------------------------------------------------------------------------!
   subroutine read_lattice(pathINPUT)
      !
      use utils_misc
      use linalg
      implicit none
      !
      character(len=*),intent(in)           :: pathINPUT
      character(len=256)                    :: path
      integer                               :: unit
      logical                               :: filexists
      !
      !
      write(*,*) "--- read_lattice ---"
      !
      !
      ! Look for LATTC
      path=pathINPUT//"LATTC"
      call inquireFile(reg(path),filexists)
      !
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="old",position="rewind",action="read")
      read(unit,*)
      read(unit,*) lat(1:3,1)
      read(unit,*) lat(1:3,2)
      read(unit,*) lat(1:3,3)
      close(unit)
      !
      rlat = lat
      call inv_sym(rlat)
      !
      rlat = 2*pi*transpose(rlat)
      !
      vol = det(lat)
      rvol = 8*pi**3 / vol
      write(*,*)"Unit cell volume=",vol
      !
      Lat_stored=.true.
      !
   end subroutine read_lattice


   !---------------------------------------------------------------------------!
   !PURPOSE: calculate lattice points inside Wigner-Seitz cell of given supercell
   !NEEDED:
   ! lat(3,3) : premitive lattice vectors
   ! nkpt : number of k-points
   ! nkpt3(3) : number of k-points along each direction
   !OUTPUT:
   ! nwig : number of points inside wigner-seitz supercell
   ! rvec(3,2*nkpt) : (INTEGER) lattice vectors
   ! nrdeg(2*nkpt) : degeneracy
   !---------------------------------------------------------------------------!
   subroutine calc_wignerseiz(nkpt,nkpt3)
      !
      use utils_misc
      implicit none
      !
      integer,intent(in)           :: nkpt,nkpt3(3)
      integer                      :: ir1,ir2,ir3,irsc1,irsc2,irsc3
      integer,parameter            :: nshell=2
      double precision             :: rtmp(3),rtmpsc(3),dr(3)
      integer                      :: i,i0
      double precision,allocatable :: dist(:)
      double precision             :: distmin
      !
      !
      write(*,*) "--- calc_wignerseiz ---"
      if(.not.Lat_stored)stop "Lattice positions not stored. Either call read_lattice(path) or read_Hk(path,Hk,kpt)"
      !
      !
      allocate(rvecwig(3,10*nkpt));rvecwig=0
      allocate(nrdegwig(10*nkpt));nrdegwig=0
      allocate(dist((2*nshell+1)**3))
      ! this i0 corresponds to irsc1=irsc2=irsc3=0
      i0=nshell*(1+(2*nshell+1)*(1+(2*nshell+1)))+1
      !
      nwig=0
      do ir1=-nkpt3(1),+nkpt3(1)
         do ir2=-nkpt3(2),+nkpt3(2)
            do ir3=-nkpt3(3),+nkpt3(3)
               rtmp(:)=matmul(lat,(/ir1,ir2,ir3/))
               i=0
               !
               do irsc1=-nshell,+nshell
                  do irsc2=-nshell,+nshell
                     do irsc3=-nshell,+nshell
                        i=i+1
                        rtmpsc(:)=matmul(lat,(/nkpt3(1)*irsc1,nkpt3(2)*irsc2,nkpt3(3)*irsc3/))
                        dr(:)=rtmp(:)-rtmpsc(:)
                        dist(i)=sum(dr(:)**2)
                     enddo ! irsc3
                  enddo ! irsc2
               enddo ! irsc1
               !
               distmin=minval(dist(:))
               if (abs(distmin-dist(i0)).le.epsWig) then
                  nwig=nwig+1
                  if (nwig.gt.10*nkpt) stop "nwig>10*nkpt"
                  rvecwig(:,nwig)=(/ir1,ir2,ir3/)
                  nrdegwig(nwig)=count(abs(distmin-dist(:)).le.epsWig)
                  !write(*,*) nwig,rvecwig(:,nwig),nrdegwig(nwig)
               endif
               !
            enddo
         enddo
      enddo
      deallocate(dist)
      !
      if (abs(sum(1d0/nrdegwig(1:nwig))-nkpt).gt.epsWig) then
         write(*,*) "Error: sum(1/nrdeg(:))=",sum(1d0/nrdegwig(1:nwig))
         stop "nrdeg"
      endif
      !
      Wig_stored=.true.
      !
   end subroutine calc_wignerseiz


   !----------------------------------------------------------------------------!
   !PURPOSE:
   !----------------------------------------------------------------------------!
   subroutine wannierinterpolation_matrix(kpt,kpt_intp,mat,mat_intp)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: kpt_orig(:,:)
      real(8),intent(in)                    :: kpt_intp(:,:)
      complex(8),intent(in)                 :: mat_orig(:,:,:,:)
      complex(8),intent(inout)              :: mat_intp(:,:,:,:)
      !
      integer                               :: Nkpt_orig,Nkpt_intp
      integer                               :: Nsize,Npoints
      integer                               :: ik,ir,id
      real(8)                               :: kR
      complex(8)                            :: cfac
      complex(8),allocatable                :: matofR(:,:,:,:)
      !
      !
      write(*,*) "--- wannierinterpolation_matrix ---"
      if (.not.lwintp_init) stop 'init'
      if (nkpt.ne.size(nrdeg)/2) stop 'nkpt'
      !
      !
      ! Size checks on Kpoint vectors
      if(size(kpt_orig,dim=1).ne.3) stop "size(kpt_orig,dim=1).ne.3"
      if(size(kpt_intp,dim=1).ne.3) stop "size(kpt_intp,dim=1).ne.3"
      Nkpt_orig = size(kpt_orig,dim=2)
      Nkpt_intp = size(kpt_intp,dim=2)
      !
      ! Size checks on Matrices
      Npoints = size(mat_orig,dim=3)
      if(size(mat_orig,dim=1).ne.size(mat_orig,dim=2)) stop "mat_orig not square."
      Nsize = size(mat_orig,dim=1)
      call assert_shape(mat_intp,[Nsize,Nsize,Npoints,Nkpt_intp],"wannierinterpolation_matrix","mat_intp")
      !




      if(size(mat_intp,dim=1).ne.size(mat_intp,dim=2)) stop "mat_intp not square."




      if(size(mat_orig,dim=1).ne.size(mat_intp,dim=1)) stop "mat_orig and mat_intp have different sizes."
      Nsize = size(mat_orig,dim=1)



















      subroutine wannierinterpolation_matrix_G(nkpt,kpt,nkpt_intp,kpt_intp,nwan,ndat,mat,mat_intp)
         implicit none
         integer,intent(in)           :: nkpt
         double precision,intent(in)  :: kpt(3,nkpt)
         integer,intent(in)           :: nkpt_intp
         double precision,intent(in)  :: kpt_intp(3,nkpt_intp)
         integer,intent(in)           :: nwan,ndat
         double complex,intent(in)    :: mat(nwan,nwan,1:ndat,nkpt)
         double complex,intent(out)   :: mat_intp(nwan,nwan,1:ndat,nkpt_intp)
         ! M(R)=\sum_{k} M(k) exp[-i*k*R]
         double complex,allocatable   :: matofR(:,:,:,:)
         integer                      :: ir,iwan1,iwan2,ik,id
         double precision,parameter   :: pi=3.14159265358979323846d0

         !
         !
         write(*,*) 'wannierinterpolation_matrix_G'
         if (.not.lwintp_init) stop 'init'
         if (nkpt.ne.size(nrdeg)/2) stop 'nkpt'
         !
         allocate(matofR(nwan,nwan,1:ndat,nwig))
         matofR(:,:,:,:)=0
         !
         ! M(R)=\sum_{k} M(k)*exp[-ik*R]
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(nkpt,nwig,ndat,nwan,kpt,rvec,mat,matofR),&
         !$OMP PRIVATE(ir,kR,cfac,id,iwan2,iwan1,ik)
         !$OMP DO
         do ir=1,nwig
            do ik=1,nkpt
               kR=2*pi*dot_product(kpt(:,ik),rvec(:,ir))
               cfac=dcmplx(cos(kR),-sin(kR))
               do id=1,ndat
                  do iwan2=1,nwan
                     do iwan1=1,nwan
                        !
                        !if (dabs(dreal(mat(iwan1,iwan2,id,ik))).gt.1.d-12.or.dabs(dimag(mat(iwan1,iwan2,id,ik))).gt.1.d-12) then
                        matofR(iwan1,iwan2,id,ir)=matofR(iwan1,iwan2,id,ir)+mat(iwan1,iwan2,id,ik)*cfac
                        !endif
                        !
                     enddo !iwan1
                  enddo ! iwan2
               enddo ! id
            enddo ! ik
         enddo ! ir
         !$OMP END DO
         !$OMP END PARALLEL
         matofR(:,:,:,:)=matofR(:,:,:,:)/nkpt
         !
         ! M(k_{intp})=\sum_{R} M(R)*exp[+ik_{intp}*R]
         mat_intp(:,:,:,:)=0
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(nkpt_intp,nwig,ndat,nwan,kpt_intp,rvec,mat_intp,matofR,nrdeg),&
         !$OMP PRIVATE(ir,kR,cfac,id,iwan2,iwan1,ik)
         !$OMP DO
         do ik=1,nkpt_intp
            do ir=1,nwig
               kR=2*pi*dot_product(kpt_intp(:,ik),rvec(:,ir))
               cfac=dcmplx(cos(kR),+sin(kR))/nrdeg(ir)
               do id=1,ndat
                  do iwan2=1,nwan
                     do iwan1=1,nwan
                        !
                        mat_intp(iwan1,iwan2,id,ik)=mat_intp(iwan1,iwan2,id,ik)+matofR(iwan1,iwan2,id,ir)*cfac
                        !
                     enddo !iwan1
                  enddo ! iwan2
               enddo ! id
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(matofR)
         !
      end subroutine wannierinterpolation_matrix_G
v


end module wannier_interpolation
