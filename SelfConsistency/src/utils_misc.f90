module utils_misc

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface str
      module procedure str_i_to_ch
      module procedure str_i_to_ch_pad
      module procedure str_r_to_ch
      module procedure str_c_to_ch
      module procedure str_l_to_ch
      module procedure str_ch_to_ch
   end interface str

   interface assert_shape
      module procedure i_assert_shape_N1
      module procedure i_assert_shape_N2
      module procedure i_assert_shape_N3
      module procedure i_assert_shape_N4
      module procedure i_assert_shape_N5
      module procedure i_assert_shape_N6
      module procedure i_assert_shape_N7
      module procedure d_assert_shape_N1
      module procedure d_assert_shape_N2
      module procedure d_assert_shape_N3
      module procedure d_assert_shape_N4
      module procedure d_assert_shape_N5
      module procedure d_assert_shape_N6
      module procedure d_assert_shape_N7
      module procedure z_assert_shape_N1
      module procedure z_assert_shape_N2
      module procedure z_assert_shape_N3
      module procedure z_assert_shape_N4
      module procedure z_assert_shape_N5
      module procedure z_assert_shape_N6
      module procedure z_assert_shape_N7
   end interface assert_shape

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),parameter,private                :: pi=3.14159265358979323846d0
   !
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: tick
   public :: tock
   public :: inquireFile
   public :: inquireDir
   public :: createDir
   public :: check_Hermiticity
   public :: check_Symmetry
   public :: assert_shape
   public :: FermionicFilon
   public :: BosonicFilon
   public :: halfbeta_symm
   !functions
   public :: FermionicFreqMesh
   public :: BosonicFreqMesh
   public :: fermidirac
   public :: diff_fermidirac
   public :: find_kpt
   public :: keq
   public :: linspace
   public :: denspace
   public :: free_unit
   public :: str
   public :: reg

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates the bosonic/fermionic Matsubara frequancy mesh
   !---------------------------------------------------------------------------!
   function FermionicFreqMesh(Beta,Nfreq) result(wmats)
      implicit none
      real(8),intent(in)                    :: Beta
      integer,intent(in)                    :: Nfreq
      real(8),dimension(Nfreq)              :: wmats
      integer                               :: iw
      !
      do iw=1,Nfreq
         wmats(iw)=(2.d0*dble(iw-1)+1.d0)*pi/Beta
      enddo
      !
   end function FermionicFreqMesh
   function BosonicFreqMesh(Beta,Nfreq) result(wmats)
      implicit none
      real(8),intent(in)                    :: Beta
      integer,intent(in)                    :: Nfreq
      real(8),dimension(Nfreq)              :: wmats
      integer                               :: iw
      !
      do iw=1,Nfreq
         wmats(iw)=2.d0*dble(iw-1)*pi/Beta
      enddo
      !
   end function BosonicFreqMesh


   !---------------------------------------------------------------------------!
   !PURPOSE: Calculates the fermi-dirac distribution f(e)
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   double precision function fermidirac(e,efermi,beta)
      implicit none
      real(8),intent(in)                    :: e,efermi,beta
      real(8)                               :: temp
      real(8)                               :: fermicut
      !
      fermicut=log(huge(1.0d0)-1e2)/2.d0 !fermicut=600.d0
      temp=(e-efermi)*beta
      if (temp.ge.fermicut) then
         fermidirac=0.d0
      elseif (temp.le.-fermicut) then
         fermidirac=1.d0
      elseif (temp.lt.0.d0.and.temp.gt.-fermicut) then
         fermidirac=1.d0/(dexp(temp)+1.d0)
      elseif (temp.ge.0.d0.and.temp.lt.fermicut) then
         fermidirac=dexp(-temp)/(dexp(-temp)+1.d0)
      endif
      !
   end function fermidirac


   !---------------------------------------------------------------------------!
   !PURPOSE: Energy derivative of fermi-dirac dist
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   double precision function diff_fermidirac(e,efermi,beta)
      implicit none
      real(8),intent(in)                    :: e,efermi,beta
      real(8)                               :: fermicut
      !
      fermicut=log(huge(1.0d0)-1e2)/2.d0 !fermicut=600.d0
      if (dabs((e-efermi)*beta).gt.fermicut) then
         diff_fermidirac=0.d0
      else
         diff_fermidirac=-dexp((e-efermi)*beta)*beta/((dexp((e-efermi)*beta)+1.d0)**2)
      endif
      !
   end function diff_fermidirac


   !---------------------------------------------------------------------------!
   !PURPOSE: Find index of K-point in list
   !---------------------------------------------------------------------------!
   function find_kpt(kvec,klist,tol) result(ikvec)
      implicit none
      real(8),intent(in)                    :: kvec(:)
      real(8),intent(in)                    :: klist(:,:)
      real(8),intent(in)                    :: tol
      integer                               :: ikvec
      integer                               :: ik,Nkpt
      Nkpt=size(klist,dim=2)
      do ik=1,Nkpt
         if((dabs(kvec(1)-klist(1,ik)).le.tol).and. &
            (dabs(kvec(2)-klist(2,ik)).le.tol).and. &
            (dabs(kvec(3)-klist(3,ik)).le.tol)) then
            ikvec=ik
            exit
         endif
      enddo
      if(ik.eq.Nkpt) then
          write(*,*)"requested k-point not found in k-point set"
          write(*,*)"kvec=",kvec
          stop "k-point for chiplot not found in k-point set"
      endif
   end function find_kpt


   !---------------------------------------------------------------------------!
   !PURPOSE: Return True if two K-points are equal
   !---------------------------------------------------------------------------!
   logical function keq(k1,k2)
      implicit none
      real(8),intent(in)                    :: k1(3),k2(3)
      real(8)                               :: dk(3),diff
      real(8),parameter                     :: eps=1e-6
      !
      dk(:)=k1(:)-k2(:)
      diff=sum(abs(dk(:)-nint(dk(:))))
      if (diff.lt.eps) then
        keq=.true.
      else
        keq=.false.
      endif
      !
   end function keq


   !---------------------------------------------------------------------------!
   !PURPOSE: analogous of python numpy.linspace
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   function linspace(start,stop,num,istart,iend,mesh) result(array)
      implicit none
      real(8),intent(in)                    :: start,stop
      integer,intent(in)                    :: num
      logical,intent(in),optional           :: istart,iend
      real(8)                               :: array(num)
      real(8),optional                      :: mesh
      !
      integer                               :: i
      real(8)                               :: step
      logical                               :: startpoint_,endpoint_
      !
      if(num<0)stop "linspace: N<0, abort."
      !
      startpoint_=.true.;if(present(istart))startpoint_=istart
      endpoint_=.true.;if(present(iend))endpoint_=iend
      !
      if(startpoint_.AND.endpoint_)then
         if(num<2)stop "linspace: N<2 with both start and end points"
         step = (stop-start)/(dble(num)-1d0)
         forall(i=1:num)array(i)=start + (dble(i)-1d0)*step
      elseif(startpoint_.AND.(.not.endpoint_))then
         step = (stop-start)/dble(num)
         forall(i=1:num)array(i)=start + (dble(i)-1d0)*step
      elseif(.not.startpoint_.AND.endpoint_)then
         step = (stop-start)/dble(num)
         forall(i=1:num)array(i)=start + dble(i)*step
      else
         step = (stop-start)/(dble(num)+1d0)
         forall(i=1:num)array(i)=start + dble(i)*step
      endif
       if(present(mesh))mesh=step
   end function linspace


   !---------------------------------------------------------------------------!
   !PURPOSE: generates imaginary time tau between 0 and beta
   ! the mesh is divided into segements where the end of the segments
   ! are distributed according to a shifted exponential mesh
   ! the mesh is symmetric about tau=beta/2, i.e., the mesh is densed
   ! around tau=0 and tau=beta
   ! Each segment is divided into two EQUAL mesh so in total there are
   ! 2*nseg+1 points.
   ! tau(2*i-1) = b * {exp[a*(i-1)] -1}, i=1, nseg+1
   ! tau(1) = 0
   ! tau(2*nseg+1) = beta  => b = beta/ {exp[a*nseg] -1}
   ! choose a = dtau*dE, where dE is about 1 eV.
   ! a test on Cu for P(iw) shows that de=1 eV gives the best result
   ! at least for T=500 and 1000K.
   ! beta and tau are in atomic unit.
   ! nseg = number of segments, must be even.
   ! nsimp fixed to 2
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   function denspace(end,num) result(array)
      implicit none
      real(8),intent(in)                    :: end
      integer,intent(in)                    :: num
      real(8)                               :: array(num)
      !
      integer                               :: nsimp=2
      integer                               :: nseg
      integer                               :: i,n,n2,nseg2
      real(8)                               :: mesh,de,a,b
      data de /1.0d0/
      !
      nseg=(num-1)/nsimp
      if (nseg .lt. 1) stop "denspace: nseg < 1"
      if (nsimp*(nseg/nsimp) .ne. nseg) stop "denspace: nseg is not a multiple of 2 or 4"
      nseg2 = nseg/2
      mesh = end/nseg
      a = mesh * de!/27.2d0
      b = (end/2.d0)/(dexp(a*nseg2)-1.d0)
      array(1) = 0.d0
      do n=1,nseg2
         n2 = nsimp * n
         array(n2+1) = b * (dexp(a*n)-1.d0)
         mesh = ( array(n2+1) - array(n2-nsimp+1) ) / dble(nsimp)
         do i=0,nsimp-2
            array(n2-i) = array(n2-i+1) - mesh
         enddo
      enddo
      do n=nseg2*nsimp+2,nsimp*nseg+1
         array(n) = end-array(nsimp*nseg+1-n+1)
      enddo
      if ( dabs(array(nsimp*nseg+1)-end).gt.1.d-9) stop "denspace: wrong endpoint"
      !
   end function denspace


   !---------------------------------------------------------------------------!
   !PURPOSE: Regularize string
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   function reg(string_in) result(string_out)
      implicit none
      character(len=*)                      :: string_in
      character(len=len_trim(trim(adjustl(trim(string_in))))) :: string_out
      string_out=trim(adjustl(trim(string_in)))
   end function reg


   !---------------------------------------------------------------------------!
   !PURPOSE: Looks for a free unit
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   function free_unit(n) result(unit_)
      implicit none
      integer,optional                      :: n
      integer                               :: unit_,ios
      logical                               :: opened
      unit_=100
      do
         unit_=unit_+1
         INQUIRE(unit=unit_,OPENED=opened,iostat=ios)
         if(.not.opened.AND.ios==0)exit
         if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
      enddo
      if(present(n))n=unit_
   end function free_unit


   !---------------------------------------------------------------------------!
   !PURPOSE: Returns time in seconds from now to time described by t
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine tick(t)
      integer, intent(out)                  :: t
      call system_clock(t)
   end subroutine tick
   real function tock(t)
      integer, intent(in)                   :: t
      integer                               :: now, clock_rate
      call system_clock(now,clock_rate)
      tock = real(now - t)/real(clock_rate)
   end function tock


   !---------------------------------------------------------------------------!
   !PURPOSE: Returns true if a file/directory exists
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine inquireFile(file,exists,hardstop,verb)
      implicit none
      character(len=*),intent(in)           :: file
      logical,intent(in),optional           :: hardstop
      logical,intent(in),optional           :: verb
      logical,intent(out)                   :: exists
      logical                               :: hardstop_,verbose_
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      verbose_=.true.
      if(present(verb))verbose_=verb
      !
      inquire(file=reg(file),exist=exists)
      if(.not.exists) then
         if(verbose_.or.hardstop_)write(*,"(A)")"Unable to find file: "//reg(file)
         if(hardstop_) stop "Stop."
      endif
      !
   end subroutine inquireFile
   subroutine inquireDir(dir,exists,hardstop,verb)
      implicit none
      character(len=*),intent(in)           :: dir
      logical,intent(in),optional           :: hardstop
      logical,intent(in),optional           :: verb
      logical,intent(out)                   :: exists
      logical                               :: hardstop_,verbose_
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      verbose_=.true.
      if(present(verb))verbose_=verb
      !
      inquire(directory=reg(dir),exist=exists)                                  !<===IFORT
      !inquire(file=reg(dir),exist=exists)                                      !<===GFORTRAN
      if(.not.exists) then
         if(verbose_.or.hardstop_)write(*,"(A)")"Unable to find directory: "//reg(dir)
         if(hardstop_) stop "Stop."
      endif
      !
   end subroutine inquireDir


   !---------------------------------------------------------------------------!
   !PURPOSE: Creat directory in path
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine createDir(dirpath,verb)
      implicit none
      character(len=*),intent(in)           :: dirpath
      logical,intent(in),optional           :: verb
      character(len=256)                    :: mkdirCmd
      logical                               :: direxists
      logical                               :: verbose_
      !
      verbose_=.true.
      if(present(verb))verbose_=verb
      !
      call inquireDir(reg(dirpath),direxists,hardstop=.false.,verb=verbose_)
      if(.not.direxists)then
         mkdirCmd = "mkdir -p "//reg(dirpath)
         if(verbose_)write(*,"(A)") "Creating new directory: "//reg(dirpath)
         if(verbose_)write(*,"(A)") reg(mkdirCmd)
         call system(reg(mkdirCmd))
         !call execute_command_line(reg(mkdirCmd))
      endif
      !
   end subroutine createDir


   !---------------------------------------------------------------------------!
   !PURPOSE: Check if matrix is Hermitian
   !TEST ON: 21-10-2020
   !---------------------------------------------------------------------------!
   subroutine check_Hermiticity(A,tol,hardstop)
      implicit none
      complex(8),intent(in)                 :: A(:,:)
      real(8),intent(in)                    :: tol
      logical,intent(in),optional           :: hardstop
      !
      complex(8)                            :: elem_ij,elem_ji
      logical                               :: hardstop_
      integer                               :: N,i,j
      !
      if(size(A,dim=1).ne.size(A,dim=2))stop "Hermitizer. Matrix not square."
      N=size(A,dim=1)
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      !
      do i=1,N
         do j=1+i,N
            !
            elem_ij = A(i,j)
            elem_ji = A(j,i)
            !
            if(abs(real(elem_ij)-real(elem_ji)).gt.tol)then
               write(*,"(2(A,1F10.5))") "Re(A_ij): ",real(elem_ij),"Re(A_ji): ",real(elem_ji)
               write(*,"(A,1F10.5)") "Real parts difference above threshold ",tol
               if(hardstop_)stop
            endif
            !
            if(abs(aimag(elem_ij)+aimag(elem_ji)).gt.tol)then
               write(*,"(2(A,1F10.5))") "Im(A_ij): ",aimag(elem_ij),"Im(A_ji): ",aimag(elem_ji)
               write(*,"(A,1F10.5)") "Imaginary parts difference above threshold ",tol
               if(hardstop_)stop
            endif
            !
         enddo
      enddo
      !
   end subroutine check_Hermiticity


   !---------------------------------------------------------------------------!
   !PURPOSE: Check if matrix is Symmetric
   !---------------------------------------------------------------------------!
   subroutine check_Symmetry(A,tol,hardstop,name)
      implicit none
      real(8),intent(in)                    :: A(:,:)
      real(8),intent(in)                    :: tol
      logical,intent(in),optional           :: hardstop
      character(len=*),intent(in),optional  :: name
      !
      real(8)                               :: elem_ij,elem_ji
      logical                               :: hardstop_
      integer                               :: N,i,j
      !
      if(present(name)) write(*,"(A)") "Matrix: "//reg(name)
      if(size(A,dim=1).ne.size(A,dim=2))stop "Symmetryzer. Matrix not square."
      N=size(A,dim=1)
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      !
      do i=1,N
         do j=1+i,N
            !
            elem_ij = A(i,j)
            elem_ji = A(j,i)
            !
            if(abs(elem_ij-elem_ji).gt.tol)then
               write(*,"(A,2I4)") "Components: ",i,j
               write(*,"(2(A,E20.12))") "Re(A_ij): ",real(elem_ij),"  Re(A_ji): ",real(elem_ji)
               write(*,"(A,E20.12)") "Real parts difference above threshold ",tol
               if(hardstop_)stop
            endif
            !
         enddo
      enddo
      !
   end subroutine check_Symmetry


   !---------------------------------------------------------------------------!
   !PURPOSE: Check if matrix is Hermitian
   !---------------------------------------------------------------------------!
   subroutine halfbeta_symm(funct,Ntau,sign)
      implicit none
      real(8),intent(inout)                 :: funct(:)
      integer,intent(in)                    :: Ntau
      real(8),intent(in)                    :: sign
      !
      integer                               :: itau
      !
      do itau=1,int(Ntau/2)
         funct(itau) = 0.5 *  funct(itau) + sign*funct(ntau-itau+1)
      enddo
      do itau=1,int(Ntau/2)
         funct(ntau-itau+1) = sign*funct(itau)
      enddo
      !
   end subroutine halfbeta_symm


   !---------------------------------------------------------------------------!
   !PURPOSE: General Filon integration
   ! I[x1,x2] dx f(x) cos(kx) and I[x1,x2] dx f(x) sin(kx)
   ! where f(x) is smooth but cos(kx) and sin(kx) can oscillate rapidly,
   ! i.e., k can be very large.
   !
   ! Divide the integration range into N segments which are NOT necessarily
   ! uniform. Each segment is divided further into two EQUAL segments of
   ! size h each. The input mesh is x(i), i=1, 2*nseg+1
   !
   ! The integral for the n-th segment centred at xn=x(2*n) is
   ! Icos(n) = I[xn-h, xn+h] dx f(x) cos(kx)
   ! = I[-h,+h] dy f(y+xn) cos(ky+ kxn)
   ! = I[-h,+h] dy g(y) [cos(ky) cos(kxn) - sin(ky) sin(kxn)]
   !
   ! Similarly
   ! Isin(n) = I[xn-h, xn+h] dx f(x) sin(kx)
   ! = I[-h,+h] dy f(y+xn) sin(ky+ kxn)
   ! = I[-h,+h] dy g(y) [sin(ky) cos(kxn) + cos(ky) sin(kxn)]
   !
   ! where y = x - x(n) and
   ! g(-h) = f(-h+xn), g(0) = f(xn), and g(h) = f(h+xn).
   !
   ! Fitting g(y) to an exponential + a square:
   ! g(y) = a  exp(b*y) + cy^2, we obtain
   ! a = g(0)
   ! B = [g(h)-g(-h)]/g(0)
   ! y+= [B + sqrt(B^2+4)] / 2 = exp(b*h) -> b*h = ln(y+)
   ! c = g(h) - a*exp(b*h)
   !   = g(-h) - a*exp(-b*h)
   !
   ! We need
   ! Ic0  = I[-h,h] dy exp(by) cos(ky)
   !      = [ (b*cos(kh) + k*sin(kh)) exp(bh)
   !         -(b*cos(kh) - k*sin(kh)) exp(-bh) ] / (b^2 + k^2)
   ! Ic2  = I[-h,h] dy y^2 cos(ky)
   !      = 2y cos(ky)/k^2 +   (y^2/k - 2/k^3) sin(ky) |y=-h,h
   !      = 4h cos(kh)/k^2 + 2 (h^2/k - 2/k^3) sin(kh)
   !
   ! Is1  = I[-h,h] dy exp(by) sin(ky)
   !      = [  (b*sin(kh) - k*cos(kh)) exp(bh)
   !         -(-b*sin(kh) - k*cos(kh)) exp(-bh) ] / (b^2 + k^2)
   !
   ! For small k:
   ! cos(ky) = 1 - (ky)^2/2 + (ky)^4/24 - (ky)^6/720 + ...
   ! sin(ky) =     (ky)     - (ky)^3/6  + (ky)^5/120 - ...
   !
   ! Ic2  = I[-h,h] dy y^2 cos(ky)
   !      = I[-h,h] dy y^2 [ 1 - (ky)^2/2 + (ky)^4/24 - (ky)^6/720 ]
   !      = 2h^3/3 - k^2 h^5/5 + k^4 h^7/84 - k^6 h^9/(9*360)
   !      = h^3 [ 2/3 - (kh)^2/5 + (kh)^4/84 - (kh)^6/(9*360) ]
   !
   ! Icos = I[-h,h] dy g(y) cos(ky)
   !      = a*Ic0 + c*Ic2
   ! Isin = I[-h,h] dy g(y) sin(ky)
   !      = a*Is1
   !
   ! Therefore
   ! Icos(n) =
   ! = I[-h,+h] dy g(y) [cos(ky) cos(kxn) - sin(ky) sin(kxn)]
   ! = cos(kxn) * Icos - sin(kxn) * Isin
   !
   ! Isin(n) =
   ! = I[-h,+h] dy g(y) [sin(ky) cos(kxn) + cos(ky) sin(kxn)]
   ! = cos(kxn) * Isin + sin(kxn) * Icos
   !
   ! Weight for Icos(n) and Isin(n)
   ! cos(kxn) * Icos - sin(kxn) * Isin
   ! = cos(kxn) * (a*Ic0 + c*Ic2) - sin(kxn) * a * Is1
   ! wcos(2*n)   = cos(kxn)*(Ic0 - Ic2/h^2)
   ! wcos(2*n-1) = cos(kxn)*Ic2/(2h^2) + sin(kxn)*Is1/(2h)
   ! wcos(2*n+1) = cos(kxn)*Ic2/(2h^2) - sin(kxn)*Is1/(2h)
   !
   ! cos(kxn) * Isin + sin(kxn) * Icos
   ! = cos(kxn) * a * Is1 + sin(kxn) * (a*Ic0 + c*Ic2)
   ! wsin(2*n)   = sin(kxn)*(Ic0 - Ic2/h^2)
   ! wsin(2*n-1) = sin(kxn)*Ic2/(2h^2) - cos(kxn)*Is1/(2h)
   ! wcos(2*n+1) = sin(kxn)*Ic2/(2h^2) + cos(kxn)*Is1/(2h)
   !
   ! 1) nseg is the number of segments.
   !    The number of mesh points MUST be odd (2*nseg+1).
   !    q = k, x = mesh with x(i) = [x(i-1) + x(i+1)]/2
   ! 2) Make sure that each segment is divided into two EQUAL segments.
   ! 3) The weights for cos and sin integration are in wcos and wsin.
   !TEST ON: 16-10-2020
   !---------------------------------------------------------------------------!
   subroutine FermionicFilon(q,x,fx,wcos,wsin)
      implicit none
      real(8),intent(in)                    :: q
      real(8),intent(in)                    :: x(:)
      real(8),intent(in)                    :: fx(:)
      real(8),intent(out)                   :: wcos,wsin
      !
      integer                               :: npoints,nseg
      integer                               :: n,n2
      real(8)                               :: oq,oq2,oq3,h,h2
      real(8)                               :: coskh,sinkh,coskx,sinkx
      real(8)                               :: c0,c2,s1,oh,oh2
      real(8)                               :: a,b,c,bb,yy
      real(8)                               :: expbh1,expbh2,scos,ssin
      !
      npoints = size(x)
      if(mod(npoints,2).eq.0) stop "FermionicFilon: npoints is even."
      nseg = (npoints-1)/2
      !
      wcos=0d0;wsin=0d0
      !
      oq=1.d0/q
      oq2=oq*oq
      oq3=oq*oq2
      !
      do n=1,nseg
         n2 = 2*n
         h = x(n2) - x(n2-1)
         !
         if(dabs(x(n2+1)-x(n2)-h) .gt. 1d-10) then
            write(*,*) "Segment= ",n
            stop "FermionicFilon: the above segment is not equally divided"
         endif
         !
         ! check that fx is not "zero"
         scos = ( fx(n2-1) + 4.d0*fx(n2) + fx(n2+1) )*h/3.d0
         if(dabs(scos).lt.1d-9) cycle !goto 1111
         h2    = h * h
         oh    = 1.d0/h
         oh2   = oh * oh
         coskx = dcos( q*x(n2) )
         sinkx = dsin( q*x(n2) )
         coskh = dcos( q*h )
         sinkh = dsin( q*h )
         !
         ! g(y) = a  exp(b*y) + cy^2, we obtain
         ! a = g(0)
         ! B = [g(h)-g(-h)]/g(0)
         ! y+= [B + sqrt(B^2+4)] / 2 = exp(b*h) -> b*h = ln(y+)
         ! c = [ g(h) - a*exp(b*h) ] / h^2 = [ g(-h) - a*exp(-b*h) ] / h^2
         a      = fx(n2)
         bb     = ( fx(n2+1) - fx(n2-1) ) / a
         yy     = 0.5d0 * ( bb + dsqrt(bb*bb+4.d0) )
         b      = dlog(yy) * oh
         expbh1 = yy
         expbh2 = 1.d0 / yy
         c      = ( fx(n2-1) - a * expbh2 ) * oh2
         !
         ! Ic0  = I[-h,h] dy exp(by) cos(ky)
         !      = [ (b*cos(kh) + k*sin(kh)) exp(bh) - (b*cos(kh) - k*sin(kh)) exp(-bh) ] / (b^2 + k^2)
         !
         ! Ic2  = I[-h,h] dy y^2 cos(ky)
         !      = 2y cos(ky)/k^2 +   (y^2/k - 2/k^3) sin(ky) |y=-h,h
         !      = 4h cos(kh)/k^2 + 2 (h^2/k - 2/k^3) sin(kh)
         !
         ! Is1  = I[-h,h] dy exp(by) sin(ky)
         !      = [  (b*sin(kh) - k*cos(kh)) exp(bh) - (-b*sin(kh) - k*cos(kh)) exp(-bh) ] / (b^2 + k^2)
         c0 = (b*coskh + q*sinkh) * expbh1 - (b*coskh - q*sinkh) * expbh2
         c0 = c0 / (b*b + q*q)
         c2 = 4.d0 * h * oq2 * coskh + 2.d0 * (oq*h2 - 2.d0*oq3) * sinkh
         s1 = (b*sinkh - q*coskh) * expbh1 + (b*sinkh + q*coskh) * expbh2
         s1 = s1 / (b*b + q*q)
         !
         ! Icos = I[-h,h] dy g(y) cos(ky) = a*Ic0 + c*Ic2
         ! Isin = I[-h,h] dy g(y) sin(ky) = a*Is1
         scos = a*c0 + c*c2
         ssin = a*s1
         !
         ! Icos(n) = I[-h,+h] dy g(y) [cos(ky) cos(kxn) - sin(ky) sin(kxn)] = cos(kxn) * Icos - sin(kxn) * Isin
         ! Isin(n) = I[-h,+h] dy g(y) [sin(ky) cos(kxn) + cos(ky) sin(kxn)] = cos(kxn) * Isin + sin(kxn) * Icos
         wcos = wcos + coskx*scos - sinkx*ssin
         wsin = wsin + coskx*ssin + sinkx*scos
         !1111   continue
      enddo
   end subroutine FermionicFilon


   !---------------------------------------------------------------------------!
   !PURPOSE: General Filon integration:
   ! I[x1,x2] dx f(x) cos(kx) and I[x1,x2] dx f(x) sin(kx)
   ! where f(x) is smooth but cos(kx) and sin(kx) can oscillate rapidly,
   ! i.e., k can be very large.
   !
   ! Divide the integration range into N segments which are NOT necessarily
   ! uniform. Each segment is divided further into two EQUAL segments of
   ! size h each. The input mesh is x(i), i=1, 2*nseg+1
   !
   ! The integral for the n-th segment centred at xn=x(2*n) is
   ! Icos(n) = I[xn-h, xn+h] dx f(x) cos(kx)
   ! = I[-h,+h] dy f(y+xn) cos(ky+ kxn)
   ! = I[-h,+h] dy g(y) [cos(ky) cos(kxn) - sin(ky) sin(kxn)]
   !
   ! Similarly
   ! Isin(n) = I[xn-h, xn+h] dx f(x) sin(kx)
   ! = I[-h,+h] dy f(y+xn) sin(ky+ kxn)
   ! = I[-h,+h] dy g(y) [sin(ky) cos(kxn) + cos(ky) sin(kxn)]
   !
   ! where y = x - x(n) and
   ! g(-h) = f(-h+xn), g(0) = f(xn), and g(h) = f(h+xn).
   !
   ! Fitting g(y) to a parabola g(y) = a + by + cy^2, we obtain
   ! a = g(0)
   ! b = [g(h)-g(-h)]/(2h)
   ! c = [g(-h)-2g(0)+g(h)] / (2h^2)
   !
   ! We need
   ! Ic0  = I[-h,h] dy cos(ky) = 2sin(kh)/k
   ! Ic2  = I[-h,h] dy y^2 cos(ky)
   !      = 2y cos(ky)/k^2 +   (y^2/k - 2/k^3) sin(ky) |y=-h,h
   !      = 4h cos(kh)/k^2 + 2 (h^2/k - 2/k^3) sin(kh)
   !
   ! Is1  = I[-h,h] dy y sin(ky)
   !      = sin(ky)/k^2 - y cos(ky)/k |y=-h,h
   !      = 2 sin(kh)/k^2 - 2h cos(kh)/k
   !
   ! For small k:
   ! cos(ky) = 1 - (ky)^2/2 + (ky)^4/24 - (ky)^6/720 + ...
   ! sin(ky) =     (ky)     - (ky)^3/6  + (ky)^5/120 - ...
   !
   ! Ic0  = I[-h,h] dy cos(ky) = 2sin(kh)/k
   !      = I[-h,h] dy [ 1 - (ky)^2/2 + (ky)^4/24 - (ky)^6/720 ]
   !      = 2h - k^2 h^3/3 + k^4 h^5/60 - k^6 h^7/(7*360)
   !      = h [ 2 - (kh)^2/3 + (kh)^4/60 - (kh)^6/(7*360) ]
   !
   ! Ic2  = I[-h,h] dy y^2 cos(ky)
   !      = I[-h,h] dy y^2 [ 1 - (ky)^2/2 + (ky)^4/24 - (ky)^6/720 ]
   !      = 2h^3/3 - k^2 h^5/5 + k^4 h^7/84 - k^6 h^9/(9*360)
   !      = h^3 [ 2/3 - (kh)^2/5 + (kh)^4/84 - (kh)^6/(9*360) ]
   !
   ! Is1  = I[-h,h] dy y sin(ky)
   !      = I[-h,h] dy y [ (ky) - (ky)^3/6  + (ky)^5/120 ]
   !      = 2k h^3/3 - k^3 h^5/15 + k^5 h^7/420
   !      = h^2 [ 2 (kh)/3 - (kh)^3/15 + (kh)^5/420 ]
   !
   ! Icos = I[-h,h] dy g(y) cos(ky)
   !      = a*Ic0 + c*Ic2
   ! Isin = I[-h,h] dy g(y) sin(ky)
   !      = b*Is1
   !
   ! Therefore
   ! Icos(n) =
   ! = I[-h,+h] dy g(y) [cos(ky) cos(kxn) - sin(ky) sin(kxn)]
   ! = cos(kxn) * Icos - sin(kxn) * Isin
   !
   ! Isin(n) =
   ! = I[-h,+h] dy g(y) [sin(ky) cos(kxn) + cos(ky) sin(kxn)]
   ! = cos(kxn) * Isin + sin(kxn) * Icos
   !
   ! Weight for Icos(n) and Isin(n)
   ! cos(kxn) * Icos - sin(kxn) * Isin
   ! = cos(kxn) * (a*Ic0 + c*Ic2) - sin(kxn) * b * Is1
   ! wcos(2*n)   = cos(kxn)*(Ic0 - Ic2/h^2)
   ! wcos(2*n-1) = cos(kxn)*Ic2/(2h^2) + sin(kxn)*Is1/(2h)
   ! wcos(2*n+1) = cos(kxn)*Ic2/(2h^2) - sin(kxn)*Is1/(2h)
   !
   ! cos(kxn) * Isin + sin(kxn) * Icos
   ! = cos(kxn) * b * Is1 + sin(kxn) * (a*Ic0 + c*Ic2)
   ! wsin(2*n)   = sin(kxn)*(Ic0 - Ic2/h^2)
   ! wsin(2*n-1) = sin(kxn)*Ic2/(2h^2) - cos(kxn)*Is1/(2h)
   ! wcos(2*n+1) = sin(kxn)*Ic2/(2h^2) + cos(kxn)*Is1/(2h)
   !
   !
   ! 1) nseg is the number of segments.
   !    The number of mesh points MUST be odd (2*nseg+1).
   !    q = k, x = mesh with x(i) = [x(i-1) + x(i+1)]/2
   ! 2) Make sure that each segment is divided into two EQUAL segments.
   ! 3) The weights for cos and sin integration are in wcos and wsin.
   !TEST ON: 21-10-2020
   !---------------------------------------------------------------------------!
   subroutine BosonicFilon(q,x,wcos,wsin)
      implicit none
      real(8),intent(in)                    :: q
      real(8),intent(in)                    :: x(:)
      real(8),intent(inout)                 :: wcos(:),wsin(:)
      !
      integer                               :: npoints,nseg
      integer                               :: n,n2
      real(8)                               :: oq,oq2,oq3,h,h2,h3
      real(8)                               :: coskh,sinkh,coskx,sinkx
      real(8)                               :: c0,c2,s1,oh,oh2,qh,qh2,qh3,qh4,qh5,qh6
      !
      npoints = size(x)
      if(mod(npoints,2).eq.0) stop "BosonicFilon: npoints is even."
      nseg = (npoints-1)/2
      !
      wcos=0d0;wsin=0d0
      !
      if(dabs(q).lt.1.d-2) then
         !
         !Small q
         do n=1,nseg
            n2 = 2 * n
            h = x(n2) - x(n2-1)
            !
            if(dabs(x(n2+1)-x(n2)-h).gt.1d-10) then
               write(*,*) "Segment= ",n
               stop "BosonicFilon: the above segment is not equally divided"
            endif
            !
            h2  = h * h
            h3  = h * h2
            oh  = 1.d0/h
            oh2 = oh * oh
            qh  = q * h
            qh2 = qh * qh
            qh3 = qh * qh2
            qh4 = qh * qh3
            qh5 = qh * qh4
            qh6 = qh * qh5
            !
            coskx = dcos( q*x(n2) )
            sinkx = dsin( q*x(n2) )
            coskh = dcos( qh )
            sinkh = dsin( qh )
            !
            ! For small k:
            ! Ic0  = h [ 2 - (kh)^2/3 + (kh)^4/60 - (kh)^6/(7*360) ]
            ! Ic2  = h^3 [ 2/3 - (kh)^2/5 + (kh)^4/84 - (kh)^6/(9*360) ]
            ! Is1  = h^2 [ 2 (kh)/3 - (kh)^3/15 + (kh)^5/420 ]
            c0 = h *  (2.d0 - qh2/3.d0 + qh4/60.d0 - qh6/2520.d0 )
            c2 = h3 * (2.d0/3.d0 - qh2/5.d0 + qh4/84.d0 - qh6/3240.d0)
            s1 = h2 * (2.d0*qh/3.d0 - qh3/15.d0 + qh5/420.d0)
            !
            ! Weight for Icos(n) and Isin(n)
            wcos(n2)   = wcos(n2)   + coskx * (c0 - oh2*c2)
            wcos(n2-1) = wcos(n2-1) + coskx * 0.5d0*oh2*c2+ sinkx * 0.5d0*oh *s1
            wcos(n2+1) = wcos(n2+1) + coskx * 0.5d0*oh2*c2- sinkx * 0.5d0*oh *s1
            !
            wsin(n2)   = wsin(n2)   + sinkx * (c0 - oh2*c2)
            wsin(n2-1) = wsin(n2-1) + sinkx * 0.5d0*oh2*c2 - coskx * 0.5d0*oh *s1
            wsin(n2+1) = wsin(n2+1) + sinkx * 0.5d0*oh2*c2 + coskx * 0.5d0*oh *s1
            !
         enddo
         !
      else
         !
         ! Not small q
         oq  = 1.d0/q
         oq2 = oq * oq
         oq3 = oq * oq2
         !
         do n=1,nseg
            n2 = 2 * n
            h = x(n2) - x(n2-1)
            !
            if(dabs(x(n2+1)-x(n2)-h).gt.1d-10) then
               write(*,*) "Segment= ",n
               stop "BosonicFilon: the above segment is not equally divided"
            endif
            !
            h2    = h * h
            oh    = 1.d0/h
            oh2   = oh * oh
            coskx = dcos( q*x(n2) )
            sinkx = dsin( q*x(n2) )
            coskh = dcos( q*h )
            sinkh = dsin( q*h )
            !
            ! Ic0  = 2sin(kh)/k
            ! Ic2  = 4h cos(kh)/k^2 + 2 (h^2/k - 2/k^3) sin(kh)
            ! Is1  = 2 sin(kh)/k^2 - 2h cos(kh)/k
            c0 = 2.d0 * oq * sinkh
            c2 = 4.d0 * h * oq2 * coskh + 2.d0 * (oq*h2 - 2.d0*oq3) * sinkh
            s1 = 2.d0 * oq2 * sinkh - 2.d0 * h * oq * coskh
            !
            ! Weight for Icos(n) and Isin(n)
            wcos(n2)   = wcos(n2)   + coskx * (c0 - oh2*c2)
            wcos(n2-1) = wcos(n2-1) + coskx * 0.5d0*oh2*c2+ sinkx * 0.5d0*oh *s1
            wcos(n2+1) = wcos(n2+1) + coskx * 0.5d0*oh2*c2- sinkx * 0.5d0*oh *s1
            !
            wsin(n2)   = wsin(n2)   + sinkx * (c0 - oh2*c2)
            wsin(n2-1) = wsin(n2-1) + sinkx * 0.5d0*oh2*c2- coskx * 0.5d0*oh *s1
            wsin(n2+1) = wsin(n2+1) + sinkx * 0.5d0*oh2*c2+ coskx * 0.5d0*oh *s1
            !
         enddo
         !
      endif
      !
    end subroutine BosonicFilon


   !---------------------------------------------------------------------------!
   !PURPOSE: Routines for the str interface
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   function str_i_to_ch(i4) result(string)
     integer                      :: i4
     character(len=:),allocatable :: string
     character(len=16)            :: string_
     call i4_to_s_left(i4,string_)
     string=trim(adjustl(trim(string_)))
   end function str_i_to_ch

   function str_i_to_ch_pad(i4,Npad) result(string)
     integer                      :: i4
     integer                      :: Npad
     character(len=:),allocatable :: string
     character(len=Npad)          :: string_pad
     call i4_to_s_zero(i4,string_pad)
     string=trim(adjustl(trim(string_pad)))
   end function str_i_to_ch_pad
   !
   function str_r_to_ch(r8,d) result(string)
     real(8)                      :: r8
     integer,optional             :: d
     integer                      :: w_,d_
     character(len=:),allocatable :: string
     character(len=:),allocatable :: string_
     d_=6 ;if(present(d))d_=d
     w_ = get_w_(r8,d_)
     allocate(character(len=w_) :: string_)
     call r8_to_s_left(r8,string_,d_)
     string=trim(adjustl(trim(string_)))
   end function str_r_to_ch
   !
   function str_c_to_ch(c,d) result(string)
     complex(8)                   :: c
     integer,optional             :: d
     integer                      :: w_,d_
     character(len=:),allocatable :: string
     character(len=:),allocatable :: sre,sim
     real(8)                      :: re,im
     d_=6 ;if(present(d))d_=d
     re=dreal(c)
     w_ = get_w_(re,d_)
     allocate(character(len=w_) :: sre)
     call r8_to_s_left(re,sre,d_)
     !
     im=dimag(c)
     w_ = get_w_(im,d_)
     allocate(character(len=w_) :: sim)
     call r8_to_s_left(im,sim,d_)
     string="("//trim(adjustl(trim(sre)))//","//trim(adjustl(trim(sim)))//")"
   end function str_c_to_ch
   !
   function str_l_to_ch(bool) result(string)
     logical          :: bool
     character(len=1) :: string
     string="F"
     if(bool)string="T"
   end function str_l_to_ch

   function str_ch_to_ch(txt) result(string)
     character(len=*)                             :: txt
     character(len=:),allocatable :: string
     string=trim(adjustl(trim(txt)))
   end function str_ch_to_ch
   !
   subroutine i4_to_s_left ( i4, s )
     !! I4_TO_S_LEFT converts an I4 to a left-justified string.
     !  Example:
     !    Assume that S is 6 characters long:
     !        I4  S
     !         1  1
     !        -1  -1
     !         0  0
     !      1952  1952
     !    123456  123456
     !   1234567  ******  <-- Not enough room!
     !  Parameters:
     !    Input, integer ( kind = 4 ) I4, an integer to be converted.
     !    Output, character ( len = * ) S, the representation of the integer.
     !    The integer will be left-justified.  If there is not enough space,
     !    the string will be filled with stars.
     character :: c
     integer   :: i
     integer   :: i4
     integer   :: idig
     integer   :: ihi
     integer   :: ilo
     integer   :: ipos
     integer   :: ival
     character(len=*) ::  s
     s = " "
     ilo = 1
     ihi = len ( s )
     if ( ihi <= 0 ) then
        return
     end if
     !  Make a copy of the integer.
     ival = i4
     !  Handle the negative sign.
     if ( ival < 0 ) then
        if ( ihi <= 1 ) then
           s(1:1) = "*"
           return
        end if
        ival = -ival
        s(1:1) = "-"
        ilo = 2
     end if
     !  The absolute value of the integer goes into S(ILO:IHI).
     ipos = ihi
     !  Find the last digit of IVAL, strip it off, and stick it into the string.
     do
        idig = mod ( ival, 10 )
        ival = ival / 10
        if ( ipos < ilo ) then
           do i = 1, ihi
              s(i:i) = "*"
           end do
           return
        end if
        call digit_to_ch ( idig, c )
        s(ipos:ipos) = c
        ipos = ipos - 1
        if ( ival == 0 ) then
           exit
        end if
     end do
     !  Shift the string to the left.
     s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
     s(ilo+ihi-ipos:ihi) = " "
   end subroutine i4_to_s_left
   !
   subroutine r8_to_s_left( r8, s, digits)
     !! R8_TO_S_LEFT writes an R8 into a left justified string.
     !    An R8 is a real ( kind = 8 ) value.
     !    A "F<len(s)>.DIGITS" format is used with a WRITE statement.
     character(len=12)   :: fmt
     !integer             :: i
     real(8)             :: r8
     character(len=*)    :: s
     integer             :: s_length,w_
     integer             :: digits
     s_length = len ( s )
     write(fmt,"(A2,I0,A1,I0,A1)")"(F",s_length,".",digits,")"
     if(r8/=0d0)then
        w_=floor(log10(abs(r8)))
        if(w_<-1)write(fmt,"(A3,I0,A1,I0,A1)")"(ES",s_length,".",digits,")"
     endif
     write ( s, fmt ) r8
     s = trim(adjustl(trim( s )))
   end subroutine r8_to_s_left
   !
   subroutine digit_to_ch(digit,ch)
     !! DIGIT_TO_CH returns the character representation of a decimal digit.
     !    Instead of CHAR, we now use the ACHAR function, which
     !    guarantees the ASCII collating sequence.
     !  Example:
     !    DIGIT   CH
     !    -----  ---
     !      0    "0"
     !      1    "1"
     !    ...    ...
     !      9    "9"
     !     17    "*"
     !  Parameters:
     !    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
     !    Output, character CH, the corresponding character.
     character :: ch
     integer   :: digit
     if ( 0 <= digit .and. digit <= 9 ) then
        ch = achar ( digit + 48 )
     else
        ch = "*"
     end if
   end subroutine digit_to_ch
   !
   subroutine i4_to_s_zero ( intval, s )
     !! I4_TO_S_ZERO converts an I4 to a string, with zero padding.
     !    An I4 is an integer ( kind = 4 ).
     !  Example:
     !    Assume that S is 6 characters long:
     !    INTVAL  S
     !         1  000001
     !        -1  -00001
     !         0  000000
     !      1952  001952
     !    123456  123456
     !   1234567  ******  <-- Not enough room!
     !  Parameters:
     !    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
     !    Output, character ( len = * ) S, the representation of the integer.
     !    The integer will be right justified, and zero padded.
     !    If there is not enough space, the string will be filled with stars.
     implicit none
     character c
     integer ( kind = 4 ) i
     integer ( kind = 4 ) idig
     integer ( kind = 4 ) ihi
     integer ( kind = 4 ) ilo
     integer ( kind = 4 ) intval
     integer ( kind = 4 ) ipos
     integer ( kind = 4 ) ival
     character ( len = * ) s
     s = " "
     ilo = 1
     ihi = len ( s )
     if ( ihi <= 0 ) then
        return
     end if
     !
     !  Make a copy of the integer.
     !
     ival = intval
     !
     !  Handle the negative sign.
     !
     if ( ival < 0 ) then
        if ( ihi <= 1 ) then
           s(1:1) = "*"
           return
        end if
        ival = -ival
        s(1:1) = "-"
        ilo = 2
     end if
     !
     !  Working from right to left, strip off the digits of the integer
     !  and place them into S(ILO:IHI).
     !
     ipos = ihi
     do while ( ival /= 0 .or. ipos == ihi )
        idig = mod ( ival, 10 )
        ival = ival / 10
        if ( ipos < ilo ) then
           do i = 1, ihi
              s(i:i) = "*"
           end do
           return
        end if
        call digit_to_ch ( idig, c )
        s(ipos:ipos) = c
        ipos = ipos - 1
     end do
     !
     !  Fill the empties with zeroes.
     !
     do i = ilo, ipos
        s(i:i) = "0"
     end do
     return
   end subroutine i4_to_s_zero
   !
   function get_w_(r8,d) result(w)
     real(8) :: r8
     integer :: d
     integer :: w
     if(r8==0d0)then
        w=d+4
     else
        w=floor(log10(abs(r8)))
        if(w < -1)then
           w = d + 4 + 4
        else
           w = w + d + 4
        endif
     endif
   end function get_w_


   !---------------------------------------------------------------------------!
   !PURPOSE: Routines for the assert_shape interface
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine i_assert_shape_N1(A,Ndim,routine,matname)
     integer,dimension(:),intent(in)          :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine i_assert_shape_N1
   subroutine i_assert_shape_N2(A,Ndim,routine,matname)
     integer,dimension(:,:),intent(in)          :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine i_assert_shape_N2
   subroutine i_assert_shape_N3(A,Ndim,routine,matname)
     integer,dimension(:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine i_assert_shape_N3
   subroutine i_assert_shape_N4(A,Ndim,routine,matname)
     integer,dimension(:,:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine i_assert_shape_N4
   subroutine i_assert_shape_N5(A,Ndim,routine,matname)
     integer,dimension(:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine i_assert_shape_N5
   subroutine i_assert_shape_N6(A,Ndim,routine,matname)
     integer,dimension(:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine i_assert_shape_N6
   subroutine i_assert_shape_N7(A,Ndim,routine,matname)
     integer,dimension(:,:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine i_assert_shape_N7
   subroutine d_assert_shape_N1(A,Ndim,routine,matname)
     real(8),dimension(:),intent(in)            :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine d_assert_shape_N1
   subroutine d_assert_shape_N2(A,Ndim,routine,matname)
     real(8),dimension(:,:),intent(in)          :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine d_assert_shape_N2
   subroutine d_assert_shape_N3(A,Ndim,routine,matname)
     real(8),dimension(:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine d_assert_shape_N3
   subroutine d_assert_shape_N4(A,Ndim,routine,matname)
     real(8),dimension(:,:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine d_assert_shape_N4
   subroutine d_assert_shape_N5(A,Ndim,routine,matname)
     real(8),dimension(:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine d_assert_shape_N5
   subroutine d_assert_shape_N6(A,Ndim,routine,matname)
     real(8),dimension(:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine d_assert_shape_N6
   subroutine d_assert_shape_N7(A,Ndim,routine,matname)
     real(8),dimension(:,:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine d_assert_shape_N7
   !
   !
   !
   subroutine z_assert_shape_N1(A,Ndim,routine,matname)
     complex(8),dimension(:),intent(in)         :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine z_assert_shape_N1
   subroutine z_assert_shape_N2(A,Ndim,routine,matname)
     complex(8),dimension(:,:),intent(in)          :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine z_assert_shape_N2
   subroutine z_assert_shape_N3(A,Ndim,routine,matname)
     complex(8),dimension(:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine z_assert_shape_N3
   subroutine z_assert_shape_N4(A,Ndim,routine,matname)
     complex(8),dimension(:,:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine z_assert_shape_N4
   subroutine z_assert_shape_N5(A,Ndim,routine,matname)
     complex(8),dimension(:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine z_assert_shape_N5
   subroutine z_assert_shape_N6(A,Ndim,routine,matname)
     complex(8),dimension(:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine z_assert_shape_N6
   subroutine z_assert_shape_N7(A,Ndim,routine,matname)
     complex(8),dimension(:,:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine z_assert_shape_N7

end module utils_misc
