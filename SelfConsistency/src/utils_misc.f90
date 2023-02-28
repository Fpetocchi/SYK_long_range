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
      module procedure l_assert_shape_N1
      module procedure l_assert_shape_N2
      module procedure l_assert_shape_N3
      module procedure l_assert_shape_N4
      module procedure l_assert_shape_N5
      module procedure l_assert_shape_N6
      module procedure l_assert_shape_N7
   end interface assert_shape

   interface sort_array
      module procedure sort_array_i
      module procedure sort_array_d
   end interface sort_array

   interface flip_array
      module procedure flip_array_i
      module procedure flip_array_d
      module procedure flip_array_z
   end interface flip_array

   interface find_vec
      module procedure find_vec_i
      module procedure find_vec_d
   end interface find_vec

   interface check_Symmetry
      module procedure check_Symmetry_d
      module procedure check_Symmetry_z
   end interface check_Symmetry

   interface halfbeta_sym
      module procedure halfbeta_sym_d
      module procedure halfbeta_sym_z
   end interface halfbeta_sym

   interface get_moments_F
      module procedure get_moments_F_d1
      module procedure get_moments_F_d2
      module procedure get_moments_F_d3
   end interface get_moments_F

   interface get_pattern
      module procedure get_pattern_i
      module procedure get_pattern_d
      module procedure get_pattern_z
   end interface get_pattern

   interface linspace
      module procedure linspace_i
      module procedure linspace_d
   end interface linspace

   interface fermidirac
      module procedure fermidirac_shifted
      module procedure fermidirac_centered
   end interface fermidirac

   interface boseeinstein
      module procedure boseeinstein_shifted
      module procedure boseeinstein_centered
   end interface boseeinstein

   interface F2Bindex
      module procedure F2Bindex_1
      module procedure F2Bindex_2
   end interface F2Bindex


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
   public :: F2Bindex
   public :: tick
   public :: tock
   public :: inquireFile
   public :: inquireDir
   public :: createDir
   public :: skip_header
   public :: check_Hermiticity
   public :: check_Symmetry
   public :: assert_shape
   public :: sort_array
   public :: FermionicFilon
   public :: BosonicFilon
   public :: halfbeta_symm
   public :: halfbeta_antisymm
   public :: halfbeta_sym
   public :: nspline
   public :: splint
   public :: get_pattern
   public :: KK_Im2Re
   public :: KK_Re2Im
   public :: get_moments_F
   public :: init_Uelements
   !functions
   public :: FermionicFreqMesh
   public :: BosonicFreqMesh
   public :: fermidirac
   public :: boseeinstein
   public :: diff_fermidirac
   public :: find_vec
   public :: keq
   public :: linspace
   public :: denspace
   public :: flip_array
   public :: free_unit
   public :: str
   public :: reg
   public :: cubic_interp
   public :: linear_interp_2y
   public :: linear_interp_2x
   public :: trapezoid_integration

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates the bosonic/fermionic Matsubara frequancy mesh
   !---------------------------------------------------------------------------!
   function FermionicFreqMesh(Beta,Nfreq,full) result(wmats)
      implicit none
      real(8),intent(in)                    :: Beta
      integer,intent(in)                    :: Nfreq
      logical,intent(in),optional           :: full
      real(8),dimension(Nfreq)              :: wmats
      integer                               :: iw,Npos
      logical                               :: full_
      !
      full_=.false.
      if(present(full))full_=full
      !
      wmats=0d0
      if(full_)then
         Npos = (Nfreq-1)/2
         do iw=1,Npos
            wmats(iw+Npos+1)=(2d0*dble(iw-1)+1d0)*pi/Beta
            wmats(Npos+1-iw)=-wmats(iw+Npos+1)
         enddo
      else
         do iw=1,Nfreq
            wmats(iw)=(2d0*dble(iw-1)+1d0)*pi/Beta
         enddo
      endif
      !
   end function FermionicFreqMesh
   !
   function BosonicFreqMesh(Beta,Nfreq,full) result(wmats)
      implicit none
      real(8),intent(in)                    :: Beta
      integer,intent(in)                    :: Nfreq
      logical,intent(in),optional           :: full
      real(8),dimension(Nfreq)              :: wmats
      integer                               :: iw,Npos
      logical                               :: full_
      !
      full_=.false.
      if(present(full))full_=full
      !
      wmats=0d0
      if(full_)then
         Npos = (Nfreq-1)/2
         do iw=1,Npos
            wmats(iw+Npos+1)=2d0*dble(iw)*pi/Beta
            wmats(Npos+1-iw)=-wmats(iw+Npos+1)
         enddo
      else
         do iw=1,Nfreq
            wmats(iw)=2d0*dble(iw-1)*pi/Beta
         enddo
      endif
      !
   end function BosonicFreqMesh


   !---------------------------------------------------------------------------!
   !PURPOSE: Universal map/ordering between Fermionic and Bosonic Wannier basis
   !---------------------------------------------------------------------------!
   subroutine F2Bindex_1(Norb,orbs,ib)
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
   end subroutine F2Bindex_1
   !
   subroutine F2Bindex_2(Norb,Lorbs,Rorbs,ib1,ib2)
      implicit none
      integer,intent(in)                    :: Lorbs(2)
      integer,intent(in)                    :: Rorbs(2)
      integer,intent(in)                    :: Norb
      integer,intent(out)                   :: ib1,ib2
      integer                               :: i,j,k,l
      !
      i = Lorbs(1)
      j = Lorbs(2)
      k = Rorbs(1)
      l = Rorbs(2)
      !
      ib1 = j + Norb*(i-1)
      ib2 = l + Norb*(k-1)
      !
   end subroutine F2Bindex_2


   !---------------------------------------------------------------------------!
   !PURPOSE: Calculates the fermi-dirac distribution f(e)
   !---------------------------------------------------------------------------!
   double precision function fermidirac_shifted(e,mu,beta)
      implicit none
      real(8),intent(in)                    :: e,mu,beta
      real(8)                               :: temp
      real(8)                               :: fermicut
      !
      fermicut=log(huge(1.0d0)-1e2)/2.d0 !fermicut=600.d0
      temp=(e-mu)*beta
      if(temp.ge.fermicut) then
         fermidirac_shifted=0.d0
      elseif(temp.le.-fermicut) then
         fermidirac_shifted=1.d0
      elseif(temp.lt.0.d0.and.temp.gt.-fermicut) then
         fermidirac_shifted=1.d0/(dexp(temp)+1.d0)
      elseif(temp.ge.0.d0.and.temp.lt.fermicut) then
         fermidirac_shifted=dexp(-temp)/(dexp(-temp)+1.d0)
      endif
      !
   end function fermidirac_shifted
   !
   double precision function fermidirac_centered(e,beta)
      implicit none
      real(8),intent(in)                    :: e,beta
      real(8)                               :: temp
      real(8)                               :: fermicut
      !
      fermicut=log(huge(1.0d0)-1e2)/2.d0 !fermicut=600.d0
      temp=e*beta
      if(temp.ge.fermicut) then
         fermidirac_centered=0.d0
      elseif(temp.le.-fermicut) then
         fermidirac_centered=1.d0
      elseif(temp.lt.0.d0.and.temp.gt.-fermicut) then
         fermidirac_centered=1.d0/(dexp(temp)+1.d0)
      elseif(temp.ge.0.d0.and.temp.lt.fermicut) then
         fermidirac_centered=dexp(-temp)/(dexp(-temp)+1.d0)
      endif
      !
   end function fermidirac_centered
   !
   double precision function boseeinstein_shifted(e,mu,beta)
      implicit none
      real(8),intent(in)                    :: e,mu,beta
      real(8)                               :: temp
      real(8)                               :: bosecut
      !
      bosecut=log(huge(1.0d0)-1e2)/2.d0 !bosecut=600.d0
      temp=(e-mu)*beta
      if(temp.ge.bosecut) then
         boseeinstein_shifted=0.d0
      elseif(temp.le.-bosecut) then
         boseeinstein_shifted=-1.d0
      elseif(abs(temp).lt.1d-6) then
         stop "boseeinstein_shifted: abs(beta*e).lt.1d-6 in Bose-Einstein distribution"
      else
         boseeinstein_shifted=1.d0/(dexp(temp)-1.d0)
      endif
      !
   end function boseeinstein_shifted
   !
   double precision function boseeinstein_centered(e,beta)
      implicit none
      real(8),intent(in)                    :: e,beta
      real(8)                               :: temp
      real(8)                               :: bosecut
      !
      bosecut=log(huge(1.0d0)-1e2)/2.d0 !bosecut=600.d0
      temp=e*beta
      if(temp.ge.bosecut) then
         boseeinstein_centered=0.d0
      elseif(temp.le.-bosecut) then
         boseeinstein_centered=-1.d0
      elseif(abs(temp).lt.1d-6) then
         stop "boseeinstein_centered: abs(beta*e).lt.1d-6 in Bose-Einstein distribution"
      else
         boseeinstein_centered=1.d0/(dexp(temp)-1.d0)
      endif
      !
   end function boseeinstein_centered


   !---------------------------------------------------------------------------!
   !PURPOSE: Energy derivative of fermi-dirac dist
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
   !PURPOSE: Find index of vector in list
   !---------------------------------------------------------------------------!
   function find_vec_i(vec,list,hardstop) result(ivec)
      implicit none
      integer,intent(in)                    :: vec(:)
      integer,intent(in)                    :: list(:,:)
      logical,intent(in),optional           :: hardstop
      logical                               :: l1,l2,l3,hardstop_
      integer                               :: i,ivec
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      if(size(list,dim=2).eq.0)stop "find_vec_i: 0 lenght list."
      !
      ivec=0
      do i=1,size(list,dim=2)
         l1 = vec(1) .eq. list(1,i)
         l2 = vec(2) .eq. list(2,i)
         l3 = vec(3) .eq. list(3,i)
         if(l1.and.l2.and.l3) then
            ivec = i
            exit
         endif
      enddo
      if((i.eq.size(list,dim=2)).or.(ivec.eq.0))then
         ivec=0
         if(hardstop_)then
            write(*,"(A,3I)")"find_vec_i: vec=",vec
            stop "find_vec_i: requested vector not found in list."
         endif
      endif
      !
   end function find_vec_i
   !
   function find_vec_d(vec,list,tol,hardstop) result(ivec)
      implicit none
      real(8),intent(in)                    :: vec(:)
      real(8),intent(in)                    :: list(:,:)
      real(8),intent(in)                    :: tol
      logical,intent(in),optional           :: hardstop
      logical                               :: l1,l2,l3,hardstop_
      integer                               :: i,ivec
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      if(size(list,dim=2).eq.0)stop "find_vec_i: 0 lenght list."
      !
      ivec=0
      do i=1,size(list,dim=2)
         l1 = dabs(vec(1)-list(1,i)) .le. tol
         l2 = dabs(vec(2)-list(2,i)) .le. tol
         l3 = dabs(vec(3)-list(3,i)) .le. tol
         if(l1.and.l2.and.l3) then
            ivec = i
            exit
         endif
      enddo
      if((i.eq.size(list,dim=2)).or.(ivec.eq.0))then
         ivec=0
         if(hardstop_)then
            write(*,"(A,3F20.12)")"find_vec_d: vec=",vec
            stop "find_vec_d: requested vector not found in list."
         endif
      endif
      !
   end function find_vec_d


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
   !---------------------------------------------------------------------------!
   function linspace_i(start,stop) result(array)
      implicit none
      integer,intent(in)                    :: start,stop
      integer                               :: array(stop-start+1)
      !
      integer                               :: i,num
      !
      num = stop-start+1
      if(num<0)stop "linspace_i: N<0, abort."
      forall(i=1:num)array(i)=start + (i-1)
      !
   end function linspace_i
   !
   function linspace_d(start,stop,num,istart,iend,mesh) result(array)
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
      if(num<0)stop "linspace_d: N<0, abort."
      !
      startpoint_=.true.;if(present(istart))startpoint_=istart
      endpoint_=.true.;if(present(iend))endpoint_=iend
      !
      if(startpoint_.AND.endpoint_)then
         if(num<2)stop "linspace_d: N<2 with both start and end points"
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
   end function linspace_d


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
   !---------------------------------------------------------------------------!
   function denspace(end,num,center) result(array)
      implicit none
      real(8),intent(in)                    :: end
      integer,intent(in)                    :: num
      logical,intent(in),optional           :: center
      real(8)                               :: array(num)
      !
      integer                               :: nsimp=2
      integer                               :: nseg
      integer                               :: i,n,n2,nseg2
      real(8)                               :: mesh,de,a,b
      logical                               :: center_
      data de /1.0d0/
      !
      nseg=(num-1)/nsimp
      if (nseg .lt. 1) stop "denspace: nseg < 1"
      if (nsimp*(nseg/nsimp) .ne. nseg) stop "denspace: nseg is not a multiple of 2 and 4"
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
      if( dabs(array(nsimp*nseg+1)-end).gt.1.d-9) stop "denspace: wrong endpoint"
      !
      center_=.false.
      if(present(center))center_=center
      if(center_)then
         n = (num+1)/2
         array(n:num) = array(1:n)
         do i=1,num/2
            array(n-i) = -array(n+i)
         enddo
      endif
      !
   end function denspace


   !---------------------------------------------------------------------------!
   !PURPOSE: Regularize string
   !---------------------------------------------------------------------------!
   function reg(string_in) result(string_out)
      implicit none
      character(len=*)                      :: string_in
      character(len=len_trim(trim(adjustl(trim(string_in))))) :: string_out
      string_out=trim(adjustl(trim(string_in)))
   end function reg


   !---------------------------------------------------------------------------!
   !PURPOSE: Looks for a free unit
   !---------------------------------------------------------------------------!
   function free_unit(n) result(unit_)
      implicit none
      integer,optional                      :: n
      integer                               :: unit_,ios
      logical                               :: opened
      unit_=100
      do
         unit_=unit_+1
         inquire(unit=unit_,OPENED=opened,iostat=ios)
         if(.not.opened.AND.ios==0)exit
         if(unit_>900) stop "free_unit: no unit free smaller than 900. Possible BUG"
      enddo
      if(present(n))n=unit_
   end function free_unit


   !---------------------------------------------------------------------------!
   !PURPOSE: Returns time in seconds from now to time described by t
   !---------------------------------------------------------------------------!
   subroutine tick(t)
      integer, intent(out)                  :: t
      call system_clock(t)
   end subroutine tick
   !
   real function tock(t)
      integer, intent(in)                   :: t
      integer                               :: now, clock_rate
      call system_clock(now,clock_rate)
      tock = real(dble(now)-dble(t))/real(clock_rate)
   end function tock


   !---------------------------------------------------------------------------!
   !PURPOSE: Returns true if a file/directory exists
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
         if(verbose_.or.hardstop_)write(*,"(A)")"     Unable to find file: "//reg(file)
         if(hardstop_) stop "inquireFile: Stop."
      endif
      !
   end subroutine inquireFile
   !
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
         if(verbose_.or.hardstop_)write(*,"(A)")"     Unable to find directory: "//reg(dir)
         if(hardstop_) stop "inquireDir: Stop."
      endif
      !
   end subroutine inquireDir


   !---------------------------------------------------------------------------!
   !PURPOSE: Creat directory in path
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
         if(verbose_)write(*,"(A)") "     Creating new directory: "//reg(dirpath)
         if(verbose_)write(*,"(A)") "     "//reg(mkdirCmd)
         call system(reg(mkdirCmd)) !call execute_command_line(reg(mkdirCmd))
      endif
      !
   end subroutine createDir


   !---------------------------------------------------------------------------!
   !PURPOSE: advance the read record of file
   !---------------------------------------------------------------------------!
   subroutine skip_header(unit,skips)
      implicit none
      integer,intent(in)                    :: unit,skips
      integer                               :: iskip
      logical                               :: itsopen=.false.
      inquire(unit=unit, opened=itsopen)
      if(itsopen)then
         do iskip=1,skips
            read(unit,*)
         enddo
      else
         stop "skip_header: requested unit does not correspond to file."
      endif
   end subroutine skip_header


   !---------------------------------------------------------------------------!
   !PURPOSE: Check if matrix is Hermitian
   !---------------------------------------------------------------------------!
   subroutine check_Hermiticity(A,tol,enforce,hardstop,name,verb,isHerm)
      implicit none
      complex(8),intent(inout)              :: A(:,:)
      real(8),intent(in)                    :: tol
      logical,intent(in),optional           :: enforce
      logical,intent(in),optional           :: hardstop
      character(len=*),intent(in),optional  :: name
      logical,intent(in),optional           :: verb
      logical,intent(out),optional          :: isHerm
      !
      real(8)                               :: ReErr,ImErr
      real(8)                               :: percReErr=0d0,percImErr=0d0
      logical                               :: hardstop_,enforce_,verb_,isHerm_
      integer                               :: N,i,j
      !
      if(size(A,dim=1).ne.size(A,dim=2))stop "check_Hermiticity: Matrix not square."
      N=size(A,dim=1)
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      enforce_=.false.
      if(present(enforce))enforce_=enforce
      verb_=.true.
      if(present(verb))verb_=verb
      !
      if(enforce_)hardstop_=.false.
      !
      do i=1,N
         do j=1+i,N
            !
            ReErr = abs(dreal(A(i,j))-dreal(A(j,i)))
            ImErr = abs(dimag(A(i,j))+dimag(A(j,i)))
            if(dreal(A(i,j)).ne.0d0) percReErr = abs(ReErr/dreal(A(i,j)))
            if(dimag(A(i,j)).ne.0d0) percImErr = abs(ImErr/dimag(A(i,j)))
            !
            if((ReErr.gt.tol).or.(ImErr.gt.tol))then
               !
               if((ReErr.gt.tol).and.verb_)then
                  if(present(name)) write(*,"(A)") "     Non-hermitian matrix: "//reg(name)
                  write(*,"(A,2I4,5(A,1E12.5))")   "     [i,j]:",i,j," Re(A_ij): ",dreal(A(i,j)) ," Re(A_ji): ",dreal(A(j,i)) ," err: ",ReErr," %err: ",percReErr," > ",tol
               endif
               if((ImErr.gt.tol).and.verb_)then
                  if(present(name)) write(*,"(A)") "     Non-hermitian matrix: "//reg(name)
                  write(*,"(A,2I4,5(A,1E12.5))")   "     [i,j]:",i,j," Im(A_ij): ",dimag(A(i,j))," Im(A_ji): ",dimag(A(j,i))," err: ",ImErr," %err: ",percImErr," > ",tol
               endif
               !
               if(hardstop_)stop  "check_Hermiticity: Matrix not Hermitian."
               isHerm_=.false.
               !
            endif
            !
            if(enforce_)then
               A(i,j) = (A(i,j)+conjg(A(j,i)))/2d0
               A(j,i) = conjg(A(i,j))
            endif
            !
         enddo
      enddo
      if(present(isHerm))isHerm=isHerm_
      !
   end subroutine check_Hermiticity


   !---------------------------------------------------------------------------!
   !PURPOSE: Check if matrix is Symmetric
   !---------------------------------------------------------------------------!
   subroutine check_Symmetry_d(A,tol,enforce,hardstop,name,verb,isSym)
      implicit none
      real(8),intent(inout)                 :: A(:,:)
      real(8),intent(in)                    :: tol
      logical,intent(in),optional           :: enforce
      logical,intent(in),optional           :: hardstop
      character(len=*),intent(in),optional  :: name
      logical,intent(in),optional           :: verb
      logical,intent(out),optional          :: isSym
      !
      real(8)                               :: ReErr,percReErr=0d0
      logical                               :: hardstop_,enforce_,verb_,isSym_
      integer                               :: N,i,j
      !
      if(size(A,dim=1).ne.size(A,dim=2))stop "check_Symmetry_d: Matrix not square."
      N=size(A,dim=1)
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      enforce_=.false.
      if(present(enforce))enforce_=enforce
      verb_=.true.
      if(present(verb))verb_=verb
      !
      if(enforce_)hardstop_=.false.
      !
      do i=1,N
         do j=1+i,N
            !
            ReErr = A(i,j)-A(j,i)
            if(A(i,j).ne.0d0) percReErr = abs(ReErr/A(i,j))
            !
            if(ReErr.gt.tol)then
               !
               if(verb_)then
                  if(present(name)) write(*,"(A)") "     Non-symmetric matrix: "//reg(name)
                  write(*,"(A,2I4,5(A,1E12.5))")   "     [i,j]:",i,j," Re(A_ij): ",A(i,j)," Re(A_ji): ",A(j,i)," err: ",ReErr," %err: ",percReErr," > ",tol
               endif
               !
               if(hardstop_)stop "check_Symmetry_d: Matrix not symmatric."
               isSym_=.false.
               !
            endif
            !
            if(enforce_)then
               A(i,j) = (A(i,j)+A(j,i))/2d0
               A(j,i) = A(i,j)
            endif
            !
         enddo
      enddo
      if(present(isSym))isSym=isSym_
      !
   end subroutine check_Symmetry_d
   !
   subroutine check_Symmetry_z(A,tol,enforce,hardstop,name,verb,isSym)
      implicit none
      complex(8),intent(inout)              :: A(:,:)
      real(8),intent(in)                    :: tol
      logical,intent(in),optional           :: enforce
      logical,intent(in),optional           :: hardstop
      character(len=*),intent(in),optional  :: name
      logical,intent(in),optional           :: verb
      logical,intent(out),optional          :: isSym
      !
      real(8)                               :: ReErr,ImErr
      real(8)                               :: percReErr=0d0,percImErr=0d0
      logical                               :: hardstop_,enforce_,verb_,isSym_
      integer                               :: N,i,j
      !
      if(size(A,dim=1).ne.size(A,dim=2))stop "check_Symmetry_z: Matrix not square."
      N=size(A,dim=1)
      !
      hardstop_=.true.
      if(present(hardstop))hardstop_=hardstop
      enforce_=.false.
      if(present(enforce))enforce_=enforce
      verb_=.true.
      if(present(verb))verb_=verb
      !
      if(enforce_)hardstop_=.false.
      !
      do i=1,N
         do j=1+i,N
            !
            ReErr = abs(dreal(A(i,j))-dreal(A(j,i)))
            ImErr = abs(dimag(A(i,j))-dimag(A(j,i)))
            if(dreal(A(i,j)).ne.0d0) percReErr = abs(ReErr/dreal(A(i,j)))
            if(dimag(A(i,j)).ne.0d0) percImErr = abs(ImErr/dimag(A(i,j)))
            !
            if((ReErr.gt.tol).or.(ImErr.gt.tol))then
               !
               if((ReErr.gt.tol).and.verb_)then
                  if(present(name)) write(*,"(A)") "     Non-symmetric matrix: "//reg(name)
                  write(*,"(A,2I4,5(A,1E12.5))")   "     [i,j]:",i,j," Re(A_ij): ",dreal(A(i,j)) ," Re(A_ji): ",dreal(A(j,i)) ," err: ",ReErr," %err: ",percReErr," > ",tol
               endif
               if((ImErr.gt.tol).and.verb_)then
                  if(present(name)) write(*,"(A)") "     Non-symmetric matrix: "//reg(name)
                  write(*,"(A,2I4,5(A,1E12.5))")   "     [i,j]:",i,j," Im(A_ij): ",dimag(A(i,j))," Im(A_ji): ",dimag(A(j,i))," err: ",ImErr," %err: ",percImErr," > ",tol
               endif
               !
               if(hardstop_)stop "check_Symmetry_z: Matrix not symmatric."
               isSym_=.false.
               !
            endif
            !
            if(enforce_)then
               A(i,j) = (A(i,j)+A(j,i))/2d0
               A(j,i) = A(i,j)
            endif
            !
         enddo
      enddo
      if(present(isSym))isSym=isSym_
      !
   end subroutine check_Symmetry_z


   !---------------------------------------------------------------------------!
   !PURPOSE: Indentify the Tensor indexes which correspond to physical (number
   !         and spin conserving) interaction elements
   !---------------------------------------------------------------------------!
   subroutine init_Uelements(Norb,Uelements)
      !
      use parameters
      implicit none
      !
      integer,intent(in)                    :: Norb
      type(physicalU),intent(inout)         :: Uelements
      !
      integer                               :: Nflavor
      integer                               :: ib1,ib2
      integer                               :: iorb,jorb,korb,lorb
      integer                               :: ispin,jspin
      !
      !
      if(verbose)write(*,"(A)") "---- init_Uelements"
      !
      !
      if(Uelements%status) write(*,"(A)") "     Warning: the Physical interaction elements container is being reinitialized."
      Nflavor = Norb*Nspin
      !
      ! Elements when the interaction is in the Norb*Nspin form
      Uelements%Flav_Size = Nflavor
      if(allocated(Uelements%Flav_Uloc))deallocate(Uelements%Flav_Uloc)
      if(allocated(Uelements%Flav_U1st))deallocate(Uelements%Flav_U1st)
      if(allocated(Uelements%Flav_U2nd))deallocate(Uelements%Flav_U2nd)
      if(allocated(Uelements%Flav_All)) deallocate(Uelements%Flav_All)
      if(allocated(Uelements%Flav_Map)) deallocate(Uelements%Flav_Map)
      allocate(Uelements%Flav_Uloc(Nflavor,Nflavor)) ;Uelements%Flav_Uloc=.false.
      allocate(Uelements%Flav_U1st(Nflavor,Nflavor)) ;Uelements%Flav_U1st=.false.
      allocate(Uelements%Flav_U2nd(Nflavor,Nflavor)) ;Uelements%Flav_U2nd=.false.
      allocate(Uelements%Flav_All(Nflavor,Nflavor))  ;Uelements%Flav_All =.false.
      allocate(Uelements%Flav_Map(Nflavor,Nflavor,4));Uelements%Flav_Map=0
      !
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            iorb = (ib1+mod(ib1,2))/2
            jorb = (ib2+mod(ib2,2))/2
            ispin = abs(mod(ib1,2)-2)
            jspin = abs(mod(ib2,2)-2)
            !
            Uelements%Flav_Uloc(ib1,ib2) = (iorb.eq.jorb).and.(ispin.ne.jspin)
            Uelements%Flav_U1st(ib1,ib2) = (iorb.ne.jorb).and.(ispin.ne.jspin)
            Uelements%Flav_U2nd(ib1,ib2) = (iorb.ne.jorb).and.(ispin.eq.jspin)
            !
            Uelements%Flav_All(ib1,ib2) = Uelements%Flav_Uloc(ib1,ib2) .or.  &
                                          Uelements%Flav_U1st(ib1,ib2) .or.  &
                                          Uelements%Flav_U2nd(ib1,ib2)
            !
            Uelements%Flav_Map(ib1,ib2,1) = iorb
            Uelements%Flav_Map(ib1,ib2,2) = jorb
            Uelements%Flav_Map(ib1,ib2,3) = ispin
            Uelements%Flav_Map(ib1,ib2,4) = jspin
            !
         enddo
      enddo
      !
      ! Elements when the interaction is in the Norb^2 form
      Uelements%Full_Size = Norb*Norb
      if(allocated(Uelements%Full_Uaa))deallocate(Uelements%Full_Uaa)
      if(allocated(Uelements%Full_Uab))deallocate(Uelements%Full_Uab)
      if(allocated(Uelements%Full_Jsf))deallocate(Uelements%Full_Jsf)
      if(allocated(Uelements%Full_Jph))deallocate(Uelements%Full_Jph)
      if(allocated(Uelements%Full_Imp))deallocate(Uelements%Full_Imp)
      if(allocated(Uelements%Full_All))deallocate(Uelements%Full_All)
      if(allocated(Uelements%Full_Map))deallocate(Uelements%Full_Map)
      allocate(Uelements%Full_Uaa(Norb*Norb,Norb*Norb))  ;Uelements%Full_Uaa=.false.
      allocate(Uelements%Full_Uab(Norb*Norb,Norb*Norb))  ;Uelements%Full_Uab=.false.
      allocate(Uelements%Full_Jsf(Norb*Norb,Norb*Norb))  ;Uelements%Full_Jsf=.false.
      allocate(Uelements%Full_Jph(Norb*Norb,Norb*Norb))  ;Uelements%Full_Jph=.false.
      allocate(Uelements%Full_Imp(Norb*Norb,Norb*Norb))  ;Uelements%Full_Imp=.false.
      allocate(Uelements%Full_All(Norb*Norb,Norb*Norb))  ;Uelements%Full_All=.false.
      allocate(Uelements%Full_Map(Norb*Norb,Norb*Norb,4));Uelements%Full_Map=0
      !
      do iorb=1,Norb
         do jorb=1,Norb
            do korb=1,Norb
               do lorb=1,Norb
                  !
                  call F2Bindex(Norb,[iorb,jorb],[korb,lorb],ib1,ib2)
                  !
                  Uelements%Full_Uaa(ib1,ib2) = (iorb.eq.jorb).and.(korb.eq.lorb).and.(iorb.eq.korb)
                  Uelements%Full_Uab(ib1,ib2) = (iorb.eq.jorb).and.(korb.eq.lorb).and.(iorb.ne.korb)
                  Uelements%Full_Jsf(ib1,ib2) = (iorb.eq.lorb).and.(jorb.eq.korb).and.(iorb.ne.jorb)
                  Uelements%Full_Jph(ib1,ib2) = (iorb.eq.korb).and.(jorb.eq.lorb).and.(iorb.ne.jorb)
                  !
                  Uelements%Full_Imp(ib1,ib2) = Uelements%Full_Uaa(ib1,ib2) .or. Uelements%Full_Uab(ib1,ib2)
                  !
                  Uelements%Full_All(ib1,ib2) = Uelements%Full_Uaa(ib1,ib2) .or.  &
                                                Uelements%Full_Uab(ib1,ib2) .or.  &
                                                Uelements%Full_Jsf(ib1,ib2) .or.  &
                                                Uelements%Full_Jph(ib1,ib2)
                  !
                  Uelements%Full_Map(ib1,ib2,1) = iorb
                  Uelements%Full_Map(ib1,ib2,2) = jorb
                  Uelements%Full_Map(ib1,ib2,3) = korb
                  Uelements%Full_Map(ib1,ib2,4) = lorb
                  !
               enddo
            enddo
         enddo
      enddo
      !
      Uelements%status=.true.
      !
   end subroutine init_Uelements


   !---------------------------------------------------------------------------!
   !PURPOSE: Deduce simmetrical indexes in a vector.
   !---------------------------------------------------------------------------!
   subroutine get_pattern_i(list,pattern,listDim,IncludeSingle)
      implicit none
      integer,allocatable,intent(out)       :: list(:,:)
      integer,intent(in)                    :: pattern(:)
      integer,allocatable,intent(out),optional :: listDim(:)
      logical,intent(in),optional           :: IncludeSingle
      integer                               :: Nelements,iel,jel,iset,jset
      integer                               :: maxLen,Neqv
      integer,allocatable                   :: actRow(:),listDim_(:)
      integer,allocatable                   :: listOld(:,:),listDimOld(:)
      logical                               :: addSet,IncludeSingle_
      !
      IncludeSingle_=.false.
      if(present(IncludeSingle))IncludeSingle_=IncludeSingle
      !
      Nelements = size(pattern)
      !
      !this is just to get the maximum row length
      maxLen = 1
      allocate(actRow(Nelements));actRow=0
      do iel=1,Nelements
         !
         actRow=0
         do jel=1,Nelements
            if(pattern(iel).eq.pattern(jel)) actRow(jel) = jel
         enddo
         Neqv = size( pack( actRow, actRow.gt.0 ) )
         if(Neqv.gt.maxLen) maxLen = Neqv
         !
      enddo
      !
      if(allocated(list)) deallocate(list)
      !
      iset=0
      do iel=1,Nelements
         !
         actRow=0
         do jel=1,Nelements
            if(pattern(iel).eq.pattern(jel)) actRow(jel) = jel
         enddo
         Neqv = size( pack( actRow, actRow.gt.0 ) )
         !
         if(Neqv.gt.1)then
            !
            if(iset.eq.0)then
               !
               iset = iset + 1
               !
               allocate(list(iset,maxLen));list=0
               allocate(listDim_(iset));listDim_=0
               !
               list(iset,1:Neqv) = pack( actRow, actRow.gt.0 )
               listDim_(iset) = Neqv
               !
            else
               !
               addSet=.true.
               do jset=1,iset
                  addSet = addSet .and. all( pack( actRow, actRow.gt.0 ) .ne. list(jset,1:listDim_(jset)) )
               enddo
               !
               if(addSet)then
                  !
                  listOld = list
                  listDimOld = listDim_
                  deallocate(list,listDim_)
                  allocate(list(iset+1,maxLen));list(1:iset,:)=listOld
                  allocate(listDim_(iset+1));listDim_(1:iset)=listDimOld
                  deallocate(listOld,listDimOld)
                  !
                  iset = iset + 1
                  !
                  list(iset,1:Neqv) = pack( actRow, actRow.gt.0 )
                  listDim_(iset) = Neqv
                  !
               endif
               !
            endif
            !
         endif
         !
      enddo
      deallocate(actRow)
      !
      if(present(listDim).and.allocated(list))then
         if(allocated(listDim))deallocate(listDim)
         listDim = listDim_
         deallocate(listDim_)
      endif
      !
   end subroutine get_pattern_i
   !
   subroutine get_pattern_d(list,pattern,tol,listDim,IncludeSingle)
      implicit none
      integer,allocatable,intent(out)       :: list(:,:)
      real(8),intent(in)                    :: pattern(:)
      real(8),intent(in)                    :: tol
      integer,allocatable,intent(out),optional :: listDim(:)
      logical,intent(in),optional           :: IncludeSingle
      integer                               :: Nelements,iel,jel,iset,jset
      integer                               :: maxLen,Neqv
      integer,allocatable                   :: actRow(:),listDim_(:)
      integer,allocatable                   :: listOld(:,:),listDimOld(:)
      logical                               :: addSet,IncludeSingle_
      !
      IncludeSingle_=.false.
      if(present(IncludeSingle))IncludeSingle_=IncludeSingle
      !
      Nelements = size(pattern)
      !
      !this is just to get the maximum row length
      maxLen = 1
      allocate(actRow(Nelements));actRow=0
      do iel=1,Nelements
         !
         actRow=0
         do jel=1,Nelements
            if(abs(pattern(iel)-pattern(jel)).lt.tol) actRow(jel) = jel
         enddo
         Neqv = size( pack( actRow, actRow.gt.0 ) )
         if(Neqv.gt.maxLen) maxLen = Neqv
         !
      enddo
      !
      if(allocated(list)) deallocate(list)
      !
      iset=0
      do iel=1,Nelements
         !
         actRow=0
         do jel=1,Nelements
            if(abs(pattern(iel)-pattern(jel)).lt.tol) actRow(jel) = jel
         enddo
         Neqv = size( pack( actRow, actRow.gt.0 ) )
         !
         if((Neqv.gt.1).or.IncludeSingle_)then
            !
            if(iset.eq.0)then
               !
               iset = iset + 1
               !
               allocate(list(iset,maxLen));list=0
               allocate(listDim_(iset));listDim_=0
               !
               list(iset,1:Neqv) = pack( actRow, actRow.gt.0 )
               listDim_(iset) = Neqv
               !
            else
               !
               addSet=.true.
               do jset=1,iset
                  addSet = addSet .and. all( pack( actRow, actRow.gt.0 ) .ne. list(jset,1:listDim_(jset)) )
               enddo
               !
               if(addSet)then
                  !
                  listOld = list
                  listDimOld = listDim_
                  deallocate(list,listDim_)
                  allocate(list(iset+1,maxLen));list(1:iset,:)=listOld
                  allocate(listDim_(iset+1));listDim_(1:iset)=listDimOld
                  deallocate(listOld,listDimOld)
                  !
                  iset = iset + 1
                  !
                  list(iset,1:Neqv) = pack( actRow, actRow.gt.0 )
                  listDim_(iset) = Neqv
                  !
               endif
               !
            endif
            !
         endif
         !
      enddo
      deallocate(actRow)
      !
      if(present(listDim).and.allocated(list))then
         if(allocated(listDim))deallocate(listDim)
         listDim = listDim_
         deallocate(listDim_)
      endif
      !
   end subroutine get_pattern_d
   !
   subroutine get_pattern_z(list,pattern,tol,listDim,IncludeSingle)
      implicit none
      integer,allocatable,intent(out)       :: list(:,:)
      complex(8),intent(in)                 :: pattern(:)
      real(8),intent(in)                    :: tol
      integer,allocatable,intent(out),optional :: listDim(:)
      logical,intent(in),optional           :: IncludeSingle
      integer                               :: Nelements,iel,jel,iset,jset
      integer                               :: maxLen,Neqv
      integer,allocatable                   :: actRow(:),listDim_(:)
      integer,allocatable                   :: listOld(:,:),listDimOld(:)
      logical                               :: addSet,IncludeSingle_
      !
      IncludeSingle_=.false.
      if(present(IncludeSingle))IncludeSingle_=IncludeSingle
      !
      Nelements = size(pattern)
      !
      !this is just to get the maximum row length
      maxLen = 1
      allocate(actRow(Nelements));actRow=0
      do iel=1,Nelements
         !
         actRow=0
         do jel=1,Nelements
            if(abs(pattern(iel)-pattern(jel)).lt.tol) actRow(jel) = jel
         enddo
         Neqv = size( pack( actRow, actRow.gt.0 ) )
         if(Neqv.gt.maxLen) maxLen = Neqv
         !
      enddo
      !
      if(allocated(list)) deallocate(list)
      !
      iset=0
      do iel=1,Nelements
         !
         actRow=0
         do jel=1,Nelements
            if(abs(pattern(iel)-pattern(jel)).lt.tol) actRow(jel) = jel
         enddo
         Neqv = size( pack( actRow, actRow.gt.0 ) )
         !
         if(Neqv.gt.1)then
            !
            if(iset.eq.0)then
               !
               iset = iset + 1
               !
               allocate(list(iset,maxLen));list=0
               allocate(listDim_(iset));listDim_=0
               !
               list(iset,1:Neqv) = pack( actRow, actRow.gt.0 )
               listDim_(iset) = Neqv
               !
            else
               !
               addSet=.true.
               do jset=1,iset
                  addSet = addSet .and. all( pack( actRow, actRow.gt.0 ) .ne. list(jset,1:listDim_(jset)) )
               enddo
               !
               if(addSet)then
                  !
                  listOld = list
                  listDimOld = listDim_
                  deallocate(list,listDim_)
                  allocate(list(iset+1,maxLen));list(1:iset,:)=listOld
                  allocate(listDim_(iset+1));listDim_(1:iset)=listDimOld
                  deallocate(listOld,listDimOld)
                  !
                  iset = iset + 1
                  !
                  list(iset,1:Neqv) = pack( actRow, actRow.gt.0 )
                  listDim_(iset) = Neqv
                  !
               endif
               !
            endif
            !
         endif
         !
      enddo
      deallocate(actRow)
      !
      if(present(listDim).and.allocated(list))then
         if(allocated(listDim))deallocate(listDim)
         listDim = listDim_
         deallocate(listDim_)
      endif
      !
   end subroutine get_pattern_z


   !---------------------------------------------------------------------------!
   !PURPOSE: Enforce symmetry/antisymmetry with respect to beta/2
   !---------------------------------------------------------------------------!
   subroutine halfbeta_symm(funct)
      implicit none
      real(8),intent(inout)                 :: funct(:)
      integer                               :: itau,Ntau
      real(8)                               :: sign=+1d0
      Ntau=size(funct)
      do itau=1,int(Ntau/2)
         funct(itau) = 0.5 * (funct(itau) + sign*funct(ntau-itau+1))
      enddo
      do itau=1,int(Ntau/2)
         funct(ntau-itau+1) = sign*funct(itau)
      enddo
   end subroutine halfbeta_symm
   !
   subroutine halfbeta_antisymm(funct)
      implicit none
      real(8),intent(inout)                 :: funct(:)
      integer                               :: itau,Ntau
      real(8)                               :: sign=-1d0
      Ntau=size(funct)
      do itau=1,int(Ntau/2)
         funct(itau) = 0.5 * (funct(itau) + sign*funct(ntau-itau+1))
      enddo
      do itau=1,int(Ntau/2)
         funct(ntau-itau+1) = sign*funct(itau)
      enddo
   end subroutine halfbeta_antisymm
   !
   subroutine halfbeta_sym_d(funct,sign)
      implicit none
      real(8),intent(inout)                 :: funct(:)
      real(8),intent(in)                    :: sign
      integer                               :: itau,Ntau
      Ntau=size(funct)
      do itau=1,int(Ntau/2)
         funct(itau) = 0.5 * (funct(itau) + sign*funct(ntau-itau+1))
      enddo
      do itau=1,int(Ntau/2)
         funct(ntau-itau+1) = sign*funct(itau)
      enddo
   end subroutine halfbeta_sym_d
   !
   subroutine halfbeta_sym_z(funct,sign)
      implicit none
      complex(8),intent(inout)              :: funct(:)
      real(8),intent(in)                    :: sign
      integer                               :: itau,Ntau
      Ntau=size(funct)
      do itau=1,int(Ntau/2)
         funct(itau) = 0.5 * (funct(itau) + sign*funct(ntau-itau+1))
      enddo
      do itau=1,int(Ntau/2)
         funct(ntau-itau+1) = sign*funct(itau)
      enddo
   end subroutine halfbeta_sym_z


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the Re/Im part of a causal function on the real frequency
   !axis using Kramers Kronig relatins
   !---------------------------------------------------------------------------!
   subroutine KK_Im2Re(RePart,ImPart,wreal,cutoff,BareVal,ZeroVal,symmetric)
      implicit none
      real(8),intent(inout)                 :: RePart(:)
      real(8),intent(in)                    :: ImPart(:)
      real(8),intent(in)                    :: wreal(:)
      real(8),intent(in)                    :: cutoff
      real(8),intent(in),optional           :: BareVal
      real(8),intent(in),optional           :: ZeroVal
      logical,intent(in),optional           :: symmetric
      !
      real(8),parameter                     :: pi=3.14159265358979323846d0
      integer                               :: iw,iwp,Nreal
      integer                               :: Nwmax,Nwmin,Nwzero
      real(8)                               :: dw,w,wp
      logical                               :: symmetric_
      !
      Nreal=size(wreal)
      call assert_shape(RePart,[Nreal],"KK_Im2Re","RePart")
      call assert_shape(ImPart,[Nreal],"KK_Im2Re","ImPart")
      !
      symmetric_=.false.
      if(present(symmetric))symmetric_=symmetric
      !
      !this is to avoid floating point error with logarithm and for the derivative boundaries
      Nwmin = 5
      Nwmax = Nreal - Nwmin
      Nwzero = minloc(abs(wreal),dim=1)
      !
      dw = abs(wreal(2)-wreal(1))/pi
      !
      RePart=0d0
      if(symmetric_)then
         !
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nwzero,Nwmax,dw,ImPart,cutoff,wreal,RePart),&
         !$OMP PRIVATE(iw,iwp,w,wp)
         !$OMP DO
         do iw=Nwzero,Nwmax
            !
            w  = wreal(iw)
            !
            do iwp=Nwzero,Nwmax
               !
               wp = wreal(iwp)
               !
               if(iw.eq.iwp)then
                  !derivative estimation
                  RePart(iw) = RePart(iw) + (ImPart(iw+1)-ImPart(iw-1))
               else
                  !KK relations
                  RePart(iw) = RePart(iw) + dw * 2d0*(ImPart(iwp)-ImPart(iw)) * w / ( wp**2 - w**2 )
               endif
               !
            enddo
            !
            RePart(iw) = RePart(iw) + ImPart(iw) * log( (cutoff-w)/(cutoff+w) )/pi
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         do iw=Nwzero,Nwmax
            if((Nwzero-iw).lt.Nwmin)exit
            RePart( Nwzero-iw ) = RePart(iw)
         enddo
         !
         !
      else
         !
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nwmin,Nwmax,dw,ImPart,cutoff,wreal,RePart),&
         !$OMP PRIVATE(iw,iwp,w,wp)
         !$OMP DO
         do iw=Nwmin,Nwmax
            !
            w  = wreal(iw)
            !
            do iwp=Nwmin,Nwmax
               !
               wp = wreal(iwp)
               !
               if(iw.eq.iwp)then
                  !derivative estimation
                  RePart(iw) = RePart(iw) + (ImPart(iw+1)-ImPart(iw-1))/2d0
               else
                  !KK relations
                  RePart(iw) = RePart(iw) + dw * (ImPart(iwp)-ImPart(iw)) / ( wp - w )
               endif
               !
            enddo
            !
            RePart(iw) = RePart(iw) + ImPart(iw) * log( (cutoff-w)/(cutoff+w) )/pi
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         !
      endif
      !
      if(present(BareVal))then
         RePart = RePart + BareVal
      elseif(present(ZeroVal))then
         RePart = RePart - RePart(Nwzero) + ZeroVal
      elseif(present(BareVal).and.present(ZeroVal))then
         write(*,"(A)")"     KK_Im2Re: Warning, both bare and screened values provided. Shift ignored."
      endif
      !
   end subroutine KK_Im2Re
   !
   subroutine KK_Re2Im(ImPart,RePart,wreal,cutoff,symmetric)
      implicit none
      real(8),intent(inout)                 :: ImPart(:)
      real(8),intent(in)                    :: RePart(:)
      real(8),intent(in)                    :: wreal(:)
      real(8),intent(in)                    :: cutoff
      logical,intent(in),optional           :: symmetric
      !
      real(8),parameter                     :: pi=3.14159265358979323846d0
      integer                               :: iw,iwp,Nreal
      integer                               :: Nwmax,Nwmin,Nwzero
      real(8)                               :: dw,w,wp
      logical                               :: symmetric_
      !
      Nreal=size(wreal)
      call assert_shape(RePart,[Nreal],"KK_Re2Im","RePart")
      call assert_shape(ImPart,[Nreal],"KK_Re2Im","ImPart")
      !
      symmetric_=.false.
      if(present(symmetric))symmetric_=symmetric
      !
      !this is to avoid floating point error with logarithm
      Nwmin = 5
      Nwmax = Nreal - Nwmin
      Nwzero = minloc(abs(wreal),dim=1)
      !
      dw = abs(wreal(2)-wreal(1))/pi
      !
      ImPart=0d0
      if(symmetric_)then
         !
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nwzero,Nwmax,dw,RePart,cutoff,wreal,ImPart),&
         !$OMP PRIVATE(iw,iwp,w,wp)
         !$OMP DO
         do iw=Nwzero,Nwmax
            !
            w  = wreal(iw)
            !
            do iwp=Nwzero,Nwmax
               !
               wp = wreal(iwp)
               !
               if(iw.eq.iwp)then
                  !derivative estimation
                  ImPart(iw) = ImPart(iw) - (RePart(iw+1)-RePart(iw-1))
               else
                  !KK relations
                  ImPart(iw) = ImPart(iw) - dw * 2d0*(RePart(iwp)-RePart(iw)) * wp / ( wp**2 - w**2 )
               endif
               !
            enddo
            !
            ImPart(iw) = ImPart(iw) - RePart(iw) * log( (cutoff-w)/(cutoff+w) )/pi
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         do iw=Nwzero,Nwmax
            if((Nwzero-iw).lt.Nwmin)exit
            ImPart( Nwzero-iw ) = -ImPart(iw)
         enddo
         !
         !
      else
         !
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nwmin,Nwmax,dw,RePart,cutoff,wreal,ImPart),&
         !$OMP PRIVATE(iw,iwp,w,wp)
         !$OMP DO
         do iw=Nwmin,Nwmax
            !
            w  = wreal(iw)
            !
            do iwp=Nwmin,Nwmax
               !
               wp = wreal(iwp)
               !
               if(iw.eq.iwp)then
                  !derivative estimation
                  ImPart(iw) = ImPart(iw) + (RePart(iw+1)-RePart(iw-1))/2d0
               else
                  !KK relations
                  ImPart(iw) = ImPart(iw) + dw * (RePart(iwp)-RePart(iw)) / ( wp - w )
               endif
               !
            enddo
            !
            ImPart(iw) = ImPart(iw) - RePart(iw) * log( (cutoff-w)/(cutoff+w) )/pi
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         !
      endif
      !
      !
   end subroutine KK_Re2Im


   !---------------------------------------------------------------------------!
   !PURPOSE: compute the moments from 0 to 4 as in mats2itau_FermionicCoeff
   !         works only for diagonal functions
   !---------------------------------------------------------------------------!
   subroutine get_moments_F_d1(moments,funct,beta,wstep,Eo)
      !
      use linalg, only : inv
      implicit none
      !
      real(8),allocatable,intent(inout)     :: moments(:)
      complex(8),intent(in)                 :: funct(:)
      real(8),intent(in)                    :: beta
      integer,intent(in),optional           :: wstep
      logical,intent(in),optional           :: Eo
      !
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: funct_e(3),funct_o(2)
      real(8),allocatable                   :: mom_e(:),tail_e(:,:)
      real(8)                               :: mom_o(2),tail_o(2,2)
      real(8)                               :: w1,w2,w3
      integer                               :: Nmats,iw,wstep_,i,j,Edim
      logical                               :: Eo_
      !
      wstep_=1
      if(present(wstep))wstep_=wstep
      !
      Eo_=.true.
      if(present(Eo))Eo_=Eo
      !
      if(allocated(moments))deallocate(moments)
      if(Eo_)then
         allocate(moments(0:4));moments=0d0
      else
         allocate(moments(1:4));moments=0d0
      endif
      !
      Nmats = size(funct)
      allocate(wmats(Nmats));wmats = FermionicFreqMesh(beta,Nmats)
      w1 = wmats(Nmats-1*wstep_)
      w2 = wmats(Nmats-2*wstep_)
      w3 = wmats(Nmats-3*wstep_)
      deallocate(wmats)
      !
      funct_e=0d0;funct_o=0d0
      do iw=1,3
         funct_e(iw) = dreal(funct(Nmats-iw*wstep_))
      enddo
      do iw=1,2
         funct_o(iw) = dimag(funct(Nmats-iw*wstep_))
      enddo
      !
      ! even tail
      if(Eo_)then
         Edim=3
         allocate(tail_e(Edim,Edim));tail_e=0d0
         allocate(mom_e(Edim));mom_e=0d0
         tail_e(1,1) = 1d0
         tail_e(1,2) = 1d0/(w1**2)
         tail_e(1,3) = 1d0/(w1**4)
         tail_e(2,1) = 1d0
         tail_e(2,2) = 1d0/(w2**2)
         tail_e(2,3) = 1d0/(w2**4)
         tail_e(3,1) = 1d0
         tail_e(3,2) = 1d0/(w3**2)
         tail_e(3,3) = 1d0/(w3**4)
      else
         Edim=2
         allocate(tail_e(Edim,Edim));tail_e=0d0
         allocate(mom_e(Edim));mom_e=0d0
         tail_e(1,1) = 1d0/(w1**2)
         tail_e(1,2) = 1d0/(w1**4)
         tail_e(2,1) = 1d0/(w2**2)
         tail_e(2,2) = 1d0/(w2**4)
      endif
      call inv(tail_e)
      mom_e=0d0
      do i=1,Edim
         do j=1,Edim
            mom_e(i) = mom_e(i) + tail_e(i,j)*funct_e(j)
         enddo
      enddo
      !
      ! odd tail
      tail_o=0d0
      tail_o(1,1) = 1d0/w1
      tail_o(1,2) = 1d0/(w1**3)
      tail_o(2,1) = 1d0/w2
      tail_o(2,2) = 1d0/(w2**3)
      call inv(tail_o)
      mom_o=0d0
      do i=1,2
         do j=1,2
            mom_o(i) = mom_o(i) + tail_o(i,j)*funct_o(j)
         enddo
      enddo
      !
      if(Eo_)then
         moments(0) = mom_e(1)
         moments(1) = mom_o(1)
         moments(2) = mom_e(2)
         moments(3) = mom_o(2)
         moments(4) = mom_e(3)
      else
         moments(1) = mom_o(1)
         moments(2) = mom_e(1)
         moments(3) = mom_o(2)
         moments(4) = mom_e(2)
      endif
      !
   end subroutine get_moments_F_d1
   !
   subroutine get_moments_F_d2(moments,funct,beta,wstep,Eo)
      !
      use linalg, only : inv
      implicit none
      !
      real(8),allocatable,intent(inout)     :: moments(:,:)
      complex(8),intent(in)                 :: funct(:,:)
      real(8),intent(in)                    :: beta
      integer,intent(in),optional           :: wstep
      logical,intent(in),optional           :: Eo
      !
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: funct_e(3),funct_o(2)
      real(8),allocatable                   :: mom_e(:),tail_e(:,:)
      real(8)                               :: mom_o(2),tail_o(2,2)
      real(8)                               :: w1,w2,w3
      integer                               :: Nmats,iw,wstep_,i,j,Edim
      logical                               :: Eo_
      integer                               :: Dim,id
      !
      Dim = size(funct,dim=1)
      !
      wstep_=1
      if(present(wstep))wstep_=wstep
      !
      Eo_=.true.
      if(present(Eo))Eo_=Eo
      !
      if(allocated(moments))deallocate(moments)
      if(Eo_)then
         allocate(moments(Dim,0:4));moments=0d0
      else
         allocate(moments(Dim,1:4));moments=0d0
      endif
      !
      Nmats = size(funct,dim=2)
      allocate(wmats(Nmats));wmats = FermionicFreqMesh(beta,Nmats)
      w1 = wmats(Nmats-1*wstep_)
      w2 = wmats(Nmats-2*wstep_)
      w3 = wmats(Nmats-3*wstep_)
      deallocate(wmats)
      !
      do id=1,Dim
         !
         funct_e=0d0;funct_o=0d0
         do iw=1,3
            funct_e(iw) = dreal(funct(id,Nmats-iw*wstep_))
         enddo
         do iw=1,2
            funct_o(iw) = dimag(funct(id,Nmats-iw*wstep_))
         enddo
         !
         ! even tail
         if(Eo_)then
            Edim=3
            allocate(tail_e(Edim,Edim));tail_e=0d0
            allocate(mom_e(Edim));mom_e=0d0
            tail_e(1,1) = 1d0
            tail_e(1,2) = 1d0/(w1**2)
            tail_e(1,3) = 1d0/(w1**4)
            tail_e(2,1) = 1d0
            tail_e(2,2) = 1d0/(w2**2)
            tail_e(2,3) = 1d0/(w2**4)
            tail_e(3,1) = 1d0
            tail_e(3,2) = 1d0/(w3**2)
            tail_e(3,3) = 1d0/(w3**4)
         else
            Edim=2
            allocate(tail_e(Edim,Edim));tail_e=0d0
            allocate(mom_e(Edim));mom_e=0d0
            tail_e(1,1) = 1d0/(w1**2)
            tail_e(1,2) = 1d0/(w1**4)
            tail_e(2,1) = 1d0/(w2**2)
            tail_e(2,2) = 1d0/(w2**4)
         endif
         call inv(tail_e)
         mom_e=0d0
         do i=1,Edim
            do j=1,Edim
               mom_e(i) = mom_e(i) + tail_e(i,j)*funct_e(j)
            enddo
         enddo
         !
         ! odd tail
         tail_o=0d0
         tail_o(1,1) = 1d0/w1
         tail_o(1,2) = 1d0/(w1**3)
         tail_o(2,1) = 1d0/w2
         tail_o(2,2) = 1d0/(w2**3)
         call inv(tail_o)
         mom_o=0d0
         do i=1,2
            do j=1,2
               mom_o(i) = mom_o(i) + tail_o(i,j)*funct_o(j)
            enddo
         enddo
         !
         if(Eo_)then
            moments(id,0) = mom_e(1)
            moments(id,1) = mom_o(1)
            moments(id,2) = mom_e(2)
            moments(id,3) = mom_o(2)
            moments(id,4) = mom_e(3)
         else
            moments(id,1) = mom_o(1)
            moments(id,2) = mom_e(1)
            moments(id,3) = mom_o(2)
            moments(id,4) = mom_e(2)
         endif
         !
         deallocate(tail_e,mom_e)
         !
      enddo
      !
   end subroutine get_moments_F_d2
   !
   subroutine get_moments_F_d3(moments,funct,beta,wstep,Eo)
      !
      use linalg, only : inv
      implicit none
      !
      real(8),allocatable,intent(inout)     :: moments(:,:,:)
      complex(8),intent(in)                 :: funct(:,:,:)
      real(8),intent(in)                    :: beta
      integer,intent(in),optional           :: wstep
      logical,intent(in),optional           :: Eo
      !
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: funct_e(3),funct_o(2)
      real(8),allocatable                   :: mom_e(:),tail_e(:,:)
      real(8)                               :: mom_o(2),tail_o(2,2)
      real(8)                               :: w1,w2,w3
      integer                               :: Nmats,iw,wstep_,i,j,Edim
      logical                               :: Eo_
      integer                               :: Dim1,id1,Dim2,id2
      !
      Dim1 = size(funct,dim=1)
      Dim2 = size(funct,dim=2)
      !
      wstep_=1
      if(present(wstep))wstep_=wstep
      !
      Eo_=.true.
      if(present(Eo))Eo_=Eo
      !
      if(allocated(moments))deallocate(moments)
      if(Eo_)then
         allocate(moments(Dim1,Dim2,0:4));moments=0d0
      else
         allocate(moments(Dim1,Dim2,1:4));moments=0d0
      endif
      !
      Nmats = size(funct,dim=2)
      allocate(wmats(Nmats));wmats = FermionicFreqMesh(beta,Nmats)
      w1 = wmats(Nmats-1*wstep_)
      w2 = wmats(Nmats-2*wstep_)
      w3 = wmats(Nmats-3*wstep_)
      deallocate(wmats)
      !
      do id1=1,Dim1
         do id2=1,Dim2
            !
            funct_e=0d0;funct_o=0d0
            do iw=1,3
               funct_e(iw) = dreal(funct(id1,id2,Nmats-iw*wstep_))
            enddo
            do iw=1,2
               funct_o(iw) = dimag(funct(id1,id2,Nmats-iw*wstep_))
            enddo
            !
            ! even tail
            if(Eo_)then
               Edim=3
               allocate(tail_e(Edim,Edim));tail_e=0d0
               allocate(mom_e(Edim));mom_e=0d0
               tail_e(1,1) = 1d0
               tail_e(1,2) = 1d0/(w1**2)
               tail_e(1,3) = 1d0/(w1**4)
               tail_e(2,1) = 1d0
               tail_e(2,2) = 1d0/(w2**2)
               tail_e(2,3) = 1d0/(w2**4)
               tail_e(3,1) = 1d0
               tail_e(3,2) = 1d0/(w3**2)
               tail_e(3,3) = 1d0/(w3**4)
            else
               Edim=2
               allocate(tail_e(Edim,Edim));tail_e=0d0
               allocate(mom_e(Edim));mom_e=0d0
               tail_e(1,1) = 1d0/(w1**2)
               tail_e(1,2) = 1d0/(w1**4)
               tail_e(2,1) = 1d0/(w2**2)
               tail_e(2,2) = 1d0/(w2**4)
            endif
            call inv(tail_e)
            mom_e=0d0
            do i=1,Edim
               do j=1,Edim
                  mom_e(i) = mom_e(i) + tail_e(i,j)*funct_e(j)
               enddo
            enddo
            !
            ! odd tail
            tail_o=0d0
            tail_o(1,1) = 1d0/w1
            tail_o(1,2) = 1d0/(w1**3)
            tail_o(2,1) = 1d0/w2
            tail_o(2,2) = 1d0/(w2**3)
            call inv(tail_o)
            mom_o=0d0
            do i=1,2
               do j=1,2
                  mom_o(i) = mom_o(i) + tail_o(i,j)*funct_o(j)
               enddo
            enddo
            !
            if(Eo_)then
               moments(id1,id2,0) = mom_e(1)
               moments(id1,id2,1) = mom_o(1)
               moments(id1,id2,2) = mom_e(2)
               moments(id1,id2,3) = mom_o(2)
               moments(id1,id2,4) = mom_e(3)
            else
               moments(id1,id2,1) = mom_o(1)
               moments(id1,id2,2) = mom_e(1)
               moments(id1,id2,3) = mom_o(2)
               moments(id1,id2,4) = mom_e(2)
            endif
            !
            deallocate(tail_e,mom_e)
            !
         enddo
      enddo
      !
   end subroutine get_moments_F_d3


   !---------------------------------------------------------------------------!
   !PURPOSE: Routines for frequency interpolation using natural splines.
   !Taken from NUMERICAL RECEPIES IN FORTRAN 77 CamUniv Press 1986-1992 p 110
   !Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.,
   !y_i = f( x_i ), with x_1 < x_2 < ... < x_N , this routine returns an
   !array y2(1:n) of length n which contains the second derivatives of the
   !interpolating function at the tabulated points x_i . The routine is signaled
   !to set the corresponding boundary condition for a natural spline, with zero
   !second derivative on that boundary.
   !Parameter: NMAX is the largest anticipated value of n .
   !---------------------------------------------------------------------------!
   subroutine nspline(x,y,y2)
      implicit none
      real(8),intent(in)                    :: x(:)
      real(8),intent(in)                    :: y(:)
      real(8),intent(inout)                 :: y2(:)
      !
      integer                               :: n,i,k
      real(8)                               :: p,qn,sig,un!yp1,ypn,
      real(8),allocatable                   :: u(:)
      !
      n = size(x)
      if(size(y).ne.n) stop "nspline: size(y).ne.size(x)."
      if(size(y2).ne.n) stop "nspline: size(y2).ne.size(x)."
      allocate(u(n));u=0d0
      y2=0d0
      !
      y2(1)=0.
      u(1)=0.
      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      qn=0.
      un=0.
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      deallocate(u)
      !
   end subroutine nspline


   !---------------------------------------------------------------------------!
   !PURPOSE: Given the arrays x_i(1:n) and y_i(1:n) of length n, which tabulate
   !a function (with the x_i is in order), and given the array y2_i(1:n),
   !which is the output from spline above, and given a value of x, this routine
   !returns a cubic-spline interpolated value y .
   !---------------------------------------------------------------------------!
   subroutine splint(x,y,y2,xp,yp)
      implicit none
      real(8),intent(in)                    :: x(:)
      real(8),intent(in)                    :: y(:)
      real(8),intent(in)                    :: y2(:)
      real(8),intent(in)                    :: xp
      real(8),intent(out)                   :: yp
      !
      integer                               :: n,k,khi,klo
      real(8)                               :: a,b,h
      !
      n = size(x)
      if(size(y).ne.n) stop "splint: size(y).ne.size(x)."
      if(size(y2).ne.n) stop "splint: size(y2).ne.size(x)."
      !
      klo=1
      khi=n
      1 if (khi-klo.gt.1) then
           k=(khi+klo)/2
           if(x(k).gt.xp)then
              khi=k
           else
         klo=k
      endif
      goto 1
      endif
      !klo and khi now bracket the input value of x.
      h=x(khi)-x(klo)
      if (h.eq.0.) stop "splint: h.eq.0." !’bad xa input in splint’ The xa’s must be distinct.
      a=(x(khi)-xp)/h
      !Cubic spline polynomial is now evaluated.
      b=(xp-x(klo))/h
      yp=a*y(klo)+b*y(khi)+((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.
      !
   end subroutine splint


   !---------------------------------------------------------------------------!
   !PURPOSE: This is a wrapper that encloses nspline and splint
   !---------------------------------------------------------------------------!
   function cubic_interp(x,y,xp) result(yp)
      implicit none
      real(8),intent(in)                    :: x(:)
      real(8),intent(in)                    :: y(:)
      real(8),intent(in)                    :: xp
      real(8)                               :: yp
      !
      integer                               :: n
      real(8),allocatable                   :: y2(:)
      !
      n=size(x)
      if(size(y).ne.n) stop "nspline: size(y).ne.size(x)."
      allocate(y2(n));y2=0d0
      call nspline(x,y,y2)
      call splint(x,y,y2,xp,yp)
      deallocate(y2)
      !
   end function cubic_interp


   !---------------------------------------------------------------------------!
   !PURPOSE: function which gives the line passing throughout 2 points
   !---------------------------------------------------------------------------!
   function linear_interp_2y(A,B,xp) result(yp)
      implicit none
      real(8),intent(in)                    :: A(2)
      real(8),intent(in)                    :: B(2)
      real(8),intent(in)                    :: xp
      real(8)                               :: yp,xA,yA,xB,yB
      !
      xA = A(1); yA = A(2)
      xB = B(1); yB = B(2)
      !
      if(xA.eq.xB)stop "linear_interp_2y: denominator is zero."
      !
      yp = (xp-xB) * (yA-yB)/(xA-xB) + yB
      !
   end function linear_interp_2y
   !
   function linear_interp_2x(A,B,yp) result(xp)
      implicit none
      real(8),intent(in)                    :: A(2)
      real(8),intent(in)                    :: B(2)
      real(8),intent(in)                    :: yp
      real(8)                               :: xp,xA,yA,xB,yB
      !
      xA = A(1); yA = A(2)
      xB = B(1); yB = B(2)
      !
      if(yA.eq.yB)stop "linear_interp_2x: denominator is zero."
      !
      xp = (yp-yB) * (xA-xB)/(yA-yB) + xB
      !
   end function linear_interp_2x


   !---------------------------------------------------------------------------!
   !PURPOSE: function which computes the integral of an array
   !---------------------------------------------------------------------------!
   function trapezoid_integration(fx,dx) result(Int)
      implicit none
      real(8),dimension(:),intent(in)       :: fx
      real(8),intent(in)                    :: dx
      real(8)                               :: Int
      integer                               :: i
      !
      if(size(fx).le.1) stop " trapezoid_integration: function with wrong dimension."
      !
      Int=0d0
      do i=2,size(fx)
         Int = Int + ( fx(i) + fx(i-1) ) * (dx/2d0)
      enddo
      !
   end function trapezoid_integration


   !---------------------------------------------------------------------------!
   !PURPOSE : flip array
   !---------------------------------------------------------------------------!
   function flip_array_i(Ain) result(Aout)
      implicit none
      integer,dimension(:),intent(in)       :: Ain(:)
      integer,dimension(size(Ain))          :: Aout
      integer                               :: i,N
      !
      N = size(Ain)
      do i=1,N
         Aout(i) = Ain(N-i+1)
      enddo
      !
   end function flip_array_i
   !
   function flip_array_d(Ain) result(Aout)
      implicit none
      real(8),dimension(:)                  :: Ain(:)
      integer,dimension(size(Ain))          :: Aout
      integer                               :: i,N
      !
      N = size(Ain)
      do i=1,N
         Aout(i) = Ain(N-i+1)
      enddo
      !
   end function flip_array_d
   !
   function flip_array_z(Ain) result(Aout)
      implicit none
      complex(8),dimension(:)               :: Ain(:)
      integer,dimension(size(Ain))          :: Aout
      integer                               :: i,N
      !
      N = size(Ain)
      do i=1,N
         Aout(i) = Ain(N-i+1)
      enddo
      !
   end function flip_array_z


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
   !---------------------------------------------------------------------------!
   subroutine FermionicFilon(q,x,fx,wcos,wsin)
      implicit none
      real(8),intent(in)                    :: q
      real(8),intent(in)                    :: x(:)
      real(8),intent(in)                    :: fx(:)
      real(8),intent(out)                   :: wcos
      real(8),intent(out)                   :: wsin
      !
      integer                               :: npoints,nseg
      integer                               :: n,n2
      real(8)                               :: oq,oq2,oq3,h,h2
      real(8)                               :: coskh,sinkh,coskx,sinkx
      real(8)                               :: c0,c2,s1,oh,oh2
      real(8)                               :: a,b,c,bb,yy
      real(8)                               :: expbh1,expbh2,scos,ssin
      logical                               :: abort
      real(8),parameter                     :: precision=1d-9
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
         !
         n2 = 2*n
         h = x(n2) - x(n2-1)
         !
         if(dabs(x(n2+1)-x(n2)-h) .gt. 1d-10) then
            write(*,"(A,I)") "Segment= ",n
            stop "FermionicFilon: the above segment is not equally divided"
         endif
         !
         ! check that fx is not "zero"
         scos = ( fx(n2-1) + 4.d0*fx(n2) + fx(n2+1) )*h/3.d0
         !
         if(dabs(scos).lt.precision) cycle
         !
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
      !
      abort=.false.
      if(wcos.ne.wcos)then
         !NaN condition
         write(*,"(A)")"FermionicFilon: wcos is NaN."
         abort=.true.
      elseif(abs(wcos).ge.huge(1d0))then
         !Infinity condition
         write(*,"(A)")"FermionicFilon: wcos is Infinity."
         abort=.true.
      endif
      if(wsin.ne.wsin)then
         !NaN condition
         write(*,"(A)")"FermionicFilon: wsin is NaN."
         abort=.true.
      elseif(abs(wsin).ge.huge(1d0))then
         !Infinity condition
         write(*,"(A)")"FermionicFilon: wsin is Infinity."
         abort=.true.
      endif
      if(abort) stop "FermionicFilon: coefficient error. Increase precision and recompile."
      !
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
   !---------------------------------------------------------------------------!
   subroutine BosonicFilon(q,x,wcos,wsin)
      implicit none
      real(8),intent(in)                    :: q
      real(8),intent(in)                    :: x(:)
      real(8),intent(inout)                 :: wcos(:)
      real(8),intent(inout)                 :: wsin(:)
      !
      integer                               :: npoints,nseg
      integer                               :: n,n2
      real(8)                               :: oq,oq2,oq3,h,h2,h3
      real(8)                               :: coskh,sinkh,coskx,sinkx
      real(8)                               :: c0,c2,s1,oh,oh2,qh,qh2,qh3,qh4,qh5,qh6
      logical                               :: abort
      real(8),parameter                     :: precision=1d-9
      !
      npoints = size(x)
      if(mod(npoints,2).eq.0) stop "BosonicFilon: npoints is even."
      nseg = (npoints-1)/2
      !
      call assert_shape(wcos,[npoints],"BosonicFilon","wcos")
      call assert_shape(wsin,[npoints],"BosonicFilon","wsin")
      !
      wcos=0d0;wsin=0d0
      !
      if(dabs(q).lt.precision)then
         !
         !Small q
         do n=1,nseg
            !
            n2 = 2 * n
            h = x(n2) - x(n2-1)
            !
            if(dabs(x(n2+1)-x(n2)-h).gt.1d-10) then
               write(*,"(A,I)") "Segment= ",n
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
               write(*,"(A,I)") "Segment= ",n
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
      abort=.false.
      do n=1,npoints
         if(wcos(n).ne.wcos(n))then
            !NaN condition
            write(*,"(A)")"BosonicFilon: wcos["//str(n)//"] is NaN."
            abort=.true.
         elseif(abs(wcos(n)).ge.huge(1d0))then
            !Infinity condition
            write(*,"(A)")"BosonicFilon: wcos["//str(n)//"] is Infinity."
            abort=.true.
         endif
         if(wsin(n).ne.wsin(n))then
            !NaN condition
            write(*,"(A)")"BosonicFilon: wsin["//str(n)//"] is NaN."
            abort=.true.
         elseif(abs(wsin(n)).ge.huge(1d0))then
            !Infinity condition
            write(*,"(A)")"BosonicFilon: wsin["//str(n)//"] is Infinity."
            abort=.true.
         endif
      enddo
      if(abort) stop "BosonicFilon: coefficient error. Increase precision and recompile."
      !
   end subroutine BosonicFilon


   !---------------------------------------------------------------------------!
   !PURPOSE : sort array using random algorithm
   !---------------------------------------------------------------------------!
   subroutine sort_array_i(array,order,replace)
      implicit none
      integer,dimension(:)                    :: array
      integer,dimension(size(array))          :: order
      logical,intent(in),optional             :: replace
      integer,dimension(size(array))          :: backup
      integer                                 :: i
      integer                                 :: lf,rg
      logical                                 :: replace_
      replace_=.false.
      if(present(replace))replace_=replace
      lf=1
      rg=size(array)
      forall(i=1:size(array))order(i)=i
      call qsort_sort(array, order,lf, rg)
      do i=1,size(array)
         backup(i)=array(order(i))
      enddo
      if(replace_)array=backup
   contains
      recursive subroutine qsort_sort( array, order, left, right )
         integer, dimension(:) :: array
         integer, dimension(:) :: order
         integer               :: left
         integer               :: right
         integer               :: i
         integer               :: last
         if ( left .ge. right ) return
         call qsort_swap( order, left, qsort_rand(left,right) )
         last = left
         do i = left+1, right
            if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
               last = last + 1
               call qsort_swap( order, last, i )
            endif
         enddo
         call qsort_swap( order, left, last )
         call qsort_sort( array, order, left, last-1 )
         call qsort_sort( array, order, last+1, right )
      end subroutine qsort_sort
      !---------------------------------------------!
      subroutine qsort_swap( order, first, second )
         integer, dimension(:) :: order
         integer               :: first, second
         integer               :: tmp
         tmp           = order(first)
         order(first)  = order(second)
         order(second) = tmp
      end subroutine qsort_swap
      !---------------------------------------------!
      integer function qsort_rand( lower, upper )
         integer            :: lower, upper
         real(8)               :: r
         call random_number(r)
         qsort_rand =  lower + nint(r * (upper-lower))
      end function qsort_rand
      !---------------------------------------------!
      function compare(f,g)
         implicit none
         integer               :: f,g
         integer               :: compare
         if(f<g) then
            compare=-1
         else
            compare=1
         endif
      end function compare
   end subroutine sort_array_i
   !
   subroutine sort_array_d(array,order,replace)
      implicit none
      real(8),dimension(:)                    :: array
      integer,dimension(size(array))          :: order
      logical,intent(in),optional             :: replace
      integer,dimension(size(array))          :: backup
      integer                                 :: i
      integer                                 :: lf,rg
      logical                                 :: replace_
      replace_=.false.
      if(present(replace))replace_=replace
      lf=1
      rg=size(array)
      forall(i=1:size(array))order(i)=i
      call qsort_sort(array, order,lf, rg)
      do i=1,size(array)
         backup(i)=array(order(i))
      enddo
      if(replace_)array=backup
   contains
      recursive subroutine qsort_sort( array, order, left, right )
         real(8), dimension(:) :: array
         integer, dimension(:) :: order
         integer               :: left
         integer               :: right
         integer               :: i
         integer               :: last
         if ( left .ge. right ) return
         call qsort_swap( order, left, qsort_rand(left,right) )
         last = left
         do i = left+1, right
             if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
                last = last + 1
             call qsort_swap( order, last, i )
             endif
         enddo
         call qsort_swap( order, left, last )
         call qsort_sort( array, order, left, last-1 )
         call qsort_sort( array, order, last+1, right )
      end subroutine qsort_sort
      !---------------------------------------------!
      subroutine qsort_swap( order, first, second )
         integer, dimension(:) :: order
         integer               :: first, second
         integer               :: tmp
         tmp           = order(first)
         order(first)  = order(second)
         order(second) = tmp
      end subroutine qsort_swap
      !---------------------------------------------!
      integer function qsort_rand( lower, upper )
         integer               :: lower, upper
         real(8)               :: r
         call random_number(r)
         qsort_rand =  lower + nint(r * (upper-lower))
      end function qsort_rand
      !---------------------------------------------!
      function compare(f,g)
         implicit none
         real(8)               :: f,g
         integer               :: compare
         if(f<g) then
            compare=-1
         else
            compare=1
         endif
      end function compare
   end subroutine sort_array_d


   !---------------------------------------------------------------------------!
   !PURPOSE: Routines for the str interface
   !---------------------------------------------------------------------------!
   function str_i_to_ch(i4) result(string)
     integer                      :: i4
     character(len=:),allocatable :: string
     character(len=16)            :: string_
     call i4_to_s_left(i4,string_)
     string=trim(adjustl(trim(string_)))
   end function str_i_to_ch
   !
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
   !
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
   !
   !
   !
   subroutine l_assert_shape_N1(A,Ndim,routine,matname)
     logical,dimension(:),intent(in)         :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine l_assert_shape_N1
   subroutine l_assert_shape_N2(A,Ndim,routine,matname)
     logical,dimension(:,:),intent(in)          :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine l_assert_shape_N2
   subroutine l_assert_shape_N3(A,Ndim,routine,matname)
     logical,dimension(:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine l_assert_shape_N3
   subroutine l_assert_shape_N4(A,Ndim,routine,matname)
     logical,dimension(:,:,:,:),intent(in)        :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine l_assert_shape_N4
   subroutine l_assert_shape_N5(A,Ndim,routine,matname)
     logical,dimension(:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine l_assert_shape_N5
   subroutine l_assert_shape_N6(A,Ndim,routine,matname)
     logical,dimension(:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine l_assert_shape_N6
   subroutine l_assert_shape_N7(A,Ndim,routine,matname)
     logical,dimension(:,:,:,:,:,:,:),intent(in)    :: A
     integer,dimension(:),intent(in)            :: Ndim
     character(len=*),optional                  :: routine, matname
     if(any(shape(A) /= Ndim)) then
        if(present(routine).AND.present(matname))&
             write(*,"(A,10I2)")trim(routine)//" error: "//trim(matname)//" has illegal shape"
        stop "assert_shape error: wrong matrix shape"
     end if
   end subroutine l_assert_shape_N7

end module utils_misc
