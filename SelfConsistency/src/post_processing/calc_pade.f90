function pade(funct_in,type,wlimit) result(funct_out)
   !
   use parameters
   use utils_misc
   use input_vars, only: Nmats, Nreal, wrealMax, Beta
   implicit none
   !
   complex(8),intent(in)                 :: funct_in(:)
   character(len=*),intent(in)           :: type
   integer,intent(in),optional           :: wlimit
   complex(8),dimension(Nreal)           :: funct_out
   !
   integer                               :: Nfreq
   real(8),allocatable                   :: wreal(:),wmats(:)
   !
   !
   if(verbose)write(*,"(A)") "---- pade"
   !
   !
   Nfreq = Nmats
   if(present(wlimit))Nfreq = wlimit
   !
   allocate(wreal(Nreal));wreal=0d0
   wreal = linspace(-wrealMax,+wrealMax,Nreal)
   allocate(wmats(Nfreq));wmats=0d0
   select case(reg(type))
      case default
         stop "Available Modes are: Fermionic, Bosonic."
      case("Fermionic")
         wmats = FermionicFreqMesh(Beta,Nfreq)
      case("Bosonic")
         wmats = BosonicFreqMesh(Beta,Nfreq)
   end select
   !
   funct_out=czero
   call padecoeff(funct_out, wreal+img*0d0, funct_in(1:Nfreq), img*wmats)
   !
end function pade

subroutine padecoeff(fwout,wout,fwin,win)
   !
   use parameters
   implicit none
   !
   complex(8),intent(out)                :: fwout(:)
   complex(8),intent(in)                 :: wout(:)
   complex(8),intent(in)                 :: fwin(:)
   complex(8),intent(in)                 :: win(:)
   !
   complex(8),allocatable                :: coeff(:,:)
   complex(8),allocatable                :: a(:),b(:)
   integer                               :: Nin,Nout,i,j
   !
   !
   if(verbose)write(*,"(A)") "---- padecoeff"
   !
   !
   Nout = size(fwout)
   if(size(wout).ne.Nout) stop "padecoeff: size(wout).ne.Nout"
   !
   Nin = size(fwin)
   if(size(win).ne.Nin) stop "padecoeff: size(win).ne.Nin"
   !
   allocate(coeff(Nin,Nin));coeff=czero
   do j=1,Nin
      coeff(1,j)=fwin(j)
   enddo
   do j=2,Nin
      do i=2,j
         if (abs(win(j)-win(i-1)).lt.1.D-10) then
            stop "pade z=0"
         endif
         if (abs(coeff(i-1,j)).lt.1.D-10) then
            write(*,*)"i,j,coeff(i-1,j)",i,j,coeff(i-1,j)
            stop "coeff=0"
         endif
         coeff(i,j) = (coeff(i-1,i-1)-coeff(i-1,j)) / (win(j)-win(i-1)) / coeff(i-1,j)
      enddo
   enddo
   !
   fwout = czero
   allocate(a(0:Nin))
   allocate(b(0:Nin))
   do j=1,Nout
      !
      a=czero
      b=czero
      !
      a(0)=0.d0
      a(1)=coeff(1,1)
      b(0)=1.d0
      b(1)=1.d0
      !
      do i=1,Nin-1
         a(i+1)=a(i)+(wout(j)-win(i))*coeff(i+1,i+1)*a(i-1)
         b(i+1)=b(i)+(wout(j)-win(i))*coeff(i+1,i+1)*b(i-1)
      enddo
      !
      fwout(j) = a(Nin)/b(Nin)
      if(fwout(j) .ne. fwout(j)) fwout(j)=czero       !remove NaN
      if(abs(fwout(j)) .ge. huge(1d0)) fwout(j)=czero !remove Infinity
      !
   enddo
   !
end subroutine padecoeff
