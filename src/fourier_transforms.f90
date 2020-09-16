module fourier_transforms

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   interface Fmats2itau_mat
      module procedure Fmats2itau_mat_Gw                                        !(beta,Gmats[Norb,Norb,Ntau],Gitau[Norb,Norb,Nmats],asympt_corr,tau_uniform)
      module procedure Fmats2itau_mat_Gwk                                       !(beta,Gmats[Norb,Norb,Ntau,Nkpt],Gitau[Norb,Norb,Nmats,Nkpt],asympt_corr,tau_uniform)
   end interface Fmats2itau_mat
   interface Fmats2itau_vec
      module procedure Fmats2itau_vec_Gw                                        !(beta,Gmats[Norb,Ntau],Gitau[Norb,Nmats],asympt_corr,tau_uniform)
      module procedure Fmats2itau_vec_Gwk                                       !(beta,Gmats[Norb,Ntau,Nkpt],Gitau[Norb,Nmats,Nkpt],asympt_corr,tau_uniform)
   end interface Fmats2itau_vec

   interface Fitau2mats_mat
      module procedure Fitau2mats_mat_Gw                                        !(beta,Gitau[Norb,Norb,Nmats],Gmats[Norb,Norb,Ntau],tau_uniform)
      module procedure Fitau2mats_mat_Gwk                                       !(beta,Gitau[Norb,Norb,Nmats,Nkpt],Gmats[Norb,Norb,Ntau,Nkpt],tau_uniform)
   end interface Fitau2mats_mat
   interface Fitau2mats_vec
      module procedure Fitau2mats_vec_Gw                                        !(beta,Gitau[Norb,Nmats],Gmats[Norb,Ntau],tau_uniform)
      module procedure Fitau2mats_vec_Gwk                                       !(beta,Gitau[Norb,Nmats,Nkpt],Gmats[Norb,Ntau,Nkpt],tau_uniform)
   end interface Fitau2mats_vec

   interface Bmats2itau
      module procedure Bmats2itau_Uw                                            !(beta,Uitau[Nbp,Nbp,Nmats],Umats[Nbp,Nbp,Ntau],asympt_corr,tau_uniform)
      module procedure Bmats2itau_Uwk                                           !(beta,Uitau[Nbp,Nbp,Nmats,Nkpt],Umats[Nbp,Nbp,Ntau,Nkpt],asympt_corr,tau_uniform)
   end interface Bmats2itau
   interface Bitau2mats
      module procedure Bitau2mats_Uw                                            !(beta,Uitau[Nbp,Nbp,Nmats],Umats[Nbp,Nbp,Ntau],tau_uniform)
      module procedure Bitau2mats_Uwk                                           !(beta,Uitau[Nbp,Nbp,Nmats,Nkpt],Umats[Nbp,Nbp,Ntau,Nkpt],tau_uniform)
   end interface Bitau2mats

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),parameter,private                :: pi=3.14159265358979323846d0
   complex(8),parameter,private             :: czero=dcmplx(0.d0,0.d0)

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: Fmats2itau_mat
   public :: Fmats2itau_vec
   public :: Fitau2mats_mat
   public :: Fitau2mats_vec
   public :: Bmats2itau
   public :: Bitau2mats

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates fermionic weights (cos,sin) and add asymptotic correction
   !---------------------------------------------------------------------------!
   ! Ge(iv) = G(iv) + G(-iv), even
   ! Go(iv) = G(iv) - G(-iv), odd
   ! G (iv) = [Ge(iv) + Go(iv)] / 2
   !
   ! G(tau) = (1/beta) S[n=-inf,inf] exp(-ivn*tau) G(ivn)
   !        = S[n=0,niv-1] [ coswt(vn) Ge(ivn) - i*sinwt(vn) Go(ivn) ]
   !
   ! coswt(vn) = cos(vn*tau)/beta
   ! sinwt(vn) = sin(vn*tau)/beta
   !
   ! cos2(tau) = (1/beta) S[n=0,niv-1] cos(vn*tau) / vn^2
   ! sin1(tau) = (1/beta) S[n=0,niv-1] sin(vn*tau) / vn
   ! vn = (2*n+1)*pi/beta
   !
   !-----------------------------
   ! How the weight is calculated:
   !-----------------------------
   ! Construct G(tau) from G(iv)
   ! G(tau) = (1/beta) S[n=-inf,inf] exp[-ivn*tau] G[ivn]
   ! vn     = (2n+1) * pi/beta
   !
   ! G(iv) decays asymptotically as
   ! Ginf(iv) = I[-inf,inf] dw' A(w')/(iv - w')
   !          = (1/iv) I[-inf,inf] dw' A(w') [1 - (w'/iv)]^-1
   !          = (1/iv) I[-inf,inf] dw' A(w') [1 + (w'/iv) + (w'/iv)^2 +...]
   !          = a1/v + a2/v^2 + a3/v^3 + ....
   ! where
   ! a1 = ar1 + i*ai1 etc.
   !
   ! The coefficients a1 etc. are obtained by fitting:
   ! A * a = Ginf  -> a = A^-1 Ginf
   ! where A(i,j)  = 1/v(niv-i)**j
   !       Ginf(i) -> Ginf(niv-i)
   !
   ! G(tau) can be rewritten
   ! G(tau) = (1/beta) S[n=0,niv-1] cos(vn*tau) [Ge(ivn)-Geinf(ivn)]
   !        + (1/beta) S[n=0,inf] cos(vn*tau) Geinf(ivn)
   !        - i * (1/beta) S[n=0,niv-1] sin(vn*tau) [Go(ivn)-Goinf(ivn)]
   !        - i * (1/beta) S[n=0,inf] sin(vn*tau) Goinf(ivn)
   !
   ! From temperature-dependent GW note:
   ! (1/beta) S[n=0,inf] sin(vn*tau) / vn    = 1/4   (beta < tau < 0)
   ! (1/beta) S[n=0,inf] cos(vn*tau) / vn^2  = beta/8 - tau/4
   ! (1/beta) S[n=0,inf] sin(vn*tau) / vn^3  = ( beta*tau - tau^2 ) / 8
   ! (1/beta) S[n=0,inf] cos(vn*tau) / vn^4
   !             = ( beta^3/96 - beta*tau^2/16 + tau^3/24 )
   !
   ! How to use the weights:
   ! ReG(tau)= S[n=0,inf] [ coswt(n) * ReGe(ivn) + sinwt(n) * ImGo(ivn) ]
   ! ImG(tau)= S[n=0,inf] [ coswt(n) * ImGe(ivn) - sinwt(n) * ReGo(ivn) ]
   subroutine mats2itau_FermionicCoeff(tau,coswt,sinwt,correct)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: tau(:)
      real(8),intent(inout)                 :: coswt(:,:)
      real(8),intent(inout)                 :: sinwt(:,:)
      logical,intent(in)                    :: correct
      !
      real(8)                               :: beta
      real(8)                               :: v1,v2,v12
      real(8)                               :: sin1,cos2,sin3,cos4
      real(8),allocatable                   :: wmats(:)
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau
      !
      !
      write(*,"(A)") "--- mats2itau_FermionicCoeff ---"
      !
      !
      Ntau = size(tau)
      Nmats = size(coswt,dim=1)
      call assert_shape(coswt,[Nmats,Ntau],"mats2itau_FermionicCoeff","coswt")
      call assert_shape(sinwt,[Nmats,Ntau],"mats2itau_FermionicCoeff","sinwt")
      beta = tau(Ntau)
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(beta,Nmats)
      v1 = wmats(Nmats-1)
      v2 = wmats(Nmats-2)
      v12 = v1*v1 - v2*v2
      !
      coswt=0d0
      sinwt=0d0
      ! loop over imaginary time tau
      do itau=1,Ntau
         sin1=0d0
         cos2=0d0
         sin3=0d0
         cos4=0d0
         !
         ! cosm(i) = (1/beta) S[n] cos(vn*tau) / vn^i
         ! sinm(i) = (1/beta) S[n] sin(vn*tau) / vn^i
         do iw=1,Nmats
            !
            coswt(iw,itau) = dcos(wmats(iw) * tau(itau)) /beta
            sinwt(iw,itau) = dsin(wmats(iw) * tau(itau)) /beta
            !
            sin1 = sin1 + sinwt(iw,itau) / wmats(iw)
            cos2 = cos2 + coswt(iw,itau) / (wmats(iw)**2)
            sin3 = sin3 + sinwt(iw,itau) / (wmats(iw)**3)
            cos4 = cos4 + coswt(iw,itau) / (wmats(iw)**4)
            !
         enddo
         !
         !
         if(correct)then
            !
            ! (1/beta) S[n=0,inf] cos(vn*tau) / vn^2  = beta/8 - tau/4
            ! (1/beta) S[n=0,inf] cos(vn*tau) / vn^4  = ( beta^3/96 - beta*tau^2/16 + tau^3/24 )
            coswt(Nmats,itau) = coswt(Nmats,itau)                                                                                     &
                              + (v1**4 / v12        ) * ( -cos2 + (beta/8.d0 - tau(itau)/ 4.d0)                                  )    &
                              + (v1**4 * v2**2 / v12) * ( +cos4 - (beta**3/96.d0 -beta*tau(itau)**2/16.d0 + tau(itau)**3/24.d0 ) )
            !
            coswt(Nmats-1,itau) = coswt(Nmats-1,itau)                                                                                 &
                              + (v2**4 / v12        ) * ( +cos2 - (beta/8.d0 - tau(itau)/ 4.d0)                                  )    &
                              + (v2**4 * v1**2 / v12) * ( -cos4 + (beta**3/96.d0 -beta*tau(itau)**2/16.d0 + tau(itau)**3/24.d0 ) )
            !
            ! (1/beta) S[n=0,inf] sin(vn*tau) / vn    = 1/4
            ! (1/beta) S[n=0,inf] sin(vn*tau) / vn^3  = ( beta*tau - tau^2 ) / 8
            sinwt(Nmats,itau) = sinwt(Nmats,itau)                                                              &
                              + (v1**3 / v12)         * ( -sin1 + 1.d0/4.d0                               )    &
                              + (v1**3 * v2**2 / v12) * ( +sin3 - (beta*tau(itau) - tau(itau)**2)/8.d0    )
            !
            sinwt(Nmats-1,itau) = sinwt(Nmats-1,itau)                                                          &
                              + (v1**3 / v12)         * ( +sin1 - 1.d0/4.d0                               )    &
                              + (v1**3 * v2**2 / v12) * ( -sin3 + (beta*tau(itau) - tau(itau)**2)/8.d0    )
            !
         endif
      enddo !itau
      deallocate(wmats)
      !
   end subroutine mats2itau_FermionicCoeff


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates bosonic weights (cos) and add asymptotic correction
   !---------------------------------------------------------------------------!
   ! Construct P(tau) from P(iw)
   ! P(tau) = T S[n=-inf,inf] exp[-iwn*tau] P[iwn]
   ! wn     = 2n * pi/beta,  T=1/beta
   !
   ! P(iw) decays asymptotically as
   ! Pinf(iw) = I[-inf,inf] dw' B(w')/(iw - w')
   !          = (1/iw) I[-inf,inf] dw' B(w') [1 - (w'/iw)]^-1
   !          = (1/iw) I[-inf,inf] dw' B(w') [1 + (w'/iw) + (w'/iw)^2 + ...]
   !          = b1/w^2 + b3/w^4 + ...
   ! where
   ! b1 = - I[-inf,inf] dw' B(w') w'
   ! b3 =   I[-inf,inf] dw' B(w') w'^3
   ! The odd terms vanish due to B being anti-symmetric.
   ! (P(iw) is even but its spectral function B(w) is odd)
   !
   ! P(tau) can be rewritten
   ! P(tau) = T*P(iw0=0)
   !        + T S[n.ne.0] exp[-iwn*tau] [P(iwn)-Pinf(iwn)]
   !        + T S[n.ne.0] exp[-iwn*tau] Pinf(iwn)
   !
   ! The second term decays much faster than P(iwn) alone and the third
   ! term can be calculated analytically.
   !
   ! From Gradshteyn and Ryzhik (1.443-3 and 1.443-4)
   ! S[n.ne.0] exp[-iwn*tau] / wn^2
   ! = S[n.ne.0] cos(wn*tau) / wn^2
   ! = 2 [beta/(2pi)]^2 [pi^2/6 - pi*x/2 + x^2/4],  x=2pi*tau/beta
   !
   ! S[n.ne.0] exp[-iwn*tau] / wn^4
   ! = S[n.ne.0] cos(wn*tau) / wn^4
   ! = 2 [beta/(2pi)]^4 [pi^4/90 - pi^2 x^2/12 + pi*x^3/12 -x^4/48]
   !
   ! To calculate b1 and b3, several options are available.
   ! iopt=1 => b1 and b3 are calculated from the last two points n and n-1.
   !           b3 = [ wn^2 P(iwn) - w(n-1)^2 P(iw(n-1)) ]
   !               / [1/wn^2 - 1/w(n-1)^2]
   !           b1 = wn^2 P(iwn) - b3/wn^2
   ! iopt=10 => b3=0
   !
   ! iopt=2 => b1 and b3 are calculated from two chosen points.
   ! iopt=20=> only b1
   !
   ! iopt=3 => b1 and b3 are calculated from a least-square fit.
   !           Minimise (P-Pinf)^2 wrt b1 and b3:
   !           S[n>N] [P(iwn) - b1/wn^2 - b3/wn^4] / wn^2 = 0
   !           S[n>N] [P(iwn) - b1/wn^2 - b3/wn^4] / wn^4 = 0
   !           N is a mesh point beyond which P is assumed to take its
   !           asymptotic form.
   !
   !           b1 = ( dp1 - bp2) / (ad - bc)
   !           b3 = (-cp1 + ap2) / (ad - bc)
   !           where
   !           a  = S[n>N] 1/wn^4; b=c= S[n>N] 1/wn^6; d = S[n>N]1/wn^8
   !           p1 = S[n>N] P(iwn)/wn^2
   !           p2 = S[n>N] P(iwn)/wn^4
   ! iopt=30=> only b1
   ! Numerical test on Cu suggests that iopt=10 with b3=0 is the best
   ! (simple is best) which is used in this routine.
   subroutine mats2itau_BosonicCoeff(tau,coswt,correct)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: tau(:)
      real(8),intent(inout)                 :: coswt(:,:)
      logical,intent(in)                    :: correct
      !
      real(8)                               :: beta
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: wb1,x
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau
      !
      !
      write(*,"(A)") "--- mats2itau_BosonicCoeff ---"
      !
      !
      Ntau = size(tau)
      Nmats = size(coswt,dim=1)
      call assert_shape(coswt,[Nmats,Ntau],"mats2itau_BosonicCoeff","coswt")
      beta = tau(Ntau)
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = BosonicFreqMesh(beta,Nmats)
      !
      ! P(tau) = (1/beta)*P(iw0=0)
      !        + (1/beta) S[n.ne.0] exp[-iwn*tau] [P(iwn)-Pinf(iwn)]
      !        + (1/beta) S[n.ne.0] exp[-iwn*tau] Pinf(iwn)
      !
      ! For any tau, the weight for n=0 is 1/beta
      ! Since Pinf(iwn) = b1/(wn*wn) and b1=w(niw)^2 P(iw(niw))
      ! Pinf in the 2nd and 3rd line contributes to the weight only for n=niw.
      do itau=1,Ntau
         !
         wb1 = 0d0
         coswt(1,itau) = 1.d0 / beta
         x = (2.d0*pi) * tau(itau) / beta
         !
         do iw=2,Nmats
            coswt(iw,itau)= 2.d0 * dcos (wmats(iw) * tau(itau)) / beta
            wb1 = wb1 - coswt(iw,itau) / (wmats(iw)*wmats(iw))
         enddo
         !
         wb1 = wb1 + 2.d0 * beta * ( 1.d0/(2.d0*pi))**2 *(pi*pi/6.d0 - 0.5d0*pi*x + 0.25d0*x*x)
         if(correct) coswt(Nmats,itau) = coswt(Nmats,itau) + wmats(Nmats)**2 * wb1
         !
      enddo
      !
   end subroutine mats2itau_BosonicCoeff


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from mats to tau of a matrix Gf
   !---------------------------------------------------------------------------!
   subroutine Fmats2itau_mat_Gw(beta,Gmats,Gitau,asympt_corr,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gmats(:,:,:)
      complex(8),intent(inout)              :: Gitau(:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: coswt(:,:),sinwt(:,:)
      real(8),allocatable                   :: tau(:)
      complex(8),allocatable                :: Ge(:,:),Go(:,:)
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau,Norb
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Fmats2itau_mat_Gw ---"
      !
      !
      Norb = size(Gmats,dim=1)
      Nmats = size(Gmats,dim=3)
      Ntau = size(Gitau,dim=3)
      if(size(Gmats,dim=1).ne.size(Gmats,dim=2)) stop "Gmats not square."
      call assert_shape(Gitau,[Norb,Norb,Ntau],"Fmats2itau_mat_Gw","Gitau")
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      allocate(coswt(Nmats,Ntau));coswt=0d0
      allocate(sinwt(Nmats,Ntau));sinwt=0d0
      call mats2itau_FermionicCoeff(tau,coswt,sinwt,asympt_corr_)
      !
      Gitau=czero
      allocate(Ge(Norb,Norb));Ge=czero
      allocate(Go(Norb,Norb));Go=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ntau,Nmats,Gmats,coswt,sinwt,Gitau),&
      !$OMP PRIVATE(itau,iw,Ge,Go)
      !$OMP DO
      do itau=1,Ntau
         do iw=1,Nmats
            !
            ! Gab(iw) = Gba*(-iwn) --> Gab(iw) = Gba*(-iwn)
            Ge = Gmats(:,:,iw) + transpose(conjg(Gmats(:,:,iw)))
            Go = Gmats(:,:,iw) - transpose(conjg(Gmats(:,:,iw)))
            !
            Gitau(:,:,itau) = Gitau(:,:,itau) + coswt(iw,itau)*Ge -dcmplx(0d0,1d0)*sinwt(iw,itau)*Go
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,sinwt,Ge,Go,tau)
      !
   end subroutine Fmats2itau_mat_Gw
   !
   subroutine Fmats2itau_mat_Gwk(beta,Gmats,Gitau,asympt_corr,tau_uniform,nkpt3,kpt)
      !
      use utils_misc
      use crystal
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in),target          :: Gmats(:,:,:,:)
      complex(8),intent(inout),target       :: Gitau(:,:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      integer,intent(in),optional           :: nkpt3(3)
      real(8),intent(in),optional           :: kpt(:,:)
      !
      real(8),allocatable                   :: coswt(:,:),sinwt(:,:)
      real(8),allocatable                   :: tau(:)
      complex(8),allocatable,target         :: Gmats_rs(:,:,:,:),Gitau_rs(:,:,:,:)
      complex(8),pointer                    :: Gft_in(:,:,:,:),Gft_out(:,:,:,:)
      complex(8),allocatable                :: Ge(:,:),Go(:,:)
      integer                               :: iw,itau,idat
      integer,pointer                       :: Ndat
      integer,target                        :: Nwig
      integer,target                        :: Nkpt
      integer                               :: Ntau,Norb,Nmats
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      logical                               :: real_space
      !
      !
      write(*,"(A)") "--- Fmats2itau_mat_Gwk ---"
      !
      !
      Norb = size(Gmats,dim=1)
      Nmats = size(Gmats,dim=3)
      Nkpt = size(Gmats,dim=4)
      Ntau = size(Gitau,dim=3)
      if(size(Gmats,dim=1).ne.size(Gmats,dim=2)) stop "Gmats not square."
      call assert_shape(Gitau,[Norb,Norb,Ntau,Nkpt],"Fmats2itau_mat_Gwk","Gitau")
      !
      real_space=.false.
      if(present(nkpt3).and.present(kpt))then
         real_space=.true.
      else
         stop "Either kpt or nkpt3 argument is missing."
      endif
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      allocate(coswt(Nmats,Ntau));coswt=0d0
      allocate(sinwt(Nmats,Ntau));sinwt=0d0
      call mats2itau_FermionicCoeff(tau,coswt,sinwt,asympt_corr_)
      !
      !connect to the wanted Gf
      if(real_space)then
         !
         call wannier_K2R(nkpt3,kpt,Gmats,Gmats_rs)
         Nwig = size(Gmats_rs,dim=4)
         allocate(Gitau_rs(Norb,Norb,Ntau,Nwig));Gitau_rs=czero
         !
         Gft_in  => Gmats_rs
         Gft_out => Gitau_rs
         Ndat => Nwig
         !
      else
         !
         Gft_in  => Gmats
         Gft_out => Gitau
         Ndat => Nkpt
         !
      endif
      !
      Gitau=czero
      allocate(Ge(Norb,Norb));Ge=czero
      allocate(Go(Norb,Norb));Go=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ndat,Ntau,Nmats,coswt,sinwt,Gft_in,Gft_out),&
      !$OMP PRIVATE(idat,itau,iw,Ge,Go)
      !$OMP DO
      do idat=1,Ndat
         do itau=1,Ntau
            do iw=1,Nmats
               !
               ! Gab(iw) = Gba*(-iwn) --> Gab(iw) = Gba*(-iwn)
               Ge = Gft_in(:,:,iw,idat) + transpose(conjg(Gft_in(:,:,iw,idat)))
               Go = Gft_in(:,:,iw,idat) - transpose(conjg(Gft_in(:,:,iw,idat)))
               !
               Gft_out(:,:,itau,idat) = Gft_out(:,:,itau,idat) + coswt(iw,itau)*Ge -dcmplx(0d0,1d0)*sinwt(iw,itau)*Go
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,sinwt,Ge,Go,tau)
      !
      if(real_space)then
         deallocate(Gmats_rs)
         call wannier_R2K(nkpt3,kpt,Gitau_rs,Gitau)
         deallocate(Gitau_rs)
      endif
      nullify(Gft_in,Gft_out,Ndat)
      !
   end subroutine Fmats2itau_mat_Gwk


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from mats to tau of a vector Gf
   !---------------------------------------------------------------------------!
   subroutine Fmats2itau_vec_Gw(beta,Gmats,Gitau,asympt_corr,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gmats(:,:)
      complex(8),intent(inout)              :: Gitau(:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: coswt(:,:),sinwt(:,:)
      real(8),allocatable                   :: tau(:)
      complex(8),allocatable                :: Ge(:),Go(:)
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau,Norb
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Fmats2itau_vec_Gw ---"
      !
      !
      Norb = size(Gmats,dim=1)
      Nmats = size(Gmats,dim=2)
      Ntau = size(Gitau,dim=2)
      call assert_shape(Gitau,[Norb,Ntau],"Fmats2itau_vec_Gw","Gitau")
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      allocate(coswt(Nmats,Ntau));coswt=0d0
      allocate(sinwt(Nmats,Ntau));sinwt=0d0
      call mats2itau_FermionicCoeff(tau,coswt,sinwt,asympt_corr_)
      !
      Gitau=czero
      allocate(Ge(Norb));Ge=czero
      allocate(Go(Norb));Go=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ntau,Nmats,Gmats,coswt,sinwt,Gitau),&
      !$OMP PRIVATE(itau,iw,Ge,Go)
      !$OMP DO
      do itau=1,Ntau
         do iw=1,Nmats
            !
            ! Gab(iw) = Gba*(-iwn) --> Gab(-iw) = Gba*(iwn)
            Ge = Gmats(:,iw) + conjg(Gmats(:,iw))
            Go = Gmats(:,iw) - conjg(Gmats(:,iw))
            !
            Gitau(:,itau) = Gitau(:,itau) + coswt(iw,itau)*Ge -dcmplx(0d0,1d0)*sinwt(iw,itau)*Go
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,sinwt,Ge,Go,tau)
      !
   end subroutine Fmats2itau_vec_Gw
   !
   subroutine Fmats2itau_vec_Gwk(beta,Gmats,Gitau,asympt_corr,tau_uniform,nkpt3,kpt)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in),target          :: Gmats(:,:,:)
      complex(8),intent(inout),target       :: Gitau(:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      integer,intent(in),optional           :: nkpt3(3)
      real(8),intent(in),optional           :: kpt(:,:)
      !
      real(8),allocatable                   :: coswt(:,:),sinwt(:,:)
      real(8),allocatable                   :: tau(:)
      complex(8),allocatable,target         :: Gmats_rs(:,:,:),Gitau_rs(:,:,:)
      complex(8),pointer                    :: Gft_in(:,:,:),Gft_out(:,:,:)
      complex(8),allocatable                :: Ge(:),Go(:)
      integer                               :: iw,itau,idat
      integer,pointer                       :: Ndat
      integer,target                        :: Nwig
      integer,target                        :: Nkpt
      integer                               :: Ntau,Norb,Nmats
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      logical                               :: real_space
      !
      !
      write(*,"(A)") "--- Fmats2itau_vec_Gwk ---"
      !
      !
      Norb = size(Gmats,dim=1)
      Nmats = size(Gmats,dim=2)
      Nkpt = size(Gmats,dim=3)
      Ntau = size(Gitau,dim=2)
      call assert_shape(Gitau,[Norb,Ntau,Nkpt],"Fmats2itau_vec_Gwk","Gitau")
      !
      real_space=.false.
      if(present(nkpt3).and.present(kpt))then
         real_space=.true.
      else
         stop "Either kpt or nkpt3 argument is missing."
      endif
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      allocate(coswt(Nmats,Ntau));coswt=0d0
      allocate(sinwt(Nmats,Ntau));sinwt=0d0
      call mats2itau_FermionicCoeff(tau,coswt,sinwt,asympt_corr_)
      !
      !connect to the wanted Gf
      if(real_space)then
         !
         call wannier_K2R(nkpt3,kpt,Gmats,Gmats_rs)
         Nwig = size(Gmats_rs,dim=3)
         allocate(Gitau_rs(Norb,Ntau,Nwig));Gitau_rs=czero
         !
         Gft_in  => Gmats_rs
         Gft_out => Gitau_rs
         Ndat => Nwig
         !
      else
         !
         Gft_in  => Gmats
         Gft_out => Gitau
         Ndat => Nkpt
         !
      endif
      !
      Gitau=czero
      allocate(Ge(Norb));Ge=czero
      allocate(Go(Norb));Go=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ndat,Ntau,Nmats,coswt,sinwt,Gft_in,Gft_out),&
      !$OMP PRIVATE(idat,itau,iw,Ge,Go)
      !$OMP DO
      do idat=1,Ndat
         do itau=1,Ntau
            do iw=1,Nmats
               !
               ! Gab(iw) = Gba*(-iwn) --> Gab(-iw) = Gba*(iwn)
               Ge = Gft_in(:,iw,idat) + conjg(Gft_in(:,iw,idat))
               Go = Gft_in(:,iw,idat) - conjg(Gft_in(:,iw,idat))
               !
               Gft_out(:,itau,idat) = Gft_out(:,itau,idat) + coswt(iw,itau)*Ge -dcmplx(0d0,1d0)*sinwt(iw,itau)*Go
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,sinwt,Ge,Go,tau)
      !
      if(real_space)then
         deallocate(Gmats_rs)
         call wannier_R2K(nkpt3,kpt,Gitau_rs,Gitau)
         deallocate(Gitau_rs)
      endif
      nullify(Gft_in,Gft_out,Ndat)
      !
   end subroutine Fmats2itau_vec_Gwk


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from tau to mats of a matrix Gf
   !---------------------------------------------------------------------------!
   subroutine Fitau2mats_mat_Gw(beta,Gitau,Gmats,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gitau(:,:,:)
      complex(8),intent(inout)              :: Gmats(:,:,:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: RealG(:),ImagG(:)
      real(8)                               :: rwcos,rwsin,cwcos,cwsin
      integer                               :: iw,iwan1,iwan2
      integer                               :: Nmats,Ntau,Norb
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Fitau2mats_mat_Gw ---"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=3)
      Nmats = size(Gmats,dim=3)
      if(size(Gitau,dim=1).ne.size(Gitau,dim=2)) stop "Gitau not square."
      call assert_shape(Gmats,[Norb,Norb,Nmats],"Fitau2mats_mat_Gw","Gmats")
      if(mod(Ntau,2).eq.0) stop "Required Filon routines are not working with odd segments."
      !
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(beta,Nmats)
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      Gmats=czero
      allocate(RealG(Ntau));RealG=0d0
      allocate(ImagG(Ntau));ImagG=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nmats,Norb,wmats,tau,Gitau,Gmats),&
      !$OMP PRIVATE(iw,iwan1,iwan2,RealG,ImagG,rwcos,rwsin,cwcos,cwsin)
      !$OMP DO
      do iw=1,Nmats
         do iwan1=1,Norb
            do iwan2=1,Norb
               !
               RealG = dreal( Gitau(iwan1,iwan2,:) )
               ImagG = dimag( Gitau(iwan1,iwan2,:) )
               !
               call FermionicFilon(wmats(iw),tau,RealG,rwcos,rwsin)
               call FermionicFilon(wmats(iw),tau,ImagG,cwcos,cwsin)
               !
               Gmats(iwan1,iwan2,iw) = dcmplx(rwcos-cwsin,0d0) + dcmplx(0d0,rwsin+cwcos)
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(RealG,ImagG,tau,wmats)
      !
   end subroutine Fitau2mats_mat_Gw
   !
   subroutine Fitau2mats_mat_Gwk(beta,Gitau,Gmats,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gitau(:,:,:,:)
      complex(8),intent(inout)              :: Gmats(:,:,:,:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: RealG(:),ImagG(:)
      real(8)                               :: rwcos,rwsin,cwcos,cwsin
      integer                               :: iw,iwan1,iwan2,ik
      integer                               :: Nmats,Ntau,Norb,Nkpt
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Fitau2mats_mat_Gwk ---"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=3)
      Nkpt = size(Gitau,dim=4)
      Nmats = size(Gmats,dim=3)
      if(size(Gitau,dim=1).ne.size(Gitau,dim=2)) stop "Gitau not square."
      call assert_shape(Gmats,[Norb,Norb,Nmats,Nkpt],"Fitau2mats_mat_Gwk","Gmats")
      if(mod(Ntau,2).eq.0) stop "Required Filon routines are not working with odd segments."
      !
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(beta,Nmats)
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      Gmats=czero
      allocate(RealG(Ntau));RealG=0d0
      allocate(ImagG(Ntau));ImagG=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Nmats,Norb,wmats,tau,Gitau,Gmats),&
      !$OMP PRIVATE(ik,iw,iwan1,iwan2,RealG,ImagG,rwcos,rwsin,cwcos,cwsin)
      !$OMP DO
      do ik=1,Nkpt
         do iw=1,Nmats
            do iwan1=1,Norb
               do iwan2=1,Norb
                  !
                  RealG = dreal( Gitau(iwan1,iwan2,:,ik) )
                  ImagG = dimag( Gitau(iwan1,iwan2,:,ik) )
                  !
                  call FermionicFilon(wmats(iw),tau,RealG,rwcos,rwsin)
                  call FermionicFilon(wmats(iw),tau,ImagG,cwcos,cwsin)
                  !
                  Gmats(iwan1,iwan2,iw,ik) = dcmplx(rwcos-cwsin,0d0) + dcmplx(0d0,rwsin+cwcos)
                  !
               enddo
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(RealG,ImagG,tau,wmats)
      !
   end subroutine Fitau2mats_mat_Gwk


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from tau to mats of a vector Gf
   !---------------------------------------------------------------------------!
   subroutine Fitau2mats_vec_Gw(beta,Gitau,Gmats,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gitau(:,:)
      complex(8),intent(inout)              :: Gmats(:,:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: RealG(:),ImagG(:)
      real(8)                               :: rwcos,rwsin,cwcos,cwsin
      integer                               :: iw,iwan1
      integer                               :: Nmats,Ntau,Norb
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Fitau2mats_vec_Gw ---"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=2)
      Nmats = size(Gmats,dim=2)
      call assert_shape(Gmats,[Norb,Nmats],"Fitau2mats_vec_Gw","Gmats")
      if(mod(Ntau,2).eq.0) stop "Required Filon routines are not working with odd segments."
      !
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(beta,Nmats)
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      Gmats=czero
      allocate(RealG(Ntau));RealG=0d0
      allocate(ImagG(Ntau));ImagG=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nmats,Norb,wmats,tau,Gitau,Gmats),&
      !$OMP PRIVATE(iw,iwan1,RealG,ImagG,rwcos,rwsin,cwcos,cwsin)
      !$OMP DO
      do iw=1,Nmats
         do iwan1=1,Norb
            !
            RealG = dreal( Gitau(iwan1,:) )
            ImagG = dimag( Gitau(iwan1,:) )
            !
            call FermionicFilon(wmats(iw),tau,RealG,rwcos,rwsin)
            call FermionicFilon(wmats(iw),tau,ImagG,cwcos,cwsin)
            !
            Gmats(iwan1,iw) = dcmplx(rwcos-cwsin,0d0) + dcmplx(0d0,rwsin+cwcos)
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(RealG,ImagG,tau,wmats)
      !
   end subroutine Fitau2mats_vec_Gw
   !
   subroutine Fitau2mats_vec_Gwk(beta,Gitau,Gmats,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gitau(:,:,:)
      complex(8),intent(inout)              :: Gmats(:,:,:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: RealG(:),ImagG(:)
      real(8)                               :: rwcos,rwsin,cwcos,cwsin
      integer                               :: iw,iwan1,ik
      integer                               :: Nmats,Ntau,Norb,Nkpt
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Fitau2mats_vec_Gwk ---"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=2)
      Nkpt = size(Gitau,dim=3)
      Nmats = size(Gmats,dim=2)
      call assert_shape(Gmats,[Norb,Nmats,Nkpt],"Fitau2mats_vec_Gwk","Gmats")
      if(mod(Ntau,2).eq.0) stop "Required Filon routines are not working with odd segments."
      !
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(beta,Nmats)
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      Gmats=czero
      allocate(RealG(Ntau));RealG=0d0
      allocate(ImagG(Ntau));ImagG=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Nmats,Norb,wmats,tau,Gitau,Gmats),&
      !$OMP PRIVATE(ik,iw,iwan1,RealG,ImagG,rwcos,rwsin,cwcos,cwsin)
      !$OMP DO
      do ik=1,Nkpt
         do iw=1,Nmats
            do iwan1=1,Norb
               !
               RealG = dreal( Gitau(iwan1,:,ik) )
               ImagG = dimag( Gitau(iwan1,:,ik) )
               !
               call FermionicFilon(wmats(iw),tau,RealG,rwcos,rwsin)
               call FermionicFilon(wmats(iw),tau,ImagG,cwcos,cwsin)
               !
               Gmats(iwan1,iw,ik) = dcmplx(rwcos-cwsin,0d0) + dcmplx(0d0,rwsin+cwcos)
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(RealG,ImagG,tau,wmats)
      !
   end subroutine Fitau2mats_vec_Gwk


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from mats to tau of a bosonic tensor
   !---------------------------------------------------------------------------!
   subroutine Bmats2itau_Uw(beta,Umats,Uitau,asympt_corr,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Umats(:,:,:)
      complex(8),intent(inout)              :: Uitau(:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: coswt(:,:)
      real(8),allocatable                   :: tau(:)
      integer                               :: iw,itau,ib1,ib2
      integer                               :: Nmats,Ntau,Nbp
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Bmats2itau_Uw ---"
      !
      !
      Nbp = size(Umats,dim=1)
      Nmats = size(Umats,dim=3)
      Ntau = size(Uitau,dim=3)
      if(size(Umats,dim=1).ne.size(Umats,dim=2)) stop "Umats not square."
      call assert_shape(Uitau,[Nbp,Nbp,Ntau],"Bmats2itau_Uw","Uitau")
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      allocate(coswt(Nmats,Ntau));coswt=0d0
      call mats2itau_BosonicCoeff(tau,coswt,asympt_corr_)
      !
      Uitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ntau,Nmats,Nbp,Umats,coswt,Uitau),&
      !$OMP PRIVATE(itau,iw,ib1,ib2)
      !$OMP DO
      do itau=1,Ntau
         do iw=1,Nmats
            !
            do ib1=1,Nbp
               do ib2=1,Nbp
                  Uitau(ib1,ib2,itau) = Uitau(ib1,ib2,itau) + coswt(iw,itau)*Umats(ib1,ib2,itau)
               enddo
            enddo
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,tau)
      !
   end subroutine Bmats2itau_Uw
   !
   subroutine Bmats2itau_Uwk(beta,Umats,Uitau,asympt_corr,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Umats(:,:,:,:)
      complex(8),intent(inout)              :: Uitau(:,:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: coswt(:,:)
      real(8),allocatable                   :: tau(:)
      integer                               :: iw,itau,ib1,ib2,ik
      integer                               :: Nmats,Ntau,Nbp,Nkpt
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Bmats2itau_Uwk ---"
      !
      !
      Nbp = size(Umats,dim=1)
      Nmats = size(Umats,dim=3)
      Nkpt = size(Umats,dim=4)
      Ntau = size(Uitau,dim=3)
      if(size(Umats,dim=1).ne.size(Umats,dim=2)) stop "Umats not square."
      call assert_shape(Uitau,[Nbp,Nbp,Ntau,Nkpt],"Bmats2itau_Uwk","Uitau")
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      allocate(coswt(Nmats,Ntau));coswt=0d0
      call mats2itau_BosonicCoeff(tau,coswt,asympt_corr_)
      !
      Uitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Ntau,Nmats,Nbp,Umats,coswt,Uitau),&
      !$OMP PRIVATE(ik,itau,iw,ib1,ib2)
      !$OMP DO
      do ik=1,Nkpt
         do itau=1,Ntau
            do iw=1,Nmats
               !
               do ib1=1,Nbp
                  do ib2=1,Nbp
                     Uitau(ib1,ib2,itau,ik) = Uitau(ib1,ib2,itau,ik) + coswt(iw,itau)*Umats(ib1,ib2,itau,ik)
                  enddo
               enddo
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,tau)
      !
   end subroutine Bmats2itau_Uwk


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from tau to mats of a bosonic tensor
   !---------------------------------------------------------------------------!
   subroutine Bitau2mats_Uw(beta,Uitau,Umats,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Uitau(:,:,:)
      complex(8),intent(inout)              :: Umats(:,:,:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: wcos(:),wsin(:)
      real(8)                               :: RealU,ImagU
      integer                               :: iw,itau,ib1,ib2
      integer                               :: Nmats,Ntau,Nbp
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Bitau2mats_Uw ---"
      !
      !
      Nbp = size(Uitau,dim=1)
      Ntau = size(Uitau,dim=3)
      Nmats = size(Umats,dim=3)
      if(size(Uitau,dim=1).ne.size(Uitau,dim=2)) stop "Uitau not square."
      call assert_shape(Umats,[Nbp,Nbp,Nmats],"Bitau2mats_Uw","Umats")
      if(mod(Ntau,2).eq.0) stop "Required Filon routines are not working with odd segments."
      !
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = BosonicFreqMesh(beta,Nmats)
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      Umats=czero
      allocate(wcos(Ntau));wcos=0d0
      allocate(wsin(Ntau));wsin=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nmats,Nbp,Ntau,wmats,tau,Uitau,Umats),&
      !$OMP PRIVATE(iw,ib1,ib2,itau,RealU,ImagU,wcos,wsin)
      !$OMP DO
      do iw=1,Nmats
         !
         call BosonicFilon(wmats(iw),tau,wcos,wsin)
         !
         do ib1=1,Nbp
            do ib2=1,Nbp
               do itau=1,Ntau
                  !
                  RealU = dreal( Uitau(ib1,ib2,itau) )
                  ImagU = dimag( Uitau(ib1,ib2,itau) )
                  !
                  Umats(ib1,ib2,iw) =   Umats(ib1,ib2,iw)                            &
                                    + ( RealU * wcos(itau) - ImagU * wsin(itau) )    &
                                    + ( ImagU * wcos(itau) + RealU * wsin(itau) )*dcmplx(0d0,1d0)
                  !
               enddo
            enddo
         enddo
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(wcos,wsin,tau,wmats)
      !
   end subroutine Bitau2mats_Uw
   !
   subroutine Bitau2mats_Uwk(beta,Uitau,Umats,tau_uniform)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Uitau(:,:,:,:)
      complex(8),intent(inout)              :: Umats(:,:,:,:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: wcos(:),wsin(:)
      real(8)                               :: RealU,ImagU
      integer                               :: iw,itau,ib1,ib2,ik
      integer                               :: Nmats,Ntau,Nbp,Nkpt
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- Bitau2mats_Uwk ---"
      !
      !
      Nbp = size(Uitau,dim=1)
      Ntau = size(Uitau,dim=3)
      Nkpt = size(Uitau,dim=4)
      Nmats = size(Umats,dim=3)
      if(size(Uitau,dim=1).ne.size(Uitau,dim=2)) stop "Uitau not square."
      call assert_shape(Umats,[Nbp,Nbp,Nmats,Nkpt],"Bitau2mats_Uwk","Umats")
      if(mod(Ntau,2).eq.0) stop "Required Filon routines are not working with odd segments."
      !
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = BosonicFreqMesh(beta,Nmats)
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      Umats=czero
      allocate(wcos(Ntau));wcos=0d0
      allocate(wsin(Ntau));wsin=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nmats,Nkpt,Nbp,Ntau,wmats,tau,Uitau,Umats),&
      !$OMP PRIVATE(iw,ik,ib1,ib2,itau,RealU,ImagU,wcos,wsin)
      !$OMP DO
      do iw=1,Nmats
         !
         call BosonicFilon(wmats(iw),tau,wcos,wsin)
         !
         do ik=1,Nkpt
            do ib1=1,Nbp
               do ib2=1,Nbp
                  do itau=1,Ntau
                     !
                     RealU = dreal( Uitau(ib1,ib2,itau,ik) )
                     ImagU = dimag( Uitau(ib1,ib2,itau,ik) )
                     !
                     Umats(ib1,ib2,iw,ik) =   Umats(ib1,ib2,iw,ik)                           &
                                            + ( RealU * wcos(itau) - ImagU * wsin(itau) )    &
                                            + ( ImagU * wcos(itau) + RealU * wsin(itau) )*dcmplx(0d0,1d0)
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(wcos,wsin,tau,wmats)
      !
   end subroutine Bitau2mats_Uwk


end module fourier_transforms
