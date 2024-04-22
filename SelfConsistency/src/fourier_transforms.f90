module fourier_transforms

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface Fmats2itau_mat
      module procedure Fmats2itau_mat_Gw                                        !(beta,Gmats[Norb,Norb,Ntau],Gitau[Norb,Norb,Nmats],asympt_corr,tau_uniform)
     !module procedure Fmats2itau_mat_Gwk                                       !(beta,Gmats[Norb,Norb,Ntau,Nkpt],Gitau[Norb,Norb,Nmats,Nkpt],asympt_corr,tau_uniform)
   end interface Fmats2itau_mat
   interface Fmats2itau_vec
      module procedure Fmats2itau_vec_Gw                                        !(beta,Gmats[Norb,Ntau],Gitau[Norb,Nmats],asympt_corr,tau_uniform)
     !module procedure Fmats2itau_vec_Gwk                                       !(beta,Gmats[Norb,Ntau,Nkpt],Gitau[Norb,Nmats,Nkpt],asympt_corr,tau_uniform)
   end interface Fmats2itau_vec

   interface Fitau2mats_mat
      module procedure Fitau2mats_mat_Gw                                        !(beta,Gitau[Norb,Norb,Ntau],Gmats[Norb,Norb,Nmats],tau_uniform)
      module procedure Fitau2mats_mat_Gwk                                       !(beta,Gitau[Norb,Norb,Ntau,Nkpt],Gmats[Norb,Norb,Nmats,Nkpt],tau_uniform)
   end interface Fitau2mats_mat
   interface Fitau2mats_vec
      module procedure Fitau2mats_vec_Gw                                        !(beta,Gitau[Norb,Ntau],Gmats[Norb,Nmats],tau_uniform)
      module procedure Fitau2mats_vec_Gwk                                       !(beta,Gitau[Norb,Ntau,Nkpt],Gmats[Norb,Nmats,Nkpt],tau_uniform)
   end interface Fitau2mats_vec

   interface Bmats2itau
      module procedure Bmats2itau_Uw_component                                  !(beta,Umats[Nmats],Uitau[Ntau],asympt_corr,tau_uniform)
      module procedure Bmats2itau_Uw                                            !(beta,Umats[Nbp,Nbp,Nmats],Uitau[Nbp,Nbp,Ntau],asympt_corr,tau_uniform)
      module procedure Bmats2itau_Uwk                                           !(beta,Umats[Nbp,Nbp,Nmats,Nkpt],Uitau[Nbp,Nbp,Ntau,Nkpt],asympt_corr,tau_uniform)
   end interface Bmats2itau
   interface Bitau2mats
      module procedure Bitau2mats_Uw_component                                  !(beta,Uitau[Ntau],Umats[Nmats],asympt_corr,tau_uniform)
      module procedure Bitau2mats_Uw                                            !(beta,Uitau[Nbp,Nbp,Ntau],Umats[Nbp,Nbp,Nmats],tau_uniform)
      module procedure Bitau2mats_Uwk                                           !(beta,Uitau[Nbp,Nbp,Ntau,Nkpt],Umats[Nbp,Nbp,Nmats,Nkpt],tau_uniform)
   end interface Bitau2mats

   interface linspace
      module procedure linspace_i
      module procedure linspace_d
   end interface linspace

   interface trapezoid_integration
      module procedure trapezoid_integration_d
      module procedure trapezoid_integration_z
   end interface trapezoid_integration

   interface cubic_interp
      module procedure cubic_interp_single
      module procedure cubic_interp_multiple
   end interface cubic_interp

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   real(8),parameter,private                :: pi=3.14159265358979323846d0
   complex(8),parameter,private             :: czero=dcmplx(0.d0,0.d0)
   real(8),private,protected                :: de_default=1d0
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
   public :: Fmats2itau_mat
   public :: Fmats2itau_vec
   public :: Fmats2itau
   public :: Fitau2mats_mat
   public :: Fitau2mats_vec
   public :: Fitau2mats
   public :: Bmats2itau
   public :: Bitau2mats
   public :: FermionicFilon
   public :: BosonicFilon
   public :: cubic_interp
   public :: inquireFile
   public :: inquireDir
   public :: createDir
   public :: removeDir
   public :: removeFile
   public :: set_powtau
   public :: tick
   public :: tock
   !functions
   public :: FermionicFreqMesh
   public :: BosonicFreqMesh
   public :: denspace
   public :: linspace
   public :: trapezoid_integration
   public :: free_unit
   public :: reg

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Creates fermionic weights (cos,sin) and add asymptotic correction
   !---------------------------------------------------------------------------!
   subroutine set_powtau(de_external)
      implicit none
      real(8),intent(in) :: de_external
      de_default = de_external
   end subroutine set_powtau


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
   !---------------------------------------------------------------------------!
   subroutine mats2itau_FermionicCoeff(tau,coswt,sinwt,correct)
      !
      ! use utils_misc
      use linalg, only : inv
      implicit none
      !
      real(8),intent(in)                    :: tau(:)
      real(8),intent(inout)                 :: coswt(:,:)
      real(8),intent(inout)                 :: sinwt(:,:)
      logical,intent(in)                    :: correct
      !
      real(8)                               :: beta
      real(8)                               :: v1,v2,v3
      real(8)                               :: sin1,cos2,sin3,cos4
      real(8)                               :: tail1,tail2,tail3,tail4
      real(8),allocatable                   :: tail(:,:),Ae(:,:),Ao(:,:)
      real(8),allocatable                   :: wmats(:)
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau
      logical                               :: abort
      !
      !
      if(verbose)write(*,"(A)") "---- mats2itau_FermionicCoeff"
      !
      !
      Ntau = size(tau)
      Nmats = size(coswt,dim=1)
      !call assert_shape(coswt,[Nmats,Ntau],"mats2itau_FermionicCoeff","coswt")
      !call assert_shape(sinwt,[Nmats,Ntau],"mats2itau_FermionicCoeff","sinwt")
      beta = tau(Ntau)
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(beta,Nmats)
      v1 = wmats(Nmats-1)
      v2 = wmats(Nmats-2)
      v3 = wmats(Nmats-3)
      !
      ! even tail
      allocate(tail(2,2));tail=0d0
      tail(1,1) = 1d0/(v1**2)
      tail(1,2) = 1d0/(v1**4)
      tail(2,1) = 1d0/(v2**2)
      tail(2,2) = 1d0/(v2**4)
      call inv(tail)
      Ae = tail
      deallocate(tail)
      ! odd tail
      allocate(tail(2,2));tail=0d0
      tail(1,1) = 1d0/v1
      tail(1,2) = 1d0/(v1**3)
      tail(2,1) = 1d0/v2
      tail(2,2) = 1d0/(v2**3)
      call inv(tail)
      Ao = tail
      deallocate(tail)
      !
      coswt=0d0
      sinwt=0d0
      do itau=1,Ntau
         !
         sin1=0d0
         cos2=0d0
         sin3=0d0
         cos4=0d0
         !
         ! cosm(i) = (1/beta) S[n] cos(vn*tau) / vn^i
         ! sinm(i) = (1/beta) S[n] sin(vn*tau) / vn^i
         do iw=1,Nmats
            !
            coswt(iw,itau) = dcos(wmats(iw) * tau(itau))/beta
            sinwt(iw,itau) = dsin(wmats(iw) * tau(itau))/beta
            !
            sin1 = sin1 + sinwt(iw,itau) / wmats(iw)
            cos2 = cos2 + coswt(iw,itau) / (wmats(iw)**2)
            sin3 = sin3 + sinwt(iw,itau) / (wmats(iw)**3)
            cos4 = cos4 + coswt(iw,itau) / (wmats(iw)**4)
            !
         enddo
         !
         if(correct)then
            !
            ! (1/beta) S[n=0,inf] cos(vn*tau) / vn^2  = beta/8 - tau/4
            ! (1/beta) S[n=0,inf] cos(vn*tau) / vn^4  = ( beta^3/96 - beta*tau^2/16 + tau^3/24 )
            !
            tail2 = Ae(1,1) * ( beta/8d0 - tau(itau)/4d0 - cos2 )
            tail4 = Ae(2,1) * ( beta**3/96d0 - beta*tau(itau)**2/16d0 + tau(itau)**3/24d0 - cos4 )
            coswt(Nmats,itau) = coswt(Nmats,itau) + tail2 + tail4
            !
            tail2 = Ae(1,2) * ( beta/8d0 - tau(itau)/4d0 - cos2 )
            tail4 = Ae(2,2) * ( beta**3/96d0 - beta*tau(itau)**2/16d0 + tau(itau)**3/24d0 - cos4 )
            coswt(Nmats-1,itau) = coswt(Nmats-1,itau) + tail2 + tail4
            !
            ! (1/beta) S[n=0,inf] sin(vn*tau) / vn    = 1/4
            ! (1/beta) S[n=0,inf] sin(vn*tau) / vn^3  = ( beta*tau - tau^2 ) / 8
            !
            tail1 = Ao(1,1) * ( 1d0/4d0 - sin1 )
            tail3 = Ao(2,1) * ( (beta*tau(itau) - tau(itau)**2)/8d0 - sin3 )
            sinwt(Nmats,itau)   = sinwt(Nmats,itau) + tail1 + tail3
            !
            tail1 = Ao(1,2) * ( 1d0/4d0 - sin1 )
            tail3 = Ao(2,2) * ( (beta*tau(itau) - tau(itau)**2)/8d0 - sin3 )
            sinwt(Nmats-1,itau) = sinwt(Nmats-1,itau) + tail1 + tail3
            !
         endif
      enddo !itau
      deallocate(wmats,Ae,Ao)
      !
      abort=.false.
      do itau=1,Ntau
         do iw=1,Nmats
            if(coswt(iw,itau).ne.coswt(iw,itau))then
               !NaN condition
               write(*,"(A)")"mats2itau_FermionicCoeff: coswt is NaN."
               abort=.true.
            elseif(abs(coswt(iw,itau)).ge.huge(1d0))then
               !Infinity condition
               write(*,"(A)")"mats2itau_FermionicCoeff: coswt is Infinity."
               abort=.true.
            endif
            if(sinwt(iw,itau).ne.sinwt(iw,itau))then
               !NaN condition
               write(*,"(A)")"mats2itau_FermionicCoeff: sinwt is NaN."
               abort=.true.
            elseif(abs(sinwt(iw,itau)).ge.huge(1d0))then
               !Infinity condition
               write(*,"(A)")"mats2itau_FermionicCoeff: sinwt is Infinity."
               abort=.true.
            endif
         enddo
      enddo
      if(abort) stop "mats2itau_FermionicCoeff: coefficient error."
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
   !---------------------------------------------------------------------------!
   subroutine mats2itau_BosonicCoeff(tau,coswt,correct)
      !
      ! use utils_misc
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
      logical                               :: abort
      !
      !
      if(verbose)write(*,"(A)") "---- mats2itau_BosonicCoeff"
      !
      !
      Ntau = size(tau)
      Nmats = size(coswt,dim=1)
      !call assert_shape(coswt,[Nmats,Ntau],"mats2itau_BosonicCoeff","coswt")
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
         if(correct) then
            !
            coswt(Nmats,itau) = coswt(Nmats,itau) + wmats(Nmats)**2 * wb1
            !
         endif
         !
      enddo
      !
      abort=.false.
      do itau=1,Ntau
         do iw=1,Nmats
            if(coswt(iw,itau).ne.coswt(iw,itau))then
               !NaN condition
               write(*,"(A)")"mats2itau_BosonicCoeff: coswt is NaN."
               abort=.true.
            elseif(abs(coswt(iw,itau)).ge.huge(1d0))then
               !Infinity condition
               write(*,"(A)")"mats2itau_BosonicCoeff: coswt is Infinity."
               abort=.true.
            endif
         enddo
      enddo
      if(abort) stop "mats2itau_BosonicCoeff: coefficient error."
      !
   end subroutine mats2itau_BosonicCoeff


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from mats to tau of a matrix Gf
   !---------------------------------------------------------------------------!
   subroutine Fmats2itau_mat_Gw(beta,Gmats,Gitau,asympt_corr,tau_uniform,atBeta)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gmats(:,:,:)
      complex(8),intent(inout)              :: Gitau(:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      logical,intent(in),optional           :: atBeta
      !
      real(8),allocatable                   :: coswt(:,:),sinwt(:,:)
      real(8),allocatable                   :: tau(:)
      complex(8),allocatable                :: Ge(:,:),Go(:,:)
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau,Norb
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      logical                               :: atBeta_
      !
      !
      if(verbose)write(*,"(A)") "---- Fmats2itau_mat_Gw"
      !
      !
      Norb = size(Gmats,dim=1)
      Nmats = size(Gmats,dim=3)
      Ntau = size(Gitau,dim=3)
      if(size(Gmats,dim=1).ne.size(Gmats,dim=2)) stop "Fmats2itau_mat_Gw: Gmats not square."
      !call assert_shape(Gitau,[Norb,Norb,Ntau],"Fmats2itau_mat_Gw","Gitau")
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      atBeta_ = .false.
      if(present(atBeta)) atBeta_ = atBeta
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
      deallocate(tau)
      !
      Gitau=czero
      allocate(Ge(Norb,Norb));Ge=czero
      allocate(Go(Norb,Norb));Go=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(atBeta_,Ntau,Nmats,Gmats,coswt,sinwt,Gitau),&
      !$OMP PRIVATE(itau,iw,Ge,Go)
      !$OMP DO
      do itau=1,Ntau
         if(atBeta_.and.(itau.ne.Ntau))cycle
         do iw=1,Nmats
            !
            !Gab(-iw) = Gba*(iwn)
            Ge = Gmats(:,:,iw) + transpose(conjg(Gmats(:,:,iw)))
            Go = Gmats(:,:,iw) - transpose(conjg(Gmats(:,:,iw)))
            !
            Gitau(:,:,itau) = Gitau(:,:,itau) + coswt(iw,itau)*Ge -dcmplx(0d0,1d0)*sinwt(iw,itau)*Go
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,sinwt,Ge,Go)
      !
   end subroutine Fmats2itau_mat_Gw


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from mats to tau of a vector Gf
   !---------------------------------------------------------------------------!
   subroutine Fmats2itau_vec_Gw(beta,Gmats,Gitau,asympt_corr,tau_uniform,atBeta)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gmats(:,:)
      complex(8),intent(inout)              :: Gitau(:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      logical,intent(in),optional           :: atBeta
      !
      real(8),allocatable                   :: coswt(:,:),sinwt(:,:)
      real(8),allocatable                   :: tau(:)
      complex(8),allocatable                :: Ge(:),Go(:)
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau,Norb
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      logical                               :: atBeta_
      !
      !
      if(verbose)write(*,"(A)") "---- Fmats2itau_vec_Gw"
      !
      !
      Norb = size(Gmats,dim=1)
      Nmats = size(Gmats,dim=2)
      Ntau = size(Gitau,dim=2)
      !call assert_shape(Gitau,[Norb,Nt:au],"Fmats2itau_vec_Gw","Gitau")
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      atBeta_ = .false.
      if(present(atBeta)) atBeta_ = atBeta
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
      deallocate(tau)
      !
      Gitau=czero
      allocate(Ge(Norb));Ge=czero
      allocate(Go(Norb));Go=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(atBeta_,Ntau,Nmats,Gmats,coswt,sinwt,Gitau),&
      !$OMP PRIVATE(itau,iw,Ge,Go)
      !$OMP DO
      do itau=1,Ntau
         if(atBeta_.and.(itau.ne.Ntau))cycle
         do iw=1,Nmats
            !
            !Gaa(-iw) = Gaa*(iwn)
            Ge = Gmats(:,iw) + conjg(Gmats(:,iw))
            Go = Gmats(:,iw) - conjg(Gmats(:,iw))
            !
            Gitau(:,itau) = Gitau(:,itau) + coswt(iw,itau)*Ge -dcmplx(0d0,1d0)*sinwt(iw,itau)*Go
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,sinwt,Ge,Go)
      !
   end subroutine Fmats2itau_vec_Gw


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from mats to tau of a Gf
   !---------------------------------------------------------------------------!
   subroutine Fmats2itau(beta,Gmats,Gitau,asympt_corr,tau_uniform,atBeta)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gmats(:)
      complex(8),intent(inout)              :: Gitau(:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      logical,intent(in),optional           :: atBeta
      !
      real(8),allocatable                   :: coswt(:,:),sinwt(:,:)
      real(8),allocatable                   :: tau(:)
      complex(8)                            :: Ge,Go
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      logical                               :: atBeta_
      !
      !
      if(verbose)write(*,"(A)") "---- Fmats2itau"
      !
      !
      Nmats = size(Gmats)
      Ntau = size(Gitau)
      !call assert_shape(Gitau,[Ntau],"Fmats2itau","Gitau")
      !
      asympt_corr_ = .true.
      if(present(asympt_corr)) asympt_corr_ = asympt_corr
      tau_uniform_ = .false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      atBeta_ = .false.
      if(present(atBeta)) atBeta_ = atBeta
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
      deallocate(tau)
      !
      Gitau=czero
      Ge=czero
      Go=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(atBeta_,Ntau,Nmats,Gmats,coswt,sinwt,Gitau),&
      !$OMP PRIVATE(itau,iw,Ge,Go)
      !$OMP DO
      do itau=1,Ntau
         if(atBeta_.and.(itau.ne.Ntau))cycle
         do iw=1,Nmats
            !
            ! Gab(iw) = Gba*(-iwn) --> Gab(-iw) = Gba*(iwn)
            Ge = Gmats(iw) + conjg(Gmats(iw))
            Go = Gmats(iw) - conjg(Gmats(iw))
            !
            Gitau(itau) = Gitau(itau) + coswt(iw,itau)*Ge -dcmplx(0d0,1d0)*sinwt(iw,itau)*Go
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt,sinwt)
      !
   end subroutine Fmats2itau


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from tau to mats of a matrix Gf
   !---------------------------------------------------------------------------!
   subroutine Fitau2mats_mat_Gw(beta,Gitau,Gmats,tau_uniform)
      !
      ! use utils_misc
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
      if(verbose)write(*,"(A)") "---- Fitau2mats_mat_Gw"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=3)
      Nmats = size(Gmats,dim=3)
      if(size(Gitau,dim=1).ne.size(Gitau,dim=2)) stop "Fitau2mats_mat_Gw: Gitau not square."
      !call assert_shape(Gmats,[Norb,Norb,Nmats],"Fitau2mats_mat_Gw","Gmats")
      if(mod(Ntau,2).eq.0) stop "Fitau2mats_mat_Gw: Required Filon routines are not working with odd segments."
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
         do iwan2=1,Norb
            do iwan1=1,Norb
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
      ! use utils_misc
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
      if(verbose)write(*,"(A)") "---- Fitau2mats_mat_Gwk"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=3)
      Nkpt = size(Gitau,dim=4)
      Nmats = size(Gmats,dim=3)
      if(size(Gitau,dim=1).ne.size(Gitau,dim=2)) stop "Fitau2mats_mat_Gwk: Gitau not square."
      !call assert_shape(Gmats,[Norb,Norb,Nmats,Nkpt],"Fitau2mats_mat_Gwk","Gmats")
      if(mod(Ntau,2).eq.0) stop "Fitau2mats_mat_Gwk: Required Filon routines are not working with odd segments."
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
            do iwan2=1,Norb
               do iwan1=1,Norb
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
      ! use utils_misc
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
      if(verbose)write(*,"(A)") "---- Fitau2mats_vec_Gw"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=2)
      Nmats = size(Gmats,dim=2)
      !call assert_shape(Gmats,[Norb,Nmats],"Fitau2mats_vec_Gw","Gmats")
      if(mod(Ntau,2).eq.0) stop "Fitau2mats_vec_Gw: Required Filon routines are not working with odd segments."
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
      ! use utils_misc
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
      if(verbose)write(*,"(A)") "---- Fitau2mats_vec_Gwk"
      !
      !
      Norb = size(Gitau,dim=1)
      Ntau = size(Gitau,dim=2)
      Nkpt = size(Gitau,dim=3)
      Nmats = size(Gmats,dim=2)
      !call assert_shape(Gmats,[Norb,Nmats,Nkpt],"Fitau2mats_vec_Gwk","Gmats")
      if(mod(Ntau,2).eq.0) stop "Fitau2mats_vec_Gwk: Required Filon routines are not working with odd segments."
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
   !PURPOSE: Perform the Fourier transform from tau to mats of a Gf
   !---------------------------------------------------------------------------!
   subroutine Fitau2mats(beta,Gitau,Gmats,tau_uniform)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Gitau(:)
      complex(8),intent(inout)              :: Gmats(:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: RealG(:),ImagG(:)
      real(8)                               :: rwcos,rwsin,cwcos,cwsin
      integer                               :: iw
      integer                               :: Nmats,Ntau
      logical                               :: tau_uniform_
      !
      !
      if(verbose)write(*,"(A)") "---- Fitau2mats"
      !
      !
      Ntau = size(Gitau)
      Nmats = size(Gmats)
      !call assert_shape(Gmats,[Nmats],"Fitau2mats","Gmats")
      if(mod(Ntau,2).eq.0) stop "Fitau2mats: Required Filon routines are not working with odd segments."
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
      !$OMP SHARED(Nmats,wmats,tau,Gitau,Gmats),&
      !$OMP PRIVATE(iw,RealG,ImagG,rwcos,rwsin,cwcos,cwsin)
      !$OMP DO
      do iw=1,Nmats
         !
         RealG = dreal( Gitau )
         ImagG = dimag( Gitau )
         !
         call FermionicFilon(wmats(iw),tau,RealG,rwcos,rwsin)
         call FermionicFilon(wmats(iw),tau,ImagG,cwcos,cwsin)
         !
         Gmats(iw) = dcmplx(rwcos-cwsin,0d0) + dcmplx(0d0,rwsin+cwcos)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(RealG,ImagG,tau,wmats)
      !
   end subroutine Fitau2mats


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from mats to tau of a bosonic tensor
   !---------------------------------------------------------------------------!
   subroutine Bmats2itau_Uw_component(beta,Umats,Uitau,asympt_corr,tau_uniform,Umats_bare)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Umats(:)
      complex(8),intent(inout)              :: Uitau(:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      complex(8),intent(in),optional        :: Umats_bare
      !
      real(8),allocatable                   :: coswt(:,:)
      real(8),allocatable                   :: tau(:)
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      !
      !
      if(verbose)write(*,"(A)") "---- Bmats2itau_Uw_component. Warning: This is only for well behaved, i.e. U(iw)=U(-iw), components."
      !
      !
      Nmats = size(Umats)
      Ntau = size(Uitau)
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
      deallocate(tau)
      !
      Uitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ntau,Nmats,Umats,coswt,Uitau,Umats_bare),&
      !$OMP PRIVATE(itau,iw)
      !$OMP DO
      do itau=1,Ntau
         do iw=1,Nmats
            !
            if(present(Umats_bare))then
               Uitau(itau) = Uitau(itau) + coswt(iw,itau) * (Umats(iw)-Umats_bare)
            else
               Uitau(itau) = Uitau(itau) + coswt(iw,itau) * Umats(iw)
            endif
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt)
      !
   end subroutine Bmats2itau_Uw_component
   !
   subroutine Bmats2itau_Uw(beta,Umats,Uitau,asympt_corr,tau_uniform,Umats_bare)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Umats(:,:,:)
      complex(8),intent(inout)              :: Uitau(:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      complex(8),intent(in),optional        :: Umats_bare(:,:)
      !
      real(8),allocatable                   :: coswt(:,:)
      real(8),allocatable                   :: tau(:)
      integer                               :: iw,itau,ib1,ib2
      integer                               :: Nmats,Ntau,Nbp
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      !
      !
      if(verbose)write(*,"(A)") "---- Bmats2itau_Uw"
      !
      !
      Nbp = size(Umats,dim=1)
      Nmats = size(Umats,dim=3)
      Ntau = size(Uitau,dim=3)
      if(size(Umats,dim=1).ne.size(Umats,dim=2)) stop "Bmats2itau_Uw: Umats not square."
      !call assert_shape(Uitau,[Nbp,Nbp,Ntau],"Bmats2itau_Uw","Uitau")
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
      deallocate(tau)
      !
      Uitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ntau,Nmats,Nbp,Umats,coswt,Uitau,Umats_bare),&
      !$OMP PRIVATE(itau,iw,ib1,ib2)
      !$OMP DO
      do itau=1,Ntau
         !
         if(present(Umats_bare))then
            !
            do ib2=1,Nbp
               do ib1=1,Nbp
                  do iw=1,Nmats
                     Uitau(ib1,ib2,itau) = Uitau(ib1,ib2,itau) + coswt(iw,itau) * (Umats(ib1,ib2,iw)-Umats_bare(ib1,ib2))
                  enddo
               enddo
            enddo
            !
         else
            !
            do ib2=1,Nbp
               do ib1=1,Nbp
                  do iw=1,Nmats
                     Uitau(ib1,ib2,itau) = Uitau(ib1,ib2,itau) + coswt(iw,itau) * Umats(ib1,ib2,iw)
                  enddo
               enddo
            enddo
            !
         endif
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt)
      !
   end subroutine Bmats2itau_Uw
   !
   subroutine Bmats2itau_Uwk(beta,Umats,Uitau,asympt_corr,tau_uniform,Umats_bare)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Umats(:,:,:,:)
      complex(8),intent(inout)              :: Uitau(:,:,:,:)
      logical,intent(in),optional           :: asympt_corr
      logical,intent(in),optional           :: tau_uniform
      complex(8),intent(in),optional        :: Umats_bare(:,:,:)
      !
      real(8),allocatable                   :: coswt(:,:)
      real(8),allocatable                   :: tau(:)
      integer                               :: iw,itau,ib1,ib2,iq
      integer                               :: Nmats,Ntau,Nbp,Nkpt
      logical                               :: asympt_corr_
      logical                               :: tau_uniform_
      !
      !
      if(verbose)write(*,"(A)") "---- Bmats2itau_Uwk"
      !
      !
      Nbp = size(Umats,dim=1)
      Nmats = size(Umats,dim=3)
      Nkpt = size(Umats,dim=4)
      Ntau = size(Uitau,dim=3)
      if(size(Umats,dim=1).ne.size(Umats,dim=2)) stop "Bmats2itau_Uwk: Umats not square."
      !call assert_shape(Uitau,[Nbp,Nbp,Ntau,Nkpt],"Bmats2itau_Uwk","Uitau")
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
      deallocate(tau)
      !
      Uitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Ntau,Nmats,Nbp,Umats,coswt,Uitau,Umats_bare),&
      !$OMP PRIVATE(iq,itau,iw,ib1,ib2)
      !$OMP DO
      do iq=1,Nkpt
         !
         if(present(Umats_bare))then
            !
            do itau=1,Ntau
               do iw=1,Nmats
                  do ib2=1,Nbp
                     do ib1=1,Nbp
                        Uitau(ib1,ib2,itau,iq) = Uitau(ib1,ib2,itau,iq) + coswt(iw,itau) * (Umats(ib1,ib2,iw,iq)-Umats_bare(ib1,ib2,iq))
                     enddo
                  enddo
               enddo
            enddo
            !
         else
            !
            do itau=1,Ntau
               do iw=1,Nmats
                  do ib2=1,Nbp
                     do ib1=1,Nbp
                        Uitau(ib1,ib2,itau,iq) = Uitau(ib1,ib2,itau,iq) + coswt(iw,itau) * Umats(ib1,ib2,iw,iq)
                     enddo
                  enddo
               enddo
            enddo
            !
         endif
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coswt)
      !
   end subroutine Bmats2itau_Uwk


   !---------------------------------------------------------------------------!
   !PURPOSE: Perform the Fourier transform from tau to mats of a bosonic tensor
   !---------------------------------------------------------------------------!
   subroutine Bitau2mats_Uw_component(beta,Uitau,Umats,tau_uniform)
      !
      ! use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: beta
      complex(8),intent(in)                 :: Uitau(:)
      complex(8),intent(inout)              :: Umats(:)
      logical,intent(in),optional           :: tau_uniform
      !
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: wcos(:),wsin(:)
      real(8)                               :: RealU,ImagU
      integer                               :: iw,itau
      integer                               :: Nmats,Ntau
      logical                               :: tau_uniform_
      !
      !
      if(verbose)write(*,"(A)") "---- Bitau2mats_Uw_component"
      !
      !
      Ntau = size(Uitau)
      Nmats = size(Umats)
      if(mod(Ntau,2).eq.0) stop "Bitau2mats_Uw_component: Required Filon routines are not working with odd segments."
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
      !$OMP SHARED(Nmats,Ntau,wmats,tau,Uitau,Umats),&
      !$OMP PRIVATE(iw,itau,RealU,ImagU,wcos,wsin)
      !$OMP DO
      do iw=1,Nmats
         !
         call BosonicFilon(wmats(iw),tau,wcos,wsin)
         !
         do itau=1,Ntau
            !
            RealU = dreal( Uitau(itau) )
            ImagU = dimag( Uitau(itau) )
            !
            Umats(iw) =   Umats(iw)                                    &
                      + ( RealU * wcos(itau) - ImagU * wsin(itau) )    &
                      + ( ImagU * wcos(itau) + RealU * wsin(itau) )*dcmplx(0d0,1d0)
            !
         enddo
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(wcos,wsin,tau,wmats)
      !
   end subroutine Bitau2mats_Uw_component
   !
   subroutine Bitau2mats_Uw(beta,Uitau,Umats,tau_uniform)
      !
      ! use utils_misc
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
      if(verbose)write(*,"(A)") "---- Bitau2mats_Uw"
      !
      !
      Nbp = size(Uitau,dim=1)
      Ntau = size(Uitau,dim=3)
      Nmats = size(Umats,dim=3)
      if(size(Uitau,dim=1).ne.size(Uitau,dim=2)) stop "Bitau2mats_Uw: Uitau not square."
      !call assert_shape(Umats,[Nbp,Nbp,Nmats],"Bitau2mats_Uw","Umats")
      if(mod(Ntau,2).eq.0) stop "Bitau2mats_Uw: Required Filon routines are not working with odd segments."
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
         do itau=1,Ntau
            do ib2=1,Nbp
               do ib1=1,Nbp
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
      ! use utils_misc
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
      integer                               :: iw,itau,ib1,ib2,iq
      integer                               :: Nmats,Ntau,Nbp,Nkpt
      logical                               :: tau_uniform_
      !
      !
      if(verbose)write(*,"(A)") "---- Bitau2mats_Uwk"
      !
      !
      Nbp = size(Uitau,dim=1)
      Ntau = size(Uitau,dim=3)
      Nkpt = size(Uitau,dim=4)
      Nmats = size(Umats,dim=3)
      if(size(Uitau,dim=1).ne.size(Uitau,dim=2)) stop "Bitau2mats_Uwk: Uitau not square."
      !call assert_shape(Umats,[Nbp,Nbp,Nmats,Nkpt],"Bitau2mats_Uwk","Umats")
      if(mod(Ntau,2).eq.0) stop "Bitau2mats_Uwk: Required Filon routines are not working with odd segments."
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
      !$OMP PRIVATE(iw,iq,ib1,ib2,itau,RealU,ImagU,wcos,wsin)
      !$OMP DO
      do iw=1,Nmats
         !
         call BosonicFilon(wmats(iw),tau,wcos,wsin)
         !
         do iq=1,Nkpt
            do itau=1,Ntau
               do ib2=1,Nbp
                  do ib1=1,Nbp
                     !
                     RealU = dreal( Uitau(ib1,ib2,itau,iq) )
                     ImagU = dimag( Uitau(ib1,ib2,itau,iq) )
                     !
                     Umats(ib1,ib2,iw,iq) =   Umats(ib1,ib2,iw,iq)                           &
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
      !call assert_shape(wcos,[npoints],"BosonicFilon","wcos")
      !call assert_shape(wsin,[npoints],"BosonicFilon","wsin")
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
            write(*,"(A)")"BosonicFilon: wcos is NaN."
            abort=.true.
         elseif(abs(wcos(n)).ge.huge(1d0))then
            !Infinity condition
            write(*,"(A)")"BosonicFilon: wcos is Infinity."
            abort=.true.
         endif
         if(wsin(n).ne.wsin(n))then
            !NaN condition
            write(*,"(A)")"BosonicFilon: wsin is NaN."
            abort=.true.
         elseif(abs(wsin(n)).ge.huge(1d0))then
            !Infinity condition
            write(*,"(A)")"BosonicFilon: wsin is Infinity."
            abort=.true.
         endif
      enddo
      if(abort) stop "BosonicFilon: coefficient error. Increase precision and recompile."
      !
   end subroutine BosonicFilon


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
      !
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
   function denspace(end,num,center,expfact) result(array)
      implicit none
      real(8),intent(in)                    :: end
      integer,intent(in)                    :: num
      logical,intent(in),optional           :: center
      real(8),intent(in),optional           :: expfact
      real(8)                               :: array(num)
      !
      integer                               :: nsimp=2
      integer                               :: nseg
      integer                               :: i,n,n2,nseg2
      real(8)                               :: mesh,de,a,b
      logical                               :: center_
      !data de /1.0d0/
      de = de_default ! 1.0d0
      if(present(expfact)) de = expfact
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
   !PURPOSE: function which computes the integral of an array
   !---------------------------------------------------------------------------!
   function trapezoid_integration_d(fx,x) result(Int)
      implicit none
      real(8),intent(in)                    :: fx(:)
      real(8),intent(in)                    :: x(:)
      real(8)                               :: Int,dx
      integer                               :: i
      !
      if(size(fx).gt.1) then
         !
         Int=0d0
         do i=2,size(fx)
            dx = x(i)-x(i-1)
            Int = Int + ( fx(i) + fx(i-1) ) * (dx/2d0)
         enddo
         !
      else
         !
         Int = fx(1)
         !
      endif
      !
   end function trapezoid_integration_d
   !
   function trapezoid_integration_z(fx,x) result(Int)
      implicit none
      complex(8),intent(in)                 :: fx(:)
      real(8),intent(in)                    :: x(:)
      complex(8)                            :: Int
      real(8)                               :: dx
      integer                               :: i
      !
      if(size(fx).gt.1) then
         !
         Int=0d0
         do i=2,size(fx)
            dx = x(i)-x(i-1)
            Int = Int + ( fx(i) + fx(i-1) ) * (dx/2d0)
         enddo
         !
      else
         !
         Int = fx(1)
         !
      endif
      !
   end function trapezoid_integration_z


   !---------------------------------------------------------------------------!
   !PURPOSE: This is a wrapper that encloses nspline and splint
   !---------------------------------------------------------------------------!
   function cubic_interp_single(x,y,xp) result(yp)
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
   end function cubic_interp_single
   !
   function cubic_interp_multiple(x,y,xp) result(yp)
      implicit none
      real(8),intent(in)                    :: x(:)
      real(8),intent(in)                    :: y(:)
      real(8),intent(in)                    :: xp(:)
      real(8)                               :: yp(size(xp))
      !
      integer                               :: n,ip
      real(8),allocatable                   :: y2(:)
      !
      n=size(x)
      if(size(y).ne.n) stop "nspline: size(y).ne.size(x)."
      allocate(y2(n));y2=0d0
      call nspline(x,y,y2)
      do ip=1,size(xp)
         call splint(x,y,y2,xp(ip),yp(ip))
      enddo
      deallocate(y2)
      !
   end function cubic_interp_multiple


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
         call system(reg(mkdirCmd))
         !call execute_command_line(reg(mkdirCmd))
      endif
      !
   end subroutine createDir
   !
   subroutine removeDir(dirpath,verb,list)
      implicit none
      character(len=*),intent(in)           :: dirpath
      logical,intent(in),optional           :: verb
      logical,intent(in),optional           :: list
      character(len=256)                    :: rmdirCmd
      logical                               :: direxists
      logical                               :: verbose_,list_
      !
      verbose_=.true.
      if(present(verb))verbose_=verb
      list_=.false.
      if(present(list))list_=list
      !
      call inquireDir(reg(dirpath),direxists,hardstop=.false.,verb=verbose_)
      if(direxists)then
         rmdirCmd = "rm -r "//reg(dirpath)
         if(verbose_)write(*,"(A)") "     Removing directory: "//reg(dirpath)
         if(verbose_)write(*,"(A)") "     "//reg(rmdirCmd)
         call system(reg(rmdirCmd))
         !call execute_command_line(reg(rmdirCmd))
      else
         if(list_)then
            rmdirCmd = "rm -r "//reg(dirpath)
            if(verbose_)write(*,"(A)") "     Removing directory: "//reg(dirpath)
            if(verbose_)write(*,"(A)") "     "//reg(rmdirCmd)
            call system(reg(rmdirCmd))
            !call execute_command_line(reg(rmdirCmd))
         endif
      endif

      !
   end subroutine removeDir
   !
   subroutine removeFile(filepath,verb,list)
      implicit none
      character(len=*),intent(in)           :: filepath
      logical,intent(in),optional           :: verb
      logical,intent(in),optional           :: list
      character(len=256)                    :: rmfileCmd
      logical                               :: filexists
      logical                               :: verbose_,list_
      !
      verbose_=.true.
      if(present(verb))verbose_=verb
      list_=.false.
      if(present(list))list_=list
      !
      call inquireFile(reg(filepath),filexists,hardstop=.false.,verb=verbose_)
      if(filexists)then
         rmfileCmd = "rm -r "//reg(filepath)
         if(verbose_)write(*,"(A)") "     Removing file: "//reg(filepath)
         if(verbose_)write(*,"(A)") "     "//reg(rmfileCmd)
         call system(reg(rmfileCmd))
         !call execute_command_line(reg(rmfileCmd))
      else
         if(list_)then
            rmfileCmd = "rm -r "//reg(filepath)
            if(verbose_)write(*,"(A)") "     Removing file: "//reg(filepath)
            if(verbose_)write(*,"(A)") "     "//reg(rmfileCmd)
            call system(reg(rmfileCmd))
            !call execute_command_line(reg(rmfileCmd))
         endif
      endif
      !
   end subroutine removeFile


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


end module fourier_transforms
