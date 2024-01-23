module interactions

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface read_U_spex
      module procedure read_U_spex_full                                         ![BosonicField,LocalOnly,save2readable,pathOUTPUT(optional to change output path),doAC(optional to override AC)]
      module procedure read_U_spex_Uloc0                                        ![Matrix,pathOUTPUT(optional to change output path)]
   end interface read_U_spex

   interface read_U_respack
      module procedure read_U_respack_full                                      ![BosonicField,LocalOnly,save2readable,pathOUTPUT(optional to change output path),doAC(optional to override AC)]
      module procedure read_U_respack_Uloc0                                     ![Matrix,pathOUTPUT(optional to change output path)]
   end interface read_U_respack

   interface read_U_vasp
      module procedure read_U_vasp_full                                         ![BosonicField,LocalOnly,save2readable,pathOUTPUT(optional to change output path),doAC(optional to override AC)]
      module procedure read_U_vasp_Uloc0                                        ![Matrix,pathOUTPUT(optional to change output path)]
   end interface read_U_vasp

   interface build_Umat
      module procedure build_Umat_singlParam                 ! (GW Format)      ![Matrix,Uaa_screened,Uab_screened,J_screened]
      module procedure build_Umat_multiParam                 ! (GW Format)      ![Matrix,Vector,Matrix,Matrix]
   end interface build_Umat

   interface build_Uret
      module procedure build_Uret_singlParam_ph              !      (GW Format) ![BosonicField,Uaa_bare,Uab_bare,J_bare,vector_g,vector_w0,LocalOnly(optional)]
      module procedure build_Uret_multiParam_ph              !      (GW Format) ![BosonicField,Vector,Matrix,Matrix,vector_g,vector_w0,LocalOnly(optional)]      !NOT USED: the input is not formatted for interactions with different matrix elements
      module procedure build_Uret_singlParam_Vn              !      (GW Format) ![BosonicField,Uaa_bare,Uab_bare,J_bare,vector_g,vector_w0,LocalOnly(optional)]
      module procedure build_Uret_multiParam_Vn              !      (GW Format) ![BosonicField,Vector,Matrix,Matrix,vector_g,vector_w0,LocalOnly(optional)]      !NOT USED: the input is not formatted for interactions with different matrix elements
   end interface build_Uret

   interface calc_QMCinteractions
      module procedure calc_QMCinteractions_static
      module procedure calc_QMCinteractions_retarded
   end interface calc_QMCinteractions

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: calc_W_full
   public :: calc_W_edmft
   public :: calc_chi
   public :: read_U_spex
   public :: read_U_respack
   public :: read_U_vasp
   public :: read_U
   public :: build_Umat
   public :: build_Uret
   public :: calc_QMCinteractions
   public :: calc_curlyU
   public :: calc_Wimp

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Lattice inversion to get fully screened interaction - GW+EDMFT
   !---------------------------------------------------------------------------!
   subroutine calc_W_full(Wmats,Umats,Pmats,Emats,Lttc,symQ)
      !
      use parameters
      use linalg, only : zeye, inv
      use utils_misc
      use utils_fields
      use file_io
      use crystal
      use input_vars, only : HandleGammaPoint, alphaGamma
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wmats
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(BosonicField),intent(inout)      :: Emats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: symQ
      !
      complex(8),allocatable                :: invW(:,:)
      complex(8),allocatable                :: epsGamma(:,:,:)
      integer,allocatable                   :: AverageList(:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw,iwU
      integer                               :: iavg,Navg
      logical                               :: Ustatic,symQ_,smear,add_iq
      !
      !
      if(verbose)write(*,"(A)") "---- calc_W_full"
      !
      !
      ! Check on the input Fields
      if(.not.Wmats%status) stop "calc_W_full: Wmats not properly initialized."
      if(.not.Umats%status) stop "calc_W_full: Umats not properly initialized."
      if(.not.Pmats%status) stop "calc_W_full: Pmats not properly initialized."
      if(Wmats%Nkpt.eq.0) stop "calc_W_full: Wmats k dependent attributes not properly initialized."
      if(Umats%Nkpt.eq.0) stop "calc_W_full: Umats k dependent attributes not properly initialized."
      if(Pmats%Nkpt.eq.0) stop "calc_W_full: Pmats k dependent attributes not properly initialized."
      if(Umats%iq_gamma.lt.0) stop "calc_W_full: Umats iq_gamma not defined."
      !
      symQ_=.false.
      if(present(symQ))symQ_=symQ
      Ustatic=.false.
      if(Umats%Npoints.eq.1) Ustatic=.true.
      if(Ustatic)write(*,"(A)")"     Static U bare."
      smear = HandleGammaPoint.gt.0
      !
      Nbp = Wmats%Nbp
      Nkpt = Wmats%Nkpt
      Beta = Wmats%Beta
      Nmats = Wmats%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "calc_W_full: Either Umats and/or Pmats have different orbital dimension with respect to Wmats."
      if(all([Umats%Nkpt-Nkpt,Pmats%Nkpt-Nkpt].ne.[0,0])) stop "calc_W_full: Either Umats and/or Pmats have different number of k-points with respect to Wmats."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "calc_W_full: Either Umats and/or Pmats have different Beta with respect to Wmats."
      if(Pmats%Npoints.ne.Nmats) stop "calc_W_full: Pmats has different number of Matsubara points with respect to Wmats."
      if(Emats%status) call assert_shape(Emats%screened,[Nbp,Nbp,Nmats,Nkpt],"calc_W_full","Emats%screened")
      !
      allocate(invW(Nbp,Nbp));invW=czero
      call clear_attributes(Wmats)
      !
      if(smear)then
         if(.not.small_ik_stored)call fill_smallk(Lttc%kpt)
         allocate(epsGamma(Nbp,Nbp,Nmats));epsGamma=czero
         AverageList = pack(small_ik(:,1),(small_ik(:,2).le.HandleGammaPoint))
         Navg = size(AverageList)
         if(verbose) write(*,"(A,100I4)")"     Averaging dielectric function over "//str(Navg)//" K-points: ",AverageList
      endif
      !
      ! Assuming that the Polarization vanishes at iw-->inf
      Wmats%bare = Umats%bare
      !
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iw,iwU,iq,invW,iavg,add_iq)
      !$OMP DO
      do iw=1,Wmats%Npoints
         !
         iwU = iw
         if(Ustatic)iwU = 1
         !
         iavg=0
         do iq=1,Wmats%Nkpt
            !
            !avoid the gamma point
            if((iq.eq.Umats%iq_gamma).and.smear)cycle
            !
            ! [ 1 - Pi*U ]
            !invW = zeye(Wmats%Nbp) - matmul(Pmats%screened(:,:,iw,iq),Umats%screened(:,:,iwU,iq))
            ! [ 1 - U*Pi ]
            invW = zeye(Wmats%Nbp) - matmul(Umats%screened(:,:,iwU,iq),Pmats%screened(:,:,iw,iq))
            if(Emats%status)Emats%screened(:,:,iw,iq) = invW
            !
            ! [ 1 - Pi*U ]^-1 or [ 1 - U*Pi ]^-1
            call inv(invW)
            !
            ! U*[ 1 - Pi*U ]^-1
            !Wmats%screened(:,:,iw,iq) = matmul(Umats%screened(:,:,iwU,iq),invW)
            ! [ 1 - U*Pi ]^-1*U
            Wmats%screened(:,:,iw,iq) = matmul(invW,Umats%screened(:,:,iwU,iq))
            !
            !Hermiticity check - print if error is bigger than 1e-3
            if(symQ_) call check_Hermiticity(Wmats%screened(:,:,iw,iq),1e7*eps,enforce=.true.,hardstop=.false.,name="Wlat_w"//str(iw)//"_q"//str(iq),verb=.true.)
            !
            !average the dielectric function around the Gamma point
            if(smear)then
               add_iq = any( AverageList .eq. iq )
               if(add_iq)then
                  epsGamma(:,:,iw) = epsGamma(:,:,iw) + invW
                  iavg = iavg + 1
               endif
            endif
            !
         enddo !iq
         !
         !Control on the number of averaged K-points close to Gamma
         if(smear.and.(iavg.ne.Navg)) stop "calc_W_full: error in the number of avergaged K-points around Gamma."
         !
      enddo !iw
      !$OMP END DO
      !$OMP END PARALLEL
      !
      ! Gamma point handling
      if(smear)then
         !
         write(*,"(A)")"     Smearing Gamma point."
         !
         !W at gamma point should be real
         epsGamma = dreal(epsGamma)/Navg
         !
         !Fill the Gamma point value - element not included in the iq loop - print if error is bigger than 1e-3
         do iw=1,Nmats
            iwU = iw
            if(Ustatic)iwU = 1
            !Wmats%screened(:,:,iw,Umats%iq_gamma) = dreal(matmul(Umats%screened(:,:,iwU,Umats%iq_gamma),epsGamma(:,:,iw)))
            Wmats%screened(:,:,iw,Umats%iq_gamma) = alphaGamma * dreal(matmul(epsGamma(:,:,iw),Umats%screened(:,:,iwU,Umats%iq_gamma)))
            call check_Symmetry(Wmats%screened(:,:,iw,Umats%iq_gamma),eps,enforce=.true.,hardstop=.false.,name="Wlat_w"//str(iw)//"_q"//str(Umats%iq_gamma),verb=.false.)
         enddo
         !
         if(Emats%status)then
            do iw=1,Nmats
               Emats%screened(:,:,iw,Umats%iq_gamma) = epsGamma(:,:,iw)
               call inv(Emats%screened(:,:,iw,Umats%iq_gamma))
            enddo
         endif
         !
         deallocate(epsGamma,AverageList)
         !
      endif
      deallocate(invW)
      !
      ! Fill the local attributes
      call BosonicKsum(Wmats)
      !
      !Check if Wlat is locally hermitian - print if relative error is bigger than 1e-3
      write(*,"(A)") "     Checking hermiticity of local Wlat (enforced)."
      do iw=1,Nmats
         call check_Hermiticity(Wmats%screened_local(:,:,iw),1e7*eps,enforce=.true.,hardstop=.false.,name="Wlat_loc_w"//str(iw),verb=.true.)
      enddo
      !
   end subroutine calc_W_full


   !---------------------------------------------------------------------------!
   !PURPOSE: Lattice inversion to get fully screened interaction - EDMFT
   !---------------------------------------------------------------------------!
   subroutine calc_W_edmft(Wmats,Umats,Pmats,Emats,Lttc,symQ,alpha)
      !
      use parameters
      use linalg, only : zeye, inv
      use utils_misc
      use utils_fields
      use file_io
      use crystal
      use input_vars, only : HandleGammaPoint, alphaGamma
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wmats
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(BosonicField),intent(inout)      :: Emats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: symQ
      real(8),intent(in),optional           :: alpha
      !
      complex(8),allocatable                :: invW(:,:),W_q(:,:)
      complex(8),allocatable                :: epsGamma(:,:,:)
      integer,allocatable                   :: AverageList(:)
      real(8)                               :: Beta,alpha_
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw,iwU
      integer                               :: iavg,Navg
      logical                               :: Ustatic,symQ_,smear,add_iq
      !
      !
      if(verbose)write(*,"(A)") "---- calc_W_edmft"
      !
      !
      ! Check on the input Fields
      if(.not.Wmats%status) stop "calc_W_edmft: Wmats not properly initialized."
      if(.not.Umats%status) stop "calc_W_edmft: Umats not properly initialized."
      if(.not.Pmats%status) stop "calc_W_edmft: Pmats not properly initialized."
      if(Umats%Nkpt.eq.0) stop "calc_W_edmft: Umats k dependent attributes not properly initialized."
      if(Wmats%Nkpt.ne.0) stop "calc_W_edmft: Wmats k dependent attributes are supposed to be unallocated."
      if(Pmats%Nkpt.ne.0) stop "calc_W_edmft: Pmats k dependent attributes are supposed to be unallocated."
      if(Umats%iq_gamma.lt.0) stop "calc_W_edmft: Umats iq_gamma not defined."
      !
      symQ_=.false.
      if(present(symQ))symQ_=symQ
      alpha_=1d0
      if(present(alpha))alpha_=alpha
      Ustatic=.false.
      if(Umats%Npoints.eq.1) Ustatic=.true.
      if(Ustatic)write(*,"(A)")"     Static U bare."
      smear = HandleGammaPoint.gt.0
      !
      Nbp = Wmats%Nbp
      Nkpt = Umats%Nkpt
      Beta = Wmats%Beta
      Nmats = Wmats%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "calc_W_edmft: Either Umats and/or Pmats have different orbital dimension with respect to Wmats."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "calc_W_edmft: Either Umats and/or Pmats have different Beta with respect to Wmats."
      if(Pmats%Npoints.ne.Nmats) stop "calc_W_edmft: Pmats has different number of Matsubara points with respect to Wmats."
      if(Emats%status) call assert_shape(Emats%screened,[Nbp,Nbp,Nmats,Nkpt],"calc_W_edmft","Emats%screened")
      !
      allocate(invW(Nbp,Nbp));invW=czero
      allocate(W_q(Nbp,Nbp));W_q=czero
      call clear_attributes(Wmats)
      !
      if(smear)then
         if(.not.small_ik_stored)call fill_smallk(Lttc%kpt)
         allocate(epsGamma(Nbp,Nbp,Nmats));epsGamma=czero
         AverageList = pack(small_ik(:,1),(small_ik(:,2).le.HandleGammaPoint))
         Navg = size(AverageList)
         if(verbose) write(*,"(A,100I4)")"     Averaging dielectric function over "//str(Navg)//" K-points: ",AverageList
      endif
      !
      ! Assuming that the Polarization vanishes at iw-->inf
      Wmats%bare_local = Umats%bare_local
      !
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iw,iwU,iq,invW,W_q,iavg,add_iq)
      !$OMP DO
      do iw=1,Wmats%Npoints
         !
         iwU = iw
         if(Ustatic)iwU = 1
         !
         iavg=0
         do iq=1,Umats%Nkpt
            !
            !avoid the gamma point
            if((iq.eq.Umats%iq_gamma).and.smear)cycle
            !
            ! [ 1 - Pi*U ]
            !invW = zeye(Umats%Nbp) - matmul(alpha_*Pmats%screened_local(:,:,iw),Umats%screened(:,:,iwU,iq))
            ! [ 1 - U*Pi ]
            invW = zeye(Umats%Nbp) - matmul(Umats%screened(:,:,iwU,iq),alpha_*Pmats%screened_local(:,:,iw))
            if(Emats%status)Emats%screened(:,:,iw,iq) = invW
            !
            ! [ 1 - Pi*U ]^-1 or [ 1 - U*Pi ]^-1
            call inv(invW)
            !
            !  U*[ 1 - U*Pi ]^-1
            !W_q = matmul(Umats%screened(:,:,iwU,iq),invW)
            !  [ 1 - U*Pi ]^-1*U
            W_q = matmul(invW,Umats%screened(:,:,iwU,iq))
            !
            !Hermiticity check - print if error is bigger than 1e-3
            if(symQ_) call check_Hermiticity(W_q,1e7*eps,enforce=.true.,hardstop=.false.,name="Wlat_w"//str(iw)//"_q"//str(iq),verb=.true.)
            !
            !Sum to local attribute
            Wmats%screened_local(:,:,iw) = Wmats%screened_local(:,:,iw) + W_q/Umats%Nkpt
            !
            !average the dielectric function around the Gamma point
            if(smear)then
               add_iq = any( AverageList .eq. iq )
               if(add_iq)then
                  epsGamma(:,:,iw) = epsGamma(:,:,iw) + invW
                  iavg = iavg + 1
               endif
            endif
            !
         enddo !iq
         !
         !Control on the number of averaged K-points close to Gamma
         if(smear.and.(iavg.ne.Navg)) stop "calc_W_edmft: error in the number of avergaged K-points around Gamma."
         !
      enddo !iw
      !$OMP END DO
      !$OMP END PARALLEL
      !
      ! Gamma point handling
      if(smear)then
         !
         write(*,"(A)")"     Smearing Gamma point."
         !
         !W at gamma point should be real
         epsGamma = dreal(epsGamma)/Navg
         !
         !Add the Gamma point value - element not summed in the iq loop - print if error is bigger than 1e-3
         do iw=1,Nmats
            iwU = iw
            if(Ustatic)iwU = 1
            !W_q = dreal(matmul(Umats%screened(:,:,iwU,Umats%iq_gamma),epsGamma(:,:,iw)))
            W_q = alphaGamma * dreal(matmul(epsGamma(:,:,iw),Umats%screened(:,:,iwU,Umats%iq_gamma)))
            call check_Symmetry(W_q,eps,enforce=.true.,hardstop=.false.,name="Wlat_w"//str(iw)//"_q"//str(Umats%iq_gamma),verb=.false.)
            Wmats%screened_local(:,:,iw) = Wmats%screened_local(:,:,iw) + W_q/Nkpt
         enddo
         !
         if(Emats%status)then
            do iw=1,Nmats
               Emats%screened(:,:,iw,Umats%iq_gamma) = epsGamma(:,:,iw)
               call inv(Emats%screened(:,:,iw,Umats%iq_gamma))
            enddo
         endif
         !
         deallocate(epsGamma,AverageList)
         !
      endif
      deallocate(invW,W_q)
      !
      !Check if Wlat is locally hermitian - print if relative error is bigger than 1e-3
      write(*,"(A)") "     Checking hermiticity of local Wlat (enforced)."
      do iw=1,Nmats
         call check_Hermiticity(Wmats%screened_local(:,:,iw),1e7*eps,enforce=.true.,hardstop=.false.,name="Wlat_loc_w"//str(iw),verb=.true.)
      enddo
      !
   end subroutine calc_W_edmft


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes [ 1 - U*Pi ]^-1 * Pi
   !---------------------------------------------------------------------------!
   subroutine calc_chi(Chi,Umats,Pmats,Lttc,pathPk)
      !
      use parameters
      use utils_misc
      use utils_fields
      use crystal
      use linalg, only : zeye, inv
      use input_vars, only : structure, Nkpt_path
      implicit none
      !
      type(BosonicField),intent(inout)      :: Chi
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(inout)           :: Lttc
      character(len=*),intent(in),optional  :: pathPk
      !
      complex(8),allocatable                :: invW(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Norb,Nkpt,Nmats
      integer                               :: iq,iw,iwU
      integer                               :: unit,iorb
      real(8),allocatable                   :: wmats(:)
      complex(8),allocatable                :: Pwk_noG(:,:,:,:),Pwk(:,:,:,:)
      logical                               :: Ustatic
      !
      !
      if(verbose)write(*,"(A)") "---- calc_chi"
      !
      !
      ! Check on the input Fields
      if(.not.Chi%status) stop "calc_chi: Chi not properly initialized."
      if(.not.Umats%status) stop "calc_chi: Umats not properly initialized."
      if(.not.Pmats%status) stop "calc_chi: Pmats not properly initialized."
      if(Chi%Nkpt.eq.0) stop "calc_chi: Chi k dependent attributes not properly initialized."
      if(Umats%Nkpt.eq.0) stop "calc_chi: Umats k dependent attributes not properly initialized."
      if(Pmats%Nkpt.eq.0) stop "calc_chi: Pmats k dependent attributes not properly initialized."
      if(Umats%iq_gamma.lt.0) stop "calc_chi: Umats iq_gamma not defined."
      !
      Ustatic=.false.
      if(Umats%Npoints.eq.1) Ustatic=.true.
      if(Ustatic)write(*,"(A)")"     Static U bare."
      !
      Norb = int(sqrt(dble(Chi%Nbp)))
      Nbp = Chi%Nbp
      Nkpt = Chi%Nkpt
      Beta = Chi%Beta
      Nmats = Chi%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "calc_chi: Either Umats and/or Pmats have different orbital dimension with respect to Chi."
      if(all([Umats%Nkpt-Nkpt,Pmats%Nkpt-Nkpt].ne.[0,0])) stop "calc_chi: Either Umats and/or Pmats have different number of k-points with respect to Chi."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "calc_chi: Either Umats and/or Pmats have different Beta with respect to Chi."
      if(Pmats%Npoints.ne.Nmats) stop "calc_chi: Pmats has different number of Matsubara points with respect to Chi."
      !
      allocate(invW(Nbp,Nbp));invW=czero
      call clear_attributes(Chi)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Pmats,Umats,Chi,Lttc,Ustatic),&
      !$OMP PRIVATE(iw,iwU,iq,invW)
      !$OMP DO
      do iw=1,Chi%Npoints
         !
         iwU = iw
         if(Ustatic)iwU = 1
         !
         do iq=1,Chi%Nkpt
            !
            ! [ 1 - U*Pi ]
            !invW = zeye(Chi%Nbp) - matmul(Umats%screened(:,:,iwU,iq),Pmats%screened(:,:,iw,iq))
            ! [ 1 - Pi*U ]
            invW = zeye(Chi%Nbp) - matmul(Pmats%screened(:,:,iw,iq),Umats%screened(:,:,iwU,iq))
            !
            ! [ 1 - U*Pi ]^-1 or [ 1 - Pi*U ]^-1
            call inv(invW)
            !
            ! [ 1 - U*Pi ]^-1 * Pi
            !Chi%screened(:,:,iw,iq) = matmul(invW,Pmats%screened(:,:,iw,iq))
            ! [ 1 - Pi*U ]^-1 * Pi
            Chi%screened(:,:,iw,iq) = matmul(invW,Pmats%screened(:,:,iw,iq))
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(invW)
      !
      call BosonicKsum(Chi)
      !
      !print along path
      if(present(pathPk).and.(reg(structure).ne."None"))then
         !
         if(.not.Lttc%pathStored)then
            write(*,"(A)") "     calc_chi: re-initializing the K-path."
            call interpolate2Path(Lttc,Nkpt_path,"Hk",store=.true.)
         endif
         !
         allocate(Pwk_noG(Nbp,Nbp,Nmats,Nkpt));Pwk_noG=czero
         do iq=1,Chi%Nkpt
            Pwk_noG(:,:,:,iq) = Chi%screened(:,:,:,iq) - Chi%screened(:,:,:,Umats%iq_gamma)
         enddo
         !
         allocate(Pwk(Nbp,Nbp,Nmats,Lttc%Nkpt_path));Pwk=czero
         call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),Pwk_noG,Pwk)
         deallocate(Pwk_noG)
         !
         !Print susceptibility along path
         call createDir(reg(pathPk),verb=verbose)
         allocate(wmats(Chi%Npoints))
         wmats = BosonicFreqMesh(Beta,Nmats)
         do iq=1,Lttc%Nkpt_path
             unit = free_unit()
             open(unit,file=reg(pathPk)//"ChiC_w_k"//str(iq)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
             do iw=1,Nmats
                write(unit,"(200E20.12)") Lttc%Kpathaxis(iq)/Lttc%Kpathaxis(Lttc%Nkpt_path),wmats(iw),&
               (dreal(  Pwk(iorb+Norb*(iorb-1),iorb+Norb*(iorb-1),iw,iq)+Chi%screened(iorb+Norb*(iorb-1),iorb+Norb*(iorb-1),iw,Umats%iq_gamma)),iorb=1,Norb)
             enddo
             close(unit)
         enddo
         deallocate(wmats,Pwk)
         !
      endif
      !
   end subroutine calc_chi


   !---------------------------------------------------------------------------!
   !PURPOSE: Read momentum and frequency dependent interaction from SPEX files.
   !---------------------------------------------------------------------------!
   subroutine read_U_spex_full(Umats,save2readable,kpt,pathOUTPUT,doAC,HartreeData,correctU0)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use input_vars, only : Nkpt3, RealPrint
      use input_vars, only : pathINPUT, UfullStructure, Uthresh, HandleGammaPoint
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      logical,intent(in)                    :: save2readable
      real(8),intent(in),optional           :: kpt(:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: doAC
      logical,intent(in),optional           :: HartreeData
      logical,intent(in),optional           :: correctU0
      !
      logical                               :: LocalOnly,filexists,ACdone,doAC_
      logical                               :: warn,smear,enforceAll,HartreeData_,correctU0_
      character(len=256)                    :: file_spex,path,pathOUTPUT_
      integer                               :: unit,Nkpt
      integer                               :: iq,iw,iqread,Nbp_spex
      integer                               :: idum,Nspin_spex,Norb_spex,Nfreq
      integer                               :: ib1,ib2,iw1,iw2
      integer                               :: iprint,Nprint
      real(8)                               :: UnitConversion
      real(8),allocatable                   :: wread(:),wmats(:)
      complex(8),allocatable                :: D1(:,:),D2(:,:),D3(:,:)
      complex(8),allocatable                :: Utmp(:,:),UR(:,:,:),URnn(:,:)
      type(BosonicField)                    :: Ureal
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_U_spex_full"
      !
      !
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      LocalOnly = .not.present(kpt)
      smear = .false.
      !
      HartreeData_=.true.
      if(present(HartreeData)) HartreeData_ = HartreeData
      UnitConversion = 1d0
      if(HartreeData_) UnitConversion = H2eV
      !
      correctU0_=.false.
      if(present(correctU0)) correctU0_ = correctU0
      !
      ! Check on the input Boson
      if(.not.Umats%status) stop "read_U_spex_full: BosonicField not properly initialized."
      if((.not.LocalOnly).and.(.not.allocated(Umats%bare))) stop "read_U_spex_full: Requested k-dependence but bare non-local attribute not allocated."
      if((.not.LocalOnly).and.(.not.allocated(Umats%screened))) stop "read_U_spex_full: Requested k-dependence but screened non-local attribute not allocated."
      if(LocalOnly.and.allocated(Umats%bare)) stop "read_U_spex_full: Bare K-dependent attributes is present but not used."
      if(LocalOnly.and.allocated(Umats%screened)) stop "read_U_spex_full: Screened K-dependent attributes is present but not used."
      !
      allocate(wmats(Umats%Npoints));wmats=0d0
      wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      !
      ! for memory demanding calculations the execute_command_line does not work
      ! so I have to create the VW_imag directory before (but I'm keeping it here).
      ! In the end instead of checking for the existence of the VW_imag folder
      ! I will check for the first Q point file within the folder.
      !
      ! Check if the data on the Matsubara axis are present
      path = reg(pathOUTPUT_)//"VW_imag"
      !call inquireDir(reg(path),ACdone,hardstop=.false.,verb=verbose)
      call inquireFile(reg(path)//"/VW.Q0001.DAT",ACdone,hardstop=.false.,verb=verbose)
      doAC_ = .not.ACdone
      if(present(doAC)) doAC_ = doAC .or. doAC_
      !
      !
      ! This has to be done before allocation of large files otherwise the "system" command will not work
      if(doAC_.and.(.not.LocalOnly))then
         call createDir(reg(trim(pathOUTPUT_)//"VW_imag"),verb=verbose)
         if(save2readable)then
            call createDir(reg(trim(pathOUTPUT_)//"VW_imag_readable"),verb=verbose)
            call createDir(reg(trim(pathINPUT)//"VW_real_readable"),verb=verbose)
         endif
      endif
      !
      call init_Uelements(int(sqrt(dble(Umats%Nbp))),PhysicalUelements)
      call clear_attributes(Umats)
      !
      !
      ! Perform cnalytical continuation on real interaction or load existing files
      if(doAC_) then
         !
         if(smear.and.LocalOnly)then
            write(*,"(A)") "     read_U_spex_full: smearing requested with LocalOnly interaction. Smearing ignored."
            smear=.false.
         endif
         if(smear.and.(size(kpt,dim=2).ne.Umats%Nkpt))stop "read_U_spex_full: K-point mesh does not correspond to Umats one."
         !
         !---------------------------------------------------------------------!
         !
         write(*,"(A)")"     Performing Analytic continuation to get UcRPA(iq,iw)."
         !
         ! Allocations from dimensions written in VW.Q0001.DAT file
         path = reg(pathINPUT)//"VW_real/VW.Q0001.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read")
         read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
         close(unit)
         !
         Nbp_spex = Norb_spex**2
         allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         write(*,"(A,I)")"     Real frequencies: ",Nfreq
         if(Nfreq.lt.5) write(*,"(A)") "     Warning: the number of real-frequency points is less than 5, AC will be skipped and UcRPA(iw)=UcRPA(w=0) for each iw."
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "read_U_spex_full: Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Nbp_spex) stop "read_U_spex_full: Size of given BosonicField and VW_real orbital space do not coincide."
         !
         ! Look for the Number of SPEX files. Which are supposed to be ordered.
         Nkpt = 0
         do iq=1,9999
            file_spex = reg(pathINPUT)//"VW_real/VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) exit
            Nkpt = Nkpt + 1
         enddo
         write(*,"(A,I)") "     The number of SPEX files (Nkpt) in VW_real is: ",Nkpt
         if((.not.LocalOnly).and.(Umats%Nkpt.ne.Nkpt)) stop "read_U_spex_full: Number of k-points of given BosonicField and number of VW_real k-points do not coincide."
         !
         ! Allocate the Bosonic field on the real axis
         if(LocalOnly)then
            call AllocateBosonicField(Ureal,Norb_spex,Nfreq,0,Nkpt=0,name="Uspex(w)",Beta=Umats%Beta,no_bare=.true.)
         else
            call AllocateBosonicField(Ureal,Norb_spex,Nfreq,0,Nkpt=Nkpt,name="Uspex(k,w)",Beta=Umats%Beta,no_bare=.true.)
         endif
         !
         ! Read VW_real accumulating local attribute and optionally storing the k-dependent part
         path = reg(pathINPUT)//"VW_real/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,verb=verbose) !redundant control
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
            if (iq.ne.iqread) stop "read_U_spex_full: iqread.ne.iq"
            !
            read(unit) wread
            wread = UnitConversion*wread
            if (dabs(wread(1)).gt.eps) stop "read_U_spex_full: wread(1) not zero"
            !
            do iw=0,Nfreq
               read(unit) Utmp
               if(iw.eq.0) then
                  !V(:,:,iq)=vwtmp(:,:)/Nkpt/Nkpt
                  !bare values on Matsubara are the same so no need to use bare attributes of Ureal
                  Umats%bare_local = Umats%bare_local + UnitConversion*Utmp/(Nkpt**3)
                  if(.not.LocalOnly) Umats%bare(:,:,iq) = UnitConversion*Utmp/(Nkpt**2)
               else
                  !Ur(:,:,iw,iq)=vwtmp(:,:)/Nkpt/Nkpt
                  Ureal%screened_local(:,:,iw) = Ureal%screened_local(:,:,iw) + UnitConversion*Utmp/(Nkpt**3)
                  if(.not.LocalOnly) Ureal%screened(:,:,iw,iq) = UnitConversion*Utmp/(Nkpt**2)
               endif
            enddo
            !
            close(unit)
            !
         enddo !iq
         !
         ! Print out the transformed Ucrpa - local
         write(*,"(A)")
         call dump_BosonicField(Ureal,reg(pathOUTPUT_),"Uloc_real.DAT",wread)
         if((.not.LocalOnly).and.save2readable) call dump_BosonicField(Ureal,reg(pathINPUT)//"VW_real_readable/",.not.save2readable,axis=wread)
         !
         !
         if(LocalOnly)then
            !
            ! Allocate the temporary quantities needed by the Analytical continuation
            allocate(D1(Nbp_spex,Nbp_spex));D1=czero
            allocate(D2(Nbp_spex,Nbp_spex));D2=czero
            allocate(D3(Nbp_spex,Nbp_spex));D3=czero
            !
            ! Analytical continuation of the local component to imag axis using spectral rep
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(Nbp_spex,wmats,wread,Nfreq,Ureal,Umats),&
            !$OMP PRIVATE(ib1,ib2,iw1,iw2,D1,D2,D3,Utmp)
            !$OMP DO
            do iw1=1,Umats%Npoints
               !
               if(Nfreq.ge.5)then
                  !
                  Utmp=czero
                  do iw2=1,Nfreq-2,2
                     !
                     ! Locally D(-w)=-D(w) for all components
                     D1 = -dimag( Ureal%screened_local(:,:,iw2)   )/pi
                     D2 = -dimag( Ureal%screened_local(:,:,iw2+1) )/pi
                     D3 = -dimag( Ureal%screened_local(:,:,iw2+2) )/pi
                     !
                     ! Integrate using Simpson method
                     if(wread(iw2).gt.0.d0) then
                        Utmp = Utmp + ( D1/(dcmplx(0.d0,wmats(iw1))-wread(iw2)  ) - D1/(dcmplx(0.d0,wmats(iw1))+wread(iw2)  ) ) *(wread(iw2+1)-wread(iw2))/3.d0
                        Utmp = Utmp + ( D2/(dcmplx(0.d0,wmats(iw1))-wread(iw2+1)) - D2/(dcmplx(0.d0,wmats(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                        Utmp = Utmp + ( D3/(dcmplx(0.d0,wmats(iw1))-wread(iw2+2)) - D3/(dcmplx(0.d0,wmats(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
                     elseif(dabs(wread(iw2)).lt.1.d-12) then
                        Utmp = Utmp + ( D2/(dcmplx(0.d0,wmats(iw1))-wread(iw2+1)) - D2/(dcmplx(0.d0,wmats(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                        Utmp = Utmp + ( D3/(dcmplx(0.d0,wmats(iw1))-wread(iw2+2)) - D3/(dcmplx(0.d0,wmats(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
                     endif
                  enddo
                  !
                  Umats%screened_local(:,:,iw1) = Utmp + Umats%bare_local
                  !
               else
                  !
                  Umats%screened_local(:,:,iw1) = Ureal%screened_local(:,:,1)
                  !
               endif
               !
            enddo !iw1
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            deallocate(D1,D2,D3)
            !
            if(correctU0_)then
               do iw1=Umats%Npoints,1,-1
                  Umats%screened_local(:,:,iw1) = Umats%screened_local(:,:,iw1) - Umats%screened_local(:,:,1) + Ureal%screened_local(:,:,1)
               enddo
            endif
            !
            ! Remove unwanted components
            if(.not.UfullStructure)then
               do ib1=1,Nbp_spex
                  do ib2=1,Nbp_spex
                     if(.not.PhysicalUelements%Full_All(ib1,ib2))then
                        Umats%bare_local(ib1,ib2) = dcmplx(0d0,0d0)
                        Umats%screened_local(ib1,ib2,:) = dcmplx(0d0,0d0)
                     endif
                  enddo
               enddo
            endif
            !
            write(*,"(A,F)") "     UcRPA(w) --> UcRPA(iw) cpu timing: ", finish-start
            !
         else
            !
            ! Allocate the temporary quantities needed by the Analytical continuation
            allocate(D1(Nbp_spex,Nbp_spex));D1=czero
            allocate(D2(Nbp_spex,Nbp_spex));D2=czero
            allocate(D3(Nbp_spex,Nbp_spex));D3=czero
            !
            ! Analytical continuation of all the K-points to imag axis using spectral rep
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(Nbp_spex,wmats,wread,Nfreq,Ureal,Umats,correctU0_,verbose),&
            !$OMP PRIVATE(iq,ib1,ib2,iw1,iw2,D1,D2,D3,Utmp)
            !$OMP DO
            do iq=1,Umats%Nkpt
               !
               !Perform analytical continuation from real to Matsubara frequency at each iq
               do iw1=1,Umats%Npoints
                  !
                  if(Nfreq.ge.5)then
                     !
                     Utmp=czero
                     do iw2=1,Nfreq-2,2
                        !
                        ! D = i/2pi * { [Uab(w)-Uab(inf)] - [Uba(w)-Uba(inf)]* }
                        D1 = img*( (Ureal%screened(:,:,iw2,iq)  -Umats%bare(:,:,iq)) - transpose(conjg(Ureal%screened(:,:,iw2,iq)  -Umats%bare(:,:,iq))) )/(2d0*pi)
                        D2 = img*( (Ureal%screened(:,:,iw2+1,iq)-Umats%bare(:,:,iq)) - transpose(conjg(Ureal%screened(:,:,iw2+1,iq)-Umats%bare(:,:,iq))) )/(2d0*pi)
                        D3 = img*( (Ureal%screened(:,:,iw2+2,iq)-Umats%bare(:,:,iq)) - transpose(conjg(Ureal%screened(:,:,iw2+2,iq)-Umats%bare(:,:,iq))) )/(2d0*pi)
                        !
                        ! Integrate using Simpson method
                        if(wread(iw2).gt.0.d0) then
                           Utmp = Utmp + (  D1/(dcmplx(0d0,wmats(iw1))-wread(iw2)  )  -  transpose(conjg(D1))/(dcmplx(0d0,wmats(iw1))+wread(iw2)  )  ) * (wread(iw2+1)-wread(iw2))/3d0
                           Utmp = Utmp + (  D2/(dcmplx(0d0,wmats(iw1))-wread(iw2+1))  -  transpose(conjg(D2))/(dcmplx(0d0,wmats(iw1))+wread(iw2+1))  ) * (wread(iw2+1)-wread(iw2))*4d0/3d0
                           Utmp = Utmp + (  D3/(dcmplx(0d0,wmats(iw1))-wread(iw2+2))  -  transpose(conjg(D3))/(dcmplx(0d0,wmats(iw1))+wread(iw2+2))  ) * (wread(iw2+1)-wread(iw2))/3d0
                        elseif(dabs(wread(iw2)).lt.1.d-12) then
                           Utmp = Utmp + (  D2/(dcmplx(0d0,wmats(iw1))-wread(iw2+1))  -  transpose(conjg(D2))/(dcmplx(0d0,wmats(iw1))+wread(iw2+1))  ) * (wread(iw2+1)-wread(iw2))*4d0/3d0
                           Utmp = Utmp + (  D3/(dcmplx(0d0,wmats(iw1))-wread(iw2+2))  -  transpose(conjg(D3))/(dcmplx(0d0,wmats(iw1))+wread(iw2+2))  ) * (wread(iw2+1)-wread(iw2))/3d0
                        endif
                     enddo
                     !
                     Umats%screened(:,:,iw1,iq) = Utmp + Umats%bare(:,:,iq)
                     !
                  else
                     !
                     Umats%screened(:,:,iw1,iq) = Ureal%screened(:,:,1,iq)
                     !
                  endif
                  !
               enddo !iw1
               !
               if(correctU0_)then
                  do iw1=Umats%Npoints,1,-1
                     Umats%screened(:,:,iw1,iq) = Umats%screened(:,:,iw1,iq) - Umats%screened(:,:,1,iq) + Ureal%screened(:,:,1,iq)
                  enddo
               endif
               !
               if(verbose)print *, "     UcRPA(q,iw) - done iq: ",iq
            enddo !iq
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            deallocate(D1,D2,D3)
            !
            write(*,"(A,F)") "     UcRPA(q,w) --> UcRPA(q,iw) cpu timing: ", finish-start
            !
            !Remove unwanted components
            if(.not.UfullStructure)then
               do ib1=1,Nbp_spex
                  do ib2=1,Nbp_spex
                     if(.not.PhysicalUelements%Full_All(ib1,ib2))then
                        Umats%bare(ib1,ib2,:) = dcmplx(0d0,0d0)
                        Umats%screened(ib1,ib2,:,:) = dcmplx(0d0,0d0)
                     endif
                  enddo
               enddo
            endif
            !
         endif !LocalOnly
         deallocate(Utmp)
         !
         !Symmetry Checks on the Ucrpa K-dependent attributes and update local ones
         write(*,"(A)")
         if(.not.LocalOnly)then
            !
            enforceAll=.true.
            write(*,"(A)") "     Symmetry check on bare Ucrpa - should be Hermitian at all iq."
            do iq=1,Umats%Nkpt
               call check_Hermiticity(Umats%bare(:,:,iq),eps,enforce=enforceAll,hardstop=.false.,name="Ucrpa_bare_q"//str(iq))
            enddo
            !
            write(*,"(A)") "     Symmetry check on screened Ucrpa - should be Hermitian at all iw and all iq."
            do iw=1,Umats%Npoints
               do iq=1,Umats%Nkpt
                  call check_Hermiticity(Umats%screened(:,:,iw,iq),eps,enforce=enforceAll,hardstop=.false.,name="Ucrpa_screened_w"//str(iw)//"_q"//str(iq))
               enddo
            enddo
            !
            !Removal of the divergence at Gamma
            if(smear)then
               write(*,"(A)") "     Smearing UcRPA."
               if(.not.small_ik_stored)call fill_smallk(kpt)
               call correct_Gamma(Umats,abs(HandleGammaPoint))
            endif
            !
            call BosonicKsum(Umats)
            !
         endif
         !
         !Symmetry Checks on the Ucrpa local attributes
         write(*,"(A)") "     Symmetry check on Uinst - should be Hermitian."
         call check_Hermiticity(Umats%screened_local(:,:,1),eps,enforce=.false.,hardstop=.false.,name="Uinst")
         !
         write(*,"(A)") "     Symmetry check on local bare Ucrpa - should be Hermitian."
         call check_Hermiticity(Umats%bare_local,eps,enforce=.false.,hardstop=.false.,name="Ucrpa_bare_local")
         !
         write(*,"(A)") "     Symmetry check on local screened Ucrpa - should be Hermitian at all iw."
         do iw=1,Umats%Npoints
            call check_Hermiticity(Umats%screened_local(:,:,iw),eps,enforce=.false.,hardstop=.false.,name="Ucrpa_screened_local_w"//str(iw))
         enddo
         !
         !Remove local components under threshold and check for inverted Im/Re symmetry
         warn=.true.
         do ib1=1,Nbp_spex
            do ib2=1,Nbp_spex
               !
               if(abs(Umats%bare_local(ib1,ib2)).lt.Uthresh)then
                  Umats%bare_local(ib1,ib2)=czero
                  Umats%screened_local(ib1,ib2,:)=czero
               endif
               !
               if(abs(dreal(Umats%bare_local(ib1,ib2))).lt.abs(dimag(Umats%bare_local(ib1,ib2))))then
                  if(warn) write(*,"(A)")new_line("A")//"     Warning: inverted Re/Im parity in local UcRPA. Check that orbital indexes belong to different sites."
                  warn=.false.
                  write(*,"(A,4I4,A,2I4)")"     Orbitals: ",PhysicalUelements%Full_Map(ib1,ib2,:)," Element:  ",ib1,ib2
                  write(*,"(2(A,F))")     "     Re[Ubare(w=inf)]: ",dreal(Umats%bare_local(ib1,ib2))," Im[Ubare(w=inf)]: ",dimag(Umats%bare_local(ib1,ib2))
                  !stop "Something wrong: Uloc cannot have inverted Re/Im parity."
               endif
               !
            enddo
         enddo
         !
         ! Print out the transformed Ucrpa
         write(*,"(A)")
         call dump_BosonicField(Umats,reg(pathOUTPUT_),"Uloc_mats.DAT")
         if(.not.LocalOnly)then
            call dump_BosonicField(Umats,reg(pathOUTPUT_)//"VW_imag/",.true.)
            if(save2readable)call dump_BosonicField(Umats,reg(pathOUTPUT_)//"VW_imag_readable/",.not.save2readable)
         endif
         !
         deallocate(wread)
         !
         if(verbose)call checkAnalyticContinuation(Umats,Ureal)
         call DeallocateBosonicField(Ureal)
         !
         !---------------------------------------------------------------------!
         !
      else
         !
         !---------------------------------------------------------------------!
         !
         write(*,"(A)")"     Reading UcRPA(q,iw) from "//reg(pathOUTPUT_)//"VW_imag/"
         !
         ! Allocations from dimensions written in W.Q0001.DAT file
         path = reg(pathOUTPUT_)//"VW_imag/VW.Q0001.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         !
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read")
         read(unit)idum,Nspin_spex,Norb_spex,Nfreq
         close(unit)
         !
         Nbp_spex = Norb_spex**2
         allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "read_U_spex_full: Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Nbp_spex) stop "read_U_spex_full: Size of given BosonicField and VW_imag orbital space do not coincide."
         if(Umats%Npoints.ne.Nfreq)then
            write(*,"(A)")"     Warning - Files grid: "//str(Nfreq)//" does not match with U mesh: "//str(Umats%Npoints)//". Reading up to the smaller."
            !stop "read_U_spex_full: Number of VW_imag Matsubara points and bosonic field mesh does not coincide."
         else
            write(*,"(A,I)")"     Matsubara frequencies: ",Nfreq
         endif
         !
         ! Look for the Number of SPEX files. Which are supposed to be ordered.
         Nkpt = 0
         do iq=1,9999
            file_spex = reg(pathOUTPUT_)//"VW_imag/VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) exit
            Nkpt = Nkpt + 1
         enddo
         write(*,"(A,I)") "     The number of SPEX files (Nkpt) in VW_imag is: ",Nkpt
         if((.not.LocalOnly).and.(Umats%Nkpt.ne.Nkpt)) stop "read_U_spex_full: Number of k-points of given BosonicField and number of VW_imag k-points do not coincide."
         !
         ! Read VW_imag accumulating local attribute and optionally storing the k-dependent part
         path = reg(pathOUTPUT_)//"VW_imag/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,verb=verbose) !redundant control
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
            if (iq.ne.iqread) stop "read_U_spex_full: iqread.ne.iq"
            !
            read(unit) wread
            wread = wread
            do iw=1,Umats%Npoints
               if(dabs(wread(iw)-wmats(iw)).gt.eps) Then
                  write(*,"(F)")dabs(wread(iw)-wmats(iw)),iw,iq
                  stop "read_U_spex_full: wread.ne.wmats"
               endif
            enddo
            !
            do iw=0,Nfreq
               read(unit) Utmp
               if(iw.le.Umats%Npoints)then
                  if(iw.eq.0) then
                     Umats%bare_local = Umats%bare_local +Utmp/(Nkpt**3)
                     if(.not.LocalOnly) Umats%bare(:,:,iq) = Utmp/(Nkpt**2)
                  else
                     Umats%screened_local(:,:,iw) = Umats%screened_local(:,:,iw) + Utmp/(Nkpt**3)
                     if(.not.LocalOnly) Umats%screened(:,:,iw,iq) = Utmp/(Nkpt**2)
                  endif
               endif
            enddo
            !
            close(unit)
            !
         enddo !iq
         !
         deallocate(wread,Utmp)
         !
         !---------------------------------------------------------------------!
         !
      endif
      !
      ! Print the nn non local interaction
      if((.not.LocalOnly).and.doAC_)then
         !
         Nprint = size(RealPrint,dim=2)
         allocate(URnn(int(sqrt(dble(Umats%Nbp))),int(sqrt(dble(Umats%Nbp)))));URnn=czero
         allocate(UR(Umats%Nbp,Umats%Nbp,Nprint));UR=czero
         call wannier_K2R_NN(RealPrint,Nkpt3,kpt,Umats%screened(:,:,1,:),UR)
         where(abs((UR))<eps) UR=czero
         do iprint=1,Nprint
            call dump_Matrix(UR(:,:,iprint),reg(trim(pathOUTPUT_)),"U_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
            call product2NN(UR(:,:,iprint),URnn)
            call dump_Matrix(URnn,reg(trim(pathOUTPUT_)),"Unn_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
         enddo
         deallocate(UR,URnn)
         !
      endif
      !
   end subroutine read_U_spex_full
   !
   subroutine read_U_spex_Uloc0(Umat,pathOUTPUT,HartreeData)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : pathINPUT
      implicit none
      !
      complex(8),allocatable,intent(inout)  :: Umat(:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: HartreeData
      !
      logical                               :: Umatsxists,Urealxists,SPEXxists,HartreeData_
      character(len=256)                    :: file_spex,path,pathOUTPUT_
      real(8)                               :: UnitConversion
      integer                               :: unit
      integer                               :: iq,iw,Nbp_spex,Nkpt
      integer                               :: iqread,Nspin_spex,Norb_spex,Nfreq
      complex(8),allocatable                :: Utmp(:,:)
      type(BosonicField)                    :: Uread
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_U_spex_Uloc0"
      !
      !
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      HartreeData_=.true.
      if(present(HartreeData)) HartreeData_ = HartreeData
      UnitConversion = 1d0
      if(HartreeData_) UnitConversion = H2eV
      !
      ! Look for Uloc_mats.DAT and Uloc_real.DAT
      path=reg(pathINPUT)//"Uloc_mats.DAT"
      call inquireFile(reg(path),Umatsxists,hardstop=.false.,verb=verbose)
      path=reg(pathINPUT)//"Uloc_real.DAT"
      call inquireFile(reg(path),Urealxists,hardstop=.false.,verb=verbose)
      path=reg(pathINPUT)//"VW_real" !/VW.Q0001.DAT"
      call inquireDir(reg(path),SPEXxists,hardstop=.false.,verb=verbose)
      !
      if(Umatsxists)then
         !
         call AllocateBosonicField(Uread,size(Umat,dim=1),1,0)
         !
         call read_BosonicField(Uread,reg(pathINPUT),"Uloc_mats.DAT")
         Umat = Uread%screened_local(:,:,1)
         call dump_matrix(Umat,reg(pathINPUT),"Umat.DAT")
         call DeallocateBosonicField(Uread)
         !
      else if(Urealxists)then
         !
         call AllocateBosonicField(Uread,size(Umat,dim=1),1,0)
         !
         call read_BosonicField(Uread,reg(pathINPUT),"Uloc_real.DAT")
         Umat = Uread%screened_local(:,:,1)
         call dump_matrix(Umat,reg(pathINPUT),"Umat.DAT")
         call DeallocateBosonicField(Uread)
         !
      else if(SPEXxists)then
         !
         Nkpt=0
         path = reg(pathINPUT)//"VW_real/"
         do iq=1,9999
            !
            file_spex = reg(path)//"VW.Q"//str(iq,4)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"VW_real/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),SPEXxists,hardstop=.false.,verb=verbose)
            !
            ! This mathod looks uglier than the previous but here I don't
            ! have to store K-dep, just sum up what's present.
            if(.not.SPEXxists)cycle
            Nkpt=Nkpt+1
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit)iqread,Nspin_spex,Norb_spex,Nfreq
            !
            Nbp_spex = Norb_spex**2
            allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
            if(iq.eq.1)call assert_shape(Umat,[Nbp_spex,Nbp_spex],"read_spex_Uloc","Umat")
            !
            read(unit) !wread
            !
            do iw=0,1
               read(unit) Utmp
               if(iw.eq.1) then
                  Umat = Umat + UnitConversion*Utmp
               endif
            enddo
            !
            close(unit)
            deallocate(Utmp)
            !
         enddo !iq
         Umat = Umat/(Nkpt**3)
         !
      else
         stop "read_U_spex_Uloc0: No useful interaction file found."
      endif
      !
      call dump_matrix(Umat,reg(pathOUTPUT_),"Umat.DAT")
      !
   end subroutine read_U_spex_Uloc0


   !---------------------------------------------------------------------------!
   !PURPOSE: Read space and real frequency dependent fully screened interaction
   !         from custom respack files, un-screen with GG
   !---------------------------------------------------------------------------!
   subroutine read_U_respack_full(Umats,LocalOnly,Lttc,pathOUTPUT)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use linalg, only : dag
      use input_vars, only : pathINPUT, U_AC, RespackRthresh
      use input_vars, only : Nkpt_path, structure
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      logical,intent(in)                    :: LocalOnly
      type(Lattice),intent(inout)           :: Lttc
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      logical                               :: filexists,exitRloop,FTdone
      character(len=256)                    :: file_respack,pathOUTPUT_
      integer                               :: unit,NRW,NRJ
      integer                               :: iR,iw,iRread,iq
      integer                               :: Nspin_respack,Norb_respack,Nfreq
      integer                               :: ib1,ib2,Nbp,iorb,jorb,Norb
      integer                               :: NwigMat,iwigMat,iwig,iwig_
      integer                               :: nx,ny,nz,Nvec(3)
      integer,allocatable                   :: wigdone(:)
      real(8),allocatable                   :: wread(:)
      complex(8),allocatable                :: Rmat(:,:),Vr(:,:,:),Ur(:,:,:,:)
      type(BosonicField)                    :: Ureal
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_U_respack_full"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      ! Check on the input Boson
      if(.not.Umats%status) stop "read_U_respack_full: BosonicField not properly initialized."
      if((.not.LocalOnly).and.(.not.allocated(Umats%bare))) stop "read_U_respack_full: Requested k-dependence but bare non-local attribute not allocated."
      if((.not.LocalOnly).and.(.not.allocated(Umats%screened))) stop "read_U_respack_full: Requested k-dependence but screened non-local attribute not allocated."
      if(LocalOnly.and.allocated(Umats%bare)) stop "read_U_respack_full: Bare K-dependent attributes is present but not used."
      if(LocalOnly.and.allocated(Umats%screened)) stop "read_U_respack_full: Screened K-dependent attributes is present but not used."
      if(.not.Lttc%status) stop "read_U_respack_full: Lattice not properly initialized."
      if(Umats%Nbp.ne.(Lttc%Norb**2)) stop "read_U_respack_full: Umats has different orbital dimension with respect to Lttc."
      !
      Norb = Lttc%Norb
      Nbp = Umats%Nbp
      !
      call inquireFile(reg(pathINPUT)//"VW_real/VW.Q0001.DAT",FTdone,hardstop=.false.,verb=verbose)
      if(.not.FTdone)then
         !
         !recover the vectors in real space
         if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
         !
         ! Look for the Number of respack files. Which are supposed to be ordered.
         NRW = 0
         do iR=1,9999
            file_respack = reg(pathINPUT)//"Respack/dir-intW_R/VW.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) cycle
            NRW = NRW + 1
         enddo
         if(NRW.eq.0) stop "read_U_respack_full: no dir-intW_R/VW.R* file found."
         write(*,"(A,I)") "     The number of RESPACK files (NRW) in dir-intW_R is: ",NRW
         !
         NRJ = 0
         do iR=1,9999
            file_respack = reg(pathINPUT)//"Respack/dir-intJ_R/XJ.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) cycle
            NRJ = NRJ + 1
         enddo
         write(*,"(A,I)") "     The number of RESPACK files (NRJ) in dir-intJ_R is: ",NRJ
         !
         !
         NwigMat = Nwig
         if(LocalOnly) NwigMat = 1
         !
         ! Read respack files to product basis form
         allocate(Rmat(Norb,Norb));Rmat=czero
         allocate(Vr(Nbp,Nbp,NwigMat));Vr=czero
         !
         ! Read density-density interaction
         exitRloop=.false.
         intWRloop: do iR=1,NRW
            !
            ! File containing density-density components
            file_respack = reg(pathINPUT)//"Respack/dir-intW_R/VW.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) cycle intWRloop
            unit = free_unit()
            open(unit,file=reg(file_respack),form="unformatted",action="read")
            !
            ! Get wigner-seitz vector
            read(unit) nx,ny,nz
            iwig = find_vec([nx,ny,nz],Nvecwig,hardstop=.false.)
            !
            ! Allocate real-frequency dependent arrays
            read(unit) iRread,Nspin_respack,Norb_respack,Nfreq
            if(iR.ne.iRread) stop "read_U_respack_full: dir-intW_R, iR.ne.iRread"
            if(Norb.ne.Norb_respack) stop "read_U_respack_full: dir-intW_R, Norb.ne.Norb_respack"
            if(iR.eq.1)then
               allocate(Ur(Nbp,Nbp,Nfreq,NwigMat));Ur=czero
               allocate(wread(Nfreq));wread=0d0
            endif
            read(unit) wread
            !
            ! Cycle conditions
            iwigMat = iwig
            if(iwig.eq.0)then
               if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"], not available."
               cycle intWRloop
            endif
            !
            if(LocalOnly)then
               if(iwig.eq.wig0)then
                  if(any([nx,ny,nz].ne.[0,0,0])) stop "read_U_respack_full: wrong index of R=0 vector (intW_R)." !redundant
                  iwigMat = 1
                  exitRloop = .true.
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, found W(R=0) in file: "//reg(file_respack)
               else
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, skipping file: "//reg(file_respack)
                  cycle intWRloop
               endif
            else
               if(radiuswig(iwigMat).gt.RespackRthresh)then
                  if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" is beyond radius threshold ("//str(radiuswig(iwigMat))//")"
                  cycle intWRloop
               else
                  write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"] iwig:"//str(iwigMat)
               endif
            endif
            !
            ! Read data from file
            do iw=0,Nfreq
               !
               Rmat=czero
               read(unit) Rmat
               !
               if(iw.eq.0) then
                  ! Fill in bare limit
                  do iorb=1,Norb
                     do jorb=1,Norb
                        call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                        Vr(ib1,ib2,iwigMat) = Rmat(iorb,jorb)
                     enddo
                  enddo
               else
                  ! Fill in screened part
                  do iorb=1,Norb
                     do jorb=1,Norb
                        call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                        Ur(ib1,ib2,iw,iwigMat) = Rmat(iorb,jorb)
                     enddo
                  enddo
               endif
               !
            enddo
            !
            if(exitRloop)then
               write(*,"(A)")"     Local density-density interaction matrix stored."
               exit intWRloop
            endif
            !
         enddo intWRloop
         !
         ! Read exchange interaction
         exitRloop=.false.
         intJRloop: do iR=1,NRJ
            !
            ! File containing exchange components
            file_respack = reg(pathINPUT)//"Respack/dir-intJ_R/XJ.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) cycle intJRloop
            unit = free_unit()
            open(unit,file=reg(file_respack),form="unformatted",action="read")
            !
            ! Fet wigner-seitz vector
            read(unit) nx,ny,nz
            iwig = find_vec([nx,ny,nz],Nvecwig,hardstop=.false.)
            !
            ! Allocate real-frequency dependent arrays
            read(unit) iRread,Nspin_respack,Norb_respack,Nfreq
            if(iR.ne.iRread) stop "read_U_respack_full: dir-intJ_R, iR.ne.iRread"
            if(Norb.ne.Norb_respack) stop "read_U_respack_full: dir-intJ_R, Norb.ne.Norb_respack"
            read(unit) wread
            !
            ! Cycle conditions
            iwigMat = iwig
            if(iwig.eq.0)then
               if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"], not available."
               cycle intJRloop
            endif
            !
            if(LocalOnly)then
               if(iwig.eq.wig0)then
                  if(any([nx,ny,nz].ne.[0,0,0])) stop "read_U_respack_full: wrong index of R=0 vector (intJ_R)." !redundant
                  iwigMat = 1
                  exitRloop = .true.
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, found J(R=0) in file: "//reg(file_respack)
               else
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, skipping file: "//reg(file_respack)
                  cycle intJRloop
               endif
            else
               if(radiuswig(iwigMat).gt.RespackRthresh)then
                  if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" is beyond radius threshold ("//str(radiuswig(iwigMat))//")"
                  cycle intJRloop
               else
                  write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"] iwig:"//str(iwigMat)
               endif
            endif
            !
            ! Read data from file
            do iw=0,Nfreq
               !
               Rmat=czero
               read(unit) Rmat
               !
               if(iw.eq.0) then
                  ! Fill in bare limit
                  do iorb=1,Norb
                     do jorb=1,Norb
                        if(iorb.ne.jorb)then
                           call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                           Vr(ib1,ib2,iwigMat) = Rmat(iorb,jorb)
                           call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                           Vr(ib1,ib2,iwigMat) = Rmat(iorb,jorb)
                        endif
                     enddo
                  enddo
               else
                  ! Fill in screened part
                  do iorb=1,Norb
                     do jorb=1,Norb
                        if(iorb.ne.jorb)then
                           call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                           Ur(ib1,ib2,iw,iwigMat) = Rmat(iorb,jorb)
                           call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                           Ur(ib1,ib2,iw,iwigMat) = Rmat(iorb,jorb)
                        endif
                     enddo
                  enddo
               endif
               !
            enddo
            !
            if(exitRloop)then
               write(*,"(A)")"     Local exchange interaction matrix stored."
               exit intJRloop
            endif
            !
         enddo intJRloop
         deallocate(Rmat)
         !
         ! FT transform to momentum space
         if(LocalOnly)then
            !
            ! Cleanup respack data: U_ab(0,w) = U_ba(0,w)
            Ur(:,:,1,1) = dcmplx(dreal(Ur(:,:,1,wig0)),0d0)
            call check_Symmetry(Vr(:,:,1),eps,enforce=.true.,hardstop=.false.,verb=.false.)
            do iw=1,Nfreq
               call check_Symmetry(Ur(:,:,iw,1),eps,enforce=.true.,hardstop=.false.,verb=.false.)
            enddo
            !
            call AllocateBosonicField(Ureal,Norb,Nfreq,1)
            !
            Ureal%screened_local = Ur(:,:,:,1)
            Ureal%bare_local = Vr(:,:,1)
            !
         else
            !
            ! Cleanup respack data: U_ab(0,w) = U_ba(0,w)
            Ur(:,:,1,wig0) = dcmplx(dreal(Ur(:,:,1,wig0)),0d0)
            call check_Symmetry(Vr(:,:,wig0),eps,enforce=.true.,hardstop=.false.,verb=.false.)
            do iw=1,Nfreq
               call check_Symmetry(Ur(:,:,iw,wig0),eps,enforce=.true.,hardstop=.false.,verb=.false.)
            enddo
            !
            ! Cleanup respack data
            allocate(wigdone(NwigMat));wigdone=-1
            fixUloop: do iwig=1,NwigMat
               !
               ! Given R retrieve -R
               Nvec = Nvecwig(:,iwig)
               iwig_ = find_vec(-1*Nvec,Nvecwig,hardstop=.false.)
               !
               if(iwig_.eq.0) cycle fixUloop
               !
               if(iwig.eq.1)then
                  !
                  ! Store -R index in wigdone
                  wigdone(iwig) = iwig_
                  !
                  ! Symmetrize
                  Vr(:,:,iwig) = dcmplx(dreal(Vr(:,:,iwig)),-dimag(Vr(:,:,iwig))) !phase fix
                  call check_Symmetry(Vr(:,:,iwig),eps,enforce=.true.,hardstop=.false.,verb=.false.) ! U_ab(R,w) = U_ba(R,w)
                  Vr(:,:,iwig_) = dag(Vr(:,:,iwig)) ! U_ab(-R,w) = U_ba*(R,w)
                  do iw=1,Nfreq
                     Ur(:,:,iw,iwig) = dcmplx(dreal(Ur(:,:,iw,iwig)),-dimag(Ur(:,:,iw,iwig))) !phase fix
                     call check_Symmetry(Ur(:,:,iw,iwig),eps,enforce=.true.,hardstop=.false.,verb=.false.) ! U_ab(R,w) = U_ba(R,w)
                     Ur(:,:,iw,iwig_) = dag(Ur(:,:,iw,iwig)) ! U_ab(-R,w) = U_ba*(R,w)
                  enddo
                  !
               else
                  !
                  ! Look if R in -R has already been treated
                  if(any(wigdone.eq.iwig))then
                     cycle fixUloop
                  else
                     !
                     ! Store -R index in wigdone
                     wigdone(iwig) = iwig_
                     !
                     ! Symmetrize
                     Vr(:,:,iwig) = dcmplx(dreal(Vr(:,:,iwig)),-dimag(Vr(:,:,iwig))) !phase fix
                     call check_Symmetry(Vr(:,:,iwig),eps,enforce=.true.,hardstop=.false.,verb=.false.) ! U_ab(R,w) = U_ba(R,w)
                     Vr(:,:,iwig_) = dag(Vr(:,:,iwig)) ! U_ab(-R,w) = U_ba*(R,w)
                     do iw=1,Nfreq
                        Ur(:,:,iw,iwig) = dcmplx(dreal(Ur(:,:,iw,iwig)),-dimag(Ur(:,:,iw,iwig))) !phase fix
                        call check_Symmetry(Ur(:,:,iw,iwig),eps,enforce=.true.,hardstop=.false.,verb=.false.) ! U_ab(R,w) = U_ba(R,w)
                        Ur(:,:,iw,iwig_) = dag(Ur(:,:,iw,iwig)) ! U_ab(-R,w) = U_ba*(R,w)
                     enddo
                     !
                  endif
                  !
               endif
               !
            enddo fixUloop
            deallocate(wigdone)
            !
            call AllocateBosonicField(Ureal,Norb,Nfreq,1,Nkpt=Lttc%Nkpt)
            !
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Vr,Ureal%bare)
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Ur,Ureal%screened)
            !
            call BosonicKsum(Ureal)
            !
            write(*,"(A)") "     Symmetry check on bare Ucrpa_respack - should be symmetric at all iq."
            do iq=1,Ureal%Nkpt
               call check_Symmetry(Ureal%bare(:,:,iq),eps,enforce=.false.,hardstop=.false.,name="Ucrpa_respack_bare_q"//str(iq))
            enddo
            !
            write(*,"(A)") "     Symmetry check on screened Ucrpa_respack - should be symmetric at all w and all iq."
            do iw=1,Ureal%Npoints
               do iq=1,Ureal%Nkpt
                  call check_Symmetry(Ureal%screened(:,:,iw,iq),eps,enforce=.false.,hardstop=.false.,name="Ucrpa_respack_screened_w"//str(iw)//"_q"//str(iq))
               enddo
            enddo
            !
         endif
         !
         write(*,"(A)") "     Symmetry check on Uinst_respack - should be symmetric."
         call check_Symmetry(Ureal%screened_local(:,:,1),eps,enforce=.false.,hardstop=.false.,name="Uinst")
         !
         write(*,"(A)") "     Symmetry check on local bare Ucrpa_respack - should be symmetric."
         call check_Symmetry(Ureal%bare_local,eps,enforce=.false.,hardstop=.false.,name="Ucrpa_respack_bare_local")
         !
         write(*,"(A)") "     Symmetry check on local screened Ucrpa_respack - should be symmetric at all iw."
         do iw=1,Ureal%Npoints
            call check_Symmetry(Ureal%screened_local(:,:,iw),eps,enforce=.false.,hardstop=.false.,name="Ucrpa_respack_screened_local_w"//str(iw))
         enddo
         !
         deallocate(Vr,Ur)
         !
         ! Store interaction on real-frequency axis
         call dump_BosonicField(Ureal,reg(pathINPUT)//"VW_real/",.true.,axis=wread)
         if(verbose)call dump_BosonicField(Ureal,reg(pathINPUT)//"VW_real_readable/",.false.,axis=wread)
         call DeallocateBosonicField(Ureal)
         deallocate(wread)
         write(*,"(A)")"     Momentum-dependent interaction on the real frequency axis written to file (SPEX format)."
         !
      endif
      !
      ! Calling spex subroutine
      call read_U_spex(Umats,save2readable=verbose,kpt=Lttc%kpt,doAC=U_AC,pathOUTPUT=reg(pathOUTPUT_),HartreeData=.false.,correctU0=.true.)
      !
      !print along path
      if(reg(structure).ne."None")then
         call interpolate2Path(Lttc,Nkpt_path,"Uk",pathOUTPUT=reg(pathINPUT),store=.false.,skipAkw=.true.,data_in=Umats%screened(:,:,1,:))
      endif
      !
   end subroutine read_U_respack_full
   !
   subroutine read_U_respack_Uloc0(Umat)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use input_vars, only : pathINPUT
      implicit none
      !
      complex(8),allocatable,intent(inout)  :: Umat(:,:)
      !
      logical                               :: filexists
      character(len=256)                    :: file_respack
      integer                               :: unit,NRW,NRJ
      integer                               :: iR,iw,iRread
      integer                               :: Nspin_respack,Norb_respack,Nfreq
      integer                               :: ib1,ib2,Nbp,iorb,jorb,Norb
      integer                               :: nx,ny,nz
      real(8),allocatable                   :: wread(:)
      complex(8),allocatable                :: Rmat(:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- read_U_respack_Uloc0"
      !
      !
      ! Check on the input field
      Nbp = size(Umat,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      !
      call assert_shape(Umat,[Nbp,Nbp],"read_U_respack_Uloc0","Umat")
      !
      ! Look for the Number of respack files. Which are supposed to be ordered.
      NRW = 0
      do iR=1,9999
         file_respack = reg(pathINPUT)//"Respack/dir-intW_R/VW.R."//str(iR,4)
         call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
         if(.not.filexists) exit
         NRW = NRW + 1
      enddo
      if(NRW.eq.0) stop "read_U_respack_Uloc0: no dir-intW_R/VW.R* file found."
      write(*,"(A,I)") "     The number of RESPACK files (NRW) in dir-intW_R is: ",NRW
      !
      NRJ = 0
      do iR=1,9999
         file_respack = reg(pathINPUT)//"Respack/dir-intJ_R/XJ.R."//str(iR,4)
         call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
         if(.not.filexists) exit
         NRJ = NRJ + 1
      enddo
      write(*,"(A,I)") "     The number of RESPACK files (NRJ) in dir-intJ_R is: ",NRJ
      !
      allocate(Rmat(Norb,Norb));Rmat=czero
      !
      ! Read density-density interaction
      do iR=1,NRW
         !
         ! File containing density-density components
         file_respack = reg(pathINPUT)//"Respack/dir-intW_R/VW.R."//str(iR,4)
         call inquireFile(reg(file_respack),filexists,verb=verbose) !redundant control
         unit = free_unit()
         open(unit,file=reg(file_respack),form="unformatted",action="read")
         !
         ! Get wigner-seitz vector
         read(unit) nx,ny,nz
         if(any([nx,ny,nz].ne.[0,0,0])) cycle
         !
         read(unit) iRread,Nspin_respack,Norb_respack,Nfreq
         !if(iR.ne.iRread) stop "read_U_respack_Uloc0: dir-intW_R, iR.ne.iRread"
         if(Norb.ne.Norb_respack) stop "read_U_respack_Uloc0: dir-intW_R, Norb.ne.Norb_respack"
         read(unit) wread
         !
         ! Read data from file
         do iw=0,1
            !
            Rmat=czero
            read(unit) Rmat
            !
            ! Take only the w=0 value
            if(iw.eq.1) then
               ! Fill in bare limit
               do iorb=1,Norb
                  do jorb=1,Norb
                     call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                     Umat(ib1,ib2) = Rmat(iorb,jorb)
                  enddo
               enddo
            endif
            !
         enddo
         !
      enddo
      !
      !
      ! Read exchange interaction
      do iR=1,NRJ
         !
         ! file containing exchange components
         file_respack = reg(pathINPUT)//"Respack/dir-intJ_R/XJ.R."//str(iR,4)
         call inquireFile(reg(file_respack),filexists,verb=verbose) !redundant control
         unit = free_unit()
         open(unit,file=reg(file_respack),form="unformatted",action="read")
         !
         ! Get wigner-seitz vector
         read(unit) nx,ny,nz
         if(any([nx,ny,nz].ne.[0,0,0])) cycle
         !
         ! Allocate real-frequency dependent arrays
         read(unit) iRread,Nspin_respack,Norb_respack,Nfreq
         !if(iR.ne.iRread) stop "read_U_respack_full: dir-intJ_R, iR.ne.iRread"
         if(Norb.ne.Norb_respack) stop "read_U_respack_full: dir-intJ_R, Norb.ne.Norb_respack"
         read(unit) wread
         !
         ! Read data from file
         do iw=0,1
            !
            Rmat=czero
            read(unit) Rmat
            !
            ! Take only the w=0 value
            if(iw.eq.1) then
               ! Fill in bare limit
               do iorb=1,Norb
                  do jorb=1,Norb
                     if(iorb.ne.jorb)then
                        call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                        Umat(ib1,ib2) = Rmat(iorb,jorb)
                        call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                        Umat(ib1,ib2) = Rmat(iorb,jorb)
                     endif
                  enddo
               enddo
            endif
            !
         enddo
         !
      enddo
      deallocate(Rmat)
      !
   end subroutine read_U_respack_Uloc0


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute momentum dependent interaction from VASP files.
   !---------------------------------------------------------------------------!
   subroutine read_U_vasp_full(Umats,Lttc)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use input_vars, only : RealPrint, Nkpt_path, structure, SiteOrbs
      use input_vars, only : pathINPUTtr, pathINPUT
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      type(Lattice),intent(inout)           :: Lttc
      !
      integer                               :: i,j,isite,jsite
      integer                               :: iorb,jorb,ib1,ib2
      integer                               :: Norb,Nbp
      integer                               :: iD,iR,iwig,ik,unit,idum
      integer                               :: iprint,Nprint
      real(8)                               :: Rdist
      real(8),allocatable                   :: Ruc(:,:),ReadLine(:),Rvec(:)
      integer,allocatable                   :: Rorder(:),Dist(:,:),DistList(:)
      real(8),allocatable                   :: Rsorted(:,:),Rsorted_bkp(:,:)
      real(8),allocatable                   :: U_iijj(:,:),U_ijji(:,:),U_ijij(:,:)
      complex(8),allocatable                :: Ur(:,:,:),Uk(:,:,:),URnn(:,:)
      character(len=255)                    :: file
      real                                  :: start,finish
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_U_vasp_full"
      !
      !
      ! Check on the input field
      if(.not.Lttc%status) stop "read_U_vasp_full: Lattice not properly initialized."
      if(.not.Umats%status) stop "read_U_vasp_full: BosonicField not properly initialized."
      if(Umats%Npoints.ne.1) stop "read_U_vasp_full: Number of matsubara points in Umats is supposed to be equal to 1."
      if(Umats%Nbp.ne.(Lttc%Norb**2)) stop "read_U_vasp_full: Umats has different orbital dimension with respect to Lttc."
      !
      !recover the vectors in real space
      if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
      call get_Ruc(Ruc)
      !
      Norb = Lttc%Norb
      Nbp = Umats%Nbp
      !
      allocate(ReadLine(Norb));ReadLine=0d0
      allocate(U_iijj(Norb,Norb));U_iijj=0d0
      allocate(U_ijji(Norb,Norb));U_ijji=0d0
      allocate(U_ijij(Norb,Norb));U_ijij=0d0
      !
      file = reg(pathINPUT)//"Vasp/U_ijkl.DAT"
      call inquireFile(reg(file),filexists,verb=verbose)
      unit = free_unit()
      open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U_vasp_full: wrong index in U_iijj."
         U_iijj(iorb,:) = ReadLine
      enddo
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U_vasp_full: wrong index in U_ijji."
         U_ijji(iorb,:) = ReadLine
      enddo
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U_vasp_full: wrong index in U_ijij."
         U_ijij(iorb,:) = ReadLine
      enddo
      deallocate(ReadLine)
      !
      if(sum(abs(U_ijji-U_ijij)).gt.eps) write(*,"(A)") "     Warning(read_U_vasp_full): non-invariant Hund's coupling."
      !
      allocate(Ur(Nbp,Nbp,Nwig));Ur=czero
      do isite=1,Lttc%Nsite
         do jsite=1,Lttc%Nsite
            !
            do i=1,SiteOrbs(isite)%Norb
               do j=1,SiteOrbs(jsite)%Norb
                  !
                  iorb = SiteOrbs(isite)%Orbs(i)
                  jorb = SiteOrbs(jsite)%Orbs(j)
                  !
                  !U_iijj
                  call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                  Ur(ib1,ib2,wig0) = dcmplx(U_iijj(iorb,jorb),0d0)
                  !
                  !U_ijji
                  call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                  Ur(ib1,ib2,wig0) = dcmplx(U_ijji(iorb,jorb),0d0)
                  !
                  !U_ijij
                  call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                  Ur(ib1,ib2,wig0) = dcmplx(U_ijij(iorb,jorb),0d0)
                  !
               enddo
            enddo
            !
         enddo
      enddo
      deallocate(U_iijj,U_ijji,U_ijij)
      !
      !build long range part of the interaction
      if(Lttc%Nsite.gt.1)then
         !
         !Get all the possible positions
         call cpu_time(start)
         allocate(Rvec(3));Rvec=0d0
         allocate(Rsorted(Nwig*Lttc%Nsite*Lttc%Nsite,4));Rsorted=0d0
         iR=0
         do iwig=1,Nwig
            do jsite=1,Lttc%Nsite
               do isite=1,Lttc%Nsite
                  !
                  Rvec = Rvecwig(:,iwig) + Ruc(:,jsite) - Ruc(:,isite)
                  Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                  !
                  iR = iR +1
                  !
                  Rsorted(iR,1) = Rdist
                  Rsorted(iR,2) = iwig
                  Rsorted(iR,3) = jsite
                  Rsorted(iR,4) = isite
                  !
               enddo
            enddo
         enddo
         deallocate(Rvec)
         !
         !Sorting the positions according to distance
         allocate(Rorder(Nwig*Lttc%Nsite*Lttc%Nsite));Rorder=0
         allocate(Rsorted_bkp(Nwig*Lttc%Nsite*Lttc%Nsite,4));Rsorted_bkp=0d0
         Rsorted_bkp = Rsorted
         call sort_array(Rsorted(:,1),Rorder)
         Rsorted=0d0
         do iR=1,Nwig*Lttc%Nsite*Lttc%Nsite
            Rsorted(iR,:) = Rsorted_bkp(Rorder(iR),:)
         enddo
         deallocate(Rsorted_bkp,Rorder)
         !
         !Regroup according to distance. The list contains the indexes of all the positions with a given distance
         call get_pattern(Dist,Rsorted(:,1),1e4*eps,listDim=DistList,IncludeSingle=.true.)
         !
         !all the possible ranges: iD=1 is the local(already done), iD=2 is the nearest neighbor (already done)
         do iD=3,size(Dist,dim=1)
            !
            !beyond becomes undefined for more than 2 sites in the uc
            if(iD.gt.3) exit
            !
            !all the indexes within that range
            do iR=1,DistList(iD)
               !
               !retrieve indexes from sorted list
               iwig = Rsorted(Dist(iD,iR),2)
               jsite = Rsorted(Dist(iD,iR),3)
               isite = Rsorted(Dist(iD,iR),4)
               !
               !sites are different because iD=2
               do i=1,SiteOrbs(isite)%Norb
                  do j=1,SiteOrbs(jsite)%Norb
                     !
                     iorb = SiteOrbs(isite)%Orbs(i)
                     jorb = SiteOrbs(jsite)%Orbs(j)
                     !
                     !U_iijj
                     call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                     Ur(ib1,ib2,iwig) = Ur(ib1,ib2,wig0)
                     !
                     !U_ijji
                     call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                     Ur(ib1,ib2,iwig) = Ur(ib1,ib2,wig0)
                     !
                     !U_ijij
                     call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                     Ur(ib1,ib2,iwig) = Ur(ib1,ib2,wig0)
                     !
                  enddo
               enddo
               !
            enddo
            !
         enddo
         deallocate(Rsorted,Dist,DistList)
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of U(R) cpu timing:", finish-start
         !
         !FT to K-space
         call cpu_time(start)
         allocate(Uk(Nbp,Nbp,Lttc%Nkpt));Uk=czero
         if(Lttc%Nkpt.gt.1)then
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Ur,Uk)
         else
            Uk(:,:,1) = Ur(:,:,wig0)
         endif
         deallocate(Ur)
         where(abs((Uk))<eps) Uk=czero
         call cpu_time(finish)
         write(*,"(A,F)") "     U(R) --> U(K) cpu timing:", finish-start
         !
         !print along path
         if(reg(structure).ne."None")then
            call interpolate2Path(Lttc,Nkpt_path,"Uk",pathOUTPUT=reg(pathINPUT),store=.false.,skipAkw=.true.,data_in=Uk)
         endif
         !
         !fill in the output
         Umats%screened(:,:,1,:) = Uk
         if(allocated(Umats%bare)) Umats%bare = Uk
         !
      else
         !
         !fill in the output
         do ik=1,Lttc%Nkpt
            Umats%screened(:,:,1,ik) = Ur(:,:,wig0)
            if(allocated(Umats%bare)) Umats%bare(:,:,ik) = Ur(:,:,wig0)
         enddo
         deallocate(Ur)
         !
      endif
      !
      call BosonicKsum(Umats)
      call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
      !
      ! Print the nn non local interaction
      if(Lttc%Nsite.gt.1)then
         !
         Nprint = size(RealPrint,dim=2)
         allocate(URnn(int(sqrt(dble(Umats%Nbp))),int(sqrt(dble(Umats%Nbp)))));URnn=czero
         allocate(UR(Umats%Nbp,Umats%Nbp,Nprint));UR=czero
         call wannier_K2R_NN(RealPrint,Lttc%Nkpt3,Lttc%kpt,Umats%screened(:,:,1,:),UR)
         where(abs((UR))<eps) UR=czero
         do iprint=1,Nprint
            call dump_Matrix(UR(:,:,iprint),reg(trim(pathINPUTtr)),"U_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
            call product2NN(UR(:,:,iprint),URnn)
            call dump_Matrix(URnn,reg(trim(pathINPUTtr)),"Unn_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
         enddo
         deallocate(UR,URnn)
         !
      endif
      !
   end subroutine read_U_vasp_full
   !
   subroutine read_U_vasp_Uloc0(Umat,Lttc)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use input_vars, only : pathINPUT, SiteOrbs
      implicit none
      !
      complex(8),allocatable,intent(inout)  :: Umat(:,:)
      type(Lattice),intent(inout)           :: Lttc
      !
      integer                               :: i,j,isite,jsite
      integer                               :: iorb,jorb,ib1,ib2
      integer                               :: Norb,Nbp
      integer                               :: unit,idum
      real(8),allocatable                   :: Ruc(:,:),ReadLine(:)
      real(8),allocatable                   :: U_iijj(:,:),U_ijji(:,:),U_ijij(:,:)
      character(len=255)                    :: file
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_U_vasp_Uloc0"
      !
      !
      ! Check on the input field
      Nbp = size(Umat,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      !
      call assert_shape(Umat,[Nbp,Nbp],"read_U_vasp_Uloc0","Umat")
      !
      !recover the vectors in real space
      if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
      call get_Ruc(Ruc)
      !
      allocate(ReadLine(Norb));ReadLine=0d0
      allocate(U_iijj(Norb,Norb));U_iijj=0d0
      allocate(U_ijji(Norb,Norb));U_ijji=0d0
      allocate(U_ijij(Norb,Norb));U_ijij=0d0
      !
      file = reg(pathINPUT)//"Vasp/U_ijkl.DAT"
      call inquireFile(reg(file),filexists,verb=verbose)
      unit = free_unit()
      open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U_vasp_Uloc0: wrong index in U_iijj."
         U_iijj(iorb,:) = ReadLine
      enddo
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U_vasp_Uloc0: wrong index in U_ijji."
         U_ijji(iorb,:) = ReadLine
      enddo
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U_vasp_Uloc0: wrong index in U_ijij."
         U_ijij(iorb,:) = ReadLine
      enddo
      deallocate(ReadLine)
      !
      if(sum(abs(U_ijji-U_ijij)).gt.eps) write(*,"(A)") "     Warning(read_U_vasp_Uloc0): non-invariant Hund's coupling."
      !
      Umat=czero
      do isite=1,Lttc%Nsite
         do jsite=1,Lttc%Nsite
            !
            do i=1,SiteOrbs(isite)%Norb
               do j=1,SiteOrbs(jsite)%Norb
                  !
                  iorb = SiteOrbs(isite)%Orbs(i)
                  jorb = SiteOrbs(jsite)%Orbs(j)
                  !
                  !U_iijj
                  call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                  Umat(ib1,ib2) = dcmplx(U_iijj(iorb,jorb),0d0)
                  !
                  !U_ijji
                  call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                  Umat(ib1,ib2) = dcmplx(U_ijji(iorb,jorb),0d0)
                  !
                  !U_ijij
                  call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                  Umat(ib1,ib2) = dcmplx(U_ijij(iorb,jorb),0d0)
                  !
               enddo
            enddo
            !
         enddo
      enddo
      deallocate(U_iijj,U_ijji,U_ijij)
      !
   end subroutine read_U_vasp_Uloc0


   !---------------------------------------------------------------------------!
   !PURPOSE: Read local interaction in product basis format.
   !---------------------------------------------------------------------------!
   subroutine read_U(Umat)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use input_vars, only : pathINPUT
      implicit none
      !
      complex(8),allocatable,intent(inout)  :: Umat(:,:)
      !
      integer                               :: iorb,jorb,ib1,ib2
      integer                               :: Norb,Nbp
      integer                               :: unit,idum
      real(8),allocatable                   :: ReadLine(:)
      real(8),allocatable                   :: U_iijj(:,:),U_ijji(:,:),U_ijij(:,:)
      character(len=255)                    :: file
      logical                               :: filexists
      !
      !
      if(verbose)write(*,"(A)") "---- read_U"
      !
      !
      ! Check on the input field
      Nbp = size(Umat,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      !
      call assert_shape(Umat,[Nbp,Nbp],"read_U","Umat")
      !
      allocate(ReadLine(Norb));ReadLine=0d0
      allocate(U_iijj(Norb,Norb));U_iijj=0d0
      allocate(U_ijji(Norb,Norb));U_ijji=0d0
      allocate(U_ijij(Norb,Norb));U_ijij=0d0
      !
      file = reg(pathINPUT)//"U_ijkl.DAT"
      call inquireFile(reg(file),filexists,verb=verbose)
      unit = free_unit()
      open(unit,file=reg(file),form="formatted",status="old",position="rewind",action="read")
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U: wrong index in U_iijj."
         U_iijj(iorb,:) = ReadLine
      enddo
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U: wrong index in U_ijji."
         U_ijji(iorb,:) = ReadLine
      enddo
      call skip_header(unit,3)
      do iorb=1,Norb
         ReadLine=0d0
         read(unit,*) idum,ReadLine
         if(idum.ne.iorb) stop "read_U: wrong index in U_ijij."
         U_ijij(iorb,:) = ReadLine
      enddo
      deallocate(ReadLine)
      !
      if(sum(abs(U_ijji-U_ijij)).gt.eps) write(*,"(A)") "     Warning(read_U): non-invariant Hund's coupling."
      !
      Umat=czero
      do iorb=1,Norb
         do jorb=1,Norb
            !
            !U_iijj
            call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
            Umat(ib1,ib2) = dcmplx(U_iijj(iorb,jorb),0d0)
            !
            !U_ijji
            call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
            Umat(ib1,ib2) = dcmplx(U_ijji(iorb,jorb),0d0)
            !
            !U_ijij
            call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
            Umat(ib1,ib2) = dcmplx(U_ijij(iorb,jorb),0d0)
            !
         enddo
      enddo
      deallocate(U_iijj,U_ijji,U_ijij)
      !
   end subroutine read_U


   !---------------------------------------------------------------------------!
   !PURPOSE: Check how the AC alters the bare and screened values
   !---------------------------------------------------------------------------!
   subroutine checkAnalyticContinuation(Umats,Ureal)
      !
      use parameters
      use utils_misc
      use input_vars, only : pathINPUTtr
      implicit none
      !
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Ureal
      character(len=256)                    :: path
      integer                               :: iq,ib1,ib2
      integer                               :: unit,Nbp,Nmats,Nreal,Nkpt
      real(8)                               :: ReErr,ImErr,thresh=1e-4
      real(8),allocatable                   :: ReErrMat(:,:),ImErrMat(:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- checkAnalyticContinuation"
      if(Umats%Nbp.ne.Ureal%Nbp) stop "checkAnalyticContinuation: Umats%Nbp.ne.Ureal%Nbp"
      Nbp = Umats%Nbp
      Nmats = Umats%Npoints
      Nreal = Ureal%Npoints
      allocate(ReErrMat(Nbp,Nbp));ReErrMat=0d0
      allocate(ImErrMat(Nbp,Nbp));ImErrMat=0d0
      !
      !
      ! Check the difference betqween bare values induced by thecutoff in the matsubara frequency
      unit = free_unit()
      path = reg(pathINPUTtr)//"ACcutoffError.DAT"
      open(unit=unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(A,1I5,A,1E14.7)")"Difference between Umats_bare value and last screened frequency: ",Nmats," thresh:",thresh
      !
      ReErrMat=0d0;ImErrMat=0d0
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            ReErr = abs(dreal(Umats%bare_local(ib1,ib2)) - dreal(Umats%screened_local(ib1,ib2,Nmats)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(dimag(Umats%bare_local(ib1,ib2)) - dimag(Umats%screened_local(ib1,ib2,Nmats)))
            if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
            !
         enddo
      enddo
      write(unit,"(A)")"Real part - local projection"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      !
      if((Umats%Nkpt.eq.Ureal%Nkpt).and.(Ureal%Nkpt.gt.0))then
         Nkpt = Umats%Nkpt
         do iq=1,Nkpt
            !
            ReErrMat=0d0;ImErrMat=0d0
            do ib1=1,Nbp
               do ib2=1,Nbp
                  !
                  ReErr = abs(dreal(Umats%bare(ib1,ib2,iq)) - dreal(Umats%screened(ib1,ib2,Nmats,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(dimag(Umats%bare(ib1,ib2,iq)) - dimag(Umats%screened(ib1,ib2,Nmats,iq)))
                  if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
                  !
               enddo
            enddo
            !
            write(unit,"(A,1I5)")"Real part - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            !
         enddo !iq
      endif
      close(unit)
      !
      !
      ! Check that the screened and bare values of Umats and Ureal are close enough
      unit = free_unit()
      path = reg(pathINPUTtr)//"ACcheck.DAT"
      open(unit=unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(A,1E14.7)")"Difference between asymptotic behaviour of Ureal and Umats. Thresh:",thresh
      !
      ReErrMat=0d0;ImErrMat=0d0
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            ReErr = abs(dreal(Umats%screened_local(ib1,ib2,1)) - dreal(Ureal%screened_local(ib1,ib2,1)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(dimag(Umats%screened_local(ib1,ib2,1)) - dimag(Ureal%screened_local(ib1,ib2,1)))
            if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
            !
         enddo
      enddo
      write(unit,"(A)")"Real part - local projection - screened limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection - screened limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      !
      ReErrMat=0d0;ImErrMat=0d0
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            ReErr = abs(dreal(Umats%screened_local(ib1,ib2,Nmats)) - dreal(Ureal%screened_local(ib1,ib2,Nreal)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(dimag(Umats%screened_local(ib1,ib2,Nmats)) - dimag(Ureal%screened_local(ib1,ib2,Nreal)))
            if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
            !
         enddo
      enddo
      write(unit,"(A)")"Real part - local projection - bare limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection - bare limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      !
      if((Umats%Nkpt.eq.Ureal%Nkpt).and.(Ureal%Nkpt.gt.0))then
         Nkpt = Umats%Nkpt
         !
         do iq=1,Nkpt
            !
            ReErrMat=0d0;ImErrMat=0d0
            do ib1=1,Nbp
               do ib2=1,Nbp
                  !
                  ReErr = abs(dreal(Umats%screened(ib1,ib2,1,iq)) - dreal(Ureal%screened(ib1,ib2,1,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(dimag(Umats%screened(ib1,ib2,1,iq)) - dimag(Ureal%screened(ib1,ib2,1,iq)))
                  if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
                  !
               enddo
            enddo
            !
            write(unit,"(A,1I5)")"Real part - screened limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - screened limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            !
            ReErrMat=0d0;ImErrMat=0d0
            do ib1=1,Nbp
               do ib2=1,Nbp
                  !
                  ReErr = abs(dreal(Umats%screened(ib1,ib2,Nmats,iq)) - dreal(Ureal%screened(ib1,ib2,Nreal,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(dimag(Umats%screened(ib1,ib2,Nmats,iq)) - dimag(Ureal%screened(ib1,ib2,Nreal,iq)))
                  if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
                  !
               enddo
            enddo
            !
            write(unit,"(A,1I5)")"Real part - bare limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - bare limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            !
         enddo !iq
      endif
      close(unit)
      !
   end subroutine checkAnalyticContinuation


   !---------------------------------------------------------------------------!
   !PURPOSE: Create the static interaction tensor from user-given parameters
   !TO DO: extend this subroutine to multi-site version with different local
   !       interactions for each site
   !---------------------------------------------------------------------------!
   subroutine build_Umat_singlParam(Umat,Uaa,Uab,Jh)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : LocalOrbs
      implicit none
      !
      complex(8),intent(inout)              :: Umat(:,:)
      real(8),intent(in)                    :: Uaa,Uab,Jh
      !
      integer                               :: Nbp,Norb,isite
      integer                               :: indx,jndx,kndx,lndx,ib1ndx,ib2ndx
      integer                               :: i,j,k,l,ib1,ib2
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- build_Umat_singlParam"
      !
      !
      ! Check on the input matrices
      Nbp = size(Umat,dim=1)
      call assert_shape(Umat,[Nbp,Nbp],"build_Umat_singlParam","Umat")
      !
      Norb = int(sqrt(dble(Nbp)))
      !
      do isite=1,size(LocalOrbs)
         call init_Uelements(LocalOrbs(isite)%Norb,PhysicalUelements)
         do indx=1,LocalOrbs(isite)%Norb
            do jndx=1,LocalOrbs(isite)%Norb
               do kndx=1,LocalOrbs(isite)%Norb
                  do lndx=1,LocalOrbs(isite)%Norb
                     !
                     call F2Bindex(LocalOrbs(isite)%Norb,[indx,jndx],[kndx,lndx],ib1ndx,ib2ndx)
                     !
                     i = LocalOrbs(isite)%Orbs(indx)
                     j = LocalOrbs(isite)%Orbs(jndx)
                     k = LocalOrbs(isite)%Orbs(kndx)
                     l = LocalOrbs(isite)%Orbs(lndx)
                     call F2Bindex(Norb,[i,j],[k,l],ib1,ib2)
                     !
                     if(PhysicalUelements%Full_Uaa(ib1ndx,ib2ndx)) Umat(ib1,ib2) = dcmplx(Uaa,0d0)
                     if(PhysicalUelements%Full_Uab(ib1ndx,ib2ndx)) Umat(ib1,ib2) = dcmplx(Uab,0d0)
                     if(PhysicalUelements%Full_J(ib1ndx,ib2ndx))   Umat(ib1,ib2) = dcmplx(Jh,0d0)
                     !
                  enddo
               enddo
            enddo
         enddo
      enddo
      !
      !
      !call init_Uelements(Norb,PhysicalUelements)
      !do ib1=1,Nbp
      !   do ib2=1,Nbp
      !      !
      !      if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uaa,0d0)
      !      if(PhysicalUelements%Full_Uab(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab,0d0)
      !      if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umat(ib1,ib2) = dcmplx(J,0d0)
      !      if(PhysicalUelements%Full_Jph(ib1,ib2)) Umat(ib1,ib2) = dcmplx(J,0d0)
      !      !
      !   enddo
      !enddo
      !
   end subroutine build_Umat_singlParam
   !
   subroutine build_Umat_multiParam(Umat,Uaa,Uab,J)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      implicit none
      !
      complex(8),intent(inout)              :: Umat(:,:)
      real(8),intent(in)                    :: Uaa(:),Uab(:,:),J(:,:)
      !
      integer                               :: Nbp,Norb
      integer                               :: ib1,ib2,iorb,jorb
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- build_Umat_multiParam"
      !
      !
      ! Check on the input matrices
      Nbp = size(Umat,dim=1)
      call assert_shape(Umat,[Nbp,Nbp],"build_Umat_singlParam","Umat")
      !
      Norb = int(sqrt(dble(Nbp)))
      call init_Uelements(Norb,PhysicalUelements)
      !
      call assert_shape(Uaa,[Norb],"build_Umat_multiParam","Uaa")
      call assert_shape(Uab,[Norb,Norb],"build_Umat_multiParam","Uab")
      call assert_shape(J,[Norb,Norb],"build_Umat_multiParam","J")
      !
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            iorb = PhysicalUelements%Flav_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Flav_Map(ib1,ib2,3)
            !
            if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uaa(iorb),0d0)
            if(PhysicalUelements%Full_Uab(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umat(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jph(ib1,ib2)) Umat(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
            !
         enddo
      enddo
      !
   end subroutine build_Umat_multiParam
   !
   subroutine build_Umat_singlParam_Flav(Umat,Uaa,Uab,J)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : Hetero
      implicit none
      !
      complex(8),intent(inout)              :: Umat(:,:)
      real(8),intent(in)                    :: Uaa,Uab,J
      !
      integer                               :: Nbp,Norb,Nflavor
      integer                               :: ib1,ib2
      integer                               :: shift,isite,Nsite
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- build_Umat_singlParam_Flav"
      !
      !
      ! Check on the input matrices
      Nbp = size(Umat,dim=1)
      if((Nspin.eq.2).and.(mod(Nbp,2).ne.0.0)) stop "build_Umat_singlParam_Flav: Wrong matrix dimension."
      call assert_shape(Umat,[Nbp,Nbp],"build_Umat_singlParam_Flav","Umat")
      !
      Norb = Nbp/Nspin
      if(Hetero%status) Norb = Hetero%Norb
      call init_Uelements(Norb,PhysicalUelements)
      !
      Nflavor = Norb*Nspin
      !
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            if(PhysicalUelements%Flav_Uloc(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uaa,0d0)
            if(PhysicalUelements%Flav_U1st(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab,0d0)
            if(PhysicalUelements%Flav_U2nd(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab-J,0d0)
            !
         enddo
      enddo
      !
      !In-plane Interaction
      if(Hetero%status)then
         Nsite = Hetero%Explicit(2)-Hetero%Explicit(1)+1
         do isite=2,Nsite
            shift = (isite-1)*Nflavor
            Umat(1+shift:Nflavor+shift,1+shift:Nflavor+shift) = Umat(1:Nflavor,1:Nflavor)
         enddo
      endif
      !
   end subroutine build_Umat_singlParam_Flav
   !
   subroutine build_Umat_multiParam_Flav(Umat,Uaa,Uab,J)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : Hetero
      implicit none
      !
      complex(8),intent(inout)              :: Umat(:,:)
      real(8),intent(in)                    :: Uaa(:),Uab(:,:),J(:,:)
      !
      integer                               :: Nbp,Norb,Nflavor
      integer                               :: ib1,ib2,iorb,jorb
      integer                               :: shift,isite,Nsite
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- build_Umat_multiParam_Flav"
      !
      !
      ! Check on the input matrices
      Nbp = size(Umat,dim=1)
      if((Nspin.eq.2).and.(mod(Nbp,2).ne.0.0)) stop "build_Umat_multiParam_Flav: Wrong matrix dimension."
      call assert_shape(Umat,[Nbp,Nbp],"build_Umat_multiParam_Flav","Umat")
      !
      Norb = Nbp/Nspin
      if(Hetero%status) Norb = Hetero%Norb
      call init_Uelements(Norb,PhysicalUelements)
      call assert_shape(Uaa,[Norb],"build_Umat_multiParam_Flav","Uaa")
      call assert_shape(Uab,[Norb,Norb],"build_Umat_multiParam_Flav","Uab")
      call assert_shape(J,[Norb,Norb],"build_Umat_multiParam_Flav","J")
      !
      Nflavor = Norb*Nspin
      !
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            iorb = PhysicalUelements%Flav_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Flav_Map(ib1,ib2,2)
            !
            if(PhysicalUelements%Flav_Uloc(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uaa(iorb),0d0)
            if(PhysicalUelements%Flav_U1st(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab(iorb,jorb),0d0)
            if(PhysicalUelements%Flav_U2nd(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab(iorb,jorb)-J(iorb,jorb),0d0)
            !
         enddo
      enddo
      !
      !In-plane Interaction
      if(Hetero%status)then
         Nsite = Hetero%Explicit(2)-Hetero%Explicit(1)+1
         do isite=2,Nsite
            shift = (isite-1)*Nflavor
            Umat(1+shift:Nflavor+shift,1+shift:Nflavor+shift) = Umat(1:Nflavor,1:Nflavor)
         enddo
      endif
      !
   end subroutine build_Umat_multiParam_Flav


   !---------------------------------------------------------------------------!
   !PURPOSE: Create the frequency dependent interaction tensor from user-given
   !         phononic modes.
   !---------------------------------------------------------------------------!
   subroutine build_Uret_singlParam_ph(Umats,Uaa,Uab,J,g_eph,wo_eph,Hetero,LocalOnly,ScreenAll)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : Nreal, wrealMax, pathINPUTtr
      implicit none
      !
      type(BosonicField),intent(inout),target :: Umats
      real(8),intent(in)                    :: Uaa,Uab,J
      complex(8),intent(in)                 :: g_eph(:)
      real(8),intent(in)                    :: wo_eph(:)
      type(Heterostructures),intent(in)     :: Hetero
      logical,intent(in),optional           :: LocalOnly
      logical,intent(in),optional           :: ScreenAll
      !
      integer                               :: Nbp,Norb,Nph
      integer                               :: ib1,ib2,ik
      integer                               :: iw,iw1,iw2,iph,iwp
      real(8)                               :: g,eta_g
      real(8)                               :: dw,RealU,ImagU
      real(8),allocatable                   :: wreal(:),wmats(:)
      complex(8),allocatable                :: D(:,:),Utmp(:,:)
      type(BosonicField)                    :: Ureal
      type(BosonicField),target             :: Umats_imp
      type(BosonicField),pointer            :: Umats_ptr
      type(physicalU)                       :: PhysicalUelements
      logical                               :: LocalOnly_,ScreenAll_
      logical                               :: Screen
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- build_Uret_singlParam_ph"
      !
      !
      ! Check on the input field
      if(.not.Umats%status) stop "build_Uret_singlParam_ph: BosonicField not properly initialized."
      !
      LocalOnly_=.true.
      if(present(LocalOnly))LocalOnly_=LocalOnly
      if(LocalOnly_.and.(Umats%Nkpt.ne.0)) stop "build_Uret_singlParam_ph: Umats k dependent attributes are supposed to be unallocated."
      !
      ScreenAll_=.false.
      if(present(ScreenAll))ScreenAll_=ScreenAll
      !
      Nph = size(g_eph)
      if(size(g_eph).ne.size(wo_eph)) stop "build_Uret_singlParam_ph: Phonon sizes does not match."
      !
      allocate(wmats(Umats%Npoints));wmats=0d0
      wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      allocate(wreal(Nreal));wreal=0d0
      wreal = linspace(0d0,+wrealMax,Nreal)
      !
      call clear_attributes(Umats)
      !
      Nbp = Umats%Nbp
      Norb = int(sqrt(dble(Nbp)))
      if(Hetero%status)then
         Norb = Hetero%Norb
         Nbp = Norb**2
         call AllocateBosonicField(Umats_imp,Norb,Umats%Npoints,Umats%iq_gamma,Nkpt=Umats%Nkpt,Beta=Umats%Beta)
         Umats_ptr => Umats_imp
      else
         Umats_ptr => Umats
      endif
      call AllocateBosonicField(Ureal,Norb,Nreal,0)
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      !setting the bare values
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uaa,0d0)
            if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uab,0d0)
            if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J,0d0)
            if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J,0d0)
            !
         enddo
      enddo
      !
      !setting the phonons
      dw = abs(wreal(3)-wreal(2))
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            Screen = ScreenAll_ .or. PhysicalUelements%Full_Uaa(ib1,ib2) .or. PhysicalUelements%Full_Uab(ib1,ib2)
            !
            do iph=1,Nph
               !
               !phonon frequency index
               iwp = minloc(abs(wreal-wo_eph(iph)),dim=1)
               !electron-phonon coupling
               g = dreal(g_eph(iph))
               !electron-phonon broadening
               eta_g = dimag(g_eph(iph))
               !
               do iw=1,Nreal
                  !
                  RealU = 2*(g**2)*wo_eph(iph) / ( (wreal(iw)**2) - (wo_eph(iph)**2) )
                  ImagU=0d0
                  !
                  !Dirac delta function
                  if(iw.eq.iwp) ImagU = -pi*(g**2)/dw
                  !Lorentzian
                  if(eta_g.ne.0d0) ImagU = -pi*(g**2) * ( eta_g**2 / ( (wreal(iw)-wreal(iwp))**2 + eta_g**2 ) )/dw
                  !
                  Ureal%screened_local(ib1,ib2,iw) = Umats_ptr%bare_local(ib1,ib2)
                  if(Screen) Ureal%screened_local(ib1,ib2,iw) = Ureal%screened_local(ib1,ib2,iw) + dcmplx(RealU,ImagU)
                  !
               enddo
            enddo
            !
         enddo
      enddo
      call dump_BosonicField(Ureal,reg(pathINPUTtr),"Uloc_real.DAT",axis=wreal)
      !
      ! Allocate the temporary quantities needed by the Analytical continuation
      allocate(Utmp(Nbp,Nbp));Utmp=czero
      allocate(D(Nbp,Nbp));D=czero
      !
      ! Analytical continuation of the local component to imag axis using spectral rep
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,wmats,wreal,Nreal,Ureal,Umats_ptr),&
      !$OMP PRIVATE(ib1,ib2,iw1,iw2,D,Utmp)
      !$OMP DO
      do iw1=1,Umats_ptr%Npoints
         !
         Utmp=czero
         !
         do iw2=1,Nreal
            !
            if((wmats(iw1).eq.0d0).and.(wreal(iw2).eq.0d0))cycle
            !
            D = -dimag( Ureal%screened_local(:,:,iw2)*abs(wreal(3)-wreal(2)) )/pi
            Utmp = Utmp +  D/( dcmplx(0.d0,wmats(iw1))-wreal(iw2) ) - D/( dcmplx(0.d0,wmats(iw1))+wreal(iw2) )
            !
         enddo
         !
         Umats_ptr%screened_local(:,:,iw1) = Utmp + Umats_ptr%bare_local
         !
      enddo !iw1
      !
      !$OMP END DO
      !$OMP END PARALLEL
      call cpu_time(finish)
      deallocate(D,Utmp,wmats,wreal)
      call DeallocateBosonicField(Ureal)
      write(*,"(A,F)") "     Ue-ph(w) --> Ue-ph(iw) cpu timing:", finish-start
      if(.not.ScreenAll_)write(*,"(A)") "     Screening only Kanamori NaNb terms."
      !
      if(.not.LocalOnly_)then
         write(*,"(A)") "     Filling the K-dependent attributes."
         do ik=1,Umats%Nkpt
            Umats_ptr%bare(:,:,ik) = Umats_ptr%bare_local
            Umats_ptr%screened(:,:,:,ik) = Umats_ptr%screened_local
         enddo
      endif
      !
      if(Hetero%status)then
         call Expand2Nsite(Umats,Umats_ptr,Hetero%Nlayer)
         call DeallocateBosonicField(Umats_imp)
      endif
      nullify(Umats_ptr)
      !
      call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
      !
   end subroutine build_Uret_singlParam_ph
   !
   subroutine build_Uret_multiParam_ph(Umats,Uaa,Uab,J,g_eph,wo_eph,Hetero,LocalOnly,ScreenAll)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : Nreal, wrealMax, pathINPUTtr
      implicit none
      !
      type(BosonicField),intent(inout),target :: Umats
      real(8),intent(in)                    :: Uaa(:),Uab(:,:),J(:,:)
      complex(8),intent(in)                 :: g_eph(:)
      real(8),intent(in)                    :: wo_eph(:)
      type(Heterostructures),intent(in)     :: Hetero
      logical,intent(in),optional           :: LocalOnly
      logical,intent(in),optional           :: ScreenAll
      !
      integer                               :: Nbp,Norb,Nph
      integer                               :: ib1,ib2,iorb,jorb,ik
      integer                               :: iw,iw1,iw2,iph,iwp
      real(8)                               :: g,eta_g
      real(8)                               :: dw,RealU,ImagU
      real(8),allocatable                   :: wreal(:),wmats(:)
      complex(8),allocatable                :: D(:,:),Utmp(:,:)
      type(BosonicField)                    :: Ureal
      type(BosonicField),target             :: Umats_imp
      type(BosonicField),pointer            :: Umats_ptr
      type(physicalU)                       :: PhysicalUelements
      logical                               :: LocalOnly_,ScreenAll_
      logical                               :: Screen
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- build_Uret_multiParam_ph"
      !
      !
      ! Check on the input field
      if(.not.Umats%status) stop "build_Uret_multiParam_ph: BosonicField not properly initialized."
      !
      LocalOnly_=.true.
      if(present(LocalOnly))LocalOnly_=LocalOnly
      if(LocalOnly_.and.(Umats%Nkpt.ne.0)) stop "build_Uret_multiParam_ph: Umats k dependent attributes are supposed to be unallocated."
      !
      ScreenAll_=.false.
      if(present(ScreenAll))ScreenAll_=ScreenAll
      !
      Nph = size(g_eph)
      if(size(g_eph).ne.size(wo_eph)) stop "build_Uret_multiParam_ph: Phonon sizes does not match."
      !
      allocate(wmats(Umats%Npoints));wmats=0d0
      wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      allocate(wreal(Nreal));wreal=0d0
      wreal = linspace(0d0,+wrealMax,Nreal)
      !
      call clear_attributes(Umats)
      !
      Nbp = Umats%Nbp
      Norb = int(sqrt(dble(Nbp)))
      if(Hetero%status)then
         Norb = Hetero%Norb
         Nbp = Norb**2
         call AllocateBosonicField(Umats_imp,Norb,Umats%Npoints,Umats%iq_gamma,Nkpt=Umats%Nkpt,Beta=Umats%Beta)
         Umats_ptr => Umats_imp
      else
         Umats_ptr => Umats
      endif
      call AllocateBosonicField(Ureal,Norb,Nreal,0)
      !
      call assert_shape(Uaa,[Norb],"build_Uret_multiParam_ph","Uaa")
      call assert_shape(Uab,[Norb,Norb],"build_Uret_multiParam_ph","Uab")
      call assert_shape(J,[Norb,Norb],"build_Uret_multiParam_ph","J")
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      !setting the bare values
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            iorb = PhysicalUelements%Full_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Full_Map(ib1,ib2,2)
            !
            if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uaa(iorb),0d0)
            if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uab(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
            !
         enddo
      enddo
      !
      !setting the phonons
      dw = abs(wreal(3)-wreal(2))
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            Screen = ScreenAll_ .or. PhysicalUelements%Full_Uaa(ib1,ib2) .or. PhysicalUelements%Full_Uab(ib1,ib2)
            !
            do iph=1,Nph
               !
               !phonon frequency index
               iwp = minloc(abs(wreal-wo_eph(iph)),dim=1)
               !electron-phonon coupling
               g = dreal(g_eph(iph))
               !electron-phonon broadening
               eta_g = dimag(g_eph(iph))
               !
               do iw=1,Nreal
                  !
                  RealU = 2*(g**2)*wo_eph(iph) / ( (wreal(iw)**2) - (wo_eph(iph)**2) )
                  ImagU=0d0
                  !
                  !Dirac delta function
                  if(iw.eq.iwp) ImagU = -pi*(g**2)/dw
                  !Lorentzian
                  if(eta_g.ne.0d0) ImagU = -pi*(g**2) * ( eta_g**2 / ( (wreal(iw)-wreal(iwp))**2 + eta_g**2 ) )/dw
                  !
                  Ureal%screened_local(ib1,ib2,iw) = Umats_ptr%bare_local(ib1,ib2)
                  if(Screen) Ureal%screened_local(ib1,ib2,iw) = Ureal%screened_local(ib1,ib2,iw) + dcmplx(RealU,ImagU)
                  !
               enddo
            enddo
            !
         enddo
      enddo
      call dump_BosonicField(Ureal,reg(pathINPUTtr),"Uloc_real.DAT",axis=wreal)
      !
      ! Allocate the temporary quantities needed by the Analytical continuation
      allocate(Utmp(Nbp,Nbp));Utmp=czero
      allocate(D(Nbp,Nbp));D=czero
      !
      ! Analytical continuation of the local component to imag axis using spectral rep
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,wmats,wreal,Nreal,Ureal,Umats_ptr),&
      !$OMP PRIVATE(ib1,ib2,iw1,iw2,D,Utmp)
      !$OMP DO
      do iw1=1,Umats_ptr%Npoints
         !
         Utmp=czero
         !
         do iw2=1,Nreal
            !
            if((wmats(iw1).eq.0d0).and.(wreal(iw2).eq.0d0))cycle
            !
            D = -dimag( Ureal%screened_local(:,:,iw2)*abs(wreal(3)-wreal(2)) )/pi
            Utmp = Utmp +  D/( dcmplx(0.d0,wmats(iw1))-wreal(iw2) ) - D/( dcmplx(0.d0,wmats(iw1))+wreal(iw2) )
            !
         enddo
         !
         Umats_ptr%screened_local(:,:,iw1) = Utmp + Umats_ptr%bare_local
         !
      enddo !iw1
      !
      !$OMP END DO
      !$OMP END PARALLEL
      call cpu_time(finish)
      deallocate(D,Utmp,wmats,wreal)
      call DeallocateBosonicField(Ureal)
      write(*,"(A,F)") "     Ue-ph(w) --> Ue-ph(iw) cpu timing:", finish-start
      if(.not.ScreenAll_)write(*,"(A)") "     Screening only Kanamori NaNb terms."
      !
      if(.not.LocalOnly_)then
         write(*,"(A)") "     Filling the K-dependent attributes."
         do ik=1,Umats%Nkpt
            Umats_ptr%bare(:,:,ik) = Umats_ptr%bare_local
            Umats_ptr%screened(:,:,:,ik) = Umats_ptr%screened_local
         enddo
      endif
      !
      if(Hetero%status)then
         call Expand2Nsite(Umats,Umats_ptr,Hetero%Nlayer)
         call DeallocateBosonicField(Umats_imp)
      endif
      nullify(Umats_ptr)
      !
      call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
      !
   end subroutine build_Uret_multiParam_ph


   !---------------------------------------------------------------------------!
   !PURPOSE: Create the K-dependent interaction tensor from user-given
   !         long-range couplings.
   !---------------------------------------------------------------------------!
   subroutine build_Uret_singlParam_Vn(Umats,Uaa,Uab,J,Vnn,Lttc,Hetero,LocalOnly)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use input_vars, only : pathINPUTtr, pathINPUT
      use input_vars, only : RealPrint, long_range, structure, Nkpt_path, attach_Coulomb
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      real(8),intent(in)                    :: Uaa,Uab,J
      real(8),intent(in)                    :: Vnn(:,:,:)
      type(Lattice),intent(inout)           :: Lttc
      type(Heterostructures),intent(in)     :: Hetero
      logical,intent(in),optional           :: LocalOnly
      !
      complex(8),allocatable                :: Ur_bulk(:,:,:),Ur(:,:,:)
      complex(8),allocatable                :: Uk(:,:,:),URnn(:,:)
      integer                               :: Norb,Nsite_bulk,Vrange
      integer                               :: ib1,ib2,ib1_l,ib2_l
      integer                               :: iorb,jorb,io,jo,io_l,jo_l
      integer                               :: iwig,iD,iR,unit
      integer                               :: iprint,Nprint
      complex(8),allocatable                :: EwaldShift(:)
      real(8)                               :: eta,den,V
      logical                               :: LocalOnly_,CoulombTail
      real                                  :: start,finish
      !multi-site
      real(8),allocatable                   :: Ruc(:,:)
      integer                               :: isite,jsite
      integer                               :: ilayer,jlayer
      integer                               :: site_i,site_j
      real(8)                               :: Rvec(3),Rdist,Rdist_last
      real(8),allocatable                   :: Rsorted(:,:),Rsorted_bkp(:,:)
      integer,allocatable                   :: Rorder(:),Dist(:,:),DistList(:)
      !
      !
      if(verbose)write(*,"(A)") "---- build_Uret_singlParam_Vn"
      !
      !
      ! Check on the input field
      if(.not.Umats%status) stop "build_Uret_singlParam_Vn: BosonicField not properly initialized."
      if(Umats%Npoints.ne.1) stop "build_Uret_singlParam_Vn: Number of matsubara points in Umats is supposed to be equal to 1."
      !
      LocalOnly_=.false.
      if(present(LocalOnly))LocalOnly_=LocalOnly
      if(LocalOnly_.and.(Umats%Nkpt.ne.0)) stop "build_Uret_singlParam_Vn: Umats k dependent attributes are supposed to be unallocated."
      if((.not.LocalOnly_).and.(Umats%Nkpt.ne.Lttc%Nkpt)) stop "build_Uret_singlParam_Vn: Umats number of K-points does not match with the lattice."
      !
      Norb = int(sqrt(dble(Umats%Nbp)))/Lttc%Nsite
      Nsite_bulk = Lttc%Nsite
      if(Hetero%status) Nsite_bulk = int(Lttc%Nsite/Hetero%Nlayer)
      !
      if(LocalOnly_)then
         !
         write(*,"(A)") "     building interaction neglecting non-local dependence."
         !
         allocate(Ur_bulk((Norb*Nsite_bulk)**2,(Norb*Nsite_bulk)**2,1));Ur_bulk=czero
         do jsite=1,Nsite_bulk
            do isite=1,Nsite_bulk
               !
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     io = iorb + Norb*(isite-1)
                     jo = jorb + Norb*(jsite-1)
                     !
                     if(isite.eq.jsite)then
                        !
                        !local Uaa and Uab
                        call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                        if(io.eq.jo) Ur_bulk(ib1,ib2,1) = dcmplx(Uaa,0d0)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,1) = dcmplx(Uab,0d0)
                        !
                        !local Jsf and Jph
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,1) = dcmplx(J,0d0)
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,1) = dcmplx(J,0d0)
                        !
                     endif
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !
         if(Hetero%status)then
            !
            allocate(Ur((Norb*Lttc%Nsite)**2,(Norb*Lttc%Nsite)**2,1));Ur=czero
            !
            !reshuffle the single-layer Ur into the heterostructured one
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  do jsite=1,Nsite_bulk
                     do iorb=1,Norb
                        do jorb=1,Norb
                           !
                           io = iorb + Norb*(isite-1)
                           jo = jorb + Norb*(jsite-1)
                           !
                           io_l = io + Norb*Nsite_bulk*(ilayer-1)
                           jo_l = jo + Norb*Nsite_bulk*(ilayer-1)
                           !
                           !Uaa and Uab
                           call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,io_l],[jo_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,1) = Ur_bulk(ib1,ib2,1)
                           !
                           !Jsf
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[jo_l,io_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,1) = Ur_bulk(ib1,ib2,1)
                           !
                           !Jph
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[io_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,1) = Ur_bulk(ib1,ib2,1)
                           !
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !
         endif
         deallocate(Ur_bulk)
         !
         !fill in the output
         Umats%screened_local(:,:,1) = Ur(:,:,1)
         if(allocated(Umats%bare)) Umats%bare_local = Ur(:,:,1)
         !
         deallocate(Ur)
         !
      else
         !
         if((reg(long_range).eq."Ewald").and.(Nsite_bulk.ne.1)) stop "build_Uret_singlParam_Vn: Ewald long-range interaction is not implemented for multi-site setups."
         !
         Vrange = size(Vnn,dim=3)
         call assert_shape(Vnn,[Norb,Norb,Vrange],"build_Uret_singlParam_Vn","Vnn")
         !
         !recover the vectors in real space
         if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
         call get_Ruc(Ruc)
         !
         !Get all the possible positions
         allocate(Rsorted(Nwig*Nsite_bulk*Nsite_bulk,4));Rsorted=0d0
         iR=0
         do iwig=1,Nwig
            do jsite=1,Nsite_bulk
               do isite=1,Nsite_bulk
                  !
                  Rvec = Rvecwig(:,iwig) + Ruc(:,jsite) - Ruc(:,isite)
                  Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                  !
                  iR = iR +1
                  !
                  Rsorted(iR,1) = Rdist
                  Rsorted(iR,2) = iwig
                  Rsorted(iR,3) = jsite
                  Rsorted(iR,4) = isite
                  !
               enddo
            enddo
         enddo
         !
         !Sorting the positions according to distance
         allocate(Rorder(Nwig*Nsite_bulk*Nsite_bulk));Rorder=0
         allocate(Rsorted_bkp(Nwig*Nsite_bulk*Nsite_bulk,4));Rsorted_bkp=0d0
         Rsorted_bkp = Rsorted
         call sort_array(Rsorted(:,1),Rorder)
         Rsorted=0d0
         do iR=1,Nwig*Nsite_bulk*Nsite_bulk
            Rsorted(iR,:) = Rsorted_bkp(Rorder(iR),:)
         enddo
         deallocate(Rsorted_bkp,Rorder)
         !
         !Regroup according to distance. The list contains the indexes of all the positions with a given distance
         call get_pattern(Dist,Rsorted(:,1),1e4*eps,listDim=DistList,IncludeSingle=.true.)
         !
         !Compute Ewald shift only if Nsite_bulk=1
         if(reg(long_range).eq."Ewald")then
            eta = Rsorted(Nwig,1)/2d0
            allocate(EwaldShift(Nwig));EwaldShift=czero
            if(any(Lttc%Nkpt3.eq.1))then
               call calc_Ewald(EwaldShift,Lttc%kpt,eta,"2D")
               den = 2d0*eta
            else
               call calc_Ewald(EwaldShift,Lttc%kpt,eta,"3D")
               den = 2d0*sqrt(eta)
            endif
         endif
         !
         !User-provided non-local interaction
         call cpu_time(start)
         allocate(Ur_bulk((Norb*Nsite_bulk)**2,(Norb*Nsite_bulk)**2,Nwig));Ur_bulk=czero
         !all the possible ranges
         do iD=1,size(Dist,dim=1)
            !
            !iD=1 is the local, iD=2 is the nearest neighbor and so on
            CoulombTail = (reg(long_range).eq."Explicit") .and. (attach_Coulomb.gt.Vrange) .and. ((iD-1).le.attach_Coulomb)
            if((reg(long_range).ne."Ewald").and.((iD-1).gt.Vrange).and.(.not.CoulombTail))exit
            !
            !all the indexes within that range
            do iR=1,DistList(iD)
               !
               if((iD.eq.1).and.(Rsorted(Dist(iD,iR),2).ne.wig0)) stop "build_Uret_singlParam_Vn: wrong index of R=0 vector."
               if((iD.eq.1).and.(Rsorted(Dist(iD,iR),1).ne.0d0)) stop "build_Uret_singlParam_Vn: wrong length of R=0 vector."
               !
               !retrieve indexes from sorted list
               iwig = Rsorted(Dist(iD,iR),2)
               jsite = Rsorted(Dist(iD,iR),3)
               isite = Rsorted(Dist(iD,iR),4)
               !
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     !site-orbital indexes of Ur
                     io = iorb + Norb*(isite-1)
                     jo = jorb + Norb*(jsite-1)
                     !
                     if((isite.eq.jsite) .and. (iwig.eq.wig0))then
                        !
                        !local Uaa and Uab
                        call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                        if(io.eq.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(Uaa,0d0)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(Uab,0d0)
                        !
                        !local Jsf and Jph
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(J,0d0)
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(J,0d0)
                        !
                     else
                        !
                        Rdist = Rsorted(Dist(iD,iR),1)
                        !
                        if(iD.eq.1)  stop "build_Uret_singlParam_Vn: wrong R=0 local index in position list."
                        if(Rdist.eq.0d0)  stop "build_Uret_singlParam_Vn: wrong R=0 distance in position list."
                        !
                        !type of long-range interaction
                        if(reg(long_range).eq."Explicit")then
                           if((iD-1).le.Vrange)then
                              V = Vnn(iorb,jorb,iD-1)
                              Rdist_last = Rdist
                           else
                              V = Vnn(iorb,jorb,Vrange)*Rdist_last/Rdist
                           endif
                        elseif(reg(long_range).eq."Coulomb")then
                           V = Vnn(iorb,jorb,1)/Rdist
                        elseif(reg(long_range).eq."Ewald")then
                           V = (Vnn(iorb,jorb,1)/Rdist)*erfc(Rdist/den) + EwaldShift(iwig)
                        elseif(reg(long_range).eq."None")then !redundant since there is also the LocalOnly flag
                           V = 0d0
                        endif
                        !
                        !long-range interction added only to desity-density components
                        call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                        Ur_bulk(ib1,ib2,iwig) =  dcmplx(V,0d0)
                        !
                     endif
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         if(allocated(EwaldShift))deallocate(EwaldShift)
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of U(R) cpu timing:", finish-start
         !
         !if(verbose)then
            unit = free_unit()
            open(unit,file=reg(pathINPUTtr)//"Ur_report.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iD=1,size(Dist,dim=1)
               write(unit,"(A)") "     Dist: "//str(iD)
               write(unit,"(A8,6A6,1A12)") "ndx" , "n1" , "n2" , "n3" , "iwig" , "jsite" , "isite" , "R"
               do iR=1,DistList(iD)
                  iwig = Rsorted(Dist(iD,iR),2)
                  jsite = Rsorted(Dist(iD,iR),3)
                  isite = Rsorted(Dist(iD,iR),4)
                  Rvec = Rvecwig(:,iwig) + Ruc(:,jsite) - Ruc(:,isite)
                  Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                  write(unit,"(I8,6I6,1F12.4)") iR,Nvecwig(:,iwig),iwig,jsite,isite,Rsorted(Dist(iD,iR),1)
                  if(abs(Rsorted(Dist(iD,iR),1)-Rdist).gt.eps)then
                     write(unit,"(A,2E20.12)") "ERROR: Rsorted(Dist(iD,iR),1).ne.Rdist",Rsorted(Dist(iD,iR),1),Rdist
                     stop "build_Uret_singlParam_Vn: check Ur_report.DAT"
                  endif
                  do iorb=1,Norb
                     do jorb=1,Norb
                        ib1 = iorb + Norb*(isite-1) + Norb*Nsite_bulk*(iorb + Norb*(isite-1)-1)
                        ib2 = jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jorb + Norb*(jsite-1)-1)
                        write(unit,"(A5,A16,1F12.4)") " ","U"//str(iorb)//","//str(jorb)//"(Ri):["//str(ib1)//","//str(ib2)//"] ",dreal(Ur_bulk(ib1,ib2,iwig))
                     enddo
                  enddo
               enddo
            enddo
         !endif
         deallocate(Rsorted,Dist,DistList)
         !
         !Set up the heterostructure
         allocate(Ur((Norb*Lttc%Nsite)**2,(Norb*Lttc%Nsite)**2,Nwig));Ur=czero
         if(Hetero%status)then
            !
            !reshuffle the single-layer Ur into the heterostructured one
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  do jsite=1,Nsite_bulk
                     do iorb=1,Norb
                        do jorb=1,Norb
                           !
                           io = iorb + Norb*(isite-1)
                           jo = jorb + Norb*(jsite-1)
                           !
                           io_l = io + Norb*Nsite_bulk*(ilayer-1)
                           jo_l = jo + Norb*Nsite_bulk*(ilayer-1)
                           !
                           !Uaa and Uab
                           call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,io_l],[jo_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,:) = Ur_bulk(ib1,ib2,:)
                           !
                           !Jsf
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[jo_l,io_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,:) = Ur_bulk(ib1,ib2,:)
                           !
                           !Jph
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[io_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,:) = Ur_bulk(ib1,ib2,:)
                           !
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !
            !adding the interaction between the layers. This implies that the Ruc are ordered correctly!
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  !
                  allocate(Rsorted(Lttc%Nsite*Nwig,4));Rsorted=0d0
                  allocate(Rsorted_bkp(Lttc%Nsite*Nwig,4));Rsorted_bkp=0d0
                  allocate(Rorder(Lttc%Nsite*Nwig));Rorder=0
                  !
                  site_i = isite + (ilayer-1)*Nsite_bulk
                  !
                  iR=0
                  do jlayer=1,Hetero%Nlayer
                     do jsite=1,Nsite_bulk
                        !
                        site_j = jsite + (jlayer-1)*Nsite_bulk
                        !
                        do iwig=1,Nwig
                           !
                           Rvec = Rvecwig(:,iwig) + Ruc(:,site_j) - Ruc(:,site_i)
                           Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                           !
                           iR = iR +1
                           !
                           Rsorted(iR,1) = Rdist
                           Rsorted(iR,2) = iwig
                           Rsorted(iR,3) = jsite
                           Rsorted(iR,4) = jlayer
                           !
                        enddo
                     enddo
                  enddo
                  !
                  !Re-ordering distances for each site in the system
                  Rorder=0
                  Rsorted_bkp = Rsorted
                  call sort_array(Rsorted(:,1),Rorder)
                  do iR=1,Nwig*Lttc%Nsite
                     Rsorted(iR,:) = Rsorted_bkp(Rorder(iR),:)
                  enddo
                  !
                  !Regrouping according to distance. The list contains the indexes of all the positions with a given distance
                  call get_pattern(Dist,Rsorted(:,1),1e4*eps,listDim=DistList,IncludeSingle=.true.)
                  !
                  !add all the inter-layer interactions
                  do iD=1,size(Dist,dim=1)
                     !
                     !iD=1 is the local, iD=2 is the nearest neighbor and so on
                     CoulombTail = (reg(long_range).eq."Explicit") .and. (attach_Coulomb.gt.Vrange) .and. ((iD-1).le.attach_Coulomb)
                     if((reg(long_range).ne."Ewald").and.((iD-1).gt.Vrange).and.(.not.CoulombTail))exit
                     !
                     !all the indexes within that range
                     do iR=1,DistList(iD)
                        !
                        !retrieve indexes from sorted list
                        Rdist = Rsorted(Dist(iD,iR),1)
                        iwig = Rsorted(Dist(iD,iR),2)
                        jsite = Rsorted(Dist(iD,iR),3)
                        jlayer = Rsorted(Dist(iD,iR),4)
                        !
                        if(ilayer.ne.jlayer)then
                           !
                           do iorb=1,Norb
                              do jorb=1,Norb
                                 !
                                 !type of long-range interaction
                                 if(reg(long_range).eq."Explicit")then
                                    if((iD-1).le.Vrange)then
                                       V = Vnn(iorb,jorb,iD-1)
                                       Rdist_last = Rdist
                                    else
                                       V = Vnn(iorb,jorb,Vrange)*Rdist_last/Rdist
                                    endif
                                 elseif(reg(long_range).eq."Coulomb")then
                                    V = Vnn(iorb,jorb,1)/Rdist
                                 elseif(reg(long_range).eq."None")then
                                    V = 0d0
                                 endif
                                 !
                                 !site-orbital indexes of Hr
                                 io = iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1)
                                 jo = jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1)
                                 !
                                 !long-range interction added only to desity-density components
                                 call F2Bindex(Norb*Lttc%Nsite,[io,io],[jo,jo],ib1,ib2)
                                 !
                                 Ur(ib1,ib2,iwig) = dcmplx(V,0d0)
                                 !
                              enddo
                           enddo
                           !
                        endif
                        !
                     enddo
                  enddo
                  !
                  !if(verbose)then
                     unit = free_unit()
                     open(unit,file=reg(pathINPUTtr)//"Ur_report_Hetero_layer"//str(site_i)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
                     do iD=1,size(Dist,dim=1)
                        write(unit,"(A)") "     Dist: "//str(iD)
                        write(unit,"(A8,8A6,20A12)") "ndx" , "n1" , "n2" , "n3" , "iwig" , "ilay" , "isite" , "jlay" , "jsite" , "R", "H(Ri)"
                        do iR=1,DistList(iD)
                           iwig = Rsorted(Dist(iD,iR),2)
                           jsite = Rsorted(Dist(iD,iR),3)
                           jlayer = Rsorted(Dist(iD,iR),4)
                           Rvec = Rvecwig(:,iwig) + Ruc(:,jsite+Nsite_bulk*(jlayer-1)) - Ruc(:,isite+Nsite_bulk*(ilayer-1))
                           Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                           write(unit,"(I8,8I6,1F12.4)") iR,Nvecwig(:,iwig),iwig,ilayer,isite,jlayer,jsite,Rsorted(Dist(iD,iR),1)
                           if(abs(Rsorted(Dist(iD,iR),1)-Rdist).gt.eps)then
                              write(unit,"(A,2E20.12)") "ERROR: Rsorted(Dist(iD,iR),1).ne.Rdist",Rsorted(Dist(iD,iR),1),Rdist
                              stop "build_Uret_multiParam_Vn:  check Ur_report_Hetero.DAT"
                           endif
                           do iorb=1,Norb
                              do jorb=1,Norb
                                 ib1 = iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1) + (iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1)-1)*Norb*Lttc%Nsite
                                 ib2 = jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1) + (jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1)-1)*Norb*Lttc%Nsite
                                 write(unit,"(A5,A16,1F12.4)") " ","U"//str(iorb)//","//str(jorb)//"(Ri):["//str(ib1)//","//str(ib2)//"] ",dreal(Ur(ib1,ib2,iwig))
                              enddo
                           enddo
                        enddo
                     enddo
                  !endif
                  !
                  deallocate(Rsorted,Rsorted_bkp,Rorder,Dist,DistList)
                  !
               enddo !isite
            enddo !ilayer
            !
         else
            !
            Ur = Ur_bulk
            deallocate(Ur_bulk)
            !
         endif
         !
         !FT to K-space
         call cpu_time(start)
         allocate(Uk((Norb*Lttc%Nsite)**2,(Norb*Lttc%Nsite)**2,Lttc%Nkpt));Uk=czero
         if(Lttc%Nkpt.gt.1)then
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Ur,Uk)
         else
            Uk(:,:,1) = Ur(:,:,wig0)
         endif
         deallocate(Ur)
         where(abs((Uk))<eps) Uk=czero
         call cpu_time(finish)
         write(*,"(A,F)") "     U(R) --> U(K) cpu timing:", finish-start
         !
         !print along path
         if(reg(structure).ne."None")then
            call interpolate2Path(Lttc,Nkpt_path,"Uk",pathOUTPUT=reg(pathINPUT),store=.false.,skipAkw=.true.,data_in=Uk)
         endif
         !
         !fill in the output
         Umats%screened(:,:,1,:) = Uk
         if(allocated(Umats%bare)) Umats%bare = Uk
         call BosonicKsum(Umats)
         !
      endif
      !
      call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
      !
      ! Print the nn non local interaction
      Nprint = size(RealPrint,dim=2)
      allocate(URnn(int(sqrt(dble(Umats%Nbp))),int(sqrt(dble(Umats%Nbp)))));URnn=czero
      allocate(UR(Umats%Nbp,Umats%Nbp,Nprint));UR=czero
      call wannier_K2R_NN(RealPrint,Lttc%Nkpt3,Lttc%kpt,Umats%screened(:,:,1,:),UR)
      where(abs((UR))<eps) UR=czero
      do iprint=1,Nprint
         call dump_Matrix(UR(:,:,iprint),reg(trim(pathINPUTtr)),"U_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
         call product2NN(UR(:,:,iprint),URnn)
         call dump_Matrix(URnn,reg(trim(pathINPUTtr)),"Unn_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
      enddo
      deallocate(UR,URnn)
      !
   end subroutine build_Uret_singlParam_Vn
   !
   subroutine build_Uret_multiParam_Vn(Umats,Uaa,Uab,J,Vnn,Lttc,Hetero,LocalOnly)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use input_vars, only : pathINPUTtr, pathINPUT
      use input_vars, only : RealPrint, long_range, structure, Nkpt_path, attach_Coulomb
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      real(8),intent(in)                    :: Uaa(:),Uab(:,:),J(:,:)
      real(8),intent(in)                    :: Vnn(:,:,:)
      type(Lattice),intent(inout)           :: Lttc
      type(Heterostructures),intent(in)     :: Hetero
      logical,intent(in),optional           :: LocalOnly
      !
      complex(8),allocatable                :: Ur_bulk(:,:,:),Ur(:,:,:)
      complex(8),allocatable                :: Uk(:,:,:),URnn(:,:)
      integer                               :: Norb,Nsite_bulk,Vrange
      integer                               :: ib1,ib2,ib1_l,ib2_l
      integer                               :: iorb,jorb,io,jo,io_l,jo_l
      integer                               :: iwig,iD,iR,unit
      integer                               :: iprint,Nprint
      complex(8),allocatable                :: EwaldShift(:)
      real(8)                               :: eta,den,V
      logical                               :: LocalOnly_,CoulombTail
      real                                  :: start,finish
      !multi-site
      real(8),allocatable                   :: Ruc(:,:)
      integer                               :: isite,jsite
      integer                               :: ilayer,jlayer
      integer                               :: site_i,site_j
      real(8)                               :: Rvec(3),Rdist,Rdist_last
      real(8),allocatable                   :: Rsorted(:,:),Rsorted_bkp(:,:)
      integer,allocatable                   :: Rorder(:),Dist(:,:),DistList(:)
      !
      !
      if(verbose)write(*,"(A)") "---- build_Uret_multiParam_Vn"
      !
      !
      ! Check on the input field
      if(.not.Umats%status) stop "build_Uret_multiParam_Vn: BosonicField not properly initialized."
      if(Umats%Npoints.ne.1) stop "build_Uret_multiParam_Vn: Number of matsubara points in Umats is supposed to be equal to 1."
      !
      LocalOnly_=.false.
      if(present(LocalOnly))LocalOnly_=LocalOnly
      if(LocalOnly_.and.(Umats%Nkpt.ne.0)) stop "build_Uret_multiParam_Vn: Umats k dependent attributes are supposed to be unallocated."
      if((.not.LocalOnly_).and.(Umats%Nkpt.ne.Lttc%Nkpt)) stop "build_Uret_multiParam_Vn: Umats number of K-points does not match with the lattice."
      !
      Norb = int(sqrt(dble(Umats%Nbp)))/Lttc%Nsite
      Nsite_bulk = Lttc%Nsite
      if(Hetero%status) Nsite_bulk = int(Lttc%Nsite/Hetero%Nlayer)
      !
      call assert_shape(Uaa,[Norb],"build_Uret_multiParam_Vn","Uaa")
      call assert_shape(Uab,[Norb,Norb],"build_Uret_multiParam_Vn","Uab")
      call assert_shape(J,[Norb,Norb],"build_Uret_multiParam_Vn","J")
      !
      if(LocalOnly_)then
         !
         write(*,"(A)") "     building interaction neglecting non-local dependence."
         !
         allocate(Ur_bulk((Norb*Nsite_bulk)**2,(Norb*Nsite_bulk)**2,1));Ur_bulk=czero
         do jsite=1,Nsite_bulk
            do isite=1,Nsite_bulk
               !
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     io = iorb + Norb*(isite-1)
                     jo = jorb + Norb*(jsite-1)
                     !
                     if(isite.eq.jsite)then
                        !
                        !local Uaa and Uab
                        call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                        if(io.eq.jo) Ur_bulk(ib1,ib2,1) = dcmplx(Uaa(iorb),0d0)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,1) = dcmplx(Uab(iorb,jorb),0d0)
                        !
                        !local Jsf and Jph
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,1) = dcmplx(J(iorb,jorb),0d0)
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,1) = dcmplx(J(iorb,jorb),0d0)
                        !
                     endif
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !
         if(Hetero%status)then
            !
            allocate(Ur((Norb*Lttc%Nsite)**2,(Norb*Lttc%Nsite)**2,1));Ur=czero
            !
            !reshuffle the single-layer Ur into the heterostructured one
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  do jsite=1,Nsite_bulk
                     do iorb=1,Norb
                        do jorb=1,Norb
                           !
                           io = iorb + Norb*(isite-1)
                           jo = jorb + Norb*(jsite-1)
                           !
                           io_l = io + Norb*Nsite_bulk*(ilayer-1)
                           jo_l = jo + Norb*Nsite_bulk*(ilayer-1)
                           !
                           !Uaa and Uab
                           call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,io_l],[jo_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,1) = Ur_bulk(ib1,ib2,1)
                           !
                           !Jsf
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[jo_l,io_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,1) = Ur_bulk(ib1,ib2,1)
                           !
                           !Jph
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[io_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,1) = Ur_bulk(ib1,ib2,1)
                           !
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !
         endif
         deallocate(Ur_bulk)
         !
         !fill in the output
         Umats%screened_local(:,:,1) = Ur(:,:,1)
         if(allocated(Umats%bare)) Umats%bare_local = Ur(:,:,1)
         !
         deallocate(Ur)
         !
      else
         !
         if((reg(long_range).eq."Ewald").and.(Nsite_bulk.ne.1)) stop "build_Uret_singlParam_Vn: Ewald long-range interaction is not implemented for multi-site setups."
         !
         Vrange = size(Vnn,dim=3)
         call assert_shape(Vnn,[Norb,Norb,Vrange],"build_Uret_singlParam_Vn","Vnn")
         !
         !recover the vectors in real space
         if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
         call get_Ruc(Ruc)
         !
         !Get all the possible positions
         allocate(Rsorted(Nwig*Nsite_bulk*Nsite_bulk,4));Rsorted=0d0
         iR=0
         do iwig=1,Nwig
            do jsite=1,Nsite_bulk
               do isite=1,Nsite_bulk
                  !
                  Rvec = Rvecwig(:,iwig) + Ruc(:,jsite) - Ruc(:,isite)
                  Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                  !
                  iR = iR +1
                  !
                  Rsorted(iR,1) = Rdist
                  Rsorted(iR,2) = iwig
                  Rsorted(iR,3) = jsite
                  Rsorted(iR,4) = isite
                  !
               enddo
            enddo
         enddo
         !
         !Sorting the positions according to distance
         allocate(Rorder(Nwig*Nsite_bulk*Nsite_bulk));Rorder=0
         allocate(Rsorted_bkp(Nwig*Nsite_bulk*Nsite_bulk,4));Rsorted_bkp=0d0
         Rsorted_bkp = Rsorted
         call sort_array(Rsorted(:,1),Rorder)
         Rsorted=0d0
         do iR=1,Nwig*Nsite_bulk*Nsite_bulk
            Rsorted(iR,:) = Rsorted_bkp(Rorder(iR),:)
         enddo
         deallocate(Rsorted_bkp,Rorder)
         !
         !Regroup according to distance. The list contains the indexes of all the positions with a given distance
         call get_pattern(Dist,Rsorted(:,1),1e4*eps,listDim=DistList,IncludeSingle=.true.)
         !
         !Compute Ewald shift only if Nsite_bulk=1
         if(reg(long_range).eq."Ewald")then
            eta = Rsorted(Nwig,1)/2d0
            allocate(EwaldShift(Nwig));EwaldShift=czero
            if(any(Lttc%Nkpt3.eq.1))then
               call calc_Ewald(EwaldShift,Lttc%kpt,eta,"2D")
               den = 2d0*eta
            else
               call calc_Ewald(EwaldShift,Lttc%kpt,eta,"3D")
               den = 2d0*sqrt(eta)
            endif
         endif
         !
         !User-provided non-local interaction
         call cpu_time(start)
         allocate(Ur_bulk((Norb*Nsite_bulk)**2,(Norb*Nsite_bulk)**2,Nwig));Ur_bulk=czero
         !all the possible ranges
         do iD=1,size(Dist,dim=1)
            !
            !iD=1 is the local, iD=2 is the nearest neighbor and so on
            CoulombTail = (reg(long_range).eq."Explicit") .and. (attach_Coulomb.gt.Vrange) .and. ((iD-1).le.attach_Coulomb)
            if((reg(long_range).ne."Ewald").and.((iD-1).gt.Vrange).and.(.not.CoulombTail))exit
            !
            !all the indexes within that range
            do iR=1,DistList(iD)
               !
               if((iD.eq.1).and.(Rsorted(Dist(iD,iR),2).ne.wig0)) stop "build_Uret_singlParam_Vn: wrong index of R=0 vector."
               if((iD.eq.1).and.(Rsorted(Dist(iD,iR),1).ne.0d0)) stop "build_Uret_singlParam_Vn: wrong length of R=0 vector."
               !
               !retrieve indexes from sorted list
               iwig = Rsorted(Dist(iD,iR),2)
               jsite = Rsorted(Dist(iD,iR),3)
               isite = Rsorted(Dist(iD,iR),4)
               !
               do iorb=1,Norb
                  do jorb=1,Norb
                     !
                     !site-orbital indexes of Ur
                     io = iorb + Norb*(isite-1)
                     jo = jorb + Norb*(jsite-1)
                     !
                     if((isite.eq.jsite) .and. (iwig.eq.wig0))then
                        !
                        !local Uaa and Uab
                        call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                        if(io.eq.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(Uaa(iorb),0d0)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(Uab(iorb,jorb),0d0)
                        !
                        !local Jsf and Jph
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(J(iorb,jorb),0d0)
                        call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                        if(io.ne.jo) Ur_bulk(ib1,ib2,iwig) = dcmplx(J(iorb,jorb),0d0)
                        !
                     else
                        !
                        Rdist = Rsorted(Dist(iD,iR),1)
                        !
                        if(iD.eq.1)  stop "build_Uret_singlParam_Vn: wrong R=0 local index in position list."
                        if(Rdist.eq.0d0)  stop "build_Uret_singlParam_Vn: wrong R=0 distance in position list."
                        !
                        !type of long-range interaction
                        if(reg(long_range).eq."Explicit")then
                           if((iD-1).le.Vrange)then
                              V = Vnn(iorb,jorb,iD-1)
                              Rdist_last = Rdist
                           else
                              V = Vnn(iorb,jorb,Vrange)*Rdist_last/Rdist
                           endif
                        elseif(reg(long_range).eq."Coulomb")then
                           V = Vnn(iorb,jorb,1)/Rdist
                        elseif(reg(long_range).eq."Ewald")then
                           V = (Vnn(iorb,jorb,1)/Rdist)*erfc(Rdist/den) + EwaldShift(iwig)
                        elseif(reg(long_range).eq."None")then !redundant since there is also the LocalOnly flag
                           V = 0d0
                        endif
                        !
                        !long-range interction added only to desity-density components
                        call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                        Ur_bulk(ib1,ib2,iwig) =  dcmplx(V,0d0)
                        !
                     endif
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         if(allocated(EwaldShift))deallocate(EwaldShift)
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of U(R) cpu timing:", finish-start
         !
         !if(verbose)then
            unit = free_unit()
            open(unit,file=reg(pathINPUTtr)//"Ur_report.DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iD=1,size(Dist,dim=1)
               write(unit,"(A)") "     Dist: "//str(iD)
               write(unit,"(A8,6A6,1A12)") "ndx" , "n1" , "n2" , "n3" , "iwig" , "jsite" , "isite" , "R"
               do iR=1,DistList(iD)
                  iwig = Rsorted(Dist(iD,iR),2)
                  jsite = Rsorted(Dist(iD,iR),3)
                  isite = Rsorted(Dist(iD,iR),4)
                  Rvec = Rvecwig(:,iwig) + Ruc(:,jsite) - Ruc(:,isite)
                  Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                  write(unit,"(I8,6I6,1F12.4)") iR,Nvecwig(:,iwig),iwig,jsite,isite,Rsorted(Dist(iD,iR),1)
                  if(abs(Rsorted(Dist(iD,iR),1)-Rdist).gt.eps)then
                     write(unit,"(A,2E20.12)") "ERROR: Rsorted(Dist(iD,iR),1).ne.Rdist",Rsorted(Dist(iD,iR),1),Rdist
                     stop "build_Uret_singlParam_Vn:  check Ur_report.DAT"
                  endif
                  do iorb=1,Norb
                     do jorb=1,Norb
                        ib1 = iorb + Norb*(isite-1) + Norb*Nsite_bulk*(iorb + Norb*(isite-1)-1)
                        ib2 = jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jorb + Norb*(jsite-1)-1)
                        write(unit,"(A5,A16,1F12.4)") " ","U"//str(iorb)//","//str(jorb)//"(Ri):["//str(ib1)//","//str(ib2)//"] ",dreal(Ur_bulk(ib1,ib2,iwig))
                     enddo
                  enddo
               enddo
            enddo
         !endif
         deallocate(Rsorted,Dist,DistList)
         !
         !Set up the heterostructure
         allocate(Ur((Norb*Lttc%Nsite)**2,(Norb*Lttc%Nsite)**2,Nwig));Ur=czero
         if(Hetero%status)then
            !
            !reshuffle the single-layer Ur into the heterostructured one
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  do jsite=1,Nsite_bulk
                     do iorb=1,Norb
                        do jorb=1,Norb
                           !
                           io = iorb + Norb*(isite-1)
                           jo = jorb + Norb*(jsite-1)
                           !
                           io_l = io + Norb*Nsite_bulk*(ilayer-1)
                           jo_l = jo + Norb*Nsite_bulk*(ilayer-1)
                           !
                           !Uaa and Uab
                           call F2Bindex(Norb*Nsite_bulk,[io,io],[jo,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,io_l],[jo_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,:) = Ur_bulk(ib1,ib2,:)
                           !
                           !Jsf
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[jo,io],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[jo_l,io_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,:) = Ur_bulk(ib1,ib2,:)
                           !
                           !Jph
                           call F2Bindex(Norb*Nsite_bulk,[io,jo],[io,jo],ib1,ib2)
                           call F2Bindex(Norb*Lttc%Nsite,[io_l,jo_l],[io_l,jo_l],ib1_l,ib2_l)
                           Ur(ib1_l,ib2_l,:) = Ur_bulk(ib1,ib2,:)
                           !
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !
            !adding the interaction between the layers. This implies that the Ruc are ordered correctly!
            do ilayer=1,Hetero%Nlayer
               do isite=1,Nsite_bulk
                  !
                  allocate(Rsorted(Lttc%Nsite*Nwig,4));Rsorted=0d0
                  allocate(Rsorted_bkp(Lttc%Nsite*Nwig,4));Rsorted_bkp=0d0
                  allocate(Rorder(Lttc%Nsite*Nwig));Rorder=0
                  !
                  site_i = isite + (ilayer-1)*Nsite_bulk
                  !
                  iR=0
                  do jlayer=1,Hetero%Nlayer
                     do jsite=1,Nsite_bulk
                        !
                        site_j = jsite + (jlayer-1)*Nsite_bulk
                        !
                        do iwig=1,Nwig
                           !
                           Rvec = Rvecwig(:,iwig) + Ruc(:,site_j) - Ruc(:,site_i)
                           Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                           !
                           iR = iR +1
                           !
                           Rsorted(iR,1) = Rdist
                           Rsorted(iR,2) = iwig
                           Rsorted(iR,3) = jsite
                           Rsorted(iR,4) = jlayer
                           !
                        enddo
                     enddo
                  enddo
                  !
                  !Re-ordering distances for each site in the system
                  Rorder=0
                  Rsorted_bkp = Rsorted
                  call sort_array(Rsorted(:,1),Rorder)
                  do iR=1,Nwig*Lttc%Nsite
                     Rsorted(iR,:) = Rsorted_bkp(Rorder(iR),:)
                  enddo
                  !
                  !Regrouping according to distance. The list contains the indexes of all the positions with a given distance
                  call get_pattern(Dist,Rsorted(:,1),1e4*eps,listDim=DistList,IncludeSingle=.true.)
                  !
                  !add the inter-layer interaction
                  do iD=1,size(Dist,dim=1)
                     !
                     !iD=1 is the local, iD=2 is the nearest neighbor and so on
                     CoulombTail = (reg(long_range).eq."Explicit") .and. (attach_Coulomb.gt.Vrange) .and. ((iD-1).le.attach_Coulomb)
                     if((reg(long_range).ne."Ewald").and.((iD-1).gt.Vrange).and.(.not.CoulombTail))exit
                     !
                     !all the indexes within that range
                     do iR=1,DistList(iD)
                        !
                        !retrieve indexes from sorted list
                        Rdist = Rsorted(Dist(iD,iR),1)
                        iwig = Rsorted(Dist(iD,iR),2)
                        jsite = Rsorted(Dist(iD,iR),3)
                        jlayer = Rsorted(Dist(iD,iR),4)
                        !
                        if(ilayer.ne.jlayer)then
                           !
                           do iorb=1,Norb
                              do jorb=1,Norb
                                 !
                                 !type of long-range interaction
                                 if(reg(long_range).eq."Explicit")then
                                    if((iD-1).le.Vrange)then
                                       V = Vnn(iorb,jorb,iD-1)
                                       Rdist_last = Rdist
                                    else
                                       V = Vnn(iorb,jorb,Vrange)*Rdist_last/Rdist
                                    endif
                                 elseif(reg(long_range).eq."Coulomb")then
                                    V = Vnn(iorb,jorb,1)/Rdist
                                 elseif(reg(long_range).eq."None")then
                                    V = 0d0
                                 endif
                                 !
                                 !site-orbital indexes of Hr
                                 io = iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1)
                                 jo = jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1)
                                 !
                                 !long-range interction added only to desity-density components
                                 call F2Bindex(Norb*Lttc%Nsite,[io,io],[jo,jo],ib1,ib2)
                                 !
                                 Ur(ib1,ib2,iwig) = dcmplx(V,0d0)
                                 !
                              enddo
                           enddo
                           !
                        endif
                        !
                     enddo
                  enddo
                  !
                  !if(verbose)then
                     unit = free_unit()
                     open(unit,file=reg(pathINPUTtr)//"Ur_report_Hetero_layer"//str(site_i)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
                     do iD=1,size(Dist,dim=1)
                        write(unit,"(A)") "     Dist: "//str(iD)
                        write(unit,"(A8,8A6,20A12)") "ndx" , "n1" , "n2" , "n3" , "iwig" , "ilay" , "isite" , "jlay" , "jsite" , "R", "H(Ri)"
                        do iR=1,DistList(iD)
                           iwig = Rsorted(Dist(iD,iR),2)
                           jsite = Rsorted(Dist(iD,iR),3)
                           jlayer = Rsorted(Dist(iD,iR),4)
                           Rvec = Rvecwig(:,iwig) + Ruc(:,jsite+Nsite_bulk*(jlayer-1)) - Ruc(:,isite+Nsite_bulk*(ilayer-1))
                           Rdist = sqrt(dble(dot_product(Rvec,Rvec)))
                           write(unit,"(I8,8I6,1F12.4)") iR,Nvecwig(:,iwig),iwig,ilayer,isite,jlayer,jsite,Rsorted(Dist(iD,iR),1)
                           if(abs(Rsorted(Dist(iD,iR),1)-Rdist).gt.eps)then
                              write(unit,"(A,2E20.12)") "ERROR: Rsorted(Dist(iD,iR),1).ne.Rdist",Rsorted(Dist(iD,iR),1),Rdist
                              stop "build_Uret_multiParam_Vn:  check Ur_report_Hetero.DAT"
                           endif
                           do iorb=1,Norb
                              do jorb=1,Norb
                                 ib1 = iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1) + (iorb + Norb*(isite-1) + Norb*Nsite_bulk*(ilayer-1)-1)*Norb*Lttc%Nsite
                                 ib2 = jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1) + (jorb + Norb*(jsite-1) + Norb*Nsite_bulk*(jlayer-1)-1)*Norb*Lttc%Nsite
                                 write(unit,"(A5,A16,1F12.4)") " ","U"//str(iorb)//","//str(jorb)//"(Ri):["//str(ib1)//","//str(ib2)//"] ",dreal(Ur(ib1,ib2,iwig))
                              enddo
                           enddo
                        enddo
                     enddo
                  !endif
                  !
                  deallocate(Rsorted,Rsorted_bkp,Rorder,Dist,DistList)
                  !
               enddo !isite
            enddo !ilayer
            !
         else
            !
            Ur = Ur_bulk
            deallocate(Ur_bulk)
            !
         endif
         !
         !FT to K-space
         call cpu_time(start)
         allocate(Uk((Norb*Lttc%Nsite)**2,(Norb*Lttc%Nsite)**2,Lttc%Nkpt));Uk=czero
         if(Lttc%Nkpt.gt.1)then
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Ur,Uk)
         else
            Uk(:,:,1) = Ur(:,:,wig0)
         endif
         deallocate(Ur)
         where(abs((Uk))<eps) Uk=czero
         call cpu_time(finish)
         write(*,"(A,F)") "     U(R) --> U(K) cpu timing:", finish-start
         !
         !print along path
         if(reg(structure).ne."None")then
            call interpolate2Path(Lttc,Nkpt_path,"Uk",pathOUTPUT=reg(pathINPUT),store=.false.,skipAkw=.true.,data_in=Uk)
         endif
         !
         !fill in the output
         Umats%screened(:,:,1,:) = Uk
         if(allocated(Umats%bare)) Umats%bare = Uk
         call BosonicKsum(Umats)
         !
      endif
      !
      call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
      !
      ! Print the nn non local interaction
      Nprint = size(RealPrint,dim=2)
      allocate(URnn(int(sqrt(dble(Umats%Nbp))),int(sqrt(dble(Umats%Nbp)))));URnn=czero
      allocate(UR(Umats%Nbp,Umats%Nbp,Nprint));UR=czero
      call wannier_K2R_NN(RealPrint,Lttc%Nkpt3,Lttc%kpt,Umats%screened(:,:,1,:),UR)
      where(abs((UR))<eps) UR=czero
      do iprint=1,Nprint
         call dump_Matrix(UR(:,:,iprint),reg(trim(pathINPUTtr)),"U_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
         call product2NN(UR(:,:,iprint),URnn)
         call dump_Matrix(URnn,reg(trim(pathINPUTtr)),"Unn_"//str(RealPrint(1,iprint))//str(RealPrint(2,iprint))//str(RealPrint(3,iprint))//".DAT")
      enddo
      deallocate(UR,URnn)
      !
   end subroutine build_Uret_multiParam_Vn


   !---------------------------------------------------------------------------!
   !PURPOSE: Given the Bosonic Field it extracts the screened interaction and
   ! retardation function.
   !---------------------------------------------------------------------------!
   subroutine calc_QMCinteractions_static(Uimp,Uinst)
      !
      use parameters
      use file_io
      use utils_misc
      use fourier_transforms
      implicit none
      !
      real(8),intent(in)                    :: Uimp(:,:)
      real(8),intent(inout)                 :: Uinst(:,:)
      !
      integer                               :: Nbp,Norb,Nflavor
      integer                               :: ib1,ib2,iorb,jorb
      integer                               :: iu1,iu2,ix1,ix2,ip1,ip2
      logical                               :: Uloc,U1st,U2nd
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- calc_QMCinteractions_static"
      !
      !
      Nbp = size(Uimp,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      Nflavor = Norb*Nspin
      !
      call assert_shape(Uinst,[Nflavor,Nflavor],"calc_QMCinteractions_static","Uinst")
      call assert_shape(Uimp,[Nbp,Nbp],"calc_QMCinteractions_static","Uimp")
      call init_Uelements(Norb,PhysicalUelements)
      !
      !computing the screened interaction
      Uinst=0d0
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            !This is just for a more compact coding
            Uloc = PhysicalUelements%Flav_Uloc(ib1,ib2)
            U1st = PhysicalUelements%Flav_U1st(ib1,ib2)
            U2nd = PhysicalUelements%Flav_U2nd(ib1,ib2)
            !
            !Orbital indexes
            iorb = PhysicalUelements%Flav_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Flav_Map(ib1,ib2,2)
            !
            ! (iorb,iorb)(jorb,jorb) indexes in the Norb^2 representaion
            call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],iu1,iu2)
            !
            ! (iorb,jorb)(jorb,iorb) indexes
            call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ix1,ix2)
            !
            ! (iorb,jorb)(iorb,jorb) indexes
            call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ip1,ip2)
            !
            if(Uloc) Uinst(ib1,ib2) = Uimp(iu1,iu2)
            if(U1st) Uinst(ib1,ib2) = Uimp(iu1,iu2)
            if(U2nd) Uinst(ib1,ib2) = Uimp(iu1,iu2) - (Uimp(ix1,ix2)+Uimp(ip1,ip2))/2d0
            !
         enddo
      enddo
      call check_Symmetry(Uinst,eps,enforce=.true.,hardstop=.false.,name="Uinst")
      !
   end subroutine calc_QMCinteractions_static
   !
   subroutine calc_QMCinteractions_retarded(Umats,Uinst,Kfunct,Ktilda,Screening,Kpfunct,sym)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use fourier_transforms
      use input_vars, only : Solver
      implicit none
      !
      type(BosonicField),intent(in)         :: Umats
      real(8),intent(inout)                 :: Uinst(:,:)
      real(8),intent(inout),optional        :: Kfunct(:,:,:)
      logical,intent(in),optional           :: Ktilda
      real(8),intent(inout),optional        :: Screening(:,:)
      real(8),intent(inout),optional        :: Kpfunct(:,:,:)
      logical,intent(in),optional           :: sym
      !
      integer                               :: Nbp,Norb,Nflavor
      integer                               :: ib1,ib2,iorb,jorb
      integer                               :: iu1,iu2,ix1,ix2,ip1,ip2
      integer                               :: iw,itau,iwlimit
      real(8),allocatable                   :: wmats(:),tau(:)
      real(8),allocatable                   :: Screening_(:,:)
      complex(8),allocatable                :: Kaux(:,:,:),Ktmp(:,:,:)
      logical                               :: Uloc,U1st,U2nd
      logical                               :: retarded,Kp_out,Scr_out
      type(physicalU)                       :: PhysicalUelements
      logical                               :: sym_,Ktilda_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_QMCinteractions_retarded"
      !
      !
      if(.not.Umats%status) stop "calc_QMCinteractions_retarded: Umats not properly initialized."
      !
      retarded=.false.
      if(present(Kfunct))retarded=.true.
      !
      Ktilda_=.true. !This is very delicate do not change.
      if(present(Ktilda).and.retarded)Ktilda_=Ktilda
      iwlimit = Umats%Npoints
      if(Ktilda_) iwlimit = 1
      !
      Kp_out=.false.
      if(present(Kpfunct).and.retarded) Kp_out=.true.
      !
      Scr_out=.false.
      if(present(Screening).and.retarded) Scr_out=.true.
      !
      sym_=.true.
      if(present(sym))sym_=sym
      !
      Nbp = Umats%Nbp
      Norb = int(sqrt(dble(Nbp)))
      Nflavor = Norb*Nspin
      !
      call assert_shape(Uinst,[Nflavor,Nflavor],"calc_QMCinteractions_retarded","Uinst")
      call init_Uelements(Norb,PhysicalUelements)
      !
      Uinst=0d0
      if(retarded)then
         call assert_shape(Kfunct,[Nflavor,Nflavor,Solver%NtauB_K],"calc_QMCinteractions_retarded","Kfunct")
         if(Kp_out)call assert_shape(Kpfunct,[Nflavor,Nflavor,Solver%NtauB_K],"calc_QMCinteractions_retarded","Kpfunct")
         if(Scr_out)call assert_shape(Screening,[Nflavor,Nflavor],"calc_QMCinteractions_retarded","Screening")
         allocate(Kaux(Nflavor,Nflavor,Umats%Npoints));Kaux=czero
         allocate(Ktmp(Nflavor,Nflavor,Solver%NtauB_K));Ktmp=czero
         allocate(Screening_(Nflavor,Nflavor));Screening_=0d0
         allocate(tau(Solver%NtauB_K));tau=0d0
         if(Solver%tau_uniform_K.eq.1)then
            tau = linspace(0d0,Umats%Beta,Solver%NtauB_K)
         else
            tau = denspace(Umats%Beta,Solver%NtauB_K)
         endif
         allocate(wmats(Umats%Npoints));wmats=0d0
         wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      endif
      !
      !computing the screened interaction
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            !This is just for a more compact coding
            Uloc = PhysicalUelements%Flav_Uloc(ib1,ib2)
            U1st = PhysicalUelements%Flav_U1st(ib1,ib2)
            U2nd = PhysicalUelements%Flav_U2nd(ib1,ib2)
            !
            !Orbital indexes
            iorb = PhysicalUelements%Flav_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Flav_Map(ib1,ib2,2)
            !
            !The maps inside PhysicalUelements contain separately the orbital
            !indexes specifially for that representation. The matching between
            !the two is not done, so I have to do it here.
            !
            ! (iorb,iorb)(jorb,jorb) indexes in the Norb^2 representaion
            call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],iu1,iu2)
            !
            ! (iorb,jorb)(jorb,iorb) indexes
            call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ix1,ix2)
            !
            ! (iorb,jorb)(iorb,jorb) indexes
            call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ip1,ip2)
            !
            if(Uloc) Uinst(ib1,ib2) = Umats%screened_local(iu1,iu2,1)
            if(U1st) Uinst(ib1,ib2) = Umats%screened_local(iu1,iu2,1)
            if(U2nd) Uinst(ib1,ib2) = Umats%screened_local(iu1,iu2,1) - (Umats%screened_local(ix1,ix2,1)+Umats%screened_local(ip1,ip2,1))/2d0
            !
            !auxiliary function to build the screening function
            if(retarded)then
               !
               if(Uloc) Kaux(ib1,ib2,:) =  Umats%screened_local(iu1,iu2,:) - Umats%screened_local(iu1,iu2,iwlimit)
               if(U1st) Kaux(ib1,ib2,:) =  Umats%screened_local(iu1,iu2,:) - Umats%screened_local(iu1,iu2,iwlimit)
               if(U2nd) Kaux(ib1,ib2,:) =  Umats%screened_local(iu1,iu2,:) - (Umats%screened_local(ix1,ix2,:)+Umats%screened_local(ip1,ip2,:))/2d0 - &
                                          (Umats%screened_local(iu1,iu2,iwlimit) - (Umats%screened_local(ix1,ix2,iwlimit)+Umats%screened_local(ip1,ip2,iwlimit))/2d0)
               !same orbital - same spin screening
               if(Uloc.and.(ib2.gt.ib1)) then
                  Kaux(ib1,ib1,:) = Kaux(ib1,ib2,:)
                  Kaux(ib2,ib2,:) = Kaux(ib1,ib2,:)
               endif
               !
               !removing numerical noise on the imaginary part
               Kaux(ib1,ib2,:) = dcmplx(dreal(Kaux(ib1,ib2,:)),0d0)
               !
            endif
            !
         enddo
      enddo
      !
      if(sym_) call check_Symmetry(Uinst,1e7*eps,enforce=.true.,hardstop=.false.,name="Uinst")
      !
      !computing the screening function and first derivative
      if(retarded)then
         !
         !save exact screening
         Screening_ = -Kaux(:,:,1)
         if(Ktilda_) Screening_ = Kaux(:,:,Umats%Npoints)
         !
         if(sym_) call check_Symmetry(Screening_,1e7*eps,enforce=.true.,hardstop=.false.,name="Screening_")
         !
         !This is the D(iw)-D(0)/iw^2 screening function
         do iw=2,Umats%Npoints
            Kaux(:,:,iw) = - (Kaux(:,:,iw)-Kaux(:,:,1)) / (wmats(iw)**2)
         enddo
         Kaux(:,:,1) = czero
         !
         Ktmp=czero
         call Bmats2itau(Umats%Beta,Kaux,Ktmp,asympt_corr=.true.,tau_uniform=(Solver%tau_uniform_K.eq.1))
         Kfunct=0d0
         do itau=2,Solver%NtauB_K-1
            Kfunct(:,:,itau) = dreal(Ktmp(:,:,itau) - Ktmp(:,:,1))
            if(sym_)call check_Symmetry(Kfunct(:,:,itau),1e7*eps,enforce=.true.,hardstop=.false.,name="K_t"//str(itau))
         enddo
         deallocate(Ktmp,Kaux)
         !
         !Enforce symmetry with respect to beta/2
         do ib1=1,Nflavor
            do ib2=1,Nflavor
               call halfbeta_sym(Kfunct(ib1,ib2,:),1d0)
            enddo
         enddo
         !
         !Mid-point derivative since the coefficients for Bmats2itau assume even bosonic functions
         allocate(Ktmp(Nflavor,Nflavor,Solver%NtauB_K));Ktmp=czero
         do itau=2,Solver%NtauB_K-1
            Ktmp(:,:,itau) = ( Kfunct(:,:,itau-1) - Kfunct(:,:,itau+1) ) / ( tau(itau-1)-tau(itau+1) )
         enddo
         Ktmp(:,:,1) = Screening_/2d0
         Ktmp(:,:,Solver%NtauB_K) = -Screening_/2d0
         !
         !Enforce anti-symmetry with respect to beta/2
         do ib1=1,Nflavor
            do ib2=1,Nflavor
               call halfbeta_sym(Ktmp(ib1,ib2,:),-1d0)
            enddo
         enddo
         !
         !pass data to output
         if(Kp_out)then
            Kpfunct = dreal(Ktmp)
            do itau=1,Solver%NtauB_K
               if(sym_)call check_Symmetry(Kpfunct(:,:,itau),1e7*eps,enforce=.true.,hardstop=.false.,name="Kp_t"//str(itau))
            enddo
            !cumbersome but exact
            Kpfunct(:,:,1) = Screening_/2d0
            Kpfunct(:,:,Solver%NtauB_K) = -Screening_/2d0
         endif
         if(Scr_out)Screening = Screening_
         !
         deallocate(Ktmp,tau,wmats,Screening_)
      endif
      !
   end subroutine calc_QMCinteractions_retarded


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the local effective interaction
   !---------------------------------------------------------------------------!
   subroutine calc_curlyU(curlyU,Wimp,Pimp,sym,curlyUcorr,mode)
      !
      use parameters
      use utils_fields
      use utils_misc
      use linalg, only : zeye, inv, inv_sym, diag, diagonal
      implicit none
      !
      type(BosonicField),intent(inout)      :: curlyU
      type(BosonicField),intent(in)         :: Wimp
      type(BosonicField),intent(in)         :: Pimp
      logical,intent(in),optional           :: sym
      type(BosonicField),intent(in),optional:: curlyUcorr
      character(len=*),intent(in),optional  :: mode
      !
      real(8),allocatable                   :: ReWimp(:,:),RePimp(:,:)
      complex(8),allocatable                :: invW(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nmats
      integer                               :: iw
      logical                               :: sym_,correctU,correctP
      character(len=10)                     :: mode_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_curlyU"
      !
      !
      ! Check on the input Fields
      if(.not.curlyU%status) stop "calc_curlyU: curlyU not properly initialized."
      if(.not.Wimp%status) stop "calc_curlyU: Wimp not properly initialized."
      if(.not.Pimp%status) stop "calc_curlyU: Pimp not properly initialized."
      if(curlyU%Nkpt.ne.0) stop "calc_curlyU: curlyU k dependent attributes are supposed to be unallocated."
      !if(Wimp%Nkpt.ne.0) stop "calc_curlyU: Wimp k dependent attributes are supposed to be unallocated."
      if(Pimp%Nkpt.ne.0) stop "calc_curlyU: Pimp k dependent attributes are supposed to be unallocated."
      !
      sym_=.true.
      if(present(sym))sym_=sym
      !
      Nbp = curlyU%Nbp
      Beta = curlyU%Beta
      Nmats = curlyU%Npoints
      !
      if(all([Wimp%Nbp-Nbp,Pimp%Nbp-Nbp].ne.[0,0])) stop "calc_curlyU: Either Wimp and/or Pimp have different orbital dimension with respect to curlyU."
      if(all([Wimp%Beta-Beta,Pimp%Beta-Beta].ne.[0d0,0d0])) stop "calc_curlyU: Either Wimp and/or Pimp have different Beta with respect to curlyU."
      if(all([Wimp%Npoints-Nmats,Pimp%Npoints-Nmats].ne.[0,0])) stop "calc_curlyU: Either Wimp and/or Pimp have different number of Matsubara points with respect to curlyU."
      !
      correctU=.false.
      correctP=.false.
      if(present(curlyUcorr))then
         if(.not.curlyUcorr%status) stop "calc_curlyU: Requested causality correction but curlyUcorr not properly initialized."
         if(Nbp.ne.curlyUcorr%Nbp) stop "calc_curlyU: curlyUcorr has different orbital dimension with respect to curlyU."
         if(Beta.ne.curlyUcorr%Beta) stop "calc_curlyU: curlyUcorr has different Beta with respect to curlyU."
         if(Nmats.ne.curlyUcorr%Npoints) stop "calc_curlyU: curlyUcorr has different number of Matsubara points with respect to curlyU."
         if(curlyUcorr%Nkpt.ne.0) stop "calc_curlyU: curlyUcorr k dependent attributes are supposed to be unallocated."
         mode_="curlyU"
         if(present(mode))mode_=reg(mode)
         select case(reg(mode_))
            case default
               stop "calc_curlyU: Available types are: curlyU, Ploc."
            case("curlyU")
               write(*,"(A)")"     Causality correction on curlyU."
               correctU=.true.
            case("Ploc")
               write(*,"(A)")"     Causality correction on Ploc."
               correctP=.true.
         end select
      endif
      !
      call clear_attributes(curlyU)
      !
      curlyU%bare_local = Wimp%bare_local
      !
      allocate(invW(Nbp,Nbp));invW=czero
      allocate(ReWimp(Nbp,Nbp));ReWimp=0d0
      allocate(RePimp(Nbp,Nbp));RePimp=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Wimp,Pimp,curlyU,correctP,curlyUcorr),&
      !$OMP PRIVATE(iw,invW,ReWimp,RePimp)
      !$OMP DO
      do iw=1,curlyU%Npoints
         !
         ReWimp = dreal(Wimp%screened_local(:,:,iw))
         RePimp = dreal(Pimp%screened_local(:,:,iw))
         !
         if(correctP)then
            invW = zeye(curlyU%Nbp) + matmul((RePimp - curlyUcorr%screened_local(:,:,iw)),ReWimp)
         else
            invW = zeye(curlyU%Nbp) + matmul(RePimp,ReWimp)
         endif
         !
         call inv(invW)
         curlyU%screened_local(:,:,iw) = matmul(ReWimp,invW)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(invW,ReWimp,RePimp)
      !
      !Causality correction on curlyU
      if(correctU)then
         do iw=1,curlyU%Npoints
            curlyU%screened_local(:,:,iw) = curlyU%screened_local(:,:,iw) - curlyUcorr%screened_local(:,:,iw)
         enddo
      endif
      call isReal(curlyU)
      !
      !Check if curlyU is locally symmetric - print if relative error is bigger than 1e-3
      if(sym_)then
         write(*,"(A)") "     Checking symmetry of curlyU (enforced)."
         do iw=1,Nmats
            call check_Symmetry(curlyU%screened_local(:,:,iw),1e7*eps,enforce=.true.,hardstop=.false.,name="curlyU_w"//str(iw),verb=.true.)
         enddo
      endif
      !
   end subroutine calc_curlyU


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the fully screened local interaction
   !---------------------------------------------------------------------------!
   subroutine calc_Wimp(Wimp,curlyU,ChiC,sym)
      !
      use parameters
      use utils_fields
      use utils_misc
      use linalg, only : zeye
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wimp
      type(BosonicField),intent(in)         :: curlyU
      type(BosonicField),intent(in)         :: ChiC
      logical,intent(in),optional           :: sym
      !
      complex(8),allocatable                :: Wtmp(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nmats
      integer                               :: iw
      logical                               :: sym_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Wimp"
      !
      !
      ! Check on the input Fields
      if(.not.Wimp%status) stop "calc_Wimp: Wimp not properly initialized."
      if(.not.curlyU%status) stop "calc_Wimp: curlyU not properly initialized."
      if(.not.ChiC%status) stop "calc_Wimp: ChiC not properly initialized."
      if(Wimp%Nkpt.ne.0) stop "calc_Wimp: Wimp k dependent attributes are supposed to be unallocated."
      if(curlyU%Nkpt.ne.0) stop "calc_Wimp: curlyU k dependent attributes are supposed to be unallocated."
      if(ChiC%Nkpt.ne.0) stop "calc_Wimp: ChiC k dependent attributes are supposed to be unallocated."
      if(allocated(ChiC%bare_local))  stop "calc_Wimp: ChiC bare_local attribute is supposed to be unallocated."
      if(allocated(ChiC%bare))  stop "calc_Wimp: ChiC bare attribute is supposed to be unallocated."
      !
      sym_=.true.
      if(present(sym))sym_=sym
      !
      Nbp = Wimp%Nbp
      Beta = Wimp%Beta
      Nmats = Wimp%Npoints
      !
      if(all([curlyU%Nbp-Nbp,ChiC%Nbp-Nbp].ne.[0,0])) stop "calc_Wimp: Either curlyU and/or ChiC have different orbital dimension with respect to Wimp."
      if(all([curlyU%Beta-Beta,ChiC%Beta-Beta].ne.[0d0,0d0])) stop "calc_Wimp: Either curlyU and/or ChiC have different Beta with respect to Wimp."
      if(all([curlyU%Npoints-Nmats,ChiC%Npoints-Nmats].ne.[0,0]))   stop "calc_Wimp: Either curlyU and/or ChiC have different number of Matsubara points with respect to Wimp."
      !
      call clear_attributes(Wimp)
      !
      Wimp%bare_local = curlyU%bare_local
      !
      allocate(Wtmp(Nbp,Nbp));Wtmp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Wimp,ChiC,curlyU),&
      !$OMP PRIVATE(iw,Wtmp)
      !$OMP DO
      do iw=1,Wimp%Npoints
         !
         Wtmp = matmul(ChiC%screened_local(:,:,iw),curlyU%screened_local(:,:,iw))
         Wimp%screened_local(:,:,iw) = curlyU%screened_local(:,:,iw) - matmul(curlyU%screened_local(:,:,iw),Wtmp)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Wtmp)
      call isReal(Wimp)
      !
      !Check if Wimp is locally symmetric - print if relative error is bigger than 1e-3
      if(sym_)then
         write(*,"(A)") "     Checking symmetry of Wimp (enforced)."
         do iw=1,Nmats
            call check_Symmetry(Wimp%screened_local(:,:,iw),1e7*eps,enforce=.true.,hardstop=.false.,name="Wimp_w"//str(iw),verb=.true.)
         enddo
      endif
      !
   end subroutine calc_Wimp


   !---------------------------------------------------------------------------!
   !PURPOSE: Correct the RPA interaction at gamma (not used)
   !---------------------------------------------------------------------------!
   subroutine correct_Gamma(U,mode)
      !
      use parameters
      use utils_misc
      use crystal
      implicit none
      !
      type(BosonicField),intent(inout)      :: U
      integer,intent(in)                    :: mode
      !
      integer                               :: Nkpt,Nbp,ib1,ib2
      integer                               :: iangle,Nangles,DimK,Navg,ik0,ik1,ik2
      integer                               :: shift,cycleCondition,idist,DistLoop
      real(8)                               :: k0,k1,k2,W1,W2
      real(8)                               :: U_ik_fact,U_gamma_fact
      logical,allocatable                   :: condition(:)
      integer,allocatable                   :: ListN(:),Korder(:)
      real(8),allocatable                   :: ListK(:),ListTheta(:),ListPhi(:)
      real(8),allocatable                   :: U_Gamma(:,:)
      type(physicalU)                       :: PhysicalUelements
      !
      !
      verbose=.true.
      if(verbose)write(*,"(A)") "---- correct_Gamma"
      !
      !
      ! Check on the input Fields
      if(.not.U%status) stop "correct_Gamma: U not properly initialized."
      if(U%Nkpt.eq.0) stop "correct_Gamma: U k dependent attributes are supposed to be allocated."
      if(U%iq_gamma.lt.0) stop "correct_Gamma: U iq_gamma not defined."
      if(.not.small_ik_stored) stop "correct_Gamma: Small k-vectors not stored."
      if(size(small_ik,dim=1)+1.ne.U%Nkpt) stop "correct_Gamma: Provided K mesh does not match the U one."
      !
      Nkpt = U%Nkpt
      Nbp = U%Nbp
      !
      call init_Uelements(int(sqrt(dble(Nbp))),PhysicalUelements)
      !
      select case(mode)
         case default
            !
            stop "correct_Gamma: Available smearing modes: 1 2 3."
            !
         case(1)
            !
            write(*,"(A)") "     Gaussian interpolation of Gamma."
            shift = 0            !points used to fit Gamma = [D1,D2]
            cycleCondition = 2   !list of agular dependend distances from Gamma must contain at least 2 elements
            DistLoop = 1         !loop over the nearest neighbors
            !
         case(2)
            !
            write(*,"(A)") "     Gaussian interpolation of Gamma and Gamma_NN."
            shift = 1            !points used to fit Gamma and D1 = [D2,D3]
            cycleCondition = 3   !list of agular dependend distances from Gamma must contain at least 3 elements
            DistLoop = 1         !loop over the nearest neighbors
            !
         case(3)
            !
            write(*,"(A)") "     Gaussian interpolation of Gamma, Gamma_NN and Gamma_NNN."
            shift = 1            !points used to fit Gamma and D1 = [D2,D3]
            cycleCondition = 3   !list of agular dependend distances from Gamma must contain at least 3 elements
            DistLoop = 2         !loop over the nearest and next-nearest neighbors
            !
      end select
      !
      U%screened(:,:,:,U%iq_gamma) = czero
      !
      Navg=0
      allocate(U_Gamma(Nbp,Nbp));U_Gamma=0d0
      allocate(condition(Nkpt-1));condition=.false.
      do idist=1,DistLoop
         !
         if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"     Scanning the K points at distance "//str(idist)//" from Gamma."
         !
         !find the K-vectors closest to Gamma and for each of them the {theta,phi} pair
         condition = small_ik(:,2) .eq. idist
         Nangles = size( pack(KvecPolar(:,1) , condition) )
         ListTheta = pack( KvecPolar(:,2) , condition )
         ListPhi = pack( KvecPolar(:,3) , condition )
         !
         !Gaussian interpolation of the Gamma and closest K-points
         condition=.false.
         do iangle=1,Nangles
            !
            !for each of the K-vectors closest to Gamma find all the others with the same {theta,phi} pair
            condition = ( KvecPolar(:,2).eq.ListTheta(iangle) ) .and. ( KvecPolar(:,3).eq.ListPhi(iangle) )
            ListK = pack( KvecPolar(:,1) , condition )
            ListN = pack( small_ik(:,1)  , condition )
            if(size(ListK).ne.size(ListN)) stop "correct_Gamma: length of K mesh for interpolation does not correspond to length of indexes array."
            DimK = size(ListK)
            !
            !I need al least two points to interpolate a gaussian
            if(DimK.lt.cycleCondition)cycle
            Navg = Navg + 1
            !
            if(verbose)then
               write(*,"(2(A,I4))") "     => iangle: ",iangle,"   iavg: ",Navg
               write(*,"(2(A,1F6.2),A,"//str(DimK)//"I6)") "        theta: ", ListTheta(iangle), "   phi: ", ListPhi(iangle),"   ik: ",ListN
               write(*,"(A,"//str(DimK)//"F8.4)") "        unsorted K: ",ListK
            endif
            !
            !Sort the K coordinate
            allocate(Korder(DimK));Korder=0
            call sort_array(ListK,Korder)
            !
            !Interpolate a Gaussian from the two K-points [D1,D2] or [D2,D3] depending on "shift" and always average the Gamma point
            k1 = ListK(Korder(1+shift)); ik1 = ListN(Korder(1+shift))
            k2 = ListK(Korder(2+shift)); ik2 = ListN(Korder(2+shift))
            if(shift.ne.0)then
               k0 = ListK(Korder(1))
               ik0 = ListN(Korder(1))
               if((k0.eq.k1).or.(k0.eq.k2)) stop "correct_Gamma: interpolated K-point correspond to interpolation one."
               if((ik0.eq.ik1).or.(ik0.eq.ik2)) stop "correct_Gamma: interpolated K-point index correspond to interpolation one."
            endif
            !
            if(verbose)then
               write(*,"(A,"//str(DimK)//"F8.4)") "        sorted K: ",ListK(Korder)
               write(*,"(2(A,I4),2(A,F8.4))") "        ik1_fit: ",ik1,"   ik2_fit: ",ik2,"   k1_fit: ",k1,"   k2_fit: ",k2
               if(shift.eq.0)then
                  write(*,"(A)") "        Adding only Gamma."
               else
                  write(*,"(A,I)") "        Adding Gamma replacing ik:",ik0
               endif
            endif
            !
            do ib1=1,Nbp
               do ib2=1,Nbp
                  !
                  if(.not.PhysicalUelements%Full_All(ib1,ib2)) cycle
                  !
                  W1 = dreal(U%bare(ib1,ib2,ik1))
                  W2 = dreal(U%bare(ib1,ib2,ik2))
                  !
                  U_Gamma(ib1,ib2) = U_Gamma(ib1,ib2) + interp_Lorentzian([k1,W1],[k2,W2],0d0)
                  if(shift.ne.0)then
                     !
                     U_ik_fact = interp_Lorentzian([k1,W1],[k2,W2],k0) / dreal(U%bare(ib1,ib2,ik0))
                     !
                     U%bare(ib1,ib2,ik0) = U%bare(ib1,ib2,ik0) * U_ik_fact
                     U%screened(ib1,ib2,:,ik0) = U%screened(ib1,ib2,:,ik0) * U_ik_fact
                     !
                  endif
                  !
               enddo
            enddo
            !
            deallocate(condition,ListK,ListN,Korder)
            !
         enddo !iangle
      enddo !idist
      U_Gamma = U_Gamma / Navg
      !
      if(verbose) write(*,"(A,I)")new_line("A")//"     Navg: ",Navg
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            if(.not.PhysicalUelements%Full_All(ib1,ib2)) cycle
            !
            U_gamma_fact = U_Gamma(ib1,ib2) / dreal(U%bare(ib1,ib2,U%iq_gamma))
            !
            U%bare(ib1,ib2,U%iq_gamma) = U%bare(ib1,ib2,U%iq_gamma) * U_gamma_fact
            U%screened(ib1,ib2,:,U%iq_gamma) = U%screened(ib1,ib2,:,U%iq_gamma) * U_gamma_fact
            !
         enddo
      enddo
      deallocate(U_Gamma)
      !
      !
      !
      contains
      !
      !
      !
      function interp_Gaussian(P1,P2,x) result(y)
         implicit none
         real(8),intent(in)                 :: P1(2)
         real(8),intent(in)                 :: P2(2)
         real(8),intent(in)                 :: x
         real(8)                            :: y
         real(8)                            :: A,B
         real(8)                            :: x1,y1,x2,y2
         !
         x1 = P1(1); y1 = P1(2)
         x2 = P2(1); y2 = P2(2)
         !
         B = (x1**2 - x2**2) / log(y2/y1)
         A = y1*exp((x1**2)/B)
         !
         y = A* exp(-(x**2)/B)
         !
      end function interp_Gaussian
      !
      function interp_Lorentzian(P1,P2,x) result(y)
         implicit none
         real(8),intent(in)                 :: P1(2)
         real(8),intent(in)                 :: P2(2)
         real(8),intent(in)                 :: x
         real(8)                            :: y
         real(8)                            :: A,B
         real(8)                            :: x1,y1,x2,y2
         !
         x1 = P1(1); y1 = P1(2)
         x2 = P2(1); y2 = P2(2)
         !
         B = abs( ( y2*(x2**2) - y1*(x1**2) ) / ( y1 - y2 ) )
         A = y1 * ( (x1**2) + B ) / B
         !
         y = A*B / ( x**2 + B )
         !
      end function interp_Lorentzian
      !
      !
      !
   end subroutine correct_Gamma


end module interactions
