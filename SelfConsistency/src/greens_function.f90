module greens_function

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface calc_density
      module procedure calc_density_loc                                         !(Gmats,n_loc[Norb,Norb,Nspin])
      module procedure calc_density_Kdep                                        !(Gmats,Lattice,n_k[Norb,Norb,Nkpt,Nspin])
   end interface calc_density

   interface set_density
      module procedure set_density_Int                                          !(Gmats,Lattice,musearch,Smats(optional))
      module procedure set_density_NonInt                                       !(mu,Beta,Lattice,musearch)
   end interface set_density

   interface calc_Gmats
      module procedure calc_Gmats_Full                                          !(Gmats,Lttc,Smats)
      module procedure calc_Gmats_Shift                                         !(Gmats,Lttc)
   end interface calc_Gmats

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
   public :: calc_density
   public :: set_density
   public :: calc_Gmats
   public :: calc_Glda
   public :: calc_G0_tau

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the density given the Green's function
   !---------------------------------------------------------------------------!
   subroutine calc_density_loc(Gmats,n_loc)
      !
      use parameters
      use utils_misc
      use fourier_transforms
      use input_vars, only : Ntau, tau_uniform, paramagnet
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats
      complex(8),allocatable,intent(inout)  :: n_loc(:,:,:)
      !
      complex(8),allocatable                :: Gitau(:,:,:)
      real(8)                               :: Beta
      integer                               :: ispin,Norb
      !
      !
      if(verbose)write(*,"(A)") "---- calc_density_loc"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "calc_density_loc: Gmats not properly initialized."
      if(Gmats%Npoints.eq.0) stop "calc_density_loc: Gmats frequency dependent attributes not properly initialized."
      !
      Norb = Gmats%Norb
      Beta = Gmats%Beta
      !
      call assert_shape(n_loc,[Norb,Norb,Nspin],"calc_density_loc","n_loc")
      !
      n_loc=czero
      allocate(Gitau(Norb,Norb,Ntau));Gitau=czero
      spinloop: do ispin=1,Nspin
         !
         call Fmats2itau_mat(Beta,Gmats%ws(:,:,:,ispin),Gitau, &
         asympt_corr=.true.,tau_uniform=tau_uniform,atBeta=.true.)
         !
         n_loc(:,:,ispin) = -Gitau(:,:,Ntau)
         if(paramagnet)then
            n_loc(:,:,Nspin) = n_loc(:,:,1)
            exit spinloop
         endif
         !
      enddo spinloop
      deallocate(Gitau)
      !
   end subroutine calc_density_loc
   !
   subroutine calc_density_Kdep(Gmats,Lttc,n_k,along_path,along_plane)
      !
      use parameters
      use utils_misc
      use fourier_transforms
      use input_vars, only : Ntau, tau_uniform, cmplxWann, paramagnet
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      complex(8),allocatable,intent(inout)  :: n_k(:,:,:,:)
      logical,intent(in),optional           :: along_path
      logical,intent(in),optional           :: along_plane
      !
      complex(8),allocatable                :: Gitau(:,:,:,:)
      real(8)                               :: Beta
      integer                               :: ispin,Norb,Nkpt
      logical                               :: along_path_,along_plane_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_density_Kdep"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "calc_density_Kdep: Gmats not properly initialized."
      if(Gmats%Npoints.eq.0) stop "calc_density_Kdep: Gmats frequency dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_density_Kdep: Gmats k dependent attributes not properly initialized."
      if(Gmats%Norb.ne.Lttc%Norb) stop "calc_density_Kdep: Lttc has different number of orbitals with respect to Gmats."
      !
      along_path_=.false.
      if(present(along_path))along_path_=along_path
      along_plane_=.false.
      if(present(along_plane))along_plane_=along_plane
      !
      if(along_path_.and.along_plane_) stop "calc_density_Kdep: cannot have along_path=T and along_plane=T simultaneously."
      !
      if(along_path_)then
         if(Gmats%Nkpt.ne.Lttc%Nkpt_path) stop "calc_density_Kdep: Lttc Kpath has different number of path k-points with respect to Gmats."
         if(.not.allocated(Lttc%kptpath)) stop "calc_density_Kdep: K-point path not allocated."
      elseif(along_plane_)then
         if(Gmats%Nkpt.ne.Lttc%Nkpt_Plane) stop "calc_density_Kdep: Lttc Kplane has different number of path k-points with respect to Gmats."
         if(.not.allocated(Lttc%kptPlane)) stop "calc_density_Kdep: K-point on plane not allocated."
      else
         if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "calc_density_Kdep: Lttc has different number of k-points with respect to Gmats."
      endif
      !
      Norb = Gmats%Norb
      Nkpt = Gmats%Nkpt
      Beta = Gmats%Beta
      !
      if(verbose)write(*,"(A,I)") "     Number of K-points:",Nkpt
      call assert_shape(n_k,[Norb,Norb,Nkpt,Nspin],"calc_density_Kdep","n_k")
      !
      n_k=czero
      allocate(Gitau(Norb,Norb,Ntau,Nkpt));Gitau=czero
      spinloop: do ispin=1,Nspin
         !
         if(cmplxWann.or.along_path_.or.along_plane_)then
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau, &
            asympt_corr=.true.,tau_uniform=tau_uniform,atBeta=.true.)
         else
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau, &
            asympt_corr=.true.,tau_uniform=tau_uniform,Nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt,atBeta=.true.)
         endif
         !
         n_k(:,:,:,ispin) = -Gitau(:,:,Ntau,:)
         if(paramagnet)then
            n_k(:,:,:,Nspin) = n_k(:,:,:,1)
            exit spinloop
         endif
         !
      enddo spinloop
      deallocate(Gitau)
      !
   end subroutine calc_density_Kdep


   !---------------------------------------------------------------------------!
   !PURPOSE: analytically compute the G0(\tau) in the diagonal basis.
   !by now only for paramagnetic G
   !---------------------------------------------------------------------------!
   subroutine calc_G0_tau(Gitau,mu,Beta,Ek,atBeta)
      !
      use parameters
      use utils_misc
      use fourier_transforms
      use input_vars, only : tau_uniform
      implicit none
      !
      complex(8),allocatable,intent(inout)  :: Gitau(:,:,:)
      real(8),intent(in)                    :: mu
      real(8),intent(in)                    :: Beta
      real(8),allocatable,intent(in)        :: Ek(:,:)
      logical,intent(in),optional           :: atBeta
      !
      real(8),allocatable                   :: tau(:)
      real(8)                               :: eu,fermicut
      integer                               :: iwan,ik,itau,Ntau
      integer                               :: Norb,Nkpt
      logical                               :: upper,lower,atBeta_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_G0_tau"
      fermicut=log(huge(1.0d0)-1e2)/2.d0
      !
      !
      ! Check on the input Fields
      Norb = size(Ek,dim=1)
      Nkpt = size(Ek,dim=2)
      Ntau = size(Gitau,dim=2)
      call assert_shape(Gitau,[Norb,Ntau,Nkpt],"calc_G0_tau","Gitau")
      atBeta_ = .false.
      if(present(atBeta)) atBeta_ = atBeta
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform)then
         tau = linspace(0d0,Beta,Ntau)
      else
         tau = denspace(Beta,Ntau)
      endif
      !
      Gitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Ntau,Norb,atBeta_,tau,Ek,mu,fermicut,Beta,Gitau),&      !VEDI SE PUOI TOGLIERE QUALCHE IF
      !$OMP PRIVATE(ik,itau,iwan,eu,upper,lower)
      !$OMP DO
      do itau=1,Ntau
         if(atBeta_.and.(itau.ne.Ntau))cycle
         do ik=1,Nkpt
            do iwan=1,Norb
               !
               eu = ( Ek(iwan,ik) - mu ) !WHY factor 0.5?
               !
               upper = (eu.gt.0d0).and.(eu*tau(itau).lt.+fermicut)
               lower = (eu.lt.0d0).and.(eu*(Beta-tau(itau)).gt.-fermicut)
               !
               if(upper) Gitau(iwan,itau,ik) = -(1.d0 - fermidirac(Ek(iwan,ik),mu,Beta)) * dexp(-eu*tau(itau))
               if(lower) Gitau(iwan,itau,ik) = -fermidirac(Ek(iwan,ik),mu,Beta)* dexp(eu*(beta-tau(itau)))
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(tau)
      !
   end subroutine calc_G0_tau


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the Matsubara Green's Function
   !---------------------------------------------------------------------------!
   subroutine calc_Gmats_Full(Gmats,Lttc,Smats,along_path,along_plane)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use input_vars, only : paramagnet, Hetero
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gmats
      type(Lattice),intent(in),target       :: Lttc
      type(FermionicField),intent(in),optional :: Smats !target
      logical,intent(in),optional           :: along_path
      logical,intent(in),optional           :: along_plane
      !
      complex(8),allocatable                :: invGf(:,:)
      !complex(8),pointer                    :: Swks(:,:,:,:,:)
      complex(8),pointer                    :: Hk(:,:,:)
      complex(8),allocatable                :: zeta(:,:,:)
      complex(8),allocatable                :: n_k(:,:,:,:)
      complex(8),allocatable                :: Potential_L(:,:,:,:,:)
      complex(8),allocatable                :: Potential_R(:,:,:,:,:)
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: Beta,mu
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iw,ik,iwan,ispin
      integer                               :: Ln(2),Rn(2),NbulkL,NbulkR
      logical                               :: along_path_,along_plane_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Gmats_Full"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "calc_Gmats_Full: Gmats not properly initialized."
      if(.not.Lttc%status) stop "calc_Gmats_Full: Lttc not properly initialized."
      if(Gmats%Npoints.eq.0) stop "calc_Gmats_Full: Gmats frequency dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_Gmats_Full: Gmats k dependent attributes not properly initialized."
      if(Gmats%Norb.ne.Lttc%Norb) stop "calc_Gmats_Full: Lttc has different number of orbitals with respect to Gmats."
      !
      along_path_=.false.
      if(present(along_path))along_path_=along_path
      along_plane_=.false.
      if(present(along_plane))along_plane_=along_plane
      !
      if(along_path_)then
         Hk => Lttc%Hk_path
         if(Gmats%Nkpt.ne.Lttc%Nkpt_path) stop "calc_Gmats_Full: Lttc Kpath has different number of path k-points with respect to Gmats."
         if(Gmats%Nkpt.ne.size(Hk,dim=3)) stop "calc_Gmats_Full: Lttc Hk_path has different number of path k-points with respect to Gmats."
         if(.not.allocated(Lttc%Hk_path)) stop "calc_Gmats_Full: H(k) along path not allocated."
         if(verbose)write(*,"(A)") "     H(k) is along a path."
      elseif(along_plane_)then
         Hk => Lttc%Hk_Plane
         if(Gmats%Nkpt.ne.Lttc%Nkpt_Plane) stop "calc_Gmats_Full: Lttc Kplane has different number of path k-points with respect to Gmats."
         if(Gmats%Nkpt.ne.size(Hk,dim=3)) stop "calc_Gmats_Full: Lttc Hk_Plane has different number of path k-points with respect to Gmats."
         if(.not.allocated(Lttc%Hk_Plane)) stop "calc_Gmats_Full: H(k) in plane not allocated."
         if(verbose)write(*,"(A)") "     H(k) is within a plane."
      else
         Hk => Lttc%Hk
         if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "calc_Gmats_Full: Lttc has different number of k-points with respect to Gmats."
         if(Gmats%Nkpt.ne.size(Hk,dim=3)) stop "calc_Gmats_Full: Lttc Hk has different number of path k-points with respect to Gmats."
         if(verbose)write(*,"(A)") "     H(k) is within the full BZ."
      endif
      !
      Norb = Gmats%Norb
      Nmats = Gmats%Npoints
      Nkpt = Gmats%Nkpt
      Beta = Gmats%Beta
      mu = Gmats%mu
      !
      if(verbose)write(*,"(A,I)") "     Number of K-points:",Nkpt
      !
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(Beta,Nmats)
      allocate(zeta(Norb,Norb,Nmats));zeta=czero
      do iwan=1,Norb
         do iw=1,Nmats
            zeta(iwan,iwan,iw) = dcmplx( mu , wmats(iw) )
         enddo
      enddo
      deallocate(wmats)
      !
      if(present(Smats))then
         if(.not.Smats%status) stop "calc_Gmats_Full: Smats not properly initialized."
         if(Smats%Npoints.ne.Nmats) stop "calc_Gmats_Full: Smats has different number of Matsubara points with respect to Gmats."
         if(Smats%Nkpt.ne.Nkpt) stop "calc_Gmats_Full: Smats has different number of k-points with respect to Gmats."
         if(verbose)write(*,"(A)") "     Interacting Green's function."
         !
         Ln=0;Rn=0
         if(Hetero%status)then
            if(Hetero%Explicit(1).ne.1)then
               !
               Ln(1) = 1
               Ln(2) = Hetero%Norb
               NbulkL = Hetero%Explicit(1)-1
               !
               allocate(Potential_L(Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin));Potential_L=czero
               call build_Potential(Potential_L,Hetero%tz(:,:,Hetero%Explicit(1)-1),NbulkL &
                                               ,zeta(Ln(1):Ln(2),Ln(1):Ln(2),:)            &
                                               ,Hk(Ln(1):Ln(2),Ln(1):Ln(2),:)              &
                                               ,Smats%wks(Ln(1):Ln(2),Ln(1):Ln(2),:,:,:) )
               write(*,"(2(A,2I4))") "     Left potential orbital indexes: ",Ln(1),Ln(2)," thickness: ",NbulkL
               !
            endif
            if(Hetero%Explicit(2).ne.Hetero%Nslab)then
               !
               Rn(1) = 1+ Norb - Hetero%Norb
               Rn(2) = Norb
               NbulkR = Hetero%Nslab-Hetero%Explicit(2)
               !
               allocate(Potential_R(Hetero%Norb,Hetero%Norb,Nmats,Nkpt,Nspin));Potential_R=czero
               call build_Potential(Potential_R,Hetero%tz(:,:,Hetero%Explicit(2)),NbulkR   &
                                               ,zeta(Rn(1):Rn(2),Rn(1):Rn(2),:)            &
                                               ,Hk(Rn(1):Rn(2),Rn(1):Rn(2),:)              &
                                               ,Smats%wks(Rn(1):Rn(2),Rn(1):Rn(2),:,:,:) )
               write(*,"(2(A,2I4))") "     Right potential orbital indexes: ",Rn(1),Rn(2)," thickness: ",NbulkR
               !
            endif
            !
         endif
         !
      else
         if(verbose)write(*,"(A)") "     LDA Green's function."
      endif
      !
      call clear_attributes(Gmats)
      allocate(invGf(Norb,Norb));invGf=czero
      spinloop: do ispin=1,Nspin
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(ispin,zeta,Gmats,Smats,Hk,Ln,Potential_L,Rn,Potential_R),&
         !$OMP PRIVATE(ik,iw,invGf)
         !$OMP DO
         do iw=1,Gmats%Npoints
            do ik=1,Gmats%Nkpt
               !
               invGf = zeta(:,:,iw) - Hk(:,:,ik)
               !
               if(present(Smats)) invGf = invGf - Smats%wks(:,:,iw,ik,ispin)
               if(allocated(Potential_L)) invGf(Ln(1):Ln(2),Ln(1):Ln(2)) = invGf(Ln(1):Ln(2),Ln(1):Ln(2)) - Potential_L(:,:,iw,ik,ispin)
               if(allocated(Potential_R)) invGf(Rn(1):Rn(2),Rn(1):Rn(2)) = invGf(Rn(1):Rn(2),Rn(1):Rn(2)) - Potential_R(:,:,iw,ik,ispin)
               !
               call inv(invGf)
               Gmats%wks(:,:,iw,ik,ispin) = invGf
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         if(paramagnet)then
            Gmats%wks(:,:,:,:,Nspin) = Gmats%wks(:,:,:,:,1)
            exit spinloop
         endif
      enddo spinloop
      deallocate(zeta,invGf)
      !if(associated(Swks))nullify(Swks)
      if(associated(Hk))nullify(Hk)
      !
      ! In the N_ks attribute is stored the k-dep occupation
      allocate(n_k(Norb,Norb,Nkpt,Nspin));n_k=czero
      if(along_path_)then
         call calc_density(Gmats,Lttc,n_k,along_path=along_path_)
      elseif((.not.along_path_).and.(.not.along_plane_))then
         call calc_density(Gmats,Lttc,n_k)
      endif
      Gmats%N_ks = n_k
      deallocate(n_k)
      !
      ! Fill the local attributes
      if(.not.(along_path_.or.along_plane_))then
         call FermionicKsum(Gmats)
         !this is more accurate than summing G(beta) at each k-point
         call calc_density_loc(Gmats,Gmats%N_s)
      endif
      !
   end subroutine calc_Gmats_Full


   !---------------------------------------------------------------------------!
   !PURPOSE: Recalculate the Matsubara Green's Function assuming that the
   !         Gmats%mu attribute has changed
   !---------------------------------------------------------------------------!
   subroutine calc_Gmats_Shift(Gmats,Lttc,mu_shift)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use input_vars, only : paramagnet
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gmats
      real(8),intent(in)                    :: mu_shift
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: invGf(:,:)
      real(8),allocatable                   :: zeta(:,:)
      complex(8),allocatable                :: n_k(:,:,:,:)
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iw,ik,ispin
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Gmats_Shift"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "calc_Gmats_Shift: Gmats not properly initialized."
      if(.not.Lttc%status) stop "calc_Gmats_Shift: Lttc not properly initialized."
      if(Gmats%Npoints.eq.0) stop "calc_Gmats_Shift: Gmats frequency dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_Gmats_Shift: Gmats k dependent attributes not properly initialized."
      if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "calc_Gmats_Shift: Lttc has different number of k-points with respect to Gmats."
      if(Gmats%Norb.ne.Lttc%Norb) stop "calc_Gmats_Shift: Lttc has different number of orbitals with respect to Gmats."
      if(mu_shift.eq.0d0) stop "calc_Gmats_Shift: Chemical potential shift is zero."
      !
      Norb = Gmats%Norb
      Nmats = Gmats%Npoints
      Nkpt = Gmats%Nkpt
      !
      if(verbose)write(*,"(A,I)") "     Number of K-points:",Nkpt
      !
      allocate(zeta(Norb,Norb));zeta=0d0
      zeta = deye(Norb)*mu_shift
      !
      allocate(invGf(Norb,Norb));invGf=czero
      spinloop: do ispin=1,Nspin
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(ispin,zeta,Gmats),&
         !$OMP PRIVATE(ik,iw,invGf)
         !$OMP DO
         do ik=1,Gmats%Nkpt
            do iw=1,Gmats%Npoints
               !
               invGf = Gmats%wks(:,:,iw,ik,ispin)
               call inv(invGf)
               !
               invGf = invGf + zeta
               !
               call inv(invGf)
               Gmats%wks(:,:,iw,ik,ispin) = invGf
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         if(paramagnet)then
            Gmats%wks(:,:,:,:,Nspin) = Gmats%wks(:,:,:,:,1)
            exit spinloop
         endif
      enddo spinloop
      deallocate(zeta,invGf)
      !
      ! In the N_ks attribute is stored the k-dep occupation
      allocate(n_k(Norb,Norb,Nkpt,Nspin));n_k=czero
      call calc_density(Gmats,Lttc,n_k)
      Gmats%N_ks = n_k
      deallocate(n_k)
      !
      ! Fill the local attributes
      call FermionicKsum(Gmats)
      !this is more accurate than summing G(beta) at each k-point
      call calc_density_loc(Gmats,Gmats%N_s)
      !
   end subroutine calc_Gmats_Shift


   !---------------------------------------------------------------------------!
   !PURPOSE: Set the density of an interacting Green's function. If the
   !         self-energy is provided the Gf is fully recalculated (v1),
   !         otherwise only a rigid shift will be applied to the inverse (v2)
   !---------------------------------------------------------------------------!
   subroutine set_density_Int(Gmats,Lttc,mu_param,Smats)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gmats
      type(Lattice),intent(in)              :: Lttc
      type(musearch),intent(in)             :: mu_param
      type(FermionicField),intent(in),optional :: Smats
      !
      real(8)                               :: n_iter,n_err,dmu
      real(8)                               :: mu_start,mu_last,mu_sign
      real(8)                               :: mu_below,mu_above
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iter
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- set_density_Int"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "set_density_Int: Gmats not properly initialized."
      if(.not.Lttc%status) stop "set_density_Int: Lttc not properly initialized."
      if(present(Smats).and.(.not.Smats%status)) stop "set_density_Int: Smats not properly initialized."
      if(Gmats%Npoints.eq.0) stop "set_density_Int: Gmats frequency dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "set_density_Int: Gmats k dependent attributes not properly initialized."
      if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "set_density_Int: Lttc has different number of k-points with respect to Gmats."
      if(Gmats%Norb.ne.Lttc%Norb) stop "set_density_Int: Lttc has different number of orbitals with respect to Gmats."
      if(mu_param%TargetDensity.eq.0d0) stop "set_density_Int: TargetDensity is set to zero."
      if(present(Smats))then
         if(Smats%Npoints.eq.0) stop "set_density_Int: Smats frequency dependent attributes not properly initialized."
         if(Smats%Nkpt.eq.0) stop "set_density_Int: Smats k dependent attributes not properly initialized."
      endif
      write(*,"(A,F6.3)") "     Target density: ",mu_param%TargetDensity
      !
      Norb = Gmats%Norb
      Nmats = Gmats%Npoints
      Nkpt = Gmats%Nkpt
      mu_start = Gmats%mu
      !
      if(present(Smats)) call calc_Gmats_Full(Gmats,Lttc,Smats)
      n_iter = get_dens()
      write(*,"(2(A,F12.5))") "     Starting density: ",n_iter,", starting mu: ",mu_start
      !
      n_err = abs(n_iter-mu_param%TargetDensity)/mu_param%TargetDensity
      if(n_err.lt.mu_param%densityRelErr)then
         write(*,"(A)") "     Density already under threshold, skipping mu scan."
         return
      endif
      !
      mu_sign = sign(1d0,mu_param%TargetDensity-n_iter)
      mu_last = mu_start
      !
      do iter=1,mu_param%muIter
         !
         Gmats%mu = mu_start + iter * mu_param%muStep * mu_sign;
         dmu = Gmats%mu - mu_last;
         !
         if(present(Smats))then
            call calc_Gmats_Full(Gmats,Lttc,Smats)
         else
            call calc_Gmats_Shift(Gmats,Lttc,dmu)
         endif
         n_iter = get_dens()
         !
         write(*,"(A,I4,3(A,F12.5))") "     Rigid shift it# ",iter,", density: ", n_iter,", mu: ", &
         Gmats%mu,", Dmu: ",dmu
         !
         if((mu_sign.gt.0.0).and.(n_iter > mu_param%TargetDensity)) exit
         if((mu_sign.lt.0.0).and.(n_iter < mu_param%TargetDensity)) exit
         !
         mu_last = Gmats%mu
         !
      enddo
      !
      mu_below = min(Gmats%mu,mu_last)
      mu_above = max(Gmats%mu,mu_last)
      !
      do iter=1,mu_param%muIter
         !
         mu_last = Gmats%mu
         Gmats%mu = (mu_below+mu_above)/2d0
         dmu = Gmats%mu - mu_last
         !
         if(present(Smats))then
            call calc_Gmats_Full(Gmats,Lttc,Smats)
         else
            call calc_Gmats_Shift(Gmats,Lttc,dmu)
         endif
         n_iter = get_dens()
         n_err = abs(n_iter-mu_param%TargetDensity)/mu_param%TargetDensity
         !
         write(*,"(A,I4,6(A,F12.5))") "     Netwon it# ",iter,", mu boundaries: ",mu_below," / ",mu_above, &
         ", dmu: ",dmu,", density: ", n_iter,", mu: ",Gmats%mu,", relative error: ",n_err
         !
         if(n_iter.gt.mu_param%TargetDensity) mu_above = Gmats%mu;
         if(n_iter.lt.mu_param%TargetDensity) mu_below = Gmats%mu;
         !
         if(n_err.lt.mu_param%densityRelErr)then
            write(*,"(A,F12.5)") "     Found correct chemical potential after "//str(iter)//" iterations: ",Gmats%mu
            exit
         elseif(iter.eq.mu_param%muIter)then
            write(*,"(A,F12.5)") "     Warning: NOT found correct chemical potential after "//str(iter)//" iterations. Last used value: ",Gmats%mu
         endif
         !
      enddo !iter
      !
   contains
      !
      !
      !
      function get_dens() result(dens)
         implicit none
         complex(8)              :: dens_C
         real(8)                 :: dens
         integer                 :: iwan,ispin
         !
         dens_C = czero
         if(allocated(mu_param%orbs))then
            do iwan=1,Norb
               if(all(mu_param%orbs.ne.iwan))cycle
               do ispin=1,Nspin
                  dens_C = dens_C + Gmats%N_s(iwan,iwan,ispin)
               enddo
            enddo
         else
            dens_C = trace(sum(Gmats%N_s,dim=3))
         endif
         !
         if(aimag(dens_C).gt.eps)then
            write(*,"(A,2F10.5)")"     Density (Int) is complex: ",real(dens_C),aimag(dens_C)
            stop "get_dens(set_density_Int)"
         endif
         dens = dreal(dens_C)
         !
      end function get_dens
      !
      !
      !
   end subroutine set_density_Int


   !---------------------------------------------------------------------------!
   !PURPOSE: Returns the chemical potential shift in order to have the wanted
   !         LDA density
   !---------------------------------------------------------------------------!
   subroutine set_density_NonInt(mu,Beta,Lttc,mu_param)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use input_vars, only : Ntau
      implicit none
      !
      real(8),intent(inout)                 :: mu
      real(8),intent(in)                    :: Beta
      type(Lattice),intent(in)              :: Lttc
      type(musearch),intent(in)             :: mu_param
      !
      real(8)                               :: n_iter,mu_start,mu_last,mu_sign
      real(8)                               :: mu_below,mu_above,n_err
      complex(8),allocatable                :: Gitau(:,:,:)
      integer                               :: Norb,Nkpt
      integer                               :: iter
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- set_density_NonInt"
      !
      !
      ! Check on the input Fields
      if(.not.Lttc%status) stop "set_density_NonInt: Lttc not properly initialized."
      if(Lttc%Nkpt.eq.0) stop "set_density_NonInt: Lttc k dependent attributes not properly initialized."
      if(mu_param%TargetDensity.eq.0d0) stop "set_density_NonInt: TargetDensity is set to zero."
      write(*,"(A,F10.5)") "     Target density: ",mu_param%TargetDensity
      !
      Norb = Lttc%Norb
      Nkpt = Lttc%Nkpt
      mu_start = mu
      !
      !Everything is paramagnetic
      allocate(Gitau(Norb,Ntau,Nkpt));Gitau=czero
      call calc_G0_tau(Gitau,mu_start,Beta,Lttc%Ek,atBeta=.true.)
      n_iter = get_dens()
      write(*,"(2(A,F12.5))") "     Starting density: ",n_iter,", starting mu: ",mu_start
      !
      n_err = abs(n_iter-mu_param%TargetDensity)/mu_param%TargetDensity
      if(n_err.lt.mu_param%densityRelErr)then
         write(*,"(A)") "     Density already under threshold, skipping mu scan."
         return
      endif
      !
      mu_sign = sign(1d0,mu_param%TargetDensity-n_iter)
      mu_last = mu_start;
      !
      do iter=1,mu_param%muIter
         !
         mu = mu_start + iter * mu_param%muStep * mu_sign;
         !
         call calc_G0_tau(Gitau,mu,Beta,Lttc%Ek,atBeta=.true.)
         n_iter = get_dens()
         !
         write(*,"(A,I4,2(A,F12.5))") "     Rigid shift it# ",iter,", density: ", n_iter,", mu: ",mu
         !
         if((mu_sign.gt.0.0).and.(n_iter > mu_param%TargetDensity)) exit
         if((mu_sign.lt.0.0).and.(n_iter < mu_param%TargetDensity)) exit
         !
         mu_last = mu;
         !
      enddo
      !
      mu_below = min(mu,mu_last)
      mu_above = max(mu,mu_last)
      !
      do iter=1,mu_param%muIter
         !
         mu = (mu_below+mu_above)/2d0
         call calc_G0_tau(Gitau,mu,Beta,Lttc%Ek,atBeta=.true.)
         n_iter = get_dens()
         n_err = abs(n_iter-mu_param%TargetDensity)/mu_param%TargetDensity
         !
         write(*,"(A,I4,5(A,F12.5))") "     Netwon it# ",iter,", mu boundaries: ",mu_below," / ",mu_above, &
         ", density: ", n_iter,", mu: ",mu,", relative error: ",n_err
         !
         if(n_iter.gt.mu_param%TargetDensity) mu_above = mu;
         if(n_iter.lt.mu_param%TargetDensity) mu_below = mu;
         !
         if(n_err.lt.mu_param%densityRelErr)then
            write(*,"(A,F12.5)") "     Found correct chemical potential after "//str(iter)//" iterations: ",mu
            exit
         elseif(iter.eq.mu_param%muIter)then
            write(*,"(A,F12.5)") "     Warning: NOT found correct chemical potential after "//str(iter)//" iterations. Last used value: ",mu
         endif
         !
      enddo !iter
      deallocate(Gitau)
      !
   contains
      !
      !
      !
      function get_dens() result(dens)
         implicit none
         complex(8)               :: dens_C
         real(8)                  :: dens
         integer                  :: iwan,ik
         complex(8),allocatable   :: n_k(:,:,:)
         complex(8),allocatable   :: n_loc(:,:)
         !
         allocate(n_k(Norb,Norb,Nkpt));n_k=czero
         do ik=1,Nkpt
            !paramagnetic
            n_k(:,:,ik) = -2.d0*diag(dreal(Gitau(:,Ntau,ik)))
            !if present, the restriction is on the orbitals in Wannier basis
            n_k(:,:,ik) = rotate(n_k(:,:,ik),transpose(conjg(Lttc%Zk(:,:,ik))))
         enddo
         !
         allocate(n_loc(Norb,Norb));n_loc=czero
         n_loc = sum(n_k,dim=3)/Nkpt
         !
         if(allocated(mu_param%orbs))then
            do iwan=1,Norb
               if(all(mu_param%orbs.ne.iwan))cycle
               dens_C = dens_C + n_loc(iwan,iwan)
            enddo
         else
            dens_C = trace(n_loc)
         endif
         deallocate(n_loc,n_k)
         !
         if(aimag(dens_C).gt.eps)then
            write(*,"(A,2F10.5)")"     Density (NonInt) is complex: ",real(dens_C),aimag(dens_C)
            stop "get_dens(set_density_NonInt)"
         endif
         dens = dreal(dens_C)
         !
      end function get_dens
      !
      !
      !
   end subroutine set_density_NonInt


   !---------------------------------------------------------------------------!
   !PURPOSE: Prints lda Gf on different axis.
   !---------------------------------------------------------------------------!
   subroutine calc_Glda(mu,Beta,Lttc,pathOUTPUT)
      !
      use parameters
      use linalg
      use utils_misc
      use file_io
      use fourier_transforms
      use input_vars, only : pathINPUT
      use input_vars, only : Nreal, wrealMax, eta
      use input_vars, only : wmatsMax
      use input_vars, only : Ntau, tau_uniform
      implicit none
      !
      real(8),intent(in)                    :: mu
      real(8),intent(in)                    :: Beta
      type(Lattice),intent(in)              :: Lttc
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      character(len=256)                    :: pathOUTPUT_
      integer                               :: iwan1,iwan2,ik,iw,itau
      integer                               :: Norb,Nkpt,Nmats
      real(8),allocatable                   :: axis(:)
      complex(8),allocatable                :: zeta(:,:,:),invGf(:,:)
      complex(8),allocatable                :: Gk_print_H(:,:,:,:),Gk_print_E(:,:,:)
      complex(8),allocatable                :: Gk_itau_H(:,:,:,:),Gk_itau_E(:,:,:)
      !deafults
      logical                               :: print_G0real=.true.
      logical                               :: print_G0itau=.true.
      logical                               :: print_G0mats=.true.
      !change before compiling
      logical                               :: printAllK_G0real=.false.
      logical                               :: printAllK_G0iatu=.false.
      logical                               :: printAllK_G0mats=.false.
      logical                               :: fullFreq_G0mats=.false.
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Glda"
      pathOUTPUT_ = reg(pathINPUT)//"LDA/"
      if(present(pathOUTPUT)) pathOUTPUT_ = reg(pathOUTPUT)//"LDA/"
      !
      !
      ! Check on the input Fields
      if(.not.Lttc%status) stop "calc_Glda: Lttc not properly initialized."
      if(Lttc%Nkpt.eq.0) stop "calc_Glda: Lttc k dependent attributes not properly initialized."
      !
      Norb = Lttc%Norb
      Nkpt = Lttc%Nkpt
      !
      !
      !
      !print G(w) in diagonal and Wannier basis---------------------------------
      if(print_G0real)then
         !
         allocate(axis(Nreal));axis=0d0
         axis = linspace(-wrealMax*1.5d0,+wrealMax*1.5d0,Nreal)
         allocate(zeta(Norb,Norb,Nreal));zeta=czero
         do iwan1=1,Norb
            do iw=1,Nreal
               zeta(iwan1,iwan1,iw) = dcmplx(  axis(iw) + mu , eta )
            enddo
         enddo
         !
         allocate(Gk_print_E(Norb,Nreal,Nkpt));Gk_print_E=czero
         allocate(Gk_print_H(Norb,Norb,Nreal,Nkpt));Gk_print_H=czero
         allocate(invGf(Norb,Norb));invGf=czero
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(iwan1,iw,ik,invGf)
         !$OMP DO
         do ik=1,Nkpt
            do iw=1,Nreal
               !
               do iwan1=1,Norb
                  Gk_print_E(iwan1,iw,ik) = 1d0/( dcmplx(  axis(iw) + mu , eta ) - Lttc%Ek(iwan1,ik) )
               enddo
               !
               invGf = zeta(:,:,iw) - Lttc%Hk(:,:,ik)
               call inv(invGf)
               Gk_print_H(:,:,iw,ik) = invGf
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(zeta,invGf)
         !
         ! Print
         call print_G("Greal",Nreal,printAllK_G0real,.true.)
         deallocate(axis,Gk_print_H,Gk_print_E)
         !
      endif
      !
      !
      !
      !print G(tau) in diagonal and Wannier basis-------------------------------
      if(print_G0itau)then
         !
         allocate(axis(Ntau));axis=0d0
         if(tau_uniform)then
            axis = linspace(0d0,Beta,Ntau)
         else
            axis = denspace(Beta,Ntau)
         endif
         !
         allocate(Gk_print_E(Norb,Ntau,Nkpt));Gk_print_E=czero
         allocate(Gk_print_H(Norb,Norb,Ntau,Nkpt));Gk_print_H=czero
         call calc_G0_tau(Gk_print_E,mu,Beta,Lttc%Ek)
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(iwan1,iw,ik)
         !$OMP DO
         do ik=1,Nkpt
            do itau=1,Ntau
               Gk_print_H(:,:,itau,ik) = rotate(diag(Gk_print_E(:,itau,ik)),conjg(transpose(Lttc%Zk(:,:,ik))))
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         ! Print
         call print_G("Gitau",Ntau,printAllK_G0iatu,.false.)
         if(print_G0mats)then
            Gk_itau_H = Gk_print_H
            Gk_itau_E = Gk_print_E
         endif
         deallocate(axis,Gk_print_H,Gk_print_E)
         !
      endif
      !
      !
      !
      !print G(iw) in diagonal and Wannier basis--------------------------------
      if(print_G0mats)then
         !
         Nmats = int(Beta*wmatsMax/(2d0*pi))
         allocate(axis(Nmats));axis=0d0
         axis = FermionicFreqMesh(Beta,Nmats)
         !
         allocate(Gk_print_E(Norb,Nmats,Nkpt));Gk_print_E=czero
         allocate(Gk_print_H(Norb,Norb,Nmats,Nkpt));Gk_print_H=czero
         !
         call Fitau2mats_vec(Beta,Gk_itau_E,Gk_print_E,tau_uniform=tau_uniform)
         call Fitau2mats_mat(Beta,Gk_itau_H,Gk_print_H,tau_uniform=tau_uniform)
         deallocate(Gk_itau_E,Gk_itau_H)
         !
         ! Print
         call print_G("Gmats",Nmats,printAllK_G0mats,.false.)
         deallocate(axis,Gk_print_H,Gk_print_E)
         !
         if(fullFreq_G0mats)then
            !
            allocate(axis(2*Nmats+1));axis=0d0
            axis = FermionicFreqMesh(Beta,2*Nmats+1,full=.true.)
            allocate(zeta(Norb,Norb,2*Nmats+1));zeta=czero
            do iwan1=1,Norb
               do iw=1,2*Nmats+1
                  zeta(iwan1,iwan1,iw) = dcmplx(  mu , axis(iw) )
               enddo
            enddo
            !
            allocate(Gk_print_E(Norb,2*Nmats+1,Nkpt));Gk_print_E=czero
            allocate(Gk_print_H(Norb,Norb,2*Nmats+1,Nkpt));Gk_print_H=czero
            allocate(invGf(Norb,Norb));invGf=czero
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(iwan1,iw,ik,invGf)
            !$OMP DO
            do ik=1,Nkpt
               do iw=1,2*Nmats+1
                  !
                  do iwan1=1,Norb
                     Gk_print_E(iwan1,iw,ik) = 1d0/( dcmplx(  mu , axis(iw) ) - Lttc%Ek(iwan1,ik) )
                  enddo
                  !
                  invGf = zeta(:,:,iw) - Lttc%Hk(:,:,ik)
                  call inv(invGf)
                  Gk_print_H(:,:,iw,ik) = invGf
                  !
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(zeta,invGf)
            !
            ! Print
            call print_G("Gmats_fullw",(2*Nmats+1),printAllK_G0mats,.false.)
            deallocate(axis,Gk_print_H,Gk_print_E)
            !
         endif
         !
      endif
      !
      !
   contains
      !
      !
      !
      subroutine print_G(name,Naxis,allK,printAkw)
         implicit none
         character(len=*),intent(in)        :: name
         integer,intent(in)                 :: Naxis
         logical,intent(in)                 :: allK
         logical,intent(in)                 :: printAkw
         real(8),allocatable                :: Akw(:)
         complex(8),allocatable             :: G_print_H(:,:,:),G_print_E(:,:)
         !
         if(printAkw)allocate(Akw(Naxis));Akw=0d0
         !
         if(allK)then
            do ik=1,Nkpt
               do iwan1=1,Norb
                  call dump_FermionicField(Gk_print_E(iwan1,:,ik),reg(pathOUTPUT_),reg(name)//"_Ek_"//str(iwan1)//"_k"//str(ik)//".lda",axis)
                  if(printAkw)then
                     Akw = dimag(G_print_E(iwan1,:))
                     Akw = Akw/(sum(Akw)*abs(axis(2)-axis(1)))
                     call dump_FermionicField(Akw,reg(pathOUTPUT_),"Akw_Ek_"//str(iwan1)//"_k"//str(ik)//".lda",axis)
                  endif
                  do iwan2=1,Norb
                     call dump_FermionicField(Gk_print_H(iwan1,iwan2,:,ik),reg(pathOUTPUT_),reg(name)//"_Hk_"//str(iwan1)//str(iwan2)//"_k"//str(ik)//".lda",axis)
                     if(printAkw)then
                        Akw = dimag(G_print_H(iwan1,iwan2,:))
                        Akw = Akw/(sum(Akw)*abs(axis(2)-axis(1)))
                        call dump_FermionicField(Akw,reg(pathOUTPUT_),"Akw_Hk_"//str(iwan1)//str(iwan2)//"_k"//str(ik)//".lda",axis)
                     endif
                  enddo
               enddo
            enddo
         endif
         !
         allocate(G_print_E(Norb,Naxis))      ; G_print_E = sum(Gk_print_E,dim=3)/Nkpt
         allocate(G_print_H(Norb,Norb,Naxis)) ; G_print_H = sum(Gk_print_H,dim=4)/Nkpt
         do iwan1=1,Norb
            call dump_FermionicField(G_print_E(iwan1,:),reg(pathOUTPUT_),reg(name)//"_Ek_"//str(iwan1)//".lda",axis)
            if(printAkw)then
               Akw = dimag(G_print_E(iwan1,:))
               Akw = Akw/(sum(Akw)*abs(axis(2)-axis(1)))
               call dump_FermionicField(Akw,reg(pathOUTPUT_),"Akw_Ek_"//str(iwan1)//".lda",axis)
            endif
            do iwan2=1,Norb
               call dump_FermionicField(G_print_H(iwan1,iwan2,:),reg(pathOUTPUT_),reg(name)//"_Hk_"//str(iwan1)//str(iwan2)//".lda",axis)
               if(printAkw)then
                  Akw = dimag(G_print_H(iwan1,iwan2,:))
                  Akw = Akw/(sum(Akw)*abs(axis(2)-axis(1)))
                  call dump_FermionicField(Akw,reg(pathOUTPUT_),"Akw_Hk_"//str(iwan1)//str(iwan2)//".lda",axis)
               endif
            enddo
         enddo
         deallocate(G_print_E,G_print_H)
         if(printAkw)deallocate(Akw)
         !
      end subroutine print_G
      !
      !
   end subroutine calc_Glda


   !---------------------------------------------------------------------------!
   !PURPOSE: Prints lda Gf on different axis.
   !---------------------------------------------------------------------------!
   subroutine build_Potential(Potential,tz,Npot,zeta,Hk,Smats)
      !
      use parameters
      use linalg, only : inv, rotate
      use utils_misc
      use input_vars, only : paramagnet
      implicit none
      !
      complex(8),intent(inout)              :: Potential(:,:,:,:,:)
      complex(8),intent(in)                 :: tz(:,:)
      integer,intent(in)                    :: Npot
      complex(8),intent(in)                 :: zeta(:,:,:)
      complex(8),intent(in)                 :: Hk(:,:,:)
      complex(8),intent(in)                 :: Smats(:,:,:,:,:)
      !
      complex(8),allocatable                :: invGbulk(:,:)
      complex(8),allocatable                :: Gbulk(:,:)
      complex(8),allocatable                :: Ptmp(:,:)
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iw,ik,ispin,ibulk
      !
      !
      if(verbose)write(*,"(A)") "---- build_Potential"
      !
      !
      Norb = size(Potential,dim=1)
      Nmats = size(zeta,dim=3)
      Nkpt = size(Hk,dim=3)
      !
      call assert_shape(Potential,[Norb,Norb,Nmats,Nkpt,Nspin],"build_Potential","Potential")
      call assert_shape(tz,[Norb,Norb],"build_Potential","tz")
      call assert_shape(zeta,[Norb,Norb,Nmats],"build_Potential","zeta")
      call assert_shape(Hk,[Norb,Norb,Nkpt],"build_Potential","Hk")
      call assert_shape(Smats,[Norb,Norb,Nmats,Nkpt,Nspin],"build_Smats","Potential")
      !
      !G and invG of each layer are constant
      allocate(invGbulk(Norb,Norb));invGbulk=czero
      allocate(Gbulk(Norb,Norb));Gbulk=czero
      allocate(Ptmp(Norb,Norb));Ptmp=czero
      spinPloop: do ispin=1,Nspin
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(ispin,Nmats,Nkpt,zeta,Hk,Smats,Npot,tz,Potential),&
         !$OMP PRIVATE(ik,iw,ibulk,invGbulk,Gbulk,Ptmp)
         !$OMP DO
         do iw=1,Nmats
            do ik=1,Nkpt
               !
               invGbulk = zeta(:,:,iw) - Hk(:,:,ik) - Smats(:,:,iw,ik,ispin)
               Gbulk = invGbulk
               call inv(Gbulk)
               !
               !first t*G*t
               Potential(:,:,iw,ik,ispin) = rotate(Gbulk,tz)
               do ibulk=2,Npot
                  Ptmp = invGbulk - Potential(:,:,iw,ik,ispin)
                  call inv(Ptmp)
                  Potential(:,:,iw,ik,ispin) = rotate(Ptmp,tz)
               enddo
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         if(paramagnet)then
            Potential(:,:,:,:,Nspin) = Potential(:,:,:,:,1)
            exit spinPloop
         endif
      enddo spinPloop
      deallocate(Gbulk,invGbulk,Ptmp)
      !
   end subroutine build_Potential


end module greens_function
