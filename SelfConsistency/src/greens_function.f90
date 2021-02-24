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

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the density given the Green's function
   !TEST ON: 16-10-2020(Kdep)
   !---------------------------------------------------------------------------!
   subroutine calc_density_loc(Gmats,n_loc)
      !
      use parameters
      use utils_misc
      use fourier_transforms
      use input_vars, only : NtauF, tau_uniform, paramagnet
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
      allocate(Gitau(Norb,Norb,NtauF));Gitau=czero
      spinloop: do ispin=1,Nspin
         !
         call Fmats2itau_mat(Beta,Gmats%ws(:,:,:,ispin),Gitau, &
         asympt_corr=.true.,tau_uniform=tau_uniform,atBeta=.true.)
         !
         n_loc(:,:,ispin) = -Gitau(:,:,NtauF)
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
   subroutine calc_density_Kdep(Gmats,Lttc,n_k,along_path)
      !
      use parameters
      use utils_misc
      use fourier_transforms
      use input_vars, only : NtauF, tau_uniform, cmplxWann, paramagnet
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      complex(8),allocatable,intent(inout)  :: n_k(:,:,:,:)
      logical,intent(in),optional           :: along_path
      !
      complex(8),allocatable                :: Gitau(:,:,:,:)
      real(8)                               :: Beta
      integer                               :: ispin,Norb,Nkpt
      logical                               :: along_path_
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
      !
      if(along_path_)then
         if(Gmats%Nkpt.ne.Lttc%Nkpt_path) stop "calc_density_Kdep: Lttc has different number of path k-points with respect to Gmats."
         if(.not.allocated(Lttc%kptpath)) stop "calc_density_Kdep: K-point path not allocated."
      else
         if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "calc_density_Kdep: Lttc has different number of k-points with respect to Gmats."
      endif
      !
      Norb = Gmats%Norb
      Nkpt = Gmats%Nkpt
      Beta = Gmats%Beta
      !
      call assert_shape(n_k,[Norb,Norb,Nkpt,Nspin],"calc_density_Kdep","n_k")
      !
      n_k=czero
      allocate(Gitau(Norb,Norb,NtauF,Nkpt));Gitau=czero
      spinloop: do ispin=1,Nspin
         !
         if(cmplxWann.or.along_path_)then
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau, &
            asympt_corr=.true.,tau_uniform=tau_uniform,atBeta=.true.)
         else
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau, &
            asympt_corr=.true.,tau_uniform=tau_uniform,Nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt,atBeta=.true.)
         endif
         !
         n_k(:,:,:,ispin) = -Gitau(:,:,NtauF,:)
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
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_G0_tau(Gitau,mu,Beta,Ek,atBeta)
      !
      use parameters
      use utils_misc
      use fourier_transforms
      use input_vars, only : NtauF, tau_uniform
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
      integer                               :: iwan,ik,itau
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
      call assert_shape(Gitau,[Norb,NtauF,Nkpt],"calc_G0_tau","Gitau")
      atBeta_ = .false.
      if(present(atBeta)) atBeta_ = atBeta
      !
      allocate(tau(NtauF));tau=0d0
      if(tau_uniform)then
         tau = linspace(0d0,Beta,NtauF)
      else
         tau = denspace(Beta,NtauF)
      endif
      !
      Gitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,NtauF,Norb,atBeta_,tau,Ek,mu,fermicut,Beta,Gitau),&      !VEDI SE PUOI TOGLIERE QUALCHE IF
      !$OMP PRIVATE(ik,itau,iwan,eu,upper,lower)
      !$OMP DO
      do itau=1,NtauF
         if(atBeta_.and.(itau.ne.NtauF))cycle
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
   !TEST ON: 16-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_Gmats_Full(Gmats,Lttc,Smats,along_path)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gmats
      type(Lattice),intent(in)              :: Lttc
      type(FermionicField),intent(in),optional,target :: Smats
      logical,intent(in),optional           :: along_path
      !
      complex(8),allocatable                :: invGf(:,:)
      complex(8),pointer                    :: Swks(:,:,:,:,:)
      complex(8),allocatable                :: zeta(:,:,:)
      complex(8),allocatable                :: n_k(:,:,:,:)
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: Beta,mu
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iw,ik,iwan,ispin
      logical                               :: along_path_
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
      !
      if(along_path_)then
         if(Gmats%Nkpt.ne.Lttc%Nkpt_path) stop "calc_Gmats_Full: Lttc has different number of path k-points with respect to Gmats."
         if(.not.allocated(Lttc%Hk_path)) stop "calc_Gmats_Full: H(k) along path not allocated."
      else
         if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "calc_Gmats_Full: Lttc has different number of k-points with respect to Gmats."
      endif
      !
      Norb = Gmats%Norb
      Nmats = Gmats%Npoints
      Nkpt = Gmats%Nkpt
      Beta = Gmats%Beta
      mu = Gmats%mu
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
         Swks => Smats%wks
         if(verbose)write(*,"(A)") "     Interacting Green's function."
      else
         if(verbose)write(*,"(A)") "     LDA Green's function."
      endif
      !
      call clear_attributes(Gmats)
      allocate(invGf(Norb,Norb));invGf=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(zeta,Lttc,Gmats,Swks,along_path_),&
      !$OMP PRIVATE(ispin,ik,iw,invGf)
      !$OMP DO
      do ispin=1,Nspin
         do ik=1,Gmats%Nkpt
            do iw=1,Gmats%Npoints
               !
               if(along_path_)then
                  invGf = zeta(:,:,iw) - Lttc%Hk_path(:,:,ik)
               else
                  invGf = zeta(:,:,iw) - Lttc%Hk(:,:,ik)
               endif
               !
               if(associated(Swks)) invGf = invGf - Swks(:,:,iw,ik,ispin)
               !
               call inv(invGf)
               Gmats%wks(:,:,iw,ik,ispin) = invGf
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(zeta,invGf)
      if(associated(Swks))nullify(Swks)
      !
      ! In the N_ks attribute is stored the k-dep occupation
      allocate(n_k(Norb,Norb,Nkpt,Nspin));n_k=czero
      call calc_density(Gmats,Lttc,n_k,along_path=along_path_)
      Gmats%N_ks = n_k
      deallocate(n_k)
      !
      ! In the N_s the local density
      if(.not.along_path_)call FermionicKsum(Gmats)
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
      allocate(zeta(Norb,Norb));zeta=0d0
      zeta = deye(Norb)*mu_shift
      !
      allocate(invGf(Norb,Norb));invGf=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(zeta,Gmats),&
      !$OMP PRIVATE(ispin,ik,iw,invGf)
      !$OMP DO
      do ispin=1,Nspin
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
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(zeta,invGf)
      !
      ! In the N_ks attribute is stored the k-dep occupation
      allocate(n_k(Norb,Norb,Nkpt,Nspin));n_k=czero
      call calc_density(Gmats,Lttc,n_k)
      Gmats%N_ks = n_k
      deallocate(n_k)
      !
      ! In the N_s the local density
      call FermionicKsum(Gmats)
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
            write(*,"(A,2F10.5)")"Density (Int) is complex: ",real(dens_C),aimag(dens_C)
            stop "get_dens(set_density_Int)"
         endif
         dens = real(dens_C)
         !
      end function get_dens
      !
      !
      !
   end subroutine set_density_Int


   !---------------------------------------------------------------------------!
   !PURPOSE: Returns the chemical potential shift in order to have the wanted
   !         LDA density
   !TEST ON: 06-11-2020
   !---------------------------------------------------------------------------!
   subroutine set_density_NonInt(mu,Beta,Lttc,mu_param)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use input_vars, only : NtauF
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
      complex(8),allocatable                :: n_loc(:,:)
      complex(8),allocatable                :: n_k(:,:,:)
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
      allocate(n_loc(Norb,Norb));n_loc=czero
      allocate(n_k(Norb,Norb,Nkpt));n_k=czero
      allocate(Gitau(Norb,NtauF,Nkpt));Gitau=czero
      call calc_G0_tau(Gitau,mu_start,Beta,Lttc%Ek,atBeta=.true.)
      n_iter = get_dens()
      write(*,"(2(A,F12.5))") "     Starting density: ",n_iter,", starting mu: ",mu_start
      !
      mu_sign = sign(1d0,mu_param%TargetDensity-n_iter)
      mu_last = mu_start;
      !
      do iter=1,mu_param%muIter
         !
         mu = mu_start + iter * 0.05 * mu_sign;
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
      deallocate(Gitau,n_loc,n_k)
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
         do ik=1,Nkpt
            !paramagnetic
            n_k(:,:,ik) = -2.d0*diag(Gitau(:,NtauF,ik))
            !if present, the restriction is on the orbitals in Wannier basis
            n_k(:,:,ik) = rotate(n_k(:,:,ik),transpose(conjg(Lttc%Zk(:,:,ik))))
         enddo
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
         !
         if(aimag(dens_C).gt.eps)then
            write(*,"(A,2F10.5)")"Density (NonInt) is complex: ",real(dens_C),aimag(dens_C)
            stop "get_dens(set_density_NonInt)"
         endif
         dens = real(dens_C)
         !
      end function get_dens
      !
      !
      !
   end subroutine set_density_NonInt


   !---------------------------------------------------------------------------!
   !PURPOSE: Prints lda Gf on different axis.
   !TEST ON: 14-10-2020
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
      use input_vars, only : NtauF, tau_uniform
      implicit none
      !
      real(8),intent(in)                    :: mu
      real(8),intent(in)                    :: Beta
      type(Lattice),intent(in)              :: Lttc
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      character(len=256)                    :: pathOUTPUT_
      real(8),allocatable                   :: axis(:),Akw(:)
      complex(8),allocatable                :: zeta(:,:,:),invGf(:,:)
      complex(8),allocatable                :: Gitau(:,:,:),Gftmp(:,:)
      complex(8),allocatable                :: Gprint_Hk(:,:,:)
      complex(8),allocatable                :: Gprint_Ek(:,:)
      integer                               :: iwan1,iwan2,ik,iw,itau
      integer                               :: Norb,Nkpt,Nmats
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
      !print G(w) in diagonal and Wannier basis---------------------------------
      allocate(axis(Nreal));axis=0d0
      axis = linspace(-wrealMax,+wrealMax,Nreal)
      allocate(zeta(Norb,Norb,Nreal));zeta=czero
      do iwan1=1,Norb
         do iw=1,Nreal
            zeta(iwan1,iwan1,iw) = dcmplx(  axis(iw) + mu , eta )
         enddo
      enddo
      !
      allocate(Gprint_Ek(Norb,Nreal));Gprint_Ek=czero
      allocate(Gprint_Hk(Norb,Norb,Nreal));Gprint_Hk=czero
      !
      allocate(invGf(Norb,Norb));invGf=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,Nreal,Nkpt,mu,zeta,eta,axis,Lttc,Gprint_Ek,Gprint_Hk),&
      !$OMP PRIVATE(iwan1,ik,iw,invGf)
      !$OMP DO
      do ik=1,Nkpt
         do iw=1,Nreal
            !
            do iwan1=1,Norb
               Gprint_Ek(iwan1,iw) = Gprint_Ek(iwan1,iw) + 1d0/( dcmplx(  axis(iw) + mu , eta ) - Lttc%Ek(iwan1,ik) )/Nkpt
            enddo
            !
            invGf = zeta(:,:,iw) - Lttc%Hk(:,:,ik)
            !
            call inv(invGf)
            Gprint_Hk(:,:,iw) = Gprint_Hk(:,:,iw) + invGf/Nkpt
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(zeta,invGf)
      !
      ! Print
      allocate(Akw(Nreal));Akw=0d0
      do iwan1=1,Norb
         call dump_FermionicField(Gprint_Ek(iwan1,:),reg(pathOUTPUT_),"Greal_Ek_"//str(iwan1)//".lda",axis)
         Akw = dimag(Gprint_Ek(iwan1,:));Akw = Akw/sum(Akw)
         call dump_FermionicField(Akw,reg(pathOUTPUT_),"Akw_Ek_"//str(iwan1)//".lda",axis)
         do iwan2=1,Norb
            call dump_FermionicField(Gprint_Hk(iwan1,iwan2,:),reg(pathOUTPUT_),"Greal_Hk_"//str(iwan1)//"_"//str(iwan2)//".lda",axis)
            Akw = dimag(Gprint_Hk(iwan1,iwan2,:));Akw = Akw/sum(Akw)
            call dump_FermionicField(Akw,reg(pathOUTPUT_),"Akw_Hk_"//str(iwan1)//"_"//str(iwan2)//".lda",axis)
         enddo
      enddo
      deallocate(axis,Gprint_Hk,Gprint_Ek,Akw)
      !
      !
      !print G(tau) in diagonal and Wannier basis-------------------------------
      allocate(axis(NtauF));axis=0d0
      if(tau_uniform)then
         axis = linspace(0d0,Beta,NtauF)
      else
         axis = denspace(Beta,NtauF)
      endif
      !
      allocate(Gitau(Norb,NtauF,Nkpt));Gitau=czero
      call calc_G0_tau(Gitau,mu,Beta,Lttc%Ek)
      !
      allocate(Gprint_Ek(Norb,NtauF));Gprint_Ek=czero
      allocate(Gprint_Hk(Norb,Norb,NtauF));Gprint_Hk=czero
      allocate(Gftmp(Norb,Norb));Gftmp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,NtauF,Nkpt,Lttc,Gitau,Gprint_Ek,Gprint_Hk),&
      !$OMP PRIVATE(iwan1,ik,iw,Gftmp)
      !$OMP DO
      do itau=1,NtauF
         do ik=1,Nkpt
            !
            do iwan1=1,Norb
               Gprint_Ek(iwan1,itau) = Gprint_Ek(iwan1,itau) + Gitau(iwan1,itau,ik)/Nkpt
            enddo
            !
            Gftmp = diag( Gitau(:,itau,ik) )/Nkpt
            !
            Gprint_Hk(:,:,itau) = Gprint_Hk(:,:,itau) + rotate(Gftmp,conjg(transpose(Lttc%Zk(:,:,ik))))
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Gitau,Gftmp)
      !
      ! Print
      do iwan1=1,Norb
         call dump_FermionicField(Gprint_Ek(iwan1,:),reg(pathOUTPUT_),"Gitau_Ek_"//str(iwan1)//".lda",axis)
         do iwan2=1,Norb
            call dump_FermionicField(Gprint_Hk(iwan1,iwan2,:),reg(pathOUTPUT_),"Gitau_Hk_"//str(iwan1)//"_"//str(iwan2)//".lda",axis)
         enddo
      enddo
      deallocate(axis)!,Gprint_Hk,Gprint_Ek)
      !
      !
      !print G(iw) in diagonal and Wannier basis--------------------------------
      Nmats = int(Beta*wmatsMax/(2d0*pi))
      allocate(axis(Nmats));axis=0d0
      axis = FermionicFreqMesh(Beta,Nmats)
      allocate(zeta(Norb,Norb,Nmats)) !Temporaty for Gmats
      !
      zeta=czero
      call Fitau2mats_mat(Beta,Gprint_Hk,zeta,tau_uniform=tau_uniform)
      do iwan1=1,Norb
         do iwan2=1,Norb
            call dump_FermionicField(zeta(iwan1,iwan2,:),reg(pathOUTPUT_),"Gmats_Hk_"//str(iwan1)//"_"//str(iwan2)//".lda",axis)
         enddo
      enddo
      !
      zeta=czero
      call Fitau2mats_vec(Beta,Gprint_Ek,zeta(1,:,:),tau_uniform=tau_uniform)
      do iwan1=1,Norb
         call dump_FermionicField(zeta(1,iwan1,:),reg(pathOUTPUT_),"Gmats_Ek_"//str(iwan1)//".lda",axis)
      enddo
      deallocate(axis,Gprint_Hk,Gprint_Ek,zeta)
      !
   end subroutine calc_Glda



end module greens_function
