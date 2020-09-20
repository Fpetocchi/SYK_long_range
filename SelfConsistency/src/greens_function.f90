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
      module procedure set_density_Int                                          !(Gmats,Lattice,n_target,Smats,orbs(optional))
      module procedure set_density_NonInt                                       !(mu,Beta,Lattice,n_target,orbs(optional))
   end interface set_density


   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: calc_density
   public :: set_density
   public :: calc_Gmats

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
      use input_vars, only : Ntau, tau_uniform
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
      write(*,"(A)") "--- calc_density_loc ---"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(Gmats%Npoints.eq.0) stop "Gmats frequency dependent attributes not properly initialized."
      !
      Norb = Gmats%Norb
      Beta = Gmats%Beta
      !
      call assert_shape(n_loc,[Norb,Norb,Nspin],"calc_density_loc","n_loc")
      !
      n_loc=czero
      allocate(Gitau(Norb,Norb,Ntau));Gitau=czero
      do ispin=1,Nspin
         !
         call Fmats2itau_mat(Beta,Gmats%ws(:,:,:,ispin),Gitau, &
         asympt_corr=.true.,tau_uniform=tau_uniform,atBeta=.true.)
         !
         n_loc(:,:,ispin) = -Gitau(:,:,Ntau)
         !
      enddo
      deallocate(Gitau)
      !
   end subroutine calc_density_loc
   !
   subroutine calc_density_Kdep(Gmats,Lttc,n_k)
      !
      use parameters
      use utils_misc
      use fourier_transforms
      use input_vars, only : Ntau, tau_uniform
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      complex(8),allocatable,intent(inout)  :: n_k(:,:,:,:)
      !
      complex(8),allocatable                :: Gitau(:,:,:,:)
      real(8)                               :: Beta
      integer                               :: ispin,Norb,Nkpt
      !
      !
      write(*,"(A)") "--- calc_density_Kdep ---"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(Gmats%Npoints.eq.0) stop "Gmats frequency dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "Gmats k dependent attributes not properly initialized."
      if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "Lttc has different number of k-points with respect to Gmats."
      !
      Norb = Gmats%Norb
      Nkpt = Gmats%Nkpt
      Beta = Gmats%Beta
      !
      call assert_shape(n_k,[Norb,Norb,Nkpt,Nspin],"calc_density_Kdep","n_k")
      !
      n_k=czero
      allocate(Gitau(Norb,Norb,Ntau,Nkpt));Gitau=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau, &
         asympt_corr=.true.,tau_uniform=tau_uniform,Nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt,atBeta=.true.)
         !
         n_k(:,:,:,ispin) = -Gitau(:,:,:,Ntau)
         !
      enddo
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
      use input_vars, only : Ntau, tau_uniform
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
      write(*,"(A)") "--- calc_G0_tau ---"
      fermicut=log(huge(1.0d0)-1e2)/2.d0
      !
      !
      ! Check on the input Fields
      Norb = size(Ek,dim=1)
      Nkpt = size(Ek,dim=2)
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
   subroutine calc_Gmats(Gmats,Lttc,Smats)
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
      !
      complex(8),allocatable                :: invGf(:,:)
      complex(8),pointer                    :: Swks(:,:,:,:,:)
      complex(8),allocatable                :: zeta(:,:,:)
      complex(8),allocatable                :: n_k(:,:,:,:)
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: Beta,mu
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iw,ik,iwan,ispin
      !
      !
      write(*,"(A)") "--- calc_Gmats ---"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(.not.Lttc%status) stop "Lttc not properly initialized."
      if(Gmats%Npoints.eq.0) stop "Gmats frequency dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "Gmats k dependent attributes not properly initialized."
      if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "Lttc has different number of k-points with respect to Gmats."
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
         if(.not.Smats%status) stop "Smats not properly initialized."
         if(Smats%Npoints.eq.0) stop "Smats frequency dependent attributes not properly initialized."
         if(Smats%Nkpt.eq.0) stop "Smats k dependent attributes not properly initialized."
         Swks => Smats%wks
      endif
      !
      call clear_attributes(Gmats)
      allocate(invGf(Norb,Norb));invGf=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(zeta,Lttc,Gmats,Swks),&
      !$OMP PRIVATE(ispin,ik,iw,invGf)
      !$OMP DO
      do ispin=1,Nspin
         do ik=1,Gmats%Nkpt
            do iw=1,Gmats%Npoints
               !
               invGf = zeta(:,:,iw) - Lttc%Hk(:,:,ik)
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
      call calc_density(Gmats,Lttc,n_k)
      Gmats%N_ks = n_k
      deallocate(n_k)
      !
      ! In the N_s the local density
      call FermionicKsum(Gmats)
      !
   end subroutine calc_Gmats



   !---------------------------------------------------------------------------!
   !PURPOSE: Set the density of a Green's function and returns the adjusted mu
   !---------------------------------------------------------------------------!
   subroutine set_density_Int(Gmats,Lttc,n_target,Smats,orbs)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use input_vars, only : densityPercErr
      implicit none
      !
      type(FermionicField),intent(inout)    :: Gmats
      type(Lattice),intent(in)              :: Lttc
      type(FermionicField),intent(in)       :: Smats
      real(8),intent(in)                    :: n_target
      real(8),intent(in),optional           :: orbs(:)
      !
      real(8)                               :: n_iter,mu_start
      real(8)                               :: emin,emax,dmu
      integer                               :: Norb,Nmats,Nkpt
      integer                               :: iter
      !
      !
      write(*,"(A)") "--- set_density_Int ---"
      write(*,"(A1,F10.5)") "Target density: ",n_target
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(.not.Lttc%status) stop "Lttc not properly initialized."
      if(.not.Smats%status) stop "Smats not properly initialized."
      if(Gmats%Npoints.eq.0) stop "Gmats frequency dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "Gmats k dependent attributes not properly initialized."
      if(Smats%Npoints.eq.0) stop "Smats frequency dependent attributes not properly initialized."
      if(Smats%Nkpt.eq.0) stop "Smats k dependent attributes not properly initialized."
      if(Gmats%Nkpt.ne.Lttc%Nkpt) stop "Lttc has different number of k-points with respect to Gmats."
      !
      Norb = Gmats%Norb
      Nmats = Gmats%Npoints
      Nkpt = Gmats%Nkpt
      mu_start = Gmats%mu
      !
      call calc_Gmats(Gmats,Lttc,Smats)
      n_iter = get_dens()
      write(*,"(A,F10.5)") "Starting density: ",n_iter
      write(*,"(A,F10.5)") "Starting mu: ",mu_start
      !
      dmu=0.d0
      emin=-30d0
      emax=+30d0
      do iter=1,100
         !
         write(*,"(A,I5)") "iteration #",iter
         !
         ! Chemical potential variation
         if (n_iter.gt.n_target) then
            emax=dmu
            dmu=(emax+emin)/2
         else
            emin=dmu
            dmu=(emax+emin)/2
         endif
         !
         Gmats%mu = mu_start + dmu
         call calc_Gmats(Gmats,Lttc,Smats)
         n_iter = get_dens()
         !
         write(*,"(2(A,F10.5))") "density: ", n_iter," Dmu: ", dmu
         !
         ! Exit condition 1: succesful search
         if (dabs((n_iter-n_target)/n_target).lt.densityPercErr) then
            write(*,"(A,I5,A)") "Density converged after ",iter," iterations."
            write(*,"(A,F10.5)") "Absolute density error: ", dabs(n_iter-n_target)
            write(*,"(A,2F10.5)") "New chemical potential: ",Gmats%mu
            write(*,"(A,F10.5)") "Old chemical potential: ", mu_start
            if(Gmats%mu.ne.(mu_start + dmu)) stop "Problem in the chemical potential."
            exit
         endif
         !
         ! Exit condition 2: un-succesful search
         if (iter.eq.100)then
            write(*,"(A)") "Chemical potential search not converged after 100 iterations."
            exit
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
         if(present(orbs))then
            do iwan=1,Norb
               if(all(orbs.ne.iwan))cycle
               do ispin=1,Nspin
                  dens_C = dens_C + Gmats%N_s(iwan,iwan,ispin)
               enddo
            enddo
         else
            dens_C = trace(sum(Gmats%N_s,dim=3))
         endif
         !
         if(aimag(dens_C).gt.eps)then
            write(*,"(A,2F10.5)")"Density is complex: ",real(dens_C),aimag(dens_C)
            stop
         endif
         dens = real(dens_C)
         !
      end function get_dens
      !
      !
      !
   end subroutine set_density_Int
   !
   subroutine set_density_NonInt(mu,Beta,Lttc,n_target,orbs)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use input_vars, only : Ntau, densityPercErr
      implicit none
      !
      real(8),intent(out)                   :: mu
      real(8),intent(in)                    :: Beta
      type(Lattice),intent(in)              :: Lttc
      real(8),intent(in)                    :: n_target
      real(8),intent(in),optional           :: orbs(:)
      !
      real(8)                               :: n_iter,mu_start
      complex(8),allocatable                :: Gitau(:,:,:)
      complex(8),allocatable                :: n_loc(:,:)
      complex(8),allocatable                :: n_k(:,:,:)
      real(8)                               :: emin,emax,dmu
      integer                               :: Norb,Nkpt
      integer                               :: iter
      !
      !
      write(*,"(A)") "--- set_density_NonInt ---"
      write(*,"(A1,F10.5)") "Target density: ",n_target
      !
      !
      ! Check on the input Fields
      if(.not.Lttc%status) stop "Lttc not properly initialized."
      if(Lttc%Nkpt.eq.0) stop "Lttc k dependent attributes not properly initialized."
      !
      Norb = Lttc%Norb
      Nkpt = Lttc%Nkpt
      mu_start = mu
      !
      !Everything is paramagnetic
      allocate(n_loc(Norb,Norb));n_loc=czero
      allocate(n_k(Norb,Norb,Nkpt));n_k=czero
      allocate(Gitau(Norb,Ntau,Nkpt));Gitau=czero
      call calc_G0_tau(Gitau,mu_start,Beta,Lttc%Ek,atBeta=.true.)
      n_iter = get_dens()
      write(*,"(A,F10.5)") "Starting density: ",n_iter
      write(*,"(A,F10.5)") "Starting mu: ",mu_start
      !
      dmu=0.d0
      emin=-30d0
      emax=+30d0
      do iter=1,100
         !
         write(*,"(A,I5)") "iteration #",iter
         !
         ! Chemical potential variation
         if (n_iter.gt.n_target) then
            emax=dmu
            dmu=(emax+emin)/2
         else
            emin=dmu
            dmu=(emax+emin)/2
         endif
         !
         mu = mu_start + dmu
         call calc_G0_tau(Gitau,mu,Beta,Lttc%Ek,atBeta=.true.)
         n_iter = get_dens()
         !
         write(*,"(2(A,F10.5))") "density: ", n_iter," Dmu: ", dmu
         !
         ! Exit condition 1: succesful search
         if (dabs((n_iter-n_target)/n_target).lt.densityPercErr) then
            write(*,"(A,I5,A)") "Density converged after ",iter," iterations."
            write(*,"(A,F10.5)") "Absolute density error: ", dabs(n_iter-n_target)
            write(*,"(A,2F10.5)") "New chemical potential: ",mu
            write(*,"(A,F10.5)") "Old chemical potential: ", mu_start
            if(mu.ne.(mu_start + dmu)) stop "Problem in the chemical potential."
            exit
         endif
         !
         ! Exit condition 2: un-succesful search
         if (iter.eq.100)then
            write(*,"(A)") "Chemical potential search not converged after 100 iterations."
            exit
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
            n_k(:,:,ik) = 2.d0*diag(Gitau(:,Ntau,ik))
            !if present, the restriction is on the orbitals in Wannier basis
            n_k(:,:,ik) = rotate(n_k(:,:,ik),transpose(conjg(Lttc%Zk(:,:,ik))))
         enddo
         n_loc = sum(n_k,dim=3)/Nkpt
         !
         if(present(orbs))then
            do iwan=1,Norb
               if(all(orbs.ne.iwan))cycle
               dens_C = dens_C + n_loc(iwan,iwan)
            enddo
         else
            dens_C = trace(n_loc)
         endif
         !
         if(aimag(dens_C).gt.eps)then
            write(*,"(A,2F10.5)")"Density is complex: ",real(dens_C),aimag(dens_C)
            stop
         endif
         dens = real(dens_C)
         !
      end function get_dens
      !
      !
      !
   end subroutine set_density_NonInt



end module greens_function
