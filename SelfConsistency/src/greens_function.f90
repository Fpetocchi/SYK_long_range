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
      module procedure calc_densityGf_loc                                       !(FermionicField,n_loc[Norb,Norb,Nspin])
      module procedure calc_densityGf_Kdep                                      !(FermionicField,Lttc,n_k[Norb,Norb,Nkpt,Nspin])
      module procedure calc_densityHk_LOC            TO DO-USING ANALYTICAL IN BETA INTERNALLY-NO CALL TO GMATS                      !(Lttc,n_k[Norb,Norb,Nkpt,Nspin])
      module procedure calc_densityHk_Kdep           TO DO-USING ANALYTICAL IN BETA INTERNALLY-NO CALL TO GMATS                      !(Lttc,n_k[Norb,Norb,Nkpt,Nspin])
   end interface calc_density

   interface set_density
      module procedure set_density_Int                                          !(Gmats,Lttc,n_target,Smats,orbs)
      module procedure set_density_NonInt      TO DO-USING calc_densityHk_LOC   !(mu,Lttc,n_target,Smats,orbs)
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
   !public :: calc_chem_pot
   public :: set_density
   public :: calc_Gmats                                                         !

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the density given the Green's function
   !---------------------------------------------------------------------------!
   subroutine calc_densityGf_loc(Gmats,n_loc)
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
      write(*,"(A)") "--- calc_densityGf_loc ---"
      !
      !
      ! Check on the input Fields
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(Gmats%Npoints.eq.0) stop "Gmats frequency dependent attributes not properly initialized."
      !
      Norb = Gmats%Norb
      Beta = Gmats%Beta
      !
      call assert_shape(n_loc,[Norb,Norb,Nspin],"calc_density","n_loc")
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
   end subroutine calc_densityGf_loc
   !
   subroutine calc_densityGf_Kdep(Gmats,Lttc,n_k)
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
      write(*,"(A)") "--- calc_densityGf_Kdep ---"
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
      call assert_shape(n_k,[Norb,Norb,Nkpt,Nspin],"calc_density","n_k")
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
   end subroutine calc_densityGf_Kdep










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
      type(FermionicField),intent(inout)       :: Gmats
      type(Lattice),intent(in)                 :: Lttc
      type(FermionicField),intent(in),optional,target :: Smats
      !
      complex(8),allocatable                   :: invGf(:,:)
      complex(8),pointer                       :: Swks(:,:,:,:,:)
      complex(8),allocatable                   :: zeta(:,:,:)
      complex(8),allocatable                   :: n_k(:,:,:,:)
      real(8)                                  :: Beta,mu
      integer                                  :: Norb,Nmats,Nkpt
      integer                                  :: ik,iwan1,iwan2,ispin
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
      do iwan1=1,Norb
         do iw=1,Nmats
            zeta(iwan1,iwan1,iw) = dcmplx( mu , wmats(iw) )
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
      allocate(n_k(Norb,Norb,Ntau,Nkpt));n_k=czero
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
      type(FermionicField),intent(inout)       :: Gmats
      type(Lattice),intent(in)                 :: Lttc
      type(FermionicField),intent(in)          :: Smats
      real(8),intent(in)                       :: n_target
      real(8),intent(in),optional              :: orbs
      !
      real(8)                                  :: mu
      real(8)                                  :: emin,emax,dmu
      integer                                  :: Norb,Nmats,Nkpt
      integer                                  :: ik,iwan1,iwan2,ispin
      integer                                  :: ik,iwan1,iwan2,ispin
      integer                                  :: iter
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
      call calc_density(Gmats,n_loc)
      n_iter = get_dens(n_loc)
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
         n_iter = get_dens(Gmats%N_s)
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
      function get_dens() result(dens)
         implicit none
         real(8)      :: dens
         complex(8)   :: dens_C
         dens_C = czero
         if(present(orbs))then
            do iwan1=1,Norb
               if(all(orbs.ne.i))cycle
               do ispin=1,Nspin
                  dens_C = dens_C + n_loc(iwan1,iwan1,ispin)
               enddo
            enddo
         else
            dens_C = trace(sum(n_loc,dim=3))
         endif
         if(aimag(dens_C).gt.eps) stop "Density is complex."
      end function linspace
      !

   end subroutine set_density_Int
   !
   subroutine set_density_NonInt(mu_out,Lttc,n_target,orbs) FINISCI: QUESTO DEVE CHIAMARE calc_density(Lttc,n_loc) in modo da interfacciarsi con calc_densityHk_LOC
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use input_vars, only : densityPercErr
      implicit none
      !
      real(8),intent(in)                       :: mu_out
      type(Lattice),intent(in)                 :: Lttc
      type(FermionicField),intent(in)          :: Smats
      real(8),intent(in)                       :: n_target
      real(8),intent(in),optional              :: orbs
      !
      real(8)                                  :: Beta,mu
      real(8)                                  :: emin,emax,dmu
      integer                                  :: Norb,Nmats,Nkpt
      integer                                  :: ik,iwan1,iwan2,ispin
      integer                                  :: ik,iwan1,iwan2,ispin
      integer                                  :: iter
      !
      !
      write(*,"(A)") "--- set_density_NonInt ---"
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
      Beta = Gmats%Beta
      mu_start = Gmats%mu
      !
      call calc_density(Gmats,n_loc)
      n_iter = get_dens(n_loc)
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
         n_iter = get_dens(Gmats%N_s)
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
      function get_dens() result(dens)
         implicit none
         real(8)      :: dens
         complex(8)   :: dens_C
         dens_C = czero
         if(present(orbs))then
            do iwan1=1,Norb
               if(all(orbs.ne.i))cycle
               do ispin=1,Nspin
                  dens_C = dens_C + n_loc(iwan1,iwan1,ispin)
               enddo
            enddo
         else
            dens_C = trace(sum(n_loc,dim=3))
         endif
         if(aimag(dens_C).gt.eps) stop "Density is complex."
      end function linspace
      !

   end subroutine set_density_NonInt











end module greens_function
