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
      module procedure calc_density_loc                                         !(FermionicField,n_loc[Norb,Norb,Nspin])
      module procedure calc_density_Kdep                                        !(FermionicField,Lttc,n_k[Norb,Norb,Nkpt,Nspin])
   end interface calc_density


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
   !public :: set_density
   !public :: calc_GreensFunction - interface between interacting and non-int

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
      call assert_shape(n_loc,[Norb,Norb,Nspin],"calc_density","n_loc")
      !
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
      call assert_shape(n_k,[Norb,Norb,Nkpt,Nspin],"calc_density","n_k")
      !
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
   !PURPOSE: Compute the chemical potential given the Green's function
   !---------------------------------------------------------------------------!




end module greens_function
