module self_energy

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface calc_VH
      module procedure calc_VH_G
      module procedure calc_VH_N
   end interface calc_VH

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
   public :: calc_sigmaGW
   public :: calc_sigmaGWdc
   public :: read_Sigma_spex
   public :: calc_VH

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the two GW self-energy components.
   !---------------------------------------------------------------------------!
   subroutine calc_sigmaGW(Smats,Gmats,Wmats,Lttc,LDAoffdiag,Smats_C,Smats_X)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use crystal
      use fourier_transforms
      use input_vars, only : Ntau, tau_uniform, cmplxWann, paramagnet
      implicit none
      !
      type(FermionicField),intent(inout)    :: Smats
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Wmats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: LDAoffdiag
      type(FermionicField),intent(out),optional :: Smats_C
      type(FermionicField),intent(out),optional :: Smats_X
      !
      type(FermionicField)                  :: Smats_C_
      type(FermionicField)                  :: Smats_X_
      complex(8),allocatable                :: Sitau(:,:,:)
      complex(8),allocatable                :: Gitau(:,:,:,:,:)
      complex(8),allocatable                :: Witau(:,:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Norb,Nmats
      integer                               :: ik1,ik2,iq,iw,itau,ispin
      integer                               :: i,j,k,l,ib1,ib2
      real                                  :: start,finish
      logical                               :: LDAoffdiag_
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_sigmaGW"
      !
      !
      ! Check on the input Fields
      if(.not.Smats%status) stop "calc_sigmaGW: Smats not properly initialized."
      if(.not.Gmats%status) stop "calc_sigmaGW: Gmats not properly initialized."
      if(.not.Wmats%status) stop "calc_sigmaGW: Wmats not properly initialized."
      if(Smats%Nkpt.eq.0) stop "calc_sigmaGW: Smats k dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_sigmaGW: Gmats k dependent attributes not properly initialized."
      if(Wmats%Nkpt.eq.0) stop "calc_sigmaGW: Wmats k dependent attributes not properly initialized."
      if(.not.allocated(Lttc%kptdif)) stop "calc_sigmaGW: kptdif not allocated."
      if((Lttc%Nkpt_irred.lt.Smats%Nkpt).and.(.not.allocated(Lttc%kptPos))) stop "calc_sigmaGW: kptPos not allocated."
      !
      Norb = Smats%Norb
      Nkpt = Smats%Nkpt
      Beta = Smats%Beta
      Nmats = Smats%Npoints
      Nbp = Norb**2
      !
      if(all([Gmats%Nkpt-Nkpt,Wmats%Nkpt-Nkpt,Lttc%Nkpt-Nkpt].ne.[0,0,0])) stop "calc_sigmaGW: Either Lattice, Gmats or Wmats have different number of k-points with respect to Smats."
      if(all([Gmats%Norb-Norb,Wmats%Nbp-Nbp].ne.[0,0])) stop "calc_sigmaGW: Either Gmats or Wmats have different orbital dimension with respect to Smats."
      if(all([Gmats%Beta-Beta,Wmats%Beta-Beta].ne.[0d0,0d0])) stop "calc_sigmaGW: Either Gmats or Wmats have different Beta with respect to Smats."
      if(all([Gmats%Npoints-Nmats,Wmats%Npoints-Nmats].ne.[0,0])) stop "calc_sigmaGW: Either Gmats or Wmats have different number of Matsubara points with respect to Smats."
      !
      call AllocateFermionicField(Smats_C_,Norb,Nmats,Nkpt=Nkpt,Nsite=Smats%Nsite,Beta=Beta)
      call AllocateFermionicField(Smats_X_,Norb,0,Nkpt=Nkpt,Nsite=Smats%Nsite,Beta=Beta)
      !
      LDAoffdiag_=.true.
      if(present(LDAoffdiag))LDAoffdiag_=LDAoffdiag
      !
      ! Compute Glat(k,tau)
      call cpu_time(start)
      allocate(Gitau(Norb,Norb,Ntau,Nkpt,Nspin));Gitau=czero
      if(cmplxWann)then
         spinloopGWc: do ispin=1,Nspin
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin), &
            asympt_corr=.true.,tau_uniform=tau_uniform)
            if(paramagnet)then
               Gitau(:,:,:,:,Nspin) = Gitau(:,:,:,:,1)
               exit spinloopGWc
            endif
         enddo spinloopGWc
      else
         spinloopGWr: do ispin=1,Nspin
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin), &
            asympt_corr=.true.,tau_uniform=tau_uniform,nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt)
            if(paramagnet)then
               Gitau(:,:,:,:,Nspin) = Gitau(:,:,:,:,1)
               exit spinloopGWr
            endif
         enddo spinloopGWr
      endif
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(k,iw) --> Glat(k,tau) cpu timing:", finish-start
      !
      ! Compute Wlat(q,tau)
      call cpu_time(start)
      allocate(Witau(Nbp,Nbp,Ntau,Nkpt));Witau=czero
      call Bmats2itau(Beta,Wmats%screened,Witau,asympt_corr=.true.,tau_uniform=tau_uniform,Umats_bare=Wmats%bare)
      call cpu_time(finish)
      write(*,"(A,F)") "     Wlat(q,iw) --> Wlat(q,tau) cpu timing:", finish-start
      !
      !Sigma_{m,n}(q,tau) = -Sum_{k,mp,np} W_{(m,mp);(n,np)}(q-k;tau)G_{mp,np}(k,tau)
      call cpu_time(start)
      call clear_attributes(Smats_C_)
      allocate(Sitau(Norb,Norb,Ntau))
      spinloopGW: do ispin=1,Nspin
         do iq=1,Lttc%Nkpt_irred
            !
            Sitau=czero
            !
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(Norb,iq,ispin,Ntau,Lttc,Nkpt,Sitau,Gitau,Witau),&
            !$OMP PRIVATE(itau,ik1,ik2,i,j,k,l,ib1,ib2)
            !$OMP DO
            do itau=1,Ntau
               do k=1,Norb
                  do i=1,Norb
                     !
                     do ik1=1,Nkpt
                        ik2=Lttc%kptdif(iq,ik1)
                        do l=1,Norb
                           do j=1,Norb
                              !
                              ib1 = i + Norb*(j-1)
                              ib2 = k + Norb*(l-1)
                              !
                              Sitau(i,k,itau) = Sitau(i,k,itau) - Gitau(j,l,itau,ik1,ispin)*Witau(ib1,ib2,itau,ik2)/Nkpt
                              !
                           enddo
                        enddo
                     enddo
                     !
                  enddo
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !
            call Fitau2mats_mat(Beta,Sitau,Smats_C_%wks(:,:,:,iq,ispin),tau_uniform=tau_uniform)
            !
         enddo
         if(paramagnet)then
            Smats_C_%wks(:,:,:,:,Nspin) = Smats_C_%wks(:,:,:,:,1)
            exit spinloopGW
         endif
      enddo spinloopGW
      deallocate(Witau,Sitau)
      call cpu_time(finish)
      write(*,"(A,F)") "     Sigma_C(k,iw) cpu timing:", finish-start
      !
      !Sigmax_nm(q) = Sum_kij V_{ni,jm}(q-k)G_ij(k,beta) <=> sigmax(r,r')=-g(r,r',tau=0-)*v(r-r')
      call cpu_time(start)
      call clear_attributes(Smats_X_)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,Ntau,Lttc,Nkpt,Smats_X_,Gitau,Wmats),&
      !$OMP PRIVATE(iq,ispin,ik1,ik2,i,j,k,l,ib1,ib2)
      !$OMP DO
      do iq=1,Lttc%Nkpt_irred
         do ispin=1,Nspin
            !
            do k=1,Norb
               do i=1,Norb
                  !
                  do ik1=1,Nkpt
                     ik2=Lttc%kptdif(iq,ik1)
                     do l=1,Norb
                        do j=1,Norb
                           !
                           ib1 = i + Norb*(j-1)
                           ib2 = k + Norb*(l-1)
                           !
                           Smats_X_%N_ks(i,k,iq,ispin) = Smats_X_%N_ks(i,k,iq,ispin) + Gitau(j,l,Ntau,ik1,ispin)*Wmats%bare(ib1,ib2,ik2)/Nkpt
                           !Smats_X_%N_ks(i,k,iq,ispin) = Smats_X_%N_ks(i,k,iq,ispin) - Gitau(j,l,1,ik1,ispin)*Wmats%bare(ib1,ib2,ik2)/Nkpt
                           !
                        enddo
                     enddo
                  enddo
                  !
               enddo
            enddo
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Gitau)
      call cpu_time(finish)
      write(*,"(A,F)") "     Sigma_X(k) cpu timing:", finish-start
      !
      if(Lttc%Nkpt_irred.lt.Nkpt) then
         !sigma(ik)=sigma(kptp(ik))
         write(*,"(A)") "     Transformation to lda eigenbasis and back."
         if(.not.LDAoffdiag_) write(*,"(A)") "     Removing off-diagonal elements in lda eigenbasis."
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nmats,Lttc,Nkpt,Smats_X_,Smats_C_,LDAoffdiag_),&
         !$OMP PRIVATE(ispin,iw,iq)
         !$OMP DO
         do ispin=1,Nspin
            !
            !rotation to lda eigenbasis
            do iq=1,Lttc%Nkpt_irred
               Smats_X_%N_ks(:,:,iq,ispin) = rotate(Smats_X_%N_ks(:,:,iq,ispin),Lttc%Zk(:,:,iq))
               do iw=1,Nmats
                  Smats_C_%wks(:,:,iw,iq,ispin) = rotate(Smats_C_%wks(:,:,iw,iq,ispin),Lttc%Zk(:,:,iq))
               enddo
            enddo
            !
            !fill up the missing Kpoints
            do iq=1,Nkpt
               Smats_X_%N_ks(:,:,iq,ispin) = Smats_X_%N_ks(:,:,Lttc%kptPos(iq),ispin)
               do iw=1,Nmats
                  Smats_C_%wks(:,:,iw,iq,ispin) = Smats_C_%wks(:,:,iw,Lttc%kptPos(iq),ispin)
               enddo
            enddo
            !
            !remove off-diagonal elements
            if(.not.LDAoffdiag_)then
               do iq=1,Nkpt
                  Smats_X_%N_ks(:,:,iq,ispin) = diag(diagonal(Smats_X_%N_ks(:,:,iq,ispin)))
                  do iw=1,Nmats
                     Smats_C_%wks(:,:,iw,iq,ispin) = diag(diagonal(Smats_C_%wks(:,:,iw,iq,ispin)))
                  enddo
               enddo
            endif
            !
            !rotate back
            do iq=1,Nkpt
               Smats_X_%N_ks(:,:,iq,ispin) = rotate(Smats_X_%N_ks(:,:,iq,ispin),transpose(conjg(Lttc%Zk(:,:,iq))))
               do iw=1,Nmats
                  Smats_C_%wks(:,:,iw,iq,ispin) = rotate(Smats_C_%wks(:,:,iw,iq,ispin),transpose(conjg(Lttc%Zk(:,:,iq))))
               enddo
            enddo
            !
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
      endif
      !
      ! Fill the local attributes
      call FermionicKsum(Smats_C_)
      call FermionicKsum(Smats_X_)
      !
      call join_SigmaCX(Smats,Smats_C_,Smats_X_)
      !
      if(present(Smats_C)) call duplicate(Smats_C,Smats_C_)
      if(present(Smats_X)) call duplicate(Smats_X,Smats_X_)
      !
      call DeallocateFermionicField(Smats_C_)
      call DeallocateFermionicField(Smats_X_)
      !
   end subroutine calc_sigmaGW


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the two local GW self-energy components as Gloc*Wloc.
   !---------------------------------------------------------------------------!
   subroutine calc_sigmaGWdc(Smats_dc,Gmats,Wmats,Smats_Cdc,Smats_Xdc)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use crystal
      use fourier_transforms
      use input_vars, only : Ntau, tau_uniform, paramagnet
      implicit none
      !
      type(FermionicField),intent(inout)    :: Smats_dc
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Wmats
      type(FermionicField),intent(out),optional :: Smats_Cdc
      type(FermionicField),intent(out),optional :: Smats_Xdc
      !
      type(FermionicField)                  :: Smats_Cdc_
      type(FermionicField)                  :: Smats_Xdc_
      complex(8),allocatable                :: Sitau_loc(:,:,:,:)
      complex(8),allocatable                :: Gitau_loc(:,:,:,:)
      complex(8),allocatable                :: Witau_loc(:,:,:)
      !complex(8),allocatable                :: Gmats_lda(:,:),Gitau_lda(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Norb,Nmats
      integer                               :: itau,ispin!,iw
      integer                               :: i,j,k,l,ib1,ib2
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- calc_sigmaGWdc"
      !
      !
      ! Check on the input Fields
      if(.not.Smats_dc%status) stop "calc_sigmaGWdc: Smats_dc not properly initialized."
      if(.not.Gmats%status) stop "calc_sigmaGWdc: Gmats not properly initialized."
      if(.not.Wmats%status) stop "calc_sigmaGWdc: Wmats not properly initialized."
      if(Smats_dc%Nkpt.ne.0) stop "calc_sigmaGWdc: Smats_dc k dependent attributes are supposed to be unallocated."
      !
      Norb = Smats_dc%Norb
      Nkpt = Smats_dc%Nkpt
      Beta = Smats_dc%Beta
      Nmats = Smats_dc%Npoints
      Nbp = Norb**2
      !
      if(all([Gmats%Norb-Norb,Wmats%Nbp-Nbp].ne.[0,0])) stop "calc_sigmaGWdc: Either Gmats or Wmats have different orbital dimension with respect to Smats_dc."
      if(all([Gmats%Beta-Beta,Wmats%Beta-Beta].ne.[0d0,0d0])) stop "calc_sigmaGWdc: Either Gmats or Wmats have different Beta with respect to Smats_dc."
      if(all([Gmats%Npoints-Nmats,Wmats%Npoints-Nmats].ne.[0,0])) stop "calc_sigmaGWdc: Either Gmats or Wmats have different number of Matsubara points with respect to Smats_dc."
      !
      call AllocateFermionicField(Smats_Cdc_,Norb,Nmats,Nsite=Smats_dc%Nsite,Beta=Beta)
      call AllocateFermionicField(Smats_Xdc_,Norb,0,Nsite=Smats_dc%Nsite,Beta=Beta)
      !
      ! Compute Glat(tau) - FT all components
      call cpu_time(start)
      allocate(Gitau_loc(Norb,Norb,Ntau,Nspin));Gitau_loc=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%ws(:,:,:,ispin),Gitau_loc(:,:,:,ispin),asympt_corr=.true.,tau_uniform=tau_uniform)
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(iw), --> Glat(tau) cpu timing:", finish-start
      !
      ! Compute Wlat(tau)
      call cpu_time(start)
      allocate(Witau_loc(Nbp,Nbp,Ntau));Witau_loc=czero
      call Bmats2itau(Beta,Wmats%screened_local,Witau_loc,asympt_corr=.true.,tau_uniform=tau_uniform,Umats_bare=Wmats%bare_local)
      call cpu_time(finish)
      write(*,"(A,F)") "     Wlat(iw) --> Wlat(tau) cpu timing:", finish-start
      !
      !Sigma_{m,n}(q,tau) = -Sum_{k,mp,np} W_{(m,mp);(n,np)}(q-k;tau)G_{mp,np}(k,tau)
      call cpu_time(start)
      allocate(Sitau_loc(Norb,Norb,Ntau,Nspin));Sitau_loc=czero
      spinloopGWdc: do ispin=1,Nspin
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Norb,Ntau,ispin,Sitau_loc,Gitau_loc,Witau_loc),&
         !$OMP PRIVATE(itau,i,j,k,l,ib1,ib2)
         !$OMP DO
         do itau=1,Ntau
            do k=1,Norb
               do i=1,Norb
                  !
                  do l=1,Norb
                     do j=1,Norb
                        !
                        ib1 = i + Norb*(j-1)
                        ib2 = k + Norb*(l-1)
                        !
                        Sitau_loc(i,k,itau,ispin) = Sitau_loc(i,k,itau,ispin) - Gitau_loc(j,l,itau,ispin)*Witau_loc(ib1,ib2,itau)
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         if(paramagnet)then
            Sitau_loc(:,:,:,Nspin) = Sitau_loc(:,:,:,1)
            exit spinloopGWdc
         endif
      enddo spinloopGWdc
      deallocate(Witau_loc)
      !
      call clear_attributes(Smats_Cdc_)
      do ispin=1,Nspin
         call Fitau2mats_mat(Beta,Sitau_loc(:,:,:,ispin),Smats_Cdc_%ws(:,:,:,ispin),tau_uniform=tau_uniform)
      enddo
      deallocate(Sitau_loc)
      call cpu_time(finish)
      write(*,"(A,F)") "     Sigma_Cdc(iw) cpu timing:", finish-start
      !
      !Sigmax_nm(q) = Sum_kij V_{ni,jm}(q-k)G_ij(k,beta) <=> sigmax(r,r')=-g(r,r',tau=0-)*v(r-r')
      call cpu_time(start)
      call clear_attributes(Smats_Xdc_)
      do ispin=1,Nspin
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Norb,Ntau,ispin,Smats_Xdc_,Gitau_loc,Wmats),&
         !$OMP PRIVATE(i,j,k,l,ib1,ib2)
         !$OMP DO
         do k=1,Norb
            do i=1,Norb
               !
               do l=1,Norb
                  do j=1,Norb
                     !
                     ib1 = i + Norb*(j-1)
                     ib2 = k + Norb*(l-1)
                     !
                     Smats_Xdc_%N_s(i,k,ispin) = Smats_Xdc_%N_s(i,k,ispin) + Gitau_loc(j,l,Ntau,ispin)*Wmats%bare_local(ib1,ib2)
                     !
                  enddo
               enddo
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
      enddo
      deallocate(Gitau_loc)
      call cpu_time(finish)
      write(*,"(A,F)") "     Sigma_Xdc cpu timing:", finish-start
      !
      call join_SigmaCX(Smats_dc,Smats_Cdc_,Smats_Xdc_)
      !
      if(present(Smats_Cdc)) call duplicate(Smats_Cdc,Smats_Cdc_)
      if(present(Smats_Xdc)) call duplicate(Smats_Xdc,Smats_Xdc_)
      !
      call DeallocateFermionicField(Smats_Cdc_)
      call DeallocateFermionicField(Smats_Xdc_)
      !
   end subroutine calc_sigmaGWdc


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute Hartree difference with G0W0
   !---------------------------------------------------------------------------!
   subroutine calc_VH_G(VH,density_LDA_spin,Gmats,Umats,sym_constrained)
      !
      use parameters
      use linalg
      use file_io
      use utils_misc
      use utils_fields
      use greens_function, only : calc_density
      use input_vars, only : Nsite, SiteNorb
      use input_vars, only : pathINPUT, VH_type, Nsite
      implicit none
      !
      complex(8),intent(inout)              :: VH(:,:)
      complex(8),intent(in)                 :: density_LDA_spin(:,:,:)
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Umats
      logical,intent(in),optional           :: sym_constrained
      !
      complex(8),allocatable                :: density_LDA(:,:)
      complex(8),allocatable                :: density(:,:),density_spin(:,:,:)
      complex(8),allocatable                :: Vgamma(:,:),Vevec(:,:)
      real(8),allocatable                   :: Veval(:)
      integer,allocatable                   :: orbs(:)
      integer                               :: Norb,Nbp
      integer                               :: isite,shift_N,shift_V
      integer                               :: ib1,ib2
      integer                               :: i,j,k,l
      integer                               :: i_N,j_N,k_N,l_N
      integer                               :: i_V,j_V,k_V,l_V
      logical                               :: filexists,sym_constrained_
      character(len=256)                    :: path
      !
      !
      if(verbose)write(*,"(A)") "---- calc_VH_G"
      !
      !
      ! Check on the input Field
      if(.not.Gmats%status) stop "calc_VH_G: Gmats not properly initialized."
      if(.not.Umats%status) stop "calc_VH_G: Umats not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_VH_G: Gmats k dependent attributes not properly initialized."
      if(Umats%Nkpt.eq.0) stop "calc_VH_G: Umats k dependent attributes not properly initialized."
      if(Umats%iq_gamma.lt.0) stop "calc_VH_G: Umats iq_gamma not defined."
      !
      sym_constrained_=.false.
      if(present(sym_constrained))sym_constrained_=sym_constrained
      !
      if(sym_constrained_.and.(reg(VH_type).eq."Ustatic")) stop "calc_VH_G: Ordering of sym_constrained is not implemented for Ustatic."
      if(sym_constrained_.and.(reg(VH_type).eq."Ubare")) stop "calc_VH_G: Ordering of sym_constrained_ is not implemented for Ubare."
      !
      Norb = Gmats%Norb
      Nbp = Umats%Nbp
      !
      if((Norb**2).ne.Nbp) stop "calc_VH_G: Umats and Gmats have different orbital space."
      call assert_shape(density_LDA_spin,[Norb,Norb,Nspin],"calc_VH_G","density_LDA_spin")
      call assert_shape(VH,[Norb,Norb],"calc_VH_G","VH")
      !
      allocate(density_spin(Norb,Norb,Nspin));density_spin=czero
      call calc_density(Gmats,density_spin)
      allocate(density(Norb,Norb));density=czero
      density = sum(density_spin,dim=3)
      deallocate(density_spin)
      !
      allocate(density_LDA(Norb,Norb));density_LDA=czero
      density_LDA = sum(density_LDA_spin,dim=3)
      !
      select case(reg(VH_type))
         case default
            stop "Available VH_type are: Ubare, Ustatic, Ubare_SPEX, Ustatic_SPEX."
         case("Ubare")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            Vgamma = Umats%bare(:,:,Umats%iq_gamma)
            !
            !rotating to diagonal basis
            allocate(Veval(Nbp));Veval=0d0
            allocate(Vevec(Nbp,Nbp));Vevec=Vgamma
            call eigh(Vevec,Veval)
            !
            !removing largest eigenvalue
            Veval(Nbp)=0d0
            !
            !rotate back
            Vgamma=czero
            Vgamma = rotate(diag(Veval),conjg(transpose(Vevec)))
            !
         case("Ustatic")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            Vgamma = Umats%screened(:,:,1,Umats%iq_gamma)
            !
            !rotating to diagonal basis
            allocate(Veval(Nbp));Veval=0d0
            allocate(Vevec(Nbp,Nbp));Vevec=Vgamma
            call eigh(Vevec,Veval)
            !
            !removing largest eigenvalue
            Veval(Nbp)=0d0
            !
            !rotate back
            Vgamma=czero
            Vgamma = rotate(diag(Veval),conjg(transpose(Vevec)))
            !
         case("Ubare_SPEX")
            !
            if(sym_constrained_) write(*,"(A)") "     Assuming [Norb*[Nsite]] Vgamma arrangement."
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            path = reg(pathINPUT)//"V_nodiv.DAT"
            call inquireFile(reg(path),filexists,verb=verbose)
            call read_Vgamma(-1)
            !
         case("Ustatic_SPEX")
            !
            if(sym_constrained_) write(*,"(A)") "     Assuming [Norb*[Nsite]] Vgamma arrangement."
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            path = reg(pathINPUT)//"V_nodiv.DAT"
            call inquireFile(reg(path),filexists,verb=verbose)
            call read_Vgamma(0)
            !
      end select
      !
      call dump_matrix(Vgamma,reg(pathINPUT),"Vgamma_"//reg(VH_type)//".DAT")
      !
      VH=czero
      if(sym_constrained_)then
         !
         allocate(orbs(SiteNorb(1)));orbs=0
         do i=1,SiteNorb(1)
            orbs(i) = (i-1)*Nsite + 1
         enddo
         write(*,"(A,"//str(SiteNorb(1))//"I4)") "     Orbs: ",orbs
         !
         do isite=1,Nsite
            !
            if(isite.gt.1)then
               if(SiteNorb(isite).ne.SiteNorb(isite-1)) stop "calc_VH_G: sym_constrained not implemented for non-identical sites."
            endif
            !
            ![Nsite*[Norb]] arrangement
            shift_N = SiteNorb(isite)*(isite-1)
            shift_V = isite-1
            !
            do i=1,SiteNorb(isite)
               do j=1,SiteNorb(isite)
                  do k=1,SiteNorb(isite)
                     do l=1,SiteNorb(isite)
                        !
                        i_N = i + shift_N
                        j_N = j + shift_N
                        k_N = k + shift_N
                        l_N = l + shift_N
                        !
                        i_V = orbs(i) + shift_V
                        j_V = orbs(j) + shift_V
                        k_V = orbs(k) + shift_V
                        l_V = orbs(l) + shift_V
                        !
                        ib1 = i_V + Norb*(j_V-1)
                        ib2 = k_V + Norb*(l_V-1)
                        !
                        VH(i_N,j_N) = VH(i_N,j_N) + dreal(density(k_N,l_N)-density_LDA(k_N,l_N)) * dreal(Vgamma(ib1,ib2))
                        !
                     enddo
                  enddo
               enddo
            enddo
            !
         enddo
         deallocate(orbs)
         !
      else
         !
         do i=1,Norb
            do j=1,Norb
               do k=1,Norb
                  do l=1,Norb
                     !
                     ib1 = i + Norb*(j-1)
                     ib2 = k + Norb*(l-1)
                     !
                     VH(i,j) = VH(i,j) + dreal(density(k,l)-density_LDA(k,l)) * dreal(Vgamma(ib1,ib2))
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      endif
      deallocate(density_LDA,density,Vgamma)
      !
      contains
         !
         subroutine read_Vgamma(Vtype)
            use utils_misc
            implicit none
            integer,intent(in)                 :: Vtype
            integer                            :: unit
            real(8)                            :: rdum1,rdum2,rdum3,rdum4
            integer                            :: idum1,idum2,idum3,idum4
            integer                            :: iwan1,iwan2,iwan3,iwan4
            integer                            :: indx1,indx2,limits
            unit = free_unit()
            open(unit,file=reg(path),position="rewind",action="read")
            read(unit,*) idum1 !,'Number of Wannier functions:'
            if(idum1.ne.Norb) stop "read_Vgamma(calc_VH_G): wrong index Norb."
            do limits=1,2
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     do iwan3=1,Norb
                       do iwan4=1,Norb
                           indx1=iwan1+Norb*(iwan2-1)
                           indx2=iwan3+Norb*(iwan4-1)
                           read(unit,*) rdum1,rdum2,idum1,idum2,idum3,idum4,rdum3,rdum4
                           if (idum1.ne.iwan1) stop "read_Vgamma(calc_VH_G): wrong index iwan1."
                           if (idum2.ne.iwan2) stop "read_Vgamma(calc_VH_G): wrong index iwan2."
                           if (idum3.ne.iwan3) stop "read_Vgamma(calc_VH_G): wrong index iwan3."
                           if (idum4.ne.iwan4) stop "read_Vgamma(calc_VH_G): wrong index iwan4."
                           if(dble(Vtype).eq.rdum1) Vgamma(indx1,indx2) = dcmplx(rdum3,rdum4) * H2eV
                       enddo
                     enddo
                  enddo
               enddo
               if(Vtype.eq.-1) exit
            enddo
         end subroutine read_Vgamma
         !
      !
   end subroutine calc_VH_G
   !
   subroutine calc_VH_N(VH,density_LDA_spin,density_spin,Umats,sym_constrained)
      !
      use parameters
      use linalg
      use file_io
      use utils_misc
      use utils_fields
      use greens_function, only : calc_density
      use input_vars, only : Nsite, SiteNorb
      use input_vars, only : pathINPUT, VH_type, Nsite
      implicit none
      !
      complex(8),intent(inout)              :: VH(:,:)
      complex(8),intent(in)                 :: density_LDA_spin(:,:,:)
      complex(8),intent(in)                 :: density_spin(:,:,:)
      type(BosonicField),intent(in)         :: Umats
      logical,intent(in),optional           :: sym_constrained
      !
      complex(8),allocatable                :: density_LDA(:,:)
      complex(8),allocatable                :: density(:,:)
      complex(8),allocatable                :: Vgamma(:,:),Vevec(:,:)
      real(8),allocatable                   :: Veval(:)
      integer,allocatable                   :: orbs(:)
      integer                               :: Norb,Nbp
      integer                               :: isite,shift_N,shift_V
      integer                               :: ib1,ib2
      integer                               :: i,j,k,l
      integer                               :: i_N,j_N,k_N,l_N
      integer                               :: i_V,j_V,k_V,l_V
      logical                               :: filexists,sym_constrained_
      character(len=256)                    :: path
      !
      !
      if(verbose)write(*,"(A)") "---- calc_VH_N"
      !
      !
      ! Check on the input Field
      if(.not.Umats%status) stop "calc_VH_N: Umats not properly initialized."
      if(Umats%Nkpt.eq.0) stop "calc_VH_N: Umats k dependent attributes not properly initialized."
      if(Umats%iq_gamma.lt.0) stop "calc_VH_N: Umats iq_gamma not defined."
      !
      sym_constrained_=.false.
      if(present(sym_constrained))sym_constrained_=sym_constrained
      !
      if(sym_constrained_.and.(reg(VH_type).eq."Ustatic")) stop "calc_VH_N: Ordering of sym_constrained is not implemented for Ustatic."
      if(sym_constrained_.and.(reg(VH_type).eq."Ubare")) stop "calc_VH_N: Ordering of sym_constrained is not implemented for Ubare."
      !
      Norb = size(density_LDA_spin,dim=1)
      Nbp = Umats%Nbp
      !
      if((Norb**2).ne.Nbp) stop "calc_VH_N: Umats and Gmats have different orbital space."
      call assert_shape(density_LDA_spin,[Norb,Norb,Nspin],"calc_VH_N","density_LDA_spin")
      call assert_shape(density_spin,[Norb,Norb,Nspin],"calc_VH_N","density_spin")
      call assert_shape(VH,[Norb,Norb],"calc_VH_N","VH")
      !
      allocate(density(Norb,Norb));density=czero
      density = sum(density_spin,dim=3)
      !
      allocate(density_LDA(Norb,Norb));density_LDA=czero
      density_LDA = sum(density_LDA_spin,dim=3)
      !
      select case(reg(VH_type))
         case default
            stop "Available VH_type are: Ubare, Ustatic, Ubare_SPEX, Ustatic_SPEX."
         case("Ubare")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            Vgamma = Umats%bare(:,:,Umats%iq_gamma)
            !
            !rotating to diagonal basis
            allocate(Veval(Nbp));Veval=0d0
            allocate(Vevec(Nbp,Nbp));Vevec=Vgamma
            call eigh(Vevec,Veval)
            !
            !removing largest eigenvalue
            Veval(Nbp)=0d0
            !
            !rotate back
            Vgamma=czero
            Vgamma = rotate(diag(Veval),conjg(transpose(Vevec)))
            !
         case("Ustatic")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            Vgamma = Umats%screened(:,:,1,Umats%iq_gamma)
            !
            !rotating to diagonal basis
            allocate(Veval(Nbp));Veval=0d0
            allocate(Vevec(Nbp,Nbp));Vevec=Vgamma
            call eigh(Vevec,Veval)
            !
            !removing largest eigenvalue
            Veval(Nbp)=0d0
            !
            !rotate back
            Vgamma=czero
            Vgamma = rotate(diag(Veval),conjg(transpose(Vevec)))
            !
         case("Ubare_SPEX")
            !
            if(sym_constrained_) write(*,"(A)") "     Assuming [Norb*[Nsite]] Vgamma arrangement."
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            path = reg(pathINPUT)//"V_nodiv.DAT"
            call inquireFile(reg(path),filexists,verb=verbose)
            call read_Vgamma(-1)
            !
         case("Ustatic_SPEX")
            !
            if(sym_constrained_) write(*,"(A)") "     Assuming [Norb*[Nsite]] Vgamma arrangement."
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            path = reg(pathINPUT)//"V_nodiv.DAT"
            call inquireFile(reg(path),filexists,verb=verbose)
            call read_Vgamma(0)
            !
      end select
      !
      call dump_matrix(Vgamma,reg(pathINPUT),"Vgamma_"//reg(VH_type)//".DAT")
      !
      allocate(orbs(SiteNorb(1)));orbs=0
      if(sym_constrained_)then
         do i=1,SiteNorb(1)
            orbs(i) = (i-1)*Nsite + 1
         enddo
      else
         do i=1,SiteNorb(1)
            orbs(i) = i
         enddo
      endif
      write(*,"(A,"//str(SiteNorb(1))//"I4)") "     Orbs: ",orbs
      !
      VH=czero
      do isite=1,Nsite
         !
         ![Nsite*[Norb]] arrangement
         shift_N = SiteNorb(isite)*(isite-1)
         shift_V = shift_N
         !
         ![Norb*[Nsite]] arrangement only for Vgamma
         if(sym_constrained_)then
            shift_V = isite-1
            if(isite.gt.1)then
               if(SiteNorb(isite).ne.SiteNorb(isite-1)) stop "calc_VH_N: sym_constrained not implemented for non-identical sites."
            endif
         endif
         !
         do i=1,SiteNorb(isite)
            do j=1,SiteNorb(isite)
               do k=1,SiteNorb(isite)
                  do l=1,SiteNorb(isite)
                     !
                     i_N = i + shift_N
                     j_N = j + shift_N
                     k_N = k + shift_N
                     l_N = l + shift_N
                     !
                     i_V = orbs(i) + shift_V
                     j_V = orbs(j) + shift_V
                     k_V = orbs(k) + shift_V
                     l_V = orbs(l) + shift_V
                     !
                     ib1 = i_V + Norb*(j_V-1)
                     ib2 = k_V + Norb*(l_V-1)
                     !
                     VH(i_N,j_N) = VH(i_N,j_N) + dreal(density(k_N,l_N)-density_LDA(k_N,l_N)) * dreal(Vgamma(ib1,ib2))
                     !
                  enddo
               enddo
            enddo
         enddo
         !
      enddo
      deallocate(orbs,density_LDA,density,Vgamma)
      !
      contains
         !
         subroutine read_Vgamma(Vtype)
            use utils_misc
            implicit none
            integer,intent(in)                 :: Vtype
            integer                            :: unit
            real(8)                            :: rdum1,rdum2,rdum3,rdum4
            integer                            :: idum1,idum2,idum3,idum4
            integer                            :: iwan1,iwan2,iwan3,iwan4
            integer                            :: indx1,indx2,limits
            unit = free_unit()
            open(unit,file=reg(path),position="rewind",action="read")
            read(unit,*) idum1 !,'Number of Wannier functions:'
            if(idum1.ne.Norb) stop "read_Vgamma(calc_VH_N): wrong index Norb."
            do limits=1,2
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     do iwan3=1,Norb
                       do iwan4=1,Norb
                           indx1=iwan1+Norb*(iwan2-1)
                           indx2=iwan3+Norb*(iwan4-1)
                           read(unit,*) rdum1,rdum2,idum1,idum2,idum3,idum4,rdum3,rdum4
                           if (idum1.ne.iwan1) stop "read_Vgamma(calc_VH_N): wrong index iwan1."
                           if (idum2.ne.iwan2) stop "read_Vgamma(calc_VH_N): wrong index iwan2."
                           if (idum3.ne.iwan3) stop "read_Vgamma(calc_VH_N): wrong index iwan3."
                           if (idum4.ne.iwan4) stop "read_Vgamma(calc_VH_N): wrong index iwan4."
                           if(dble(Vtype).eq.rdum1) Vgamma(indx1,indx2) = dcmplx(rdum3,rdum4) * H2eV
                       enddo
                     enddo
                  enddo
               enddo
               if(Vtype.eq.-1) exit
            enddo
         end subroutine read_Vgamma
         !
      !
   end subroutine calc_VH_N


   !---------------------------------------------------------------------------!
   !PURPOSE: Read self-energy from SPEX files.
   !---------------------------------------------------------------------------!
   subroutine read_Sigma_spex(mode,Smats_GoWo,Lttc,save2readable,Vxc_out,pathOUTPUT,recompute)
      !
      use parameters
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: mode
      type(FermionicField),intent(inout)    :: Smats_GoWo
      type(Lattice),intent(inout)           :: Lttc
      logical,intent(in)                    :: save2readable
      complex(8),intent(inout),optional     :: Vxc_out(:,:,:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: recompute
      !
      select case(reg(mode))
         case default
            stop "read_Sigma_spex: Available modes are: Julich, Lund."
         case("Julich")
            !
            call read_Sigma_spex_Julich(Smats_GoWo,Lttc,save2readable,Vxc_out,pathOUTPUT,recompute)
            !
         case("Lund")
            !
            call read_Sigma_spex_Lund(Smats_GoWo,Lttc,save2readable,Vxc_out,pathOUTPUT,recompute)
            !
      end select
      !
   end subroutine read_Sigma_spex
   !
   subroutine read_Sigma_spex_Lund(Smats_GoWo,Lttc,save2readable,Vxc_out,pathOUTPUT,recompute)
      !
      use linalg
      use parameters
      use file_io
      use utils_misc
      use crystal
      use utils_fields
      use input_vars, only : pathINPUT,UseXepsKorder,paramagneticSPEX,XEPSisread
      implicit none
      !
      type(FermionicField),intent(inout)    :: Smats_GoWo
      type(Lattice),intent(inout)           :: Lttc
      logical,intent(in)                    :: save2readable
      complex(8),intent(inout),optional     :: Vxc_out(:,:,:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: recompute
      !
      logical                               :: filexists,ACdone,recompute_,Vxcdone,doVxc
      character(len=256)                    :: path,pathOUTPUT_
      integer                               :: iseg,SigmaSegments
      integer                               :: iq,ik,iw,iw2,ispin,iwan1,iwan2,unit
      integer                               :: Nkpt,Norb,Nmats,Nfreq
      real(8),allocatable                   :: wread(:),wmats(:)
      real                                  :: start,finish
      !Uwan
      integer                               :: Nspin_Uwan,Nkpt_Uwan
      integer                               :: ib_Uwan1,ib_Uwan2,Norb_Uwan
      complex(8),allocatable                :: Uwan(:,:,:,:)
      !spex
      integer                               :: ik_spex,Nspin_spex,Nkpt_irred_spex
      integer                               :: ib,ib_sigma1,ib_sigma2,ispin_spex
      integer                               :: Nspin_spex_old,Nkpt_irred_spex_old
      integer                               :: ib_sigma1_old,ib_sigma2_old
      integer                               :: iseg_old,Nfreq_old,Nkpt_file,Nkpt_file_old
      logical                               :: ldum,lexonly
      real(8)                               :: wS_old,wE_old
      integer,allocatable                   :: NfreqSeg(:)
      real(8)                               :: gamma1,gamma2
      complex(8)                            :: trap
      real(8),allocatable                   :: SigmaX_seg(:,:,:),SigmaX_tmp(:,:)
      complex(8),allocatable                :: SigmaC_seg(:,:,:,:),SigmaC_tmp(:,:,:)
      complex(8),allocatable                :: SigmaC_diag(:,:,:,:)
      complex(8),allocatable                :: Vxc(:,:,:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_Sigma_spex_Lund"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Check on the input Fields
      if(.not.Smats_GoWo%status) stop "read_Sigma_spex_Lund: FermionicField not properly initialized."
      if(.not.Lttc%status) stop "read_Sigma_spex_Lund: Lattice container not properly initialized."
      if(Lttc%Nkpt.ne.Smats_GoWo%Nkpt) stop "read_Sigma_spex_Lund: Lattice has different number of k-points with respect to Smats_GoWo."
      !
      Norb = Smats_GoWo%Norb
      Nkpt = Smats_GoWo%Nkpt
      Nmats = Smats_GoWo%Npoints
      !
      allocate(Vxc(Norb,Norb,Nkpt,Nspin));Vxc=czero
      !
      allocate(wmats(Smats_GoWo%Npoints));wmats=0d0
      wmats = FermionicFreqMesh(Smats_GoWo%Beta,Smats_GoWo%Npoints)
      !
      ! Read XEPS data
      if(.not.XEPSisread)then
         path = reg(pathINPUT)//"XEPS.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         call read_xeps(reg(path),Lttc%kpt,Lttc%Nkpt3,UseXepsKorder, &
         Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,Lttc%iq_gamma,paramagneticSPEX)
      endif
      !
      ! Check if the data on the Matsubara axis are present if(.not.paramagneticSPEX) look also for spin2
      path = reg(pathOUTPUT_)//"SGoWo_w_k_s1.DAT"
      call inquireFile(reg(path),ACdone,hardstop=.false.,verb=verbose)
      recompute_ = .not.ACdone
      if(present(recompute)) recompute_ = recompute .or. recompute_
      !
      ! Check if the Vxc_wann is present
      path = reg(pathINPUT)//"Vxc_k_s1.DAT"
      call inquireFile(reg(path),Vxcdone,hardstop=.false.,verb=verbose)
      doVxc = .not.Vxcdone
      !
      !
      if(doVxc.and.(.not.recompute_))then
         write(*,"(A)")"     Sorry but I can't produce Vxc_wann without reading the self-energy (internal band indexes needed)."
         write(*,"(A)")"     Analytic continuation will be perforemd anyway."
         recompute_ = .true.
      endif
      !
      !
      ! Perform cnalytical continuation on the self-energy on the real axis
      if(recompute_)then
         !
         !---------------------------------------------------------------------!
         !
         write(*,"(A)")"     Performing Analytic continuation to get SigmaGoWo(ik,iw)."
         !
         ! Read UWAN file
         if(Lttc%UseDisentangledBS)then
            path = reg(pathINPUT)//"UWAN_NEW.DAT"
         else
            path = reg(pathINPUT)//"UWAN.DAT"
         endif
         call inquireFile(reg(path),filexists,verb=verbose)
         write(*,"(A)") "     Opening "//reg(path)
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
         read(unit) Nspin_Uwan,Nkpt_Uwan,ib_Uwan1,ib_Uwan2,Norb_Uwan
         if(paramagneticSPEX.and.(Nspin_Uwan.ne.1)) stop "read_Sigma_spex_Lund: UWAN file is not paramagnetic."
         if(Nkpt_Uwan.ne.Nkpt) stop "read_Sigma_spex_Lund: UWAN file has wrong number of k-points (not irreducible)."
         if(Norb_Uwan.ne.Norb) stop "read_Sigma_spex_Lund: UWAN file has wrong orbital dimension."
         write(*,"(A,2I4)") "     The band indexes in the "//reg(path)//" rotation are: ",ib_Uwan1,ib_Uwan2
         allocate(Uwan(ib_Uwan1:ib_Uwan2,Norb,Nkpt,Nspin_Uwan))
         do ispin=1,Nspin_Uwan
            do ik=1,Nkpt
               read(unit) Uwan(:,:,ik,ispin)
            enddo
         enddo
         close(unit)
         !
         !In disentangled case only UWAN_NEW(ib_wan1:ib_wan1+nbasis-1,ibwan1:ibwan1+nbasis-1) is non-zero
         if(Lttc%UseDisentangledBS) ib_Uwan2 = ib_Uwan1 + Norb-1
         !
         ! Look for the Number of Sigma segments.
         SigmaSegments=0
         do iseg=1,99
            path = reg(pathINPUT)//"Sigma_real_"//str(iseg,2)
            call inquireDir(reg(path),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) exit
            SigmaSegments = SigmaSegments + 1
         enddo
         write(*,"(A,1I6)") "     The number of segment of the SPEX self-energy is: ",SigmaSegments
         allocate(NfreqSeg(SigmaSegments));NfreqSeg=0
         !
         ! Look for the Number of Kpoints in each segment (supposed to be the same of Lttc%Nkpt_irred)
         do iseg=1,SigmaSegments
            Nkpt_file = 0
            do ik=1,2000
               path = reg(pathINPUT)//"Sigma_real_"//str(iseg,2)//"/SIGMA.Q"//str(ik,4)//".DAT"
               call inquireFile(reg(path),filexists,hardstop=.false.,verb=verbose)
               if(.not.filexists) exit
               Nkpt_file = Nkpt_file + 1
            enddo
            write(*,"(A,1I4,A,1I6)") "     The number k-points in the segment number: ",iseg," is: ",Nkpt_file
            Nkpt_file_old = Nkpt_file
            if((iseg.gt.1).and.(Nkpt_file.ne.Nkpt_file_old)) stop "read_Sigma_spex_Lund: Number of K-points does not match among segments."
            if(Nkpt_file.ne.Lttc%Nkpt_irred) stop "read_Sigma_spex_Lund: Number of K-points does not match with Nkpt_irred readed from XEPS."
         enddo
         !
         ! Check that all the Sigma parameters are cosistent and that the segments match
         do iseg=1,SigmaSegments
            do ik=1,Lttc%Nkpt_irred
               !
               path = reg(pathINPUT)//"Sigma_real_"//str(iseg,2)//"/SIGMA.Q"//str(ik,4)//".DAT"
               if(verbose)write(*,"(A)") "     Checking "//reg(path)
               call inquireFile(reg(path),filexists,verb=verbose)!redundant control
               !
               unit = free_unit()
               open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
               read(unit) ik_spex,Nspin_spex,Nkpt_irred_spex,ib_sigma1,ib_sigma2,ldum,lexonly,Nfreq
               allocate(wread(Nfreq));wread=0d0
               read(unit) wread
               close(unit)
               if(lexonly) stop "read_Sigma_spex_Lund: lexonly"
               !
               iseg_old = iseg
               Nspin_spex_old = Nspin_spex
               Nkpt_irred_spex_old = Nkpt_irred_spex
               ib_sigma1_old = ib_sigma1
               ib_sigma2_old = ib_sigma2
               Nfreq_old = Nfreq
               wS_old = wread(1)
               wE_old = wread(Nfreq)
               !
               ! Each k-point controls
               if(ib_sigma1.gt.ib_Uwan1) stop "read_Sigma_spex_Lund: ib_sigma1>ib_Uwan1"
               if(ib_sigma2.lt.ib_Uwan2) stop "read_Sigma_spex_Lund: ib_sigma2<ib_Uwan2"
               if(paramagneticSPEX.and.(Nspin_spex.ne.1)) stop "read_Sigma_spex_Lund: Spex self-energy file is not paramagnetic."
               if(ik.ne.ik_spex) stop "read_Sigma_spex_Lund: K-point index in SPEX not match the expected index."
               if(ik.gt.1)then
                  if(ib_sigma1_old.ne.ib_sigma1) stop "read_Sigma_spex_Lund: ib_sigma1 does not match with previous file."
                  if(Nspin_spex_old.ne.Nspin_spex) stop "read_Sigma_spex_Lund: Nspin_spex does not match with previous file."
                  if(Nkpt_irred_spex_old.ne.Nkpt_irred_spex) stop "read_Sigma_spex_Lund: Nkpt_irred_spex does not match with previous file."
                  if(ib_sigma2_old.ne.ib_sigma2) stop "read_Sigma_spex_Lund: ib_sigma2 does not match with previous file."
                  if(iseg.eq.iseg_old)then
                     if(Nfreq_old.ne.Nfreq) stop "read_Sigma_spex_Lund: Nfreq does not match among different k-points same segment."
                     if(wS_old.ne.wread(1)) stop "read_Sigma_spex_Lund: First freq does not match among different k-points same segment."
                     if(wE_old.ne.wread(Nfreq)) stop "read_Sigma_spex_Lund: Last freq does not match among different k-points same segment."
                  else
                     if(dabs(wread(1)-wE_old).gt.eps)then
                        write(*,"(A,2F10.5)") "w_old, w(1):",wE_old,wread(1)
                        write(*,"(A,2I5)") "iq,iseg:",ik,iseg
                        stop "read_Sigma_spex_Lund: Segments in sigma do not match."
                     endif
                  endif
               endif
               deallocate(wread)
               !
            enddo !ik
            !
            !saving the number of frequency points in each segment
            NfreqSeg(iseg) = Nfreq
            !
         enddo !iseg
         write(*,"(A,2I4)") "     The band indexes in the SPEX self-energy are: ",ib_sigma1,ib_sigma2
         !
         allocate(SigmaC_diag(ib_sigma1:ib_sigma2,Nmats,Lttc%Nkpt_irred,Nspin_spex));SigmaC_diag=czero
         !
         ! Perform the transformation on each fraction of the real-axis
         call cpu_time(start)
         do iseg=1,SigmaSegments
            !
            allocate(SigmaX_seg(ib_sigma1:ib_sigma2,Lttc%Nkpt_irred,Nspin_spex));SigmaX_seg=0d0
            allocate(SigmaC_seg(NfreqSeg(iseg),ib_sigma1:ib_sigma2,Lttc%Nkpt_irred,Nspin_spex));SigmaC_seg=czero
            !
            allocate(SigmaX_tmp(3,ib_sigma1:ib_sigma2));SigmaX_tmp=0d0
            allocate(SigmaC_tmp(2,NfreqSeg(iseg),ib_sigma1:ib_sigma2));SigmaC_tmp=czero
            !
            allocate(wread(NfreqSeg(iseg)));wread=0d0
            !
            do iq=1,Lttc%Nkpt_irred
               !
               path = reg(pathINPUT)//"Sigma_real_"//str(iseg,2)//"/SIGMA.Q"//str(iq,4)//".DAT"
               if(verbose)write(*,"(A)") "     Opening "//reg(path)
               call inquireFile(reg(path),filexists,verb=verbose) !redundant control
               !
               unit = free_unit()
               open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
               read(unit) ik_spex,Nspin_spex,Nkpt_irred_spex,ib_sigma1,ib_sigma2,ldum,lexonly,Nfreq
               read(unit) wread
               !
               do ispin=1,Nspin_spex
                  !
                  ! Every iq file contains all the different k component of the self-energy
                  ! so that here ik=1 is the sum over all the internal momentum q
                  do ik=1,Lttc%Nkpt_irred
                     !
                     SigmaX_tmp=0d0;SigmaC_tmp=czero
                     read(unit) SigmaX_tmp
                     read(unit) SigmaC_tmp
                     !
                     do ib=ib_sigma1,ib_sigma2
                        SigmaX_seg(ib,ik,ispin)=SigmaX_seg(ib,ik,ispin)+sum(SigmaX_tmp(:,ib))
                        do iw=1,NfreqSeg(iseg)
                           SigmaC_seg(iw,ib,ik,ispin)=SigmaC_seg(iw,ib,ik,ispin)+sum(SigmaC_tmp(:,iw,ib))
                        enddo !iw
                     enddo !ib
                     !
                  enddo !ik
               enddo !ispin
               close(unit)
               !
            enddo !iq
            deallocate(SigmaX_tmp,SigmaC_tmp)
            !
            !
            do ik=1,Lttc%Nkpt_irred
               do ib=ib_sigma1,ib_sigma2
                  !
                  ! Check that the GoWo self-energy is vanishing at w --> +/-inf
                  if (iseg.eq.1.and.dabs(dimag(SigmaC_seg(1,ib,ik,1))).gt.1.d-3) then
                     write(*,"(A,2E20.12)") "     Warning: ImSigmaC_spex("//str(ik)//",1) orb "//str(ib)//" is > 1.d-4: ",SigmaC_seg(1,ib,ik,1)
                  endif
                  if (iseg.eq.SigmaSegments.and.dabs(dimag(SigmaC_seg(NfreqSeg(iseg),ib,ik,1))).gt.1.d-3) then
                     write(*,"(A,2E20.12)") "     Warning: ImSigmaC_spex("//str(ik)//",Nw) orb "//str(ib)//" is > 1.d-4: ",SigmaC_seg(NfreqSeg(iseg),ib,ik,1)
                  endif
                  !
                  !Calc Sigma along imag axis using
                  !Sigmac(w)=int dw'Gamma(w')/(w-w')
                  !Gamma= - 1/pi Im(Sigmac) Sign(w-mu)
                  !$OMP PARALLEL DEFAULT(NONE),&
                  !$OMP SHARED(Nspin_spex,ib,ik,iseg,Nmats,NfreqSeg,wread,wmats,SigmaC_seg,SigmaC_diag),&
                  !$OMP PRIVATE(ispin,iw,iw2,gamma1,gamma2,trap)
                  !$OMP DO
                  do iw=1,Nmats
                     do iw2=1,NfreqSeg(iseg)-1
                        do ispin=1,Nspin_spex
                           !
                           !try trapetziodal method
                           if(wread(iw2).lt.0.d0) then
                              gamma1 = (+1.d0/pi) * dimag( SigmaC_seg(iw2,ib,ik,ispin)   )
                           else
                              gamma1 = (-1.d0/pi) * dimag( SigmaC_seg(iw2,ib,ik,ispin)   )
                           endif
                           if(wread(iw2+1).lt.0.d0) then
                              gamma2 = (+1.d0/pi) * dimag( SigmaC_seg(iw2+1,ib,ik,ispin) )
                           else
                              gamma2 = (-1.d0/pi) * dimag( SigmaC_seg(iw2+1,ib,ik,ispin) )
                           endif
                           !
                           trap = ( dcmplx(gamma1,0.d0) / (dcmplx(0.d0,wmats(iw)/H2eV) - dcmplx(wread(iw2),0.d0)  )  &
                                  + dcmplx(gamma2,0.d0) / (dcmplx(0.d0,wmats(iw)/H2eV) - dcmplx(wread(iw2+1),0.d0))  )/2.d0
                           !trap=(-dcmplx(gamma*w(iw2),gamma*w_F(iw))/(w_F(iw)**2 + w(iw2)**2) - dcmplx(gamma*w(iw2+1),gamma*w_F(iw))/(w_F(iw)**2 + w(iw2+1)**2))/2.d0
                           SigmaC_diag(ib,iw,ik,ispin) = SigmaC_diag(ib,iw,ik,ispin) + dcmplx(wread(iw2+1)-wread(iw2),0.d0)*trap
                           !
                        enddo !ispin
                     enddo !iw2
                  enddo !iw
                  !$OMP END DO
                  !$OMP END PARALLEL
               enddo !ib
            enddo !ik
            !
            deallocate(wread)
            deallocate(SigmaC_seg)
            if(iseg.lt.SigmaSegments) deallocate(SigmaX_seg)
            !
         enddo !iseg
         !
         !
         call clear_attributes(Smats_GoWo)
         !Sigma=sigmax+sigmac and transform to Wannier basis
         do ispin_spex=1,Nspin_spex
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(ispin_spex,Nkpt,Norb,ib_Uwan1,ib_Uwan2,Nmats,Lttc,Uwan,SigmaC_diag,SigmaX_seg,Smats_GoWo),&
            !$OMP PRIVATE(ik,iwan1,iwan2,iw)
            !$OMP DO
            do ik=1,Nkpt
               do iwan2=1,Norb
                  do iwan1=1,Norb
                     do iw=1,Nmats
                        !
                        Smats_GoWo%wks(iwan1,iwan2,iw,ik,ispin_spex) = &
                        + sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * SigmaC_diag(ib_Uwan1:ib_Uwan2,iw,Lttc%kptPos(ik),ispin_spex) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))  &
                        + sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * SigmaX_seg(ib_Uwan1:ib_Uwan2,Lttc%kptPos(ik),ispin_spex) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))
                        !
                     enddo
                  enddo
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
         enddo
         if(paramagneticSPEX) Smats_GoWo%wks(:,:,:,:,Nspin) = Smats_GoWo%wks(:,:,:,:,1)
         deallocate(Uwan,SigmaC_diag,SigmaX_seg)
         !
         Smats_GoWo%wks = Smats_GoWo%wks * H2eV
         !
         call FermionicKsum(Smats_GoWo)
         !
         call cpu_time(finish)
         write(*,"(A,F)") "     Sigma_GoWo(k,w) --> Sigma_GoWo(k,iw) cpu timing:", finish-start
         !
         ! Print out the transformed stuff
         call dump_FermionicField(Smats_GoWo,reg(pathOUTPUT_),"SGoWo_w",.true.,Lttc%kpt,paramagneticSPEX)
         if(save2readable)call dump_FermionicField(Smats_GoWo,reg(pathOUTPUT_)//"Sigma_imag/","SGoWo_w",.false.,Lttc%kpt,paramagneticSPEX)
         !
         ! Read the Vxc and print it out
         if(doVxc)then
            call read_Vxc_Lund(Vxc,Lttc,ib_sigma1,ib_sigma2,save2readable)
         else
            write(*,"(A)")"     Reading Vxc(k) from "//reg(pathINPUT)
            call read_matrix(Vxc(:,:,:,1),reg(pathINPUT)//"Vxc_k_s1.DAT")
            if(paramagneticSPEX)then
               Vxc(:,:,:,2) = Vxc(:,:,:,1)
            else
               call read_matrix(Vxc(:,:,:,2),reg(pathINPUT)//"Vxc_k_s2.DAT")
            endif
         endif
         !
         !---------------------------------------------------------------------!
         !
      else
         !
         !---------------------------------------------------------------------!
         !
         ! Just read all
         write(*,"(A)")"     Reading SigmaGoWo(k,iw) from "//reg(pathOUTPUT_)//"SGoWo_w_k_s[1,2].DAT"
         call clear_attributes(Smats_GoWo)
         call read_FermionicField(Smats_GoWo,reg(pathOUTPUT_),"SGoWo_w",Lttc%kpt)
         call FermionicKsum(Smats_GoWo)
         !
         write(*,"(A)")"     Reading Vxc(k) from "//reg(pathINPUT)
         call read_matrix(Vxc(:,:,:,1),reg(pathINPUT)//"Vxc_k_s1.DAT")
         if(paramagneticSPEX)then
            Vxc(:,:,:,2) = Vxc(:,:,:,1)
         else
            call read_matrix(Vxc(:,:,:,2),reg(pathINPUT)//"Vxc_k_s2.DAT")
         endif
         !
         !
      endif
      !
      ! Provide output
      if(present(Vxc_out))then
         !
         write(*,"(A)")"     SigmaGoWo(k,iw) and  Vxc(k) are provided separately."
         Vxc_out = Vxc
         !
      else
         !
         write(*,"(A)")"     SigmaGoWo(k,iw)-Vxc(k) is provided in the Fermionic Field."
         do iw=1,Nmats
            Smats_GoWo%wks(:,:,iw,:,:) = Smats_GoWo%wks(:,:,iw,:,:) - Vxc
         enddo
         call FermionicKsum(Smats_GoWo)
         !
      endif
      !
   end subroutine read_Sigma_spex_Lund
   !
   subroutine read_Sigma_spex_Julich(Smats_GoWo,Lttc,save2readable,Vxc_out,pathOUTPUT,recompute)
      !
      use linalg
      use parameters
      use file_io
      use utils_misc
      use crystal
      use utils_fields
      use input_vars, only : pathINPUT,UseXepsKorder,paramagneticSPEX,XEPSisread
      implicit none
      !
      type(FermionicField),intent(inout)    :: Smats_GoWo
      type(Lattice),intent(inout)           :: Lttc
      logical,intent(in)                    :: save2readable
      complex(8),intent(inout),optional     :: Vxc_out(:,:,:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: recompute
      !
      logical                               :: filexists,interpdone,recompute_,Vxcdone,doVxc
      character(len=256)                    :: path,pathOUTPUT_
      integer                               :: ifile,ierr
      integer                               :: ik,iw,ispin,iwan1,iwan2,unit
      integer                               :: Nkpt,Norb,Nmats,Nfreq
      real(8),allocatable                   :: wread(:),wmats(:)

      !Uwan
      integer                               :: Nspin_Uwan,Nkpt_Uwan,Nwan
      integer                               :: ib_Uwan1,ib_Uwan2,Norb_Uwan
      complex(8),allocatable                :: Uwan(:,:,:,:)
      !spex
      type(FermionicField)                  :: Smats_GoWo_C
      real(8)                               :: axispoint,RealS,ImagS
      real(8)                               :: axispoint_prev
      integer                               :: ik_spex,Nspin_spex
      integer                               :: ib_sigma1,ib_sigma2,ispin_spex
      integer                               :: Nfreq_old,ib_read,ik_read
      character(len=10)                     :: commentS
      character(len=12)                     :: commentV
      type(OldBeta)                         :: Beta_Match
      real(8),allocatable                   :: Vxc_read(:,:,:)
      real(8),allocatable                   :: Smats_GoWo_X(:,:,:)
      complex(8),allocatable                :: Vxc(:,:,:,:),Vxc_loc(:,:,:)
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_Sigma_spex_Julich"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Check on the input Fields
      if(.not.Smats_GoWo%status) stop "read_Sigma_spex_Julich: FermionicField not properly initialized."
      if(.not.Lttc%status) stop "read_Sigma_spex_Julich: Lattice container not properly initialized."
      if(Lttc%Nkpt.ne.Smats_GoWo%Nkpt) stop "read_Sigma_spex_Julich: Lattice has different number of k-points with respect to Smats_GoWo."
      !
      Norb = Smats_GoWo%Norb
      Nkpt = Smats_GoWo%Nkpt
      Nmats = Smats_GoWo%Npoints
      !
      allocate(Vxc(Norb,Norb,Nkpt,Nspin));Vxc=czero
      !
      allocate(wmats(Smats_GoWo%Npoints));wmats=0d0
      wmats = FermionicFreqMesh(Smats_GoWo%Beta,Smats_GoWo%Npoints)
      !
      ! Read XEPS data
      if(.not.XEPSisread)then
         path = reg(pathINPUT)//"XEPS.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         call read_xeps(reg(path),Lttc%kpt,Lttc%Nkpt3,UseXepsKorder, &
         Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,Lttc%iq_gamma,paramagneticSPEX)
      endif
      !
      ! Check if the data on the Matsubara axis are present if(.not.paramagneticSPEX) look also for spin2
      path = reg(pathOUTPUT_)//"SGoWo_w_k_s1.DAT"
      call inquireFile(reg(path),interpdone,hardstop=.false.,verb=verbose)
      recompute_ = .not.interpdone
      if(present(recompute)) recompute_ = recompute .or. recompute_
      !
      ! Check if the Vxc_wann is present
      path = reg(pathINPUT)//"Vxc_k_s1.DAT"
      call inquireFile(reg(path),Vxcdone,hardstop=.false.,verb=verbose)
      doVxc = .not.Vxcdone
      !
      !
      if(doVxc.and.(.not.recompute_))then
         write(*,"(A)")"     Sorry but I can't produce Vxc_wann without reading the self-energy (internal band indexes needed)."
         write(*,"(A)")"     Analytic continuation will be perforemd anyway."
         recompute_ = .true.
      endif
      !
      !
      ! Perform cnalytical continuation on the self-energy on the real axis
      if(recompute_)then
         !
         !---------------------------------------------------------------------!
         !
         write(*,"(A)")"     Reading SPEX self-energy on the imaginary axis."
         !
         ! Read UWAN file
         path = reg(pathINPUT)//"UWAN.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         write(*,"(A)") "     Opening "//reg(path)
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
         read(unit) Nspin_Uwan,Nkpt_Uwan,ib_Uwan1,ib_Uwan2,Norb_Uwan
         if(paramagneticSPEX.and.(Nspin_Uwan.ne.1)) stop "read_Sigma_spex_Julich: UWAN file is not paramagnetic."
         if(Nkpt_Uwan.ne.Nkpt) stop "read_Sigma_spex_Julich: UWAN file has wrong number of k-points (not irreducible)."
         if(Norb_Uwan.ne.Norb) stop "read_Sigma_spex_Julich: UWAN file has wrong orbital dimension."
         write(*,"(A,2I4)") "     The band indexes in the "//reg(path)//" rotation are: ",ib_Uwan1,ib_Uwan2
         allocate(Uwan(ib_Uwan1:ib_Uwan2,Norb,Nkpt,Nspin_Uwan))
         do ispin=1,Nspin_Uwan
            do ik=1,Nkpt
               read(unit) Uwan(:,:,ik,ispin)
            enddo
         enddo
         close(unit)
         !
         ! Bands in the energy window used for the disentangling procedure
         Nwan = abs(ib_Uwan2-ib_Uwan1+1)
         !
         ! Look for the Number of Sigma files and internal consistency.
         Nspin_spex = Nspin
         if(paramagneticSPEX) Nspin_spex = 1
         do ifile=1,Nwan*Lttc%Nkpt_irred*Nspin_spex
            !
            path = reg(pathINPUT)//"Sigma_imag/spex.sew."//str(ifile-1,3)
            if(verbose)write(*,"(A)") "     Checking "//reg(path)
            call inquireFile(reg(path),filexists,verb=verbose)
            !
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",action="read",position="rewind")
            call skip_header(unit,5)
            !
            ierr=0
            Nfreq=0
            do while (ierr.eq.0)
               read(unit,*,iostat=ierr) axispoint,RealS,ImagS
               if((axispoint.ge.0d0).and.(axispoint.ne.axispoint_prev)) Nfreq = Nfreq + 1
               axispoint_prev = axispoint
            enddo
            close(unit)
            !
            if(verbose.or.(ifile.eq.1))then
               write(*,"(A,1I6)") "     The number frequency points in the file number "//str(ifile-1,3)//" is: ",Nfreq
            endif
            !
            if(ifile.gt.1.and.Nfreq.ne.Nfreq_old) stop "read_Sigma_spex_Julich: Nfreq does not match among different files."
            Nfreq_old = Nfreq
            !
         enddo
         write(*,"(A)")"     All the SPEX self-energy files have been found."
         !
         !reading the self-energy
         allocate(wread(Nfreq));wread=0d0
         call AllocateFermionicField(Smats_GoWo_C,Nwan,Nfreq,Nkpt=Lttc%Nkpt_irred)
         ifile=0
         do ispin=1,Nspin_spex
            do ik=1,Lttc%Nkpt_irred
               do iwan1=1,Nwan
                  !
                  path = reg(pathINPUT)//"Sigma_imag/spex.sew."//str(ifile,3)
                  if(verbose)write(*,"(A)") "     Reading "//reg(path)
                  call inquireFile(reg(path),filexists,verb=verbose)
                  !
                  unit = free_unit()
                  open(unit,file=reg(path),form="formatted",action="read",position="rewind")
                  read(unit,*) !# Expectation value of the self-energy on the imaginary frequency axis.
                  read(unit,*) !#
                  read(unit,"(A,2I5)") commentS,ik_spex,ispin_spex
                  read(unit,"(A,I)") commentS,ib_read
                  read(unit,*)
                  !
                  if(ifile.eq.0) ib_sigma1 = ib_read
                  if(ifile.eq.(Lttc%Nkpt_irred*Nwan-1)) ib_sigma2 = ib_read
                  if(ik.ne.ik_spex) stop "read_Sigma_spex_Julich: K-point index in SPEX not match the expected index."
                  if(ispin.ne.ispin_spex) stop "read_Sigma_spex_Julich: spin index in SPEX not match the expected index."
                  if(paramagneticSPEX.and.(ispin_spex.gt.1)) stop "read_Sigma_spex_Julich: Spex self-energy file is not paramagnetic."
                  !
                  ierr=0
                  iw=1
                  do while (ierr.eq.0)
                     read(unit,*,iostat=ierr) axispoint,RealS,ImagS
                     if((axispoint.ge.0d0).and.(axispoint.ne.axispoint_prev))then
                        wread(iw) = axispoint * H2eV
                        Smats_GoWo_C%wks(iwan1,iwan1,iw,ik,ispin) = dcmplx(RealS,ImagS) * H2eV
                        iw = iw + 1
                     endif
                     axispoint_prev = axispoint
                  enddo
                  close(unit)
                  !
                  ifile = ifile + 1
                  !
               enddo
            enddo
            if(ib_sigma1.gt.ib_Uwan1) stop "read_Sigma_spex_Julich: ib_sigma1>ib_Uwan1"
            if(ib_sigma2.lt.ib_Uwan2) stop "read_Sigma_spex_Julich: ib_sigma2<ib_Uwan2"
            if(abs(ib_sigma2-ib_sigma1+1).ne.Nwan) stop "read_Sigma_spex_Julich: wannier dimension in SPEX not match the expected one from UWAN.DAT."
            write(*,"(A,2I4)") "     The band indexes in the SPEX self-energy are: ",ib_sigma1,ib_sigma2
         enddo
         if(paramagneticSPEX) Smats_GoWo_C%wks(:,:,:,:,Nspin) = Smats_GoWo_C%wks(:,:,:,:,1)
         !
         !Interpolate to the target beta
         write(*,"(A)")"     Interpolating to current beta."
         Beta_Match%Nmats_old = Nfreq
         Beta_Match%Nmats_new = Nmats
         Beta_Match%Beta_new = Smats_GoWo%Beta
         call interpolate2Beta_Fermionic(Smats_GoWo_C,Beta_Match,"lat",.false.,wmats_in=wread)
         deallocate(wread)
         !
         !read Smats_GoWo_X and Vxc
         allocate(Vxc_read(Nwan,Lttc%Nkpt_irred,Nspin_spex));Vxc_read=0d0
         allocate(Smats_GoWo_X(Nwan,Lttc%Nkpt_irred,Nspin_spex));Smats_GoWo_X=0d0
         !
         path = reg(pathINPUT)//"Sigma_imag/spex.out"
         if(verbose)write(*,"(A)") "     Reading "//reg(path)
         call inquireFile(reg(path),filexists,verb=verbose)
         !
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",action="read",position="rewind")
         ik=0
         do while (ierr.eq.0)
            read(unit,"(1A12,1I,1A10)",iostat=ierr) commentV,ik_read,commentS
            if(reg(commentV).eq."### K POINT:")then
               ik = ik + 1
               if(ik.ne.ik_read) stop "read_Sigma_spex_Julich: K-point index in spex.out not match the expected index."
               call skip_header(unit,6)
               do iwan1=1,Nwan
                  read(unit,*) ib_read, Vxc_read(iwan1,ik,1), Smats_GoWo_X(iwan1,ik,1)
                  read(unit,*)
                  if(ib_read.ne.ib_sigma1+(iwan1-1)) stop "read_Sigma_spex_Julich:orbital index in spex.out not match the expected index."
               enddo
            endif
         enddo
         close(unit)
         !
         !It should be already in eV
         !Vxc_read = Vxc_read * H2eV
         !Smats_GoWo_X = Smats_GoWo_X * H2eV
         !
         !Transform to Wannier basis
         call clear_attributes(Smats_GoWo)
         !Sigma=sigmax+sigmac and transform to Wannier basis
         do ispin_spex=1,Nspin_spex
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(ispin_spex,Nkpt,Norb,ib_Uwan1,ib_Uwan2,Nmats,Lttc,Uwan,Smats_GoWo_C,Smats_GoWo_X,Smats_GoWo,&
            !$OMP doVxc,Vxc,Vxc_read),&
            !$OMP PRIVATE(ik,iwan1,iwan2,iw)
            !$OMP DO
            do ik=1,Nkpt
               do iwan2=1,Norb
                  do iwan1=1,Norb
                     do iw=1,Nmats
                        !
                        Smats_GoWo%wks(iwan1,iwan2,iw,ik,ispin_spex) = &
                        + sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * diagonal(Smats_GoWo_C%wks(:,:,iw,Lttc%kptPos(ik),ispin_spex)) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))  &
                        + sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * Smats_GoWo_X(:,Lttc%kptPos(ik),ispin_spex) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))
                        !
                     enddo
                     !
                     if(doVxc) Vxc(iwan1,iwan2,ik,ispin_spex) = sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * Vxc_read(ib_Uwan1:ib_Uwan2,Lttc%kptPos(ik),ispin_spex) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))
                     !
                  enddo
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
         enddo
         deallocate(Uwan,Vxc_read,Smats_GoWo_X)
         call DeallocateFermionicField(Smats_GoWo_C)
         if(paramagneticSPEX)then
            Smats_GoWo%wks(:,:,:,:,Nspin) = Smats_GoWo%wks(:,:,:,:,1)
            Vxc(:,:,:,Nspin) = Vxc(:,:,:,1)
         endif
         !
         call FermionicKsum(Smats_GoWo)
         !
         ! Print out the transformed stuff
         call dump_FermionicField(Smats_GoWo,reg(pathOUTPUT_),"SGoWo_w",.true.,Lttc%kpt,paramagneticSPEX)
         if(save2readable)call dump_FermionicField(Smats_GoWo,reg(pathOUTPUT_)//"Sigma_imag/","SGoWo_w",.false.,Lttc%kpt,paramagneticSPEX)
         !
         ! Read or print print Vxc out
         if(doVxc)then
            do ispin=1,Nspin_Uwan
               call dump_matrix(Vxc(:,:,:,ispin),.true.,reg(pathINPUT),"Vxc",ispin=ispin)
               if(save2readable)call dump_matrix(Vxc(:,:,:,ispin),.false.,reg(pathINPUT)//"Vxc/","Vxc",ispin=ispin)
            enddo
            !
            allocate(Vxc_loc(Norb,Norb,Nspin));Vxc_loc=czero
            do ik=1,Nkpt
               Vxc_loc = Vxc_loc + Vxc(:,:,ik,:)/Nkpt
            enddo
            call dump_matrix(Vxc_loc,reg(pathINPUT),"Vxc",(Nspin_Uwan.eq.2))
            deallocate(Vxc_loc)
         else
            write(*,"(A)")"     Reading Vxc(k) from "//reg(pathINPUT)
            call read_matrix(Vxc(:,:,:,1),reg(pathINPUT)//"Vxc_k_s1.DAT")
            if(paramagneticSPEX)then
               Vxc(:,:,:,2) = Vxc(:,:,:,1)
            else
               call read_matrix(Vxc(:,:,:,2),reg(pathINPUT)//"Vxc_k_s2.DAT")
            endif
         endif
         !
         !---------------------------------------------------------------------!
         !
      else
         !
         !---------------------------------------------------------------------!
         !
         ! Just read all
         write(*,"(A)")"     Reading SigmaGoWo(k,iw) from "//reg(pathOUTPUT_)//"SGoWo_w_k_s[1,2].DAT"
         call clear_attributes(Smats_GoWo)
         call read_FermionicField(Smats_GoWo,reg(pathOUTPUT_),"SGoWo_w",Lttc%kpt)
         call FermionicKsum(Smats_GoWo)
         !
         write(*,"(A)")"     Reading Vxc(k) from "//reg(pathINPUT)
         call read_matrix(Vxc(:,:,:,1),reg(pathINPUT)//"Vxc_k_s1.DAT")
         if(paramagneticSPEX)then
            Vxc(:,:,:,2) = Vxc(:,:,:,1)
         else
            call read_matrix(Vxc(:,:,:,2),reg(pathINPUT)//"Vxc_k_s2.DAT")
         endif
         !
         !
      endif
      !
      ! Provide output
      if(present(Vxc_out))then
         !
         write(*,"(A)")"     SigmaGoWo(k,iw) and  Vxc(k) are provided separately."
         Vxc_out = Vxc
         !
      else
         !
         write(*,"(A)")"     SigmaGoWo(k,iw)-Vxc(k) is provided in the Fermionic Field."
         do iw=1,Nmats
            Smats_GoWo%wks(:,:,iw,:,:) = Smats_GoWo%wks(:,:,iw,:,:) - Vxc
         enddo
         call FermionicKsum(Smats_GoWo)
         !
      endif
      !
   end subroutine read_Sigma_spex_Julich


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the exchange potential from gwa file.
   !---------------------------------------------------------------------------!
   subroutine read_Vxc_Lund(Vxc,Lttc,ib_sigma1,ib_sigma2,save2readable)
      !
      use linalg
      use parameters
      use file_io
      use utils_misc
      use crystal
      use input_vars, only : pathINPUT,UseXepsKorder,paramagneticSPEX,XEPSisread
      implicit none
      !
      complex(8),intent(inout)              :: Vxc(:,:,:,:)
      type(Lattice),intent(inout)           :: Lttc
      integer,intent(in)                    :: ib_sigma1,ib_sigma2
      logical,intent(in)                    :: save2readable
      !
      logical                               :: filexists
      character(len=256)                    :: path
      integer                               :: ik,ispin,iwan1,iwan2,unit
      integer                               :: Nkpt,Norb,ispin_spex
      !Uwan
      integer                               :: Nspin_Uwan,Nkpt_Uwan
      integer                               :: ib_Uwan1,ib_Uwan2,Norb_Uwan
      complex(8),allocatable                :: Uwan(:,:,:,:)
      !gwa
      integer                               :: idum,l,i,j
      integer                               :: lcutd,ntypd,nlod,neigd,ncent
      integer                               :: n_stride,n_start,n_size,n_rank
      real(8)                               :: rdum,latpar,lat2(3,3)
      logical                               :: invs,l_soc
      integer(8)                            :: irecl
      complex(8),allocatable                :: vxcmat(:,:)
      !DISENT_EVEC
      integer                               :: unit_dis
      integer                               :: Nspin_disent,Nkpt_irred_disent,Norb_disent
      integer                               :: ib_Dwan1,ib_Dwan2
      complex(8),allocatable                :: dis_evec(:,:),cmat(:,:)
      !eig and vxcfull
      integer                               :: unit_eig,unit_vxc
      integer                               :: nband,nrec
      real(8),allocatable                   :: vxc_diag(:,:,:)
      complex(8),allocatable                :: Vxc_loc(:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- read_Vxc_Lund"
      !
      !
      ! Check on the input Vxc matrix
      Norb = size(Vxc,dim=1)
      if(Norb.ne.size(Vxc,dim=2)) stop "read_Vxc_Lund: Vxc is not a square matrix."
      Nkpt = size(Vxc,dim=3)
      if(.not.Lttc%status) stop "read_Vxc_Lund: Lattice container not properly initialized."
      if(Lttc%Nkpt.ne.Nkpt) stop "read_Vxc_Lund: Lattice has different number of k-points with respect to Vxc."
      !
      ! Read XEPS data
      if(.not.XEPSisread)then
         path = reg(pathINPUT)//"XEPS.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         call read_xeps(reg(path),Lttc%kpt,Lttc%Nkpt3,UseXepsKorder, &
         Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,Lttc%iq_gamma,paramagneticSPEX)
      endif
      !
      ! Read UWAN file
      if(Lttc%UseDisentangledBS)then
         path = reg(pathINPUT)//"UWAN_NEW.DAT"
      else
         path = reg(pathINPUT)//"UWAN.DAT"
      endif
      call inquireFile(reg(path),filexists,verb=verbose)
      write(*,"(A)") "     Opening "//reg(path)
      unit = free_unit()
      open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
      read(unit) Nspin_Uwan,Nkpt_Uwan,ib_Uwan1,ib_Uwan2,Norb_Uwan
      if(paramagneticSPEX.and.(Nspin_Uwan.ne.1)) stop "read_Vxc_Lund: UWAN file is not paramagnetic."
      if(Nkpt_Uwan.ne.Nkpt) stop "read_Vxc_Lund: UWAN file has wrong number of k-points (not irreducible)."
      if(Norb_Uwan.ne.Norb) stop "read_Vxc_Lund: UWAN file has wrong orbital dimension."
      write(*,"(A,2I4)") "     The band indexes in the "//reg(path)//" rotation are: ",ib_Uwan1,ib_Uwan2
      allocate(Uwan(ib_Uwan1:ib_Uwan2,Norb,Nkpt,Nspin_Uwan))
      do ispin=1,Nspin_Uwan
         do ik=1,Nkpt
            read(unit) Uwan(:,:,ik,ispin)
         enddo
      enddo
      close(unit)
      !
      !In disentangled case only UWAN_NEW(ib_wan1:ib_wan1+nbasis-1,ibwan1:ibwan1+nbasis-1) is non-zero
      if(Lttc%UseDisentangledBS) ib_Uwan2 = ib_Uwan1 + Norb-1
      !
      ! Read gwa file
      path = reg(pathINPUT)//"gwa.DAT"
      call inquireFile(reg(path),filexists,verb=verbose)
      write(*,"(A)") "     Opening "//reg(path)
      unit = free_unit()
      open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
      read(unit) idum,ncent,ntypd,l,nlod
      read(unit) (idum,i=1,ncent),(idum,i=1,ntypd),(idum,i=1,ntypd*(l+1)),&
           (rdum,i=1,ntypd),(rdum,i=1,3*ncent),&
           latpar,lat2,rdum,neigd,lcutd
      read(unit) invs,l_soc
      read(unit,iostat=i) irecl,n_start,n_stride,n_rank,n_size
      close(unit)
      allocate(vxcmat(neigd,neigd));vxcmat=czero
      !
      ! Read DISENT_EVEC file
      if(Lttc%UseDisentangledBS) then
         path = reg(pathINPUT)//"DISENT_EVEC.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         write(*,"(A)") "     Opening "//reg(path)
         unit_dis = free_unit()
         open(unit_dis,file=reg(path),form="unformatted",action="read",position="rewind")
         read(unit_dis) Nspin_disent,Nkpt_irred_disent,Norb_disent,ib_Dwan1,ib_Dwan2
         if (Nspin_disent.ne.Nspin_Uwan) stop "DISENT_EVEC.DAT: wrong index nspin."
         if (Nkpt_irred_disent.ne.Lttc%Nkpt_irred) stop "read_Vxc_Lund: DISENT_EVEC.DAT: wrong index nkpt1."
         if (Norb_disent.ne.Norb) stop "read_Vxc_Lund: DISENT_EVEC.DAT: wrong index nwan."
         if (ib_Dwan1.ne.ib_Uwan1) stop "read_Vxc_Lund: DISENT_EVEC.DAT: wrong index ib_wan1."
         if (ib_Dwan2.ne.ib_Uwan2) then
            write(*,"(A,2I4)") "     ib_Uwan2,ib_Dwan2: ",ib_Uwan2, ib_Dwan2
         endif
         allocate(dis_evec(ib_Uwan1:ib_Dwan2,ib_Uwan1:ib_Dwan2));dis_evec=czero
         allocate(cmat(ib_Uwan1:ib_Dwan2,ib_Uwan1:ib_Dwan2));cmat=czero
      endif
      !
      ! Read eig and vxcfull files
      path = reg(pathINPUT)//"eig.DAT"
      call inquireFile(reg(path),filexists,verb=verbose)
      write(*,"(A)") "     Opening "//reg(path)
      unit_eig = free_unit()
      open(unit_eig,file=reg(path),form='unformatted',access='direct',action='read',recl=irecl)
      !
      path = reg(pathINPUT)//"vxcfull.DAT"
      call inquireFile(reg(path),filexists,verb=verbose)
      write(*,"(A)") "     Opening "//reg(path)
      unit_vxc = free_unit()
      open(unit_vxc,file=reg(path),form='unformatted',action='read')
      !
      ! Read diagonal Vxc
      allocate(vxc_diag(ib_sigma1:ib_sigma2,Lttc%Nkpt_irred,Nspin_Uwan));vxc_diag=0d0
      do ispin=1,Nspin_Uwan
         do ik=1,Lttc%Nkpt_irred
            nrec = Lttc%Nkpt_irred*(ispin-1)+(ik-1)*n_stride+n_start   ! according to
            nrec = n_size*(nrec-1) + n_rank + 1                        ! eigen.f (Fleur)
            read(unit_eig,rec=nrec) ((rdum,l=0,lcutd),i=1,ntypd),rdum,rdum, &
                                    ((rdum,l=1,nlod),i=1,ntypd),(rdum,i=1,3),rdum,nband
            !write(*,*) 'ispin,ik,neigd,nband=',ispin,ik,neigd,nband
            if (nband.lt.1.or.nband.gt.neigd) then
               write(*,"(A,4I5)") "ispin,ik,neigd,nband: ",ispin,ik,neigd,nband
               stop "read_Vxc_Lund: nband is out of bound."
            endif
            read(unit_vxc) ((vxcmat(i,j),i=1,j),j=1,nband)
            do j=1,nband
               do i=1,j-1
                  vxcmat(j,i)=conjg(vxcmat(i,j))
               enddo
            enddo
            !
            do i=ib_sigma1,ib_sigma2
               vxc_diag(i,ik,ispin)=dble(vxcmat(i,i))
            enddo
            if(Lttc%UseDisentangledBS) then
               read(unit_dis) dis_evec(:,:)
               cmat(:,:)=matmul(vxcmat(ib_Uwan1:ib_Dwan2,ib_Uwan1:ib_Dwan2),dis_evec)
               do i=ib_Uwan1,ib_Uwan2
                  vxc_diag(i,ik,ispin)=dble(dot_product(dis_evec(:,i),cmat(:,i)))
                  !F.N:
                  !write(*,*) 'vxc i=',i,'ik=',ik
                  !write(*,*) vxc(i,ik,ispin)
               enddo
            endif
            !
         enddo ! ik
      enddo ! ispin
      close(unit_eig)
      close(unit_vxc)
      !
      if(Lttc%UseDisentangledBS) then
         close(unit_dis)
         deallocate(cmat,dis_evec)
      endif
      deallocate(vxcmat)
      !
      !
      Vxc=czero
      ! Transform to Wannier basis
      !Sigma=sigmax+sigmac and transform to Wannier basis
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(paramagneticSPEX,Norb,Nkpt,ib_Uwan1,ib_Uwan2,Lttc,Uwan,vxc_diag,Vxc),&
      !$OMP PRIVATE(ispin,ispin_spex,iwan1,iwan2)
      !$OMP DO
      do ispin=1,Nspin
         do ik=1,Nkpt
            do iwan2=1,Norb
               do iwan1=1,Norb
                  !
                  ispin_spex=1
                  if(.not.paramagneticSPEX)ispin_spex=ispin
                  !
                  Vxc(iwan1,iwan2,ik,ispin) = sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * vxc_diag(ib_Uwan1:ib_Uwan2,Lttc%kptPos(ik),ispin_spex) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))
                  !
                  Vxc(iwan1,iwan2,ik,ispin) = Vxc(iwan1,iwan2,ik,ispin) * H2eV
                  !
               enddo
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(vxc_diag,Uwan)
      !
      do ispin=1,Nspin_Uwan
         call dump_matrix(Vxc(:,:,:,ispin),.true.,reg(pathINPUT),"Vxc",ispin=ispin)
         if(save2readable)call dump_matrix(Vxc(:,:,:,ispin),.false.,reg(pathINPUT)//"Vxc/","Vxc",ispin=ispin)
      enddo
      !
      allocate(Vxc_loc(Norb,Norb,Nspin));Vxc_loc=czero
      do ik=1,Nkpt
         Vxc_loc = Vxc_loc + Vxc(:,:,ik,:)/Nkpt
      enddo
      call dump_matrix(Vxc_loc,reg(pathINPUT),"Vxc",(Nspin_Uwan.eq.2))
      deallocate(Vxc_loc)
      !
   end subroutine read_Vxc_Lund


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate to a new frequency mesh a Fermionic field
   !         copied from post_processing.f90
   !---------------------------------------------------------------------------!
   subroutine interpolate2Beta_Fermionic(G,Beta_Match,mode,offDiag,wmats_in)
      !
      use parameters
      use utils_misc
      use utils_fields
      use interactions
      use input_vars, only: Nsite, SiteNorb, SiteOrbs
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      type(OldBeta),intent(in)              :: Beta_Match
      character(len=*),intent(in)           :: mode
      logical,intent(in)                    :: offDiag
      real(8),intent(in),optional           :: wmats_in(:)
      !
      type(FermionicField)                  :: G_old
      integer                               :: isite,ik,iw
      integer                               :: i,j,ispin
      integer                               :: i_lat,j_lat
      logical                               :: LocalOnly,replace
      real(8),allocatable                   :: wmats_new(:),wmats_old(:)
      real(8),allocatable                   :: ReGf(:),ImGf(:)
      real(8),allocatable                   :: ReD(:),ImD(:)
      !
      !
      if(verbose)write(*,"(A)") "---- interpolate2Beta_Fermionic"
      !
      !
      if(.not.G%status) stop "interpolate2Beta_Fermionic: FermionicField not properly initialized."
      if(G%Beta.ne.Beta_Match%Beta_old) stop "interpolate2Beta_Fermionic: Fermionic field Beta is different from the expected one."
      if(G%Npoints.ne.Beta_Match%Nmats_old) stop "interpolate2Beta_Fermionic: Fermionic field Npoints is different from the expected one."
      !
      LocalOnly=.true.
      if(G%Nkpt.ne.0)LocalOnly=.false.
      !
      allocate(wmats_new(Beta_Match%Nmats_new)); wmats_new=BosonicFreqMesh(Beta_Match%Beta_new,Beta_Match%Nmats_new)
      if(present(wmats_in))then
         wmats_old = wmats_in
      else
         allocate(wmats_old(Beta_Match%Nmats_old)); wmats_old=BosonicFreqMesh(Beta_Match%Beta_old,Beta_Match%Nmats_old)
      endif
      !
      call duplicate(G_old,G)
      call DeallocateFermionicField(G)
      call AllocateFermionicField(G,G_old%Norb,Beta_Match%Nmats_new,Nkpt=G_old%Nkpt,Nsite=G_old%Nsite,Beta=Beta_Match%Beta_new,mu=G_old%mu)
      !
      select case(reg(mode))
         case default
            !
            stop "interpolate2Beta_Fermionic: Available Modes are: imp, lat."
            !
         case("imp")
            !
            allocate(ReD(Beta_Match%Nmats_old)) ;allocate(ImD(Beta_Match%Nmats_old))
            allocate(ReGf(Beta_Match%Nmats_new));allocate(ImGf(Beta_Match%Nmats_new))
            !
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(G,G_old,Nsite,SiteNorb,SiteOrbs),&
            !$OMP SHARED(Beta_Match,offDiag,wmats_new,wmats_old),&
            !$OMP PRIVATE(ispin,isite,i,j,i_lat,j_lat,replace,iw,ReD,ReGf,ImD,ImGf)
            !$OMP DO
            do isite=1,Nsite
               !
               do i=1,SiteNorb(isite)
                  do j=1,SiteNorb(isite)
                     !
                     i_lat = SiteOrbs(isite,i)
                     j_lat = SiteOrbs(isite,j)
                     !
                     replace = i_lat.eq.j_lat
                     if(offDiag) replace = .true.
                     !
                     if(replace)then
                        !
                        do ispin=1,Nspin
                           ReD=0d0;ImD=0d0;ReGf=0d0;ImGf=0d0
                           call nspline(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD)
                           call nspline(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD)
                           do iw=1,Beta_Match%Nmats_new
                              call splint(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD,wmats_new(iw),ReGf(iw))
                              call splint(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD,wmats_new(iw),ImGf(iw))
                           enddo
                           do iw=1,Beta_Match%Nmats_new
                              G%ws(i_lat,j_lat,iw,ispin) = dcmplx(ReGf(iw),ImGf(iw))
                           enddo
                        enddo
                        !
                     endif
                     !
                  enddo
               enddo
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(ReD,ImD,ReGf,ImGf)
            call DeallocateFermionicField(G_old)
            !
         case("lat")
            !
            allocate(ReD(Beta_Match%Nmats_old)) ;allocate(ImD(Beta_Match%Nmats_old))
            allocate(ReGf(Beta_Match%Nmats_new));allocate(ImGf(Beta_Match%Nmats_new))
            !
            !
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(G,G_old,LocalOnly),&
            !$OMP SHARED(Beta_Match,offDiag,wmats_new,wmats_old),&
            !$OMP PRIVATE(ispin,i,j,i_lat,j_lat,replace,iw,ik,ReD,ReGf,ImD,ImGf)
            !$OMP DO
            do i_lat=1,G%Norb
               do j_lat=1,G%Norb
                  !
                  replace = i_lat.eq.j_lat
                  if(offDiag) replace = .true.
                  !
                  if(replace)then
                     !
                     do ispin=1,Nspin
                        !
                        if(LocalOnly)then
                           !
                           ReD=0d0;ImD=0d0;ReGf=0d0;ImGf=0d0
                           call nspline(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD)
                           call nspline(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD)
                           do iw=1,Beta_Match%Nmats_new
                              call splint(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD,wmats_new(iw),ReGf(iw))
                              call splint(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD,wmats_new(iw),ImGf(iw))
                           enddo
                           do iw=1,Beta_Match%Nmats_new
                              G%ws(i_lat,j_lat,iw,ispin) = dcmplx(ReGf(iw),ImGf(iw))
                           enddo
                           !
                        else
                           !
                           do ik=1,G%Nkpt
                              ReD=0d0;ImD=0d0;ReGf=0d0;ImGf=0d0
                              call nspline(wmats_old, real(G_old%wks(i_lat,j_lat,:,ik,ispin)),ReD)
                              call nspline(wmats_old,aimag(G_old%wks(i_lat,j_lat,:,ik,ispin)),ImD)
                              do iw=1,Beta_Match%Nmats_new
                                 call splint(wmats_old, real(G_old%wks(i_lat,j_lat,:,ik,ispin)),ReD,wmats_new(iw),ReGf(iw))
                                 call splint(wmats_old,aimag(G_old%wks(i_lat,j_lat,:,ik,ispin)),ImD,wmats_new(iw),ImGf(iw))
                              enddo
                              do iw=1,Beta_Match%Nmats_new
                                 G%wks(i_lat,j_lat,iw,ik,ispin) = dcmplx(ReGf(iw),ImGf(iw))
                              enddo
                           enddo
                           !
                        endif
                        !
                     enddo
                     !
                  endif
                  !
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(ReD,ImD,ReGf,ImGf)
            call DeallocateFermionicField(G_old)
            if(.not.LocalOnly) call FermionicKsum(G)
            !
      end select
      !
   end subroutine interpolate2Beta_Fermionic

end module self_energy
