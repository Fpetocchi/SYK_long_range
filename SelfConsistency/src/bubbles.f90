module bubbles

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface calc_PiGG
      module procedure calc_Pi_GoGo                                             ![BosonicField,Lattice]
      module procedure calc_Pi_scGG                                             ![BosonicField,FermionicField,Lattice,tau_output(optional)]
   end interface calc_PiGG

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
   public :: calc_PiGG
   public :: calc_PiGGdc
   public :: calc_Pimp

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes analytically the non-interacting polarization bubble
   !COMMENT: For some reason this gives a bit more wiggling W.
   !         Using calc_Pi_scGkGk from the Glda seems to work better.
   !---------------------------------------------------------------------------!
   subroutine calc_Pi_GoGo(Pmats,Lttc)
      !
      use parameters
      use utils_misc
      use utils_fields
      use crystal
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: cprod(:,:,:,:,:)
      complex(8),allocatable                :: alpha(:)
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats,Norb
      integer                               :: ik1,ik2,iq,iw
      integer                               :: iwan1,iwan2,iwan3,iwan4,ib1,ib2
      logical                               :: clear
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Pi_GoGo"
      call cpu_time(start)
      !
      !
      ! Check on the input Fields
      if(.not.Pmats%status) stop "calc_Pi_GoGo: Pmats not properly initialized."
      if(Pmats%Nkpt.eq.0) stop "calc_Pi_GoGo: Pmats k dependent attributes not properly initialized."
      !
      Nbp = Pmats%Nbp
      Nkpt = Pmats%Nkpt
      Beta = Pmats%Beta
      Nmats = Pmats%Npoints
      Norb = int(sqrt(dble(Nbp)))
      if(Lttc%Norb.ne.Norb) stop "calc_Pi_GoGo: Pmats and Lattice have different orbital dimension."
      if(.not.allocated(Lttc%kptsum)) stop "calc_Pi_GoGo: kptsum not allocated."
      if(.not.allocated(Lttc%Zk)) stop "calc_Pi_GoGo: Zk not allocated."
      if(.not.allocated(Lttc%Ek)) stop "calc_Pi_GoGo: Ek not allocated."
      allocate(wmats(Nmats));wmats=0d0
      wmats = BosonicFreqMesh(Beta,Nmats)
      !
      !cprod(alpha,i,n,ik)= < B_q,alpha Psi_kn |Psi_q+k,i>
      allocate(cprod(Nbp,Norb,Norb,Nkpt,Nkpt));cprod=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Nmats,Norb,Lttc,cprod,beta),&
      !$OMP PRIVATE(iq,ik1,ik2,ib1,iwan1,iwan2,iwan3,iwan4)
      !$OMP DO
      do iq=1,Nkpt
         do ik1=1,Nkpt
            !
            ik2 = Lttc%kptsum(ik1,iq)
            !
            do iwan4=1,Norb
               do iwan3=1,Norb
                  ib1=0
                  do iwan2=1,Norb
                     do iwan1=1,Norb
                        !
                        ib1=ib1+1
                        cprod(ib1,iwan3,iwan4,ik1,iq) = dconjg(Lttc%Zk(iwan1,iwan4,ik1))*Lttc%Zk(iwan2,iwan3,ik2)
                        !
                     enddo
                  enddo
               enddo
            enddo
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      !
      allocate(alpha(Nmats));alpha=czero
      call clear_attributes(Pmats)
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iq,ik1,ik2,iwan1,iwan2,alpha)
      !$OMP DO
      do iq=1,Nkpt
         alpha=czero
         do iw=1,Nmats
            !
            do ik1=1,Nkpt
               ik2 = Lttc%kptsum(ik1,iq)
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     !
                     if (dabs(-Lttc%Ek(iwan1,ik1)+Lttc%Ek(iwan2,ik2)).lt.1d-6.and.iw.eq.1) then
                        !
                        !lim_{E'->E} (n(E)-n(E'))/(E'-E) = (n(E) - (n(E) + n'(E)*(E'-E)))/(E'-E) -> -n'(E)
                        alpha(iw) = +2.d0 * diff_fermidirac(Lttc%Ek(iwan1,ik1),Lttc%mu,Pmats%Beta) / Nkpt
                        !
                     else
                        !
                        alpha(iw) = -2.d0 * ( fermidirac(Lttc%Ek(iwan1,ik1),Lttc%mu,Pmats%Beta) - fermidirac(Lttc%Ek(iwan2,ik2),Lttc%mu,Pmats%Beta) ) &
                                          / ( dcmplx(0d0,1d0) * wmats(iw) - Lttc%Ek(iwan1,ik1) + Lttc%Ek(iwan2,ik2) ) / nkpt
                        !
                     endif
                     !
                     !alpha(iw)=alpha(iw)!/nkpt/nkpt !to account for the fact that cprod basis = 1/nkpt * prod basis of U in spex
                     call ZHER('U',Nbp,alpha(iw),cprod(:,iwan2,iwan1,ik1,iq),1,Pmats%screened(:,:,iw,iq),Nbp)
                     !
                  enddo !iwan2
               enddo !iwan1
            enddo !ik1
         enddo !iw
         if(verbose)print *, "     PiGG(q,iw) - done iq: ",iq
      enddo !iq
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(cprod,alpha,wmats)
      !
      !Clean up numerical noise
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iq,ib1,ib2,clear)
      !$OMP DO
      do iq=1,Nkpt
         do ib1=1,Nbp
            do ib2=ib1,Nbp
               clear = sum(abs(Pmats%screened(ib1,ib2,:,iq))) .lt. eps
               if(clear)then
                  Pmats%screened(ib1,ib2,:,iq)=czero
                  Pmats%screened(ib2,ib1,:,iq)=czero
               endif
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      call BosonicKsum(Pmats)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     PiGG cpu timing: ", finish-start
      !
   end subroutine calc_Pi_GoGo


   !---------------------------------------------------------------------------!
   !PURPOSE: Interface between the two possible ways to compute the bubble.
   !         Real-space implementation never used
   !---------------------------------------------------------------------------!
   subroutine calc_Pi_scGG(Pout,Gmats,Lttc,tau_output)
      !
      use parameters
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pout
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: tau_output
      logical                               :: tau_output_
      !
      logical                               :: cmplxHyb=.false.
      !
      tau_output_=.false.
      if(present(tau_output)) tau_output_ = tau_output
      !
      if(cmplxHyb)then
         call calc_Pi_scGrGr(Pout,Gmats,Lttc,tau_output=tau_output_)
      else
         call calc_Pi_scGkGk(Pout,Gmats,Lttc,tau_output=tau_output_)
      endif
      !
   end subroutine calc_Pi_scGG


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes polarization bubble from the interacting Gf in momentum
   !         space. if tau_output=T P(K,tau) is returned
   !---------------------------------------------------------------------------!
   subroutine calc_Pi_scGkGk(Pout,Gmats,Lttc,tau_output)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use fourier_transforms
      use crystal
      use fourier_transforms
      use file_io
      use input_vars, only : Ntau, tau_uniform, paramagnet
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pout
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: tau_output
      !
      complex(8),allocatable                :: Gitau(:,:,:,:,:)
      complex(8),allocatable                :: Pq_tau(:,:,:)
      real(8),allocatable                   :: tau(:)
      real(8)                               :: Beta,tau2
      integer                               :: Nbp,Nkpt,Norb,Ntau_,NaxisB
      integer                               :: ik1,ik2,iq,itau,ispin,ip
      integer                               :: i,j,k,l,ib1,ib2
      logical                               :: tau_output_,OD
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish,start_in,finish_in
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Pi_scGkGk"
      call cpu_time(start)
      !
      !
      ! Check on the input Fields
      if(.not.Pout%status) stop "calc_Pi_scGkGk: Pout not properly initialized."
      if(.not.Gmats%status) stop "calc_Pi_scGkGk: Green's function not properly initialized."
      if(Pout%Nkpt.eq.0) stop "calc_Pi_scGkGk: Pout k dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_Pi_scGkGk: Green's function k dependent attributes not properly initialized."
      if(Pout%Beta.ne.Gmats%Beta) stop "calc_Pi_scGkGk: Pout and Green's have different Beta."
      !
      Nbp = Pout%Nbp
      Nkpt = Pout%Nkpt
      Beta = Pout%Beta
      NaxisB = Pout%Npoints
      Norb = int(sqrt(dble(Nbp)))
      if(Gmats%Norb.ne.Norb) stop "calc_Pi_scGkGk: Pout and Green's function have different orbital dimension."
      if(Lttc%Norb.ne.Norb) stop "calc_Pi_scGkGk: Pout and Lttc have different orbital dimension."
      if(Gmats%Nkpt.ne.Nkpt) stop "calc_Pi_scGkGk: Pout and Green's function have different number of Kpoints."
      if(Lttc%Nkpt.ne.Nkpt) stop "calc_Pi_scGkGk: Pout and Lttc have different number of Kpoints."
      if(.not.allocated(Lttc%kptdif)) stop "calc_Pi_scGkGk: kptdif not allocated."
      !
      tau_output_=.false.
      if(present(tau_output)) tau_output_ = tau_output
      Ntau_ = Ntau
      if(tau_output_) Ntau_ = NaxisB
      !
      allocate(tau(Ntau_));tau=0d0
      if(tau_uniform)then
         tau = linspace(0d0,Beta,Ntau_)
      else
         tau = denspace(beta,Ntau_)
      endif
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      ! Compute Glat(k,tau)
      call cpu_time(start)
      allocate(Gitau(Norb,Norb,Ntau_,Nkpt,Nspin));Gitau=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin),asympt_corr=.true.,tau_uniform=tau_uniform)
        !if(.not.Gtau_K) call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin),asympt_corr=.true.,tau_uniform=tau_uniform,nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt)
         if(paramagnet)then
            Gitau(:,:,:,:,Nspin) = Gitau(:,:,:,:,1)
            exit
         endif
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(k,iw) --> Glat(k,tau) cpu timing:", finish-start
      !
      !Hermiticity check
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(ip,iq)
      !$OMP DO
      do ip=1,Ntau_
         do iq=1,Nkpt
            call check_Hermiticity(Gitau(:,:,ip,iq,1),eps,enforce=.false.,hardstop=.false.,name="Glat_t"//str(ip)//"_q"//str(iq)//"_s1",verb=.true.)
            if(.not.paramagnet)call check_Hermiticity(Gitau(:,:,ip,iq,Nspin),eps,enforce=.false.,hardstop=.false.,name="Glat_t"//str(ip)//"_q"//str(iq)//"_s2",verb=.true.)
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      call cpu_time(finish)
      write(*,"(A,F)") "     Hermiticity check on Glat(k,tau) cpu timing:", finish-start
      !
      !Compute the bubble in momentum space and tau axis
      allocate(Pq_tau(Nbp,Nbp,Ntau_))
      call clear_attributes(Pout)
      do iq=1,Nkpt
         !
         Pq_tau=czero
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(itau,tau2,ispin,ik1,ik2,i,j,k,l,ib1,ib2,OD)
         !$OMP DO
         do itau=1,Ntau_
            !
            tau2=tau(Ntau_)-tau(itau)
            if (dabs(tau2-tau(Ntau_-itau+1)).gt.eps) stop "calc_Pi_scGkGk: itau2 not found."
            !
            do ik1=1,Nkpt
               ik2=Lttc%kptdif(ik1,iq)
               !
               do ib1=1,Nbp
                  do ib2=ib1,Nbp
                     !
                     OD = ib1.ne.ib2
                     !
                     i = PhysicalUelements%Full_Map(ib1,ib2,1)
                     j = PhysicalUelements%Full_Map(ib1,ib2,2)
                     k = PhysicalUelements%Full_Map(ib1,ib2,3)
                     l = PhysicalUelements%Full_Map(ib1,ib2,4)
                     !
                     do ispin=1,Nspin
                        !
                        !Pi_(i,j)(k,l)(q,tau) = - sum_k G_ik(k,tau) * G_lj(q-k,beta-tau)
                        Pq_tau(ib1,ib2,itau) = Pq_tau(ib1,ib2,itau) - ( Gitau(i,k,itau,ik1,ispin) * Gitau(l,j,Ntau_-itau+1,ik2,ispin) )/Nkpt
                        !
                        !Pi_(k,l)(i,j)(q,tau) = - sum_k G_ki(k,tau) * G_jl(q-k,beta-tau)
                        if(OD)Pq_tau(ib2,ib1,itau) = Pq_tau(ib2,ib1,itau) - ( Gitau(k,i,itau,ik1,ispin) * Gitau(j,l,Ntau_-itau+1,ik2,ispin) )/Nkpt
                        !
                     enddo
                     !
                  enddo
               enddo
               !
            enddo !ik1
            !
         enddo !itau
         !$OMP END DO
         !$OMP END PARALLEL
         !
         !FT to matsubara all the q points
         if(tau_output_)then
            Pout%screened(:,:,:,iq) = Pq_tau
         else
            call cpu_time(start_in)
            call Bitau2mats(Beta,Pq_tau,Pout%screened(:,:,:,iq),tau_uniform=tau_uniform)
            call cpu_time(finish_in)
            if(verbose)write(*,"(A,F)") "     PiGGsc(q_"//str(iq)//",tau) --> PiGGsc(q_"//str(iq)//",iw) cpu timing:", finish_in-start_in
         endif
         !
      enddo !iq
      deallocate(tau,Gitau,Pq_tau)
      !
      ! Fill the local attributes
      call BosonicKsum(Pout)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     PiGGsc[K] cpu timing: ", finish-start
      !
      !call dump_BosonicField(Pout,"./Plat_readable/",.false.)
      !
   end subroutine calc_Pi_scGkGk


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes polarization bubble from the interacting Gf in real space
   !         if tau_output=T P(R,tau) is returned
   !         Not used, implemented just for testing
   !---------------------------------------------------------------------------!
   subroutine calc_Pi_scGrGr(Pout,Gmats,Lttc,tau_output)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use fourier_transforms
      use crystal
      use fourier_transforms
      use file_io
      use input_vars, only : Ntau, tau_uniform, paramagnet
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pout
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: tau_output
      !
      complex(8),allocatable                :: Gitau(:,:,:,:,:),Gitau_k(:,:,:,:)
      complex(8),allocatable                :: Pr_itau(:,:,:,:),Pr_mats(:,:,:,:)
      real(8),allocatable                   :: tau(:)
      real(8)                               :: Beta,tau2
      integer                               :: Nbp,Nkpt,Norb,Ntau_,NaxisB
      integer                               :: iw,itau,ispin,ir,ir2
      integer                               :: i,j,k,l,ib1,ib2
      logical                               :: tau_output_,OD,Rcond
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish,start_in,finish_in
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Pi_scGrGr"
      call cpu_time(start)
      !
      !
      ! Check on the input Fields
      if(.not.Pout%status) stop "calc_Pi_scGrGr: Pout not properly initialized."
      if(.not.Gmats%status) stop "calc_Pi_scGrGr: Green's function not properly initialized."
      if(Pout%Nkpt.eq.0) stop "calc_Pi_scGrGr: Pout k dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_Pi_scGrGr: Green's function k dependent attributes not properly initialized."
      if(Pout%Beta.ne.Gmats%Beta) stop "calc_Pi_scGrGr: Pout and Green's have different Beta."
      !
      Nbp = Pout%Nbp
      Nkpt = Pout%Nkpt
      Beta = Pout%Beta
      NaxisB = Pout%Npoints
      Norb = int(sqrt(dble(Nbp)))
      if(Gmats%Norb.ne.Norb) stop "calc_Pi_scGrGr: Pout and Green's function have different orbital dimension."
      if(Lttc%Norb.ne.Norb) stop "calc_Pi_scGrGr: Pout and Lttc have different orbital dimension."
      if(Gmats%Nkpt.ne.Nkpt) stop "calc_Pi_scGrGr: Pout and Green's function have different number of Kpoints."
      if(Lttc%Nkpt.ne.Nkpt) stop "calc_Pi_scGrGr: Pout and Lttc have different number of Kpoints."
      if(.not.allocated(Lttc%kptdif)) stop "calc_Pi_scGrGr: kptdif not allocated."
      !
      if(.not.Wig_stored) call calc_wignerseiz(Lttc%Nkpt3)
      !
      tau_output_=.false.
      if(present(tau_output)) tau_output_ = tau_output
      Ntau_ = Ntau
      if(tau_output_)then
         Ntau_ = NaxisB
         if(NaxisB.ne.Nwig) stop  "calc_Pi_scGrGr: P(R,tau) requested but Npoints attribute does not match with the number of real space vectors."
      endif
      !
      allocate(tau(Ntau_));tau=0d0
      if(tau_uniform)then
         tau = linspace(0d0,Beta,Ntau_)
      else
         tau = denspace(beta,Ntau_)
      endif
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      !FT to tau axis and then FT to real space
      call cpu_time(start)
      allocate(Gitau(Norb,Norb,Ntau_,Nwig,Nspin));Gitau=czero
      do ispin=1,Nspin
         !
         allocate(Gitau_k(Norb,Norb,Ntau_,Nkpt));Gitau_k=czero
         call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau_k,asympt_corr=.true.,tau_uniform=tau_uniform)
         !
         call wannier_K2R(Lttc%Nkpt3,Lttc%kpt,Gitau_k,Gitau(:,:,:,:,ispin))
         deallocate(Gitau_k)
         !
         if(paramagnet)then
            Gitau(:,:,:,:,Nspin) = Gitau(:,:,:,:,1)
            exit
         endif
         !
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(k,iw) --> Glat(R,tau) cpu timing:", finish-start
      !
      !Compute the bubble in real space and tau axis
      call cpu_time(start)
      allocate(Pr_itau(Nbp,Nbp,Ntau_,Nwig));Pr_itau=czero
      if(.not.tau_output_) allocate(Pr_mats(Nbp,Nbp,NaxisB,Nwig));Pr_mats=czero
      loopR: do ir=1,Nwig
         !
         Rcond = (Nvecwig(3,ir).lt.0) .or. ((Nvecwig(3,ir).eq.0).and.(Nvecwig(1,ir).lt.0))
         !
         ir2=0
         if(all(Nvecwig(:,ir).eq.[0,0,0]))then
            !
            !R=0 --> R = -R
            ir2 = ir
            if(ir.ne.wig0) stop "calc_Pi_scGrGr: ir.ne.wig0 something is wrong in the real-space vectors indexing."
            !
         elseif(Rcond)then
            !
            !exploiting: Pi_(ab)(cd)(R,tau) = Pi*_(cd)(ab)(-R,tau)
            cycle loopR
            !
         else
            !
            !find index of -R if not found returns 0
            ir2 = find_vec(-Nvecwig(:,ir),Nvecwig,hardstop=.false.)
            !
         endif
         !
         if(verbose)write(*,"(2(A,1I5),A,3I4,A,1I5,A,3I4,A,1I4)") "     Nwig: ",Nwig," iwig: ",ir," -> ",Nvecwig(:,ir),"  -  ",ir2," : ",-Nvecwig(:,ir)," deg:",nrdegwig(ir)
         !
         !index of -R if not found
         if(ir2.eq.0)cycle loopR
         !
         call cpu_time(start_in)
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(itau,tau2,ispin,i,j,k,l,ib1,ib2,OD)
         !$OMP DO
         do itau=1,Ntau_
            !
            tau2=tau(Ntau_)-tau(itau)
            if (dabs(tau2-tau(Ntau_-itau+1)).gt.eps) stop "calc_Pi_scGrGr: itau2 not found."
            !
            do ib1=1,Nbp
               do ib2=ib1,Nbp
                  !
                  OD = ib1.ne.ib2
                  !
                  i = PhysicalUelements%Full_Map(ib1,ib2,1)
                  j = PhysicalUelements%Full_Map(ib1,ib2,2)
                  k = PhysicalUelements%Full_Map(ib1,ib2,3)
                  l = PhysicalUelements%Full_Map(ib1,ib2,4)
                  !
                  !Filling the R>=0 directions
                  do ispin=1,Nspin
                     !
                     !Pi_(i,j)(k,l)(R,tau) = - G_ik(R,tau) * G_lj(-R,beta-tau)
                     Pr_itau(ib1,ib2,itau,ir) = Pr_itau(ib1,ib2,itau,ir) - ( Gitau(i,k,itau,ir,ispin) * Gitau(l,j,Ntau_-itau+1,ir2,ispin) )!*nrdegwig(ir)
                     !
                     !Pi_(k,l)(i,j)(R,tau) = - G_ki(R,tau) * G_jl(-R,beta-tau)
                     if(OD)Pr_itau(ib2,ib1,itau,ir) = Pr_itau(ib1,ib2,itau,ir) - ( Gitau(k,i,itau,ir,ispin) * Gitau(j,l,Ntau_-itau+1,ir2,ispin) )!*nrdegwig(ir)
                     !
                  enddo
                  !
                  !Filling the R<0 directions
                  if(ir2.ne.ir)then
                     !
                     !Pi_(i,j)(k,l)(-R,tau) = Pi*_(k,l)(i,j)(R,tau)
                     Pr_itau(ib1,ib2,itau,ir2) = conjg(Pr_itau(ib2,ib1,itau,ir))
                     !
                     !Pi_(k,l)(i,j)(-R,tau) = Pi*_(i,j)(k,l)(R,tau)
                     Pr_itau(ib2,ib1,itau,ir2) = conjg(Pr_itau(ib1,ib2,itau,ir))
                     !
                  endif
                  !
               enddo
            enddo
            !
         enddo !itau
         !$OMP END DO
         !$OMP END PARALLEL
         call cpu_time(finish_in)
         write(*,"(A,F)") new_line("A")//"     PiGGsc(R_"//str(ir)//",tau) = Glat(R_"//str(ir)//",tau)*Glat(-R_"//str(ir)//",-tau) cpu timing:", finish_in-start_in
         !
         !FT to matsubara each R vector
         if(tau_output_)then
            Pout%screened(:,:,:,ir) = Pr_itau(:,:,:,ir)
         else
            call cpu_time(start_in)
            call Bitau2mats(Beta,Pr_itau(:,:,:,ir),Pr_mats(:,:,:,ir),tau_uniform=tau_uniform)
            call cpu_time(finish_in)
            write(*,"(A,F)") "     PiGGsc(R_"//str(ir)//",tau) --> PiGGsc(R_"//str(ir)//",iw) cpu timing:", finish_in-start_in
         endif
         !
      enddo loopR
      deallocate(Pr_itau)
      !
      !FT to momentum space
      if(.not.tau_output_)then
         call cpu_time(start_in)
         do iw=1,NaxisB
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Pr_mats(:,:,iw,:),Pout%screened(:,:,iw,:))
         enddo
         deallocate(Pr_mats)
         call cpu_time(finish_in)
         write(*,"(A,F)") "     PiGGsc(R,iw) --> PiGGsc(K,iw) cpu timing:", finish_in-start_in
      endif
      !
      ! Fill the local attributes
      call BosonicKsum(Pout)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     PiGGsc[R] cpu timing: ", finish-start
      !
   end subroutine calc_Pi_scGrGr


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes polarization bubble double counting as Gloc*Gloc
   !---------------------------------------------------------------------------!
   subroutine calc_PiGGdc(Pout,Gmats)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use fourier_transforms
      use crystal
      use fourier_transforms
      use file_io
      use input_vars, only : Ntau, tau_uniform, paramagnet
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pout
      type(FermionicField),intent(in)       :: Gmats
      !
      complex(8),allocatable                :: Gitau(:,:,:,:)
      complex(8),allocatable                :: Ptau_dc(:,:,:)
      real(8),allocatable                   :: tau(:)
      real(8)                               :: Beta,tau2
      integer                               :: Nbp,Norb,itau,ispin,ip
      integer                               :: i,j,k,l,ib1,ib2
      logical                               :: OD
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_PiGGdc"
      call cpu_time(start)
      !
      !
      ! Check on the input Fields
      if(.not.Pout%status) stop "calc_PiGGdc: Pout not properly initialized."
      if(.not.Gmats%status) stop "calc_PiGGdc: Green's function not properly initialized."
      if(Pout%Nkpt.ne.0) stop "calc_PiGGdc: Pout k dependent attributes are supposed to be unallocated."
      if(Pout%Beta.ne.Gmats%Beta) stop "calc_PiGGdc: Pout and Green's have different Beta."
      !
      Nbp = Pout%Nbp
      Beta = Pout%Beta
      Norb = int(sqrt(dble(Nbp)))
      if(Gmats%Norb.ne.Norb) stop "calc_PiGGdc: Pout and Green's function have different orbital dimension."
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform)then
         tau = linspace(0d0,Beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      ! Compute Glat(tau) - FT all components
      call cpu_time(start)
      allocate(Gitau(Norb,Norb,Ntau,Nspin));Gitau=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%ws(:,:,:,ispin),Gitau(:,:,:,ispin),asympt_corr=.true.,tau_uniform=tau_uniform)
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(iw), --> Glat(tau) cpu timing:", finish-start
      !
      !Hermiticity check
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(ip)
      !$OMP DO
      do ip=1,Ntau
         call check_Hermiticity(Gitau(:,:,ip,1),eps,enforce=.false.,hardstop=.false.,name="Glat_t"//str(ip)//"_s1",verb=.true.)
         if(.not.paramagnet)call check_Hermiticity(Gitau(:,:,ip,Nspin),eps,enforce=.false.,hardstop=.false.,name="Glat_t"//str(ip)//"_s2",verb=.true.)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      call cpu_time(finish)
      write(*,"(A,F)") "     Hermiticity check on Glat(tau) cpu timing:", finish-start
      !
      !Compute the bubble double counting and tau axis
      allocate(Ptau_dc(Nbp,Nbp,Ntau))
      call clear_attributes(Pout)
      !
      Ptau_dc=czero
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(itau,tau2,ispin,i,j,k,l,ib1,ib2,OD)
      !$OMP DO
      do itau=1,Ntau
         !
         tau2=tau(Ntau)-tau(itau)
         if (dabs(tau2-tau(Ntau-itau+1)).gt.eps) stop "calc_PiGGdc: itau2 not found."
         !
         do ib1=1,Nbp
            do ib2=ib1,Nbp
               !
               OD = ib1.ne.ib2
               !
               i = PhysicalUelements%Full_Map(ib1,ib2,1)
               j = PhysicalUelements%Full_Map(ib1,ib2,2)
               k = PhysicalUelements%Full_Map(ib1,ib2,3)
               l = PhysicalUelements%Full_Map(ib1,ib2,4)
               !
               do ispin=1,Nspin
                  !
                  !Pi_(i,j)(k,l)(q,tau) = - sum_k G_ik(k,tau) * G_lj(q-k,beta-tau)
                  Ptau_dc(ib1,ib2,itau) = Ptau_dc(ib1,ib2,itau) - ( Gitau(i,k,itau,ispin) * Gitau(l,j,Ntau-itau+1,ispin) )
                  !
                  !Pi_(k,l)(i,j)(q,tau) = - sum_k G_ki(k,tau) * G_jl(q-k,beta-tau)
                  if(OD)Ptau_dc(ib2,ib1,itau) = Ptau_dc(ib2,ib1,itau) - ( Gitau(k,i,itau,ispin) * Gitau(j,l,Ntau-itau+1,ispin) )
                  !
               enddo
               !
            enddo
         enddo
         !
      enddo !itau
      !$OMP END DO
      !$OMP END PARALLEL
      !
      call Bitau2mats(Beta,Ptau_dc,Pout%screened_local,tau_uniform=tau_uniform)
      !
      deallocate(tau,Gitau,Ptau_dc)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     PiGGdc cpu timing: ", finish-start
      !
   end subroutine calc_PiGGdc


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the local polarization vertex
   !---------------------------------------------------------------------------!
   subroutine calc_Pimp(Pimp,curlyU,ChiC,sym,NaNb)
      !
      use parameters
      use utils_fields
      use utils_misc
      use linalg, only : zeye, inv
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pimp
      type(BosonicField),intent(in)         :: curlyU
      type(BosonicField),intent(in)         :: ChiC
      logical,intent(in),optional           :: sym
      logical,intent(in),optional           :: NaNb
      !
      real(8),allocatable                   :: RecurlyU(:,:),ReChiC(:,:)
      complex(8),allocatable                :: invP(:,:)
      type(physicalU)                       :: PhysicalUelements
      real(8)                               :: Beta
      integer                               :: Nbp,Nmats,ib1,ib2
      integer                               :: iw
      logical                               :: sym_,NaNb_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Pimp"
      !
      !
      ! Check on the input Fields
      if(.not.Pimp%status) stop "calc_Pimp: Pimp not properly initialized."
      if(.not.curlyU%status) stop "calc_Pimp: curlyU not properly initialized."
      if(.not.ChiC%status) stop "calc_Pimp: ChiC not properly initialized."
      if(Pimp%Nkpt.ne.0) stop "calc_Pimp: Pimp k dependent attributes are supposed to be unallocated."
      if(curlyU%Nkpt.ne.0) stop "calc_Pimp: curlyU k dependent attributes are supposed to be unallocated."
      if(ChiC%Nkpt.ne.0) stop "calc_Pimp: ChiC k dependent attributes are supposed to be unallocated."
      if(allocated(Pimp%bare_local))  stop "calc_Pimp: Pimp bare_local attribute is supposed to be unallocated."
      if(allocated(Pimp%bare))  stop "calc_Pimp: Pimp bare attribute is supposed to be unallocated."
      if(allocated(ChiC%bare_local))  stop "calc_Pimp: ChiC bare_local attribute is supposed to be unallocated."
      if(allocated(ChiC%bare))  stop "calc_Pimp: ChiC bare attribute is supposed to be unallocated."
      !
      sym_=.true.
      if(present(sym))sym_=sym
      NaNb_=.true.
      if(present(NaNb))NaNb_=NaNb
      !
      Nbp = Pimp%Nbp
      Beta = Pimp%Beta
      Nmats = Pimp%Npoints
      !
      if(all([curlyU%Nbp-Nbp,ChiC%Nbp-Nbp].ne.[0,0])) stop "calc_Pimp: Either curlyU and/or ChiC have different orbital dimension with respect to Pimp."
      if(all([curlyU%Beta-Beta,ChiC%Beta-Beta].ne.[0d0,0d0])) stop "calc_Pimp: Either curlyU and/or ChiC have different Beta with respect to Pimp."
      if(all([curlyU%Npoints-Nmats,ChiC%Npoints-Nmats].ne.[0,0])) stop "calc_Pimp: Either curlyU and/or ChiC have different number of Matsubara points with respect to Pimp."
      !
      if(NaNb_)call init_Uelements(int(sqrt(dble(Nbp))),PhysicalUelements)
      !
      call clear_attributes(Pimp)
      !
      allocate(invP(Nbp,Nbp));invP=czero
      allocate(RecurlyU(Nbp,Nbp));RecurlyU=0d0
      allocate(ReChiC(Nbp,Nbp));ReChiC=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iw,invP,RecurlyU,ReChiC,ib1,ib2)
      !$OMP DO
      do iw=1,Pimp%Npoints
         !
         ReChiC = dreal(ChiC%screened_local(:,:,iw))
         RecurlyU = dreal(curlyU%screened_local(:,:,iw))
         !
         ! keeping only curlyU_(aa)(bb) to compute Pimp
         if(NaNb_)then
            do ib1=1,Pimp%Nbp
               do ib2=ib1,Pimp%Nbp
                  if(.not.PhysicalUelements%Full_Imp(ib1,ib2))then
                     RecurlyU(ib1,ib2) = 0d0
                     RecurlyU(ib2,ib1) = 0d0
                  endif
               enddo
            enddo
         endif
         !
         ! [ curlyU*ChiC - 1 ]
         invP = matmul(RecurlyU,ReChiC) - zeye(Pimp%Nbp)
         !
         ! [ curlyU*ChiC - 1 ]^-1
         call inv(invP)
         !
         ! ChiC*[ curlyU*ChiC - 1 ]^-1
         Pimp%screened_local(:,:,iw) = matmul(ReChiC,invP)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(invP,RecurlyU,ReChiC)
      call isReal(Pimp)
      !
      !Check if Pimp is locally symmetric - print if relative error is bigger than 1e-3
      if(sym_)then
         write(*,"(A)") "     Checking symmetry of Pimp (enforced)."
         do iw=1,Nmats
            call check_Symmetry(Pimp%screened_local(:,:,iw),1e7*eps,enforce=.true.,hardstop=.false.,name="Pimp_w"//str(iw),verb=.true.)
         enddo
      endif
      !
   end subroutine calc_Pimp


end module bubbles



!
!old implementation
!do l=1,Norb
!   do k=1,Norb
!      do j=1,Norb
!         do i=1,Norb
!            !
!            !(i,j)(k,l). Second index varying faster: (1,1),(1,2),(1,3),...
!            ib1 = j + Norb*(i-1)
!            ib2 = l + Norb*(k-1)
!            !
!            do ispin=1,Nspin
!               Pq_tau(ib1,ib2,itau) = Pq_tau(ib1,ib2,itau) - ( Gitau(i,k,itau,ik1,ispin) * Gitau(l,j,Ntau_-itau+1,ik2,ispin) )/Nkpt
!            enddo
!            !
!         enddo
!      enddo
!   enddo
!enddo
