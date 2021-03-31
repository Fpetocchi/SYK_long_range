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
   interface calc_Pi
      module procedure calc_Pi_GoGo                                             ![BosonicField,Lattice]
      module procedure calc_Pi_scGG                                             ![BosonicField,FermionicField,Lattice,tau_output(optional)]
   end interface calc_Pi

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
   public :: calc_Pi
   public :: calc_Pimp

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes analytically the non-interacting polarization bubble
   !TEST ON: 21-10-2020
   !COMMENT: For some reason this gives a bit more wiggling W.
   !         Using calc_Pi_scGG from the Glda seems to work better.
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
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Norb,wmats,cprod,Lttc,Pmats,verbose),&
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
            !
            do ib2=1,Nbp
               do ib1=ib2+1,Nbp
                  if(abs(Pmats%screened(ib2,ib1,iw,iq)).lt.eps)Pmats%screened(ib2,ib1,iw,iq)=czero
                  Pmats%screened(ib1,ib2,iw,iq)=conjg(Pmats%screened(ib2,ib1,iw,iq))
               enddo
            enddo
            if(verbose)call check_Hermiticity(Pmats%screened(:,:,iw,iq),eps,hardstop=.true.)
            !
         enddo !iw
         if(verbose)print *, "     PiGG(q,iw) - done iq: ",iq
      enddo !iq
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(cprod,alpha,wmats)
      !
      call BosonicKsum(Pmats)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     PiGG cpu timing: ", finish-start
      !
   end subroutine calc_Pi_GoGo


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes polarization bubble from the interacting Gf
   !TEST ON: 21-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_Pi_scGG(Pout,Gmats,Lttc,tau_output,sym)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use fourier_transforms
      use crystal
      use fourier_transforms
      use file_io
      use interactions, only: init_Uelements
      use input_vars, only : NtauB, tau_uniform, cmplxWann, UfullStructure, paramagnet
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pout
      type(FermionicField),intent(in)       :: Gmats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: tau_output
      logical,intent(in),optional           :: sym
      !
      type(physicalU)                       :: PhysicalUelements
      complex(8),allocatable                :: Gitau(:,:,:,:,:)
      complex(8),allocatable                :: Pq_tau(:,:,:)
      real(8),allocatable                   :: tau(:)
      real(8)                               :: Beta,tau2
      integer                               :: Nbp,Nkpt,Norb,Ntau_,NaxisB
      integer                               :: ik1,ik2,iq,itau,ispin,ip
      integer                               :: i,j,k,l,ib1,ib2
      logical                               :: tau_output_,sym_
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Pi_scGG"
      call cpu_time(start)
      !
      !
      ! Check on the input Fields
      if(.not.Pout%status) stop "calc_Pi_scGG: Pout not properly initialized."
      if(.not.Gmats%status) stop "calc_Pi_scGG: Green's function not properly initialized."
      if(Pout%Nkpt.eq.0) stop "calc_Pi_scGG: Pout k dependent attributes not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "calc_Pi_scGG: Green's function k dependent attributes not properly initialized."
      if(Pout%Beta.ne.Gmats%Beta) stop "calc_Pi_scGG: Pout and Green's have different Beta."
      !
      Nbp = Pout%Nbp
      Nkpt = Pout%Nkpt
      Beta = Pout%Beta
      NaxisB = Pout%Npoints
      Norb = int(sqrt(dble(Nbp)))
      if(Gmats%Norb.ne.Norb) stop "calc_Pi_scGG: Pout and Green's function have different orbital dimension."
      if(Gmats%Nkpt.ne.Nkpt) stop "calc_Pi_scGG: Pout and Green's function have different number of Kpoints."
      if(.not.allocated(Lttc%kptdif)) stop "calc_Pi_scGG: kptdif not allocated."
      !
      sym_=.true.
      if(present(sym))sym_=sym
      !
      tau_output_=.false.
      if(present(tau_output)) tau_output_ = tau_output
      Ntau_ = NtauB
      if(tau_output_) Ntau_ = NaxisB
      !
      call init_Uelements(Norb,PhysicalUelements)
      !
      allocate(tau(Ntau_));tau=0d0
      if(tau_uniform)then
         tau = linspace(0d0,Beta,Ntau_)
      else
         tau = denspace(beta,Ntau_)
      endif
      !
      ! Compute Glat(k,tau)
      call cpu_time(start)
      allocate(Gitau(Norb,Norb,Ntau_,Nkpt,Nspin));Gitau=czero
      if(cmplxWann)then
         spinloopC: do ispin=1,Nspin
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin), &
            asympt_corr=.true.,tau_uniform=tau_uniform)
            if(paramagnet)then
               Gitau(:,:,:,:,Nspin) = Gitau(:,:,:,:,1)
               exit spinloopC
            endif
         enddo spinloopC
      else
         spinloopR: do ispin=1,Nspin
            call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin), &
            asympt_corr=.true.,tau_uniform=tau_uniform,nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt)
            if(paramagnet)then
               Gitau(:,:,:,:,Nspin) = Gitau(:,:,:,:,1)
               exit spinloopR
            endif
         enddo spinloopR
      endif
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(k,iw) --> Glat(k,tau) cpu timing:", finish-start
      !
      !Hermiticity check
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Ntau_,Nkpt,Gitau,paramagnet),&
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
      ! Compute the bubble
      allocate(Pq_tau(Nbp,Nbp,Ntau_))
      call clear_attributes(Pout)
      do iq=1,Nkpt
         !
         Pq_tau=czero
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(iq,Ntau_,Nkpt,Norb,tau,Lttc,Gitau,Pq_tau,UfullStructure,PhysicalUelements),&
         !$OMP PRIVATE(itau,tau2,ispin,ik1,ik2,i,j,k,l,ib1,ib2)
         !$OMP DO
         do itau=1,Ntau_
            !
            tau2=tau(Ntau_)-tau(itau)
            if (dabs(tau2-tau(Ntau_-itau+1)).gt.eps) stop "calc_Pi_scGG: itau2"
            !
            do ik1=1,Nkpt
               ik2=Lttc%kptdif(ik1,iq)
               !
               do j=1,Norb
                  do l=1,Norb
                     do k=1,Norb
                        do i=1,Norb
                           !
                           ib1 = i + Norb*(j-1)
                           ib2 = k + Norb*(l-1)
                           !
                           !Avoid the unwanted bubble components
                           if((.not.UfullStructure).and.(.not.PhysicalUelements%Full_All(ib1,ib2)))cycle
                           !
                           do ispin=1,Nspin
                              Pq_tau(ib1,ib2,itau) = Pq_tau(ib1,ib2,itau) - ( Gitau(i,k,itau,ik1,ispin) * Gitau(l,j,Ntau_-itau+1,ik2,ispin) )/Nkpt
                           enddo
                           !
                        enddo
                     enddo
                  enddo
               enddo
               !
            enddo !ik1
            !
         enddo !itau
         !$OMP END DO
         !$OMP END PARALLEL
         !
         !FT to matsubara
         if(tau_output_)then
            Pout%screened(:,:,:,iq) = Pq_tau
         else
            call Bitau2mats(Beta,Pq_tau,Pout%screened(:,:,:,iq),tau_uniform=tau_uniform)
         endif
         !
         !Hermiticity check - generally it is not
         !if(sym_)then
         !   do ip=1,NaxisB
         !      call check_Hermiticity(Pout%screened(:,:,ip,iq),eps,enforce=.false.,hardstop=.false.,name="Plat_w"//str(ip)//"_q"//str(iq),verb=(mod(iq,251).eq.0))
         !   enddo
         !endif
         !
         if(verbose) write(*,"(A,I)") "     PiGGsc(q,iw) - done iq: ",iq
         !
      enddo !iq
      deallocate(tau,Gitau,Pq_tau)
      !
      ! Fill the local attributes
      call BosonicKsum(Pout)
      !
      !Hermiticity check - generally it is not
      !if(sym_)then
      !   write(*,"(A)") "     Checking hermiticity of local Plat."
      !   do ip=1,NaxisB
      !      call check_Hermiticity(Pout%screened_local(:,:,ip),eps,enforce=.false.,hardstop=.false.,name="Plat_loc_w"//str(ip))
      !   enddo
      !endif
      !
      call cpu_time(finish)
      write(*,"(A,F)") "     PiGGsc cpu timing: ", finish-start
      !call dump_BosonicField(Pout,"./Plat_readable/",.false.)
      !
   end subroutine calc_Pi_scGG


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the local polarization vertex
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine calc_Pimp(Pimp,curlyU,ChiC,sym)
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
      !
      complex(8),allocatable                :: invP(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nmats
      integer                               :: iw
      logical                               :: sym_
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
      !
      Nbp = Pimp%Nbp
      Beta = Pimp%Beta
      Nmats = Pimp%Npoints
      !
      if(all([curlyU%Nbp-Nbp,ChiC%Nbp-Nbp].ne.[0,0])) stop "calc_Pimp: Either curlyU and/or ChiC have different orbital dimension with respect to Pimp."
      if(all([curlyU%Beta-Beta,ChiC%Beta-Beta].ne.[0d0,0d0])) stop "calc_Pimp: Either curlyU and/or ChiC have different Beta with respect to Pimp."
      if(all([curlyU%Npoints-Nmats,ChiC%Npoints-Nmats].ne.[0,0]))   stop "calc_Pimp: Either curlyU and/or ChiC have different number of Matsubara points with respect to Pimp."
      !
      call clear_attributes(Pimp)
      !
      allocate(invP(Nbp,Nbp));invP=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Pimp,ChiC,curlyU),&
      !$OMP PRIVATE(iw,invP)
      !$OMP DO
      do iw=1,Pimp%Npoints
         !
         ! [ curlyU*ChiC - 1 ]
         invP = matmul(curlyU%screened_local(:,:,iw),ChiC%screened_local(:,:,iw)) - zeye(Pimp%Nbp)
         !
         ! [ curlyU*ChiC - 1 ]^-1
         call inv(invP)
         !
         ! ChiC*[ curlyU*ChiC - 1 ]^-1
         Pimp%screened_local(:,:,iw) = matmul(ChiC%screened_local(:,:,iw),invP)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(invP)
      call isReal(Pimp)
      !
      if(sym_)then
         write(*,"(A)") "     Checking symmetry of Pimp (enforced)"
         do iw=1,Nmats
            call check_Symmetry(Pimp%screened_local(:,:,iw),eps,enforce=.true.,hardstop=.false.,name="Pimp_w"//str(iw))
         enddo
      endif
      !
   end subroutine calc_Pimp


end module bubbles
