module bubbles

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   interface calc_Pi
      module procedure calc_Pi_GkGk                                             ![BosonicField,LocalOnly,save2bin,pathOUTPUT(optional to change output path),doAC(optional to override AC)]
      module procedure calc_Pi_selfcons                                         ![Matrix,pathOUTPUT(optional to change output path)]
   end interface calc_Pi

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: calc_Pi
   !public :: calc_Optcond
   !public :: calc_Hall

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes analytically the non-interacting polarization bubble
   !---------------------------------------------------------------------------!
   subroutine calc_Pi_GkGk(Pmats,Lttc)
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
      complex(8),allocatable                :: cprod(:,:,:,:)
      complex(8),allocatable                :: alpha(:)
      real(8),allocatable                   :: wmats(:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats,Norb
      integer                               :: ik1,ik2,iq,iw
      integer                               :: iwan1,iwan2,ib1,ib2
      !
      !
      write(*,*) "--- calc_Pi_GkGk ---"
      !
      !
      ! Check on the input Bosons
      if(.not.Pmats%status) stop "Pmats not properly initialized."
      if(Pmats%Nkpt.eq.0) stop "Pmats k dependent attributes not properly initialized."
      !
      Nbp = Pmats%Nbp
      Nkpt = Pmats%Nkpt
      Beta = Pmats%Beta
      Nmats = Pmats%Npoints
      Norb = int(sqrt(dble(Nbp)))
      if(Lttc%Norb.ne.Norb) stop "Pmats and Lattice have different orbital dimension."
      if(.not.allocated(Lttc%kptsum)) stop "kptsum not allocated."
      if(.not.allocated(Lttc%Zk)) stop "Zk not allocated."
      if(.not.allocated(Lttc%Ek)) stop "Ek not allocated."
      allocate(wmats(Nmats));wmats=0d0
      wmats = BosonicFreqMesh(Beta,Nmats)
      !
      !cprod(alpha,i,n,ik)= < B_q,alpha Psi_kn |Psi_q+k,i>
      allocate(cprod(Nbp,Norb,Norb,Nkpt));cprod=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Nmats,Norb,Lttc,cprod),&
      !$OMP PRIVATE(iq,ik1,ik2,ib1,ib2,ipb1,iwan1,iwan2)
      !$OMP DO
      do iq=1,Nkpt
         do ik1=1,Nkpt
            ik2 = Lttc%kptsum(ik1,iq)
            !
            do ib1=1,Norb
               do ib2=1,Norb
                  ipb1=0
                  do iwan2=1,Norb
                     do iwan1=1,Norb
                        !
                        ipb1=ipb1+1
                        cprod(ipb1,ib2,ib1,ik1) = dconjg(Lttc%Zk(iwan1,ib1,ik1))*Lttc%Zk(iwan2,ib2,ik2)
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
      Pmats%bare=czero; Pmats%bare_local=czero
      Pmats%screened=czero; Pmats%screened_local=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Norb,wmats,cprod,alpha,Lttc,Pmats),&
      !$OMP PRIVATE(iw,iq,ik1,ik2,iwan1,iwan2)
      !$OMP DO
      do iw=1,Nmats
         do iq=1,Nkpt
            do ik1=1,Nkpt
               ik2 = Lttc%kptsum(ik1,iq)
               !
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     !
                     if (dabs(-Lttc%Ek(iwan1,ik1)+Lttc%Ek(iwan2,ik2)).lt.eps.and.iw.eq.1) then
                        !
                        !lim_{E'->E} (n(E)-n(E'))/(E'-E) = (n(E) - (n(E) + n'(E)*(E'-E)))/(E'-E) -> -n'(E)
                        alpha(iw) = +2.d0 * diff_fermidirac(Lttc%Ek(iwan1,ik1),Lttc%mu) / Nkpt
                        !
                     else
                        !
                        alpha(iw) = -2.d0 * ( fermidirac(Lttc%Ek(iwan1,ik1),Lttc%mu) - fermidirac(Lttc%Ek(iwan2,ik2),Lttc%mu) ) &
                                          / ( dcmplx(0d0,1d0) * wmats(iw) - Lttc%Ek(iwan1,ik1) + Lttc%Ek(iwan2,ik2) ) / nkpt
                        !
                     endif
                     !
                     !alpha(iw)=alpha(iw)!/nkpt/nkpt !to account for the fact that cprod basis = 1/nkpt * prod basis of U in spex
                     call ZHER('U',Nbp,alpha(iw),cprod(:,iwan2,iwan1,ik1),1,Pmats%screened(:,:,iw,iq),Nbp)
                     !
                  enddo
               enddo
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(cprod,alpha,wmats)
      !
      call BosonicKsum(Pmats)
      !
   end subroutine calc_Pi_GkGk


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes polarization bubble from the interacting Gf
   !---------------------------------------------------------------------------!
   subroutine calc_Pi_selfcons(Pout,Gin,Ntau,tau_uniform,tau_output)
      !
      use parameters
      use utils_misc
      use utils_fields
      use fourier_transforms
      implicit none
      !
      type(BosonicField),intent(inout)      :: Pout
      type(FermionicField),intent(in)       :: G
      integer,intent(in)                    :: Ntau
      logical,intent(in),optional           :: tau_uniform
      logical,intent(in),optional           :: tau_output
      !
      real(8),allocatable                   :: tau(:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,NaxisB
      integer                               :: Norb,NmatsF,Ntau_
      integer                               :: ik1,ik2,iq,iw,ispin
      integer                               :: iwan1,iwan2,ib1,ib2
      logical                               :: tau_uniform_
      logical                               :: tau_output_
      !
      !
      write(*,*) "--- calc_Pi_selfcons ---"
      !
      !
      ! Check on the input Fields
      if(.not.Pout%status) stop "Pout not properly initialized."
      if(.not.G%status) stop "Green's function not properly initialized."
      if(Pout%Nkpt.eq.0) stop "Pout k dependent attributes not properly initialized."
      if(Gin%Nkpt.eq.0) stop "Green's function k dependent attributes not properly initialized."
      !
      Nbp = Pout%Nbp
      Nkpt = Pout%Nkpt
      Beta = Pout%Beta
      NaxisB = Pout%Npoints
      Norb = int(sqrt(dble(Nbp)))
      NmatsF = Gin%Npoints
      if(Gin%Norb.ne.Norb) stop "Pout and Green's function have different orbital dimension."
      if(Gin%Nkpt.ne.Nkpt) stop "Pout and Green's function have different number of Kpoints."
      !
      tau_uniform_=.false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      tau_output_=.false.
      if(present(tau_output)) tau_output_ = tau_output
      Ntau_ = Ntau
      if(tau_output_) Ntau_ = NaxisB
      !
      allocate(tau(Ntau_));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau_)
      else
         tau = denspace(beta,Ntau_)
      endif
      !
      !
      Pout%bare=czero; Pout%bare_local=czero
      Pout%screened=czero; Pout%screened_local=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(nkpt,kptdiff,iq,nbasis,gtau_wann,Pq_tau,ntau,npb,tau,nspin_DMFT),&
      !$OMP PRIVATE(ik1,ik2,m,n,mp,np,ipb1,ipb2,iwan1,iwan2,itau,tau2,ispin)
      !$OMP DO
      do iq=1,Nkpt
         do itau=1,Ntau_
            !
            tau2=tau(Ntau_)-tau(itau)
            if (dabs(tau2-tau(Ntau_-itau+1)).gt.eps) stop "itau2"
            !
            do ispin=1,Nspin
               do ik1=1,Nkpt
                  ik2=kptdiff(ik1,iq)
                  !
                  do m=1,Norb
                     do n=1,Norb
                        do mp=1,Norb
                           do np=1,Norb
                              !
                              ib1=mp+nbasis*(m-1)
                              ib2=np+nbasis*(n-1)
                              !
                              Pq_tau(ipb1,ipb2,itau) = Pq_tau(ipb1,ipb2,itau)    &
                              -gtau_wann(m,n,itau,ik1,ispin) * gtau_wann(np,mp,ntau-itau+1,ik2,ispin)/nkpt








        !
        Pq_tau(:,:,itau)=0.d0


                          !

                          if (nspin_DMFT.eq.1) then
                             Pq_tau(ipb1,ipb2,itau)=Pq_tau(ipb1,ipb2,itau)-2.d0*gtau_wann(m,n,itau,ik1,1)*gtau_wann(np,mp,ntau-itau+1,ik2,1)/nkpt
                          else

                          endif
                          !
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
        !
     enddo
     !$OMP END DO
     !$OMP END PARALLEL



















      endif
      !
      if(.not.tau_output_)then

      endif
      !







         do ispin=1,Nspin
            call Fmats2itau_mat()
            call Fmats2itau_mat()









      if(all([Umats%Nbp-Nbp,Pout%Nbp-Nbp].ne.[0,0])) stop "Either Umats and/or Pmats have different orbital dimension with respect to Wmats."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "Either Umats and/or Pmats have different Beta with respect to Wmats."
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))  write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Wmats%Npoints,Umats%Npoints,Pmats%Npoints])




      allocate(wmats(Nmats));wmats=0d0
      wmats = BosonicFreqMesh(Beta,Nmats)
      !
      !cprod(alpha,i,n,ik)= < B_q,alpha Psi_kn |Psi_q+k,i>
      allocate(cprod(Nbp,Norb,Norb,Nkpt));cprod=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Nmats,Norb,Lttc,cprod),&
      !$OMP PRIVATE(iq,ik1,ik2,ib1,ib2,ipb1,iwan1,iwan2)
      !$OMP DO
      do iq=1,Nkpt
         do ik1=1,Nkpt
            ik2 = Lttc%kptsum(ik1,iq)
            !
            do ib1=1,Norb
               do ib2=1,Norb
                  ipb1=0
                  do iwan2=1,Norb
                     do iwan1=1,Norb
                        !
                        ipb1=ipb1+1
                        cprod(ipb1,ib2,ib1,ik1) = dconjg(Lttc%Zk(iwan1,ib1,ik1))*Lttc%Zk(iwan2,ib2,ik2)
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
      Wmats%bare=czero; Wmats%bare_local=czero
      Wmats%screened=czero; Wmats%screened_local=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Norb,wmats,cprod,alpha,Lttc,Pmats),&
      !$OMP PRIVATE(iw,iq,ik1,ik2,iwan1,iwan2)
      !$OMP DO
      do iw=1,Nmats
         do iq=1,Nkpt
            do ik1=1,Nkpt
               ik2 = Lttc%kptsum(ik1,iq)
               !
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     !
                     if (dabs(-Lttc%Ek(iwan1,ik1)+Lttc%Ek(iwan2,ik2)).lt.eps.and.iw.eq.1) then
                        !
                        !lim_{E'->E} (n(E)-n(E'))/(E'-E) = (n(E) - (n(E) + n'(E)*(E'-E)))/(E'-E) -> -n'(E)
                        alpha(iw) = +2.d0 * diff_fermidirac(Lttc%Ek(iwan1,ik1),Lttc%mu) / Nkpt
                        !
                     else
                        !
                        alpha(iw) = -2.d0 * ( fermidirac(Lttc%Ek(iwan1,ik1),Lttc%mu) - fermidirac(Lttc%Ek(iwan2,ik2),Lttc%mu) ) &
                                          / ( dcmplx(0d0,1d0) * wmats(iw) - Lttc%Ek(iwan1,ik1) + Lttc%Ek(iwan2,ik2) ) / nkpt
                        !
                     endif
                     !
                     !alpha(iw)=alpha(iw)!/nkpt/nkpt !to account for the fact that cprod basis = 1/nkpt * prod basis of U in spex
                     call ZHER('U',Nbp,alpha(iw),cprod(:,iwan2,iwan1,ik1),1,Pmats%screened(:,:,iw,iq),Nbp)
                     !
                  enddo
               enddo
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(cprod,alpha,wmats)
      !
   end subroutine calc_Pi_selfcons





end module bubbles
