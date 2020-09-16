module self_energy

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !interface calc_Pi
   !   module procedure calc_Pi_GkGk
   !   module procedure calc_Pi_selfcons
   !end interface calc_Pi

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: calc_sigmaGW

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute the two GW self-energy components.
   !---------------------------------------------------------------------------!
   subroutine calc_sigmaGW(Smats_C,Smats_X,Gmats,Wmats,Lttc,tau_uniform)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use crystal
      use fourier_transforms
      use global_vars, only : Ntau
      implicit none
      !
      type(FermionicField),intent(inout)    :: Smats_C
      type(FermionicField),intent(inout)    :: Smats_X
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Wmats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in),optional           :: tau_uniform
      !
      complex(8),allocatable                :: Sitau(:,:,:,:,:)
      complex(8),allocatable                :: Gitau(:,:,:,:,:)
      complex(8),allocatable                :: Witau(:,:,:,:)
      real(8),allocatable                   :: tau(:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Norb,Nmats
      integer                               :: ik1,ik2,iq,iw,itau,ispin
      integer                               :: m,n,mp,np,ib1,ib2
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- calc_sigmaGW ---"
      !
      !
      ! Check on the input Fields
      if(.not.Smats_C%status) stop "Smats_C not properly initialized."
      if(.not.Smats_X%status) stop "Smats_X not properly initialized."
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(.not.Wmats%status) stop "Wmats not properly initialized."
      if(Smats_C%Nkpt.eq.0) stop "Smats_C k dependent attributes not properly initialized."
      if(Smats_X%Nkpt.eq.0) stop "Smats_X k dependent attributes not properly initialized."
      if(Smats_X%Npoints.ne.0) stop "Smats_X frequency dependent attributes are supposed to be unallocated."
      if(Gmats%Nkpt.eq.0) stop "Gmats k dependent attributes not properly initialized."
      if(Wmats%Nkpt.eq.0) stop "Wmats k dependent attributes not properly initialized."
      if(.not.allocated(Lttc%kptdif)) stop "kptdif not allocated."
      if(.not.allocated(Lttc%kptPos)) stop "kptPos not allocated."
      !
      Norb = Smats_C%Norb
      Nkpt = Smats_C%Nkpt
      Beta = Smats_C%Beta
      Nmats = Smats_C%Npoints
      Nbp = Norb**2
      !
      if(all([Smats_X%Norb-Norb,Gmats%Norb-Norb,Wmats%Nbp-Nbp].ne.[0,0,0])) stop "Either Smats_X, Gmats or Wmats have different orbital dimension with respect to Smats_C."
      if(all([Smats_X%Beta-Beta,Gmats%Beta-Beta,Wmats%Beta-Beta].ne.[0d0,0d0,0d0])) stop "Either Smats_X, Gmats or Wmats have different Beta with respect to Smats_C."
      if(all([Gmats%Npoints-Nmats,Wmats%Npoints-Nmats].ne.[0,0]))  write(*,"(A)") "Warning: Either Smats_C, Gmats or Wmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Smats_C%Npoints,Wmats%Npoints,Wmats%Npoints])
      !
      tau_uniform_=.false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(tau(Ntau));tau=0d0
      if(tau_uniform_)then
         tau = linspace(0d0,beta,Ntau)
      else
         tau = denspace(beta,Ntau)
      endif
      !
      allocate(Gitau(Norb,Norb,Ntau,Nkpt,Nspin));Gitau=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin), &
         asympt_corr=.true.,tau_uniform=tau_uniform_,Nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt)
      enddo
      !
      allocate(Witau(Nbp,Nbp,Ntau,Nkpt));Witau=czero
      call Bmats2itau(Beta,Wmats%screened,Witau,asympt_corr=.true.,tau_uniform=tau_uniform_)
      do itau=1,Ntau
         Witau(:,:,Ntau,:) = Witau(:,:,Ntau,:) - Wmats%bare
      enddo
      !
      write(*,"(A)") "Sigma_C(tau)"
      !Sigma_{m,n}(q,tau) = -Sum_{k,mp,np} W_{(m,mp);(n,np)}(q-k;tau)G_{mp,np}(k,tau)
      allocate(Sitau(Norb,Norb,Ntau,Nkpt,Nspin));Sitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,Ntau,Lttc,Nkpt,Sitau,Gitau,Witau),&
      !$OMP PRIVATE(m,n,itau,iq,ispin,ik1,ik2,mp,np,ib1,ib2)
      !$OMP DO
      do m=1,Norb
         do n=1,Norb
            do itau=1,Ntau
               do iq=1,Lttc%Nkpt_irred
                  do ispin=1,Nspin
                     !
                     do ik1=1,Nkpt
                        ik2=Lttc%kptdif(iq,ik1)
                        do mp=1,Norb
                           do np=1,Norb
                              !
                              ib1 = mp + Norb*(m-1)
                              ib2 = np + Norb*(n-1)
                              !
                              Sitau(m,n,itau,iq,ispin) = Sitau(m,n,itau,iq,ispin) - Gitau(mp,np,itau,ik1,ispin)*Witau(ib1,ib2,itau,ik2)/Nkpt
                              !
                           enddo
                        enddo
                     enddo
                     !
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Gitau,Witau)
      !
      write(*,"(A)") "Sigma_C(tau)-->Sigma_C(iw)"
      call clear_attributes(Smats_C)
      do ispin=1,Nspin
         do iq=1,Lttc%Nkpt_irred
            call Fitau2mats_mat(Beta,Sitau(:,:,:,iq,ispin),Smats_C%wks(:,:,:,iq,ispin),tau_uniform=tau_uniform_)
         enddo
      enddo
      deallocate(Sitau)
      !
      write(*,"(A)") "Sigma_X"
      call clear_attributes(Smats_X)
      !sigmax(r,r')=-g(r,r',tau=0-)*v(r-r')
      !Sigmax_nm(q) = Sum_kij V_{ni,jm}(q-k)G_ij(k,beta)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,Ntau,Lttc,Nkpt,Smats_X,Gitau,Wmats),&
      !$OMP PRIVATE(m,n,iq,ispin,ik1,ik2,mp,np,ib1,ib2)
      !$OMP DO
      do m=1,Norb
         do n=1,Norb
            do iq=1,Lttc%Nkpt_irred
               do ispin=1,Nspin
                  !
                  do ik1=1,Nkpt
                     ik2=Lttc%kptdif(iq,ik1)
                     do mp=1,Norb
                        do np=1,Norb
                           !
                           ib1 = mp + Norb*(m-1)
                           ib2 = np + Norb*(n-1)
                           !
                           Smats_X%ks(m,n,iq,ispin) = Smats_X%ks(m,n,iq,ispin) + Gitau(mp,np,Ntau,ik1,ispin)*Wmats%bare(ib1,ib2,ik2)/Nkpt
                           !
                        enddo
                     enddo
                  enddo
                  !
               enddo
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      if(Lttc%Nkpt_irred.lt.Nkpt) then
         !sigma(ik)=sigma(kptp(ik))
         write(*,"(A)") "transformation to lda eigenbasis and back"
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nmats,Lttc,Nkpt,Smats_X,Smats_C),&
         !$OMP PRIVATE(ispin,iw,iq)
         !$OMP DO
         do ispin=1,Nspin
            !
            !rotation to lda eigenbasis
            do iq=1,Lttc%Nkpt_irred
               Smats_X%ks(:,:,iq,ispin) = rotate(Smats_X%ks(:,:,iq,ispin),Lttc%Zk(:,:,iq))
               do iw=1,Nmats
                  Smats_C%wks(:,:,iw,iq,ispin) = rotate(Smats_C%wks(:,:,iw,iq,ispin),Lttc%Zk(:,:,iq))
               enddo
            enddo
            !
            !fill up the missing Kpoints
            do iq=1,nkpt
               Smats_X%ks(:,:,iq,ispin) = Smats_X%ks(:,:,Lttc%kptPos(iq),ispin)
               do iw=1,Nmats
                  Smats_C%wks(:,:,iw,iq,ispin) = Smats_C%wks(:,:,iw,Lttc%kptPos(iq),ispin)
               enddo
            enddo
            !
            !rotate back
            do iq=1,nkpt
               Smats_X%ks(:,:,iq,ispin) = rotate(Smats_X%ks(:,:,iq,ispin),transpose(conjg(Lttc%Zk(:,:,iq))))
               do iw=1,Nmats
                  Smats_C%wks(:,:,iw,iq,ispin) = rotate(Smats_C%wks(:,:,iw,iq,ispin),transpose(conjg(Lttc%Zk(:,:,iq))))
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
      call FermionicKsum(Smats_C)
      call FermionicKsum(Smats_X)
      !
   end subroutine calc_sigmaGW


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute Hartree difference with G0W0
   !---------------------------------------------------------------------------!
   subroutine calc_VH(density_GoWo,Gmats,Umats,VH)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      !use greens_functions
      use global_vars, only :  pathINPUT, VH_type
      implicit none
      !
      complex(8),intent(in)                 :: density_GoWo(:,:)
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Umats
      complex(8),intent(inout)              :: VH(:,:)
      !
      complex(8),allocatable                :: density(:,:)
      complex(8),allocatable                :: Vgamma(:,:),Vevec(:,:)
      real(8),allocatable                   :: Veval(:)
      integer                               :: Norb,Nbp
      integer                               :: ib1,ib2,ib3
      integer                               :: m,n,mp,np
      logical                               :: filexists
      character(len=256)                    :: path
      !
      !
      write(*,"(A)") "--- calc_VH ---"
      !
      !
      ! Check on the input Field
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(.not.Umats%status) stop "Umats not properly initialized."
      if(Gmats%Nkpt.eq.0) stop "Gmats k dependent attributes not properly initialized."
      if(Umats%Nkpt.eq.0) stop "Umats k dependent attributes not properly initialized."
      if(Umats%iq_gamma.lt.0) stop "Umats iq_gamma not defined."
      Norb = Gmats%Norb
      Nbp = Umats%Nbp
      !
      allocate(density(Norb,Norb))
      !call calc_density(density,Gmats)
      !
      select case(VH_type)
         case default
            stop "Available VH_type are: Ubare, Ustatic, Ubare_SPEX, Ustatic_SPEX."
         case("Ubare")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            Vgamma = Umats%bare(:,:,Umats%iq_gamma)
            !
            allocate(Veval(Nbp));Veval=0d0
            allocate(Vevec(Nbp,Nbp));Vevec=Vgamma
            call eigh(Vevec,Veval)
            !
            do ib1=1,Nbp
               do ib2=1,Nbp
                  do ib3=1,Nbp
                     Vgamma(ib1,ib2) = Vgamma(ib1,ib2) + Vevec(ib1,ib3)*Veval(ib3)*conjg(Vevec(ib2,ib3))
                  enddo
               enddo
            enddo
            !
         case("Ustatic")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            Vgamma = Umats%screened(:,:,1,Umats%iq_gamma)
            !
            allocate(Veval(Nbp));Veval=0d0
            allocate(Vevec(Nbp,Nbp));Vevec=Vgamma
            call eigh(Vevec,Veval)
            !
            do ib1=1,Nbp
               do ib2=1,Nbp
                  do ib3=1,Nbp
                     Vgamma(ib1,ib2) = Vgamma(ib1,ib2) + Vevec(ib1,ib3)*Veval(ib3)*conjg(Vevec(ib2,ib3))
                  enddo
               enddo
            enddo
            !
         case("Ubare_SPEX")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            path = pathINPUT//"V_nodiv.DAT"
            call inquireFile(reg(path),filexists)
            call read_Vgamma(-1)
            !
         case("Ustatic_SPEX")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            path = pathINPUT//"V_nodiv.DAT"
            call inquireFile(reg(path),filexists)
            call read_Vgamma(0)
            !
      end select
      !
      VH=czero
      do mp=1,Norb
         do  m=1,Norb
           do np=1,Norb
               do n=1,Norb
                 !
                 ib1 = mp + Norb*(m-1)
                 ib2 = np + Norb*(n-1)
                 !
                 VH(mp,m) = VH(mp,m) + real( (density(np,n)-density_GoWo(np,n))*Vgamma(ib1,ib2) )
                 !
               enddo
           enddo
         enddo
      enddo
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
            if(idum1.ne.Norb) stop "read_Vgamma Norb"
            do limits=1,2
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     do iwan3=1,Norb
                       do iwan4=1,Norb
                           indx1=iwan1+Norb*(iwan2-1)
                           indx2=iwan3+Norb*(iwan4-1)
                           read(unit,*) rdum1,rdum2,idum1,idum2,idum3,idum4,rdum3,rdum4
                           if (idum1.ne.iwan1) stop "read_Vgamma: iwan1"
                           if (idum2.ne.iwan2) stop "read_Vgamma: iwan2"
                           if (idum3.ne.iwan3) stop "read_Vgamma: iwan3"
                           if (idum4.ne.iwan4) stop "read_Vgamma: iwan4"
                           if(dble(Vtype).eq.rdum1) Vgamma(indx1,indx2)=dcmplx(rdum3,rdum4)
                       enddo
                     enddo
                  enddo
               enddo
               if(Vtype.eq.-1) exit
            enddo
         end subroutine read_Vgamma
      !
   end subroutine calc_VH
   !
end module self_energy
