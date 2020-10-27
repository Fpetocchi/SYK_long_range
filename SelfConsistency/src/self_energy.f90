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
   !

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
   !TEST ON: 27-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_sigmaGW(Smats_C,Smats_X,Gmats,Wmats,Lttc)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use crystal
      use fourier_transforms
      use input_vars, only : NtauB, tau_uniform
      implicit none
      !
      type(FermionicField),intent(inout)    :: Smats_C
      type(FermionicField),intent(inout)    :: Smats_X
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Wmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: Sitau(:,:,:,:,:)
      complex(8),allocatable                :: Gitau(:,:,:,:,:)
      complex(8),allocatable                :: Witau(:,:,:,:)
      complex(8),allocatable                :: WmatsC(:,:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Norb,Nmats
      integer                               :: ik1,ik2,iq,iw,itau,ispin
      integer                               :: m,n,mp,np,ib1,ib2
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- calc_sigmaGW"
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
      if(all([Lttc%Nkpt-Nkpt,Gmats%Nkpt-Nkpt,Wmats%Nkpt-Nkpt].ne.[0,0,0])) stop "Either Lattice, Gmats or Wmats have different number of k-points with respect to Smats_C."
      if(all([Smats_X%Norb-Norb,Gmats%Norb-Norb,Wmats%Nbp-Nbp].ne.[0,0,0])) stop "Either Smats_X, Gmats or Wmats have different orbital dimension with respect to Smats_C."
      if(all([Smats_X%Beta-Beta,Gmats%Beta-Beta,Wmats%Beta-Beta].ne.[0d0,0d0,0d0])) stop "Either Smats_X, Gmats or Wmats have different Beta with respect to Smats_C."
      if(all([Gmats%Npoints-Nmats,Wmats%Npoints-Nmats].ne.[0,0])) write(*,"(A)") "Warning: Either Smats_C, Gmats or Wmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Smats_C%Npoints,Wmats%Npoints,Wmats%Npoints])
      !
      call cpu_time(start)
      !
      allocate(Gitau(Norb,Norb,NtauB,Nkpt,Nspin));Gitau=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin), &
         asympt_corr=.true.,tau_uniform=tau_uniform,Nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt)
      enddo
      !
      allocate(Witau(Nbp,Nbp,NtauB,Nkpt));Witau=czero
      allocate(WmatsC(Nbp,Nbp,Nmats,Nkpt));WmatsC=czero
      do iw=1,Nmats
         WmatsC(:,:,iw,:) = Wmats%screened(:,:,iw,:) - Wmats%bare
      enddo
      call Bmats2itau(Beta,WmatsC,Witau,asympt_corr=.true.,tau_uniform=tau_uniform)
      deallocate(WmatsC)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "Glat(ik,iw),Wlat(iq,iw) --> Glat(ik,itau),Wlat(iq,itau) cpu timing:", finish-start
      !
      call cpu_time(start)
      !Sigma_{m,n}(q,tau) = -Sum_{k,mp,np} W_{(m,mp);(n,np)}(q-k;tau)G_{mp,np}(k,tau)
      allocate(Sitau(Norb,Norb,NtauB,Lttc%Nkpt_irred,Nspin));Sitau=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,NtauB,Lttc,Nkpt,Sitau,Gitau,Witau),&
      !$OMP PRIVATE(m,n,itau,iq,ispin,ik1,ik2,mp,np,ib1,ib2)
      !$OMP DO
      do m=1,Norb
         do n=1,Norb
            do itau=1,NtauB
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
      deallocate(Witau)
      !
      call clear_attributes(Smats_C)
      do ispin=1,Nspin
         do iq=1,Lttc%Nkpt_irred
            call Fitau2mats_mat(Beta,Sitau(:,:,:,iq,ispin),Smats_C%wks(:,:,:,iq,ispin),tau_uniform=tau_uniform)
         enddo
      enddo
      deallocate(Sitau)
      call cpu_time(finish)
      write(*,"(A,F)") "Sigma_C(ik,iw) cpu timing:", finish-start
      !
      call cpu_time(start)
      call clear_attributes(Smats_X)
      !sigmax(r,r')=-g(r,r',tau=0-)*v(r-r')
      !Sigmax_nm(q) = Sum_kij V_{ni,jm}(q-k)G_ij(k,beta)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,NtauB,Lttc,Nkpt,Smats_X,Gitau,Wmats),&
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
                           Smats_X%N_ks(m,n,iq,ispin) = Smats_X%N_ks(m,n,iq,ispin) + Gitau(mp,np,NtauB,ik1,ispin)*Wmats%bare(ib1,ib2,ik2)/Nkpt
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
      deallocate(Gitau)
      call cpu_time(finish)
      write(*,"(A,F)") "Sigma_X(ik) cpu timing:", finish-start
      !
      if(Lttc%Nkpt_irred.lt.Nkpt) then
         !sigma(ik)=sigma(kptp(ik))
         write(*,"(A)") "Transformation to lda eigenbasis and back."
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nmats,Lttc,Nkpt,Smats_X,Smats_C),&
         !$OMP PRIVATE(ispin,iw,iq)
         !$OMP DO
         do ispin=1,Nspin
            !
            !rotation to lda eigenbasis
            do iq=1,Lttc%Nkpt_irred
               Smats_X%N_ks(:,:,iq,ispin) = rotate(Smats_X%N_ks(:,:,iq,ispin),Lttc%Zk(:,:,iq))
               do iw=1,Nmats
                  Smats_C%wks(:,:,iw,iq,ispin) = rotate(Smats_C%wks(:,:,iw,iq,ispin),Lttc%Zk(:,:,iq))
               enddo
            enddo
            !
            !fill up the missing Kpoints
            do iq=1,nkpt
               Smats_X%N_ks(:,:,iq,ispin) = Smats_X%N_ks(:,:,Lttc%kptPos(iq),ispin)
               do iw=1,Nmats
                  Smats_C%wks(:,:,iw,iq,ispin) = Smats_C%wks(:,:,iw,Lttc%kptPos(iq),ispin)
               enddo
            enddo
            !
            !rotate back
            do iq=1,nkpt
               Smats_X%N_ks(:,:,iq,ispin) = rotate(Smats_X%N_ks(:,:,iq,ispin),transpose(conjg(Lttc%Zk(:,:,iq))))
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
   !PURPOSE: Compute the two local GW self-energy components as Gloc*Wloc.
   !TEST ON: 27-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_sigmaGWdc(Smats_Cdc,Smats_Xdc,Gmats,Wmats)
      !
      use parameters
      use linalg
      use utils_misc
      use utils_fields
      use crystal
      use fourier_transforms
      use input_vars, only : NtauB, tau_uniform
      implicit none
      !
      type(FermionicField),intent(inout)    :: Smats_Cdc
      type(FermionicField),intent(inout)    :: Smats_Xdc
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Wmats
      !
      complex(8),allocatable                :: Sitau_loc(:,:,:,:)
      complex(8),allocatable                :: Gitau_loc(:,:,:,:)
      complex(8),allocatable                :: Witau_loc(:,:,:)
      complex(8),allocatable                :: WmatsC_loc(:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Norb,Nmats
      integer                               :: iw,itau,ispin
      integer                               :: m,n,mp,np,ib1,ib2
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- calc_sigmaGWdc"
      !
      !
      ! Check on the input Fields
      if(.not.Smats_Cdc%status) stop "Smats_Cdc not properly initialized."
      if(.not.Smats_Xdc%status) stop "Smats_Xdc not properly initialized."
      if(.not.Gmats%status) stop "Gmats not properly initialized."
      if(.not.Wmats%status) stop "Wmats not properly initialized."
      if(Smats_Cdc%Nkpt.ne.0) stop "Smats_Cdc k dependent attributes are supposed to be unallocated."
      if(Smats_Xdc%Nkpt.ne.0) stop "Smats_Cdc k dependent attributes are supposed to be unallocated."
      if(Smats_Xdc%Npoints.ne.0) stop "Smats_Xdc frequency dependent attributes are supposed to be unallocated."
      !
      Norb = Smats_Cdc%Norb
      Nkpt = Smats_Cdc%Nkpt
      Beta = Smats_Cdc%Beta
      Nmats = Smats_Cdc%Npoints
      Nbp = Norb**2
      !
      if(all([Smats_Xdc%Norb-Norb,Gmats%Norb-Norb,Wmats%Nbp-Nbp].ne.[0,0,0])) stop "Either Smats_Xdc, Gmats or Wmats have different orbital dimension with respect to Smats_Cdc."
      if(all([Smats_Xdc%Beta-Beta,Gmats%Beta-Beta,Wmats%Beta-Beta].ne.[0d0,0d0,0d0])) stop "Either Smats_Xdc, Gmats or Wmats have different Beta with respect to Smats_Cdc."
      if(all([Gmats%Npoints-Nmats,Wmats%Npoints-Nmats].ne.[0,0]))  write(*,"(A)") "Warning: Either Smats_Cdc, Gmats or Wmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Smats_Cdc%Npoints,Wmats%Npoints,Wmats%Npoints])
      !
      call cpu_time(start)
      !
      allocate(Gitau_loc(Norb,Norb,NtauB,Nspin));Gitau_loc=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%ws(:,:,:,ispin),Gitau_loc(:,:,:,ispin), &
                             asympt_corr=.true.,tau_uniform=tau_uniform)
      enddo
      !
      allocate(Witau_loc(Nbp,Nbp,NtauB));Witau_loc=czero
      allocate(WmatsC_loc(Nbp,Nbp,Nmats));WmatsC_loc=czero
      do iw=1,Nmats
         WmatsC_loc(:,:,iw) = Wmats%screened_local(:,:,iw) - Wmats%bare_local
      enddo
      call Bmats2itau(Beta,WmatsC_loc,Witau_loc,asympt_corr=.true.,tau_uniform=tau_uniform)
      deallocate(WmatsC_loc)
      !
      call cpu_time(finish)
      write(*,"(A,F)") "Glat(iw),Wlat(iw) --> Glat(itau),Wlat(itau) cpu timing:", finish-start
      !
      call cpu_time(start)
      !Sigma_{m,n}(q,tau) = -Sum_{k,mp,np} W_{(m,mp);(n,np)}(q-k;tau)G_{mp,np}(k,tau)
      allocate(Sitau_loc(Norb,Norb,NtauB,Nspin));Sitau_loc=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,NtauB,Sitau_loc,Gitau_loc,Witau_loc),&
      !$OMP PRIVATE(m,n,itau,ispin,mp,np,ib1,ib2)
      !$OMP DO
      do m=1,Norb
         do n=1,Norb
            do itau=1,NtauB
               do ispin=1,Nspin
                  !
                  do mp=1,Norb
                     do np=1,Norb
                        !
                        ib1 = mp + Norb*(m-1)
                        ib2 = np + Norb*(n-1)
                        !
                        Sitau_loc(m,n,itau,ispin) = Sitau_loc(m,n,itau,ispin) - Gitau_loc(mp,np,itau,ispin)*Witau_loc(ib1,ib2,itau)
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Witau_loc)
      !
      call clear_attributes(Smats_Cdc)
      do ispin=1,Nspin
         call Fitau2mats_mat(Beta,Sitau_loc(:,:,:,ispin),Smats_Cdc%ws(:,:,:,ispin),tau_uniform=tau_uniform)
      enddo
      deallocate(Sitau_loc)
      call cpu_time(finish)
      write(*,"(A,F)") "Sigma_Cdc(iw) cpu timing:", finish-start
      !
      call cpu_time(start)
      call clear_attributes(Smats_Xdc)
      !sigmax(r,r')=-g(r,r',tau=0-)*v(r-r')
      !Sigmax_nm(q) = Sum_kij V_{ni,jm}(q-k)G_ij(k,beta)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,NtauB,Smats_Xdc,Gitau_loc,Wmats),&
      !$OMP PRIVATE(m,n,ispin,mp,np,ib1,ib2)
      !$OMP DO
      do m=1,Norb
         do n=1,Norb
            do ispin=1,Nspin
               !
               do mp=1,Norb
                  do np=1,Norb
                     !
                     ib1 = mp + Norb*(m-1)
                     ib2 = np + Norb*(n-1)
                     !
                     Smats_Xdc%N_s(m,n,ispin) = Smats_Xdc%N_s(m,n,ispin) + Gitau_loc(mp,np,NtauB,ispin)*Wmats%bare_local(ib1,ib2)
                     !
                  enddo
               enddo
               !
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Gitau_loc)
      call cpu_time(finish)
      write(*,"(A,F)") "Sigma_Xdc cpu timing:", finish-start
      !
   end subroutine calc_sigmaGWdc


   !---------------------------------------------------------------------------!
   !PURPOSE: Compute Hartree difference with G0W0
   !TEST ON: 27-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_VH(density_LDA,Gmats,Umats,VH)
      !
      use parameters
      use linalg
      use file_io
      use utils_misc
      use utils_fields
      use greens_function, only : calc_density
      use input_vars, only :  pathINPUT, VH_type
      implicit none
      !
      complex(8),intent(in)                 :: density_LDA(:,:)
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Umats
      complex(8),intent(inout)              :: VH(:,:)
      !
      complex(8),allocatable                :: density(:,:),density_spin(:,:,:)
      complex(8),allocatable                :: Vgamma(:,:),Vevec(:,:)
      real(8),allocatable                   :: Veval(:)
      integer                               :: Norb,Nbp
      integer                               :: ib1,ib2,ib3
      integer                               :: m,n,mp,np
      logical                               :: filexists
      character(len=256)                    :: path
      !
      !
      if(verbose)write(*,"(A)") "---- calc_VH"
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
      if((Norb**2).ne.Nbp) stop "Umats and Gmats have different orbital space."
      call assert_shape(density_LDA,[Norb,Norb],"calc_VH","density_LDA")
      call assert_shape(VH,[Norb,Norb],"calc_VH","VH")
      !
      allocate(density(Norb,Norb));density=0d0
      allocate(density_spin(Norb,Norb,Nspin));density_spin=0d0
      call calc_density(Gmats,density_spin)
      density = sum(density_spin,dim=3)
      deallocate(density_spin)
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
            path = reg(pathINPUT)//"V_nodiv.DAT"
            call inquireFile(reg(path),filexists)
            call read_Vgamma(-1)
            !
         case("Ustatic_SPEX")
            !
            allocate(Vgamma(Nbp,Nbp));Vgamma=czero
            path = reg(pathINPUT)//"V_nodiv.DAT"
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
                 VH(mp,m) = VH(mp,m) + real( (density(np,n)-density_LDA(np,n))*Vgamma(ib1,ib2) )
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
      !
   end subroutine calc_VH



   !---------------------------------------------------------------------------!
   !PURPOSE: Read self-energy from SPEX files.
   !TEST ON: 27-10-2020(both write and read)
   !---------------------------------------------------------------------------!
   subroutine read_Sigma_spex(Smats_GoWo,Lttc,save2readable,Vxc,pathOUTPUT,doAC)
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
      complex(8),intent(inout),optional     :: Vxc(:,:,:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: doAC
      !
      logical                               :: filexists,ACdone,doAC_,Vxcdone,doVxc
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
      !
      !
      if(verbose)write(*,"(A)") "---- read_Sigma_spex"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Check on the input Fields
      if(.not.Smats_GoWo%status) stop "FermionicField not properly initialized."
      if(.not.Lttc%status) stop "Lattice container not properly initialized."
      if(Lttc%Nkpt.ne.Smats_GoWo%Nkpt) stop "Lattice has different number of k-points with respect to Smats_GoWo."
      !
      Norb = Smats_GoWo%Norb
      Nkpt = Smats_GoWo%Nkpt
      Nmats = Smats_GoWo%Npoints
      !
      allocate(wmats(Smats_GoWo%Npoints));wmats=0d0
      wmats = FermionicFreqMesh(Smats_GoWo%Beta,Smats_GoWo%Npoints)
      !
      ! Read XEPS data
      if(.not.XEPSisread)then
         path = reg(pathINPUT)//"XEPS.DAT"
         call inquireFile(reg(path),filexists)
         call read_xeps(reg(path),Lttc%kpt,Lttc%Nkpt3,UseXepsKorder, &
         Lttc%kptPos,Lttc%Nkpt_irred,Lttc%UseDisentangledBS,Lttc%iq_gamma,paramagneticSPEX)
      endif
      !
      ! Check if the data on the Matsubara axis are present if(.not.paramagneticSPEX) look also for spin2
      path = reg(pathINPUT)//"Sigma_GoWo_k_up.DAT"
      call inquireFile(reg(path),ACdone,hardstop=.false.)
      doAC_ = .not.ACdone
      if(present(doAC)) doAC_ = doAC
      !
      ! Check if the Vxc_wann is present
      path = reg(pathINPUT)//"Vxc_wann_k_up.DAT"
      call inquireFile(reg(path),Vxcdone,hardstop=.false.)
      doVxc = .not.Vxcdone .and. present(Vxc)
      if(doVxc.and.(.not.doAC_))then
         call assert_shape(Vxc,[Norb,Norb,Nkpt,Nspin],"read_Sigma_spex","Vxc")
         write(*,"(A)")"Sorry but I can't produce Vxc_wann without reading the self-energy."
         write(*,"(A)")"Analytic continuation will be perforemd anyway."
         doAC_ = .true.
      endif
      !
      !
      ! Perform cnalytical continuation on the self-energy on the real axis
      if(doAC_)then
         !
         !---------------------------------------------------------------------!
         !
         ! Read UWAN file
         if(Lttc%UseDisentangledBS)then
            path = reg(pathINPUT)//"UWAN_NEW.DAT"
         else
            path = reg(pathINPUT)//"UWAN.DAT"
         endif
         call inquireFile(reg(path),filexists)
         if(verbose)write(*,"(A)") "Opening "//reg(path)
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
         read(unit) Nspin_Uwan,Nkpt_Uwan,ib_Uwan1,ib_Uwan2,Norb_Uwan
         if(paramagneticSPEX.and.(Nspin_Uwan.ne.1)) stop "UWAN file is not paramagnetic."
         if(Nkpt_Uwan.ne.Nkpt) stop "UWAN file has wrong number of k-points (not irreducible)."
         if(Norb_Uwan.ne.Norb) stop "UWAN file has wrong orbital dimension."
         write(*,"(A,2I4)") "The band indexes in the UWAN rotation are: ",ib_Uwan1,ib_Uwan2
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
            call inquireDir(reg(path),filexists,hardstop=.false.)
            if(.not.filexists) exit
            SigmaSegments = SigmaSegments + 1
         enddo
         write(*,"(A,1I6)") "The number of segment of the SPEX self-energy is: ",SigmaSegments
         allocate(NfreqSeg(SigmaSegments));NfreqSeg=0
         !
         ! Look for the Number of Kpoints in each segment (supposed to be the same of Lttc%Nkpt_irred)
         do iseg=1,SigmaSegments
            Nkpt_file = 0
            do ik=1,2000
               path = reg(pathINPUT)//"Sigma_real_"//str(iseg,2)//"/SIGMA.Q"//str(ik,4)//".DAT"
               call inquireFile(reg(path),filexists,hardstop=.false.)
               if(.not.filexists) exit
               Nkpt_file = Nkpt_file + 1
            enddo
            write(*,"(A,1I4,A,1I6)") "The number k-points in the segment number: ",iseg," is: ",Nkpt_file
            Nkpt_file_old = Nkpt_file
            if((iseg.gt.1).and.(Nkpt_file.ne.Nkpt_file_old)) stop "Number of K-points does not match among segments."
            if(Nkpt_file.ne.Lttc%Nkpt_irred) stop "Number of K-points does not match with Nkpt_irred readed from XEPS."
         enddo
         !
         ! Check that all the Sigma parameters are cosistent and that the segments match
         do iseg=1,SigmaSegments
            do ik=1,Lttc%Nkpt_irred
               !
               path = reg(pathINPUT)//"Sigma_real_"//str(iseg,2)//"/SIGMA.Q"//str(ik,4)//".DAT"
               write(*,"(A)") "Checking "//reg(path)
               call inquireFile(reg(path),filexists) !redundant control
               !
               unit = free_unit()
               open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
               read(unit) ik_spex,Nspin_spex,Nkpt_irred_spex,ib_sigma1,ib_sigma2,ldum,lexonly,Nfreq
               allocate(wread(Nfreq));wread=0d0
               read(unit) wread
               close(unit)
               if(lexonly) stop "lexonly"
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
               if(ib_sigma1.gt.ib_Uwan1) stop "ib_sigma1>ib_Uwan1"
               if(ib_sigma2.lt.ib_Uwan2) stop "ib_sigma2<ib_Uwan2"
               if(paramagneticSPEX.and.(Nspin_spex.ne.1)) stop "Spex self-energy file is not paramagnetic."
               if(ik.ne.ik_spex) stop "K-point index in SPEX not match the expected index."
               if(ik.gt.1)then
                  if(ib_sigma1_old.ne.ib_sigma1) stop "ib_sigma1 does not match with previous file."
                  if(Nspin_spex_old.ne.Nspin_spex) stop "Nspin_spex does not match with previous file."
                  if(Nkpt_irred_spex_old.ne.Nkpt_irred_spex) stop "Nkpt_irred_spex does not match with previous file."
                  if(ib_sigma2_old.ne.ib_sigma2) stop "ib_sigma2 does not match with previous file."
                  if(iseg.eq.iseg_old)then
                     if(Nfreq_old.ne.Nfreq) stop "Nfreq does not match among different k-points same segment."
                     if(wS_old.ne.wread(1)) stop "First freq does not match among different k-points same segment."
                     if(wE_old.ne.wread(Nfreq)) stop "Last freq does not match among different k-points same segment."
                  else
                     if(dabs(wread(1)-wE_old).gt.eps)then
                        write(*,"(A,2F10.5)") "w_old, w(1):",wE_old,wread(1)
                        write(*,"(A,2I5)") "iq,iseg:",ik,iseg
                        stop "Segments in sigma do not match."
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
         write(*,"(A,2I4)") "The band indexes in the SPEX self-energy are: ",ib_sigma1,ib_sigma2
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
               if(verbose)write(*,"(A)") "Opening "//reg(path)
               call inquireFile(reg(path),filexists) !redundant control
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
                  if (iseg.eq.1.and.dabs(dimag(SigmaC_seg(1,ib,ik,1))).gt.1.d-6) then
                     write(*,"(A,2E20.12)") "Warning: ImSigmaC_spex("//str(ik)//",1) orb "//str(ib)//" is > 1.d-6: ",SigmaC_seg(1,ib,ik,1)
                  endif
                  if (iseg.eq.SigmaSegments.and.dabs(dimag(SigmaC_seg(NfreqSeg(iseg),ib,ik,1))).gt.1.d-6) then
                     write(*,"(A,2E20.12)") "Warning: ImSigmaC_spex("//str(ik)//",Nw) orb "//str(ib)//" is > 1.d-6: ",SigmaC_seg(NfreqSeg(iseg),ib,ik,1)
                  endif
                  !
                  !Calc Sigma along imag axis using
                  !Sigmac(w)=int dw'Gamma(w')/(w-w')
                  !Gamma= - 1/pi Im(Sigmac) Sign(w-mu)
                  !$OMP PARALLEL DEFAULT(NONE),&
                  !$OMP SHARED(Nspin_spex,ib,ik,iseg,Nmats,NfreqSeg,wread,wmats,SigmaC_seg,SigmaC_diag),&
                  !$OMP PRIVATE(ispin,iw,iw2,gamma1,gamma2,trap)
                  !$OMP DO
                  do ispin=1,Nspin_spex
                     do iw=1,Nmats
                        do iw2=1,NfreqSeg(iseg)-1
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
                        enddo !iw2
                     enddo !iw
                  enddo !ispin
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
         !Sigma=sigmax+sigmac and transform to wannier basis
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(paramagneticSPEX,Nkpt,Norb,ib_Uwan1,ib_Uwan2,Nmats,Lttc,Uwan,SigmaC_diag,SigmaX_seg,Smats_GoWo),&
         !$OMP PRIVATE(ispin,ispin_spex,ik,iwan1,iwan2,iw)
         !$OMP DO
         do ispin=1,Nspin
            do ik=1,Nkpt
               do iwan2=1,Norb
                  do iwan1=1,Norb
                     do iw=1,Nmats
                        !
                        ispin_spex=1
                        if(.not.paramagneticSPEX)ispin_spex=ispin
                        !
                        Smats_GoWo%wks(iwan1,iwan2,iw,ik,ispin) = &
                        + sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * SigmaC_diag(ib_Uwan1:ib_Uwan2,iw,Lttc%kptPos(ik),ispin_spex) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))  &
                        + sum(conjg(Uwan(ib_Uwan1:ib_Uwan2,iwan1,ik,ispin_spex)) * SigmaX_seg(ib_Uwan1:ib_Uwan2,Lttc%kptPos(ik),ispin_spex) * Uwan(ib_Uwan1:ib_Uwan2,iwan2,ik,ispin_spex))
                        !
                        !Smats_GoWo%wks(iwan1,iwan2,iw,ik,ispin) = Smats_GoWo%wks(iwan1,iwan2,iw,ik,ispin) !* H2eV
                        !
                     enddo
                  enddo
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(Uwan,SigmaC_diag,SigmaX_seg)
         !
         Smats_GoWo%wks = Smats_GoWo%wks * H2eV
         !
         call cpu_time(finish)
         write(*,"(A,F)") "Sigma_GoWo(k,w) --> Sigma_GoWo(k,iw) cpu timing:", finish-start
         !
         call FermionicKsum(Smats_GoWo)
         !
         ! Print out the transformed stuff
         ! local
         call dump_FermionicField(Smats_GoWo,1,reg(pathOUTPUT_),"Sigma_GoWo_loc_up.DAT")
         if(Nspin_spex.eq.2)call dump_FermionicField(Smats_GoWo,2,reg(pathOUTPUT_),"Sigma_GoWo_loc_dw.DAT")
         ! k-dependent
         call dump_FermionicField(Smats_GoWo,reg(pathOUTPUT_),"Sigma_GoWo",.true.,Lttc%kpt)
         if(save2readable)call dump_FermionicField(Smats_GoWo,reg(pathOUTPUT_)//"Sigma_imag/","Sigma_GoWo",.false.,Lttc%kpt)
         !
         !
         if(doVxc)call read_Vxc(Vxc,Lttc,ib_sigma1,ib_sigma2,save2readable)
         !
         !---------------------------------------------------------------------!
         !
      else
         !
         !---------------------------------------------------------------------!
         !
         ! Just read all
         call clear_attributes(Smats_GoWo)
         call read_FermionicField(Smats_GoWo,reg(pathINPUT),"Sigma_GoWo")
         !
         if(present(Vxc))then
            call assert_shape(Vxc,[Norb,Norb,Nkpt,Nspin],"read_Sigma_spex","Vxc")
            Vxc=czero
            call read_matrix(Vxc(:,:,:,1),reg(pathINPUT),"Vxc_wann_k_up.DAT")
            if(paramagneticSPEX)then
               Vxc(:,:,:,2) = Vxc(:,:,:,1)
            else
               call read_matrix(Vxc(:,:,:,2),reg(pathINPUT),"Vxc_wann_k_dw.DAT")
            endif
         endif
         !
         !
      endif
      !
   end subroutine read_Sigma_spex


   !---------------------------------------------------------------------------!
   !PURPOSE: Read the exchange potential from gwa file.
   !TEST ON: 27-10-2020
   !---------------------------------------------------------------------------!
   subroutine read_Vxc(Vxc,Lttc,ib_sigma1,ib_sigma2,save2readable)
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
      if(verbose)write(*,"(A)") "---- read_Vxc"
      !
      !
      ! Check on the input Vxc matrix
      Norb = size(Vxc,dim=1)
      if(Norb.ne.size(Vxc,dim=2)) stop "Vxc is not a square matrix."
      Nkpt = size(Vxc,dim=3)
      if(.not.Lttc%status) stop "Lattice container not properly initialized."
      if(Lttc%Nkpt.ne.Nkpt) stop "Lattice has different number of k-points with respect to Vxc."
      !
      ! Read XEPS data
      if(.not.XEPSisread)then
         path = reg(pathINPUT)//"XEPS.DAT"
         call inquireFile(reg(path),filexists)
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
      call inquireFile(reg(path),filexists)
      write(*,"(A)") "Opening "//reg(path)
      unit = free_unit()
      open(unit,file=reg(path),form="unformatted",action="read",position="rewind")
      read(unit) Nspin_Uwan,Nkpt_Uwan,ib_Uwan1,ib_Uwan2,Norb_Uwan
      if(paramagneticSPEX.and.(Nspin_Uwan.ne.1)) stop "UWAN file is not paramagnetic."
      if(Nkpt_Uwan.ne.Nkpt) stop "UWAN file has wrong number of k-points (not irreducible)."
      if(Norb_Uwan.ne.Norb) stop "UWAN file has wrong orbital dimension."
      write(*,"(A,2I4)") "The band indexes in the UWAN rotation are: ",ib_Uwan1,ib_Uwan2
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
      call inquireFile(reg(path),filexists)
      write(*,"(A)") "Opening "//reg(path)
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
         call inquireFile(reg(path),filexists)
         write(*,"(A)") "Opening "//reg(path)
         unit_dis = free_unit()
         open(unit_dis,file=reg(path),form="unformatted",action="read",position="rewind")
         read(unit_dis) Nspin_disent,Nkpt_irred_disent,Norb_disent,ib_Dwan1,ib_Dwan2
         if (Nspin_disent.ne.Nspin_Uwan) stop 'DISENT_EVEC.DAT: nspin'
         if (Nkpt_irred_disent.ne.Lttc%Nkpt_irred) stop 'DISENT_EVEC.DAT: nkpt1'
         if (Norb_disent.ne.Norb) stop 'DISENT_EVEC.DAT: nwan'
         if (ib_Dwan1.ne.ib_Uwan1) stop 'DISENT_EVEC.DAT: ib_wan1'
         if (ib_Dwan2.ne.ib_Uwan2) then
            write(*,"(A2I4)") 'ib_Uwan2,ib_Dwan2',ib_Uwan2, ib_Dwan2
         endif
         allocate(dis_evec(ib_Uwan1:ib_Dwan2,ib_Uwan1:ib_Dwan2));dis_evec=czero
         allocate(cmat(ib_Uwan1:ib_Dwan2,ib_Uwan2:ib_Dwan2));cmat=czero
      endif
      !
      ! Read eig and vxcfull files
      path = reg(pathINPUT)//"eig.DAT"
      call inquireFile(reg(path),filexists)
      write(*,"(A)") "Opening "//reg(path)
      unit_eig = free_unit()
      open(unit_eig,file=reg(path),form='unformatted',access='direct',action='read',recl=irecl)
      !
      path = reg(pathINPUT)//"vxcfull.DAT"
      call inquireFile(reg(path),filexists)
      write(*,"(A)") "Opening "//reg(path)
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
               write(*,*) 'ispin,ik,neigd,nband=',ispin,ik,neigd,nband
               stop 'read_vxc:nband'
            endif
            read(unit_vxc) ((vxcmat(i,j),i=1,j),j=1,nband)
            do j=1,nband
               do i=1,j-1
                  vxcmat(j,i)=conjg(vxcmat(i,j))
               enddo
            enddo
            !do i=1,nbandmax
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
      ! Transform to wannier basis
      !Sigma=sigmax+sigmac and transform to wannier basis
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
         call dump_matrix(Vxc(:,:,:,ispin),reg(pathINPUT),"Vxc_wann",.true.,ispin=ispin)
         if(save2readable)call dump_matrix(Vxc(:,:,:,ispin),reg(pathINPUT)//"Vxc_wann/","Vxc_wann",.false.,ispin=ispin)
      enddo
      !
      allocate(Vxc_loc(Norb,Norb,Nspin));Vxc_loc=czero
      do ik=1,Nkpt
         Vxc_loc = Vxc_loc + Vxc(:,:,ik,:)/Nkpt
      enddo
      call dump_matrix(Vxc_loc(:,:,1),reg(pathINPUT)//"Vxc_wann_loc_up.DAT")
      if(Nspin_Uwan.eq.2)call dump_matrix(Vxc_loc(:,:,2),reg(pathINPUT)//"Vxc_wann_loc_dw.DAT")
      deallocate(Vxc_loc)
      !
    end subroutine read_vxc

end module self_energy
