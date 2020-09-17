module self_energy

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: calc_sigmaGW
   public :: calc_sigmaGW_DC
   public :: read_Sigma_spex

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
      complex(8),allocatable                :: WmatsC(:,:,:,:)
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
      allocate(Gitau(Norb,Norb,Ntau,Nkpt,Nspin));Gitau=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%wks(:,:,:,:,ispin),Gitau(:,:,:,:,ispin), &
         asympt_corr=.true.,tau_uniform=tau_uniform_,Nkpt3=Lttc%Nkpt3,kpt=Lttc%kpt)
      enddo
      !
      allocate(Witau(Nbp,Nbp,Ntau,Nkpt));Witau=czero
      allocate(WmatsC(Nbp,Nbp,Nmats,Nkpt));WmatsC=czero
      do iw=1,Nmats
         WmatsC(:,:,iw,:) = Wmats%screened(:,:,iw,:) - Wmats%bare
      enddo
      call Bmats2itau(Beta,WmatsC,Witau,asympt_corr=.true.,tau_uniform=tau_uniform_)
      deallocate(WmatsC)
      !
      write(*,"(A)") "Sigma_C(tau)"
      !Sigma_{m,n}(q,tau) = -Sum_{k,mp,np} W_{(m,mp);(n,np)}(q-k;tau)G_{mp,np}(k,tau)
      allocate(Sitau(Norb,Norb,Ntau,Lttc%Nkpt_irred,Nspin));Sitau=czero
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
      deallocate(Witau)
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
      deallocate(Gitau)
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
   !PURPOSE: Compute the two GW self-energy components.
   !---------------------------------------------------------------------------!
   subroutine calc_sigmaGW_DC(Smats_Cdc,Smats_Xdc,Gmats,Wmats,tau_uniform)
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
      type(FermionicField),intent(inout)    :: Smats_Cdc
      type(FermionicField),intent(inout)    :: Smats_Xdc
      type(FermionicField),intent(in)       :: Gmats
      type(BosonicField),intent(in)         :: Wmats
      logical,intent(in),optional           :: tau_uniform
      !
      complex(8),allocatable                :: Sitau_loc(:,:,:,:)
      complex(8),allocatable                :: Gitau_loc(:,:,:,:)
      complex(8),allocatable                :: Witau_loc(:,:,:)
      complex(8),allocatable                :: WmatsC_loc(:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Norb,Nmats
      integer                               :: iw,itau,ispin
      integer                               :: m,n,mp,np,ib1,ib2
      logical                               :: tau_uniform_
      !
      !
      write(*,"(A)") "--- calc_sigmaGW_DC ---"
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
      tau_uniform_=.false.
      if(present(tau_uniform)) tau_uniform_ = tau_uniform
      !
      allocate(Gitau_loc(Norb,Norb,Ntau,Nspin));Gitau_loc=czero
      do ispin=1,Nspin
         call Fmats2itau_mat(Beta,Gmats%ws(:,:,:,ispin),Gitau_loc(:,:,:,ispin), &
                             asympt_corr=.true.,tau_uniform=tau_uniform_)
      enddo
      !
      allocate(Witau_loc(Nbp,Nbp,Ntau));Witau_loc=czero
      allocate(WmatsC_loc(Nbp,Nbp,Nmats));WmatsC_loc=czero
      do iw=1,Nmats
         WmatsC_loc(:,:,iw) = Wmats%screened_local(:,:,iw) - Wmats%bare_local
      enddo
      call Bmats2itau(Beta,WmatsC_loc,Witau_loc,asympt_corr=.true.,tau_uniform=tau_uniform_)
      deallocate(WmatsC_loc)
      !
      write(*,"(A)") "Sigma_Cdc(tau)"
      !Sigma_{m,n}(q,tau) = -Sum_{k,mp,np} W_{(m,mp);(n,np)}(q-k;tau)G_{mp,np}(k,tau)
      allocate(Sitau_loc(Norb,Norb,Ntau,Nspin));Sitau_loc=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,Ntau,Sitau_loc,Gitau_loc,Witau_loc),&
      !$OMP PRIVATE(m,n,itau,ispin,mp,np,ib1,ib2)
      !$OMP DO
      do m=1,Norb
         do n=1,Norb
            do itau=1,Ntau
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
      write(*,"(A)") "Sigma_Cdc(tau)-->Sigma_Cdc(iw)"
      call clear_attributes(Smats_Cdc)
      do ispin=1,Nspin
         call Fitau2mats_mat(Beta,Sitau_loc(:,:,:,ispin),Smats_Cdc%ws(:,:,:,ispin),tau_uniform=tau_uniform_)
      enddo
      deallocate(Sitau_loc)
      !
      write(*,"(A)") "Sigma_X"
      call clear_attributes(Smats_Xdc)
      !sigmax(r,r')=-g(r,r',tau=0-)*v(r-r')
      !Sigmax_nm(q) = Sum_kij V_{ni,jm}(q-k)G_ij(k,beta)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Norb,Ntau,Smats_Xdc,Gitau_loc,Wmats),&
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
                     Smats_Xdc%s(m,n,ispin) = Smats_Xdc%s(m,n,ispin) + Gitau_loc(mp,np,Ntau,ispin)*Wmats%bare_local(ib1,ib2)
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
      !
   end subroutine calc_sigmaGW_DC


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



   !---------------------------------------------------------------------------!
   !PURPOSE: Read self-energy from SPEX files.
   !---------------------------------------------------------------------------!
   subroutine read_Sigma_spex(Smats,Lttc,save2bin,pathOUTPUT,doAC)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use global_vars, only :  pathINPUT,UseXepsKorder
      implicit none
      !
      type(BosonicField),intent(inout)      :: Smats
      type(Lattice),intent(in)              :: Lttc
      logical,intent(in)                    :: save2bin
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: doAC
      !
      logical                               :: filexists,ACdone,doAC_
      logical                               :: UseDisentangledBS
      character(len=256)                    :: file_spex,path,pathOUTPUT_
      integer                               :: iseg,SigmaSegments
      integer                               :: idum1,idum2,idum3,idum4
      logical                               :: ldum,lexonly





      integer                               :: unit,Nkpt
      integer                               :: iq,iw,iqread,Nbp_spex

      integer                               :: Nspin_spex,Norb_spex,Nfreq
      integer                               :: ib1,ib2,iw1,iw2
      real(8),allocatable                   :: wread(:),wmats(:)
      complex(8),allocatable                :: D1(:,:),D2(:,:),D3(:,:),imgFact(:,:,:)
      complex(8),allocatable                :: Utmp(:,:)
      type(BosonicField)                    :: Ureal
      real                                  :: start,finish
      !
      !
      write(*,"(A)") "--- read_Sigma_spex ---"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Check on the input Boson
      if(.not.Smats%status) stop "FermionicField not properly initialized."
      if(.not.Lttc%status) stop "Lattice container not properly initialized."
      if(Lttc%mu.eq.0d0) stop "Chemical potwential in Lattice container not properly initialized (0d0)."
      allocate(wmats(Smats%Npoints));wmats=0d0
      wmats = FermionicFreqMesh(Smats%Beta,Smats%Npoints)
      !
      !
      ! Read XEPS data
      path = pathINPUT//"XEPS.DAT"
      call inquireFile(reg(path),filexists)
      allocate(kptPos_xeps(size(Lttc%kpt,dim=2)));kptPos_xeps=0
      call read_xeps(reg(path),Lttc%kpt,Lttc%Nkpt3,UseXepsKorder,Lttc%kptPos,Lttc%Nkpt_irred,UseDisentangledBS)
      !
      !
      ! Read UWAN file



      !
      !
      ! Check if the data on the Matsubara axis are present
      path = pathINPUT//"Sigma_imag" !/SIGMA.Q0001.DAT"
      call inquireDir(reg(path),ACdone,hardstop=.false.)
      doAC_ = .not.ACdone
      if(present(doAC)) doAC_ = doAC
      !
      !
      ! Perform cnalytical continuation on the self-energy on the real axis
      if(doAC_) then
         !
         !---------------------------------------------------------------------!
         !
         ! Look for the Number of Sigma segments.
         SigmaSegments=0
         do iseg=1,99
            file_spex = reg(path)//"Sigma_real_"//str(iq,1)
            call inquireDir(reg(file_spex),filexists,hardstop=.false.)
            if(.not.filexists) exit
            SigmaSegments = SigmaSegments + 1
         enddo
         write(*,"(A,1I6)") "The number of segment of the SPEX self-energy is: ",SigmaSegments
         allocate(NfreqSeg(SigmaSegments));NfreqSeg=0
         !
         ! Look for the Number of Kpoints in each segment (supposed to be the same of Lttc%Nkpt_irred)
         do iseg=1,SigmaSegments
            Nkpt = 0
            do ik=1,2000
               file_spex = reg(path)//"Sigma_real_"//str(iq,1)//"/SIGMA.Q"//str(iq,4)//".DAT"
               call inquireFile(reg(file_spex),filexists,hardstop=.false.)
               if(.not.filexists) exit
               Nkpt = Nkpt + 1
            enddo
            write(*,"(A,1I4,A,1I6)") "The number k-points in the segment number: ",iseg," is: ",Nkpt
            Nkpt_old = Nkpt
            if((iseg.gt.1).and.(Nkpt.ne.Nkpt_old)) stop "Number of K-points does not match among segments."
            if(Nkpt.ne.Lttc%Nkpt_irred) stop "Number of K-points does not match with Nkpt_irred readed from XEPS."
         enddo
         !
         ! Check that all the Sigma parameters are cosistent and that the segments match
         do iseg=1,SigmaSegments
            do ik=1,Lttc%Nkpt_irred
               !
               file_spex = reg(path)//"Sigma_real_"//str(iseg,1)//"/SIGMA.Q"//str(ik,4)//".DAT"
               write(*,"(A)") "Checking "//reg(file_spex)
               call inquireFile(reg(file_spex),filexists) !redundant control
               !
               unit = free_unit()
               open(unit,file=reg(path),form="unformatted",action="read")
               read(unit) ik_spex,Nspin_spex,Nkpt_irred_spex,ib_sigma1,ib_sigma2,ldum,lexonly,Nfreq
               allocate(wread(Nfreq));wred=0d0
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
               ! Eack k-point controls
               if(Nspin_spex.ne.1) stop "Nspin_spex is expected to be 1."
               if(ik.ne.ik_spex) stop "ik_spex mismatch."
               if(ik.gt.1)then
                  if(ib_sigma1_old.ne.ib_sigma1) stop "ib_sigma1 does not match with previous file."
                  if(Nspin_spex_old.ne.Nspin_spex) stop "Nspin_spex does not match with previous file."
                  if(Nkpt_irred_spex_old.ne.Nkpt_irred_spex) stop "Nkpt_irred_spex does not match with previous file."
                  if(ib_sigma2_old.ne.ib_sigma2) stop "ib_sigma2 does not match with previous file."
                  if(iseg.eq.iseg_old)then
                     if(Nfreq_old.ne.Nfreq) stop "Nfreq does not match among different k-points same segment."
                     if(wS_old.ne.wread(1))) stop "First freq does not match among different k-points same segment."
                     if(wE_old.ne.wread(Nfreq))) stop "Last freq does not match among different k-points same segment."
                  else
                     if(dabs(wread(1)-wE_old).gt.eps)then
                        write(*,"(A,2F10.5)") "w_old, w(1):",w_old,w(1)
                        write(*,"(A,2I5)") "iq,iseg:",ik,iseg
                        stop "Segments in sigma do not match."
                    endif
               endif
               deallocate(wread)
               !
            enddo
            NfreqSeg(iseg) = Nfreq
         enddo
         write(*,"(A,2I4)") "The band indexes in the SPEX self-energy are: ",ib_sigma1,ib_sigma2
         !
         do iseg=1,SigmaSegments
            !
            allocate(SigmaX_GW(ib_sigma1:ib_sigma2,Lttc%Nkpt_irred,Nspin_spex))
            allocate(SigmaC_GW(NfreqSeg(iseg),ib_sigma1:ib_sigma2,Lttc%Nkpt_irred,Nspin_spex))
            !
            allocate(wread(NfreqSeg(iseg)));wread=0d0
            allocate(SigmaX_tmp(3,ib_sigma1:ib_sigma2))
            allocate(SigmaC_tmp(2,NfreqSeg(iseg),ib_sigma1:ib_sigma2))
            do iq=1,Lttc%Nkpt_irred
               !
               file_spex = reg(path)//"Sigma_real_"//str(iq,1)//"/SIGMA.Q"//str(iq,4)//".DAT"
               write(*,"(A)") "Opening "//reg(file_spex)
               call inquireFile(reg(file_spex),filexists) !redundant control
               !
               unit = free_unit()
               open(unit,file=reg(path),form="unformatted",action="read")
               read(unit) idum1,idum2,idum3,ib_sigma1,ib_sigma2,ldum,lexonly,Nfreq
               read(unit) wread
               !
               do ispin=1,Nspin_spex
                  !
                  ! Every iq file contains all the different k component of the self-energy
                  ! so that here ik=1 is the sum over all the interanl momentum q
                  do ik=1,Lttc%Nkpt_irred
                     !
                     SigmaX_tmp=0d0;SigmaC_tmp=czero
                     read(ifile) SigmaX_tmp
                     read(ifile) SigmaC_tmp
                     !
                     do ib=ib_sigma1,ib_sigma2
                        SigmaX_GW(ib,ik,ispin)=SigmaX_GW(ib,ik,ispin)+sum(SigmaX_tmp(:,ib))
                        do iw=1,NfreqSeg(iseg)
                           SigmaC_GW(iw,ib,ik,ispin)=SigmaC_GW(iw,ib,ik,ispin)+sum(SigmaC_tmp(:,iw,ib))
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


























         !
         ! Allocate the temporary quantities needed by the Analytical continuation
         allocate(D1(Nbp_spex,Nbp_spex));D1=czero
         allocate(D2(Nbp_spex,Nbp_spex));D2=czero
         allocate(D3(Nbp_spex,Nbp_spex));D3=czero
         !
         !
         ! Analytical continuation of the local component to imag axis using spectral rep
         call cpu_time(start)
         !
         ! Check if any local Urpa component has inverted Im/Re symmetry
         do ib1=1,Nbp_spex
            do ib2=1,Nbp_spex
               if( abs(real(Ureal%bare_local(ib1,ib2))).lt.abs(aimag(Ureal%bare_local(ib1,ib2))))then
                  write(*,"(A,2I5)")"Element: ",ib1,ib2
                  write(*,"(A,E14.7)")"Re[Ubare(w=inf)]: ",real(Ureal%bare_local(ib1,ib2))
                  write(*,"(A,E14.7)")"Im[Ubare(w=inf)]: ",aimag(Ureal%bare_local(ib1,ib2))
                  stop "Something wrong: Uloc cannot have inverted Re/Im parity."
               endif
            enddo
         enddo
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Nbp_spex,wmats,wread,Nfreq,Ureal,Umats),&
         !$OMP PRIVATE(ib1,ib2,iw1,iw2,D1,D2,D3,Utmp)
         !$OMP DO
         do iw1=1,Umats%Npoints
            Utmp=czero
            do iw2=1,Nfreq-2,2
               !
               do ib1=1,Nbp_spex
                  do ib2=1,Nbp_spex
                     D1(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2)   )/pi
                     D2(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+1) )/pi
                     D3(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+2) )/pi
                  enddo
               enddo
               !
               !D(-w)=-D(w), integrate using Simpson method
               if(wread(iw2).gt.0.d0) then
                  Utmp(:,:) = Utmp(:,:) + ( D1(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2)  ) - D1(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2)  ) ) *(wread(iw2+1)-wread(iw2))/3.d0
                  Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                  Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
               elseif(dabs(wread(iw2)).lt.1.d-12) then
                  Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                  Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
               endif
            enddo
            !
            do ib1=1,Nbp_spex
               do ib2=1,Nbp_spex
                  Umats%screened_local(ib1,ib2,iw1) = Utmp(ib1,ib2) + Umats%bare_local(ib1,ib2)
               enddo
            enddo
            !
         enddo !iw1
         !
         !$OMP END DO
         !$OMP END PARALLEL
         call cpu_time(finish)
         deallocate(D1,D2,D3)
         write(*,"(A)") "UcRPA(w) --> UcRPA(iw) cpu timing:", finish-start
         !
         !
         if(.not.LocalOnly)then
            !
            ! Allocate the temporary quantities needed by the Analytical continuation
            allocate(D1(Nbp_spex,Nbp_spex));D1=czero
            allocate(D2(Nbp_spex,Nbp_spex));D2=czero
            allocate(D3(Nbp_spex,Nbp_spex));D3=czero
            !
            !
            ! Analytical continuation of all the K-points to imag axis using spectral rep
            allocate(imgFact(Nbp_spex,Nbp_spex,2));imgFact=cone
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(Nbp_spex,wmats,wread,Nfreq,Ureal,Umats,UfullStructure),&
            !$OMP PRIVATE(iq,ib1,ib2,iw1,iw2,D1,D2,D3,Utmp,imgFact)
            !$OMP DO
            do iq=1,Umats%Nkpt
               !
               ! Some elelments of U, usually the k-dependent one, have inverted Im/Re symmetry
               imgFact=cone
               if(UfullStructure)then
                  do ib1=1,Nbp_spex
                     do ib2=1,Nbp_spex
                        if( abs(real(Ureal%bare(ib1,ib2,iq))).lt.abs(aimag(Ureal%bare(ib1,ib2,iq))))then
                           imgFact(ib1,ib2,1) = -img !this correspond to dividing by I
                           imgFact(ib1,ib2,2) = +img !this correspond to multiplying by I
                        endif
                     enddo
                  enddo
               endif
               !
               do iw1=1,Umats%Npoints
                  Utmp=czero
                  do iw2=1,Nfreq-2,2
                     !
                     do ib1=1,Nbp_spex
                        do ib2=1,Nbp_spex
                           D1(ib1,ib2) = -dimag( ( imgFact(ib1,ib2,1) * Ureal%screened(ib1,ib2,iw2,iq)   ) )/pi
                           D2(ib1,ib2) = -dimag( ( imgFact(ib1,ib2,1) * Ureal%screened(ib1,ib2,iw2+1,iq) ) )/pi
                           D3(ib1,ib2) = -dimag( ( imgFact(ib1,ib2,1) * Ureal%screened(ib1,ib2,iw2+2,iq) ) )/pi
                        enddo
                     enddo
                     !
                     !D(-w)=-D(w), integrate using Simpson method
                     if(wread(iw2).gt.0.d0) then
                        Utmp(:,:) = Utmp(:,:) + ( D1(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2)  ) - D1(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2)  ) ) *(wread(iw2+1)-wread(iw2))/3.d0
                        Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                        Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
                     elseif(dabs(wread(iw2)).lt.1.d-12) then
                        Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                        Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
                     endif
                  enddo
                  !
                  do ib1=1,Nbp_spex
                     do ib2=1,Nbp_spex
                        Umats%screened(ib1,ib2,iw1,iq) = imgFact(ib1,ib2,2)*Utmp(ib1,ib2) + Umats%bare(ib1,ib2,iq)
                     enddo
                  enddo
                  !
               enddo !iw1
               !
            enddo !iq
            !
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            deallocate(D1,D2,D3)
            write(*,"(A)") "UcRPA(q,w) --> UcRPA(q,iw) cpu timing:", finish-start
            !
         endif !LocalOnly
         call checkAnalyticContinuation(Umats,Ureal)
         deallocate(Utmp)
         !
         ! Print out the transformed stuff - local
         call dump_BosonicField(Umats,reg(pathOUTPUT_),"Uloc_mats.DAT")
         call dump_BosonicField(Ureal,reg(pathOUTPUT_),"Uloc_real.DAT",wread)
         !
         ! Print out the transformed stuff - Kdep
         call dump_BosonicField(Umats,reg(pathOUTPUT_//"VW_imag/"),.true.)
         if(.not.save2bin)then
            call dump_BosonicField(Umats,reg(pathOUTPUT_//"VW_imag_readable/"),save2bin)
            call dump_BosonicField(Ureal,reg(pathOUTPUT_//"VW_real_readable/"),save2bin,axis=wread)
         endif
         !
         deallocate(wread)
         call DeallocateBosonicField(Ureal)
         !
         !---------------------------------------------------------------------!
         !
      else
         !
         !---------------------------------------------------------------------!
         !
         ! Allocations from dimensions written in W.Q0001.DAT file
         path = pathINPUT//"VW_imag/VW.Q0001.DAT"
         call inquireFile(reg(path),filexists)
         !
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read")
         read(unit)idum,Nspin_spex,Norb_spex,Nfreq
         close(unit)
         !
         Nbp_spex = Norb_spex**2
         allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         write(*,"(A,I5)")"Matsubara frequencies: ",Nfreq
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Nbp_spex) stop "Size of given BosonicField and VW_imag orbital space do not coincide."
         if(Umats%Npoints.ne.Nfreq) stop "Number of VW_imag Matsubara points and bosonic field mesh does not coincide."
         !
         ! Look for the Number of SPEX files. Which are supposed to be ordered.
         Nkpt = 0
         do iq=1,2000
            file_spex = reg(path)//"VW_imag/VW.Q"//str(iq)//".DAT"
            call inquireFile(reg(file_spex),filexists,hardstop=.false.)
            if(.not.filexists) exit
            Nkpt = Nkpt + 1
         enddo
         write(*,"(A1,1I6)") "The number of SPEX files (Nkpt) in VW_imag is: ",Nkpt
         if((.not.LocalOnly).and.(Umats%Nkpt.ne.Nkpt)) stop "Number of k-points of given BosonicField and number of VW_imag k-points do not coincide."
         !
         ! Read VW_imag accumulating local attribute and optionally storing the k-dependent part
         path = pathINPUT//"VW_imag/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"VW_imag/VW.Q"//str(iq)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"VW_imag/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),filexists) !redundant control
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
            write(*,"(A)")"read iq",iq   !!!!>>>>>TEST<<<<<!!!!
            if (iq.ne.iqread) stop "iqread.ne.iq"
            !
            read(unit) wread
            wread = H2eV*wread
            write(*,"(A)")"read wread",wread   !!!!>>>>>TEST<<<<<!!!!
            do iw=1,Nfreq
               if (dabs(wread(iw)-wmats(iw)).gt.eps) stop "wread.ne.wmats"
            enddo
            !
            do iw=0,Nfreq
               read(unit) Utmp
               if(iw.eq.0) then
                  Umats%bare_local = Umats%bare_local + H2eV*Utmp/(Nkpt**3)
                  if(.not.LocalOnly) Umats%bare(:,:,iq) = H2eV*Utmp/(Nkpt**2)
               else
                  Umats%screened_local(:,:,iw) = Umats%screened_local(:,:,iw) + H2eV*Utmp/(Nkpt**3)
                  if(.not.LocalOnly) Umats%screened(:,:,iw,iq) = H2eV*Utmp/(Nkpt**2)
               endif
            enddo
            !
            close(unit)
            !
         enddo !iq
         !
         deallocate(wread,Utmp)
         !
         !---------------------------------------------------------------------!
         !
      endif
      !
      ! Remove elements with inverted parity from the k-dependent fields.
      if(allocated(Umats%screened).and.allocated(Umats%bare).and.(.not.UfullStructure))then
         do iq=1,Nkpt
            do ib1=1,Nbp_spex
               do ib2=1,Nbp_spex
                  if (dabs(dimag(Umats%bare(ib1,ib2,iq))).gt.1.d-6) then
                     write(*,"(A,2I5)") "Warning Umats%bare imaginary. Set matrix element to static value",ib1,ib2
                     Umats%bare(ib1,ib2,iq) = Umats%screened(ib1,ib2,1,iq)
                     Umats%screened(ib1,ib2,:,iq) = Umats%screened(ib1,ib2,1,iq)
                  endif
               enddo
            enddo
         enddo
      endif
      !
   end subroutine read_Sigma_spex



end module self_energy
