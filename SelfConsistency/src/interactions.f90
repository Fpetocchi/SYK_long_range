module interactions

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface read_U_spex
      module procedure read_U_spex_full                                         ![BosonicField,LocalOnly,save2readable,pathOUTPUT(optional to change output path),doAC(optional to override AC)]
      module procedure read_U_spex_Uloc0                                        ![Matrix,pathOUTPUT(optional to change output path)]
   end interface read_U_spex

   interface build_Umat
      module procedure build_Umat_singlParam                  ! (Solver Format) ![Matrix,Uaa_screened,Uab_screened,J_screened]
      module procedure build_Umat_multiParam                  ! (Solver Format) ![Matrix,Vector,Matrix,Matrix]
   end interface build_Umat

   interface build_Uret
      module procedure build_Uretloc_singlParam                   !   (GW Format)  ![BosonicField,Uaa_bare,Uab_bare,J_bare,vector_g,vector_w0]
      module procedure build_Uretloc_multiParam                   !   (GW Format)  ![BosonicField,Vector,Matrix,Matrix,vector_g,vector_w0]
   end interface build_Uret

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif
   complex(8),allocatable,private           :: Ugamma(:,:,:)

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: init_Uelements
   public :: calc_W_full
   public :: calc_W_edmft
   public :: calc_chi_full
   public :: calc_chi_edmft
   public :: read_U_spex
   public :: build_Umat
   public :: build_Uret
   public :: calc_QMCinteractions
   public :: calc_curlyU
   public :: calc_Wimp
   !public :: rescale_interaction

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Indentify the Tensor indexes which correspond to physical (number
   !         and spin conserving) interaction elements
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine init_Uelements(Norb,Uelements)
      !
      use parameters
      implicit none
      !
      integer,intent(in)                    :: Norb
      type(physicalU),intent(inout)         :: Uelements
      !
      integer                               :: Nflavor
      integer                               :: ib1,ib2
      integer                               :: iorb,jorb,korb,lorb
      integer                               :: ispin,jspin
      !
      !
      if(verbose)write(*,"(A)") "---- init_Uelements"
      !
      !
      if(Uelements%status) write(*,"(A)") "Warning: the Physical interaction elements container is being reinitialized."
      Nflavor = Norb*Nspin
      !
      ! Elements when the interaction is in the Norb*Nspin form
      Uelements%Flav_Size = Norb
      if(allocated(Uelements%Flav_Uloc))deallocate(Uelements%Flav_Uloc)
      if(allocated(Uelements%Flav_U1st))deallocate(Uelements%Flav_U1st)
      if(allocated(Uelements%Flav_U2nd))deallocate(Uelements%Flav_U2nd)
      if(allocated(Uelements%Flav_All)) deallocate(Uelements%Flav_All)
      if(allocated(Uelements%Flav_Map)) deallocate(Uelements%Flav_Map)
      allocate(Uelements%Flav_Uloc(Nflavor,Nflavor)) ;Uelements%Flav_Uloc=.false.
      allocate(Uelements%Flav_U1st(Nflavor,Nflavor)) ;Uelements%Flav_U1st=.false.
      allocate(Uelements%Flav_U2nd(Nflavor,Nflavor)) ;Uelements%Flav_U2nd=.false.
      allocate(Uelements%Flav_All(Nflavor,Nflavor))  ;Uelements%Flav_All =.false.
      allocate(Uelements%Flav_Map(Nflavor,Nflavor,4));Uelements%Flav_Map=0
      !
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            iorb = (ib1+mod(ib1,2))/2
            jorb = (ib2+mod(ib2,2))/2
            ispin = abs(mod(ib1,2)-2)
            jspin = abs(mod(ib2,2)-2)
            !
            Uelements%Flav_Uloc(ib1,ib2) = (iorb.eq.jorb).and.(ispin.ne.jspin)
            Uelements%Flav_U1st(ib1,ib2) = (iorb.ne.jorb).and.(ispin.ne.jspin)
            Uelements%Flav_U2nd(ib1,ib2) = (iorb.ne.jorb).and.(ispin.eq.jspin)
            !
            Uelements%Flav_All(ib1,ib2) = Uelements%Flav_Uloc(ib1,ib2) .or.  &
                                          Uelements%Flav_U1st(ib1,ib2) .or.  &
                                          Uelements%Flav_U2nd(ib1,ib2)
            !
            Uelements%Flav_Map(ib1,ib2,1) = iorb
            Uelements%Flav_Map(ib1,ib2,2) = jorb
            Uelements%Flav_Map(ib1,ib2,3) = ispin
            Uelements%Flav_Map(ib1,ib2,4) = jspin
            !
         enddo
      enddo
      !
      ! Elements when the interaction is in the Norb^2 form
      Uelements%Flav_Size = Norb*Norb
      if(allocated(Uelements%Full_Uaa))deallocate(Uelements%Full_Uaa)
      if(allocated(Uelements%Full_Uab))deallocate(Uelements%Full_Uab)
      if(allocated(Uelements%Full_Jsf))deallocate(Uelements%Full_Jsf)
      if(allocated(Uelements%Full_Jph))deallocate(Uelements%Full_Jph)
      if(allocated(Uelements%Full_All))deallocate(Uelements%Full_All)
      if(allocated(Uelements%Full_Map))deallocate(Uelements%Full_Map)
      allocate(Uelements%Full_Uaa(Norb*Norb,Norb*Norb))  ;Uelements%Full_Uaa=.false.
      allocate(Uelements%Full_Uab(Norb*Norb,Norb*Norb))  ;Uelements%Full_Uab=.false.
      allocate(Uelements%Full_Jsf(Norb*Norb,Norb*Norb))  ;Uelements%Full_Jsf=.false.
      allocate(Uelements%Full_Jph(Norb*Norb,Norb*Norb))  ;Uelements%Full_Jph=.false.
      allocate(Uelements%Full_All(Norb*Norb,Norb*Norb))  ;Uelements%Full_All=.false.
      allocate(Uelements%Full_Map(Norb*Norb,Norb*Norb,4));Uelements%Full_Map=0
      !
      do iorb=1,Norb
         do jorb=1,Norb
            do korb=1,Norb
               do lorb=1,Norb
                  !
                  ib1 = iorb + Norb*(jorb-1)
                  ib2 = korb + Norb*(lorb-1)
                  !
                  Uelements%Full_Uaa(ib1,ib2) = (iorb.eq.jorb).and.(korb.eq.lorb).and.(iorb.eq.korb)
                  Uelements%Full_Uab(ib1,ib2) = (iorb.eq.jorb).and.(korb.eq.lorb).and.(iorb.ne.korb)
                  Uelements%Full_Jsf(ib1,ib2) = (iorb.eq.lorb).and.(jorb.eq.korb).and.(iorb.ne.jorb)
                  Uelements%Full_Jph(ib1,ib2) = (iorb.eq.korb).and.(jorb.eq.lorb).and.(iorb.ne.jorb)
                  !
                  Uelements%Full_All(ib1,ib2) = Uelements%Full_Uaa(ib1,ib2) .or.  &
                                                Uelements%Full_Uab(ib1,ib2) .or.  &
                                                Uelements%Full_Jsf(ib1,ib2) .or.  &
                                                Uelements%Full_Jph(ib1,ib2)
                  !
                  Uelements%Full_Map(ib1,ib2,1) = iorb
                  Uelements%Full_Map(ib1,ib2,2) = jorb
                  Uelements%Full_Map(ib1,ib2,3) = korb
                  Uelements%Full_Map(ib1,ib2,4) = lorb
                  !
               enddo
            enddo
         enddo
      enddo
      !
   end subroutine init_Uelements


   !---------------------------------------------------------------------------!
   !PURPOSE: Lattice inversion to get fully screened interaction - GW+EDMFT
   !TEST ON: 21-10-2020
   !---------------------------------------------------------------------------!
   subroutine calc_W_full(Wmats,Umats,Pmats,Lttc)
      !
      use parameters
      use utils_misc
      use utils_fields
      use linalg, only : zeye, inv  !_her !or sym??
      use input_vars, only : HandleGammaPoint
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wmats
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: invW(:,:)
      complex(8),allocatable                :: den_smallk(:,:,:,:)
      complex(8),allocatable                :: den_smallk_avrg(:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      integer                               :: ismall,num_k
      !
      !
      if(verbose)write(*,"(A)") "---- calc_W_full"
      !
      !
      ! Check on the input Fields
      if(.not.Wmats%status) stop "Wmats not properly initialized."
      if(.not.Umats%status) stop "Umats not properly initialized."
      if(.not.Pmats%status) stop "Pmats not properly initialized."
      if(Wmats%Nkpt.eq.0) stop "Wmats k dependent attributes not properly initialized."
      if(Umats%Nkpt.eq.0) stop "Umats k dependent attributes not properly initialized."
      if(Pmats%Nkpt.eq.0) stop "Pmats k dependent attributes not properly initialized."
      if(Umats%iq_gamma.lt.0) stop "Umats iq_gamma not defined."
      if(.not.allocated(Lttc%small_ik)) stop "Kpoints near Gamma not stored. W divergence cannot be cured."
      !
      Nbp = Wmats%Nbp
      Nkpt = Wmats%Nkpt
      Beta = Wmats%Beta
      Nmats = Wmats%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "Either Umats and/or Pmats have different orbital dimension with respect to Wmats."
      if(all([Umats%Nkpt-Nkpt,Pmats%Nkpt-Nkpt].ne.[0,0])) stop "Either Umats and/or Pmats have different number of k-points with respect to Wmats."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "Either Umats and/or Pmats have different Beta with respect to Wmats."
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))then
         Nmats = minval([Wmats%Npoints,Umats%Npoints,Pmats%Npoints])
         write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller: "//str(Nmats)
      endif
      !
      allocate(invW(Nbp,Nbp));invW=czero
      call clear_attributes(Wmats)
      !
      if(HandleGammaPoint)then
         allocate(den_smallk(Nbp,Nbp,Nmats,12));den_smallk=czero
         allocate(den_smallk_avrg(Nbp,Nbp,Nmats));den_smallk_avrg=czero
      endif
      !
      ! Assuming that the Polarization vanishes at iw-->inf
      Wmats%bare = Umats%bare
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Wmats,den_smallk,Lttc,HandleGammaPoint),&
      !$OMP PRIVATE(iq,iw,invW,ismall)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if((iq.eq.Umats%iq_gamma).and.HandleGammaPoint)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            invW = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened(:,:,iw,iq))
            !
            ! [ 1 - U*Pi ]^-1
            call inv(invW)
            !
            ! [ 1 - U*Pi ]^-1 * U
            Wmats%screened(:,:,iw,iq) = matmul(invW,Umats%screened(:,:,iw,iq))
            !
            !store the dielectric function around the Gamma point
            if(HandleGammaPoint)then
               do ismall=1,12
                  if (Lttc%small_ik(ismall,1).eq.iq) den_smallk(:,:,iw,ismall) = invW
               enddo
            endif
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      ! Gamma point handling
      if(HandleGammaPoint)then
         !
         num_k=1
         do ismall=2,12
            if (Lttc%small_ik(ismall,2).eq.1) num_k=num_k+1
         enddo
         !
         !if num_k only is one include also next nearset points
         if(num_k.eq.1)then
            num_k=1
            do ismall=2,12
               if (Lttc%small_ik(ismall,2).le.2) num_k=num_k+1
            enddo
         endif
         !
         !calc average epsinv for small k
         do ismall=1,num_k
            den_smallk_avrg = den_smallk_avrg + den_smallk(:,:,:,ismall)/num_k
         enddo
         !
         !W at gamma point should be real
         den_smallk_avrg = real(den_smallk_avrg)
         !
         !Fill the Gamma point value - element not included in the iq loop
         do iw=1,Nmats
            Wmats%screened(:,:,iw,Umats%iq_gamma) = matmul(den_smallk_avrg(:,:,iw),Umats%screened(:,:,iw,Umats%iq_gamma))
         enddo
         !
         deallocate(den_smallk,den_smallk_avrg)
         !
      endif
      deallocate(invW)
      !
      ! Fill the local attributes
      call BosonicKsum(Wmats)
      !
   end subroutine calc_W_full


   !---------------------------------------------------------------------------!
   !PURPOSE: Lattice inversion to get fully screened interaction - EDMFT
   !---------------------------------------------------------------------------!
   subroutine calc_W_edmft(Wmats,Umats,Pmats,Lttc)
      !
      use parameters
      use utils_misc
      use utils_fields
      use linalg, only : zeye, inv
      use input_vars, only : HandleGammaPoint
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wmats
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: invW(:,:)
      complex(8),allocatable                :: den_smallk(:,:,:,:)
      complex(8),allocatable                :: den_smallk_avrg(:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      integer                               :: ismall,num_k
      !
      !
      if(verbose)write(*,"(A)") "---- calc_W_edmft"
      !
      !
      ! Check on the input Fields
      if(.not.Wmats%status) stop "Wmats not properly initialized."
      if(.not.Umats%status) stop "Umats not properly initialized."
      if(.not.Pmats%status) stop "Pmats not properly initialized."
      if(Umats%Nkpt.eq.0) stop "Umats k dependent attributes not properly initialized."
      if(Wmats%Nkpt.ne.0) stop "Wmats k dependent attributes are supposed to be unallocated."
      if(Pmats%Nkpt.ne.0) stop "Pmats k dependent attributes are supposed to be unallocated."
      if(Umats%iq_gamma.lt.0) stop "Umats iq_gamma not defined."
      if(.not.allocated(Lttc%small_ik)) stop "Kpoints near Gamma not stored. W divergence cannot be cured."
      !
      Nbp = Wmats%Nbp
      Nkpt = Umats%Nkpt
      Beta = Wmats%Beta
      Nmats = Wmats%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "Either Umats and/or Pmats have different orbital dimension with respect to Wmats."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "Either Umats and/or Pmats have different Beta with respect to Wmats."
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))then
         Nmats = minval([Wmats%Npoints,Umats%Npoints,Pmats%Npoints])
         write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller: "//str(Nmats)
      endif
      !
      allocate(invW(Nbp,Nbp));invW=czero
      call clear_attributes(Wmats)
      !
      if(HandleGammaPoint)then
         allocate(den_smallk(Nbp,Nbp,Nmats,12));den_smallk=czero
         allocate(den_smallk_avrg(Nbp,Nbp,Nmats));den_smallk_avrg=czero
      endif
      !
      ! Assuming that the Polarization vanishes at iw-->inf
      Wmats%bare_local = Umats%bare_local
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Wmats,den_smallk,Lttc,HandleGammaPoint),&
      !$OMP PRIVATE(iq,iw,invW,ismall)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if((iq.eq.Umats%iq_gamma).and.HandleGammaPoint)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            invW = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened_local(:,:,iw))
            !
            ! [ 1 - U*Pi ]^-1
            call inv(invW)
            !
            ! [ 1 - U*Pi ]^-1 * U
            Wmats%screened_local(:,:,iw) = Wmats%screened_local(:,:,iw) + matmul(invW,Umats%screened(:,:,iw,iq))/Nkpt
            !
            !store the dielectric function around the Gamma point
            if(HandleGammaPoint)then
               do ismall=1,12
                  if (Lttc%small_ik(ismall,1).eq.iq) den_smallk(:,:,iw,ismall) = invW
               enddo
            endif
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      ! Gamma point handling
      if(HandleGammaPoint)then
         !
         num_k=1
         do ismall=2,12
            if (Lttc%small_ik(ismall,2).eq.1) num_k=num_k+1
         enddo
         !
         !if num_k only is one include also next nearset points
         if(num_k.eq.1)then
            num_k=1
            do ismall=2,12
               if (Lttc%small_ik(ismall,2).le.2) num_k=num_k+1
            enddo
         endif
         !
         !calc average epsinv for small k
         do ismall=1,num_k
            den_smallk_avrg = den_smallk_avrg + den_smallk(:,:,:,ismall)/num_k
         enddo
         !
         !W at gamma point should be real
         den_smallk_avrg = real(den_smallk_avrg)
         !
         !Add the Gamma point value - element not summed in the iq loop
         do iw=1,Nmats
            Wmats%screened_local(:,:,iw) = Wmats%screened_local(:,:,iw) + matmul(den_smallk_avrg(:,:,iw),Umats%screened(:,:,iw,Umats%iq_gamma))/Nkpt
         enddo
         !
         deallocate(den_smallk,den_smallk_avrg)
         !
      endif
      deallocate(invW)
      !
   end subroutine calc_W_edmft


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes [ 1 - U*Pi ]^-1 * Pi - GW+EDMFT
   !---------------------------------------------------------------------------!
   subroutine calc_chi_full(Chi,Umats,Pmats,Lttc)
      !
      use parameters
      use utils_misc
      use utils_fields
      use linalg, only : zeye, inv
      implicit none
      !
      type(BosonicField),intent(inout)      :: Chi
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: invW(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      !
      !
      if(verbose)write(*,"(A)") "---- calc_chi_full"
      !
      !
      ! Check on the input Fields
      if(.not.Chi%status) stop "Chi not properly initialized."
      if(.not.Umats%status) stop "Umats not properly initialized."
      if(.not.Pmats%status) stop "Pmats not properly initialized."
      if(Chi%Nkpt.eq.0) stop "Chi k dependent attributes not properly initialized."
      if(Umats%Nkpt.eq.0) stop "Umats k dependent attributes not properly initialized."
      if(Pmats%Nkpt.eq.0) stop "Pmats k dependent attributes not properly initialized."
      if(Umats%iq_gamma.lt.0) stop "Umats iq_gamma not defined."
      !
      Nbp = Chi%Nbp
      Nkpt = Chi%Nkpt
      Beta = Chi%Beta
      Nmats = Chi%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "Either Umats and/or Pmats have different orbital dimension with respect to Chi."
      if(all([Umats%Nkpt-Nkpt,Pmats%Nkpt-Nkpt].ne.[0,0])) stop "Either Umats and/or Pmats have different number of k-points with respect to Chi."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "Either Umats and/or Pmats have different Beta with respect to Chi."
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))then
         Nmats = minval([Chi%Npoints,Umats%Npoints,Pmats%Npoints])
         write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller: "//str(Nmats)
      endif
      !
      allocate(invW(Nbp,Nbp));invW=czero
      call clear_attributes(Chi)
      Chi%bare = Umats%bare
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Chi,Lttc),&
      !$OMP PRIVATE(iq,iw,invW)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if(iq.eq.Umats%iq_gamma)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            invW = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened(:,:,iw,iq))
            !
            ! [ 1 - U*Pi ]^-1
            call inv(invW)
            !
            ! [ 1 - U*Pi ]^-1 * Pi
            Chi%screened(:,:,iw,iq) = matmul(invW,Pmats%screened(:,:,iw,iq))
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(invW)
      !
      call BosonicKsum(Chi)
      !
   end subroutine calc_chi_full


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes [ 1 - U*Pi ]^-1 * Pi - EDMFT
   !---------------------------------------------------------------------------!
   subroutine calc_chi_edmft(Chi,Umats,Pmats,Lttc)
      !
      use parameters
      use utils_misc
      use utils_fields
      use linalg, only : zeye, inv
      implicit none
      !
      type(BosonicField),intent(inout)      :: Chi
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: invW(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      !
      !
      if(verbose)write(*,"(A)") "---- calc_chi_edmft"
      !
      !
      ! Check on the input Fields
      if(.not.Chi%status) stop "Chi not properly initialized."
      if(.not.Umats%status) stop "Umats not properly initialized."
      if(.not.Pmats%status) stop "Pmats not properly initialized."
      if(Chi%Nkpt.ne.0) stop "Chi k dependent attributes are supposed to be unallocated."
      if(Umats%Nkpt.eq.0) stop "Umats k dependent attributes not properly initialized."
      if(Pmats%Nkpt.ne.0) stop "Pmats k dependent attributes are supposed to be unallocated."
      if(Umats%iq_gamma.lt.0) stop "Umats iq_gamma not defined."
      !
      Nbp = Chi%Nbp
      Nkpt = Umats%Nkpt
      Beta = Chi%Beta
      Nmats = Chi%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "Either Umats and/or Pmats have different orbital dimension with respect to Chi."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "Either Umats and/or Pmats have different Beta with respect to Chi."
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))then
         Nmats = minval([Chi%Npoints,Umats%Npoints,Pmats%Npoints])
         write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller: "//str(Nmats)
      endif
      !
      allocate(invW(Nbp,Nbp));invW=czero
      call clear_attributes(Chi)
      Chi%bare_local = Umats%bare_local
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Chi,Lttc),&
      !$OMP PRIVATE(iq,iw,invW)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if(iq.eq.Umats%iq_gamma)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            invW = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened_local(:,:,iw))
            !
            ! [ 1 - U*Pi ]^-1
            call inv(invW)
            !
            ! [ 1 - U*Pi ]^-1 * Pi
            Chi%screened_local(:,:,iw) = Chi%screened_local(:,:,iw) + matmul(invW,Pmats%screened_local(:,:,iw))/Nkpt
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(invW)
      !
   end subroutine calc_chi_edmft


   !---------------------------------------------------------------------------!
   !PURPOSE: Read frequancy dependent interactions from SPEX files.
   !TEST ON: 20-10-2020 (both doAC)
   !---------------------------------------------------------------------------!
   subroutine read_U_spex_full(Umats,save2readable,LocalOnly,pathOUTPUT,doAC)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only :  pathINPUT, UfullStructure, Uthresh
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      logical,intent(in)                    :: save2readable
      logical,intent(in)                    :: LocalOnly
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: doAC
      !
      logical                               :: filexists,ACdone,doAC_
      character(len=256)                    :: file_spex,path,pathOUTPUT_
      integer                               :: unit,Nkpt
      integer                               :: iq,iw,iqread,Nbp_spex
      integer                               :: idum,Nspin_spex,Norb_spex,Nfreq
      integer                               :: ib1,ib2,iw1,iw2
      real(8),allocatable                   :: wread(:),wmats(:)
      complex(8),allocatable                :: D1(:,:),D2(:,:),D3(:,:),imgFact(:,:,:)
      complex(8),allocatable                :: Utmp(:,:)
      logical,allocatable                   :: RevSym(:,:,:)
      type(BosonicField)                    :: Ureal
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_U_spex_full"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Check on the input Boson
      if(.not.Umats%status) stop "BosonicField not properly initialized."
      if((.not.LocalOnly).and.(.not.allocated(Umats%bare))) stop "Requested k-dependence but bare attribute not allocated."
      if((.not.LocalOnly).and.(.not.allocated(Umats%screened))) stop "Requested k-dependence but screened attribute not allocated."
      if(LocalOnly.and.allocated(Umats%bare)) stop "Bare K-dependent attributes is present but not used."
      if(LocalOnly.and.allocated(Umats%screened)) stop "Screened K-dependent attributes is present but not used."
      allocate(wmats(Umats%Npoints));wmats=0d0
      wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      !
      !
      ! for memory demanding calculations the execute_command_line does not work
      ! so I have to create the VW_imag directory before (but I'm keeping it here).
      ! In the end instead of checking for the existence of the VW_imag folder
      ! I will check for the first Q point file within the folder.
      !
      !
      ! Check if the data on the Matsubara axis are present
      path = reg(pathOUTPUT_)//"VW_imag"
      !call inquireDir(reg(path),ACdone,hardstop=.false.,verb=verbose)
      call inquireFile(reg(path)//"/VW.Q0001.DAT",ACdone,hardstop=.false.,verb=verbose)
      doAC_ = .not.ACdone
      if(present(doAC)) doAC_ = doAC .or. doAC_
      !
      !
      ! This has to be done before allocation of large files otherwise the "system" command will not work
      if(doAC_.and.(.not.LocalOnly))then
         call createDir(reg(trim(pathOUTPUT_)//"VW_imag"),verb=verbose)
         if(save2readable)then
            call createDir(reg(trim(pathOUTPUT_)//"VW_imag_readable"),verb=verbose)
            call createDir(reg(trim(pathOUTPUT_)//"VW_real_readable"),verb=verbose)
         endif
      endif
      !
      call init_Uelements(int(sqrt(dble(Umats%Nbp))),PhysicalUelements)
      !
      !
      ! Perform cnalytical continuation on real interaction or load existing files
      if(doAC_) then
         !
         !---------------------------------------------------------------------!
         !
         write(*,"(A)")"     Performing Analytic continuation to get UcRPA(iq,iw)."
         !
         ! Allocations from dimensions written in VW.Q0001.DAT file
         path = reg(pathINPUT)//"VW_real/VW.Q0001.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read")
         read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
         close(unit)
         !
         Nbp_spex = Norb_spex**2
         allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         write(*,"(A,I)")"     Real frequencies: ",Nfreq
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Nbp_spex) stop "Size of given BosonicField and VW_real orbital space do not coincide."
         !
         ! Look for the Number of SPEX files. Which are supposed to be ordered.
         Nkpt = 0
         do iq=1,2000
            file_spex = reg(pathINPUT)//"VW_real/VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) exit
            Nkpt = Nkpt + 1
         enddo
         write(*,"(A,I)") "     The number of SPEX files (Nkpt) in VW_real is: ",Nkpt
         if((.not.LocalOnly).and.(Umats%Nkpt.ne.Nkpt)) stop "Number of k-points of given BosonicField and number of VW_real k-points do not coincide."
         !
         ! Allocate the Bosonic field on the real axis
         if(LocalOnly)then
            call AllocateBosonicField(Ureal,Norb_spex,Nfreq,0,Nkpt=0,name="Uspex(w)",Beta=Umats%Beta,no_bare=.true.)
         else
            call AllocateBosonicField(Ureal,Norb_spex,Nfreq,0,Nkpt=Nkpt,name="Uspex(k,w)",Beta=Umats%Beta,no_bare=.true.)
         endif
         !
         ! Read VW_real accumulating local attribute and optionally storing the k-dependent part
         path = reg(pathINPUT)//"VW_real/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,verb=verbose) !redundant control
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
            if (iq.ne.iqread) stop "iqread.ne.iq"
            !
            read(unit) wread
            wread = H2eV*wread
            if (dabs(wread(1)).gt.eps) stop "wread(1) not zero"
            !
            do iw=0,Nfreq
               read(unit) Utmp
               if(iw.eq.0) then
                  !V(:,:,iq)=vwtmp(:,:)/Nkpt/Nkpt
                  !bare values on Matsubara are the same so no need to use bare attributes of Ureal
                  Umats%bare_local = Umats%bare_local + H2eV*Utmp/(Nkpt**3)
                  if(.not.LocalOnly) Umats%bare(:,:,iq) = H2eV*Utmp/(Nkpt**2)
               else
                  !Ur(:,:,iw,iq)=vwtmp(:,:)/Nkpt/Nkpt
                  Ureal%screened_local(:,:,iw) = Ureal%screened_local(:,:,iw) + H2eV*Utmp/(Nkpt**3)
                  if(.not.LocalOnly) Ureal%screened(:,:,iw,iq) = H2eV*Utmp/(Nkpt**2)
               endif
            enddo
            !
            close(unit)
            !
         enddo !iq
         !
         ! Allocate the temporary quantities needed by the Analytical continuation
         allocate(D1(Nbp_spex,Nbp_spex));D1=czero
         allocate(D2(Nbp_spex,Nbp_spex));D2=czero
         allocate(D3(Nbp_spex,Nbp_spex));D3=czero
         !
         ! Check if any local Urpa component has inverted Im/Re symmetry
         do ib1=1,Nbp_spex
            do ib2=1,Nbp_spex
               !
               if(abs(Umats%bare_local(ib1,ib2)).lt.Uthresh)then
                  Umats%bare_local(ib1,ib2)=czero ; Ureal%screened_local(ib1,ib2,:)=czero
                  Umats%bare_local(ib2,ib1)=czero ; Ureal%screened_local(ib2,ib1,:)=czero
               endif
               !
               if(PhysicalUelements%Full_All(ib1,ib2).and.abs(real(Umats%bare_local(ib1,ib2))).lt.abs(aimag(Umats%bare_local(ib1,ib2))))then
                  write(*,"(A,2I)")"Element: ",ib1,ib2
                  write(*,"(A,F)")"Re[Ubare(w=inf)]: ",real(Umats%bare_local(ib1,ib2))
                  write(*,"(A,F)")"Im[Ubare(w=inf)]: ",aimag(Umats%bare_local(ib1,ib2))
                  stop "Something wrong: Uloc cannot have inverted Re/Im parity."
               endif
               !
            enddo
         enddo
         !
         ! Analytical continuation of the local component to imag axis using spectral rep
         call cpu_time(start)
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
         !$OMP END DO
         !$OMP END PARALLEL
         call cpu_time(finish)
         deallocate(D1,D2,D3)
         write(*,"(A,F)") "     UcRPA(w) --> UcRPA(iw) cpu timing: ", finish-start
         !
         !
         if(.not.LocalOnly)then
            !
            ! Allocate the temporary quantities needed by the Analytical continuation
            allocate(D1(Nbp_spex,Nbp_spex));D1=czero
            allocate(D2(Nbp_spex,Nbp_spex));D2=czero
            allocate(D3(Nbp_spex,Nbp_spex));D3=czero
            !
            ! Analytical continuation of all the K-points to imag axis using spectral rep
            allocate(imgFact(Nbp_spex,Nbp_spex,2));imgFact=cone
            allocate(RevSym(Nbp_spex,Nbp_spex,Umats%Nkpt));RevSym=.false.
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(Nbp_spex,wmats,wread,Nfreq,Ureal,Umats,UfullStructure,verbose,RevSym,Uthresh),&
            !$OMP PRIVATE(iq,ib1,ib2,iw1,iw2,D1,D2,D3,Utmp,imgFact)
            !$OMP DO
            do iq=1,Umats%Nkpt
               !
               !Some elements of U, usually the k-dependent ones, might have inverted Im/Re symmetry
               imgFact=cone
               do ib1=1,Nbp_spex
                  do ib2=1,Nbp_spex
                     !
                     if(abs(Umats%bare(ib1,ib2,iq)).lt.Uthresh)then
                        Umats%bare(ib1,ib2,iq)=czero ; Ureal%screened(ib1,ib2,:,iq)=czero
                        Umats%bare(ib2,ib1,iq)=czero ; Ureal%screened(ib2,ib1,:,iq)=czero
                     endif
                     !
                     if(UfullStructure.and.(abs(real(Umats%bare(ib1,ib2,iq))).lt.abs(aimag(Umats%bare(ib1,ib2,iq)))) )then
                        imgFact(ib1,ib2,1) = -img !this correspond to dividing by I
                        imgFact(ib1,ib2,2) = +img !this correspond to multiplying by I
                        !
                        RevSym(ib1,ib2,iq) = .true.
                        !
                     endif
                     !
                  enddo
               enddo
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
               if(verbose)print *, "     UcRPA(q,iw) - done iq: ",iq
            enddo !iq
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            !
            !Deal with elements with inverted Im/Re symmetry
            if(UfullStructure)then
               !
               unit = free_unit()
               path = reg(pathOUTPUT_)//"ReversedSymmetry.DAT"
               open(unit=unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               write(unit,"(3A5)")"iq","ib1","ib2"
               do iq=1,Umats%Nkpt
                  do ib1=1,Nbp_spex
                     do ib2=ib1,Nbp_spex
                        if(RevSym(ib1,ib2,iq).and.RevSym(ib2,ib1,iq))then
                           write(unit,"(5I5)")iq,PhysicalUelements%Full_Map(ib1,ib2,:)
                        elseif(RevSym(ib1,ib2,iq).and.(.not.RevSym(ib2,ib1,iq)))then
                           write(unit,"(5I5,A)")iq,PhysicalUelements%Full_Map(ib1,ib2,:),"   WARNING - non symmetrical with orbital index exchange."
                        endif
                     enddo
                  enddo
               enddo
               close(unit)
               !
            else
               !
               unit = free_unit()
               path = reg(pathOUTPUT_)//"ReversedSymmetry_Removed.DAT"
               open(unit=unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               write(unit,"(3A5)")"iq","ib1","ib2"
               do iq=1,Umats%Nkpt
                  do ib1=1,Nbp_spex
                     do ib2=1,Nbp_spex
                        if(abs(real(Umats%bare(ib1,ib2,iq))).lt.abs(aimag(Umats%bare(ib1,ib2,iq))))then
                           write(unit,"(5I5)")iq,PhysicalUelements%Full_Map(ib1,ib2,:)
                           Umats%bare(ib1,ib2,iq) = dcmplx(0d0,0d0)
                           Umats%screened(ib1,ib2,:,iq) = dcmplx(0d0,0d0)
                        endif
                     enddo
                  enddo
               enddo
               close(unit)
               !
            endif
            !
            deallocate(D1,D2,D3,imgFact,RevSym)
            write(*,"(A,F)") "     UcRPA(q,w) --> UcRPA(q,iw) cpu timing: ", finish-start
            !
            call BosonicKsum(Umats)
            !
         endif !LocalOnly
         deallocate(Utmp)
         !
         ! Print out the transformed stuff - local
         call dump_BosonicField(Umats,reg(pathOUTPUT_),"Uloc_mats.DAT")
         call dump_BosonicField(Ureal,reg(pathOUTPUT_),"Uloc_real.DAT",wread)
         !
         ! Print out the transformed stuff - Kdep
         if(.not.LocalOnly)then
            call dump_BosonicField(Umats,reg(pathOUTPUT_)//"VW_imag/",.true.)
            if(save2readable)then
               call dump_BosonicField(Umats,reg(pathOUTPUT_)//"VW_imag_readable/",.not.save2readable)
               call dump_BosonicField(Ureal,reg(pathOUTPUT_)//"VW_real_readable/",.not.save2readable,axis=wread)
            endif
         endif
         !
         deallocate(wread)
         !
         if(verbose)call checkAnalyticContinuation(Umats,Ureal)
         call DeallocateBosonicField(Ureal)
         !
         !---------------------------------------------------------------------!
         !
      else
         !
         !---------------------------------------------------------------------!
         !
         write(*,"(A)")"     Reading UcRPA(q,iw) from "//reg(pathOUTPUT_)//"VW_imag/"
         !
         ! Allocations from dimensions written in W.Q0001.DAT file
         path = reg(pathOUTPUT_)//"VW_imag/VW.Q0001.DAT"
         call inquireFile(reg(path),filexists,verb=verbose)
         !
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read")
         read(unit)idum,Nspin_spex,Norb_spex,Nfreq
         close(unit)
         !
         Nbp_spex = Norb_spex**2
         allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         write(*,"(A,I)")"     Matsubara frequencies: ",Nfreq
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Nbp_spex) stop "Size of given BosonicField and VW_imag orbital space do not coincide."
         if(Umats%Npoints.ne.Nfreq) stop "Number of VW_imag Matsubara points and bosonic field mesh does not coincide."
         !
         ! Look for the Number of SPEX files. Which are supposed to be ordered.
         Nkpt = 0
         do iq=1,2000
            file_spex = reg(pathOUTPUT_)//"VW_imag/VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) exit
            Nkpt = Nkpt + 1
         enddo
         write(*,"(A,I)") "     The number of SPEX files (Nkpt) in VW_imag is: ",Nkpt
         if((.not.LocalOnly).and.(Umats%Nkpt.ne.Nkpt)) stop "Number of k-points of given BosonicField and number of VW_imag k-points do not coincide."
         !
         ! Read VW_imag accumulating local attribute and optionally storing the k-dependent part
         path = reg(pathOUTPUT_)//"VW_imag/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,verb=verbose) !redundant control
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
            if (iq.ne.iqread) stop "iqread.ne.iq"
            !
            read(unit) wread
            !wread = H2eV*wread
            do iw=1,Nfreq
               if(dabs(wread(iw)-wmats(iw)).gt.eps) Then
                  write(*,"(F)")dabs(wread(iw)-wmats(iw)),iw,iq
                  stop "wread.ne.wmats"
               endif
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
   end subroutine read_U_spex_full


   !---------------------------------------------------------------------------!
   !PURPOSE: Read screened interaction tensor from Ucrpa(0)
   !---------------------------------------------------------------------------!
   subroutine read_U_spex_Uloc0(Umat,pathOUTPUT)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only :  pathINPUT
      implicit none
      !
      real(8),allocatable,intent(inout)     :: Umat(:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      logical                               :: Umatsxists,Urealxists,SPEXxists
      character(len=256)                    :: file_spex,path,pathOUTPUT_
      integer                               :: unit
      integer                               :: iq,iw,Nbp_spex,Nkpt
      integer                               :: iqread,Nspin_spex,Norb_spex,Nfreq
      complex(8),allocatable                :: Utmp(:,:)
      type(BosonicField)                    :: Uread
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_U_spex_Uloc0"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Look for Uloc_mats.DAT and Uloc_real.DAT
      path=reg(pathINPUT)//"Uloc_mats.DAT"
      call inquireFile(reg(path),Umatsxists,hardstop=.false.,verb=verbose)
      path=reg(pathINPUT)//"Uloc_real.DAT"
      call inquireFile(reg(path),Urealxists,hardstop=.false.,verb=verbose)
      path=reg(pathINPUT)//"VW_real" !/VW.Q0001.DAT"
      call inquireDir(reg(path),SPEXxists,hardstop=.false.,verb=verbose)
      !
      if(Umatsxists)then
         !
         call AllocateBosonicField(Uread,size(Umat,dim=1),1,0)
         !
         call read_BosonicField(Uread,reg(pathINPUT),"Uloc_mats.DAT")
         Umat = Uread%screened_local(:,:,1)
         call dump_matrix(Umat,reg(pathINPUT//"Umat.DAT"))
         call DeallocateBosonicField(Uread)
         return
         !
      else if(Urealxists)then
         !
         call AllocateBosonicField(Uread,size(Umat,dim=1),1,0)
         !
         call read_BosonicField(Uread,reg(pathINPUT),"Uloc_real.DAT")
         Umat = Uread%screened_local(:,:,1)
         call dump_matrix(Umat,reg(pathINPUT//"Umat.DAT"))
         call DeallocateBosonicField(Uread)
         return
         !
      else if(SPEXxists)then
         !
         Nkpt=0
         path = reg(pathINPUT)//"VW_real/"
         do iq=1,2000
            !
            file_spex = reg(path)//"VW_real/VW.Q"//str(iq,4)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"VW_real/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),SPEXxists,hardstop=.false.,verb=verbose)
            !
            ! This mathod looks uglier than the previous but here I don't
            ! have to store K-dep, just sum up what's present.
            if(.not.SPEXxists)cycle
            Nkpt=Nkpt+1
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit)iqread,Nspin_spex,Norb_spex,Nfreq
            !
            Nbp_spex = Norb_spex**2
            allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
            if(iq.eq.1)call assert_shape(Umat,[Nbp_spex,Nbp_spex],"read_spex_Uloc","Umat")
            !
            read(unit) !wread
            !
            do iw=0,1
               read(unit) Utmp
               if(iw.eq.1) then
                  Umat = Umat + H2eV*Utmp
               endif
            enddo
            !
            close(unit)
            deallocate(Utmp)
            !
         enddo !iq
         Umat = Umat/(Nkpt**3)
         call dump_matrix(Umat,reg(reg(pathOUTPUT_)//"Umat.DAT"))
         return
         !
      else
         stop "No useful interaction file found."
      endif
      !
   end subroutine read_U_spex_Uloc0


   !---------------------------------------------------------------------------!
   !PURPOSE: Check if the AC alters the bare and screened values
   !TEST ON: 21-10-2020
   !---------------------------------------------------------------------------!
   subroutine checkAnalyticContinuation(Umats,Ureal)
      !
      use parameters
      use utils_misc
      use input_vars, only :  pathINPUTtr
      implicit none
      !
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Ureal
      character(len=256)                    :: path
      integer                               :: iq,ib1,ib2
      integer                               :: unit,Nbp,Nmats,Nreal,Nkpt
      real(8)                               :: ReErr,ImErr,thresh=1e-4
      real(8),allocatable                   :: ReErrMat(:,:),ImErrMat(:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- checkAnalyticContinuation"
      if(Umats%Nbp.ne.Ureal%Nbp) stop "Umats%Nbp.ne.Ureal%Nbp"
      Nbp = Umats%Nbp
      Nmats = Umats%Npoints
      Nreal = Ureal%Npoints
      allocate(ReErrMat(Nbp,Nbp));ReErrMat=0d0
      allocate(ImErrMat(Nbp,Nbp));ImErrMat=0d0
      !
      !
      ! Check the difference betqween bare values induced by thecutoff in the matsubara frequency
      unit = free_unit()
      path = reg(pathINPUTtr)//"ACcutoffError.DAT"
      open(unit=unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(A,1I5,A,1E14.7)")"Difference between Umats_bare value and last screened frequency: ",Nmats," thresh:",thresh
      !
      ReErrMat=0d0;ImErrMat=0d0
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            ReErr = abs(real(Umats%bare_local(ib1,ib2)) - real(Umats%screened_local(ib1,ib2,Nmats)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(aimag(Umats%bare_local(ib1,ib2)) - aimag(Umats%screened_local(ib1,ib2,Nmats)))
            if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
            !
         enddo
      enddo
      write(unit,"(A)")"Real part - local projection"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      !
      if((Umats%Nkpt.eq.Ureal%Nkpt).and.(Ureal%Nkpt.gt.0))then
         Nkpt = Umats%Nkpt
         do iq=1,Nkpt
            !
            ReErrMat=0d0;ImErrMat=0d0
            do ib1=1,Nbp
               do ib2=1,Nbp
                  !
                  ReErr = abs(real(Umats%bare(ib1,ib2,iq)) - real(Umats%screened(ib1,ib2,Nmats,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(aimag(Umats%bare(ib1,ib2,iq)) - aimag(Umats%screened(ib1,ib2,Nmats,iq)))
                  if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
                  !
               enddo
            enddo
            !
            write(unit,"(A,1I5)")"Real part - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            !
         enddo !iq
      endif
      close(unit)
      !
      !
      ! Check that the screened and bare values of Umats and Ureal are close enough
      unit = free_unit()
      path = reg(pathINPUTtr)//"ACcheck.DAT"
      open(unit=unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(A,1E14.7)")"Difference between asymptotic behaviour of Ureal and Umats. Thresh:",thresh
      !
      ReErrMat=0d0;ImErrMat=0d0
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            ReErr = abs(real(Umats%screened_local(ib1,ib2,1)) - real(Ureal%screened_local(ib1,ib2,1)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(aimag(Umats%screened_local(ib1,ib2,1)) - aimag(Ureal%screened_local(ib1,ib2,1)))
            if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
            !
         enddo
      enddo
      write(unit,"(A)")"Real part - local projection - screened limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection - screened limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      !
      ReErrMat=0d0;ImErrMat=0d0
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            ReErr = abs(real(Umats%screened_local(ib1,ib2,Nmats)) - real(Ureal%screened_local(ib1,ib2,Nreal)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(aimag(Umats%screened_local(ib1,ib2,Nmats)) - aimag(Ureal%screened_local(ib1,ib2,Nreal)))
            if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
            !
         enddo
      enddo
      write(unit,"(A)")"Real part - local projection - bare limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection - bare limit"
      do ib1=1,Nbp
         write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      !
      if((Umats%Nkpt.eq.Ureal%Nkpt).and.(Ureal%Nkpt.gt.0))then
         Nkpt = Umats%Nkpt
         !
         do iq=1,Nkpt
            !
            ReErrMat=0d0;ImErrMat=0d0
            do ib1=1,Nbp
               do ib2=1,Nbp
                  !
                  ReErr = abs(real(Umats%screened(ib1,ib2,1,iq)) - real(Ureal%screened(ib1,ib2,1,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(aimag(Umats%screened(ib1,ib2,1,iq)) - aimag(Ureal%screened(ib1,ib2,1,iq)))
                  if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
                  !
               enddo
            enddo
            !
            write(unit,"(A,1I5)")"Real part - screened limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - screened limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            !
            ReErrMat=0d0;ImErrMat=0d0
            do ib1=1,Nbp
               do ib2=1,Nbp
                  !
                  ReErr = abs(real(Umats%screened(ib1,ib2,Nmats,iq)) - real(Ureal%screened(ib1,ib2,Nreal,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(aimag(Umats%screened(ib1,ib2,Nmats,iq)) - aimag(Ureal%screened(ib1,ib2,Nreal,iq)))
                  if(ImErr.gt.thresh) ImErrMat(ib1,ib2) = ImErr
                  !
               enddo
            enddo
            !
            write(unit,"(A,1I5)")"Real part - bare limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - bare limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E14.7)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            !
         enddo !iq
      endif
      close(unit)
      !
   end subroutine checkAnalyticContinuation


   !---------------------------------------------------------------------------!
   !PURPOSE: Create the static interaction tensor from user-given parameters
   !---------------------------------------------------------------------------!
   subroutine build_Umat_singlParam(Umat,Uaa,Uab,J)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      implicit none
      !
      real(8),allocatable,intent(inout)     :: Umat(:,:)
      real(8),intent(in)                    :: Uaa,Uab,J
      !
      integer                               :: Nbp,Norb,Nflavor
      integer                               :: ib1,ib2
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- build_Umat_singlParam"
      !
      !
      ! Check on the input matrices
      Nbp = size(Umat,dim=1)
      Norb = Nbp/Nspin
      Nflavor = Norb*Nspin
      if((Nspin.eq.2).and.(mod(Nbp,2).ne.0.0)) stop "Wrong matrix dimension."
      call init_Uelements(Norb,PhysicalUelements)
      call assert_shape(Umat,[Nbp,Nbp],"build_Umat_singlParam","Umat")
      !
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            if(PhysicalUelements%Flav_Uloc(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uaa,0d0)
            if(PhysicalUelements%Flav_U1st(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab,0d0)
            if(PhysicalUelements%Flav_U2nd(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab-J,0d0)
            !
         enddo
      enddo
      !
   end subroutine build_Umat_singlParam
   !
   !
   subroutine build_Umat_multiParam(Umat,Uaa,Uab,J)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      implicit none
      !
      real(8),allocatable,intent(inout)     :: Umat(:,:)
      real(8),allocatable,intent(in)        :: Uaa(:),Uab(:,:),J(:,:)
      !
      integer                               :: Nbp,Norb,Nflavor
      integer                               :: ib1,ib2,iorb,jorb
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- build_Umat_multiParam"
      !
      !
      ! Check on the input matrices
      Nbp = size(Umat,dim=1)
      Norb = Nbp/Nspin
      Nflavor = Norb*Nspin
      if((Nspin.eq.2).and.(mod(Nbp,2).ne.0.0)) stop "Wrong matrix dimension."
      call init_Uelements(Norb,PhysicalUelements)
      call assert_shape(Umat,[Nbp,Nbp],"build_Umat_multiParam","Umat")
      call assert_shape(Uaa,[Norb],"build_Umat_multiParam","Uaa")
      call assert_shape(Uab,[Norb,Norb],"build_Umat_multiParam","Uab")
      call assert_shape(J,[Norb,Norb],"build_Umat_multiParam","J")
      !
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            iorb = PhysicalUelements%Flav_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Flav_Map(ib1,ib2,2)
            !
            if(PhysicalUelements%Flav_Uloc(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uaa(iorb),0d0)
            if(PhysicalUelements%Flav_U1st(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab(iorb,jorb),0d0)
            if(PhysicalUelements%Flav_U2nd(ib1,ib2)) Umat(ib1,ib2) = dcmplx(Uab(iorb,jorb)-J(iorb,jorb),0d0)
            !
         enddo
      enddo
      !
   end subroutine build_Umat_multiParam


   !---------------------------------------------------------------------------!
   !PURPOSE: Create the freq. dependent interaction tensor from user-given parameters
   !---------------------------------------------------------------------------!
   subroutine build_Uretloc_singlParam(Umats,Uaa,Uab,J,g_eph,wo_eph)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : Nreal, wrealMax
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      real(8),intent(in)                    :: Uaa,Uab,J
      real(8),allocatable,intent(in)        :: g_eph(:),wo_eph(:)
      !
      integer                               :: Nbp,Norb,Nph
      integer                               :: ib1,ib2
      integer                               :: iw,iw1,iw2,iph,iwp
      real(8)                               :: RealU,ImagU
      real(8),allocatable                   :: wreal(:),wmats(:)
      complex(8),allocatable                :: D1(:,:),D2(:,:),D3(:,:)
      complex(8),allocatable                :: Utmp(:,:)
      type(BosonicField)                    :: Ureal
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- build_Uretloc_singlParam"
      !
      !
      ! Check on the input field
      if(.not.Umats%status) stop "BosonicField not properly initialized."
      !
      Nbp = Umats%Nbp
      Norb = int(sqrt(dble(Nbp)))
      Nph = size(g_eph)
      !
      call init_Uelements(Norb,PhysicalUelements)
      if(size(g_eph).ne.size(wo_eph)) stop "Phonon sizes does not match."
      !
      allocate(wmats(Umats%Npoints));wmats=0d0
      wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      allocate(wreal(Nreal));wreal=0d0
      wreal = linspace(0d0,+wrealMax,Nreal)
      !
      call AllocateBosonicField(Ureal,Nbp,Nreal,0)
      !
      !setting the bare values
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(Uaa,0d0)
            if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(Uab,0d0)
            if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(J,0d0)
            if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(J,0d0)
            !
         enddo
      enddo
      !
      !setting the phonons
      do ib1=1,Nbp
         do ib2=1,Nbp
            do iph=1,Nph
               iwp=minloc(wreal-wo_eph(iph),dim=1)
               do iw=1,Nreal
                  !
                  RealU = 2*(g_eph(iph)**2)*wo_eph(iph) / ( (wreal(iw)**2) - (wo_eph(iph)**2) )
                  ImagU=0d0
                  if(iw.eq.iwp) ImagU = -pi*(g_eph(iph)**2)
                  !
                  Ureal%screened_local(ib1,ib2,iw) = Umats%bare_local(ib1,ib2) + dcmplx(RealU,ImagU)
                  !
               enddo
            enddo
         enddo
      enddo
      !
      ! Allocate the temporary quantities needed by the Analytical continuation
      allocate(Utmp(Nbp,Nbp));Utmp=czero
      allocate(D1(Nbp,Nbp));D1=czero
      allocate(D2(Nbp,Nbp));D2=czero
      allocate(D3(Nbp,Nbp));D3=czero
      !
      ! Analytical continuation of the local component to imag axis using spectral rep
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,wmats,wreal,Nreal,Ureal,Umats),&
      !$OMP PRIVATE(ib1,ib2,iw1,iw2,D1,D2,D3,Utmp)
      !$OMP DO
      do iw1=1,Umats%Npoints
         Utmp=czero
         do iw2=1,Nreal-2,2
            !
            do ib1=1,Nbp
               do ib2=1,Nbp
                  D1(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2)   )/pi
                  D2(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+1) )/pi
                  D3(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+2) )/pi
               enddo
            enddo
            !
            !D(-w)=-D(w), integrate using Simpson method
            if(wreal(iw2).gt.0.d0) then
               Utmp(:,:) = Utmp(:,:) + ( D1(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2)  ) - D1(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2)  ) ) *(wreal(iw2+1)-wreal(iw2))/3.d0
               Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+1)) ) *(wreal(iw2+1)-wreal(iw2))*4.d0/3.d0
               Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+2)) ) *(wreal(iw2+1)-wreal(iw2))/3.d0
            elseif(dabs(wreal(iw2)).lt.1.d-12) then
               Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+1)) ) *(wreal(iw2+1)-wreal(iw2))*4.d0/3.d0
               Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+2)) ) *(wreal(iw2+1)-wreal(iw2))/3.d0
            endif
         enddo
         !
         do ib1=1,Nbp
            do ib2=1,Nbp
               Umats%screened_local(ib1,ib2,iw1) = Utmp(ib1,ib2) + Umats%bare_local(ib1,ib2)
            enddo
         enddo
         !
      enddo !iw1
      !
      !$OMP END DO
      !$OMP END PARALLEL
      call cpu_time(finish)
      deallocate(D1,D2,D3,Utmp,wmats,wreal)
      call DeallocateBosonicField(Ureal)
      write(*,"(A,F)") "Ue-ph(w) --> Ue-ph(iw) cpu timing:", finish-start
      !
   end subroutine build_Uretloc_singlParam
   !
   subroutine build_Uretloc_multiParam(Umats,Uaa,Uab,J,g_eph,wo_eph)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : Nreal, wrealMax
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      real(8),allocatable,intent(in)        :: Uaa(:),Uab(:,:),J(:,:)
      real(8),allocatable,intent(in)        :: g_eph(:),wo_eph(:)
      !
      integer                               :: Nbp,Norb,Nph
      integer                               :: iorb,jorb,ib1,ib2
      integer                               :: iw,iw1,iw2,iph,iwp
      real(8)                               :: RealU,ImagU
      real(8),allocatable                   :: wreal(:),wmats(:)
      complex(8),allocatable                :: D1(:,:),D2(:,:),D3(:,:)
      complex(8),allocatable                :: Utmp(:,:)
      type(BosonicField)                    :: Ureal
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      !
      !
      if(verbose)write(*,"(A)") "---- build_Uretloc_multiParam"
      !
      !
      ! Check on the input field
      if(.not.Umats%status) stop "BosonicField not properly initialized."
      !
      Nbp = Umats%Nbp
      Norb = int(sqrt(dble(Nbp)))
      Nph = size(g_eph)
      !
      call init_Uelements(Norb,PhysicalUelements)
      if(size(g_eph).ne.size(wo_eph)) stop "Phonon sizes does not match."
      !
      call assert_shape(Uaa,[Norb],"build_Uretloc_multiParam","Uaa")
      call assert_shape(Uab,[Norb,Norb],"build_Uretloc_multiParam","Uab")
      call assert_shape(J,[Norb,Norb],"build_Uretloc_multiParam","J")
      !
      allocate(wmats(Umats%Npoints));wmats=0d0
      wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      allocate(wreal(Nreal));wreal=0d0
      wreal = linspace(0d0,+wrealMax,Nreal)
      !
      call AllocateBosonicField(Ureal,Nbp,Nreal,0)
      !
      !setting the bare values
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            iorb = PhysicalUelements%Full_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Full_Map(ib1,ib2,2)
            !
            if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(Uaa(iorb),0d0)
            if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(Uab(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats%bare_local(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
            !
         enddo
      enddo
      !
      !setting the phonons
      do ib1=1,Nbp
         do ib2=1,Nbp
            do iph=1,Nph
               iwp=minloc(wreal-wo_eph(iph),dim=1)
               do iw=1,Nreal
                  !
                  RealU = 2*(g_eph(iph)**2)*wo_eph(iph) / ( (wreal(iw)**2) - (wo_eph(iph)**2) )
                  ImagU=0d0
                  if(iw.eq.iwp) ImagU = -pi*(g_eph(iph)**2)
                  !
                  Ureal%screened_local(ib1,ib2,iw) = Umats%bare_local(ib1,ib2) + dcmplx(RealU,ImagU)
                  !
               enddo
            enddo
         enddo
      enddo
      !
      ! Allocate the temporary quantities needed by the Analytical continuation
      allocate(Utmp(Nbp,Nbp));Utmp=czero
      allocate(D1(Nbp,Nbp));D1=czero
      allocate(D2(Nbp,Nbp));D2=czero
      allocate(D3(Nbp,Nbp));D3=czero
      !
      ! Analytical continuation of the local component to imag axis using spectral rep
      call cpu_time(start)
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,wmats,wreal,Nreal,Ureal,Umats),&
      !$OMP PRIVATE(ib1,ib2,iw1,iw2,D1,D2,D3,Utmp)
      !$OMP DO
      do iw1=1,Umats%Npoints
         Utmp=czero
         do iw2=1,Nreal-2,2
            !
            do ib1=1,Nbp
               do ib2=1,Nbp
                  D1(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2)   )/pi
                  D2(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+1) )/pi
                  D3(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+2) )/pi
               enddo
            enddo
            !
            !D(-w)=-D(w), integrate using Simpson method
            if(wreal(iw2).gt.0.d0) then
               Utmp(:,:) = Utmp(:,:) + ( D1(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2)  ) - D1(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2)  ) ) *(wreal(iw2+1)-wreal(iw2))/3.d0
               Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+1)) ) *(wreal(iw2+1)-wreal(iw2))*4.d0/3.d0
               Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+2)) ) *(wreal(iw2+1)-wreal(iw2))/3.d0
            elseif(dabs(wreal(iw2)).lt.1.d-12) then
               Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+1)) ) *(wreal(iw2+1)-wreal(iw2))*4.d0/3.d0
               Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wmats(iw1))-wreal(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wmats(iw1))+wreal(iw2+2)) ) *(wreal(iw2+1)-wreal(iw2))/3.d0
            endif
         enddo
         !
         do ib1=1,Nbp
            do ib2=1,Nbp
               Umats%screened_local(ib1,ib2,iw1) = Utmp(ib1,ib2) + Umats%bare_local(ib1,ib2)
            enddo
         enddo
         !
      enddo !iw1
      !
      !$OMP END DO
      !$OMP END PARALLEL
      call cpu_time(finish)
      deallocate(D1,D2,D3,Utmp,wmats,wreal)
      call DeallocateBosonicField(Ureal)
      write(*,"(A,F)") "Ue-ph(w) --> Ue-ph(iw) cpu timing:", finish-start
      !
   end subroutine build_Uretloc_multiParam


   !---------------------------------------------------------------------------!
   !PURPOSE: Given the Bosonic Field it extracts the screened interaction and
   ! retardation function.
   !---------------------------------------------------------------------------!
   subroutine calc_QMCinteractions(Umats,Uinst,Kfunct,sym)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only : Solver
      implicit none
      !
      type(BosonicField),intent(in)         :: Umats
      real(8),intent(inout)                 :: Uinst(:,:)
      real(8),intent(inout),optional        :: Kfunct(:,:,:)
      logical,intent(in),optional           :: sym
      !
      integer                               :: Nbp,Norb,Nflavor
      integer                               :: ib1,ib2,iorb,jorb
      integer                               :: iu1,iu2,ix1,ix2,ip1,ip2
      integer                               :: iw,itau
      real(8),allocatable                   :: wmats(:),tau(:)
      complex(8),allocatable                :: Kaux(:,:,:)
      logical                               :: Uloc,U1st,U2nd,retarded
      type(physicalU)                       :: PhysicalUelements
      logical                               :: sym_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_QMCinteractions"
      !
      !
      ! Check on the input field
      if(.not.Umats%status) stop "Umats not properly initialized."
      retarded=.false.
      if(present(Kfunct))retarded=.true.
      !
      sym_=.true.
      if(present(sym))sym_=sym
      !
      Nbp = Umats%Nbp
      Norb = int(sqrt(dble(Nbp)))
      Nflavor = Norb*Nspin
      !
      call init_Uelements(Norb,PhysicalUelements)
      call assert_shape(Uinst,[Nflavor,Nflavor],"calc_QMCinteractions","Uinst")
      Uinst=0d0
      if(retarded)then
         call assert_shape(Kfunct,[Nflavor,Nflavor,Solver%NtauB],"calc_QMCinteractions","Kfunct")
         allocate(Kaux(Nflavor,Nflavor,Umats%Npoints));Kaux=czero
         allocate(tau(Solver%NtauB));tau=0d0
         tau = linspace(0d0,Umats%Beta,Solver%NtauB)
         allocate(wmats(Umats%Npoints));wmats=0d0
         wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      endif
      !
      !setting the istantaneous values
      do ib1=1,Nflavor
         do ib2=1,Nflavor
            !
            !This is just for a more compact code
            Uloc = PhysicalUelements%Flav_Uloc(ib1,ib2)
            U1st = PhysicalUelements%Flav_U1st(ib1,ib2)
            U2nd = PhysicalUelements%Flav_U2nd(ib1,ib2)
            !
            !Orbital indexes
            iorb = PhysicalUelements%Flav_Map(ib1,ib2,1)
            jorb = PhysicalUelements%Flav_Map(ib1,ib2,2)
            !
            !The maps inside PhysicalUelements contain separately the orbital
            !indexes specifially for that representation. The matching between
            !the two is not done, so I have to do it here.
            !
            ! (aa)(bb) indexes in the Norb^2 representaion
            iu1 = iorb + Norb*(iorb-1)
            iu2 = jorb + Norb*(jorb-1)
            !
            ! (ab)(ba) indexes
            ix1 = iorb + Norb*(jorb-1)
            ix2 = jorb + Norb*(iorb-1)
            !
            ! (ab)(ab) indexes
            ip1 = iorb + Norb*(jorb-1)
            ip2 = iorb + Norb*(jorb-1)
            !
            if(Uloc) Uinst(ib1,ib2) = Umats%screened_local(iu1,iu2,1)
            if(U1st) Uinst(ib1,ib2) = Umats%screened_local(iu1,iu2,1)
            if(U2nd) Uinst(ib1,ib2) = Umats%screened_local(iu1,iu2,1) - (Umats%screened_local(ix1,ix2,1)+Umats%screened_local(ip1,ip2,1))/2d0
            !
            if(retarded)then
               !
               if(Uloc) Kaux(ib1,ib2,:) =  Umats%screened_local(iu1,iu2,:) - Umats%screened_local(iu1,iu2,1)
               if(U1st) Kaux(ib1,ib2,:) =  Umats%screened_local(iu1,iu2,:) - Umats%screened_local(iu1,iu2,1)
               if(U2nd) Kaux(ib1,ib2,:) =  Umats%screened_local(iu1,iu2,:) - (Umats%screened_local(ix1,ix2,:)+Umats%screened_local(ip1,ip2,:))/2d0 - &
                                          (Umats%screened_local(iu1,iu2,1) - (Umats%screened_local(ix1,ix2,1)+Umats%screened_local(ip1,ip2,1))/2d0)
               !same orbital - same spin screening
               if(Uloc.and.(ib2.gt.ib1)) then
                  Kaux(ib1,ib1,:) = Kaux(ib1,ib2,:)
                  Kaux(ib2,ib2,:) = Kaux(ib1,ib2,:)
               endif
               !
            endif
            !
         enddo
      enddo
      if(sym_)call check_Symmetry(Uinst,eps,enforce=.true.,hardstop=.false.,name="Uinst")
      !
      !setting the reterdation function
      if(retarded)then
         Kfunct=0d0
         do itau=2,Solver%NtauB-1
            !
            do iw=2,Umats%Npoints
               Kfunct(:,:,itau) = Kfunct(:,:,itau) - 2d0*Kaux(:,:,iw) * ( cos(wmats(iw)*tau(itau)) - 1d0 ) / ( Umats%Beta*wmats(iw)**2 )
            enddo
            !
            if(sym_)call check_Symmetry(Kfunct(:,:,itau),eps,enforce=.true.,hardstop=.false.,name="Kfunct_t"//str(itau))
            !
         enddo
         deallocate(Kaux,tau,wmats)
      endif
      !
   end subroutine calc_QMCinteractions


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the local effective interaction
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine calc_curlyU(curlyU,Wimp,Pimp,sym)
      !
      use parameters
      use utils_fields
      use utils_misc
      use linalg, only : zeye, inv, inv_sym
      implicit none
      !
      type(BosonicField),intent(inout)      :: curlyU
      type(BosonicField),intent(in)         :: Wimp
      type(BosonicField),intent(in)         :: Pimp
      logical,intent(in),optional           :: sym
      !
      complex(8),allocatable                :: invW(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nmats
      integer                               :: iw
      logical                               :: sym_
      !
      !
      if(verbose)write(*,"(A)") "---- calc_curlyU"
      !
      !
      ! Check on the input Fields
      if(.not.curlyU%status) stop "curlyU not properly initialized."
      if(.not.Wimp%status) stop "Wimp not properly initialized."
      if(.not.Pimp%status) stop "Pimp not properly initialized."
      if(curlyU%Nkpt.ne.0) stop "curlyU k dependent attributes are supposed to be unallocated."
      if(Wimp%Nkpt.ne.0) stop "Wimp k dependent attributes are supposed to be unallocated."
      if(Pimp%Nkpt.ne.0) stop "Pimp k dependent attributes are supposed to be unallocated."
      !
      sym_=.true.
      if(present(sym))sym_=sym
      !
      Nbp = curlyU%Nbp
      Beta = curlyU%Beta
      Nmats = curlyU%Npoints
      !
      if(all([Wimp%Nbp-Nbp,Pimp%Nbp-Nbp].ne.[0,0])) stop "Either Wimp and/or Pimp have different orbital dimension with respect to curlyU."
      if(all([Wimp%Beta-Beta,Pimp%Beta-Beta].ne.[0d0,0d0])) stop "Either Wimp and/or Pimp have different Beta with respect to curlyU."
      if(all([Wimp%Npoints-Nmats,Pimp%Npoints-Nmats].ne.[0,0]))   stop "Either Wimp and/or Pimp have different number of Matsubara points with respect to curlyU."
      !
      call clear_attributes(curlyU)
      !
      curlyU%bare_local = Wimp%bare_local
      !
      allocate(invW(Nbp,Nbp));invW=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nmats,Wimp,Pimp,curlyU),&
      !$OMP PRIVATE(iw,invW)
      !$OMP DO
      do iw=1,Nmats
         !
         invW = zeye(Nbp) + matmul(Pimp%screened_local(:,:,iw),Wimp%screened_local(:,:,iw))
         !
         call inv(invW)
         !call inv_sym(invW)
         !
         curlyU%screened_local(:,:,iw) = matmul(Wimp%screened_local(:,:,iw),invW)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(invW)
      call isReal(curlyU)
      !
      if(sym_)then
         do iw=1,Nmats
            call check_Symmetry(curlyU%screened_local(:,:,iw),eps,enforce=.true.,hardstop=.false.,name="curlyU_"//str(iw))
         enddo
      endif
      !
   end subroutine calc_curlyU


   !---------------------------------------------------------------------------!
   !PURPOSE: Computes the fully screened local interaction
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine calc_Wimp(Wimp,curlyU,ChiC)
      !
      use parameters
      use utils_fields
      use linalg, only : zeye
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wimp
      type(BosonicField),intent(in)         :: curlyU
      type(BosonicField),intent(in)         :: ChiC
      !
      complex(8),allocatable                :: Wtmp(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nmats
      integer                               :: iw
      !
      !
      if(verbose)write(*,"(A)") "---- calc_Wimp"
      !
      !
      ! Check on the input Fields
      if(.not.Wimp%status) stop "Wimp not properly initialized."
      if(.not.curlyU%status) stop "curlyU not properly initialized."
      if(.not.ChiC%status) stop "ChiC not properly initialized."
      if(Wimp%Nkpt.ne.0) stop "Wimp k dependent attributes are supposed to be unallocated."
      if(curlyU%Nkpt.ne.0) stop "curlyU k dependent attributes are supposed to be unallocated."
      if(ChiC%Nkpt.ne.0) stop "ChiC k dependent attributes are supposed to be unallocated."
      if(allocated(ChiC%bare_local))  stop "ChiC bare_local attribute is supposed to be unallocated."
      if(allocated(ChiC%bare))  stop "ChiC bare attribute is supposed to be unallocated."
      !
      Nbp = Wimp%Nbp
      Beta = Wimp%Beta
      Nmats = Wimp%Npoints
      !
      if(all([curlyU%Nbp-Nbp,ChiC%Nbp-Nbp].ne.[0,0])) stop "Either curlyU and/or ChiC have different orbital dimension with respect to Wimp."
      if(all([curlyU%Beta-Beta,ChiC%Beta-Beta].ne.[0d0,0d0])) stop "Either curlyU and/or ChiC have different Beta with respect to Wimp."
      if(all([curlyU%Npoints-Nmats,ChiC%Npoints-Nmats].ne.[0,0]))   stop "Either curlyU and/or ChiC have different number of Matsubara points with respect to Wimp."
      !
      call clear_attributes(Wimp)
      !
      Wimp%bare_local = curlyU%bare_local
      !
      allocate(Wtmp(Nbp,Nbp));Wtmp=czero
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nmats,Wimp,ChiC,curlyU),&
      !$OMP PRIVATE(iw,Wtmp)
      !$OMP DO
      do iw=1,Nmats
         !
         Wtmp = matmul(ChiC%screened_local(:,:,iw),curlyU%screened_local(:,:,iw))
         Wimp%screened_local(:,:,iw) = curlyU%screened_local(:,:,iw) - matmul(curlyU%screened_local(:,:,iw),Wtmp)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Wtmp)
      call isReal(Wimp)
      !
   end subroutine calc_Wimp


   !---------------------------------------------------------------------------!
   !PURPOSE: Correct the polarization at gamma
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine correct_Ugamma(Ulat)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(BosonicField),intent(inout)      :: Ulat
      !
      !
      if(verbose)write(*,"(A)") "---- correct_Ugamma"
      !
      !
      ! Check on the input Fields
      if(.not.Ulat%status) stop "Plat not properly initialized."
      if(Ulat%Nkpt.eq.0) stop "Plat k dependent attributes are supposed to be allocated."
      if(.not.allocated(Ugamma)) stop "Ugamma is not allocated."
      !
      call assert_shape(Ugamma,[Ulat%Nbp,Ulat%Nbp,Ulat%Npoints],"correct_Ugamma","Ugamma")
      !
      Ulat%screened(:,:,:,Ulat%iq_gamma) = Ugamma
      !
      call BosonicKsum(Ulat)
      !
   end subroutine correct_Ugamma


end module interactions
