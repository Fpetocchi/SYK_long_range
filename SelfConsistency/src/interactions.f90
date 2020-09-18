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
      module procedure read_U_spex_full                                         ![BosonicField,LocalOnly,save2bin,pathOUTPUT(optional to change output path),doAC(optional to override AC)]
      module procedure read_U_spex_Uloc0                                        ![Matrix,pathOUTPUT(optional to change output path)]
   end interface read_U_spex

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: calc_W_full
   public :: calc_W_edmft
   public :: calc_chi_full
   public :: calc_chi_edmft
   public :: read_U_spex
   !public :: build_Umatrix
   !public :: calc_QMCinteractions
   !public :: rescale_interaction

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Lattice inversion to get fully screened interaction - GW+EDMFT
   !---------------------------------------------------------------------------!
   subroutine calc_W_full(Wmats,Umats,Pmats,Lttc)
      !
      use parameters
      use utils_misc
      use utils_fields
      use linalg, only : zeye, inv_her !or sym??
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wmats
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: den(:,:)
      complex(8),allocatable                :: den_smallk(:,:,:,:)
      complex(8),allocatable                :: den_smallk_avrg(:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      integer                               :: ismall,num_k
      !
      !
      write(*,*) "--- calc_W_full ---"
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
      if(allocated(Lttc%small_ik)) stop "Kpoints near Gamma not stored. W divergence cannot be cured."
      !
      Nbp = Wmats%Nbp
      Nkpt = Wmats%Nkpt
      Beta = Wmats%Beta
      Nmats = Wmats%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "Either Umats and/or Pmats have different orbital dimension with respect to Wmats."
      if(all([Umats%Nkpt-Nkpt,Pmats%Nkpt-Nkpt].ne.[0,0])) stop "Either Umats and/or Pmats have different number of k-points with respect to Wmats."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "Either Umats and/or Pmats have different Beta with respect to Wmats."
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))  write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Wmats%Npoints,Umats%Npoints,Pmats%Npoints])
      !
      allocate(den_smallk(Nbp,Nbp,Nmats,12))
      allocate(den_smallk_avrg(Nbp,Nbp,Nmats))
      allocate(den(Nbp,Nbp))
      call clear_attributes(Wmats)
      !
      ! Assuming that the Polarization vanishes at iw__>inf
      Wmats%bare = Umats%bare
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Wmats,den_smallk,Lttc),&
      !$OMP PRIVATE(iq,iw,den,ismall)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if(iq.eq.Umats%iq_gamma)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            den=dcmplx(0d0,0d0)
            den = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened(:,:,iw,iq))
            !
            ! [ 1 - U*Pi ]^-1
            call inv_her(den)
            !
            ! [ 1 - U*Pi ]^-1 * U
            Wmats%screened(:,:,iw,iq) = matmul(den,Umats%screened(:,:,iw,iq))
            !
            do ismall=1,12
               if (Lttc%small_ik(ismall,1).eq.iq) den_smallk(:,:,iw,ismall) = den
            enddo
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(den)
      write(*,"(A,1I7)") "Gamma point is at kpoint list index: ",Umats%iq_gamma
      !
      !
      ! Gamma point handling
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
      !Replace the Gamma point value
      do iw=1,Nmats
         Wmats%screened(:,:,iw,Umats%iq_gamma) = matmul(den_smallk_avrg(:,:,iw),Umats%screened(:,:,iw,Umats%iq_gamma))
      enddo
      !
      deallocate(den_smallk,den_smallk_avrg)
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
      use linalg, only : zeye, inv_her !or sym??
      implicit none
      !
      type(BosonicField),intent(inout)      :: Wmats
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: den(:,:)
      complex(8),allocatable                :: den_smallk(:,:,:,:)
      complex(8),allocatable                :: den_smallk_avrg(:,:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      integer                               :: ismall,num_k
      !
      !
      write(*,*) "--- calc_W_edmft ---"
      !
      !
      ! Check on the input Fields
      if(.not.Wmats%status) stop "Wmats not properly initialized."
      if(.not.Umats%status) stop "Umats not properly initialized."
      if(.not.Pmats%status) stop "Pmats not properly initialized."
      if(Wmats%Nkpt.ne.0) stop "Wmats k dependent attributes are supposed to be unallocated."
      if(Umats%Nkpt.eq.0) stop "Umats k dependent attributes not properly initialized."
      if(Pmats%Nkpt.ne.0) stop "Pmats k dependent attributes are supposed to be unallocated."
      if(Umats%iq_gamma.lt.0) stop "Umats iq_gamma not defined."
      if(allocated(Lttc%small_ik)) stop "Kpoints near Gamma not stored. W divergence cannot be cured."
      !
      Nbp = Wmats%Nbp
      Nkpt = Umats%Nkpt
      Beta = Wmats%Beta
      Nmats = Wmats%Npoints
      !
      if(all([Umats%Nbp-Nbp,Pmats%Nbp-Nbp].ne.[0,0])) stop "Either Umats and/or Pmats have different orbital dimension with respect to Wmats."
      if(all([Umats%Beta-Beta,Pmats%Beta-Beta].ne.[0d0,0d0])) stop "Either Umats and/or Pmats have different Beta with respect to Wmats."
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))  write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Wmats%Npoints,Umats%Npoints,Pmats%Npoints])
      !
      allocate(den_smallk(Nbp,Nbp,Nmats,12))
      allocate(den_smallk_avrg(Nbp,Nbp,Nmats))
      allocate(den(Nbp,Nbp))
      call clear_attributes(Wmats)
      !
      ! Assuming that the Polarization vanishes at iw__>inf
      Wmats%bare_local = Umats%bare_local
      !
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Wmats,den_smallk,Lttc),&
      !$OMP PRIVATE(iq,iw,den,ismall)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if(iq.eq.Umats%iq_gamma)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            den=dcmplx(0d0,0d0)
            den = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened_local(:,:,iw))
            !
            ! [ 1 - U*Pi ]^-1
            call inv_her(den)
            !
            ! [ 1 - U*Pi ]^-1 * U
            Wmats%screened_local(:,:,iw) = Wmats%screened_local(:,:,iw) + matmul(den,Umats%screened(:,:,iw,iq))/Nkpt
            !
            do ismall=1,12
               if (Lttc%small_ik(ismall,1).eq.iq) den_smallk(:,:,iw,ismall) = den
            enddo
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(den)
      write(*,"(A,1I7)") "Gamma point is at kpoint list index: ",Umats%iq_gamma
      !
      !
      ! Gamma point handling
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
      !Add the Gamma point value
      do iw=1,Nmats
         Wmats%screened_local(:,:,iw) = Wmats%screened_local(:,:,iw) + matmul(den_smallk_avrg(:,:,iw),Umats%screened(:,:,iw,Umats%iq_gamma))/Nkpt
      enddo
      !
      deallocate(den_smallk,den_smallk_avrg)
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
      use linalg, only : zeye, inv_her !or sym??
      implicit none
      !
      type(BosonicField),intent(inout)      :: Chi
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: den(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      !
      !
      write(*,*) "--- calc_chi_full ---"
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
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))  write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Chi%Npoints,Umats%Npoints,Pmats%Npoints])
      !
      allocate(den(Nbp,Nbp))
      call clear_attributes(Chi)
      Chi%bare = Umats%bare
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Chi,Lttc),&
      !$OMP PRIVATE(iq,iw,den)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if(iq.eq.Umats%iq_gamma)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            den=dcmplx(0d0,0d0)
            den = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened(:,:,iw,iq))
            !
            ! [ 1 - U*Pi ]^-1
            call inv_her(den)
            !
            ! [ 1 - U*Pi ]^-1 * Pi
            Chi%screened(:,:,iw,iq) = matmul(den,Pmats%screened(:,:,iw,iq))
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
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
      use linalg, only : zeye, inv_her !or sym??
      implicit none
      !
      type(BosonicField),intent(inout)      :: Chi
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Pmats
      type(Lattice),intent(in)              :: Lttc
      !
      complex(8),allocatable                :: den(:,:)
      real(8)                               :: Beta
      integer                               :: Nbp,Nkpt,Nmats
      integer                               :: iq,iw
      !
      !
      write(*,*) "--- calc_chi_edmft ---"
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
      if(all([Umats%Npoints-Nmats,Pmats%Npoints-Nmats].ne.[0,0]))  write(*,"(A)") "Warning: Either Umats and/or Pmats have different number of Matsubara points. Computing up to the smaller."
      Nmats = minval([Chi%Npoints,Umats%Npoints,Pmats%Npoints])
      !
      allocate(den(Nbp,Nbp))
      call clear_attributes(Chi)
      Chi%bare_local = Umats%bare_local
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nbp,Nkpt,Nmats,Pmats,Umats,Chi,Lttc),&
      !$OMP PRIVATE(iq,iw,den)
      !$OMP DO
      do iq=1,Nkpt
         !
         !avoid the gamma point
         if(iq.eq.Umats%iq_gamma)cycle
         !
         do iw=1,Nmats
            !
            ! [ 1 - U*Pi ]
            den=dcmplx(0d0,0d0)
            den = zeye(Nbp) - matmul(Umats%screened(:,:,iw,iq),Pmats%screened_local(:,:,iw))
            !
            ! [ 1 - U*Pi ]^-1
            call inv_her(den)
            !
            ! [ 1 - U*Pi ]^-1 * Pi
            Chi%screened_local(:,:,iw) = Chi%screened_local(:,:,iw) + matmul(den,Pmats%screened_local(:,:,iw))/Nkpt
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine calc_chi_edmft


   !---------------------------------------------------------------------------!
   !PURPOSE: Read frequancy dependent interactions from SPEX files.
   !---------------------------------------------------------------------------!
   subroutine read_U_spex_full(Umats,save2bin,LocalOnly,pathOUTPUT,doAC)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use input_vars, only :  pathINPUT, UfullStructure
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      logical,intent(in)                    :: save2bin
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
      type(BosonicField)                    :: Ureal
      real                                  :: start,finish
      !
      !
      write(*,*) "--- read_U_spex_full ---"
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
      ! Read XEPS data
      !path = pathINPUT//"XEPS.DAT"
      !call inquireFile(reg(path),filexists)
      !call read_XEPS(reg(path))
      !
      !
      ! Check if the data on the Matsubara axis are present
      path = pathINPUT//"VW_imag"
      call inquireDir(reg(path),ACdone,hardstop=.false.)
      doAC_ = .not.ACdone
      if(present(doAC)) doAC_ = doAC
      !
      !
      ! Perform cnalytical continuation on real interaction or load existing files
      if(doAC_) then
         !
         !---------------------------------------------------------------------!
         !
         ! Allocations from dimensions written in W.Q0001.DAT file
         path = pathINPUT//"VW_real/VW.Q0001.DAT"
         call inquireFile(reg(path),filexists)
         unit = free_unit()
         open(unit,file=reg(path),form="unformatted",action="read")
         read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
         close(unit)
         !
         Nbp_spex = Norb_spex**2
         allocate(Utmp(Nbp_spex,Nbp_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         write(*,"(A,I5)")"Real frequencies: ",Nfreq
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Nbp_spex) stop "Size of given BosonicField and VW_real orbital space do not coincide."
         !
         ! Look for the Number of SPEX files. Which are supposed to be ordered.
         Nkpt = 0
         do iq=1,2000
            file_spex = reg(path)//"VW_real/VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists,hardstop=.false.)
            if(.not.filexists) exit
            Nkpt = Nkpt + 1
         enddo
         write(*,"(A1,1I6)") "The number of SPEX files (Nkpt) in VW_real is: ",Nkpt
         if((.not.LocalOnly).and.(Umats%Nkpt.ne.Nkpt)) stop "Number of k-points of given BosonicField and number of VW_real k-points do not coincide."
         !
         ! Allocate the Bosonic field on the real axis
         call AllocateBosonicField(Ureal,Nbp_spex,Nfreq,Nkpt,"Uspex(w)")
         !
         ! Read VW_real accumulating local attribute and optionally storing the k-dependent part
         path = pathINPUT//"VW_real/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"VW_real/VW.Q"//str(iq,4)//".DAT"
            call inquireFile(reg(file_spex),filexists) !redundant control
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
            write(*,*)"read iq",iq   !!!!>>>>>TEST<<<<<!!!!
            if (iq.ne.iqread) stop "iqread.ne.iq"
            !
            read(unit) wread
            wread = H2eV*wread
            write(*,*)"read wread",wread   !!!!>>>>>TEST<<<<<!!!!
            if (dabs(wread(1)).gt.eps) stop "wread(1) not zero"
            !
            do iw=0,Nfreq
               read(unit) Utmp
               if(iw.eq.0) then
                  !V(:,:,iq)=vwtmp(:,:)/Nkpt/Nkpt
                  Ureal%bare_local = Ureal%bare_local + H2eV*Utmp/(Nkpt**3)
                  if(.not.LocalOnly) Ureal%bare(:,:,iq) = H2eV*Utmp/(Nkpt**2)
                  !bare values on Matsubara are the same
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
         write(*,*) "UcRPA(w) --> UcRPA(iw) cpu timing:", finish-start
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
            write(*,*) "UcRPA(q,w) --> UcRPA(q,iw) cpu timing:", finish-start
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
            file_spex = reg(path)//"VW_imag/VW.Q"//str(iq,4)//".DAT"
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
            file_spex = reg(path)//"VW_imag/VW.Q"//str(iq,4)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"VW_imag/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),filexists) !redundant control
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit) iqread,Nspin_spex,Norb_spex,Nfreq
            write(*,*)"read iq",iq   !!!!>>>>>TEST<<<<<!!!!
            if (iq.ne.iqread) stop "iqread.ne.iq"
            !
            read(unit) wread
            wread = H2eV*wread
            write(*,*)"read wread",wread   !!!!>>>>>TEST<<<<<!!!!
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
      complex(8),allocatable,intent(inout)  :: Umat(:,:)
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
      write(*,*) "--- read_U_spex_Uloc0 ---"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Look for Uloc_mats.DAT and Uloc_real.DAT
      path=pathINPUT//"Uloc_mats.DAT"
      call inquireFile(reg(path),Umatsxists,hardstop=.false.)
      path=pathINPUT//"Uloc_real.DAT"
      call inquireFile(reg(path),Urealxists,hardstop=.false.)
      path=pathINPUT//"VW_real" !/VW.Q0001.DAT"
      call inquireDir(reg(path),SPEXxists,hardstop=.false.)      !
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
         path = pathINPUT//"VW_real/"
         do iq=1,2000
            !
            file_spex = reg(path)//"VW_real/VW.Q"//str(iq,4)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"VW_real/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),SPEXxists,hardstop=.false.)
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
   !---------------------------------------------------------------------------!
   subroutine checkAnalyticContinuation(Umats,Ureal)
      !
      use parameters
      use utils_misc
      use input_vars, only :  pathINPUT
      implicit none
      !
      type(BosonicField),intent(in)         :: Umats
      type(BosonicField),intent(in)         :: Ureal
      integer                               :: iq,ib1,ib2
      integer                               :: unit,Nbp,Nfreq,Nkpt
      real(8)                               :: ReErr,ImErr,thresh=1e-4
      !real(8)                               :: RUm_0,IUm_0,RUr_0,IUr_0
      !real(8)                               :: RUm_inf,IUm_inf,RUr_inf,IUr_inf
      real(8),allocatable                   :: ReErrMat(:,:),ImErrMat(:,:)
      !
      !
      write(*,*) "--- checkAnalyticContinuation ---"
      if(Umats%Nbp.ne.Ureal%Nbp) stop "Umats%Nbp.ne.Ureal%Nbp"
      if(Umats%Npoints.ne.Ureal%Npoints) stop "Umats%Npoints.ne.Ureal%Npoints"
      Nbp = Umats%Nbp
      Nfreq = Umats%Npoints
      allocate(ReErrMat(Nbp,Nbp));ReErrMat=0d0
      allocate(ImErrMat(Nbp,Nbp));ImErrMat=0d0
      !
      !
      ! Check the difference betqween bare values induced by thecutoff in the matsubara frequency
      unit = free_unit()
      open(unit=unit,file=reg(pathINPUT//"ACcutoffError.DAT"),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(A,1I5,A,1E14.7)")"Difference between Umats_bare value and last screened frequency: ",Nfreq," thresh:",thresh
      !
      ReErrMat=0d0;ImErrMat=0d0
      do ib1=1,Nbp
         do ib2=1,Nbp
            !
            ReErr = abs(real(Umats%bare_local(ib1,ib2)) - real(Umats%screened_local(ib1,ib2,Nfreq)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(aimag(Umats%bare_local(ib1,ib2)) - aimag(Umats%screened_local(ib1,ib2,Nfreq)))
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
                  ReErr = abs(real(Umats%bare(ib1,ib2,iq)) - real(Umats%screened(ib1,ib2,Nfreq,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(aimag(Umats%bare(ib1,ib2,iq)) - aimag(Umats%screened(ib1,ib2,Nfreq,iq)))
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
      open(unit=unit,file=reg(pathINPUT//"ACcheck.DAT"),form="formatted",status="unknown",position="rewind",action="write")
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
            ReErr = abs(real(Umats%screened_local(ib1,ib2,Nfreq)) - real(Ureal%screened_local(ib1,ib2,Nfreq)))
            if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
            !
            ImErr = abs(aimag(Umats%screened_local(ib1,ib2,Nfreq)) - aimag(Ureal%screened_local(ib1,ib2,Nfreq)))
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
                  ReErr = abs(real(Umats%screened(ib1,ib2,Nfreq,iq)) - real(Ureal%screened(ib1,ib2,Nfreq,iq)))
                  if(ReErr.gt.thresh) ReErrMat(ib1,ib2) = ReErr
                  !
                  ImErr = abs(aimag(Umats%screened(ib1,ib2,Nfreq,iq)) - aimag(Ureal%screened(ib1,ib2,Nfreq,iq)))
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
   !subroutine build_Umatrix(Umat,Uaa,Uab,Jsf,Jph)
   !end subroutine build_Umatrix



end module interactions
