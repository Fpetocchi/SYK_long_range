module interactions

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   interface read_spex
      module procedure read_spex_full                                           ![BosonicField,LocalOnly,save2bin,pathOUTPUT(optional to change output path),doAC(optional to override AC)]
      module procedure read_spex_Uloc0                                          ![Matrix,pathOUTPUT(optional to change output path)]
   end interface read_spex

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   integer,private                          :: Nspin_xeps
   integer,private                          :: Nkpt3_xeps(3)
   integer,private                          :: Nkpt_xeps
   integer,private                          :: Nkpt_xeps_irred
   integer,private                          :: Nband_xeps
   real(8),private                          :: Efermi_xeps
   logical,private                          :: UseDisentangledBS
   real(8),allocatable,private              :: Ene_xeps(:,:,:)
   real(8),allocatable,private              :: kpt_xeps(:,:)
   integer,allocatable,private              :: kptPos_xeps(:)
   integer,allocatable,private              :: kptPos(:)

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: read_spex
   !public :: build_Umatrix
   !public :: rescale_interaction

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Read XEPS.DAT file
   !---------------------------------------------------------------------------!
   subroutine read_xeps(path2xeps)
      !
      use global_vars,                 only :  Nkpt,Nkpt3,kpt,Nkpt_irred
      use global_vars,                 only :  UseXepsKorder
      use utils_misc
      implicit none
      !
      character(len=*),intent(in),optional  :: path2xeps
      integer                               :: ik,unit
      logical                               :: dumlogical
      !
      !
      write(*,*) "--- read_xeps ---"
      !
      !
      unit = free_unit()
      open(unit,file=reg(path2xeps),form="unformatted",action="read")
      read(unit) Nspin_xeps
      read(unit) Nkpt3_xeps
      read(unit) Nkpt_xeps
      read(unit) Nkpt_xeps_irred
      read(unit) Nband_xeps
      read(unit) Efermi_xeps
      read(unit) dumlogical
      read(unit) UseDisentangledBS
      !
      allocate(kpt_xeps(3,Nkpt_xeps));kpt_xeps=0d0
      allocate(kptPos_xeps(Nkpt_xeps));kptPos_xeps=0
      allocate(Ene_xeps(Nband_xeps,Nkpt_xeps_irred,Nspin_xeps));Ene_xeps=0d0
      !
      read(unit) kpt_xeps
      read(unit) kptPos_xeps
      read(unit) Ene_xeps
      !
      close(unit)
      !
      Nkpt_irred = Nkpt
      if(UseXepsKorder) Nkpt_irred = Nkpt_xeps_irred
      !
      ! Global checks
      if(Nspin_xeps.ne.1)           stop "Nspin_xeps.ne.1 in XEPS.DAT"
      if(Nkpt_xeps.ne.Nkpt)         stop "Nkpt_xeps.ne.Nkpt in XEPS.DAT"
      if(Nkpt3_xeps(1).ne.Nkpt3(1)) stop "Nkpt(1)_xeps.ne.Nkpt(1) in XEPS.DAT"
      if(Nkpt3_xeps(2).ne.Nkpt3(2)) stop "Nkpt(2)_xeps.ne.Nkpt(2) in XEPS.DAT"
      if(Nkpt3_xeps(3).ne.Nkpt3(3)) stop "Nkpt(3)_xeps.ne.Nkpt(3) in XEPS.DAT"
      !
      ! Check of the K-point ordering
      do ik=1,Nkpt
         if (.not.keq(kpt_xeps(:,ik),kpt(:,ik))) then
            write(*,*)"ik=",ik,"kpt(:,ik)=",kpt(:,ik),"kpt_loc(:,ik=)",kpt_xeps(:,ik)
            !write(*,*) "kptp(ik)=",kptPos(ik),"kptp_loc(ik)=",kptPos_xeps(ik)
            stop "K-points grid does not match"
         endif
      enddo
      !
      allocate(kptPos(Nkpt));kptPos=0
      if(UseXepsKorder)then
         kptPos=kptPos_xeps
      else
         do ik=1,Nkpt
            kptPos(ik)=ik
         enddo
      endif
      !
   end subroutine read_xeps


   !---------------------------------------------------------------------------!
   !PURPOSE: Read frequancy dependent interactions
   !---------------------------------------------------------------------------!
   subroutine read_spex_full(Umats,LocalOnly,save2bin,pathOUTPUT,doAC)
      !
      use parameters
      use file_io
      use global_vars,                 only :  pathINPUT,UfullStructure
      use utils_misc
      use utils_fields
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      logical,intent(in)                    :: LocalOnly
      logical,intent(in)                    :: save2bin
      character(len=*),intent(in),optional  :: pathOUTPUT
      logical,intent(in),optional           :: doAC
      !
      logical                               :: filexists,ACdone,doAC_
      character(len=256)                    :: file_spex,path,pathOUTPUT_
      integer                               :: unit,Nkpt
      integer                               :: iq,iw,iqread,Npb_spex
      integer                               :: idum,Nspin_spex,Norb_spex,Nfreq
      integer                               :: ib1,ib2,iw1,iw2
      real(8),allocatable                   :: wread(:),wm_B(:)
      complex(8),allocatable                :: D1(:,:),D2(:,:),D3(:,:),imgFact(:,:,:)
      complex(8),allocatable                :: Utmp(:,:)
      type(BosonicField)                    :: Ureal
      real                                  :: start,finish
      !
      !
      write(*,*) "--- read_spex_full ---"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      !
      ! Check on the input Boson
      if(.not.Umats%status) stop "BosonicField not properly initialized."
      allocate(wm_B(Umats%Npoints));wm_B=0d0
      wm_B = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
      if(.not.LocalOnly)Nkpt = Umats%Nkpt
      !
      !
      ! Read XEPS data
      path = pathINPUT//"XEPS.DAT"
      call inquireFile(reg(path),filexists)
      call read_XEPS(reg(path))
      !
      !
      ! Check if the data on the Matsubara axis are present
      path = pathINPUT//"VW_imag" !/VW.Q0001.DAT"
      call inquireDir(reg(path),ACdone)
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
         Npb_spex = Norb_spex**2
         allocate(Utmp(Npb_spex,Npb_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         write(*,"(A,I5)")"Real frequencies: ",Nfreq
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Npb_spex) stop "Size of given BosonicField and VW_real orbital space do not coincide."
         if((.not.LocalOnly).and.(.not.allocated(Umats%bare))) stop "Requested K-dependence but bare attribute not allocated."
         if((.not.LocalOnly).and.(.not.allocated(Umats%screened))) stop "Requested K-dependence but screened attribute not allocated."
         !
         ! Allocate the Bosonic field on the real axis
         Ureal%Nbp = Npb_spex
         Ureal%Npoints = Nfreq
         if(.not.LocalOnly)Ureal%Nkpt = Nkpt
         call selfAllocateBosonicField(Ureal)
         !
         ! Read VW_real accumulating local attribute and optionally storing the k-dependent part
         path = pathINPUT//"VW_real/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"GWinput/VW_real/VW.Q"//str(iq)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"GWinput/VW_real/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),filexists)
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
            if (dabs(wread(1)).gt.1.d-9) stop "wread(1) not zero"
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
         allocate(D1(Npb_spex,Npb_spex));D1=czero
         allocate(D2(Npb_spex,Npb_spex));D2=czero
         allocate(D3(Npb_spex,Npb_spex));D3=czero
         !
         !
         ! Analytical continuation of the local component to imag axis using spectral rep
         call cpu_time(start)
         !
         ! Check if any local Urpa component has inverted Im/Re symmetry
         do ib1=1,Npb_spex
            do ib2=1,Npb_spex
               if( abs(real(Ureal%bare_local(ib1,ib2))).lt.abs(aimag(Ureal%bare_local(ib1,ib2))))then
                  write(*,"(A,2I5)")"Element: ",ib1,ib2
                  write(*,"(A,E10.4)")"Re[Ubare(w=inf)]: ",real(Ureal%bare_local(ib1,ib2))
                  write(*,"(A,E10.4)")"Im[Ubare(w=inf)]: ",aimag(Ureal%bare_local(ib1,ib2))
                  stop "Something wrong: Uloc cannot have inverted Re/Im parity."
               endif
            enddo
         enddo
         !
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Npb_spex,wm_B,Nw_B,wread,Nfreq,Ureal,Umats),&
         !$OMP PRIVATE(ib1,ib2,iw1,iw2,D1,D2,D3,Utmp)
         !$OMP DO
         do iw1=1,Umats%Npoints
            Utmp=czero
            do iw2=1,Nfreq-2,2
               !
               do ib1=1,Npb_spex
                  do ib2=1,Npb_spex
                     D1(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2)   )/pi
                     D2(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+1) )/pi
                     D3(ib1,ib2) = -dimag( Ureal%screened_local(ib1,ib2,iw2+2) )/pi
                  enddo
               enddo
               !
               !D(-w)=-D(w), integrate using Simpson method
               if(wread(iw2).gt.0.d0) then
                  Utmp(:,:) = Utmp(:,:) + ( D1(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2)  ) - D1(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2)  ) ) *(wread(iw2+1)-wread(iw2))/3.d0
                  Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                  Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
               elseif(dabs(wread(iw2)).lt.1.d-12) then
                  Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                  Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
               endif
            enddo
            !
            do ib1=1,Npb_spex
               do ib2=1,Npb_spex
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
            allocate(D1(Npb_spex,Npb_spex));D1=czero
            allocate(D2(Npb_spex,Npb_spex));D2=czero
            allocate(D3(Npb_spex,Npb_spex));D3=czero
            !
            !
            ! Analytical continuation of all the K-points to imag axis using spectral rep
            allocate(imgFact(Npb_spex,Npb_spex,2));imgFact=cone
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(Nk,Npb_spex,wm_B,wread,Nfreq,Ureal,Umats,UfullStructure),&
            !$OMP PRIVATE(iq,ib1,ib2,iw1,iw2,D1,D2,D3,Utmp,imgFact)
            !$OMP DO
            do iq=1,Umats%Nkpt
               !
               ! Some elelments of U, usually the k-dependent one, have inverted Im/Re symmetry
               imgFact=cone
               if(UfullStructure)then
                  do ib1=1,Npb_spex
                     do ib2=1,Npb_spex
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
                     do ib1=1,Npb_spex
                        do ib2=1,Npb_spex
                           D1(ib1,ib2) = -dimag( ( imgFact(ib1,ib2,1) * Ureal%screened(ib1,ib2,iw2,iq)   ) )/pi
                           D2(ib1,ib2) = -dimag( ( imgFact(ib1,ib2,1) * Ureal%screened(ib1,ib2,iw2+1,iq) ) )/pi
                           D3(ib1,ib2) = -dimag( ( imgFact(ib1,ib2,1) * Ureal%screened(ib1,ib2,iw2+2,iq) ) )/pi
                        enddo
                     enddo
                     !
                     !D(-w)=-D(w), integrate using Simpson method
                     if(wread(iw2).gt.0.d0) then
                        Utmp(:,:) = Utmp(:,:) + ( D1(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2)  ) - D1(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2)  ) ) *(wread(iw2+1)-wread(iw2))/3.d0
                        Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                        Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
                     elseif(dabs(wread(iw2)).lt.1.d-12) then
                        Utmp(:,:) = Utmp(:,:) + ( D2(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+1)) - D2(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+1)) ) *(wread(iw2+1)-wread(iw2))*4.d0/3.d0
                        Utmp(:,:) = Utmp(:,:) + ( D3(:,:)/(dcmplx(0.d0,wm_B(iw1))-wread(iw2+2)) - D3(:,:)/(dcmplx(0.d0,wm_B(iw1))+wread(iw2+2)) ) *(wread(iw2+1)-wread(iw2))/3.d0
                     endif
                  enddo
                  !
                  do ib1=1,Npb_spex
                     do ib2=1,Npb_spex
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
         call selfDeallocateBosonicField(Ureal)
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
         Npb_spex = Norb_spex**2
         allocate(Utmp(Npb_spex,Npb_spex));Utmp=czero
         allocate(wread(Nfreq));wread=0d0
         write(*,"(A,I5)")"Matsubara frequencies: ",Nfreq
         !
         ! Few checks
         if(Nspin_spex.ne.1) stop "Nspin_spex.ne.1"
         if(Umats%Nbp.ne.Npb_spex) stop "Size of given BosonicField and VW_imag orbital space do not coincide."
         if((.not.LocalOnly).and.(.not.allocated(Umats%bare))) stop "Requested K-dependence but bare attribute not allocated."
         if((.not.LocalOnly).and.(.not.allocated(Umats%screened))) stop "Requested K-dependence but screened attribute not allocated."
         if(Umats%Npoints.ne.Nfreq) stop "Number of Matsubara points and bosonic field mesh does not coincide."
         !
         ! Read VW_real accumulating local attribute and optionally storing the k-dependent part
         path = pathINPUT//"VW_imag/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"GWinput/VW_imag/VW.Q"//str(iq)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"GWinput/VW_imag/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),filexists)
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
               if (dabs(wread(iw)-wm_B(iw)).gt.1.d-6) stop "wread.ne.wm_B"
            enddo
            !
            do iw=0,Nfreq
               read(unit) Utmp
               if(iw.eq.0) then
                  Umats%bare_local = Umats%bare_local + H2eV*Utmp/(Nkpt**3)
                  Umats%bare(:,:,iq) = H2eV*Utmp/(Nkpt**2)
               else
                  Umats%screened_local(:,:,iw) = Umats%screened_local(:,:,iw) + H2eV*Utmp/(Nkpt**3)
                  Umats%screened(:,:,iw,iq) = H2eV*Utmp/(Nkpt**2)
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
            do ib1=1,Npb_spex
               do ib2=1,Npb_spex
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
   end subroutine read_spex_full


   !---------------------------------------------------------------------------!
   !PURPOSE: Read screened interaction tensor from Ucrpa(0)
   !---------------------------------------------------------------------------!
   subroutine read_spex_Uloc0(Umat,pathOUTPUT)
      !
      use parameters
      use file_io
      use global_vars,                 only :  Nkpt
      use global_vars,                 only :  pathINPUT
      use utils_misc
      use utils_fields
      implicit none
      !
      complex(8),allocatable,intent(inout)  :: Umat(:,:)
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      logical                               :: Umatsxists,Urealxists,XEPSxists
      character(len=256)                    :: file_spex,path,pathOUTPUT_
      integer                               :: unit
      integer                               :: iq,iw,Npb_spex
      integer                               :: iqread,Nspin_spex,Norb_spex,Nfreq
      complex(8),allocatable                :: Utmp(:,:)
      type(BosonicField)                    :: Uread
      !
      !
      write(*,*) "--- read_spex_Uloc0 ---"
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
      call inquireDir(reg(path),XEPSxists,hardstop=.false.)      !
      !
      if(Umatsxists)then
         !
         Uread%status = .true.
         Uread%Nbp = size(Umat,dim=1)
         Uread%Npoints = 2
         call selfAllocateBosonicField(Uread)
         !
         call read_BosonicField(Uread,reg(pathINPUT),"Uloc_mats.DAT")
         Umat = Uread%screened_local(:,:,1)
         call dump_matrix(Umat,reg(pathINPUT//"Umat.DAT"))
         call selfDeallocateBosonicField(Uread)
         return
         !
      else if(Urealxists)then
         !
         Uread%status = .true.
         Uread%Nbp = size(Umat,dim=1)
         Uread%Npoints = 2
         call selfAllocateBosonicField(Uread)
         !
         call read_BosonicField(Uread,reg(pathINPUT),"Uloc_real.DAT")
         Umat = Uread%screened_local(:,:,1)
         call dump_matrix(Umat,reg(pathINPUT//"Umat.DAT"))
         call selfDeallocateBosonicField(Uread)
         return
         !
      else if(XEPSxists)then
         !
         path = pathINPUT//"VW_real/"
         do iq=1,Nkpt
            !
            file_spex = reg(path)//"GWinput/VW_real/VW.Q"//str(iq)//".DAT"        !write(fn,"(a,a,i4.4,a)") reg(path),"GWinput/VW_real/VW.Q",iq,".DAT"
            call inquireFile(reg(file_spex),XEPSxists)
            !
            unit = free_unit()
            open(unit,file=reg(file_spex),form="unformatted",action="read")
            read(unit)iqread,Nspin_spex,Norb_spex,Nfreq
            !
            Npb_spex = Norb_spex**2
            allocate(Utmp(Npb_spex,Npb_spex));Utmp=czero
            if(iq.eq.1)call assert_shape(Umat,[Npb_spex,Npb_spex],"read_spex_Uloc","Umat")
            !
            read(unit) !wread
            !
            do iw=0,1
               read(unit) Utmp
               if(iw.eq.1) then
                  Umat = Umat + H2eV*Utmp/(Nkpt**3)
               endif
            enddo
            !
            close(unit)
            deallocate(Utmp)
            !
         enddo !iq
         call dump_matrix(Umat,reg(reg(pathOUTPUT_)//"Umat.DAT"))
         return
         !
      else
         stop "No useful interaction file found."
      endif
      !
   end subroutine read_spex_Uloc0


   !---------------------------------------------------------------------------!
   !PURPOSE: Check if the AC alters the bare and screened values
   !---------------------------------------------------------------------------!
   subroutine checkAnalyticContinuation(Umats,Ureal)
      !
      use parameters
      use global_vars,                 only :  pathINPUT
      use utils_misc
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
      write(unit,"(A,1I5,A,1E6.2)")"Difference between Umats_bare value and last screened frequency: ",Nfreq," thresh:",thresh
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
         write(unit,"(999E10.5)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection"
      do ib1=1,Nbp
         write(unit,"(999E10.5)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
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
               write(unit,"(999E10.5)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E10.5)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
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
      write(unit,"(A,1E6.2)")"Difference between asymptotic behaviour of Ureal and Umats. Thresh:",thresh
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
         write(unit,"(999E10.5)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection - screened limit"
      do ib1=1,Nbp
         write(unit,"(999E10.5)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
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
         write(unit,"(999E10.5)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
      enddo
      write(unit,"(A)")"Imag part - local projection - bare limit"
      do ib1=1,Nbp
         write(unit,"(999E10.5)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
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
               write(unit,"(999E10.5)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - screened limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E10.5)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
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
               write(unit,"(999E10.5)")(ReErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            write(unit,"(A,1I5)")"Imag part - bare limit - iq: ",iq
            do ib1=1,Nbp
               write(unit,"(999E10.5)")(ImErrMat(ib1,ib2),ib2=1,Nbp)
            enddo
            !
         enddo !iq
      endif
      close(unit)
      !
   end subroutine checkAnalyticContinuation


   !---------------------------------------------------------------------------!
   !PURPOSE: Create the interaction tensor from user-given values of U and J
   !---------------------------------------------------------------------------!
   !subroutine build_Umatrix(Umat,Uaa,Uab,Jsf,Jph)
   !end subroutine build_Umatrix



end module interactions
