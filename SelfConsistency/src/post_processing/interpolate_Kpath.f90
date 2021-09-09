subroutine interpolateG2Path(Sfull,Lttc,pathOUTPUT)
   !
   use parameters
   use utils_misc
   use utils_fields
   use linalg, only : eigh, inv, zeye, det
   use crystal
   use file_io
   use greens_function, only : calc_Gmats
   use fourier_transforms
   use input_vars, only : structure, path_funct, Nkpt_path, FermiSurf, Nkpt_Fermi
   use input_vars, only : paramagnet, CalculationType, Hetero
   implicit none
   !
   type(FermionicField),intent(inout)    :: Sfull
   type(Lattice),intent(inout)           :: Lttc
   character(len=*),intent(in)           :: pathOUTPUT
   !
   type(FermionicField)                  :: Spath
   type(FermionicField)                  :: Sfermi
   type(FermionicField)                  :: Sloc
   type(FermionicField)                  :: Gpath
   type(FermionicField)                  :: Gfull
   type(FermionicField)                  :: Gfermi
   !
   real(8),allocatable                   :: Zk(:,:)
   complex(8),allocatable                :: correction(:,:,:)
   integer                               :: Norb,Nmats,unit!,Ntau
   integer                               :: ik,iw,itau,ispin,iorb
   integer                               :: ikx,iky
   real(8)                               :: kx,ky,Bvec(3),Blat(3,3)
   character(len=256)                    :: path
   logical                               :: Kdependence
   real                                  :: start,finish
   !
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- interpolateG2Path"
   !
   !
   if(.not.Lttc%status) stop "interpolateG2Path: Lttc not properly initialized."
   if(.not.Sfull%status) stop "interpolateG2Path: Sfull not properly initialized."
   !
   if(Sfull%Norb.ne.Lttc%Norb) stop "interpolateG2Path: Lttc has different number of orbitals with respect to Sfull."
   if(Sfull%Nkpt.ne.Lttc%Nkpt) stop "interpolateG2Path: Lttc has different number of k-points with respect to Sfull."
   Norb = Sfull%Norb
   Nmats = Sfull%Npoints
   !
   if(FermiSurf)call get_Blat(Blat)
   !
   !
   !---------------------- path along high-symmetry points --------------------!
   !
   !
   !non-interacting data (Bands, spectral function, Fermi-surface)
   call interpolateHk2Path(Lttc,structure,Nkpt_path,reg(pathOUTPUT)//"K_resolved/",doplane=FermiSurf,hetero=Hetero%status)
   !
   !Dump MaxEnt data for the local self-energy (always done for every setup)
   call calc_MaxEnt_on_Sigma_imp(Sfull)
   !
   !
   !Interpolate the slef-energy along the path if its K-dependent otherwise duplicate the local one
   select case(reg(CalculationType))
      case default
         !
         stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
         !
      case("G0W0","scGW","GW+EDMFT")
         !
         !
         !
         Kdependence = .true.
         !
         !
         !
         !--------------- Green's function in the full BZ ---------------!
         if(scan(reg(path_funct),"S").gt.0)then
            call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
         endif
         !
         !
         !
         !--------------- Green's function along the path ---------------!
         !
         !Interpolate the self-energy along the path
         call cpu_time(start)
         call AllocateFermionicField(Spath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
         do ispin=1,Nspin
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),Sfull%wks(:,:,:,:,ispin),Spath%wks(:,:,:,:,ispin))
            if(paramagnet)then
               Spath%wks(:,:,:,:,Nspin) = Spath%wks(:,:,:,:,1)
               exit
            endif
         enddo
         call cpu_time(finish)
         write(*,"(A,F)") new_line("A")//new_line("A")//"     Sigma(fullBZ,iw) --> Sigma(Kpath,iw) cpu timing:", finish-start
         !
         !Compute the quasiparticle weight along the path
         do ispin=1,Nspin
            allocate(Zk(Lttc%Nkpt_path,Norb));Zk=0d0
            do ik=1,Lttc%Nkpt_path
               do iorb=1,Norb
                  Zk(ik,iorb) = 1d0 / (1d0 + abs(dimag(Spath%wks(iorb,iorb,1,ik,ispin)))*Spath%Beta/pi)
               enddo
            enddo
            path = reg(pathOUTPUT)//"K_resolved/Zk_path_s"//str(ispin)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,Lttc%Nkpt_path
               write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik)/Lttc%Kpathaxis(Lttc%Nkpt_path),(Zk(ik,iorb),iorb=1,Norb)
            enddo
            close(unit)
            deallocate(Zk)
            if(paramagnet)exit
         enddo
         !
         !Recompute the Green's function
         call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
         call calc_Gmats(Gpath,Lttc,Smats=Spath,along_path=.true.)
         call DeallocateFermionicField(Spath)
         !
         !
         !
         !---------- Green's function along the {kx,ky} plane -----------!
         if(FermiSurf)then
            !
            !
            !Interpolate the self-energy along the plane
            call cpu_time(start)
            call AllocateFermionicField(Sfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            do ispin=1,Nspin
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,Sfull%wks(:,:,:,:,ispin),Sfermi%wks(:,:,:,:,ispin))
               if(paramagnet)then
                  Sfermi%wks(:,:,:,:,Nspin) = Sfermi%wks(:,:,:,:,1)
                  exit
               endif
            enddo
            call cpu_time(finish)
            write(*,"(A,F)") "     Sigma(fullBZ,iw) --> Sigma(kx,ky,iw) cpu timing:", finish-start
            !
            !Compute the quasiparticle weight along the plane - this part can be optimized by avoiding to interpolate the full frequency range of S_Full
            do ispin=1,Nspin
               allocate(Zk(Lttc%Nkpt_Plane,Norb));Zk=0d0
               do ik=1,Lttc%Nkpt_Plane
                  do iorb=1,Norb
                     Zk(ik,iorb) = 1d0 / (1d0 + abs(dimag(Sfermi%wks(iorb,iorb,1,ik,ispin)))*Sfermi%Beta/pi)
                  enddo
               enddo
               path = reg(pathOUTPUT)//"K_resolved/Zk_plane_s"//str(ispin)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_Plane
                  ikx = int(ik/(Nkpt_Fermi+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Fermi-1) - 0.5d0
                  iky = ik - (ikx-1)*Nkpt_Fermi      ; ky = (iky-1)/dble(Nkpt_Fermi-1) - 0.5d0
                  Bvec = kx*Blat(:,1) + ky*Blat(:,2)
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Zk(ik,iorb),iorb=1,Norb)
                  if(iky.eq.Nkpt_Fermi)write(unit,*)
               enddo
               close(unit)
               deallocate(Zk)
               if(paramagnet)exit
            enddo
            !
            !Recompute the Green's function
            call AllocateFermionicField(Gfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gfermi,Lttc,Smats=Sfermi,along_plane=.true.)
            call DeallocateFermionicField(Sfermi)
            !
         endif
         !
         !
         !
      case("DMFT+statU","DMFT+dynU","EDMFT")
         !
         !
         Kdependence = .false.
         !
         !
         !
         call AllocateFermionicField(Sloc,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
         do ik=1,Lttc%Nkpt_path
            Sloc%wks(:,:,:,ik,:) = Sfull%ws
         enddo
         !
         !
         !
         !--------------- Green's function along the path ---------------!
         call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
         call calc_Gmats(Gpath,Lttc,Smats=Sloc,along_path=.true.)
         !
         !
         !
         !---------- Green's function along the {kx,ky} plane -----------!
         if(FermiSurf)then
            call AllocateFermionicField(Gfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gfermi,Lttc,Smats=Sloc,along_plane=.true.)
         endif
         call DeallocateFermionicField(Sloc)
         !
         !
         !
   end select
   !
   !correction to LDA given by the real part of the local self-energy in iw=0
   allocate(correction(Norb,Norb,Sfull%Nkpt));correction=czero
   do ispin=1,Nspin
      !
      do ik=1,Sfull%Nkpt
         correction(:,:,ik) = zeye(Norb)*Sfull%mu - dreal(Sfull%ws(:,:,1,ispin))
      enddo
      call interpolateHk2Path(Lttc,structure,Nkpt_path,reg(pathOUTPUT)//"K_resolved/",corrname="dmft_s"//str(ispin),correction=correction,doplane=FermiSurf,hetero=Hetero%status)
      if(paramagnet) exit
      !
   enddo
   deallocate(correction)
   !
   !correction to LDA given by the non-local self-energy in iw=0 with scattering rates removed
   if(Kdependence)then
      allocate(correction(Norb,Norb,Sfull%Nkpt));correction=czero
      do ispin=1,Nspin
         !
         do ik=1,Sfull%Nkpt
            correction(:,:,ik) = zeye(Norb)*Sfull%mu - Sfull%wks(:,:,1,ik,ispin)
         enddo
         do iorb=1,Norb
            correction(iorb,iorb,:) = dcmplx(dreal(correction(iorb,iorb,:)),0d0)
         enddo
         call interpolateHk2Path(Lttc,structure,Nkpt_path,reg(pathOUTPUT)//"K_resolved/",corrname="qpsc_s"//str(ispin),correction=correction,doplane=FermiSurf,hetero=Hetero%status)
         if(paramagnet) exit
         !
      enddo
      deallocate(correction)
   endif
   !
   !
   if(scan(reg(path_funct),"G").gt.0)then
     !
     !Dump K-resolved MaxEnt data along the path
     call calc_MaxEnt_on_G_K(Gpath,"path")
     !
     !Dump K-resolved MaxEnt data along the {kx,ky} plane
     if(FermiSurf)call calc_MaxEnt_on_G_K(Gfermi,"plane")
     !
   endif
   !
   !
   if(scan(reg(path_funct),"S").gt.0)then
     !
     if(Kdependence)then
         !
         !Dump K-resolved MaxEnt data in the full BZ
         call calc_MaxEnt_on_Sigma_K(Gfull,"full")
         !
         !Dump K-resolved MaxEnt data along the path
         call calc_MaxEnt_on_Sigma_K(Gpath,"path")
         !
         !Dump K-resolved MaxEnt data along the {kx,ky} plane
         call calc_MaxEnt_on_Sigma_K(Gfermi,"plane")
         !
     endif
     !
   endif
   !
   call DeallocateFermionicField(Gpath)
   call DeallocateFermionicField(Gfull)
   call DeallocateFermionicField(Gfermi)
   !
   !
contains
   !
   !
   !
   subroutine calc_MaxEnt_on_G_K(Gmats_in,mode)
      !
      use input_vars, only : Ntau
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats_in
      character(len=*),intent(in)           :: mode
      !
      complex(8),allocatable                :: Gmats_diag(:,:,:,:),Gitau_diag(:,:,:,:)
      real(8),allocatable                   :: Ak(:,:)
      real(8),allocatable                   :: tau(:)
      integer                               :: Nkpt!Ntau,Nmats_cutoff
      integer                               :: ikx,iky
      !
      !
      if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- calc_MaxEnt_on_G_K"
      !
      !
      if(.not.Gmats_in%status) stop "calc_MaxEnt_on_G_K: Gmats_in not properly allocated."
      select case(reg(mode))
         case default
            !
            stop "calc_MaxEnt_on_G_K: Available Modes are: path, full, plane."
            !
         case("path")
            !
            Nkpt = Lttc%Nkpt_path
            write(*,"(A,I)") "     G path. Total number of K-points along path:",Nkpt
            !
         case("plane")
            !
            Nkpt = Lttc%Nkpt_Plane
            write(*,"(A,I)") "     G plane. Total number of K-points in the {kx,ky} sheet:",Nkpt
            !
      end select
      !
      !Extract the diagonal of the Green's function
      allocate(Gmats_diag(Norb,Nmats,Nkpt,Nspin));Gmats_diag=czero
      do ispin=1,Nspin
         do ik=1,Nkpt
            do iw=1,Nmats
               do iorb=1,Norb
                  Gmats_diag(iorb,iw,ik,ispin) = Gmats_in%wks(iorb,iorb,iw,ik,ispin)
               enddo
            enddo
         enddo
      enddo
      !
      !Fourier transform the diagonal of the Green's function
      !Ntau = NtauF  <-- should I decrease this guy to ease MaxEnt?
      call cpu_time(start)
      allocate(Gitau_diag(Norb,Ntau,Nkpt,Nspin));Gitau_diag=czero
      do ispin=1,Nspin
         call Fmats2itau_vec(Sfull%Beta,Gmats_diag(:,:,:,ispin),Gitau_diag(:,:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
         if(paramagnet)then
            Gitau_diag(:,:,:,Nspin) = Gitau_diag(:,:,:,1)
            exit
         endif
      enddo
      deallocate(Gmats_diag)
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(K"//reg(mode)//",iw) --> Glat(K"//reg(mode)//",tau) cpu timing:", finish-start
      !
      !Print data for K-resolved MaxEnt
      allocate(tau(Ntau));tau = linspace(0d0,Sfull%Beta,Ntau)
      do ispin=1,Nspin
        do ik=1,Nkpt
            !
            path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Gk_"//reg(mode)//"_t_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do itau=1,Ntau
                !write(unit,"(200E20.12)") tau(itau),(dreal(Gitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
                write(unit,"(200E20.12)") tau(itau),(-abs(Gitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
            enddo
            close(unit)
            !
        enddo
        if(paramagnet)exit
      enddo
      deallocate(tau)
      !
      !Compute the spectral weight at Fermi along the path. See arxiv:0805.3778 Eq.(5)
      do ispin=1,Nspin
         !
         allocate(Ak(Nkpt,Norb));Ak=0d0
         do ik=1,Nkpt
            do iorb=1,Norb
               Ak(ik,iorb) = -dreal(Gitau_diag(iorb,int(Ntau/2),ik,ispin))*Gmats_in%Beta
            enddo
         enddo
         !
         path = reg(pathOUTPUT)//"K_resolved/Ak_"//reg(mode)//"_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         if(reg(mode).eq."path")then
            do ik=1,Nkpt
               write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik)/Lttc%Kpathaxis(Lttc%Nkpt_path),(Ak(ik,iorb),iorb=1,Norb)
            enddo
         elseif(reg(mode).eq."plane")then
            do ik=1,Nkpt
               ikx = int(ik/(Nkpt_Fermi+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Fermi-1) - 0.5d0
               iky = ik - (ikx-1)*Nkpt_Fermi      ; ky = (iky-1)/dble(Nkpt_Fermi-1) - 0.5d0
               Bvec = kx*Blat(:,1) + ky*Blat(:,2)
               write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Ak(ik,iorb),iorb=1,Norb)
               if(iky.eq.Nkpt_Fermi)write(unit,*)
            enddo
         endif
         close(unit)
         deallocate(Ak)
         !
         if(paramagnet)exit
      enddo
      deallocate(Gitau_diag)
      !





      if(Hetero%status)then
         !
 !        Norb_layer = Hetero%Norb
 !        if(Norb_layer.ne.int(Lttc%Norb/Lttc%Nsite)) stop "calc_MaxEnt_on_G_K: wrong hetero orbital dimension."
 !        !
 !        allocate(Gmats_kpkz(Norb_layer,Norb_layer,Nmats,Nkpt,Nkpt_path,Nspin));Gmats_kpkz=czero
 !        do ispin=1,Nspin
 !           do iw=1,Nmats
 !              call fill_Gamma_A(Gmats_in%wks(:,:,iw,:,ispin),Gmats_kpkz(:,:,iw,:,:,ispin))
 !           enddo
 !           if(paramagnet)then
 !              Gmats_kpkz(:,:,:,:,:,Nspin) = Gmats_kpkz(:,:,:,:,:,1)
 !              exit
 !           endif
 !        enddo
 !        !



      endif
      !
   end subroutine calc_MaxEnt_on_G_K
   !
   !
   !
   subroutine fill_Gamma_A(data_in,data_out)
      !
      implicit none
      !
      complex(8),intent(in)              :: data_in(:,:,:)
      complex(8),intent(out)             :: data_out(:,:,:,:)
      !
      integer                            :: ra,rb,ca,cb,ikz
      integer                            :: isite,jsite
      integer                            :: Nkpt_layer,Norb_layer
      real(8)                            :: kR
      complex(8)                         :: cfac
      !
      Nkpt_layer = size(data_in,dim=3)
      Norb_layer = size(data_out,dim=1)
      !
      !$OMP PARALLEL DEFAULT(PRIVATE),&
      !$OMP SHARED(Lttc,Nkpt_layer,Nkpt_path,Norb_layer,data_out,data_in)
      !$OMP DO
      do ik=1,Nkpt_layer
         do ikz=1,Nkpt_path
            do isite=1,Lttc%Nsite
               do jsite=1,Lttc%Nsite
                  !
                  ra = 1+(isite-1)*Norb_layer ; rb = ra + Norb_layer-1
                  ca = 1+(jsite-1)*Norb_layer ; cb = ca + Norb_layer-1
                  !
                  kR = 2*pi * Lttc%kptpath(3,Lttc%Nkpt_path+ikz) * (isite-jsite)
                  cfac = dcmplx(cos(kR),+sin(kR))
                  !
                  data_out(:,:,ik,ikz) = data_out(:,:,ik,ikz) + data_in(ra:rb,ca:cb,ik)*cfac / Lttc%Nsite
                  !
               enddo
            enddo
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
   end subroutine fill_Gamma_A
   !
   !
   !
   subroutine calc_MaxEnt_on_Sigma_K(Gmats_in,mode)
      !
      use linalg, only : diagonal, rotate
      use input_vars, only : ReplaceTail_Simp
      use input_vars, only : Ntau
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats_in
      character(len=*),intent(in)           :: mode
      !
      real(8),allocatable                   :: wmats(:),Sparams(:,:,:,:)
      complex(8),allocatable                :: Gmats_rho(:,:,:,:)
      complex(8),allocatable                :: Smats_diag(:,:,:,:),Sitau_diag(:,:,:,:)
      complex(8),allocatable                :: RotN(:,:,:,:)
      real(8),allocatable                   :: EigN(:,:,:)
      real(8),allocatable                   :: tau(:)
      character(len=256)                    :: ParaFile
      real(8)                               :: M0,M1,M2,M3,M4
      real(8)                               :: x1,x2,x3,y1,y2,y3
      integer                               :: unit,wndx,iw1,iw2,iw3
      integer                               :: Nkpt
      !
      !
      if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- calc_MaxEnt_on_Sigma_K"
      !
      !
      if(.not.Gmats_in%status) stop "calc_MaxEnt_on_Sigma_K: Gmats_in not properly allocated."
      select case(reg(mode))
         case default
            !
            stop "calc_MaxEnt_on_Sigma_K: Available Modes are: path, full, plane."
            !
         case("full")
            !
            Nkpt = Lttc%Nkpt
            write(*,"(A,I)") "     Sigma full. Total number of K-points in the BZ:",Nkpt
            !
         case("path")
            !
            Nkpt = Lttc%Nkpt_path
            write(*,"(A,I)") "     Sigma path. Total number of K-points along path:",Nkpt
            !
         case("plane")
            !
            Nkpt = Lttc%Nkpt_Plane
            write(*,"(A,I)") "     Sigma plane. Total number of K-points in the {kx,ky} sheet:",Nkpt
            !
      end select
      !
      allocate(Sparams(Norb,Nkpt,Nspin,2));Sparams=0d0
      allocate(wmats(Nmats));wmats=FermionicFreqMesh(Sfull%Beta,Nmats)
      wndx = minloc(abs(wmats-ReplaceTail_Simp),dim=1)
      !
      !Find the basis where the K-dependent density matrix is diagoal on the path
      allocate(RotN(Norb,Norb,Nkpt,Nspin));RotN=czero
      allocate(EigN(Norb,Nkpt,Nspin));EigN=0d0
      do ik=1,Nkpt
         do ispin=1,Nspin
            RotN(:,:,ik,ispin) = dreal(Gmats_in%N_ks(:,:,ik,ispin))
            call eigh(RotN(:,:,ik,ispin),EigN(:,ik,ispin))
         enddo
      enddo
      !
      !Operations at each K-point
      allocate(Smats_diag(Norb,Nmats,Nkpt,Nspin));Smats_diag=czero
      do ik=1,Nkpt
         !
         !
         !Bring the Green's function on the path to basis where the K-dependent density matrix is diagonal
         allocate(Gmats_rho(Norb,Norb,Nmats,Nspin));Gmats_rho=czero
         do ispin=1,Nspin
            do iw=1,Nmats
               Gmats_rho(:,:,iw,ispin) = rotate(Gmats_in%wks(:,:,iw,ik,ispin),RotN(:,:,ik,ispin))
            enddo
         enddo
         if(paramagnet) Gmats_rho(:,:,:,Nspin) = Gmats_rho(:,:,:,1)
         !
         !Check Print - This is to check if the rotation has brought the Gf to a diagonal basis
         do ispin=1,Nspin
            path = reg(pathOUTPUT)//"K_resolved/Gk_"//reg(mode)//"_wm_s"//str(ispin)//"/Gk_wm_rho_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nmats
               write(unit,"(200E20.12)") wmats(iw),(Gmats_rho(iorb,iorb,iw,ispin),iorb=1,Norb),(Gmats_rho(1,iorb,iw,ispin),iorb=2,Norb)
            enddo
            close(unit)
            if(paramagnet)exit
         enddo
         !
         !Reference frequencies
         iw1 = wndx-30   ; x1 = wmats(iw1)
         iw2 = wndx-1    ; x2 = wmats(iw2)
         iw3 = Nmats     ; x3 = wmats(iw3) !Never really used
         !
         !Replace the tail of the diagonal Green's function
         do ispin=1,Nspin
            do iorb=1,Norb
               !
               !Real part asymptotic behaviour M0=0, M2, M4
               y1 = dreal(Gmats_rho(iorb,iorb,iw1,ispin))
               y2 = dreal(Gmats_rho(iorb,iorb,iw2,ispin))
               M4 = (y1*x1**2 - y2*x2**2) * ( x1**2 * x2**2 )/( x2**2 - x1**2 )
               M4 = 0d0
               M2 = y2*x2**2 - M4/x2**2
               !
               !Imaginary part asymptotic behaviour M1=1, M3>0 by definition
               y2 = dimag(Gmats_rho(iorb,iorb,iw2,ispin))
               M3 = y2*x2**3 + x2**2
               !
               !Replace the tail of the diagoal Green's function
               do iw=iw2,Nmats
                  Gmats_rho(iorb,iorb,iw,ispin) = dcmplx( M2/(wmats(iw)**2)+M4/(wmats(iw)**4) , -1d0/wmats(iw)+M3/(wmats(iw)**3) )
               enddo
               !
               !The corresponding M0 and M1 of the self-energy
               Sparams(iorb,ik,ispin,1) = Sfull%mu - M2
               Sparams(iorb,ik,ispin,2) = M2**2 - M3
               !
            enddo
         enddo
         !
         !Check Print - This is to check if the tail replacement has any problem
         do ispin=1,Nspin
            path = reg(pathOUTPUT)//"K_resolved/Gk_"//reg(mode)//"_wm_s"//str(ispin)//"/Gk_wm_rho_tail_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nmats
               write(unit,"(200E20.12)") wmats(iw),(Gmats_rho(iorb,iorb,iw,ispin),iorb=1,Norb)
            enddo
            close(unit)
            if(paramagnet)exit
         enddo
         !
         !Get the self-energy on the path in the basis where the K-dependent density matrix is diagonal
         !Here H(k) is enclosed inside the self-energy
         do ispin=1,Nspin
            do iorb=1,Norb
               do iw=1,Nmats
                  Smats_diag(iorb,iw,ik,ispin) = img*wmats(iw) + Sfull%mu - 1d0/Gmats_rho(iorb,iorb,iw,ispin)
               enddo
            enddo
         enddo
         deallocate(Gmats_rho)
         !
         !Check Print - This is to check the self-energy in the diagonal basis
         do ispin=1,Nspin
            path = reg(pathOUTPUT)//"K_resolved/Sk_"//reg(mode)//"_wm_s"//str(ispin)//"/Sk_wm_rho_tail_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nmats
               write(unit,"(200E20.12)") wmats(iw),(Smats_diag(iorb,iw,ik,ispin),iorb=1,Norb)
            enddo
            close(unit)
            if(paramagnet)exit
         enddo
         !
         !Extract the rescaling parameters from the Self-energy just for comparison
         do ispin=1,Nspin
            do iorb=1,Norb
               !
               !Real part asymptotic behaviour M0, M2
               y1 = dreal(Smats_diag(iorb,iw1,ik,ispin))
               y2 = dreal(Smats_diag(iorb,iw2,ik,ispin))
               y3 = dreal(Smats_diag(iorb,iw3,ik,ispin))
               !M0 = get_Meven(y1,y2,y3,x1,x2,x3,0)
               !M2 = get_Meven(y1,y2,y3,x1,x2,x3,2)
               !M4 = get_Meven(y1,y2,y3,x1,x2,x3,4)
               M2 = (y1-y2)*(x1**2 * x2**2)/(x2**2 - x1**2)
               M0 = y2 - M2/x2**2
               !
               !Imaginary part asymptotic behaviour M1
               y2 = dimag(Smats_diag(iorb,iw2,ik,ispin))
               !M1 = get_Modd(y1,y2,x1,x2,1)
               !M3 = get_Modd(y1,y2,x1,x2,3)
               M1 = y2*x2
               !
               !Print comparison
               if(verbose)then
                  write(*,"(5X,2(A8,I5))")"ik=",ik,"  iorb=",iorb
                  write(*,"(5X,2(A12,1F20.8))")"M0(1/G): ",Sparams(iorb,ik,ispin,1),"M0(S): ",M0
                  write(*,"(5X,2(A12,1F20.8))")"M1(1/G): ",Sparams(iorb,ik,ispin,2),"M1(S): ",M1
               endif
               !
               !Remove the bare limit M0
               Smats_diag(iorb,:,ik,ispin) = Smats_diag(iorb,:,ik,ispin) - Sparams(iorb,ik,ispin,1)
               !
               !Rescale so as to have -ImS~1/iw
               Smats_diag(iorb,:,ik,ispin) = Smats_diag(iorb,:,ik,ispin)/abs(Sparams(iorb,ik,ispin,2))
               !
               !Revert the real part
               Smats_diag(iorb,:,ik,ispin) = -conjg(Smats_diag(iorb,:,ik,ispin))
               !
            enddo
         enddo
         !
         !Check Print - This is to check the self-energy in the diagonal basis afer the removal of M0 and rescaling
         do ispin=1,Nspin
            path = reg(pathOUTPUT)//"K_resolved/Sk_"//reg(mode)//"_wm_s"//str(ispin)//"/Sk_wm_rho_tail_rescaled_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nmats
               write(unit,"(200E20.12)") wmats(iw),(Smats_diag(iorb,iw,ik,ispin),iorb=1,Norb)
            enddo
            close(unit)
            if(paramagnet)exit
         enddo
         !
         !
      enddo
      deallocate(wmats)
      !
      !Fourier transform
      !Ntau = 300
      call cpu_time(start)
      allocate(Sitau_diag(Norb,Ntau,Nkpt,Nspin));Sitau_diag=czero
      do ispin=1,Nspin
         call Fmats2itau_vec(Sfull%Beta,Smats_diag(:,:,:,ispin),Sitau_diag(:,:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
         if(paramagnet)then
            Sitau_diag(:,:,:,Nspin) = Sitau_diag(:,:,:,1)
            exit
         endif
      enddo
      deallocate(Smats_diag)
      call cpu_time(finish)
      write(*,"(A,F)") "     Slat(K"//reg(mode)//",iw) --> Slat(K"//reg(mode)//",tau) cpu timing:", finish-start
      !
      !Final correction
      do ik=1,Nkpt
         do iorb=1,Norb
            if(dreal(Sitau_diag(iorb,1,ik,1)).gt.0d0) then
               write(*,"(A)")"     Warning: orbital# "//str(iorb)//" of K-point # "//str(ik)//" is positive in tau=0"
            elseif(dreal(Sitau_diag(iorb,Ntau,ik,1)).gt.0d0)then
               write(*,"(A)")"     Warning: orbital# "//str(iorb)//" of K-point # "//str(ik)//" is positive in tau=beta"
            else
               write(*,"(A)")"     Orbital# "//str(iorb)//" of K-point # "//str(ik)//" is fine. Reverting positive noise."
               do ispin=1,Nspin
                  Sitau_diag(iorb,:,ik,ispin) = -abs(Sitau_diag(iorb,:,ik,ispin))
               enddo
            endif
         enddo
      enddo
      !
      !Print data for K-resolved MaxEnt
      allocate(tau(Ntau));tau = linspace(0d0,Sfull%Beta,Ntau)
      do ispin=1,Nspin
         do ik=1,Nkpt
            !
            path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Sk_"//reg(mode)//"_t_s"//str(ispin)//"/Sk_t_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do itau=1,Ntau
               write(unit,"(200E20.12)") tau(itau),(dreal(Sitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
            enddo
            close(unit)
            !
         enddo
         if(paramagnet)exit
      enddo
      deallocate(Sitau_diag,tau)
      !
      !Print data needed to reconstruct the Self-energy
      ParaFile = reg(pathOUTPUT)//"K_resolved/Sigma_vars/S"//reg(mode)//"_Params.DAT"
      unit = free_unit()
      open(unit,file=reg(ParaFile),form="formatted",status="unknown",position="rewind",action="write")
      do ik=1,Nkpt
         !
         !rotation
         do ispin=1,Nspin
            path = reg(pathOUTPUT)//"K_resolved/Sigma_vars/"
            call dump_Matrix(RotN(:,:,ik,ispin),reg(path),"S"//reg(mode)//"_Rot_k"//str(ik)//"_s"//str(ispin)//".DAT")
            if(paramagnet)exit
         enddo
         !
         !Parameters
         write(unit,"(200E20.12)") (Sparams(iorb,ik,1,1),iorb=1,Norb),(Sparams(iorb,ik,1,2),iorb=1,Norb),         &
                                   (Sparams(iorb,ik,Nspin,1),iorb=1,Norb),(Sparams(iorb,ik,Nspin,2),iorb=1,Norb)
         !
      enddo
      deallocate(Sparams,RotN,EigN)
      close(unit)
      !
   end subroutine calc_MaxEnt_on_Sigma_K
   !
   !
   !
   subroutine calc_MaxEnt_on_Sigma_imp(Smats_in)
      !
      use input_vars, only : ReplaceTail_Simp, PadeWlimit
      use input_vars, only : SiteNorb, SiteOrbs, SiteName, Nsite, EqvGWndx
      use input_vars, only : OlocSite, OlocRot, OlocRotDag, OlocEig
      use input_vars, only : RotateHloc, ExpandImpurity, AFMselfcons, wrealMax
      use input_vars, only : Ntau
      implicit none
      !
      type(FermionicField),intent(in)       :: Smats_in
      !
      type(FermionicField)                  :: Simp
      real(8),allocatable                   :: wmats(:),Sparams(:,:,:)
      complex(8),allocatable                :: Smats_diag(:,:,:)
      complex(8),allocatable                :: Sitau_diag(:,:,:)
      complex(8),allocatable                :: Rot(:,:)
      integer,allocatable                   :: Orbs(:)
      real(8)                               :: M0,M1,M2,M3,M4
      real(8)                               :: x1,x2,x3,y1,y2,y3
      integer                               :: unit,wndx,iw1,iw2,iw3
      integer                               :: isite
      character(len=256)                    :: ParaFile
      !
      !
      if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- calc_MaxEnt_on_Sigma_imp"
      !
      !
      if(.not.Smats_in%status) stop "calc_MaxEnt_on_Sigma_imp: Smats_in not properly allocated."
      !
      allocate(wmats(Nmats));wmats=FermionicFreqMesh(Smats_in%Beta,Nmats)
      if(ReplaceTail_Simp.gt.0d0)then
         wndx = minloc(abs(wmats-ReplaceTail_Simp),dim=1)
      else
         wndx = minloc(abs(wmats-wrealMax),dim=1)
      endif
      !
      !Operations at each site
      do isite=1,Nsite
         !
         allocate(Orbs(SiteNorb(isite)))
         Orbs = SiteOrbs(isite,1:SiteNorb(isite))
         !
         allocate(Smats_diag(SiteNorb(isite),Nmats,Nspin));Smats_diag=czero
         allocate(Sparams(SiteNorb(isite),Nspin,2));Sparams=0d0
         !
         !Get the irreducible local self-energy
         call AllocateFermionicField(Simp,SiteNorb(isite),Nmats,Beta=Smats_in%Beta)
         if(RotateHloc)then
            !
            allocate(Rot(SiteNorb(isite),SiteNorb(isite)))
            Rot=OlocRot(1:SiteNorb(isite),1:SiteNorb(isite),isite)
            call loc2imp(Simp,Smats_in,Orbs,U=Rot)
            !
         else
            !
            call loc2imp(Simp,Smats_in,Orbs)
            !
         endif
         !
         !I need a standard array for the FT
         do ispin=1,Nspin
            do iorb=1,SiteNorb(isite)
               Smats_diag(iorb,:,ispin) = Simp%ws(iorb,iorb,:,ispin)
            enddo
         enddo
         !
         !Check Print - This has to be identical to the one in the Solver_* folder
         call dump_MaxEnt(Smats_diag,"mats",reg(pathOUTPUT)//"Convergence/","Sqmc_"//reg(SiteName(isite)))
         !
         !Reference frequencies
         if(wndx.gt.4)  iw1 = wndx-1
         if(wndx.gt.10) iw1 = wndx-5
         if(wndx.gt.20) iw1 = wndx-10
         if(wndx.gt.40) iw1 = wndx-20
         if(wndx.gt.60) iw1 = wndx-30
                         ; x1 = wmats(iw1)
         iw2 = wndx-1    ; x2 = wmats(iw2)
         iw3 = Nmats     ; x3 = wmats(iw3) !Never really used
         !
         !Extract the rescaling parameters from the Self-energy
         do ispin=1,Nspin
            do iorb=1,SiteNorb(isite)
               !
               !Real part asymptotic behaviour M0, M2
               y1 = dreal(Smats_diag(iorb,iw1,ispin))
               y2 = dreal(Smats_diag(iorb,iw2,ispin))
               y3 = dreal(Smats_diag(iorb,iw3,ispin))
               M0 = get_Meven(y1,y2,y3,x1,x2,x3,0)
               !M2 = get_Meven(y1,y2,y3,x1,x2,x3,2)
               !M4 = get_Meven(y1,y2,y3,x1,x2,x3,4)
               M2 = (y1-y2)*(x1**2 * x2**2)/(x2**2 - x1**2)
               !M0 = y2 - M2/x2**2
               !
               !Imaginary part asymptotic behaviour M1
               y2 = dimag(Smats_diag(iorb,iw2,ispin))
               !M1 = get_Modd(y1,y2,x1,x2,1)
               !M3 = get_Modd(y1,y2,x1,x2,3)
               M1 = y2*x2
               !
               Sparams(iorb,ispin,1) = M0
               Sparams(iorb,ispin,2) = M1
               !
               !------------------------------------------------------------!
               !SPbca - EDMFT correction
               if(iorb.eq.1)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.197
               if(iorb.eq.2)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.1
               if(iorb.eq.3)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.1
               !LPbca - EDMFT correction
               !if(iorb.eq.1)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.213
               !if(iorb.eq.2)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.15
               !if(iorb.eq.3)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.25
               !------------------------------------------------------------!
               !
               !Print comparison
               !if(verbose)then
                  write(*,"(5X,A8,I5,A)")"isite=",isite,"  Element: "//reg(SiteName(isite))
                  write(*,"(5X,2(A8,I5))")"iorb=",iorb,"ispin=",ispin
                  write(*,"(5X,A12,1F20.8)")"M0(S): ",Sparams(iorb,ispin,1)
                  write(*,"(5X,A12,1F20.8)")"M1(S): ",Sparams(iorb,ispin,2)
               !endif
               !
               !Remove the bare limit M0
               Smats_diag(iorb,:,ispin) = Smats_diag(iorb,:,ispin) - Sparams(iorb,ispin,1)
               !
               !Rescale so as to have -ImS~1/iw
               Smats_diag(iorb,:,ispin) = Smats_diag(iorb,:,ispin)/abs(Sparams(iorb,ispin,2))
               !
               !Revert the real part
               Smats_diag(iorb,:,ispin) = -conjg(Smats_diag(iorb,:,ispin))
               !
            enddo
         enddo
         !
         !Check Print - This is to check the self-energy in the diagonal basis after the removal of M0 and rescaling of M1
         call dump_MaxEnt(Smats_diag,"mats",reg(pathOUTPUT)//"Convergence/","Sqmc_rescaled_"//reg(SiteName(isite)))
         !
         !Fourier transform
         !Ntau = 300
         call cpu_time(start)
         allocate(Sitau_diag(SiteNorb(isite),Ntau,Nspin));Sitau_diag=czero
         do ispin=1,Nspin
            call Fmats2itau_vec(Sfull%Beta,Smats_diag(:,:,ispin),Sitau_diag(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
            if(paramagnet)then
               Sitau_diag(:,:,Nspin) = Sitau_diag(:,:,1)
               exit
            endif
         enddo
         deallocate(Smats_diag)
         call cpu_time(finish)
         write(*,"(A,F)") "     Sqmc_"//reg(SiteName(isite))//"(iw) --> Sqmc_"//reg(SiteName(isite))//"(tau) cpu timing:", finish-start
         !
         !Final correction
         do iorb=1,SiteNorb(isite)
            !
            if(dreal(Sitau_diag(iorb,1,1)).gt.0d0) then
               write(*,"(A)")"     Warning: orbital# "//str(iorb)//" is positive in tau=0"
            else
               write(*,"(A)")"     Orbital# "//str(iorb)//" is fine."
               write(*,"(A)")"     Reverting positive noise."
               do ispin=1,Nspin
                  Sitau_diag(iorb,:,ispin) = -abs(Sitau_diag(iorb,:,ispin))
                  Sitau_diag(iorb,Ntau,ispin) = -1d0-Sitau_diag(iorb,1,ispin)
               enddo
            endif
            !
         enddo
         !
         !Print data for MaxEnt
         call dump_MaxEnt(Sitau_diag,"itau",reg(pathOUTPUT)//"Convergence/","Sqmc_"//reg(SiteName(isite)))
         deallocate(Sitau_diag)
         !
         !Print data needed to reconstruct the Self-energy
         ParaFile = reg(pathOUTPUT)//"K_resolved/Sigma_vars/Sqmc_"//reg(SiteName(isite))//"_Params.DAT"
         unit = free_unit()
         open(unit,file=reg(ParaFile),form="formatted",status="unknown",position="rewind",action="write")
         write(unit,"(200E20.12)") (Sparams(iorb,1,1),iorb=1,SiteNorb(isite)),(Sparams(iorb,1,2),iorb=1,SiteNorb(isite)),         & !spin=1
                                   (Sparams(iorb,Nspin,1),iorb=1,SiteNorb(isite)),(Sparams(iorb,Nspin,2),iorb=1,SiteNorb(isite))    !spin=2
         close(unit)
         !
         deallocate(Orbs,Sparams)
         call DeallocateFermionicField(Simp)
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      deallocate(wmats)
      !
   end subroutine calc_MaxEnt_on_Sigma_imp
   !
   !
   !
   function get_Meven(b1,b2,b3,w1,w2,w3,coeff) result(M)
      use linalg, only : det
      implicit none
      real(8),intent(in)      :: b1,b2,b3,w1,w2,w3
      integer,intent(in)      :: coeff
      real(8)                 :: M
      !
      real(8)                 :: num(3,3),den(3,3)
      !
      den(:,1)=1d0
      den(:,2)=[1d0/w1**2,1d0/w2**2,1d0/w3**2]
      den(:,3)=[1d0/w1**4,1d0/w2**4,1d0/w3**4]
      if(det(den).eq.0d0)stop"get_Meven"
      if(coeff.gt.4)stop"wrong even coeff"
      !
      num=den
      num(:,1+coeff/2)=[b1,b2,b3]
      !
      M = det(num) / det(den)
      !
   end function get_Meven
   !
   function get_Modd(b1,b2,w1,w2,coeff) result(M)
      use linalg, only : det
      implicit none
      real(8),intent(in)      :: b1,b2,w1,w2
      integer,intent(in)      :: coeff
      real(8)                 :: M
      !
      real(8)                 :: num(2,2),den(2,2)
      !
      den(:,1)=[1d0/w1   ,1d0/w2   ]
      den(:,2)=[1d0/w1**3,1d0/w2**3]
      if(det(den).eq.0d0)stop"get_Modd"
      if(coeff.gt.3)stop"wrong odd coeff"
      !
      num=den
      if(coeff.eq.1)num(:,1)=[b1,b2]
      if(coeff.eq.3)num(:,2)=[b1,b2]
      !num(:,1+(coeff-1)/2)=[b1,b2]
      !
      M = det(num) / det(den)
      !
   end function get_Modd
   !
   !
   !
end subroutine interpolateG2Path
