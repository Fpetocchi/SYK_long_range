subroutine interpolate2kpath_Fermionic(Sfull,Lttc,pathOUTPUT)
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
   use input_vars, only : wmatsMax, dampStail
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
   !integer                               :: tzl,tzr,ilayer
   real(8)                               :: kx,ky,Bvec(3),Blat(3,3)
   character(len=256)                    :: path
   logical                               :: Kdependence
   real                                  :: start,finish
   !
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- interpolate2kpath_Fermionic"
   !
   !
   if(.not.Lttc%status) stop "interpolate2kpath_Fermionic: Lttc not properly initialized."
   if(.not.Sfull%status) stop "interpolate2kpath_Fermionic: Sfull not properly initialized."
   !
   if(Sfull%Norb.ne.Lttc%Norb) stop "interpolate2kpath_Fermionic: Lttc has different number of orbitals with respect to Sfull."
   if(Sfull%Nkpt.ne.Lttc%Nkpt) stop "interpolate2kpath_Fermionic: Lttc has different number of k-points with respect to Sfull."
   !
   write(*,"(A,F)")"     Chemical potential: ",Sfull%mu
   Norb = Sfull%Norb
   Nmats = Sfull%Npoints
   !
   Kdependence=.false.
   if(scan(reg(CalculationType),"W").gt.0)Kdependence=.true.
   !
   if(FermiSurf)call get_Blat(Blat)
   !
   !
   !
   !---------------------- LDA Hamiltonian and corrections --------------------!
   !
   !non-interacting data (Bands, spectral function, Fermi-surface)
   call interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT=reg(pathOUTPUT)//"K_resolved/",doplane=FermiSurf,hetero=Hetero)
   !
   !correction to LDA given by the real part of the local self-energy in iw=0
   allocate(correction(Norb,Norb,Sfull%Nkpt));correction=czero
   do ispin=1,Nspin
      !
      do ik=1,Sfull%Nkpt
         correction(:,:,ik) = dreal(Sfull%ws(:,:,1,ispin)) - zeye(Norb)*Sfull%mu
      enddo
      call interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT=reg(pathOUTPUT)//"K_resolved/",corrname="dmft_s"//str(ispin),correction=correction,doplane=FermiSurf,hetero=Hetero)
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
            correction(:,:,ik) = Sfull%wks(:,:,1,ik,ispin) - zeye(Norb)*Sfull%mu
         enddo
         do iorb=1,Norb
            correction(iorb,iorb,:) = dcmplx(dreal(correction(iorb,iorb,:)),0d0)
         enddo
         call interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT=reg(pathOUTPUT)//"K_resolved/",corrname="qpsc_s"//str(ispin),correction=correction,doplane=FermiSurf,hetero=Hetero)
         if(paramagnet) exit
         !
      enddo
      deallocate(correction)
   endif
   !
   !
   !
   !------------------- Interpolation of interacting solutuon -----------------!
   !
   !
   !recalculate the internal K-meshes just for Nkpt_Fermi different from default
   if(FermiSurf.and.(Nkpt_Fermi.ne.201))call interpolateHk2Path(Lttc,structure,Nkpt_path,doplane=FermiSurf,Nkpt_Kside=Nkpt_Fermi,hetero=Hetero)
   !
   !Dump MaxEnt data for the local self-energy (done for every CalculationType)
   if(scan(reg(path_funct),"S").gt.0)call calc_MaxEnt_on_Sigma_imp(Sfull)
   !
   !
   !Interpolate the slef-energy along the path if its K-dependent otherwise duplicate the local one
   select case(reg(CalculationType))
      case default
         !
         stop "interpolate2kpath_Fermionic: Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
         !
      case("G0W0","scGW","GW+EDMFT")
         !
         !
         !
         !--------------- Green's function in the full BZ ---------------!
         call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
         call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
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
         !Useless and potentially harmful - dampStail is a testflag inactive by dafault
         if((dampStail.gt.0d0).and.(dampStail.lt.wmatsMax)) call smooth_ImSigma_tail(Spath)
         call dump_FermionicField(Spath,reg(pathOUTPUT)//"K_resolved/Sk_path/","Sk_w",.false.,Lttc%kptpath(:,1:Lttc%Nkpt_path),paramagnet)
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
            !Useless and potentially harmful - dampStail is a testflag inactive by dafault
            if((dampStail.gt.0d0).and.(dampStail.lt.wmatsMax)) call smooth_ImSigma_tail(Sfermi)
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
   !
   !
   if(scan(reg(path_funct),"G").gt.0)then
      !
     !Dump K-resolved MaxEnt data in the full BZ
     call dump_MaxEnt_on_G_K(Gfull,"full")
     !
     !Dump K-resolved MaxEnt data along the path
     call dump_MaxEnt_on_G_K(Gpath,"path")
     !
     !Dump K-resolved MaxEnt data along the {kx,ky} plane
     if(FermiSurf)call dump_MaxEnt_on_G_K(Gfermi,"plane")
     !
   endif
   !
   !
   if(scan(reg(path_funct),"S").gt.0)then
     !
     if(Kdependence)then
         !
         !Dump K-resolved MaxEnt data in the full BZ
         call dump_MaxEnt_on_Sigma_K(Gfull,"full")
         !
         !Dump K-resolved MaxEnt data along the path
         call dump_MaxEnt_on_Sigma_K(Gpath,"path")
         !
         !Dump K-resolved MaxEnt data along the {kx,ky} plane
         call dump_MaxEnt_on_Sigma_K(Gfermi,"plane")
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
   subroutine dump_MaxEnt_on_G_K(Gmats_in,mode)
      !
      use input_vars, only : Ntau_MaxEnt, Nmats_MaxEnt, Solver
      use utils_misc
      implicit none
      !
      type(FermionicField),intent(in)       :: Gmats_in
      character(len=*),intent(in)           :: mode
      !
      complex(8),allocatable                :: Gmats_diag(:,:,:,:),Gmats_diag_tmp(:,:,:,:)
      complex(8),allocatable                :: Gitau_diag(:,:,:,:)
      real(8),allocatable                   :: Ak(:,:),tau(:),wmats(:),Moments(:,:)
      real(8)                               :: Gmax,ReGtail,ImGtail
      integer                               :: Nkpt,NtauFT
      integer                               :: ikx,iky,iw
      !Hetero
      integer                               :: Norb_layer,ikz
      complex(8),allocatable                :: Gmats_kz(:,:,:,:,:,:),Gmats_kz_diag(:,:,:,:,:)
      complex(8),allocatable                :: Gitau_kpkz_diag(:,:,:,:,:)
      !
      !
      if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- dump_MaxEnt_on_G_K"
      !
      !
      NtauFT = Ntau_MaxEnt
      !
      !
      if(.not.Gmats_in%status) stop "dump_MaxEnt_on_G_K: Gmats_in not properly allocated."
      select case(reg(mode))
         case default
            !
            stop "dump_MaxEnt_on_G_K: Available Modes are: path, full, plane."
            !
         case("full")
            !
            Nkpt = Lttc%Nkpt
            write(*,"(A,I)") "     G full. Total number of K-points in the BZ:",Nkpt
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
      do iorb=1,Norb
         Gmats_diag(iorb,:,:,:) = Gmats_in%wks(iorb,iorb,:,:,:)
      enddo
      !
      !Attach the tail before the FT to tau axis
      allocate(wmats(Nmats));wmats=FermionicFreqMesh(Sfull%Beta,Nmats)
      if(Nmats_MaxEnt.gt.Nmats)then
         !
         deallocate(wmats)
         allocate(wmats(Nmats_MaxEnt));wmats=FermionicFreqMesh(Sfull%Beta,Nmats_MaxEnt)
         !
         write(*,"(A)") "     Attaching tail from wm: "//str(wmats(Nmats+1),5)//" to wm: "//str(wmats(Nmats_MaxEnt),5)
         !
         allocate(Gmats_diag_tmp(Norb,Nmats,Nkpt,Nspin));Gmats_diag_tmp=czero
         Gmats_diag_tmp = Gmats_diag
         deallocate(Gmats_diag)
         allocate(Gmats_diag(Norb,Nmats_MaxEnt,Nkpt,Nspin));Gmats_diag=czero
         Gmats_diag(:,1:Nmats,:,:) = Gmats_diag_tmp
         !
         do ispin=1,Nspin
            do ik=1,Nkpt
               !
               call get_moments_F(Moments,Gmats_diag_tmp(:,:,ik,ispin),Gmats_in%Beta,wstep=10,Eo=.false.)
               !
               do iw=Nmats+1,Nmats_MaxEnt
                  do iorb=1,Norb
                     ReGtail = Moments(iorb,2)/(wmats(iw)**2) + Moments(iorb,4)/(wmats(iw)**4)
                     ImGtail =            -1d0/(wmats(iw)**1) + Moments(iorb,3)/(wmats(iw)**3)
                     Gmats_diag(iorb,iw,ik,ispin) = dcmplx(ReGtail,ImGtail)
                  enddo
               enddo
               !
            enddo
            !
            if(paramagnet)then
               Gmats_diag(:,:,:,Nspin) = Gmats_diag(:,:,:,1)
               exit
            endif
            !
         enddo
         deallocate(Gmats_diag_tmp)
         !
      else
         !
         Nmats_MaxEnt = Nmats
         !
      endif
      !
      !Fourier transform the diagonal of the Green's function
      call cpu_time(start)
      allocate(Gitau_diag(Norb,NtauFT,Nkpt,Nspin));Gitau_diag=czero
      do ispin=1,Nspin
         call Fmats2itau_vec(Sfull%Beta,Gmats_diag(:,:,:,ispin),Gitau_diag(:,:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
         if(paramagnet)then
            Gitau_diag(:,:,:,Nspin) = Gitau_diag(:,:,:,1)
            exit
         endif
      enddo
      call cpu_time(finish)
      write(*,"(A,F)") "     Glat(K"//reg(mode)//",iw) --> Glat(K"//reg(mode)//",tau) cpu timing:", finish-start
      !
      !Print data for K-resolved MaxEnt
      allocate(tau(NtauFT));tau = linspace(0d0,Sfull%Beta,NtauFT)
      do ispin=1,Nspin
         do ik=1,Nkpt
            !
            !end-point correction if -G(0)-G(beta) != -1
            do iorb=1,Norb
               Gmax = - (dreal(Gitau_diag(iorb,1,ik,ispin)) + dreal(Gitau_diag(iorb,NtauFT,ik,ispin)))
               Gitau_diag(iorb,1,ik,ispin) = Gitau_diag(iorb,1,ik,ispin)/abs(Gmax)
               Gitau_diag(iorb,NtauFT,ik,ispin) = Gitau_diag(iorb,NtauFT,ik,ispin)/abs(Gmax)
            enddo
            !
            !print on imaginary time axis
            path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do itau=1,NtauFT
                write(unit,"(200E20.12)") tau(itau),(dreal(Gitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
            enddo
            close(unit)
            !
            !print on imaginary frequency axis
            path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_w_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nmats_MaxEnt
                write(unit,"(200E20.12)") wmats(iw),(Gmats_diag(iorb,iw,ik,ispin),iorb=1,Norb)
            enddo
            close(unit)
            !
         enddo
         if(paramagnet)exit
      enddo
      deallocate(tau,wmats,Gmats_diag)
      !
      !Compute the spectral weight at Fermi along the path. See arxiv:0805.3778 Eq.(5)
      do ispin=1,Nspin
         !
         allocate(Ak(Nkpt,Norb));Ak=0d0
         do ik=1,Nkpt
            do iorb=1,Norb
               Ak(ik,iorb) = -dreal(Gitau_diag(iorb,int(NtauFT/2),ik,ispin))*Gmats_in%Beta
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
      !Add the dispersion along the Gamma-A direction
      if(Hetero%status)then
         !
         Norb_layer = Hetero%Norb
         if(Norb_layer.ne.int(Lttc%Norb/Lttc%Nsite)) stop "dump_MaxEnt_on_G_K: wrong hetero orbital dimension."
         !
         !Fill in Gamma-A direction
         allocate(Gmats_kz(Norb_layer,Norb_layer,Nmats,Nkpt,0:Nkpt_path,Nspin));Gmats_kz=czero
         do ispin=1,Nspin
            do iw=1,Nmats
               call fill_Gamma_A(Gmats_in%wks(:,:,iw,:,ispin),Gmats_kz(:,:,iw,:,:,ispin))
            enddo
            if(paramagnet)then
               Gmats_kz(:,:,:,:,:,Nspin) = Gmats_kz(:,:,:,:,:,1)
               exit
            endif
         enddo
         !
         !Extract the diagonal of the Green's function
         allocate(Gmats_kz_diag(Norb_layer,Nmats,Nkpt,0:Nkpt_path,Nspin));Gmats_kz_diag=czero
         do iorb=1,Norb
            Gmats_kz_diag(iorb,:,:,:,:) = Gmats_kz(iorb,iorb,:,:,:,:)
         enddo
         deallocate(Gmats_kz)
         !
         !Fourier transform the diagonal of the Green's function
         call cpu_time(start)
         allocate(Gitau_kpkz_diag(Norb_layer,NtauFT,Nkpt,0:Nkpt_path,Nspin));Gitau_kpkz_diag=czero
         do ispin=1,Nspin
            do ikz=0,Nkpt_path
               call Fmats2itau_vec(Sfull%Beta,Gmats_kz_diag(:,:,:,ikz,ispin),Gitau_kpkz_diag(:,:,:,ikz,ispin),asympt_corr=.true.,tau_uniform=.true.)
            enddo
            if(paramagnet)then
               Gitau_kpkz_diag(:,:,:,:,Nspin) = Gitau_kpkz_diag(:,:,:,:,1)
               exit
            endif
         enddo
         deallocate(Gmats_kz_diag)
         call cpu_time(finish)
         write(*,"(A,F)") "     Glat(K"//reg(mode)//",kz,iw) --> Glat(K"//reg(mode)//",kz,tau) cpu timing:", finish-start
         !
         !Print data for K-resolved MaxEnt
         allocate(tau(NtauFT));tau = linspace(0d0,Sfull%Beta,NtauFT)
         do ispin=1,Nspin
            do ik=1,Nkpt
               !
               !end-point correction if G(0)+G(beta) != -1
               do iorb=1,Norb
                  Gmax = - (dreal(Gitau_kpkz_diag(iorb,1,ik,0,ispin)) + dreal(Gitau_kpkz_diag(iorb,NtauFT,ik,0,ispin)))
                  Gitau_kpkz_diag(iorb,1,ik,0,ispin) = Gitau_kpkz_diag(iorb,1,ik,0,ispin)/abs(Gmax)
                  Gitau_kpkz_diag(iorb,NtauFT,ik,0,ispin) = Gitau_kpkz_diag(iorb,NtauFT,ik,0,ispin)/abs(Gmax)
               enddo
               !
               !print on imaginary time axis
               path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//"_Hetero.DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do itau=1,NtauFT
                  write(unit,"(200E20.12)") tau(itau),(dreal(Gitau_kpkz_diag(iorb,itau,ik,0,ispin)),iorb=1,Norb_layer)
               enddo
               close(unit)
               !
            enddo
            do ik=Nkpt+1,Nkpt+Nkpt_path
               !
               !end-point correction if G(0)+G(beta) != -1
               do iorb=1,Norb
                  Gmax = - (dreal(Gitau_kpkz_diag(iorb,1,Lttc%iq_gamma,ik-Nkpt,ispin)) + dreal(Gitau_kpkz_diag(iorb,NtauFT,Lttc%iq_gamma,ik-Nkpt,ispin)))
                  Gitau_kpkz_diag(iorb,1,Lttc%iq_gamma,ik-Nkpt,ispin) = Gitau_kpkz_diag(iorb,1,Lttc%iq_gamma,ik-Nkpt,ispin)/abs(Gmax)
                  Gitau_kpkz_diag(iorb,NtauFT,Lttc%iq_gamma,ik-Nkpt,ispin) = Gitau_kpkz_diag(iorb,NtauFT,Lttc%iq_gamma,ik-Nkpt,ispin)/abs(Gmax)
               enddo
               !
               !print on imaginary time axis
               path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//"_Hetero.DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do itau=1,NtauFT
                  write(unit,"(200E20.12)") tau(itau),(dreal(Gitau_kpkz_diag(iorb,itau,Lttc%iq_gamma,ik-Nkpt,ispin)),iorb=1,Norb_layer)
               enddo
               close(unit)
               !
            enddo
           if(paramagnet)exit
         enddo
         deallocate(tau,Gitau_kpkz_diag)
      endif
      !
   end subroutine dump_MaxEnt_on_G_K
   !
   !
   !
   subroutine fill_Gamma_A(data_in,data_out)
      !
      implicit none
      !
      complex(8),intent(in)              :: data_in(:,:,:)
      complex(8),intent(inout)           :: data_out(:,:,:,0:)
      !
      integer                            :: ra,rb,ca,cb,ikz
      integer                            :: isite,jsite
      integer                            :: Nkpt_layer
      real(8)                            :: kR
      complex(8)                         :: cfac
      !
      Nkpt_layer = size(data_in,dim=3)
      if(Nkpt_layer.ne.size(data_out,dim=3)) stop "fill_Gamma_A: planar K-mesh does not coincide between layer-resolved and kz integrated."
      !
      data_out=czero
      !$OMP PARALLEL DEFAULT(PRIVATE),&
      !$OMP SHARED(Lttc,Nkpt_layer,Nkpt_path,Hetero,data_out,data_in)
      !$OMP DO
      do ik=1,Nkpt_layer
         do ikz=0,Nkpt_path
            do isite=1,Lttc%Nsite
               do jsite=1,Lttc%Nsite
                  !
                  ra = 1+(isite-1)*Hetero%Norb ; rb = ra + Hetero%Norb-1
                  ca = 1+(jsite-1)*Hetero%Norb ; cb = ca + Hetero%Norb-1
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
   subroutine dump_MaxEnt_on_Sigma_K(Gmats_in,mode)
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
      if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- dump_MaxEnt_on_Sigma_K"
      !
      !
      if(.not.Gmats_in%status) stop "dump_MaxEnt_on_Sigma_K: Gmats_in not properly allocated."
      select case(reg(mode))
         case default
            !
            stop "dump_MaxEnt_on_Sigma_K: Available Modes are: path, full, plane."
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
            path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Sk_"//reg(mode)//"_s"//str(ispin)//"/Sk_t_k"//str(ik)//".DAT"
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
   end subroutine dump_MaxEnt_on_Sigma_K
   !
   !
   !
   subroutine calc_MaxEnt_on_Sigma_imp(Smats_in)
      !
      use input_vars, only : ReplaceTail_Simp, PadeWlimit
      use input_vars, only : LocalOrbs, EqvGWndx
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
      real(8)                               :: M0,M1,M2,M3,M4
      real(8)                               :: x1,x2,x3,y1,y2,y3
      integer                               :: unit,wndx,iw1,iw2,iw3
      integer                               :: isite,Nsite
      character(len=256)                    :: ParaFile
      !
      !
      if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- calc_MaxEnt_on_Sigma_imp"
      !
      !
      if(.not.Smats_in%status) stop "calc_MaxEnt_on_Sigma_imp: Smats_in not properly allocated."
      if(.not.allocated(LocalOrbs)) stop "calc_MaxEnt_on_Sigma_imp: LocalOrbs not properly initialized."
      !
      Nsite = size(LocalOrbs)
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
         allocate(Smats_diag(LocalOrbs(isite)%Norb,Nmats,Nspin));Smats_diag=czero
         allocate(Sparams(LocalOrbs(isite)%Norb,Nspin,2));Sparams=0d0
         !
         !Get the irreducible local self-energy
         call AllocateFermionicField(Simp,LocalOrbs(isite)%Norb,Nmats,Beta=Smats_in%Beta)
         if(RotateHloc)then
            !
            call loc2imp(Simp,Smats_in,LocalOrbs(isite)%Orbs,U=LocalOrbs(isite)%Rot)
            !
         else
            !
            call loc2imp(Simp,Smats_in,LocalOrbs(isite)%Orbs)
            !
         endif
         !
         !I need a standard array for the FT
         do ispin=1,Nspin
            do iorb=1,LocalOrbs(isite)%Norb
               Smats_diag(iorb,:,ispin) = Simp%ws(iorb,iorb,:,ispin)
            enddo
         enddo
         !
         !Check Print - This has to be identical to the one in the Solver_* folder
         call dump_MaxEnt(Smats_diag,"mats",reg(pathOUTPUT)//"Convergence/","Sqmc_"//reg(LocalOrbs(isite)%Name))
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
            do iorb=1,LocalOrbs(isite)%Norb
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
                  write(*,"(5X,A8,I5,A)")"isite=",isite,"  Element: "//reg(LocalOrbs(isite)%Name)
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
         call dump_MaxEnt(Smats_diag,"mats",reg(pathOUTPUT)//"Convergence/","Sqmc_rescaled_"//reg(LocalOrbs(isite)%Name))
         !
         !Fourier transform
         !Ntau = 300
         call cpu_time(start)
         allocate(Sitau_diag(LocalOrbs(isite)%Norb,Ntau,Nspin));Sitau_diag=czero
         do ispin=1,Nspin
            call Fmats2itau_vec(Sfull%Beta,Smats_diag(:,:,ispin),Sitau_diag(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
            if(paramagnet)then
               Sitau_diag(:,:,Nspin) = Sitau_diag(:,:,1)
               exit
            endif
         enddo
         deallocate(Smats_diag)
         call cpu_time(finish)
         write(*,"(A,F)") "     Sqmc_"//reg(LocalOrbs(isite)%Name)//"(iw) --> Sqmc_"//reg(LocalOrbs(isite)%Name)//"(tau) cpu timing:", finish-start
         !
         !Final correction
         do iorb=1,LocalOrbs(isite)%Norb
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
         call dump_MaxEnt(Sitau_diag,"itau",reg(pathOUTPUT)//"Convergence/","Sqmc_"//reg(LocalOrbs(isite)%Name))
         deallocate(Sitau_diag)
         !
         !Print data needed to reconstruct the Self-energy
         ParaFile = reg(pathOUTPUT)//"K_resolved/Sigma_vars/Sqmc_"//reg(LocalOrbs(isite)%Name)//"_Params.DAT"
         unit = free_unit()
         open(unit,file=reg(ParaFile),form="formatted",status="unknown",position="rewind",action="write")
         write(unit,"(200E20.12)") (Sparams(iorb,1,1),iorb=1,LocalOrbs(isite)%Norb),(Sparams(iorb,1,2),iorb=1,LocalOrbs(isite)%Norb),         & !spin=1
                                   (Sparams(iorb,Nspin,1),iorb=1,LocalOrbs(isite)%Norb),(Sparams(iorb,Nspin,2),iorb=1,LocalOrbs(isite)%Norb)    !spin=2
         close(unit)
         !
         deallocate(Sparams)
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
end subroutine interpolate2kpath_Fermionic

subroutine interpolate2kpath_Bosonic(Wfull,Lttc,pathOUTPUT,name,remove_Gamma,NaNb)
   !
   use parameters
   use utils_misc
   use utils_fields
   use linalg, only : eigh, inv, zeye, det
   use crystal
   use file_io
   use greens_function, only : calc_Gmats
   use fourier_transforms
   use input_vars, only : structure, Nkpt_path
   implicit none
   !
   type(BosonicField),intent(in)         :: Wfull
   type(Lattice),intent(inout)           :: Lttc
   character(len=*),intent(in)           :: pathOUTPUT
   character(len=*),intent(in)           :: name
   logical,intent(in),optional           :: remove_Gamma
   logical,intent(in),optional           :: NaNb
   !
   complex(8),allocatable                :: W_orig(:,:,:,:)
   complex(8),allocatable                :: W_intp(:,:,:,:)
   complex(8),allocatable                :: Wgamma(:,:,:)
   real(8),allocatable                   :: wmats(:)
   !
   integer                               :: Norb,Nmats,unit,Wdim
   integer                               :: iq,iw,iorb,jorb,ib1,ib2
   character(len=256)                    :: path
   logical                               :: remove_Gamma_,NaNb_
   real                                  :: start,finish
   !
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- interpolate2kpath_Bosonic"
   !
   !
   if(.not.Lttc%status) stop "interpolate2kpath_Bosonic: Lttc not properly initialized."
   if(.not.Wfull%status) stop "interpolate2kpath_Bosonic: Wfull not properly initialized."
   !
   Norb = int(sqrt(dble(Wfull%Nbp)))
   if(Norb.ne.Lttc%Norb) stop "interpolate2kpath_Bosonic: Lttc has different number of orbitals with respect to Wfull."
   if(Wfull%Nkpt.ne.Lttc%Nkpt) stop "interpolate2kpath_Bosonic: Lttc has different number of k-points with respect to Wfull."
   !
   Nmats = Wfull%Npoints
   !
   remove_Gamma_ = .true.
   if(present(remove_Gamma)) remove_Gamma_ = remove_Gamma
   !
   NaNb_ = .true.
   if(present(NaNb)) NaNb_ = NaNb
   !
   Wdim = Wfull%Nbp
   if(NaNb_) Wdim = Norb
   !
   !
   !
   !------------------- Interpolation of interacting solutuon -----------------!
   !
   !
   !recalculate the internal K-meshes just for Nkpt_Fermi different from default
   if(.not.Lttc%pathStored) call interpolateHk2Path(Lttc,structure,Nkpt_path,doplane=.false.,store=.false.)
   !
   !
   allocate(W_orig(Wdim,Wdim,Nmats,Lttc%Nkpt));W_orig=czero
   allocate(Wgamma(Wdim,Wdim,Nmats));Wgamma=czero
   if(NaNb_)then
      !
      do iorb=1,Norb
         do jorb=1,Norb
            !
            call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
            !
            W_orig(iorb,jorb,:,:) = Wfull%screened(ib1,ib2,:,:)
            !
            if(remove_Gamma_) then
               Wgamma(iorb,jorb,:) = Wfull%screened(ib1,ib2,:,Wfull%iq_gamma)
               do iq=1,Wfull%Nkpt
                  W_orig(iorb,jorb,:,iq) = W_orig(iorb,jorb,:,iq) - Wgamma(iorb,jorb,:)
               enddo
            endif
            !
         enddo
      enddo
      !
   else
      !
      W_orig = Wfull%screened
      !
      if(remove_Gamma_) then
         Wgamma = Wfull%screened(:,:,:,Wfull%iq_gamma)
         do iq=1,Wfull%Nkpt
            W_orig(:,:,:,iq) = W_orig(:,:,:,iq) - Wgamma
         enddo
      endif
      !
   endif
   !
   !
   !Interpolate along the path
   call cpu_time(start)
   allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_path));W_intp=czero
   call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),W_orig,W_intp)
   deallocate(W_orig)
   call cpu_time(finish)
   write(*,"(A,F)") new_line("A")//new_line("A")//"     "//reg(name)//"(fullBZ,iw) --> "//reg(name)//"(Kpath,iw) cpu timing:", finish-start
   !
   !
   !TEST>>>
   !
   !Print along path - I'm printing one orbital per file because my MaxEnt is shitty for bosons
   path = reg(pathOUTPUT)//"K_resolved/MaxEnt_"//reg(name)//"k_path/"
   call createDir(reg(path),verb=verbose)
   allocate(wmats(Nmats))
   wmats = BosonicFreqMesh(Wfull%Beta,Nmats)
   do iq=1,Lttc%Nkpt_path
      unit = free_unit()
      open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
      do iw=1,Nmats-1
          write(unit,"(200E20.12)") wmats(iw),(dreal(W_intp(iorb,iorb,iw,iq)+Wgamma(iorb,iorb,iw)),iorb=1,Wdim)
      enddo
      close(unit)
   enddo
   !
   !temporary for the time-being I have a shitty maxent
   do iorb=1,Wdim
      path = reg(pathOUTPUT)//"K_resolved/MaxEnt_"//reg(name)//"k_path/"
      call createDir(reg(path),verb=verbose)
      do iq=1,Lttc%Nkpt_path
         unit = free_unit()
         open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_o"//str(iorb)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nmats-1
             write(unit,"(200E20.12)") wmats(iw),dreal(W_intp(iorb,iorb,iw,iq)+Wgamma(iorb,iorb,iw))
         enddo
         close(unit)
      enddo
   enddo
   !>>>TEST
   deallocate(wmats,W_intp,Wgamma)
   !
   !
end subroutine interpolate2kpath_Bosonic

!---------------------------------------------------------------------------!
!PURPOSE: Smooth the imaginary part of a K-dependent self-energy on the
!         Imaginary frequency axis. This is done to ease maxent,
!         absolutely NOT during self-consistency.
!---------------------------------------------------------------------------!
subroutine smooth_ImSigma_tail(Smats)
   !
   use parameters
   use utils_misc
   use input_vars, only : paramagnet, dampStail
   implicit none
   !
   type(FermionicField),intent(inout)    :: Smats
   !
   integer                               :: Nmats,Nkpt,Norb,wndx
   integer                               :: ispin,iw,ik,iorb,jorb
   real(8)                               :: Beta,wncut,Moment,ReS,ImS
   real(8)                               :: limit=0.999
   real(8),allocatable                   :: wmats(:)
   integer                               :: TailPower=1
   !
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- smooth_ImSigma_tail"
   !
   !
   if(.not.Smats%status) stop "smooth_ImSigma_tail: Self-energy not properly initialized."
   !
   Nkpt = Smats%Nkpt
   Nmats = Smats%Npoints
   Norb = Smats%Norb
   Beta = Smats%Beta
   !
   allocate(wmats(Nmats));wmats=FermionicFreqMesh(Beta,Nmats)
   !
   wndx = minloc(abs(wmats-dampStail),dim=1)
   wncut = wmats(wndx) + log(limit/(1d0-limit))/Beta
   write(*,"(A)")"     Center of the Fermi function: "//str(wncut,5)
   !
   !$OMP PARALLEL DEFAULT(SHARED),&
   !$OMP PRIVATE(ik,ispin,iorb,jorb,iw,ReS,ImS,Moment)
   !$OMP DO
   do ik=1,Nkpt
      do ispin=1,Nspin
         do iorb=1,Norb
            do jorb=1,Norb
               !
               Moment = dimag(Smats%wks(iorb,jorb,wndx,ik,ispin))*(wmats(wndx)**TailPower)
               do iw = wndx,Nmats
                 !ReS = dreal(Smats%wks(iorb,jorb,wndx,ik,ispin)) ! constant real part
                  ReS = dreal(Smats%wks(iorb,jorb,iw,ik,ispin))   ! pristine real part
                  ImS = dimag(Smats%wks(iorb,jorb,iw,ik,ispin))
                  ImS = ImS + ( Moment/(wmats(iw)**TailPower) - ImS )*(1d0-fermidirac(wmats(iw),wncut,Beta))
                  Smats%wks(iorb,jorb,iw,ik,ispin) = dcmplx(ReS,ImS)
               enddo
               !
            enddo
         enddo
         if(paramagnet)then
            Smats%wks(:,:,wndx:Nmats,ik,Nspin) = Smats%wks(:,:,wndx:Nmats,ik,1)
            cycle
         endif
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
   deallocate(wmats)
   !
end subroutine smooth_ImSigma_tail





!call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
!call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
!write(*,*)"Gfull(K) computed"
!!
!allocate(Gfull_R(Norb,Norb,Nmats,Nwig));Gfull_R=czero
!call wannier_K2R(Lttc%Nkpt3,Lttc%kpt,Gfull%wks(:,:,:,:,1),Gfull_R)
!call DeallocateField(Gfull)
!write(*,*)"Gfull(R) computed"
!!
!allocate(Gfull_R_red(2,Nmats,Nwig));Gfull_R_red=czero
!do iwig=1,Nwig
!   if(Nvecwig(3,iwig).eq.0) then
!      Gfull_R_red(1,:,iwig) = Gfull_R(1,1,:,iwig)
!      Gfull_R_red(2,:,iwig) = Gfull_R(2,2,:,iwig)
!   endif
!enddo
!write(*,*)"Gfull_red(R) plane set"
!!
!write(*,*)"Max nx: ",maxval(Nvecwig(1,:))
!write(*,*)"Min nx: ",minval(Nvecwig(1,:))
!write(*,*)"Max ny: ",maxval(Nvecwig(2,:))
!write(*,*)"Min ny: ",minval(Nvecwig(2,:))
!write(*,*)"Max nz: ",maxval(Nvecwig(3,:))
!write(*,*)"Min nz: ",minval(Nvecwig(3,:))
!!
!do nz=minval(Nvecwig(3,:)),maxval(Nvecwig(3,:)),1
!   !
!   if(nz.ne.0)then
!      do nx=minval(Nvecwig(1,:)),maxval(Nvecwig(1,:)),1
!         do ny=minval(Nvecwig(2,:)),maxval(Nvecwig(2,:)),1
!            !
!            iwig = find_vec([nx,ny,nz],Nvecwig,hardstop=.false.)
!            if(iwig.eq.0)cycle
!            !
!            ! effective single band of site 1 ----------------------
!            nz_old = ceiling(nz/2d0)
!            iwig_old = find_vec([nx,ny,nz_old],Nvecwig,hardstop=.true.)
!            !
!            comp = 2
!            if(mod(nz,2).eq.0) comp = 1
!            !
!            Gfull_R_red(1,:,iwig) = Gfull_R(1,comp,:,iwig_old)
!            if((nx.eq.0).and.(ny.eq.0))then
!               write(*,"(A,10I)")"upper- ",nz,nz_old,1,comp
!            endif
!            !
!            ! effective single band of site 2 ----------------------
!            nz_old = floor(nz/2d0)
!            iwig_old = find_vec([nx,ny,nz_old],Nvecwig,hardstop=.true.)
!            !
!            comp = 1
!            if(mod(nz,2).eq.0) comp = 2
!            !
!            Gfull_R_red(2,:,iwig) = Gfull_R(2,comp,:,iwig_old)
!            if((nx.eq.0).and.(ny.eq.0))then
!               write(*,"(A,10I)")"lower- ",nz,nz_old,2,comp
!            endif
!            !
!         enddo
!      enddo
!   endif
!   !
!enddo
!deallocate(Gfull_R)
!write(*,*)"Gfull_red(R) out-of-plane set"
!!
!do iorb=1,2
!   call wannier_R2K(Lttc%Nkpt3,Lttc%kptpath(:,1:Lttc%Nkpt_path),Gfull_R_red(iorb,:,:),Gpath%wks(iorb,iorb,:,:,1))
!   write(*,*)"Gfull_red(K) FT orb:",iorb
!enddo
!Gpath%wks(:,:,:,:,2) = Gpath%wks(:,:,:,:,1)
!deallocate(Gfull_R_red)
!
