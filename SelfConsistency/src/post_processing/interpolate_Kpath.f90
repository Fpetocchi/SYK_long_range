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
   use input_vars, only : structure, Nkpt_path, Nkpt_plane
   use input_vars, only : print_path_G, print_full_G , print_path_S, print_plane_G
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
   integer                               :: ik,ispin,iorb
   integer                               :: ikx,iky
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
   if(print_plane_G)call get_Blat(Blat)
   !
   call createDir(reg(pathOUTPUT),verb=verbose)
   !
   !
   !
   !---------------------- LDA Hamiltonian and corrections --------------------!
   !
   !non-interacting data (Bands, spectral function, Fermi-surface)
   call interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT=reg(pathOUTPUT),doplane=.true.,hetero=Hetero,store=.true.)
   !
   !correction to LDA given by the real part of the local self-energy in iw=0
   allocate(correction(Norb,Norb,Sfull%Nkpt));correction=czero
   do ispin=1,Nspin
      !
      do ik=1,Sfull%Nkpt
         correction(:,:,ik) = dreal(Sfull%ws(:,:,1,ispin)) - zeye(Norb)*Sfull%mu
      enddo
      call interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT=reg(pathOUTPUT),corrname="dmft_s"//str(ispin),correction=correction,doplane=.true.,hetero=Hetero)
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
         call interpolateHk2Path(Lttc,structure,Nkpt_path,pathOUTPUT=reg(pathOUTPUT),corrname="qpsc_s"//str(ispin),correction=correction,doplane=.true.,hetero=Hetero)
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
   !recalculate the internal K-meshes just for Nkpt_plane different from default
   if(print_plane_G.and.(Nkpt_plane.ne.201))call interpolateHk2Path(Lttc,structure,Nkpt_path,doplane=print_plane_G,Nkpt_Kside=Nkpt_plane,hetero=Hetero,store=.true.)
   !
   !
   ! print the Hamiltonian on the path used
   if(print_path_G)call dump_Hk(Lttc%Hk_path,Lttc%kptpath,reg(pathOUTPUT),"Hk_path.DAT")
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
         !------------------ Green's function in the full BZ ------------------!
         if(print_full_G)then
            call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
         endif
         !
         !
         !
         !------------------ Green's function along the path ------------------!
         if(print_path_G)then
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
            !print the self-energy on the path - matsubara, not tau
            if(print_path_S) call dump_FermionicField(Spath,reg(pathOUTPUT)//"Sk_path/","Sk_w",.false.,Lttc%kptpath(:,1:Lttc%Nkpt_path),paramagnet)
            !
            !Compute the quasiparticle weight along the path
            do ispin=1,Nspin
               allocate(Zk(Lttc%Nkpt_path,Norb));Zk=0d0
               do ik=1,Lttc%Nkpt_path
                  do iorb=1,Norb
                     Zk(ik,iorb) = 1d0 / (1d0 + abs(dimag(Spath%wks(iorb,iorb,1,ik,ispin)))*Spath%Beta/pi)
                  enddo
               enddo
               path = reg(pathOUTPUT)//"Zk_path_s"//str(ispin)//".DAT"
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
         endif
         !
         !
         !
         !-------------- Green's function along the {kx,ky} plane -------------!
         if(print_plane_G)then
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
               path = reg(pathOUTPUT)//"Zk_plane_s"//str(ispin)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do ik=1,Lttc%Nkpt_Plane
                  ikx = int(ik/(Nkpt_plane+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_plane-1) - 0.5d0
                  iky = ik - (ikx-1)*Nkpt_plane      ; ky = (iky-1)/dble(Nkpt_plane-1) - 0.5d0
                  Bvec = kx*Blat(:,1) + ky*Blat(:,2)
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Zk(ik,iorb),iorb=1,Norb)
                  if(iky.eq.Nkpt_plane)write(unit,*)
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
         ! store a local self-energy
         call AllocateFermionicField(Sloc,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
         do ik=1,Lttc%Nkpt_path
            Sloc%wks(:,:,:,ik,:) = Sfull%ws
         enddo
         !
         !
         !
         !------------------ Green's function in the full BZ ------------------!
         if(print_full_G)then
            call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gfull,Lttc,Smats=Sloc,along_path=.false.)
         endif
         !
         !
         !
         !------------------ Green's function along the path ------------------!
         if(print_path_G)then
            call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gpath,Lttc,Smats=Sloc,along_path=.true.)
         endif
         !
         !
         !
         !-------------- Green's function along the {kx,ky} plane -------------!
         if(print_plane_G)then
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
   !Dump K-resolved MaxEnt data in the full BZ
   if(print_full_G)call dump_MaxEnt_on_G_K(Gfull,"full")
   !Dump K-resolved MaxEnt data along the path
   if(print_path_G)call dump_MaxEnt_on_G_K(Gpath,"path")
   !Dump K-resolved MaxEnt data along the {kx,ky} plane
   if(print_plane_G)call dump_MaxEnt_on_G_K(Gfermi,"plane")
   !
   !
   !
   call DeallocateFermionicField(Gpath)
   call DeallocateFermionicField(Gfull)
   call DeallocateFermionicField(Gfermi)
   !
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
      complex(8),allocatable                :: Gmats_diag(:,:,:,:),Gmats_diag_tmp(:,:,:,:),Gmats_trace(:,:,:)
      complex(8),allocatable                :: Gitau_diag(:,:,:,:),Gitau_trace(:,:,:)
      real(8),allocatable                   :: Ak(:,:),tau(:),wmats(:),Moments(:,:)
      real(8)                               :: Gmax,ReGtail,ImGtail
      integer                               :: Nkpt,NtauFT
      integer                               :: ikx,iky,iw,itau
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
      do ispin=1,Nspin
         call createDir(reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin),verb=verbose)
         if(paramagnet)exit
      enddo
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
      !Compute the trace
      if(Norb.gt.1)then
         allocate(Gmats_trace(Nkpt,Nmats_MaxEnt,Nspin));Gmats_trace=czero
         do ik=1,Nkpt
            do iorb=1,Norb
               Gmats_trace(ik,:,:) = Gmats_trace(ik,:,:) + Gmats_diag(iorb,:,ik,:)/Norb
            enddo
         enddo
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
      if(Norb.gt.1)then
         call cpu_time(start)
         allocate(Gitau_trace(Nkpt,NtauFT,Nspin));Gitau_trace=czero
         do ispin=1,Nspin
            call Fmats2itau_vec(Sfull%Beta,Gmats_trace(:,:,ispin),Gitau_trace(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
            if(paramagnet)then
               Gitau_trace(:,:,Nspin) = Gitau_trace(:,:,1)
               exit
            endif
         enddo
         call cpu_time(finish)
         write(*,"(A,F)") "     TrGlat(K"//reg(mode)//",iw) --> TrGlat(K"//reg(mode)//",tau) cpu timing:", finish-start
      endif
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
            path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do itau=1,NtauFT
                write(unit,"(200E20.12)") tau(itau),(dreal(Gitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
            enddo
            close(unit)
            !
            !print on imaginary frequency axis
            path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_w_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nmats_MaxEnt
                write(unit,"(200E20.12)") wmats(iw),(Gmats_diag(iorb,iw,ik,ispin),iorb=1,Norb)
            enddo
            close(unit)
            !
            if(Norb.gt.1)then
               !
               Gmax = - (dreal(Gitau_trace(ik,1,ispin)) + dreal(Gitau_trace(ik,NtauFT,ispin)))
               Gitau_trace(ik,1,ispin) = Gitau_trace(ik,1,ispin)/abs(Gmax)
               Gitau_trace(ik,NtauFT,ispin) = Gitau_trace(ik,NtauFT,ispin)/abs(Gmax)
               !
               path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//"_Tr.DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do itau=1,NtauFT
                   write(unit,"(200E20.12)") tau(itau),dreal(Gitau_trace(ik,itau,ispin))
               enddo
               close(unit)
               !
               path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_w_k"//str(ik)//"_Tr.DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats_MaxEnt
                   write(unit,"(200E20.12)") wmats(iw),Gmats_trace(ik,iw,ispin)
               enddo
               close(unit)
               !
            endif
            !
         enddo
         if(paramagnet)exit
      enddo
      if(Norb.gt.1)deallocate(Gmats_trace,Gitau_trace)
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
         path = reg(pathOUTPUT)//"Ak_"//reg(mode)//"_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         if(reg(mode).eq."path")then
            do ik=1,Nkpt
               write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik)/Lttc%Kpathaxis(Lttc%Nkpt_path),(Ak(ik,iorb),iorb=1,Norb)
            enddo
         elseif(reg(mode).eq."plane")then
            do ik=1,Nkpt
               ikx = int(ik/(Nkpt_plane+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_plane-1) - 0.5d0
               iky = ik - (ikx-1)*Nkpt_plane      ; ky = (iky-1)/dble(Nkpt_plane-1) - 0.5d0
               Bvec = kx*Blat(:,1) + ky*Blat(:,2)
               write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Ak(ik,iorb),iorb=1,Norb)
               if(iky.eq.Nkpt_plane)write(unit,*)
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
               path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//"_Hetero.DAT"
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
               path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//"_Hetero.DAT"
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
   use input_vars, only : structure, Nkpt_path, Nkpt_plane
   use input_vars, only : print_path_Chi, print_path_W, print_plane_W
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
   !I forgot why is better to remove the Gamma point
   remove_Gamma_ = .true.
   if(present(remove_Gamma)) remove_Gamma_ = remove_Gamma
   !
   !this is to restrict the printing to the density-density terms only
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
   !recalculate the internal K-meshes just for Nkpt_plane different from default
   if(print_plane_W.and.(Nkpt_plane.ne.201))call interpolateHk2Path(Lttc,structure,Nkpt_path,doplane=print_plane_W,Nkpt_Kside=Nkpt_plane,store=.true.)
   !
   !
   !Optionally keep only the density-density terms
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
   !
   !-------------------------- Boson along the path ---------------------------!
   if(print_path_W.or.print_path_Chi)then
      !
      !Interpolate the boson along the path
      call cpu_time(start)
      allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_path));W_intp=czero
      call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),W_orig,W_intp)
      call cpu_time(finish)
      write(*,"(A,F)") new_line("A")//new_line("A")//"     "//reg(name)//"(fullBZ,iw) --> "//reg(name)//"(Kpath,iw) cpu timing:", finish-start
      !
      !Print along path
      path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_path/"
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
      !TEST>>>
      !temporary for the time-being I have a shitty maxent
      do iorb=1,Wdim
         path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_path/"
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
      deallocate(wmats,W_intp)
      !
   endif
   !
   !
   !
   !---------------------- Boson along the {kx,ky} plane ----------------------!
   if(print_plane_W)then
      !
      !Interpolate the boson along the plane
      call cpu_time(start)
      allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_Plane));W_intp=czero
      call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,W_orig,W_intp)
      call cpu_time(finish)
      write(*,"(A,F)") new_line("A")//new_line("A")//"     "//reg(name)//"(fullBZ,iw) --> "//reg(name)//"(kx,ky,iw) cpu timing:", finish-start
      !
      !Print along plane
      path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_plane/"
      call createDir(reg(path),verb=verbose)
      allocate(wmats(Nmats))
      wmats = BosonicFreqMesh(Wfull%Beta,Nmats)
      do iq=1,Lttc%Nkpt_Plane
         unit = free_unit()
         open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nmats-1
             write(unit,"(200E20.12)") wmats(iw),(dreal(W_intp(iorb,iorb,iw,iq)+Wgamma(iorb,iorb,iw)),iorb=1,Wdim)
         enddo
         close(unit)
      enddo
      !
      !TEST>>>
      !temporary for the time-being I have a shitty maxent
      do iorb=1,Wdim
         path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_plane/"
         call createDir(reg(path),verb=verbose)
         do iq=1,Lttc%Nkpt_Plane
            unit = free_unit()
            open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_o"//str(iorb)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nmats-1
                write(unit,"(200E20.12)") wmats(iw),dreal(W_intp(iorb,iorb,iw,iq)+Wgamma(iorb,iorb,iw))
            enddo
            close(unit)
         enddo
      enddo
      !>>>TEST
      deallocate(wmats,W_intp)
      !
   endif
   !
   deallocate(W_orig,Wgamma)
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
