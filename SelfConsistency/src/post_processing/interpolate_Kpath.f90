subroutine interpolate2kpath_Fermionic(Sfull,Lttc,pathOUTPUT)
   !
   use parameters
   use utils_misc
   use utils_fields
   use linalg, only : eigh, inv, zeye, det, rotate
   use crystal
   use file_io
   use greens_function, only : calc_Gmats
   use fourier_transforms
   use input_vars, only : Nkpt_path_default, Nkpt_plane_default
   use input_vars, only : Nkpt_path, Nkpt_plane
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
   integer                               :: ik,ispin,iorb,iw
   integer                               :: ikx,iky
   real(8)                               :: k1,k2,k3,Bx,By,Bz,Bx_old,Blat(3,3)
   character(len=256)                    :: path
   logical                               :: Kdependence
   logical                               :: WannInterpSigma=.true.
   real                                  :: start,finish
   !
   !TEST>>>
   logical                               :: integrate_kz=.false.
   logical                               :: dump_BAB=.false.
   logical                               :: dump_unfolded=.false.,print_AB=.false.
   integer                               :: Nkp,Nkz,Nkz_,ikz,iwig
   integer                               :: nx,ny,nz,isite,jsite,iwig_old
   real(8)                               :: kz,kzfact,Bz_thresh,KR,cfac
   real(8),allocatable                   :: kptpath_auxkz(:,:)
   complex(8),allocatable                :: Gfull_R(:,:,:,:),Gfull_R_halfkz(:,:,:,:)
   type(FermionicField)                  :: Gpath_auxkz,Gfull_aux
   !>>>TEST
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
   !non-interacting data (Bands, spectral function, Fermi-surface). This has the wrong chemical potential stored in Lttc. Kee the data stored in GWinput
   call interpolate2Path(Lttc,Nkpt_path_default,"Hk",pathOUTPUT=reg(pathOUTPUT),store=.false.,skipAkw=.false.)
   call interpolate2Plane(Lttc,Nkpt_plane_default,"Hk",pathOUTPUT=reg(pathOUTPUT),store=.false.,skipFk=.false.)
   !
   !correction to LDA given by the real part of the local self-energy in iw=0.
   allocate(correction(Norb,Norb,Sfull%Nkpt));correction=czero
   do ispin=1,Nspin
      !
      do ik=1,Sfull%Nkpt
         correction(:,:,ik) = Lttc%Hk(:,:,ik) + dreal(Sfull%ws(:,:,1,ispin)) - zeye(Norb)*Sfull%mu
      enddo
      !
      call interpolate2Path(Lttc,Nkpt_path_default,"Hk_dmft_s"//str(ispin),pathOUTPUT=reg(pathOUTPUT),data_in=correction,store=.false.,skipAkw=.true.)   !Spectral function not computed as the correction is valid only at w=0
      call interpolate2Plane(Lttc,Nkpt_plane_default,"Hk_dmft_s"//str(ispin),pathOUTPUT=reg(pathOUTPUT),data_in=correction,store=.false.,skipFk=.false.) !Fermi surface is computed as the correction is valid only at w=0
      !
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
            correction(:,:,ik) = Lttc%Hk(:,:,ik) + Sfull%wks(:,:,1,ik,ispin) - zeye(Norb)*Sfull%mu
         enddo
         do iorb=1,Norb
            correction(iorb,iorb,:) = dcmplx(dreal(correction(iorb,iorb,:)),0d0)
         enddo
         !
         call interpolate2Path(Lttc,Nkpt_path_default,"Hk_qpsc_s"//str(ispin),pathOUTPUT=reg(pathOUTPUT),data_in=correction,store=.false.,skipAkw=.true.)   !Spectral function not computed as the correction is valid only at w=0
         call interpolate2Plane(Lttc,Nkpt_plane_default,"Hk_qpsc_s"//str(ispin),pathOUTPUT=reg(pathOUTPUT),data_in=correction,store=.false.,skipFk=.false.) !Fermi surface is computed as the correction is valid only at w=0
         !
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
   !recalculate the internal K-meshes
   if(print_path_G)then
      call interpolate2Path(Lttc,Nkpt_path,"Hk",store=.true.)
      call dump_Hk(Lttc%Hk_path,Lttc%kptpath,reg(pathOUTPUT),"Hk_path.DAT")
   endif
   !
   if(print_plane_G)call interpolate2Plane(Lttc,Nkpt_plane,"Hk",store=.true.)
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
            !
            !TEST>>>
            do ispin=1,Nspin
               !
               do ik=1,Lttc%Nkpt
                  do iw=1,Nmats
                     Gfull%wks(:,:,iw,ik,ispin) = rotate(Gfull%wks(:,:,iw,ik,ispin),Lttc%Zk(:,:,ik))
                  enddo
               enddo
               if(paramagnet)then
                  Gfull%wks(:,:,:,:,Nspin) = Gfull%wks(:,:,:,:,1)
                  exit
               endif
            enddo
            !>>>TEST
            !
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
            call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            if(WannInterpSigma)then
               !
               !Recompute the Green's function using interpolated Sigma and Hk
               call calc_Gmats(Gpath,Lttc,Smats=Spath,along_path=.true.)
               call DeallocateFermionicField(Spath)
               !
            else
               !
               !Interpolate directly the Green's function from the BZ to the path
               call DeallocateFermionicField(Spath)
               !
               call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
               !
               call cpu_time(start)
               do ispin=1,Nspin
                  call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),Gfull%wks(:,:,:,:,ispin),Gpath%wks(:,:,:,:,ispin))
                  if(paramagnet)then
                     Gpath%wks(:,:,:,:,Nspin) = Gpath%wks(:,:,:,:,1)
                     exit
                  endif
               enddo
               call DeallocateFermionicField(Gfull)
               call cpu_time(finish)
               write(*,"(A,F)") new_line("A")//new_line("A")//"     G(fullBZ,iw) --> G(Kpath,iw) cpu timing:", finish-start
               !
            endif
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
                  !
                  ikx = int(ik/(Nkpt_plane+0.001))+1
                  iky = ik - (ikx-1)*Nkpt_plane
                  !
                  k1 = Lttc%kptPlane(1,ik)
                  k2 = Lttc%kptPlane(2,ik)
                  k3 = Lttc%kptPlane(3,ik)
                  !
                  Bx = k1*Blat(1,1) + k2*Blat(1,2) + k3*Blat(1,3) ; if(ik.eq.1) Bx_old = Bx
                  By = k1*Blat(2,1) + k2*Blat(2,2) + k3*Blat(2,3)
                  Bz = k1*Blat(3,1) + k2*Blat(3,2) + k3*Blat(3,3)
                  !
                  write(unit,"(3I10,200E20.12)") ik,ikx,iky,k1,k2,k3,Bx,By,Bz,(Zk(ik,iorb),iorb=1,Norb)
                  if(iky.eq.Nkpt_plane)write(unit,*)
                  !
               enddo
               close(unit)
               deallocate(Zk)
               if(paramagnet)exit
            enddo
            !
            call AllocateFermionicField(Gfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            if(WannInterpSigma)then
               !
               !Recompute the Green's function using interpolated Sigma and Hk
               call calc_Gmats(Gfermi,Lttc,Smats=Sfermi,along_plane=.true.)
               call DeallocateFermionicField(Sfermi)
               !
            else
               !
               !Interpolate directly the Green's function from the BZ to the plane
               call DeallocateFermionicField(Sfermi)
               !
               call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
               !
               call cpu_time(start)
               do ispin=1,Nspin
                  call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,Gfull%wks(:,:,:,:,ispin),Gfermi%wks(:,:,:,:,ispin))
                  if(paramagnet)then
                     Gfermi%wks(:,:,:,:,Nspin) = Gfermi%wks(:,:,:,:,1)
                     exit
                  endif
               enddo
               call DeallocateFermionicField(Gfull)
               call cpu_time(finish)
               write(*,"(A,F)") "     G(fullBZ,iw) --> G(kx,ky,iw) cpu timing:", finish-start
               !
            endif
            !
         endif
         !
         !
         !
      case("DMFT+statU","DMFT+dynU","EDMFT")
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
            call AllocateFermionicField(Sloc,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            do ik=1,Lttc%Nkpt_path
               Sloc%wks(:,:,:,ik,:) = Sfull%ws
            enddo
            !
            call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gpath,Lttc,Smats=Sloc,along_path=.true.)
            !
            !TEST>>>
            !!if((.not.integrate_kz).and.(.not.dump_BAB).and.(.not.dump_unfolded)) call calc_Gmats(Gpath,Lttc,Smats=Sloc,along_path=.true.)
            !>>>TEST
            !
            !TEST>>> This portion is to integrate over a given thickness in kz
            !!if(integrate_kz)then
            !!   !
            !!   call get_Blat(Blat)
            !!   allocate(kptpath_auxkz(3,Lttc%Nkpt_path));kptpath_auxkz=0d0
            !!   !
            !!   Bz_thresh = 5.1
            !!   Nkz = 100
            !!   Nkz_ = Nkz
            !!   !
            !!   do ikz=1,Nkz
            !!      !
            !!      kptpath_auxkz = Lttc%kptpath(:,1:Lttc%Nkpt_path)
            !!      !
            !!      kz = (ikz-1)*1d0/Nkz - 0.5d0
            !!      kptpath_auxkz(3,:) = kptpath_auxkz(3,:) + kz
            !!      !
            !!      unit = free_unit()
            !!      open(unit,file="Ktest_"//str(ikz)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
            !!      do ik=1,Lttc%Nkpt_path
            !!         !
            !!         k1 = kptpath_auxkz(1,ik)
            !!         k2 = kptpath_auxkz(2,ik)
            !!         k3 = kptpath_auxkz(3,ik)
            !!         !
            !!         Bx = k1*Blat(1,1) + k2*Blat(1,2) + k3*Blat(1,3)
            !!         By = k1*Blat(2,1) + k2*Blat(2,2) + k3*Blat(2,3)
            !!         Bz = k1*Blat(3,1) + k2*Blat(3,2) + k3*Blat(3,3)
            !!         !
            !!         write(unit,"(1I5,200E20.12)") ik,0d0,kptpath_auxkz(:,ik),Bx,By,Bz
            !!         !
            !!      enddo
            !!      close(unit)
            !!      !
            !!      kzfact = 1d0
            !!      if(abs(Bz).gt.Bz_thresh)then
            !!         kzfact = 0d0
            !!         Nkz_ = Nkz_ - 1
            !!      endif
            !!      write(*,*)ikz,"Bz",Bz,abs(Bz),kzfact,Nkz_
            !!      !
            !!      call cpu_time(start)
            !!      Lttc%Hk_path = czero
            !!      call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,kptpath_auxkz(:,1:Lttc%Nkpt_path),Lttc%Hk,Lttc%Hk_path)
            !!      call cpu_time(finish)
            !!      !write(*,"(A,F)") "     (fullBZ) --> (Kpath,kz="//str(ikz)//") cpu timing:", finish-start
            !!      !
            !!      call AllocateFermionicField(Gpath_auxkz,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            !!      call calc_Gmats(Gpath_auxkz,Lttc,Smats=Sloc,along_path=.true.)
            !!      !
            !!      Gpath%wks = Gpath%wks + kzfact*Gpath_auxkz%wks
            !!      !
            !!      call DeallocateFermionicField(Gpath_auxkz)
            !!      !
            !!   enddo
            !!   Gpath%wks = Gpath%wks/Nkz_
            !!   !
            !!endif
            !>>>TEST
            !
            !TEST>>> This portion is to pass from orbital basis to B/AB basis directly on Gpath
            !!if(dump_BAB)then
            !!   !
            !!   call calc_Gmats(Gpath,Lttc,Smats=Sloc,along_path=.true.)
            !!   !
            !!   call duplicate(Gpath_auxkz,Gpath)
            !!   call DeallocateField(Gpath)
            !!   !
            !!   call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            !!   !
            !!   Gpath%wks(1,1,:,:,:) = (Gpath_auxkz%wks(1,1,:,:,:)+Gpath_auxkz%wks(2,2,:,:,:)+Gpath_auxkz%wks(1,2,:,:,:)+Gpath_auxkz%wks(2,1,:,:,:))/2d0
            !!   Gpath%wks(2,2,:,:,:) = (Gpath_auxkz%wks(1,1,:,:,:)+Gpath_auxkz%wks(2,2,:,:,:)-Gpath_auxkz%wks(1,2,:,:,:)-Gpath_auxkz%wks(2,1,:,:,:))/2d0
            !!   !
            !!   call DeallocateFermionicField(Gpath_auxkz)
            !!   !
            !!endif
            !>>>TEST
            !
            !TEST>>> This portion is to pass from orbital basis to B/AB basis directly on Gpath
            !!if(dump_unfolded)then
            !!   !
            !!   call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            !!   call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
            !!   write(*,*)"Gfull(K) computed"
            !!   !
            !!   if(print_AB)then
            !!      Gfull%wks(1,2,:,:,:) = -Gfull%wks(1,2,:,:,:)
            !!      Gfull%wks(2,1,:,:,:) = -Gfull%wks(2,1,:,:,:)
            !!   endif
            !!   !
            !!   allocate(Gfull_R(Norb,Norb,Nmats,Nwig));Gfull_R=czero
            !!   call wannier_K2R(Lttc%Nkpt3,Lttc%kpt,Gfull%wks(:,:,:,:,1),Gfull_R)
            !!   call DeallocateField(Gfull)
            !!   write(*,*)"Gfull_AB(R) computed"
            !!   !
            !!   allocate(Gfull_R_halfkz(3*Norb,3*Norb,Nmats,Nwig));Gfull_R_halfkz=czero
            !!   !
            !!   do nx=minval(Nvecwig(1,:)),maxval(Nvecwig(1,:)),1
            !!      do ny=minval(Nvecwig(2,:)),maxval(Nvecwig(2,:)),1
            !!         !
            !!         !new iwig
            !!         iwig = find_vec([nx,ny,0],Nvecwig,hardstop=.false.)
            !!         if(iwig.eq.0)cycle
            !!         !
            !!         do nz=minval(Nvecwig(3,:)),maxval(Nvecwig(3,:)),1
            !!            !
            !!            if(nz.eq.0)then
            !!               !diagonals
            !!               do ik=1,3
            !!                  Gfull_R_halfkz(1+(ik-1)*Norb:ik*Norb,1+(ik-1)*Norb:ik*Norb,:,iwig) = Gfull_R(:,:,:,iwig)
            !!               enddo
            !!            elseif(nz.eq.+1)then
            !!               !distance +1
            !!               iwig_old = find_vec([nx,ny,+1],Nvecwig,hardstop=.false.)
            !!               if(iwig_old.eq.0)cycle
            !!               !
            !!               Gfull_R_halfkz(1+Norb:2*Norb,1:Norb,:,iwig) = Gfull_R(:,:,:,iwig_old)
            !!               Gfull_R_halfkz(1+2*Norb:3*Norb,1+Norb:2*Norb,:,iwig) = Gfull_R(:,:,:,iwig_old)
            !!            elseif(nz.eq.-1)then
            !!               !distance -1
            !!               iwig_old = find_vec([nx,ny,-1],Nvecwig,hardstop=.false.)
            !!               if(iwig_old.eq.0)cycle
            !!               !
            !!               Gfull_R_halfkz(1:Norb,1+Norb:2*Norb,:,iwig) = Gfull_R(:,:,:,iwig_old)
            !!               Gfull_R_halfkz(1+Norb:2*Norb,1+2*Norb:3*Norb,:,iwig) = Gfull_R(:,:,:,iwig_old)
            !!            elseif(nz.eq.+2)then
            !!               !distance +2
            !!               iwig_old = find_vec([nx,ny,+2],Nvecwig,hardstop=.false.)
            !!               if(iwig_old.eq.0)cycle
            !!               !
            !!               Gfull_R_halfkz(1+2*Norb:3*Norb,1:Norb,:,iwig) = Gfull_R(:,:,:,iwig_old)
            !!            elseif(nz.eq.-2)then
            !!               !distance -2
            !!               iwig_old = find_vec([nx,ny,+2],Nvecwig,hardstop=.false.)
            !!               if(iwig_old.eq.0)cycle
            !!               !
            !!               Gfull_R_halfkz(1:Norb,1+2*Norb:3*Norb,:,iwig) = Gfull_R(:,:,:,iwig_old)
            !!            endif
            !!            !
            !!         enddo
            !!      enddo
            !!   enddo
            !!   deallocate(Gfull_R)
            !!   write(*,*)"extraction to hetero-like"
            !!   !
            !!   Nkp = 151
            !!   Nkz = 50
            !!   !
            !!   allocate(Gfull_R(3*Norb,3*Norb,Nmats,Lttc%Nkpt_path));Gfull_R=czero
            !!   call wannier_R2K(Lttc%Nkpt3,Lttc%kptpath(:,1:Nkp),Gfull_R_halfkz,Gfull_R(:,:,:,1:Nkp))
            !!   deallocate(Gfull_R_halfkz)
            !!   write(*,*)"path-interpolation"
            !!   !
            !!   allocate(Gfull_R_halfkz(1,Nmats,Nkp,0:Nkz));Gfull_R_halfkz=czero
            !!   do ik=1,Nkp
            !!      do ikz=0,Nkz
            !!         do isite=1,6
            !!            do jsite=1,6
            !!               !
            !!               kR = 2*pi * Lttc%kptpath(3,Nkp+ikz) * (isite-jsite)
            !!               cfac = dcmplx(cos(kR),+sin(kR))
            !!               !
            !!               Gfull_R_halfkz(1,:,ik,ikz) = Gfull_R_halfkz(1,:,ik,ikz) + Gfull_R(isite,jsite,:,ik)*cfac / 6
            !!               !
            !!            enddo
            !!         enddo
            !!      enddo
            !!   enddo
            !!   deallocate(Gfull_R)
            !!   write(*,*)"Gamma-A filled-1"
            !!   !
            !!   do ispin=1,Nspin
            !!      do iorb=1,2
            !!         do ik=1,Nkp
            !!            Gpath%wks(iorb,iorb,:,ik,ispin) = Gfull_R_halfkz(1,:,ik,0)
            !!         enddo
            !!         do ikz=0,Nkz
            !!            Gpath%wks(iorb,iorb,:,Nkp+ikz,ispin) = Gfull_R_halfkz(1,:,Nkp,ikz)
            !!         enddo
            !!      enddo
            !!   enddo
            !!   deallocate(Gfull_R_halfkz)
            !!   write(*,*)"Gamma-A filled-2"
            !!   !
            !!endif
            !
            !>>>TEST
            !
            call DeallocateFermionicField(Sloc)
            !
         endif
         !
         !
         !
         !-------------- Green's function along the {kx,ky} plane -------------!
         if(print_plane_G)then
            !
            call AllocateFermionicField(Sloc,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            do ik=1,Lttc%Nkpt_Plane
               Sloc%wks(:,:,:,ik,:) = Sfull%ws
            enddo
            !
            call AllocateFermionicField(Gfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
            call calc_Gmats(Gfermi,Lttc,Smats=Sloc,along_plane=.true.)
            !
            call DeallocateFermionicField(Sloc)
            !
         endif
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
      use input_vars, only : PadeWlimit, Nreal, wrealMax
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
      !Pade
      integer                               :: PadeWlimit_ndx
      real(8),allocatable                   :: wreal(:)
      complex(8),allocatable                :: Gpade(:,:)
      !
      logical                               :: printGmats,doPade=.false.
      !
      !
      if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- dump_MaxEnt_on_G_K"
      !
      !
      NtauFT = Ntau_MaxEnt
      !
      printGmats = Nmats_MaxEnt.gt.0
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
            write(*,"(A,I)") "     dump_MaxEnt_on_G_K: G full. Total number of K-points in the BZ:",Nkpt
            !
         case("path")
            !
            Nkpt = Lttc%Nkpt_path
            write(*,"(A,I)") "     dump_MaxEnt_on_G_K: G path. Total number of K-points along path:",Nkpt
            !
         case("plane")
            !
            Nkpt = Lttc%Nkpt_Plane
            write(*,"(A,I)") "     dump_MaxEnt_on_G_K: G plane. Total number of K-points in the {kx,ky} sheet:",Nkpt
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
            !print all diagonal elements on imaginary time axis
            path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do itau=1,NtauFT
                write(unit,"(200E20.12)") tau(itau),(dreal(Gitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
            enddo
            close(unit)
            !
            !print all diagonal elements on imaginary frequency axis
            if(printGmats)then
               path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_w_k"//str(ik)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats_MaxEnt
                   write(unit,"(200E20.12)") wmats(iw),(Gmats_diag(iorb,iw,ik,ispin),iorb=1,Norb)
               enddo
               close(unit)
            endif
            !
            !print all diagonal elements on real frequency axis
            if((PadeWlimit.gt.0d0).and.doPade)then
               !
               allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
               PadeWlimit_ndx = minloc(abs(wmats-PadeWlimit),dim=1)
               !
               allocate(Gpade(Norb,Nreal));Gpade=czero
               do iorb=1,Norb
                  Gpade(iorb,:) = pade(Gmats_diag(iorb,:,ik,ispin),"Fermionic",wlimit=PadeWlimit_ndx)
               enddo
               !
               path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_w_k"//str(ik)//".DAT_pade.dat"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nreal
                   write(unit,"(200E20.12)") wreal(iw),(abs(dimag(Gpade(iorb,iw))),iorb=1,Norb)
               enddo
               close(unit)
               !
               deallocate(wreal,Gpade)
               !
            endif
            !
            if(Norb.gt.1)then
               !
               !end-point correction if -G(0)-G(beta) != -1
               Gmax = - (dreal(Gitau_trace(ik,1,ispin)) + dreal(Gitau_trace(ik,NtauFT,ispin)))
               Gitau_trace(ik,1,ispin) = Gitau_trace(ik,1,ispin)/abs(Gmax)
               Gitau_trace(ik,NtauFT,ispin) = Gitau_trace(ik,NtauFT,ispin)/abs(Gmax)
               !
               !print trace on imaginary time axis
               path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//"_Tr.DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do itau=1,NtauFT
                   write(unit,"(200E20.12)") tau(itau),dreal(Gitau_trace(ik,itau,ispin))
               enddo
               close(unit)
               !
               !print trace on imaginary frequency axis
               if(printGmats)then
                  path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_w_k"//str(ik)//"_Tr.DAT"
                  unit = free_unit()
                  open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
                  do iw=1,Nmats_MaxEnt
                      write(unit,"(200E20.12)") wmats(iw),Gmats_trace(ik,iw,ispin)
                  enddo
                  close(unit)
               endif
               !
               !print all diagonal elements on real frequency axis
               if((PadeWlimit.gt.0d0).and.doPade)then
                  !
                  allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
                  PadeWlimit_ndx = minloc(abs(wmats-PadeWlimit),dim=1)
                  !
                  allocate(Gpade(1,Nreal));Gpade=czero
                  Gpade(1,:) = pade(Gmats_trace(ik,:,ispin),"Fermionic",wlimit=PadeWlimit_ndx)
                  !
                  path = reg(pathOUTPUT)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_w_k"//str(ik)//"_Tr.DAT_pade.dat"
                  unit = free_unit()
                  open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
                  do iw=1,Nreal
                      write(unit,"(200E20.12)") wreal(iw),abs(dimag(Gpade(1,iw)))
                  enddo
                  close(unit)
                  !
                  deallocate(wreal,Gpade)
                  !
               endif
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
               !
               ikx = int(ik/(Nkpt_plane+0.001))+1
               iky = ik - (ikx-1)*Nkpt_plane
               !
               k1 = Lttc%kptPlane(1,ik)
               k2 = Lttc%kptPlane(2,ik)
               k3 = Lttc%kptPlane(3,ik)
               !
               Bx = k1*Blat(1,1) + k2*Blat(1,2) + k3*Blat(1,3) ; if(ik.eq.1) Bx_old = Bx
               By = k1*Blat(2,1) + k2*Blat(2,2) + k3*Blat(2,3)
               Bz = k1*Blat(3,1) + k2*Blat(3,2) + k3*Blat(3,3)
               !
               write(unit,"(3I10,200E20.12)") ik,ikx,iky,k1,k2,k3,Bx,By,Bz,(Ak(ik,iorb),iorb=1,Norb)
               if(iky.eq.Nkpt_plane)write(unit,*)
               !
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
      if(Hetero%status.and.Hetero%fill_Gamma_A)then
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

subroutine interpolate2kpath_Bosonic(Wfull,Lttc,pathOUTPUT,name,mode,remove_Gamma)
   !
   use parameters
   use utils_misc
   use utils_fields
   use linalg, only : eigh, inv, zeye, det, trace
   use crystal
   use file_io
   use greens_function, only : calc_Gmats
   use fourier_transforms
   use input_vars, only : Nkpt_path, Nkpt_plane
   use input_vars, only : print_path_Chi, print_path_W, print_path_E
   implicit none
   !
   type(BosonicField),intent(in)         :: Wfull
   type(Lattice),intent(inout)           :: Lttc
   character(len=*),intent(in)           :: pathOUTPUT
   character(len=*),intent(in)           :: name
   character(len=*),intent(in)           :: mode
   logical,intent(in),optional           :: remove_Gamma
   !
   complex(8),allocatable                :: W_orig(:,:,:,:),W_intp(:,:,:,:)
   complex(8),allocatable                :: Wgamma(:,:,:)
   complex(8),allocatable                :: TrW_orig(:,:),TrW_intp(:,:)
   complex(8),allocatable                :: TrWgamma(:)
   real(8),allocatable                   :: wmats(:)
   !
   integer                               :: Norb,Nmats,unit,Wdim
   integer                               :: iq,iw,iorb,jorb,ib1,ib2
   character(len=256)                    :: path
   logical                               :: print_path,print_plane
   logical                               :: remove_Gamma_
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
   print_path = print_path_Chi .or. print_path_W .or. print_path_E
   !
   !flag removed from input_vars
   print_plane = .false. !print_plane_W
   !
   call createDir(reg(pathOUTPUT),verb=verbose)
   !
   !recalculate the internal K-meshes if not already stored
   if(print_path.and.(.not.Lttc%pathStored))then
      call interpolate2Path(Lttc,Nkpt_path,"Hk",store=.true.)
      call dump_Hk(Lttc%Hk_path,Lttc%kptpath,reg(pathOUTPUT),"Hk_path.DAT")
   endif
   if(print_plane.and.(.not.Lttc%planeStored))then
      call interpolate2Plane(Lttc,Nkpt_plane,"Hk",store=.true.)
   endif
   !
   ! Different modes explanation
   ! NaNa:          prints the interpolation of all the density-density elements, namely W(aa)(aa)
   ! Avg_NaNa:      prints the average of all the density-density elements, namely sum_a W(aa)(aa)/Norb
   ! Na:            prints for each index (aa) the average of all the off-diagonal density-density elements, namely sum_b W(aa)(bb)/Norb
   ! EigvProd_NaNb: prints the product of the eigenvalues of the (aa)(bb) elements, namely det{ W(aa)(bb) }
   ! EigvProd:      prints the product of the eigenvalues of the (ab)(cd) elements, namely det{ W(ab)(cd) }. Only for small orbital spaces.
   !
   select case(reg(mode))
      case default
         stop "interpolate2kpath_Bosonic: Available modes are: NaNa, Avg_NaNa, Na, EigvProd."
      case("NaNa","Avg_NaNa","Na","EigvProd_NaNb")
         !
         Wdim = Norb
         !
         allocate(W_orig(Wdim,Wdim,Nmats,Lttc%Nkpt));W_orig=czero
         allocate(Wgamma(Wdim,Wdim,Nmats));Wgamma=czero
         !
         do iorb=1,Norb
            do jorb=1,Norb
               call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
               W_orig(iorb,jorb,:,:) = Wfull%screened(ib1,ib2,:,:)
            enddo
         enddo
         !
         if(remove_Gamma_) then
            do iorb=1,Norb
               do jorb=1,Norb
                  call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                  Wgamma(iorb,jorb,:) = Wfull%screened(ib1,ib2,:,Wfull%iq_gamma)
                  do iq=1,Wfull%Nkpt
                     W_orig(iorb,jorb,:,iq) = W_orig(iorb,jorb,:,iq) - Wgamma(iorb,jorb,:)
                  enddo
               enddo
            enddo
         endif
         !
      case("EigvProd")
         !
         Wdim = Wfull%Nbp
         !
         allocate(W_orig(Wdim,Wdim,Nmats,Lttc%Nkpt));W_orig=czero
         allocate(Wgamma(Wdim,Wdim,Nmats));Wgamma=czero
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
   end select
   !
   !
   !
   !-------------------------- Boson along the path ---------------------------!
   if(print_path)then
      !
      allocate(wmats(Nmats));wmats = BosonicFreqMesh(Wfull%Beta,Nmats)
      !
      select case(reg(mode))
         !
         case("NaNa")
            !
            !Interpolate the boson in NaNb format along the path
            call cpu_time(start)
            allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_path));W_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),W_orig,W_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_NaNb(fullBZ,iw) --> "//reg(name)//"_NaNb(Kpath,iw) cpu timing:", finish-start
            !
            !Print along path all NaNb components
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_path/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_path
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_NaNb.DAT",form="formatted",status="unknown",position="rewind",action="write")
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
            !
            deallocate(W_intp)
            !
         case("Avg_NaNa")
            !
            !Print along path only the trace of the NaNa (diagonal) elements
            allocate(TrW_orig(Nmats,Wfull%Nkpt));TrW_orig=czero
            allocate(TrWgamma(Nmats));TrWgamma=czero
            do iw=1,Nmats
               do iq=1,Wfull%Nkpt
                  do iorb=1,Norb
                     TrW_orig(iw,iq) = TrW_orig(iw,iq) + W_orig(iorb,iorb,iw,iq)/Norb
                  enddo
               enddo
               do iorb=1,Norb
                  TrWgamma(iw) = TrWgamma(iw) + Wgamma(iorb,iorb,iw)/Norb
               enddo
            enddo
            !
            !Interpolate the boson trace along the path
            call cpu_time(start)
            allocate(TrW_intp(Nmats,Lttc%Nkpt_path));TrW_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),TrW_orig,TrW_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_Avg_NaNa(fullBZ,iw) --> "//reg(name)//"_Avg_NaNa(Kpath,iw) cpu timing:", finish-start
            !
            !Print along path the components of the trace
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_path/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_path
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_Avg_NaNa.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats-1
                  write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)+TrWgamma(iw)),dimag(TrW_intp(iw,iq)+TrWgamma(iw))
               enddo
               close(unit)
            enddo
            deallocate(TrW_orig,TrW_intp,TrWgamma)
            !
         case("Na")
            !
            !Print along path the sum of all the off-diag components for each orbial
            do jorb=1,Norb
               !
               allocate(TrW_orig(Nmats,Wfull%Nkpt));TrW_orig=czero
               allocate(TrWgamma(Nmats));TrWgamma=czero
               do iw=1,Nmats
                  do iq=1,Wfull%Nkpt
                     do iorb=1,Norb
                        TrW_orig(iw,iq) = TrW_orig(iw,iq) + W_orig(jorb,iorb,iw,iq)/Norb
                     enddo
                  enddo
                  do iorb=1,Norb
                     TrWgamma(iw) = TrWgamma(iw) + Wgamma(jorb,iorb,iw)/Norb
                  enddo
               enddo
               !
               !Interpolate the boson trace along the path
               call cpu_time(start)
               allocate(TrW_intp(Nmats,Lttc%Nkpt_path));TrW_intp=czero
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),TrW_orig,TrW_intp)
               call cpu_time(finish)
               write(*,"(A,F)") "     "//reg(name)//"_Na_o"//str(jorb)//"(fullBZ,iw) --> "//reg(name)//"_Na_o"//str(jorb)//"(Kpath,iw) cpu timing:", finish-start
               !
               !Print along path the components of the trace
               path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_path/"
               call createDir(reg(path),verb=verbose)
               do iq=1,Lttc%Nkpt_path
                  unit = free_unit()
                  open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_o"//str(jorb)//"_Na.DAT",form="formatted",status="unknown",position="rewind",action="write")
                  do iw=1,Nmats-1
                     write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)+TrWgamma(iw)),dimag(TrW_intp(iw,iq)+TrWgamma(iw))
                  enddo
                  close(unit)
               enddo
               deallocate(TrW_orig,TrW_intp,TrWgamma)
               !
            enddo
            !
         case("EigvProd_NaNb")
            !
            !Interpolate the full boson along the path
            call cpu_time(start)
            allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_path));W_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),W_orig,W_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_EigvProd_NaNb(fullBZ,iw) --> "//reg(name)//"_EigvProd_NaNb(Kpath,iw) ["//str(ib1)//","//str(ib2)//"] cpu timing:", finish-start
            !
            !extract the eigenvalues along the path
            allocate(TrW_intp(Nmats,Lttc%Nkpt_path));TrW_intp=czero
            do iw=1,Nmats
               do iq=1,Lttc%Nkpt_path
                  !
                  !get the eigenvalues of the full Field with also Gamma
                  TrW_intp(iw,iq) = det( W_intp(:,:,iw,iq) + Wgamma(:,:,iw) )
                  !
               enddo
            enddo
            !
            !Print along path the components of the EigvProd function
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_path/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_path
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_EigvProd_NaNb.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats-1
                  write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)),dimag(TrW_intp(iw,iq))
               enddo
               close(unit)
            enddo
            deallocate(TrW_intp)
            !
         case("EigvProd")
            !
            !Interpolate the full boson along the path
            call cpu_time(start)
            allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_path));W_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath(:,1:Lttc%Nkpt_path),W_orig,W_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_EigvProd(fullBZ,iw) --> "//reg(name)//"_EigvProd(Kpath,iw) ["//str(ib1)//","//str(ib2)//"] cpu timing:", finish-start
            !
            !extract the eigenvalues along the path
            allocate(TrW_intp(Nmats,Lttc%Nkpt_path));TrW_intp=czero
            do iw=1,Nmats
               do iq=1,Lttc%Nkpt_path
                  !
                  !get the eigenvalues of the full Field with also Gamma
                  TrW_intp(iw,iq) = det( W_intp(:,:,iw,iq) + Wgamma(:,:,iw) )
                  !
               enddo
            enddo
            !
            !Print along path the components of the EigvProd function
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_path/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_path
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_EigvProd.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats-1
                  write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)),dimag(TrW_intp(iw,iq))
               enddo
               close(unit)
            enddo
            deallocate(TrW_intp)
            !
            !
            !
      end select
      !
      deallocate(wmats)
      !
   endif
   !
   !
   !
   !---------------------- Boson along the {kx,ky} plane ----------------------!
   if(print_plane)then
      !
      allocate(wmats(Nmats));wmats = BosonicFreqMesh(Wfull%Beta,Nmats)
      !
      select case(reg(mode))
         !
         case("NaNa")
            !
            !Interpolate the boson in NaNb format along the path
            call cpu_time(start)
            allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_Plane));W_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,W_orig,W_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_NaNb(fullBZ,iw) --> "//reg(name)//"_NaNb(kx,ky,iw) cpu timing:", finish-start
            !
            !Print along path all NaNb components
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_plane/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_Plane
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_NaNb.DAT",form="formatted",status="unknown",position="rewind",action="write")
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
            !
            deallocate(W_intp)
            !
         case("Avg_NaNa")
            !
            !Print along path only the trace of the NaNa (diagonal) elements
            allocate(TrW_orig(Nmats,Wfull%Nkpt));TrW_orig=czero
            allocate(TrWgamma(Nmats));TrWgamma=czero
            do iw=1,Nmats
               do iq=1,Wfull%Nkpt
                  do iorb=1,Norb
                     TrW_orig(iw,iq) = TrW_orig(iw,iq) + W_orig(iorb,iorb,iw,iq)/Norb
                  enddo
               enddo
               do iorb=1,Norb
                  TrWgamma(iw) = TrWgamma(iw) + Wgamma(iorb,iorb,iw)/Norb
               enddo
            enddo
            !
            !Interpolate the boson trace along the path
            call cpu_time(start)
            allocate(TrW_intp(Nmats,Lttc%Nkpt_Plane));TrW_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,TrW_orig,TrW_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_Avg_NaNa(fullBZ,iw) --> "//reg(name)//"_Avg_NaNa(kx,ky,iw) cpu timing:", finish-start
            !
            !Print along path the components of the trace
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_plane/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_Plane
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_Avg_NaNa.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats-1
                  write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)+TrWgamma(iw)),dimag(TrW_intp(iw,iq)+TrWgamma(iw))
               enddo
               close(unit)
            enddo
            deallocate(TrW_orig,TrW_intp,TrWgamma)
            !
         case("Na")
            !
            !Print along path the sum of all the off-diag components for each orbial
            do jorb=1,Norb
               !
               allocate(TrW_orig(Nmats,Wfull%Nkpt));TrW_orig=czero
               allocate(TrWgamma(Nmats));TrWgamma=czero
               do iw=1,Nmats
                  do iq=1,Wfull%Nkpt
                     do iorb=1,Norb
                        TrW_orig(iw,iq) = TrW_orig(iw,iq) + W_orig(jorb,iorb,iw,iq)/Norb
                     enddo
                  enddo
                  do iorb=1,Norb
                     TrWgamma(iw) = TrWgamma(iw) + Wgamma(jorb,iorb,iw)/Norb
                  enddo
               enddo
               !
               !Interpolate the boson trace along the path
               call cpu_time(start)
               allocate(TrW_intp(Nmats,Lttc%Nkpt_Plane));TrW_intp=czero
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,TrW_orig,TrW_intp)
               call cpu_time(finish)
               write(*,"(A,F)") "     "//reg(name)//"_Na_o"//str(jorb)//"(fullBZ,iw) --> "//reg(name)//"_Na_o"//str(jorb)//"(kx,ky,iw) cpu timing:", finish-start
               !
               !Print along path the components of the trace
               path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_plane/"
               call createDir(reg(path),verb=verbose)
               do iq=1,Lttc%Nkpt_Plane
                  unit = free_unit()
                  open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_o"//str(jorb)//"_Na.DAT",form="formatted",status="unknown",position="rewind",action="write")
                  do iw=1,Nmats-1
                     write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)+TrWgamma(iw)),dimag(TrW_intp(iw,iq)+TrWgamma(iw))
                  enddo
                  close(unit)
               enddo
               deallocate(TrW_orig,TrW_intp,TrWgamma)
               !
            enddo
            !
         case("EigvProd_NaNb")
            !
            !Interpolate the full boson along the path
            call cpu_time(start)
            allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_Plane));W_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,W_orig,W_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_EigvProd_NaNb(fullBZ,iw) --> "//reg(name)//"_EigvProd_NaNb(kx,ky,iw) ["//str(ib1)//","//str(ib2)//"] cpu timing:", finish-start
            !
            !extract the eigenvalues along the path
            allocate(TrW_intp(Nmats,Lttc%Nkpt_Plane));TrW_intp=czero
            do iw=1,Nmats
               do iq=1,Lttc%Nkpt_Plane
                  !
                  !get the eigenvalues of the full Field with also Gamma
                  TrW_intp(iw,iq) = det( W_intp(:,:,iw,iq) + Wgamma(:,:,iw) )
                  !
               enddo
            enddo
            !
            !Print along path the components of the EigvProd function
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_plane/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_Plane
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_EigvProd_NaNb.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats-1
                  write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)),dimag(TrW_intp(iw,iq))
               enddo
               close(unit)
            enddo
            deallocate(TrW_intp)
            !
         case("EigvProd")
            !
            !Interpolate the full boson along the path
            call cpu_time(start)
            allocate(W_intp(Wdim,Wdim,Nmats,Lttc%Nkpt_Plane));W_intp=czero
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,W_orig,W_intp)
            call cpu_time(finish)
            write(*,"(A,F)") "     "//reg(name)//"_EigvProd(fullBZ,iw) --> "//reg(name)//"_EigvProd(kx,ky,iw) ["//str(ib1)//","//str(ib2)//"] cpu timing:", finish-start
            !
            !extract the eigenvalues along the path
            allocate(TrW_intp(Nmats,Lttc%Nkpt_Plane));TrW_intp=czero
            do iw=1,Nmats
               do iq=1,Lttc%Nkpt_Plane
                  !
                  !get the eigenvalues of the full Field with also Gamma
                  TrW_intp(iw,iq) = det( W_intp(:,:,iw,iq) + Wgamma(:,:,iw) )
                  !
               enddo
            enddo
            !
            !Print along path the components of the EigvProd function
            path = reg(pathOUTPUT)//"MaxEnt_"//reg(name)//"k_plane/"
            call createDir(reg(path),verb=verbose)
            do iq=1,Lttc%Nkpt_Plane
               unit = free_unit()
               open(unit,file=reg(path)//reg(name)//"k_w_k"//str(iq)//"_EigvProd.DAT",form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats-1
                  write(unit,"(200E20.12)") wmats(iw),dreal(TrW_intp(iw,iq)),dimag(TrW_intp(iw,iq))
               enddo
               close(unit)
            enddo
            deallocate(TrW_intp)
            !
            !
            !
      end select
      !
      deallocate(wmats)
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
