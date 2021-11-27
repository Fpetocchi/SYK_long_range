Vec_test = [4d0,5d0,4d0,1d0,2d0,1d0,3d0,4d0,9d0,4d0,2d0,3d0,3d0,4d0]
call get_pattern(Dist,Vec_test,eps,listDim=DistList,IncludeSingle=.true.)

write(*,"(A,100I3)")  "     Vec: ",int(Vec_test)
do iD=1,size(Dist,dim=1)
   write(*,"(A,1I3)")  new_line("A")//"     list row: ",iD
   write(*,"(A,10I3)") "     list col: ",(Dist(iD,iR),iR=1,DistList(iD))
   write(*,"(A,10I3)") "     list els: ",(int(Vec_test(Dist(iD,iR))),iR=1,DistList(iD))
enddo





subroutine build_Uret_singlParam_ph(Umats,Uaa,Uab,J,g_eph,wo_eph,LocalOnly)
   !
   use parameters
   use file_io
   use utils_misc
   use utils_fields
   use input_vars, only : Nreal, wrealMax, pathINPUTtr
   use input_vars, only : Hetero
   !TEST>>>
   use input_vars, only : tau_uniform, Ntau, Test_flag_3
   use fourier_transforms
   !>>>TEST
   implicit none
   !
   type(BosonicField),intent(inout),target :: Umats
   real(8),intent(in)                    :: Uaa,Uab,J
   real(8),intent(in)                    :: g_eph(:),wo_eph(:)
   logical,intent(in),optional           :: LocalOnly
   !
   integer                               :: Nbp,Norb,Nph
   integer                               :: ib1,ib2,ik
   integer                               :: iw,iw1,iw2,iph,iwp
   integer                               :: Nsite
   real(8)                               :: RealU,ImagU
   real(8),allocatable                   :: wreal(:),wmats(:)
   complex(8),allocatable                :: D(:,:),Utmp(:,:)
   type(BosonicField)                    :: Ureal
   type(BosonicField),target             :: Umats_imp
   type(BosonicField),pointer            :: Umats_ptr
   type(physicalU)                       :: PhysicalUelements
   logical                               :: LocalOnly_,Screen
   real                                  :: start,finish
   !TEST>>>
   type(BosonicField)                    :: Utmp1,Utmp2
   real(8)                               :: g,w,b
   integer                               :: itau
   real(8),allocatable                   :: tau(:)
   !>>>TEST
   !
   !
   if(verbose)write(*,"(A)") "---- build_Uret_singlParam_ph"
   !
   !
   ! Check on the input field
   if(.not.Umats%status) stop "build_Uret_singlParam_ph: BosonicField not properly initialized."
   !
   LocalOnly_=.true.
   if(present(LocalOnly))LocalOnly_=LocalOnly
   if(LocalOnly_.and.(Umats%Nkpt.ne.0)) stop "build_Uret_singlParam_ph: Umats k dependent attributes are supposed to be unallocated."
   !
   Nph = size(g_eph)
   if(size(g_eph).ne.size(wo_eph)) stop "build_Uret_singlParam_ph: Phonon sizes does not match."
   !
   allocate(wmats(Umats%Npoints));wmats=0d0
   wmats = BosonicFreqMesh(Umats%Beta,Umats%Npoints)
   allocate(wreal(Nreal));wreal=0d0
   wreal = linspace(0d0,+wrealMax,Nreal)
   !
   call clear_attributes(Umats)
   !
   Nbp = Umats%Nbp
   Norb = int(sqrt(dble(Nbp)))
   if(Hetero%status)then
      Norb = Hetero%Norb
      Nbp = Norb**2
      call AllocateBosonicField(Umats_imp,Norb,Umats%Npoints,Umats%iq_gamma,Nkpt=Umats%Nkpt,Beta=Umats%Beta)
      Umats_ptr => Umats_imp
   else
      Umats_ptr => Umats
   endif
   call AllocateBosonicField(Ureal,Norb,Nreal,0)
   !
   call init_Uelements(Norb,PhysicalUelements)
   !
   !setting the bare values
   do ib1=1,Nbp
      do ib2=1,Nbp
         !
         if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uaa,0d0)

         if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uab,0d0)
         if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J,0d0)
         if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J,0d0)
         !
         if(PhysicalUelements%Full_Uaa(ib1,ib2).and.(ib1.eq.1)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uaa,0d0) + 0.5d0
         !
      enddo
   enddo
   !
   !setting the phonons
   do ib1=1,Nbp
      do ib2=1,Nbp
         !
         Screen = .not. (PhysicalUelements%Full_Jsf(ib1,ib2).or.PhysicalUelements%Full_Jph(ib1,ib2))
         !
         do iph=1,Nph
            iwp=minloc(abs(wreal-wo_eph(iph)),dim=1)
            do iw=1,Nreal
               !
               RealU = 2*(g_eph(iph)**2)*wo_eph(iph) / ( (wreal(iw)**2) - (wo_eph(iph)**2) )
               ImagU=0d0
               if(iw.eq.iwp) ImagU = -pi*(g_eph(iph)**2)/abs(wreal(3)-wreal(2))
               !
               !TEST>>>
               if(Test_flag_3)then
                  g = 1d0
                  if(PhysicalUelements%Full_Uaa(ib1,ib2).and.(ib1.eq.4)) g = 1.185d0
                  if(iw.eq.iwp) ImagU = -pi*(g**2)/abs(wreal(3)-wreal(2))
               endif
               !>>>TEST
               !
               Ureal%screened_local(ib1,ib2,iw) = Umats_ptr%bare_local(ib1,ib2)
               if(Screen) Ureal%screened_local(ib1,ib2,iw) = Ureal%screened_local(ib1,ib2,iw) + dcmplx(RealU,ImagU)
               !
            enddo
         enddo
      enddo
   enddo
   !
   ! Allocate the temporary quantities needed by the Analytical continuation
   allocate(Utmp(Nbp,Nbp));Utmp=czero
   allocate(D(Nbp,Nbp));D=czero
   !
   ! Analytical continuation of the local component to imag axis using spectral rep
   call cpu_time(start)
   !$OMP PARALLEL DEFAULT(NONE),&
   !$OMP SHARED(Nbp,wmats,wreal,Nreal,Ureal,Umats_ptr),&
   !$OMP PRIVATE(ib1,ib2,iw1,iw2,D,Utmp)
   !$OMP DO
   do iw1=1,Umats_ptr%Npoints
      !
      Utmp=czero
      !
      do iw2=1,Nreal
         !
         if((wmats(iw1).eq.0d0).and.(wreal(iw2).eq.0d0))cycle
         !
         D = -dimag( Ureal%screened_local(:,:,iw2)*abs(wreal(3)-wreal(2)) )/pi
         Utmp = Utmp +  D/( dcmplx(0.d0,wmats(iw1))-wreal(iw2) ) - D/( dcmplx(0.d0,wmats(iw1))+wreal(iw2) )
         !
      enddo
      !
      Umats_ptr%screened_local(:,:,iw1) = Utmp + Umats_ptr%bare_local
      !
   enddo !iw1
   !
   !$OMP END DO
   !$OMP END PARALLEL
   call cpu_time(finish)
   deallocate(D,Utmp,wmats,wreal)
   call DeallocateBosonicField(Ureal)
   write(*,"(A,F)") "     Ue-ph(w) --> Ue-ph(iw) cpu timing:", finish-start
   if(Norb.gt.1)write(*,"(A)") "     Screening not considered for Hund coupling."
   !
   if(.not.LocalOnly_)then
      write(*,"(A)") "     Filling the K-dependent attributes."
      do ik=1,Umats%Nkpt
         Umats_ptr%bare(:,:,ik) = Umats_ptr%bare_local
         Umats_ptr%screened(:,:,:,ik) = Umats_ptr%screened_local
      enddo
   endif
   !
   if(Hetero%status)then
      Nsite = Hetero%Explicit(2)-Hetero%Explicit(1)+1
      call Expand2Nsite(Umats,Umats_ptr,Nsite)
      call DeallocateBosonicField(Umats_imp)
   endif
   nullify(Umats_ptr)
   !
   call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
   !
   !TEST>>>
   allocate(tau(Ntau));tau=0d0
   if(tau_uniform)tau = linspace(0d0,Umats%Beta,Ntau)
   if(.not.tau_uniform)tau = denspace(Umats%Beta,Ntau)
   !
   !FT from U(iw) computed with AC to numerical U(tau)
   call AllocateBosonicField(Utmp1,Norb,Ntau,0,Beta=Umats%Beta)
   call Bmats2itau(Umats%Beta,Umats%screened_local,Utmp1%screened_local,asympt_corr=.true.,tau_uniform=tau_uniform,Umats_bare=Umats%bare_local)
   call dump_BosonicField(Utmp1,reg(pathINPUTtr),"Uloc_tau_FTofAC.DAT",axis=tau)
   !
   analytical U(tau)
   g = g_eph(1)
   w = wo_eph(1)
   b = Umats%Beta/2
   call clear_attributes(Utmp1)
   do itau=1,Ntau
      Utmp1%screened_local(:,:,itau) = -g**2 * cosh(w*(b-tau(itau))) / ( sinh(w*b) )
   enddo
   call dump_BosonicField(Utmp1,reg(pathINPUTtr),"Uloc_tau_analytic.DAT",axis=tau)
   !
   !FT from analytical U(tau) to U(iw)
   call AllocateBosonicField(Utmp2,Norb,Umats%Npoints,0,Beta=Umats%Beta)
   call Bitau2mats(Umats%Beta,Utmp1%screened_local,Utmp2%screened_local,tau_uniform=tau_uniform)
   do iw=1,Utmp2%Npoints
      Utmp2%screened_local(:,:,iw) = Utmp2%screened_local(:,:,iw) + Umats%bare_local
   enddo
   Utmp2%bare_local = Umats%bare_local
   call dump_BosonicField(Utmp2,reg(pathINPUTtr),"Uloc_mats_FTofAnalytical.DAT")
   !
   call duplicate(Umats,Utmp2)
   !
   call DeallocateField(Utmp1)
   call DeallocateField(Utmp2)
   !>>>TEST
   !
end subroutine build_Uret_singlParam_ph

subroutine calc_QMCinteractions(Umats,Uinst,Kfunct,Kpfunct,Screening,sym)
   !
   use parameters
   use file_io
   use utils_misc
   use utils_fields
   use input_vars, only : Solver
   !TEST>>>
   use input_vars, only : g_eph,wo_eph,Test_flag_2
   !>>>TEST
   implicit none
   !
   type(BosonicField),intent(in)         :: Umats
   real(8),intent(inout)                 :: Uinst(:,:)
   real(8),intent(inout),optional        :: Kfunct(:,:,:)
   real(8),intent(inout),optional        :: Kpfunct(:,:,:)
   real(8),intent(inout),optional        :: Screening(:,:)
   logical,intent(in),optional           :: sym
   !
   integer                               :: Nbp,Norb,Nflavor
   integer                               :: ib1,ib2,iorb,jorb
   integer                               :: iu1,iu2,ix1,ix2,ip1,ip2
   integer                               :: iw,itau
   real(8),allocatable                   :: wmats(:),tau(:)
   complex(8),allocatable                :: Kaux(:,:,:)
   logical                               :: Uloc,U1st,U2nd,retarded,Kp,Scr
   type(physicalU)                       :: PhysicalUelements
   logical                               :: sym_
   !TEST>>>
   real(8)                               :: g,w,b
   !>>>TEST
   !
   !
   if(verbose)write(*,"(A)") "---- calc_QMCinteractions"
   !
   !
   if(.not.Umats%status) stop "calc_QMCinteractions: Umats not properly initialized."
   !
   retarded=.false.
   if(present(Kfunct))retarded=.true.
   !
   Kp=.false.
   if(present(Kpfunct).and.retarded)Kp=.true.
   !
   Scr=.false.
   if(present(Screening).and.retarded)Scr=.true.
   !
   sym_=.true.
   if(present(sym))sym_=sym
   !
   Nbp = Umats%Nbp
   Norb = int(sqrt(dble(Nbp)))
   Nflavor = Norb*Nspin
   !
   call init_Uelements(Norb,PhysicalUelements)
   !
   call assert_shape(Uinst,[Nflavor,Nflavor],"calc_QMCinteractions","Uinst")
   Uinst=0d0
   if(retarded)then
      call assert_shape(Kfunct,[Nflavor,Nflavor,Solver%NtauB],"calc_QMCinteractions","Kfunct")
      if(Kp)call assert_shape(Kpfunct,[Nflavor,Nflavor,Solver%NtauB],"calc_QMCinteractions","Kpfunct")
      if(Scr)call assert_shape(Screening,[Nflavor,Nflavor],"calc_QMCinteractions","Screening")
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
         ! (iorb,iorb)(jorb,jorb) indexes in the Norb^2 representaion
         call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],iu1,iu2)
         !
         ! (iorb,jorb)(jorb,iorb) indexes
         call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ix1,ix2)
         !
         ! (iorb,jorb)(iorb,jorb) indexes
         call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ip1,ip2)
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
         if(Scr)then
            !
            if(Uloc) Screening(ib1,ib2) =  Umats%bare_local(iu1,iu2) - Umats%screened_local(iu1,iu2,1)
            if(U1st) Screening(ib1,ib2) =  Umats%bare_local(iu1,iu2) - Umats%screened_local(iu1,iu2,1)
            if(U2nd) Screening(ib1,ib2) =  Umats%bare_local(iu1,iu2) - (Umats%bare_local(ix1,ix2)+Umats%bare_local(ip1,ip2))/2d0 - &
                                       (Umats%screened_local(iu1,iu2,1) - (Umats%screened_local(ix1,ix2,1)+Umats%screened_local(ip1,ip2,1))/2d0)
            !same orbital - same spin screening
            if(Uloc.and.(ib2.gt.ib1)) then
               Screening(ib1,ib1) = Screening(ib1,ib2)
               Screening(ib2,ib2) = Screening(ib1,ib2)
            endif
            !
         endif
         !
      enddo
   enddo
   if(sym_)call check_Symmetry(Uinst,eps,enforce=.true.,hardstop=.false.,name="Uinst")
   !
   !computing the retarded function
   if(retarded)then
      Kfunct=0d0
      do itau=2,Solver%NtauB-1
         do iw=2,Umats%Npoints
            Kfunct(:,:,itau) = Kfunct(:,:,itau) - 2d0*Kaux(:,:,iw) * ( cos(wmats(iw)*tau(itau)) - 1d0 ) / ( Umats%Beta*wmats(iw)**2 )
         enddo
         if(sym_)call check_Symmetry(Kfunct(:,:,itau),eps,enforce=.true.,hardstop=.false.,name="Kfunct_t"//str(itau))
      enddo
   endif
   !
   !computing the first derivative of retarded function
   if(retarded.and.kp)then
      Kpfunct=0d0
      do itau=2,Solver%NtauB-1
         do iw=2,Umats%Npoints
            Kpfunct(:,:,itau) = Kpfunct(:,:,itau) + 2d0*Kaux(:,:,iw) * sin(wmats(iw)*tau(itau)) / ( Umats%Beta*wmats(iw) )
         enddo
         if(sym_)call check_Symmetry(Kpfunct(:,:,itau),eps,enforce=.true.,hardstop=.false.,name="Kpfunct_t"//str(itau))
      enddo
   endif
   !
   !TEST>>>
   if(retarded.and.Test_flag_2)then
      write(*,"(A)") new_line("A")//new_line("A")//"---- calc_QMCinteractions: Analytical screening function."
      Kfunct=0d0
      g = g_eph(1)
      w = wo_eph(1)
      b = Umats%Beta/2
      do itau=2,Solver%NtauB-1
         Kfunct(:,:,itau) = -( g**2/w ) * ( cosh(w*(b-tau(itau))) - cosh(w*b) ) / ( sinh(w*b) )
      enddo
   endif
   !>>>TEST
   !
   if(retarded)deallocate(Kaux,tau,wmats)
   !
end subroutine calc_QMCinteractions

subroutine build_Hk(Norb,hopping,Nkpt3,alphaHk,readHr,Hetero,Hk,kpt,Ek,Zk,Hloc,iq_gamma,pathOUTPUT)
   !
   use utils_misc
   use parameters, only : Heterostructures !WHY IS THIS WORKING?
   use linalg, only : zeye, diagonal, diag, eigh, dag
   implicit none
   !
   integer,intent(in)                    :: Norb
   real(8),intent(in)                    :: hopping(:)
   integer,intent(in)                    :: Nkpt3(3)
   real(8),intent(in)                    :: alphaHk
   logical,intent(in)                    :: readHr
   type(Heterostructures),intent(inout)  :: Hetero
   complex(8),allocatable,intent(out)    :: Hk(:,:,:)
   real(8),allocatable,intent(out)       :: kpt(:,:)
   real(8),allocatable,intent(out)       :: Ek(:,:)
   complex(8),allocatable,intent(out)    :: Zk(:,:,:)
   complex(8),allocatable,intent(out)    :: Hloc(:,:)
   integer,intent(out),optional          :: iq_gamma
   character(len=*),intent(in),optional  :: pathOUTPUT
   !
   !User
   integer                               :: unit,Nkpt
   integer                               :: iwan1,iwan2,ik
   integer                               :: Trange,idist,iwig
   !W90
   integer,parameter                     :: W90NumCol=15
   integer                               :: Num_wann,Nrpts
   integer                               :: Qst,Rst,i,j,ir
   integer                               :: nx,ny,nz
   integer,allocatable                   :: Ndegen(:)
   real(8)                               :: ReHr,ImHr
   character(len=256)                    :: path
   logical                               :: filexists,Tcond
   !Hetero
   integer                               :: isite,Nsite,na,nb,ilayer
   logical,allocatable                   :: inHomo(:)
   real(8)                               :: tzRatio,angle,Rvec(3)
   real(8),allocatable                   :: Rsorted(:)
   integer,allocatable                   :: Rorder(:),itz(:)
   complex(8),allocatable                :: Hr(:,:,:),Hk_single(:,:,:),Hk_single_offdiag(:,:,:)
   !
   !
   if(verbose)write(*,"(A)") "---- build_Hk"
   !
   !
   if(.not.Lat_stored) stop "build_Hk: Lattice vectors not stored."
   if(readHr.and.(.not.present(pathOUTPUT))) stop "build_Hk: reading of Hr.DAT requested but missing path."
   !
   Nkpt = Nkpt3(1)*Nkpt3(2)*Nkpt3(3)
   call assert_shape(hopping,[Norb],"build_Hk","hopping")
   !
   if(allocated(Hk))deallocate(Hk)
   allocate(Hk(Norb,Norb,Nkpt));Hk=czero
   !
   if(allocated(kpt))deallocate(kpt)
   allocate(kpt(3,Nkpt));kpt=0d0
   !
   call build_kpt(Nkpt3,kpt,pathOUTPUT=reg(pathOUTPUT))
   !
   !recover the vectors in real space and allocate hopping in real space
   if(.not.Wig_stored)call calc_wignerseiz(Nkpt3)
   allocate(Rsorted(Nwig));Rsorted = radiuswig
   allocate(Rorder(Nwig))
   call sort_array(Rsorted,Rorder)
   allocate(Hr(Norb,Norb,Nwig));Hr=czero
   !
   if(readHr)then
      !
      ! Look for Hk.DAT
      path=reg(pathOUTPUT)//"Hr.DAT"
      call inquireFile(reg(path),filexists)
      !
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
      read(unit,*)                      !skip first line
      read(unit,*) Num_wann !Number of Wannier orbitals
      read(unit,*) Nrpts    !Number of Wigner-Seitz vectors
      !
      if(Num_wann.ne.Norb) stop "build_Hk: number of Wannier orbital in Hr.DAT and model orbital space does not coincide."
      !
      Qst = int(Nrpts/W90NumCol)
      Rst = mod(Nrpts,W90NumCol)
      !
      allocate(Ndegen(Nrpts));Ndegen=0
      do i=1,Qst
         read(unit,*)(Ndegen(j+(i-1)*W90NumCol),j=1,W90NumCol)
      enddo
      if(Rst.ne.0)read(unit,*)(Ndegen(j+Qst*W90NumCol),j=1,Rst)
      !
      !Read W90 TB hoppings in real space. Assumed paramagnetic
      do ir=1,Nrpts
         do i=1,Num_wann
            do j=1,Num_wann
               !
               read(unit,*) nx, ny, nz, iwan1, iwan2, ReHr, ImHr
               !
               iwig = find_vec([nx,ny,nz],Nvecwig)
               !
               Hr(iwan1,iwan2,iwig) = dcmplx(ReHr,ImHr)/Ndegen(ir)
               !nrdegwig(iwig) = Ndegen(ir) <-- this would mess-up things in the FT
               !
            enddo
         enddo
      enddo
      close(unit)
      deallocate(Ndegen)
      !
   else
      !
      !User-provided hopping is only nearest neighbor by now
      Trange=1
      !
      !loop over the sorted Wigner Seiz positions
      idist=1
      loopwigD:do iwig=1,Nwig
         !
         !setting the local energy
         if(Rsorted(Rorder(iwig)).eq.0d0)then
            if(Rorder(iwig).ne.wig0)stop "build_Hk: wrong index of R=0 vector."
            cycle
         endif
         !
         !increasing range
         if(iwig.gt.2)then
            if((Rsorted(Rorder(iwig))-Rsorted(Rorder(iwig-1))).gt.1e-5) idist=idist+1  !if(Rsorted(Rorder(iwig)).gt.Rsorted(Rorder(iwig-1))) idist=idist+1
            if(idist.gt.Trange) exit loopwigD
         endif
         !
         !setting matrix element
         do iwan1=1,Norb
            Hr(iwan1,iwan1,Rorder(iwig)) = -dcmplx(hopping(iwan1),0d0)
         enddo
         !
      enddo loopwigD
      !
   endif
   !
   if(verbose)then
      if(readHr)then
         write(*,'(1A)')        "     H_W90:"
         write(*,'(A,I6)')      "     Number of Wannier functions:   ",Num_wann
         write(*,'(A,I6)')      "     Number of Wigner-Seitz vectors:",Nrpts
         write(*,'(A,I6,A,I6)') "     Deg rows:",Qst," N last row   :",Rst
      endif
      write(*,'(1A)')"     Real-space hopping elements:"
      write(*,"(A6,3A12,1A4)") "  i  ","  Ri  ","  H(Ri)  "," [n1,n2,n3] "," Ndeg "
      do iwig=1,Nwig
         write(*,"(1I6,2F12.4,5I4)")Rorder(iwig),Rsorted(Rorder(iwig)),real(Hr(1,1,Rorder(iwig))),Nvecwig(:,Rorder(iwig)),nrdegwig(Rorder(iwig))
      enddo
   endif
   !
   !FT Hr-->Hk
   call wannier_R2K(Nkpt3,kpt,Hr,Hk)
   deallocate(Hr)
   !
   do ik=1,nkpt
      do iwan1=1,Norb
         Hk(iwan1,iwan1,ik) = dcmplx(dreal(Hk(iwan1,iwan1,ik)),0d0)
      enddo
      if(Norb.gt.1)call check_Hermiticity(Hk(:,:,ik),eps)
   enddo
   !
   !Build up the Heterostructure Hamiltonian
   Nsite = 1
   if(Hetero%status)then
      !
      !this should be already been checked in input_vars
      Nsite = Hetero%Explicit(2)-Hetero%Explicit(1)+1
      !
      !Setting up off-diagonal dispersion if requested
      if(Hetero%offDiagEk)then
         !
         allocate(Hk_single_offdiag(Norb,Norb,Nkpt));Hk_single_offdiag=czero
         allocate(Hr(Norb,Norb,Nwig));Hr=czero
         !
         !User-provided hopping is only nearest neighbor by now
         Trange=1
         !
         !loop over the sorted Wigner Seiz positions
         idist=1
         loopwigOD:do iwig=1,Nwig
            !
            !setting the local energy
            if(Rsorted(Rorder(iwig)).eq.0d0)then
               if(Rorder(iwig).ne.wig0)stop "build_Hk: wrong index of R=0 vector."
               cycle
            endif
            !
            !increasing range
            if(iwig.gt.2)then
               if((Rsorted(Rorder(iwig))-Rsorted(Rorder(iwig-1))).gt.1e-5) idist=idist+1
               if(idist.gt.Trange) exit loopwigOD
            endif
            !
            !setting matrix element
            !PROJECT SPECIFIC (TaS2)>>> The vertical hopping has only three next neighbor
            !do iwan1=1,Norb
            !   Hr(iwan1,iwan1,Rorder(iwig)) = dcmplx(1d0,0d0)
            !enddo
            Rvec = Nvecwig(1,Rorder(iwig))*Rlat(:,1) + Nvecwig(2,Rorder(iwig))*Rlat(:,2) + Nvecwig(3,Rorder(iwig))*Rlat(:,3)
            angle = atan2(Rvec(2),Rvec(1))
            if(angle.lt.0d0) angle = angle + 2d0*pi
            Tcond = (mod(nint(angle*180/pi)/60,2)-1) .eq. 0
            if(Tcond)then
               !write(*,*)angle,angle*180/pi,nint(angle*180/pi),mod(nint(angle*180/pi)/60,2),(mod(nint(angle*180/pi)/60,2)-1)
               !write(*,*)Nvecwig(:,Rorder(iwig))
               do iwan1=1,Norb
                  Hr(iwan1,iwan1,Rorder(iwig)) = dcmplx(1d0,0d0)
               enddo
            endif
            !>>>PROJECT SPECIFIC (TaS2)
            !
         enddo loopwigOD
         !
         call wannier_R2K(Nkpt3,kpt,Hr,Hk_single_offdiag)
         deallocate(Hr)
         !
      endif
      !
      !Setting up the out-of-plane hopping array
      Hetero%tzIndex(1) = Hetero%Explicit(1)
      Hetero%tzIndex(2) = Hetero%Explicit(2) - 1
      if(Hetero%Explicit(1).ne.1) Hetero%tzIndex(1) = Hetero%tzIndex(1) - 1              ! hopping to the left potential
      if(Hetero%Explicit(2).ne.Hetero%Nslab) Hetero%tzIndex(2) = Hetero%tzIndex(2) + 1   ! hopping to the right potential
      !
      allocate(Hetero%tz(Norb,Norb,Nkpt,Hetero%tzIndex(1):Hetero%tzIndex(2)));Hetero%tz=czero
      allocate(inHomo(Hetero%tzIndex(1):Hetero%tzIndex(2)));inHomo=.false.
      write(*,"(A)")new_line("A")//"     Hetero:"
      do ilayer = Hetero%tzIndex(1),Hetero%tzIndex(2)
         !
         inHomo(ilayer) = (Hetero%NtzExplicit.gt.0) !.and. any(Hetero%ExplicitTzPos.eq.ilayer)
         if(inHomo(ilayer)) inHomo(ilayer) = inHomo(ilayer) .and. any(Hetero%ExplicitTzPos.eq.ilayer)
         !
         tzRatio = 1d0
         if(inHomo(ilayer))then
            allocate(itz(Hetero%NtzExplicit));itz=0
            itz = findloc(Hetero%ExplicitTzPos,value=ilayer)
            if(itz(1).eq.0) stop "build_Hk: something wrong with the Hetero%ExplicitTzPos"
            tzRatio = Hetero%ExplicitTzRatios(itz(1))
            deallocate(itz)
         else
            tzRatio = Hetero%GlobalTzRatio
         endif
         !
         do ik=1,Nkpt
            !PROJECT SPECIFIC (TaS2)>>> The ihomogeneous vertical hopping is also without dispersion
            !if(Hetero%offDiagEk)then
            if(Hetero%offDiagEk.and.(.not.inHomo(ilayer)))then
            !>>>PROJECT SPECIFIC (TaS2)
               Hetero%tz(:,:,ik,ilayer) = matmul(diag(hopping)*tzRatio,Hk_single_offdiag(:,:,ik))
            else
               Hetero%tz(:,:,ik,ilayer) = diag(hopping)*tzRatio
            endif
         enddo
         !
         write(*,"(A,F)")"     tz/tplane ["//str(ilayer)//"-"//str(ilayer+1)//"]:",tzRatio
         !
      enddo
      !
      !Setting up multi-site H(k)
      allocate(Hk_single(Norb,Norb,Nkpt));Hk_single=czero
      Hk_single = Hk
      deallocate(Hk)
      allocate(Hk(Norb*Nsite,Norb*Nsite,Nkpt));Hk=czero
      !
      !adding non-diaongonal part
      do isite=1,Nsite
         !
         !Index of the layer inside the slab - Needed because Hetero%tz has a different indexing
         ilayer = Hetero%Explicit(1) + (isite-1)
         !
         !In-plane orbital block
         na = 1+(isite-1)*Norb
         nb = isite*Norb
         !
         !In-plane Hk
         Hk(na:nb,na:nb,:) = Hk_single
         !
         !Out-of-plane hopping
         if(isite.ne.Nsite)then
            do ik=1,Nkpt
               Hk(na:nb,na+Norb:nb+Norb,ik) = Hetero%tz(:,:,ik,ilayer)
               Hk(na+Norb:nb+Norb,na:nb,ik) = dag(Hk(na:nb,na+Norb:nb+Norb,ik))
            enddo
         endif
         !
      enddo
      deallocate(Hk_single,inHomo)
      !
   endif
   deallocate(Rorder,Rsorted)
   if(Hetero%offDiagEk)deallocate(Hk_single_offdiag)
   !
   Hk = Hk*alphaHk
   !
   if(allocated(Ek))deallocate(Ek)
   if(allocated(Zk))deallocate(Zk)
   if(allocated(Hloc))deallocate(Hloc)
   allocate(Ek(Norb*Nsite,Nkpt));Ek=0d0
   allocate(Zk(Norb*Nsite,Norb*Nsite,Nkpt));Zk=czero
   allocate(Hloc(Norb*Nsite,Norb*Nsite));Hloc=czero
   !
   do ik=1,Nkpt
      !
      call check_Hermiticity(Hk(:,:,ik),eps)
      !
      Ek(:,ik) = 0d0
      Zk(:,:,ik) = Hk(:,:,ik)
      call eigh(Zk(:,:,ik),Ek(:,ik))
      !
   enddo
   Hloc = sum(Hk,dim=3)/Nkpt
   !
   if(present(iq_gamma))iq_gamma = find_vec([0d0,0d0,0d0],kpt,eps)
   write(*,"(A)")"     Gamma point index: "//str(iq_gamma)
   !
   if(present(pathOUTPUT))then
      !
      unit = free_unit()
      open(unit,file=reg(pathOUTPUT)//"Hk.DAT",form="formatted",status="unknown",position="rewind",action="write")
      write(unit,("(3I10)")) 1,Nkpt,Norb
      do ik=1,Nkpt
         write(unit,("(3F14.8)")) kpt(:,ik)
         do iwan1=1,Norb*Nsite
            do iwan2=1,Norb*Nsite
               write(unit,("(2I4,2E20.12)")) iwan1,iwan2,dreal(Hk(iwan1,iwan2,ik)),dimag(Hk(iwan1,iwan2,ik))
            enddo
         enddo
      enddo
      !
   endif
   !
end subroutine build_Hk

subroutine build_Uret_singlParam_Vn(Umats,Uaa,Uab,J,Vnn,Lttc,LocalOnly)
   !
   use parameters
   use file_io
   use utils_misc
   use utils_fields
   use crystal
   use input_vars, only : pathINPUTtr, pathINPUT
   use input_vars, only : long_range, structure, Nkpt_path
   use input_vars, only : Hetero
   implicit none
   !
   type(BosonicField),intent(inout),target :: Umats
   real(8),intent(in)                    :: Uaa,Uab,J
   real(8),intent(in)                    :: Vnn(:,:)
   type(Lattice),intent(inout)           :: Lttc
   logical,intent(in),optional           :: LocalOnly
   !
   complex(8),allocatable                :: U_K(:,:,:)
   complex(8),allocatable                :: U_R(:,:,:)
   integer                               :: Nbp,Norb,Vrange
   integer                               :: ib1,ib2,iorb
   integer                               :: iwig,idist,Nsite
   real(8),allocatable                   :: Rsorted(:)
   integer,allocatable                   :: Rorder(:)
   type(physicalU)                       :: PhysicalUelements
   type(BosonicField),target             :: Umats_imp
   type(BosonicField),pointer            :: Umats_ptr
   complex(8),allocatable                :: EwaldShift(:)
   real(8),allocatable                   :: V(:)
   real(8)                               :: eta,den
   logical                               :: LocalOnly_,RealSpace
   real                                  :: start,finish
   !
   !
   if(verbose)write(*,"(A)") "---- build_Uret_singlParam_Vn"
   !
   !
   ! Check on the input field
   if(.not.Umats%status) stop "build_Uret_singlParam_Vn: BosonicField not properly initialized."
   if(Umats%Npoints.ne.1) stop "build_Uret_singlParam_Vn: Number of matsubara points in Umats is supposed to be equal to 1."
   !
   LocalOnly_=.false.
   if(present(LocalOnly))LocalOnly_=LocalOnly
   if(LocalOnly_.and.(Umats%Nkpt.ne.0)) stop "build_Uret_singlParam_Vn: Umats k dependent attributes are supposed to be unallocated."
   if((.not.LocalOnly_).and.(Umats%Nkpt.ne.Lttc%Nkpt)) stop "build_Uret_singlParam_Vn: Umats number of K-points does not match with the lattice."
   !
   Nbp = Umats%Nbp
   Norb = int(sqrt(dble(Nbp)))
   RealSpace = Lttc%Nsite.gt.1
   !
   if(RealSpace)then
      Norb = Lttc%Norb/Lttc%Nsite
      if(Hetero%status) Norb = Hetero%Norb
      Nbp = Norb**2
      call AllocateBosonicField(Umats_imp,Norb,Umats%Npoints,Umats%iq_gamma,Nkpt=Umats%Nkpt,Beta=Umats%Beta)
      Umats_ptr => Umats_imp
   else
      Umats_ptr => Umats
   endif
   !
   Vrange = size(Vnn,dim=2)
   call assert_shape(Vnn,[Norb,Vrange],"build_Uret_singlParam_Vn","Vnn")
   !
   call init_Uelements(Norb,PhysicalUelements)
   !
   call cpu_time(start)
   !
   !recover the vectors in real space and allocate interaction in real space
   if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
   allocate(Rsorted(Nwig));Rsorted = radiuswig
   allocate(Rorder(Nwig))
   call sort_array(Rsorted,Rorder)
   allocate(U_R(Nbp,Nbp,Nwig));U_R=czero
   allocate(V(Norb));V=czero
   if(reg(long_range).eq."Ewald")then
      eta = Rsorted(Rorder(Nwig))/2d0
      allocate(EwaldShift(Nwig));EwaldShift=czero
      if(any(Lttc%Nkpt3.eq.1))then
         call calc_Ewald(EwaldShift,Lttc%kpt,eta,"2D")
      else
         call calc_Ewald(EwaldShift,Lttc%kpt,eta,"3D")
      endif
   endif
   !
   !loop over the sorted Wigner Seiz positions
   idist=1
   loopwig:do iwig=1,Nwig
      !
      !setting the local interaction
      if(Rsorted(Rorder(iwig)).eq.0d0)then
         if(Rorder(iwig).ne.wig0)stop "build_Uret_singlParam_Vn: wrong index of R=0 vector."
         do ib1=1,Nbp
            do ib2=1,Nbp
               !
               if(PhysicalUelements%Full_Uaa(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(Uaa,0d0)
               if(PhysicalUelements%Full_Uab(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(Uab,0d0)
               if(PhysicalUelements%Full_Jsf(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(J,0d0)
               if(PhysicalUelements%Full_Jph(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(J,0d0)
               !
            enddo
         enddo
         !
         cycle
      endif
      !
      !increasing range
      if(iwig.gt.2)then
         if((Rsorted(Rorder(iwig))-Rsorted(Rorder(iwig-1))).gt.1e-5) idist=idist+1
         if((idist.gt.Vrange).and.(reg(long_range).ne."Ewald")) exit loopwig
      endif
      !
      !setting the R dependence
      if(reg(long_range).eq."Explicit")then
         V = Vnn(:,idist)
      elseif(reg(long_range).eq."Coulomb")then
         V = Vnn(:,1)/Rsorted(Rorder(iwig))
      elseif(reg(long_range).eq."Ewald")then
         den = 2d0*sqrt(eta)
         if(any(Lttc%Nkpt3.eq.1)) den = 2d0*eta
         V = (Vnn(:,1)/Rsorted(Rorder(iwig)))*erfc(Rsorted(Rorder(iwig))/den) + EwaldShift(Rorder(iwig))
      else
         stop "build_Uret_singlParam_Vn: the long_range varibale is not set."
      endif
      !
      !setting matrix element
      do ib1=1,Nbp
         iorb = PhysicalUelements%Full_Map(ib1,ib1,1)
         if(PhysicalUelements%Full_Uaa(ib1,ib1)) U_R(ib1,ib1,Rorder(iwig)) = dcmplx(V(iorb),0d0)
      enddo
      !
   enddo loopwig
   deallocate(V)
   if(allocated(EwaldShift))deallocate(EwaldShift)
   !
   if(verbose)then
      write(*,*)"     Real-space interaction elements:"
      write(*,"(A6,3A12)") "    i","    Ri","    H(Ri)","  [n1,n2,n3]"
      do iwig=1,Nwig
         write(*,"(1I6,2F12.4,3I4)")Rorder(iwig),Rsorted(Rorder(iwig)),real(U_R(1,1,Rorder(iwig))),Nvecwig(:,Rorder(iwig))
      enddo
   endif
   !
   !FT to K-space
   allocate(U_K(Nbp,Nbp,Lttc%Nkpt));U_K=czero
   if(Lttc%Nkpt.gt.1)then
      call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,U_R,U_K)
   else
      U_K(:,:,1) = U_R(:,:,wig0)
   endif
   deallocate(U_R,Rorder,Rsorted)
   !
   call cpu_time(finish)
   write(*,"(A,F)") "     Unn(R) --> Unn(K) cpu timing:", finish-start
   !
   if(reg(structure).ne."None")then
      call interpolateHk2Path(Lttc,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),filename="Uk",data_in=U_K)
   endif
   !
   !fill in the output
   do ib1=1,Nbp
      do ib2=1,Nbp
         !
         if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(Uaa,0d0)
         if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(Uab,0d0)
         if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(J,0d0)
         if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(J,0d0)
         !
         if(allocated(Umats_ptr%bare_local))then
            if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uaa,0d0)
            if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uab,0d0)
            if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J,0d0)
            if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J,0d0)
         endif
         !
      enddo
   enddo
   !
   if(.not.LocalOnly_)then
      !
      Umats_ptr%screened(:,:,1,:) = U_K
      !
      if(allocated(Umats_ptr%bare))then
         Umats_ptr%bare = U_K
      endif
      !
   endif
   deallocate(U_K)
   !
   if(RealSpace)then
      Nsite = Lttc%Nsite
      if(Hetero%status) Nsite = Hetero%Explicit(2)-Hetero%Explicit(1)+1
      call Expand2Nsite(Umats,Umats_ptr,Nsite)
      call DeallocateBosonicField(Umats_imp)
   endif
   nullify(Umats_ptr)
   !
   call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats_nosum.DAT")
   call BosonicKsum(Umats)
   call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
   !
end subroutine build_Uret_singlParam_Vn
!
subroutine build_Uret_multiParam_Vn(Umats,Uaa,Uab,J,Vnn,Lttc,LocalOnly)
   !
   use parameters
   use file_io
   use utils_misc
   use utils_fields
   use crystal
   use input_vars, only : pathINPUTtr, pathINPUT
   use input_vars, only : long_range, structure, Nkpt_path
   use input_vars, only : Hetero
   implicit none
   !
   type(BosonicField),intent(inout),target :: Umats
   real(8),intent(in)                    :: Uaa(:),Uab(:,:),J(:,:)
   real(8),intent(in)                    :: Vnn(:,:)
   type(Lattice),intent(inout)           :: Lttc
   logical,intent(in),optional           :: LocalOnly
   !
   complex(8),allocatable                :: U_K(:,:,:)
   complex(8),allocatable                :: U_R(:,:,:)
   integer                               :: Nbp,Norb,Vrange
   integer                               :: ib1,ib2,iorb,jorb
   integer                               :: iwig,idist,Nsite
   real(8),allocatable                   :: Rsorted(:)
   integer,allocatable                   :: Rorder(:)
   type(physicalU)                       :: PhysicalUelements
   type(BosonicField),target             :: Umats_imp
   type(BosonicField),pointer            :: Umats_ptr
   complex(8),allocatable                :: EwaldShift(:)
   real(8),allocatable                   :: V(:)
   real(8)                               :: eta,den
   logical                               :: LocalOnly_,RealSpace
   real                                  :: start,finish
   !
   !
   if(verbose)write(*,"(A)") "---- build_Uret_multiParam_Vn"
   !
   !
   ! Check on the input field
   if(.not.Umats%status) stop "build_Uret_multiParam_Vn: BosonicField not properly initialized."
   if(Umats%Npoints.ne.1) stop "build_Uret_multiParam_Vn: Number of matsubara points in Umats is supposed to be equal to 1."
   !
   LocalOnly_=.false.
   if(present(LocalOnly))LocalOnly_=LocalOnly
   if(LocalOnly_.and.(Umats%Nkpt.ne.0)) stop "build_Uret_multiParam_Vn: Umats k dependent attributes are supposed to be unallocated."
   if((.not.LocalOnly_).and.(Umats%Nkpt.ne.Lttc%Nkpt)) stop "build_Uret_multiParam_Vn: Umats number of K-points does not match with the lattice."
   !
   Nbp = Umats%Nbp
   Norb = int(sqrt(dble(Nbp)))
   RealSpace = Lttc%Nsite.gt.1
   !
   if(RealSpace)then
      Norb = Lttc%Norb/Lttc%Nsite
      if(Hetero%status) Norb = Hetero%Norb
      Nbp = Norb**2
      call AllocateBosonicField(Umats_imp,Norb,Umats%Npoints,Umats%iq_gamma,Nkpt=Umats%Nkpt,Beta=Umats%Beta)
      Umats_ptr => Umats_imp
   else
      Umats_ptr => Umats
   endif
   call assert_shape(Uaa,[Norb],"build_Uret_multiParam_Vn","Uaa")
   call assert_shape(Uab,[Norb,Norb],"build_Uret_multiParam_Vn","Uab")
   call assert_shape(J,[Norb,Norb],"build_Uret_multiParam_Vn","J")
   !
   Vrange = size(Vnn,dim=2)
   call assert_shape(Vnn,[Norb,Vrange],"build_Uret_multiParam_Vn","Vnn")
   !
   call init_Uelements(Norb,PhysicalUelements)
   !
   call cpu_time(start)
   !
   !recover the vectors in real space and allocate interaction in real space
   if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
   allocate(Rsorted(Nwig));Rsorted = radiuswig
   allocate(Rorder(Nwig))
   call sort_array(Rsorted,Rorder)
   allocate(U_R(Nbp,Nbp,Nwig));U_R=czero
   allocate(V(Norb));V=czero
   if(reg(long_range).eq."Ewald")then
      eta = Rsorted(Rorder(Nwig))/2d0
      allocate(EwaldShift(Nwig));EwaldShift=czero
      if(any(Lttc%Nkpt3.eq.1))then
         call calc_Ewald(EwaldShift,Lttc%kpt,eta,"2D")
      else
         call calc_Ewald(EwaldShift,Lttc%kpt,eta,"3D")
      endif
   endif
   !
   !loop over the sorted Wigner Seiz positions
   idist=1
   loopwig:do iwig=1,Nwig
      !
      !setting the local interaction
      if(Rsorted(Rorder(iwig)).eq.0d0)then
         if(Rorder(iwig).ne.wig0)stop "build_Uret_multiParam_Vn: wrong index of R=0 vector."
         do ib1=1,Nbp
            do ib2=1,Nbp
               !
               iorb = PhysicalUelements%Full_Map(ib1,ib2,1)
               jorb = PhysicalUelements%Full_Map(ib1,ib2,2)
               !
               if(PhysicalUelements%Full_Uaa(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(Uaa(iorb),0d0)
               if(PhysicalUelements%Full_Uab(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(Uab(iorb,jorb),0d0)
               if(PhysicalUelements%Full_Jsf(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(J(iorb,jorb),0d0)
               if(PhysicalUelements%Full_Jph(ib1,ib2)) U_R(ib1,ib2,Rorder(iwig)) = dcmplx(J(iorb,jorb),0d0)
               !
            enddo
         enddo
         !
         cycle
      endif
      !
      !increasing range
      if(iwig.gt.2)then
         if((Rsorted(Rorder(iwig))-Rsorted(Rorder(iwig-1))).gt.1e-5) idist=idist+1
         if((idist.gt.Vrange).and.(reg(long_range).ne."Ewald")) exit loopwig
      endif
      !
      !setting the R dependence
      if(reg(long_range).eq."Explicit")then
         V = Vnn(:,idist)
      elseif(reg(long_range).eq."Coulomb")then
         V = Vnn(:,1)/Rsorted(Rorder(iwig))
      elseif(reg(long_range).eq."Ewald")then
         den = 2d0*sqrt(eta)
         if(any(Lttc%Nkpt3.eq.1)) den = 2d0*eta
         V = (Vnn(:,1)/Rsorted(Rorder(iwig)))*erfc(Rsorted(Rorder(iwig))/den) + EwaldShift(Rorder(iwig))
      else
         stop "build_Uret_singlParam_Vn: the long_range varibale is not set."
      endif
      !
      !setting matrix element
      do ib1=1,Nbp
         iorb = PhysicalUelements%Full_Map(ib1,ib1,1)
         if(PhysicalUelements%Full_Uaa(ib1,ib1)) U_R(ib1,ib1,Rorder(iwig)) = dcmplx(V(iorb),0d0)
      enddo
      !
   enddo loopwig
   deallocate(V)
   if(allocated(EwaldShift))deallocate(EwaldShift)
   !
   if(verbose)then
      write(*,*)"     Real-space interaction elements:"
      write(*,"(A6,3A12)") "    i","    Ri","    H(Ri)","  [n1,n2,n3]"
      do iwig=1,Nwig
         write(*,"(1I6,2F12.4,3I4)")Rorder(iwig),Rsorted(Rorder(iwig)),real(U_R(1,1,Rorder(iwig))),Nvecwig(:,Rorder(iwig))
      enddo
   endif
   !
   !FT to K-space
   allocate(U_K(Nbp,Nbp,Lttc%Nkpt));U_K=czero
   if(Lttc%Nkpt.gt.1)then
      call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,U_R,U_K)
   else
      U_K(:,:,1) = U_R(:,:,wig0)
   endif
   deallocate(U_R,Rorder,Rsorted)
   !
   call cpu_time(finish)
   write(*,"(A,F)") "     Unn(R) --> Unn(K) cpu timing:", finish-start
   !
   if(reg(structure).ne."None")then
      call interpolateHk2Path(Lttc,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),filename="Uk",data_in=U_K)
   endif
   !
   !fill in the output
   do ib1=1,Nbp
      do ib2=1,Nbp
         !
         iorb = PhysicalUelements%Full_Map(ib1,ib2,1)
         jorb = PhysicalUelements%Full_Map(ib1,ib2,2)
         !
         if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(Uaa(iorb),0d0)
         if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(Uab(iorb,jorb),0d0)
         if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(J(iorb,jorb),0d0)
         if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats_ptr%screened_local(ib1,ib2,1) = dcmplx(J(iorb,jorb),0d0)
         !
         if(allocated(Umats_ptr%bare_local))then
            if(PhysicalUelements%Full_Uaa(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uaa(iorb),0d0)
            if(PhysicalUelements%Full_Uab(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(Uab(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jsf(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
            if(PhysicalUelements%Full_Jph(ib1,ib2)) Umats_ptr%bare_local(ib1,ib2) = dcmplx(J(iorb,jorb),0d0)
         endif
         !
      enddo
   enddo
   !
   if(.not.LocalOnly_)then
      !
      Umats_ptr%screened(:,:,1,:) = U_K
      !
      if(allocated(Umats_ptr%bare))then
         Umats_ptr%bare = U_K
      endif
      !
   endif
   deallocate(U_K)
   !
   if(RealSpace)then
      Nsite = Lttc%Nsite
      if(Hetero%status) Nsite = Hetero%Explicit(2)-Hetero%Explicit(1)+1
      call Expand2Nsite(Umats,Umats_ptr,Nsite)
      call DeallocateBosonicField(Umats_imp)
   endif
   nullify(Umats_ptr)
   !
   call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats_nosum.DAT")
   call BosonicKsum(Umats)
   call dump_BosonicField(Umats,reg(pathINPUTtr),"Uloc_mats.DAT")
   !
end subroutine build_Uret_multiParam_Vn
