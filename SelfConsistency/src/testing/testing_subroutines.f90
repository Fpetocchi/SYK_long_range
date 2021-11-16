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
