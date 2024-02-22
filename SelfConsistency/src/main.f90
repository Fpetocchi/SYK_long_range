program Syk
   !
   use fourier_transforms
   use omp_lib
   implicit none
   !
   !parameters
   real(8),parameter                        :: eps=1e-12
   real(8),parameter                        :: pi=3.14159265358979323846d0
   complex(8),parameter                     :: czero=dcmplx(0.d0,0.d0)
   complex(8),parameter                     :: img=dcmplx(0.d0,1.d0)
   !
   !input file variables
   character(len=20)                        :: InputFile="input.in"
   integer                                  :: Ntau,NR,Nk,NE,Qpower,Nloops
   real(8)                                  :: Beta,alpha,wmatsMax
   real(8)                                  :: Emin,Emax,error_thr
   logical                                  :: logtau,verbose
   !
   !generic variables
   integer                                  :: Nthread,TimeStart,unit
   character(len=1024)                      :: path,filename,alphapad
   character(len=1)                         :: gnupad
   logical                                  :: filexists
   !
   !DoS variables
   integer                                  :: iE,iR,ik
   integer                                  :: NR_,Nk_,NE_
   real(8)                                  :: L_alpha,tk,Norm
   real(8)                                  :: alpha_,Emin_,Emax_
   real(8),allocatable                      :: Egrid(:),DoS(:)
   !
   !Bands variables
   integer                                  :: Nkpath
   real(8),allocatable                      :: dispersion(:)
   !
   !Fields variables
   integer                                  :: in,Nmats
   real(8),allocatable                      :: wmats(:)
   complex(8),allocatable                   :: Gmats(:),Smats(:)
   !
   !self-consistency variables
   integer                                  :: iloop
   character(len=1024)                      :: iloopad
   real(8)                                  :: error
   logical                                  :: converged
   complex(8),allocatable                   :: Gmats_old(:)
   !
   !Free energy variables
   complex(8),allocatable                   :: KE(:)
   complex(8)                               :: K,U,Etot
   !
   !
   !
   !---------------------------------------------------------------------------!
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   !
   Nthread = omp_get_max_threads()
   write(*,"(A,1I4)") new_line("A")//"Setting Nthread:",Nthread
   !
   path = reg("./"//InputFile)
   !
   call inquireFile(path,filexists,verb=.true.)
   if(filexists)then
      !
      write(*,"(A)") "Reading InputFile: "//reg(path)//new_line("A")
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) !"VERBOSE"
      read(unit,*) verbose
      read(unit,*) !"BETA  MAXWMATS  NTAU"
      read(unit,*) Beta, wmatsMax, Ntau
      read(unit,*) !"ALPHA    NR     NK"
      read(unit,*) alpha, NR, Nk
      read(unit,*) !"EMIN     EMAX   NE"
      read(unit,*) Emin, Emax, NE
      read(unit,*) !"QPOWER   LOGTAU"
      read(unit,*) Qpower, logtau
      read(unit,*) !"NLOOPS  ERROR_THR"
      read(unit,*) Nloops, error_thr
      close(unit)
      !
   else
      !
      stop "unable to find input file"
      !
   endif
   !
   if(logtau)then
      if(mod(Ntau,2).eq.0)Ntau=Ntau+1
      if(mod(Ntau-1,4).ne.0)Ntau=Ntau+mod(Ntau-1,4)
   endif
   if(mod(Nk,2).eq.0)then
      write(*,"(A)") "Adding 1 to Nk in order to make it odd"
      Nk = Nk + 1
   endif
   Nmats = int(Beta*wmatsMax/(2d0*pi))
   write(*,"(A,I)") "Number of Matsubara frequencies:",Nmats
   !
   !
   !
   !
   !---------------------------------------------------------------------------!
   !                               BUILD/READ DOS                              !
   !---------------------------------------------------------------------------!
   !
   allocate(Egrid(NE));Egrid=0d0
   allocate(DoS(NE));DoS=0d0
   !
   !DoS
   write(alphapad,"(1F8.2)") alpha
   write(filename,"(1A30)") "DoS_alpha"//reg(alphapad)//".DAT"
   call inquireFile(reg(filename),filexists,verb=.true.,hardstop=.false.)
   if(filexists)then
      !
      write(*,"(A)") "Reading DoS from: "//reg(filename)
      unit = free_unit()
      open(unit,file=reg(filename),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) !"-------------------"
      read(unit,*) !"ALPHA    NR     NK"
      read(unit,*) gnupad, alpha_, NR_, Nk_
      read(unit,*) !"EMIN     EMAX   NE"
      read(unit,*) gnupad, Emin_, Emax_, NE_
      read(unit,*) !"-------------------"
      if(alpha_.ne.alpha) stop "DoS from file has the wrong alpha"
      if(NR_.ne.NR) stop "DoS from file has the wrong NR"
      if(Nk_.ne.Nk) stop "DoS from file has the wrong Nk"
      if(Emin_.ne.Emin) stop "DoS from file has the wrong Emin"
      if(Emax_.ne.Emax) stop "DoS from file has the wrong Emax"
      if(NE_.ne.NE) stop "DoS from file has the wrong NE"
      !
      do iE=1,NE
         read(unit,"(2E20.12)") Egrid(iE),DoS(iE)
      enddo
      close(unit)
      !
   else
      !
      write(*,"(A)") "Building DoS from scratch"
      call tick(TimeStart)
      !
      !Computing the DoS energy grid
      Egrid = linspace(Emin,Emax,NE,istart=.true.,iend=.true.)
      !
      !Computing the normalization factor
      L_alpha=0d0
      do iR=1,NR
         L_alpha = L_alpha + 1 / ( iR**alpha)
      enddo
      write(*,"(A,F)") "Dispersion normalization factor:",L_alpha
      !
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(ik,iR,tk,iE)
      !$OMP DO
      do ik=-floor(Nk/2d0),+floor(Nk/2d0)
         !
         tk = 0d0
         do iR=1,NR
            tk =  tk + cos( (2*pi*ik/Nk) * iR ) / ( iR**alpha )
         enddo
         tk = tk/L_alpha
         !
         iE = minloc(abs(  Egrid - ( 1d0 - tk ) ),dim=1)
         DoS(iE) =  DoS(iE) + 1
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      !Normalization
      Norm = trapezoid_integration(DoS,Egrid)
      write(*,"(A,F)") "DoS normalization factor:,",Norm
      DoS = DoS/Norm
      Norm = trapezoid_integration(DoS,Egrid)
      write(*,"(A,F)") "DoS normalization factor:,",Norm
      !
      !Write to file
      unit = free_unit()
      open(unit,file=reg(filename),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,*) "#-------------------"
      write(unit,*) "# ALPHA    NR     NK"
      write(unit,"(A2,1F10.5,2I10)") " #",alpha, NR, Nk
      write(unit,*) "#EMIN     EMAX   NE   LOGRID"
      write(unit,"(A2,2F10.5,1I10,1L)") " #",Emin, Emax, NE
      write(unit,*) "#-------------------"
      do iE=1,NE
         write(unit,"(2E20.12)") Egrid(iE),DoS(iE)
      enddo
      close(unit)
      !
      write(*,"(A,F)") new_line("A")//"DoS construction finished. Total timing (s): ",tock(TimeStart)
      !
   endif
   !
   !Bands, this is just for fun
   if(verbose)then
      !
      !Computing the normalization factor
      L_alpha=0d0
      do iR=1,NR
         L_alpha = L_alpha + 1 / ( iR**alpha)
      enddo
      write(*,"(A,F)") "Dispersion normalization factor:",L_alpha
      !
      Nkpath=201
      allocate(dispersion(Nkpath));dispersion=0d0
      !
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(ik,iR,tk,iE)
      !$OMP DO
      do ik=-floor(Nkpath/2d0),+floor(Nkpath/2d0)
         !
         tk = 0d0
         do iR=1,NR
            tk =  tk + cos( (2*pi*ik/Nkpath) * iR ) / ( iR**alpha )
         enddo
         dispersion(ik+1+floor(Nkpath/2d0)) = 1d0 - tk/L_alpha
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      write(filename,"(1A30)") "Ek_alpha"//reg(alphapad)//".DAT"
      unit = free_unit()
      open(unit,file=reg(filename),form="formatted",status="unknown",position="rewind",action="write")
      do ik=-floor(Nkpath/2d0),+floor(Nkpath/2d0)
         write(unit,"(2E20.12)") 2*pi*ik/Nkpath,dispersion(ik+1+floor(Nkpath/2d0))
      enddo
      close(unit)
      !
   endif
   !
   !
   !
   !
   !---------------------------------------------------------------------------!
   !                             INITIALIZE FIELDS                             !
   !---------------------------------------------------------------------------!
   !
   allocate(wmats(Nmats));wmats=0d0
   wmats = FermionicFreqMesh(Beta,Nmats)
   !
   !Non interacting Green's function
   allocate(Gmats(Nmats));Gmats=czero
   call calcGmats()
   call dumpField(Gmats,"./G",pad="it0")
   !
   !Zeroth iteration self-energy
   allocate(Smats(Nmats));Smats=czero
   call calcSmats()
   call dumpField(Smats,"./S",pad="it0")
   !
   !
   !
   !
   !---------------------------------------------------------------------------!
   !                             SELF-CONSISTENCY LOOP                         !
   !---------------------------------------------------------------------------!
   !
   allocate(Gmats_old(Nmats));Gmats_old=czero
   converged = .false.
   SCloops: do iloop=1,Nloops
      !
      write(*,"(A,1I5,A)")new_line("A")//new_line("A")//"---- Loop #",iloop," ----"
      write(iloopad,"(1I100)") iloop
      !
      !Store old Green's function for convergence check
      Gmats_old = Gmats
      !
      call calcGmats(Sigma=Smats)
      if(verbose) call dumpField(Gmats,"./G",pad="it"//reg(iloopad))
      !
      error = check_error(Gmats,Gmats_old)
      if(error.gt.error_thr)then
         write(*,"(2(A,1E10.3))")"Error: ",error," > ",error_thr
      else
         write(*,"(2(A,1E10.3),A)")"Error: ",error," < ",error_thr," Converged!"
         converged = .true.
         exit SCloops
      endif
      if(iloop.gt.Nloops)write(*,"(A)")"WARNING: self-consistency cylce not converged, increase NLOOP."
      !
      !Compute self-energy for next iteration
      call calcSmats()
      if(verbose) call dumpField(Smats,"./S",pad="it"//reg(iloopad))
      !
   enddo SCloops
   deallocate(Gmats_old)
   !
   if(converged)then
      call dumpField(Gmats,"./G",pad="converged")
      call dumpField(Smats,"./S",pad="converged")
   endif
   !
   !
   !
   !
   !---------------------------------------------------------------------------!
   !                                FREE ENERGY                                !
   !---------------------------------------------------------------------------!
   !
   !Kinetic term
   allocate(KE(NE));KE=czero
   !$OMP PARALLEL DEFAULT(SHARED),&
   !$OMP PRIVATE(iE,in)
   !$OMP DO
   do iE=1,NE
      do in=1,Nmats
         KE(iE) = KE(iE) - zlog( -img*wmats(in) + Egrid(iE) + Smats(in) ) * DoS(iE) / Beta
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
   K = czero
   K = trapezoid_integration(KE,Egrid)
   deallocate(KE)
   !
   !Potential term
   U=czero
   do in=1,Nmats
      U = U - ( (2*Qpower-1d0)/(2*Qpower) ) * Smats(in)*Gmats(in) / Beta
   enddo
   !
   !Total energy
   Etot = K + U
   !
   !Report values
   write(*,"(A,2E20.12)")"Ekin: ", dreal(K), dimag(K)
   write(*,"(A,2E20.12)")"Epot: ", dreal(U), dimag(U)
   write(*,"(A,2E20.12)")"Epot: ", dreal(Etot), dimag(Etot)
   !
   !
   !
   contains
   !
   !
   !
   subroutine calcGmats(Sigma)
      !
      implicit none
      !
      complex(8),intent(in),optional        :: Sigma(:)
      integer                               :: iw
      complex(8),allocatable                :: GwE(:)
      complex(8),allocatable                :: zeta(:)
      !
      Gmats=czero
      !
      allocate(GwE(NE));GwE=czero
      allocate(zeta(Nmats));zeta=czero
      !
      zeta = img*wmats
      if(present(Sigma)) zeta =  zeta - Sigma
      !
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iw,iE,GwE)
      !$OMP DO
      do iw=1,Nmats
         !
         GwE = czero
         do iE=1,NE
            GwE(iE) = DoS(iE) / ( zeta(iw) - Egrid(iE) )
         enddo
         !
         Gmats(iw) = trapezoid_integration(GwE,Egrid)
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(GwE,zeta)
      !
   end subroutine calcGmats
   !
   subroutine calcSmats()
      !
      implicit none
      !
      integer                               :: itau
      real(8)                               :: tau2
      real(8),allocatable                   :: tau(:)
      complex(8),allocatable                :: Gitau(:),Sitau(:)
      !
      Smats=czero
      !
      allocate(tau(Ntau));tau=0d0
      if(logtau)then
         tau = denspace(beta,Ntau)
      else
         tau = linspace(0d0,Beta,Ntau)
      endif
      !
      !Compute G(tau)
      call tick(TimeStart)
      allocate(Gitau(Ntau));Gitau=czero
      call Fmats2itau(Beta,Gmats,Gitau,asympt_corr=.true.,tau_uniform=logtau)
      if(verbose)write(*,"(A,F)") new_line("A")//"G(iw) --> G(tau). Total timing (s): ",tock(TimeStart)
      !
      !Compute S(tau)
      call tick(TimeStart)
      allocate(Sitau(Ntau));Sitau=czero
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(itau,tau2)
      !$OMP DO
      do itau=1,Ntau
         !
         tau2 = tau(Ntau)-tau(itau)
         if (dabs(tau2-tau(Ntau-itau+1)).gt.eps) stop "calcSmats: itau2 not found."
         !
         !Note that G(-tau) = -G(beta-tau)
         Sitau(itau) = (-1)**(Qpower+1) * (Gitau(itau)**Qpower) * (-Gitau(Ntau-itau+1)**(Qpower-1))
         !
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Gitau)
      if(verbose)write(*,"(A,F)") "S(tau) calculation. Total timing (s): ",tock(TimeStart)
      !
      !Compute S(iw)
      call tick(TimeStart)
      call Fitau2mats(Beta,Sitau,Smats,tau_uniform=logtau)
      if(verbose)write(*,"(A,F)") "G(iw) --> G(tau). Total timing (s): ",tock(TimeStart)
      deallocate(tau,Sitau)
      !
   end subroutine calcSmats
   !
   subroutine dumpField(Field,header,pad)
      !
      implicit none
      !
      complex(8),intent(in)                 :: Field(:)
      character(len=*),intent(in)           :: header
      character(len=*),intent(in),optional  :: pad
      integer                               :: iw
      character(len=1024)                   :: fname,betapad
      !
      write(betapad,"(1F8.2)") beta
      write(fname,"(1A100)") reg(header)//"mats_alpha"//reg(alphapad)//"_beta"//reg(betapad)
      if(present(pad)) fname = reg(fname)//"_"//reg(pad)
      fname = reg(fname)//".DAT"
      !
      unit = free_unit()
      open(unit,file=reg(fname),form="formatted",status="unknown",position="rewind",action="write")
      do iw=1,Nmats
         write(unit,"(3E20.12)") wmats(iw),dreal(Field(iw)),dimag(Field(iw))
      enddo
      close(unit)
      !
   end subroutine dumpField
   !
   function check_error(fnew,fold) result(err)
      implicit none
      complex(8),intent(in)                 :: fnew(:)
      complex(8),intent(in)                 :: fold(:)
      real(8)                               :: err
      real(8)                               :: S,M
      integer                               :: i
      !
      S=0.d0 ; M=0.d0
      do i=1,size(fnew)
         M = M + abs(fnew(i)-fold(i))
         S = S + abs(fnew(i))
      enddo
      !
      err = M / S
      !
   end function check_error
   !
   !
   !
end program Syk
