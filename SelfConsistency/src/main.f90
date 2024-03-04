program Syk
   !
   use fourier_transforms
   use omp_lib
   implicit none
   !
   !parameters
   real(8),parameter                        :: eps=1e-12
   real(8),parameter                        :: pi=3.14159265358979323846d0
   real(8),parameter                        :: K2eV= 1d0 ! 8.617333262d-5 Temperature directly in eV
   complex(8),parameter                     :: czero=dcmplx(0.d0,0.d0)
   complex(8),parameter                     :: img=dcmplx(0.d0,1.d0)
   !
   !input file variables
   character(len=20)                        :: InputFile="input.in"
   integer                                  :: Ntau_in,NR,Nk,NE
   integer                                  :: Qpower,Nloops
   integer                                  :: NT,NT_inf,NT_intp
   integer                                  :: min_Nmats,max_Nmats
   real(8)                                  :: alpha,wmatsMax,mu,Wband
   real(8)                                  :: Emin,Emax,error_thr
   real(8)                                  :: thop,Uloc
   real(8)                                  :: Tmin,Tmax
   logical                                  :: logtau,verbose
   !
   !generic variables
   integer                                  :: Nthread,TimeStart,unit
   character(len=1024)                      :: path,filename
   character(len=1024)                      :: alphapad,betapad
   character(len=1024)                      :: iloopad,Tpad
   character(len=1)                         :: gnupad
   logical                                  :: filexists
   !
   !DoS variables
   integer                                  :: iE,iR,ik
   integer                                  :: NR_,Nk_,NE_
   real(8)                                  :: L_alpha,tk,Norm
   real(8)                                  :: alpha_,Emin_,Emax_,thop_
   real(8),allocatable                      :: Egrid(:),DoS(:)
   !
   !Bands variables
   integer                                  :: Nkpath
   real(8),allocatable                      :: dispersion(:)
   !
   !Temperature variables
   integer                                  :: iT
   real(8)                                  :: T,dT,Beta
   real(8),allocatable                      :: Ts(:),Ts_intp(:)
   !
   !Fields variables
   integer                                  :: in,Nmats,Ntau,itau
   real(8),allocatable                      :: wmats(:),tau(:)
   complex(8),allocatable                   :: Gmats(:),Smats(:)
   complex(8),allocatable                   :: Gitau(:),Sitau(:)
   !
   !self-consistency variables
   integer                                  :: iloop
   real(8)                                  :: error,mixing
   logical                                  :: converged,Smatsexists,Gmatsexists
   complex(8),allocatable                   :: Gmats_old(:),Smats_old(:)
   real(8),allocatable                      :: ReSmats(:),ImSmats(:)
   real(8),allocatable                      :: wmats_old(:)
   !
   !Free energy variables
   integer                                  :: fact
   real(8)                                  :: Cv_intp_last
   complex(8)                               :: logarg,K,U
   complex(8),allocatable                   :: KE(:),Ut(:),Energy(:,:)
   real(8),allocatable                      :: Energy_intp(:,:)
   real(8),allocatable                      :: Cv_intp(:)
   real(8),allocatable                      :: Ds_intp(:)
   !
   !
   !
   !
   !###########################################################################!
   !#    READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE     #!
   !###########################################################################!
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
      read(unit,*) !#VERBOSE
      read(unit,*) verbose
      read(unit,*) !# ALPHA     NR       NK      WBAND
      read(unit,*) alpha, NR, Nk, Wband
      read(unit,*) !# EMIN      EMAX     NE
      read(unit,*) Emin, Emax, NE
      read(unit,*) !# QPOWER    LOGTAU   MU
      read(unit,*) Qpower, logtau, mu
      read(unit,*) !# THOP      ULOC
      read(unit,*) thop, Uloc
      read(unit,*) !# MAXWMATS  NTAU     MIN_NMATS   MAX_NAMTS
      read(unit,*) wmatsMax, Ntau_in, min_Nmats, max_Nmats
      read(unit,*) !# NLOOPS    ERROR    MIXING
      read(unit,*) Nloops, error_thr, mixing
      read(unit,*) !# TMIN      TMAX     NT     NT_inf    NT_intp
      read(unit,*) Tmin, Tmax, NT, NT_inf, NT_intp
      close(unit)
      !
   else
      !
      stop "unable to find input file"
      !
   endif
   !
   if(mod(Nk,2).eq.0)then
      write(*,"(A)") "Adding 1 to Nk in order to make it odd"
      Nk = Nk + 1
   endif
   !
   !
   !
   !
   !###########################################################################!
   !#                             BUILD/READ DOS                              #!
   !###########################################################################!
   !
   allocate(Egrid(NE));Egrid=0d0
   allocate(DoS(NE));DoS=0d0
   !
   !DoS
   write(alphapad,"(1F8.2)") alpha
   write(filename,"(1A30)") "DoS_alpha"//reg(alphapad)//".DAT"
   call inquireFile(reg(filename),filexists,verb=verbose,hardstop=.false.)
   if(filexists)then
      !
      write(*,"(A)") "Reading DoS from: "//reg(filename)
      unit = free_unit()
      open(unit,file=reg(filename),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) !"#-------------------------------------------#"
      read(unit,*) !"ALPHA    NR     NK"
      read(unit,*) gnupad, alpha_, NR_, Nk_
      read(unit,*) !"EMIN     EMAX   NE    THOP"
      read(unit,*) gnupad, Emin_, Emax_, NE_, thop_
      read(unit,*) !"#-------------------------------------------#"
      if(alpha_.ne.alpha) stop "DoS from file has the wrong alpha"
      if(NR_.ne.NR) stop "DoS from file has the wrong NR"
      if(Nk_.ne.Nk) stop "DoS from file has the wrong Nk"
      if(Emin_.ne.Emin) stop "DoS from file has the wrong Emin"
      if(Emax_.ne.Emax) stop "DoS from file has the wrong Emax"
      if(NE_.ne.NE) stop "DoS from file has the wrong NE"
      if(thop_.ne.thop) stop "DoS from file has the wrong thop"
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
      write(unit,*) "#------------------------------------------------#"
      write(unit,*) "# ALPHA    NR     NK"
      write(unit,"(A2,1F10.5,2I10)") " #",alpha, NR, Nk
      write(unit,*) "# EMIN     EMAX   NE    THOP"
      write(unit,"(A2,2F10.5,1I10,1F10.5)") " #",Emin, Emax, NE, thop
      write(unit,*) "#------------------------------------------------#"
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
         tk =  tk + thop*cos( (2*pi*ik/Nkpath) * iR ) / ( iR**alpha )
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
   !
   !
   !
   !###########################################################################!
   !#                             TEMPERATURE LOOP                            #!
   !###########################################################################!
   !
   allocate(Ts(NT+NT_inf));Ts=0d0
   allocate(Energy(3,NT+NT_inf));Energy=czero
   do iT=1,NT+NT_inf
      !
      !adaptive T increment
      dT=0d0
      if(NT.gt.1)then
         if(iT.le.NT)then
            dT = (iT-1)*abs(Tmax-Tmin)/dble(NT-1)
         else
            dT = abs(Tmax-Tmin) + (iT-NT)*abs(Wband/K2eV-Tmax)/dble(NT_inf)
         endif
      endif
      !
      !adaptive Matsubara mesh
      T = Tmin + dT
      Beta = 1d0 / (T*K2eV)
      Nmats = int(Beta*wmatsMax/(2d0*pi))
      if(Nmats.gt.max_Nmats)then
         write(*,"(2(A,I),A,1F12.4,A)") "Nmats: ",Nmats,"   cutoffed to: ",max_Nmats, "   corresponding to: ",max_Nmats*pi/beta," eV "
         Nmats = max_Nmats
      endif
      if(Nmats.lt.min_Nmats)then
         write(*,"(2(A,I),A,1F12.4,A)") "Nmats: ",Nmats,"   increased to: ",max_Nmats, "   corresponding to: ",min_Nmats*pi/beta," eV "
         Nmats = min_Nmats
      endif
      Ts(iT) = T
      !
      !adaptive tau mesh
      Ntau = Ntau_in
      if(Ntau_in.eq.0) Ntau = int(Nmats)/2
      if(logtau)then
         if(mod(Ntau,2).eq.0)Ntau=Ntau+1
         if(mod(Ntau-1,4).ne.0)Ntau=Ntau+mod(Ntau-1,4)
      endif
      !
      write(Tpad,"(1F30.8)") T
      write(betapad,"(1F8.2)") beta
      write(path,"(1A30)") "./loops_T"//reg(Tpad)//"/"
      call createDir(reg(path),verb=verbose)
      !
      write(*,"(A)") new_line("A")//new_line("A")//".................................................."//new_line("A")
      if(K2eV.eq.1d0)then
         write(*,"(1A6,1F12.5,2(A16,1F12.5),A16,I)") "T(K): ",T/8.617333262d-5,"T(eV): ",1d0/Beta,"Beta(1/eV): ",Beta,"Nmats: ",Nmats
      else
         write(*,"(1A6,1F12.5,2(A16,1F12.5),A16,I)") "T(K): ",T,"T(eV): ",1d0/Beta,"Beta(1/eV): ",Beta,"Nmats: ",Nmats
      endif
      write(*,"(A)")"Data sored in: "//reg(path)
      !
      !------------------------------------------------------------------------!
      !                            INITIALIZE FIELDS                           !
      !------------------------------------------------------------------------!
      !
      if(allocated(wmats))deallocate(wmats)
      allocate(wmats(Nmats));wmats=0d0
      wmats = FermionicFreqMesh(Beta,Nmats)
      !
      !Non interacting Green's function
      if(allocated(Gmats))deallocate(Gmats)
      allocate(Gmats(Nmats));Gmats=czero
      call calcGmats()
      call dumpField(Gmats,reg(path)//"G",pad="it0")
      !
      !Zeroth iteration self-energy
      if(allocated(Smats))deallocate(Smats)
      allocate(Smats(Nmats));Smats=czero
      if(iT.eq.1)then
         !look for user-provided initial guess
         write(filename,"(1A100)") "Smats_alpha"//reg(alphapad)//"_beta"//reg(betapad)//"_init.DAT"
         call inquireFile(reg(path)//"../"//reg(filename),filexists,verb=verbose,hardstop=.false.)
         if(filexists)then
            write(*,"(A)")"Reading starting Sigma."
            call readField(Smats,reg(path)//"../S",pad="init")
         else
            write(*,"(A)")"Initializing Sigma from bare G."
            call calcSmats()
         endif
         call dumpField(Smats,reg(path)//"S",pad="it0")
      else
         allocate(ReSmats(Nmats));ReSmats=0d0
         allocate(ImSmats(Nmats));ImSmats=0d0
         ReSmats = cubic_interp( wmats_old, dreal(Smats_old), wmats )
         ImSmats = cubic_interp( wmats_old, dimag(Smats_old), wmats )
         Smats = dcmplx(ReSmats,ImSmats)
         deallocate(ReSmats,ImSmats)
      endif
      !
      !------------------------------------------------------------------------!
      !                          SELF-CONSISTENCY LOOP                         !
      !------------------------------------------------------------------------!
      !
      if(allocated(Smats_old))deallocate(Smats_old)
      allocate(Smats_old(Nmats));Smats_old=czero
      allocate(Gmats_old(Nmats));Gmats_old=czero
      converged = .false.
      SCloops: do iloop=1,Nloops
         !
         write(*,"(A,1I5,A)")new_line("A")//"---- Loop #",iloop," ----"
         write(iloopad,"(1I1000)") iloop
         !
         !Check if converged fields are already present
         write(filename,"(1A100)") "Smats_alpha"//reg(alphapad)//"_beta"//reg(betapad)//"_converged.DAT"
         call inquireFile(reg(path)//reg(filename),Smatsexists,verb=verbose,hardstop=.false.)
         write(filename,"(1A100)") "Gmats_alpha"//reg(alphapad)//"_beta"//reg(betapad)//"_converged.DAT"
         call inquireFile(reg(path)//reg(filename),Gmatsexists,verb=verbose,hardstop=.false.)
         if(Smatsexists.and.Gmatsexists)then
            !
            write(filename,"(1A100)") "Smats_alpha"//reg(alphapad)//"_beta"//reg(betapad)//"_converged.DAT"
            write(*,"(A)")"Reading converged Smats from: "//reg(path)//reg(filename)
            call readField(Smats,reg(path)//"S",pad="converged")
            Smats_old = Smats
            !
            write(filename,"(1A100)") "Gmats_alpha"//reg(alphapad)//"_beta"//reg(betapad)//"_converged.DAT"
            write(*,"(A)")"Reading converged Gmats from: "//reg(path)//reg(filename)
            call readField(Gmats,reg(path)//"S",pad="converged")
            Gmats_old = Gmats
            !
            write(*,"(A)")"Skipping self-consistency."
            exit SCloops
            !
         endif
         !
         !Store old Green's function for convergence check
         Gmats_old = Gmats
         !
         call calcGmats(Sigma=Smats)
         if(verbose) call dumpField(Gmats,reg(path)//"G",pad="it"//reg(iloopad))
         !
         !Store old self-energy for mixing
         Smats_old = Smats
         !
         error = check_error(Gmats,Gmats_old) !the first error is always with respect to G0
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
         Smats = (1d0-mixing)*Smats + mixing*Smats_old
         if(verbose) call dumpField(Smats,reg(path)//"/S",pad="it"//reg(iloopad))
         !
      enddo SCloops
      deallocate(Gmats_old)
      !
      if(converged)then
         call dumpField(Gmats,reg(path)//"G",pad="converged")
         call dumpField(Smats,reg(path)//"S",pad="converged")
      endif
      !
      !Store for nex Temperature interpolation
      if(allocated(wmats_old))deallocate(wmats_old)
      allocate(wmats_old(Nmats));wmats_old=0d0
      wmats_old = wmats
      Smats_old = Smats
      !
      !------------------------------------------------------------------------!
      !                         FREE ENERGY CALCULATION                        !
      !------------------------------------------------------------------------!
      !
      !Kinetic term
      allocate(KE(NE));KE=czero
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(iE,in,logarg)
      !$OMP DO
      do iE=1,NE
         do in=1,Nmats
            !
            !Eq.35
            !logarg = ( -img*wmats(in) + Egrid(iE) + Smats(in) )
            !
            !Eq.37
            logarg = ( img*wmats(in) + mu - Egrid(iE) - Smats(in) )/( img*wmats(in) + mu )
            KE(iE) = KE(iE) - zlog( logarg ) * DoS(iE) / Beta
            !
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      K = czero
      K = trapezoid_integration(KE,Egrid) - log( 1d0 + exp(Beta*mu) )/Beta
      deallocate(KE)
      !
      !Potential term
      !
      !Eq.35
      !U=czero
      !do in=1,Nmats
      !   U = U - ( (2*Qpower-1d0)/(2*Qpower) ) * Smats(in)*Gmats(in) / Beta
      !enddo
      !
      !Eq.37
      if(allocated(Gitau))deallocate(Gitau)
      allocate(Gitau(Ntau));Gitau=czero
      call Fmats2itau(Beta,Gmats,Gitau,asympt_corr=.true.,tau_uniform=logtau)
      if(allocated(Sitau))deallocate(Sitau)
      allocate(Sitau(Ntau));Sitau=czero
      call Fmats2itau(Beta,Smats,Sitau,asympt_corr=.true.,tau_uniform=logtau)
      !
      allocate(Ut(Ntau));Ut=czero
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(itau)
      !$OMP DO
      do itau=1,Ntau
         Ut(itau) = Ut(itau) - (Gitau(Ntau-itau+1)**Qpower)*(Gitau(itau)**Qpower) + Sitau(itau)*Gitau(Ntau-itau+1)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(Gitau,Sitau)
      !
      if(allocated(tau))deallocate(tau)
      allocate(tau(Ntau));tau=0d0
      if(logtau)then
         tau = denspace(beta,Ntau)
      else
         tau = linspace(0d0,Beta,Ntau)
      endif
      !
      U = czero
      U = trapezoid_integration(Ut,tau) * (Uloc**2)/(2*Qpower)
      deallocate(Ut,tau)
      !
      !Store energy components for a given temperature
      Energy(1,iT) = 2*K   !the 2 comes from the fact that the intergral is only for iw>0
      Energy(2,iT) = U     !this is a tau convolution so no factors
      Energy(3,iT) = K + U
      !
   enddo ! end of Tloop
   !
   !
   !
   !
   !###########################################################################!
   !#                              SPECIFIC HEAT                              #!
   !###########################################################################!
   !
   !Cubic interpolation for better derivative
   allocate(Ts_intp(NT_intp));Ts_intp=0d0
   allocate(Energy_intp(3,NT_intp));Energy_intp=0d0
   !
   if(NT_inf.eq.0)then
      Ts_intp = linspace(Tmin,Tmax,NT_intp,istart=.true.,iend=.true.)
   else
      Ts_intp = linspace(Tmin,Wband/K2eV,NT_intp,istart=.true.,iend=.true.)
   endif
   !interpolate energy components (Real parts only)
   Energy_intp(1,:) = cubic_interp( Ts, dreal(Energy(1,:)), Ts_intp )
   Energy_intp(2,:) = cubic_interp( Ts, dreal(Energy(2,:)), Ts_intp )
   Energy_intp(3,:) = cubic_interp( Ts, dreal(Energy(3,:)), Ts_intp )
   !
   !Specific Heat obtained by central difference
   allocate(Cv_intp(NT_intp));Cv_intp=0d0
   do iT=2,NT_intp-1
      Cv_intp(iT) = ( Energy_intp(3,iT+1) - Energy_intp(3,iT-1) ) / ( Ts_intp(iT+1) - Ts_intp(iT-1) )
   enddo
   Cv_intp(1) = cubic_interp( Ts_intp(2:NT_intp-1), Cv_intp(2:NT_intp-1), Tmin )
   Cv_intp_last = Tmax
   if(NT_inf.gt.0)Cv_intp_last = Wband/K2eV
   Cv_intp(NT_intp) = cubic_interp( Ts_intp(2:NT_intp-1), Cv_intp(2:NT_intp-1), Cv_intp_last )
   !
   !Entropy difference: S_f(T) - S_i(T) = int ^f _i Cv/T dT
   allocate(Ds_intp(NT_intp));Ds_intp=0d0
   do iT=NT_intp-1,1,-1
      Ds_intp(iT) = trapezoid_integration( Cv_intp(iT:NT_intp)/Ts_intp(iT:NT_intp) ,Ts_intp(iT:NT_intp))
   enddo
   !
   !Print results
   write(filename,"(1A30)") "DeltaS_alpha"//reg(alphapad)//".DAT"
   unit = free_unit()
   open(unit,file=reg(filename),form="formatted",status="unknown",position="rewind",action="write")
   write(unit,"(6A20)") "# Temp"," Energy"," Specific Heat"," Entropy Diff", " Re Ek", " Re Ep"
   do iT=1,NT_intp
      fact=0
      write(unit,"(6E20.12)") Ts_intp(iT)       , &
                              Energy_intp(3,iT) , &
                              Cv_intp(iT)       , &
                              Ds_intp(iT)       , &
                              Energy_intp(1,iT) , &
                              Energy_intp(2,iT)
   enddo
   close(unit)
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
      real(8)                               :: tau2
      !
      Smats=czero
      !
      if(allocated(tau))deallocate(tau)
      allocate(tau(Ntau));tau=0d0
      if(logtau)then
         tau = denspace(beta,Ntau)
      else
         tau = linspace(0d0,Beta,Ntau)
      endif
      !
      !Compute G(tau)
      call tick(TimeStart)
      if(allocated(Gitau))deallocate(Gitau)
      allocate(Gitau(Ntau));Gitau=czero
      call Fmats2itau(Beta,Gmats,Gitau,asympt_corr=.true.,tau_uniform=logtau)
      if(verbose)write(*,"(A,F)") new_line("A")//"G(iw) --> G(tau). Total timing (s): ",tock(TimeStart)
      !
      !Compute S(tau)
      call tick(TimeStart)
      if(allocated(Sitau))deallocate(Sitau)
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
         Sitau(itau) = (-1)**(Qpower+1) * (Uloc**2) * (Gitau(itau)**Qpower) * ((-Gitau(Ntau-itau+1))**(Qpower-1))
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
      if(verbose)write(*,"(A,F)") "S(tau) --> S(iw). Total timing (s): ",tock(TimeStart)
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
      character(len=1024)                   :: fname,bpad
      !
      write(bpad,"(1F8.2)") beta
      write(fname,"(1A100)") reg(header)//"mats_alpha"//reg(alphapad)//"_beta"//reg(bpad)
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
   subroutine readField(Field,header,pad)
      !
      implicit none
      !
      complex(8),intent(inout)              :: Field(:)
      character(len=*),intent(in)           :: header
      character(len=*),intent(in),optional  :: pad
      real(8)                               :: wmats_read,ReF,ImF
      integer                               :: iw
      character(len=1024)                   :: fname,bpad
      !
      write(bpad,"(1F8.2)") beta
      write(fname,"(1A100)") reg(header)//"mats_alpha"//reg(alphapad)//"_beta"//reg(bpad)
      if(present(pad)) fname = reg(fname)//"_"//reg(pad)
      fname = reg(fname)//".DAT"
      !
      Field=czero
      !
      unit = free_unit()
      open(unit,file=reg(fname),form="formatted",status="unknown",position="rewind",action="read")
      do iw=1,Nmats
         read(unit,*) wmats_read,ReF,ImF
         Field(iw) = dcmplx(ReF,ImF)
         !if(wmats_read.eq.wmats(iw))then
         !   Field(iw) = dcmplx(ReF,ImF)
         !else
         !   stop "readField: wrong imaginary frequency point."
         !endif
      enddo
      close(unit)
      !
   end subroutine readField
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
