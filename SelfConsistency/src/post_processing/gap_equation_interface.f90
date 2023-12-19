subroutine calc_Tc(pathOUTPUT,Inputs,Lttc,Wlat)
   !
   use parameters
   use linalg, only : zeye
   use utils_misc
   use utils_fields
   use crystal
   use file_io
   use gap_equation
   use input_vars, only : pathINPUT, pathINPUTtr, eta
   implicit none
   !
   character(len=*),intent(in)           :: pathOUTPUT
   type(SCDFT),intent(in)                :: Inputs
   type(Lattice),intent(in)              :: Lttc
   type(BosonicField),intent(in)         :: Wlat
   !
   real(8)                               :: dE,dT
   integer                               :: iE,iE1,iE2,EE_dim,Ngrid
   integer                               :: Norb,Nkpt,unit,Ngrid_read
   complex(8)                            :: Kint
   real(8),allocatable                   :: Egrid_print(:),DoS_Model_interp(:)
   real(8),allocatable                   :: Tlist(:),Delta_T(:)
   real(8),allocatable                   :: Zph(:),Kph(:)
   complex(8),allocatable                :: Kel(:)
   complex(8),allocatable                :: Hk_used(:,:,:)
   character(len=255)                    :: printpath,printpath_T
   !
   integer                               :: iloop,iT,Q_p(2),Q_f(2)
   real(8)                               :: Temp,Beta,Beta_DFT,errDelta
   real(8)                               :: dumE,ReD,ImD,tanh_b,tanh_f
   complex(8)                            :: Kel_b,Kel_f
   real(8),allocatable                   :: EDsq(:),kpt_QP(:,:)
   complex(8),allocatable                :: Delta(:),oldDelta(:),newDelta(:)
   logical                               :: converged,filexists
   !
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Tc"
   !
   !
   !grid meshes and dimension of the UT representation of the Kernels
   Ngrid = Inputs%Ngrid
   EE_dim = Ngrid * (Ngrid+1)/2
   !
   Norb = Lttc%Norb
   Nkpt = Lttc%Nkpt
   !
   if(.not.Wlat%status) stop "calc_Tc: Wlat not properly initialized."
   if(Wlat%Nkpt.ne.Nkpt) stop "calc_Tc: Wlat has different number of k-points with respect to Lttc."
   if(Wlat%Nbp.ne.(Norb**2)) stop "calc_Tc: Wlat has different orbital dimension with respect to Lttc."
   !
   allocate(Hk_used(Norb,Norb,Nkpt));Hk_used=czero
   if(Inputs%HkRenorm)then
      !
      printpath=reg(pathOUTPUT)//"Gap_Equation_Renorm_Hk/"
      if(.not.allocated(Lttc%Hk_qp))then
         write(*,"(A)")"     calc_Tc: QP bandstructure not allocated."
         call inquireFile(reg(pathINPUTtr)//"G0W0plots/Hk_qp_s1.DAT",filexists,hardstop=.true.,verb=.true.)
         if(filexists) call read_Hk(Hk_used,kpt_QP,reg(pathINPUTtr)//"G0W0plots/","Hk_qp_s1.DAT")
         if(.not.all(kpt_QP.eq.Lttc%kpt)) stop "calc_Tc: Lttc k-grid and Hk_qp_s1.DAT one do not coincide."
         deallocate(kpt_QP)
      else
         Hk_used = Lttc%Hk_qp(:,:,:,1)
      endif
      !
   else
      !
      printpath=reg(pathOUTPUT)//"Gap_Equation/"
      Hk_used = Lttc%Hk
      !
   endif
   call createDir(reg(printpath),verb=verbose)
   !
   call Initialize_inputs(reg(pathINPUT),Inputs,Lttc,Hk_used)
   deallocate(Hk_used)
   !
   !recompute the DoS accounting for the imaginary part of the G0W0 self-energy, verbatim of print_G0W0_dispersion
   if(calc_Kel)then
      !
      if(Inputs%G0W0Renorm)then
         call overload_G0W0(reg(pathINPUT),reg(pathINPUTtr),Lttc,eta,Inputs%DoSthresh)
         printpath=reg(pathOUTPUT)//"Gap_Equation_Renorm_G0W0/"
      endif
      !
      if(Inputs%DMFTRenorm)then
         call overload_DMFT(reg(pathOUTPUT)//"K_resolved/",Lttc,Inputs%DoSthresh)
         printpath=reg(pathOUTPUT)//"Gap_Equation_Renorm_DMFT/"
      endif
      !
      call dump_Field_component(DoS_Model,reg(printpath),"DoS_Model.DAT",Egrid_Model)
      !
      allocate(DoS_Model_interp(Ngrid));DoS_Model_interp=0d0
      if(calc_phonons)then
         DoS_Model_interp = cubic_interp( Egrid_Model, DoS_Model, Egrid_Phonons )
         call dump_Field_component(DoS_Model_interp,reg(printpath),"DoS_Model_interp.DAT",Egrid_Phonons)
      else
         DoS_Model_interp = DoS_Model
      endif
      !
   endif
   call dump_Field_component(DoS_DFT,reg(printpath),"DoS_DFT.DAT",Egrid_Phonons)
   !
   !Store inside the module the required energy averages
   if(calc_Kel)call calc_energy_averages(Wlat%screened,Lttc,Wlat%Beta,Inputs%Wk_cutoff,reg(printpath),reg(Inputs%printmode_el))
   !
   if(Inputs%calc_Tc)then
      !
      !initialization of Delta
      allocate(Delta(Ngrid));Delta=czero
      call inquireFile(reg(printpath)//"Delta_init.DAT",filexists,hardstop=.false.)
      if(filexists)then
         !read file
         unit = free_unit()
         open(unit,file=reg(printpath)//"Delta_init.DAT",form="formatted",status="old",position="rewind",action="read")
         read(unit,*)Ngrid_read
         if(Ngrid_read.ne.Ngrid)stop "calc_Tc: the provided Delta_init.DAT haw wrong size."
         do iE=1,Ngrid
            read(unit,*) dumE,ReD,ImD
            Delta(iE) = dcmplx(ReD,ImD)
         enddo
      else
         !gaussian initialization
         do iE=1,Ngrid
            Delta(iE) = exp( -0.5d0*(Egrid_Phonons(iE)/Inputs%DeltaInit_B)**2 )
         enddo
         Delta = Delta * Inputs%DeltaInit_M * eV2DFTgrid / Delta(minloc(abs(Egrid_Phonons),dim=1))
      endif
      !
      where(DoS_Model.lt.Inputs%DoSthresh) Delta=czero
      !
      allocate(oldDelta(Ngrid));oldDelta=Delta
      allocate(Tlist(Inputs%Tsteps));Tlist=0d0
      allocate(Delta_T(Inputs%Tsteps));Delta_T=0d0
      !
      !
      !============================ TEMPERATURE LOOP =============================!
      !
      write(*,"(A)") new_line("A")//"---- Starting temperature scan"
      do iT=1,Inputs%Tsteps
         !
         dT=0d0
         if(Inputs%Tsteps.gt.1) dT = (iT-1)*abs(Inputs%Tbounds(2)-Inputs%Tbounds(1))/dble(Inputs%Tsteps-1)
         !
         Temp = Inputs%Tbounds(1) + dT
         Beta = 1d0 / (Temp*K2eV)
         Beta_DFT = 1d0 / (Temp*K2eV*eV2DFTgrid)
         !
         printpath_T = reg(printpath)//"loops_"//str(iT)//"_T"//str(Temp,2)//"/"
         !
         write(*,"(A)") new_line("A")//"     ................................................"//new_line("A")
         write(*,"(3(A,1F20.5))") "     T(K): ",Temp,"    Beta(1/eV): ",Beta,"    Beta(1/"//DFTgrid//"): ",Beta_DFT
         !
         write(*,"(A)") new_line("A")//"     Computing Kernels."
         !
         !phononic Kernels
         if(calc_phonons)then
            !
            allocate(Zph(Ngrid));Zph=0d0
            call calc_Zph_e(Beta_DFT,Zph,reg(Inputs%mode_Zph),reg(printpath_T))
            !
            allocate(Kph(EE_dim));Kph=0d0
            call calc_Kph_e(Beta_DFT,Kph,reg(Inputs%printmode_ph),reg(printpath_T))
            !
         endif
         !
         !get the electronic Kernel depending on the input specifications
         if(calc_Kel)then
            !
            allocate(Kel(EE_dim));Kel=0d0
            call get_Kel(Kel,Beta_DFT,reg(Inputs%printmode_el),reg(printpath_T))
            !
         endif
         !
         !Convergence loop over Delta(e)
         write(*,"(A)") new_line("A")//"     Solving gap equation."
         call dump_Field_component(oldDelta,reg(printpath_T),"0_Delta.DAT",Egrid_Phonons)
         allocate(EDsq(Ngrid)); EDsq = czero
         allocate(newDelta(Ngrid)); newDelta = czero
         converged=.false.
         !
         SCloop: do iloop=1,Inputs%loops
            !
            do iE=1,Ngrid
               EDsq(iE) = sqrt( Egrid_Phonons(iE)**2 + conjg(Delta(iE))*Delta(iE) )
            enddo
            !
            newDelta = czero
            !$OMP PARALLEL DEFAULT(PRIVATE),&
            !$OMP SHARED(Ngrid,Beta_DFT,EDsq,Delta,newDelta),&
            !$OMP SHARED(calc_phonons,Egrid_Phonons,DoS_DFT,Kph,Zph),&
            !$OMP SHARED(calc_Kel,Egrid_Model,DoS_Model_interp,Kel,bilinear_map)
            !$OMP DO
            do iE1=1,Ngrid
               !
               Kint=czero
               do iE2=2,Ngrid
                  !
                  tanh_b = Beta_DFT/2d0
                  tanh_f = Beta_DFT/2d0
                  if(EDsq(iE2-1).ne.0d0) tanh_b = (tanh(Beta_DFT/2d0*EDsq(iE2-1))/EDsq(iE2-1))
                  if(EDsq(iE2).ne.0d0)   tanh_f = (tanh(Beta_DFT/2d0*EDsq(iE2))  /EDsq(iE2)  )
                  !
                  dE = Egrid_Phonons(iE2)-Egrid_Phonons(iE2-1)
                  !
                  !Integral of the phonon Kernel done with the DFT DoS - logarithmic mesh
                  if(calc_phonons)then
                     !
                     Kint = Kint + ( DoS_DFT(iE2-1) * Kph(rc2ut(iE1,iE2-1,Ngrid)) * tanh_b * Delta(iE2-1) + &
                                     DoS_DFT(iE2)   * Kph(rc2ut(iE1,iE2,Ngrid))   * tanh_f * Delta(iE2)   ) * (dE/2d0)
                     !
                  endif
                  !
                  !Integral of the electronic Kernel done with the Model DoS - bilinear interpolation to logarithmic mesh
                  if(calc_Kel)then
                     !
                     !default values
                     Kel_b = Kel(rc2ut(iE1,iE2-1,Ngrid))
                     Kel_f = Kel(rc2ut(iE1,iE2,Ngrid))
                     !
                     !the phonons are on a logarithmic mesh so bilinear interpolation is required
                     if(calc_phonons)then
                        !
                        Q_p = bilinear_map(iE1,iE2-1,:)
                        Kel_b = bilinear_interp( Egrid_Model(Q_p(1))             , & !x1 = y1
                                                 Egrid_Model(Q_p(2))             , & !x2
                                                 Egrid_Model(Q_p(1))             , & !y1
                                                 Egrid_Model(Q_p(2))             , & !y2
                                                 Kel(rc2ut(Q_p(1),Q_p(1),Ngrid)) , & !f(x1,y1)
                                                 Kel(rc2ut(Q_p(1),Q_p(2),Ngrid)) , & !f(x1,y2)
                                                 Kel(rc2ut(Q_p(2),Q_p(2),Ngrid)) , & !f(x2,y2)
                                                 Kel(rc2ut(Q_p(2),Q_p(1),Ngrid)) , & !f(x2,y1)
                                                 Egrid_Phonons(iE2)              , & !x
                                                 Egrid_Phonons(iE2-1)              ) !y
                        !
                        Q_f = bilinear_map(iE1,iE2,:)
                        Kel_f = bilinear_interp( Egrid_Model(Q_f(1))             , & !x1 = y1
                                                 Egrid_Model(Q_f(2))             , & !x2
                                                 Egrid_Model(Q_f(1))             , & !y1
                                                 Egrid_Model(Q_f(2))             , & !y2
                                                 Kel(rc2ut(Q_f(1),Q_f(1),Ngrid)) , & !f(x1,y1)
                                                 Kel(rc2ut(Q_f(1),Q_f(2),Ngrid)) , & !f(x1,y2)
                                                 Kel(rc2ut(Q_f(2),Q_f(2),Ngrid)) , & !f(x2,y2)
                                                 Kel(rc2ut(Q_f(2),Q_f(1),Ngrid)) , & !f(x2,y1)
                                                 Egrid_Phonons(iE2)              , & !x
                                                 Egrid_Phonons(iE2)                ) !y
                        !
                     endif
                     !
                     Kint = Kint + ( DoS_Model_interp(iE2-1) * Kel_b * tanh_b * Delta(iE2-1) + &
                                     DoS_Model_interp(iE2)   * Kel_f * tanh_f * Delta(iE2)   ) * (dE/2d0)
                     !
                  endif
                  !
               enddo
               !
               newDelta(iE1) = - 0.5d0*dreal(Kint)
               if(calc_phonons) newDelta(iE1) = newDelta(iE1) - Zph(iE1)*Delta(iE1)
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            !
            !This is to fix the phase of Delta
            newDelta = dcmplx(dreal(newDelta),0d0)
            if(calc_Kel) where(DoS_Model_interp.lt.Inputs%DoSthresh) newDelta=czero
            !
            !Convergence check
            !errDelta = maxval(abs(Delta-newDelta))
            errDelta = min(maxval(abs(Delta-newDelta)),maxval(abs(Delta+newDelta)))
            !
            if(errDelta.lt.Inputs%DeltaErr)then
               if(verbose) write(*,"(2(A,1E20.10),A3,1E20.10)")"     loop #"//str(iloop)//" Delta(0): ",abs(newDelta(1+minloc(abs(Egrid_Phonons),dim=1))),"   error: ",errDelta," < ",Inputs%DeltaErr
               write(*,"(A)")"     Delta at T "//str(Temp,2)//"K is converged in "//str(iloop)//" loops. Moving to next Temperature."
               converged=.true.
               exit SCloop
            else
               if(verbose) write(*,"(2(A,1E20.10),A3,1E20.10)")"     loop #"//str(iloop)//" Delta(0): ",abs(newDelta(1+minloc(abs(Egrid_Phonons),dim=1))),"   error: ",errDelta," > ",Inputs%DeltaErr
               Delta = (1d0-Inputs%DeltaMix)*newDelta + Inputs%DeltaMix*oldDelta
            endif
            oldDelta = Delta
            !
            call dump_Field_component(Delta,reg(printpath_T),str(iloop)//"_Delta.DAT",Egrid_Phonons)
            !
         enddo SCloop !iloop
         !
         deallocate(EDsq,newDelta)
         if(allocated(Zph))deallocate(Zph)
         if(allocated(Kph))deallocate(Kph)
         if(allocated(Kel))deallocate(Kel)
         !
         if(.not.converged)write(*,"(A)")"     Warning: Delta is not converged after "//str(iloop)//" loops. Moving to next Temperature."
         Delta = Delta*DFTgrid2eV
         allocate(Egrid_print(Ngrid));Egrid_print = Egrid_Phonons*DFTgrid2eV
         call dump_Field_component(Delta,reg(printpath_T),"Delta_e_T"//str(Temp,2)//".DAT",Egrid_print)
         Tlist(iT) = Temp
         Delta_T(iT) = abs(Delta(1+minloc(abs(Egrid_print),dim=1)))
         deallocate(Egrid_print)
         !
      enddo !iT
      !
      call dump_Field_component(Delta_T,reg(printpath),"Delta_T.DAT",Tlist)
      deallocate(Tlist,Delta_T,DoS_Model_interp)
      !
   endif
   !
end subroutine calc_Tc
