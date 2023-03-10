subroutine calc_Tc(pathOUTPUT,Inputs,Lttc,Wlat)
   !
   use parameters
   use linalg, only : zeye
   use utils_misc
   use utils_fields
   use crystal
   use file_io
   use gap_equation
   use input_vars, only : pathINPUT, pathINPUTtr
   implicit none
   !
   character(len=*),intent(in)           :: pathOUTPUT
   type(SCDFT),intent(in)                :: Inputs
   type(Lattice),intent(in)              :: Lttc
   type(BosonicField),intent(in)         :: Wlat
   !
   real(8)                               :: dE
   integer                               :: iE,iE1,iE2
   integer                               :: Norb,Nkpt,Ngrid
   complex(8)                            :: Kint
   real(8),allocatable                   :: Zph(:),Kph(:,:)
   complex(8),allocatable,target         :: Kel_stat(:,:),Kel_dyn(:,:)
   complex(8),pointer                    :: K(:,:)
   complex(8),allocatable                :: Hk_used(:,:,:)
   character(len=255)                    :: printpath
   !
   integer                               :: iloop,iT
   real(8)                               :: Temp,Beta,Beta_DFT,errDelta
   real(8)                               :: tanh_b,tanh_f
   real(8),allocatable                   :: EDsq(:),kpt_QP(:,:)
   complex(8),allocatable                :: Delta(:),oldDelta(:),newDelta(:)
   logical                               :: converged,filexists
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Tc"
   !
   !
   printpath=reg(pathOUTPUT)//"Gap_Equation/"
   call createDir(reg(printpath),verb=verbose)
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
      Hk_used = Lttc%Hk
      !
   endif
   !
   call Initialize_inputs(reg(pathINPUT),Inputs,Lttc,Hk_used)
   deallocate(Hk_used)
   !
   call dump_Field_component(DoS_DFT,reg(printpath),"DoS_DFT.DAT",Egrid)
   call dump_Field_component(DoS_Model,reg(printpath),"DoS_Model.DAT",Egrid)
   !
   if (calc_Int_static.or.calc_Int_dynamic) then
      call store_Wk4gap(Wlat%screened,Lttc,Wlat%Beta,Inputs%Wk_cutoff,reg(printpath))
   endif
   !
   if(Inputs%calc_Tc)then
      !
      !Allocating delta here so as to imply annealing between temperatures
      Ngrid = size(Egrid)
      allocate(Delta(Ngrid));Delta = Inputs%DeltaInit*eV2DFTgrid
      do iE=1,Ngrid
         if(abs(Egrid(iE)).gt.(0.1d0*Inputs%wrealMax*eV2DFTgrid))Delta(iE)=czero
      enddo
      allocate(EDsq(Ngrid)); EDsq = czero
      !
      !
      !============================ TEMPERATURE LOOP =============================!
      !
      write(*,"(A)") new_line("A")//"---- Starting temperature scan"
      do iT=1,Inputs%Tsteps
         !
         Temp = Inputs%Tbounds(2) - (iT-1)*abs(Inputs%Tbounds(2)-Inputs%Tbounds(1))/dble(Inputs%Tsteps-1)
         Beta = 1d0 / (Temp*K2eV)
         Beta_DFT = 1d0 / (Temp*K2eV*eV2DFTgrid)
         !
         write(*,"(A)") new_line("A")//"     ...................................."//new_line("A")
         write(*,"(A,1F15.3)") "     T:   ",Temp
         write(*,"(2(A,1F15.3))") "     Beta:",Beta," Beta(DFT):",Beta_DFT
         !
         write(*,"(A)") new_line("A")//"     Computing Kernels."
         !
         if(calc_phonons)then
            !
            allocate(Zph(Ngrid));Zph=0d0
            call calc_Zph_e(Beta_DFT,Zph,mode=reg(Inputs%mode_Zph),printZpath=reg(printpath))
            !
            allocate(Kph(Ngrid,Ngrid));Kph=0d0
            call calc_Kph_e(Beta_DFT,Kph,printmode=reg(Inputs%printmode_ph),printKpath=reg(printpath))
            !
         endif
         !
         if(calc_Int_static.and.calc_Int_dynamic)then
            !
            allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=czero
            call calc_Kel_stat_e(Beta_DFT,Kel_stat,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
            !
            allocate(Kel_dyn(Ngrid,Ngrid));Kel_dyn=czero
            call calc_Kel_dyn_e(Beta_DFT,Kel_dyn,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
            !
            Kel_dyn = Kel_stat + Kel_dyn
            K => Kel_dyn
            !
         elseif(calc_Int_static) then
            !
            allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=czero
            call calc_Kel_stat_e(Beta_DFT,Kel_stat,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
            K => Kel_stat
            !
         elseif(calc_Int_dynamic) then
            !
            allocate(Kel_dyn(Ngrid,Ngrid));Kel_dyn=czero
            call calc_Kel_dyn_e(Beta_DFT,Kel_dyn,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
            K => Kel_dyn
            !
         endif
         !
         !Convergence loop over Delta(e)
         write(*,"(A)") new_line("A")//"     Solving gap equation."
         !
         allocate(newDelta(Ngrid));newDelta = czero
         allocate(oldDelta(Ngrid));oldDelta = Delta
         converged=.false.
         SCloop: do iloop=1,Inputs%loops
            !
            do iE=1,Ngrid
              EDsq(iE) = sqrt( Egrid(iE)**2 + conjg(Delta(iE))*Delta(iE) )
            enddo
            !
            !$OMP PARALLEL DEFAULT(PRIVATE),&
            !$OMP SHARED(Ngrid,Egrid,Beta_DFT,DoS_Model,DoS_DFT,Zph,Kph,K,EDsq,Delta,newDelta)
            newDelta = czero
            !$OMP DO SCHEDULE(DYNAMIC)
            do iE1=1,Ngrid
               !
               !integral over E2
               Kint=czero
               do iE2=2,Ngrid
                 !
                 dE=Egrid(iE2)-Egrid(iE2-1)
                 !
                 tanh_b = Beta_DFT/2d0
                 tanh_f = Beta_DFT/2d0
                 if(EDsq(iE2-1).ne.0d0) tanh_b = (tanh(Beta_DFT/2d0*EDsq(iE2-1))/EDsq(iE2-1))
                 if(EDsq(iE2).ne.0d0)   tanh_f = (tanh(Beta_DFT/2d0*EDsq(iE2))  /EDsq(iE2)  )
                 !
                 !Integral of the phonon Kernel done with the DFT DoS
                 if(calc_phonons)then
                    Kint = Kint + ( DoS_DFT(iE2-1) * Kph(iE1,iE2-1) * tanh_b * Delta(iE2-1) + &
                                    DoS_DFT(iE2)   * Kph(iE1,iE2)   * tanh_f * Delta(iE2)   ) * (dE/2d0)
                 endif
                 !
                 !Integral of the electronic Kernel done with the custom DoS
                 if(calc_Int_static.or.calc_Int_dynamic)then
                    Kint = Kint + ( DoS_Model(iE2-1) * K(iE1,iE2-1) * tanh_b * Delta(iE2-1) + &
                                    DoS_Model(iE2)   * K(iE1,iE2)   * tanh_f * Delta(iE2)   ) * (dE/2d0)
                 endif
                 !
               enddo
               !
               newDelta(iE1) = - 0.5d0*Kint
               if(calc_phonons) newDelta(iE1) = newDelta(iE1) - Zph(iE1)*Delta(iE1)
               !
            enddo
            !$OMP END DO
            !$OMP BARRIER
            !$OMP END PARALLEL
            !
            !Convergence check
            errDelta = maxval(abs(Delta-newDelta))
            if(errDelta.lt.Inputs%DeltaErr)then
               write(*,"(2(A,1E20.10),A3,1E20.10)")"     loop #"//str(iloop)//" Delta(0): ",abs(newDelta(minloc(abs(Egrid),dim=1))),"   error: ",errDelta," < ",Inputs%DeltaErr
               write(*,"(A)")"     Delta is converged moving to next Temperature."
               converged=.true.
               oldDelta = Delta
               exit SCloop
            else
               write(*,"(2(A,1E20.10),A3,1E20.10)")"     loop #"//str(iloop)//" Delta(0): ",abs(newDelta(minloc(abs(Egrid),dim=1))),"   error: ",errDelta," > ",Inputs%DeltaErr
            endif
            !
            Delta = (1d0-Inputs%DeltaMix)*newDelta + Inputs%DeltaMix*oldDelta
            oldDelta = Delta
            !
         enddo SCloop !iloop
         deallocate(EDsq,newDelta,oldDelta)
         if(calc_phonons)deallocate(Zph,Kph)
         if(associated(K))nullify(K)
         !
         if(.not.converged)write(*,"(A)")"     Warning: Delta is not converged."
         call dump_Field_component(Delta,reg(pathOUTPUT)//"Gap_Equation/","Delta_T"//str(Temp,2)//".DAT",Egrid)
         !
      enddo !iT
      !
   endif
   !
end subroutine calc_Tc
