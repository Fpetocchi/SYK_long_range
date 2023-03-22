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
   real(8)                               :: dE,dT
   integer                               :: iE,iE1,iE2
   integer                               :: Norb,Nkpt,Ngrid
   complex(8)                            :: Kint
   real(8),allocatable                   :: Egrid_print(:)
   real(8),allocatable                   :: Tlist(:),Delta_T(:)
   real(8),allocatable                   :: Zph(:),Kph(:,:)
   complex(8),allocatable,target         :: Kel_stat(:,:),Kel_dyn(:,:)
   complex(8),pointer                    :: Kel(:,:)
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
      printpath=reg(pathOUTPUT)//"Gap_Equation_QP/"
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
   call dump_Field_component(DoS_DFT,reg(printpath),"DoS_DFT.DAT",Egrid)
   call dump_Field_component(DoS_Model,reg(printpath),"DoS_Model.DAT",Egrid)
   !
   if (calc_Int_static.or.calc_Int_dynamic) then
      call store_Wk4gap(Wlat%screened,Lttc,Wlat%Beta,Inputs%Wk_cutoff,reg(printpath),.true.)
      !QUI CAMBIA TUTTO INVECE CHE STORARE WK STORA I DUE INTEGRALI SU E ED EPRIMO TANTO QUELLO CHE DIPENDE DALLA TEMPERATURA É SOLO KELDYN
      !QUI DEVI CAMBIARE PROPRIO TUTTO
   endif
   !
   if(Inputs%calc_Tc)then
      !
      Ngrid = size(Egrid)
      allocate(Delta(Ngrid));Delta = czero
      do iE=1,Ngrid
         Delta(iE) = -diff_fermidirac(Egrid(iE),0d0,(50d0/(Inputs%wrealMax*eV2DFTgrid)))
      enddo
      Delta = Delta * Inputs%DeltaInit * eV2DFTgrid / Delta(minloc(abs(Egrid),dim=1))
      allocate(oldDelta(Ngrid));oldDelta=Delta
      !
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
         write(*,"(A)") new_line("A")//"     ................................................"//new_line("A")
         write(*,"(3(A,1F12.5))") "     T(K): ",Temp,"    Beta(1/eV): ",Beta,"    Beta(1/"//DFTgrid//"): ",Beta_DFT
         !
         write(*,"(A)") new_line("A")//"     Computing Kernels."
         !
         if(calc_phonons)then
            !
            allocate(Zph(Ngrid));Zph=0d0
            call calc_Zph_e(Beta_DFT,Zph,reg(Inputs%mode_Zph),reg(printpath))
            !
            allocate(Kph(Ngrid,Ngrid));Kph=0d0
            call calc_Kph_e(Beta_DFT,Kph,reg(Inputs%printmode_ph),reg(printpath))
            !
         endif
         !
         if(calc_Int_static.and.calc_Int_dynamic)then
            !
            allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=czero
            call calc_Kel_stat_e(Beta_DFT,Kel_stat,reg(Inputs%printmode_el),reg(printpath))
            !
            allocate(Kel_dyn(Ngrid,Ngrid));Kel_dyn=czero
            call calc_Kel_dyn_e(Beta_DFT,Kel_dyn,reg(Inputs%printmode_el),reg(printpath))
            !
            Kel_dyn = Kel_stat + Kel_dyn
            Kel => Kel_dyn
            !
         elseif(calc_Int_static) then
            !
            allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=czero
            call calc_Kel_stat_e(Beta_DFT,Kel_stat,reg(Inputs%printmode_el),reg(printpath))
            Kel => Kel_stat
            !
         elseif(calc_Int_dynamic) then
            !
            allocate(Kel_dyn(Ngrid,Ngrid));Kel_dyn=czero
            call calc_Kel_dyn_e(Beta_DFT,Kel_dyn,reg(Inputs%printmode_el),reg(printpath))
            Kel => Kel_dyn
            !
         endif
         !
         !Convergence loop over Delta(e)
         write(*,"(A)") new_line("A")//"     Solving gap equation."
         !
         call dump_Field_component(oldDelta,reg(printpath)//"loops_T"//str(Temp,2)//"/","0_Delta.DAT",Egrid)
         !
         allocate(EDsq(Ngrid)); EDsq = czero
         allocate(newDelta(Ngrid));newDelta = czero
         converged=.false.
         SCloop: do iloop=1,Inputs%loops
            !
            do iE=1,Ngrid
              EDsq(iE) = sqrt( Egrid(iE)**2 + conjg(Delta(iE))*Delta(iE) )
            enddo
            !
            newDelta = czero
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
                    Kint = Kint + ( DoS_Model(iE2-1) * Kel(iE1,iE2-1) * tanh_b * Delta(iE2-1) + &
                                    DoS_Model(iE2)   * Kel(iE1,iE2)   * tanh_f * Delta(iE2)   ) * (dE/2d0)
                 endif
                 !
               enddo
               !
               newDelta(iE1) = - 0.5d0*Kint
               if(calc_phonons) newDelta(iE1) = newDelta(iE1) - Zph(iE1)*Delta(iE1)
               !
            enddo
            !
            !Convergence check
            errDelta = maxval(abs(Delta-newDelta))
            if(errDelta.lt.Inputs%DeltaErr)then
               write(*,"(2(A,1E20.10),A3,1E20.10)")"     loop #"//str(iloop)//" Delta(0): ",abs(newDelta(minloc(abs(Egrid),dim=1))),"   error: ",errDelta," < ",Inputs%DeltaErr
               write(*,"(A)")"     Delta at T "//str(Temp,2)//"K is converged. Moving to next Temperature."
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
            call dump_Field_component(Delta,reg(printpath)//"loops_T"//str(Temp,2)//"/",str(iloop)//"_Delta.DAT",Egrid)
            !
         enddo SCloop !iloop
         deallocate(EDsq,newDelta)
         if(associated(Kel))nullify(Kel)
         if(allocated(Zph))deallocate(Zph)
         if(allocated(Kph))deallocate(Kph)
         if(allocated(Kel_stat))deallocate(Kel_stat)
         if(allocated(Kel_dyn))deallocate(Kel_dyn)
         !
         if(.not.converged)write(*,"(A)")"     Warning: Delta is not converged."
         Delta = Delta*DFTgrid2eV
         allocate(Egrid_print(Ngrid));Egrid_print = Egrid*DFTgrid2eV
         call dump_Field_component(Delta,reg(printpath),"Delta_e_T"//str(Temp,2)//".DAT",Egrid_print)
         Tlist(iT) = Temp
         Delta_T(iT) = abs(Delta(minloc(abs(Egrid_print),dim=1)))
         deallocate(Egrid_print)
         !
      enddo !iT
      !
      call dump_Field_component(Delta_T,reg(printpath),"Delta_T.DAT",Tlist)
      deallocate(Tlist,Delta_T)
      !
   endif
   !
end subroutine calc_Tc
