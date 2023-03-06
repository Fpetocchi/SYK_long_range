subroutine calc_Tc(pathOUTPUT,Inputs,Lttc,Wlat)
   !
   use parameters
   use linalg, only : zeye
   use utils_misc
   use utils_fields
   use crystal
   use file_io
   use gap_equation
   use input_vars, only : Nreal, wrealMax, pathINPUT, pathINPUTtr
   implicit none
   !
   character(len=*),intent(in)           :: pathOUTPUT
   type(SCDFT),intent(in)                :: Inputs
   type(Lattice),intent(in)              :: Lttc
   type(BosonicField),intent(in),optional   :: Wlat
   !
   real(8)                               :: Beta_input,dE
   integer                               :: iorb,ik,iE,iE1,iE2
   integer                               :: Norb,Nkpt,Ngrid
   complex(8)                            :: Kint
   real(8),allocatable                   :: Zph(:),Kph(:,:)
   complex(8),allocatable                :: Kel_stat(:,:),Kel_dyn(:,:)
   complex(8),allocatable                :: Z(:),K(:,:)
   complex(8),allocatable                :: Hk_used(:,:,:),SfullCorr(:,:,:)
   character(len=255)                    :: printpath
   !
   integer                               :: iloop,iT
   real(8)                               :: Temp,Beta,errDelta
   real(8)                               :: tanh_b,tanh_f
   real(8),allocatable                   :: ED(:),kpt_QP(:,:)
   complex(8),allocatable                :: Delta(:),oldDelta(:),newDelta(:)
   logical                               :: converged,filexists
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Tc"
   !
   !
   printpath=reg(pathOUTPUT)//"Gap_Equation/"
   call createDir(reg(printpath),verb=verbose)
   !
   !Dimensionla checks
   Norb = Lttc%Norb
   Nkpt = Lttc%Nkpt
   write(*,"(A,1F15.3)") "     W:   ",wrealMax
   !
   allocate(Hk_used(Norb,Norb,Nkpt));Hk_used=czero
   Hk_used = Lttc%Hk
   !
   !Only needed if one wants to include electronic Kernels
   if(present(Wlat))then
      !
      Beta_input = Wlat%Beta


      Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)


      !
      if(.not.Wlat%status) stop "calc_Tc: Wlat not properly initialized."
      if(Wlat%Nkpt.eq.Nkpt) stop "calc_Tc: Wlat has different number of k-points with respect to Lttc."
      if(Wlat%Nbp.ne.(Norb**2)) stop "calc_Tc: Wlat has different orbital dimension with respect to Lttc."
      !
      !Use the QP bandstructure
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
      endif
      !
   else
      !
      if (calc_Int_static.or.calc_Int_full) stop "calc_Tc: electronic Kernels requested but Wlat not provided."
      !
   endif
   !
   call Initialize_inputs(reg(pathINPUT),reg(Inputs%mode_ph),reg(Inputs%mode_el), &
                          Lttc%Nkpt3,Lttc%kpt,Hk_used,Nreal,wrealMax,             &
                          Nkpt3_intp_Hk=Inputs%Nkpt3_intp_Hk,Nkpt3_intp_Wk=Inputs%Nkpt3_intp_Wk)
   deallocate(Hk_used)
   !
   Ngrid = size(Egrid)
   if(.not.associated(DoS_Hk))stop"calc_Tc: DoS_Hk pointer not associated."
   call dump_Field_component(DoS_Hk,reg(printpath),"DoS_Hk.DAT",Egrid)
   if (calc_Int_static.or.calc_Int_full) call dump_Field_component(DoS_Wk,reg(pathOUTPUT)//"Gap_Equation","DoS_Wk.DAT",Egrid)
   !
   !Allocating delta here so as to imply annealing between temperatures
   allocate(Delta(Ngrid));Delta = Inputs%DeltaInit
   do iE=1,Ngrid
      if(abs(Egrid(iE)).gt.(0.1d0*wrealMax))Delta(iE)=czero
   enddo
   !
   !TEMPERATURE LOOP
   write(*,"(A)") new_line("A")//"---- Starting temperature scan"
   do iT=1,Inputs%Tsteps
      !
      Temp = Inputs%Tbounds(1) + (iT-1)*abs(Inputs%Tbounds(2)-Inputs%Tbounds(1))/dble(Inputs%Tsteps-1)
      Beta = 1d0 / (Temp*K2eV)
      !
      write(*,"(A)") new_line("A")//"     ...................................."//new_line("A")
      write(*,"(A,1F15.3)") "     T:   ",Temp
      write(*,"(A,1F15.3)") "     Beta:",Beta
      !
      write(*,"(A)") new_line("A")//"     Computing Kernels."
      allocate(Z(Ngrid));Z=czero
      allocate(K(Ngrid,Ngrid));K=czero
      !
      if(calc_phonons)then
         !
         allocate(Zph(Ngrid));Zph=0d0
         allocate(Kph(Ngrid,Ngrid));Kph=0d0
         !
         call calc_Zph_e(Beta,Egrid,DoS_Hk,Zph,mode=reg(Inputs%mode_Zph),printZpath=reg(printpath))
         call calc_Kph_e(Beta,Egrid,DoS_Hk,Kph,printKpath=reg(printpath),printmode=reg(Inputs%printmode_ph))
         !
         Z = Z + Zph
         K = K + Kph
         !
         deallocate(Zph,Kph)
         !
      endif
      !
      if(calc_Int_static)then
         !
         if(.not.associated(weights_Wk))stop"calc_Tc: weights_Wk pointer not associated."
         allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=0d0
         !
         call store_Wk(Wlat%screened,Beta,eps)
         call calc_Kel_stat_e(Beta,Egrid,weights_Wk,Kel_stat,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
         !
         K = K + Kel_stat
         !
         deallocate(Kel_stat)
         !
      elseif(calc_Int_full)then
         !
         if(.not.associated(weights_Wk))stop"calc_Tc: weights_Wk pointer not associated."
         allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=0d0
         allocate(Kel_dyn(Ngrid,Ngrid));Kel_dyn=0d0
         !
         call store_Wk(Wlat%screened,Beta,Inputs%Wk_cutoff)
         call calc_Kel_stat_e(Beta,Egrid,weights_Wk,Kel_stat,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
         call calc_Kel_dyn_e(Beta,Egrid,weights_Wk,Kel_dyn,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
         !
         K = K + Kel_stat + Kel_dyn
         !
         deallocate(Kel_stat,Kel_dyn)
         !
      endif
      !
      write(*,"(A)") new_line("A")//"     Solving gap equation."
      allocate(ED(Ngrid));ED = czero
      allocate(newDelta(Ngrid));newDelta = czero
      allocate(oldDelta(Ngrid));oldDelta = Delta
      converged=.false.
      SCloop: do iloop=1,Inputs%loops
         !
         !$OMP PARALLEL DEFAULT(PRIVATE),&
         !$OMP SHARED(Ngrid,Egrid,Beta,DoS_Hk,Z,K,ED,Delta,newDelta)
         !
         !$OMP DO SCHEDULE(STATIC)
         do iE=1,Ngrid
           ED(iE) = sqrt( Egrid(iE)**2 + conjg(Delta(iE))*Delta(iE) )
         enddo
         !$OMP END DO
         !
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
              tanh_b = Beta/2d0
              tanh_f = Beta/2d0
              if(ED(iE2-1).ne.0d0) tanh_b = (tanh(Beta/2d0*ED(iE2-1))/ED(iE2-1))
              if(ED(iE2).ne.0d0)   tanh_f = (tanh(Beta/2d0*ED(iE2))  /ED(iE2)  )
              !
              SPLIT IT UP IN 2 SUMS ONE FOR THE KPHONON WHERE THE DOS FROM DFT IS USED AND THE ELECTRONINC ONE WHERE I CAN USE MINE
              AND



              Kint = Kint + ( DoS_Hk(iE2-1) * K(iE1,iE2-1) * tanh_b * Delta(iE2-1) + &
                              DoS_Hk(iE2)   * K(iE1,iE2)   * tanh_f * Delta(iE2)   ) * (dE/2d0)
              !
            enddo
            !
            newDelta(iE1) = -Z(iE1)*Delta(iE1) - 0.5d0*Kint
            !
         enddo
         !$OMP END DO
         !$OMP BARRIER
         !$OMP END PARALLEL
         !
         !Convergence check
         errDelta = maxval(abs(Delta-newDelta))
         if(errDelta.lt.Inputs%DeltaErr)then
            write(*,"(2(A,1E20.10),A3,1E20.10)")"     loop #"//str(iloop)//" Delta: ",maxval(abs(dble(newDelta)),abs(Egrid).lt.1d-2),"   error: ",errDelta," < ",Inputs%DeltaErr
            write(*,"(A)")"     Delta is converged moving to next Temperature."
            converged=.true.
            Delta = (1d0-Inputs%DeltaMix)*newDelta + Inputs%DeltaMix*oldDelta
            oldDelta = Delta
            exit SCloop
         else
            write(*,"(2(A,1E20.10),A3,1E20.10)")"     loop #"//str(iloop)//" Delta: ",maxval(abs(dble(newDelta)),abs(Egrid).lt.1d-2),"   error: ",errDelta," > ",Inputs%DeltaErr
         endif
         !
         Delta = (1d0-Inputs%DeltaMix)*newDelta + Inputs%DeltaMix*oldDelta
         oldDelta = Delta
         !
      enddo SCloop !iloop
      deallocate(ED,Z,K,newDelta,oldDelta)
      !
      if(.not.converged)write(*,"(A)")"     Warning: Delta is not converged."
      call dump_Field_component(Delta,reg(pathOUTPUT)//"Gap_Equation/","Delta_T"//str(Temp,2)//".DAT",Egrid)
      !
   enddo !iT
   !
end subroutine calc_Tc
