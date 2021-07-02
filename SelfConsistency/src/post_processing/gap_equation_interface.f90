subroutine calc_Tc(pathOUTPUT,Inputs,Lttc,Wlat,Sfull)
   !
   use parameters
   use linalg, only : zeye
   use utils_misc
   use utils_fields
   use crystal
   use file_io
   use gap_equation
   use input_vars, only : Nreal, wrealMax, pathINPUT
   implicit none
   !
   character(len=*),intent(in)           :: pathOUTPUT
   type(SCDFT),intent(in)                :: Inputs
   type(Lattice),intent(in)              :: Lttc
   type(BosonicField),intent(in),optional   :: Wlat
   type(FermionicField),intent(in),optional :: Sfull
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
   real(8),allocatable                   :: ED(:)
   complex(8),allocatable                :: Delta(:),newDelta(:)
   logical                               :: converged
   !
   write(*,"(A)") new_line("A")//new_line("A")//"---- calc_Tc"
   !
   wrealMax = 17d0
   print *,"wrealMax",wrealMax
   !
   printpath=reg(pathOUTPUT)//"Gap_Equation/"
   call createDir(reg(printpath),verb=verbose)
   !
   !Dimensionla checks
   Norb = Lttc%Norb
   Nkpt = Lttc%Nkpt
   !
   allocate(Hk_used(Norb,Norb,Nkpt));Hk_used=czero
   Hk_used = Lttc%Hk
   !
   !Only needed if one wants to include electronic Kernels
   if(present(Wlat))then
      !
      Beta_input = Wlat%Beta
      !
      if(.not.Wlat%status) stop "calc_Tc: Wlat not properly initialized."
      if(Wlat%Nkpt.eq.Nkpt) stop "calc_Tc: Wlat has different number of k-points with respect to Lttc."
      if(Wlat%Nbp.ne.(Norb**2)) stop "calc_Tc: Wlat has different orbital dimension with respect to Lttc."
      !
      !Adding the renormalization from the Self-energy at iw=0
      if(present(Sfull))then
         !
         if(.not.Sfull%status) stop "calc_Tc: Sfull not properly initialized."
         if(Sfull%Nkpt.eq.Nkpt) stop "calc_Tc: Sfull has different number of k-points with respect to Lttc."
         if(Sfull%Norb.ne.Norb) stop "calc_Tc: Sfull has different orbital dimension with respect to Lttc."
         if(Sfull%Beta.ne.Beta_input) stop "calc_Tc: Sfull has different orbital dimension with respect to Wlat."
         !
         allocate(SfullCorr(Norb,Norb,Nkpt));SfullCorr=czero
         SfullCorr = Sfull%wks(:,:,1,:,1)
         do iorb=1,Norb
            SfullCorr(iorb,iorb,:) = dcmplx(dreal(SfullCorr(iorb,iorb,:)),0d0)
         enddo
         !
         do ik=1,Nkpt
            Hk_used(:,:,ik) = Lttc%Hk(:,:,ik) + SfullCorr(:,:,ik) - Sfull%mu*zeye(Norb)
         enddo
         deallocate(SfullCorr)
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
      allocate(Z(Ngrid));Z=0d0
      allocate(K(Ngrid,Ngrid));K=0d0
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
      !stop
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
         call calc_Kel_dyn_e(Beta,Egrid,weights_Wk,Inputs%wstep,Kel_dyn,printKpath=reg(printpath),printmode=reg(Inputs%printmode_el))
         !
         K = K + Kel_stat + Kel_dyn
         !
         deallocate(Kel_stat,Kel_dyn)
         !
      endif
      !
      write(*,"(A)") new_line("A")//"     Solving gap equation."
      allocate(ED(Ngrid));ED = czero
      converged=.false.
      SCloop: do iloop=1,Inputs%loops
         !
         write(*,"(A)")"     self-consistency loop #"//str(iloop)//" starting parallel."
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
              if(ED(iE2-1).ne.0d0) tanh_f = (tanh(Beta/2d0*ED(iE2))  /ED(iE2)  )
              !
              Kint = Kint - dE * ( DoS_Hk(iE2-1) * K(iE1,iE2-1) * tanh_b * Delta(iE2-1) + &
                                   DoS_Hk(iE2)   * K(iE1,iE2)   * tanh_f * Delta(iE2)   )
              !
            enddo
            Kint = Kint/2d0 !From Trapezoid integration
            !
            newDelta(iE1) = -Z(iE1)*Delta(iE1) - 0.5d0*Kint
            !
         enddo
         !$OMP END DO
         !$OMP BARRIER
         !$OMP END PARALLEL
         !
         write(*,"(A,1E20.10)")"     Delta: ",maxval(abs(dble(Delta)),abs(Egrid).lt.1d-2)
         !
         errDelta = maxval(abs(Delta-newDelta))
         if(errDelta.lt.Inputs%DeltaErr)then
            write(*,"(A,1E20.10,A3,1E20.10)")"     error on Delta: ",errDelta," < ",Inputs%DeltaErr
            write(*,"(A)")"     Delta is converged moving to next Temperature."
            converged=.true.
            exit SCloop
         else
            write(*,"(A,1E20.10,A3,1E20.10)")"     error on Delta: ",errDelta," > ",Inputs%DeltaErr
         endif
         !
      enddo SCloop !iloop
      deallocate(ED,Z,K)
      !
      if(.not.converged)write(*,"(A)")"     Warning: Delta is not converged."
      call dump_Field_component(newDelta,reg(pathOUTPUT)//"Gap_Equation","Delta_T"//str(Temp)//".DAT",Egrid)
      !
   enddo !iT
   !
end subroutine calc_Tc
