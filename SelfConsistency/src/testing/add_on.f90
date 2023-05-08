!
!TEST>>>
path = reg(pathOUTPUT)//"Sigma_interp_ik1.DAT"
unit = free_unit()
open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
write(unit,"(I5,A)") Norb," Number of Wannier functions"
write(unit,"(I5,A)") Nreal_sigma," Number of grid points"
write(unit,"(1E20.12,A)") Lttc%mu," chemical potential"
write(unit,"(A)") "Wannier-projected fermionic components:"
do iw=1,Nreal_sigma
   do iorb=1,Norb
      do jorb=1,Norb
         write(unit,"(1E20.12,2I4,2E20.12)") Sigma_axis(iw),iorb,jorb,dreal(Sigma_intp(iorb,jorb,iw,1)),dimag(Sigma_intp(iorb,jorb,iw,1))
      enddo
   enddo
enddo
close(unit)
!
unit = free_unit()
open(unit,file=reg(pathOUTPUT)//"Hk_interp_ik1.DAT",form="formatted",status="unknown",position="rewind",action="write")
write(unit,("(3I10)")) 1,Lttc%Nkpt_path,Norb
write(unit,("(2I6,3F14.8)")) 1,1,Lttc%kptpath(:,1)
do iorb=1,Norb
   do jorb=1,Norb
      write(unit,("(2I4,2E20.12)")) iorb,jorb,dreal(Lttc%Hk_path(iorb,jorb,1)),dimag(Lttc%Hk_path(iorb,jorb,1))
   enddo
enddo
close(unit)
!>>>TEST
!








case("Hartree_lat_Nimp","Hartree_lat_Nlat")   ! DEPRECATED
   !
   !try to see if the SPEX Hartree is present otherwise use curlyU(0)
   call inquireFile(reg(PrevItFolder)//"Hartree_lat.DAT",filexists,hardstop=.false.,verb=verbose)
   if(filexists)then
      !
      if(allocated(Hartree_lat))deallocate(Hartree_lat)
      allocate(Hartree_lat(Crystal%Norb,Crystal%Norb));Hartree_lat=czero
      call read_Matrix(Hartree_lat,reg(PrevItFolder)//"Hartree_lat.DAT")
      do ispin=1,Nspin
         if(RotateHloc)then
            call loc2imp(Simp(isite)%N_s(:,:,ispin),Hartree_lat,LocalOrbs(isite)%Orbs,U=LocalOrbs(isite)%Rot)
         else
            call loc2imp(Simp(isite)%N_s(:,:,ispin),Hartree_lat,LocalOrbs(isite)%Orbs)
         endif
      enddo
      deallocate(Hartree_lat)
      !
   else
      !
      allocate(rho(LocalOrbs(isite)%Norb,LocalOrbs(isite)%Norb,Nspin));rho=0d0
      if(reg(DC_type).eq."Hartree_lat_Nimp")then
         rho = LocalOrbs(isite)%rho_OrbSpin
      elseif(reg(DC_type).eq."Hartree_lat_Nlat")then
         call read_Matrix(rho,reg(PrevItFolder)//"Solver_"//reg(LocalOrbs(isite)%Name)//"/Nloc_"//reg(LocalOrbs(isite)%Name),paramagnet)
      endif
      !
      Simp(isite)%N_s = czero
      do ispin=1,Nspin
         do iorb=1,LocalOrbs(isite)%Norb
            do jorb=1,LocalOrbs(isite)%Norb
               do korb=1,LocalOrbs(isite)%Norb
                  do lorb=1,LocalOrbs(isite)%Norb
                     !
                     call F2Bindex(LocalOrbs(isite)%Norb,[iorb,jorb],[korb,lorb],ib1,ib2)
                     Simp(isite)%N_s(iorb,jorb,ispin) = Simp(isite)%N_s(iorb,jorb,ispin) + curlyU%screened_local(ib1,ib2,1)*rho(korb,lorb,ispin)
                     !
                  enddo
               enddo
            enddo
         enddo
      enddo
      deallocate(rho)
      !
      !The magnetization will be given only by the self-energy beyond Hartree
      Simp(isite)%N_s(:,:,1) = (Simp(isite)%N_s(:,:,1)+Simp(isite)%N_s(:,:,Nspin))
      Simp(isite)%N_s(:,:,Nspin) = Simp(isite)%N_s(:,:,1)
      !
   endif
   !
   !the self-energy in the solver basis is always diagonal
   do ispin=1,Nspin
      Simp(isite)%N_s(:,:,ispin) = diag(diagonal(Simp(isite)%N_s(:,:,ispin)))
   enddo



















!do iorb=1,Lttc%Norb
!   do jorb=1,Lttc%Norb
!      do ispin=1,Nspin
!         !
!         ReS = cubic_interp( wreal_read, dreal(S_G0W0%wks(iorb,jorb,1:Nreal_read,ik,ispin)), wreal(iw) )
!         ImS = cubic_interp( wreal_read, dimag(S_G0W0%wks(iorb,jorb,1:Nreal_read,ik,ispin)), wreal(iw) )
!         S_G0W0_interp%wks(iorb,jorb,iw,ik,ispin) = dcmplx(ReS,ImS)
!         !
!         if(paramagnet)then
!            S_G0W0_interp%wks(iorb,jorb,iw,ik,Nspin) = S_G0W0_interp%wks(iorb,jorb,iw,ik,1)
!            cycle
!         endif
!         !
!      enddo
!   enddo
!enddo


!do iorb=1,Lttc%Norb
!   do jorb=1,Lttc%Norb
!      do ispin=1,Nspin
!         !
!         ReS = cubic_interp( wreal_read, dreal(S_G0W0dc%wks(iorb,jorb,1:Nreal_read,ik,ispin)), wreal(iw) )
!         ImS = cubic_interp( wreal_read, dimag(S_G0W0dc%wks(iorb,jorb,1:Nreal_read,ik,ispin)), wreal(iw) )
!         S_G0W0dc_interp%wks(iorb,jorb,iw,ik,ispin) = dcmplx(ReS,ImS)
!         !
!         if(paramagnet)then
!            S_G0W0dc_interp%wks(iorb,jorb,iw,ik,Nspin) = S_G0W0dc_interp%wks(iorb,jorb,iw,ik,1)
!            cycle
!         endif
!         !
!      enddo
!   enddo
!enddo
