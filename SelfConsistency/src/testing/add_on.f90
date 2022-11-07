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
