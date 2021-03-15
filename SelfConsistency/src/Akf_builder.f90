program Akw_builder
   !
   use module_container
   use utils_main
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: ItStart,Itend
   !
   !
   integer                                  :: ispin,ik,iw,iorb
   integer                                  :: unit,ierr
   integer                                  :: Nreal_old
   integer                                  :: Nreal_min,Nreal_max,Nreal_read
   logical,allocatable                      :: Kmask(:)
   real(8),allocatable                      :: wreal(:),wreal_read(:)
   real(8),allocatable                      :: ImSigma_read(:,:,:),ImG_read(:,:,:)
   real(8),allocatable                      :: ImSigma(:,:,:),ReSigma(:)
   real(8),allocatable                      :: Sparams(:,:,:,:)
   complex(8),allocatable                   :: Sigma_rho(:,:),RotN(:,:,:)
   complex(8),allocatable                   :: Greal_rho(:,:),Greal_orb(:,:,:)
   real(8),allocatable                      :: Akw_rho(:,:,:),Akw_orb(:,:,:)
   real(8)                                  :: dw
   character(len=256)                       :: path
   !
   !
   !
   !---------------------------------------------------------------------------!
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   call tick(TimeStart)
   call read_InputFile("input.in")
   write(*,"(A,1I4)") "Setting Nthread:",Nthread
   call omp_set_num_threads(Nthread)
   call printHeader()
   call initialize_DataStructure(ItStart,Itend)
   call initialize_Lattice(Crystal,ItStart)
   !
   !
   !
   !---------------------------------------------------------------------------!
   !                         ADDING ADDITIONAL FOLDERS                         !
   !---------------------------------------------------------------------------!
   do ispin=1,Nspin
      if(reg(path_funct).eq."G")then
         call createDir(reg(MaxEnt_K)//"Akw_G_s"//str(ispin))
      elseif(reg(path_funct).eq."S")then
         call createDir(reg(MaxEnt_K)//"Akw_S_s"//str(ispin))
         call createDir(reg(MaxEnt_K)//"Sk_wr_s"//str(ispin))
      elseif(reg(path_funct).eq."GS")then
         call createDir(reg(MaxEnt_K)//"Akw_G_s"//str(ispin))
         call createDir(reg(MaxEnt_K)//"Akw_S_s"//str(ispin))
         call createDir(reg(MaxEnt_K)//"Sk_wr_s"//str(ispin))
      endif
   enddo
   !
   !
   !
   !---------------------------------------------------------------------------!
   !             COLLECTING RESULTS FROM MAXENT ON THE SELF-ENERGY             !
   !---------------------------------------------------------------------------!
   !
   !
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,reg(pathINPUT))
   Crystal%Nkpt_path = Crystal%Nkpt_path-1
   !
   !
   do ispin=1,Nspin
      !
      !
      if(scan(reg(path_funct),"G").gt.0)then
         !
         allocate(Kmask(Crystal%Nkpt_path));Kmask=.false.
         !
         !First check that all the files contains the same number fo real frequecies
         do ik=1,Crystal%Nkpt_path
            !
            path = reg(MaxEnt_K)//"MaxEnt_Gk_t_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT_dos.dat"
            call inquireFile(reg(path),Kmask(ik),hardstop=.false.,verb=.true.)
            if(.not.Kmask(ik)) cycle
            !
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            !
            Nreal_read=0
            ierr=0
            do while (ierr.eq.0)
               Nreal_read = Nreal_read + 1
               read(unit,*,iostat=ierr)
            enddo
            close(unit)
            !
            !MaxEnt parameters written i the last line
            Nreal_read = Nreal_read - 2
            !
            write(*,"(A,1I5)") "     The file "//reg(path)//" contains "//str(Nreal_read)//" real frequencies."
            if((ik.gt.1).and.(Nreal_read.ne.Nreal_old))then
               write(*,"(A,1I5)") "     Stop. Real frequency mesh is not consistent among K-points."
               stop
            endif
            Nreal_old=Nreal_read
            !
         enddo
         !
         allocate(wreal_read(Nreal_read));wreal=0d0
         allocate(ImG_read(Crystal%Norb,Nreal_read,Crystal%Nkpt_path));ImG_read=0d0
         do ik=1,Crystal%Nkpt_path
            !
            if(.not.Kmask(ik)) cycle
            !
            path = reg(MaxEnt_K)//"MaxEnt_Gk_t_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT_dos.dat"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            do iw=1,Nreal_read
               read(unit,*) wreal_read(iw),(ImG_read(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            close(unit)
            !
         enddo
         dw = abs(wreal_read(10)-wreal_read(9))
         write(*,"(A)") "     MaxEnt output on Green's function is read."
         !
         !Define a smaller frequency array
         if(wreal_read(Nreal_read).gt.KKcutoff)then
            Nreal_min = minloc(abs(wreal_read+KKcutoff),dim=1)
            Nreal_max = minloc(abs(wreal_read-KKcutoff),dim=1)
            Nreal = size(wreal_read(Nreal_min:Nreal_max))
            allocate(wreal(Nreal));wreal=wreal_read(Nreal_min:Nreal_max)
            write(*,"(A,F)") "     Reduced real frequency mesh (old_min): iw_["//str(Nreal_min)//"]=",wreal_read(Nreal_min)
            write(*,"(A,F)") "     Reduced real frequency mesh (old_max): iw_["//str(Nreal_max)//"]=",wreal_read(Nreal_max)
            write(*,"(A,F)") "     Reduced real frequency mesh (new_min): iw_["//str(1)//"]=",wreal(1)
            write(*,"(A,F)") "     Reduced real frequency mesh (new_max): iw_["//str(Nreal)//"]=",wreal(Nreal)
            write(*,"(A,F)") "     dw(old): "//str(dw)
            write(*,"(A,F)") "     dw(new): "//str(abs(wreal(2)-wreal(1)))
         else
            Nreal_min = 1
            Nreal_max = Nreal_read
            Nreal = Nreal_read
            allocate(wreal(Nreal));wreal=wreal_read
         endif
         deallocate(wreal_read)
         !
         !Manipulate MaxEnt output
         allocate(Akw_orb(Crystal%Norb,Nreal,Crystal%Nkpt_path));Akw_orb=0d0
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Crystal,wreal,Nreal,Nreal_max,Nreal_min,dw,ispin,KKcutoff,Akw_orb,ImG_read),&
         !$OMP PRIVATE(ik,iorb,iw)
         !$OMP DO
         do ik=1,Crystal%Nkpt_path
            do iorb=1,Crystal%Norb
               !
               !Chunk the data to the smaller frequency array and revert the poles
               do iw=1,Nreal
                  Akw_orb(iorb,iw,ik) = +abs(ImG_read(iorb,Nreal_min+(iw-1),ik))
                  if((Nreal_max-(iw-1)).lt.Nreal_min)stop "chunking issue."
               enddo
               !
               !Fix normalization
               Akw_orb(iorb,:,ik) = Akw_orb(iorb,:,ik) / abs(sum(Akw_orb(iorb,:,ik))*dw)
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(ImG_read)
         where(abs((Akw_orb))<1.d-12)Akw_orb=0d0
         write(*,"(A)") "     MaxEnt output is Normalized."
         !
         !Print
         do ik=1,Crystal%Nkpt_path
            path = reg(MaxEnt_K)//"Akw_G_s"//str(ispin)//"/Akw_orb_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
               write(unit,"(200E20.12)") wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            close(unit)
         enddo
         !
         !Write down spectral function on the path in usual format
         path = reg(MaxEnt_K)//"Akw_G_orb_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Crystal%Nkpt_path
            do iw=1,Nreal
               write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik),wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
         deallocate(Kmask,Akw_orb,wreal)
         !
      endif
      !
      !
      if(scan(reg(path_funct),"S").gt.0)then
         !
         !
         !This is only to fetch the chemical potential
         call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
         write(*,"(A)")"     Lattice chemical potential is: "//str(S_Full%mu)
         !
         allocate(Kmask(Crystal%Nkpt_path));Kmask=.false.
         !
         !First check that all the files contains the same number fo real frequecies
         do ik=1,Crystal%Nkpt_path
            !
            path = reg(MaxEnt_K)//"MaxEnt_Sk_t_s"//str(ispin)//"/Sk_t_k"//str(ik)//".DAT_dos.dat"
            call inquireFile(reg(path),Kmask(ik),hardstop=.false.,verb=.true.)
            if(.not.Kmask(ik)) cycle
            !
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            !
            Nreal_read=0
            ierr=0
            do while (ierr.eq.0)
               Nreal_read = Nreal_read + 1
               read(unit,*,iostat=ierr)
            enddo
            close(unit)
            !
            !MaxEnt parameters written i the last line
            Nreal_read = Nreal_read - 2
            !
            write(*,"(A,1I5)") "     The file "//reg(path)//" contains "//str(Nreal_read)//" real frequencies."
            if((ik.gt.1).and.(Nreal_read.ne.Nreal_old))then
               write(*,"(A,1I5)") "     Stop. Real frequency mesh is not consistent among K-points."
               stop
            endif
            Nreal_old=Nreal_read
            !
         enddo
         !
         allocate(wreal_read(Nreal_read));wreal_read=0d0
         allocate(ImSigma_read(Crystal%Norb,Nreal_read,Crystal%Nkpt_path));ImSigma_read=0d0
         do ik=1,Crystal%Nkpt_path
            !
            if(.not.Kmask(ik)) cycle
            !
            path = reg(MaxEnt_K)//"MaxEnt_Sk_t_s"//str(ispin)//"/Sk_t_k"//str(ik)//".DAT_dos.dat"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            do iw=1,Nreal_read
               read(unit,*) wreal_read(iw),(ImSigma_read(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            close(unit)
            !
            if(KKcutoff.gt.wreal_read(Nreal_read))then
               write(*,"(A)") "     Real frequency cutoff is too high. Decrease KK_CUTOFF below "//str(wreal_read(Nreal_read))//". Aborting."
               stop
            endif
            !
         enddo
         dw = abs(wreal_read(10)-wreal_read(9))
         write(*,"(A)") "     MaxEnt output on self-energy is read."
         !
         !Read the parameters for the Self-energy rescaling
         allocate(Sparams(Crystal%Norb,Crystal%Nkpt_path,Nspin,2));Sparams=0d0
         path = reg(MaxEnt_K)//"Spath_vars/Spath_Params.DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
         do ik=1,Crystal%Nkpt_path
            read(unit,*) (Sparams(iorb,ik,1,1),iorb=1,Crystal%Norb),(Sparams(iorb,ik,1,2),iorb=1,Crystal%Norb), &
                         (Sparams(iorb,ik,Nspin,1),iorb=1,Crystal%Norb),(Sparams(iorb,ik,Nspin,2),iorb=1,Crystal%Norb)
         enddo
         close(unit)
         write(*,"(A)") "     Self-energy parameters are read."
         !
         !Read the rotations of the K-dependent self-energy matrix
         allocate(RotN(Crystal%Norb,Crystal%Norb,Crystal%Nkpt_path));RotN=czero
         do ik=1,Crystal%Nkpt_path
            path = reg(MaxEnt_K)//"Spath_vars/Spath_Rot_k"//str(ik)//"_s"//str(ispin)//".DAT"
            call read_Matrix(RotN(:,:,ik),reg(path))
         enddo
         write(*,"(A)") "     Rotations are read."
         !
         !Define a smaller frequency array
         Nreal_max = minloc(abs(wreal_read-KKcutoff),dim=1)
         Nreal_min = minloc(abs(wreal_read+KKcutoff),dim=1)
         Nreal = size(wreal_read(Nreal_min:Nreal_max))
         allocate(wreal(Nreal));wreal=wreal_read(Nreal_min:Nreal_max)
         write(*,"(A,F)") "     Reduced real frequency mesh (old_min): iw_["//str(Nreal_min)//"]=",wreal_read(Nreal_min)
         write(*,"(A,F)") "     Reduced real frequency mesh (old_max): iw_["//str(Nreal_max)//"]=",wreal_read(Nreal_max)
         write(*,"(A,F)") "     Reduced real frequency mesh (new_min): iw_["//str(1)//"]=",wreal(1)
         write(*,"(A,F)") "     Reduced real frequency mesh (new_max): iw_["//str(Nreal)//"]=",wreal(Nreal)
         write(*,"(A,F)") "     dw(old): "//str(dw)
         write(*,"(A,F)") "     dw(new): "//str(abs(wreal(2)-wreal(1)))
         !
         !Manipulate MaxEnt output
         allocate(ImSigma(Crystal%Norb,Nreal,Crystal%Nkpt_path));ImSigma=0d0
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(Crystal,wreal,Nreal,Nreal_max,Nreal_min,dw,ispin,KKcutoff,ImSigma,ImSigma_read,Sparams),&
         !$OMP PRIVATE(ik,iorb,iw)
         !$OMP DO
         do ik=1,Crystal%Nkpt_path
            do iorb=1,Crystal%Norb
               !
               !Chunk the data to the smaller frequency array and revert the poles
               do iw=1,Nreal
                  ImSigma(iorb,iw,ik) = +abs(ImSigma_read(iorb,Nreal_max-(iw-1),ik))
                  if((Nreal_max-(iw-1)).lt.Nreal_min)stop "chunking issue."
               enddo
               !
               !Remove high-frequency bias
               do iw=1,Nreal
                  if(abs(wreal(iw)).gt.0.9*KKcutoff) ImSigma(iorb,iw,ik) = 0d0
               enddo
               !
               !Fix sign (negative)
               ImSigma(iorb,:,ik) = -abs(ImSigma(iorb,:,ik))
               !
               !Fix normalization
               ImSigma(iorb,:,ik) = abs(Sparams(iorb,ik,ispin,2))*pi * ImSigma(iorb,:,ik) / abs(sum(ImSigma(iorb,:,ik))*dw)
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(ImSigma_read,wreal_read)
         write(*,"(A)") "     MaxEnt output is cooked."
         !
         !Operations at each K-point
         allocate(Akw_rho(Crystal%Norb,Nreal,Crystal%Nkpt_path));Akw_rho=0d0
         allocate(Akw_orb(Crystal%Norb,Nreal,Crystal%Nkpt_path));Akw_orb=0d0
         do ik=1,Crystal%Nkpt_path
            !
            if(.not.Kmask(ik)) cycle
            !
            !Compute the real part and build the full self-energy
            allocate(Sigma_rho(Crystal%Norb,Nreal));Sigma_rho=czero
            allocate(ReSigma(Nreal));ReSigma=0d0
            do iorb=1,Crystal%Norb
               !
               call KK_Im2Re(ReSigma,ImSigma(iorb,:,ik),wreal,KKcutoff,BareVal=Sparams(iorb,ik,ispin,1),symmetric=.false.)
               !
               Sigma_rho(iorb,:) = dcmplx(ReSigma,ImSigma(iorb,:,ik))
               !
            enddo
            deallocate(ReSigma)
            write(*,"(A)") "     Real Part of the self-energy computed for K-point #"//str(ik)
            !
            !Check Print
            path = reg(MaxEnt_K)//"Sk_wr_s"//str(ispin)//"/Sk_wr_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
              write(unit,"(200E20.12)") wreal(iw),(Sigma_rho(iorb,iw),iorb=1,Crystal%Norb)
            enddo
            close(unit)
            !
            !Rebuild the Green's function in the basis where the K-dependent density matrix is diagonal
            allocate(Greal_rho(Crystal%Norb,Nreal));Greal_rho=czero
            do iorb=1,Crystal%Norb
               do iw=1,Nreal
                  Greal_rho(iorb,iw) = 1d0 / (  dcmplx(wreal(iw),eta) + S_Full%mu - Sigma_rho(iorb,iw) )
               enddo
               !
               Akw_rho(iorb,:,ik) = dimag(Greal_rho(iorb,:))
               Akw_rho(iorb,:,ik) = Akw_rho(iorb,:,ik) / (sum(Akw_rho(iorb,:,ik))*dw)
               !
            enddo
            deallocate(Sigma_rho)
            where(abs((Akw_rho(:,:,ik)))<1.d-12)Akw_rho(:,:,ik)=0d0
            !
            !Print
            path = reg(MaxEnt_K)//"Akw_S_s"//str(ispin)//"/Akw_rho_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
               write(unit,"(200E20.12)") wreal(iw),(Akw_rho(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            close(unit)
            !
            !Rotate the Green's function to the orbital basis
            allocate(Greal_orb(Crystal%Norb,Crystal%Norb,Nreal));Greal_orb=czero
            do iw=1,Nreal
               Greal_orb(:,:,iw) = rotate(diag(Greal_rho(:,iw)),transpose(conjg(RotN(:,:,ik))))
            enddo
            do iorb=1,Crystal%Norb
               Akw_orb(iorb,:,ik) = dimag(Greal_orb(iorb,iorb,:))
               Akw_orb(iorb,:,ik) = Akw_orb(iorb,:,ik) / (sum(Akw_orb(iorb,:,ik))*dw)
            enddo
            deallocate(Greal_rho,Greal_orb)
            where(abs((Akw_orb(:,:,ik)))<1.d-12)Akw_orb(:,:,ik)=0d0
            !
            !Print
            path = reg(MaxEnt_K)//"Akw_S_s"//str(ispin)//"/Akw_orb_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
               write(unit,"(200E20.12)") wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            close(unit)
            !
         enddo
         !
         !Write down spectral function on the path in usual format   Skw_s1.DAT
         path = reg(MaxEnt_K)//"Akw_S_rho_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Crystal%Nkpt_path
            do iw=1,Nreal
               write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik),wreal(iw),(Akw_rho(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
         !Write down spectral function on the path in usual format   Skw_s1.DAT
         path = reg(MaxEnt_K)//"Akw_S_orb_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Crystal%Nkpt_path
            do iw=1,Nreal
               write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik),wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
         deallocate(Kmask,wreal,Akw_orb,Akw_rho)
         !
      endif
      !
      if(paramagnet)exit
      !
   enddo
   !
   write(*,"(A)") "     Green's function on the K-path rebuilt from the self-energy. Exiting script."
   stop
   !
end program Akw_builder
