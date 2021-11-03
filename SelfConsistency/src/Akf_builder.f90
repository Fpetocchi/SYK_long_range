program Akw_builder
   !
   use module_container
   use utils_main
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: ItStart,Itend
   integer                                  :: ispin

   !
   !

   !
   !
   !
   !---------------------------------------------------------------------------!
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   call tick(TimeStart)
   call read_InputFile("input.in.akw")
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
         !
         call createDir(reg(MaxEnt_K)//"Akw_Gk_path_s"//str(ispin))  ! Spectral functions from MaxEnt on G along the K-path
         if(FermiSurf)call createDir(reg(MaxEnt_K)//"Akw_Gk_plane_s"//str(ispin)) ! Spectral functions from MaxEnt on G in the full BZ to get Fermi surface
         !
      elseif(reg(path_funct).eq."S")then
         !
         if((reg(CalculationType).eq."G0W0").or.(reg(CalculationType).eq."scGW").or.(reg(CalculationType).eq."GW+EDMFT"))then
            call createDir(reg(MaxEnt_K)//"Sk_path_wr_s"//str(ispin)) ! Rebuilt self-energy from MaxEnt on S along the K-path
            call createDir(reg(MaxEnt_K)//"Sk_full_wr_s"//str(ispin)) ! Rebuilt self-energy from MaxEnt on S in the full BZ to get Gloc
            if(FermiSurf)call createDir(reg(MaxEnt_K)//"Sk_plane_wr_s"//str(ispin)) ! Rebuilt self-energy from MaxEnt on S in the full BZ to get Fermi surface
         endif
         !
         call createDir(reg(MaxEnt_K)//"Akw_Sk_path_s"//str(ispin)) ! Spectral functions from MaxEnt on S along the K-path
         call createDir(reg(MaxEnt_K)//"Akw_Sk_full_s"//str(ispin)) ! Spectral functions from MaxEnt on S in the full BZ to get Gloc
         if(FermiSurf)call createDir(reg(MaxEnt_K)//"Akw_Sk_plane_s"//str(ispin)) ! Spectral functions from MaxEnt on S in the full BZ to get Fermi surface
         !
      elseif(reg(path_funct).eq."GS")then
         !
         call createDir(reg(MaxEnt_K)//"Akw_Gk_path_s"//str(ispin)) ! Spectral functions from MaxEnt on G along the K-path
         if(FermiSurf)call createDir(reg(MaxEnt_K)//"Akw_Gk_plane_s"//str(ispin)) ! Spectral functions from MaxEnt on G in the full BZ to get Fermi surface
         !
         if((reg(CalculationType).eq."G0W0").or.(reg(CalculationType).eq."scGW").or.(reg(CalculationType).eq."GW+EDMFT"))then
            call createDir(reg(MaxEnt_K)//"Sk_path_wr_s"//str(ispin)) ! Rebuilt self-energy from MaxEnt on S along the K-path
            call createDir(reg(MaxEnt_K)//"Sk_full_wr_s"//str(ispin)) ! Rebuilt self-energy from MaxEnt on S in the full BZ to get Gloc
            if(FermiSurf)call createDir(reg(MaxEnt_K)//"Sk_plane_wr_s"//str(ispin)) ! Rebuilt self-energy from MaxEnt on S in the full BZ to get Fermi surface
         endif
         !
         call createDir(reg(MaxEnt_K)//"Akw_Sk_path_s"//str(ispin)) ! Spectral functions from MaxEnt on S along the K-path
         call createDir(reg(MaxEnt_K)//"Akw_Sk_full_s"//str(ispin)) ! Spectral functions from MaxEnt on S in the full BZ to get Gloc
         if(FermiSurf)call createDir(reg(MaxEnt_K)//"Akw_Sk_plane_s"//str(ispin)) ! Spectral functions from MaxEnt on S in the full BZ to get Fermi surface
         !
      endif
      if(paramagnet)exit
   enddo
   !
   !
   !
   !---------------------------------------------------------------------------!
   !                       COLLECTING RESULTS FROM MAXENT                      !
   !---------------------------------------------------------------------------!
   !
   !
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,doplane=FermiSurf,Nkpt_Kside=Nkpt_Fermi,hetero=Hetero)
   Crystal%Nkpt_path = Crystal%Nkpt_path-1
   !
   !
   do ispin=1,Nspin
      !
      !Dump MaxEnt data for the local self-energy (always done for every setup)
      call rebuild_Sigma_imp("full",justSigma=.true.)
      !
      if(scan(reg(path_funct),"G").gt.0)then
         !
         !
         write(*,"(A)") new_line("A")//new_line("A")//"---- Collecting results from K-resolved MaxEnt on the Green's function."
         !
         !Collect the spectral function from MaxEnt on G
         call rebuild_G("path")
         if(FermiSurf)call rebuild_G("plane")
         if(Hetero%status)call rebuild_G("path",suffix="Hetero")
         !
      endif
      !
      if(scan(reg(path_funct),"S").gt.0)then
         !
         !
         write(*,"(A)") new_line("A")//new_line("A")//"---- Collecting results from K-resolved MaxEnt on the self-energy."
         !
         !
         !This is only to fetch the chemical potential
         call AllocateFermionicField(S_Full,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
         call read_FermionicField(S_Full,reg(ItFolder),"Sfull_w",Crystal%kpt)
         write(*,"(A)")"     Lattice chemical potential is: "//str(S_Full%mu)
         !
         select case(reg(CalculationType))
            case default
               !
               stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
               !
            case("G0W0","scGW","GW+EDMFT")
               !
               !Collect the spectral function from MaxEnt on S
               call rebuild_Sigma_K("path")
               call rebuild_Sigma_K("full")
               if(FermiSurf)call rebuild_Sigma_K("plane")
               !
            case("DMFT+statU","DMFT+dynU","EDMFT")
               !
               !Collect the spectral function from MaxEnt on Simp
               call rebuild_Sigma_imp("path")
               call rebuild_Sigma_imp("full")
               if(FermiSurf)call rebuild_Sigma_imp("plane")
               !
         end  select
         !
      endif
      !
      !
      if(paramagnet)exit
   enddo
   !!
   write(*,"(A)") "     Exiting script."
   stop
   !
   !
   !
contains
   !
   !
   !
   subroutine rebuild_G(mode,suffix)
      implicit none
      character(len=*),intent(in)           :: mode
      character(len=*),intent(in),optional  :: suffix
      integer                               :: ik,Nkpt,Nkpt_Kside,Norb
      integer                               :: iw,iorb,ikx,iky,wndx_cut
      integer                               :: unit,ierr
      integer                               :: Nreal_min,Nreal_max,Nreal_read,Nreal_old,ik1st
      logical,allocatable                   :: Kmask(:)
      real(8),allocatable                   :: wreal(:),wreal_read(:)
      real(8),allocatable                   :: Akw_orb(:,:,:)
      real(8)                               :: dw
      character(len=256)                    :: path,suffix_
      logical                               :: ik1st_read
      !
      real(8),allocatable                   :: ImG_read(:,:,:)
      !
      select case(reg(mode))
         case default
            !
            stop "rebuild_G: Available Modes are: path, full."
            !
         case("path")
            !
            Nkpt = Crystal%Nkpt_path
            Norb = Crystal%Norb
            suffix_=" "
            if(present(suffix))then
               suffix_="_"//reg(suffix)
               if(reg(suffix_).eq."_Hetero")then
                  Nkpt = Crystal%Nkpt_path + Nkpt_path
                  Norb = Hetero%Norb
               endif
            endif
            write(*,"(A,I)") new_line("A")//new_line("A")//"     Gpath. Total number of K-points along path:",Nkpt
            !
         case("plane")
            !
            Nkpt = Crystal%Nkpt_Plane
            Nkpt_Kside = int(sqrt(dble(Nkpt)))
            Norb = Crystal%Norb
            suffix_=" "
            write(*,"(2(A,I))") "     G plane. Total number of K-points in the {kx,ky} sheet:",Nkpt," number of K-points per dimension:",Nkpt_Kside
            !
      end select
      !
      !
      if(allocated(Kmask))deallocate(Kmask)
      allocate(Kmask(Nkpt));Kmask=.false.
      !
      !First check that all the files contains the same number fo real frequecies
      ik1st_read=.false.
      do ik=1,Nkpt
         !
         path = reg(MaxEnt_K)//"MaxEnt_Gk_"//reg(mode)//"_t_s"//str(ispin)//"/Gk_t_k"//str(ik)//reg(suffix_)//".DAT_dos.dat"
         !
         call inquireFile(reg(path),Kmask(ik),hardstop=.false.,verb=.true.)
         if(.not.Kmask(ik))then
            cycle
         elseif(.not.ik1st_read)then
            ik1st_read=.true.
            ik1st = ik
         endif
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
         if(Kmask(ik).and.(ik.gt.ik1st).and.(Nreal_read.ne.Nreal_old))then
            write(*,"(A,1I5)") "     Aborting. Real frequency mesh is not consistent among K-points."
            return
         endif
         Nreal_old=Nreal_read
         !
      enddo
      !
      if(all(Kmask.eq..false.)) then
         write(*,"(A)") "     Warning: the MaxEnt_Gk_"//reg(mode)//"_t_s"//str(ispin)//" folder is empty."
         return
      endif
      !
      allocate(wreal_read(Nreal_read));wreal=0d0
      allocate(ImG_read(Norb,Nreal_read,Nkpt));ImG_read=0d0
      do ik=1,Nkpt
         !
         if(.not.Kmask(ik)) cycle
         !
         path = reg(MaxEnt_K)//"MaxEnt_Gk_"//reg(mode)//"_t_s"//str(ispin)//"/Gk_t_k"//str(ik)//reg(suffix_)//".DAT_dos.dat"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
         do iw=1,Nreal_read
            read(unit,*) wreal_read(iw),(ImG_read(iorb,iw,ik),iorb=1,Norb)
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
      allocate(Akw_orb(Norb,Nreal,Nkpt));Akw_orb=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Crystal,Norb,wreal,Nreal,Nreal_max,Nreal_min,dw,ispin,KKcutoff,Akw_orb,ImG_read),&
      !$OMP PRIVATE(ik,iorb,iw)
      !$OMP DO
      do ik=1,Nkpt
         do iorb=1,Norb
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
      !
      !Print
      do ik=1,Nkpt
         path = reg(MaxEnt_K)//"Akw_Gk_"//reg(mode)//"_s"//str(ispin)//"/Akw"//reg(suffix_)//"_k"//str(ik)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nreal
            write(unit,"(200E20.12)") wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Norb)
         enddo
         close(unit)
      enddo
      !
      !
      !
      if(reg(mode).eq."path")then
         !
         !Write down spectral function on the path in usual format - orbital basis
         path = reg(MaxEnt_K)//"Akw_Gk"//reg(suffix_)//"_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Nkpt
            do iw=1,Nreal
               if(abs(wreal(iw)).gt.0.5*KKcutoff)cycle
               write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik),wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
      elseif(reg(mode).eq."plane")then
         !
         !Print cut at energy EcutSheet in the {kx,ky} plane
         wndx_cut = minloc(abs(wreal-EcutSheet),dim=1)
         path = reg(MaxEnt_K)//"Fk_Gk"//reg(suffix_)//"_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Nkpt
            ikx = int(ik/(Nkpt_Kside+0.001))+1
            iky = ik - (ikx-1)*Nkpt_Kside
            write(unit,"(3I5,200E20.12)") ik,ikx,iky,(Akw_orb(iorb,wndx_cut,ik),iorb=1,Norb)
            if(iky.eq.Nkpt_Kside)write(unit,*)
         enddo
         close(unit)
         !
      endif
      deallocate(Kmask,Akw_orb,wreal)
      !
   end subroutine rebuild_G
   !
   !
   !
   subroutine rebuild_Sigma_K(mode)
      implicit none
      character(len=*),intent(in)           :: mode
      !
      integer                               :: ik,Nkpt,Nkpt_Kside
      integer                               :: iw,iorb,ikx,iky,wndx_cut
      integer                               :: unit,ierr
      integer                               :: Nreal_min,Nreal_max,Nreal_read,Nreal_old
      logical,allocatable                   :: Kmask(:)
      real(8),allocatable                   :: wreal(:),wreal_read(:)
      real(8),allocatable                   :: Akw_orb(:,:,:)
      real(8)                               :: dw
      character(len=256)                    :: path
      !
      real(8)                               :: etafact
      real(8),allocatable                   :: ImSigma_read(:,:,:)
      real(8),allocatable                   :: ImSigma(:,:,:),ReSigma(:)
      real(8),allocatable                   :: Sparams(:,:,:,:)
      complex(8),allocatable                :: Sigma_rho(:,:),RotN(:,:,:),Srot(:,:)
      complex(8),allocatable                :: Sigma_orb(:,:,:),Greal_orb(:,:,:)
      complex(8),allocatable                :: Hk(:,:,:),Hkrot(:,:,:)
      !
      select case(reg(mode))
         case default
            !
            stop "rebuild_Sigma_K: Available Modes are: path, full, plane."
            !
         case("path")
            !
            etafact=1d0
            Nkpt = Crystal%Nkpt_path
            allocate(Hk(Crystal%Norb,Crystal%Norb,Nkpt))
            Hk = Crystal%Hk_path
            write(*,"(A,I)") new_line("A")//new_line("A")//"     Sigma path. Total number of K-points along path:",Nkpt
            write(*,"(A,F)") "     Used broadening:",eta*etafact
            !
         case("full")
            !
            etafact=1d0
            Nkpt = Crystal%Nkpt
            allocate(Hk(Crystal%Norb,Crystal%Norb,Nkpt))
            Hk = Crystal%Hk
            write(*,"(A,I)") new_line("A")//new_line("A")//"     Sigma full. Total number of K-points in the BZ:",Nkpt
            write(*,"(A,F)") "     Used broadening:",eta*etafact
            !
         case("plane")
            !
            etafact=1d0
            Nkpt = Crystal%Nkpt_Plane
            Nkpt_Kside = int(sqrt(dble(Nkpt)))
            allocate(Hk(Crystal%Norb,Crystal%Norb,Nkpt))
            Hk = Crystal%Hk_Plane
            write(*,"(2(A,I))") new_line("A")//new_line("A")//"     Sigma plane. Total number of K-points in the {kx,ky} sheet:",Nkpt," number of K-points per dimension:",Nkpt_Kside
            write(*,"(A,F)") "     Used broadening:",eta*etafact
            !
      end select
      !
      !
      if(allocated(Kmask))deallocate(Kmask)
      allocate(Kmask(Nkpt));Kmask=.false.
      !
      !First check that all the files contains the same number fo real frequecies
      do ik=1,Nkpt
         !
         path = reg(MaxEnt_K)//"MaxEnt_Sk_"//reg(mode)//"_t_s"//str(ispin)//"/Sk_t_k"//str(ik)//".DAT_dos.dat"
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
         if(Kmask(ik).and.(ik.gt.1).and.(Nreal_read.ne.Nreal_old))then
            write(*,"(A,1I5)") "     Aborting. Real frequency mesh is not consistent among K-points."
            return
         endif
         Nreal_old=Nreal_read
         !
      enddo
      !
      if(all(Kmask.eq..false.)) then
         write(*,"(A)") "     Warning: the MaxEnt_Sk_"//reg(mode)//"_t_s"//str(ispin)//" folder is empty."
         return
      endif
      !
      allocate(wreal_read(Nreal_read));wreal_read=0d0
      allocate(ImSigma_read(Crystal%Norb,Nreal_read,Nkpt));ImSigma_read=0d0
      do ik=1,Nkpt
         !
         if(.not.Kmask(ik)) cycle
         !
         path = reg(MaxEnt_K)//"MaxEnt_Sk_"//reg(mode)//"_t_s"//str(ispin)//"/Sk_t_k"//str(ik)//".DAT_dos.dat"
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
      write(*,"(A)") "     MaxEnt output on self-energy is read(K"//reg(mode)//")."
      !
      !Read the parameters for the Self-energy rescaling
      allocate(Sparams(Crystal%Norb,Nkpt,Nspin,2));Sparams=0d0
      path = reg(MaxEnt_K)//"Sigma_vars/S"//reg(mode)//"_Params.DAT"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
      do ik=1,Nkpt
         read(unit,*) (Sparams(iorb,ik,1,1),iorb=1,Crystal%Norb),(Sparams(iorb,ik,1,2),iorb=1,Crystal%Norb), &
                      (Sparams(iorb,ik,Nspin,1),iorb=1,Crystal%Norb),(Sparams(iorb,ik,Nspin,2),iorb=1,Crystal%Norb)
      enddo
      close(unit)
      write(*,"(A)") "     Self-energy parameters are read(K"//reg(mode)//")."
      !
      !Read the rotations of the K-dependent self-energy matrix
      allocate(RotN(Crystal%Norb,Crystal%Norb,Nkpt));RotN=czero
      allocate(Hkrot(Crystal%Norb,Crystal%Norb,Nkpt));Hkrot=czero
      do ik=1,Nkpt
         !
         path = reg(MaxEnt_K)//"Sigma_vars/S"//reg(mode)//"_Rot_k"//str(ik)//"_s"//str(ispin)//".DAT"
         call read_Matrix(RotN(:,:,ik),reg(path))
         !
         Hkrot(:,:,ik) = rotate(Hk(:,:,ik),RotN(:,:,ik))
         !
      enddo
      write(*,"(A)") "     Rotations are read(K"//reg(mode)//")."
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
      allocate(ImSigma(Crystal%Norb,Nreal,Nkpt));ImSigma=0d0
      !$OMP PARALLEL DEFAULT(NONE),&
      !$OMP SHARED(Nkpt,Crystal,wreal,Nreal,Nreal_max,Nreal_min,dw,ispin,KKcutoff,ImSigma,ImSigma_read,Sparams),&
      !$OMP PRIVATE(ik,iorb,iw)
      !$OMP DO
      do ik=1,Nkpt
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
               if(abs(wreal(iw)).gt.0.5*KKcutoff) ImSigma(iorb,iw,ik) = 0d0
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
      write(*,"(A)") "     MaxEnt output is cooked(K"//reg(mode)//")."
      !
      !Operations at each K-point
      allocate(Akw_orb(Crystal%Norb,Nreal,Nkpt));Akw_orb=0d0
      do ik=1,Nkpt
         !
         if(.not.Kmask(ik)) cycle
         !
         !Compute the real part and build the full self-energy
         allocate(Sigma_rho(Crystal%Norb,Nreal));Sigma_rho=czero
         allocate(ReSigma(Nreal));ReSigma=0d0
         do iorb=1,Crystal%Norb
            !
            call KK_Im2Re(ReSigma,ImSigma(iorb,:,ik),wreal,KKcutoff,BareVal=Sparams(iorb,ik,ispin,1)-dreal(Hkrot(iorb,iorb,ik)),symmetric=.false.)
            !
            Sigma_rho(iorb,:) = dcmplx(ReSigma,ImSigma(iorb,:,ik))
            !
         enddo
         deallocate(ReSigma)
         write(*,"(A)") "     Real Part of the self-energy computed for K-point #"//str(ik)
         !
         !Check Print
         path = reg(MaxEnt_K)//"Sk_"//reg(mode)//"_wr_s"//str(ispin)//"/Sk_wr_k"//str(ik)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nreal
           write(unit,"(200E20.12)") wreal(iw),(Sigma_rho(iorb,iw),iorb=1,Crystal%Norb)
         enddo
         close(unit)
         !
         !Rotate the self-energy to the orbital basis
         allocate(Sigma_orb(Crystal%Norb,Crystal%Norb,Nreal));Sigma_orb=czero
         allocate(Srot(Crystal%Norb,Crystal%Norb));Srot=czero
         do iw=1,Nreal
            Srot = diag(Sigma_rho(:,iw)) !- ( Hkrot(:,:,ik)-diag(diagonal(Hkrot(:,:,ik))) )
            Sigma_orb(:,:,iw) = rotate(Srot,transpose(conjg(RotN(:,:,ik))))
         enddo
         deallocate(Sigma_rho,Srot)
         !
         !Rebuild the Green's function in the orbital basis
         allocate(Greal_orb(Crystal%Norb,Crystal%Norb,Nreal));Greal_orb=czero
         do iw=1,Nreal
            Greal_orb(:,:,iw) = zeye(Crystal%Norb)*dcmplx(wreal(iw)+S_Full%mu,eta*etafact) - Hk(:,:,ik) - Sigma_orb(:,:,iw)
            call inv(Greal_orb(:,:,iw))
         enddo
         deallocate(Sigma_orb)
         !
         !Get the spectral function in the orbital basis
         do iorb=1,Crystal%Norb
            Akw_orb(iorb,:,ik) = dimag(Greal_orb(iorb,iorb,:))
            Akw_orb(iorb,:,ik) = Akw_orb(iorb,:,ik) / (sum(Akw_orb(iorb,:,ik))*dw)
         enddo
         deallocate(Greal_orb)
         where(abs((Akw_orb(:,:,ik)))<1.d-12)Akw_orb(:,:,ik)=0d0
         !
         !Print
         path = reg(MaxEnt_K)//"Akw_S_"//reg(mode)//"_s"//str(ispin)//"/Akw_orb_k"//str(ik)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nreal
            write(unit,"(200E20.12)") wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
         enddo
         close(unit)
         !
      enddo
      !
      !
      if(reg(mode).eq."path")then
         !
         !Write down spectral function on the path in usual format - orbital basis
         path = reg(MaxEnt_K)//"Akw_S_orb_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Nkpt
            do iw=1,Nreal
               if(abs(wreal(iw)).gt.0.5*KKcutoff)cycle
               write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik),wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
      elseif(reg(mode).eq."full")then
         !
         !Write down the local spectral function - orbital basis
         path = reg(MaxEnt_K)//"Aw_S_orb_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nreal
            write(unit,"(200E20.12)") wreal(iw),(sum(Akw_orb(iorb,iw,:)/Nkpt),iorb=1,Crystal%Norb)
         enddo
         close(unit)
         !
      elseif(reg(mode).eq."plane")then
         !
         !Print cut at energy EcutSheet in the {kx,ky} plane
         wndx_cut = minloc(abs(wreal-EcutSheet),dim=1)
         path = reg(MaxEnt_K)//"Fk_S_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Nkpt
            ikx = int(ik/(Nkpt_Kside+0.001))+1
            iky = ik - (ikx-1)*Nkpt_Kside
            write(unit,"(3I5,200E20.12)") ik,ikx,iky,(Akw_orb(iorb,wndx_cut,ik),iorb=1,Crystal%Norb)
            if(iky.eq.Nkpt_Kside)write(unit,*)
         enddo
         close(unit)
         !
      endif
      deallocate(Kmask,Sparams,RotN,Hkrot,ImSigma,wreal,Akw_orb)
      !
   end subroutine rebuild_Sigma_K
   !
   !
   !
   subroutine rebuild_Sigma_imp(mode,justSigma)
      implicit none
      character(len=*),intent(in)           :: mode
      logical,intent(in),optional           :: justSigma
      !
      integer                               :: ik,Nkpt,Nkpt_Kside
      integer                               :: iw,iorb,ikx,iky,wndx_cut
      integer                               :: unit,ierr
      integer                               :: Nreal_min,Nreal_max,Nreal_read,Nreal_old
      logical,allocatable                   :: Kmask(:)
      real(8),allocatable                   :: wreal(:),wreal_read(:)
      real(8),allocatable                   :: Akw_orb(:,:,:)
      real(8)                               :: dw
      character(len=256)                    :: path
      !
      integer                               :: isite
      real(8)                               :: etafact
      real(8),allocatable                   :: ImSigma_read(:,:)
      real(8),allocatable                   :: ImSigma(:,:),ReSigma(:)
      real(8),allocatable                   :: Sparams(:,:,:)
      complex(8),allocatable                :: Sigma_rho(:,:,:)
      complex(8),allocatable                :: Sigma_orb(:,:,:),Greal_orb(:,:,:)
      complex(8),allocatable                :: Hk(:,:,:)
      integer,allocatable                   :: Orbs(:)
      logical                               :: justSigma_
      !
      justSigma_=.false.
      if(present(justSigma))justSigma_=justSigma
      !
      select case(reg(mode))
         case default
            !
            stop "rebuild_Sigma_imp: Available Modes are: path, full, plane."
            !
         case("path")
            !
            etafact=1d0
            Nkpt = Crystal%Nkpt_path
            allocate(Hk(Crystal%Norb,Crystal%Norb,Nkpt))
            Hk = Crystal%Hk_path
            write(*,"(A,I)") new_line("A")//new_line("A")//"     Sigma path. Total number of K-points along path:",Nkpt
            write(*,"(A,F)") "     Used broadening:",eta*etafact
            !
         case("full")
            !
            etafact=1d0
            Nkpt = Crystal%Nkpt
            allocate(Hk(Crystal%Norb,Crystal%Norb,Nkpt))
            Hk = Crystal%Hk
            write(*,"(A,I)") new_line("A")//new_line("A")//"     Sigma full. Total number of K-points in the BZ:",Nkpt
            write(*,"(A,F)") "     Used broadening:",eta*etafact
            !
         case("plane")
            !
            etafact=1d0
            Nkpt = Crystal%Nkpt_Plane
            Nkpt_Kside = int(sqrt(dble(Nkpt)))
            allocate(Hk(Crystal%Norb,Crystal%Norb,Nkpt))
            Hk = Crystal%Hk_Plane
            write(*,"(2(A,I))") new_line("A")//new_line("A")//"     Sigma plane. Total number of K-points in the {kx,ky} sheet:",Nkpt," number of K-points per dimension:",Nkpt_Kside
            write(*,"(A,F)") "     Used broadening:",eta*etafact
            !
      end select
      !
      !
      do isite=1,Nsite
         !
         allocate(Orbs(SiteNorb(isite)))
         Orbs = SiteOrbs(isite,1:SiteNorb(isite))
         !
         if(allocated(Kmask))deallocate(Kmask)
         allocate(Kmask(SiteNorb(isite)));Kmask=.false.
         !
         !First check that all the files contains the same number fo real frequecies
         do iorb=1,SiteNorb(isite)
            !
            path = reg(ItFolder)//"Convergence/MaxEnt_Sqmc_"//reg(SiteName(isite))//"_o"//str(iorb)//"_s1/Sqmc_"//reg(SiteName(isite))//"_t_o"//str(iorb)//"_s"//str(ispin)//".DAT_dos.dat"
            call inquireFile(reg(path),Kmask(iorb),hardstop=.false.,verb=.true.)
            if(.not.Kmask(iorb)) then
               write(*,"(A)") "     Warning: the MaxEnt_Sqmc_"//reg(SiteName(isite))//"_o"//str(iorb)//"_s"//str(ispin)//" folder is missing one orbital."
               return
            endif
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
            if(Kmask(iorb).and.(iorb.gt.1).and.(Nreal_read.ne.Nreal_old))then
               write(*,"(A,1I5)") "     Aborting. Real frequency mesh is not consistent among orbitals."
               return
            endif
            Nreal_old=Nreal_read
            !
         enddo
         !
         allocate(wreal_read(Nreal_read));wreal_read=0d0
         allocate(ImSigma_read(SiteNorb(isite),Nreal_read));ImSigma_read=0d0
         do iorb=1,SiteNorb(isite)
            !
            path = reg(ItFolder)//"Convergence/MaxEnt_Sqmc_"//reg(SiteName(isite))//"_o"//str(iorb)//"_s1/Sqmc_"//reg(SiteName(isite))//"_t_o"//str(iorb)//"_s"//str(ispin)//".DAT_dos.dat"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            do iw=1,Nreal_read
               read(unit,*) wreal_read(iw),ImSigma_read(iorb,iw)
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
         allocate(Sparams(SiteNorb(isite),Nspin,2));Sparams=0d0
         path = reg(MaxEnt_K)//"Sigma_vars/Sqmc_"//reg(SiteName(isite))//"_Params.DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
         read(unit,*) (Sparams(iorb,1,1),iorb=1,SiteNorb(isite)),(Sparams(iorb,1,2),iorb=1,SiteNorb(isite)), &
                      (Sparams(iorb,Nspin,1),iorb=1,SiteNorb(isite)),(Sparams(iorb,Nspin,2),iorb=1,SiteNorb(isite))
         close(unit)
         write(*,"(A)") "     Self-energy parameters are read."
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
         if(isite.eq.1)then
            allocate(Akw_orb(Crystal%Norb,Nreal,Nkpt));Akw_orb=0d0
            allocate(Sigma_orb(Crystal%Norb,Crystal%Norb,Nreal));Sigma_orb=czero
         endif
         !
         !Manipulate MaxEnt output
         allocate(ImSigma(SiteNorb(isite),Nreal));ImSigma=0d0
         do iorb=1,SiteNorb(isite)
            !
            !Chunk the data to the smaller frequency array and revert the poles
            do iw=1,Nreal
               ImSigma(iorb,iw) = +abs(ImSigma_read(iorb,Nreal_max-(iw-1)))
               if((Nreal_max-(iw-1)).lt.Nreal_min)stop "chunking issue."
            enddo
            !
            !Remove high-frequency bias
            do iw=1,Nreal
               if(abs(wreal(iw)).gt.0.5*KKcutoff) ImSigma(iorb,iw) = 0d0
            enddo
            !
            !Fix sign (negative)
            ImSigma(iorb,:) = -abs(ImSigma(iorb,:))
            !
            !Fix normalization
            ImSigma(iorb,:) = abs(Sparams(iorb,ispin,2))*pi * ImSigma(iorb,:) / abs(sum(ImSigma(iorb,:))*dw)
            !
         enddo
         deallocate(ImSigma_read,wreal_read)
         write(*,"(A)") "     MaxEnt output is cooked."
         !
         !
         !Compute the real part and build the full self-energy
         allocate(Sigma_rho(SiteNorb(isite),SiteNorb(isite),Nreal));Sigma_rho=czero
         allocate(ReSigma(Nreal));ReSigma=0d0
         do iorb=1,SiteNorb(isite)
            !
            call KK_Im2Re(ReSigma,ImSigma(iorb,:),wreal,KKcutoff,BareVal=Sparams(iorb,ispin,1),symmetric=.false.)
            !
            Sigma_rho(iorb,iorb,:) = dcmplx(ReSigma,ImSigma(iorb,:))
            !
         enddo
         deallocate(ReSigma,ImSigma)
         write(*,"(A)") "     Real Part of the self-energy computed."
         !
         !Check Print
         do iorb=1,SiteNorb(isite)
            path = reg(ItFolder)//"Convergence/MaxEnt_Sqmc_"//reg(SiteName(isite))//"_o"//str(iorb)//"_s1/Sqmc_rebuilt_"//reg(SiteName(isite))//"_t_o"//str(iorb)//"_s"//str(ispin)//".DAT_dos.dat"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
              write(unit,"(200E20.12)") wreal(iw),dreal(Sigma_rho(iorb,iorb,iw)),dimag(Sigma_rho(iorb,iorb,iw))
            enddo
            close(unit)
         enddo
         !
         !Expand to the Lattice basis
         if(RotateHloc)then
            do iw=1,Nreal
               call imp2loc(Sigma_orb(:,:,iw),Sigma_rho(:,:,iw),isite,Orbs,ExpandImpurity,U=OlocRotDag)
            enddo
         else
            do iw=1,Nreal
               call imp2loc(Sigma_orb(:,:,iw),Sigma_rho(:,:,iw),isite,Orbs,ExpandImpurity)
            enddo
         endif
         !
         deallocate(Orbs)
         if(ExpandImpurity.or.AFMselfcons)exit
         !
      enddo
      !
      !Operations at each K-point
      if(.not.justSigma_)then
         !
         do ik=1,Nkpt
            !
            !Rebuild the Green's function in the orbital basis
            allocate(Greal_orb(Crystal%Norb,Crystal%Norb,Nreal));Greal_orb=czero
            do iw=1,Nreal
               Greal_orb(:,:,iw) = zeye(Crystal%Norb)*dcmplx(wreal(iw)+S_Full%mu,eta*etafact) - Hk(:,:,ik) - Sigma_orb(:,:,iw)
               call inv(Greal_orb(:,:,iw))
            enddo
            !
            !Get the spectral function in the orbital basis
            do iorb=1,Crystal%Norb
               Akw_orb(iorb,:,ik) = dimag(Greal_orb(iorb,iorb,:))
               Akw_orb(iorb,:,ik) = Akw_orb(iorb,:,ik) / (sum(Akw_orb(iorb,:,ik))*dw)
            enddo
            deallocate(Greal_orb)
            where(abs((Akw_orb(:,:,ik)))<1.d-12)Akw_orb(:,:,ik)=0d0
            !
            !Print
            path = reg(MaxEnt_K)//"Akw_S_"//reg(mode)//"_s"//str(ispin)//"/Akw_orb_k"//str(ik)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
               write(unit,"(200E20.12)") wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
            enddo
            close(unit)
            !
         enddo
         deallocate(Sigma_orb)
         !
         !
         if(reg(mode).eq."path")then
            !
            !Write down spectral function on the path in usual format - orbital basis
            path = reg(MaxEnt_K)//"Akw_S_orb_s"//str(ispin)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,Nkpt
               do iw=1,Nreal
                  if(abs(wreal(iw)).gt.0.5*KKcutoff)cycle
                  write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik),wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Crystal%Norb)
               enddo
               write(unit,*)
            enddo
            close(unit)
            !
         elseif(reg(mode).eq."full")then
            !
            !Write down the local spectral function - orbital basis
            path = reg(MaxEnt_K)//"Aw_S_orb_s"//str(ispin)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
               write(unit,"(200E20.12)") wreal(iw),(sum(Akw_orb(iorb,iw,:)/Nkpt),iorb=1,Crystal%Norb)
            enddo
            close(unit)
            !
         elseif(reg(mode).eq."plane")then
            !
            !Print cut at energy EcutSheet in the {kx,ky} plane
            wndx_cut = minloc(abs(wreal-EcutSheet),dim=1)
            path = reg(MaxEnt_K)//"Fk_S_s"//str(ispin)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,Nkpt
               ikx = int(ik/(Nkpt_Kside+0.001))+1
               iky = ik - (ikx-1)*Nkpt_Kside
               write(unit,"(3I5,200E20.12)") ik,ikx,iky,(Akw_orb(iorb,wndx_cut,ik),iorb=1,Crystal%Norb)
               if(iky.eq.Nkpt_Kside)write(unit,*)
            enddo
            close(unit)
            !
         endif
         !
      endif
      deallocate(Kmask,Sparams,wreal,Akw_orb)
      !
   end subroutine rebuild_Sigma_imp
   !
   !
   !
end program Akw_builder
