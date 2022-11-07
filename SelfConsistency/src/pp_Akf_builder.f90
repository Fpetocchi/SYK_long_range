program Akw_builder
   !
   use module_container
   use utils_main
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: ItStart,Itend
   integer                                  :: ispin
   logical                                  :: FermiSurf
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
   !                       COLLECTING RESULTS FROM MAXENT                      !
   !---------------------------------------------------------------------------!
   !
   !
   FermiSurf = print_plane_G .or. (print_plane_W.and.print_path_W) .or. (print_plane_W.and.print_path_Chi)
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(MaxEnt_K),doplane=FermiSurf,Nkpt_Kside=Nkpt_plane,hetero=Hetero)
   Crystal%Nkpt_path = Crystal%Nkpt_path !-1
   !
   !
   do ispin=1,Nspin
      !
      if(print_path_G)then
         call rebuild_G("path")
         if(Hetero%status)call rebuild_G("path",suffix="Hetero")
      endif
      !
      if(print_full_G) call rebuild_G("full")
      !
      if(print_plane_G) call rebuild_G("plane")
      !
      if(paramagnet)exit
   enddo
   !
   !
   if(print_path_Chi)then
      call rebuild_W("C","path")
      if(print_plane_W) call rebuild_W("C","plane")
   endif
   !
   if(print_path_W)then
      call rebuild_W("W","path")
      if(print_plane_W) call rebuild_W("W","plane")
   endif
   !
   !
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
      !
      use input_vars, only : Nreal, wrealMax
      use utils_misc
      implicit none
      !
      character(len=*),intent(in)           :: mode
      character(len=*),intent(in),optional  :: suffix
      integer                               :: ik,Nkpt,Nkpt_Kside,Norb
      integer                               :: iw,iorb,ikx,iky,wndx_cut
      integer                               :: unit,ierr
      integer                               :: Nreal_read,Nreal_old,ik1st
      integer                               :: Nreal_MaxEnt,Nreal_min,Nreal_max
      logical,allocatable                   :: Kmask(:)
      real(8),allocatable                   :: wreal(:),wreal_read(:),wreal_interp(:)
      real(8),allocatable                   :: Akw_orb(:,:,:),Akw_orb_interp(:,:,:)
      real(8)                               :: dw,fact
      real(8)                               :: kx,ky,Bvec(3),Blat(3,3)
      character(len=256)                    :: path,suffix_
      logical                               :: ik1st_read
      logical                               :: wr_interp=.false.
      !
      real(8),allocatable                   :: ImG_read(:,:,:)
      !
      select case(reg(mode))
         case default
            !
            stop "rebuild_G: Available Modes are: full, path, plane."
            !
         case("full")
            !
            Nkpt = Crystal%Nkpt
            Norb = Crystal%Norb
            suffix_=" "
            write(*,"(A,I)") new_line("A")//new_line("A")//"     Gfull. Total number of K-points in the BZ:",Nkpt
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
            write(*,"(2(A,I))") new_line("A")//new_line("A")//"     Gplane. Total number of K-points in the {kx,ky} sheet:",Nkpt," number of K-points per dimension:",Nkpt_Kside
            call get_Blat(Blat)
            !
      end select
      !
      !
      if(allocated(Kmask))deallocate(Kmask)
      allocate(Kmask(Nkpt));Kmask=.false.
      !
      !First check that all the files contains the same number of real frequecies
      ik1st_read=.false.
      do ik=1,Nkpt
         !
         path = reg(MaxEnt_K)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//reg(suffix_)//".DAT_dos.dat"
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
         write(*,"(A)") "     Warning: the MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//" folder is empty."
         return
      endif
      !
      allocate(wreal_read(Nreal_read));wreal_read=0d0
      allocate(ImG_read(Norb,Nreal_read,Nkpt));ImG_read=0d0
      do ik=1,Nkpt
         !
         if(.not.Kmask(ik)) cycle
         !
         path = reg(MaxEnt_K)//"MaxEnt_Gk_"//reg(mode)//"_s"//str(ispin)//"/Gk_t_k"//str(ik)//reg(suffix_)//".DAT_dos.dat"
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
         Nreal_MaxEnt = size(wreal_read(Nreal_min:Nreal_max))
         allocate(wreal(Nreal_MaxEnt));wreal=wreal_read(Nreal_min:Nreal_max)
         write(*,"(A,F)") "     Reduced real frequency mesh (old_min): iw_["//str(Nreal_min)//"]=",wreal_read(Nreal_min)
         write(*,"(A,F)") "     Reduced real frequency mesh (old_max): iw_["//str(Nreal_max)//"]=",wreal_read(Nreal_max)
         write(*,"(A,F)") "     Reduced real frequency mesh (new_min): iw_["//str(1)//"]=",wreal(1)
         write(*,"(A,F)") "     Reduced real frequency mesh (new_max): iw_["//str(Nreal_MaxEnt)//"]=",wreal(Nreal_MaxEnt)
         write(*,"(A,F)") "     dw(old): "//str(dw)
         write(*,"(A,F)") "     dw(new): "//str(abs(wreal(2)-wreal(1)))
      else
         Nreal_min = 1
         Nreal_max = Nreal_read
         Nreal_MaxEnt = Nreal_read
         allocate(wreal(Nreal_MaxEnt));wreal=wreal_read
      endif
      deallocate(wreal_read,Kmask)
      !
      !Manipulate MaxEnt output
      allocate(Akw_orb(Norb,Nreal_MaxEnt,Nkpt));Akw_orb=0d0
      !$OMP PARALLEL DEFAULT(SHARED),&
      !$OMP PRIVATE(ik,iorb,iw)
      !$OMP DO
      do ik=1,Nkpt
         do iorb=1,Norb
            !
            !Chunk the data to the smaller frequency array and revert the poles
            do iw=1,Nreal_MaxEnt
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
      !Print each K-point separately
      call createDir(reg(MaxEnt_K)//"Akw_Gk_"//reg(mode)//"_s"//str(ispin),verb=verbose)
      do ik=1,Nkpt
         path = reg(MaxEnt_K)//"Akw_Gk_"//reg(mode)//"_s"//str(ispin)//"/Akw"//reg(suffix_)//"_k"//str(ik)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nreal_MaxEnt
            write(unit,"(200E20.12)") wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Norb)
         enddo
         close(unit)
      enddo
      !
      !
      !Print in gnuplot pm3d map format
      if(reg(mode).eq."path")then
         !
         !Write down spectral function on the path in usual format - orbital basis
         path = reg(MaxEnt_K)//"Akw_Gk"//reg(suffix_)//"_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         fact=Crystal%Kpathaxis(Crystal%Nkpt_path)
         if(present(suffix).and.(reg(suffix_).eq."_Hetero"))fact=1d0
         do ik=1,Nkpt
            do iw=1,Nreal_MaxEnt
               if(abs(wreal(iw)).gt.0.5*KKcutoff)cycle
               write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik)/fact,wreal(iw),(Akw_orb(iorb,iw,ik),iorb=1,Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
      elseif(reg(mode).eq."plane")then
         !
         !Print cut at energy FermiCut in the {kx,ky} plane
         wndx_cut = minloc(abs(wreal-FermiCut),dim=1)
         write(*,"(A,F)") "     Cutting Fermi surface at w_["//str(wndx_cut)//"]=",wreal(wndx_cut)
         path = reg(MaxEnt_K)//"Fk_Gk"//reg(suffix_)//"_s"//str(ispin)//"_E"//str(FermiCut,3)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Nkpt
            ikx = int(ik/(Nkpt_Kside+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside-1) - 0.5d0
            iky = ik - (ikx-1)*Nkpt_Kside      ; ky = (iky-1)/dble(Nkpt_Kside-1) - 0.5d0
            Bvec = kx*Blat(:,1) + ky*Blat(:,2)
            write(unit,"(3I5,200E20.12)") ik,ikx,iky,Bvec(1),Bvec(2),(Akw_orb(iorb,wndx_cut,ik),iorb=1,Norb)
            if(iky.eq.Nkpt_Kside)write(unit,*)
         enddo
         close(unit)
         !
      endif
      !
      !
      !Print interpolate on the user-provided real frequency mesh
      if(wr_interp.and.(Nreal.ne.Nreal_MaxEnt).and.(reg(mode).eq."path"))then
         !
         write(*,"(A)") "     Interpolating on the input real frequency mesh."
         !
         allocate(wreal_interp(Nreal));wreal_interp=linspace(-wrealMax,+wrealMax,Nreal)
         allocate(Akw_orb_interp(Norb,Nreal,Nkpt));Akw_orb_interp=0d0
         !
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(ik,iorb,iw)
         !$OMP DO
         do ik=1,Nkpt
            do iorb=1,Norb
               do iw=1,Nreal
                  !
                  Akw_orb_interp(iorb,iw,ik) = cubic_interp( wreal, Akw_orb(iorb,:,ik), wreal_interp(iw) )
                  !
               enddo
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         !
         !Print each K-point separately
         do ik=1,Nkpt
            path = reg(MaxEnt_K)//"Akw_Gk_path_s"//str(ispin)//"/Akw"//reg(suffix_)//"_k"//str(ik)//"_interp.DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do iw=1,Nreal
               write(unit,"(200E20.12)") wreal_interp(iw),(Akw_orb_interp(iorb,iw,ik),iorb=1,Norb)
            enddo
            close(unit)
         enddo
         !
         !Print in gnuplot pm3d map format
         path = reg(MaxEnt_K)//"Akw_Gk"//reg(suffix_)//"_s"//str(ispin)//"_interp.DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         fact=Crystal%Kpathaxis(Crystal%Nkpt_path)
         if(present(suffix).and.(reg(suffix_).eq."_Hetero"))fact=1d0
         do ik=1,Nkpt
            do iw=1,Nreal
               write(unit,"(1I5,200E20.12)") ik,Crystal%Kpathaxis(ik)/fact,wreal_interp(iw),(Akw_orb_interp(iorb,iw,ik),iorb=1,Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
         deallocate(wreal_interp,Akw_orb_interp)
      endif
      deallocate(wreal,Akw_orb)
      !
   end subroutine rebuild_G
   !
   !
   !
   subroutine rebuild_W(name,mode,orbsep)
      implicit none
      character(len=*),intent(in)           :: name
      character(len=*),intent(in)           :: mode
      logical,intent(in),optional           :: orbsep
      integer                               :: iq,Nkpt,Nkpt_Kside,Norb
      integer                               :: iw,iorb,ikx,iky,wndx_cut
      integer                               :: unit,ierr
      integer                               :: Nreal_min,Nreal_max,Nreal_read,Nreal_old,ik1st
      logical,allocatable                   :: Kmask(:)
      real(8),allocatable                   :: wreal(:),wreal_read(:)
      real(8),allocatable                   :: Akw_orb(:,:,:)
      real(8)                               :: dw,fact
      real(8)                               :: kx,ky,Bvec(3),Blat(3,3)
      character(len=256)                    :: path
      logical                               :: ik1st_read,orbsep_
      !
      real(8),allocatable                   :: ImW_read(:,:,:)
      !
      select case(reg(mode))
         case default
            !
            stop "rebuild_W: Available Modes are: full, path, plane."
            !
         case("full")
            !
            Nkpt = Crystal%Nkpt
            Norb = Crystal%Norb
            write(*,"(A,I)") new_line("A")//new_line("A")//"     "//reg(name)//"full. Total number of K-points in the BZ:",Nkpt
            !
         case("path")
            !
            Nkpt = Crystal%Nkpt_path
            Norb = Crystal%Norb
            write(*,"(A,I)") new_line("A")//new_line("A")//"     "//reg(name)//"path. Total number of K-points along path:",Nkpt
            !
         case("plane")
            !
            Nkpt = Crystal%Nkpt_Plane
            Nkpt_Kside = int(sqrt(dble(Nkpt)))
            Norb = Crystal%Norb
            write(*,"(2(A,I))") new_line("A")//new_line("A")//"     "//reg(name)//"plane. Total number of K-points in the {kx,ky} sheet:",Nkpt," number of K-points per dimension:",Nkpt_Kside
            call get_Blat(Blat)
            !
      end select
      !
      ! this will be removed as soon as I find a better MaxEnt procedure
      orbsep_ = .true.
      if(present(orbsep))orbsep_=orbsep
      !
      if(orbsep_)then
         !
         do iorb=1,Norb
            !
            if(allocated(Kmask))deallocate(Kmask)
            allocate(Kmask(Nkpt));Kmask=.false.
            !
            !First check that all the files contains the same number of real frequecies
            ik1st_read=.false.
            do iq=1,Nkpt
               !
               path = reg(MaxEnt_K)//"MaxEnt_"//reg(name)//"k_"//reg(mode)//"/"//reg(name)//"k_w_k"//str(iq)//"_o"//str(iorb)//".DAT_dos.dat"
               !
               call inquireFile(reg(path),Kmask(iq),hardstop=.false.,verb=.true.)
               if(.not.Kmask(iq))then
                  cycle
               elseif(.not.ik1st_read)then
                  ik1st_read=.true.
                  ik1st = iq
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
               if(Kmask(iq).and.(iq.gt.ik1st).and.(Nreal_read.ne.Nreal_old))then
                  write(*,"(A,1I5)") "     Aborting. Real frequency mesh is not consistent among K-points."
                  return
               endif
               Nreal_old=Nreal_read
               !
            enddo
            !
            if(all(Kmask.eq..false.)) then
               write(*,"(A)") "     Warning: the MaxEnt_"//reg(name)//"k_path folder is empty for orbital "//str(iorb)
               return
            endif
            !
            allocate(wreal_read(Nreal_read));wreal_read=0d0
            allocate(ImW_read(1,Nreal_read,Nkpt));ImW_read=0d0
            !
            do iq=1,Nkpt
               !
               if(.not.Kmask(iq)) cycle
               !
               path = reg(MaxEnt_K)//"MaxEnt_"//reg(name)//"k_"//reg(mode)//"/"//reg(name)//"k_w_k"//str(iq)//"_o"//str(iorb)//".DAT_dos.dat"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
               do iw=1,Nreal_read
                  read(unit,*) wreal_read(iw),ImW_read(1,iw,iq)
               enddo
               close(unit)
               !
            enddo
            dw = abs(wreal_read(10)-wreal_read(9))
            write(*,"(A)") "     MaxEnt output on "//reg(name)//" function is read."
            !
            !Define a smaller frequency array
            if(wreal_read(Nreal_read).gt.KKcutoff)then
               Nreal_min = minloc(abs(wreal_read+KKcutoff),dim=1)
               Nreal_max = minloc(abs(wreal_read-KKcutoff),dim=1)
               Nreal = size(wreal_read(Nreal_min:Nreal_max))
               if(iorb.eq.1)allocate(wreal(Nreal))
               wreal=wreal_read(Nreal_min:Nreal_max)
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
               if(iorb.eq.1)allocate(wreal(Nreal))
               wreal=wreal_read
            endif
            deallocate(wreal_read)
            !
            !Manipulate MaxEnt output
            if(iorb.eq.1) allocate(Akw_orb(Norb,Nreal,Nkpt))
            !
            do iq=1,Nkpt
               !
               !Chunk the data to the smaller frequency array and revert the poles
               do iw=1,Nreal
                  Akw_orb(iorb,iw,iq) = ImW_read(1,Nreal_min+(iw-1),iq)
                  if((Nreal_max-(iw-1)).lt.Nreal_min)stop "chunking issue."
               enddo
               !
               !Fix normalization
               Akw_orb(iorb,:,iq) = Akw_orb(iorb,:,iq) / abs(sum(Akw_orb(iorb,:,iq))*dw)
               !
            enddo
            !
            deallocate(ImW_read)
            !
         enddo
         where(abs((Akw_orb))<1.d-12)Akw_orb=0d0
         write(*,"(A)") "     MaxEnt output is Normalized."
         !
      else
         !
         if(allocated(Kmask))deallocate(Kmask)
         allocate(Kmask(Nkpt));Kmask=.false.
         !
         !First check that all the files contains the same number of real frequecies
         ik1st_read=.false.
         do iq=1,Nkpt
            !
            path = reg(MaxEnt_K)//"MaxEnt_"//reg(name)//"k_"//reg(mode)//"/"//reg(name)//"k_w_k"//str(iq)//".DAT_dos.dat"
            !
            call inquireFile(reg(path),Kmask(iq),hardstop=.false.,verb=.true.)
            if(.not.Kmask(iq))then
               cycle
            elseif(.not.ik1st_read)then
               ik1st_read=.true.
               ik1st = iq
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
            if(Kmask(iq).and.(iq.gt.ik1st).and.(Nreal_read.ne.Nreal_old))then
               write(*,"(A,1I5)") "     Aborting. Real frequency mesh is not consistent among K-points."
               return
            endif
            Nreal_old=Nreal_read
            !
         enddo
         !
         if(all(Kmask.eq..false.)) then
            write(*,"(A)") "     Warning: the MaxEnt_"//reg(name)//"k_path folder is empty."
            return
         endif
         !
         allocate(wreal_read(Nreal_read));wreal_read=0d0
         allocate(ImW_read(Norb,Nreal_read,Nkpt));ImW_read=0d0
         do iq=1,Nkpt
            !
            if(.not.Kmask(iq)) cycle
            !
            path = reg(MaxEnt_K)//"MaxEnt_"//reg(name)//"k_"//reg(mode)//"/"//reg(name)//"k_w_k"//str(iq)//".DAT_dos.dat"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="read")
            do iw=1,Nreal_read
               read(unit,*) wreal_read(iw),(ImW_read(iorb,iw,iq),iorb=1,Norb)
            enddo
            close(unit)
            !
         enddo
         dw = abs(wreal_read(10)-wreal_read(9))
         write(*,"(A)") "     MaxEnt output on "//reg(name)//" function is read."
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
         !$OMP SHARED(Nkpt,Crystal,Norb,wreal,Nreal,Nreal_max,Nreal_min,dw,KKcutoff,Akw_orb,ImW_read),&
         !$OMP PRIVATE(iq,iorb,iw)
         !$OMP DO
         do iq=1,Nkpt
            do iorb=1,Norb
               !
               !Chunk the data to the smaller frequency array and revert the poles
               do iw=1,Nreal
                  Akw_orb(iorb,iw,iq) = ImW_read(iorb,Nreal_min+(iw-1),iq)
                  if((Nreal_max-(iw-1)).lt.Nreal_min)stop "chunking issue."
               enddo
               !
               !Fix normalization
               Akw_orb(iorb,:,iq) = Akw_orb(iorb,:,iq) / abs(sum(Akw_orb(iorb,:,iq))*dw)
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(ImW_read)
         where(abs((Akw_orb))<1.d-12)Akw_orb=0d0
         write(*,"(A)") "     MaxEnt output is Normalized."
         !
      endif
      !
      !
      !Print each K-point separately
      call createDir(reg(MaxEnt_K)//"Akw_"//reg(name)//"k_"//reg(mode),verb=verbose)
      do iq=1,Nkpt
         path = reg(MaxEnt_K)//"Akw_"//reg(name)//"k_"//reg(mode)//"/Akw_k"//str(iq)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do iw=1,Nreal
            write(unit,"(200E20.12)") wreal(iw),(Akw_orb(iorb,iw,iq),iorb=1,Norb)
         enddo
         close(unit)
      enddo
      !
      !
      !Print in gnuplot pm3d map format
      if(reg(mode).eq."path")then
         !
         path = reg(MaxEnt_K)//"Akw_"//reg(name)//"k.DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         fact=Crystal%Kpathaxis(Crystal%Nkpt_path)
         do iq=1,Nkpt
            do iw=1,Nreal
               if(abs(wreal(iw)).gt.0.5*KKcutoff)cycle
               write(unit,"(1I5,200E20.12)") iq,Crystal%Kpathaxis(iq)/fact,wreal(iw),(Akw_orb(iorb,iw,iq),iorb=1,Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         !
      elseif(reg(mode).eq."plane")then
         !
         !Print cut at energy FermiCut in the {kx,ky} plane
         wndx_cut = minloc(abs(wreal-FermiCut),dim=1)
         write(*,"(A,F)") "     Cutting Fermi surface at w_["//str(wndx_cut)//"]=",wreal(wndx_cut)
         path = reg(MaxEnt_K)//"Fk_"//reg(name)//"k_E"//str(FermiCut,3)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do iq=1,Nkpt
            ikx = int(iq/(Nkpt_Kside+0.001))+1 ; kx = (ikx-1)/dble(Nkpt_Kside-1) - 0.5d0
            iky = iq - (ikx-1)*Nkpt_Kside      ; ky = (iky-1)/dble(Nkpt_Kside-1) - 0.5d0
            Bvec = kx*Blat(:,1) + ky*Blat(:,2)
            write(unit,"(3I5,200E20.12)") iq,ikx,iky,Bvec(1),Bvec(2),(Akw_orb(iorb,wndx_cut,iq),iorb=1,Norb)
            if(iky.eq.Nkpt_Kside)write(unit,*)
         enddo
         close(unit)
         !
      endif
      deallocate(Kmask,Akw_orb,wreal)
      !
   end subroutine rebuild_W
   !
   !
   !
end program Akw_builder
