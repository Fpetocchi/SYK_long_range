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












   !
   subroutine read_U_respack_full_TESTING(Umats,LocalOnly,Lttc,pathOUTPUT)
      !
      use parameters
      use file_io
      use utils_misc
      use utils_fields
      use crystal
      use greens_function, only : set_density, calc_Gmats
      use bubbles, only : calc_PiGG
      use linalg, only : zeye, inv
      use bubbles, only : calc_PiGG
      use input_vars, only : pathINPUT, U_AC, RespackRthresh
      use input_vars, only : Nkpt_path, structure, look4dens
      implicit none
      !
      type(BosonicField),intent(inout)      :: Umats
      logical,intent(in)                    :: LocalOnly
      type(Lattice),intent(inout)           :: Lttc
      character(len=*),intent(in),optional  :: pathOUTPUT
      !
      logical                               :: filexists,exitRloop,FTdone
      character(len=256)                    :: file_respack,pathOUTPUT_
      integer                               :: unit,NRW,NRJ
      integer                               :: iR,iw,iRread
      integer                               :: Nspin_respack,Norb_respack,Nfreq,Nmats
      integer                               :: ib1,ib2,Nbp,iorb,jorb,Norb,iq
      integer                               :: NwigMat,iwigMat,iwig,nx,ny,nz
      real(8),allocatable                   :: wread(:)
      complex(8),allocatable                :: invW(:,:),Rmat(:,:)
      complex(8),allocatable                :: Vr(:,:,:),Ur(:,:,:,:)
      type(BosonicField)                    :: Wreal,Wmats,Pmats
      type(FermionicField)                  :: Gmats
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- read_U_respack_full"
      pathOUTPUT_ = pathINPUT
      if(present(pathOUTPUT)) pathOUTPUT_ = pathOUTPUT
      !
      ! Check on the input Boson
      if(.not.Umats%status) stop "read_U_respack_full: BosonicField not properly initialized."
      if((.not.LocalOnly).and.(.not.allocated(Umats%bare))) stop "read_U_respack_full: Requested k-dependence but bare non-local attribute not allocated."
      if((.not.LocalOnly).and.(.not.allocated(Umats%screened))) stop "read_U_respack_full: Requested k-dependence but screened non-local attribute not allocated."
      if(LocalOnly.and.allocated(Umats%bare)) stop "read_U_respack_full: Bare K-dependent attributes is present but not used."
      if(LocalOnly.and.allocated(Umats%screened)) stop "read_U_respack_full: Screened K-dependent attributes is present but not used."
      if(.not.Lttc%status) stop "read_U_respack_full: Lattice not properly initialized."
      if(Umats%Nbp.ne.(Lttc%Norb**2)) stop "read_U_respack_full: Umats has different orbital dimension with respect to Lttc."
      !
      Norb = Lttc%Norb
      Nbp = Umats%Nbp
      Nmats = Umats%Npoints
      !
      call inquireFile(reg(pathOUTPUT_)//"VW_imag/VW.Q0001.DAT",FTdone,hardstop=.false.,verb=verbose)
      if(.not.FTdone)then
         !
         !recover the vectors in real space
         if(.not.Wig_stored)call calc_wignerseiz(Lttc%Nkpt3)
         !
         ! Look for the Number of respack files. Which are supposed to be ordered.
         NRW = 0
         do iR=1,9999
            file_respack = reg(pathINPUT)//"Respack/dir-intW_R/VW.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) exit
            NRW = NRW + 1
         enddo
         if(NRW.eq.0) stop "read_U_respack_full: no dir-intW_R/VW.R* file found."
         write(*,"(A,I)") "     The number of RESPACK files (NRW) in dir-intW_R is: ",NRW
         !
         NRJ = 0
         do iR=1,9999
            file_respack = reg(pathINPUT)//"Respack/dir-intJ_R/XJ.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,hardstop=.false.,verb=verbose)
            if(.not.filexists) exit
            NRJ = NRJ + 1
         enddo
         write(*,"(A,I)") "     The number of RESPACK files (NRJ) in dir-intJ_R is: ",NRJ
         !
         !
         NwigMat = Nwig
         if(LocalOnly) NwigMat = 1
         !
         ! Read respack files to product basis form
         allocate(Rmat(Norb,Norb));Rmat=czero
         allocate(Vr(Nbp,Nbp,NwigMat));Vr=czero
         !
         !
         ! Read density-density interaction
         exitRloop=.false.
         intWRloop: do iR=1,NRW
            !
            ! File containing density-density components
            file_respack = reg(pathINPUT)//"Respack/dir-intW_R/VW.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,verb=verbose) !redundant control
            unit = free_unit()
            open(unit,file=reg(file_respack),form="unformatted",action="read")
            !
            ! Get wigner-seitz vector
            read(unit) nx,ny,nz
            iwig = find_vec([nx,ny,nz],Nvecwig,hardstop=.false.)
            !
            ! Allocate real-frequency dependent arrays
            read(unit) iRread,Nspin_respack,Norb_respack,Nfreq
            if(iR.ne.iRread) stop "read_U_respack_full: dir-intW_R, iR.ne.iRread"
            if(Norb.ne.Norb_respack) stop "read_U_respack_full: dir-intW_R, Norb.ne.Norb_respack"
            if(iR.eq.1)then
               allocate(Ur(Nbp,Nbp,Nfreq,NwigMat));Ur=czero
               allocate(wread(Nfreq));wread=0d0
            endif
            read(unit) wread
            !
            ! Cycle conditions
            iwigMat = iwig
            if(iwig.eq.0)then
               if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"], not available."
               cycle intWRloop
            endif
            !
            if(LocalOnly)then
               if(iwig.eq.wig0)then
                  if(any([nx,ny,nz].ne.[0,0,0])) stop "read_U_respack_full: wrong index of R=0 vector (intW_R)." !redundant
                  iwigMat = 1
                  exitRloop = .true.
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, found W(R=0) in file: "//reg(file_respack)
               else
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, skipping file: "//reg(file_respack)
                  cycle intWRloop
               endif
            else
               if(radiuswig(iwigMat).gt.RespackRthresh)then
                  if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" is beyond radius threshold ("//str(radiuswig(iwigMat))//")"
                  cycle intWRloop
               else
                  write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"] iwig:"//str(iwigMat)
               endif
            endif
            !
            ! Read data from file
            do iw=0,Nfreq
               !
               Rmat=czero
               read(unit) Rmat
               !
               if(iw.eq.0) then
                  ! Fill in bare limit
                  do iorb=1,Norb
                     do jorb=1,Norb
                        call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                        Vr(ib1,ib2,iwigMat) = Rmat(iorb,jorb)
                     enddo
                  enddo
               else
                  ! Fill in screened part
                  do iorb=1,Norb
                     do jorb=1,Norb
                        call F2Bindex(Norb,[iorb,iorb],[jorb,jorb],ib1,ib2)
                        Ur(ib1,ib2,iw,iwigMat) = Rmat(iorb,jorb)
                     enddo
                  enddo
               endif
               !
            enddo
            !
            !TEST>>>
            if(iwigMat.eq.wig0)then
               call dump_Field_component(Ur(1,1,:,iwigMat),reg(pathOUTPUT_),"Ur_11_R0.DAT",wread)
            else
               call dump_Field_component(Ur(1,1,:,iwigMat),reg(pathOUTPUT_),"Ur_11_R"//str(iwigMat)//".DAT",wread)
            endif
            !>>>TEST
            !
            if(exitRloop)then
               write(*,"(A)")"     Density-density interaction matrix read."
               exit intWRloop
            endif
            !
         enddo intWRloop
         !
         !
         ! Read exchange interaction
         exitRloop=.false.
         intJRloop: do iR=1,NRJ
            !
            ! File containing exchange components
            file_respack = reg(pathINPUT)//"Respack/dir-intJ_R/XJ.R."//str(iR,4)
            call inquireFile(reg(file_respack),filexists,verb=verbose) !redundant control
            unit = free_unit()
            open(unit,file=reg(file_respack),form="unformatted",action="read")
            !
            ! Fet wigner-seitz vector
            read(unit) nx,ny,nz
            iwig = find_vec([nx,ny,nz],Nvecwig,hardstop=.false.)
            !
            ! Allocate real-frequency dependent arrays
            read(unit) iRread,Nspin_respack,Norb_respack,Nfreq
            if(iR.ne.iRread) stop "read_U_respack_full: dir-intJ_R, iR.ne.iRread"
            if(Norb.ne.Norb_respack) stop "read_U_respack_full: dir-intJ_R, Norb.ne.Norb_respack"
            read(unit) wread
            !
            ! Cycle conditions
            iwigMat = iwig
            if(iwig.eq.0)then
               if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"], not available."
               cycle intJRloop
            endif
            !
            if(LocalOnly)then
               if(iwig.eq.wig0)then
                  if(any([nx,ny,nz].ne.[0,0,0])) stop "read_U_respack_full: wrong index of R=0 vector (intJ_R)." !redundant
                  iwigMat = 1
                  exitRloop = .true.
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, found J(R=0) in file: "//reg(file_respack)
               else
                  write(*,"(A)")"     read_U_respack_full in LocalOnly mode, skipping file: "//reg(file_respack)
                  cycle intJRloop
               endif
            else
               if(radiuswig(iwigMat).gt.RespackRthresh)then
                  if(verbose)write(*,"(A)")"     file: "//reg(file_respack)//" is beyond radius threshold ("//str(radiuswig(iwigMat))//")"
                  cycle intJRloop
               else
                  write(*,"(A)")"     file: "//reg(file_respack)//" correspond to: ["//str(nx)//","//str(ny)//","//str(nz)//"] iwig:"//str(iwigMat)
               endif
            endif
            !
            ! Read data from file
            do iw=0,Nfreq
               !
               Rmat=czero
               read(unit) Rmat
               !
               if(iw.eq.0) then
                  ! Fill in bare limit
                  do iorb=1,Norb
                     do jorb=1,Norb
                        if(iorb.ne.jorb)then
                           call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                           Vr(ib1,ib2,iwigMat) = Rmat(iorb,jorb)
                           call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                           Vr(ib1,ib2,iwigMat) = Rmat(iorb,jorb)
                        endif
                     enddo
                  enddo
               else
                  ! Fill in screened part
                  do iorb=1,Norb
                     do jorb=1,Norb
                        if(iorb.ne.jorb)then
                           call F2Bindex(Norb,[iorb,jorb],[iorb,jorb],ib1,ib2)
                           Ur(ib1,ib2,iw,iwigMat) = Rmat(iorb,jorb)
                           call F2Bindex(Norb,[iorb,jorb],[jorb,iorb],ib1,ib2)
                           Ur(ib1,ib2,iw,iwigMat) = Rmat(iorb,jorb)
                        endif
                     enddo
                  enddo
               endif
               !
            enddo
            !
            if(exitRloop)then
               write(*,"(A)")"     Exchange interaction matrix read."
               exit intJRloop
            endif
            !
         enddo intJRloop
         deallocate(Rmat)
         !
         ! FT transform to momentum space
         if(LocalOnly)then
            !
            call AllocateBosonicField(Wreal,Norb,Nfreq,1)
            !
            Wreal%screened = Ur
            Wreal%bare = Vr
            !
         else
            !
            call AllocateBosonicField(Wreal,Norb,Nfreq,1,Nkpt=Lttc%Nkpt)
            !
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Vr,Wreal%bare)
            call wannier_R2K(Lttc%Nkpt3,Lttc%kpt,Ur,Wreal%screened)
            !
         endif
         !
         deallocate(Vr,Ur)
         !
         call BosonicKsum(Wreal)
         !
         Umats%bare = Wreal%bare
         Umats%bare_local = Wreal%bare_local
         !
         ! Store fully screened interaction on real-frequency axis
         call dump_BosonicField(Wreal,reg(pathINPUT)//"Respack/","Wloc_real.DAT",axis=wread)
         call dump_BosonicField(Wreal,reg(pathINPUT)//"Respack/W_real_readable/",.false.,axis=wread) ! to be put if(verbose) in front
         ! This is going to be deleted later on because of wrong convention name
         call dump_BosonicField(Wreal,reg(pathINPUT)//"VW_real/",.true.,axis=wread)
         call DeallocateField(Wreal)
         deallocate(wread)
         !
         ! Calling spex subroutine to do the analytic continuation on Wreal written in the VW_real folder
         call duplicate(Wmats,Umats)
         call read_U_spex(Wmats,save2readable=.false.,kpt=Lttc%kpt,doAC=U_AC,pathOUTPUT=reg(pathOUTPUT_),HartreeData=.false.)
         call dump_BosonicField(Wmats,reg(pathOUTPUT_),"Wloc_mats.DAT")
         !
         ! Auxiliray files printed by read_U_spex corresponding to the fully screened interaction
         call removeDir(reg(pathINPUT)//"VW_real")
         call removeFile(reg(pathOUTPUT_)//"Uloc_real.DAT")
         call removeFile(reg(pathOUTPUT_)//"Uloc_mats.DAT")
         call removeFile(reg(pathOUTPUT_)//"U_*",list=.true.)
         call removeFile(reg(pathOUTPUT_)//"Unn_*",list=.true.)
         !
         ! Compute the non-interacting polarization of the model
         call AllocateFermionicField(Gmats,Lttc%Norb,Nmats,Nkpt=Lttc%Nkpt,Beta=Umats%Beta)
         Gmats%mu = look4dens%mu
         call set_density(Gmats%mu,Umats%Beta,Lttc,look4dens)
         call calc_Gmats(Gmats,Lttc)
         !
         call AllocateBosonicField(Pmats,Norb,Nmats,1,Nkpt=Lttc%Nkpt,no_bare=.true.,Beta=Umats%Beta)
         call calc_PiGG(Pmats,Gmats,Lttc)
         call dump_BosonicField(Pmats,reg(pathOUTPUT_),"Ploc_mats.DAT")
         call DeallocateField(Gmats)
         !
         ! Upscreen the model
         allocate(invW(Nbp,Nbp));invW=czero
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(iw,iq,invW)
         !$OMP DO
         do iw=1,Umats%Npoints
            do iq=1,Umats%Nkpt
               !
               ! [ 1 + W*Pi ]
               invW = zeye(Umats%Nbp) + matmul(Wmats%screened(:,:,iw,iq),Pmats%screened(:,:,iw,iq))
               !
               ! [ 1 + W*Pi ]^-1
               call inv(invW)
               !
               ! [ 1  W*Pi ]^-1*W
               Umats%screened(:,:,iw,iq) = matmul(invW,Wmats%screened(:,:,iw,iq))
               !
            enddo
            !
         enddo !iw
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(invW)
         call DeallocateField(Wmats)
         call DeallocateField(Pmats)
         !
         call BosonicKsum(Umats)
         !
         ! Overwrite what's printed by spex previously
         call dump_BosonicField(Umats,reg(pathOUTPUT_),"Uloc_mats_NEW.DAT")
         call dump_BosonicField(Umats,reg(pathOUTPUT_)//"VW_imag/",.true.)
         call dump_BosonicField(Umats,reg(pathOUTPUT_)//"VW_imag_readable/",.false.) ! to be put if(verbose) in front
         !
         write(*,"(A)")"     Momentum-dependent interaction on the real frequency axis written to file (SPEX format)."
         !
      endif
      !
      !print along path
      if(reg(structure).ne."None")then
         call interpolate2Path(Lttc,Nkpt_path,"Uk_wm",pathOUTPUT=reg(pathINPUT),store=.false.,skipAkw=.true.,data_in=Umats%screened(:,:,1,:))
      endif
      !
   end subroutine read_U_respack_full_TESTING





   subroutine solve_Hk_along_BZpath(hk_model,Nlso,kpath,Nkpath,colors_name,points_name,file,iproject)
     interface
        function hk_model(kpoint,N)
          real(8),dimension(:)      :: kpoint
          integer                   :: N
          complex(8),dimension(N,N) :: hk_model
        end function hk_model
     end interface
     integer                                   :: Nlso
     real(8),dimension(:,:)                    :: kpath
     integer                                   :: Nkpath
     type(rgb_color),dimension(Nlso)           :: colors_name
     character(len=*),dimension(size(kpath,1)) :: points_name
     character(len=*),optional                 :: file
     logical,optional                          :: iproject
     character(len=256)                        :: file_
     logical                                   :: iproject_
     character(len=256)                        :: xtics
     integer                                   :: Npts,Ndim,Nktot
     integer                                   :: ipts,ik,ic,unit,iorb
     real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff
     real(8)                                   :: eval(Nlso),coeff(Nlso),klen,ktics(size(Kpath,1))
     complex(8)                                :: h(Nlso,Nlso)
     type(rgb_color)                           :: corb(Nlso),c(Nlso)
     character(len=10)                         :: chpoint
     character(len=32)                         :: fmt
     real(8),allocatable                       :: kseg(:),Ekval(:,:)
     integer,allocatable                       :: Ekcol(:,:)
     !
     file_    = "Eigenbands.tb";if(present(file))file_=file
     iproject_= TB_w90%status
     if(TB_w90%status)write(*,*)"Using iproject=.TRUE. in W90 interface. Disable it explicitly using iproject=.false. "
     if(present(iproject))iproject_=iproject
     !
     Npts = size(kpath,1)
     Ndim = size(kpath,2)
     Nktot= (Npts-1)*Nkpath

!variable definition
type(rgb_color) :: corb(Nlso),c(Nlso)
!here ou initialize corb to some inital color different for each iorb
do iorb=1,Nlso
   corb(iorb) = colors_name(iorb)
enddo


allocate(Ekval(Nkpath,Norb)) !real
allocate(Ekcol(Nkpath,Norb)) !integers
!loop over the points along the path
do ik=1,Nkpath
   Zk = Hk(:,:,ik)
   !diagonalize Hamiltonian at each kpoint
   call eigh(Zk,Ek)
   do iorb=1,Norb
      coeff(:)=Zk(:,iorb)*conjg(Zk(:,iorb)) !<-- coeff(:) is a real vector of dimension Norb
      c(iorb) = coeff.dot.corb              !<-- c(:) is a new rgb vector type that renormalizes the initial corb via the coeff vector
      Ekval(ik,iorb) = Eval(iorb)           !<-- hamiltonian eigenvalue
      !gnucol is an integer number interpreted as a color by gnuplot
      gnucol = int(c(iorb)%r)*65536 + int(c(iorb)%g)*256 + int(c(iorb)%b)
      Ekcol(ik,iorb) = gnucol               !<-- 2D array of integers
   enddo
enddo

!write in Norb blocks of length Nkpath
do iorb=1,Norb
  do ik=1,Nkpath
     write(unit,*)kpoint(ik),Ekval(ik,iorb),Ekcol(ik,iorb)
  enddo
  write(unit,*)""
enddo

function rgb(c) result(num)
  type(rgb_color),intent(in) :: c
  integer :: num
  num = int(c%r)*65536 + int(c%g)*256 + int(c%b)
end function rgb

     !
     if(.not.set_bkvec)stop "solve_w90hk_along_BZpath ERROR: bk vectors not set!"
     !


     if(iproject_)then
        select case(Ndim)
        case (1)
           forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
        case(2)
           forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
        case (3)
           forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
        end select
     endif
     !
     !
     write(*,*)"Solving model along the path:"
     write(fmt,"(A3,I0,A)")"(A,",size(kpath,2),"F7.4,A1)"
     do ipts=1,Npts
        write(*,fmt)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
     enddo
     !
     ic = 0
     allocate(kseg(Nktot))
     allocate(ekval(Nktot,Nlso))
     allocate(ekcol(Nktot,Nlso))
     klen=0d0
     do ipts=1,Npts-1
        kstart = kpath(ipts,:)
        kstop  = kpath(ipts+1,:)
        kdiff  = (kstop-kstart)/dble(Nkpath)
        ktics(ipts)  = klen

     enddo
     ktics(Npts) = kseg(ic-1)
     !
     unit=free_unit()
     open(unit,file=reg(file_))
     do iorb=1,Nlso
        do ic=1,Nktot
           write(unit,*)kseg(ic),Ekval(ic,iorb),Ekcol(ic,iorb)
        enddo
        write(unit,*)""
     enddo
     close(unit)
     !
     xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
     do ipts=2,Npts-1
        xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
     enddo
     xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
     !
     open(unit,file=reg(file_)//".gp")
     write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
     write(unit,*)"#set out '"//reg(file_)//".png'"
     write(unit,*)""
     write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
     write(unit,*)"#set out '"//reg(file_)//".svg'"
     write(unit,*)""
     write(unit,*)"#set term postscript eps enhanced color 'Times'"
     write(unit,*)"#set output '|ps2pdf -dEPSCrop - "//reg(file_)//".pdf'"
     write(unit,*)"unset key"
     write(unit,*)"set xtics ("//reg(xtics)//")"
     write(unit,*)"set grid noytics xtics"
     !
     do iorb=1,Nlso
        chpoint=str(0.95d0-(iorb-1)*0.05d0)
        write(unit,"(A)")str("#set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
             " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
     enddo
     !
     write(unit,*)"plot '"//reg(file_)//"' every :::0 u 1:2:3 w l lw 3 lc rgb variable"
     write(unit,*)"# to print from the i-th to the j-th block use: every :::i::j"
     !
     close(unit)
     call system("chmod +x "//reg(file_)//".gp")
   end subroutine solve_Hk_along_BZpath






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





   if(dump_unfolded_v0)then
      !
      call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
      call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
      write(*,*)"Gfull(K) computed"
      !
      allocate(Gfull_R(Norb,Norb,Nmats,Nwig,2));Gfull_R=czero
      call wannier_K2R(Lttc%Nkpt3,Lttc%kpt,Gfull%wks(:,:,:,:,1),Gfull_R(:,:,:,:,1))
      Gfull%wks(1,2,:,:,:) = -Gfull%wks(1,2,:,:,:)
      Gfull%wks(2,1,:,:,:) = -Gfull%wks(2,1,:,:,:)
      call wannier_K2R(Lttc%Nkpt3,Lttc%kpt,Gfull%wks(:,:,:,:,1),Gfull_R(:,:,:,:,2))
      call DeallocateField(Gfull)
      write(*,*)"Gfull_AB(R) computed"
      !
      allocate(Gfull_R_halfkz(3*Norb,3*Norb,Nmats,Nwig,2));Gfull_R_halfkz=czero
      !
      do nx=minval(Nvecwig(1,:)),maxval(Nvecwig(1,:)),1
         do ny=minval(Nvecwig(2,:)),maxval(Nvecwig(2,:)),1
            !
            !new iwig
            iwig = find_vec([nx,ny,0],Nvecwig,hardstop=.false.)
            if(iwig.eq.0)cycle
            !
            do nz=minval(Nvecwig(3,:)),maxval(Nvecwig(3,:)),1
               !
               if(nz.eq.0)then
                  !diagonals
                  do ik=1,3
                     !B
                     Gfull_R_halfkz(1+(ik-1)*Norb:ik*Norb,1+(ik-1)*Norb:ik*Norb,:,iwig,1) = Gfull_R(:,:,:,iwig,1)
                     !AB
                     Gfull_R_halfkz(1+(ik-1)*Norb:ik*Norb,1+(ik-1)*Norb:ik*Norb,:,iwig,2) = Gfull_R(:,:,:,iwig,2)
                  enddo
               elseif(nz.eq.+1)then
                  !distance +1
                  iwig_old = find_vec([nx,ny,+1],Nvecwig,hardstop=.false.)
                  if(iwig_old.eq.0)cycle
                  !
                  !B
                  Gfull_R_halfkz(1+Norb:2*Norb,1:Norb,:,iwig,1) = Gfull_R(:,:,:,iwig_old,1)
                  Gfull_R_halfkz(1+2*Norb:3*Norb,1+Norb:2*Norb,:,iwig,1) = Gfull_R(:,:,:,iwig_old,1)
                  !AB
                  Gfull_R_halfkz(1+Norb:2*Norb,1:Norb,:,iwig,2) = Gfull_R(:,:,:,iwig_old,2)
                  Gfull_R_halfkz(1+2*Norb:3*Norb,1+Norb:2*Norb,:,iwig,2) = Gfull_R(:,:,:,iwig_old,2)
               elseif(nz.eq.-1)then
                  !distance -1
                  iwig_old = find_vec([nx,ny,-1],Nvecwig,hardstop=.false.)
                  if(iwig_old.eq.0)cycle
                  !
                  !B
                  Gfull_R_halfkz(1:Norb,1+Norb:2*Norb,:,iwig,1) = Gfull_R(:,:,:,iwig_old,1)
                  Gfull_R_halfkz(1+Norb:2*Norb,1+2*Norb:3*Norb,:,iwig,1) = Gfull_R(:,:,:,iwig_old,1)
                  !AB
                  Gfull_R_halfkz(1:Norb,1+Norb:2*Norb,:,iwig,2) = Gfull_R(:,:,:,iwig_old,2)
                  Gfull_R_halfkz(1+Norb:2*Norb,1+2*Norb:3*Norb,:,iwig,2) = Gfull_R(:,:,:,iwig_old,2)
               elseif(nz.eq.+2)then
                  !distance +2
                  iwig_old = find_vec([nx,ny,+2],Nvecwig,hardstop=.false.)
                  if(iwig_old.eq.0)cycle
                  !
                  !B
                  Gfull_R_halfkz(1+2*Norb:3*Norb,1:Norb,:,iwig,1) = Gfull_R(:,:,:,iwig_old,1)
                  !AB
                  Gfull_R_halfkz(1+2*Norb:3*Norb,1:Norb,:,iwig,2) = Gfull_R(:,:,:,iwig_old,2)
               elseif(nz.eq.-2)then
                  !distance -2
                  iwig_old = find_vec([nx,ny,+2],Nvecwig,hardstop=.false.)
                  if(iwig_old.eq.0)cycle
                  !
                  !B
                  Gfull_R_halfkz(1:Norb,1+2*Norb:3*Norb,:,iwig,1) = Gfull_R(:,:,:,iwig_old,1)
                  !AB
                  Gfull_R_halfkz(1:Norb,1+2*Norb:3*Norb,:,iwig,2) = Gfull_R(:,:,:,iwig_old,2)
               endif
               !
            enddo
         enddo
      enddo
      deallocate(Gfull_R)
      write(*,*)"extraction to hetero-like"
      !
      Nkp = 151
      Nkz = 50
      !
      allocate(Gfull_R(3*Norb,3*Norb,Nmats,Lttc%Nkpt_path,2));Gfull_R=czero
      call wannier_R2K(Lttc%Nkpt3,Lttc%kptpath(:,1:Nkp),Gfull_R_halfkz(:,:,:,:,1),Gfull_R(:,:,:,1:Nkp,1))
      call wannier_R2K(Lttc%Nkpt3,Lttc%kptpath(:,1:Nkp),Gfull_R_halfkz(:,:,:,:,2),Gfull_R(:,:,:,1:Nkp,2))
      deallocate(Gfull_R_halfkz)
      write(*,*)"path-interpolation"
      !
      allocate(Gkw(2,Nmats,Nkp,0:Nkz));Gkw=czero
      do ik=1,Nkp
         do ikz=0,Nkz
            do isite=1,6
               do jsite=1,6
                  !
                  kR = 2*pi * Lttc%kptpath(3,Nkp+ikz) * (isite-jsite)
                  cfac = dcmplx(cos(kR),+sin(kR))
                  !
                  !B
                  Gkw(1,:,ik,ikz) = Gkw(1,:,ik,ikz) + Gfull_R(isite,jsite,:,ik,1)*cfac / 6
                  !AB
                  Gkw(2,:,ik,ikz) = Gkw(2,:,ik,ikz) + Gfull_R(isite,jsite,:,ik,2)*cfac / 6
                  !
               enddo
            enddo
         enddo
      enddo
      deallocate(Gfull_R)
      write(*,*)"Gamma-A filled-1"
      !
      do ispin=1,Nspin
         do iorb=1,2
            do ik=1,Nkp
               Gpath%wks(iorb,iorb,:,ik,ispin) = Gkw(iorb,:,ik,0)
            enddo
            do ikz=0,Nkz
               Gpath%wks(iorb,iorb,:,Nkp+ikz,ispin) = Gkw(iorb,:,Nkp,ikz)
            enddo
         enddo
      enddo
      deallocate(Gkw)
      write(*,*)"Gamma-A filled-2"
      !
   endif










   interface calc_energy_averages
      module procedure calc_energy_averages_integral
      module procedure calc_energy_averages_list
   end interface calc_energy_averages



      !
   subroutine calc_energy_averages_list(Inputs,Wk_orig,Lttc,Beta,pathOUTPUT,printmode)
      !
      use parameters
      use utils_misc
      use utils_fields
      use file_io
      use omp_lib
      use linalg, only : tensor_transform
      use crystal, only : fill_ksumkdiff, wannierinterpolation
      implicit none
      !
      type(SCDFT),intent(in)                :: Inputs
      complex(8),intent(in),target          :: Wk_orig(:,:,:,:)
      type(Lattice),intent(in)              :: Lttc
      real(8),intent(in)                    :: Beta
      character(len=*),intent(in)           :: pathOUTPUT
      character(len=*),intent(in)           :: printmode
      !
      integer                               :: iw,Nmats
      integer                               :: Ngrid,Norb,Nbp
      integer                               :: iorb,jorb,ib1,ib2,a,b,c,d
      integer                               :: Wk_dim,row,col,ndx,ndx1,ndx2
      integer                               :: ik1,ik2,iq,Nkpt
      integer                               :: iweig,jweig,iE1,iE2,iT,Efermi_ndx
      integer                               :: ithread,Nthread
      integer                               :: Ngrid_y=100
      real(8)                               :: E1,E2,dT,tanhs!DoS0
      real(8)                               :: DE_m,DE_p,nF_m,nF_p
      real(8)                               :: DosWeights,cutoff
      integer,allocatable                   :: kptsum(:,:),kptdif(:,:),map(:,:)
      real(8),allocatable                   :: Beta_DFT_list(:)
      complex(8),allocatable,target         :: Wk_interp(:,:,:,:)
      complex(8),pointer                    :: Wk_used(:,:,:,:)
      complex(8),allocatable                :: Wk_full(:,:)
      type(physicalU)                       :: PhysicalUelements
      real                                  :: start,finish
      logical                               :: Kstat_exists,Kdyn_exists
      logical                               :: calc_Kel_stat,calc_Wee_dyn
      logical                               :: Wrot_exists
      !
      !
      write(*,"(A)") new_line("A")//"---- calc_energy_averages_list"
      !
      !
      if(.not.initialized)stop "calc_energy_averages_list: input meshes not initialized. Call Initialize_inputs."
      !
      !Various checks
      Ngrid = size(Egrid)
      Nbp = size(Wk_orig,dim=1)
      Norb = int(sqrt(dble(Nbp)))
      Nmats = size(Wk_orig,dim=3)
      Nkpt = size(kpt_Model,dim=2)
      cutoff = Inputs%Wk_cutoff
      call assert_shape(Wk_orig,[Nbp,Nbp,Nmats,Lttc%Nkpt],"calc_energy_averages_list","Wk_orig")
      call assert_shape(kpt_Model,[3,Nkpt],"calc_energy_averages_list","kpt_Model")
      !
      Efermi_ndx = minloc(abs(Egrid),dim=1)
      !DoS0 = DoS_Model(Efermi_ndx)
      !write(*,"(A,F10.5)") new_line("A")//"     get_Kel_integral: Model DoS at the Fermi level:",DoS0
      if(Egrid(Efermi_ndx).ne.0d0) stop "get_Kel_integral: the energy grid requires the E=0 point."
      !
      if(Interpolate2Model)then
         if(Nkpt.eq.Lttc%Nkpt)stop "calc_energy_averages_list: something is wrong with the K-point dimension (interpolation)."
      else
         if(Nkpt.ne.Lttc%Nkpt)stop "calc_energy_averages_list: something is wrong with the K-point dimension."
      endif
      !
      allocate(wmats_orig(Nmats));wmats_orig=BosonicFreqMesh(Beta,Nmats)
      wmax_ndx = minloc(abs(wmats_orig-cutoff),dim=1)
      write(*,"(A)") "     Interaction frequency cut at iw_["//str(wmax_ndx)//"]="//str(wmats_orig(wmax_ndx),5)//" eV -> "//str(wmats_orig(wmax_ndx)*eV2DFTgrid,5)//" "//DFTgrid
      write(*,"(A)") "     Bosonic frequency step="//str(abs(wmats_orig(2)-wmats_orig(1)),5)//" eV -> "//str(abs(wmats_orig(2)-wmats_orig(1))*eV2DFTgrid,5)//" "//DFTgrid
      wmats_orig = wmats_orig * eV2DFTgrid
      wmax = wmats_orig(wmax_ndx)
      MatsStep = abs(wmats_orig(2)-wmats_orig(1))
      !
      !check if any of the needed kernels is already printed
      calc_Kel_stat=.false.
      calc_Wee_dyn=.false.
      !
      if(calc_Int_static)then
         !
         allocate(Kel_stat(Ngrid,Ngrid));Kel_stat=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Kel_stat.DAT",Kstat_exists,hardstop=.false.,verb=verbose)
         if(Kstat_exists)then
            call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","read")
         else
            calc_Kel_stat=.true.
            deallocate(Kel_stat)
         endif
         !
      endif
      !
      if(calc_Int_dynamic)then
         !
         allocate(Kel_dyn_list(Inputs%Tsteps,Ngrid,Ngrid));Kel_dyn_list=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Kel_dyn_list.DAT",Kdyn_exists,hardstop=.false.,verb=verbose)
         if(Kdyn_exists)then
            call io_Kel(Kel_dyn_list,reg(pathOUTPUT)//"Kel_dyn_list.DAT","read")
         else
            calc_Wee_dyn=.true.
            deallocate(Kel_dyn_list)
            !
            allocate(Beta_DFT_list(Inputs%Tsteps));Beta_DFT_list=0d0
            do iT=1,Inputs%Tsteps
               dT=0d0
               if(Inputs%Tsteps.gt.1) dT = (iT-1)*abs(Inputs%Tbounds(2)-Inputs%Tbounds(1))/dble(Inputs%Tsteps-1)
               Beta_DFT_list(iT) = 1d0 / ((Inputs%Tbounds(1) + dT)*K2eV*eV2DFTgrid)
            enddo
         endif
         !
      endif
      !
      if(calc_Kel_stat.or.calc_Wee_dyn)then
         !
         write(*,"(A)") "     One of the required Kernels is missing. Starting operations on screened interacion in Wannier basis."
         !
         if(Interpolate2Model)then
            !
            allocate(Wk_interp(Nbp,Nbp,wmax_ndx,size(kpt_Model,dim=2)));Wk_interp=czero
            call cpu_time(start)
            do iw=1,wmax_ndx
               call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,kpt_Model,Wk_orig(:,:,iw,:),Wk_interp(:,:,iw,:))
            enddo
            call cpu_time(finish)
            write(*,"(A,F)") "     Interpolating Wk to the new k grid cpu timing:", finish-start
            !
            Wk_used => Wk_interp
            !
         else
            !
            write(*,"(A)") "     Interpolation skipped, linking to original Wk."
            Wk_used => Wk_orig
            !
         endif
         !
         call init_Uelements(Norb,PhysicalUelements)
         !
         !rotation done externally and stored
         Wk_dim = (Nkpt*Norb) * (Nkpt*Norb+1)/2
         allocate(Wk_full(wmax_ndx,Wk_dim));Wk_full=czero
         !
         call inquireFile(reg(pathOUTPUT)//"Wrot.DAT",Wrot_exists,hardstop=.false.,verb=verbose)
         if(Wrot_exists)then
            !
            call cpu_time(start)
            call io_Kel(Wk_full,reg(pathOUTPUT)//"Wrot.DAT","read")
            call cpu_time(finish)
            write(*,"(A,F)") "     Read interaction in band basis cpu timing:", finish-start
            !
         else
            !
            !map between upper triangular and row,col
            call cpu_time(start)
            allocate(map(2,Wk_dim));map=0
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(row,col,ndx)
            !$OMP DO
            do row=1,Nkpt*Norb
               do col=row,Nkpt*Norb
                  !
                  ! upper triangular map (fixed)
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  map(1,ndx) = row
                  map(2,ndx) = col
                  !
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            write(*,"(A,F)") "     Upper triangular Wmap(dim="//str(Wk_dim)//") stored cpu timing:", finish-start
            !
            call fill_ksumkdiff(kpt_Model,kptsum,kptdif)
            deallocate(kptsum)
            !
            !rotation of the interaction W_(ab)(cd)(q) --> W_(i.k1,j.k2)(j.k2,i.k1)
            !W_(i.k1,j.k2)(j.k2,i.k1) = sum_abcd Zdag_ia(k1) * Z_bj(k2) * Zdag_jc(k2) * Z_di(k1) W_(ab)(cd)(q=k1-k2)
            !NOTE: only the upper triangular (c>=r) part of W_(i.k1,j.k2)(j.k2,i.k1) is computed
            call cpu_time(start)
            !$OMP PARALLEL DEFAULT(SHARED),&
            !$OMP PRIVATE(ndx,row,col,iorb,jorb,ik1,ik2,iq),&
            !$OMP PRIVATE(a,b,c,d,ib1,ib2,ithread)
            Nthread = omp_get_num_threads()
            !$OMP DO
            do ndx=1,Wk_dim
               !
               ithread = omp_get_thread_num()
               print *, "thread", ithread, " / ", Nthread, " ndx: ", ndx, " over: ", Wk_dim
               !
               row = map(1,ndx)
               col = map(2,ndx)
               !
               !row = ik1 + (iorb-1)*Nkpt
               iorb = floor((row-0.01)/Nkpt)+1
               ik1  = row - (iorb-1)*Nkpt
               !col = ik2 + (jorb-1)*Nkpt
               jorb = floor((col-0.01)/Nkpt)+1
               ik2  = col - (jorb-1)*Nkpt
               !
               iq = kptdif(ik1,ik2)
               !
               do ib1=1,Nbp
                  !
                  !diagonal elements
                  a = PhysicalUelements%Full_Map(ib1,ib1,1)
                  b = PhysicalUelements%Full_Map(ib1,ib1,2)
                  c = PhysicalUelements%Full_Map(ib1,ib1,3)
                  d = PhysicalUelements%Full_Map(ib1,ib1,4)
                  !
                  Wk_full(:,ndx) = Wk_full(:,ndx)                                                                      &
                                 + Wk_used(ib1,ib1,1:wmax_ndx,iq) * conjg(Zk_Model(a,iorb,ik1)) * Zk_Model(b,jorb,ik2) &
                                                                  * conjg(Zk_Model(c,jorb,ik2)) * Zk_Model(d,iorb,ik1)
                  !
                  do ib2=ib1+1,Nbp
                     !
                     !off-diagonal elements
                     a = PhysicalUelements%Full_Map(ib1,ib2,1)
                     b = PhysicalUelements%Full_Map(ib1,ib2,2)
                     c = PhysicalUelements%Full_Map(ib1,ib2,3)
                     d = PhysicalUelements%Full_Map(ib1,ib2,4)
                     !
                     Wk_full(:,ndx) = Wk_full(:,ndx)                                                                       &
                                    + Wk_used(ib1,ib2,1:wmax_ndx,iq) * conjg(Zk_Model(a,iorb,ik1)) * Zk_Model(b,jorb,ik2)  &
                                                                     * conjg(Zk_Model(c,jorb,ik2)) * Zk_Model(d,iorb,ik1)  &
                                    + Wk_used(ib2,ib1,1:wmax_ndx,iq) * conjg(Zk_Model(c,iorb,ik1)) * Zk_Model(d,jorb,ik2)  &
                                                                     * conjg(Zk_Model(a,jorb,ik2)) * Zk_Model(b,iorb,ik1)
                     !
                  enddo
               enddo
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            call cpu_time(finish)
            deallocate(map,kptdif)
            !
            !print the interaction in band basis depending explicitly of the two kpoints
            call io_Kel(Wk_full,reg(pathOUTPUT)//"Wrot.DAT","write")
            !
            write(*,"(A,F)") "     Rotation of the interaction to band basis cpu timing:", finish-start
            !
         endif
         if(allocated(Wk_interp))deallocate(Wk_interp)
         nullify(Wk_used)
         !
         !Kernels calculation
         call cpu_time(start)
         if(calc_Kel_stat)then
            allocate(Kel_stat(Ngrid,Ngrid))
            Kel_stat=czero
         endif
         if(calc_Wee_dyn)then
            allocate(Kel_dyn_list(Inputs%Tsteps,Ngrid,Ngrid))
            Kel_dyn_list=czero
         endif
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(iweig,jweig,iE1,iE2,DosWeights),&
         !$OMP PRIVATE(ndx,ndx1,ndx2,row,col),&
         !$OMP PRIVATE(iorb,jorb,ik1,ik2,ithread),&
         !$OMP PRIVATE(E1,E2,DE_m,DE_p,nF_m,nF_p,iT,tanhs)
         Nthread = omp_get_num_threads()
         !$OMP DO
         do jweig=1,size(finite_weights_Model,dim=1)
            !
            iE2 = finite_weights_Model(jweig,1)
            jorb = finite_weights_Model(jweig,2)
            ik2 = finite_weights_Model(jweig,3)
            !
            ithread = omp_get_thread_num()
            print *, "thread", ithread, " / ", Nthread, " jweig: ", jweig, " over: ", size(finite_weights_Model,dim=1)
            !
            do iweig=1,size(finite_weights_Model,dim=1)
               !
               iE1 = finite_weights_Model(iweig,1)
               iorb = finite_weights_Model(iweig,2)
               ik1 = finite_weights_Model(iweig,3)
               !
               DosWeights = (weights_Model(iE1,iorb,ik1)/DoS_Model(iE1)) * (weights_Model(iE2,jorb,ik2)/DoS_Model(iE2))
               !
               E1 = Ek_Model(iorb,ik1)
               E2 = Ek_Model(jorb,ik2)
               !
               DE_m = E1-E2
               nF_m = fermidirac(+E1,Beta_DFT_list(iT))-fermidirac(+E2,Beta_DFT_list(iT))
               DE_p = E1+E2
               nF_p = fermidirac(-E1,Beta_DFT_list(iT))-fermidirac(+E2,Beta_DFT_list(iT))
               !
               !product basis map, the indexes spanned by iweig,jweig cover all the possible
               !(ik1,iorb) pairs, so the whole Wk_full matrix, both LT and UT.
               ndx1 = ik1 + (iorb-1)*Nkpt
               ndx2 = ik2 + (jorb-1)*Nkpt
               !
               if(ndx2.ge.ndx1)then
                  !
                  !I'm looking for an element in the UT. ndx gives me the position
                  row = ndx1 !this is the row
                  col = ndx2 !this is the col
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  !static Kernel
                  if(calc_Kel_stat) Kel_stat(iE1,iE2) = Kel_stat(iE1,iE2) + Wk_full(1,ndx) * eV2DFTgrid * DosWeights
                  !
                  !dynamic Kernel stored for all the required temeperatures
                  if(calc_Wee_dyn)then
                     !
                     do iT=1,Inputs%Tsteps
                        !
                        tanhs = 0d0
                        if((E1*E2).ne.0d0) tanhs = 1d0/(tanh(Beta_DFT_list(iT)*E1/2d0)*tanh(Beta_DFT_list(iT)*E2/2d0))
                        !
                        Kel_dyn_list(iT,iE1,iE2) = Kel_dyn_list(iT,iE1,iE2) + DosWeights * eV2DFTgrid * (2d0/pi) * tanhs *  &
                        (                                                                                                                                                             &
                           nF_m * ( aux_integral(DE_m,(Wk_full(:,ndx)-Wk_full(1,ndx))) + (Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_m)) ) +                         &
                           nF_p * ( aux_integral(DE_p,(Wk_full(:,ndx)-Wk_full(1,ndx))) + (Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_p)) )                           &
                        )
                        !
                     enddo
                     !
                  endif
                  !
               else
                  !
                  !I'm looking for an element in the LT. I look via ndx his complex conjg in the UT
                  row = ndx2 !this is the row
                  col = ndx1 !this is the col
                  ndx = (Nkpt*Norb)*(row-1) - (row-1)*row/2 + col
                  !
                  !static Kernel
                  if(calc_Kel_stat) Kel_stat(iE1,iE2) = Kel_stat(iE1,iE2) + conjg(Wk_full(1,ndx)) * eV2DFTgrid * DosWeights
                  !
                  !dynamic Kernel stored for all the required temeperatures
                  if(calc_Wee_dyn)then
                     !
                     do iT=1,Inputs%Tsteps
                        !
                        tanhs = 0d0
                        if((E1*E2).ne.0d0) tanhs = 1d0/(tanh(Beta_DFT_list(iT)*E1/2d0)*tanh(Beta_DFT_list(iT)*E2/2d0))
                        !
                        Kel_dyn_list(iT,iE1,iE2) = Kel_dyn_list(iT,iE1,iE2) + DosWeights * eV2DFTgrid * (2d0/pi) * tanhs *  &
                        (                                                                                                                                                             &
                           nF_m * ( aux_integral(DE_m,conjg(Wk_full(:,ndx)-Wk_full(1,ndx))) + conjg(Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_m)) ) +               &
                           nF_p * ( aux_integral(DE_p,conjg(Wk_full(:,ndx)-Wk_full(1,ndx))) + conjg(Wk_full(wmax_ndx,ndx)-Wk_full(1,ndx))*(pi/2d0-atan2(wmax,DE_p)) )                 &
                        )
                        !
                     enddo
                     !
                  endif
                  !
               endif
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
         deallocate(Wk_full)
         call cpu_time(finish)
         write(*,"(A,F)") "     Calculation of static electronic Kernel cpu timing:", finish-start
         !
         !Filling the Fermi lines
         do iT=1,Inputs%Tsteps
            call interpFermi(Kel_dyn_list(iT,:,:),Egrid,Egrid,Efermi_ndx,Efermi_ndx)
         enddo
         !
         if(calc_Kel_stat) call io_Kel(Kel_stat,reg(pathOUTPUT)//"Kel_stat.DAT","write")
         if(calc_Wee_dyn) call io_Kel(Kel_dyn_list,reg(pathOUTPUT)//"Kel_dyn_list.DAT","write")
         !
      endif
      deallocate(weights_Model,finite_weights_Model)
      !
      if(calc_Int_static.and.(reg(printmode).ne."None"))then
         call print_Kernel("electronic",reg(printmode),reg(pathOUTPUT),"Kel_stat",Egrid,Egrid,Kel_stat)
      endif
      !
      Kernels_stored = .true.
      !
      !
      !
   contains
      !
      !
      !
      function ygrid(stop,num,ndx) result(yval)
         implicit none
         real(8),intent(in)                    :: stop
         integer,intent(in)                    :: num,ndx
         real(8)                               :: yval
         real(8)                               :: start,step
         !
         if(num.lt.0)stop "ygrid: N<0, abort."
         if(ndx.le.0)stop "ygrid: ndx<=0, abort."
         start = -1d0
         step = (stop-start)/(dble(num)+1d0)
         yval = start + dble(ndx)*step
         !
      end function ygrid
      !
      function aux_integral(DE,W) result(Integral)
         implicit none
         real(8),intent(in)                    :: DE
         complex(8),intent(in)                 :: W(:)
         complex(8)                            :: Integral
         integer                               :: iy,wndx_a,wndx_b
         real(8)                               :: ymax,dy,y_i,y_j,wm
         real(8)                               :: ReW_wm_intp,ImW_wm_intp
         complex(8)                            :: W_wm_i,W_wm_j,Int_i,Int_j
         !
         ymax = (wmax-abs(DE)) / (wmax+abs(DE))
         dy = abs(ygrid(ymax,Ngrid_y,2)-ygrid(ymax,Ngrid_y,1))
         !
         Integral=czero
         if((DE.ne.0d0).and.(ymax.ne.1d0))then
            !
            do iy=2,Ngrid_y
               !
               !continous frequency correspnding to iy
               y_i = ygrid(ymax,Ngrid_y,iy)
               wm = abs(DE) * (1+y_i)/(1-y_i)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               if(wndx_b.gt.wmax_ndx) stop"aux_integral (DE)_i: the frequency index is beyond the cutoff."
               ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(W(wndx_a))] , [wmats_orig(wndx_b),dreal(W(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(W(wndx_a))] , [wmats_orig(wndx_b),dimag(W(wndx_b))] , wm )
               W_wm_i = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !continous frequency correspnding to iy-1
               y_j = ygrid(ymax,Ngrid_y,iy-1)
               wm = abs(DE) * (1+y_j)/(1-y_j)
               !linear interpolation of Wee between the two points on the matsubara grid enclosing wm
               wndx_a = floor(wm/MatsStep) + 1
               wndx_b = wndx_a + 1
               if(wndx_b.gt.wmax_ndx) stop"aux_integral (DE)_j: the frequency index is beyond the cutoff."
               ReW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dreal(W(wndx_a))] , [wmats_orig(wndx_b),dreal(W(wndx_b))] , wm )
               ImW_wm_intp = linear_interp_2y( [wmats_orig(wndx_a),dimag(W(wndx_a))] , [wmats_orig(wndx_b),dimag(W(wndx_b))] , wm )
               W_wm_j = dcmplx(ReW_wm_intp,ImW_wm_intp)
               !
               !integrand for iy and iy-1
               Int_i = W_wm_i / ( 1d0 + y_i**2 )
               Int_j = W_wm_j / ( 1d0 + y_j**2 )
               !
               !trapezoidal integration
               Integral = Integral + sign(1d0,DE)*(Int_i+Int_j)*(dy/2d0)
               !
            enddo
            !
         endif
         !
      end function aux_integral
      !
      !
      !
   end subroutine calc_energy_averages_list

!Store inside the module the required energy averages
   if(calc_Kel)then
      select case(reg(Inputs%mode_avg))
         case default
            stop "Available entries for MODE_AVG: integral, list."
         case("integral")
            call calc_energy_averages(Wlat%screened,Lttc,Wlat%Beta,Inputs%Wk_cutoff,reg(printpath),reg(Inputs%printmode_el))
         case("list")
            call calc_energy_averages(Inputs,Wlat%screened,Lttc,Wlat%Beta,reg(printpath),reg(Inputs%printmode_el))
      end select
   endif


   if(calc_Kel)then
      allocate(Kel(Ngrid,Ngrid));Kel=0d0
      select case(reg(Inputs%mode_avg))
         case default
            stop "Available entries for MODE_AVG: integral, list."
         case("integral")
            call get_Kel(Kel,Beta_DFT,reg(Inputs%printmode_el),reg(printpath_T))
         case("list")
            call get_Kel(Kel,iT,Beta_DFT,reg(Inputs%printmode_el),reg(printpath_T))
      end select
   endif


   !
   subroutine get_Kel_list(Kel,iT,Beta,printmode,printKpath)
      !
      use parameters
      use utils_misc
      implicit none
      !
      complex(8),intent(out)                :: Kel(:,:)
      integer,intent(in)                    :: iT
      real(8),intent(in)                    :: Beta
      character(len=*),intent(in)           :: printmode
      character(len=*),intent(in)           :: printKpath
      integer                               :: Ngrid
      real(8)                               :: Temp
      !
      !
      if(verbose)write(*,"(A)") "---- get_Kel_list"
      !
      !
      if(.not.calc_Kel)stop "get_Kel_list: inputs not initialized. call Initialize_inputs."
      if(.not.Kernels_stored)stop "get_Kel_list: fully screened interaction not stored. call calc_energy_averages."
      if(calc_Int_static.and.(.not.allocated(Kel_stat)))stop "get_Kel_list: strange Kel_stat should be allocated."
      if(calc_Int_dynamic.and.(.not.allocated(Kel_dyn_list)))stop "get_Kel_list: strange Kel_dyn_list should be allocated."
      !
      Ngrid = size(Egrid)
      call assert_shape(Kel,[Ngrid,Ngrid],"get_Kel_list","Kel")
      !
      Kel=czero
      if(calc_Int_dynamic)then
         !
         Kel = Kel_dyn_list(iT,:,:)
         !
         if(reg(printmode).ne."None")then
            Temp = 1d0 / (K2eV*eV2DFTgrid*Beta)
            call print_Kernel("electronic",reg(printmode),reg(printKpath),"Kel_dyn_T"//str(Temp,2),Egrid,Egrid,Kel)
         endif
         !
      endif
      !
      if(calc_Int_static) Kel = Kel + Kel_stat
      !
   end subroutine get_Kel_list