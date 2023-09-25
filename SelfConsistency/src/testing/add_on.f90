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
