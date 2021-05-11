module post_processing

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   !
   !

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface dump_MaxEnt
      module procedure :: dump_MaxEnt_Gfunct                                    ![Array(Norb,Np,Nspin),Beta,mode,dirpath,filename,iorb(optinal)]
      module procedure :: dump_MaxEnt_Gfield                                    ![FermionicField,mode,dirpath,filename,orbset(:,:)]
      module procedure :: dump_MaxEnt_Wfunct                                    ![Array(Np),Beta,mode,dirpath,filename,iorb(2)]
      module procedure :: dump_MaxEnt_Wfield                                    ![BosonicField,mode,dirpath,filename,orbset(:,:)]
   end interface dump_MaxEnt

   interface interpolate2Beta
      module procedure :: interpolate2Beta_Fermionic                            ![FermionicField,OldBeta,offDiag]
      module procedure :: interpolate2Beta_Bosonic                              ![BosonicField,OldBeta,offDiag]
   end interface interpolate2Beta

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !functions
   public :: pade
   !subroutines
   public :: dump_MaxEnt
   public :: remove_CDW
   public :: interpolate2Beta
   public :: interpolateG2Path

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine remove_CDW(W,mode,site)
      !
      use parameters
      use utils_misc
      use utils_fields
      use interactions
      use linalg
      use input_vars, only: Nsite, SiteNorb, SiteOrbs
      implicit none
      !
      type(BosonicField),intent(inout)      :: W
      character(len=*),intent(in)           :: mode
      integer,intent(in),optional           :: site
      !
      integer                               :: Norb,Nmats,iw
      integer                               :: isite,ib1,ib2
      integer                               :: i,j,k,l
      integer                               :: i_lat,j_lat,k_lat,l_lat
      logical                               :: replace
      real(8)                               :: ReW,ImW
      real(8),allocatable                   :: wmats(:),Wdiag(:,:)
      complex(8),allocatable                :: Rot(:,:)
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- remove_CDW"
      !
      !
      if(.not.W%status) stop "remove_CDW: BosonicField not properly initialized."
      Norb = int(sqrt(dble(W%Nbp)))
      Nmats = W%Npoints
      !
      allocate(wmats(Nmats))
      wmats=BosonicFreqMesh(W%Beta,Nmats)
      !
      select case(reg(mode))
         case default
            !
            stop "remove_CDW: Available Modes are: imp, imp_exp, diag, lat."
            !
         case("imp")
            !
            !Removing CDW from all the (aa)(aa) and (aa)(bb) input indexes
            !
            call init_Uelements(Norb,PhysicalUelements)
            !
            do ib1=1,W%Nbp
               do ib2=1,W%Nbp
                  !
                  replace = PhysicalUelements%Full_Uaa(ib1,ib2) .or. PhysicalUelements%Full_Uab(ib1,ib2)
                  !
                  if(replace)then
                     ReW = cubic_interp(wmats(2:Nmats),real(W%screened_local(ib1,ib2,2:Nmats)),0d0)
                     W%screened_local(ib1,ib2,1) = dcmplx(ReW,0d0)
                  endif
                  !
               enddo
            enddo
            !
         case("imp_exp")
            !
            !Removing CDW from the (aa)(aa) and (aa)(bb) input indexes belonging to all or given sites
            !
            call init_Uelements(Norb,PhysicalUelements)
            !
            do isite=1,Nsite
               !
               if(present(site).and.(isite.ne.site))cycle
               !
               do i=1,SiteNorb(isite)
                  do j=1,SiteNorb(isite)
                     do k=1,SiteNorb(isite)
                        do l=1,SiteNorb(isite)
                           !
                           ! mapping
                           i_lat = SiteOrbs(isite,i)
                           j_lat = SiteOrbs(isite,j)
                           k_lat = SiteOrbs(isite,k)
                           l_lat = SiteOrbs(isite,l)
                           !
                           if(any([i_lat,j_lat,k_lat,l_lat].gt.Norb)) stop "remove_CDW: the input field is not in the lattice space."
                           !
                           ! bosonic indexes on the lattice
                           ib1 = i_lat + Norb*(j_lat-1)
                           ib2 = k_lat + Norb*(l_lat-1)
                           !
                           replace = PhysicalUelements%Full_Uaa(ib1,ib2) .or. PhysicalUelements%Full_Uab(ib1,ib2)
                           !
                           if(replace)then
                              ReW = cubic_interp(wmats(2:Nmats),real(W%screened_local(ib1,ib2,2:Nmats)),0d0)
                              W%screened_local(ib1,ib2,1) = dcmplx(ReW,0d0)
                           endif
                           !
                        enddo
                     enddo
                  enddo
               enddo
               !
            enddo
            !
         case("diag")
            !
            !Removing CDW from all the input indexes in the diagonal basis
            !
            allocate(Wdiag(W%Nbp,Nmats));Wdiag=0d0
            !
            !the rotation is the one given by the second freq
            allocate(Rot(W%Nbp,W%Nbp))
            Rot=W%screened_local(:,:,2)
            call eigh(Rot,Wdiag(:,2))
            !
            !rotate the orbital space at all freq
            do iw=2,Nmats
               Wdiag(:,iw) = diagonal(rotate(W%screened_local(:,:,iw),Rot))
            enddo
            !
            !fit the missing point in the diagonal basis
            do ib1=1,W%Nbp
               Wdiag(ib1,1) = cubic_interp(wmats(2:Nmats),Wdiag(ib1,2:Nmats),0d0)
            enddo
            !
            !bring back to the original basis the fitted point
            W%screened_local(:,:,1) = rotate(diag(Wdiag(:,1)),transpose(conjg(Rot)))
            where(abs(W%screened_local(:,:,1))<1e-6)W%screened_local(:,:,1)=czero
            !
         case("lat")
            !
            !Removing CDW from all the input indexes
            !
            do ib1=1,W%Nbp
               do ib2=1,W%Nbp
                  !
                  ReW = cubic_interp(wmats(2:Nmats),real(W%screened_local(ib1,ib2,2:Nmats)),0d0)
                  ImW = cubic_interp(wmats(2:Nmats),aimag(W%screened_local(ib1,ib2,2:Nmats)),0d0)
                  W%screened_local(ib1,ib2,1) = dcmplx(ReW,ImW)
                  !
               enddo
            enddo
            !
      end select
      !
   end subroutine remove_CDW


   !---------------------------------------------------------------------------!
   !PURPOSE:
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine dump_MaxEnt_Gfunct(G,mode,dirpath,filename,iorb,WmaxPade)
      !
      use parameters
      use utils_misc
      use input_vars, only: Nmats, Solver, Beta, Nreal, wrealMax
      use fourier_transforms
      use file_io
      implicit none
      !
      complex(8),intent(in)                 :: G(:,:,:)
      character(len=*),intent(in)           :: mode
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      integer,intent(in),optional           :: iorb
      integer,intent(in),optional           :: WmaxPade
      !
      complex(8),allocatable                :: Gft(:,:,:)
      complex(8),allocatable                :: Gpade(:)
      integer                               :: iwan,ispin
      integer                               :: Norb,Npoints
      real(8),allocatable                   :: tau(:),wmats(:),wreal(:)
      !
      !
      if(verbose)write(*,"(A)") "---- dump_MaxEnt_Gfunct"
      !
      !
      Norb = size(G,dim=1)
      Npoints = size(G,dim=2)
      if(Nspin.ne.size(G,dim=3)) stop "dump_MaxEnt_Gfunct: wrong Nspin dimension."
      !
      select case(reg(mode))
         case default
            !
            stop "Available Modes are: itau, mats, itau2mats, mats2itau."
            !
         case("itau")
            !
            if(Npoints.ne.Solver%NtauF) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct number fo imaginary time points differ from input."
            if(present(WmaxPade).and.(WmaxPade.gt.0)) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct pade from itau not done."
            allocate(tau(Npoints));tau=0d0
            tau = linspace(0d0,Beta,Npoints)
            !
            do iwan=1,Norb
               if(present(iorb).and.(iwan.ne.iorb))cycle
               do ispin=1,Nspin
                  call dump_Field_component(real(G(iwan,:,ispin)),reg(dirpath),reg(filename)//"_t_o"//str(iwan)//"_s"//str(ispin)//".DAT",tau)
               enddo
            enddo
            !
         case("mats")
            !
            if(Npoints.ne.Nmats) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct number fo Matsubara points differ from input."
            allocate(wmats(Npoints));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Npoints)
            !
            do iwan=1,Norb
               if(present(iorb).and.(iwan.ne.iorb))cycle
               do ispin=1,Nspin
                  !
                  call dump_Field_component(G(iwan,:,ispin),reg(dirpath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//".DAT",wmats)
                  !
                  if(present(WmaxPade).and.(WmaxPade.gt.0))then
                     allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
                     allocate(Gpade(Nreal));Gpade=czero
                     Gpade = pade(G(iwan,:,ispin),"Fermionic",wlimit=WmaxPade)
                     call dump_Field_component(Gpade,reg(dirpath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//"_pade.DAT",wreal)
                     deallocate(Gpade,wreal)
                  endif
                  !
               enddo
            enddo
            !
         case("itau2mats")
            !
            if(Npoints.ne.Solver%NtauF) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct number fo imaginary time points differ from input."
            allocate(wmats(Nmats));wmats=0d0
            wmats = FermionicFreqMesh(Beta,Nmats)
            !
            allocate(Gft(Norb,Nmats,Nspin));Gft=czero
            do ispin=1,Nspin
               call Fitau2mats_vec(Beta,G(:,:,ispin),Gft(:,:,ispin),tau_uniform=.true.)
            enddo
            !
            do iwan=1,Norb
               if(present(iorb).and.(iwan.ne.iorb))cycle
               do ispin=1,Nspin
                  !
                  call dump_Field_component(Gft(iwan,:,ispin),reg(dirpath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//".DAT",wmats)
                  !
                  if(present(WmaxPade).and.(WmaxPade.gt.0))then
                     allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
                     allocate(Gpade(Nreal));Gpade=czero
                     Gpade = pade(G(iwan,:,ispin),"Fermionic",wlimit=WmaxPade)
                     call dump_Field_component(Gpade,reg(dirpath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//"_pade.DAT",wreal)
                     deallocate(Gpade,wreal)
                  endif
                  !
               enddo
            enddo
            deallocate(Gft)
            !
         case("mats2itau")
            !
            if(Npoints.ne.Nmats) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct number fo Matsubara points differ from input."
            if(present(WmaxPade).and.(WmaxPade.gt.0)) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct pade from mats2itau not done."
            allocate(tau(Solver%NtauF));tau=0d0
            tau = linspace(0d0,Beta,Solver%NtauF)
            !
            allocate(Gft(Norb,Solver%NtauF,Nspin));Gft=czero
            do ispin=1,Nspin
               call Fmats2itau_vec(Beta,G(:,:,ispin),Gft(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
            enddo
            !
            do iwan=1,Norb
               if(present(iorb).and.(iwan.ne.iorb))cycle
               do ispin=1,Nspin
                  call dump_Field_component(real(Gft(iwan,:,ispin)),reg(dirpath),reg(filename)//"_t_o"//str(iwan)//"_s"//str(ispin)//".DAT",tau)
               enddo
            enddo
            deallocate(Gft)
            !
      end select
      !
   end subroutine dump_MaxEnt_Gfunct
   !
   subroutine dump_MaxEnt_Gfield(G,mode,dirpath,filename,Orbs,WmaxPade)
      !
      use parameters
      use utils_misc
      use input_vars, only: Beta
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(in)       :: G
      character(len=*),intent(in)           :: mode
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      integer,allocatable,intent(in)        :: Orbs(:,:)
      integer,intent(in),optional           :: WmaxPade
      !
      integer                               :: iwan,iset
      complex(8),allocatable                :: Gprint(:,:,:)
      !
      !
      if(verbose)write(*,"(A)") "---- dump_MaxEnt_Gfield"
      !
      !
      if(.not.G%status) stop "dump_MaxEnt_Gfield: Fermionic field not properly initialized."
      if(G%Beta.ne.Beta) stop "dump_MaxEnt_Gfield: Beta attribute of the field does not match with data from input file."
      !
      allocate(Gprint(G%Norb,G%Npoints,Nspin));Gprint=czero
      do iwan=1,G%Norb
         Gprint(iwan,:,:) = G%ws(iwan,iwan,:,:)
      enddo
      !
      if(allocated(Orbs))then
         do iset=1,size(Orbs,dim=1)
            call dump_MaxEnt_Gfunct(Gprint,reg(mode),reg(dirpath),reg(filename),iorb=Orbs(iset,1))
            if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Gfunct(Gprint,reg(mode),reg(dirpath),reg(filename),iorb=Orbs(iset,1),WmaxPade=WmaxPade)
         enddo
      else
         call dump_MaxEnt_Gfunct(Gprint,reg(mode),reg(dirpath),reg(filename))
         if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Gfunct(Gprint,reg(mode),reg(dirpath),reg(filename),WmaxPade=WmaxPade)
      endif
      !
   end subroutine dump_MaxEnt_Gfield
   !
   subroutine dump_MaxEnt_Wfunct(W,mode,dirpath,filename,iorb,type,WmaxPade)
      !
      use parameters
      use utils_misc
      use input_vars, only: Nmats, Solver, Beta, Nreal, wrealMax
      use fourier_transforms
      use file_io
      implicit none
      !
      complex(8),intent(in)                 :: W(:)
      character(len=*),intent(in)           :: mode
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      integer,intent(in)                    :: iorb(2)
      character(len=*),intent(in)           :: type
      integer,intent(in),optional           :: WmaxPade
      !
      integer                               :: ndx(4)
      complex(8),allocatable                :: Wft(:),Wpade(:)
      integer                               :: Npoints
      real(8),allocatable                   :: tau(:),wmats(:),wreal(:)
      !
      !
      if(verbose)write(*,"(A)") "---- dump_MaxEnt_Wfunct"
      !
      !
      Npoints = size(W)
      !
      select case(reg(type))
         case default
            stop "Available types are: Uaa, Uab, J."
         case("Uaa")
            if(iorb(1).ne.iorb(2)) stop "dump_MaxEnt_Wfunct: wrong indexes provided for type Uaa."
            ndx(1)=iorb(1)
            ndx(2)=iorb(2)
            ndx(3)=iorb(2)
            ndx(4)=iorb(1)
         case("J")
            if(iorb(1).eq.iorb(2)) stop "dump_MaxEnt_Wfunct: wrong indexes provided for type J."
            ndx(1)=iorb(1)
            ndx(2)=iorb(2)
            ndx(3)=iorb(2)
            ndx(4)=iorb(1)
         case("Uab")
            if(iorb(1).eq.iorb(2)) stop "dump_MaxEnt_Wfunct: wrong indexes provided for type Uab."
            ndx(1)=iorb(1)
            ndx(2)=iorb(1)
            ndx(3)=iorb(2)
            ndx(4)=iorb(2)
      end select
      !
      select case(reg(mode))
         case default
            !
            stop "Available Modes are: itau, mats, itau2mats, mats2itau."
            !
         case("itau")
            !
            if(Npoints.ne.Solver%NtauB) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct number fo imaginary time points differ from input."
            if(present(WmaxPade).and.(WmaxPade.gt.0)) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct pade from itau not done."
            allocate(tau(Npoints));tau=0d0
            tau = linspace(0d0,Beta,Npoints)
            !
            call dump_Field_component(real(W),reg(dirpath),reg(filename)//"_t_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",tau)
            !
         case("mats")
            !
            if(Npoints.ne.Nmats) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct number fo Matsubara points differ from input."
            allocate(wmats(Npoints));wmats=0d0
            wmats = BosonicFreqMesh(Beta,Npoints)
            !
            call dump_Field_component(real(W),reg(dirpath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",wmats)
            !
            if(present(WmaxPade).and.(WmaxPade.gt.0))then
               allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
               allocate(Wpade(Nreal));Wpade=czero
               Wpade = pade(W,"Bosonic",wlimit=WmaxPade)
               call dump_Field_component(real(Wpade),reg(dirpath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//")_pade.DAT",wreal)
               deallocate(Wpade,wreal)
            endif
            !
         case("itau2mats")
            !
            if(Npoints.ne.Solver%NtauB) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct number fo imaginary time points differ from input."
            allocate(wmats(Nmats));wmats=0d0
            wmats = BosonicFreqMesh(Beta,Nmats)
            !
            allocate(Wft(Nmats));Wft=czero
            call Bitau2mats(Beta,W,Wft,tau_uniform=.true.)
            !
            call dump_Field_component(real(Wft),reg(dirpath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",wmats)
            !
            if(present(WmaxPade).and.(WmaxPade.gt.0))then
               allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
               allocate(Wpade(Nreal));Wpade=czero
               Wpade = pade(Wft,"Bosonic",wlimit=WmaxPade)
               call dump_Field_component(real(Wpade),reg(dirpath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//")_pade.DAT",wreal)
               deallocate(Wpade,wreal)
            endif
            !
            deallocate(Wft)
            !
         case("mats2itau")
            !
            if(Npoints.ne.Nmats) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct number fo Matsubara points differ from input."
            if(present(WmaxPade).and.(WmaxPade.gt.0)) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct pade from mats2itau not done."
            allocate(tau(Solver%NtauB));tau=0d0
            tau = linspace(0d0,Beta,Solver%NtauB)
            !
            allocate(Wft(Solver%NtauB));Wft=czero
            call Bmats2itau(Beta,W,Wft,asympt_corr=.true.,tau_uniform=.true.,Umats_bare=W(Npoints))
            Wft = Wft + W(Npoints)
            !
            call dump_Field_component(real(Wft),reg(dirpath),reg(filename)//"_t_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",tau)
            deallocate(Wft)
            !
      end select
      !
   end subroutine dump_MaxEnt_Wfunct
   !
   subroutine dump_MaxEnt_Wfield(W,mode,dirpath,filename,Orbs,WmaxPade)
      !
      use parameters
      use utils_misc
      use input_vars, only: Beta
      use utils_fields
      use file_io
      implicit none
      !
      type(BosonicField),intent(in)         :: W
      character(len=*),intent(in)           :: mode
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      integer,allocatable,intent(in)        :: Orbs(:,:)
      integer,intent(in),optional           :: WmaxPade
      !
      integer                               :: iwan1,iwan2,iJ1,iJ2,iU1,iU2
      integer                               :: Norb,iset
      !
      !
      if(verbose)write(*,"(A)") "---- dump_MaxEnt_Wfield"
      !
      !
      if(.not.W%status) stop "dump_MaxEnt_Wfield: Bosonic field not properly initialized."
      if(W%Beta.ne.Beta) stop "dump_MaxEnt_Wfield: Beta attribute of the field does not match with data from input file."
      Norb = int(sqrt(dble(W%Nbp)))
      !
      if(allocated(Orbs))then
         do iset=1,size(Orbs,dim=1)
            iwan1 = Orbs(iset,1)
            iwan2 = Orbs(iset,2)
            iU1 = iwan1+Norb*(iwan1-1)
            call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU1,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan1],"Uaa")
            if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU1,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan1],"Uaa",WmaxPade=WmaxPade)
            !
            iU2 = iwan2+Norb*(iwan2-1)
            iJ1 = iwan1+Norb*(iwan2-1)
            iJ2 = iwan2+Norb*(iwan1-1)
            if(iwan2.ne.0)then
               call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"Uab")
               if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"Uab",WmaxPade=WmaxPade)
               call dump_MaxEnt_Wfunct(W%screened_local(iJ1,iJ2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"J")
               if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Wfunct(W%screened_local(iJ1,iJ2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"J",WmaxPade=WmaxPade)
            endif
         enddo
      else
         do iwan1=1,Norb
            iU1 = iwan1+Norb*(iwan1-1)
            call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU1,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan1],"Uaa")
            if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU1,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan1],"Uaa",WmaxPade=WmaxPade)
            do iwan2=1+iwan1,Norb
               iU2 = iwan2+Norb*(iwan2-1)
               iJ1 = iwan1+Norb*(iwan2-1)
               iJ2 = iwan2+Norb*(iwan1-1)
               call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"Uab")
               if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Wfunct(W%screened_local(iU1,iU2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"Uab",WmaxPade=WmaxPade)
               call dump_MaxEnt_Wfunct(W%screened_local(iJ1,iJ2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"J")
               if(present(WmaxPade).and.(WmaxPade.gt.0))call dump_MaxEnt_Wfunct(W%screened_local(iJ1,iJ2,:),reg(mode),reg(dirpath),reg(filename),[iwan1,iwan2],"J",WmaxPade=WmaxPade)
            enddo
         enddo
      endif
      !
   end subroutine dump_MaxEnt_Wfield


   !---------------------------------------------------------------------------!
   !PURPOSE: Pade analytic continuation
   !TEST ON:
   !---------------------------------------------------------------------------!
   function pade(funct_in,type,wlimit) result(funct_out)
      !
      use parameters
      use utils_misc
      use input_vars, only: Nmats, Nreal, wrealMax, eta, Beta
      implicit none
      !
      complex(8),intent(in)                 :: funct_in(:)
      character(len=*),intent(in)           :: type
      integer,intent(in),optional           :: wlimit
      complex(8),dimension(Nreal)           :: funct_out
      !
      integer                               :: Nfreq
      real(8),allocatable                   :: wreal(:),wmats(:)
      !
      !
      if(verbose)write(*,"(A)") "---- pade"
      !
      !
      Nfreq = Nmats
      if(present(wlimit))Nfreq = wlimit
      !
      allocate(wreal(Nreal));wreal=0d0
      wreal = linspace(-wrealMax,+wrealMax,Nreal)
      allocate(wmats(Nfreq));wmats=0d0
      select case(reg(type))
         case default
            stop "Available Modes are: Fermionic, Bosonic."
         case("Fermionic")
            wmats = FermionicFreqMesh(Beta,Nfreq)
         case("Bosonic")
            wmats = BosonicFreqMesh(Beta,Nfreq)
      end select
      !
      funct_out=czero
      call padecoeff(funct_out, wreal+img*eta, funct_in(1:Nfreq), img*wmats)
      !
   end function pade


   !---------------------------------------------------------------------------!
   !PURPOSE: Recursion for Pade coefficient (J.Serene)
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine padecoeff(fwout,wout,fwin,win)
      !
      use parameters
      implicit none
      !
      complex(8),intent(out)                :: fwout(:)
      complex(8),intent(in)                 :: wout(:)
      complex(8),intent(in)                 :: fwin(:)
      complex(8),intent(in)                 :: win(:)
      !
      complex(8),allocatable                :: coeff(:,:)
      complex(8),allocatable                :: a(:),b(:)
      integer                               :: Nin,Nout,i,j
      !
      !
      if(verbose)write(*,"(A)") "---- padecoeff"
      !
      !
      Nout = size(fwout)
      if(size(wout).ne.Nout) stop "padecoeff: size(wout).ne.Nout"
      !
      Nin = size(fwin)
      if(size(win).ne.Nin) stop "padecoeff: size(win).ne.Nin"
      !
      allocate(coeff(Nin,Nin));coeff=czero
      do j=1,Nin
         coeff(1,j)=fwin(j)
      enddo
      do j=2,Nin
         do i=2,j
            if (abs(win(j)-win(i-1)).lt.1.D-10) then
               stop "pade z=0"
            endif
            if (abs(coeff(i-1,j)).lt.1.D-10) then
               write(*,*)"i,j,coeff(i-1,j)",i,j,coeff(i-1,j)
               stop "coeff=0"
            endif
            coeff(i,j) = (coeff(i-1,i-1)-coeff(i-1,j)) / (win(j)-win(i-1)) / coeff(i-1,j)
         enddo
      enddo
      !
      allocate(a(0:Nin))
      allocate(b(0:Nin))
      do j=1,Nout
         !
         a=czero
         b=czero
         !
         a(0)=0.d0
         a(1)=coeff(1,1)
         b(0)=1.d0
         b(1)=1.d0
         !
         do i=1,Nin-1
            a(i+1)=a(i)+(wout(j)-win(i))*coeff(i+1,i+1)*a(i-1)
            b(i+1)=b(i)+(wout(j)-win(i))*coeff(i+1,i+1)*b(i-1)
         enddo
         !
         fwout(j) = a(Nin)/b(Nin)
         !
      enddo
      !
   end subroutine padecoeff


   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate to a new frequency mesh a Fermionic field
   !TEST ON:
   !---------------------------------------------------------------------------!
   subroutine interpolate2Beta_Fermionic(G,Beta_Match,mode,offDiag)
      !
      use parameters
      use utils_misc
      use utils_fields
      use interactions
      use input_vars, only: Nsite, SiteNorb, SiteOrbs
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      type(OldBeta),intent(in)              :: Beta_Match
      character(len=*),intent(in)           :: mode
      logical,intent(in)                    :: offDiag
      !
      type(FermionicField)                  :: G_old
      integer                               :: isite,ik,iw
      integer                               :: i,j,ispin
      integer                               :: i_lat,j_lat
      logical                               :: LocalOnly,replace
      real(8),allocatable                   :: wmats_new(:),wmats_old(:)
      real(8),allocatable                   :: ReGf(:),ImGf(:)
      real(8),allocatable                   :: ReD(:),ImD(:)
      !
      !
      if(verbose)write(*,"(A)") "---- interpolate2Beta_Fermionic"
      !
      !
      if(.not.G%status) stop "interpolate2Beta_Fermionic: FermionicField not properly initialized."
      if(G%Beta.ne.Beta_Match%Beta_old) stop "interpolate2Beta_Fermionic: Fermionic field Beta is different from the expected one."
      if(G%Npoints.ne.Beta_Match%Nmats_old) stop "interpolate2Beta_Fermionic: Fermionic field Npoints is different from the expected one."
      !
      LocalOnly=.true.
      if(G%Nkpt.ne.0)LocalOnly=.false.
      !
      allocate(wmats_new(Beta_Match%Nmats_new)); wmats_new=BosonicFreqMesh(Beta_Match%Beta_new,Beta_Match%Nmats_new)
      allocate(wmats_old(Beta_Match%Nmats_old)); wmats_old=BosonicFreqMesh(Beta_Match%Beta_old,Beta_Match%Nmats_old)
      !
      call duplicate(G_old,G)
      call DeallocateFermionicField(G)
      call AllocateFermionicField(G,G_old%Norb,Beta_Match%Nmats_new,Nkpt=G_old%Nkpt,Nsite=G_old%Nsite,Beta=Beta_Match%Beta_new,mu=G_old%mu)
      !
      select case(reg(mode))
         case default
            !
            stop "interpolate2Beta_Fermionic: Available Modes are: imp, lat."
            !
         case("imp")
            !
            allocate(ReD(Beta_Match%Nmats_old)) ;allocate(ImD(Beta_Match%Nmats_old))
            allocate(ReGf(Beta_Match%Nmats_new));allocate(ImGf(Beta_Match%Nmats_new))
            !
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(G,G_old,Nsite,SiteNorb,SiteOrbs),&
            !$OMP SHARED(Beta_Match,offDiag,wmats_new,wmats_old),&
            !$OMP PRIVATE(ispin,isite,i,j,i_lat,j_lat,replace,iw,ReD,ReGf,ImD,ImGf)
            !$OMP DO
            do isite=1,Nsite
               !
               do i=1,SiteNorb(isite)
                  do j=1,SiteNorb(isite)
                     !
                     i_lat = SiteOrbs(isite,i)
                     j_lat = SiteOrbs(isite,j)
                     !
                     replace = i_lat.eq.j_lat
                     if(offDiag) replace = .true.
                     !
                     if(replace)then
                        !
                        do ispin=1,Nspin
                           ReD=0d0;ImD=0d0;ReGf=0d0;ImGf=0d0
                           call nspline(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD)
                           call nspline(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD)
                           do iw=1,Beta_Match%Nmats_new
                              call splint(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD,wmats_new(iw),ReGf(iw))
                              call splint(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD,wmats_new(iw),ImGf(iw))
                           enddo
                           do iw=1,Beta_Match%Nmats_new
                              G%ws(i_lat,j_lat,iw,ispin) = dcmplx(ReGf(iw),ImGf(iw))
                           enddo
                        enddo
                        !
                     endif
                     !
                  enddo
               enddo
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(ReD,ImD,ReGf,ImGf)
            call DeallocateFermionicField(G_old)
            !
         case("lat")
            !
            allocate(ReD(Beta_Match%Nmats_old)) ;allocate(ImD(Beta_Match%Nmats_old))
            allocate(ReGf(Beta_Match%Nmats_new));allocate(ImGf(Beta_Match%Nmats_new))
            !
            !
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(G,G_old,LocalOnly),&
            !$OMP SHARED(Beta_Match,offDiag,wmats_new,wmats_old),&
            !$OMP PRIVATE(ispin,i,j,i_lat,j_lat,replace,iw,ik,ReD,ReGf,ImD,ImGf)
            !$OMP DO
            do i_lat=1,G%Norb
               do j_lat=1,G%Norb
                  !
                  replace = i_lat.eq.j_lat
                  if(offDiag) replace = .true.
                  !
                  if(replace)then
                     !
                     do ispin=1,Nspin
                        !
                        if(LocalOnly)then
                           !
                           ReD=0d0;ImD=0d0;ReGf=0d0;ImGf=0d0
                           call nspline(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD)
                           call nspline(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD)
                           do iw=1,Beta_Match%Nmats_new
                              call splint(wmats_old, real(G_old%ws(i_lat,j_lat,:,ispin)),ReD,wmats_new(iw),ReGf(iw))
                              call splint(wmats_old,aimag(G_old%ws(i_lat,j_lat,:,ispin)),ImD,wmats_new(iw),ImGf(iw))
                           enddo
                           do iw=1,Beta_Match%Nmats_new
                              G%ws(i_lat,j_lat,iw,ispin) = dcmplx(ReGf(iw),ImGf(iw))
                           enddo
                           !
                        else
                           !
                           do ik=1,G%Nkpt
                              ReD=0d0;ImD=0d0;ReGf=0d0;ImGf=0d0
                              call nspline(wmats_old, real(G_old%wks(i_lat,j_lat,:,ik,ispin)),ReD)
                              call nspline(wmats_old,aimag(G_old%wks(i_lat,j_lat,:,ik,ispin)),ImD)
                              do iw=1,Beta_Match%Nmats_new
                                 call splint(wmats_old, real(G_old%wks(i_lat,j_lat,:,ik,ispin)),ReD,wmats_new(iw),ReGf(iw))
                                 call splint(wmats_old,aimag(G_old%wks(i_lat,j_lat,:,ik,ispin)),ImD,wmats_new(iw),ImGf(iw))
                              enddo
                              do iw=1,Beta_Match%Nmats_new
                                 G%wks(i_lat,j_lat,iw,ik,ispin) = dcmplx(ReGf(iw),ImGf(iw))
                              enddo
                           enddo
                           !
                        endif
                        !
                     enddo
                     !
                  endif
                  !
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(ReD,ImD,ReGf,ImGf)
            call DeallocateFermionicField(G_old)
            if(.not.LocalOnly) call FermionicKsum(G)
            !
      end select
      !
   end subroutine interpolate2Beta_Fermionic
   !
   subroutine interpolate2Beta_Bosonic(W,Beta_Match,mode,offDiag)
      !
      use parameters
      use utils_misc
      use utils_fields
      use interactions
      use input_vars, only: Nsite, SiteNorb, SiteOrbs, UfullStructure
      implicit none
      !
      type(BosonicField),intent(inout)      :: W
      type(OldBeta),intent(in)              :: Beta_Match
      character(len=*),intent(in)           :: mode
      logical,intent(in)                    :: offDiag
      !
      type(BosonicField)                    :: W_old
      integer                               :: Norb
      integer                               :: isite,iq,ib1,ib2,iw
      integer                               :: i,j,k,l
      integer                               :: i_lat,j_lat,k_lat,l_lat
      logical                               :: LocalOnly,replace
      real(8),allocatable                   :: wmats_new(:),wmats_old(:)
      real(8),allocatable                   :: ReW(:),ImW(:)
      real(8),allocatable                   :: ReD(:),ImD(:)
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- interpolate2Beta_Bosonic"
      !
      !
      if(.not.W%status) stop "interpolate2Beta_Bosonic: BosonicField not properly initialized."
      if(W%Beta.ne.Beta_Match%Beta_old) stop "interpolate2Beta_Bosonic: Bosonic field Beta is different from the expected one."
      if(W%Npoints.ne.Beta_Match%Nmats_old) stop "interpolate2Beta_Bosonic: Bosonic field Npoints is different from the expected one."
      !
      Norb = int(sqrt(dble(W%Nbp)))
      !
      LocalOnly=.true.
      if(W%Nkpt.ne.0)LocalOnly=.false.
      !
      allocate(wmats_new(Beta_Match%Nmats_new)); wmats_new=BosonicFreqMesh(Beta_Match%Beta_new,Beta_Match%Nmats_new)
      allocate(wmats_old(Beta_Match%Nmats_old)); wmats_old=BosonicFreqMesh(Beta_Match%Beta_old,Beta_Match%Nmats_old)
      !
      call duplicate(W_old,W)
      call DeallocateBosonicField(W)
      call AllocateBosonicField(W,Norb,Beta_Match%Nmats_new,W_old%iq_gamma,Nkpt=W_old%Nkpt,Nsite=W_old%Nsite,no_bare=(.not.allocated(W_old%bare)),Beta=Beta_Match%Beta_new)
      !
      select case(reg(mode))
         case default
            !
            stop "interpolate2Beta_Bosonic: Available Modes are: imp, lat."
            !
         case("imp")
            !
            allocate(ReD(Beta_Match%Nmats_old))
            allocate(ReW(Beta_Match%Nmats_new))
            !
            call init_Uelements(Norb,PhysicalUelements)
            !
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(W,W_old,Nsite,SiteNorb,SiteOrbs,Norb),&
            !$OMP SHARED(Beta_Match,offDiag,UfullStructure,PhysicalUelements,wmats_new,wmats_old),&
            !$OMP PRIVATE(isite,i,j,k,l,i_lat,j_lat,k_lat,l_lat,ib1,ib2,replace,iw,ReD,ReW)
            !$OMP DO
            do isite=1,W%Nsite
               !
               do i=1,SiteNorb(isite)
                  do j=1,SiteNorb(isite)
                     do k=1,SiteNorb(isite)
                        do l=1,SiteNorb(isite)
                           !
                           ! mapping
                           i_lat = SiteOrbs(isite,i)
                           j_lat = SiteOrbs(isite,j)
                           k_lat = SiteOrbs(isite,k)
                           l_lat = SiteOrbs(isite,l)
                           !
                           ! bosonic indexes on the lattice
                           ib1 = i_lat + Norb*(j_lat-1)
                           ib2 = k_lat + Norb*(l_lat-1)
                           !
                           replace = PhysicalUelements%Full_Uaa(ib1,ib2) .or. PhysicalUelements%Full_Uab(ib1,ib2)
                           if(offDiag)then
                              replace = .true.
                              if(.not.UfullStructure) replace = PhysicalUelements%Full_All(ib1,ib2)
                           endif
                           !
                           if(replace)then
                              !
                              ReD=0d0;ReW=0d0
                              call nspline(wmats_old,real(W_old%screened_local(ib1,ib2,:)),ReD)
                              do iw=1,Beta_Match%Nmats_new
                                 call splint(wmats_old,real(W_old%screened_local(ib1,ib2,:)),ReD,wmats_new(iw),ReW(iw))
                              enddo
                              do iw=1,Beta_Match%Nmats_new
                                 W%screened_local(ib1,ib2,iw) = dcmplx(ReW(iw),0d0)
                              enddo
                              !
                           endif
                           !
                        enddo
                     enddo
                  enddo
               enddo
               !
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
            deallocate(ReW,ReD)
            call DeallocateBosonicField(W_old)
            !
         case("lat")
            !
            allocate(ReD(Beta_Match%Nmats_old));allocate(ImD(Beta_Match%Nmats_old))
            allocate(ReW(Beta_Match%Nmats_new));allocate(ImW(Beta_Match%Nmats_new))
            !
            !$OMP PARALLEL DEFAULT(NONE),&
            !$OMP SHARED(W,W_old,LocalOnly),&
            !$OMP SHARED(Beta_Match,offDiag,UfullStructure,PhysicalUelements,wmats_new,wmats_old),&
            !$OMP PRIVATE(ib1,ib2,replace,iw,iq,ReD,ReW,ImD,ImW)
            !$OMP DO
            do ib1=1,W%Nbp
               do ib2=1,W%Nbp
                  !
                  replace = PhysicalUelements%Full_Uaa(ib1,ib2) .or. PhysicalUelements%Full_Uab(ib1,ib2)
                  if(offDiag)then
                     replace = .true.
                     if(.not.UfullStructure) replace = PhysicalUelements%Full_All(ib1,ib2)
                  endif
                  !
                  if(replace)then
                     !
                     if(LocalOnly)then
                        !
                        ReD=0d0;ImD=0d0;ReW=0d0;ImW=0d0
                        call nspline(wmats_old, real(W_old%screened_local(ib1,ib2,:)),ReD)
                        call nspline(wmats_old,aimag(W_old%screened_local(ib1,ib2,:)),ImD)
                        do iw=1,Beta_Match%Nmats_new
                           call splint(wmats_old, real(W_old%screened_local(ib1,ib2,:)),ReD,wmats_new(iw),ReW(iw))
                           call splint(wmats_old,aimag(W_old%screened_local(ib1,ib2,:)),ImD,wmats_new(iw),ImW(iw))
                        enddo
                        do iw=1,Beta_Match%Nmats_new
                           W%screened_local(ib1,ib2,iw) = dcmplx(ReW(iw),ImW(iw))
                        enddo
                        !
                     else
                        !
                        do iq=1,W%Nkpt
                           ReD=0d0;ImD=0d0;ReW=0d0;ImW=0d0
                           call nspline(wmats_old, real(W_old%screened(ib1,ib2,:,iq)),ReD)
                           call nspline(wmats_old,aimag(W_old%screened(ib1,ib2,:,iq)),ImD)
                           do iw=1,Beta_Match%Nmats_new
                              call splint(wmats_old, real(W_old%screened(ib1,ib2,:,iq)),ReD,wmats_new(iw),ReW(iw))
                              call splint(wmats_old,aimag(W_old%screened(ib1,ib2,:,iq)),ImD,wmats_new(iw),ImW(iw))
                           enddo
                           do iw=1,Beta_Match%Nmats_new
                              W%screened(ib1,ib2,iw,iq) = dcmplx(ReW(iw),ImW(iw))
                           enddo
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
            deallocate(ReW,ReD,ImW,ImD)
            call DeallocateBosonicField(W_old)
            if(.not.LocalOnly) call BosonicKsum(W)
            !
      end select
      !
   end subroutine interpolate2Beta_Bosonic



   !---------------------------------------------------------------------------!
   !PURPOSE: Interpolate to a user provided K-point path a Fermionic field
   !---------------------------------------------------------------------------!
   subroutine interpolateG2Path(Sfull,Lttc,structure,Nkpt_path,pathOUTPUT)
      !
      use parameters
      use utils_misc
      use utils_fields
      use linalg, only : eigh, inv, zeye
      use crystal
      use file_io
      use greens_function, only : calc_Gmats
      use fourier_transforms
      use input_vars, only : Nreal, wrealMax, eta, path_funct, FermiSurf
      use input_vars, only : paramagnet, CalculationType
      implicit none
      !
      type(FermionicField),intent(inout)    :: Sfull
      type(Lattice),intent(inout)           :: Lttc
      character(len=*),intent(in)           :: structure
      integer,intent(in)                    :: Nkpt_path
      character(len=*),intent(in)           :: pathOUTPUT
      !
      type(FermionicField)                  :: Spath
      type(FermionicField)                  :: Sfermi
      type(FermionicField)                  :: Gpath
      type(FermionicField)                  :: Gfull
      type(FermionicField)                  :: Gfermi
      !
      real(8),allocatable                   :: wreal(:)
      complex(8),allocatable                :: zeta(:,:,:),invGf(:,:)
      real(8),allocatable                   :: Akw(:,:,:,:),Fk(:,:,:)
      real(8),allocatable                   :: Zk(:,:)
      integer                               :: Norb,Nmats,Ntau,unit
      integer                               :: ik,iw,itau,ispin,iorb
      integer                               :: ikx,iky
      integer                               :: Nkpt_Kside,wndx_cut
      character(len=256)                    :: path
      real                                  :: start,finish
      !
      !
      write(*,"(A)") new_line("A")//new_line("A")//"---- interpolateG2Path"
      !
      !
      ! Check on the input Fields
      if(.not.Lttc%status) stop "interpolateG2Path: Lttc not properly initialized."
      if(Sfull%status)then
         if(Sfull%Norb.ne.Lttc%Norb) stop "interpolateG2Path: Lttc has different number of orbitals with respect to Sfull."
         if(Sfull%Nkpt.ne.Lttc%Nkpt) stop "interpolateG2Path: Lttc has different number of k-points with respect to Sfull."
         Norb = Sfull%Norb
         Nmats = Sfull%Npoints
      else
         Norb = Lttc%Norb
      endif
      !
      !
      !-------------------- path along high-symmetry points -------------------!
      !
      !
      if(.not.Lttc%pathStored)then
         !
         !Create K-points along high-symmetry points
         if(allocated(Lttc%kptpath))deallocate(Lttc%kptpath)
         call calc_Kpath(Lttc%kptpath,reg(structure),Nkpt_path,Kaxis=Lttc%Kpathaxis,KaxisPoints=Lttc%KpathaxisPoints)
         Lttc%Nkpt_path = size(Lttc%kptpath,dim=2)
         !
         !Fill in Hk along points
         if(allocated(Lttc%Hk_path))deallocate(Lttc%Hk_path)
         allocate(Lttc%Hk_path(Norb,Norb,Lttc%Nkpt_path));Lttc%Hk_path=czero
         call cpu_time(start)
         call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath,Lttc%Hk,Lttc%Hk_path)
         call cpu_time(finish)
         write(*,"(A,F)") "     H(fullBZ) --> H(Kpath) cpu timing:", finish-start
         !
         !Fill in Ek along points
         if(allocated(Lttc%Ek_path))deallocate(Lttc%Ek_path)
         allocate(Lttc%Ek_path(Norb,Lttc%Nkpt_path));Lttc%Ek_path=0d0
         if(allocated(Lttc%Zk_path))deallocate(Lttc%Zk_path)
         allocate(Lttc%Zk_path(Norb,Norb,Lttc%Nkpt_path));Lttc%Zk_path=czero
         do ik=1,Lttc%Nkpt_path
            Lttc%Zk_path(:,:,ik) = Lttc%Hk_path(:,:,ik)
            call eigh(Lttc%Zk_path(:,:,ik),Lttc%Ek_path(:,ik))
         enddo
         !
      endif
      !
      !Re-Print bands in the same folder where the function is
      path = reg(pathOUTPUT)//"K_resolved/Bands.DAT"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      do ik=1,Lttc%Nkpt_path
         write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),(Lttc%Ek_path(:,ik),iorb=1,Norb)
      enddo
      close(unit)
      write(*,"(A,I)") "     Total number of K-points along path:",Lttc%Nkpt_path
      !
      !Re-Print position of High-symmetry points in the same folder where the function is
      path = reg(pathOUTPUT)//"K_resolved/Kpoints_path.DAT"
      unit = free_unit()
      open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
      do ik=1,size(Lttc%KpathaxisPoints,dim=1)
         write(unit,"(1I5,200E20.12)") ik,Lttc%KpathaxisPoints(ik)
      enddo
      close(unit)
      write(*,"(A,I)") "     Total number of High symmetry points:",size(Lttc%KpathaxisPoints,dim=1)
      !
      !Compute non-interacting spectral function in the Wannier basis
      allocate(wreal(Nreal));wreal=0d0
      wreal = linspace(-wrealMax,+wrealMax,Nreal)
      wndx_cut = minloc(abs(wreal-EcutSheet),dim=1)
      !
      allocate(zeta(Norb,Norb,Nreal));zeta=czero
      do iorb=1,Norb
         do iw=1,Nreal
            zeta(iorb,iorb,iw) = dcmplx(wreal(iw),eta)
         enddo
      enddo
      !
      allocate(Akw(Norb,Nreal,Lttc%Nkpt_path,Nspin));Akw=0d0
      allocate(invGf(Norb,Norb));invGf=czero
      do ispin=1,Nspin
         !$OMP PARALLEL DEFAULT(NONE),&
         !$OMP SHARED(ispin,Nreal,Norb,zeta,Lttc,Akw),&
         !$OMP PRIVATE(ik,iw,iorb,invGf)
         !$OMP DO
         do ik=1,Lttc%Nkpt_path
            do iw=1,Nreal
               !
               invGf = zeta(:,:,iw) - Lttc%Hk_path(:,:,ik)
               call inv(invGf)
               do iorb=1,Norb
                  Akw(iorb,iw,ik,ispin) = dimag(invGf(iorb,iorb))
               enddo
               !
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL
      enddo
      deallocate(zeta,invGf)
      !
      !Normalize spectral function
      do ispin=1,Nspin
         do ik=1,Lttc%Nkpt_path
            do iorb=1,Norb
               Akw(iorb,:,ik,ispin) = Akw(iorb,:,ik,ispin)/(sum(Akw(iorb,:,ik,ispin))*abs(wreal(2)-wreal(1)))
            enddo
         enddo
      enddo
      !
      !Print spectral function
      do ispin=1,Nspin
         path = reg(pathOUTPUT)//"K_resolved/Akw_nonInt_s"//str(ispin)//".DAT"
         unit = free_unit()
         open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Lttc%Nkpt_path
            do iw=1,Nreal
                write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),wreal(iw),(Akw(iorb,iw,ik,ispin),iorb=1,Norb)
            enddo
            write(unit,*)
         enddo
         close(unit)
         if(paramagnet)exit
      enddo
      deallocate(Akw,wreal)
      !
      !
      !------------------- path along planar sheet on kx,ky -------------------!
      !
      !
      if(FermiSurf)then
         !
         Nkpt_Kside = int(Nkpt_path/2)
         !
         if(.not.Lttc%planeStored)then
            !
            !Create K-points along high-symmetry points
            if(allocated(Lttc%kptPlane))deallocate(Lttc%kptPlane)
            call calc_Kplane(Lttc%kptPlane,Nkpt_Kside)
            Lttc%Nkpt_Plane = size(Lttc%kptPlane,dim=2)
            !
            !Fill in Hk in the points on the kx,ky plane - stored row-wise
            if(allocated(Lttc%Hk_Plane))deallocate(Lttc%Hk_Plane)
            allocate(Lttc%Hk_Plane(Norb,Norb,Lttc%Nkpt_Plane));Lttc%Hk_Plane=czero
            call cpu_time(start)
            call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,Lttc%Hk,Lttc%Hk_Plane)
            call cpu_time(finish)
            write(*,"(A,F)") "     H(fullBZ) --> H(kx,ky) cpu timing:", finish-start
            !
            Lttc%planeStored=.true.
            !
         endif
         write(*,"(A,I)") "     Total number of K-points along {kx,ky} plane:",Lttc%Nkpt_Plane
         !
         !Compute non-interacting Fermisurface
         allocate(Fk(Norb,Lttc%Nkpt_Plane,Nspin));Fk=0d0
         allocate(invGf(Norb,Norb));invGf=czero
         do ispin=1,Nspin
            do ik=1,Lttc%Nkpt_Plane
               !
               invGf = zeye(Norb)*dcmplx(EcutSheet,eta) - Lttc%Hk_Plane(:,:,ik)
               !
               call inv(invGf)
               do iorb=1,Norb
                  Fk(iorb,ik,ispin) = -dimag(invGf(iorb,iorb))
               enddo
               !
            enddo
         enddo
         deallocate(invGf)
         !
         !Print Fermi surface
         do ispin=1,Nspin
            path = reg(pathOUTPUT)//"K_resolved/Fk_nonInt_s"//str(ispin)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            do ik=1,Lttc%Nkpt_Plane
               ikx = int(ik/(Nkpt_Kside+0.001))+1
               iky = ik - (ikx-1)*Nkpt_Kside
               write(unit,"(3I5,200E20.12)") ik,ikx,iky,(Fk(iorb,ik,ispin),iorb=1,Norb)
               if(iky.eq.Nkpt_Kside)write(unit,*)
            enddo
            close(unit)
            if(paramagnet)exit
         enddo
         deallocate(Fk)
         !
      endif
      !
      !
      if(Sfull%status)then
         !
         !Interpolate the slef-energy along the path if its K-dependent otherwise duplicate the local one
         select case(reg(CalculationType))
            case default
               !
               stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
               !
            case("G0W0","scGW","GW+EDMFT")
               !
               !
               !
               !--------------- Green's function in the full BZ ---------------!
               call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
               !
               !
               !
               !--------------- Green's function along the path ---------------!
               !
               !Interpolate the self-energy along the path
               call cpu_time(start)
               call AllocateFermionicField(Spath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               spinloopSpath: do ispin=1,Nspin
                  call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptpath,Sfull%wks(:,:,:,:,ispin),Spath%wks(:,:,:,:,ispin))
                  if(paramagnet)then
                     Spath%wks(:,:,:,:,Nspin) = Spath%wks(:,:,:,:,1)
                     exit spinloopSpath
                  endif
               enddo spinloopSpath
               call cpu_time(finish)
               write(*,"(A,F)") new_line("A")//new_line("A")//"     Sigma(fullBZ,iw) --> Sigma(Kpath,iw) cpu timing:", finish-start
               !
               !Compute the quasiparticle weight along the path
               do ispin=1,Nspin
                  !
                  allocate(Zk(Lttc%Nkpt_path,Norb));Zk=0d0
                  do ik=1,Lttc%Nkpt_path
                     do iorb=1,Norb
                        Zk(ik,iorb) = 1d0 / (1d0 + abs(dimag(Spath%wks(iorb,iorb,1,ik,ispin)))*Spath%Beta/pi)
                     enddo
                  enddo
                  !
                  path = reg(pathOUTPUT)//"K_resolved/Zk_path_s"//str(ispin)//".DAT"
                  unit = free_unit()
                  open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
                  do ik=1,Lttc%Nkpt_path
                     write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),(Zk(ik,iorb),iorb=1,Norb)
                  enddo
                  close(unit)
                  deallocate(Zk)
                  !
                  if(paramagnet)exit
               enddo
               !
               !Recompute the Green's function
               call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               call calc_Gmats(Gpath,Lttc,Smats=Spath,along_path=.true.)
               call DeallocateFermionicField(Spath)
               !
               !
               !
               !---------- Green's function along the {kx,ky} plane -----------!
               if(FermiSurf)then
                  !
                  !Interpolate the self-energy along the plane
                  call cpu_time(start)
                  call AllocateFermionicField(Sfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
                  spinloopSplane: do ispin=1,Nspin
                     call wannierinterpolation(Lttc%Nkpt3,Lttc%kpt,Lttc%kptPlane,Sfull%wks(:,:,:,:,ispin),Sfermi%wks(:,:,:,:,ispin))
                     if(paramagnet)then
                        Sfermi%wks(:,:,:,:,Nspin) = Sfermi%wks(:,:,:,:,1)
                        exit spinloopSplane
                     endif
                  enddo spinloopSplane
                  call cpu_time(finish)
                  write(*,"(A,F)") "     Sigma(fullBZ,iw) --> Sigma(kx,ky,iw) cpu timing:", finish-start
                  !
                  !Compute the quasiparticle weight along the plane
                  do ispin=1,Nspin
                     !
                     allocate(Zk(Lttc%Nkpt_Plane,Norb));Zk=0d0
                     do ik=1,Lttc%Nkpt_Plane
                        do iorb=1,Norb
                           Zk(ik,iorb) = 1d0 / (1d0 + abs(dimag(Sfermi%wks(iorb,iorb,1,ik,ispin)))*Sfermi%Beta/pi)
                        enddo
                     enddo
                     !
                     path = reg(pathOUTPUT)//"K_resolved/Zk_plane_s"//str(ispin)//".DAT"
                     unit = free_unit()
                     open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
                     do ik=1,Lttc%Nkpt_Plane
                        ikx = int(ik/(Nkpt_Kside+0.001))+1
                        iky = ik - (ikx-1)*Nkpt_Kside
                        write(unit,"(3I5,200E20.12)") ik,ikx,iky,(Zk(ik,iorb),iorb=1,Norb)
                        if(iky.eq.Nkpt_Kside)write(unit,*)
                     enddo
                     close(unit)
                     deallocate(Zk)
                     !
                     if(paramagnet)exit
                  enddo
                  !
                  !Recompute the Green's function
                  call AllocateFermionicField(Gfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
                  call calc_Gmats(Gfermi,Lttc,Smats=Sfermi,along_plane=.true.)
                  call DeallocateFermionicField(Sfermi)
                  !
               endif
               !
               !
               !
            case("DMFT+statU","DMFT+dynU","EDMFT")
               !
               !
               !
               !--------------- Green's function in the full BZ ---------------!
               call AllocateFermionicField(Gfull,Norb,Nmats,Nkpt=Lttc%Nkpt,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               call calc_Gmats(Gfull,Lttc,Smats=Sfull,along_path=.false.)
               !
               !
               !
               !--------------- Green's function along the path ---------------!
               call AllocateFermionicField(Spath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               call AllocateFermionicField(Gpath,Norb,Nmats,Nkpt=Lttc%Nkpt_path,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
               do ik=1,Lttc%Nkpt_path
                  Spath%wks(:,:,:,ik,:) = Sfull%ws
               enddo
               call calc_Gmats(Gpath,Lttc,Smats=Spath,along_path=.true.)
               call DeallocateFermionicField(Spath)
               !
               !
               !
               !---------- Green's function along the {kx,ky} plane -----------!
               if(FermiSurf)then
                  call AllocateFermionicField(Sfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
                  call AllocateFermionicField(Gfermi,Norb,Nmats,Nkpt=Lttc%Nkpt_Plane,Nsite=Sfull%Nsite,Beta=Sfull%Beta,mu=Sfull%mu)
                  do ik=1,Lttc%Nkpt_Plane
                     Sfermi%wks(:,:,:,ik,:) = Sfull%ws
                  enddo
                  call calc_Gmats(Gfermi,Lttc,Smats=Sfermi,along_plane=.true.)
                  call DeallocateFermionicField(Sfermi)
               endif

               !
               !
               !
         end  select
         !
         !
         if(scan(reg(path_funct),"G").gt.0)then
            !
            !Dump K-resolved MaxEnt data along the path
            call calc_MaxEnt_on_G_K(Gpath,"path")
            !
            !Dump K-resolved MaxEnt data along the {kx,ky} plane
            if(FermiSurf)call calc_MaxEnt_on_G_K(Gfermi,"plane")
            !
         endif
         !
         !
         if(scan(reg(path_funct),"S").gt.0)then
            !
            select case(reg(CalculationType))
               case default
                  !
                  stop "Available Calculation types are: G0W0, scGW, DMFT+statU, DMFT+dynU, EDMFT, GW+EDMFT."
                  !
               case("G0W0","scGW","GW+EDMFT")
                  !
                  !Dump K-resolved MaxEnt data in the full BZ
                  call calc_MaxEnt_on_Sigma_K(Gfull,"full")
                  !
                  !Dump K-resolved MaxEnt data along the path
                  call calc_MaxEnt_on_Sigma_K(Gpath,"path")
                  !
                  !Dump K-resolved MaxEnt data along the {kx,ky} plane
                  call calc_MaxEnt_on_Sigma_K(Gfermi,"plane")
                  !
               case("DMFT+statU","DMFT+dynU","EDMFT")
                  !
                  !Dump MaxEnt data for the local self-energy
                  call calc_MaxEnt_on_Sigma_imp(Sfull)
                  !
            end  select
            !
         endif
         !
         call DeallocateFermionicField(Gpath)
         call DeallocateFermionicField(Gfull)
         call DeallocateFermionicField(Gfermi)
         !
      endif
      !
      !
   contains
      !
      !
      !
      subroutine calc_MaxEnt_on_G_K(Gmats_in,mode)
         !
         use input_vars, only : Solver
         implicit none
         !
         type(FermionicField),intent(in)       :: Gmats_in
         character(len=*),intent(in)           :: mode
         !
         complex(8),allocatable                :: Gmats_diag(:,:,:,:),Gitau_diag(:,:,:,:)
         real(8),allocatable                   :: Ak(:,:)
         real(8),allocatable                   :: tau(:)
         integer                               :: Ntau,Nkpt
         integer                               :: ikx,iky
         !
         !
         if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- calc_MaxEnt_on_G_K"
         !
         !
         if(.not.Gmats_in%status) stop "calc_MaxEnt_on_G_K: Gmats_in not properly allocated."
         select case(reg(mode))
            case default
               !
               stop "calc_MaxEnt_on_G_K: Available Modes are: path, full, plane."
               !
            case("path")
               !
               Nkpt = Lttc%Nkpt_path
               write(*,"(A,I)") "     G path. Total number of K-points along path:",Nkpt
               !
            case("plane")
               !
               Nkpt = Lttc%Nkpt_Plane
               write(*,"(A,I)") "     G plane. Total number of K-points in the {kx,ky} sheet:",Nkpt
               !
         end select
         !
         !Extract the diagonal of the Green's function
         allocate(Gmats_diag(Norb,Nmats,Nkpt,Nspin));Gmats_diag=czero
         do ispin=1,Nspin
            do ik=1,Nkpt
               do iw=1,Nmats
                  do iorb=1,Norb
                     Gmats_diag(iorb,iw,ik,ispin) = Gmats_in%wks(iorb,iorb,iw,ik,ispin)
                  enddo
               enddo
            enddo
         enddo
         !
         !Fourier transform the diagonal of the Green's function
         Ntau = Solver%NtauF
         call cpu_time(start)
         allocate(Gitau_diag(Norb,Ntau,Nkpt,Nspin));Gitau_diag=czero
         spinloopGftP: do ispin=1,Nspin
            call Fmats2itau_vec(Sfull%Beta,Gmats_diag(:,:,:,ispin),Gitau_diag(:,:,:,ispin), &
            asympt_corr=.true.,tau_uniform=.true.)
            if(paramagnet)then
               Gitau_diag(:,:,:,Nspin) = Gitau_diag(:,:,:,1)
               exit spinloopGftP
            endif
         enddo spinloopGftP
         deallocate(Gmats_diag)
         call cpu_time(finish)
         write(*,"(A,F)") "     Glat(K"//reg(mode)//",iw) --> Glat(K"//reg(mode)//",tau) cpu timing:", finish-start
         !
         !Print data for K-resolved MaxEnt
         allocate(tau(Ntau));tau = linspace(0d0,Sfull%Beta,Ntau)
         do ispin=1,Nspin
           do ik=1,Nkpt
               !
               path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Gk_"//reg(mode)//"_t_s"//str(ispin)//"/Gk_t_k"//str(ik)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do itau=1,Ntau
                   write(unit,"(200E20.12)") tau(itau),(dreal(Gitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
               enddo
               close(unit)
               !
           enddo
           if(paramagnet)exit
         enddo
         deallocate(tau)
         !
         !Compute the spectral weight at Fermi along the path. See arxiv:0805.3778 Eq.(5)
         do ispin=1,Nspin
            !
            allocate(Ak(Nkpt,Norb));Ak=0d0
            do ik=1,Nkpt
               do iorb=1,Norb
                  Ak(ik,iorb) = -dreal(Gitau_diag(iorb,int(Ntau/2),ik,ispin))*Gmats_in%Beta
               enddo
            enddo
            !
            path = reg(pathOUTPUT)//"K_resolved/Ak_"//reg(mode)//"_s"//str(ispin)//".DAT"
            unit = free_unit()
            open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
            if(reg(mode).eq."path")then
               do ik=1,Nkpt
                  write(unit,"(1I5,200E20.12)") ik,Lttc%Kpathaxis(ik),(Ak(ik,iorb),iorb=1,Norb)
               enddo
            elseif(reg(mode).eq."plane")then
               do ik=1,Nkpt
                  ikx = int(ik/(Nkpt_Kside+0.001))+1
                  iky = ik - (ikx-1)*Nkpt_Kside
                  write(unit,"(3I5,200E20.12)") ik,ikx,iky,(Ak(ik,iorb),iorb=1,Norb)
                  if(iky.eq.Nkpt_Kside)write(unit,*)
               enddo
            endif
            close(unit)
            deallocate(Ak)
            !
            if(paramagnet)exit
         enddo
         deallocate(Gitau_diag)
         !
      end subroutine calc_MaxEnt_on_G_K
      !
      !
      !
      subroutine calc_MaxEnt_on_Sigma_K(Gmats_in,mode)
         !
         use linalg, only : diagonal, rotate
         use input_vars, only : Solver, ReplaceTail_Simp
         implicit none
         !
         type(FermionicField),intent(in)       :: Gmats_in
         character(len=*),intent(in)           :: mode
         !
         real(8),allocatable                   :: wmats(:),Sparams(:,:,:,:)
         complex(8),allocatable                :: Gmats_rho(:,:,:,:)
         complex(8),allocatable                :: Smats_diag(:,:,:,:),Sitau_diag(:,:,:,:)
         complex(8),allocatable                :: RotN(:,:,:,:)
         real(8),allocatable                   :: EigN(:,:,:)
         real(8),allocatable                   :: tau(:)
         character(len=256)                    :: ParaFile
         real(8)                               :: M0,M1,M2,M3,M4
         real(8)                               :: x1,x2,x3,y1,y2,y3
         integer                               :: unit,wndx,iw1,iw2,iw3
         integer                               :: Nkpt
         !
         !
         if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- calc_MaxEnt_on_Sigma_K"
         !
         !
         if(.not.Gmats_in%status) stop "calc_MaxEnt_on_Sigma_K: Gmats_in not properly allocated."
         select case(reg(mode))
            case default
               !
               stop "calc_MaxEnt_on_Sigma_K: Available Modes are: path, full, plane."
               !
            case("full")
               !
               Nkpt = Lttc%Nkpt
               write(*,"(A,I)") "     Sigma full. Total number of K-points in the BZ:",Nkpt
               !
            case("path")
               !
               Nkpt = Lttc%Nkpt_path
               write(*,"(A,I)") "     Sigma path. Total number of K-points along path:",Nkpt
               !
            case("plane")
               !
               Nkpt = Lttc%Nkpt_Plane
               write(*,"(A,I)") "     Sigma plane. Total number of K-points in the {kx,ky} sheet:",Nkpt
               !
         end select
         !
         allocate(Sparams(Norb,Nkpt,Nspin,2));Sparams=0d0
         allocate(wmats(Nmats));wmats=FermionicFreqMesh(Sfull%Beta,Nmats)
         wndx = minloc(abs(wmats-ReplaceTail_Simp),dim=1)
         !
         !Find the basis where the K-dependent density matrix is diagoal on the path
         allocate(RotN(Norb,Norb,Nkpt,Nspin));RotN=czero
         allocate(EigN(Norb,Nkpt,Nspin));EigN=0d0
         do ik=1,Nkpt
            do ispin=1,Nspin
               RotN(:,:,ik,ispin) = dreal(Gmats_in%N_ks(:,:,ik,ispin))
               call eigh(RotN(:,:,ik,ispin),EigN(:,ik,ispin))
            enddo
         enddo
         !
         !Operations at each K-point
         allocate(Smats_diag(Norb,Nmats,Nkpt,Nspin));Smats_diag=czero
         do ik=1,Nkpt
            !
            !
            !Bring the Green's function on the path to basis where the K-dependent density matrix is diagonal
            allocate(Gmats_rho(Norb,Norb,Nmats,Nspin));Gmats_rho=czero
            do ispin=1,Nspin
               do iw=1,Nmats
                  Gmats_rho(:,:,iw,ispin) = rotate(Gmats_in%wks(:,:,iw,ik,ispin),RotN(:,:,ik,ispin))
               enddo
            enddo
            if(paramagnet) Gmats_rho(:,:,:,Nspin) = Gmats_rho(:,:,:,1)
            !
            !Check Print - This is to check if the rotation has brought the Gf to a diagonal basis
            do ispin=1,Nspin
               path = reg(pathOUTPUT)//"K_resolved/Gk_"//reg(mode)//"_wm_s"//str(ispin)//"/Gk_wm_rho_k"//str(ik)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats
                  write(unit,"(200E20.12)") wmats(iw),(Gmats_rho(iorb,iorb,iw,ispin),iorb=1,Norb),(Gmats_rho(1,iorb,iw,ispin),iorb=2,Norb)
               enddo
               close(unit)
               if(paramagnet)exit
            enddo
            !
            !Reference frequencies
            iw1 = wndx-30   ; x1 = wmats(iw1)
            iw2 = wndx-1    ; x2 = wmats(iw2)
            iw3 = Nmats     ; x3 = wmats(iw3) !Never really used
            !
            !Replace the tail of the diagonal Green's function
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  !Real part asymptotic behaviour M0=0, M2, M4
                  y1 = dreal(Gmats_rho(iorb,iorb,iw1,ispin))
                  y2 = dreal(Gmats_rho(iorb,iorb,iw2,ispin))
                  M4 = (y1*x1**2 - y2*x2**2) * ( x1**2 * x2**2 )/( x2**2 - x1**2 )
                  M4 = 0d0
                  M2 = y2*x2**2 - M4/x2**2
                  !
                  !Imaginary part asymptotic behaviour M1=1, M3>0 by definition
                  y2 = dimag(Gmats_rho(iorb,iorb,iw2,ispin))
                  M3 = y2*x2**3 + x2**2
                  !
                  !Replace the tail of the diagoal Green's function
                  do iw=iw2,Nmats
                     Gmats_rho(iorb,iorb,iw,ispin) = dcmplx( M2/(wmats(iw)**2)+M4/(wmats(iw)**4) , -1d0/wmats(iw)+M3/(wmats(iw)**3) )
                  enddo
                  !
                  !The corresponding M0 and M1 of the self-energy
                  Sparams(iorb,ik,ispin,1) = Sfull%mu - M2
                  Sparams(iorb,ik,ispin,2) = M2**2 - M3
                  !
               enddo
            enddo
            !
            !Check Print - This is to check if the tail replacement has any problem
            do ispin=1,Nspin
               path = reg(pathOUTPUT)//"K_resolved/Gk_"//reg(mode)//"_wm_s"//str(ispin)//"/Gk_wm_rho_tail_k"//str(ik)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats
                  write(unit,"(200E20.12)") wmats(iw),(Gmats_rho(iorb,iorb,iw,ispin),iorb=1,Norb)
               enddo
               close(unit)
               if(paramagnet)exit
            enddo
            !
            !Get the self-energy on the path in the basis where the K-dependent density matrix is diagonal
            !Here H(k) is enclosed inside the self-energy
            do ispin=1,Nspin
               do iorb=1,Norb
                  do iw=1,Nmats
                     Smats_diag(iorb,iw,ik,ispin) = img*wmats(iw) + Sfull%mu - 1d0/Gmats_rho(iorb,iorb,iw,ispin)
                  enddo
               enddo
            enddo
            deallocate(Gmats_rho)
            !
            !Check Print - This is to check the self-energy in the diagonal basis
            do ispin=1,Nspin
               path = reg(pathOUTPUT)//"K_resolved/Sk_"//reg(mode)//"_wm_s"//str(ispin)//"/Sk_wm_rho_tail_k"//str(ik)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats
                  write(unit,"(200E20.12)") wmats(iw),(Smats_diag(iorb,iw,ik,ispin),iorb=1,Norb)
               enddo
               close(unit)
               if(paramagnet)exit
            enddo
            !
            !Extract the rescaling parameters from the Self-energy just for comparison
            do ispin=1,Nspin
               do iorb=1,Norb
                  !
                  !Real part asymptotic behaviour M0, M2
                  y1 = dreal(Smats_diag(iorb,iw1,ik,ispin))
                  y2 = dreal(Smats_diag(iorb,iw2,ik,ispin))
                  y3 = dreal(Smats_diag(iorb,iw3,ik,ispin))
                  !M0 = get_Meven(y1,y2,y3,x1,x2,x3,0)
                  !M2 = get_Meven(y1,y2,y3,x1,x2,x3,2)
                  !M4 = get_Meven(y1,y2,y3,x1,x2,x3,4)
                  M2 = (y1-y2)*(x1**2 * x2**2)/(x2**2 - x1**2)
                  M0 = y2 - M2/x2**2
                  !
                  !Imaginary part asymptotic behaviour M1
                  y2 = dimag(Smats_diag(iorb,iw2,ik,ispin))
                  !M1 = get_Modd(y1,y2,x1,x2,1)
                  !M3 = get_Modd(y1,y2,x1,x2,3)
                  M1 = y2*x2
                  !
                  !Print comparison
                  if(verbose)then
                     write(*,"(5X,2(A8,I5))")"ik=",ik,"  iorb=",iorb
                     write(*,"(5X,2(A12,1F20.8))")"M0(1/G): ",Sparams(iorb,ik,ispin,1),"M0(S): ",M0
                     write(*,"(5X,2(A12,1F20.8))")"M1(1/G): ",Sparams(iorb,ik,ispin,2),"M1(S): ",M1
                  endif
                  !
                  !Remove the bare limit M0
                  Smats_diag(iorb,:,ik,ispin) = Smats_diag(iorb,:,ik,ispin) - Sparams(iorb,ik,ispin,1)
                  !
                  !Rescale so as to have -ImS~1/iw
                  Smats_diag(iorb,:,ik,ispin) = Smats_diag(iorb,:,ik,ispin)/abs(Sparams(iorb,ik,ispin,2))
                  !
                  !Revert the real part
                  Smats_diag(iorb,:,ik,ispin) = -conjg(Smats_diag(iorb,:,ik,ispin))
                  !
               enddo
            enddo
            !
            !Check Print - This is to check the self-energy in the diagonal basis afer the removal of M0 and rescaling
            do ispin=1,Nspin
               path = reg(pathOUTPUT)//"K_resolved/Sk_"//reg(mode)//"_wm_s"//str(ispin)//"/Sk_wm_rho_tail_rescaled_k"//str(ik)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do iw=1,Nmats
                  write(unit,"(200E20.12)") wmats(iw),(Smats_diag(iorb,iw,ik,ispin),iorb=1,Norb)
               enddo
               close(unit)
               if(paramagnet)exit
            enddo
            !
            !
         enddo
         deallocate(wmats)
         !
         !Fourier transform
         Ntau = 300
         call cpu_time(start)
         allocate(Sitau_diag(Norb,Ntau,Nkpt,Nspin));Sitau_diag=czero
         spinloopSft: do ispin=1,Nspin
            call Fmats2itau_vec(Sfull%Beta,Smats_diag(:,:,:,ispin),Sitau_diag(:,:,:,ispin), &
            asympt_corr=.true.,tau_uniform=.true.)
            if(paramagnet)then
               Sitau_diag(:,:,:,Nspin) = Sitau_diag(:,:,:,1)
               exit spinloopSft
            endif
         enddo spinloopSft
         deallocate(Smats_diag)
         call cpu_time(finish)
         write(*,"(A,F)") "     Slat(K"//reg(mode)//",iw) --> Slat(K"//reg(mode)//",tau) cpu timing:", finish-start
         !
         !Final correction
         do ik=1,Nkpt
            do iorb=1,Norb
               if(dreal(Sitau_diag(iorb,1,ik,1)).gt.0d0) then
                  write(*,"(A)")"     Warning: orbital# "//str(iorb)//" of K-point # "//str(ik)//" is positive in tau=0"
               elseif(dreal(Sitau_diag(iorb,Ntau,ik,1)).gt.0d0)then
                  write(*,"(A)")"     Warning: orbital# "//str(iorb)//" of K-point # "//str(ik)//" is positive in tau=beta"
               else
                  write(*,"(A)")"     Orbital# "//str(iorb)//" of K-point # "//str(ik)//" is fine. Reverting positive noise."
                  do ispin=1,Nspin
                     Sitau_diag(iorb,:,ik,ispin) = -abs(Sitau_diag(iorb,:,ik,ispin))
                  enddo
               endif
            enddo
         enddo
         !
         !Print data for K-resolved MaxEnt
         allocate(tau(Ntau));tau = linspace(0d0,Sfull%Beta,Ntau)
         do ispin=1,Nspin
            do ik=1,Nkpt
               !
               path = reg(pathOUTPUT)//"K_resolved/MaxEnt_Sk_"//reg(mode)//"_t_s"//str(ispin)//"/Sk_t_k"//str(ik)//".DAT"
               unit = free_unit()
               open(unit,file=reg(path),form="formatted",status="unknown",position="rewind",action="write")
               do itau=1,Ntau
                  write(unit,"(200E20.12)") tau(itau),(dreal(Sitau_diag(iorb,itau,ik,ispin)),iorb=1,Norb)
               enddo
               close(unit)
               !
            enddo
            if(paramagnet)exit
         enddo
         deallocate(Sitau_diag,tau)
         !
         !Print data needed to reconstruct the Self-energy
         ParaFile = reg(pathOUTPUT)//"K_resolved/Sigma_vars/S"//reg(mode)//"_Params.DAT"
         unit = free_unit()
         open(unit,file=reg(ParaFile),form="formatted",status="unknown",position="rewind",action="write")
         do ik=1,Nkpt
            !
            !rotation
            do ispin=1,Nspin
               path = reg(pathOUTPUT)//"K_resolved/Sigma_vars/S"//reg(mode)//"_Rot_k"//str(ik)//"_s"//str(ispin)//".DAT"
               call dump_Matrix(RotN(:,:,ik,ispin),reg(path))
               if(paramagnet)exit
            enddo
            !
            !Parameters
            write(unit,"(200E20.12)") (Sparams(iorb,ik,1,1),iorb=1,Norb),(Sparams(iorb,ik,1,2),iorb=1,Norb),         &
                                      (Sparams(iorb,ik,Nspin,1),iorb=1,Norb),(Sparams(iorb,ik,Nspin,2),iorb=1,Norb)
            !
         enddo
         deallocate(Sparams,RotN,EigN)
         close(unit)
         !
      end subroutine calc_MaxEnt_on_Sigma_K
      !
      !
      !
      subroutine calc_MaxEnt_on_Sigma_imp(Smats_in)
         !
         use input_vars, only : ReplaceTail_Simp, PadeWlimit, Solver
         use input_vars, only : SiteNorb, SiteOrbs, SiteName, Nsite, EqvGWndx
         use input_vars, only : OlocSite, OlocRot, OlocRotDag, OlocEig
         use input_vars, only : RotateHloc, ExpandImpurity, AFMselfcons
         implicit none
         !
         type(FermionicField),intent(in)       :: Smats_in
         !
         type(FermionicField)                  :: Simp
         real(8),allocatable                   :: wmats(:),Sparams(:,:,:)
         complex(8),allocatable                :: Smats_diag(:,:,:)
         complex(8),allocatable                :: Sitau_diag(:,:,:)
         complex(8),allocatable                :: Rot(:,:)
         integer,allocatable                   :: Orbs(:)
         real(8)                               :: M0,M1,M2,M3,M4
         real(8)                               :: x1,x2,x3,y1,y2,y3
         integer                               :: unit,wndx,iw1,iw2,iw3
         integer                               :: isite
         character(len=256)                    :: ParaFile
         !
         !
         if(verbose) write(*,"(A)") new_line("A")//new_line("A")//"---- calc_MaxEnt_on_Sigma_imp"
         !
         !
         if(.not.Smats_in%status) stop "calc_MaxEnt_on_Sigma_imp: Smats_in not properly allocated."
         !
         allocate(wmats(Nmats));wmats=FermionicFreqMesh(Smats_in%Beta,Nmats)
         wndx = minloc(abs(wmats-ReplaceTail_Simp),dim=1)
         !
         !Operations at each site
         do isite=1,Nsite
            !
            allocate(Orbs(SiteNorb(isite)))
            Orbs = SiteOrbs(isite,1:SiteNorb(isite))
            !
            allocate(Smats_diag(SiteNorb(isite),Nmats,Nspin));Smats_diag=czero
            allocate(Sparams(SiteNorb(isite),Nspin,2));Sparams=0d0
            !
            !Get the irreducible local self-energy
            call AllocateFermionicField(Simp,SiteNorb(isite),Nmats,Beta=Smats_in%Beta)
            if(RotateHloc)then
               !
               allocate(Rot(SiteNorb(isite),SiteNorb(isite)))
               Rot=OlocRot(1:SiteNorb(isite),1:SiteNorb(isite),isite)
               call loc2imp(Simp,Smats_in,Orbs,U=Rot)
               !
            else
               !
               call loc2imp(Simp,Smats_in,Orbs)
               !
            endif
            !
            !I need a standard array for the FT
            do ispin=1,Nspin
               do iorb=1,SiteNorb(isite)
                  Smats_diag(iorb,:,ispin) = Simp%ws(iorb,iorb,:,ispin)
               enddo
            enddo
            !
            !Check Print - This has to be identical to the one in the Solver_* folder
            call dump_MaxEnt(Smats_diag,"mats",reg(pathOUTPUT)//"Convergence/","Sqmc_"//reg(SiteName(isite)))
            !
            !Reference frequencies
            iw1 = wndx-30   ; x1 = wmats(iw1)
            iw2 = wndx-1    ; x2 = wmats(iw2)
            iw3 = Nmats     ; x3 = wmats(iw3) !Never really used
            !
            !Extract the rescaling parameters from the Self-energy
            do ispin=1,Nspin
               do iorb=1,SiteNorb(isite)
                  !
                  !Real part asymptotic behaviour M0, M2
                  y1 = dreal(Smats_diag(iorb,iw1,ispin))
                  y2 = dreal(Smats_diag(iorb,iw2,ispin))
                  y3 = dreal(Smats_diag(iorb,iw3,ispin))
                  M0 = get_Meven(y1,y2,y3,x1,x2,x3,0)
                  !M2 = get_Meven(y1,y2,y3,x1,x2,x3,2)
                  !M4 = get_Meven(y1,y2,y3,x1,x2,x3,4)
                  M2 = (y1-y2)*(x1**2 * x2**2)/(x2**2 - x1**2)
                  !M0 = y2 - M2/x2**2
                  !
                  !Imaginary part asymptotic behaviour M1
                  y2 = dimag(Smats_diag(iorb,iw2,ispin))
                  !M1 = get_Modd(y1,y2,x1,x2,1)
                  !M3 = get_Modd(y1,y2,x1,x2,3)
                  M1 = y2*x2
                  !
                  Sparams(iorb,ispin,1) = M0
                  Sparams(iorb,ispin,2) = M1
                  !
                  !------------------------------------------------------------!
                  !SPbca - EDMFT correction
                  if(iorb.eq.1)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.197
                  if(iorb.eq.2)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.1
                  if(iorb.eq.3)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.1
                  !LPbca - EDMFT correction
                  !if(iorb.eq.1)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.213
                  !if(iorb.eq.2)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.15
                  !if(iorb.eq.3)Sparams(iorb,ispin,1)=Sparams(iorb,ispin,1)-0.25
                  !------------------------------------------------------------!
                  !
                  !Print comparison
                  !if(verbose)then
                     write(*,"(5X,A8,I5,A)")"isite=",isite,"  Element: "//reg(SiteName(isite))
                     write(*,"(5X,2(A8,I5))")"iorb=",iorb,"ispin=",ispin
                     write(*,"(5X,A12,1F20.8)")"M0(S): ",Sparams(iorb,ispin,1)
                     write(*,"(5X,A12,1F20.8)")"M1(S): ",Sparams(iorb,ispin,2)
                  !endif
                  !
                  !Remove the bare limit M0
                  Smats_diag(iorb,:,ispin) = Smats_diag(iorb,:,ispin) - Sparams(iorb,ispin,1)
                  !
                  !Rescale so as to have -ImS~1/iw
                  Smats_diag(iorb,:,ispin) = Smats_diag(iorb,:,ispin)/abs(Sparams(iorb,ispin,2))
                  !
                  !Revert the real part
                  Smats_diag(iorb,:,ispin) = -conjg(Smats_diag(iorb,:,ispin))
                  !
               enddo
            enddo
            !
            !Check Print - This is to check the self-energy in the diagonal basis afer the removal of M0 and rescaling
            call dump_MaxEnt(Smats_diag,"mats",reg(pathOUTPUT)//"Convergence/","Sqmc_rescaled_"//reg(SiteName(isite)))
            !
            deallocate(wmats)
            !
            !Fourier transform
            Ntau = 300
            call cpu_time(start)
            allocate(Sitau_diag(SiteNorb(isite),Ntau,Nspin));Sitau_diag=czero
            spinloopSft: do ispin=1,Nspin
               call Fmats2itau_vec(Sfull%Beta,Smats_diag(:,:,ispin),Sitau_diag(:,:,ispin), &
               asympt_corr=.true.,tau_uniform=.true.)
               if(paramagnet)then
                  Sitau_diag(:,:,Nspin) = Sitau_diag(:,:,1)
                  exit spinloopSft
               endif
            enddo spinloopSft
            deallocate(Smats_diag)
            call cpu_time(finish)
            write(*,"(A,F)") "     Sqmc_"//reg(SiteName(isite))//"(iw) --> Sqmc_"//reg(SiteName(isite))//"(tau) cpu timing:", finish-start
            !
            !Final correction
            do iorb=1,SiteNorb(isite)
               !
               if(dreal(Sitau_diag(iorb,1,1)).gt.0d0) then
                  write(*,"(A)")"     Warning: orbital# "//str(iorb)//" is positive in tau=0"
               else
                  write(*,"(A)")"     Orbital# "//str(iorb)//" is fine."
                  write(*,"(A)")"     Reverting positive noise."
                  do ispin=1,Nspin
                     Sitau_diag(iorb,:,ispin) = -abs(Sitau_diag(iorb,:,ispin))
                     Sitau_diag(iorb,Ntau,ispin) = -1d0-Sitau_diag(iorb,1,ispin)
                  enddo
               endif
               !
            enddo
            !
            !Print data for MaxEnt
            call dump_MaxEnt(Sitau_diag,"itau",reg(pathOUTPUT)//"Convergence/","Sqmc_"//reg(SiteName(isite)))
            deallocate(Sitau_diag)
            !
            !Print data needed to reconstruct the Self-energy
            ParaFile = reg(pathOUTPUT)//"K_resolved/Sigma_vars/Sqmc_"//reg(SiteName(isite))//"_Params.DAT"
            unit = free_unit()
            open(unit,file=reg(ParaFile),form="formatted",status="unknown",position="rewind",action="write")
            write(unit,"(200E20.12)") (Sparams(iorb,1,1),iorb=1,SiteNorb(isite)),(Sparams(iorb,1,2),iorb=1,SiteNorb(isite)),         & !spin=1
                                      (Sparams(iorb,Nspin,1),iorb=1,SiteNorb(isite)),(Sparams(iorb,Nspin,2),iorb=1,SiteNorb(isite))    !spin=2
            deallocate(Sparams)
            close(unit)
            !
            deallocate(Orbs)
            if(ExpandImpurity.or.AFMselfcons)exit
            !
         enddo
         !
      end subroutine calc_MaxEnt_on_Sigma_imp
      !
      !
      !
      function get_Meven(b1,b2,b3,w1,w2,w3,coeff) result(M)
         use linalg, only : det
         implicit none
         real(8),intent(in)      :: b1,b2,b3,w1,w2,w3
         integer,intent(in)      :: coeff
         real(8)                 :: M
         !
         real(8)                 :: num(3,3),den(3,3)
         !
         den(:,1)=1d0
         den(:,2)=[1d0/w1**2,1d0/w2**2,1d0/w3**2]
         den(:,3)=[1d0/w1**4,1d0/w2**4,1d0/w3**4]
         if(det(den).eq.0d0)stop"get_Meven"
         if(coeff.gt.4)stop"wrong even coeff"
         !
         num=den
         num(:,1+coeff/2)=[b1,b2,b3]
         !
         M = det(num) / det(den)
         !
      end function get_Meven
      !
      function get_Modd(b1,b2,w1,w2,coeff) result(M)
         use linalg, only : det
         implicit none
         real(8),intent(in)      :: b1,b2,w1,w2
         integer,intent(in)      :: coeff
         real(8)                 :: M
         !
         real(8)                 :: num(2,2),den(2,2)
         !
         den(:,1)=[1d0/w1   ,1d0/w2   ]
         den(:,2)=[1d0/w1**3,1d0/w2**3]
         if(det(den).eq.0d0)stop"get_Modd"
         if(coeff.gt.3)stop"wrong odd coeff"
         !
         num=den
         if(coeff.eq.1)num(:,1)=[b1,b2]
         if(coeff.eq.3)num(:,2)=[b1,b2]
         !num(:,1+(coeff-1)/2)=[b1,b2]
         !
         M = det(num) / det(den)
         !
      end function get_Modd
      !
      !
      !
   end subroutine interpolateG2Path



end module post_processing
