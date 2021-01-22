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
      use input_vars, only: Nsite, SiteNorb, SiteOrbs
      implicit none
      !
      type(BosonicField),intent(inout)      :: W
      character(len=*),intent(in)           :: mode
      integer,intent(in),optional           :: site
      !
      integer                               :: Norb,Nmats
      integer                               :: isite,ib1,ib2
      integer                               :: i,j,k,l
      integer                               :: i_lat,j_lat,k_lat,l_lat
      logical                               :: LocalOnly,replace
      real(8)                               :: ReW,ImW
      real(8),allocatable                   :: wmats(:)
      type(physicalU)                       :: PhysicalUelements
      !
      !
      if(verbose)write(*,"(A)") "---- remove_CDW"
      !
      !
      if(.not.W%status) stop "remove_CDW. BosonicField not properly initialized."
      Norb = int(sqrt(dble(W%Nbp)))
      Nmats = W%Npoints
      !
      LocalOnly=.true.
      if(W%Nkpt.ne.0)LocalOnly=.false.
      !
      allocate(wmats(Nmats))
      wmats=BosonicFreqMesh(W%Beta,Nmats)
      !
      select case(reg(mode))
         case default
            !
            stop "remove_CDW. Available Modes are: imp, lat."
            !
         case("imp")
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
         case("lat")
            !
            do ib1=1,W%Nbp
               do ib2=1,W%Nbp
                  !
                  ReW = 0d0 ; ImW = 0d0
                  ReW = cubic_interp(wmats(2:Nmats),real(W%screened_local(ib1,ib2,2:Nmats)),0d0)
                  ImW = cubic_interp(wmats(2:Nmats),aimag(W%screened_local(ib1,ib2,2:Nmats)),0d0)
                  W%screened_local(ib1,ib2,1) = dcmplx(ReW,ImW)
                  !
                  !if(.not.LocalOnly)then
                  !   do iq=1,W%Nkpt
                  !      ReW = 0d0 ; ImW = 0d0
                  !      ReW = cubic_interp(wmats(2:Nmats),real(W%screened(ib1,ib2,2:Nmats,iq)),0d0)
                  !      ImW = cubic_interp(wmats(2:Nmats),aimag(W%screened(ib1,ib2,2:Nmats,iq)),0d0)
                  !      W%screened(ib1,ib2,1,iq) = dcmplx(ReW,ImW)
                  !   enddo
                  !endif
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
      if(.not.G%status) stop "interpolate2Beta_Fermionic. FermionicField not properly initialized."
      if(G%Beta.ne.Beta_Match%Beta_old) stop "interpolate2Beta_Fermionic. Fermionic field Beta is different from the expected one."
      if(G%Npoints.ne.Beta_Match%Nmats_old) stop "interpolate2Beta_Fermionic. Fermionic field Npoints is different from the expected one."
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
            stop "interpolate2Beta_Fermionic. Available Modes are: imp, lat."
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
      if(.not.W%status) stop "interpolate2Beta_Bosonic. BosonicField not properly initialized."
      if(W%Beta.ne.Beta_Match%Beta_old) stop "interpolate2Beta_Bosonic. Bosonic field Beta is different from the expected one."
      if(W%Npoints.ne.Beta_Match%Nmats_old) stop "interpolate2Beta_Bosonic. Bosonic field Npoints is different from the expected one."
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
            stop "interpolate2Beta_Bosonic. Available Modes are: imp, lat."
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

end module post_processing
