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
   integer                               :: Norb,Npoints,Ntau
   real(8),allocatable                   :: tau(:),wmats(:),wreal(:)
   character(len=255)                    :: filepath
   !
   !
   if(verbose)write(*,"(A)") "---- dump_MaxEnt_Gfunct"
   !
   !
   Norb = size(G,dim=1)
   Npoints = size(G,dim=2)
   if(Nspin.ne.size(G,dim=3)) stop "dump_MaxEnt_Gfunct: wrong Nspin dimension."
   !
   filepath = reg(dirpath)//reg(filename)//"/"
   call createDir(reg(filepath),verb=verbose)
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
               call dump_Field_component(real(G(iwan,:,ispin)),reg(filepath),reg(filename)//"_t_o"//str(iwan)//"_s"//str(ispin)//".DAT",tau)
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
               call dump_Field_component(G(iwan,:,ispin),reg(filepath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//".DAT",wmats)
               !
               if(present(WmaxPade).and.(WmaxPade.gt.0))then
                  allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
                  allocate(Gpade(Nreal));Gpade=czero
                  Gpade = pade(G(iwan,:,ispin),"Fermionic",wlimit=WmaxPade)
                  call dump_Field_component(Gpade,reg(filepath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//"_pade.DAT",wreal)
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
               call dump_Field_component(Gft(iwan,:,ispin),reg(filepath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//".DAT",wmats)
               !
               if(present(WmaxPade).and.(WmaxPade.gt.0))then
                  allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
                  allocate(Gpade(Nreal));Gpade=czero
                  Gpade = pade(G(iwan,:,ispin),"Fermionic",wlimit=WmaxPade)
                  call dump_Field_component(Gpade,reg(filepath),reg(filename)//"_w_o"//str(iwan)//"_s"//str(ispin)//"_pade.DAT",wreal)
                  deallocate(Gpade,wreal)
               endif
               !
            enddo
         enddo
         deallocate(Gft)
         !
      case("mats2itau")
         !
         Ntau = int(2d0*pi*Nmats)
         if(Npoints.ne.Nmats) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct number fo Matsubara points differ from input."
         if(present(WmaxPade).and.(WmaxPade.gt.0)) write(*,"(A)")"     Warning: dump_MaxEnt_Gfunct pade from mats2itau not done."
         allocate(tau(Ntau));tau=0d0
         tau = linspace(0d0,Beta,Ntau)
         !
         allocate(Gft(Norb,Ntau,Nspin));Gft=czero
         do ispin=1,Nspin
            call Fmats2itau_vec(Beta,G(:,:,ispin),Gft(:,:,ispin),asympt_corr=.true.,tau_uniform=.true.)
         enddo
         !
         do iwan=1,Norb
            if(present(iorb).and.(iwan.ne.iorb))cycle
            do ispin=1,Nspin
               call dump_Field_component(real(Gft(iwan,:,ispin)),reg(filepath),reg(filename)//"_t_o"//str(iwan)//"_s"//str(ispin)//".DAT",tau)
            enddo
         enddo
         deallocate(Gft,tau)
         !
   end select
   !
end subroutine dump_MaxEnt_Gfunct

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
   integer                               :: Npoints,Ntau
   real(8),allocatable                   :: tau(:),wmats(:),wreal(:)
   character(len=255)                    :: filepath
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
   filepath = reg(dirpath)//reg(filename)//"/"
   call createDir(reg(filepath),verb=verbose)
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
         call dump_Field_component(real(W),reg(filepath),reg(filename)//"_t_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",tau)
         !
      case("mats")
         !
         if(Npoints.ne.Nmats) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct number fo Matsubara points differ from input."
         allocate(wmats(Npoints));wmats=0d0
         wmats = BosonicFreqMesh(Beta,Npoints)
         !
         call dump_Field_component(real(W),reg(filepath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",wmats)
         !
         if(present(WmaxPade).and.(WmaxPade.gt.0))then
            allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
            allocate(Wpade(Nreal));Wpade=czero
            Wpade = pade(W,"Bosonic",wlimit=WmaxPade)
            call dump_Field_component(real(Wpade),reg(filepath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//")_pade.DAT",wreal)
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
         call dump_Field_component(real(Wft),reg(filepath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",wmats)
         !
         if(present(WmaxPade).and.(WmaxPade.gt.0))then
            allocate(wreal(Nreal));wreal=linspace(-wrealMax,+wrealMax,Nreal)
            allocate(Wpade(Nreal));Wpade=czero
            Wpade = pade(Wft,"Bosonic",wlimit=WmaxPade)
            call dump_Field_component(real(Wpade),reg(filepath),reg(filename)//"_w_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//")_pade.DAT",wreal)
            deallocate(Wpade,wreal)
         endif
         !
         deallocate(Wft)
         !
      case("mats2itau")
         !
         Ntau = int(2d0*pi*Nmats)
         if(Npoints.ne.Nmats) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct number fo Matsubara points differ from input."
         if(present(WmaxPade).and.(WmaxPade.gt.0)) write(*,"(A)")"     Warning: dump_MaxEnt_Wfunct pade from mats2itau not done."
         allocate(tau(Ntau));tau=0d0
         tau = linspace(0d0,Beta,Ntau)
         !
         allocate(Wft(Ntau));Wft=czero
         call Bmats2itau(Beta,W,Wft,asympt_corr=.true.,tau_uniform=.true.,Umats_bare=W(Npoints))
         Wft(1) = Wft(1) + W(Npoints)
         Wft(Ntau) = Wft(Ntau) + W(Npoints)
         !
         call dump_Field_component(real(Wft),reg(filepath),reg(filename)//"_t_("//str(ndx(1))//","//str(ndx(2))//")("//str(ndx(3))//","//str(ndx(4))//").DAT",tau)
         deallocate(Wft,tau)
         !
   end select
   !
end subroutine dump_MaxEnt_Wfunct

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
