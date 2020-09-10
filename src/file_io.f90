module file_io

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   ! LO SALVO KPT NEL FILE? SE SI COME LO PASSO IL VETTORE?
   ! all dirpaths must end with /

   interface dump_FermionicField
      module procedure :: dump_FermionicField_local                             ![FermionicField,dirpath,filename,ispin,axis(optional for real freq. or tau)]     (Only local projection. Writes only to formatted output.)
      module procedure :: dump_FermionicField_Kdep                              ![FermionicField,dirpath,filename,binfmt,axis(optional for real freq. or tau))]   (Fixed format: [filename]k.dat.[1,2]. Writes to format chosen by binfmt.)
      module procedure :: dump_Field_component                                  ![Array,dirpath,filename,axis]
   end interface dump_FermionicField

   interface read_FermionicField
      module procedure :: read_FermionicField_local                             ![FermionicField,dirpath,filename,ispin,axis(optional for real freq. or tau)]     (Only local projection. Reads only from formatted input.)
      module procedure :: read_FermionicField_Kdep                              ![FermionicField,dirpath,filename,axis(optional for real freq. or tau))]          (Fixed format: [filename]k.dat.[1,2]. Reads only from unformatted input.)
   end interface read_FermionicField

   interface dump_BosonicField
      module procedure :: dump_BosonicField_local                               ![BosonicField,dirpath,filename,axis(optional for real freq. or tau)]             (Only local projection)
      module procedure :: dump_BosonicField_Kdep_SPEXlike                       ![BosonicField,dirpath,binfmt,axis(optional for real freq. or tau))]              (Fixed format: VW_{real,imag}/VW.Q****.DAT )
      module procedure :: dump_Field_component                                  ![Array,dirpath,filename,axis]
   end interface dump_BosonicField

   interface read_BosonicField
      module procedure :: read_BosonicField_local                               ![BosonicField,dirpath,filename,axis(optional for real freq. or tau))]            (Only local projection. Reads only from formatted input.)
   end interface read_BosonicField

   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   integer,private                          :: Naxis
   real(8),allocatable,private              :: axis_(:)

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: dump_Matrix
   public :: dump_FermionicField
   public :: read_FermionicField
   public :: dump_BosonicField
   public :: read_BosonicField

   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file a generic matrix
   !---------------------------------------------------------------------------!
   subroutine dump_matrix(Umat,printpath)
      !
      use utils_misc
      implicit none
      !
      complex(8),intent(in)                 :: Umat(:,:)
      character(len=*),intent(in)           :: printpath
      !
      integer                               :: unit
      integer                               :: i,j
      !
      !
      write(*,*) "--- dump_matrix ---"
      write(*,"(A)") "Open "//reg(printpath)
      !
      !
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do i=1,size(Umat,dim=1)
         write(unit,"(999E20.12)") (real(Umat(i,j)),j=1,size(Umat,dim=2))
      enddo
      write(unit,*)
      do i=1,size(Umat,dim=1)
         write(unit,"(999E20.12)") (aimag(Umat(i,j)),j=1,size(Umat,dim=2))
      enddo
      close(unit)
      !
   end subroutine dump_matrix


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file a single Fermionic/Bosonic component
   !---------------------------------------------------------------------------!
   subroutine dump_Field_component(Fcomp,dirpath,filename,axis)
      !
      use utils_misc
      implicit none
      !
      complex(8),intent(in)                 :: Fcomp(:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),allocatable,intent(in)        :: axis(:)
      !
      integer                               :: unit,iaxis
      character(len=256)                    :: printpath
      !
      !
      write(*,*) "--- dump_Field_component ---"
      printpath = dirpath//filename
      write(*,"(A)") "Open "//reg(printpath)
      !
      !
      ! Create directory
      call createDir(reg(dirpath))
      !
      ! Write to file
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do iaxis=1,size(axis)
         write(unit,"(3E20.12)") axis(iaxis), real(Fcomp(iaxis)), aimag(Fcomp(iaxis))
      enddo
      close(unit)
      !
   end subroutine dump_Field_component


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file the local attributes of a Fermionic field
   !---------------------------------------------------------------------------!
   subroutine dump_FermionicField_local(G,ispin,dirpath,filename,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(in)       :: G
      integer,intent(in)                    :: ispin
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),intent(in),optional           :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb
      integer                               :: iwan1,iwan2
      character(len=256)                    :: printpath
      !
      !
      write(*,*) "--- dump_FermionicField_local ---"
      printpath = dirpath//filename
      write(*,"(A)") "Open "//reg(printpath)
      !
      !
      ! Check on the input Boson
      if(.not.G%status) stop "Field not properly initialized."
      if(present(axis))then
         if(size(axis).ne.G%Npoints)write(*,"(A)")"Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Writing up to the smaller."
         Naxis = min(size(axis),G%Npoints)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = axis(1:Naxis)
      else
         Naxis = G%Npoints
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = FermionicFreqMesh(G%Beta,G%Npoints)
      endif
      Norb = G%Norb
      !
      ! Create directory
      call createDir(reg(dirpath))
      !
      ! Write to file
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(I5,A)") Norb," Number of Wannier functions"
      write(unit,"(I5,A)") G%Npoints," Number of grid points"
      write(unit,"(/A)") "Wannier-projected fermionic components:"
      do iaxis=1,G%Npoints
         do iwan1=1,Norb
            do iwan2=1,Norb
               write(unit,"(1E20.12,2I4,2E20.12)") axis_(iaxis),iwan1,iwan2,real(G%w(iwan1,iwan2,iaxis,ispin)),aimag(G%w(iwan1,iwan2,iaxis,ispin))
            enddo
         enddo
      enddo
      close(unit)
      !
   end subroutine dump_FermionicField_local


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file the K-dependent Fermionic field
   !---------------------------------------------------------------------------!
   subroutine dump_FermionicField_Kdep(G,dirpath,filename,binfmt,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(in)       :: G
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      logical,intent(in)                    :: binfmt
      real(8),intent(in),optional           :: axis(:)
      !
      integer                               :: unit,ik,iaxis,Norb
      integer                               :: ispin,iwan1,iwan2
      character(len=256)                    :: printpath
      !
      !
      write(*,*) "--- dump_FermionicField_Kdep ---"
      !
      !
      ! Check on the input Boson
      if(.not.G%status) stop "Field not properly initialized."
      if(G%Nkpt.eq.0) stop "K-dependent part not allocated."
      if(present(axis))then
         if(size(axis).ne.G%Npoints)write(*,"(A)")"Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Writing up to the smaller."
         Naxis = min(size(axis),G%Npoints)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = axis(1:Naxis)
      else
         Naxis = G%Npoints
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = FermionicFreqMesh(G%Beta,G%Npoints)
      endif
      Norb = G%Norb
      !
      ! Create directory
      call createDir(reg(dirpath))
      !
      ! Write to file
      do ispin=1,Nspin
         !
         if(binfmt) then
            !
            printpath = reg(dirpath)//reg(filename)//"k.dat."//str(ispin)
            write(*,"(A)") "Open "//reg(printpath)
            !
            unit = free_unit()
            open(unit,file=reg(printpath),form="unformatted",status="unknown",position="rewind",action="write")
            !
            write(unit) ispin,G%Nkpt,Norb,Naxis,G%mu
            do ik=1,G%Nkpt
               write(unit) ispin,ik!,kpt(1:3,ik)
               do iaxis=1,Naxis
                  write(unit) iaxis,axis_(iaxis)
                  do iwan1=1,Norb
                     do iwan2=1,Norb
                        write(unit) iwan1,iwan2,real(G%wk(iwan1,iwan2,iaxis,ik,ispin)),aimag(G%wk(iwan1,iwan2,iaxis,ik,ispin))
                     enddo
                  enddo
               enddo !iaxis
            enddo !ik
            close(unit)
            !
         else
            !
            do ik=1,G%Nkpt
               !
               printpath = reg(dirpath)//"Gk_"//str(ik)//".dat."//str(ispin)
               write(*,"(A)") "Open "//reg(printpath)
               unit = free_unit()
               open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
               !
               write(unit,"(1I3,1I7)") ispin,ik!,kpt(1:3,ik)
               do iaxis=1,Naxis
                  write(unit,"(1I7,1E20.12)") iaxis,axis_(iaxis)
                  do iwan1=1,Norb
                     do iwan2=1,Norb
                        write(unit,"(2I4,2E20.12)") iwan1,iwan2,real(G%wk(iwan1,iwan2,iaxis,ik,ispin)),aimag(G%wk(iwan1,iwan2,iaxis,ik,ispin))
                     enddo
                  enddo
               enddo !iaxis
               !
               close(unit)
               !
            enddo !ik
            !
         endif
         !
      enddo !ispin
      !
   end subroutine dump_FermionicField_Kdep


   !---------------------------------------------------------------------------!
   !PURPOSE: Read from file the local attributes of a Fermionic field
   !---------------------------------------------------------------------------!
   subroutine read_FermionicField_local(G,dirpath,filename,ispin,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(out)      :: G
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      integer,intent(in)                    :: ispin
      real(8),intent(inout),optional        :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb
      integer                               :: iwan1,iwan2
      integer                               :: idum1,idum2
      integer                               :: Norb_read,Naxis_read
      real(8)                               :: axispoint,RealG,ImagG
      logical                               :: filexists
      character(len=256)                    :: readpath
      !
      !
      write(*,*) "--- read_FermionicField_local ---"
      readpath = dirpath//filename
      write(*,"(A)") "Open "//reg(readpath)
      !
      !
      ! Check on the input Boson
      if(.not.G%status) stop "Field not properly initialized."
      Norb = G%Norb
      !
      ! Check file existence
      call inquireFile(reg(readpath),filexists)
      !
      ! Read file
      unit = free_unit()
      open(unit,file=reg(readpath),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) Norb_read !," Number of Wannier functions"
      read(unit,*) Naxis_read !," Number of grid points"
      !
      if(Norb_read.ne.Norb) stop "File with wrong number of Wannier functions."
      if(present(axis))then
         if(size(axis).ne.G%Npoints)write(*,"(A)")"Warning! Axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
         Naxis = min(size(axis),G%Npoints)
         if(Naxis.ne.Naxis_read)write(*,"(A)")"Warning! Expected grid: "//str(Naxis)//" does not match with files grid: "//str(Naxis_read)//". Reading up to the smaller."
         Naxis = min(Naxis,Naxis_read)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
      else
         if(Naxis_read.ne.G%Npoints)write(*,"(A)")"Warning! Files grid: "//str(Naxis_read)//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
         Naxis = min(Naxis_read,G%Npoints)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = FermionicFreqMesh(G%Beta,G%Npoints)
      endif
      !
      read(unit,*) !"Wannier-projected fermionic components:"
      do iaxis=1,Naxis
         do iwan1=1,Norb
            do iwan2=1,Norb
               !
               read(unit,"(1F20.10,2I4,2E20.12)") axispoint,idum1,idum2,RealG,ImagG
               if (idum1.ne.iwan1) stop "iwan1 (screened) does not match"
               if (idum2.ne.iwan2) stop "iwan2 (screened) does not match"
               G%w(iwan1,iwan2,iaxis,ispin) = dcmplx(RealG,ImagG)
               !
            enddo
         enddo
         if((.not.present(axis)).and.(axispoint.ne.axis_(iaxis)))stop "axispoint does not match with expected grid."
         axis_(iaxis) = axispoint
      enddo
      close(unit)
      !
      if(present(axis)) axis(1:Naxis)=axis_
      !
   end subroutine read_FermionicField_local


   !---------------------------------------------------------------------------!
   !PURPOSE: Read from file the local attributes of a Fermionic field
   !---------------------------------------------------------------------------!
   subroutine read_FermionicField_Kdep(G,dirpath,filename,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(out)      :: G
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),intent(inout),optional        :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb,ispin,ik
      integer                               :: iwan1,iwan2
      integer                               :: idum1,idum2
      integer                               :: ispin_read,Nkpt_read,Norb_read
      integer                               :: Naxis_read,mu_read
      real(8)                               :: axispoint,RealG,ImagG
      logical                               :: filexists
      character(len=256)                    :: readpath
      !
      !
      write(*,*) "--- read_FermionicField_Kdep ---"
      !
      !
      ! Check on the input Boson
      if(.not.G%status) stop "Field not properly initialized."
      Norb = G%Norb
      !
      ! Check file existence
      call inquireFile(reg(readpath),filexists)
      !
      ! Read file
      do ispin=1,Nspin
         !
         readpath = reg(dirpath)//reg(filename)//"k.dat."//str(ispin)
         write(*,"(A)") "Open "//reg(readpath)
         !
         unit = free_unit()
         open(unit,file=reg(readpath),form="unformatted",status="old",position="rewind",action="read")
         !
         read(unit) ispin_read,Nkpt_read,Norb_read,Naxis_read,mu_read
         !
         if(ispin_read.ne.ispin) stop "File with wrong spin index."
         if(Nkpt_read.ne.G%Nkpt) stop "File with wrong number of K-points."
         if(Norb_read.ne.Norb) stop "File with wrong number of Wannier functions."
         if(present(axis))then
            if(size(axis).ne.G%Npoints)write(*,"(A)")"Warning! Axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
            Naxis = min(size(axis),G%Npoints)
            if(Naxis.ne.Naxis_read)write(*,"(A)")"Warning! Expected grid: "//str(Naxis)//" does not match with files grid: "//str(Naxis_read)//". Reading up to the smaller."
            Naxis = min(Naxis,Naxis_read)
            if(allocated(axis_))deallocate(axis_)
            allocate(axis_(Naxis));axis_=0d0
         else
            if(Naxis_read.ne.G%Npoints)write(*,"(A)")"Warning! Files grid: "//str(Naxis_read)//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
            Naxis = min(Naxis_read,G%Npoints)
            if(allocated(axis_))deallocate(axis_)
            allocate(axis_(Naxis));axis_=0d0
            axis_ = FermionicFreqMesh(G%Beta,G%Npoints)
         endif
         !
         do ik=1,G%Nkpt
            !
            read(unit) idum1,idum2!,kpt(1:3,ik)
            if (idum1.ne.ispin) stop "ispin does not match"
            if (idum2.ne.ik) stop "ik does not match"
            !
            do iaxis=1,Naxis
               !
               read(unit) idum1,axispoint
               if (idum1.ne.iaxis) stop "iaxis does not match"
               if((.not.present(axis)).and.(axispoint.ne.axis_(iaxis)))stop "axispoint does not match with expected grid."
               axis_(iaxis) = axispoint
               !
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     !
                     read(unit) idum1,idum2,RealG,ImagG
                     if (idum1.ne.iwan1) stop "iwan1 does not match"
                     if (idum2.ne.iwan2) stop "iwan2 does not match"
                     G%wk(iwan1,iwan2,iaxis,ik,ispin) = dcmplx(RealG,ImagG)
                     !
                  enddo
               enddo
               !
            enddo !iaxis
            !
         enddo !ik
         close(unit)
         !
      enddo !ispin
      !
   end subroutine read_FermionicField_Kdep


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file the local attributes of a Bosonic field
   !---------------------------------------------------------------------------!
   subroutine dump_BosonicField_local(U,dirpath,filename,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(BosonicField),intent(in)         :: U
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),intent(in),optional           :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb
      integer                               :: iwan1,iwan2,iwan3,iwan4
      integer                               :: ib1,ib2
      character(len=256)                    :: printpath
      !
      !
      write(*,*) "--- dump_BosonicField_local ---"
      printpath = dirpath//filename
      write(*,"(A)") "Open "//reg(printpath)
      !
      !
      ! Check on the input Boson
      if(.not.U%status) stop "Field not properly initialized."
      if(present(axis))then
         if(size(axis).ne.U%Npoints)write(*,"(A)")"Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(U%Npoints)//". Writing up to the smaller."
         Naxis = min(size(axis),U%Npoints)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = axis(1:Naxis)
      else
         Naxis = U%Npoints
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = BosonicFreqMesh(U%Beta,U%Npoints)
      endif
      Norb = int(sqrt(dble(U%Nbp)))
      !
      ! Create directory
      call createDir(reg(dirpath))
      !
      ! Write to file
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(I5,A)") Norb," Number of Wannier functions"
      write(unit,"(I5,A)") U%Npoints," Number of grid points"
      write(unit,"(/A)") "Wannier-projected bare limit:"
      do iwan1=1,Norb
         do iwan2=1,Norb
            do iwan3=1,Norb
               do iwan4=1,Norb
                  !
                  ib1 = iwan1+Norb*(iwan2-1)
                  ib2 = iwan3+Norb*(iwan4-1)
                  !
                  write(unit,"(4I4,2E20.12)") iwan1,iwan2,iwan3,iwan4,real(U%bare_local(ib1,ib2)),aimag(U%bare_local(ib1,ib2))
                  !
               enddo
            enddo
         enddo
      enddo
      write(unit,"(/A)") "Wannier-projected screening dependence:"
      do iaxis=1,U%Npoints
         do iwan1=1,Norb
            do iwan2=1,Norb
               do iwan3=1,Norb
                  do iwan4=1,Norb
                     !
                     ib1 = iwan1+Norb*(iwan2-1)
                     ib2 = iwan3+Norb*(iwan4-1)
                     !
                     write(unit,"(1E20.12,4I4,2E20.12)") axis_(iaxis),iwan1,iwan2,iwan3,iwan4,real(U%screened_local(ib1,ib2,iaxis)),aimag(U%screened_local(ib1,ib2,iaxis))
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
      !
   end subroutine dump_BosonicField_local


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file the K-dependent Bosonic field - SPEX format & units
   !---------------------------------------------------------------------------!
   subroutine dump_BosonicField_Kdep_SPEXlike(U,dirpath,binfmt,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(BosonicField),intent(in)         :: U
      character(len=*),intent(in)           :: dirpath
      logical,intent(in)                    :: binfmt
      real(8),intent(in),optional           :: axis(:)
      !
      integer                               :: unit,iq,iw,Norb
      character(len=256)                    :: printpath
      !
      !
      write(*,*) "--- dump_BosonicField_Kdep_SPEXlike ---"
      !
      !
      ! Check on the input Boson
      if(.not.U%status) stop "Field not properly initialized."
      if(U%Nkpt.eq.0) stop "K-dependent part not allocated."
      if(present(axis))then
         if(size(axis).ne.U%Npoints)write(*,"(A)")"Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(U%Npoints)//". Writing up to the smaller."
         Naxis = min(size(axis),U%Npoints)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = axis(1:Naxis)
      else
         Naxis = U%Npoints
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = BosonicFreqMesh(U%Beta,U%Npoints)
      endif
      Norb = int(sqrt(dble(U%Nbp)))
      !
      ! Create directory
      call createDir(reg(dirpath))
      !
      ! Write to file
      do iq=1,U%Nkpt
         !
         if(binfmt) then
            !
            printpath = reg(dirpath)//"VW.Q"//str(iq)//".DAT"
            write(*,"(A)") "Open "//reg(printpath)
            !
            unit = free_unit()
            open(unit,file=reg(printpath),form="unformatted",status="unknown",position="rewind",action="write")
            write(unit) iq,1,Norb,Naxis,.true.
            write(unit) axis_
            write(unit) U%bare(:,:,iq)*U%Nkpt*U%Nkpt/H2eV
            do iw=1,Naxis
               write(unit) U%screened(:,:,iw,iq)*U%Nkpt*U%Nkpt/H2eV
            enddo
            close(unit)
            !
         else
            !
            call dump_BosonicField_local(U,reg(dirpath),reg("VW.Q"//str(iq)//".DAT"),axis_)
            !
         endif
         !
      enddo
      !
   end subroutine dump_BosonicField_Kdep_SPEXlike


   !---------------------------------------------------------------------------!
   !PURPOSE: Read from file the local attributes of a Bosonic field
   !---------------------------------------------------------------------------!
   subroutine read_BosonicField_local(U,dirpath,filename,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(BosonicField),intent(out)        :: U
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),intent(inout),optional        :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb
      integer                               :: iwan1,iwan2,iwan3,iwan4
      integer                               :: idum1,idum2,idum3,idum4
      integer                               :: ib1,ib2
      integer                               :: Norb_read,Naxis_read
      real(8)                               :: axispoint,RealU,ImagU
      logical                               :: filexists
      character(len=256)                    :: readpath
      !
      !
      write(*,*) "--- read_BosonicField_local ---"
      readpath = dirpath//filename
      write(*,"(A)") "Open "//reg(readpath)
      !
      !
      ! Check on the input Boson
      if(.not.U%status) stop "Field not properly initialized."
      Norb = int(sqrt(dble(U%Nbp)))
      !
      ! Check file existence
      call inquireFile(reg(readpath),filexists)
      !
      ! Read file
      unit = free_unit()
      open(unit,file=reg(readpath),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) Norb_read !," Number of Wannier functions"
      read(unit,*) Naxis_read !," Number of grid points"
      !
      if(Norb_read.ne.Norb) stop "File with wrong number of Wannier functions."
      if(present(axis))then
         if(size(axis).ne.U%Npoints)write(*,"(A)")"Warning! Axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(U%Npoints)//". Reading up to the smaller."
         Naxis = min(size(axis),U%Npoints)
         if(Naxis.ne.Naxis_read)write(*,"(A)")"Warning! Expected grid: "//str(Naxis)//" does not match with files grid: "//str(Naxis_read)//". Reading up to the smaller."
         Naxis = min(Naxis,Naxis_read)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
      else
         if(Naxis_read.ne.U%Npoints)write(*,"(A)")"Warning! Files grid: "//str(Naxis_read)//" does not match with field mesh: "//str(U%Npoints)//". Reading up to the smaller."
         Naxis = min(Naxis_read,U%Npoints)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
         axis_ = BosonicFreqMesh(U%Beta,U%Npoints)
      endif
      !
      read(unit,*) !"Wannier-projected bare limit:"
      do iwan1=1,Norb
         do iwan2=1,Norb
            do iwan3=1,Norb
               do iwan4=1,Norb
                  !
                  ib1 = iwan1+Norb*(iwan2-1)
                  ib2 = iwan3+Norb*(iwan4-1)
                  !
                  read(unit,"(4I4,2E20.12)") idum1,idum2,idum3,idum4,RealU,ImagU
                  if (idum1.ne.iwan1) stop "iwan1 (bare) does not match"
                  if (idum2.ne.iwan2) stop "iwan2 (bare) does not match"
                  if (idum3.ne.iwan3) stop "iwan3 (bare) does not match"
                  if (idum4.ne.iwan4) stop "iwan4 (bare) does not match"
                  U%bare_local(ib1,ib2) = dcmplx(RealU,ImagU)
                  !
               enddo
            enddo
         enddo
      enddo
      read(unit,"(/A)") !"Wannier-projected screening dependence:"
      do iaxis=1,Naxis
         do iwan1=1,Norb
            do iwan2=1,Norb
               do iwan3=1,Norb
                  do iwan4=1,Norb
                     !
                     ib1 = iwan1+Norb*(iwan2-1)
                     ib2 = iwan3+Norb*(iwan4-1)
                     !
                     read(unit,"(1F20.10,4I4,2E20.12)") axispoint,idum1,idum2,idum3,idum4,RealU,ImagU
                     if (idum1.ne.iwan1) stop "iwan1 (screened) does not match"
                     if (idum2.ne.iwan2) stop "iwan2 (screened) does not match"
                     if (idum3.ne.iwan3) stop "iwan3 (screened) does not match"
                     if (idum4.ne.iwan4) stop "iwan4 (screened) does not match"
                     U%screened_local(ib1,ib2,iaxis) = dcmplx(RealU,ImagU)
                     !
                  enddo
               enddo
            enddo
         enddo
         if((.not.present(axis)).and.(axispoint.ne.axis_(iaxis)))stop "axispoint does not match with expected grid."
         axis_(iaxis) = axispoint
      enddo
      close(unit)
      !
      if(present(axis)) axis(1:Naxis)=axis_
      !
   end subroutine read_BosonicField_local




end module file_io
