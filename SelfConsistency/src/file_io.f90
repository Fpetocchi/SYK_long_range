module file_io

   implicit none
   private

   !===========================================================================!

   ! COMMENTS:
   ! all dirpaths must end with /
   ! Here anything can be written in binfmt or in standard format however
   ! the idea is that the main code as default writes in binary and only optionally
   ! stores data in readable format. That's why several routines to read data
   ! from standard format are missing.
   ! read_BosonicField_Kdep is missing because implemnted ad hoc by now

   !---------------------------------------------------------------------------!
   !PURPOSE: Module interfaces
   !---------------------------------------------------------------------------!
   interface dump_Matrix
      module procedure :: dump_Matrix_local_d                                   ![Umat(Norb,Norb),printpath]............................................( Writes only to formatted output )
      module procedure :: dump_Matrix_local_z                                   ![Umat(Norb,Norb),printpath]............................................( Writes only to formatted output )
      module procedure :: dump_Matrix_Kdep_d                                    ![Umat(Norb,Norb,Nkpt),dirpath,filename,binfmt].........................( Writes output depending on binfmt )
      module procedure :: dump_Matrix_Kdep_z                                    ![Umat(Norb,Norb,Nkpt),dirpath,filename,binfmt].........................( Writes output depending on binfmt )
   end interface dump_Matrix

   interface read_Matrix
      module procedure :: read_Matrix_local_d                                   ![Umat(Norb,Norb),readpath].............................................( Reads only from formatted input )
      module procedure :: read_Matrix_local_z                                   ![Umat(Norb,Norb),readpath].............................................( Reads only from formatted input )
      module procedure :: read_Matrix_Kdep_d                                    ![Umat(Norb,Norb,Nkpt),dirpath,filename]................................( Reads only from unformatted input )
      module procedure :: read_Matrix_Kdep_z                                    ![Umat(Norb,Norb,Nkpt),dirpath,filename]................................( Reads only from unformatted input )
   end interface read_Matrix

   interface dump_Field_component
      module procedure :: dump_Field_component_d                                ![Array(Nfreq),dirpath,filename,axis]...................................( Writes only to formatted output )
      module procedure :: dump_Field_component_z                                ![Array(Nfreq),dirpath,filename,axis]...................................( Writes only to formatted output )
   end interface dump_Field_component

   interface dump_FermionicField
      module procedure :: dump_Field_component_d                                ![Array(Nfreq),dirpath,filename,axis]...................................( Writes only to formatted output )
      module procedure :: dump_Field_component_z                                ![Array(Nfreq),dirpath,filename,axis]...................................( Writes only to formatted output )
      module procedure :: dump_FermionicField_local                             ![FermionicField,dirpath,filename,axis(optional)].......................( Writes only to formatted output )
      module procedure :: dump_FermionicField_Kdep                              ![FermionicField,dirpath,filename,binfmt,kpt(3,Nkpt),axis(optional)]....( Writes output depending on binfmt )
   end interface dump_FermionicField

   interface read_FermionicField
      module procedure :: read_FermionicField_local                             ![FermionicField,dirpath,filename,axis(optional)].......................( Reads only from formatted input )
      module procedure :: read_FermionicField_Kdep                              ![FermionicField,dirpath,filename,kpt(3,Nkpt),axis(optional)]...........( Reads only from unformatted input )
   end interface read_FermionicField

   interface dump_BosonicField
      module procedure :: dump_Field_component_d                                ![Array(Nfreq),dirpath,filename,axis]...................................( Writes only to formatted output )
      module procedure :: dump_Field_component_z                                ![Array(Nfreq),dirpath,filename,axis]...................................( Writes only to formatted output )v
      module procedure :: dump_BosonicField_local                               ![BosonicField,dirpath,filename,axis(optional)].........................( Writes only to formatted output )
      module procedure :: dump_BosonicField_Kdep_SPEXlike                       ![BosonicField,dirpath,binfmt,axis(optional)]...........................( Writes output depending on binfmt )
   end interface dump_BosonicField

   interface read_BosonicField
      module procedure :: read_BosonicField_local                               ![BosonicField,dirpath,filename,axis(optional)].........................( Reads only from unformatted input )
   end interface read_BosonicField


   !---------------------------------------------------------------------------!
   !PURPOSE: Module variables
   !---------------------------------------------------------------------------!
   integer,private                          :: Naxis
   real(8),allocatable,private              :: axis_(:)
   !
#ifdef _verb
   logical,private                          :: verbose=.true.
#else
   logical,private                          :: verbose=.false.
#endif

   !---------------------------------------------------------------------------!
   !PURPOSE: Rutines available for the user. Description only for interfaces.
   !---------------------------------------------------------------------------!
   !subroutines
   public :: dump_Matrix
   public :: read_Matrix
   public :: dump_Field_component
   public :: dump_FermionicField
   public :: read_FermionicField
   public :: dump_BosonicField
   public :: read_BosonicField


   !===========================================================================!

contains


   !---------------------------------------------------------------------------!
   !PURPOSE: Write square matrices to file
   !TEST ON: 17-11-2020
   !---------------------------------------------------------------------------!
   subroutine dump_Matrix_local_z(Umat,printpath)
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
      if(verbose)write(*,"(A)") "---- dump_Matrix_local_z"
      write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
      !
      !
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do i=1,size(Umat,dim=1)
         write(unit,"(999E20.12)") (real(Umat(i,j)),j=1,size(Umat,dim=2))
      enddo
      write(unit,*)
      do i=1,size(Umat,dim=1)
         write(unit,"(999E20.12)") (dimag(Umat(i,j)),j=1,size(Umat,dim=2))
      enddo
      close(unit)
      !
   end subroutine dump_Matrix_local_z
   !
   subroutine dump_Matrix_local_d(Umat,printpath)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Umat(:,:)
      character(len=*),intent(in)           :: printpath
      !
      integer                               :: unit
      integer                               :: i,j
      !
      !
      if(verbose)write(*,"(A)") "---- dump_Matrix_local_d"
      write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
      !
      !
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do i=1,size(Umat,dim=1)
         write(unit,"(999E20.12)") (real(Umat(i,j)),j=1,size(Umat,dim=2))
      enddo
      close(unit)
      !
   end subroutine dump_Matrix_local_d
   !
   subroutine dump_Matrix_Kdep_z(Umat,dirpath,filename,binfmt,ispin)
      !
      use utils_misc
      use utils_fields
      implicit none
      !
      complex(8),intent(in)                 :: Umat(:,:,:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      logical,intent(in)                    :: binfmt
      integer,intent(in),optional           :: ispin
      !
      integer                               :: unit,ik,Norb,Nkpt
      integer                               :: iwan1,iwan2
      character(len=256)                    :: printpath
      !
      !
      if(verbose)write(*,"(A)") "---- dump_Matrix_Kdep_z"
      !
      !
      ! Check on the input matrix
      Nkpt = size(Umat,dim=3)
      Norb = size(Umat,dim=1)
      if(Norb.ne.size(Umat,dim=2)) stop "dump_Matrix_Kdep_z: The provided matrix is not square."
      !
      ! Create directory
      call createDir(reg(dirpath),verb=verbose)
      !
      ! Write to file
      if(binfmt) then
         !
         if(present(ispin))then
            printpath = reg(dirpath)//reg(filename)//"_k_s"//str(ispin)//".DAT"
         else
            printpath = reg(dirpath)//reg(filename)//"_k.DAT."
         endif
         write(*,"(A)") "     Dump "//reg(printpath)//" (binary)"
         !
         unit = free_unit()
         open(unit,file=reg(printpath),form="unformatted",status="unknown",position="rewind",action="write")
         !
         write(unit) Nkpt,Norb
         do ik=1,Nkpt
            write(unit) ik
            do iwan1=1,Norb
               do iwan2=1,Norb
                  write(unit) iwan1,iwan2,dreal(Umat(iwan1,iwan2,ik)),dimag(Umat(iwan1,iwan2,ik))
               enddo
            enddo
         enddo !ik
         close(unit)
         !
      else
         !
         do ik=1,Nkpt
            !
            if(present(ispin))then
               printpath = reg(dirpath)//reg(filename)//"_ik"//str(ik)//"_s"//str(ispin)//".DAT"
            else
               printpath = reg(dirpath)//reg(filename)//"_ik"//str(ik)//".DAT"
            endif
            write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
            unit = free_unit()
            open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
            !
            write(unit,"(2(A,1I7))") "ik",ik,"Norb",Norb
            do iwan1=1,Norb
               do iwan2=1,Norb
                  write(unit,"(2I4,2E20.12)") iwan1,iwan2,dreal(Umat(iwan1,iwan2,ik)),dimag(Umat(iwan1,iwan2,ik))
               enddo
            enddo
            !
            close(unit)
            !
         enddo !ik
         !
      endif
      !
   end subroutine dump_Matrix_Kdep_z
   !
   subroutine dump_Matrix_Kdep_d(Umat,dirpath,filename,binfmt,ispin)
      !
      use utils_misc
      use utils_fields
      implicit none
      !
      real(8),intent(in)                    :: Umat(:,:,:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      logical,intent(in)                    :: binfmt
      integer,intent(in),optional           :: ispin
      !
      integer                               :: unit,ik,Norb,Nkpt
      integer                               :: iwan1,iwan2
      character(len=256)                    :: printpath
      !
      !
      if(verbose)write(*,"(A)") "---- dump_Matrix_Kdep_d"
      !
      !
      ! Check on the input matrix
      Nkpt = size(Umat,dim=3)
      Norb = size(Umat,dim=1)
      if(Norb.ne.size(Umat,dim=2)) stop "dump_Matrix_Kdep_d: The provided matrix is not square."
      !
      ! Create directory
      call createDir(reg(dirpath),verb=verbose)
      !
      ! Write to file
      if(binfmt) then
         !
         if(present(ispin))then
            printpath = reg(dirpath)//reg(filename)//"_k_s"//str(ispin)//".DAT"
         else
            printpath = reg(dirpath)//reg(filename)//"_k.DAT."
         endif
         write(*,"(A)") "     Dump "//reg(printpath)//" (binary)"
         !
         unit = free_unit()
         open(unit,file=reg(printpath),form="unformatted",status="unknown",position="rewind",action="write")
         !
         write(unit) Nkpt,Norb
         do ik=1,Nkpt
            write(unit) ik
            do iwan1=1,Norb
               do iwan2=1,Norb
                  write(unit) iwan1,iwan2,Umat(iwan1,iwan2,ik)
               enddo
            enddo
         enddo !ik
         close(unit)
         !
      else
         !
         do ik=1,Nkpt
            !
            if(present(ispin))then
               printpath = reg(dirpath)//reg(filename)//"_ik"//str(ik)//"_s"//str(ispin)//".DAT"
            else
               printpath = reg(dirpath)//reg(filename)//"_ik"//str(ik)//".DAT"
            endif
            write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
            unit = free_unit()
            open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
            !
            write(unit,"(2(A,1I7))") "ik",ik,"Norb",Norb
            do iwan1=1,Norb
               do iwan2=1,Norb
                  write(unit,"(2I4,2E20.12)") iwan1,iwan2,Umat(iwan1,iwan2,ik)
               enddo
            enddo
            !
            close(unit)
            !
         enddo !ik
         !
      endif
      !
   end subroutine dump_Matrix_Kdep_d


   !---------------------------------------------------------------------------!
   !PURPOSE: Read square matrices from file
   !TEST ON: 17-11-2020
   !---------------------------------------------------------------------------!
   subroutine read_Matrix_local_z(Umat,readpath)
      !
      use utils_misc
      implicit none
      !
      complex(8),intent(inout)              :: Umat(:,:)
      character(len=*),intent(in)           :: readpath
      !
      logical                               :: filexists
      real(8),allocatable                   :: RealM(:,:)
      real(8),allocatable                   :: ImagM(:,:)
      integer                               :: unit
      integer                               :: i,j
      !
      !
      if(verbose)write(*,"(A)") "---- read_Matrix_local_z"
      write(*,"(A)") "     Read "//reg(readpath)
      !
      call inquireFile(reg(readpath),filexists,verb=verbose)
      allocate(RealM(size(Umat,dim=1),size(Umat,dim=2)));RealM=0d0
      allocate(ImagM(size(Umat,dim=1),size(Umat,dim=2)));ImagM=0d0
      !
      unit = free_unit()
      open(unit,file=reg(readpath),form="formatted",status="old",position="rewind",action="read")
      do i=1,size(Umat,dim=1)
         read(unit,"(999E20.12)") (RealM(i,j),j=1,size(Umat,dim=2))
      enddo
      read(unit,*)
      do i=1,size(Umat,dim=1)
         read(unit,"(999E20.12)") (ImagM(i,j),j=1,size(Umat,dim=2))
      enddo
      close(unit)
      !
      Umat = dcmplx(0d0,0d0)
      Umat = RealM + dcmplx(0d0,1d0)*ImagM
      !
   end subroutine read_Matrix_local_z
   !
   subroutine read_Matrix_local_d(Umat,readpath)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(inout)                 :: Umat(:,:)
      character(len=*),intent(in)           :: readpath
      !
      logical                               :: filexists
      real(8),allocatable                   :: RealM(:,:)
      integer                               :: unit
      integer                               :: i,j
      !
      !
      if(verbose)write(*,"(A)") "---- read_Matrix_local_d"
      write(*,"(A)") "     Read "//reg(readpath)
      !
      call inquireFile(reg(readpath),filexists,verb=verbose)
      allocate(RealM(size(Umat,dim=1),size(Umat,dim=2)));RealM=0d0
      !
      unit = free_unit()
      open(unit,file=reg(readpath),form="formatted",status="old",position="rewind",action="read")
      do i=1,size(Umat,dim=1)
         read(unit,"(999E20.12)") (RealM(i,j),j=1,size(Umat,dim=2))
      enddo
      close(unit)
      !
      Umat = 0d0
      Umat = RealM
      !
   end subroutine read_Matrix_local_d
   !
   subroutine read_Matrix_Kdep_z(Umat,dirpath,filename)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      complex(8),intent(inout)              :: Umat(:,:,:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      !
      integer                               :: unit
      integer                               :: Norb,ik,Nkpt
      integer                               :: iwan1,iwan2
      integer                               :: idum1,idum2
      integer                               :: Nkpt_read,Norb_read
      real(8)                               :: RealM,ImagM
      logical                               :: filexists
      character(len=256)                    :: readpath
      !
      !
      if(verbose)write(*,"(A)") "---- read_Matrix_Kdep_z"
      readpath = reg(dirpath)//filename
      write(*,"(A)") "     Read "//reg(readpath)
      !
      ! Check on the input Matrix
      Nkpt = size(Umat,dim=3)
      Norb = size(Umat,dim=1)
      if(Norb.ne.size(Umat,dim=2)) stop "read_Matrix_Kdep_z:The provided matrix is not square."
      !
      ! Check file existence
      call inquireFile(reg(readpath),filexists,verb=verbose)
      !
      !
      unit = free_unit()
      open(unit,file=reg(readpath),form="unformatted",status="old",position="rewind",action="read")
      !
      read(unit) Nkpt_read,Norb_read
      !
      if(Nkpt_read.ne.Nkpt) stop "read_Matrix_Kdep_z: File with wrong number of K-points."
      if(Norb_read.ne.Norb) stop "read_Matrix_Kdep_z: File with wrong number of Wannier functions."
      !
      do ik=1,Nkpt
         !
         read(unit) idum1
         if (idum1.ne.ik) stop "read_Matrix_Kdep_z: ik does not match."
         !
         do iwan1=1,Norb
            do iwan2=1,Norb
               !
               read(unit) idum1,idum2,RealM,ImagM
               if (idum1.ne.iwan1) stop "read_Matrix_Kdep_z: iwan1 does not match."
               if (idum2.ne.iwan2) stop "read_Matrix_Kdep_z: iwan2 does not match."
               Umat(iwan1,iwan2,ik) = dcmplx(RealM,ImagM)
               !
            enddo
         enddo
         !
      enddo !ik
      close(unit)
      !
   end subroutine read_Matrix_Kdep_z
   !
   subroutine read_Matrix_Kdep_d(Umat,dirpath,filename)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      real(8),intent(inout)                 :: Umat(:,:,:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      !
      integer                               :: unit
      integer                               :: Norb,ik,Nkpt
      integer                               :: iwan1,iwan2
      integer                               :: idum1,idum2
      integer                               :: Nkpt_read,Norb_read
      real(8)                               :: RealM
      logical                               :: filexists
      character(len=256)                    :: readpath
      !
      !
      if(verbose)write(*,"(A)") "---- read_Matrix_Kdep_d"
      readpath = reg(dirpath)//filename
      write(*,"(A)") "     Read "//reg(readpath)
      !
      ! Check on the input Matrix
      Nkpt = size(Umat,dim=3)
      Norb = size(Umat,dim=1)
      if(Norb.ne.size(Umat,dim=2)) stop "read_Matrix_Kdep_d: The provided matrix is not square."
      !
      ! Check file existence
      call inquireFile(reg(readpath),filexists,verb=verbose)
      !
      !
      unit = free_unit()
      open(unit,file=reg(readpath),form="unformatted",status="old",position="rewind",action="read")
      !
      read(unit) Nkpt_read,Norb_read
      !
      if(Nkpt_read.ne.Nkpt) stop "read_Matrix_Kdep_d: File with wrong number of K-points."
      if(Norb_read.ne.Norb) stop "read_Matrix_Kdep_d: File with wrong number of Wannier functions."
      !
      do ik=1,Nkpt
         !
         read(unit) idum1
         if (idum2.ne.ik) stop "read_Matrix_Kdep_d: ik does not match."
         !
         do iwan1=1,Norb
            do iwan2=1,Norb
               !
               read(unit) idum1,idum2,RealM
               if (idum1.ne.iwan1) stop "read_Matrix_Kdep_d: iwan1 does not match."
               if (idum2.ne.iwan2) stop "read_Matrix_Kdep_d: iwan2 does not match."
               Umat(iwan1,iwan2,ik) = RealM
               !
            enddo
         enddo
         !
      enddo !ik
      close(unit)
      !
   end subroutine read_Matrix_Kdep_d


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file a single Fermionic/Bosonic component
   !TEST ON: 14-10-2020
   !---------------------------------------------------------------------------!
   subroutine dump_Field_component_d(Fcomp,dirpath,filename,axis)
      !
      use utils_misc
      implicit none
      !
      real(8),intent(in)                    :: Fcomp(:)
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),allocatable,intent(in)        :: axis(:)
      !
      integer                               :: unit,iaxis
      character(len=256)                    :: printpath
      !
      !
      if(verbose)write(*,"(A)") "---- dump_Field_component_d"
      printpath = reg(dirpath)//filename
      write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
      !
      !
      ! Create directory
      call createDir(reg(dirpath),verb=verbose)
      !
      ! Write to file
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do iaxis=1,size(axis)
         write(unit,"(2E20.12)") axis(iaxis), Fcomp(iaxis)
      enddo
      close(unit)
      !
   end subroutine dump_Field_component_d
   !
   subroutine dump_Field_component_z(Fcomp,dirpath,filename,axis)
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
      if(verbose)write(*,"(A)") "---- dump_Field_component_z"
      printpath = reg(dirpath)//filename
      write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
      !
      !
      ! Create directory
      call createDir(reg(dirpath),verb=verbose)
      !
      ! Write to file
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      do iaxis=1,size(axis)
         write(unit,"(3E20.12)") axis(iaxis), real(Fcomp(iaxis)), dimag(Fcomp(iaxis))
      enddo
      close(unit)
      !
   end subroutine dump_Field_component_z


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file Fermionic fields
   !TEST ON: 16-10-2020
   !---------------------------------------------------------------------------!
   subroutine dump_FermionicField_local(G,dirpath,filename,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(in)       :: G
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),intent(in),optional           :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb
      integer                               :: iwan1,iwan2,ispin
      character(len=256)                    :: printpath
      !
      !
      if(verbose)write(*,"(A)") "---- dump_FermionicField_local"
      !
      !
      ! Check on the input Field
      if(.not.G%status) stop "dump_FermionicField_local: Field not properly initialized."
      if(present(axis))then
         if(size(axis).ne.G%Npoints)write(*,"(A)")"     Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Writing up to the smaller."
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
      call createDir(reg(dirpath),verb=verbose)
      !
      ! Write to file
      do ispin=1,Nspin
         !
         printpath = reg(dirpath)//filename//"_s"//str(ispin)//".DAT"
         write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
         !
         unit = free_unit()
         open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
         write(unit,"(I5,A)") Norb," Number of Wannier functions"
         write(unit,"(I5,A)") G%Npoints," Number of grid points"
         write(unit,"(1E20.12,A)") G%mu," chemical potential"
         write(unit,"(A)") "Wannier-projected fermionic components:"
         do iaxis=1,G%Npoints
            do iwan1=1,Norb
               do iwan2=1,Norb
                  write(unit,"(1E20.12,2I4,2E20.12)") axis_(iaxis),iwan1,iwan2,dreal(G%ws(iwan1,iwan2,iaxis,ispin)),dimag(G%ws(iwan1,iwan2,iaxis,ispin))
               enddo
            enddo
         enddo
         close(unit)
         !
      enddo
      !
   end subroutine dump_FermionicField_local
   !
   subroutine dump_FermionicField_Kdep(G,dirpath,filename,binfmt,kpt,axis)
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
      real(8),intent(in)                    :: kpt(:,:)
      real(8),intent(in),optional           :: axis(:)
      !
      integer                               :: unit,ik,iaxis,Norb
      integer                               :: ispin,iwan1,iwan2
      character(len=256)                    :: printpath
      !
      !
      if(verbose)write(*,"(A)") "---- dump_FermionicField_Kdep"
      !
      !
      ! Check on the input Field
      if(.not.G%status) stop "dump_FermionicField_Kdep: Field not properly initialized."
      if(G%Nkpt.eq.0) stop "dump_FermionicField_Kdep: K-dependent part not allocated."
      call assert_shape(kpt,[3,G%Nkpt],"dump_FermionicField_Kdep","kpt")
      if(present(axis))then
         if(size(axis).ne.G%Npoints)write(*,"(A)")"     Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Writing up to the smaller."
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
      call createDir(reg(dirpath),verb=verbose)
      !
      ! Write to file
      do ispin=1,Nspin
         !
         if(binfmt) then
            !
            printpath = reg(dirpath)//reg(filename)//"_k_s"//str(ispin)//".DAT"
            write(*,"(A)") "     Dump "//reg(printpath)//" (binary)"
            !
            unit = free_unit()
            open(unit,file=reg(printpath),form="unformatted",status="unknown",position="rewind",action="write")
            !
            write(unit) ispin,G%Nkpt,Norb,Naxis,G%mu
            do ik=1,G%Nkpt
               write(unit) ispin,ik,kpt(:,ik)
               do iaxis=1,Naxis
                  write(unit) iaxis,axis_(iaxis)
                  do iwan1=1,Norb
                     do iwan2=1,Norb
                        write(unit) iwan1,iwan2,dreal(G%wks(iwan1,iwan2,iaxis,ik,ispin)),dimag(G%wks(iwan1,iwan2,iaxis,ik,ispin))
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
               printpath = reg(dirpath)//reg(filename)//"_ik"//str(ik)//"_s"//str(ispin)//".DAT"
               write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
               !
               unit = free_unit()
               open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
               !
               write(unit,"(1I3,1I7,3E20.12)") ispin,ik,kpt(:,ik)
               do iaxis=1,Naxis
                  write(unit,"(1I7,1E20.12)") iaxis,axis_(iaxis)
                  do iwan1=1,Norb
                     do iwan2=1,Norb
                        write(unit,"(2I4,2E20.12)") iwan1,iwan2,dreal(G%wks(iwan1,iwan2,iaxis,ik,ispin)),dimag(G%wks(iwan1,iwan2,iaxis,ik,ispin))
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
   !PURPOSE: Read from file Fermionic fields
   !TEST ON: 16-10-2020
   !---------------------------------------------------------------------------!
   subroutine read_FermionicField_local(G,dirpath,filename,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),intent(inout),optional        :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb,ispin
      integer                               :: iwan1,iwan2
      integer                               :: idum1,idum2
      integer                               :: Norb_read,Naxis_read
      real(8)                               :: mu_read(2)
      real(8)                               :: axispoint,RealG,ImagG
      logical                               :: filexists
      character(len=256)                    :: readpath
      !
      !
      if(verbose)write(*,"(A)") "---- read_FermionicField_local"
      !
      !
      ! Check on the input Field
      if(.not.G%status) stop "read_FermionicField_local: Field not properly initialized."
      Norb = G%Norb
      !
      do ispin=1,Nspin
         !
         readpath = reg(dirpath)//filename//"_s"//str(ispin)//".DAT"
         write(*,"(A)") "     Read "//reg(readpath)
         call inquireFile(reg(readpath),filexists,verb=verbose)
         !
         ! Read file
         unit = free_unit()
         open(unit,file=reg(readpath),form="formatted",status="old",position="rewind",action="read")
         !
         read(unit,*) Norb_read !," Number of Wannier functions"
         read(unit,*) Naxis_read !," Number of grid points"
         read(unit,*) mu_read(ispin) !," chemical potential"
         !
         if(Norb_read.ne.Norb) stop "read_FermionicField_local: File with wrong number of Wannier functions."
         if(present(axis))then
            if(size(axis).ne.G%Npoints)write(*,"(A)")"     Warning: Axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
            Naxis = min(size(axis),G%Npoints)
            if(Naxis.ne.Naxis_read)write(*,"(A)")"     Warning: Expected grid: "//str(Naxis)//" does not match with files grid: "//str(Naxis_read)//". Reading up to the smaller."
            Naxis = min(Naxis,Naxis_read)
            if(allocated(axis_))deallocate(axis_)
            allocate(axis_(Naxis));axis_=0d0
         else
            if(Naxis_read.ne.G%Npoints)write(*,"(A)")"     Warning: Files grid: "//str(Naxis_read)//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
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
                  if (idum1.ne.iwan1) stop "read_FermionicField_local: iwan1 does not match."
                  if (idum2.ne.iwan2) stop "read_FermionicField_local: iwan2 does not match."
                  G%ws(iwan1,iwan2,iaxis,ispin) = dcmplx(RealG,ImagG)
                  !
               enddo
            enddo
            if((.not.present(axis)).and.(abs(axispoint-axis_(iaxis)).gt.eps)) stop "read_FermionicField_local: axispoint does not match with expected grid."
            axis_(iaxis) = axispoint
         enddo
         close(unit)
         !
      enddo
      !
      if((Nspin.eq.2).and.(mu_read(1).ne.mu_read(2))) stop "read_FermionicField_local: Chemical potential is different between up and down."
      G%mu=mu_read(1)
      if(present(axis)) axis(1:Naxis)=axis_
      !
   end subroutine read_FermionicField_local
   !
   subroutine read_FermionicField_Kdep(G,dirpath,filename,kpt,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(FermionicField),intent(inout)    :: G
      character(len=*),intent(in)           :: dirpath
      character(len=*),intent(in)           :: filename
      real(8),intent(in)                    :: kpt(:,:)
      real(8),intent(inout),optional        :: axis(:)
      !
      integer                               :: unit
      integer                               :: iaxis,Norb,ispin,ik
      integer                               :: iwan1,iwan2
      integer                               :: idum1,idum2
      integer                               :: ispin_read,Nkpt_read,Norb_read,Naxis_read
      real(8)                               :: mu_read(2)
      real(8)                               :: axispoint,RealG,ImagG
      real(8)                               :: kvec(3)
      logical                               :: filexists
      character(len=256)                    :: readpath
      !
      !
      if(verbose)write(*,"(A)") "---- read_FermionicField_Kdep"
      !
      !
      ! Check on the input Field
      if(.not.G%status) stop "read_FermionicField_Kdep: Field not properly initialized."
      if(G%Nkpt.eq.0) stop "read_FermionicField_Kdep: K-dependent part not allocated."
      Norb = G%Norb
      call assert_shape(kpt,[3,G%Nkpt],"read_FermionicField_Kdep","kpt")
      !
      !
      ! Read file
      do ispin=1,Nspin
         !
         readpath = reg(dirpath)//reg(filename)//"_k_s"//str(ispin)//".DAT"
         write(*,"(A)") "     Read "//reg(readpath)
         call inquireFile(reg(readpath),filexists,verb=verbose)
         !
         unit = free_unit()
         open(unit,file=reg(readpath),form="unformatted",status="old",position="rewind",action="read")
         !
         read(unit) ispin_read,Nkpt_read,Norb_read,Naxis_read,mu_read(ispin)
         !
         if(ispin_read.ne.ispin) stop "read_FermionicField_Kdep: File with wrong spin index."
         if(Nkpt_read.ne.G%Nkpt) stop "read_FermionicField_Kdep: File with wrong number of K-points."
         if(Norb_read.ne.Norb) stop "read_FermionicField_Kdep: File with wrong number of Wannier functions."
         if(present(axis))then
            if(size(axis).ne.G%Npoints)write(*,"(A)")"     Warning: Axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
            Naxis = min(size(axis),G%Npoints)
            if(Naxis.ne.Naxis_read)write(*,"(A)")"     Warning: Expected grid: "//str(Naxis)//" does not match with files grid: "//str(Naxis_read)//". Reading up to the smaller."
            Naxis = min(Naxis,Naxis_read)
            if(allocated(axis_))deallocate(axis_)
            allocate(axis_(Naxis));axis_=0d0
         else
            if(Naxis_read.ne.G%Npoints)write(*,"(A)")"     Warning: Files grid: "//str(Naxis_read)//" does not match with field mesh: "//str(G%Npoints)//". Reading up to the smaller."
            Naxis = min(Naxis_read,G%Npoints)
            if(allocated(axis_))deallocate(axis_)
            allocate(axis_(Naxis));axis_=0d0
            axis_ = FermionicFreqMesh(G%Beta,G%Npoints)
         endif
         !
         do ik=1,G%Nkpt
            !
            read(unit) idum1,idum2,kvec
            if (idum1.ne.ispin) stop "read_FermionicField_Kdep: ispin does not match."
            if (idum2.ne.ik) stop "read_FermionicField_Kdep: ik does not match."
            if(abs(kvec(1)-kpt(1,ik)).gt.eps)stop "read_FermionicField_Kdep: kvec(1) does not match."
            if(abs(kvec(2)-kpt(2,ik)).gt.eps)stop "read_FermionicField_Kdep: kvec(1) does not match."
            if(abs(kvec(3)-kpt(3,ik)).gt.eps)stop "read_FermionicField_Kdep: kvec(1) does not match."
            !
            do iaxis=1,Naxis
               !
               read(unit) idum1,axispoint
               if (idum1.ne.iaxis) stop "read_FermionicField_Kdep: iaxis does not match"
               if((.not.present(axis)).and.(abs(axispoint-axis_(iaxis)).gt.eps)) stop "read_FermionicField_Kdep: axispoint does not match with expected grid."
               axis_(iaxis) = axispoint
               !
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     !
                     read(unit) idum1,idum2,RealG,ImagG
                     if (idum1.ne.iwan1) stop "read_FermionicField_Kdep: iwan1 does not match."
                     if (idum2.ne.iwan2) stop "read_FermionicField_Kdep: iwan2 does not match."
                     G%wks(iwan1,iwan2,iaxis,ik,ispin) = dcmplx(RealG,ImagG)
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
      call FermionicKsum(G)
      if((Nspin.eq.2).and.(mu_read(1).ne.mu_read(2))) stop "read_FermionicField_Kdep: Chemical potential is different between up and down."
      G%mu=mu_read(1)
      !
   end subroutine read_FermionicField_Kdep


   !---------------------------------------------------------------------------!
   !PURPOSE: Write to file Bosonic fields
   !TEST ON: 20-10-2020(with and without axis)
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
      if(verbose)write(*,"(A)") "---- dump_BosonicField_local"
      printpath = reg(dirpath)//filename
      write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
      !
      !
      ! Check on the input Field
      if(.not.U%status) stop "dump_BosonicField_local: Field not properly initialized."
      if(present(axis))then
         if(size(axis).ne.U%Npoints)write(*,"(A)")"     Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(U%Npoints)//". Writing up to the smaller."
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
      call createDir(reg(dirpath),verb=verbose)
      !
      ! Write to file
      unit = free_unit()
      open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
      write(unit,"(I5,A)") Norb," Number of Wannier functions"
      write(unit,"(I5,A)") U%Npoints," Number of grid points"
      write(unit,"(A)") "Wannier-projected bare limit:"
      do iwan1=1,Norb
         do iwan2=1,Norb
            do iwan3=1,Norb
               do iwan4=1,Norb
                  !
                  ib1 = iwan1+Norb*(iwan2-1)
                  ib2 = iwan3+Norb*(iwan4-1)
                  !
                  if(allocated(U%bare_local))then
                     write(unit,"(4I4,2E20.12)") iwan1,iwan2,iwan3,iwan4,dreal(U%bare_local(ib1,ib2)),dimag(U%bare_local(ib1,ib2))
                  else
                     write(unit,"(4I4,2E20.12)") iwan1,iwan2,iwan3,iwan4,0d0,0d0
                  endif
                  !
               enddo
            enddo
         enddo
      enddo
      write(unit,"(/A)") "Wannier-projected screening dependence:"
      do iaxis=1,Naxis
         do iwan1=1,Norb
            do iwan2=1,Norb
               do iwan3=1,Norb
                  do iwan4=1,Norb
                     !
                     ib1 = iwan1+Norb*(iwan2-1)
                     ib2 = iwan3+Norb*(iwan4-1)
                     !
                     write(unit,"(1E20.12,4I4,2E20.12)") axis_(iaxis),iwan1,iwan2,iwan3,iwan4,dreal(U%screened_local(ib1,ib2,iaxis)),dimag(U%screened_local(ib1,ib2,iaxis))
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
      !
   end subroutine dump_BosonicField_local
   !
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
      integer                               :: unit,iq
      integer                               :: iaxis,Norb
      integer                               :: iwan1,iwan2,iwan3,iwan4
      integer                               :: ib1,ib2
      complex(8),allocatable                :: Utmp(:,:)
      character(len=256)                    :: printpath
      !
      !
      if(verbose)write(*,"(A)") "---- dump_BosonicField_Kdep_SPEXlike"
      !
      !
      ! Check on the input Field
      if(.not.U%status) stop "dump_BosonicField_Kdep_SPEXlike: Field not properly initialized."
      if(U%Nkpt.eq.0) stop "dump_BosonicField_Kdep_SPEXlike: K-dependent part not allocated."
      if(present(axis))then
         if(size(axis).ne.U%Npoints)write(*,"(A)")"     Warning: axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(U%Npoints)//". Writing up to the smaller."
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
      call createDir(reg(dirpath),verb=verbose)
      !
      !Just in case I'm printing a polarization
      if(.not.allocated(U%bare))then
         allocate(Utmp(U%Nbp,U%Nbp))
         Utmp=czero
      endif
      !
      ! Write to file
      do iq=1,U%Nkpt
         !
         if(binfmt) then
            !
            printpath = reg(dirpath)//"VW.Q"//str(iq,4)//".DAT"
            write(*,"(A)") "     Dump "//reg(printpath)//" (binary)"
            !
            unit = free_unit()
            open(unit,file=reg(printpath),form="unformatted",status="unknown",position="rewind",action="write")
            write(unit) iq,1,Norb,Naxis,.true.
            write(unit) axis_
            if(allocated(U%bare))then
               write(unit) U%bare(:,:,iq)*U%Nkpt*U%Nkpt/H2eV
            else
               write(unit) Utmp
            endif
            do iaxis=1,Naxis
               write(unit) U%screened(:,:,iaxis,iq)*U%Nkpt*U%Nkpt/H2eV
            enddo
            close(unit)
            !
         else
            !
            if((iq.gt.1).and.mod(iq,251).ne.0)cycle
            !
            printpath = reg(dirpath)//"VW.Q"//str(iq,4)//".DAT"
            write(*,"(A)") "     Dump "//reg(printpath)//" (readable)"
            !
            unit = free_unit()
            open(unit,file=reg(printpath),form="formatted",status="unknown",position="rewind",action="write")
            write(unit,"(I5,A)") Norb," Number of Wannier functions"
            write(unit,"(I5,A)") U%Npoints," Number of grid points"
            write(unit,"(A)") "Wannier-projected bare limit:"
            do iwan1=1,Norb
               do iwan2=1,Norb
                  do iwan3=1,Norb
                     do iwan4=1,Norb
                        !
                        ib1 = iwan1+Norb*(iwan2-1)
                        ib2 = iwan3+Norb*(iwan4-1)
                        !
                        if(allocated(U%bare_local))then
                           write(unit,"(4I4,2E20.12)") iwan1,iwan2,iwan3,iwan4,dreal(U%bare(ib1,ib2,iq)),dimag(U%bare(ib1,ib2,iq))
                        else
                           write(unit,"(4I4,2E20.12)") iwan1,iwan2,iwan3,iwan4,0d0,0d0
                        endif
                        !
                     enddo
                  enddo
               enddo
            enddo
            write(unit,"(/A)") "Wannier-projected screening dependence:"
            do iaxis=1,Naxis
               do iwan1=1,Norb
                  do iwan2=1,Norb
                     do iwan3=1,Norb
                        do iwan4=1,Norb
                           !
                           ib1 = iwan1+Norb*(iwan2-1)
                           ib2 = iwan3+Norb*(iwan4-1)
                           !
                           write(unit,"(1E20.12,4I4,2E20.12)") axis_(iaxis),iwan1,iwan2,iwan3,iwan4,dreal(U%screened(ib1,ib2,iaxis,iq)),dimag(U%screened(ib1,ib2,iaxis,iq))
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            close(unit)
            !
         endif
         !
      enddo
      if(allocated(Utmp))deallocate(Utmp)
      !
   end subroutine dump_BosonicField_Kdep_SPEXlike


   !---------------------------------------------------------------------------!
   !PURPOSE: Read from file the local attributes of a Bosonic field
   !TEST ON: 23-10-2020
   !---------------------------------------------------------------------------!
   subroutine read_BosonicField_local(U,dirpath,filename,axis)
      !
      use parameters
      use utils_misc
      use utils_fields
      implicit none
      !
      type(BosonicField),intent(inout)      :: U
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
      if(verbose)write(*,"(A)") "---- read_BosonicField_local"
      readpath = reg(dirpath)//filename
      write(*,"(A)") "     Read "//reg(readpath)
      !
      !
      ! Check on the input Field
      if(.not.U%status) stop "read_BosonicField_local: Field not properly initialized."
      Norb = int(sqrt(dble(U%Nbp)))
      !
      ! Check file existence
      call inquireFile(reg(readpath),filexists,verb=verbose)
      !
      ! Read file
      unit = free_unit()
      open(unit,file=reg(readpath),form="formatted",status="old",position="rewind",action="read")
      read(unit,*) Norb_read !," Number of Wannier functions"
      read(unit,*) Naxis_read !," Number of grid points"
      !
      if(Norb_read.ne.Norb) stop "read_BosonicField_local: File with wrong number of Wannier functions."
      if(present(axis))then
         if(size(axis).ne.U%Npoints)write(*,"(A)")"     Warning: Axis provided but its length: "//str(size(axis))//" does not match with field mesh: "//str(U%Npoints)//". Reading up to the smaller."
         Naxis = min(size(axis),U%Npoints)
         if(Naxis.ne.Naxis_read)write(*,"(A)")"     Warning: Expected grid: "//str(Naxis)//" does not match with files grid: "//str(Naxis_read)//". Reading up to the smaller."
         Naxis = min(Naxis,Naxis_read)
         if(allocated(axis_))deallocate(axis_)
         allocate(axis_(Naxis));axis_=0d0
      else
         if(Naxis_read.ne.U%Npoints)write(*,"(A)")"     Warning: Files grid: "//str(Naxis_read)//" does not match with field mesh: "//str(U%Npoints)//". Reading up to the smaller."
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
                  if (idum1.ne.iwan1) stop "read_BosonicField_local: iwan1 (bare) does not match."
                  if (idum2.ne.iwan2) stop "read_BosonicField_local: iwan2 (bare) does not match."
                  if (idum3.ne.iwan3) stop "read_BosonicField_local: iwan3 (bare) does not match."
                  if (idum4.ne.iwan4) stop "read_BosonicField_local: iwan4 (bare) does not match."
                  if(allocated(U%bare_local))U%bare_local(ib1,ib2) = dcmplx(RealU,ImagU)
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
                     if (idum1.ne.iwan1) stop "read_BosonicField_local: iwan1 (screened) does not match."
                     if (idum2.ne.iwan2) stop "read_BosonicField_local: iwan2 (screened) does not match."
                     if (idum3.ne.iwan3) stop "read_BosonicField_local: iwan3 (screened) does not match."
                     if (idum4.ne.iwan4) stop "read_BosonicField_local: iwan4 (screened) does not match."
                     U%screened_local(ib1,ib2,iaxis) = dcmplx(RealU,ImagU)
                     !
                  enddo
               enddo
            enddo
         enddo
         if((.not.present(axis)).and.(abs(axispoint-axis_(iaxis)).gt.eps)) stop "read_BosonicField_local: axispoint does not match with expected grid."
         axis_(iaxis) = axispoint
      enddo
      close(unit)
      !
      if(present(axis)) axis(1:Naxis)=axis_
      !
   end subroutine read_BosonicField_local


end module file_io
