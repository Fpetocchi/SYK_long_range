
!
! USAGE: Gaux2Sigma flag=F parameters=Nw,-wmin,wmax,eta filename=***.DAT_dos.dat
!

program Integrtor
    !
    use omp_lib
    use utils_misc 
    implicit none
    !
    real(8),parameter       :: pi=3.14159265358979323846d0
    !
    !var parsing
    real(8)                 :: parameters(4)=0d0
    logical                 :: flag
    character(1024)         :: filename="./"
    character(1024)         :: command_line
    !
    logical                 :: exists
    integer                 :: Nthread
    integer                 :: ierr,Nreal_read,Nreal_cut,unit=999
    integer                 :: iw,Nreal,iw1,iw2
    real(8)                 :: dw_read,Norm
    real(8)                 :: dw,Wmin,Wmax,eta
    real(8)                 :: dw1,dw2,dw3,D1,D2,D3
    complex(8)              :: Gw
    real(8),allocatable     :: wreal_read(:),ImGaux(:)
    real(8),allocatable     :: wreal(:)
    complex(8),allocatable  :: Gaux(:)
    complex(8),allocatable  :: Sigma(:)
    !
    !reading data from command line
    call read_command_line
    call parse_command_line
    !
    print *, "-------------------------"
    print *, flag
    print *, "'",reg(filename),"'"
    print *, parameters
    print *, "-------------------------"
    print *, " "
    !
    Nthread = omp_get_max_threads()
    write(*,"(A,1I4)") new_line("A")//"Setting Nthread:",Nthread
    !
    Nreal_cut = 20
    !
    !opening file and read how many frequency are present in
    inquire(file=reg(filename),exist=exists)
    if(.not.exists) stop "file not found"
    open(unit,file=reg(filename),form="formatted",status="unknown",position="rewind",action="read")
    !
    Nreal_read=0
    ierr=0
    do while (ierr.eq.0)
       Nreal_read = Nreal_read + 1
       read(unit,*,iostat=ierr)
    enddo
    close(unit)
    !
    Nreal_read = Nreal_read - 2 - 2*Nreal_cut
    write(*,"(A,1I7)") "number of frequency points in MaxEnt file", Nreal_read
    if(Nreal_read.lt.10) stop "Nreal_read.lt.10)"
    !
    !read the data
    allocate(wreal_read(Nreal_read));wreal_read=0d0
    allocate(ImGaux(Nreal_read));ImGaux=0d0
    open(unit,file=reg(filename),form="formatted",status="unknown",position="rewind",action="read")
    do iw=1,Nreal_cut
        read(unit,*)
    enddo
    do iw=1,Nreal_read
       read(unit,*) wreal_read(iw),ImGaux(iw)
    enddo
    close(unit)
    dw_read = abs(wreal_read(10)-wreal_read(9))
    Norm = sum(ImGaux)*dw_read
    write(*,"(A,2F10.6)") "frequency boundaries in MaxEnt file", wreal_read(1),wreal_read(Nreal_read)
    write(*,"(2(A,F10.6))") "normalization integral", Norm, "  dw= ",dw_read
    ImGaux = ImGaux/Norm
    !
    !data of the new axis
    Nreal = Nreal_read
    if(parameters(1).ne.0d0) Nreal = int(parameters(1)) !parameters(1) => new number of real frequancy points
    Wmin = wreal_read(1)
    if(parameters(2).ne.0d0) Wmin = parameters(2) !parameters(2) => lower frequancy boundary
    Wmax = wreal_read(Nreal_read)
    if(parameters(3).ne.0d0) Wmax = parameters(3) !parameters(3) => upper frequancy boundary
    if(parameters(4).ne.0d0) eta = parameters(4)
    !
    allocate(wreal(Nreal));wreal=0d0
    wreal = linspace(Wmin,Wmax,Nreal,mesh=dw)
    write(*,"(A,2F10.6)") "frequency boundaries in output file",  wreal(1),wreal(Nreal)
    write(*,"(2(A,F10.6))") "eta= ", eta, "  dw= ",dw_read
    !
    !KK integral
    allocate(Gaux(Nreal));Gaux=dcmplx(0d0,0d0)
    !$OMP PARALLEL DEFAULT(PRIVATE),&
    !$OMP SHARED(Gaux,Nreal,Nreal_read,ImGaux,wreal_read,wreal,eta)
    !$OMP DO
    do iw1=1,Nreal
        Gw=dcmplx(0d0,0d0)
        do iw2=1,Nreal_read-2,2
            !
            D1 = ImGaux(iw2)   ; dw1 = (wreal_read(iw2+1)-wreal_read(iw2))/3.d0
            D2 = ImGaux(iw2+1) ; dw2 = (wreal_read(iw2+1)-wreal_read(iw2))*4.d0/3.d0
            D3 = ImGaux(iw2+2) ; dw3 = (wreal_read(iw2+1)-wreal_read(iw2))/3.d0
            !
            !Integrate using Simpson method
            !if(wreal_read(iw2).gt.0.d0) then
            if(dabs(wreal_read(iw2)).ge.1.d-12) then
                !
                Gw = Gw + ( D1/(dcmplx(wreal(iw1),eta)-wreal_read(iw2)  ) ) * dw1
                Gw = Gw + ( D2/(dcmplx(wreal(iw1),eta)-wreal_read(iw2+1)) ) * dw2
                Gw = Gw + ( D3/(dcmplx(wreal(iw1),eta)-wreal_read(iw2+2)) ) * dw3
                !
            else! if(dabs(wreal_read(iw2)).lt.1.d-9) then
                !
                Gw = Gw + ( D2/(dcmplx(wreal(iw1),eta)-wreal_read(iw2+1)) ) * dw2
                Gw = Gw + ( D3/(dcmplx(wreal(iw1),eta)-wreal_read(iw2+2)) ) * dw3
                !
            endif
            !
        enddo
        Gaux(iw1) = Gw
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !
    !print auxiliary Gf
    filename="Gaux.DAT"
    open(unit,file=reg(filename),form="formatted",status="unknown",position="rewind",action="write")
    do iw=1,Nreal
       write(unit,"(3E20.12)") wreal(iw),dreal(Gaux(iw))/pi,dimag(Gaux(iw))/pi
    enddo
    close(unit)
    !
    !compute Sigma
    allocate(Sigma(Nreal));Sigma=dcmplx(0d0,0d0)
    do iw=1,Nreal
        Sigma(iw) = dcmplx(wreal(iw),0d0) - 1d0/Gaux(iw)
    enddo
    !
    !print Sigma
    !write(filename,'(A9,1E10.3,A4)') "Sigma_eta",eta,".DAT"
    filename="Sigma.DAT"
    open(unit,file=reg(filename),form="formatted",status="unknown",position="rewind",action="write")
    do iw=1,Nreal
       write(unit,"(3E20.12)") wreal(iw),dreal(Sigma(iw)),dimag(Sigma(iw))
    enddo
    close(unit)
    !
    !
    !
contains
    !
    !
    !
    subroutine read_command_line
        integer :: exenamelength
        integer :: io, io2
        command_line = ""
        call get_command(command = command_line,status = io)
        !
        if(io==0)then
            call get_command_argument(0,length=exenamelength,status=io2)
            if(io2==0) then
                command_line = "&cmd "//adjustl(reg(command_line(exenamelength+1:)))//" /"
            else
                command_line = "&cmd "//adjustl(reg(command_line))//" /"
            endif
        else
          write(*,*) io,"Error getting command line."
        endif
        !
    end subroutine read_command_line
    !
    subroutine parse_command_line
        character(256) :: msg
        namelist /cmd/ filename, parameters, flag
        integer :: io
        !
        if(len_trim(command_line)>0)then
            msg = ''
            read(command_line,nml = cmd,iostat = io,iomsg = msg)
            if(io/=0)then
                print *,"Error parsing the command line or cmd.conf " // msg
                stop 
            endif
        endif
        !
    end subroutine parse_command_line
    !
    !
    !
end program Integrtor
