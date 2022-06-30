program pp_dispersion_fit
   !
   use module_container
   use utils_main
   use crystal, only : Nwig, Nvecwig
   implicit none
   !
   integer                                  :: TimeStart
   integer                                  :: unit,ndx,ik,ir,iorb,io,jo,Norb_dft,Npoints
   complex(8),allocatable                   :: Hr_mat(:,:,:),dumCorr(:,:,:)
   real(8),allocatable                      :: Hr_vec(:)
   integer,allocatable                      :: map(:,:),maplen(:)
   real(8),allocatable                      :: Ek_dft(:,:),Ek_dft_read(:,:),Ek_dft_refined(:,:),Ek_dft_read_unsort(:)
   real(8),allocatable                      :: Kpathaxis_dft_read(:),Kpathaxis_dft_refined(:)
   integer,allocatable                      :: order(:),used_ndx(:)
   integer                                  :: Niter!,Kpathaxis_refined
   real(8)                                  :: chi2,Kpathaxis_frac
   real(8)                                  :: Kvec(3),kx,ky,kz,tk,Ek_1,Ek_2
   real(8)                                  :: t(2,6),Ekmodel(2),Blat(3,3),csi,neta
   !
   !
   !
   !---------------------------------------------------------------------------!
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   call tick(TimeStart)
   call read_InputFile("input.in.fit")
   if(.not.Hmodel) stop
   if(reg(structure).eq."None") stop
   if(Hetero%status) stop
   !
   !
   ! Compute the starting dispersion from the Hr in the GWinput folder
   call initialize_Lattice(Crystal,1000)
   allocate(dumCorr(Crystal%Norb,Crystal%Norb,Crystal%Nkpt));dumCorr=czero
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="natcom",correction=dumCorr)
   !
   write(*,*) "Done1"
   !
   t(1,:) = [ -45.1, 157.8, 203.2,  25.7,  0.02, 0.48  ]/1000
   t(2,:) = [ 501.3, -15.9, 557.1, -72.0, -13.9, 12.2  ]/1000
   call get_Blat(Blat)
   !
   do ik=1,Crystal%Nkpt
      !
      Kvec =  Crystal%kpt(1,ik) * Blat(:,1) + Crystal%kpt(2,ik) * Blat(:,2) + Crystal%kpt(3,ik) * Blat(:,3)
      kx = Kvec(1)
      ky = Kvec(2)
      kz = Kvec(3)
      !
      csi = 0.5d0 * Crystal%kpt(1,ik) *2*pi
      neta = (sqrt(3d0)/2d0) * Crystal%kpt(2,ik) *2*pi
      !
      do io=1,2
         Ekmodel(io) = t(io,1) + t(io,2) * ( 2*cos(csi)*cos(neta)     + cos(2*csi) )   &
                               + t(io,3) * ( 2*cos(3*csi)*cos(neta)   + cos(2*neta) )   &
                               + t(io,4) * ( 2*cos(2*csi)*cos(2*neta) + cos(4*csi) )   &
                               + t(io,5) * ( cos(csi)*cos(3*neta) + cos(5*csi)*cos(neta) + cos(4*csi)*cos(2*neta) )   &
                               + t(io,6) * ( cos(3*csi)*cos(3*neta) + cos(6*csi) )
      enddo
      !
      Ek_1 = maxval(Ekmodel)
      Ek_2 = minval(Ekmodel)
      !
      tk = (Ek_1-Ek_2)/2d0
      !
      Crystal%Hk(1,1,ik) = Ek_1 - tk
      Crystal%Hk(2,2,ik) = Ek_2 + tk
      Crystal%Hk(1,2,ik) = tk
      Crystal%Hk(2,1,ik) = tk
      !
   enddo
   write(*,*) "Done2"
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="folded",correction=dumCorr,doplane=.true.)
   write(*,*) "Done3"
   allocate(Hr_mat(Crystal%Norb,Crystal%Norb,Nwig));Hr_mat=czero
   call wannier_K2R(Crystal%Nkpt3,Crystal%kpt,Crystal%Hk,Hr_mat)
   write(*,*) "Done4"
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Hr_folded.DAT",form="formatted",status="unknown",position="rewind",action="write")
   write(unit,*)                    !skip first line
   write(unit,"(1I5)") Crystal%Norb !Number of Wannier orbitals
   write(unit,"(1I5)") Nwig         !Number of Wigner-Seitz vectors
   do ir=1,Nwig
      do io=1,Crystal%Norb
         do jo=1,Crystal%Norb
            if(abs(Hr_mat(io,jo,iR)).gt.0d0)write(unit,"(5I5,2F12.6)") Nvecwig(:,iR), io, jo, dreal(Hr_mat(io,jo,iR)), dimag(Hr_mat(io,jo,iR))
         enddo
      enddo
   enddo
   close(unit)
   deallocate(Hr_mat)

   write(*,*) "Done5"


   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"model.DAT",form="formatted",status="unknown",position="rewind",action="write")
   do ik=1,Crystal%Nkpt_path
      !
      Kvec =  Crystal%kptpath(1,ik) * Blat(:,1) + Crystal%kptpath(2,ik) * Blat(:,2) + Crystal%kptpath(3,ik) * Blat(:,3)
      kx = Kvec(1)
      ky = Kvec(2)
      kz = Kvec(3)
      !
      csi = 0.5d0 * Crystal%kptpath(1,ik) *2*pi
      neta = (sqrt(3d0)/2d0) * Crystal%kptpath(2,ik) *2*pi
      !
      do io=1,2
         Ekmodel(io) = t(io,1) + t(io,2) * ( 2*cos(csi)*cos(neta)     + cos(2*csi) )   &
                               + t(io,3) * ( 2*cos(3*csi)*cos(neta)   + cos(2*neta) )   &
                               + t(io,4) * ( 2*cos(2*csi)*cos(2*neta) + cos(4*csi) )   &
                               + t(io,5) * ( cos(csi)*cos(3*neta) + cos(5*csi)*cos(neta) + cos(4*csi)*cos(2*neta) )   &
                               + t(io,6) * ( cos(3*csi)*cos(3*neta) + cos(6*csi) )
      enddo
      !
      write(unit,"(1I5,10000E20.12)") ik,Crystal%Kpathaxis(ik)/Crystal%Kpathaxis(Crystal%Nkpt_path),(Ekmodel(iorb),iorb=1,Crystal%Norb)
   enddo
   close(unit)


   stop












   allocate(dumCorr(Crystal%Norb,Crystal%Norb,Crystal%Nkpt));dumCorr=czero
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="estimate",correction=dumCorr)
   allocate(Hr_mat(Crystal%Norb,Crystal%Norb,Nwig));Hr_mat=czero
   call wannier_K2R(Crystal%Nkpt3,Crystal%kpt,Crystal%Hk,Hr_mat)
   call mat2vec(Hr_mat,Hr_vec,map,maplen) !the map is global and ever updated
   !
   !
   ! Read the benchmark dispersion
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Bands_input.DAT",form="formatted",status="old",position="rewind",action="read")
   read(unit,*) Npoints,Norb_dft
   if(Norb_dft.ne.Crystal%Norb) stop
   allocate(Ek_dft_read_unsort(Crystal%Norb),order(Crystal%Norb));Ek_dft_read_unsort=0d0
   allocate(Ek_dft_read(Crystal%Norb,Npoints));Ek_dft_read=0d0
   allocate(Kpathaxis_dft_read(Npoints));Kpathaxis_dft_read=-1d0
   allocate(used_ndx(Npoints));used_ndx=0d0
   do ik=1,Npoints
      !
      read(unit,*) Kpathaxis_dft_read(ik),(Ek_dft_read_unsort(iorb),iorb=1,Crystal%Norb)
      !
      call sort_array(Ek_dft_read_unsort,order)
      Ek_dft_read(:,ik) = Ek_dft_read_unsort(order)
      !
   enddo
   close(unit)
   deallocate(Ek_dft_read_unsort)
   !
   !
   ! Interpolate the benchmark dispersion to the same grid of Hr
   allocate(Ek_dft(Crystal%Norb,Crystal%Nkpt_path));Ek_dft=0d0
   do ik=1,Crystal%Nkpt_path
      do iorb=1,Crystal%Norb
         Ek_dft(iorb,ik) = cubic_interp( Kpathaxis_dft_read*Crystal%Kpathaxis(Crystal%Nkpt_path), Ek_dft_read(iorb,:), Crystal%Kpathaxis(ik) ) !/Crystal%Kpathaxis(Crystal%Nkpt_path)
      enddo
   enddo
   deallocate(Kpathaxis_dft_read,Ek_dft_read)
   !
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Bands_input_interp.DAT",form="formatted",status="unknown",position="rewind",action="write")
   do ik=1,Crystal%Nkpt_path
      write(unit,"(1I5,10000E20.12)") ik,Crystal%Kpathaxis(ik)/Crystal%Kpathaxis(Crystal%Nkpt_path),(Ek_dft(iorb,ik),iorb=1,Crystal%Norb)
   enddo
   close(unit)
   !
   !
   call fit_wrapper(chi2_Dispersion,Hr_vec,chi2,Niter,ftol=1d-15)
   write(*,"(A,I,A,F)")"     Iterations: ",Niter," Chi^2: ",chi2
   !
   !
   call vec2mat(Hr_vec,map,maplen,Hr_mat)
   call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat,Crystal%Hk)
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="fitted",correction=dumCorr)
   deallocate(Hr_vec,map,maplen)
   !
   !
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Hr_fitted.DAT",form="formatted",status="unknown",position="rewind",action="write")
   write(unit,*)                    !skip first line
   write(unit,"(1I5)") Crystal%Norb !Number of Wannier orbitals
   write(unit,"(1I5)") Nwig         !Number of Wigner-Seitz vectors
   do ir=1,Nwig
      do io=1,Crystal%Norb
         do jo=1,Crystal%Norb
            if(abs(Hr_mat(io,jo,iR)).gt.0d0)write(unit,"(5I5,2F12.6)") Nvecwig(:,iR), io, jo, dreal(Hr_mat(io,jo,iR)), dimag(Hr_mat(io,jo,iR))
         enddo
      enddo
   enddo
   close(unit)
   deallocate(Hr_mat)
   !
   !
   !
contains
   !
   !
   !
   subroutine mat2vec(Hr_mat_,Hr_vec_,map_,maplen_)
      implicit none
      complex(8),intent(in)                 :: Hr_mat_(:,:,:)
      real(8),allocatable,intent(out)       :: Hr_vec_(:)
      integer,allocatable,intent(out)       :: map_(:,:)
      integer,allocatable,intent(out)       :: maplen_(:)
      real(8),allocatable                   :: Hr_mat_flattened(:)
      integer                               :: Nelements,io,jo,iR,ndx
      !
      call assert_shape(Hr_mat_,[Crystal%Norb,Crystal%Norb,Nwig],"mat2vec","Hr_mat_")
      Nelements = Crystal%Norb*Crystal%Norb*Nwig
      allocate(Hr_mat_flattened(Nelements));Hr_mat_flattened=0d0
      !
      ndx=1
      do io=1,Crystal%Norb
         do jo=1,Crystal%Norb
            do iR=1,Nwig
               Hr_mat_flattened(ndx) = dreal(Hr_mat_(io,jo,iR))
               ndx = ndx + 1
            enddo
         enddo
      enddo
      !
      if(allocated(map_))deallocate(map_)
      if(allocated(maplen_))deallocate(maplen_)
      call get_pattern(map_,Hr_mat_flattened,1e4*eps,listDim=maplen_,IncludeSingle=.false.)
      !
      if(allocated(Hr_vec_))deallocate(Hr_vec_)
      allocate(Hr_vec_(size(maplen_)));Hr_vec_=0d0
      do io=1,size(maplen_)
         Hr_vec_(io) = Hr_mat_flattened(map_(io,1))
      enddo
      deallocate(Hr_mat_flattened)
      !
   end subroutine mat2vec
   !
   subroutine vec2mat(Hr_vec_,map_,maplen_,Hr_mat_)
      implicit none
      real(8),intent(in)                    :: Hr_vec_(:)
      integer,allocatable,intent(in)        :: map_(:,:)
      integer,allocatable,intent(in)        :: maplen_(:)
      complex(8),allocatable,intent(out)    :: Hr_mat_(:,:,:)
      real(8),allocatable                   :: Hr_mat_flattened(:)
      integer                               :: Nelements,io,jo,iR,ndx
      !
      if(.not.allocated(map_)) stop
      if(size(Hr_vec_).ne.size(map_,dim=1)) stop
      if(size(Hr_vec_).ne.size(maplen_)) stop
      !
      Nelements = Crystal%Norb*Crystal%Norb*Nwig
      allocate(Hr_mat_flattened(Nelements));Hr_mat_flattened=0d0
      !
      do io=1,size(maplen_)
         do jo=1,maplen_(io)
            Hr_mat_flattened(map_(io,jo)) = Hr_vec_(io)
         enddo
      enddo
      !
      allocate(Hr_mat_(Crystal%Norb,Crystal%Norb,Nwig));Hr_mat_=czero
      ndx=1
      do io=1,Crystal%Norb
         do jo=1,Crystal%Norb
            do iR=1,Nwig
               Hr_mat_(io,jo,iR) = dcmplx(Hr_mat_flattened(ndx),0d0)
               ndx = ndx + 1
            enddo
         enddo
      enddo
      !
      deallocate(Hr_mat_flattened)
      !
   end subroutine vec2mat
   !
   !
   subroutine chi2_Dispersion(Npara,Hr_vec_,chi2)
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: Hr_vec_
      real(8)                               :: chi2
      integer                               :: io,ik
      complex(8),allocatable                :: Hr_mat_(:,:,:)
      real(8),allocatable                   :: Ek(:,:)
      !
      ! Recompute dispersion with new Hr here
      write(*,*)Hr_vec_
      call vec2mat(Hr_vec_,map,maplen,Hr_mat_)
      call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat_,Crystal%Hk)
      call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.false.)!,corrname="fitted",correction=dumCorr)
      if(allocated(Ek))deallocate(Ek)
      Ek = Crystal%Ek
      !
      ! Estimate difference with interpolation
      chi2=0d0
      do io=1,Crystal%Norb
         do ik=1,Crystal%Nkpt_path
            chi2 = chi2 + abs(Ek(io,ik)-Ek_dft(io,ik))**2
         enddo
      enddo
      deallocate(Ek)
      !
   end subroutine chi2_Dispersion
   !
   !
   !
end program pp_dispersion_fit
