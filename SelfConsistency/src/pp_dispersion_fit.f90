program pp_dispersion_fit
   !
   use module_container
   use utils_main
   use crystal, only : Nwig
   use crystal, only : Rvecwig
   use crystal, only : radiuswig
   use crystal, only : Nvecwig
   use crystal, only : nrdegwig
   use crystal, only : calc_wignerseiz
   implicit none
   !
   type maptype
      integer,allocatable                   :: map(:,:)
      integer,allocatable                   :: maplen(:)
   end type maptype
   !
   integer                                  :: iorb,jorb,io,jo,Norb_dft,Npoints
   integer                                  :: unit,ik,iR,iwig,idist,idum,ndx,isite,jsite
   integer                                  :: ihop,maxdist,tirred,outerit,iouterit
   complex(8),allocatable                   :: Hr_mat(:,:,:),dumCorr(:,:,:)
   complex(8),allocatable                   :: cosinesEk_mat(:,:,:,:,:),cosinesHk_mat(:,:,:,:,:)
   real(8),allocatable                      :: Hr_vec(:),Hr_vec_old(:),thopping(:,:)
   integer,allocatable                      :: map(:,:),maplen(:)
   real(8),allocatable                      :: Ek_model(:,:),Ek_model_BZ(:,:)
   real(8),allocatable                      :: Eo_model_BZ(:),Rot(:,:)
   real(8),allocatable                      :: Ek_dft(:,:),Ek_dft_read(:,:)
   real(8),allocatable                      :: Kpathaxis_dft_read(:),radiuswig_uc(:),radiuswig_sorted(:)
   real(8),allocatable                      :: cosinesEk(:,:),cosinesHk(:,:),Rvec(:,:)
   complex(8),allocatable                   :: Hr_mat_layered(:,:,:),Hk_mat_layered(:,:,:)
   integer,allocatable                      :: orderE(:),orderK(:),orderR(:),Zrad(:),ZradPlane(:)
   integer                                  :: Niter,Qst,Rst,i,Nwig_,zerondx,zrange,Npoints_limit
   integer                                  :: Nwig_plane,Nkpt_plane,Nkpt3_plane(3),Nbilayer,ibilayer
   integer                                  :: ikx,iky,ikz,ik1,ik2
   complex(8),allocatable                   :: Hk_intralayer(:,:),Hk_interlayer(:,:)
   real(8)                                  :: chi2,hoppFactor,est,fto,kR,Kfact
   real(8)                                  :: Kvec(3),kx,ky,kz,tk,Ek_1,Ek_2
   real(8)                                  :: t(2,6),Ekmodel(2),Blat(3,3),Bvec(3),csi,neta
   logical                                  :: dofolded,dofit,domanip,filexists,justprint,verbs,optimize,Hkft
   integer                                  :: ndxmanip(10)
   real(8)                                  :: valmanip(10)
   integer,parameter                        :: W90NumCol=15
   complex(8)                               :: Mat(2,2)
   real(8)                                  :: Emat(2)
   !
   type(maptype),allocatable                :: mapOrb(:)
   !
   !
   dofolded=.false.
   justprint=.true.
   domanip=.false.
   !
   verbs=.false.
   !
   !
   !---------------------------------------------------------------------------!
   !     READING INPUT FILE, INITIALIZING OMP, AND CHECK FOLDER STRUCTURE      !
   !---------------------------------------------------------------------------!
   call read_InputFile("input.in.fit")
   if(.not.Hmodel) stop
   if(reg(structure).eq."None") stop
   if(Hetero%status) stop
   !
   call parse_input_variable(hoppFactor, "HOPPFACT", "input.in.fit", default=1d0 )
   call parse_input_variable(hoppFactor, "HOPPFACT", "input.in.fit", default=1d0 )
   call parse_input_variable(ndxmanip, "NDXMANIP", "input.in.fit", default=[0,0,0,0,0,0,0,0,0,0] )
   call parse_input_variable(valmanip, "VALMANIP", "input.in.fit", default=[0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0] )
   call parse_input_variable(est, "EST", "input.in.fit", default=0.5d0 )
   call parse_input_variable(fto, "FTO", "input.in.fit", default=1d-15 )
   call parse_input_variable(dofit, "DOFIT", "input.in.fit", default=.false. )
   call parse_input_variable(outerit, "OUTIT", "input.in.fit", default=1 )
   call parse_input_variable(optimize, "OPT", "input.in.fit", default=.false. )
   call parse_input_variable(Hkft, "HKFT", "input.in.fit", default=.false. )
   call parse_input_variable(Nbilayer, "NBILAYER", "input.in.fit", default=3 )
   !
   !
   call set_lattice(LatticeVec,ucVec)
   Crystal%Nkpt3 = Nkpt3
   Crystal%Nkpt = product(Nkpt3)
   Crystal%Norb = Norb_model*Nsite
   Crystal%Nsite = Nsite
   call build_kpt(Crystal%Nkpt3,Crystal%kpt,pathOUTPUT=reg(pathINPUT))
   if(reg(structure).eq."User") call set_UserPath(UserPath)
   call calc_Kpath(Crystal%kptpath,reg(structure),Nkpt_path,Crystal%Kpathaxis,Crystal%KpathaxisPoints)
   Crystal%Nkpt_path = size(Crystal%kptpath,dim=2)
   allocate(Crystal%Hk(Crystal%Norb,Crystal%Norb,Crystal%Nkpt));Crystal%Hk=czero
   allocate(Crystal%Ek(Crystal%Norb,Crystal%Nkpt));Crystal%Hk=czero
   allocate(Crystal%Hloc(Crystal%Norb,Crystal%Norb));Crystal%Hloc=czero
   Crystal%status=.true.
   call get_Blat(Blat)
   !
   !
   !---------------------------------------------------------------------------!
   !                       READ-IN BENCHMARK DISPERSION                        !
   !---------------------------------------------------------------------------!
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Bands_input.DAT",form="formatted",status="old",position="rewind",action="read")
   read(unit,*) Npoints,Norb_dft,Npoints_limit,Kfact
   if(Crystal%Norb.ne.Norb_dft) stop "Crystal%Norb.ne.Norb_dft (Bands_input.DAT)"
   allocate(orderE(Crystal%Norb),orderK(Npoints))
   allocate(Ek_dft_read(Crystal%Norb,Npoints));Ek_dft_read=0d0
   allocate(Kpathaxis_dft_read(Npoints));Kpathaxis_dft_read=-1d0
   !
   !
   !sorting along energy and kpoints
   do ik=1,Npoints
      !
      read(unit,*) Kpathaxis_dft_read(ik),(Ek_dft_read(iorb,ik),iorb=1,Crystal%Norb)
      !
      if(ik.le.Npoints_limit) Kpathaxis_dft_read(ik) = Kpathaxis_dft_read(ik)*Kfact
      Ek_dft_read(:,ik) = Ek_dft_read(:,ik) - look4dens%mu
      !
      call sort_array(Ek_dft_read(:,ik),orderE)
      Ek_dft_read(:,ik) = Ek_dft_read(orderE,ik)
      !
   enddo
   close(unit)
   call sort_array(Kpathaxis_dft_read,orderK)
   Kpathaxis_dft_read=Kpathaxis_dft_read(orderK)
   Ek_dft_read=Ek_dft_read(:,orderK)
   deallocate(orderK,orderE)
   !
   open(unit,file=reg(pathINPUT)//"Bands_input_ordered.DAT",form="formatted",status="unknown",position="rewind",action="write")
   do ik=1,Npoints
      write(unit,"(10F)") Kpathaxis_dft_read(ik),(Ek_dft_read(iorb,ik),iorb=1,Crystal%Norb)
   enddo
   close(unit)
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
   !---------------------------------------------------------------------------!
   !                   CREATE A DISPERSION MODEL SUM OF COSINES                !
   !---------------------------------------------------------------------------!
   !read in the starting hopping t0,t1,...
   unit = free_unit()
   open(unit,file="./hopping.DAT",form="formatted",status="old",position="rewind",action="read")
   read(unit,*) maxdist, Norb_dft
   if(Crystal%Norb.ne.Norb_dft) stop "Crystal%Norb.ne.Norb_dft (hopping.DAT)"
   allocate(thopping(maxdist, Crystal%Norb));thopping=0d0
   do idist=1,maxdist
      read(unit,*) idum,(thopping(idist,iorb),iorb=1,Crystal%Norb)
      if(idum.ne.idist) stop "idum.ne.idist (hopping.DAT)"
   enddo
   !
   !
   !create real-space lattice and divide it into shells with fixed distance
   call calc_wignerseiz(Crystal%Nkpt3)
   !
   allocate(orderR(Nwig))
   call sort_array(radiuswig,orderR)
   radiuswig_sorted=radiuswig(orderR)
   call get_pattern(map,radiuswig_sorted,1e4*eps,listDim=maplen,IncludeSingle=.true.)
   if(maxdist.gt.size(maplen)) stop "maxdist.gt.size(maplen) increase kpoints"
   !
   !
   !create cosines functions along the path on each distance
   allocate(cosinesEk(maxdist,Crystal%Nkpt_path));cosinesEk=0d0
   allocate(cosinesEk_mat(maxdist,Crystal%Nkpt,2,2,0:2));cosinesEk_mat=czero
   do ik=1,Crystal%Nkpt_path
      do idist=1,maxdist
         if(ik.eq.1) write(*,"(A,I)") "distance: ",idist
         do iR=1,maplen(idist)
            iwig = map(idist,iR)
            !if(ik.eq.1) write(*,"(2(A,I4),A,2F12.6,2(A,1F12.6),A,3I4)") "iR: ",iR," radius index: ",iwig," hopping: ",thopping(idist,:),"    radius value: ",radiuswig_sorted(iwig),"    Rlat value: ",radiuswig(Zrad(orderR(iwig)))," Nvecwig: ",Nvecwig(:,Zrad(orderR(iwig)))
            if(ik.eq.1) write(*,"(2(A,I4),A,2F12.6,1(A,1F12.6),A,3I4)") "iR: ",iR," radius index: ",iwig," hopping: ",thopping(idist,:),"    radius value: ",radiuswig_sorted(iwig)," Nvecwig: ",Nvecwig(:,orderR(iwig))
            kR = 2*pi * dot_product(Crystal%kptpath(:,ik),Nvecwig(:,orderR(iwig)))
            cosinesEk(idist,ik) = cosinesEk(idist,ik) + cos(kR)
         enddo
      enddo
   enddo
   cosinesEk(1,:) = 1d0
   !
   !
   !create cosines functions in the full BZ on each distance
   allocate(cosinesHk(maxdist,Crystal%Nkpt));cosinesHk=0d0
   allocate(cosinesHk_mat(maxdist,Crystal%Nkpt,2,2,0:2));cosinesHk_mat=czero
   do ik=1,Crystal%Nkpt
      do idist=1,maxdist
         do iR=1,maplen(idist)
            iwig = map(idist,iR)
            kR = 2*pi * dot_product(Crystal%kpt(:,ik),Nvecwig(:,orderR(iwig)))
            cosinesHk(idist,ik) = cosinesHk(idist,ik) + cos(kR)
         enddo
      enddo
   enddo
   cosinesHk(1,:) = 1d0
   deallocate(orderR,radiuswig_sorted,map,maplen)
   !
   !
   !create parameter vector with only non-vanishing hoppings separately for each band
   allocate(mapOrb(Crystal%Norb))
   tirred=0
   do iorb=1,Crystal%Norb
      !
      call get_pattern(mapOrb(iorb)%map,thopping(:,iorb),1e4*eps,listDim=mapOrb(iorb)%maplen,IncludeSingle=.true.)
      !check on the map
      write(*,"(A)") new_line("(A)")//" orbital #"//str(iorb)//" hopping read from file"
      do idist=1,maxdist
         write(*,"(10000E20.12)") thopping(idist,iorb)
      enddo
      write(*,"(A)") " map"
      do idist=1,size(mapOrb(iorb)%maplen)
         write(*,"(10000E20.12)") (thopping(mapOrb(iorb)%map(idist,ihop),iorb), ihop=1,mapOrb(iorb)%maplen(idist))
      enddo
      !
      tirred = tirred + size(mapOrb(iorb)%maplen)
      !
   enddo
   write(*,"(A)") " tirred: "//str(tirred)
   allocate(Hr_vec(tirred));Hr_vec=0d0
   ndx=1
   do iorb=1,Crystal%Norb
      do idist=1,size(mapOrb(iorb)%maplen)
         Hr_vec(ndx) = thopping(mapOrb(iorb)%map(idist,1),iorb)
         write(*,*) Hr_vec(ndx)
         ndx = ndx + 1
      enddo
   enddo
   !
   !
   call calc_dispersion(Ek_model,Hr_vec,"initial",.true.)
   !
   !
   ! run the minimization
   chi2=0d0
   if(dofit)then
      !
      do iouterit=1,outerit
         call fit_wrapper(chi2_DispersionModel,Hr_vec,chi2,Niter,ftol=fto,estm=est)
         write(*,"(A,I,A,F)")str(iouterit)//" Iterations: ",Niter," Chi^2: ",chi2
      enddo
      !
   endif
   call calc_dispersion(Ek_model,Hr_vec,"final",.true.)
   !
   !
   ! get updated thopping and build the Ek in the full BZ
   call vec2vec(Hr_vec,thopping)
   allocate(Ek_model_BZ(Crystal%Norb,Crystal%Nkpt));Ek_model_BZ=0d0
   allocate(Eo_model_BZ(Crystal%Norb));Eo_model_BZ=0d0
   do ik=1,Crystal%Nkpt
      do iorb=1,Crystal%Norb
         do idist=1,maxdist
            Ek_model_BZ(iorb,ik) = Ek_model_BZ(iorb,ik) + thopping(idist,iorb) * cosinesHk(idist,ik) !<=====
         enddo
         Crystal%Ek(iorb,ik) = Ek_model_BZ(iorb,ik)
      enddo
   enddo
   do iorb=1,Crystal%Norb
      Eo_model_BZ(iorb) = sum(Ek_model_BZ(iorb,:))/Crystal%Nkpt
      write(*,"(A,F)")" Eo_"//str(iorb)//"= ",Eo_model_BZ(iorb)
   enddo
   !
   !
   !rotate with B-AB basis
   allocate(Rot(Crystal%Norb,Crystal%Norb));Rot=1d0/sqrt(2d0)
   Rot(2,2)=-1d0/sqrt(2d0)
   do ik=1,Crystal%Nkpt
      Crystal%Hk(:,:,ik) = rotate(diag(Ek_model_BZ(:,ik)),Rot)
      Crystal%Hloc = Crystal%Hloc + Crystal%Hk(:,:,ik)/Crystal%Nkpt
   enddo
   do iorb=1,Crystal%Norb
      write(*,"(100F)")(dreal(Crystal%Hloc(iorb,jorb)),jorb=1,Crystal%Norb)
   enddo
   write(*,*)
   do iorb=1,Crystal%Norb
      write(*,"(100F)")(dimag(Crystal%Hloc(iorb,jorb)),jorb=1,Crystal%Norb)
   enddo
   !
   !
   ! compute chemical potential
   call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
   Glat%mu=0d0
   if(look4dens%TargetDensity.ne.0d0)then
      look4dens%mu=0d0
      call calc_Gmats(Glat,Crystal)
      call set_density(Glat,Crystal,look4dens)
      call calc_density(Glat,Glat%N_s)
      densityLDA = Glat%N_s
      call dump_Matrix(densityLDA,reg(pathINPUT),"Nlda",paramagnet)
   endif
   !
   !
   ! reprint the bands coming from the fitting
   if(allocated(dumCorr))deallocate(dumCorr)
   allocate(dumCorr(Crystal%Norb,Crystal%Norb,Crystal%Nkpt));dumCorr=czero
   dumCorr=czero
   allocate(Hr_mat(Crystal%Norb,Crystal%Norb,Nwig));Hr_mat=czero
   call wannier_K2R(Crystal%Nkpt3,Crystal%kpt,Crystal%Hk,Hr_mat)
   !
   !
   !cleanup spurious complex stuff
   where(abs(Hr_mat)<eps)Hr_mat=czero
   do iwig=1,Nwig
      do iorb=1,Crystal%Norb
         do jorb=1,Crystal%Norb
            if(dimag(Hr_mat(iorb,jorb,iwig)).lt.eps) Hr_mat(iorb,jorb,iwig) = dcmplx(dreal(Hr_mat(iorb,jorb,iwig)),0d0)
         enddo
      enddo
   enddo
   !
   !
   !go back to k space
   call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat,Crystal%Hk)
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="model_final",correction=dumCorr,doplane=FermiSurf)
   !
   !
   ! Ft along Kz
   if(Hkft)then
      Nkpt3_plane = [Nkpt3(1),Nkpt3(2),1]
      Nkpt_plane = product(Nkpt3_plane)
      allocate(Hk_mat_layered(Crystal%Norb*Nbilayer,Crystal%Norb*Nbilayer,Nkpt_plane));Hk_mat_layered=czero
      allocate(Hk_interlayer(Crystal%Norb,Crystal%Norb));Hk_interlayer=czero
      allocate(Hk_intralayer(Crystal%Norb,Crystal%Norb));Hk_intralayer=czero
      ik1=1
      ik2=1
      do ikx=1,Nkpt3(1)
         do iky=1,Nkpt3(2)
            !
            Hk_interlayer = czero
            Hk_intralayer = czero
            do ikz=1,Nkpt3(3)
               !
               kR = 2*pi * Crystal%kpt(3,ikz) * 0d0
               Hk_intralayer = Hk_intralayer + Crystal%Hk(:,:,ik1) * dcmplx(cos(kR),+sin(kR)) / Nkpt3(3)
               !
               kR = 2*pi * Crystal%kpt(3,ikz) * 1d0
               Hk_interlayer = Hk_interlayer + Crystal%Hk(:,:,ik1) * dcmplx(cos(kR),+sin(kR)) / Nkpt3(3)
               !
               ik1=ik1+1
               !
            enddo
            !
            Hk_mat_layered(1:Crystal%Norb,1:Crystal%Norb,ik2) = Hk_intralayer
            Hk_mat_layered(1:Crystal%Norb,1+Crystal%Norb:2*Crystal%Norb,ik2) = Hk_interlayer
            Hk_mat_layered(1+Crystal%Norb:2*Crystal%Norb,1:Crystal%Norb,ik2) = Hk_interlayer
            ik2=ik2+1
            !
         enddo
      enddo
      !
      !
      do ibilayer=2,Nbilayer
         !
         io = 1 + (ibilayer-1)*Crystal%Norb
         jo = ibilayer*Crystal%Norb
         Hk_mat_layered(io:jo,io:jo,:) = Hk_mat_layered(1:Crystal%Norb,1:Crystal%Norb,:)
         !
         Hk_mat_layered(io:jo,io-Crystal%Norb:jo-Crystal%Norb,:) = Hk_mat_layered(1:Crystal%Norb,1+Crystal%Norb:2*Crystal%Norb,:)
         Hk_mat_layered(io-Crystal%Norb:jo-Crystal%Norb,io:jo,:) = Hk_mat_layered(1+Crystal%Norb:2*Crystal%Norb,1:Crystal%Norb,:)
         !
      enddo
      !
      !
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"Hk_built.DAT",form="formatted",status="unknown",position="rewind",action="write")
      write(unit,("(3I10)")) 1,Nkpt_plane,Crystal%Norb*Nbilayer
      do ik=1,Nkpt_plane
         write(unit,("(2I6,3F14.8)")) 1,ik,Crystal%kpt(:,1+(ik-1)*Nkpt3(3))
         do iorb=1,Crystal%Norb*Nbilayer
            do jorb=1,Crystal%Norb*Nbilayer
               write(unit,("(2I4,2E20.12)")) iorb,jorb,dreal(Hk_mat_layered(iorb,jorb,ik)),dimag(Hk_mat_layered(iorb,jorb,ik))
            enddo
         enddo
      enddo
      close(unit)
      !
      stop "message"
      !
   endif
   !
   !
   !optimize the tuned real-space structure after the modifications
   if(optimize)then
      !
      ! remove unwanted hoppings
      do iwig=1,Nwig
         if(Nvecwig(3,iwig).eq.1)then
            Hr_mat(1,1,iwig)=czero;Hr_mat(2,2,iwig)=czero
            Hr_mat(2,1,iwig)=czero
         elseif(Nvecwig(3,iwig).eq.-1)then
            Hr_mat(1,1,iwig)=czero;Hr_mat(2,2,iwig)=czero
            Hr_mat(1,2,iwig)=czero
         endif
      enddo
      !
      call mat2vec(Hr_mat,Hr_vec,map,maplen) !the map is global and never updated
      !
      call inquireFile("./Hr_vec.DAT",filexists,hardstop=.false.)
      if(filexists)then
         !
         unit = free_unit()
         open(unit,file="./Hr_vec.DAT",form="formatted",status="old",position="rewind",action="read")
         do iorb=1,size(Hr_vec)
            read(unit,"(1E20.12)") Hr_vec(iorb)
         enddo
         close(unit)
         !
         write(*,*)"updated from file Hr_vec old -> updated"
         !
      endif
      !
      write(*,*)"starting Hr_vec"
      do iorb=1,size(Hr_vec)
         write(*,"(1E20.12)") Hr_vec(iorb)
      enddo
      !
      chi2=0d0
      if(dofit)then
         call fit_wrapper(chi2_Dispersion,Hr_vec,chi2,Niter,ftol=fto,estm=est)
         write(*,"(A,I,A,F)")"     Iterations: ",Niter," Chi^2: ",chi2
      endif
      !
      write(*,*)"final Hr_vec"
      do iorb=1,size(Hr_vec)
         write(*,"(1E20.12)") Hr_vec(iorb)
      enddo
      !
      Hr_mat=czero
      call vec2mat(Hr_vec,map,maplen,Hr_mat)
      where(abs(Hr_mat)<eps)Hr_mat=czero
      do iwig=1,Nwig
         do iorb=1,Crystal%Norb
            do jorb=1,Crystal%Norb
               if(dimag(Hr_mat(iorb,jorb,iwig)).lt.eps) Hr_mat(iorb,jorb,iwig) = dcmplx(dreal(Hr_mat(iorb,jorb,iwig)),0d0)
            enddo
         enddo
         if(Nvecwig(3,iwig).eq.1)then
            Hr_mat(1,1,iwig)=czero;Hr_mat(2,2,iwig)=czero
            Hr_mat(2,1,iwig)=czero
         elseif(Nvecwig(3,iwig).eq.-1)then
            Hr_mat(1,1,iwig)=czero;Hr_mat(2,2,iwig)=czero
            Hr_mat(1,2,iwig)=czero
         endif
      enddo
      !
      call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat,Crystal%Hk)
      !
      if(allocated(dumCorr))deallocate(dumCorr)
      allocate(dumCorr(Crystal%Norb,Crystal%Norb,Crystal%Nkpt));dumCorr=czero
      call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="optimized_final",correction=dumCorr,doplane=FermiSurf)
      !call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,doplane=FermiSurf)
      !
   endif
   !
   !
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,doplane=FermiSurf)
   !
   !
   !extract real space hamiltonian for the single layer - 11,12
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Hr_fitted_11.DAT",form="formatted",status="unknown",position="rewind",action="write")
   write(unit,*)                    !skip first line
   write(unit,"(1I5)") Crystal%Norb !Number of Wannier orbitals
   write(unit,"(1I5)") Nwig         !Number of Wigner-Seitz vectors
   Qst = int(Nwig/W90NumCol)
   Rst = mod(Nwig,W90NumCol)
   do i=1,Qst
     write(unit,"("//str(W90NumCol)//"I5)")(1,j=1,W90NumCol)
   enddo
   if(Rst.ne.0)write(unit,"("//str(Rst)//"I5)")(1,j=1,Rst)
   Nwig_=0
   do ir=1,Nwig
      if(Nvecwig(3,iR).eq.0)then
          write(unit,"(5I5,2E20.12)") Nvecwig(:,iR), 1, 1, dreal(Hr_mat(1,1,iR)), dimag(Hr_mat(1,1,iR))
          Nwig_=Nwig_+1
          if(abs(Hr_mat(1,2,iR)).gt.eps)then
             if(.not.allocated(Zrad)) then
                Zrad = [ iR ]
             else
                Zrad = [ Zrad, iR ]
             endif
          endif
      else
          if(abs(Hr_mat(1,2,iR)).gt.eps)write(*,"(5I5,3E20.12)") Nvecwig(:,iR), 1, 2, radiuswig(iR), dreal(Hr_mat(1,2,iR)), dimag(Hr_mat(1,2,iR))
      endif
   enddo
   close(unit)
   write(*,*)"Nwig_",Nwig_,int(Nwig_/W90NumCol)
   write(*,*)
   write(*,*)
   write(*,*)
   !
   call get_pattern(map,radiuswig(Zrad),1e4*eps,listDim=maplen,IncludeSingle=.true.)
   Zrad = [ Zrad(map(:,1)) ]
   zrange = size(Zrad)
   !
   allocate(orderR(zrange));orderR=0
   call sort_array(radiuswig(Zrad),orderR)
   !
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Hr_fitted_12_b.DAT",form="formatted",status="unknown",position="rewind",action="write")
   do idist=1,zrange
      iR = Zrad(orderR(idist))
      write(unit,"(3I5,5E20.12)")  Nvecwig(:,iR), dreal(Hr_mat(1,2,iR)), radiuswig(iR)
   enddo
   close(unit)
   deallocate(orderR,Zrad)
   !
   !
   ! allocate real-space Hamiltonian for the multilayer. HARD CODED
   allocate(Hr_mat_layered(Crystal%Norb*Nbilayer,Crystal%Norb*Nbilayer,Nwig));Hr_mat_layered=czero
   !
   !
   !extract real space hamiltonian for the bilayer
   !first check how many Nwig to use
   Nwig_=0
   do iR=1,Nwig
      if(Nvecwig(3,iR).eq.0)then
         orb_p_search: do iorb=1,Crystal%Norb
            do jorb=1,Crystal%Norb
               if(abs(Hr_mat(iorb,jorb,iR)).gt.eps)then
                  if(.not.allocated(Zrad)) then
                     Zrad = [ iR ]
                     Nwig_ = Nwig_ + 1
                     exit orb_p_search
                  else
                     Zrad = [ Zrad, iR ]
                     Nwig_ = Nwig_ + 1
                     exit orb_p_search
                  endif
               endif
            enddo
         enddo orb_p_search
      endif
   enddo
   write(*,*)"Nwig_ plane",Nwig_,size(Zrad)
   Nwig_plane = Nwig_
   !
   !then write to file
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Hr_fitted_11_22.DAT",form="formatted",status="unknown",position="rewind",action="write")
   write(unit,*)                    !skip first line
   write(unit,"(1I5)") Crystal%Norb !Number of Wannier orbitals
   write(unit,"(1I5)") Nwig_         !Number of Wigner-Seitz vectors
   Qst = int(Nwig_/W90NumCol)
   Rst = mod(Nwig_,W90NumCol)
   do i=1,Qst
     write(unit,"("//str(W90NumCol)//"I5)")(nrdegwig(j+(i-1)*W90NumCol),j=1,W90NumCol)
   enddo
   if(Rst.ne.0)write(unit,"("//str(Rst)//"I5)")(nrdegwig(j+Qst*W90NumCol),j=1,Rst)
   do iR=1,Nwig_
      iwig = Zrad(iR)
      write(*,*)iwig
      do iorb=1,Crystal%Norb
         do jorb=1,Crystal%Norb
            write(unit,"(5I5,2E20.12)") Nvecwig(:,iwig), iorb, jorb, dreal(Hr_mat(iorb,jorb,iwig)), dimag(Hr_mat(iorb,jorb,iwig))
            !
            Hr_mat_layered(iorb,jorb,iwig) = Hr_mat(iorb,jorb,iwig)
            do ibilayer=2,Nbilayer
               io = iorb+(ibilayer-1)*2
               jo = jorb+(ibilayer-1)*2
               write(*,*)io,jo, Hr_mat(iorb,jorb,iwig)
               Hr_mat_layered(io,jo,iwig) = Hr_mat(iorb,jorb,iwig)
            enddo
            !
         enddo
      enddo
   enddo
   close(unit)
   ZradPlane = Zrad
   deallocate(Zrad)
   !
   !
   !extract real space connection between bilayers
   Nwig_=0
   do iR=1,Nwig
      if(Nvecwig(3,iR).ne.0)then
         orb_v_search: do iorb=1,Crystal%Norb
            do jorb=1,Crystal%Norb
               if(abs(Hr_mat(iorb,jorb,iR)).gt.eps)then
                  if(.not.allocated(Zrad)) then
                     Zrad = [ iR ]
                     Nwig_ = Nwig_ + 1
                     exit orb_v_search
                  else
                     Zrad = [ Zrad, iR ]
                     Nwig_ = Nwig_ + 1
                     exit orb_v_search
                  endif
               endif
            enddo
         enddo orb_v_search
      endif
   enddo
   write(*,*)"Nwig_ vert",Nwig_,size(Zrad)
   !
   call get_pattern(map,radiuswig(Zrad),1e4*eps,listDim=maplen,IncludeSingle=.true.)
   !Zrad = [ Zrad(map(:,1)) ]
   zrange = size(maplen)
   !
   allocate(orderR(zrange));orderR=0
   call sort_array(radiuswig(Zrad(map(:,1))),orderR)
   !
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Hr_fitted_12.DAT",form="formatted",status="unknown",position="rewind",action="write")
   do idist=1,zrange
      iR = Zrad(map(orderR(idist),1))
      write(unit,"(3I5,50E20.12)")  Nvecwig(:,iR), dreal(Hr_mat(1,2,iR)), dreal(Hr_mat(2,1,iR)), radiuswig(iR)
   enddo
   write(unit,*)
   write(unit,*)"---"
   write(unit,*)
   do iR=1,Nwig_
      do iorb=1,Crystal%Norb
         do jorb=1,Crystal%Norb
            write(unit,"(5I5,2E20.12)") Nvecwig(:,Zrad(iR)), iorb, jorb, dreal(Hr_mat(iorb,jorb,Zrad(iR))), dimag(Hr_mat(iorb,jorb,Zrad(iR)))
            !
            ikz = Nvecwig(3,Zrad(iR))
            iwig = find_vec([ Nvecwig(1,Zrad(iR)), Nvecwig(2,Zrad(iR)), 0 ],Nvecwig,hardstop=.true.)
            if(optimize)then
               if(abs(Hr_mat(iorb,jorb,Zrad(iR))).gt.eps)then
                  if(iorb.eq.jorb) stop " iorb.eq.jorb "
                  do ibilayer=1,Nbilayer-1
                     io = (iorb+1)+(ibilayer-1)*2
                     jo = (jorb+1)+(ibilayer-1)*2
                     write(*,*)io,jo, Hr_mat(iorb,jorb,Zrad(iR))
                     Hr_mat_layered(io,jo,iwig) = Hr_mat(iorb,jorb,Zrad(iR))
                  enddo
               endif
            else
               if(ikz.eq.0) stop " ikz.eq.0 "
               do ibilayer=2,Nbilayer-1
                  io = iorb + (ibilayer-1)*2
                  jo = jorb + (ibilayer-1)*2 + ikz*2
                  write(*,*)io,jo, Hr_mat(iorb,jorb,Zrad(iR))
                  Hr_mat_layered(io,jo,iwig) = Hr_mat(iorb,jorb,Zrad(iR))
                  Hr_mat_layered(jo,io,iwig) = Hr_mat(iorb,jorb,Zrad(iR))
               enddo
            endif
            !
         enddo
      enddo
   enddo
   close(unit)
   deallocate(orderR,Zrad)
   !
   !
   !
   !then write to file
   unit = free_unit()
   open(unit,file=reg(pathINPUT)//"Hr_multilayer.DAT",form="formatted",status="unknown",position="rewind",action="write")
   write(unit,*)                    !skip first line
   write(unit,"(1I5)") Crystal%Norb*Nbilayer
   write(unit,"(1I5)") Nwig_plane         !Number of Wigner-Seitz vectors
   Qst = int(Nwig_plane/W90NumCol)
   Rst = mod(Nwig_plane,W90NumCol)
   do i=1,Qst
     write(unit,"("//str(W90NumCol)//"I5)")(1,j=1,W90NumCol)
   enddo
   if(Rst.ne.0)write(unit,"("//str(Rst)//"I5)")(1,j=1,Rst)

   do iR=1,Nwig_plane
      do iorb=1,Crystal%Norb*Nbilayer
         do jorb=1,Crystal%Norb*Nbilayer
            write(unit,"(5I5,2E20.12)") Nvecwig(:,ZradPlane(iR)), iorb, jorb, dreal(Hr_mat_layered(iorb,jorb,ZradPlane(iR))), dimag(Hr_mat_layered(iorb,jorb,ZradPlane(iR)))
         enddo
      enddo
   enddo
   close(unit)
   !






























   stop








































































   stop

   !
   !
   !---------------------------------------------------------------------------!
   !              DISPERSION COMING FROM PRB 97, 115118 (2018)                 !
   !---------------------------------------------------------------------------!
   if(dofolded)then
      call initialize_Lattice(Crystal,1000)
      allocate(dumCorr(Crystal%Norb,Crystal%Norb,Crystal%Nkpt));dumCorr=czero
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
      dumCorr=czero
      call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.false.,corrname="folded",correction=dumCorr,doplane=.true.)
      write(*,*) "printed folded Dispersion Ek -> Hr -> Ek "
      !
      !
      allocate(Hr_mat(Crystal%Norb,Crystal%Norb,Nwig));Hr_mat=czero
      call wannier_K2R(Crystal%Nkpt3,Crystal%kpt,Crystal%Hk,Hr_mat)
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
      write(*,*) "printed folded Hr "
      !
      !
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"model_Ek_folded.DAT",form="formatted",status="unknown",position="rewind",action="write")
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
      write(*,*) "printed folded Ek model dispersion "
   endif
   !
   !
   !
   !---------------------------------------------------------------------------!
   !                  OPTIMIZING DISPERSION FROM ncomms11043                   !
   !---------------------------------------------------------------------------!
   if(dofolded)call DeallocateLattice(Crystal)
   call initialize_Lattice(Crystal,1000)
   if(allocated(dumCorr))deallocate(dumCorr)
   allocate(dumCorr(Crystal%Norb,Crystal%Norb,Crystal%Nkpt));dumCorr=czero
   dumCorr=czero
   call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="initialized",correction=dumCorr,doplane=FermiSurf)
   write(*,*) "initialized Hr "
   do ik=1,Crystal%Nkpt_path
      call check_Hermiticity(Crystal%Hk(:,:,ik),eps,enforce=.true.,hardstop=.true.,name="Hk_k"//str(ik),verb=.true.)
   enddo
   !
   !
   ! Read the benchmark dispersion

   !
   !
   ! rotate the interpolated Hk with the rotation of natcom
   if(domanip)then
      !
      ! creating a map of the irreducible elements from Hr_mat to Hr_vec
      allocate(Hr_mat(Crystal%Norb,Crystal%Norb,Nwig));Hr_mat=czero
      call wannier_K2R(Crystal%Nkpt3,Crystal%kpt,Crystal%Hk,Hr_mat)
      call mat2vec(Hr_mat,Hr_vec,map,maplen) !the map is global and never updated
      !
      !
      !manipulate elementwise
      write(*,"(10F)")Hr_vec
      do iorb=1,6
         if(ndxmanip(iorb).eq.1)Hr_vec(iorb) = valmanip(iorb)
      enddo
      write(*,"(10F)")Hr_vec
      !
      Hr_vec = Hr_vec * hoppFactor
      write(*,"(10F)")Hr_vec
      !
      ! get back Hk
      call vec2mat(Hr_vec,map,maplen,Hr_mat)
      call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat,Crystal%Hk)
      call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="manip",correction=dumCorr)
      deallocate(Hr_vec,map,maplen)
      !
   endif
   !
   !
   ! run the minimization
   if(dofit)then
      !
      ! creating a map of the irreducible elements from Hr_mat to Hr_vec
      allocate(Hr_mat(Crystal%Norb,Crystal%Norb,Nwig));Hr_mat=czero
      call wannier_K2R(Crystal%Nkpt3,Crystal%kpt,Crystal%Hk,Hr_mat)
      call mat2vec(Hr_mat,Hr_vec,map,maplen) !the map is global and never updated
      !
      !
      call inquireFile("./Hr_vec.DAT",filexists,hardstop=.false.)
      if(filexists)then
         !
         Hr_vec_old = Hr_vec
         !
         unit = free_unit()
         open(unit,file="./Hr_vec.DAT",form="formatted",status="unknown",position="rewind",action="read")
         do iorb=1,size(Hr_vec)
            read(unit,"(1E20.12)") Hr_vec(iorb)
         enddo
         close(unit)
         !
         write(*,*)"updated from file Hr_vec old -> updated"
         do iorb=1,size(Hr_vec)
            write(*,"(2E20.12)") Hr_vec_old(iorb),Hr_vec(iorb)
         enddo
         !
         Hr_mat=czero
         call vec2mat(Hr_vec,map,maplen,Hr_mat)
         call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat,Crystal%Hk)
         !
      else
         !
         write(*,*)"starting Hr_vec"
         do iorb=1,size(Hr_vec)
            write(*,"(1E20.12)") Hr_vec(iorb)
         enddo
         !
      endif
      !
      !
      do ik=1,Crystal%Nkpt_path
         call check_Hermiticity(Crystal%Hk(:,:,ik),eps,enforce=.true.,hardstop=.true.,name="Hk_k"//str(ik),verb=.true.)
      enddo
      call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="inital",correction=dumCorr,doplane=FermiSurf)
      !
      !
      if(.not.justprint)then
         call fit_wrapper(chi2_Dispersion,Hr_vec,chi2,Niter,ftol=1d-15,estm=est)
         write(*,"(A,I,A,F)")"     Iterations: ",Niter," Chi^2: ",chi2
      endif
      !
      !
      write(*,*)"final Hr_vec"
      unit = free_unit()
      open(unit,file="./Hr_vec.DAT",form="formatted",status="unknown",position="rewind",action="write")
      do iorb=1,size(Hr_vec)
         write(*,"(1E20.12)") Hr_vec(iorb)
         write(unit,"(1E20.12)") Hr_vec(iorb)
      enddo
      close(unit)
      !
      !
      call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="final",correction=dumCorr,doplane=FermiSurf)
      Hr_mat=czero
      call vec2mat(Hr_vec,map,maplen,Hr_mat)
      call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat,Crystal%Hk)
      deallocate(Hr_vec,map,maplen)
      !
      !
      unit = free_unit()
      open(unit,file=reg(pathINPUT)//"Hr_fitted.DAT",form="formatted",status="unknown",position="rewind",action="write")
      write(unit,*)                    !skip first line
      write(unit,"(1I5)") Crystal%Norb !Number of Wannier orbitals
      write(unit,"(1I5)") Nwig         !Number of Wigner-Seitz vectors
      Qst = int(Nwig/W90NumCol)
      Rst = mod(Nwig,W90NumCol)
      do i=1,Qst
         write(unit,"("//str(W90NumCol)//"I5)")(nrdegwig(j+(i-1)*W90NumCol),j=1,W90NumCol)
      enddo
      if(Rst.ne.0)write(unit,"("//str(Rst)//"I5)")(nrdegwig(j+Qst*W90NumCol),j=1,Rst)
      Nwig_=0
      do ir=1,Nwig
        !do io=1,Crystal%Norb
        !    do jo=1,Crystal%Norb
        !       if(abs(Hr_mat(io,jo,iR)).gt.0d0)write(unit,"(5I5,2F12.6)") Nvecwig(:,iR), io, jo, dreal(Hr_mat(io,jo,iR)), dimag(Hr_mat(io,jo,iR))
        !    enddo
        !enddo
        if(Nvecwig(3,iR).eq.0)then
           write(unit,"(5I5,2E20.12)") Nvecwig(:,iR), 1, 1, dreal(Hr_mat(1,1,iR)), dimag(Hr_mat(1,1,iR))
           Nwig_=Nwig_+1
        else
           if(abs(Hr_mat(1,2,iR)).gt.0d0)write(*,"(3I3,2E20.12)") Nvecwig(:,iR),dreal(Hr_mat(1,2,iR)),radiuswig(iR)
        endif
        !
      enddo
      close(unit)
      deallocate(Hr_mat)
   endif
   write(*,*)"Nwig_ ",Nwig_
   !
   !
   call AllocateFermionicField(Glat,Crystal%Norb,Nmats,Nkpt=Crystal%Nkpt,Nsite=Nsite,Beta=Beta)
   call set_density(Glat%mu,Beta,Crystal,look4dens)
   call calc_Gmats(Glat,Crystal)
   call calc_density(Glat,Glat%N_s)
   densityLDA = Glat%N_s
   !
   call dump_Matrix(densityLDA,reg(pathINPUT),"Nlda",paramagnet)
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
      integer                               :: Nelements,io,jo,iR,ndx,ndx_shift
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
      call get_pattern(map_,Hr_mat_flattened,eps,listDim=maplen_,IncludeSingle=.true.)
      !
      if(allocated(Hr_vec_))deallocate(Hr_vec_)
      allocate(Hr_vec_(size(maplen_)));Hr_vec_=0d0
      ndx_shift=0
      do io=1,size(maplen_)
         Hr_vec_(io) = Hr_mat_flattened(map_(io,1))
      enddo
      Hr_vec_ = [Hr_vec_(2:size(maplen_))]
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
      real(8),allocatable                   :: Hr_vec_used(:)
      real(8),allocatable                   :: Hr_mat_flattened(:)
      integer                               :: Nelements,io,jo,iR,ndx,ndx_shift
      !
      if(.not.allocated(map_)) stop ".not.allocated(map_)"
      if((size(Hr_vec_)+1).ne.size(map_,dim=1)) stop "size(Hr_vec_).ne.size(map_,dim=1)"
      if((size(Hr_vec_)+1).ne.size(maplen_)) stop "size(Hr_vec_).ne.size(maplen_)"
      !
      Nelements = Crystal%Norb*Crystal%Norb*Nwig
      allocate(Hr_mat_flattened(Nelements));Hr_mat_flattened=0d0
      !
      Hr_vec_used = [ 0d0, Hr_vec_ ]
      !
      do io=1,size(maplen_)
         do jo=1,maplen_(io)
            Hr_mat_flattened(map_(io,jo)) = Hr_vec_used(io)
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
      deallocate(Hr_mat_flattened,Hr_vec_used)
      !
   end subroutine vec2mat
   !
   !
   subroutine chi2_Dispersion(Npara,Hr_vec_,chi2)
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: Hr_vec_
      real(8)                               :: chi2
      integer                               :: io,ik,unit_
      complex(8),allocatable                :: Hr_mat_(:,:,:)
      !
      ! Recompute dispersion with new Hr here
      write(*,*)"Internal iteration"
      unit_ = free_unit()
      open(unit_,file="./Hr_vec.DAT",form="formatted",status="unknown",position="rewind",action="write")
      do io=1,size(Hr_vec_)
         write(*,"(1E20.12)") Hr_vec_(io)
         write(unit_,"(1E20.12)") Hr_vec_(io)
      enddo
      close(unit_)
      !
      call vec2mat(Hr_vec_,map,maplen,Hr_mat_)
      call wannier_R2K(Crystal%Nkpt3,Crystal%kpt,Hr_mat_,Crystal%Hk)
      call interpolateHk2Path(Crystal,reg(structure),Nkpt_path,pathOUTPUT=reg(pathINPUT),store=.true.,corrname="optimized_fitting",correction=dumCorr)
      !
      chi2=0d0
      do io=1,Crystal%Norb
         do ik=1,Crystal%Nkpt_path
            chi2 = chi2 + abs(Crystal%Ek_path(io,ik)-Ek_dft(io,ik))**2
         enddo
      enddo
      !
   end subroutine chi2_Dispersion
   !
   !
   subroutine calc_dispersion(Ek_model_,Hr_vec_,name,printEk)
      implicit none
      real(8),allocatable,intent(inout)     :: Ek_model_(:,:)
      real(8),intent(in)                    :: Hr_vec_(:)
      character(len=*),intent(in)           :: name
      logical,intent(in)                    :: printEk
      real(8),allocatable                   :: thopping_(:,:)
      integer                               :: iorb_,ik_,idist_,unit_
      !
      !recreate the original thopping_ matrix from the vector
      call vec2vec(Hr_vec_,thopping_)
      !
      ! print to file
      if(verbs)write(*,"(A)")" call to calc_dispersion at stage: "//reg(name)
      unit_ = free_unit()
      open(unit_,file="./hopping_"//reg(name)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
      write(unit_,*) maxdist, Crystal%Norb
      do idist_=1,maxdist
         if(verbs)write(*,*) idist_,(thopping_(idist_,iorb_),iorb_=1,Crystal%Norb)
         write(unit_,*) idist_,(thopping_(idist_,iorb_),iorb_=1,Crystal%Norb)
      enddo
      close(unit_)
      !
      !compute dispersion
      if(allocated(Ek_model_))deallocate(Ek_model_)
      allocate(Ek_model_(Crystal%Norb,Crystal%Nkpt_path));Ek_model_=0d0
      do ik_=1,Crystal%Nkpt_path
         do iorb_=1,Crystal%Norb
            do idist_=1,maxdist
               Ek_model_(iorb_,ik_) = Ek_model_(iorb_,ik_) + thopping_(idist_,iorb_) * cosinesEk(idist_,ik_) !<=====
            enddo
         enddo
      enddo
      !
      deallocate(thopping_)
      !
      if(printEk)then
         unit_ = free_unit()
         open(unit_,file=reg(pathINPUT)//"Bands_model_"//reg(name)//".DAT",form="formatted",status="unknown",position="rewind",action="write")
         do ik_=1,Crystal%Nkpt_path
            write(unit_,"(1I5,10000E20.12)") ik_,Crystal%Kpathaxis(ik_)/Crystal%Kpathaxis(Crystal%Nkpt_path),(Ek_model_(iorb_,ik_),iorb_=1,Crystal%Norb)
         enddo
         close(unit_)
      endif
      !
   end subroutine calc_dispersion
   !
   !
   subroutine vec2vec(Hr_vec_,thopping_)
      implicit none
      real(8),intent(in)                    :: Hr_vec_(:)
      real(8),allocatable,intent(inout)     :: thopping_(:,:)
      integer                               :: ndx_,iorb_,idist_
      !
      if(allocated(thopping_))deallocate(thopping_)
      allocate(thopping_(maxdist, Crystal%Norb));thopping_=0d0
      ndx_=1
      do iorb_=1,Crystal%Norb
         do idist_=1,size(mapOrb(iorb_)%maplen)
            do ihop=1,mapOrb(iorb_)%maplen(idist_)
               thopping_(mapOrb(iorb_)%map(idist_,ihop),iorb_) = Hr_vec_(ndx_)
            enddo
            ndx_ = ndx_ + 1
         enddo
      enddo
      !
   end subroutine vec2vec
   !
   !
   subroutine chi2_DispersionModel(Npara,Hr_vec_,chi2)
      implicit none
      integer,intent(in)                    :: Npara
      real(8),dimension(Npara),intent(in)   :: Hr_vec_
      real(8)                               :: chi2
      integer                               :: iorb_,ik_
      real(8),allocatable                   :: Ek_model_(:,:)
      !
      ! Recompute dispersion with new Hr here
      call calc_dispersion(Ek_model_,Hr_vec_,"fitting",.false.)
      !
      chi2=0d0
      do iorb_=1,Crystal%Norb
         do ik_=1,Crystal%Nkpt_path
            chi2 = chi2 + abs(Ek_model_(iorb_,ik_)-Ek_dft(iorb_,ik_))**2
         enddo
      enddo
      !
   end subroutine chi2_DispersionModel
   !
   !
   function fcosinesEk(idist_,ik_) result(cosEk_id)
     integer,intent(in)           :: idist_,ik_
     real(8)                      :: cosEk(6),Kvec_(3)
     real(8)                      :: csi_,neta_,cosEk_id
     real(8)                      :: kx_,ky_,kz_
     !
     Kvec_ =  Crystal%kptpath(1,ik_) * Blat(:,1) + Crystal%kptpath(2,ik_) * Blat(:,2) + Crystal%kptpath(3,ik_) * Blat(:,3)
     kx_ = Kvec_(1)
     ky_ = Kvec_(2)
     kz_ = Kvec_(3)
     !
     csi_ = kx_ !0.5d0 * Crystal%kptpath(1,ik_) *2*pi
     neta_ = ky_ !(sqrt(3d0)/2d0) * Crystal%kptpath(2,ik_) *2*pi
     !
     cosEk(1) = 1d0
     cosEk(2) = ( 2*cos(csi_)*cos(neta_)     + cos(2*csi_) )
     cosEk(3) = ( 2*cos(3*csi_)*cos(neta_)   + cos(2*neta_) )
     cosEk(4) = ( 2*cos(2*csi_)*cos(2*neta_) + cos(4*csi_) )
     cosEk(5) = ( cos(csi_)*cos(3*neta_) + cos(5*csi_)*cos(neta_) + cos(4*csi_)*cos(2*neta_) )
     cosEk(6) = ( cos(3*csi_)*cos(3*neta_) + cos(6*csi_) )
     !
     cosEk_id = cosEk(idist_)
     !
  end function fcosinesEk
  function fcosinesHk(idist_,ik_) result(cosEk_id)
    integer,intent(in)           :: idist_,ik_
    real(8)                      :: cosEk(6),Kvec_(3)
    real(8)                      :: csi_,neta_,cosEk_id
    real(8)                      :: kx_,ky_,kz_
    !
    Kvec_ =  Crystal%kpt(1,ik_) * Blat(:,1) + Crystal%kpt(2,ik_) * Blat(:,2) + Crystal%kpt(3,ik_) * Blat(:,3)
    kx_ = Kvec_(1)
    ky_ = Kvec_(2)
    kz_ = Kvec_(3)
    !
    csi_ = kx_ !0.5d0 * Crystal%kptpath(1,ik_) *2*pi
    neta_ = ky_ !(sqrt(3d0)/2d0) * Crystal%kptpath(2,ik_) *2*pi
    !
    cosEk(1) = 1d0
    cosEk(2) = ( 2*cos(csi_)*cos(neta_)     + cos(2*csi_) )
    cosEk(3) = ( 2*cos(3*csi_)*cos(neta_)   + cos(2*neta_) )
    cosEk(4) = ( 2*cos(2*csi_)*cos(2*neta_) + cos(4*csi_) )
    cosEk(5) = ( cos(csi_)*cos(3*neta_) + cos(5*csi_)*cos(neta_) + cos(4*csi_)*cos(2*neta_) )
    cosEk(6) = ( cos(3*csi_)*cos(3*neta_) + cos(6*csi_) )
    !
    cosEk_id = cosEk(idist_)
    !
end function fcosinesHk
   !
   !
   !
end program pp_dispersion_fit
