subroutine interpolate2Beta_Fermionic(G,Beta_Match,mode,offDiag,wmats_in)
   !
   use parameters
   use utils_misc
   use utils_fields
   use interactions
   use input_vars, only: LocalOrbs
   implicit none
   !
   type(FermionicField),intent(inout)    :: G
   type(OldBeta),intent(in)              :: Beta_Match
   character(len=*),intent(in)           :: mode
   logical,intent(in)                    :: offDiag
   real(8),intent(in),optional           :: wmats_in(:)
   !
   type(FermionicField)                  :: G_old
   integer                               :: Nsite,isite,ik,iw
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
   if(.not.allocated(LocalOrbs)) stop "interpolate2Beta_Fermionic: LocalOrbs not properly initialized."
   !
   Nsite = size(LocalOrbs)
   !
   LocalOnly=.true.
   if(G%Nkpt.ne.0)LocalOnly=.false.
   !
   allocate(wmats_new(Beta_Match%Nmats_new)); wmats_new=FermionicFreqMesh(Beta_Match%Beta_new,Beta_Match%Nmats_new)
   if(present(wmats_in))then
      wmats_old = wmats_in
   else
      allocate(wmats_old(Beta_Match%Nmats_old)); wmats_old=FermionicFreqMesh(Beta_Match%Beta_old,Beta_Match%Nmats_old)
   endif
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
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(ispin,isite,i,j,i_lat,j_lat,replace,iw,ReD,ReGf,ImD,ImGf)
         !$OMP DO
         do isite=1,Nsite
            !
            do i=1,LocalOrbs(isite)%Norb
               do j=1,LocalOrbs(isite)%Norb
                  !
                  i_lat = LocalOrbs(isite)%Orbs(i)
                  j_lat = LocalOrbs(isite)%Orbs(j)
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

subroutine interpolate2Beta_Bosonic(W,Beta_Match,mode,offDiag,wmats_in)
   !
   use parameters
   use utils_misc
   use utils_fields
   use interactions
   use input_vars, only: LocalOrbs, UfullStructure
   implicit none
   !
   type(BosonicField),intent(inout)      :: W
   type(OldBeta),intent(in)              :: Beta_Match
   character(len=*),intent(in)           :: mode
   logical,intent(in)                    :: offDiag
   real(8),intent(in),optional           :: wmats_in(:)
   !
   type(BosonicField)                    :: W_old
   integer                               :: Norb
   integer                               :: Nsite,isite,iq,ib1,ib2,iw
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
   if(.not.allocated(LocalOrbs)) stop "interpolate2Beta_Bosonic: LocalOrbs not properly initialized."
   !
   Nsite = size(LocalOrbs)
   Norb = int(sqrt(dble(W%Nbp)))
   !
   LocalOnly=.true.
   if(W%Nkpt.ne.0)LocalOnly=.false.
   !
   allocate(wmats_new(Beta_Match%Nmats_new)); wmats_new=BosonicFreqMesh(Beta_Match%Beta_new,Beta_Match%Nmats_new)
   if(present(wmats_in))then
      wmats_old = wmats_in
   else
      allocate(wmats_old(Beta_Match%Nmats_old)); wmats_old=BosonicFreqMesh(Beta_Match%Beta_old,Beta_Match%Nmats_old)
   endif
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
         !$OMP PARALLEL DEFAULT(SHARED),&
         !$OMP PRIVATE(isite,i,j,k,l,i_lat,j_lat,k_lat,l_lat,ib1,ib2,replace,iw,ReD,ReW)
         !$OMP DO
         do isite=1,W%Nsite
            !
            do i=1,LocalOrbs(isite)%Norb
               do j=1,LocalOrbs(isite)%Norb
                  do k=1,LocalOrbs(isite)%Norb
                     do l=1,LocalOrbs(isite)%Norb
                        !
                        ! mapping
                        i_lat = LocalOrbs(isite)%Orbs(i)
                        j_lat = LocalOrbs(isite)%Orbs(j)
                        k_lat = LocalOrbs(isite)%Orbs(k)
                        l_lat = LocalOrbs(isite)%Orbs(l)
                        !
                        ! bosonic indexes on the lattice
                        call F2Bindex(Norb,[i_lat,j_lat],[k_lat,l_lat],ib1,ib2)
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
