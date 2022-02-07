subroutine remove_CDW(W,mode,site)
   !
   use parameters
   use utils_misc
   use utils_fields
   use interactions
   use linalg
   use input_vars, only: LocalOrbs
   implicit none
   !
   type(BosonicField),intent(inout)      :: W
   character(len=*),intent(in)           :: mode
   integer,intent(in),optional           :: site
   !
   integer                               :: Norb,Nmats,iw
   integer                               :: Nsite,isite,ib1,ib2
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
   if(.not.allocated(LocalOrbs)) stop "remove_CDW: LocalOrbs not properly initialized."
   !
   Norb = int(sqrt(dble(W%Nbp)))
   Nmats = W%Npoints
   Nsite = size(LocalOrbs)
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
                        if(any([i_lat,j_lat,k_lat,l_lat].gt.Norb)) stop "remove_CDW: the input field is not in the lattice space."
                        !
                        ! bosonic indexes on the lattice
                        call F2Bindex(Norb,[i_lat,j_lat],[k_lat,l_lat],ib1,ib2)
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
