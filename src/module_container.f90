module module_container


   ! COMMENTS:
   !
   !


   !(a) Standalone module. Can be used as generic library.
   use linalg, only:                         &
   eig, eigh                                ,&
   inv, inv_sym, inv_her                    ,&
   det,trace                                ,&
   deye, zeye, zeros, ones                  ,&
   diag, diagonal, rotate


   !(b) Standalone module. Can be used as generic library.
   use utils_misc, only:                     &
   find_kpt, keq                            ,&
   FermionicFreqMesh,BosonicFreqMesh        ,&
   fermidirac,diff_fermidirac               ,&
   linspace,denspace                        ,&
   FermionicFilon,BosonicFilon              ,&
   reg,str                                  ,&
   free_unit                                ,&
   tick,tock                                ,&
   inquireFile,inquireDir,createDir         ,&
   check_Hermiticity                        ,&
   assert_shape


   !(c) Lattice related quantities. Depends on (a),(b) not on specific types.
   use crystal, only:                        &
   read_lattice                             ,&
   read_xeps                                ,&
   read_Hk                                  ,&
   fill_ksumkdiff                           ,&
   fill_smallk                              ,&
   wannierinterpolation                     ,&
   wannier_K2R,wannier_R2K,wannier_K2R_NN   !,&
   !add_crystalfields


   !(d) Fermionic and Bosonic Fourier transforms. Depends on (b),(c) not on specific types.
   use fourier_transforms, only:             &
   Fmats2itau_mat,Fmats2itau_vec            ,&
   Fitau2mats_mat,Fitau2mats_vec            ,&
   Bmats2itau,Bitau2mats


   !(d) Contains units and definitions of lattice, fermionic and bosnoic types.
   use parameters


   !(e) Standalone module. Contains input variables.
   use global_vars


   !(f) Container attributes manipulations. Depends on parameters
   use utils_fields, only:                   &
   FermionicKsum                            ,&
   BosonicKsum                              ,&
   selfAllocateFermionicField               ,&
   selfDeallocateFermionicField             ,&
   selfAllocateBosonicField                 ,&
   selfDeallocateBosonicField               ,&
   clear_attributes


   !(g) Input/Output routines. Depends on (utils_misc),(parameters),(utils_fields)
   use file_io, only:                        &
   dump_Matrix                              ,&
   dump_FermionicField                      ,&
   read_FermionicField                      ,&
   dump_BosonicField                        ,&
   read_BosonicField


   !(h) Interactions container. Depends on (utils_misc),(crystal),(parameters),(global_vars),(utils_fields),(file_io)
   use interactions, only:                   &
   read_spex                                ,&
   calc_W_full                              ,&
   calc_W_edmft                             ,&
   calc_chi_full                            ,&
   calc_chi_edmft !build_Umatrix,rescale_interaction


   !(i) Bubble diagram container. Depends on (utils_misc),(crystal),(parameters),(global_vars),(utils_fields),(fourier_transforms)
   use bubbles, only:                        &
   calc_Pi !calc_Optcond,calc_Hall


   !(i) Self-energy container. Depends on (linalg),(utils_misc),(crystal),(parameters),(global_vars),(utils_fields),(fourier_transforms)
   use self_energy, only:                    &
   calc_sigmaGW !calc_Optcond,calc_Hall





end module module_container
