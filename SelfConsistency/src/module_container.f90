module module_container


   ! COMMENTS:
   !
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
   wannier_K2R                              ,&
   wannier_R2K                              ,&
   wannier_K2R_NN!                          ,&
   !add_crystalfields


   !(d) Fermionic and Bosonic Fourier transforms. Depends on (b),(c) not on specific types.
   use fourier_transforms, only:             &
   Fmats2itau_mat,Fmats2itau_vec            ,&
   Fitau2mats_mat,Fitau2mats_vec            ,&
   Bmats2itau,Bitau2mats


   !(e) Contains units and definitions of lattice, fermionic and bosnoic types.
   use parameters


   !(f) Standalone module. Contains input variables.
   use input_vars!                           &
   !read_inputfile


   !(g) Container attributes manipulations. Depends on (e)
   use utils_fields, only:                   &
   FermionicKsum                            ,&
   BosonicKsum                              ,&
   AllocateFermionicField                   ,&
   DeallocateFermionicField                 ,&
   AllocateLattice                          ,&
   DeallocateLattice                        ,&
   AllocateBosonicField                     ,&
   DeallocateBosonicField                   ,&
   clear_attributes!                        ,&
   !reshuffle_imp2loc                       ,&
   !reshuffle_loc2imp


   !(h) Input/Output routines. Depends on (b),(e),(g)
   use file_io, only:                        &
   dump_Matrix                              ,&
   dump_FermionicField                      ,&
   read_FermionicField                      ,&
   dump_BosonicField                        ,&
   read_BosonicField


   !(i) Interactions container. Depends on (b),(c),(e),(f),(g),(h)
   use interactions, only:                   &
   read_U_spex                              ,&
   calc_W_full                              ,&
   calc_W_edmft                             ,&
   calc_chi_full                            ,&
   calc_chi_edmft !                         ,&
   !build_Umatrix                           ,&
   !rescale_interaction


   !(j) Bubble diagram container. Depends on (b),(c),(d),(e),(f),(g)
   use bubbles, only:                        &
   calc_Pi!                                 ,&
   !calc_Optcond                            ,&
   !calc_Hall


   !(k) greens_function container. Depends on (a),(b),(c),(d),(e),(f),(g),(h)
   use greens_function, only:                &
   calc_density                             ,&
   set_density                              ,&
   calc_Gmats                               ,&
   calc_Glda


   !(l) Self-energy container. Depends on (a),(b),(c),(d),(e),(f),(g),(h),(k)
   use self_energy, only:                    &
   calc_sigmaGW                             ,&
   calc_sigmaGW_DC                          ,&
   read_Sigma_spex


end module module_container
