module module_container

   !(a) Standalone module. Can be used as generic library.
   use linalg, only:                         &
   eig, eigh                                ,&
   inv, inv_sym, inv_her                    ,&
   det,trace                                ,&
   deye, zeye, zeros, ones                  ,&
   diag, diagonal, rotate, diag_factor      ,&
   kronecker_product, outerprod             ,&
   cross_product, s3_product                ,&
   tensor_transform


   !(b) Standalone module. Can be used as generic library.
   use utils_misc, only:                     &
   find_kpt, keq                            ,&
   FermionicFreqMesh,BosonicFreqMesh        ,&
   fermidirac,diff_fermidirac               ,&
   linspace,denspace                        ,&
   FermionicFilon,BosonicFilon              ,&
   get_Tier_occupation                      ,&
   reg,str                                  ,&
   free_unit                                ,&
   tick,tock                                ,&
   inquireFile,inquireDir,createDir         ,&
   check_Hermiticity                        ,&
   check_Symmetry                           ,&
   KK_Im2Re,KK_Re2Im                        ,&
   halfbeta_symm                            ,&
   halfbeta_antisymm                        ,&
   assert_shape                             ,&
   splint,nspline,cubic_interp              ,&
   get_pattern,sort_array


   !(c) Lattice related quantities. Depends on (a),(b) not on specific types.
   use crystal, only:                        &
   read_lattice                             ,&
   read_xeps                                ,&
   read_Hk                                  ,&
   build_Hk                                 ,&
   fill_ksumkdiff                           ,&
   fill_smallk                              ,&
   wannierinterpolation                     ,&
   wannier_K2R                              ,&
   wannier_R2K                              ,&
   wannier_K2R_NN                           ,&
   calc_path                                ,&
   interpolateHk2Path


   !(d) Fermionic and Bosonic Fourier transforms. Depends on (b),(c) not on specific types.
   use fourier_transforms, only:             &
   Fmats2itau_mat,Fmats2itau_vec,Fmats2itau ,&
   Fitau2mats_mat,Fitau2mats_vec,Fitau2mats ,&
   Bmats2itau,Bitau2mats


   !(e) Contains units and definitions of lattice, fermionic and bosnoic types.
   use parameters


   !(f) Variables and read inputifile. Depends on (b) and (e)
   use input_vars


   !(g) Container attributes manipulations. Depends on (a),(b) and (e)
   use utils_fields, only:                   &
   FermionicKsum                            ,&
   BosonicKsum                              ,&
   AllocateLattice                          ,&
   DeallocateLattice                        ,&
   AllocateFermionicField                   ,&
   AllocateBosonicField                     ,&
   DeallocateFermionicField                 ,&
   DeallocateBosonicField                   ,&
   DeallocateField                          ,&
   TransformBosonicField                    ,&
   clear_attributes                         ,&
   clear_MatrixElements                     ,&
   isReal                                   ,&
   duplicate                                ,&
   loc2imp,imp2loc                          ,&
   symmetrize_GW                            ,&
   symmetrize_imp                           ,&
   MergeFields                              ,&
   join_SigmaCX


   !(h) Input/Output routines. Depends on (b),(e),(g)
   use file_io, only:                        &
   dump_Matrix                              ,&
   read_Matrix                              ,&
   dump_Field_component                     ,&
   dump_FermionicField                      ,&
   read_FermionicField                      ,&
   dump_BosonicField                        ,&
   read_BosonicField


   !(i) minimization routines. Depends on (b),(e),(h)
   use fit, only:                            &
   fit_moments                              ,&
   fit_delta                                ,&
   G_Moments                                ,&
   S_Moments                                ,&
   W_Moments


   !(l) Interactions container. Depends on (b),(c),(e),(f),(g),(h)
   use interactions, only:                   &
   init_Uelements                           ,&
   read_U_spex                              ,&
   calc_W_full                              ,&
   calc_W_edmft                             ,&
   calc_chi_full                            ,&
   calc_chi_edmft                           ,&
   build_Umat                               ,&
   build_Uret                               ,&
   calc_curlyU                              ,&
   calc_Wimp                                ,&
   calc_QMCinteractions


   !(m) Bubble diagram container. Depends on (a),(b),(c),(d),(e),(f),(g),(h),(l)
   use bubbles, only:                        &
   calc_Pi                                  ,&
   calc_Pimp


   !(n) greens_function container. Depends on (a),(b),(c),(d),(e),(f),(g),(h)
   use greens_function, only:                &
   calc_density                             ,&
   set_density                              ,&
   calc_G0_tau                              ,&
   calc_Gmats                               ,&
   calc_Glda


   !(o) Self-energy container. Depends on (a),(b),(c),(d),(e),(f),(g),(h),(n)
   use self_energy, only:                    &
   calc_sigmaGW                             ,&
   calc_sigmaGWdc                           ,&
   read_Sigma_spex                          ,&
   calc_VH


   !(p) post-processing container. Depends on (a),(b),(c),(d),(e),(f),(g),(h)
   use post_processing, only:                &
   dump_MaxEnt                              ,&
   pade                                     ,&
   remove_CDW                               ,&
   interpolate2Beta                         ,&
   interpolateG2Path
   !calc_OptCond,calc_HallCond


end module module_container
