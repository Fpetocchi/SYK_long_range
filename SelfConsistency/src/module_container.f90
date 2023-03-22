module module_container

   !(a) Standalone module. Can be used as generic library.
   use linalg, only:                         &
   eig, eigh                                ,&
   inv, inv_sym, inv_her                    ,&
   det, det3, trace, Cramer                 ,&
   deye, zeye, zeros, ones, dag             ,&
   diag, diagonal, rotate, diag_factor      ,&
   kronecker_product, outerprod             ,&
   cross_product, s3_product, RotSO3        ,&
   tensor_transform


   !(b) Standalone module. Can be used as generic library.
   use utils_misc, only:                     &
   F2Bindex                                 ,&
   find_vec, keq                            ,&
   init_Uelements                           ,&
   FermionicFreqMesh,BosonicFreqMesh        ,&
   fermidirac,diff_fermidirac               ,&
   linspace,denspace                        ,&
   FermionicFilon,BosonicFilon              ,&
   reg,str,free_unit                        ,&
   tick,tock                                ,&
   inquireFile,inquireDir                   ,&
   createDir,skip_header                    ,&
   check_Hermiticity,check_Symmetry         ,&
   KK_Im2Re,KK_Re2Im                        ,&
   halfbeta_sym                             ,&
   get_moments_F                            ,&
   assert_shape                             ,&
   splint,nspline                           ,&
   cubic_interp                             ,&
   linear_interp_2y,linear_interp_2x        ,&
   trapezoid_integration                    ,&
   get_pattern,sort_array,flip_array


   !(c) Lattice related quantities. Depends on (a),(b) not on specific types.
   use crystal, only:                        &
   read_lattice                             ,&
   set_lattice                              ,&
   read_xeps                                ,&
   read_Hk                                  ,&
   dump_Hk                                  ,&
   build_kpt                                ,&
   build_Hk                                 ,&
   build_Hk_from_Hr                         ,&
   fill_ksumkdiff                           ,&
   fill_smallk                              ,&
   wannierinterpolation                     ,&
   wannier_K2R                              ,&
   wannier_R2K                              ,&
   wannier_K2R_NN                           ,&
   calc_Kpath                               ,&
   calc_Kplane                              ,&
   calc_Ewald                               ,&
   set_UserPath                             ,&
   get_Ruc,get_Rlat,get_Blat                ,&
   interpolate2Path                         ,&
   interpolate2Plane                        ,&
   tetrahedron_integration


   !(d) Fermionic and Bosonic Fourier transforms. Depends on (a),(b),(c) not on specific types.
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
   FieldKsum                                ,&
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
   TransformMatrix                          ,&
   clear_attributes                         ,&
   clear_MatrixElements                     ,&
   isReal                                   ,&
   duplicate                                ,&
   loc2imp,imp2loc,Expand2Nsite             ,&
   symmetrize_GW                            ,&
   symmetrize_imp                           ,&
   MergeFields                              ,&
   join_SigmaCX                             ,&
   product2NN,NN2product                    ,&
   calc_Ek,calc_Ep                          ,&
   check_QP_poles                           ,&
   build_Potential


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
   fit_wrapper                              ,&
   fit_moments                              ,&
   fit_delta                                ,&
   G_Moments                                ,&
   S_Moments                                ,&
   W_Moments


   !(l) Interactions container. Depends on (b),(c),(d),(e),(f),(g),(h)
   use interactions, only:                   &
   read_U_spex                              ,&
   read_U_vasp                              ,&
   calc_W_full                              ,&
   calc_W_edmft                             ,&
   calc_chi                                 ,&
   build_Umat                               ,&
   build_Uret                               ,&
   calc_curlyU                              ,&
   calc_Wimp                                ,&
   calc_QMCinteractions


   !(m) Bubble diagram container. Depends on (a),(b),(c),(d),(e),(f),(g),(h)
   use bubbles, only:                        &
   calc_PiGG                                ,&
   calc_PiGGdc                              ,&
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
   calc_VH                                  ,&
   print_G0W0_dispersion,Uwan_stored


   !(p) post-processing container. Depends on (a),(b),(c),(d),(e),(f),(g),(h)
   use post_processing, only:                &
   dump_MaxEnt                              ,&
   pade                                     ,&
   remove_CDW                               ,&
   interpolate2Beta                         ,&
   interpolate2kpath                        ,&
   calc_Tc
   !calc_OptCond,calc_HallCond


   !------------------------ Modules Obscured to User -------------------------!
   !gap_equation container. Depends only on (a),(b),(c),(e).
   !Interface with main inside post_processing module.


end module module_container
