// -----------------------------------------------------------------------------

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <string>

// -----------------------------------------------------------------------------

#ifndef CNTR_USE_OMP
#define CNTR_USE_OMP
#endif

#ifndef CNTR_USE_MPI
#define CNTR_USE_MPI
#endif

#define NCA_SOLVER_ASSERT_0 0
#define NCA_SOLVER_ASSERT_1 0

#include <cntr/cntr.hpp>
#include <cntr/utils/read_inputfile.hpp>

// -----------------------------------------------------------------------------

#include "./ppsc/ppsc.hpp"
#include "./ppsc/solver.hpp"

#include "./ppsc/hilbert_spaces/three_band_fermi_diag.hpp"
#include "./ppsc/hilbert_spaces/three_band_fermi_densdens.hpp"
#include "./ppsc/hamiltonians/three_band_hubbard_solar.hpp"

#include "./ppsc/baths/non_int_boson_propagator.hpp"

// -----------------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------------

typedef ppsc::operator_type operator_type;
typedef ppsc::mam::dynamic_matrix_type matrix_type;


////////////////////////////////////////////////////////////////////////////////


void integrate_time(std::vector<double> &int_ft, std::vector<double> &ft,
   double h,int kt)
{
   // int_ft[t+1] = int_0^{t} ds ft(s), where ft(t)=ft[t+1], int_ft[0]=0
   int nt=ft.size()-2,tstp;
   int kt1=(nt==-1 || nt>=kt ? kt : nt),n,n1;
   double tmp;
   if(nt<-1) throw("integrate_time: nt too small");
   int_ft.assign(nt+2,0);
   for(tstp=1;tstp<=nt;tstp++)
   {
      tmp=0.0;
      n1=(tstp<kt1 ? kt1 : tstp);
      for(n=0;n<=n1;n++)tmp+=integration::I<double>(kt1).gregory_weights(tstp,n)*ft[n+1];
      int_ft[tstp+1]=h*tmp;
   }
};


////////////////////////////////////////////////////////////////////////////////


std::vector<std::complex<double>> get_S_star_Gles_timestep(int tstp,
                                                  ppsc::gf_tstps_type & Spot,
                                                  ppsc::gf_tstps_type & Gloc,
                                                  double beta, double h, int kt)
{
   //
   std::vector<std::complex<double>> Fvec;
   //
   for(int i=0; i < Gloc.size(); i++)
   {
      std::complex<double> Fval;
      cntr::herm_matrix<double> G(tstp, Gloc[0].ntau(), 1, -1);
      cntr::herm_matrix<double> S(tstp, Gloc[0].ntau(), 1, -1);
      cntr::herm_matrix<double> F(tstp, Gloc[0].ntau(), 1, -1);
      G.set_timestep(tstp,Gloc[i]);
      S.set_timestep(tstp,Spot[i]);
      ::cntr::convolution_timestep_les<double,cntr::herm_matrix<double>,1>
      (tstp,F,G,G,S,S,integration::I<double>(kt),beta,h);
      //G.print_to_file("Gpot",16);
      //S.print_to_file("Spot",16);
      //F.print_to_file("Fpot",16);
      F.get_les(tstp,tstp,Fval);
      Fvec.push_back(Fval);
   }
   return Fvec;
};


////////////////////////////////////////////////////////////////////////////////


class solver_type
{
public:
   solver_type(int nt, int ntau, double beta, double h, int kt, int nomp
      ,ppsc::hilbert_spaces::three_band_fermi_diag & hilb, int order):
      nt_(nt),ntau_(ntau),kt_(kt),beta_(beta),h_(h),
      imp(nt, ntau, beta, h, kt, nomp, hilb, order)//,
      //nt(nt),ntau(ntau),beta(beta),h(h),kt(kt)
      {}
      int nt(void)      const {return nt_;}
      int ntau(void)    const {return ntau_;}
      int kt(void)      const {return kt_;}
      double h(void)    const {return h_;}
      double beta(void) const {return beta_;}
      solver_type & operator=( solver_type& Sin) { return Sin; }
      //
      //int nt, ntau, kt, order, nomp;
      //double beta, h;
      typedef ppsc::hilbert_spaces::three_band_fermi_diag hilbert_space_type;
      hilbert_space_type hilb;
      typedef ppsc::hamiltonians::three_band_hubbard_solar<hilbert_space_type,
      ppsc::hamiltonians::interaction_type::kanamori> hamiltonian_type;
      ppsc::solver<hamiltonian_type> imp;
      //
      int nt_, ntau_, kt_;
      double beta_, h_;
};


////////////////////////////////////////////////////////////////////////////////


#include "./LVOheaders/SpinUNResolvedFields.hpp"


////////////////////////////////////////////////////////////////////////////////



int main(int argc, char *argv[])
{
   int nt, ntau, kt, tstp, order, nomp, store_pp;
   int iter_equil, iter_warmup, itermax,iter_rt, iter_rtime;
   int read_eq_sym, read_rt_sym, read_state_from_file, sym_flag;
   int setup;
   int EndVRamp, equil, photon_ndx, j_flag;
   int af_bethe_flag, Ti_phflag, LVO_phflag, put_nonint_leads;

   double beta, h, sym_tol, errmax, linear_mixing;
   double HinplaneXY, HtoBaths_Z;
   double gradient, gradientBulk, Eloc_top, Eloc_bot, Vbias, AdjYTO;
   double lead_beta_bot, lead_W_bot, lead_beta_top, lead_W_top;

   double omega;
   std::vector<double> g;

   bool matsubara_converged = false;
   bool warmup_converged = false;

   int Nspin=1;
   int Norb=3;

   AdjYTO=0.0;

   try
   {
      // Read general input - no microscopic parameters-------------------------
      {
        if (argc < 2) throw("COMMAND LINE ARGUMENT MISSING");
        find_param(argv[1], "__nt="                  , nt                   );
        find_param(argv[1], "__ntau="                , ntau                 );
        find_param(argv[1], "__beta="                , beta                 );
        find_param(argv[1], "__h="                   , h                    );
        find_param(argv[1], "__itermax="             , itermax              );
        find_param(argv[1], "__errmax="              , errmax               );
        find_param(argv[1], "__iter_rtime="          , iter_rtime           );
        find_param(argv[1], "__kt="                  , kt                   );
        find_param(argv[1], "__order="               , order                );
        find_param(argv[1], "__store_pp="            , store_pp             );
        find_param(argv[1], "__linear_mixing="       , linear_mixing        );
        find_param(argv[1], "__sym_flag="            , sym_flag             );
        find_param(argv[1], "__sym_tol="             , sym_tol              );
        find_param(argv[1], "__read_eq_sym="         , read_eq_sym          );
        find_param(argv[1], "__read_rt_sym="         , read_rt_sym          );
        find_param(argv[1], "__read_state_from_file=", read_state_from_file );
        find_param(argv[1], "__nomp="                , nomp                 );
        find_param(argv[1], "__equil="               , equil                );
        //
        find_param(argv[1], "__setup="               , setup                );
        //
        find_param(argv[1], "__af_bethe_flag="       , af_bethe_flag        );
        //
        find_param(argv[1], "__put_nonint_leads="    , put_nonint_leads     );
        //
        find_param(argv[1], "__layer_phn_flag="      , LVO_phflag           );
        find_param(argv[1], "__lead_phn_flag="       , Ti_phflag            );
        find_param_tvector(argv[1], "__g="           , g, nt                );
        find_param(argv[1], "__omega="               , omega                );
        //
        find_param(argv[1], "__photon_ndx="          , photon_ndx           );
        find_param(argv[1], "__gradient="            , gradient             );
        find_param(argv[1], "__gradientBulk="        , gradientBulk         );
        find_param(argv[1], "__Vbias="               , Vbias                );
        find_param(argv[1], "__EndVRamp="            , EndVRamp             );
        find_param(argv[1], "__HinplaneXY="          , HinplaneXY           );
        find_param(argv[1], "__HtoBaths_Z="          , HtoBaths_Z           );
        find_param(argv[1], "__j_flag="              , j_flag               );
        if(setup==3)find_param(argv[1], "__AdjYTO=", AdjYTO);
      }

      //dimensione della Hamiltoniana
      int Nlat, Nlayer;
      // BULK SETUP
      if(setup==1)
      {
         Nlat=4;
         Nlayer=1;
      }
      // LVO SETUP
      if(setup==2)
      {
         Nlat=12;
         Nlayer=6;
      }
      // YTO SETUP
      if(setup==3)
      {
         Nlat=12;
         Nlayer=6;
      }

      // Setup potential bias---------------------------------------------------
      double Vbiasv[nt+2];
      double Vstep=Vbias/EndVRamp;
      for (int it=0; it<=EndVRamp; it++ )     Vbiasv[it]=Vstep*it;
      for (int it=EndVRamp+1; it<=nt+1; it++ )Vbiasv[it]=Vbias;
      FILE * VbiasFile; VbiasFile = fopen ("Vbias.out","w");
      for (int it=0; it <nt+2; it++)fprintf (VbiasFile,"%.12e\n",Vbiasv[it]);
      fclose(VbiasFile);

      // Setup boson bath-------------------------------------------------------
      cntr::herm_matrix<double> D0(nt, ntau, 1, +1);
      ppsc::boson_utils::green_from_eps_phonon(beta, D0, omega, h);
      cntr::function<double> gfunc(nt);
      for(int tstp=-1;tstp<=nt;tstp++)
      {
         cdmatrix tmp(1,1);
         tmp(0,0)=g[tstp+1];
         gfunc.set_value(tstp,tmp);
         //std::cout<<g[tstp+1]<<std::endl;
      }
      for(int tstp=-1; tstp <= nt; tstp++)
      {
         D0.left_multiply(tstp, gfunc);
         D0.right_multiply(tstp, gfunc);
      }
      D0.print_to_file("Bosonprop.out",16);

      // Setup bulk fake bias---------------------------------------------------
      cntr::function<double> psi(nt),psidag(nt);
      if(gradientBulk!=0.0)
      {
         double Estep=gradientBulk/20.0;
         std::vector<double> intdE(nt+2,0),dE(nt+2,0);
         for (int it=0; it<=20; it++ )dE[it]=Estep*it;
         for (int it=20+1; it<=nt+1; it++ )dE[it]=gradientBulk;
         integrate_time(intdE,dE,h,kt);
         for(int tstp=-1;tstp<=nt;tstp++)
         {
            double sde=sin(intdE[tstp+1]);
            double cde=cos(intdE[tstp+1]);
            psi[tstp]= std::complex<double>(cde,-sde);
            psidag[tstp]=conj(psi[tstp]);
         }
      }

      // Setup DMFT fields------------------------------------------------------
      std::vector<single_particle_greensfunction_type> Glocs;
      for(int i=0; i < Nlayer; i++)Glocs.push_back(single_particle_greensfunction_type(nt, ntau));
      if(put_nonint_leads==1)
      {
         Glocs.push_back(single_particle_greensfunction_type(nt, ntau)); Glocs[Nlayer].null()  ;//TOP BATH-->FILLED
         Glocs.push_back(single_particle_greensfunction_type(nt, ntau)); Glocs[Nlayer+1].null();//BOT BATH-->EMPTY
      }
      //
      std::vector<hybridization_function_type> Deltas;
      for(int i=0; i < Nlayer; i++)Deltas.push_back(hybridization_function_type(nt,ntau));

      //Setup pp calculator-----------------------------------------------------
      ppsc::hilbert_spaces::three_band_fermi_diag hilbert_space;
      hilbert_space.init(); //std::cout<<hilbert_space;
      solver_type *Solvers;
      Solvers = new solver_type[Nlayer](nt, ntau, beta, h, kt, nomp, hilbert_space, order);

      //Read Hamiltonian parameters-----------------------------------------------
      int readlimit;
      if(Nlayer==6)readlimit=Nlayer-2;
      if(Nlayer==1)readlimit=Nlayer;
      for(int i=0; i < readlimit; i++)
      {
         find_param(argv[1]        , "__mu=" , Solvers[i].imp.hamiltonian.mu     );
         find_param_tvector(argv[1], "__U="  , Solvers[i].imp.hamiltonian.U  , nt);
         find_param_tvector(argv[1], "__J="  , Solvers[i].imp.hamiltonian.J  , nt);
         find_param_tvector(argv[1], "__Bz=" , Solvers[i].imp.hamiltonian.Bz , nt);
         Solvers[i].imp.hamiltonian.site = i;
      }
      if(Nlayer==6)
      {
         // TOP INTERFACE --> FILLED
         find_param(argv[1]        , "__mu="   , Solvers[4].imp.hamiltonian.mu     );
         find_param_tvector(argv[1], "__Utop=" , Solvers[4].imp.hamiltonian.U  , nt);
         find_param_tvector(argv[1], "__Jtop=" , Solvers[4].imp.hamiltonian.J  , nt);
         find_param(argv[1],         "__Eloc_top="     , Eloc_top                  );
         find_param(argv[1],         "__lead_beta_top=", lead_beta_top             );
         find_param(argv[1],         "__lead_W_top="   , lead_W_top                );
         Solvers[4].imp.hamiltonian.site = 4;
         // BOTTOM INTERFACE --> EMPTY
         find_param(argv[1]        , "__mu="   , Solvers[5].imp.hamiltonian.mu     );
         find_param_tvector(argv[1], "__Ubot=" , Solvers[5].imp.hamiltonian.U  , nt);
         find_param_tvector(argv[1], "__Jbot=" , Solvers[5].imp.hamiltonian.J  , nt);
         find_param(argv[1],         "__Eloc_bot="     , Eloc_bot                  );
         find_param(argv[1],         "__lead_beta_bot=", lead_beta_bot             );
         find_param(argv[1],         "__lead_W_bot="   , lead_W_bot                );
         Solvers[5].imp.hamiltonian.site = 5;
      }

      // LOCAL HAMILTONIAN STRUCTURE--------------------------------------------
      // time dependent stuff that I'm gonna keep
      // 1) Null hopping and template
      cntr::function<double> t0(nt,Norb);
      // 2) vertical hybridizations
      std::vector<cntr::function<double>> tva,tvb; //ultimi due da togliere
      // 3) in-plane nearest neighbor hoppings
      std::vector<cntr::function<double>> t_nn_Uc,t_nn_Rd,t_nn_Ry,t_nn_Rx;
      // 4) in-plane next nearest neighbor hoppings
      std::vector<cntr::function<double>> tnnnRxp,tnnnRyp,tnnnRxm,tnnnRym;
      //template initialization
      for(int i=0; i < Nlayer; i++)
      {
         tva.push_back(t0); tvb.push_back(t0);
         t_nn_Uc.push_back(t0); t_nn_Rd.push_back(t0); t_nn_Ry.push_back(t0); t_nn_Rx.push_back(t0);
         tnnnRxp.push_back(t0); tnnnRyp.push_back(t0); tnnnRxm.push_back(t0); tnnnRym.push_back(t0);
      }
      // time dependent stuff that I'm gonna update while reading
      {
         //full H matrix at each tstp
         Eigen::MatrixXcd hloc(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopXp(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopXm(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopYp(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopYm(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopZp(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopZm(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopDL(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         Eigen::MatrixXcd hopDR(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
         //null hopping and template
         Eigen::MatrixXcd tnull(Norb,Norb);
         tnull.setZero(Norb,Norb);
         t0.set_constant(tnull);
         //dum var for reading
         double locdumR,locdumI;
         double hopXp_dumR,hopXp_dumI,hopXm_dumR,hopXm_dumI;
         double hopYp_dumR,hopYp_dumI,hopYm_dumR,hopYm_dumI;
         double hopZp_dumR,hopZp_dumI,hopZm_dumR,hopZm_dumI;
         double hopDL_dumR,hopDL_dumI,hopDR_dumR,hopDR_dumI;
         std::string hdr, ph;
         //matrix at each timestep
         std::vector<Eigen::MatrixXcd> hva,hvb;
         std::vector<Eigen::MatrixXcd> h_nn_Uc,h_nn_Rd,h_nn_Ry,h_nn_Rx;
         std::vector<Eigen::MatrixXcd> hnnnRxp,hnnnRyp,hnnnRxm,hnnnRym;
         //matrix initialization
         for(int i=0; i < Nlayer; i++)
         {
            hva.push_back(tnull); hvb.push_back(tnull);
            h_nn_Uc.push_back(tnull); h_nn_Rd.push_back(tnull); h_nn_Ry.push_back(tnull); h_nn_Rx.push_back(tnull);
            hnnnRxp.push_back(tnull); hnnnRyp.push_back(tnull); hnnnRxm.push_back(tnull); hnnnRym.push_back(tnull);
         }
         // read Hamilt---------------------------------------------------------
         std::cout << "read Hamiltonians " << endl;

         hdr="RE_Hloct_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_Hloct(hdr.append(ph).append(".dat"));
         hdr="IM_Hloct_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_Hloct(hdr.append(ph).append(".dat"));

         hdr="RE_HopptXp_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppXp(hdr.append(ph).append(".dat"));
         hdr="IM_HopptXp_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppXp(hdr.append(ph).append(".dat"));
         hdr="RE_HopptXm_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppXm(hdr.append(ph).append(".dat"));
         hdr="IM_HopptXm_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppXm(hdr.append(ph).append(".dat"));

         hdr="RE_HopptYp_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppYp(hdr.append(ph).append(".dat"));
         hdr="IM_HopptYp_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppYp(hdr.append(ph).append(".dat"));
         hdr="RE_HopptYm_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppYm(hdr.append(ph).append(".dat"));
         hdr="IM_HopptYm_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppYm(hdr.append(ph).append(".dat"));

         hdr="RE_HopptZp_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppZp(hdr.append(ph).append(".dat"));
         hdr="IM_HopptZp_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppZp(hdr.append(ph).append(".dat"));
         hdr="RE_HopptZm_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppZm(hdr.append(ph).append(".dat"));
         hdr="IM_HopptZm_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppZm(hdr.append(ph).append(".dat"));

         hdr="RE_HopptDL_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppDL(hdr.append(ph).append(".dat"));
         hdr="IM_HopptDL_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppDL(hdr.append(ph).append(".dat"));
         hdr="RE_HopptDR_ph"; ph = std::to_string(photon_ndx);
         ifstream RE_HoppDR(hdr.append(ph).append(".dat"));
         hdr="IM_HopptDR_ph"; ph = std::to_string(photon_ndx);
         ifstream IM_HoppDR(hdr.append(ph).append(".dat"));
         /*
         FILE * REFile;
         FILE * IMFile;
         REFile = fopen ("REtest_used.out","w");
         IMFile = fopen ("IMtest_used.out","w");
         */

         for (int it= -1; it <= nt; it++)
         {
            // reading----------------------------------------------------------
            //equil==1: read only the first timestep once
            //equil==0:  read each timestep
            if((equil==1 && it==-1)||(equil==0))
            {
               hloc.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopXp.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopXm.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopYp.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopYm.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopZp.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopZm.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopDL.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               hopDR.setZero(Nlat*Nspin*Norb,Nlat*Nspin*Norb);
               //
               for (int i= 0; i < Nlat*Nspin*Norb; i++)
               {
                  for (int j = 0; j < Nlat*Nspin*Norb; j++)
                  {
                     //reading values
                     RE_Hloct  >> locdumR    ; IM_Hloct  >> locdumI;
                     RE_HoppXp >> hopXp_dumR ; IM_HoppXp >> hopXp_dumI;
                     RE_HoppXm >> hopXm_dumR ; IM_HoppXm >> hopXm_dumI;
                     RE_HoppYp >> hopYp_dumR ; IM_HoppYp >> hopYp_dumI;
                     RE_HoppYm >> hopYm_dumR ; IM_HoppYm >> hopYm_dumI;
                     RE_HoppZp >> hopZp_dumR ; IM_HoppZp >> hopZp_dumI;
                     RE_HoppZm >> hopZm_dumR ; IM_HoppZm >> hopZm_dumI;
                     RE_HoppDL >> hopDL_dumR ; IM_HoppDL >> hopDL_dumI;
                     RE_HoppDR >> hopDR_dumR ; IM_HoppDR >> hopDR_dumI;
                     //assignment
                      hloc(i,j) = std::complex<double>(locdumR,locdumI);
                     hopXp(i,j) = std::complex<double>(hopXp_dumR,hopXp_dumI);
                     hopXm(i,j) = std::complex<double>(hopXm_dumR,hopXm_dumI);
                     hopYp(i,j) = std::complex<double>(hopYp_dumR,hopYp_dumI);
                     hopYm(i,j) = std::complex<double>(hopYm_dumR,hopYm_dumI);
                     hopZp(i,j) = std::complex<double>(hopZp_dumR,hopZp_dumI);
                     hopZm(i,j) = std::complex<double>(hopZm_dumR,hopZm_dumI);
                     hopDL(i,j) = std::complex<double>(hopDL_dumR,hopDL_dumI);
                     hopDR(i,j) = std::complex<double>(hopDR_dumR,hopDR_dumI);
                     //fprintf (REFile, "%12.6f",hloc(i,j).real());
                     //fprintf (IMFile, "%12.6f",hloc(i,j).imag());
                  }
                  //fprintf (REFile, "\n");
                  //fprintf (IMFile, "\n");
               }
            }

            // extract hybridizations-------------------------------------------
            // BULK LVO
            if(setup==1)
            {
               //vertical
               hva[0] =hopZp.block( 0*Nspin, 6*Nspin, 3, 3);
               hvb[0] = hloc.block( 0*Nspin, 6*Nspin, 3, 3);
               for(int i=0; i < Nlayer; i++)
               {
                  //in-plane nearest neighbor
                  h_nn_Uc[i] = hloc.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Rd[i] =hopDL.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Ry[i] =hopYp.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Rx[i] =hopXm.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  //in-plane next-nearest neighbor
                  hnnnRxp[i]=hopXp.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRyp[i]=hopYp.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRxm[i]=hopXm.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRym[i]=hopYm.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
               }

            }
            // HETERO LVO
            if(setup==2)
            {
               //vertical
               // SITE 9 - R
               hva[4] = hloc.block(24*Nspin,30*Nspin, 3, 3)*HtoBaths_Z;
               hvb[4] =hopZm.block(24*Nspin,18*Nspin, 3, 3);
               // SITE 7
               hva[3] =hopZp.block(18*Nspin,24*Nspin, 3, 3);
               hvb[3] = hloc.block(18*Nspin,12*Nspin, 3, 3);
               // SITE 5
               hva[2] = hloc.block(12*Nspin,18*Nspin, 3, 3);
               hvb[2] = hloc.block(12*Nspin, 6*Nspin, 3, 3);
               // SITE 3
               hva[1] = hloc.block( 6*Nspin,12*Nspin, 3, 3);
               hvb[1] = hloc.block( 6*Nspin, 0*Nspin, 3, 3);
               // SITE 1
               hva[0] = hloc.block( 0*Nspin, 6*Nspin, 3, 3);
               hvb[0] = hloc.block( 0*Nspin,30*Nspin, 3, 3);
               // SITE 11 - L
               hva[5] = hloc.block(30*Nspin, 0*Nspin, 3, 3);
               hvb[5] = hloc.block(30*Nspin,24*Nspin, 3, 3)*HtoBaths_Z;
               //
               for(int i=0; i < Nlayer; i++)
               {
                  //in-plane nearest neighbor
                  h_nn_Uc[i] = hloc.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Rd[i] =hopDR.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Ry[i] =hopYm.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Rx[i] =hopXp.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  //in-plane next-nearest neighbor
                  hnnnRxp[i]=hopXp.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRyp[i]=hopYp.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRxm[i]=hopXm.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRym[i]=hopYm.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
               }
            }
            // HETERO YTO
            if(setup==3)
            {
               //vertical
               // SITE 9 - R
               hva[4] = hloc.block(24*Nspin,30*Nspin, 3, 3)*HtoBaths_Z;
               hvb[4] =hopZm.block(24*Nspin,18*Nspin, 3, 3);
               // SITE 7
               hva[3] =hopZp.block(18*Nspin,24*Nspin, 3, 3);
               hvb[3] = hloc.block(18*Nspin,12*Nspin, 3, 3);
               // SITE 5
               hva[2] = hloc.block(12*Nspin,18*Nspin, 3, 3);
               hvb[2] = hloc.block(12*Nspin, 6*Nspin, 3, 3);
               // SITE 3
               hva[1] = hloc.block( 6*Nspin,12*Nspin, 3, 3);
               hvb[1] = hloc.block( 6*Nspin, 0*Nspin, 3, 3);
               // SITE 1
               hva[0] = hloc.block( 0*Nspin, 6*Nspin, 3, 3);
               hvb[0] = hloc.block( 0*Nspin,30*Nspin, 3, 3);
               // SITE 11 - L
               hva[5] = hloc.block(30*Nspin, 0*Nspin, 3, 3);
               hvb[5] = hloc.block(30*Nspin,24*Nspin, 3, 3)*HtoBaths_Z;
               //
               for(int i=0; i < Nlayer; i++)
               {
                  //in-plane nearest neighbor
                  h_nn_Uc[i] = hloc.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Rd[i] =hopDL.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Ry[i] =hopYp.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  h_nn_Rx[i] =hopXm.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3)*HinplaneXY;
                  //in-plane next-nearest neighbor
                  hnnnRxp[i]=hopXp.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRyp[i]=hopYp.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRxm[i]=hopXm.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
                  hnnnRym[i]=hopYm.block((6*i)*Nspin,(6*i)*Nspin, 3, 3)*HinplaneXY;
               }
            }

            //assignment to cntr functions--------------------------------------
            for(int i=0; i < Nlayer; i++)
            {
               //vertical
               tva[i].set_value(it,hva[i]);
               tvb[i].set_value(it,hvb[i]);
               //in-plane nearest neighbor
               t_nn_Uc[i].set_value(it,h_nn_Uc[i]);
               t_nn_Rd[i].set_value(it,h_nn_Rd[i]);
               t_nn_Ry[i].set_value(it,h_nn_Ry[i]);
               t_nn_Rx[i].set_value(it,h_nn_Rx[i]);
               //in-plane next-nearest neighbor
               tnnnRxp[i].set_value(it,hnnnRxp[i]);
               tnnnRyp[i].set_value(it,hnnnRyp[i]);
               tnnnRxm[i].set_value(it,hnnnRxm[i]);
               tnnnRym[i].set_value(it,hnnnRym[i]);
            }

            // extract local energies-------------------------------------------
            for(int i=0; i < Nlayer; i++)
            {
               int orb1, orb2, orb3;
               orb1 = 0 + (2*i)*Norb*Nspin;
               orb2 = 1 + (2*i)*Norb*Nspin;
               orb3 = 2 + (2*i)*Norb*Nspin;
               Solvers[i].imp.hamiltonian.E1[it+1] = hloc( orb1 , orb1 ).real();
               Solvers[i].imp.hamiltonian.E2[it+1] = hloc( orb2 , orb2 ).real();
               Solvers[i].imp.hamiltonian.E3[it+1] = hloc( orb3 , orb3 ).real();
            }

            // gradient adjustment and Vbias---------------------->(!CUSTOMIZE!)
            if(Nlayer==6)
            {
               Solvers[4].imp.hamiltonian.E1[it+1] += Eloc_top;
               Solvers[4].imp.hamiltonian.E2[it+1] += Eloc_top;
               Solvers[4].imp.hamiltonian.E3[it+1] += Eloc_top;
               //
               Solvers[5].imp.hamiltonian.E1[it+1] += Eloc_bot;
               Solvers[5].imp.hamiltonian.E2[it+1] += Eloc_bot;
               Solvers[5].imp.hamiltonian.E3[it+1] += Eloc_bot;
               //
               if(AdjYTO!=0.0)
               {
                  Solvers[1].imp.hamiltonian.E1[it+1] += AdjYTO;
                  Solvers[1].imp.hamiltonian.E2[it+1] += AdjYTO;
                  Solvers[1].imp.hamiltonian.E3[it+1] += AdjYTO;
               }
               //
               if(gradient!=0.0 || Vbias!=0.0)
               {
                  Solvers[0].imp.hamiltonian.E1[it+1] += + 0*gradient + Vbiasv[it+1]/2;
                  Solvers[0].imp.hamiltonian.E2[it+1] += + 0*gradient + Vbiasv[it+1]/2;
                  Solvers[0].imp.hamiltonian.E3[it+1] += + 0*gradient + Vbiasv[it+1]/2;

                  Solvers[1].imp.hamiltonian.E1[it+1] += + 1*gradient + Vbiasv[it+1]*0.166666666666;
                  Solvers[1].imp.hamiltonian.E2[it+1] += + 1*gradient + Vbiasv[it+1]*0.166666666666;
                  Solvers[1].imp.hamiltonian.E3[it+1] += + 1*gradient + Vbiasv[it+1]*0.166666666666;

                  Solvers[2].imp.hamiltonian.E1[it+1] += + 2*gradient - Vbiasv[it+1]*0.166666666666;
                  Solvers[2].imp.hamiltonian.E2[it+1] += + 2*gradient - Vbiasv[it+1]*0.166666666666;
                  Solvers[2].imp.hamiltonian.E3[it+1] += + 2*gradient - Vbiasv[it+1]*0.166666666666;

                  Solvers[3].imp.hamiltonian.E1[it+1] += + 3*gradient - Vbiasv[it+1]/2;
                  Solvers[3].imp.hamiltonian.E2[it+1] += + 3*gradient - Vbiasv[it+1]/2;
                  Solvers[3].imp.hamiltonian.E3[it+1] += + 3*gradient - Vbiasv[it+1]/2;

                  Solvers[4].imp.hamiltonian.E1[it+1] += - Vbiasv[it+1]/2;
                  Solvers[4].imp.hamiltonian.E2[it+1] += - Vbiasv[it+1]/2;
                  Solvers[4].imp.hamiltonian.E3[it+1] += - Vbiasv[it+1]/2;

                  Solvers[5].imp.hamiltonian.E1[it+1] += + Vbiasv[it+1]/2;
                  Solvers[5].imp.hamiltonian.E2[it+1] += + Vbiasv[it+1]/2;
                  Solvers[5].imp.hamiltonian.E3[it+1] += + Vbiasv[it+1]/2;
               }
            }
            /*
            if(it==300)
            {
               int i;
               i=0;
               std::ofstream file("./hoptest.out");
               Eigen::IOFormat myFmt(2, 0, ", ", ";\n", "", "", "[", "]");
               //
               file << '\n'<< '\n'<< '\n';
               file << "##########   in-plane nearest neighbor  ##########" << '\n';
               file << "RE Ryp" << '\n';
               file << hopYp.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3).real().format(myFmt) << '\n'<< '\n';
               file << "IM Ryp" << '\n';
               file << hopYp.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3).imag().format(myFmt) << '\n'<< '\n';
               file << "RE Rym" << '\n';
               file << hopYm.block((6*i+3)*Nspin,(6*i)*Nspin, 3, 3).real().format(myFmt) << '\n'<< '\n';
               file << "IM Rym" << '\n';
               file << hopYm.block((6*i+3)*Nspin,(6*i)*Nspin, 3, 3).imag().format(myFmt) << '\n'<< '\n';
               file << "RE Ryp.adj" << '\n';
               file << hopYp.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3).adjoint().eval().real().format(myFmt) << '\n'<< '\n';
               file << "IM Ryp.adj" << '\n';
               file << hopYp.block((6*i)*Nspin,(6*i+3)*Nspin, 3, 3).adjoint().eval().imag().format(myFmt) << '\n'<< '\n';
               file.close();
            }
            */
         }//ENDloop on time
         RE_Hloct.close() ; IM_Hloct.close();
         RE_HoppXp.close();IM_HoppXp.close();RE_HoppXm.close();IM_HoppXm.close();
         RE_HoppYp.close();IM_HoppYp.close();RE_HoppYm.close();IM_HoppYm.close();
         RE_HoppZp.close();IM_HoppZp.close();RE_HoppZm.close();IM_HoppZm.close();
         RE_HoppDL.close();IM_HoppDL.close();RE_HoppDR.close();IM_HoppDR.close();
         //fclose (REFile);
         //fclose (IMFile);
      }//END scope hopping assignment


      // GLOC CONNECTIONS STRUCTURE------------------------------->(!CUSTOMIZE!)
      int map[Nlayer][4];
      if(Nlayer==6)
      {
         //PERODIC SLAB WITH BATHS (that can be also vanishing)
         //      Gloc         Gplane           Gv_a           Gv_b
         map[0][0]= 0 ; map[0][1]= 0 ; map[0][2]= 1 ; map[0][3]= 5 ;
         map[1][0]= 1 ; map[1][1]= 1 ; map[1][2]= 2 ; map[1][3]= 0 ;
         map[2][0]= 2 ; map[2][1]= 2 ; map[2][2]= 3 ; map[2][3]= 1 ;
         map[3][0]= 3 ; map[3][1]= 3 ; map[3][2]= 4 ; map[3][3]= 2 ;
         map[4][0]= 4 ; map[4][1]= 4 ; map[4][2]= 6 ; map[4][3]= 3 ;
         map[5][0]= 5 ; map[5][1]= 5 ; map[5][2]= 0 ; map[5][3]= 7 ;
      }
      if(Nlayer==1)for(int j=0; j < 4; j++) map[0][j]=0;

      // FLAT BATH SETUP--------------------------------------------------------
      if(put_nonint_leads==1)
      {
         double shift;
         // TOP BATH-->FILLED <--imp[4]
         shift=3.22283+0.3+0.2;
         buildBath(Glocs[Nlayer]  , Solvers[4], lead_W_top, Vbias, shift);
         {
            std::string ph = std::to_string(photon_ndx);
            std::string filename = "ppsc_L_topB_ph_"; filename.append(ph).append(".out");
            std::cout << filename << std::endl;
            hid_t file_id = open_hdf5_file(filename);
            std::cout << "opened" << std::endl;
            Glocs[Nlayer].store(file_id);
            std::cout << "Gloc stored" << std::endl;
            close_hdf5_file(file_id);
            std::cout << "closed" << std::endl;
         }
         // BOT BATH-->EMPTY <--imp[5]
         shift=1.31717-0.3;
         buildBath(Glocs[Nlayer+1], Solvers[5], lead_W_bot, Vbias, shift);
         {
            std::string ph = std::to_string(photon_ndx);
            std::string filename = "ppsc_L_botB_ph_"; filename.append(ph).append(".out");
            std::cout << filename << std::endl;
            hid_t file_id = open_hdf5_file(filename);
            std::cout << "opened" << std::endl;
            Glocs[Nlayer+1].store(file_id);
            std::cout << "Gloc stored" << std::endl;
            close_hdf5_file(file_id);
            std::cout << "closed" << std::endl;
         }
      }

      // UPDATE HAMILTONIAN-----------------------------------------------------
      for(int i=0; i < Nlayer; i++)Solvers[i].imp.update_hamiltonian();
      FILE * ElocFile; ElocFile = fopen ("LocE.out","w");
      for (int it=0; it <nt+2; it++)
      {
         for(int i=0; i < Nlayer; i++)fprintf (ElocFile,"%.12e\t%.12e\t%.12e\t",
         Solvers[i].imp.hamiltonian.E1[it],Solvers[i].imp.hamiltonian.E2[it],Solvers[i].imp.hamiltonian.E3[it]);
         fprintf (ElocFile,"\n");
      }
      fclose(ElocFile);

      // MATSUBARA PART - EQUILIBRIUM INITIAL STATE-----------------------------
      {
         // -- Setup vars for error check
         std::vector<cntr::herm_matrix<double>> Gerr;
         double dmfterr_equil[Nlayer];
         for(int i=0; i < Nlayer; i++)
         {
            cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);
            Gerr.push_back(gtmp);
         }
         // -- DMFT loops
         for (iter_equil = 1; iter_equil <= itermax; iter_equil++)
         {
            // -- Seed for magnetic calculations
            if(af_bethe_flag==1)
            {
               double Bz_seed=0.01;
               if(iter_equil>10 && iter_equil<=20 )
               {
                  std::cout << "--> Applying Bz_seed = " << Bz_seed << std::endl;
                  for(int i=0; i < Nlayer-2; i++)
                  {
                     Solvers[i].imp.hamiltonian.Bz[0] = -Bz_seed;
                     Solvers[i].imp.update_hamiltonian(-1);
                  }
               }
               else
               {
                  for(int i=0; i < Nlayer-2; i++)
                  {
                     Solvers[i].imp.hamiltonian.Bz[0] = 0.0;
                     Solvers[i].imp.update_hamiltonian(-1);
                  }
               }
            }
            // -- Initialize ppsc
            if(iter_equil==1)for(int i=0; i < Nlayer; i++)Solvers[i].imp.solve_atomic();
            // -- loop over sites
            for(int i=0; i < Nlayer; i++)
            {
               // -- Construct interactions
               std::cout << "--> Building G "<<i+1<< std::endl;
               ppsc::pp_ints_type pp_ints = get_pp_ints(Deltas[i], Solvers[i].imp.hamiltonian);
               ppsc::gf_verts_type gf_verts = get_gf_verts(Solvers[i].imp.hamiltonian);
               Solvers[i].imp.update_diagrams(pp_ints, gf_verts);
               if(iter_equil == 1 && read_eq_sym && !Solvers[i].imp.has_symmetries())
               {
                  std::string site = std::to_string(i+1);
                  std::string filename = "sym_eq_site_";
                  filename.append(site).append(".out");
                  Solvers[i].imp.read_symmetries(filename);
               }
               // -- Solve pseudo particle problem
               Solvers[i].imp.pp_step(-1);
               // -- Get spgf and mix
               ppsc::gf_tstps_type gf_tstps = Solvers[i].imp.get_spgf(-1);
               Glocs[i].update(-1, gf_tstps, linear_mixing, af_bethe_flag);
               // -- Compute error
               dmfterr_equil[i] = cntr::distance_norm2(-1, Gerr[i], Glocs[i].o1up);
               Gerr[i].set_timestep(-1, Glocs[i].o1up);
            }
            // -- Update Hybridizations
            for(int i=0; i < Nlayer; i++)
            {
               int mag_flag=af_bethe_flag;
               int ph_flag=LVO_phflag;
               if(i>=4) mag_flag=0;
               if(i>=4) ph_flag=Ti_phflag;
               Deltas[i].update( -1 ,  i , Glocs[map[i][0]]   // Gloc
                                         , Glocs[map[i][1]]   // Gplane
                                         , Glocs[map[i][2]]   // Gv_a
                                         , Glocs[map[i][3]]   // Gv_b
                                         , t_nn_Uc[i]
                                         , t_nn_Rd[i]
                                         , t_nn_Ry[i]
                                         , t_nn_Rx[i]
                                         , tnnnRxp[i]
                                         , tnnnRyp[i]
                                         , tnnnRxm[i]
                                         , tnnnRym[i]
                                         , tva[i]
                                         , tvb[i]
                                         , psi, psidag
                                         , D0                 // phonon bubble
                                         , mag_flag           // AFM
                                         , ph_flag            // put phonon
                                         , 0                  // compute current
                                         , beta, h , kt
                                         , gradientBulk, Nlayer);
            }
            // -- Error check
            double dmfterr=0.0;
            for(int i=0; i < Nlayer; i++) dmfterr+= dmfterr_equil[i]/((double) Nlayer);
            cout << "Eq. iter:  " << iter_equil<< " err: " << dmfterr << endl;
            if( dmfterr < errmax)
            {
               matsubara_converged = true;
               cout << "MATSUBARA CONVERGED" << endl;
               break;
            }
         } // END DMFT loops
         if (iter_equil > itermax)
         {
            cerr << "WARNING: Matsubara not converged  after " << itermax
            << "steps ... abort" << endl;
            cerr << "skip real-time calculation " << endl;
         }
      } // END MATSUBARA

      // save symmetries--------------------------------------------------------
      if(sym_flag)
      {
         for(int i=0; i < Nlayer; i++)
         {
            if(!Solvers[i].imp.has_symmetries())
            {
               std::string site = std::to_string(i+1);
               std::string filename = "sym_eq_site_";
               filename.append(site).append(".out");
               Solvers[i].imp.symmetry_reduction(-1);
               Solvers[i].imp.write_symmetries(filename);
            }
         }
      }

      // REAL TIME PART - BOOTSTRAP---------------------------------------------
      if (nt > 0 && matsubara_converged == true)
      {
         if(sym_flag && nt > 0)
         {
            for(int i=0; i < Nlayer; i++)Solvers[i].imp.clear_symmetries();
            if(read_rt_sym)
            {
               for(int i=0; i < Nlayer; i++)
               {
                  std::string site = std::to_string(i+1);
                  std::string filename = "sym_rt_site_";
                  filename.append(site).append(".out");
                  Solvers[i].imp.read_symmetries(filename);
               }
            }
         }
         warmup_converged = false;
         // -- Setup vars for error check
         std::vector<cntr::herm_matrix<double>> Gerr;
         double dmfterr_warm[Nlayer];
         for(int i=0; i < Nlayer; i++)
         {
            cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
            Gerr.push_back(gtmp);
         }
         // -- Initialize ppsc
         if(iter_warmup==1)for(int i=0; i < Nlayer; i++)Solvers[i].imp.init_real_time();
         // -- DMFT loops
         for (iter_warmup = 1; iter_warmup <= itermax; iter_warmup++)
         {
            // -- loop over sites
            for(int i=0; i < Nlayer; i++)
            {
               // -- Construct interactions
               std::cout << "--> Building G "<<i+1<< std::endl;
               ppsc::pp_ints_type pp_ints = get_pp_ints(Deltas[i], Solvers[i].imp.hamiltonian);
               ppsc::gf_verts_type gf_verts = get_gf_verts(Solvers[i].imp.hamiltonian);
               Solvers[i].imp.update_diagrams(pp_ints, gf_verts);
               // -- Solve pseudo particle problem
               Solvers[i].imp.pp_step(kt);
               // -- Get spgf
               for (int n = 0; n <= kt; n++)
               {
                  ppsc::gf_tstps_type gf_tstps = Solvers[i].imp.get_spgf(n);
                  Glocs[i].update(n, gf_tstps, af_bethe_flag);
               }
               // -- Compute error
               dmfterr_warm[i] = cntr::distance_norm2(kt, Gerr[i], Glocs[i].o1up);
               Gerr[i].set_timestep(kt, Glocs[i].o1up);
            }
            // -- Update Hybridizations
            for (int n = 0; n <= kt; n++)
            {
               for(int i=0; i < Nlayer; i++)
               {
                  int mag_flag=af_bethe_flag;
                  int ph_flag=LVO_phflag;
                  if(i>=4) mag_flag=0;
                  if(i>=4) ph_flag=Ti_phflag;
                  Deltas[i].update(  n ,  i , Glocs[map[i][0]]   // Gloc
                                            , Glocs[map[i][1]]   // Gplane
                                            , Glocs[map[i][2]]   // Gv_a
                                            , Glocs[map[i][3]]   // Gv_b
                                            , t_nn_Uc[i]
                                            , t_nn_Rd[i]
                                            , t_nn_Ry[i]
                                            , t_nn_Rx[i]
                                            , tnnnRxp[i]
                                            , tnnnRyp[i]
                                            , tnnnRxm[i]
                                            , tnnnRym[i]
                                            , tva[i]
                                            , tvb[i]
                                            , psi, psidag
                                            , D0                 // phonon bubble
                                            , mag_flag           // AFM
                                            , ph_flag            // put phonon
                                            , 0                  // compute current
                                            , beta, h , kt
                                            , gradientBulk, Nlayer);
               }
            }
            // -- Error check
            double dmfterr=0.0;
            for(int i=0; i < Nlayer; i++) dmfterr+= dmfterr_warm[i]/((double) Nlayer);
            cout << "WARMUP: iter:  " << iter_warmup << " err: " << dmfterr << endl;
            if (dmfterr < errmax)
            {
               warmup_converged = true;
               break;
            }
         } // END DMFT-WARMUP loops
      } // END BOOTSTRAP

      // save symmetries--------------------------------------------------------
      if(sym_flag && matsubara_converged && nt > 0 )
      {
         for(int i=0; i < Nlayer; i++)
         {
            if(!Solvers[i].imp.has_symmetries())
            {
               std::string site = std::to_string(i+1);
               std::string filename = "sym_rt_site_";
               filename.append(site).append(".out");
               Solvers[i].imp.symmetry_reduction(-1);
               Solvers[i].imp.write_symmetries(filename);
            }
         }
      }

      // REAL TIME PART - NON EQUILIBRIUM STATE---------------------------------
      if(matsubara_converged && warmup_converged)

      // Real time steps
      for (tstp = kt + 1; tstp <= nt; tstp++)
      {
         // -- Initialize ppsc
         for(int i=0; i < Nlayer; i++)Solvers[i].imp.extrapolate_timestep(tstp - 1);
         // DMFT loops
         for (iter_rt = 1; iter_rt <= iter_rtime; iter_rt++)
         {
            std::cout << "iter_rtime:  "<< iter_rt << endl;
            // loop over sites
            for(int i=0; i < Nlayer; i++)
            {
               // -- Construct interactions
               std::cout << "--> Building G "<<i+1<< std::endl;
               ppsc::pp_ints_type pp_ints = get_pp_ints(Deltas[i], Solvers[i].imp.hamiltonian);
               ppsc::gf_verts_type gf_verts = get_gf_verts(Solvers[i].imp.hamiltonian);
               Solvers[i].imp.update_diagrams(pp_ints, gf_verts);
               // -- Solve pseudo particle problem
               Solvers[i].imp.pp_step(tstp);
               // -- Get spgf and mix
               ppsc::gf_tstps_type gf_tstps = Solvers[i].imp.get_spgf(tstp);
               Glocs[i].update(tstp, gf_tstps, linear_mixing, af_bethe_flag);
            }
            // -- Update Hybridizations
            int jf=0;
            if(iter_rt==iter_rtime) jf=j_flag;
            for(int i=0; i < Nlayer; i++)
            {
               int mag_flag=af_bethe_flag;
               int ph_flag=LVO_phflag;
               if(i>=4) mag_flag=0;
               if(i>=4) ph_flag=Ti_phflag;
               Deltas[i].update(tstp,  i , Glocs[map[i][0]]   // Gloc
                                         , Glocs[map[i][1]]   // Gplane
                                         , Glocs[map[i][2]]   // Gv_a
                                         , Glocs[map[i][3]]   // Gv_b
                                         , t_nn_Uc[i]
                                         , t_nn_Rd[i]
                                         , t_nn_Ry[i]
                                         , t_nn_Rx[i]
                                         , tnnnRxp[i]
                                         , tnnnRyp[i]
                                         , tnnnRxm[i]
                                         , tnnnRym[i]
                                         , tva[i]
                                         , tvb[i]
                                         , psi, psidag
                                         , D0                 // phonon bubble
                                         , mag_flag           // AFM
                                         , ph_flag            // put phonon
                                         , jf                 // compute current
                                         , beta, h , kt
                                         , gradientBulk, Nlayer);

            }
         } // END DMFT loops
      }// END time steps

      // STORE------------------------------------------------------------------
      std::string ph = std::to_string(photon_ndx);
      for(int i=0; i < Nlayer; i++)
      {
        std::string filename = "ppsc_L_";
        std::string site = std::to_string(i+1);
        filename.append(site).append("_ph_").append(ph).append(".out");
        std::cout << "site: "<< i << std::endl;
        std::cout << filename << std::endl;
        hid_t file_id = open_hdf5_file(filename);
        std::cout << "opened" << std::endl;
        Solvers[i].imp.store(file_id, store_pp);
        std::cout << "imp stored" << std::endl;
        Glocs[i].store(file_id);
        std::cout << "Gloc stored" << std::endl;
        Deltas[i].store(file_id);
        std::cout << "Delta stored" << std::endl;
        close_hdf5_file(file_id);
        std::cout << "closed" << std::endl;
     }
     std::cout << "  end of try  "<<std::endl;
  } // END of global try
  catch (char *message)
  {
    cerr << "exception\n**** " << message << " ****" << endl;
    cerr << "CDMFT input_file [ --test ]\n" << endl;
  } catch (...)
  {
    cerr << "unspecified exception " << endl;
    cerr << "\nCDMFT input_file [ --test ]\n" << endl;
  }
  std::cout << "  END OF ALL  "<<std::endl;
  return 0;
}



// ENERGY CONSERVATION------------------------------------------------------
/*
if( tstp>=(kt+2) && j_flag==1 )
{
   std::cout << "  photon used for energy cons " << photon_ndx  << std::endl;
   double Et[15][2*nt];
   ifstream Efld("E_t.dat");
   for (int it= 1; it <=nt; it++)
   {
      for (int j = 0; j < 12; j++)
      {
         Efld >> Et[j][it];
      }
   }
   std::cout << "  Efield read  " << std::endl;
   Efld.close();
   //
   double jintra[Nlayer][Nspin*Norb][2*nt];
   double dum;
   for(int ilayer=1; ilayer<=(Nlat-1); ilayer=ilayer+2)
   {
      std::string hdrp = "F_plane_L_"; std::string layndx = std::to_string(ilayer);
      hdrp.append(layndx).append(".out"); ifstream Fplane(hdrp);
      std::cout << "  read current:  " << ilayer << "  " << (ilayer-1)/2 << " file:  " << hdrp << std::endl;
      for (int it=kt+2; it <=tstp; it++)
      {
         for (int j = 0; j < 2*Nspin*Norb+2; j++)
         {
            dum=0;
            Fplane >> dum;
            if(j>=Nspin*Norb+2) jintra[(ilayer-1)/2][j-(Nspin*Norb+2)][it] = -dum;
         }
      }
      Fplane.close();
   }
   //
   double Einj[Nlayer][2*nt];
   double dum1up, dum1do;
   double dum2up, dum2do;
   double dum3up, dum3do;
   for(int ilayer=0; ilayer<Nlayer; ilayer++)
   {
      for (int it1=kt+2; it1 <=tstp; it1++)
      {
         dum=0;
         dum1up = 0; dum2up = 0; dum3up = 0;
         dum1do = 0; dum2do = 0; dum3do = 0;
         for (int it2=kt+2; it2 <=it1; it2++)
         {
            dum1up += Et[photon_ndx][it2] * jintra[ilayer][0][it2] * h;
            dum1do += Et[photon_ndx][it2] * jintra[ilayer][1][it2] * h;
            dum2up += Et[photon_ndx][it2] * jintra[ilayer][2][it2] * h;
            dum2do += Et[photon_ndx][it2] * jintra[ilayer][3][it2] * h;
            dum3up += Et[photon_ndx][it2] * jintra[ilayer][4][it2] * h;
            dum3do += Et[photon_ndx][it2] * jintra[ilayer][5][it2] * h;
         }
         //std::cout << " ilayer it1 dum1up dum1do dum2up dum2do dum3up dum3do  " << "  " << ilayer << "  "  << it1 << "  "  << dum1up << "  " << dum1do << "  " << dum2up << "  " << dum2do << "  " << dum3up << "  " << dum3do << std::endl;
         Einj[ilayer][it1] = dum1up + dum1do + dum2up + dum2do + dum3up + dum3do;
      }
   }
   FILE * EinjFile; EinjFile = fopen ("Einj.out","w");
   for (int it1=kt+2; it1 <=tstp; it1++)
   {
      fprintf (EinjFile, "%5i",it1);
      for(int ilayer=0; ilayer<Nlayer; ilayer++)
      {
         fprintf (EinjFile, "\t%.12e",Einj[ilayer][it1]);
      }
      fprintf (EinjFile, "\n");
   }
   fclose(EinjFile);
}
*/
