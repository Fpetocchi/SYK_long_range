#ifndef ___IMPURITY___
#define ___IMPURITY___

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <complex>
#include <cmath>
#include <math.h>
#include <time.h>
#include <string>
//
#include <set>
#include <valarray>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
//
#include "times.hpp"
#include "file_io.hpp"


//============================================================================//


class ct_hyb
{
public:
   ct_hyb(double mu, double beta, int Nspin, int Norb, int Ntau, int MaxTime, int Therm, int Nmeas, bool retarded ):
          mu(mu), Beta(beta),Nspin(Nspin), Norb(Norb), Ntau(Ntau),
          MaxTime(MaxTime), Therm(Therm), Nmeas(Nmeas), retarded(retarded) {}
   double get_mu(void)const {return mu;}
   double get_Beta(void)const {return Beta;}
   int get_Nspin(void)const {return Nspin;}
   int get_Norb(void)const {return Norb;}
   int get_Nflavor(void)const {return Nspin*Norb;}
   int get_Ntau(void)const {return Ntau;}
   int get_MaxTime(void)const {return MaxTime;}
   int get_Nmeas(void)const {return Nmeas;}
   bool get_retarded(void)const {return retarded;}
   //
   void init(path);
   void reset_mu(double mu_new) {mu = mu_new; mu_is_reset=true;}
   //
   void dostep(bool);
   //
   void print_results(path);
   bool is_thermalized() const;
   double work_done() const;

private:
   double                              mu;
   double                              Beta;
   int                                 Nspin;
   int                                 Norb;
   int                                 Ntau;
   int                                 MaxTime;
   int                                 Therm;
   int                                 Nmeas;
   bool                                retarded;
   //
   bool                                mu_is_reset=false;
   bool                                initialized=false;
   int                                 Nflavor=Nspin*Norb;                      // spin-orbitla flavours
   int                                 sweeps;                                  // sweeps done
   int                                 thermalization_sweeps=Therm;             // sweeps to be done for equilibration
   vector_t                            Eloc;                                    // mu-<\epsilon>
   dense_matrix                        U;                                       // U matrix
   std::vector<segment_container_t >   segments;                                // stores configurations with 0,1,... segments (but not full line)
   std::vector<int>                    full_line;                               // if 1 means that particle occupies full time-line
   std::vector<double>                 sign;                                    // sign of Z_n_up [actually not needed for density-density]
   std::vector<dense_matrix>           M;                                       // inverse hybridization matrix
   std::valarray<double>               G_Nmeas;                                 // Nmeasured GF
   hybridization_container_t           F;                                       // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu)
   Kfunct_container_t                  K_table;                                 // K function matrix for retarded interactions


   //---------------------------------------------------------------------------


   void init(path &inputDir)
   {
      //
      //set the local energies
      Eloc.resize(Nflavor);
      read_vector(inputDir+"/Eloc.DAT", Eloc);

      //
      //read the hybridization function
      F.resize(Ntau+1,vector_t(Nflavor));
      read_hybridization(inputDir+"/Delta.DAT", F);

      //
      //read the screening function
      if(retarded==true)
      {
         K_table.resize(Ntau+1,dense_matrix(Nflavor,Nflavor));
         read_Kfunct(inputDir+"/K.DAT", K_table);
      }

      //
      std::cout << " Solver is initialized." << std::endl;
      initialized=true;

   }


   //-------------------------------------------------------------------------


};


#endif
