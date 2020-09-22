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


//============================================================================//


class ct_hyb
{
public:
   ct_hyb(int Nspin, int Norb, int Ntau, int MaxTime, int Therm, int Nmeas, bool retarded, path inputFolder ):
          Nspin(Nspin), Norb(Norb), Ntau(Ntau),
          MaxTime(MaxTime), Therm(Therm), Nmeas(Nmeas),
          retarded(retarded), inputFolder(inputFolder) {}
   int get_Nspin(void)const {return Nspin;}
   int get_Norb(void)const {return Norb;}
   int get_Flavors(void)const {return Nspin*Norb;}
   int get_Ntau(void)const {return Ntau;}
   int get_MaxTime(void)const {return MaxTime;}
   int get_Nmeas(void)const {return Nmeas;}
   bool get_retarded(void)const {return retarded;}
   path get_inputFolder(void)const {return inputFolder;}
   //
   void read_Eloc();
   void read_Hybridization();
   void read_Kfunct();
   void dostep();
   bool is_thermalized() const;
   double work_done() const;

private:
   int                                 Nspin;
   int                                 Norb;
   int                                 Ntau;
   int                                 MaxTime;
   int                                 Therm;
   int                                 Nmeas;
   bool                                retarded;
   path                                inputFolder;
   //
   int                                 Flavors=Nspin*Norb;                      // spin-orbitla flavours
   int                                 sweeps;                                  // sweeps done
   int                                 thermalization_sweeps=Therm;             // sweeps to be done for equilibration
   double                              mu;                                      // chemical potential
   vector_t                            mu_e;                                    // mu-<\epsilon>
   dense_matrix                        u;                                       // U matrix
   double                              t;                                       // bandwidth=4t (for semi-circle)
   std::vector<segment_container_t >   segments;                                // stores configurations with 0,1,... segments (but not full line)
   std::vector<int>                    full_line;                               // if 1 means that particle occupies full time-line
   std::vector<double>                 sign;                                    // sign of Z_n_up [actually not needed for density-density]
   std::vector<dense_matrix>           M;                                       // inverse hybridization matrix
   std::valarray<double>               G_Nmeas;                                 // Nmeasured GF
   hybridization_container_t           F;                                       // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu)
   Kfunct_container_t                  K_table;                                 // K function matrix for retarded interactions
   //-------------------------------------------------------------------------






   //-------------------------------------------------------------------------


};


#endif
