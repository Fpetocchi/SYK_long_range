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
   ct_hyb(double mu, double beta, int Nspin, int Norb, int Ntau,
          int MaxTime, int CheckTime, int Therm, int Nmeas, bool retarded ):
          mu(mu), Beta(beta),Nspin(Nspin), Norb(Norb), Ntau(Ntau),
          MaxTime(MaxTime), CheckTime(CheckTime), Therm(Therm), Nmeas(Nmeas),
          retarded(retarded) {}
   //
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
   int                                 CheckTime;
   int                                 Therm;
   int                                 Nmeas;
   bool                                retarded;
   //
   duration                            Tstart,Tnow;
   bool                                mu_is_reset=false;
   bool                                initialized=false;
   int                                 Nflavor=Nspin*Norb;                      // spin-orbitla flavours
   int                                 sweeps;                                  // sweeps done
   int                                 thermalization_sweeps=Therm;             // sweeps to be done for equilibration
   vector_t                            Eloc;                                    // mu-<\epsilon>
   dense_matrix                        Uloc;                                       // U matrix
   std::vector<segment_container_t >   segments;                                // stores configurations with 0,1,... segments (but not full line)
   std::vector<int>                    full_line;                               // if 1 means that particle occupies full time-line
   std::vector<double>                 sign;                                    // sign of Z_n_up [actually not needed for density-density]
   std::vector<dense_matrix>           M;                                       // inverse hybridization matrix
   std::valarray<double>               G_Nmeas;                                 // Nmeasured GF
   hybridization_container_t           F,G;                                     // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu)
   Kfunct_container_t                  K_table;                                 // K function matrix for retarded interactions


   //---------------------------------------------------------------------------


   void start_timer(void)
   {
      Tstart = std::chrono::system_clock::now();
   }

   bool check_timer(void)
   {
      bool stop=false;
      Tnow = std::chrono::system_clock::now();
      std::chrono::duration<double> runtime_seconds = Tnow-Tstart;
      if(runtime_seconds.count() > MaxTime) stop=true;
      return stop;
   }


   //---------------------------------------------------------------------------


   void init(path &inputDir)
   {
      //
      //set the local energies
      Eloc.resize(Nflavor);
      read_vector(inputDir+"/Eloc.DAT", Eloc);
      for (int ifl=0; ifl<Nflavor; ifl++) Eloc(ifl) = mu - Eloc(ifl);

      //
      //read the hybridization function (std::vector<Eigen::VectorXd>)
      F.resize(Ntau+1,vector_t(Nflavor));
      read_hybridization(inputDir+"/Delta.DAT", F);

      //
      //read the istantaneous interaction (Eigen::MatrixXd)
      Uloc.resize(Nflavor,Nflavor);
      read_Umat(inputDir+"/Uloc.DAT", Uloc);

      //
      //read the screening function (std::vector<Eigen::MatrixXd>)
      if(retarded==true)
      {
         K_table.resize(Ntau+1,dense_matrix(Nflavor,Nflavor));
         read_Kfunct(inputDir+"/K.DAT", K_table);
      }

      //
      //initialize the Green's function (std::vector<Eigen::VectorXd>)
      G.resize(Ntau+1,vector_t(Nflavor));

      //
      //initialize segment container (std::vector<std::set<times>>)
      segments.resize(Nflavor);
      for (int ifl=0; ifl<Nflavor; ifl++)segments[ifl].clear();

      //
      //initialize full_line (std::vector<int>)
      full_line.resize(Nflavor,0);

      //
      // initialize seed aggiungi qualcosa che dipenda dal task
      unsigned long seed = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::system_clock::now().time_since_epoch()).count()-1566563000000;
      generator.seed(abs(seed));

      //
      //start the timer for the calculation
      start_timer();

      //
      std::cout << " Solver is initialized." << std::endl;
      initialized=true;

   }


   //---------------------------------------------------------------------------


   void dostep(bool &atBeta)
   {
      //
      if(!initialized)
      {
         std::cout << " Solver is not initialized - Exiting." << std::endl;
         exit(1);
      }

      //
      bool finalize=false;
      do
      {
         // The measurments I'm going to do anyway
         for (int imeas=0; imeas<Nmeas; imeas++)
         {
            //Flavor loop
             for (int ifl=0; ifl<FLAVORS; ifl++)
             {
                //
                // insert or remove full line
                 if (segments[ifl].size() == 0) insert_remove_full_line(rndm, Eloc(ifl), Umat, Beta, full_line[ifl], segments, full_line, ifl);

                 insert_remove_antisegment(rndm, Beta*rndm(), N, Beta, mu_e[ifl], u, F[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);

                 if (!full_line[ifl])
                 {
                     // local update
                     insert_remove_segment(rndm, Beta*rndm(), N, Beta, mu_e[ifl], u, F[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);

                     // shift segments
                     for (int k=0; k<N_shift; k++) shift_segment(rndm, segments[ifl], N, Beta, mu_e[ifl], u, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);

                     /*
                      // flip segment
                      for (int i=0; i<N_flip; i++)
                      flip_segment(rndm, segments_up, N, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
                     */
                 }

                 // measure perturbation order

                 if (segments[ifl].size()<N_order)
                     order_meas[ifl*N_order+segments[ifl].size()] += 1;


                 // measure Green functions

                 if (segments[ifl].size()>0) {
                     for (int i=0; i<M[ifl].size1(); i++) {
                         (i==0 ? it1 = segments[ifl].begin() : it1++);
                         for (int k=0; k<M[ifl].size1(); k++) {
                             (k==0 ? it2 = segments[ifl].begin() : it2++);
                             if (M[ifl](k,i)!=0) {
                                 double argument = it1->t_end()-it2->t_start();
                                 double bubble_sign=1;
                                 if (argument > 0) {
                                     bubble_sign = 1;
                                 }
                                 else {
                                     bubble_sign = -1;
                                     argument += Beta;
                                 }

                                 int index = argument/Beta*N+0.5;
                                 G_meas[ifl*(N+1)+index] += M[ifl](k,i)*bubble_sign/(Beta*Beta);
                             }
                         }
                     }
                 }

                 s *= sign[ifl];

                 n_meas[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/Beta;

             }

             sign_meas += s;

             // swap updates may be needed in case of trapping in a symmetry-broken state
             /*
              if (i%N_swap==1) {
              int orbital=rndm()*(FLAVORS/2);
              if (M[2*orbital].size1()!=0 && M[2*orbital+1].size1()!=0)
              swap_segments(rndm, Beta, F[orbital], F[orbital+1], segments[orbital], segments[orbital+1], full_line[orbital], full_line[orbital+1], sign[orbital], sign[orbital+1], M[orbital], M[orbital+1]);
              }
             */

             /*
              // global spin flip (does not require overlap calculation)
              if (i%N_swap==1){
              //check if doable
              bool ok=true;
              for(int j=0; j<FLAVORS; j++) if(M[j].size1()==0) ok=false;
              if(ok) swap_spins(rndm, Beta, FLAVORS, F, segments, full_line, sign, M);
              }
             */


         }

      }while( !finalize );
   }


};


#endif
