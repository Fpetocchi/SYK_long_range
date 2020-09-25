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
#include "observables.hpp"


//============================================================================//


class ct_hyb
{
public:
   ct_hyb(double mu, double beta, int Nspin, int Norb, int Ntau, int Norder,
          int MaxTime, int Therm, int Nmeas, bool retarded ):
          mu(mu), Beta(beta),Nspin(Nspin), Norb(Norb), Ntau(Ntau), Norder(Norder),
          MaxTime(MaxTime), Therm(Therm), Nmeas(Nmeas), retarded(retarded) {}
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
   void solve(bool);
   //
   void print_results(path);
   bool is_thermalized() const;
   double work_done() const;

private:
   // Pointer to the internal function
   void(ct_hyb::*dostep)(void) = NULL;
   // Input vars
   double                              mu;                                      // chemical potential
   double                              Beta;                                    // inverse temperature
   int                                 Nspin;
   int                                 Norb;
   int                                 Ntau;
   int                                 Norder;                                  // Max perturbation order
   int                                 MaxTime;                                 // Max runtime
   int                                 Therm;                                   // Thernalization sweeps
   int                                 Nmeas;                                   // Number of sweeps between expensive measurments
   bool                                retarded;                                // flag to switch on retarded interactions
   // Internal vars
   duration                            Tstart,Tnow;
   bool                                mu_is_reset=false;
   bool                                initialized=false;
   int                                 Nflavor=Nspin*Norb;                      // Spin-orbital flavours
   int                                 Ntau_p1=Ntau+1;                          // Spin-orbital flavours
   int                                 sweeps;                                  // Sweeps done
   int                                 thermalization_sweeps=Therm;             // REMOVE
   Vec                                 Eloc;                                    // mu-<\epsilon>
   Mat                                 Uloc;                                    // Istantaneous U matrix
   std::vector<segment_container_t>    segments;                                // Stores configurations with 0,1,... segments (but not full line)
   std::vector<int>                    full_line;                               // if 1 means that particle occupies full time-line
   VecMat                              M;                                       // Inverse hybridization matrix
   VecVec                              F;                                       // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu)
   VecVecVec                           K_table;                                 // K function matrix for retarded interactions
   // Observables
   double                              sign_meas;                               // Total Sign
   Vec                                 sign;                                    // Sign per flavor
   Vec                                 Nloc;                                    // Density per flavor
   Vec                                 Pert;                                    // Perturbation order
   Vec                                 Nhist;                                   //
   Vec                                 Szhist;                                  //
   VecVec                              G;                                       // Impurity Green's function
   VecVec                              nt;                                      // n(\tau)
   VecVec                              nnt;                                     // Charge susceptibility


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


   void init(path &inputDir, bool &atBeta)
   {
      //
      //set the local energies ( std::vector<double> )
      read_Vec(inputDir+"/Eloc.DAT", Eloc, Nflavor);
      for (int ifl=0; ifl<Nflavor; ifl++) Eloc[ifl] = mu - Eloc[ifl];

      //
      //read the hybridization function ( std::vector<std::vector<double>> )
      read_VecVec(inputDir+"/Delta.DAT", F, Nflavor, Ntau_p1, true);

      //
      //read the istantaneous interaction ( Eigen::MatrixXd )
      read_EigenMat(inputDir+"/Uloc.DAT", Uloc, Nflavor, Nflavor);

      //read the screening function ( std::vector<std::vector<std::vector<double>>> )
      if(retarded==true) read_VecVecVec(inputDir+"/K.DAT", K_table, Nflavor, Nflavor, Ntau_p1);

      //
      //initialize segment container ( std::vector<std::set<times>> )
      segments.resize(Nflavor);
      for (int ifl=0; ifl<Nflavor; ifl++)segments[ifl].clear();

      //
      //initialize full_line ( std::vector<int> )
      full_line.resize(Nflavor,0);

      //
      //initialize empty Inverse hybridization matrix M ( std::vector<Eigen::MatrixXd> )
      M.resize(Nflavor,Mat::Zero(Nflavor,Nflavor));

      //
      //initialize observables
      sign.resize(Nflavor,0.0);                                                 // ( std::vector<double> )
      Nloc.resize(Nflavor,0.0);                                                 // ( std::vector<double> )
      Nhist.resize(Nflavor+1,0.0);                                              // ( std::vector<double> )
      Szhist.resize(Nflavor/2+1,0.0);                                           // ( std::vector<double> )
      Pert.resize(Norder*Nflavor,0.0);                                          // ( std::vector<double> )
      G.resize(Nflavor,std::vector<double>(Ntau_p1,0.0));                       // ( std::vector<std::vector<double>> )
      nt.resize(Nflavor,std::vector<double>(Ntau_p1,0.0));                      // ( std::vector<std::vector<double>> )
      nnt.resize(Nflavor*(Nflavor+1)/2,std::vector<double>(Ntau_p1,0.0));       // ( std::vector<std::vector<double>> )
      sign_meas=0.0;
      //
      // initialize seed - aggiungi qualcosa che dipenda dal task
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


   void solve(bool &atBeta)
   {
      //
      bool finalize=false;

      //
      if(!initialized)
      {
         std::cout << " Solver is not initialized - Exiting." << std::endl;
         exit(1);
      }
      else
      {
         std::cout << " Solver is about to start." << std::endl;
         //
         if((retarded==true)&&(atBeta==false))  dostep = &ct_hyb::dostep_retarded_full;
         if((retarded==true)&&(atBeta==true))   dostep = &ct_hyb::dostep_retarded_atBeta;
         if((retarded==false)&&(atBeta==false)) dostep = &ct_hyb::dostep_istantaneous_full;
         if((retarded==false)&&(atBeta==true))  dostep = &ct_hyb::dostep_istantaneous_atBeta;
         //
         do
         {
            dostep();

            //
            //accumulate_observables()

            //
            // Global time condition
            finalize = check_timer();

         }while( !finalize );

      }

   }


   //---------------------------------------------------------------------------


   void dostep_retarded_full()
   {
         //
         double s=1;

         // The measurments I'm going to do regardless from the time
         for (int imeas=0; imeas<Nmeas; imeas++)
         {
            //
            //Flavor loop
            for (int ifl=0; ifl<Nflavor; ifl++)
            {
                //
                // insert or remove full line
//                if (segments[ifl].size() == 0) insert_remove_full_line(rndm(), Eloc[ifl], Umat, Beta, full_line[ifl], segments, full_line, ifl);
//                insert_remove_antisegment(rndm(), Beta*rndm(), Ntau, Beta, Eloc[ifl], Umat, F[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);
                //
                if (!full_line[ifl])
                {
                   //
                   // local update
//                   insert_remove_segment(rndm(), Beta*rndm(), Ntau, Beta, Eloc[ifl], Umat, F[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);
                   //
                   // shift segments
//                   for (int k=0; k<N_shift; k++) shift_segment(rndm(), segments[ifl], Ntau, Beta, Eloc[ifl], Umat, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);
                   //
                   // flip segment
                   //for (int i=0; i<N_flip; i++)flip_segment(rndm(), segments_up, Ntau, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
                }
                //
                // measure perturbation order
                if (segments[ifl].size()<Norder) Pert[ifl*Norder+segments[ifl].size()] += 1;
                //
                // measure Green's functions
                if (segments[ifl].size()>0) G[ifl] = measure_G( segments[ifl], M[ifl], Ntau_p1, Beta );
                //
                // collect the sign among the segments
                s *= sign[ifl];
                //
                // Measure the flavour occupation
//                Nloc[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/Beta;
            } // end Flavor loop
            //
            //
            sign_meas += s;
         }  // end Nmeas loop

         //
         // measure n_a(\tau)
         nt = measure_nt( segments, full_line, Ntau_p1, Beta );

         //
         // measure n_a(\tau)n_b(0)
         nnt = measure_nnt( nt );

         //
         // measure density histogram
         Nhist = measure_Nhist( nt );

         //
         // measure spin histogram
         Szhist = measure_Szhist( nt );

   }


   //---------------------------------------------------------------------------


   void dostep_retarded_atBeta()
   {
      //
      if(!initialized)
      {
         std::cout << " Solver is not initialized - Exiting." << std::endl;
         exit(1);
      }

      //
      bool finalize=false;

   }


   //---------------------------------------------------------------------------


   void dostep_istantaneous_full()
   {
      //
      if(!initialized)
      {
         std::cout << " Solver is not initialized - Exiting." << std::endl;
         exit(1);
      }

      //
      bool finalize=false;

   }


   //---------------------------------------------------------------------------


   void dostep_istantaneous_atBeta()
   {
      //
      if(!initialized)
      {
         std::cout << " Solver is not initialized - Exiting." << std::endl;
         exit(1);
      }

      //
      bool finalize=false;

   }


};


#endif
