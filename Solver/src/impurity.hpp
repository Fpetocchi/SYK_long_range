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
#include "segments.hpp"
#include "file_io.hpp"
#include "observables.hpp"
#include "moves.hpp"


//============================================================================//


class ct_hyb
{
public:
   ct_hyb(double mu, double beta, int Nspin, int Norb, int Ntau, int Norder,
          int Nmeas, int Nshift, int PrintTime, int binlength, bool retarded ):
   mu(mu), Beta(beta),Nspin(Nspin), Norb(Norb), Ntau(Ntau), Norder(Norder),
   Nmeas(Nmeas), Nshift(Nshift), PrintTime(PrintTime),
   binlength(binlength),retarded(retarded)
   {}
   //
   double get_mu(void)const {return mu;}
   double get_Beta(void)const {return Beta;}
   int get_Nspin(void)const {return Nspin;}
   int get_Norb(void)const {return Norb;}
   int get_Nflavor(void)const {return Nspin*Norb;}
   int get_Ntau(void)const {return Ntau;}
   int get_Nmeas(void)const {return Nmeas;}
   Vec get_Nloc(void)const {return Nloc;}
   double get_Density(void)const {return std::accumulate(Nloc.begin(), Nloc.end(), 0.0);}
   void reset_mu(double mu_new) {mu = mu_new; mu_is_reset=true;}


   //---------------------------------------------------------------------------


   void init(path &inputDir)
   {
      //
      // This can be made optional
      resultsDir = inputDir;

      //
      std::vector<path> mandatoryFiles{"/Eloc.DAT", "/Delta.DAT", "/Uloc.DAT"};
      if(retarded==true) mandatoryFiles.push_back("/K.DAT");
      for(int ifile=0; ifile < (int)mandatoryFiles.size(); ifile++)
      {
         path filepath=inputDir+mandatoryFiles[ifile];
         if(!PathExist(strcpy(new char[filepath.length() + 1], filepath.c_str())))
         {
            std::cout << filepath+" (Not Found) - Exiting." << std::endl;
            exit(1);
         }
      }

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
      if(retarded==true) read_VecVecVec(inputDir+"/K.DAT", K_table, Norb, Ntau_p1); // read_VecVecVec_extended(inputDir+"/K.DAT", K_table, Nflavor, Nflavor, Ntau_p1);

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
      Gerr.resize(Nflavor,std::vector<double>(Ntau_p1,0.0));                    // ( std::vector<std::vector<double>> )
      nt.resize(Nflavor,std::vector<double>(Ntau_p1,0.0));                      // ( std::vector<std::vector<double>> )
      nnt.resize(Nflavor*(Nflavor+1)/2,std::vector<double>(Ntau_p1,0.0));       // ( std::vector<std::vector<double>> )
      sign_meas=0.0;
      sweeps=0;

      //
      K_table_ptr = NULL;
      if(retarded==true) K_table_ptr = &K_table;

      //
      // initialize seed - aggiungi qualcosa che dipenda dal task
      unsigned long seed = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::system_clock::now().time_since_epoch()).count()-1566563000000;
      generator.seed(abs(seed));

      //
      //start the timer for the calculation
      start_timer();

      //
      initialized=true;
      std::cout << " Solver is initialized and timer has started." << std::endl;
   }


   //---------------------------------------------------------------------------


   void solve(int MaxTime, bool atBeta=false)
   {
      //
      if(!initialized)
      {
         std::cout << " Solver is not initialized - Exiting." << std::endl;
         exit(1);
      }
      else
      {
         //
         if(atBeta==false) dostep = &ct_hyb::dostep_full;
         if(atBeta==true)  dostep = &ct_hyb::dostep_atBeta;

         //
         std::cout << " Solver is about to start." << std::endl;

         //
         while( !check_timer(MaxTime) )
         {
            //
            (this->*dostep)();
            sweeps++;

            //
            if( check_timer(PrintTime) )
            {
               std::cout << " Printing observables." << std::endl;
               std::cout << " Partial sweeps: " << sweeps << std::endl;
               std::cout << " Average sign: " << sign_meas/sweeps << std::endl;
               print_observables(sweeps);
            }

         }

         //
         std::cout << " Solver is done." << std::endl;
         std::cout << " Printing observables." << std::endl;
         std::cout << " Total sweeps: " << sweeps << std::endl;
         std::cout << " Average sign: " << sign_meas/sweeps << std::endl;
         print_observables(sweeps);

      }

   }

private:
   // Pointer to the internal function defining the solution mode
   void(ct_hyb::*dostep)(void) = NULL;
   VecVecVec *K_table_ptr = NULL;
   // Input vars
   double                              mu;                                      // chemical potential
   double                              Beta;                                    // inverse temperature
   int                                 Nspin;
   int                                 Norb;
   int                                 Ntau;
   int                                 Norder;                                  // Max perturbation order
   int                                 Therm;                                   // Thernalization sweeps
   int                                 Nmeas;                                   // Number of sweeps between expensive measurments
   int                                 Nshift;
   int                                 PrintTime;
   int                                 binlength;
   bool                                retarded;                                // flag to switch on retarded interactions
   // Internal vars
   path                                resultsDir;
   duration                            Tstart,Tnow;
   bool                                mu_is_reset=false;
   bool                                initialized=false;
   int                                 Nflavor=Nspin*Norb;                      // Spin-orbital flavours
   int                                 Ntau_p1=Ntau+1;                          // Spin-orbital flavours
   unsigned long long int              sweeps;                                  // Sweeps done
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
   VecVec                              Gerr;                                    // Impurity Green's function
   VecVec                              nt;                                      // n(\tau)
   VecVec                              nnt;                                     // Charge susceptibility


   //---------------------------------------------------------------------------


   void start_timer(void)
   {
      Tstart = std::chrono::system_clock::now();
   }

   inline bool check_timer(int MaxTime)
   {
      bool stop=false;
      Tnow = std::chrono::system_clock::now();
      std::chrono::duration<double> runtime_seconds = Tnow-Tstart;
      if(runtime_seconds.count() > MaxTime) stop=true;
      return stop;
   }


   //---------------------------------------------------------------------------


   void dostep_full()
   {
      //
      times full_segment(0,Beta);
      double s=1;
      VecVec G_tmp(Nflavor,Vec(Ntau_p1,0.0));
      VecVec Gerr_tmp(Nflavor,Vec(Ntau_p1,0.0));

      // The measurments I'm going to do regardless from the time
      for (int imeas=0; imeas<Nmeas; imeas++)
      {
         for (int ifl=0; ifl<Nflavor; ifl++)
         {

            // insert or remove full line
            if (segments[ifl].size() == 0) insert_remove_full_line( Eloc[ifl], Uloc, Beta, full_line[ifl], segments, full_line, ifl );
            insert_remove_antisegment( Beta*rndm(), Beta, Eloc[ifl], Uloc, F[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table_ptr );

            //
            if (!full_line[ifl])
            {
               // local update
               insert_remove_segment( Beta*rndm(), Beta, Eloc[ifl], Uloc, F[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table_ptr);
               // shift segments
               for (int k=0; k<Nshift; k++) shift_segment( segments[ifl], Beta, Eloc[ifl], Uloc, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table_ptr);
               // flip segment
               //for (int i=0; i<N_flip; i++)flip_segment( segments_up, Ntau, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
            }

            //.........................Cheap measurments........................
            // perturbation order
            if (segments[ifl].size()<(unsigned int)Norder) Pert[ifl*Norder+segments[ifl].size()] += (1./Nmeas);
            // Green's functions
            if (segments[ifl].size()>0)  G_tmp[ifl] = measure_G( segments[ifl], M[ifl], Ntau_p1, Beta, (double)(Ntau/Nmeas)); //accumulate_G( G[ifl], segments[ifl], M[ifl], Ntau_p1, Beta, (double)(Ntau/Nmeas));
            // sign among the segments
            s *= sign[ifl];
            // flavour occupation
            Nloc[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/(Beta*Nmeas);

         }

         //
         sign_meas += s/Nmeas;

      }

      //.........................Expensive measurments..........................
      // n_a(\tau)
      nt = measure_nt( segments, full_line, Ntau_p1, Beta );
      // correct G(0^+)and G(beta^-) and estimate the error
      correct_G( Nloc, binlength, G_tmp, Gerr_tmp );
      if(binlength>0)accumulate_G( Gerr, Gerr_tmp );
      accumulate_G( G, G_tmp );
      // n_a(\tau)n_b(0)
      accumulate_nnt( nnt, nt );
      // density histogram
      accumulate_Nhist( Nhist, nt );
      // spin histogram
      accumulate_Szhist( Szhist, nt );

   }


   //---------------------------------------------------------------------------


   void dostep_atBeta()
   {
      //
      times full_segment(0,Beta);

      // The measurments I'm going to do regardless from the time
      for (int imeas=0; imeas<Nmeas; imeas++)
      {
         //
         //Flavor loop
         for (int ifl=0; ifl<Nflavor; ifl++)
         {

            // insert or remove full line
            if (segments[ifl].size() == 0) insert_remove_full_line( Eloc[ifl], Uloc, Beta, full_line[ifl], segments, full_line, ifl );
            insert_remove_antisegment( Beta*rndm(), Beta, Eloc[ifl], Uloc, F[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table_ptr );

            //
            if (!full_line[ifl])
            {
               // local update
               insert_remove_segment( Beta*rndm(), Beta, Eloc[ifl], Uloc, F[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table_ptr);
               // shift segments
               for (int k=0; k<Nshift; k++) shift_segment( segments[ifl], Beta, Eloc[ifl], Uloc, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table_ptr);
               // flip segment
               //for (int i=0; i<N_flip; i++)flip_segment( segments_up, Ntau, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
            }

            //.........................Cheap measurments........................
            // flavour occupation
            Nloc[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/(Beta*Nmeas);

         }
      }
   }


   //---------------------------------------------------------------------------


   void print_observables(unsigned long long int &iterations)
   {
      print_Vec(resultsDir+"/Nloc.DAT", Nloc, iterations);
      print_Vec(resultsDir+"/PertOrder.DAT", Pert, iterations);
      print_VecVec(resultsDir+"/Gimp.DAT", G, iterations);
      print_VecVec(resultsDir+"/nnt.DAT", nnt, iterations);
   }


};


#endif
