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
#include <cstdlib>
#include <mpi.h>
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
#include "mpi.hpp"


//============================================================================//


class ct_hyb
{
   public:

      //----------------------------------------------------------------------//

      ct_hyb( path SiteName, double beta, int Nspin, int Norb, int NtauF, int NtauB,
              int Norder, int Nmeas, int Ntherm, int Nshift,
              bool paramagnet, bool retarded, bool nnt_meas,
              int printTime, std::vector<int> bins,
              CustomMPI &mpi, bool testing_mpi=false):
      SiteName(SiteName), Beta(beta), Nspin(Nspin), Norb(Norb), NtauF(NtauF), NtauB(NtauB),
      Norder(Norder), Nmeas(Nmeas), Ntherm(Ntherm), Nshift(Nshift),
      paramagnet(paramagnet), retarded(retarded), nnt_meas(nnt_meas),
      printTime(printTime), bins(bins), mpi(mpi), testing_mpi(testing_mpi)
      {}

      //----------------------------------------------------------------------//

      void reset_mu(double mu_new) {mu = mu_new; mu_is_reset=true;}
      double get_mu(void) {return mu;}

      Vec get_Nloc(bool EachRank=false)
      {
         //
         Vec NormNloc(Nflavor,0.0);
         for (int ifl=0; ifl<Nflavor; ifl++) NormNloc[ifl] = Nloc[ifl] / RankSweeps;

         //
         if(EachRank)
         {
            // return the spin-orbital occupation of the rank
            return NormNloc;
         }
         else
         {
            // return the spin-orbital occupation averaged over the ranks
            Vec WorldNloc(Nflavor,0.0);
            mpi.allreduce(NormNloc, WorldNloc, true);
            return WorldNloc;
         }
      }

      double get_Density(bool EachRank=false)
      {
         Vec WorldNloc = get_Nloc(EachRank);
         return std::accumulate(WorldNloc.begin(), WorldNloc.end(), 0.0);
      }

      //----------------------------------------------------------------------//

      void init(path &inputDir)
      {
         //
         // This can be made optional or changed
         resultsDir = inputDir+"/resultsQMC/";
         if(!PathExist(strcpy(new char[resultsDir.length() + 1], resultsDir.c_str())) && mpi.is_master())
         {
            int check = mkdir(resultsDir.c_str(),0777);
            mpi.report(" Folder = "+resultsDir+" Created. "+str(check));
         }

         //
         // Check if mandatory files are present
         std::vector<path> mandatoryFiles{"/Eloc.DAT", "/Delta_t.DAT", "/Umat.DAT"};
         if(retarded) mandatoryFiles.push_back("/K_t.DAT");
         for(int ifile=0; ifile < (int)mandatoryFiles.size(); ifile++)
         {
            path filepath=inputDir+mandatoryFiles[ifile];
            if(!PathExist(strcpy(new char[filepath.length() + 1], filepath.c_str())))
            {
               mpi.StopError(filepath+" (Not Found) - Exiting.");
            }
         }

         //
         //read the istantaneous interaction ( Eigen::MatrixXd )
         read_EigenMat(inputDir+"/Umat.DAT", Uloc, Nflavor, Nflavor);
         mu_correction.resize(Norb);
         for(int iorb=0; iorb < Norb; iorb++)
         {
            mu_correction[iorb] += Uloc(2*iorb,2*iorb+1)/2.0;
            mpi.report(" Orbital "+str(iorb)+" correction of the local level: "+str(mu_correction[iorb],4));
         }

         //
         //set the local energies ( std::vector<double> )
         read_Vec(inputDir+"/Eloc.DAT", Eloc, Nflavor, mu);
         set_levels();
         mpi.report(" Chemical potential: "+str(mu,4));

         //
         //read the hybridization function ( std::vector<std::vector<double>> )
         read_VecVec(inputDir+"/Delta_t.DAT", F, Nflavor, NtauF, true, true);  // last flag is to reverse the tau index

         //
         //read the screening function ( std::vector<std::vector<std::vector<double>>> )
         if(retarded)
         {
            read_VecVecVec(inputDir+"/K_t.DAT", K_table, Nflavor, true); // read_VecVecVec(inputDir+"/K.DAT", K_table, Norb, Ntau, true);
            //
            for(int ifl=0; ifl < Norb; ifl++)
            {
               for(int jfl=0; jfl <= ifl; jfl++)
               {
                  int Ntau_K = K_table[ifl][jfl].size();
                  path Kcomp = "K["+str(ifl)+"]["+str(jfl)+"]";
                  if(Ntau_K!=NtauB) mpi.report(" The Number of tau points in "+Kcomp+" is: "+str(Ntau_K)+" different from NtauB: "+str(NtauB));
                  //
                  double K_zero =  K_table[ifl][jfl].front();
                  double K_beta =  K_table[ifl][jfl].back();
                  if(K_zero!=0.0) mpi.StopError( "->"+Kcomp+" at tau=0 is not vanishing - Exiting.");
                  if(K_beta!=0.0) mpi.StopError( "->"+Kcomp+" at tau=beta is not vanishing - Exiting.");
               }
            }
         }

         //
         if(testing_mpi && mpi.is_master())
         {
            print_Vec(inputDir+"/used.Eloc.DAT", Eloc, mu);
            print_Vec(inputDir+"/used.Levels.DAT", Levels);
            print_VecVec(inputDir+"/used.Delta_t.DAT", F);
            if(retarded)print_VecVecVec(inputDir+"/used.K_t.DAT", K_table);
         }

         //
         //initialize segment container ( std::vector<std::set<times>> )
         segments.resize(Nflavor);
         for (int ifl=0; ifl<Nflavor; ifl++)segments[ifl].clear();

         //
         //initialize full_line ( std::vector<int> )
         full_line.resize(Nflavor,0);

         //
         //initialize empty Inverse hybridization matrix M ( std::vector<Eigen::MatrixXd> )
         M.resize(Nflavor);

         //
         //initialize observables
         sign_meas=0.0;
         sign.resize(Nflavor,1.0);                                              // ( std::vector<double> )
         Nloc.resize(Nflavor,0.0);                                              // ( std::vector<double> )
         Nhist.resize(Nflavor+1,0.0);                                           // ( std::vector<double> )
         Szhist.resize(Nflavor/2+1,0.0);                                        // ( std::vector<double> )
         Pert.resize(Nflavor,std::vector<double>(Norder,0.0));                  // ( std::vector<std::vector<double>> )
         G.resize(Nflavor,std::vector<double>(NtauF,0.0));                      // ( std::vector<std::vector<double>> )
         Gerr.resize(Nflavor,std::vector<double>(NtauF,0.0));                   // ( std::vector<std::vector<double>> )
         nt.resize(Nflavor,std::vector<double>(NtauB,0.0));                     // ( std::vector<std::vector<double>> )
         nnt.resize(Nflavor*(Nflavor+1)/2,std::vector<double>(NtauB,0.0));      // ( std::vector<std::vector<double>> )
         //
         RankSign=0.0;
         WorldSign=0.0;
         RankSweeps=0;
         WorldSweeps=0;

         //
         // initialize seed
         //unsigned long seed = std::chrono::duration_cast<std::chrono::minutes>(std::chrono::system_clock::now().time_since_epoch()).count()-4321*mpi.rank();
         unsigned long seed = random_seed();
         srand((unsigned)time(NULL)+mpi.rank()*mpi.size() + 7);

         //
         generator.seed(fabs(seed));
         if(testing_mpi) mpi.report(" Seed: "+str(seed));
         //
         initialized=true;
         mpi.report(" Solver is initialized.");

      }

      //----------------------------------------------------------------------//

      void solve(double MaxTime, bool atBeta=false)
      {
         //
         if(!initialized)
         {
            mpi.StopError(" Solver is not initialized - Exiting.");
         }
         else
         {
            // Execution mode
            if(!atBeta ) dostep = &ct_hyb::dostep_full;
            if( atBeta ) dostep = &ct_hyb::dostep_atBeta;

            //
            if(mu_is_reset)
            {
               // Internal vars reset
               set_levels();
               for (int ifl=0; ifl<Nflavor; ifl++)
               {
                  segments[ifl].clear();
                  M[ifl].resize(0,0);
                  Nloc[ifl]=0.0;
                  sign[ifl]=1.0;
               }
               full_line.resize(Nflavor,0);
               RankSign=0.0;
               WorldSign=0.0;
               RankSweeps=0;
               WorldSweeps=0;
               mpi.barrier();
            }

            //thermalization pre-steps - I'm using always dostep_atBeta either or not for quick
            dostep_atBeta(((atBeta == true) ? 10 : Ntherm)*Nmeas,true);
            mpi.report(" Thermalization is done.");

            //
            start_timer();
            if(!atBeta)
            {
               print_line_space(1,mpi.is_master());
               print_line_equal(80,mpi.is_master());
               mpi.report(" Solver for impurity "+SiteName+" has started.");
            }

            //full solution
            while( !check_timer_global(MaxTime) )
            {
               //
               (this->*dostep)(Nmeas,false);
               //
               RankSweeps++;
               mpi.allreduce(RankSweeps,WorldSweeps);
               ////printf(" [Rank %d] - RankSweeps: %d - WorldSweeps: %d \n",mpi.rank(), RankSweeps,WorldSweeps);
               //
               RankSign = sign_meas / RankSweeps;
               mpi.allreduce(RankSign,WorldSign,true);
               //
               if( check_timer_print(printTime) && (!atBeta) )
               {
                  print_line_space(1,mpi.is_master());
                  print_timestamp("min");
                  mpi.report(" Partial sweeps: "+str( (testing_mpi == true) ? RankSweeps : WorldSweeps ));
                  mpi.report(" Average sign: "+str( (testing_mpi == true) ? RankSign : WorldSign ));
                  mpi.report(" Density: "+str(get_Density(testing_mpi)));
                  //
                  path pad = "_T"+str(printTime*(TimeStamp-1))+".DAT";
                  print_observables(pad,bins);
                  //
               }
            }
            mu_is_reset=false;

            // testing_mpi>>>
            //mpi.report(" RankSweeps: "+str(RankSweeps),true);
            //mpi.report(" Dens: "+str(std::accumulate(Nloc.begin(), Nloc.end(), 0.0)/RankSweeps),true);
            // >>>testing_mpi

            //
            if(!atBeta)
            {
               print_line_space(2,mpi.is_master());
               print_line_minus(80,mpi.is_master());
               print_timestamp("min");
               mpi.report(" Solver for impurity "+SiteName+" is done. Results written to: "+resultsDir);
               mpi.report(" Total sweeps: "+str( (testing_mpi == true) ? RankSweeps : WorldSweeps ));
               mpi.report(" Average sign: "+str( (testing_mpi == true) ? RankSign : WorldSign ));
               mpi.report(" Density: "+str(get_Density(testing_mpi)));
               //
               path pad = ".DAT";
               print_observables(pad,bins);
               print_line_minus(80,mpi.is_master());
            }
            else
            {
               print_timestamp("sec");
               mpi.report(" Total sweeps: "+str( (testing_mpi == true) ? RankSweeps : WorldSweeps ));
               mpi.report(" Solver for impurity "+SiteName+" is done.");
            }
            //
            mpi.barrier();
         }
      }

      //----------------------------------------------------------------------//

   private:

      //----------------------------------------------------------------------//
      // Input vars
      path                                SiteName;
      double                              Beta;                                 // inverse temperature
      int                                 Nspin;
      int                                 Norb;
      int                                 NtauF;                                // Number of POINTS in Fermionic tau grid
      int                                 NtauB;                                // Number of POINTS in Bosonic tau grid
      int                                 Norder;                               // Max perturbation order
      int                                 Nmeas;                                // Number of sweeps between expensive measurments
      int                                 Ntherm;
      int                                 Nshift;
      bool                                paramagnet;
      bool                                retarded;
      bool                                nnt_meas;
      int                                 printTime;
      std::vector<int>                    bins;                                 // 2D vec Contains binlength and binstart in [0] and [1] respectively
      CustomMPI                           mpi;
      bool                                testing_mpi;
      // Internal vars
      void(ct_hyb::*dostep)(int,bool) = NULL;                                   // Pointer to the internal function defining the solution mode
      path                                resultsDir;
      double                              mu;                                   // chemical potential
      bool                                mu_is_reset=false;
      bool                                initialized=false;
      duration                            Tstart;
      int                                 TimeStamp;
      int                                 Nflavor=Nspin*Norb;                   // Spin-orbital flavours
      int                                 NtauF_m1=NtauF-1;                     // Number of SEGMENTS in Fermionic tau grid
      int                                 NtauB_m1=NtauB-1;                     // Number of SEGMENTS in Bosonic tau grid
      unsigned long long int              RankSweeps,WorldSweeps;               // Sweeps done
      // Input data
      Vec                                 Levels;                               // mu-<\epsilon>
      Vec                                 mu_correction;                        // chemical potential correction due to the interaction shift
      Vec                                 Eloc;                                 // <\epsilon>
      Mat                                 Uloc;                                 // Istantaneous U matrix
      VecVec                              F;                                    // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu)
      VecVecVec                           K_table;                              // K function matrix for retarded interactions
      // Solver vars
      VecMat                              M;                                    // Inverse hybridization matrix
      std::vector<int>                    full_line;                            // if 1 means that particle occupies full time-line
      std::vector<segment_container_t>    segments;                             // Stores configurations with 0,1,... segments (but not full line)
      // Observables
      double                              sign_meas;
      double                              RankSign,WorldSign;                   // Total Sign
      Vec                                 sign;                                 // Sign per flavor
      Vec                                 Nloc;                                 // Density per flavor
      VecVec                              Pert;                                 // Perturbation order
      Vec                                 Nhist;                                //
      Vec                                 Szhist;                               //
      VecVec                              G;                                    // Impurity Green's function
      VecVec                              Gerr;                                 // Impurity Green's function
      VecVec                              nt;                                   // n(\tau)
      VecVec                              nnt;                                  // Charge susceptibility

      //----------------------------------------------------------------------//

      void start_timer(void)
      {
         TimeStamp=1;
         Tstart = std::chrono::system_clock::now();
      }

      inline bool check_timer_global(double DeltaMin)
      {
         bool stop=false;
         duration Tnow = std::chrono::system_clock::now();
         std::chrono::duration<double> runtime_seconds = Tnow-Tstart;
         if(runtime_seconds.count() > 60*DeltaMin) stop=true;
         mpi.broadcast(stop,mpi.master());
         return stop;
      }

      inline bool check_timer_print(double DeltaMin)
      {
         bool stop=false;
         duration Tnow = std::chrono::system_clock::now();
         std::chrono::duration<double> runtime_seconds = Tnow-Tstart;
         if(runtime_seconds.count() > 60*(DeltaMin*TimeStamp)) stop=true;
         mpi.broadcast(stop,mpi.master());
         if(stop) TimeStamp++;
         return stop;
      }

      void print_timestamp(path kind="sec")
      {
         duration Tnow = std::chrono::system_clock::now();
         std::chrono::duration<double> runtime_seconds = Tnow-Tstart;
         if(kind=="sec")mpi.report(" Timestamp(sec): "+str(runtime_seconds.count(),3));
         if(kind=="min")mpi.report(" Timestamp(min): "+str(runtime_seconds.count()/60,3));
         if(kind=="hrs")mpi.report(" Timestamp(hrs): "+str(runtime_seconds.count()/3600,3));
      }

      //----------------------------------------------------------------------//

      void dostep_full(int Nmeas_, bool fulldry)
      {
         //
         times full_segment(0,Beta);
         double s=1;
         VecVec G_tmp(Nflavor,Vec(NtauF,0.0));
         VecVec n_tau(Nflavor,Vec(NtauB,0.0));

         // The measurments I'm going to do regardless from the time
         for (int imeas=0; imeas<Nmeas_; imeas++)
         {
            for (int ifl=0; ifl<Nflavor; ifl++)
            {
               // insert or remove full line
               if (segments[ifl].size() == 0) insert_remove_full_line( Levels[ifl], Uloc, Beta, full_line[ifl], segments, full_line, ifl );
               insert_remove_antisegment( Beta*rndm(), Beta, Levels[ifl], Uloc, F[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );

               //
               if (!full_line[ifl])
               {
                  // local update
                  insert_remove_segment( Beta*rndm(), Beta, Levels[ifl], Uloc, F[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);
                  // shift segments
                  for (int k=0; k<Nshift; k++) shift_segment( segments[ifl], Beta, Levels[ifl], Uloc, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );
                  // flip segment
                  //for (int i=0; i<N_flip; i++)flip_segment( segments_up, NtauF_m1, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
               }

               //
               //.........................Cheap measurments..................... add perturbation order
               // Green's functions
               if (segments[ifl].size()>0) measure_G( G_tmp[ifl], segments[ifl], M[ifl], NtauF, Beta );
               // sign among the segments
               s *= sign[ifl];
               // flavour occupation
               for (int i=0; i<G_tmp[ifl].size(); i++)G_tmp[ifl][i]*=(1.*NtauF_m1)/Nmeas_;
               Nloc[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/(Beta*Nmeas_);
               //...............................................................
            }
            //
            sign_meas += s/(double)Nmeas_;
         }

         //
         //.........................Expensive measurments.......................
         //
         if(paramagnet) spin_symm(G_tmp);
         accumulate_G( G, G_tmp );
         //
         if(nnt_meas)
         {
            n_tau = measure_nt( segments, full_line, NtauB, Beta );
            accumulate_nt( nt, n_tau );
            accumulate_nnt( nnt, n_tau );
            accumulate_Szhist( Szhist, n_tau );
            accumulate_Nhist( Nhist, n_tau );
         }
         //.....................................................................

      }

      //----------------------------------------------------------------------//

      void dostep_atBeta(int Nmeas_, bool fulldry)
      {
         //
         times full_segment(0,Beta);

         // The measurments I'm going to do regardless from the time
         for (int imeas=0; imeas<Nmeas_; imeas++)
         {
            for (int ifl=0; ifl<Nflavor; ifl++)
            {
               // insert or remove full line
               if (segments[ifl].size() == 0) insert_remove_full_line( Levels[ifl], Uloc, Beta, full_line[ifl], segments, full_line, ifl );
               insert_remove_antisegment( Beta*rndm(), Beta, Levels[ifl], Uloc, F[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );

               //
               if (!full_line[ifl])
               {
                  // local update
                  insert_remove_segment( Beta*rndm(), Beta, Levels[ifl], Uloc, F[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );
                  // shift segments
                  for (int k=0; k<Nshift; k++) shift_segment( segments[ifl], Beta, Levels[ifl], Uloc, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );
                  // flip segment
                  //for (int i=0; i<N_flip; i++)flip_segment( segments_up, NtauF_m1, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
               }

               //.........................Cheap measurments.....................
               // flavour occupation
               if(!fulldry) Nloc[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/(Beta*Nmeas_);
            }
         }
      }

      //----------------------------------------------------------------------//

      void set_levels()
      {
         Levels.resize(Nflavor,0.0);
         for (int ifl=0; ifl<Nflavor; ifl++) Levels[ifl] = mu - Eloc[ifl] + mu_correction[(int)(ifl/2)];
      }

      //----------------------------------------------------------------------//

      void print_observables(path pad, std::vector<int> bins)
      {
         if(testing_mpi) // All ranks print their own observables
         {
            //
            mpi.report(" Rank #"+str(mpi.rank())+" is printing observables.");
            //
            // density
            Vec PrintNloc = get_Nloc(testing_mpi); // provides spin-orbital occupation already normalized
            print_Vec(resultsDir+"/Nqmc_rank"+str(mpi.rank())+pad, PrintNloc, mu);
            mpi.report(" Nqmc_rank"+str(mpi.rank())+pad+" is printed.");
            //
            // error estimate and binning
            if(bins[0]>0)
            {
               VecVec Gerr(Nflavor,Vec(NtauF,0.0));
               binAverageVec( bins, G, Gerr );
               print_VecVec(resultsDir+"/Gerr_rank"+str(mpi.rank())+pad, Gerr, Beta, (double)RankSweeps);
               mpi.report(" Gerr_rank"+str(mpi.rank())+pad+" is printed.");
            }
            //
            // Green's function
            Vec Ntmp = Nloc;
            if(paramagnet) spin_symm(Ntmp);
            for (int ifl=0; ifl<Nflavor; ifl++)
            {
               G[ifl].front() = +Ntmp[ifl]-(long double)RankSweeps ;
               G[ifl].back()  = -Ntmp[ifl];
            }
            print_VecVec(resultsDir+"/Gimp_t_rank"+str(mpi.rank())+pad, G, Beta, (double)RankSweeps);
            mpi.report(" Gimp_t_rank"+str(mpi.rank())+pad+" is printed.");
            //
            // observables derived from n(tau)n(0)
            if(nnt_meas)
            {
               //
               print_Vec(resultsDir+"/Nhist_rank"+str(mpi.rank())+pad, Nhist );
               mpi.report(" Nhist_rank"+str(mpi.rank())+pad+" is printed.");
               //
               print_Vec(resultsDir+"/Szhist_rank"+str(mpi.rank())+pad, Szhist );
               mpi.report(" Szhist_rank"+str(mpi.rank())+pad+" is printed.");

               print_VecVec(resultsDir+"/n_t_rank"+str(mpi.rank())+pad, nt, Beta, (double)RankSweeps);
               mpi.report(" n_t_rank"+str(mpi.rank())+pad+" is printed.");
               //
               print_VecVec(resultsDir+"/nn_t_rank"+str(mpi.rank())+pad, nnt, Beta, (double)RankSweeps);
               mpi.report(" nn_t_rank"+str(mpi.rank())+pad+" is printed.");
            }
         }

         //
         // In any case at the ened Master prints after the collapse of the observables
         //
         mpi.report(" Master (Rank #"+str(mpi.master())+") is printing observables.");
         //
         // density
         Vec PrintNloc = get_Nloc(); // provides spin-orbital occupation already normalized
         if(mpi.is_master()) print_Vec(resultsDir+"/Nqmc"+pad, PrintNloc, mu);
         mpi.report(" Nqmc"+pad+" is printed.");
         //
         // error estimate and binning
         if(bins[0]>0)
         {
            VecVec Gerr(Nflavor,Vec(NtauF,0.0));
            binAverageVec( bins, G, Gerr );
            Gerr = normalize_VecVec(Gerr, RankSweeps);
            VecVec PrintGerr(Nflavor,std::vector<double>(NtauF,0.0));
            mpi.allreduce(Gerr, PrintGerr, true);
            if(mpi.is_master()) print_VecVec(resultsDir+"/Gerr"+pad, PrintGerr, Beta);
            mpi.report(" Gerr"+pad+" is printed.");
         }
         //
         // Green's function
         Vec Ntmp = Nloc;
         if(paramagnet) spin_symm(Ntmp);
         for (int ifl=0; ifl<Nflavor; ifl++)
         {
            G[ifl].front() = +Ntmp[ifl]-(long double)RankSweeps ;
            G[ifl].back()  = -Ntmp[ifl];
         }
         VecVec NormG = normalize_VecVec(G, RankSweeps);
         VecVec PrintG(Nflavor,std::vector<double>(NtauF,0.0));
         mpi.allreduce(NormG, PrintG, true);
         if(mpi.is_master()) print_VecVec(resultsDir+"/Gimp_t"+pad, PrintG, Beta);
         mpi.report(" Gimp_t"+pad+" is printed.");
         //
         // observables derived from n(tau)n(0)
         if(nnt_meas)
         {
            //
            Vec NormNhist = normalize_Vec(Nhist, RankSweeps);
            Vec PrintNhist(Nflavor+1,0.0);
            mpi.allreduce(NormNhist, PrintNhist, true);
            if(mpi.is_master()) print_Vec(resultsDir+"/Nhist"+pad, PrintNhist);
            mpi.report(" Nhist"+pad+" is printed.");
            //
            Vec NormSzhist = normalize_Vec(Szhist, RankSweeps);
            Vec PrintSzhist(Nflavor/2+1,0.0);
            mpi.allreduce(NormSzhist, PrintSzhist, true);
            if(mpi.is_master()) print_Vec(resultsDir+"/Szhist"+pad, PrintSzhist);
            mpi.report(" Szhist"+pad+" is printed.");
            //
            VecVec Normnt = normalize_VecVec(nt, RankSweeps);
            VecVec Printnt(Nflavor,std::vector<double>(NtauB,0.0));
            mpi.allreduce(Normnt, Printnt, true);
            if(mpi.is_master()) print_VecVec(resultsDir+"/n_t"+pad, Printnt, Beta);
            mpi.report(" n_t"+pad+" is printed.");
            //
            VecVec Normnnt = normalize_VecVec(nnt, RankSweeps);
            VecVec Printnnt(Nflavor*(Nflavor+1)/2,std::vector<double>(NtauB,0.0));
            mpi.allreduce(Normnnt, Printnnt, true);
            if(mpi.is_master()) print_VecVec(resultsDir+"/nn_t"+pad, Printnnt, Beta);
            mpi.report(" nn_t"+pad+" is printed.");
         }
      }

      //----------------------------------------------------------------------//

};


//============================================================================//


#endif
