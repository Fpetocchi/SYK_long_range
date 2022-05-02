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

      ct_hyb( path inputDir, path SiteName, double beta, int Nspin, int Norb, int NtauF, int NtauB,
              int Norder, bool Gexp, int Nmeas, int Ntherm, int NsegShift, int NspinSwap, int NnntMeas,
              bool removeUhalf, bool paramagnet, bool retarded, std::vector<int> SetsNorb,
              int printTime, std::vector<int> bins, CustomMPI &mpi):
              inputDir(inputDir),
              SiteName(SiteName),
              Beta(beta),
              Nspin(Nspin),
              Norb(Norb),
              NtauF(NtauF),
              NtauB(NtauB),
              Norder(Norder),
              Gexp(Gexp),
              Nmeas(Nmeas),
              Ntherm(Ntherm),
              NsegShift(NsegShift),
              NspinSwap(NspinSwap),
              NnntMeas(NnntMeas),
              removeUhalf(removeUhalf),
              paramagnet(paramagnet),
              retarded(retarded),
              SetsNorb(SetsNorb),
              printTime(printTime),
              bins(bins),
              mpi(mpi)
              {}

      //----------------------------------------------------------------------//

      void reset_mu( double mu_new ) { mu = mu_new; mu_is_reset=true; }

      double get_mu(void) { return mu; }

      Vec get_Nloc( bool EachRank=false )
      {
         //
         Vec NormNloc(Nflavor,0.0);
         for( int ifl=0; ifl<Nflavor; ifl++ ) NormNloc[ifl] = Nloc[ifl] / RankSweeps;

         //
         if(EachRank)
         {
            // return the spin-orbital occupation of the rank
            return NormNloc;
         }
         else
         {
            // return the spin-orbital occupation averaged over the ranks
            Vec WorldNloc( Nflavor, 0.0 );
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

      void init(void)
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
         // Check if Orbital symmetrization is required
         OrbSym=false;
         if(SetsNorb.size()>0)OrbSym=true;

         //
         // Check if mandatory files are present
         std::vector<path> mandatoryFiles{"/Eloc.DAT", "/Delta_t.DAT", "/Umat.DAT"};
         if(retarded) mandatoryFiles.insert( mandatoryFiles.end(), { "/K_t.DAT", "/Screening.DAT" } );
         if(OrbSym) mandatoryFiles.push_back("/Eqv.DAT");
         for(int ifile=0; ifile < (int)mandatoryFiles.size(); ifile++)
         {
            path filepath=inputDir+mandatoryFiles[ifile];
            if(!PathExist(strcpy(new char[filepath.length() + 1], filepath.c_str())))
            {
               mpi.StopError(filepath+" (Not Found) - Exiting.");
            }
         }

         //
         //read the local energies
         read_Vec(inputDir+"/Eloc.DAT", Eloc, Nflavor, mu);
         mpi.report(" Chemical potential: "+str(mu,4));

         //
         //read the istantaneous interaction
         read_EigenMat(inputDir+"/Umat.DAT", Uloc, Nflavor, Nflavor);

         //
         //read the hybridization function
         read_VecVec(inputDir+"/Delta_t.DAT", F, Nflavor, NtauF, true, true);  // last flag is to reverse the tau index
         for(int ifl=0; ifl < Nflavor; ifl++)
         {
            for(int itau=0; itau < F[ifl].size(); itau++)
            {
               path Dcomp = "Delta["+str(ifl)+"]";
               if(F[ifl][itau]<0.0) mpi.StopError( " ->"+Dcomp+" at tau "+str(itau)+" is positive - Exiting.");
            }
         }

         //
         //read the screening function
         Screening_shift.resize(Norb,0.0);
         if(retarded)
         {
            //
            //Read the screening from file
            if(!removeUhalf)
            {
               mpi.report(" Reading screening file.");
               read_EigenMat(inputDir+"/Screening.DAT", Screening_Mat, Nflavor, Nflavor);
               for(int iorb=0; iorb < Norb; iorb++)
               {
                  // big number
                  //for(int ifl=0; ifl < Nflavor; ifl++) Screening_shift[iorb] += (2*iorb!=ifl) ? Screening_Mat(2*iorb,ifl)/2.0 : 0.0 ;
                  // small number
                  Screening_shift[iorb] += Screening_Mat(2*iorb,2*iorb)/2.0;
                  mpi.report(" Orbital "+str(iorb)+" - screening shift = S/2: "+str(Screening_shift[iorb],6));
               }
            }

            //
            read_VecVecVec(inputDir+"/K_t.DAT", K_table, Nflavor, true);
            for(int ifl=0; ifl < Nflavor; ifl++)
            {
               for(int jfl=0; jfl <= ifl; jfl++)
               {
                  int Ntau_K = K_table[ifl][jfl].size();
                  path Kcomp = "K["+str(ifl)+"]["+str(jfl)+"]";
                  if(Ntau_K!=NtauB) mpi.report(" The Number of tau points in "+Kcomp+" is: "+str(Ntau_K)+" different from NtauB: "+str(NtauB));
                  //
                  double K_zero =  K_table[ifl][jfl].front();
                  double K_beta =  K_table[ifl][jfl].back();
                  if(K_zero!=0.0) mpi.StopError( " ->"+Kcomp+" at tau=0 is not vanishing - Exiting.");
                  if(K_beta!=0.0) mpi.StopError( " ->"+Kcomp+" at tau=beta is not vanishing - Exiting.");
               }
            }
         }

         //
         //rescale the chemical potential rescaling with the half-filling one
         Hartree_shift.resize(Norb,0.0);
         if(removeUhalf)
         {
            for(int iorb=0; iorb < Norb; iorb++)
            {
               for(int ifl=0; ifl < Nflavor; ifl++) Hartree_shift[iorb] += (2*iorb!=ifl) ? ( Uloc(2*iorb,ifl) )/2.0 : 0.0 ;
               mpi.report(" Orbital "+str(iorb)+" - Hartree term = Uscr/2: "+str(Hartree_shift[iorb],6));
            }
         }

         //
         //set the levels for the impurity solver
         set_levels();

         //
         //read in the indexes contained to each equivalent set
         if(OrbSym)
         {
            read_list(inputDir+"/Eqv.DAT",SetsNorb,SetsOrbs);
            for(int iset=0; iset < SetsNorb.size(); iset++)
            {
               mpi.report(" "+SiteName+" Eqv set #"+str(iset+1)+" indexes:");
               for(int iorb=0; iorb < SetsNorb[iset]; iorb++) mpi.report(" "+str(SetsOrbs[iset][iorb]));
            }
         }

         //
         if(mpi.is_master()) //debug &&
         {
            print_Vec(inputDir+"/used.Eloc.DAT", Eloc, mu);
            print_VecVec(inputDir+"/used.Delta_t.DAT", F);
            if(retarded)print_VecVecVec(inputDir+"/used.K_t.DAT", K_table);
         }

         //
         //initialize segment container
         segments.resize(Nflavor);
         for (int ifl=0; ifl<Nflavor; ifl++)segments[ifl].clear();

         //
         //initialize full_line
         full_line.resize(Nflavor,0);

         //
         //initialize empty Inverse hybridization matrix M
         M.resize(Nflavor);

         //
         //initialize observables
         sign_meas=0.0;
         sign.resize( Nflavor, 1.0 );
         Order.resize( Nflavor, Vec( Norder, 0.0 ) );
         Nloc.resize( Nflavor, 0.0 );
         Nhist.resize( Nflavor+1, 0.0 );
         Szhist.resize( Nflavor/2+1, 0.0 );
         G.resize( Nflavor, Vec( NtauF, 0.0 ) );
         Gerr.resize( Nflavor, Vec( NtauF, 0.0 ) );
         nt.resize( Nflavor, Vec( NtauB, 0.0 ) );
         nnt.resize( Nflavor*(Nflavor+1)/2, Vec( NtauB, 0.0 ) );
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
         if(debug) mpi.report(" Seed: "+str(seed));
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
                  mpi.report(" Partial sweeps: "+str( (debug == true) ? RankSweeps : WorldSweeps ));
                  mpi.report(" Average sign: "+str( (debug == true) ? RankSign : WorldSign ));
                  mpi.report(" Density: "+str(get_Density(debug))+" mu: "+str(mu,6) );
                  //
                  path pad = "_T"+str(printTime*(TimeStamp-1))+".DAT";
                  print_observables(pad,bins);
                  //
               }
            }
            mu_is_reset=false;

            //
            if(!atBeta)
            {
               print_line_space(2,mpi.is_master());
               print_line_minus(80,mpi.is_master());
               print_timestamp("min");
               mpi.report(" Solver for impurity "+SiteName+" is done. Results written to: "+resultsDir);
               mpi.report(" Total sweeps: "+str( (debug == true) ? RankSweeps : WorldSweeps ));
               mpi.report(" Average sign: "+str( (debug == true) ? RankSign : WorldSign ));
               mpi.report(" Density: "+str(get_Density(debug))+" mu: "+str(mu,6) );
               //
               path pad = ".DAT";
               print_observables(pad,bins);
               print_line_minus(80,mpi.is_master());
            }
            else
            {
               print_timestamp("sec");
               mpi.report(" Total sweeps: "+str( (debug == true) ? RankSweeps : WorldSweeps ));
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
      path                                inputDir;
      path                                SiteName;
      double                              Beta;                                 // inverse temperature
      int                                 Nspin;
      int                                 Norb;
      int                                 NtauF;                                // Number of POINTS in Fermionic tau grid
      int                                 NtauB;                                // Number of POINTS in Bosonic tau grid
      int                                 Norder;                               // Max perturbation order
      bool                                Gexp;
      int                                 Nmeas;                                // Number of sweeps between expensive measurments
      int                                 Ntherm;
      int                                 NsegShift;
      int                                 NspinSwap;
      int                                 NnntMeas;
      bool                                removeUhalf;
      bool                                screenshift;
      bool                                paramagnet;
      bool                                retarded;
      int                                 printTime;
      std::vector<int>                    SetsNorb;
      std::vector<int>                    bins;                                 // 2D vec Contains binlength and binstart in [0] and [1] respectively
      CustomMPI                           mpi;
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
   #ifdef _verb
      bool debug=true;
   #else
      bool debug=false;
   #endif
      std::vector<std::vector<int>>       SetsOrbs;
      bool                                OrbSym;
      // Input data
      Vec                                 Levels;                               // mu - <\epsilon> + mu_correction
      Vec                                 Hartree_shift;                        // chemical potential shift
      Vec                                 Screening_shift;                      // chemical potential correction due to shift encoded in the screening function
      Vec                                 Eloc;                                 // <\epsilon>
      Mat                                 Uloc;                                 // Istantaneous U matrix
      Mat                                 Screening_Mat;                        // Screening matrix
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
      VecVec                              Order;                                // Perturbation order
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
         Vec N_tmp(Nflavor,0.0);
         VecVec G_tmp(Nflavor,Vec(NtauF,0.0));

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
                  for (int ishift=0; ishift<NsegShift; ishift++) // if(imeas%NsegShift==1)
                     shift_segment( segments[ifl], Beta, Levels[ifl], Uloc, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );
                  // flip segment - NOT IMPLEMENTED
                  //for (int i=0; i<N_flip; i++)flip_segment( segments_up, NtauF_m1, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
               }
               //
               //....................Observables Measurments....................
               //
               // sign among the segments
               s *= sign[ifl];
               //
               // pertrbation order
               if( segments[ifl].size()<Norder ) Order[ifl][(int)segments[ifl].size()] += 1/(double)Nmeas_;
               //
               // flavour occupation - measurment averaged
               N_tmp[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/(Beta*Nmeas_);
               //
               // Green's functions - measurment averaged
               if( !Gexp && segments[ifl].size()>0) measure_G( G_tmp[ifl], segments[ifl], M[ifl], Beta );
               //...............................................................
            }
            //
            //  sign among the segments - measurment averaged - step summed
            sign_meas += s/(double)Nmeas_;

            //
            //global spin flip (does not require overlap calculation)
            for (int iswap=0; iswap<NspinSwap; iswap++) // if(imeas%NspinSwap==1)
            {
               bool SpinSwap;
               for (int ifl=0; ifl<Nflavor; ifl++) SpinSwap = ( M[ifl].rows() == 0 ) ? false : true;
               if(SpinSwap) swap_spins( Beta, F, segments, full_line, sign, M );
            }

         }

         //
         //........................Observables sweep sums.......................
         //
         // flavour occupation - symmetrization - step sum
         if(OrbSym) orb_symm( N_tmp, SetsOrbs );
         accumulate_Vec( Nloc, N_tmp );
         //
         // Green's functions - symmetrization - step sum
         if(Gexp)
         {
            for (int ifl=0; ifl<Nflavor; ifl++)
            {
               if (segments[ifl].size()>0) measure_G( G_tmp[ifl], segments[ifl], M[ifl], Beta );
            }
         }
         double Gnorm = 1.0/(NtauF-1);
         if(!Gexp) Gnorm *= Nmeas_;
         G_tmp = normalize_VecVec( G_tmp, Gnorm );
         //
         if(paramagnet) spin_symm( G_tmp );
         if(OrbSym) orb_symm( G_tmp, SetsOrbs );
         accumulate_VecVec( G, G_tmp );
         //
         // n(tau) - computed every Nmeas - symmetrization - step sum
         VecVec n_tau(Nflavor,Vec(NtauB,0.0));
         n_tau = measure_nt( segments, full_line, NtauB, Beta );
         if(OrbSym) orb_symm( n_tau, SetsOrbs ); // TEST
         accumulate_VecVec( nt, n_tau );
         //
         // Histograms - computed every Nmeas - step sum
         accumulate_Nhist( Nhist, n_tau );
         accumulate_Szhist( Szhist, n_tau );
         //
         // n(tau)n(0) - computed every Nmeas - step sum
         if(NnntMeas>1) accumulate_nnt( nnt, n_tau, s, NnntMeas );
         if(NnntMeas==1) accumulate_nnt( nnt, n_tau );
         //
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
                  insert_remove_segment( Beta*rndm(), Beta, Levels[ifl], Uloc, F[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);
                  // shift segments
                  for (int ishift=0; ishift<NsegShift; ishift++) // if(imeas%NsegShift==1)
                     shift_segment( segments[ifl], Beta, Levels[ifl], Uloc, F[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );
                  // flip segment - NOT IMPLEMENTED
                  //for (int i=0; i<N_flip; i++)flip_segment( segments_up, NtauF_m1, Beta, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
               }
               //
               //....................Observables Measurments....................
               //
               // flavour occupation - measurment averaged - step summed
               if(!fulldry) Nloc[ifl] += compute_overlap(full_segment, segments[ifl], full_line[ifl], Beta)/(Beta*Nmeas_);
               //...............................................................
            }
            //
            //global spin flip (does not require overlap calculation)
            for (int iswap=0; iswap<NspinSwap; iswap++) // if(imeas%NspinSwap==1)
            {
               bool SpinSwap;
               for (int ifl=0; ifl<Nflavor; ifl++) SpinSwap = ( M[ifl].rows() == 0 ) ? false : true;
               if(SpinSwap) swap_spins( Beta, F, segments, full_line, sign, M );
            }

         }
      }

      //----------------------------------------------------------------------//

      void set_levels(void)
      {
         Levels.resize(Nflavor,0.0);
         for (int ifl=0; ifl<Nflavor; ifl++) Levels[ifl] = mu - Eloc[ifl] + Hartree_shift[(int)(ifl/2)] + Screening_shift[(int)(ifl/2)];
         print_Vec(inputDir+"/used.Levels.DAT", Levels);
      }

      //----------------------------------------------------------------------//

      void print_observables(path pad, std::vector<int> bins)
      {
         if(debug) // All ranks print their own observables - NOT orbitally symmetrized
         {
            //
            mpi.report(" Rank #"+str(mpi.rank())+" is printing observables.");
            //
            // density
            Vec PrintNloc = get_Nloc(debug); // provides spin-orbital occupation already normalized
            print_Vec(resultsDir+"/Nqmc_rank"+str(mpi.rank())+pad, PrintNloc, mu);
            mpi.report(" Nqmc_rank"+str(mpi.rank())+pad+" is printed.");
            //
            // perturbation order
            print_VecVec(resultsDir+"/Order_rank"+str(mpi.rank())+pad, Order, 0.0, (double)RankSweeps);
            mpi.report(" Order_rank"+str(mpi.rank())+pad+" is printed.");
            //
            // error estimate and binning
            if(bins[0]>0)
            {
               VecVec Gerr(Nflavor,Vec(NtauF,0.0));
               binAverageVecVec( bins, G, Gerr );
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
            // observables derived from n(tau)
            print_VecVec(resultsDir+"/n_t_rank"+str(mpi.rank())+pad, nt, Beta, (double)RankSweeps);
            mpi.report(" n_t_rank"+str(mpi.rank())+pad+" is printed.");
            //
            print_Vec(resultsDir+"/Nhist_rank"+str(mpi.rank())+pad, Nhist );
            mpi.report(" Nhist_rank"+str(mpi.rank())+pad+" is printed.");
            //
            print_Vec(resultsDir+"/Szhist_rank"+str(mpi.rank())+pad, Szhist );
            mpi.report(" Szhist_rank"+str(mpi.rank())+pad+" is printed.");
            //
            // n(tau)n(0)
            if(NnntMeas>0)
            {
               // error estimate and binning
               if(bins[0]>0)
               {
                  VecVec Nerr(Nflavor*(Nflavor+1)/2,Vec(NtauB,0.0));
                  binAverageVecVec( bins, nnt, Nerr );
                  print_VecVec(resultsDir+"/nn_err_rank"+str(mpi.rank())+pad, Nerr, Beta, (double)RankSweeps);
                  mpi.report(" nn_err_rank"+str(mpi.rank())+pad+" is printed.");
               }
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
         // perturbation order
         VecVec NormOrder = normalize_VecVec(Order, RankSweeps); // done before allreduce to overcome overfloat error for large RankSweeps
         VecVec PrintOrder(Nflavor,Vec(Norder,0.0));
         mpi.allreduce(NormOrder, PrintOrder, true);
         if(mpi.is_master()) print_VecVec(resultsDir+"/Order"+pad, PrintOrder);
         mpi.report(" Order"+pad+" is printed.");
         //
         // error estimate and binning
         if(bins[0]>0)
         {
            VecVec Gerr(Nflavor,Vec(NtauF,0.0));
            binAverageVecVec( bins, G, Gerr );
            Gerr = normalize_VecVec(Gerr, RankSweeps);
            VecVec PrintGerr(Nflavor,Vec(NtauF,0.0));
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
         VecVec NormG = normalize_VecVec(G, RankSweeps); // done before allreduce to overcome overfloat error for large RankSweeps
         VecVec PrintG(Nflavor,Vec(NtauF,0.0));
         mpi.allreduce(NormG, PrintG, true);
         if(mpi.is_master()) print_VecVec(resultsDir+"/Gimp_t"+pad, PrintG, Beta);
         mpi.report(" Gimp_t"+pad+" is printed.");
         //
         // observables derived from n(tau)
         VecVec Normnt = normalize_VecVec(nt, RankSweeps);
         VecVec Printnt(Nflavor,Vec(NtauB,0.0));
         mpi.allreduce(Normnt, Printnt, true);
         if(mpi.is_master()) print_VecVec(resultsDir+"/n_t"+pad, Printnt, Beta);
         mpi.report(" n_t"+pad+" is printed.");
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
         // n(tau)n(0)
         if(NnntMeas>0)
         {
            // error estimate and binning
            if(bins[0]>0)
            {
               VecVec Nerr(Nflavor*(Nflavor+1)/2,Vec(NtauB,0.0));
               binAverageVecVec( bins, nnt, Nerr );
               Nerr = normalize_VecVec(Nerr, RankSweeps);
               VecVec PrintNerr(Nflavor,Vec(NtauF,0.0));
               mpi.allreduce(Nerr, PrintNerr, true);
               if(mpi.is_master()) print_VecVec(resultsDir+"/nn_err"+pad, PrintNerr, Beta);
               mpi.report(" nn_err"+pad+" is printed.");
            }
            VecVec Normnnt = normalize_VecVec(nnt, RankSweeps);
            VecVec Printnnt(Nflavor*(Nflavor+1)/2,Vec(NtauB,0.0));
            mpi.allreduce(Normnnt, Printnnt, true);
            if(mpi.is_master())
            {
               if(OrbSym) orb_symm( Printnnt, SetsOrbs, Norb );
               print_VecVec(resultsDir+"/nn_t"+pad, Printnnt, Beta);
            }
            mpi.report(" nn_t"+pad+" is printed.");
         }
      }

      //----------------------------------------------------------------------//

};


//============================================================================//


#endif
