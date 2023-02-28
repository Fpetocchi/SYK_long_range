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
              bool Improved_F, bool Improved_B,
              bool removeUhalf, bool paramagnet, bool retarded,
              std::vector<int> SetsNorb, bool full_ntOrbSym,
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
              Improved_F(Improved_F),
              Improved_B(Improved_B),
              removeUhalf(removeUhalf),
              paramagnet(paramagnet),
              retarded(retarded),
              SetsNorb(SetsNorb),
              full_ntOrbSym(full_ntOrbSym),
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
         if(retarded&&Improved_F) mandatoryFiles.insert( mandatoryFiles.end(), { "/Kp_t.DAT" } );
         if(OrbSym) mandatoryFiles.push_back("/Eqv.DAT");
         for(int ifile=0; ifile < (int)mandatoryFiles.size(); ifile++)
         {
            path filepath=inputDir+mandatoryFiles[ifile];
            if(!PathExist(strcpy(new char[filepath.length() + 1], filepath.c_str())))
            {
               mpi.StopError(filepath+" (Not Found) - Exiting.");
            }
            else
            {
               mpi.report(" Found file: "+filepath);
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
         //read_VecVec(inputDir+"/Delta_t.DAT", Delta, Nflavor, true, true);  // last flag is to reverse the tau index
         tau_uniform_D = read_Delta( inputDir+"/Delta_t.DAT", Delta, Nflavor );
         for(int ifl=0; ifl < Nflavor; ifl++)
         {
            path Dcomp = "Delta["+str(ifl)+"]";
            int Ntau_D =  Delta[ifl].size();
            //if(Ntau_D!=NtauF) mpi.report(" The Number of tau points in "+Dcomp+" is: "+str(Ntau_D)+" different from NtauF: "+str(NtauF));
            for(int itau=0; itau < Ntau_D; itau++)
            {
               if(Delta[ifl][itau]<0.0) mpi.StopError( " ->"+Dcomp+" at tau "+str(itau)+" is positive - Exiting.");
            }
         }

         //
         //read the screening function
         Screening_shift.resize(Norb,0.0);
         if(retarded)
         {
            //
            //Read the screening from file
            mpi.report(" Reading screening file.");
            read_EigenMat(inputDir+"/Screening.DAT", Screening_Mat, Nflavor, Nflavor);
            if(!removeUhalf)
            {
               for(int iorb=0; iorb < Norb; iorb++)
               {
                  // Hartree-like screening
                  // for(int ifl=0; ifl < Nflavor; ifl++) Screening_shift[iorb] += (2*iorb!=ifl) ? Screening_Mat(2*iorb,ifl)/2.0 : 0.0 ;
                  // diagonal screening
                  Screening_shift[iorb] = Screening_Mat(2*iorb,2*iorb)/2.0;
                  mpi.report(" Orbital "+str(iorb)+" - screening shift = S/2: "+str(Screening_shift[iorb],6));
               }
            }

            //
            //read_VecVecVec(inputDir+"/K_t.DAT", K_table, Nflavor, true);
            tau_uniform_K = read_K( inputDir+"/K_t.DAT", K_table, Nflavor );
            for(int ifl=0; ifl < Nflavor; ifl++)
            {
               for(int jfl=0; jfl <= ifl; jfl++)
               {
                  int Ntau_K = K_table[ifl][jfl].size();
                  path Kcomp = "K["+str(ifl)+"]["+str(jfl)+"]";
                  //if(Ntau_K!=NtauB) mpi.report(" The Number of tau points in "+Kcomp+" is: "+str(Ntau_K)+" different from NtauB: "+str(NtauB));
                  //
                  double K_zero =  K_table[ifl][jfl].front();
                  double K_beta =  K_table[ifl][jfl].back();
                  if(K_zero!=0.0) mpi.StopError( " ->"+Kcomp+" at tau=0 is not vanishing - Exiting.");
                  if(K_beta!=0.0) mpi.StopError( " ->"+Kcomp+" at tau=beta is not vanishing - Exiting.");
               }
            }
            //
            if(Improved_F)
            {
               //read_VecVecVec(inputDir+"/Kp_t.DAT", Kp_table, Nflavor, true);
               tau_uniform_Kp = read_K( inputDir+"/Kp_t.DAT", Kp_table, Nflavor );
               if( tau_uniform_K != tau_uniform_Kp ) mpi.StopError( " K_t.DAT and Kp_t.DAT are on different tau meshes - Exiting.");
               for(int ifl=0; ifl < Nflavor; ifl++)
               {
                  for(int jfl=0; jfl <= ifl; jfl++)
                  {
                     int Ntau_K = Kp_table[ifl][jfl].size();
                     path Kcomp = "Kp["+str(ifl)+"]["+str(jfl)+"]";
                     //if(Ntau_K!=NtauB) mpi.report(" The Number of tau points in "+Kcomp+" is: "+str(Ntau_K)+" different from NtauB: "+str(NtauB));
                     //
                     double Kp_zero =  Kp_table[ifl][jfl].front();
                     double Kp_beta =  Kp_table[ifl][jfl].back();
                     if((Kp_zero+Kp_beta)!=0.0) mpi.StopError( " ->"+Kcomp+" is not symmetric with respect to beta/2 - Exiting.");
                     //if((Kp_zero+Kp_beta)!=0.0) mpi.report( " "+Kcomp+" is not symmetric with respect to beta/2 - Warning.");
                     if(removeUhalf)
                     {
                        Kp_table[ifl][ifl][0] = 0.0;
                        Kp_table[ifl][ifl][Ntau_K-1] = 0.0;
                     }
                  }
               }
            }
         }

         //
         //rescale the chemical potential with the half-filling one
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
            print_VecVec(inputDir+"/used.Delta_t.DAT", Delta);
            print_EigenMat(inputDir+"/used.Umat.DAT", Uloc);
            if(retarded)
            {
               print_VecVecVec(inputDir+"/used.K_t.DAT", K_table);
               print_EigenMat(inputDir+"/used.Screening.DAT", Screening_Mat);
               if(Improved_F)print_VecVecVec(inputDir+"/used.Kp_t.DAT", Kp_table);
            }
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
         if(Improved_F) F_S.resize( Nflavor, Vec( NtauF, 0.0 ) );
         if(Improved_F&&retarded) F_R.resize( Nflavor, Vec( NtauF, 0.0 ) );
         Gerr.resize( Nflavor, Vec( NtauF, 0.0 ) );
         nt.resize( Nflavor, Vec( NtauB, 0.0 ) );
         nnt.resize( Nflavor*(Nflavor+1)/2, Vec( NtauB, 0.0 ) );
         Dimp.resize( Nflavor, Vec( Nflavor, 0.0 ) );
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
                  for (int ifl=0; ifl<Nflavor; ifl++)
                     if (full_line[ifl]==1) mpi.report(" Flavor #"+str(ifl)+" has a full segment.");
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
      bool                                Improved_F;
      bool                                Improved_B;
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
      bool                                full_ntOrbSym;
      bool                                OrbSym;
      // Input data
      Vec                                 Levels;                               // mu - <\epsilon> + mu_correction
      Vec                                 Hartree_shift;                        // chemical potential shift
      Vec                                 Screening_shift;                      // chemical potential correction due to shift encoded in the screening function
      Vec                                 Eloc;                                 // <\epsilon>
      Mat                                 Uloc;                                 // Istantaneous U matrix
      Mat                                 Screening_Mat;                        // Screening matrix
      VecVec                              Delta;                                // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu)
      VecVecVec                           K_table;                              // K function matrix for retarded interactions
      VecVecVec                           Kp_table;
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
      VecVec                              Dimp;                                 // Double occupation matrix
      Vec                                 Nhist;                                //
      Vec                                 Szhist;                               //
      VecVec                              G;                                    // Impurity Green's function
      VecVec                              F_S, F_R;                             // Impurity Green's function Improved estimator
      VecVec                              Gerr;                                 // Binning error on Impurity Green's function
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
         VecVec F_tmp_S(Nflavor,Vec(NtauF,0.0));
         VecVec F_tmp_R(Nflavor,Vec(NtauF,0.0));

         // The measurments I'm going to do regardless from the time
         for (int imeas=0; imeas<Nmeas_; imeas++)
         {
            for (int ifl=0; ifl<Nflavor; ifl++)
            {
               // insert or remove full line
               if (segments[ifl].size() == 0) insert_remove_full_line( Levels[ifl], Uloc, Beta, full_line[ifl], segments, full_line, ifl );
               insert_remove_antisegment( Beta*rndm(), Beta, Levels[ifl], Uloc, Delta[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );

               //
               if (!full_line[ifl])
               {
                  // local update
                  insert_remove_segment( Beta*rndm(), Beta, Levels[ifl], Uloc, Delta[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);
                  // shift segments
                  for (int ishift=0; ishift<NsegShift; ishift++) // if(imeas%NsegShift==1)
                     shift_segment( segments[ifl], Beta, Levels[ifl], Uloc, Delta[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );
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
               //...............................................................
            }
            //
            // sign among the segments - measurment averaged - step sum
            sign_meas += s/(double)Nmeas_;

            //
            // Green's functions - measurment averaged - step sum
            if(!Gexp)
            {
               if(Improved_F)
               {
                  measure_GF( G_tmp, F_tmp_S, F_tmp_R, segments, full_line, M, Beta, Uloc, Kp_table, removeUhalf );
               }
               else
               {
                  measure_G( G_tmp, segments, M, Beta );
               }
            }

            //
            // global spin flip (does not require overlap calculation)
            for (int iswap=0; iswap<NspinSwap; iswap++) // if(imeas%NspinSwap==1)
            {
               bool SpinSwap;
               for (int ifl=0; ifl<Nflavor; ifl++) SpinSwap = ( M[ifl].rows() == 0 ) ? false : true;
               if(SpinSwap) swap_spins( Beta, Delta, segments, full_line, sign, M );
            }

         }

         //
         //........................Observables sweep sums.......................
         //
         // flavour occupation - symmetrization - step sum
         accumulate_Vec( Nloc, N_tmp );
         //
         // Green's functions - symmetrization - step sum
         if(Gexp)
         {
            if(Improved_F)
            {
               measure_GF( G_tmp, F_tmp_S, F_tmp_R, segments, full_line, M, Beta, Uloc, Kp_table, removeUhalf );
            }
            else
            {
               measure_G( G_tmp, segments, M, Beta );
            }
         }
         double Gnorm = 1.0/(NtauF-1);
         if(!Gexp) Gnorm *= Nmeas_;
         G_tmp = normalize_VecVec( G_tmp, Gnorm );
         F_tmp_S = normalize_VecVec( F_tmp_S, Gnorm );
         F_tmp_R = normalize_VecVec( F_tmp_R, Gnorm );
         //
         accumulate_VecVec( G, G_tmp );
         if(Improved_F)
         {
            accumulate_VecVec( F_S, F_tmp_S );
            if(retarded) accumulate_VecVec( F_R, F_tmp_R );
         }
         //
         // n(tau) - computed every Nmeas - symmetrization - step sum
         VecVec n_tau(Nflavor,Vec(NtauB,0.0));
         n_tau = measure_nt( segments, full_line, NtauB, Beta );
         //
         // this gives NaNa = NaNb if a is equivalent to b: is a more strict constraints
         if(full_ntOrbSym&&OrbSym) orb_symm( n_tau, SetsOrbs );
         //
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
               insert_remove_antisegment( Beta*rndm(), Beta, Levels[ifl], Uloc, Delta[ifl], full_line[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );

               //
               if (!full_line[ifl])
               {
                  // local update
                  insert_remove_segment( Beta*rndm(), Beta, Levels[ifl], Uloc, Delta[ifl], segments[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table);
                  // shift segments
                  for (int ishift=0; ishift<NsegShift; ishift++) // if(imeas%NsegShift==1)
                     shift_segment( segments[ifl], Beta, Levels[ifl], Uloc, Delta[ifl], M[ifl], sign[ifl], segments, full_line, ifl, K_table );
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
               if(SpinSwap) swap_spins( Beta, Delta, segments, full_line, sign, M );
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
            // n(tau)
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
         //
         // density - note that its never spin symmetrized because nn_t is not either
         Vec PrintNloc = get_Nloc();
         if(mpi.is_master())
         {
            if(OrbSym) orb_symm( PrintNloc, SetsOrbs );
            print_Vec(resultsDir+"/Nqmc"+pad, PrintNloc, mu);
         }
         mpi.report(" Nqmc"+pad+" is printed.");
         //
         //
         // perturbation order--------------------------------------------------
         dump_data( Order, "Order", pad, Nflavor, Norder );
         //
         //
         // charge configuration Histogram--------------------------------------
         dump_data( Nhist, "Nhist", pad, Nflavor+1 );
         //
         //
         // spin configuration Histogram----------------------------------------
         dump_data( Szhist, "Szhist", pad, Nflavor/2+1 );
         //
         //
         // n(tau)--------------------------------------------------------------
         dump_data( nt, "n_t", pad, Nflavor, NtauB, Beta, ( (!full_ntOrbSym && OrbSym) ? 0 : -1 ) );
         //
         //
         // Green's function----------------------------------------------------
         if(bins[0]>0)
         {
            // binning and related error estimate
            VecVec Err(Nflavor,Vec(NtauF,0.0));
            binAverageVecVec( bins, G, Err );
            dump_data( Err, "Gerr", pad, Nflavor, NtauF, Beta );
         }
         //
         // fix endpoints, note that Ntmp has no orbital symmetrizations same as G
         Vec Ntmp = Nloc;
         if(paramagnet)
         {
            spin_symm(Ntmp);
            spin_symm(G);
         }
         for (int ifl=0; ifl<Nflavor; ifl++)
         {
            G[ifl].front() = +Ntmp[ifl]-(long double)RankSweeps ;
            G[ifl].back()  = -Ntmp[ifl];
         }
         //
         // print
         dump_data( G, "Gimp_t", pad, Nflavor, NtauF, Beta, ( OrbSym ? 0 : -1 ) );
         //
         //
         // n(tau)n(0)----------------------------------------------------------
         if(NnntMeas>0)
         {
            // error estimate and binning
            if(bins[0]>0)
            {
               VecVec Err(Nflavor*(Nflavor+1)/2,Vec(NtauB,0.0));
               binAverageVecVec( bins, nnt, Err );
               dump_data( Err, "nn_err", pad, Nflavor*(Nflavor+1)/2, NtauB, Beta );
            }
            //
            // print
            dump_data( nnt, "nn_t", pad, Nflavor*(Nflavor+1)/2, NtauB, Beta, ( (!full_ntOrbSym && OrbSym) ? Norb : -1 ) );
            //
         }
         //
         //
         // double occupations: n(0)n(0)----------------------------------------
         if(NnntMeas>0)
         {
            VecVec nntD = nnt;
            if(OrbSym) orb_symm( nntD, SetsOrbs, Norb );
            int ndx=0;
            for (int ifl=0; ifl<Nflavor; ++ifl)
            {
               for (int jfl=0; jfl<=ifl; ++jfl)
               {
                  Dimp[ifl][jfl] = nntD[ndx][0];
                  Dimp[jfl][ifl] = nntD[ndx][0];
                  ndx++;
               }
            }
            dump_data( Dimp, "Dimp", pad, Nflavor, Nflavor );
         }
         //
         //
         // Improved estimator for the self-energy------------------------------
         if(Improved_F)
         {
            //
            if(bins[0]>0)
            {
               // binning and related error estimate
               {
                  VecVec Err(Nflavor,Vec(NtauF,0.0));
                  binAverageVecVec( bins, F_S, Err );
                  dump_data( Err, "Ferr_S", pad, Nflavor, NtauF, Beta );
               }
               if(retarded)
               {
                  VecVec Err(Nflavor,Vec(NtauF,0.0));
                  binAverageVecVec( bins, F_R, Err );
                  dump_data( Err, "Ferr_S", pad, Nflavor, NtauF, Beta );
               }
            }
            //
            // fix endpoints
            for (int ifl=0; ifl<Nflavor; ++ifl)
            {
               F_S[ifl].front() *=2;
               F_S[ifl].back() *=2;
               if(retarded)
               {
                  F_R[ifl].front() *=2;
                  F_R[ifl].back() *=2;
               }
            }
            //
            if(paramagnet)
            {
               spin_symm(F_S);
               if(retarded) spin_symm(F_R);
            }
            //
            // print
            dump_data( F_S, "Fimp_S_t", pad, Nflavor, NtauF, Beta, ( OrbSym ? 0 : -1 ) );
            if(retarded) dump_data( F_R, "Fimp_R_t", pad, Nflavor, NtauF, Beta, ( OrbSym ? 0 : -1 ) );
            //
         }
      }

      //----------------------------------------------------------------------//

      void dump_data( Vec &Data, path name, path pad, int dim )
      {
         Vec NormData = normalize_Vec(Data, RankSweeps);
         Vec PrintData(dim,0.0);
         mpi.allreduce(NormData, PrintData, true);
         if(mpi.is_master()) print_Vec(resultsDir+"/"+name+pad, PrintData);
         mpi.report(" "+name+pad+" is printed.");
      }

      void dump_data( VecVec &Data, path name, path pad, int dim1, int dim2, double Beta=0.0, int sym=-1)
      {
         VecVec NormData = normalize_VecVec(Data, RankSweeps);
         VecVec PrintData(dim1,Vec(dim2,0.0));
         mpi.allreduce(NormData, PrintData, true);
         if(mpi.is_master())
         {
            if(sym!=-1)
            {
               orb_symm( PrintData, SetsOrbs, sym );
               mpi.report(" "+name+pad+" is orbital symmetrized.");
            }
            print_VecVec(resultsDir+"/"+name+pad, PrintData, Beta );
         }
         mpi.report(" "+name+pad+" is printed.");
      }

      //----------------------------------------------------------------------//

      bool read_Delta( std::string path, std::vector<std::vector<double>> &D, int &dim )
      {
         //
         mpi.report( " Reading "+path );

         //
         // Read axis
         std::vector<double> tau_axis;
         std::string line;
         ifstream file( path );
         while( std::getline(file, line) ) // read one line from ifs
         {
            std::istringstream iss(line); // access line as a stream
            double tau_point;
            iss >> tau_point;
            tau_axis.push_back(tau_point);
         }
         file.close();

         //
         // Read Delta
         read_VecVec( path, D, dim, true, true);  // last flag is to reverse the tau index

         // consistency check
         int Ntau_axis = tau_axis.size();
         int Ntau_field = D[0].size();
         if(Ntau_axis!=Ntau_field) mpi.StopError( " Something is wrong with the tau mesh - Exiting.");

         //
         double dt1 = fabs( tau_axis[1]-tau_axis[0] );
         double dt2 = fabs( tau_axis[3]-tau_axis[2] );
         bool tau_uniform = fabs(dt2-dt1) < 1e-12;
         if( !tau_uniform )
         {
            int nseg = (Ntau_axis-1)/2;
            if( 2*(nseg/2) != nseg ) mpi.StopError( " Dense tau mesh is not a multiple of 2 and 4 - Exiting.");
            mpi.report( " tau mesh is found to be not uniform." );
         }

         //
         return tau_uniform;
      }

      bool read_K( std::string path, std::vector<std::vector<std::vector<double>>> &K, int &dim )
      {
         //
         mpi.report( " Reading "+path );

         //
         // Read axis
         std::vector<double> tau_axis;
         std::string line;
         ifstream file( path );
         while( std::getline(file, line) ) // read one line from ifs
         {
            std::istringstream iss(line); // access line as a stream
            double tau_point;
            iss >> tau_point;
            tau_axis.push_back(tau_point);
         }
         file.close();

         //
         // Read K_table
         read_VecVecVec( path, K, dim, true);

         // consistency check
         int Ntau_axis = tau_axis.size();
         int Ntau_field = K[0][0].size();
         if(Ntau_axis!=Ntau_field) mpi.StopError( " Something is wrong with the tau mesh - Exiting.");

         //
         double dt1 = fabs( tau_axis[1]-tau_axis[0] );
         double dt2 = fabs( tau_axis[3]-tau_axis[2] );
         bool tau_uniform = fabs(dt2-dt1) < 1e-12;
         if( !tau_uniform )
         {
            int nseg = (Ntau_axis-1)/2;
            if( 2*(nseg/2) != nseg ) mpi.StopError( " Dense tau mesh is not a multiple of 2 and 4 - Exiting.");
            mpi.report( " tau mesh is found to be not uniform." );
         }

         //
         return tau_uniform;
      }


};


//============================================================================//


#endif
