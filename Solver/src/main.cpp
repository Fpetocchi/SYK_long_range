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
#include "file_io.hpp"
#include "impurity.hpp"
#include "mpi.hpp"
//
using namespace std;


//============================================================================//


int main(int argc, char *argv[])
{
   //.........................................................................//
   //                           MPI Environment                               //
   //.........................................................................//
   CustomMPI mpi;
   mpi.init();


   //.........................................................................//
   //                           input variables                               //
   //.........................................................................//
   // Global Vars
   double Beta;
   int Nspin,NtauF,NtauB,Norder,Nmeas,Ntherm,Nshift,printTime;
   //logical flags and compatibility typo fix
   bool paramagnet,retarded,nnt_meas,quickloops,dichotomy,OrbSym;
   int para_read,ret_read,nnt_read,quick_read,sym_read;
   // Post-processing of the Green's function
   int binlength,binstart;
   // Density lookup algorithm (dichotomy)
   int muIter;
   double density,muStep,muTime,muErr;
   // Site Dependent Vars
   int Nimp;
   std::vector<int> SiteTime;
   std::vector<int> SiteNorb;
   std::vector<std::string> SiteName;
   std::vector<std::string> SiteDir;
   char* IterationDir;
   //Equivalent manifold per site
   std::vector<std::vector<int>> SiteSetsNorb;

#ifdef _verb
   bool debug=true;
#else
   bool debug=false;
#endif

   //.........................................................................//
   //                         start global timer                              //
   //.........................................................................//
   std::chrono::time_point<std::chrono::system_clock> Tstart, Tend;
   Tstart = std::chrono::system_clock::now();


   try
   {
      //......................................................................//
      //                           read inputfile                             //
      //......................................................................//
      if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");
      //
      // Iteration folder
      IterationDir = argv[2];
      // Global Vars
      find_param(argv[1], "BETA"       , Beta      );
      find_param(argv[1], "NSPIN"      , Nspin     );
      find_param(argv[1], "NTAU_F_IMP" , NtauF     );
      find_param(argv[1], "NTAU_B_IMP" , NtauB     );
      find_param(argv[1], "NORDER"     , Norder    );
      find_param(argv[1], "NMEAS"      , Nmeas     );
      find_param(argv[1], "NTHERM"     , Ntherm    );
      find_param(argv[1], "NSHIFT"     , Nshift    );
      find_param(argv[1], "PRINT_TIME" , printTime );
      find_param(argv[1], "PARAMAGNET" , para_read ); paramagnet = (para_read == 1) ? true : false;
      find_param(argv[1], "RETARDED"   , ret_read  ); retarded = (ret_read == 1) ? true : false;
      find_param(argv[1], "NNT_MEAS"   , nnt_read  ); nnt_meas = (nnt_read == 1) ? true : false;
      // Post-processing of the Green's function
      find_param(argv[1], "BINLENGTH"  , binlength );
      find_param(argv[1], "BINSTART"   , binstart  );
      // Density lookup algorithm (dichotomy)
      find_param(argv[1], "N_READ_IMP" , density   );
      find_param(argv[1], "MU_STEP"    , muStep    );
      find_param(argv[1], "MU_ITER"    , muIter    );
      find_param(argv[1], "MU_TIME"    , muTime    );
      find_param(argv[1], "N_ERR"      , muErr     );
      find_param(argv[1], "N_QUICK"    , quick_read); quickloops = (quick_read == 1) ? true : false;
      //Symmetrization type
      find_param(argv[1], "SYM_MODE"   , sym_read  ); OrbSym = (sym_read > 1) ? true : false;
      // Site Dependent Vars
      find_param(argv[1], "NIMP"       , Nimp      );
      //
      if(mpi.is_master()) //debug &&
      {
         mpi.report(" beta= "+str(Beta));
         mpi.report(" Nspin= "+str(Nspin));
         mpi.report(" NtauF= "+str(NtauF));
         mpi.report(" NtauB= "+str(NtauB));
         mpi.report(" Norder= "+str(Norder));
         mpi.report(" Nmeas= "+str(Nmeas));
         mpi.report(" Ntherm= "+str(Ntherm));
         mpi.report(" Nshift= "+str(Nshift));
         mpi.report(" printTime= "+str(printTime)+"min");
         mpi.report(" retarded= "+str(retarded));
         mpi.report(" quickloops= "+str(quickloops));
         mpi.report(" paramagnet= "+str(paramagnet));
         mpi.report(" nnt_meas= "+str(nnt_meas));
         mpi.report(" OrbSym= "+str(OrbSym));
         mpi.report(" debug= "+str(debug));
         if(binlength>0)
         {
            mpi.report(" binlength= "+str(binlength));
            mpi.report(" binstart= "+str(binstart));
         }
         mpi.report(" Nimp= "+str(Nimp));
         if((density>0.0)&&quickloops)
         {
            mpi.report(" density= "+str(density));
            mpi.report(" muStep= "+str(muStep));
            mpi.report(" muIter= "+str(muIter));
            mpi.report(" muErr= "+str(muErr));
            mpi.report(" muTime= "+str(muTime)+"min");
         }
      }

      //
      if(OrbSym)SiteSetsNorb.resize(Nimp);

      //
      for(int isite=0; isite < Nimp; isite++)
      {
         std::string ss=str(isite+1);
         const char * s  = ss.c_str();
         char lineN[5];
         char lineT[5];
         char lineO[5];
         char folder[999];
         char element[2];
         int minutes,Norb;
         //
         strcpy(lineT,"TIME_");
         strcat(lineT,s);
         find_param(argv[1], lineT, minutes );
         SiteTime.push_back(minutes);
         //
         strcpy(lineO,"NORB_");
         strcat(lineO,s);
         find_param(argv[1], lineO, Norb );
         SiteNorb.push_back(Norb);
         //
         strcpy(lineN,"NAME_");
         strcat(lineN,s);
         find_param(argv[1], lineN, element );
         SiteName.push_back(element);
         //
         strcpy(folder,IterationDir);
         strcat(folder,"Solver_");
         strcat(folder,element);
         SiteDir.push_back(folder);
         //

         //
         if(OrbSym)
         {
            //
            char lineSets[10];
            char lineNorb[12];
            int Sets,Norb;

            //
            strcpy(lineSets,"EQV_");
            strcat(lineSets,s);
            strcat(lineSets,"_SETS");
            find_param(argv[1], lineSets, Sets );
            if(Sets>0)
            {
               //
               SiteSetsNorb[isite].resize(Sets);

               //
               for(int iset=0; iset < Sets; iset++)
               {
                  std::string sst=str(iset+1);
                  const char * st  = sst.c_str();

                  //
                  strcpy(lineNorb,"EQV_");
                  strcat(lineNorb,s);
                  strcat(lineNorb,"_NORB_");
                  strcat(lineNorb,st);
                  find_param(argv[1], lineNorb, Norb );
                  SiteSetsNorb[isite][iset] = Norb;

               }
            }
         }


      }


      //......................................................................//
      //                             Solver setup                             //
      //......................................................................//
      print_line_space(1,mpi.is_master());
      std::vector<ct_hyb> ImpurityList;
      for(int isite=0; isite < Nimp; isite++)
      {
         mpi.report(" isite = "+str(isite));
         mpi.report(" Name = "+SiteName[isite]);
         mpi.report(" Norb = "+str(SiteNorb[isite]));
         mpi.report(" Time = "+str(SiteTime[isite])+"min");
         if(OrbSym)
         {
            for(int iset=0; iset < SiteSetsNorb[isite].size(); iset++) mpi.report(" Norb in Eqv set #"+str(iset+1)+" = "+str(SiteSetsNorb[isite][iset]));
         }

         //
         if(PathExist(strcpy(new char[SiteDir[isite].length() + 1], SiteDir[isite].c_str())))
         {
            mpi.report(" Folder = "+SiteDir[isite]+" (Found).");
            ImpurityList.push_back(ct_hyb( SiteName[isite], Beta, Nspin, SiteNorb[isite],
                                           NtauF, NtauB, Norder, Nmeas, Ntherm, Nshift,
                                           paramagnet, retarded, nnt_meas, SiteSetsNorb[isite],
                                           printTime, std::vector<int> { binlength,binstart }, mpi ));
            ImpurityList[isite].init( SiteDir[isite]);
         }
         else
         {
            mpi.StopError(" Folder = "+SiteDir[isite]+" (Not Found) - Exiting. \n");
         }
         print_line_space(1,mpi.is_master());
      }


      //......................................................................//
      //                     Chemical potential Lookup                        //
      //......................................................................//
      if((density>0.0)&&quickloops)
      {
         print_line_space(1,mpi.is_master());
         double trial_density;
         double mu_start=ImpurityList[0].get_mu();
         double mu_new,mu_last=ImpurityList[0].get_mu();
         std::vector<double>Ntmp(Nimp,0.0);

         //
         dichotomy=true;

         //
         print_line_minus(80,mpi.is_master());
         mpi.report(" Starting chemical potential: "+str(mu_start));
         for(int isite=0; isite < Nimp; isite++)
         {
            mpi.report(" Quick solution ("+str((int)(muTime*60))+"sec) of site "+SiteName[isite]);
            ImpurityList[isite].solve( muTime, true );
            Ntmp[isite]=ImpurityList[isite].get_Density();
            mpi.report(" Site density: "+str(Ntmp[isite]));
         }
         trial_density = std::accumulate(Ntmp.begin(), Ntmp.end(), 0.0);
         mpi.report(" Total density: "+str(trial_density));
         print_line_space(1,mpi.is_master());

         //
         double muSign = (trial_density < density) ? +1.0 : -1.0;
         for (int imu=1; imu < muIter; imu++)
         {
            //
            mu_new = mu_start + imu*muStep*muSign;
            mpi.report(" Setting chemical potential: "+str(mu_new));
            //
            for(int isite=0; isite < Nimp; isite++)
            {
               mpi.report(" Quick solution ("+str((int)(muTime*60))+"sec) of site "+SiteName[isite]);
               ImpurityList[isite].reset_mu( mu_new );
               ImpurityList[isite].solve( muTime, true );
               Ntmp[isite]=ImpurityList[isite].get_Density();
               mpi.report(" Site density: "+str(Ntmp[isite]));
            }
            trial_density = std::accumulate(Ntmp.begin(), Ntmp.end(), 0.0);
            mpi.report(" Total density: "+str(trial_density));
            print_line_space(1,mpi.is_master());
            //
            if((muSign>0.0)&&(trial_density > density)) break;
            if((muSign<0.0)&&(trial_density < density)) break;
            if(imu == muIter-1)
            {
               mpi.report(" *NOT* found chemical potential after "+str(muIter)+" rigid shifts. Last value: "+str(mu_new));
               mpi.report(" (User should try to increase either MU_ITER or MU_STEP.)");
               dichotomy=false;
               break;
            }
            //
            mu_last=mu_new;
         }

         //
         if(dichotomy)
         {
            print_line_space(3,mpi.is_master());
            double mu_below=std::min(mu_new,mu_last);
            double mu_above=std::max(mu_new,mu_last);
            for (int imu=0; imu < muIter; imu++)
            {
               //
               mpi.report(" Chemical potential boundaries: "+str(mu_below)+" / "+str(mu_above));
               mu_new = (mu_below+mu_above)/2.0;
               mpi.report(" Setting chemical potential: "+str(mu_new));
               //
               for(int isite=0; isite < Nimp; isite++)
               {
                  mpi.report(" Quick solution ("+str((int)(muTime*60))+"sec) of site "+SiteName[isite]);
                  ImpurityList[isite].reset_mu( mu_new );
                  ImpurityList[isite].solve( muTime, true );
                  Ntmp[isite]=ImpurityList[isite].get_Density();
                  mpi.report(" Site density: "+str(Ntmp[isite]));
               }
               trial_density = std::accumulate(Ntmp.begin(), Ntmp.end(), 0.0);
               mpi.report(" Total density: "+str(trial_density)+" relative error: "+str(fabs(trial_density-density)/density));
               print_line_space(1,mpi.is_master());
               //
               if(trial_density > density) mu_above=mu_new;
               if(trial_density < density) mu_below=mu_new;
               if((fabs(trial_density-density)/density)<muErr)
               {
                  mpi.report(" Found correct chemical potential after "+str(imu)+" iterations: "+str(mu_new));
                  break;
               }
               else if(imu == muIter-1)
               {
                  mpi.report(" *NOT* found chemical potential after "+str(muIter)+" iterations. Last value: "+str(mu_new));
                  break;
               }
            }
         }

         //
         for(int isite=0; isite < Nimp; isite++) ImpurityList[isite].reset_mu( mu_new );

      }


      //..................................................
      //               Solve the impurities
      //..................................................
      for(int isite=0; isite < Nimp; isite++) ImpurityList[isite].solve( SiteTime[isite] );



      //..................................................
      //                 stop global timer
      //..................................................
      Tend = std::chrono::system_clock::now();
      std::chrono::duration<double> runtime_seconds = Tend-Tstart;
      cout << endl;
      mpi.report(" Time [total] = "+str(runtime_seconds.count()));
      print_line_star(60,mpi.is_master());
      mpi.finalize();
   } // try
   catch (char *message)
   {
      cerr << "eception\n**** " << message << " ****" << endl;
      cerr << "input_file [ --test ]\n" << endl;
   }
   catch (...)
   {
      cerr << "unspecified exception " << endl;
      cerr << "\n input_file [ --test ]\n" << endl;
   }
   return 0;
}


//============================================================================//
