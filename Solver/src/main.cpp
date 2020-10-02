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
   int Nsite,Nspin,Ntau;
   int Norder,Nmeas,Nshift;
   int PrintTime,binlength,verbosity,muIter;
   double Beta,mu,muStep,muTime,muErr,density;
   bool retarded;
   std::vector<int> SiteTime;
   std::vector<int> SiteNorb;
   std::vector<std::string> SiteName;
   std::vector<std::string> SiteDir;
   char* IterationDir;


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
      IterationDir = argv[2]; // iteration folder
      //
      find_param(argv[1], "mu"         , mu        );
      find_param(argv[1], "beta"       , Beta      );
      find_param(argv[1], "Nsite"      , Nsite     );
      find_param(argv[1], "Nspin"      , Nspin     );
      find_param(argv[1], "Ntau"       , Ntau      );
      find_param(argv[1], "Norder"     , Norder    );
      find_param(argv[1], "Nmeas"      , Nmeas     );
      find_param(argv[1], "Nshift"     , Nshift    );
      find_param(argv[1], "PrintTime"  , PrintTime );
      find_param(argv[1], "verbosity"  , verbosity );
      find_param(argv[1], "binlength"  , binlength );
      find_param(argv[1], "retarded"   , retarded  );
      find_param(argv[1], "density"    , density   );
      find_param(argv[1], "muStep"     , muStep    );
      find_param(argv[1], "muIter"     , muIter    );
      find_param(argv[1], "muTime"     , muTime    );
      find_param(argv[1], "muErr"      , muErr     );
      //
      if(verbosity>=1 && mpi.is_master())
      {
         mpi.report(" mu= "+str(mu));
         mpi.report(" beta= "+str(Beta));
         mpi.report(" Nsite= "+str(Nsite));
         mpi.report(" Nspin= "+str(Nspin));
         mpi.report(" Ntau= "+str(Ntau));
         mpi.report(" Norder= "+str(Norder));
         mpi.report(" Nmeas= "+str(Nmeas));
         mpi.report(" Nshift= "+str(Nshift));
         mpi.report(" PrintTime= "+str(PrintTime));
         mpi.report(" binlength= "+str(binlength));
         mpi.report(" retarded= "+str(retarded));
         if(density>0.0)
         {
            mpi.report(" density= "+str(density));
            mpi.report(" muStep= "+str(muStep));
            mpi.report(" muIter= "+str(muIter));
            mpi.report(" muErr= "+str(muErr));
            mpi.report(" muTime= "+str(muTime));
         }
      }

      //
      for(int isite=0; isite < Nsite; isite++)
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
         strcpy(lineT,"Time");
         strcat(lineT,s);
         find_param(argv[1], lineT, minutes );
         SiteTime.push_back(minutes);
         //
         strcpy(lineO,"Norb");
         strcat(lineO,s);
         find_param(argv[1], lineO, Norb );
         SiteNorb.push_back(Norb);
         //
         strcpy(lineN,"Name");
         strcat(lineN,s);
         find_param(argv[1], lineN, element );
         SiteName.push_back(element);
         //
         strcpy(folder,IterationDir);
         strcat(folder,"Solver_");
         strcat(folder,element);
         SiteDir.push_back(folder);
         //
      }


      //......................................................................//
      //                             Solver setup                             //
      //......................................................................//
      print_line_space(1,mpi.is_master());
      std::vector<ct_hyb> ImpurityList;
      for(int isite=0; isite < Nsite; isite++)
      {
         mpi.report(" isite = "+str(isite));
         mpi.report(" Name = "+SiteName[isite]);
         mpi.report(" Norb = "+str(SiteNorb[isite]));
         mpi.report(" Time = "+str(SiteTime[isite]));

         //
         if(PathExist(strcpy(new char[SiteDir[isite].length() + 1], SiteDir[isite].c_str())))
         {
            mpi.report(" Folder = "+SiteDir[isite]+" (Found).");
            ImpurityList.push_back(ct_hyb( SiteName[isite], mu, Beta, Nspin, SiteNorb[isite],
                                           Ntau, Norder, Nmeas, Nshift, retarded,
                                           PrintTime, binlength, mpi, true));
            ImpurityList[isite].init( SiteDir[isite] );
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
      if(density>0.0)
      {
         print_line_space(1,mpi.is_master());
         double trial_density;
         double mu_start=mu;
         double mu_new,mu_last;
         std::vector<double>Ntmp(Nsite,0.0);

         //
         print_line_minus(80,mpi.is_master());
         mpi.report(" Starting chemical potential: "+str(mu_start));
         for(int isite=0; isite < Nsite; isite++)
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
            for(int isite=0; isite < Nsite; isite++)
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
            mu_last=mu_new;
         }

         //
         print_line_space(1,mpi.is_master());
         double mu_below=std::min(mu_new,mu_last);
         double mu_above=std::max(mu_new,mu_last);
         for (int imu=0; imu < muIter; imu++)
         {
            //
            mpi.report(" Chemical potential boundaries: "+str(mu_below)+" / "+str(mu_above));
            mu_new = (mu_below+mu_above)/2.0;
            mpi.report(" Setting chemical potential: "+str(mu_new));
            //
            for(int isite=0; isite < Nsite; isite++)
            {
               mpi.report(" Quick solution ("+str((int)(muTime*60))+"sec) of site "+SiteName[isite]);
               ImpurityList[isite].reset_mu( mu_new );
               ImpurityList[isite].solve( muTime, true );
               Ntmp[isite]=ImpurityList[isite].get_Density();
               mpi.report(" Site density: "+str(Ntmp[isite]));
            }
            trial_density = std::accumulate(Ntmp.begin(), Ntmp.end(), 0.0);
            mpi.report(" Total density: "+str(trial_density)+" Error: "+str(fabs(trial_density-density)));
            print_line_space(1,mpi.is_master());
            //
            if(trial_density > density) mu_above=mu_new;
            if(trial_density < density) mu_below=mu_new;
            if(fabs(trial_density-density)<muErr)
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

         //
         for(int isite=0; isite < Nsite; isite++) ImpurityList[isite].reset_mu( mu_new );

      }


      //..................................................
      //               Solve the impurities
      //..................................................
      for(int isite=0; isite < Nsite; isite++) ImpurityList[isite].solve( SiteTime[isite] );



      //..................................................
      //                 stop global timer
      //..................................................
      Tend = std::chrono::system_clock::now();
      std::chrono::duration<double> runtime_seconds = Tend-Tstart;
      cout << endl;
      cout << "Time [total] = " << runtime_seconds.count() << "s\n";
      print_line_star(60);
      exit(1);
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
