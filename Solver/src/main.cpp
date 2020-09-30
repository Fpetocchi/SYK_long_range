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
//
#include "file_io.hpp"
#include "impurity.hpp"
//
using namespace std;


//============================================================================//



//std::random_device generator;
//std::default_random_engine generator;


int main(int argc, char *argv[])
{
   //..................................................
   //               input variables
   //..................................................
   int Nsite,Nspin,Ntau;
   int Norder,Nmeas,Nshift;
   int PrintTime,binlength,verbosity,muIter;
   double Beta,mu,muStep,muErr,density;
   bool retarded;
   std::vector<int> SiteTime;
   std::vector<int> SiteNorb;
   std::vector<std::string> SiteName;
   std::vector<std::string> SiteDir;
   char* IterationDir;


   //..................................................
   //                 start global timer
   //..................................................
   std::chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;
   print_line_star(80);
   start_tot = std::chrono::system_clock::now();


   //..................................................
   //                      main code
   //..................................................
   try
   {
      //..................................................
      //                read inputfile
      //..................................................
      if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");
      //
      // iteration folder
      IterationDir = argv[2];
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
      find_param(argv[1], "muErr"      , muErr     );
      //
      if(verbosity>=1)
      {
         std::cout << "mu= " << mu << std::endl;
         std::cout << "beta= " << Beta << std::endl;
         std::cout << "Nsite= " << Nsite << std::endl;
         std::cout << "Nspin= " << Nspin << std::endl;
         std::cout << "Ntau= " << Ntau << std::endl;
         std::cout << "Norder= " << Norder << std::endl;
         std::cout << "Nmeas= " << Nmeas << std::endl;
         std::cout << "Nshift= " << Nshift << std::endl;
         std::cout << "Nshift= " << Nshift << std::endl;
         std::cout << "PrintTime= " << PrintTime << std::endl;
         std::cout << "binlength= " << binlength << std::endl;
         std::cout << "retarded= " << retarded << std::endl;
         std::cout << "muStep= " << muStep << std::endl;
         std::cout << "muIter= " << muIter << std::endl;
         std::cout << "muErr= " << muErr << std::endl;
         std::cout << "density= " << density << std::endl;
      }

      //
      for(int isite=0; isite < Nsite; isite++)
      {
         //
         std::string ss=std::to_string(isite+1);
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



      //..................................................
      //                    Solver setup
      //..................................................
      std::vector<ct_hyb> ImpurityList;
      for(int isite=0; isite < Nsite; isite++)
      {
         std::cout << " \nisite = " << isite << std::endl;
         std::cout << " Name = " << SiteName[isite] << std::endl;
         std::cout << " Norb = " << SiteNorb[isite] << std::endl;
         std::cout << " Time = " << SiteTime[isite] << std::endl;
         std::cout << " Folder = " << SiteDir[isite];

         if(PathExist(strcpy(new char[SiteDir[isite].length() + 1], SiteDir[isite].c_str())))
         {
            std::cout << " (Found) - ";
            ImpurityList.push_back(ct_hyb( mu, Beta, Nspin, SiteNorb[isite], Ntau, Norder, Nmeas, Nshift, PrintTime, binlength, retarded ));
            ImpurityList[isite].init( SiteDir[isite] );
         }
         else
         {
            std::cout << " (Not Found) - Exiting." << std::endl;
            exit(1);
         }
         std::cout << std::endl;
      }
      print_line_space(1);



      //..................................................
      //             Chemical potential Lookup
      //..................................................
      if(density>0.0)
      {
         int QuickTime=120;
         double trial_density;
         double mu_start=mu;
         double mu_new;
         std::vector<double>Ntmp(Nsite,0.0);

         //
         print_line_minus(80);
         std::cout << " Starting chemical potential: " << mu_start << std::endl;
         for(int isite=0; isite < Nsite; isite++)
         {
            std::cout << " Quick solution (120sec) of site: " << SiteName[isite] << std::endl;
            ImpurityList[isite].solve( QuickTime, true );
            Ntmp[isite]=ImpurityList[isite].get_Density();
            std::cout << " Site density: " << Ntmp[isite] << std::endl;
         }
         trial_density = std::accumulate(Ntmp.begin(), Ntmp.end(), 0.0);
         std::cout << " Total density: " << trial_density << std::endl;
         print_line_space(1);

         //
         print_line_minus(80);
         double muSign = (trial_density < density) ? +1.0 : -1.0;
         for (int imu=0; imu < muIter; imu++)
         {
            //
            mu_new = mu_start + imu*muStep*muSign;
            std::cout << " Setting chemical potential: " << mu_new << std::endl;
            //
            for(int isite=0; isite < Nsite; isite++)
            {
               std::cout << " Quick solution (120sec) of site: " << SiteName[isite] << std::endl;
               ImpurityList[isite].reset_mu( mu_new );
               ImpurityList[isite].solve( QuickTime, true );
               Ntmp[isite]=ImpurityList[isite].get_Density();
               std::cout << " Site density: " << Ntmp[isite] << std::endl;
            }
            trial_density = std::accumulate(Ntmp.begin(), Ntmp.end(), 0.0);
            std::cout << " Total density: " << trial_density << std::endl;
            print_line_space(1);
            //
            if((muSign>0.0)&&(trial_density > density)) break;
            if((muSign<0.0)&&(trial_density < density)) break;
         }

         //
         print_line_minus(80);
         double mu_below=std::min(mu_new,mu_start);
         double mu_above=std::max(mu_new,mu_start);
         for (int imu=0; imu < muIter; imu++)
         {
            //
            std::cout << " Chemical potential boundaries: " << mu_below  << " - "  << mu_above << std::endl;
            mu_new = (mu_below+mu_above)/2.0;
            std::cout << " Setting chemical potential: " << mu_new << std::endl;
            //
            for(int isite=0; isite < Nsite; isite++)
            {
               std::cout << " Quick solution (120sec) of site: " << SiteName[isite] << std::endl;
               ImpurityList[isite].reset_mu( mu_new );
               ImpurityList[isite].solve( QuickTime, true );
               Ntmp[isite]=ImpurityList[isite].get_Density();
               std::cout << " Site density: " << Ntmp[isite] << std::endl;
            }
            trial_density = std::accumulate(Ntmp.begin(), Ntmp.end(), 0.0);
            std::cout << " Total density: " << trial_density << " Error: " << fabs(trial_density-density) << std::endl;
            print_line_space(1);
            //
            if(trial_density > density) mu_above=mu_new;
            if(trial_density < density) mu_below=mu_new;
            if(fabs(trial_density-density)<muErr)
            {
               std::cout << " Found correct chemical potential: " << mu_new << std::endl;
               break;
            }
            else if(imu == muIter-1)
            {
               std::cout << " *NOT* found correct chemical potential. Last value: " << mu_new << std::endl;
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
      end_tot = std::chrono::system_clock::now();
      std::chrono::duration<double> runtime_seconds = end_tot-start_tot;
      cout << endl;
      cout << "Time [total] = " << runtime_seconds.count() << "s\n";
      print_line_star(60);

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
