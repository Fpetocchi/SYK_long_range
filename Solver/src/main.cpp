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
   int Ntau;
   int Norb;
   int Nspin;
   int Nsite;
   std::vector<int> SiteTime;
   std::vector<std::string> SiteName;
   std::vector<std::string> SiteDir;
   char* IterationDir;


   //..................................................
   //                 start global timer
   //..................................................
   std::chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;
   print_line_star(60);
   cout << "            Test progam: QMC ct-hyb" << endl;
   print_line_star(60);
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
      find_param(argv[1], "Nspin"  , Nspin  );
      find_param(argv[1], "Norb"   , Norb   );
      find_param(argv[1], "Ntau"   , Ntau   );
      find_param(argv[1], "Nsite"  , Nsite  );
      //
      std::cout << "Nspin= " << Nspin << std::endl;
      std::cout << "Norb= " << Norb << std::endl;
      std::cout << "Ntau= " << Ntau << std::endl;
      std::cout << "Nsite= " << Nsite << std::endl;
      //
      for(int isite=0; isite < Nsite; isite++)
      {
         //
         std::string ss=std::to_string(isite+1);
         const char * s  = ss.c_str();
         char lineN[5];
         char lineT[5];
         char folder[999];
         char element[2];
         int minutes;
         //
         strcpy(lineT,"Time");
         strcat(lineT,s);
         find_param(argv[1], lineT, minutes );
         SiteTime.push_back(minutes);
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

      for(int isite=0; isite < Nsite; isite++)
      {
         std::cout << " \nisite = " << isite << std::endl;
         std::cout << " Name = " << SiteName[isite] << std::endl;
         std::cout << " Time = " << SiteTime[isite] << std::endl;
         std::cout << " Folder = " << SiteDir[isite];

         if(PathExist(strcpy(new char[SiteDir[isite].length() + 1], SiteDir[isite].c_str())))
         {
            std::cout << " (Found) ";
         }
         else
         {
            std::cout << " (Not Found) - Exiting." << std::endl;
            exit(1);
         }
         std::cout << std::endl;
      }




      //..................................................
      //           SETUP THE RANDOM GENERATOR SEED
      //..................................................
      unsigned long seed = std::chrono::duration_cast<std::chrono::milliseconds>
      (std::chrono::system_clock::now().time_since_epoch()).count()-1566563000000;
      generator.seed(abs(seed));




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
