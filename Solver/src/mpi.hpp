#ifndef ___MPI___
#define ___MPI___
//
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
#include <unistd.h>
#include "file_io.hpp"


//============================================================================//


class CustomMPI
{
   public:

      //----------------------------------------------------------------------//

      int rank(void)const {return MPIrank;}
      int size(void)const {return MPIsize;}
      int master(void)const {return MPImaster;}
      bool is_master(void)const {return (MPIrank == MPImaster) ? true : false;}
      void finalize(void)const { MPI_Finalize(); exit(1); }

      //----------------------------------------------------------------------//

      void init(int master_=0)
      {
         MPI_Init(NULL, NULL);
         MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
         MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
         if(master_!=0) MPImaster=master_;
         //
         print_line_star(80,is_master());
         sleep(1);
         printf(" Rank %d out of %d available processors is alive.\n",MPIrank, MPIsize);
         if(is_master())printf(" Master is %d. \n",MPImaster);
         if((is_master())&&(MPIverbose))printf(" Verbosity ON.\n");
         sleep(1);
         print_line_star(80,is_master());
         print_line_space(1,is_master());
         //
         MPI_is_init=true;
      }

      //----------------------------------------------------------------------//

      void report(std::string message)
      {
         unsigned int microseconds=10000*MPIrank;
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         if(!MPIverbose)
         {
            if(MPIrank==MPImaster) printf("%s\n", message.c_str());
         }
         else
         {
            printf(" [Rank %d]: %s\n",MPIrank, message.c_str());
            usleep(microseconds);
         }

      }

      //----------------------------------------------------------------------//

      void barrier(void)
      {
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         MPI_Barrier(MPI_COMM_WORLD);
      }

      //----------------------------------------------------------------------//

      void StopError(std::string message)
      {
         report(message);
         finalize();
      }

      //----------------------------------------------------------------------//

      void allreduce( unsigned long long int &RankInt, unsigned long long int &WorldInt, bool average=false, path location="none")
      {
         if(location!="none")report("allreduce(int) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize==1)
         {
            WorldInt = RankInt;
         }
         else
         {
            MPI_Allreduce(&RankInt, &WorldInt, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            if(average)WorldInt/=MPIsize;
            MPI_Barrier(MPI_COMM_WORLD);
         }
      }

      void allreduce( double &RankFloat, double &WorldFloat, bool average=false, path location="none")
      {
         if(location!="none")report("allreduce(float) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize==1)
         {
            WorldFloat = RankFloat;
         }
         else
         {
            MPI_Allreduce(&RankFloat, &WorldFloat, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if(average) WorldFloat/=MPIsize;
            MPI_Barrier(MPI_COMM_WORLD);
         }
      }

      void allreduce( std::vector<double> &RankVec, std::vector<double> &WorldVec, bool average=false, path location="none")
      {
         if(location!="none")report("allreduce(Vec) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize==1)
         {
            WorldVec = RankVec;
         }
         else
         {
            int iDim = WorldVec.size();
            for(int i=0; i < iDim; i++)
            {
               MPI_Allreduce(&RankVec[i], &WorldVec[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
               if(average) WorldVec[i]/=MPIsize;
               MPI_Barrier(MPI_COMM_WORLD);
            }
         }
      }

      void allreduce( std::vector<std::vector<double>> &RankVecVec, std::vector<std::vector<double>> &WorldVecVec, bool average=false, path location="none")
      {
         if(location!="none")report("allreduce(VecVec) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize==1)
         {
            WorldVecVec = RankVecVec;
         }
         else
         {
            int iDim = WorldVecVec.size();
            int jDim = WorldVecVec[0].size();
            for(int i=0; i < iDim; i++)
            {
               for(int j=0; j < jDim; j++)
               {
                  MPI_Allreduce(&RankVecVec[i][j], &WorldVecVec[i][j], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  if(average) WorldVecVec[i][j]/=MPIsize;
                  MPI_Barrier(MPI_COMM_WORLD);
               }
            }
         }
      }

      //----------------------------------------------------------------------//

      void broadcast(bool &RankBool, int SendingRank, path location="none")
      {
         if(location!="none")report("broadcast(bool) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize>=2)
         {
            MPI_Bcast( &RankBool, 1, MPI_C_BOOL, SendingRank, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
         }
      }

      void broadcast(int &RankInt, int SendingRank, path location="none")
      {
         if(location!="none")report("broadcast(int) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize>=2)
         {
            MPI_Bcast( &RankInt, 1, MPI_INT, SendingRank, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
         }
      }

      void broadcast(double &RankFloat, int SendingRank, path location="none")
      {
         if(location!="none")report("broadcast(double) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize>=2)
         {
            MPI_Bcast( &RankFloat, 1, MPI_DOUBLE, SendingRank, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
         }
      }

      void broadcast(std::vector<int> &Vec, int SendingRank, path location="none")
      {
         if(location!="none")report("broadcast(Vec<int>) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize>=2)
         {
            int iDim = Vec.size();
            for(int i=0; i < iDim; i++)
            {
               MPI_Bcast( &Vec[i], 1, MPI_INT, SendingRank, MPI_COMM_WORLD);
               MPI_Barrier(MPI_COMM_WORLD);
            }
         }
      }

      void broadcast(std::vector<double> &Vec, int SendingRank, path location="none")
      {
         if(location!="none")report("broadcast(Vec<double>) called in: "+location);
         if(!MPI_is_init) StopError(" MPI environment is not initialized. Exiting.");
         //
         if(MPIsize>=2)
         {
            int iDim = Vec.size();
            for(int i=0; i < iDim; i++)
            {
               MPI_Bcast( &Vec[i], 1, MPI_DOUBLE, SendingRank, MPI_COMM_WORLD);
               MPI_Barrier(MPI_COMM_WORLD);
            }
         }
      }

      //----------------------------------------------------------------------//

   //private:

      //----------------------------------------------------------------------//

      int                                 MPIsize;
      int                                 MPIrank;
      int                                 MPImaster=0;
      bool                                MPI_is_init=false;
      #ifdef _verb
         bool MPIverbose=true;
      #else
         bool MPIverbose=false;
      #endif

      //----------------------------------------------------------------------//

};


//============================================================================//


#endif
