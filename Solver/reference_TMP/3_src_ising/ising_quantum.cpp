#include <ctime>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <mpi.h>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "quantumMC.hpp"

using namespace std;

//-------------------------------------------------------------------------

float randfloat()
{
  return rand()/(float(RAND_MAX)+1);
}

//-------------------------------------------------------------------------

template <typename T>
std::string to_string_p(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

//-------------------------------------------------------------------------

std::vector<std::complex<double>> invert_motionEq( std::complex<double> Eloc
                                                 , std::vector<double> expC
                                                 , std::vector<std::complex<double>> expEC
                                                 , std::vector<std::vector<double>> expCC
                                                 , bool verbose )
{
   //
   int directions = expC.size();
   std::vector<std::complex<double>> alphadot_vec(expC.size(),std::complex<double> (0.0,0.0));

   //
   Eigen::MatrixXcd A(directions,directions); A.setZero(directions,directions);
   Eigen::VectorXcd B(directions); B.setZero(directions);
   Eigen::VectorXcd adot(directions); adot.setZero(directions);

   //
   int idirndx=-1;
   for (int idir = 0; idir < expC.size(); idir++)
   {
      //
      //if ((idir!=2))//&&(idir!=4))
      //{
      //   idirndx++;
      //   B(idirndx) = expEC[idir] - Eloc*expC[idir];
      //}
      B(idir) = expEC[idir] - Eloc*expC[idir];
      //int jdirndx=-1;
      for (int jdir = 0; jdir < expC.size(); jdir++)
      {
         A(idir,jdir) = expCC[idir][jdir] - expC[idir]*expC[jdir];
         //if ((jdir!=2))//&&(jdir!=4))
         //{
         //   jdirndx++;
         //   A(idirndx,jdirndx) = expCC[idir][jdir] - expC[idir]*expC[jdir];
         //   //std::cout << " " << idirndx << " " << jdirndx << std::endl;
         //}
      }
   }

   //
   if( verbose == true )
   {
      printf("\n--------------\n");
      std::cout << "Here is B:\n" << B << std::endl;
      std::cout << "Here is A:\n" << A << std::endl;
      std::cout << "Here is A^-1:\n" << A.inverse() << std::endl;
      std::cout << "Here is A*A^-1:\n" << A*A.inverse() << std::endl;
      std::cout << "Here is A^-1*B:\n" << A.inverse()* B << std::endl;
      printf("--------------\n\n");
   }

   //
   adot = A.inverse() * B;

   //
   //alphadot_vec[0]=adot[0];
   //alphadot_vec[1]=adot[1];
   //alphadot_vec[2]=adot[2];
   //alphadot_vec[3]=adot[2];
   //alphadot_vec[4]=adot[2];
   for (int idir = 0; idir < expC.size(); idir++)alphadot_vec[idir] = adot[idir];
   return alphadot_vec;
}

//-------------------------------------------------------------------------

std::vector<std::complex<double>> ode_RK4( double h, std::vector<std::complex<double>> alphadot
                                                   , std::vector<std::complex<double>> alpha_old )
{
   //
   int directions = alphadot.size();
   std::vector<std::complex<double>> alpha_new(directions,std::complex<double> (0.0,0.0));

   // just not to modify too much
   int n = directions;
   std::vector<std::complex<double>> y0(n),dy(n);
   std::vector<std::complex<double>> y1(n),k1(n),k2(n),k3(n),k4(n);
   std::complex<double> img(0.0,1.0);

   //
   for (int idir = 0; idir < directions; idir++)
   {
      y0[idir] = alpha_old[idir];
      dy[idir] = -img*alphadot[idir];
   }

   //
   // K1 = h * f( tn , Sn )
   //dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<n;i++)
   {
      k1[i] = h*dy[i];
      y1[i] = y0[i] + 0.5*k1[i];
   }
   //dy = ode.deriv(t+0.5*h,y1,Jvec);
   for(int i=0;i<n;i++)
   {
      k2[i] = h*dy[i];
      y1[i] = y0[i] + 0.5*k2[i];
   }
   //dy = ode.deriv(t+0.5*h,y1,Jvec);
   for(int i=0;i<n;i++)
   {
      k3[i] = h*dy[i];
      y1[i] = y0[i] + k3[i];
   }
   //dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<n;i++)
   {
      k4[i] = h*dy[i];
      y0[i] = y0[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
   }

   //
   for (int idir = 0; idir < directions; idir++) alpha_new[idir] = y0[idir];
   return alpha_new;
}

//-------------------------------------------------------------------------

std::vector<std::complex<double>> ode_Eul( double h, std::vector<std::complex<double>> alphadot
                                                   , std::vector<std::complex<double>> alpha_old )
{
   //
   int directions = alphadot.size();
   std::vector<std::complex<double>> alpha_new(directions,std::complex<double> (0.0,0.0));
   std::complex<double> img(0.0,1.0);
   //
   for (int idir = 0; idir < directions; idir++) alpha_new[idir] = -img*alphadot[idir]*h + alpha_old[idir];
   return alpha_new;
}

//-------------------------------------------------------------------------

double Jramp(double t, double J_ini, double J_fin, double T_o, double speed)
{
   double J;
   double offset = (J_fin - J_ini)/2.0;
   J = J_ini + offset * (tanh(2.0*speed*( t-T_o-3.0/(2.0*speed) ))+1.0);
   return J;
}


//===========================================================================
//===========================================================================


int main()
{
   // Initialize the MPI environment
   MPI_Init(NULL, NULL);
   // Get the number of processes
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   printf("rank %d out of %d processors\n",world_rank, world_size);

   // input vars
   int nx=4, ny=4;

   // time evolution vars
   int it;
   int tmax=1000;
   double h=0.005;

   // Hamiltonian parameters
   double Hfield=5.0;
   double Jstart=1.0;
   double Jend=1.0;
   double Jspeed=4.0;
   int Jswitch=1005;

   bool testing=false;
   bool verbose=true;
   int progress_resolution = 25;

   // relaxation after initialization
   unsigned long int max_MEASsteps=30000;
   unsigned long int max_THRMsteps=10000;
   unsigned long int max_VMCsteps=max_MEASsteps+max_THRMsteps;

   // random generator setup //srand((unsigned)time(0));
   srand((unsigned)time(NULL)+world_rank*world_size + 7);

   // setup class
   QuantumIsing_2D spinsys(nx, ny, verbose, testing);
   spinsys.set_directions();
   MPI_Barrier(MPI_COMM_WORLD);

   // directions and conserved operators
   int directions = spinsys.get_directions();
   std::vector<std::vector<double>> Coperators = spinsys.gen_operators();

   // Hilbert space various maps
   std::vector<std::vector<int>> Hilbert_bin = spinsys.gen_Hilbert_bin(false);
   std::vector<std::vector<int>> Hilbert_Sz = spinsys.gen_Hilbert_bin(true);
   std::vector<std::vector<int>> Flip_spin = spinsys.gen_Flips();
   std::vector<double> Diag_H = spinsys.gen_Hamiltonian_diag();
   //std::vector<std::vector<double>> Diag_H = spinsys.gen_Hamiltonian_diag();
   //std::vector<std::vector<double>> OffDiag_H = spinsys.gen_Hamiltonian_offdiag(Flip_spin);
   std::vector<double> Sz = spinsys.gen_Sz();

   // Jastrow factors
   std::vector<std::complex<double>> coeff;
   std::vector<std::vector<std::complex<double>>> alpha_vec( tmax , std::vector<std::complex<double>> (directions, std::complex<double> (0.0,0.0)));
   std::vector<std::vector<std::complex<double>>> alphadot_vec( tmax , std::vector<std::complex<double>> (directions, std::complex<double> (0.0,0.0)));

   // TESTING
   if (testing==true && world_rank==0)
   {
       for (int istate = 0; istate < 20 /*spinsys.get_Hil()*/; istate++)
       {
           printf(" state: %i  bin:  ",istate);
           for (int ilatt = 0; ilatt < spinsys.get_dim(); ilatt++) printf("  %i  ",Hilbert_bin[istate][ilatt]);
           printf(" spin:    ");
           for (int ilatt = 0; ilatt < spinsys.get_dim(); ilatt++) printf("  % i  ",Hilbert_Sz[istate][ilatt]);
           printf(" <Sz>:  % 2.4f  ",Sz[istate]);
           printf(" Oper:    ");
           for (int idir = 0; idir < directions; idir++) printf("  % 2.4f  ",Coperators[istate][idir]);
           //printf(" conn:    ");
           for (int ilatt = 0; ilatt < spinsys.get_dim(); ilatt++) printf("  %i  ",Flip_spin[istate][ilatt]);
           printf("\n");
       }
       for (int istate = 0; istate < spinsys.get_Hil(); istate++)
       {
           printf(" neig of state : %i \n",istate);
           for (int ilatt = 0; ilatt < spinsys.get_dim(); ilatt++)
           {
               std::vector<int> neigh = spinsys.get_neigh(ilatt, Hilbert_Sz[istate]);
               printf(" ilat: %i  ",spinsys.get_dim()-1-ilatt);
               for (int ineig = 0; ineig < 4; ineig++) printf("  % i  ",neigh[ineig]);
               printf("\n");
           }
       }
   }

   // Time evolution
   for (int it = 0; it < tmax; it++)
   {
       double Jfield=Jstart; // = Jramp( it*h, Jstart, Jend, h*Jswitch, Jspeed);
       if (it>=Jswitch) Jfield=Jend;

       if (world_rank==0) printf("\n\n--------- it: %5i - Jfield: %2.2f - Hfield: %2.2f ---------\n",it,Jfield,Hfield);

       // Progress variables
       int progress = progress_resolution;

       // Measured quantities
       int accepted = 0;
       std::complex<double> Eloc_p = 0;
       std::complex<double> Eloc = 0;
       double Mloc_p = 0;
       double Mloc = 0;
       double MlocAv = 0;
       std::vector<double> expC_p(directions,0.0), expC(directions,0.0);
       std::vector<std::complex<double>> expEC_p(directions,0.0), expEC(directions,0.0);
       std::vector<std::vector<double>> expCC_p( directions , std::vector<double> (directions, 0.0));
       std::vector<std::vector<double>> expCC( directions , std::vector<double> (directions, 0.0));

       // Initialize all coefficient given the operators(fixed for each configuration) and the alpha_vec of the particular timestep
       coeff = spinsys.init_alpha_coeff(alpha_vec[it], Coperators );

       // Initial state
       //int istate = std::distance(coeff.begin(),std::max_element(coeff.begin(), coeff.end())); //170//23130
       int istate = 5; //(int) (spinsys.get_Hil()/2);

       // Montercarlo steps
       for (int k = 0; k < (max_THRMsteps+max_MEASsteps) ; k++)
       {
          //progress
          if ((int)(100*(double)k/max_VMCsteps)>progress)
          {
             printf("Processor: %i - MC Progress - %i %% \n",world_rank,progress);
             progress += progress_resolution;
          }

          if(k>=max_THRMsteps)
          {
             // 1) Local energy evaluation (diagonal + non-diagonal)
             std::complex<double> Etmp = Jfield * Diag_H[istate];
             for (int ilatt = 0; ilatt < spinsys.get_dim(); ilatt++) Etmp -= coeff[Flip_spin[istate][ilatt]] * Hfield / coeff[istate];
             Eloc_p += Etmp / ((double)(nx*ny));

             // 2) <C>, <Eloc,C>, <CC>
             for (int idir = 0; idir < directions; idir++)
             {
                expC_p[idir]  += Coperators[istate][idir];
                expEC_p[idir] += Eloc_p*Coperators[istate][idir];
                for (int jdir = 0; jdir < directions; jdir++) expCC_p[idir][jdir]  += Coperators[istate][idir]*Coperators[istate][jdir];
             }

             // 3) Magnetization(diagonal)
             Mloc_p += Sz[istate];
          }


          // random choice of the spin to flip
          int ilatt = rand()%spinsys.get_dim();

          // new coefficient and state given by the flip
          int istate_new = Flip_spin[istate][ilatt];
          std::complex<double> coeff_used = coeff[istate];
          std::complex<double> coeff_new  = coeff[istate_new];

          // Metropolis
          double rf = randfloat();
          double acceptance = 0.0; //pow(coeff_new,2)/pow(coeff[istate],2)
          for (int idir = 0; idir < directions; idir++) acceptance += 2.0*real(alpha_vec[it][idir])*(Coperators[istate_new][idir]-Coperators[istate][idir]);
          if( rf < min( 1. , exp(acceptance) ) )
          {
              istate = istate_new;
              accepted++;
          }
       }

       //Eloc_p = 2 + Hfield + Eloc_p / (double)max_VMCsteps;
       Eloc_p = Eloc_p / (double)max_MEASsteps;
       Mloc_p = Mloc_p / (double)max_MEASsteps;
       for (int idir = 0; idir < directions; idir++)
       {
          expC_p[idir]  /= (double)max_MEASsteps;
          expEC_p[idir] /= (double)max_MEASsteps;
          for (int jdir = 0; jdir < directions; jdir++) expCC_p[idir][jdir]  /= (double)max_MEASsteps;
       }

       if (world_rank==0)printf("\n");
       MPI_Barrier(MPI_COMM_WORLD);
       printf("Processor: %i - Re(Eloc): %4.7f - Im(Eloc): %4.7f - Magn: %4.7f - accepted: %10i\n",world_rank,real(Eloc_p),imag(Eloc_p),Mloc_p,accepted);

       // Gather the Elocs from all the procs and average
       MPI_Allreduce(&Eloc_p, &Eloc, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
       Eloc = Eloc / (double) world_size;

       // Gather the Mlocs from all the procs and average
       MPI_Allreduce(&Mloc_p, &Mloc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
       Mloc = Mloc / (double) world_size;
       MlocAv += sqrt(Mloc*Mloc)*h;

       // Gather the directional operators and average
       for (int idir = 0; idir < directions; idir++)
       {
          MPI_Allreduce(&expC_p[idir], &expC[idir], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          expC[idir] = expC[idir] / (double) world_size;
          MPI_Allreduce(&expEC_p[idir], &expEC[idir], 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
          expEC[idir] = expEC[idir] / (double) world_size;
          for (int jdir = 0; jdir < directions; jdir++)
          {
             MPI_Allreduce(&expCC_p[idir][jdir], &expCC[idir][jdir], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
             expCC[idir][jdir] = expCC[idir][jdir] / (double) world_size;
          }
       }

       // Build the alphadot vector
       alphadot_vec[it] = invert_motionEq(Eloc,expC,expEC,expCC,verbose&&(world_rank==0));

       // Solve eq of motion with standard integration
       alpha_vec[it+1] = ode_RK4(h,alphadot_vec[it],alpha_vec[it]);
       //alpha_vec[it+1] = ode_Eul(h,alphadot_vec[it],alpha_vec[it]);

       //
       for (int idir = 0; idir < directions; idir++)
       {
          printf("Processor: %i - idir: %i, alphadot: ( %f , %f ), alpha_new: ( %f , %f )\n",world_rank,idir
          ,real(alphadot_vec[it][idir]),imag(alphadot_vec[it][idir]),real(alpha_vec[it+1][idir]),imag(alpha_vec[it+1][idir]));
       }
       MPI_Barrier(MPI_COMM_WORLD);

       //
       MPI_Barrier(MPI_COMM_WORLD);
       printf("\nUsed  - Re(Eloc): %4.7f - Im(Eloc): %4.7f - Magn: %4.7f  - AvMagn: %4.7f \n",real(Eloc),imag(Eloc),Mloc,MlocAv);
       if (world_rank==0)
       {
           FILE *fobservables=fopen("Observables.dat","a");
           fprintf (fobservables , "%6i  ",it);
           fprintf (fobservables , "% 4.4f  % 4.4f  ",Hfield,Jfield);
           fprintf (fobservables , "% 4.4f  % 4.4f  ",real(Eloc),imag(Eloc));
           fprintf (fobservables , "% 4.4f  % 4.4f  ",Mloc,MlocAv);
           for (int idir = 0; idir < directions; idir++)fprintf (fobservables , "% 4.4f  ",real(alpha_vec[it][idir]));
           for (int idir = 0; idir < directions; idir++)fprintf (fobservables , "% 4.4f  ",imag(alpha_vec[it][idir]));
           for (int idir = 0; idir < directions; idir++)fprintf (fobservables , "% 4.4f  ",real(alphadot_vec[it][idir]));
           for (int idir = 0; idir < directions; idir++)fprintf (fobservables , "% 4.4f  ",imag(alphadot_vec[it][idir]));
           for (int idir = 0; idir < directions; idir++)fprintf (fobservables , "% 4.4f  ",expC[idir]);
           for (int idir = 0; idir < directions; idir++)fprintf (fobservables , "% 4.4f  ",real(expEC[idir]));
           for (int idir = 0; idir < directions; idir++)fprintf (fobservables , "% 4.4f  ",imag(expEC[idir]));
           fprintf (fobservables , "\n");
           fclose(fobservables);
       }
       MPI_Barrier(MPI_COMM_WORLD);

   } // Loop on time

   // Finalize the MPI environment.
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   //
   exit(1);
   return 0;

}
