#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <math.h>
#include <time.h>
#include <fstream>
#include <limits>



#include "odesolvers.hpp"
#include "formats.hpp"
#include "read_inputfile.hpp"

#define PI 3.14159265

using namespace std;

//-------------------------------------------------------------------------

//std::random_device generator;
std::default_random_engine generator;


//-------------------------------------------------------------------------


std::fstream& GotoLine(std::fstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i)
    {
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


//-------------------------------------------------------------------------

double compute_polar(int i, int j, vector<double> spins)
{
   double component;
   std::vector<double> Polar(3);
   double rho = sqrt(pow(spins[3*i+0],2) + pow(spins[3*i+1],2) + pow(spins[3*i+2],2));
   double theta = acos(spins[3*i+2]/rho);
   double phi = atan2 (spins[3*i+1],spins[3*i+0]);
   //
   Polar[0] = rho;
   Polar[1] = theta;
   Polar[2] = phi;
   //
   component = Polar[j];
   //
   return component;
}

//-------------------------------------------------------------------------

double gen_Temp(double t, double Temp_ini, double Temp_fin, double To, double speed)
{
   double Temp;
   double offset = (Temp_fin - Temp_ini)/2.0;
   Temp = Temp_ini + offset * (tanh(2.0*speed*( t-To-3.0/(2.0*speed) ))+1.0);
   return Temp;
}

double gen_Jex(double t, double Jex_ini, double Jex_fin, double J_o, double speed)
{
   double Jex;
   double offset = (Jex_fin - Jex_ini)/2.0;
   Jex = Jex_ini + offset * (tanh(2.0*speed*( t-J_o-3.0/(2.0*speed) ))+1.0);
   return Jex;
}

double gen_Jey(double t, double Jey_ini, double Jey_fin, double J_o, double speed)
{
   double Jey;
   double offset = (Jey_fin - Jey_ini)/2.0;
   Jey = Jey_ini + offset * (tanh(2.0*speed*( t-J_o-3.0/(2.0*speed) ))+1.0);
   return Jey;
}

double gen_Jez(double t, double Jez_ini, double Jez_fin, double J_o, double speed)
{
   double Jez;
   double offset = (Jez_fin - Jez_ini)/2.0;
   Jez = Jez_ini + offset * (tanh(2.0*speed*( t-J_o-3.0/(2.0*speed) ))+1.0);
   return Jez;
}

//===========================================================================

class spin2D_t
{
public:
   spin2D_t (int nx, int ny, double damping):
             nx(nx), ny(ny), damping(damping){}
   int get_nx(void)const {return nx;}
   int get_ny(void)const {return ny;}
   int get_dim(void)const {return nx*ny;}
   double get_damping(void)const {return abs(damping);}

   //-------------------------------------------------------------------------

   void gen_latt()
   {
      coords.resize(2*nx*ny);
      int i=0;
      for(int ix=0;ix<nx;ix++)
      {
         for(int iy=0;iy<ny;iy++)
         {
            coords[2*i] = ix+1;
            coords[2*i+1] = iy+1;
            i++;
         }
      }
   }

   //-------------------------------------------------------------------------

   void gen_neigh()
   {
      neigh.resize(4*nx*ny);
      int i=0;
      for (int i=0; i<nx*ny; i++)
      {
         int k=coords[2*i];
         int k1=coords[2*i]-1;
         int k2=coords[2*i]+1;
         int l=coords[2*i+1];
         int l1=coords[2*i+1]-1;
         int l2=coords[2*i+1]+1;
         //
         if (k1 < 1)  k1 = nx;
         if (k2>nx)   k2=1;
         if (l1 < 1)  l1=ny;
         if (l2 > ny) l2=1;
         //
         // -X
         neigh[4*i+0]=ny*(k1-1)+l-1;
         // -Y
         neigh[4*i+1]=ny*(k-1)+l1-1;
         // +X
         neigh[4*i+2]=ny*(k-1)+l2-1;
         // +Y
         neigh[4*i+3]=ny*(k2-1)+l-1;
         //
      }
   }

   //-------------------------------------------------------------------------

   std::vector<double> get_spin(int i, vector<double> spins)
   {
      std::vector<double> spin_i(3);
      for (int j=0; j<3; j++) spin_i[j] = spins[3*i+j];
      return spin_i;
   }

   //-------------------------------------------------------------------------

   std::vector<double> gen_beff(int i, vector<double> spins, vector<double> Jvec)
   {
      std::vector<double> beff(3);
      std::vector<double> s1(3),s2(3),s3(3),s4(3);
      //
      //loop over the three directions
      for (int j=0; j<3; j++)
      {
         s1[j] = spins[3*neigh[4*i+0]+j];
         s2[j] = spins[3*neigh[4*i+1]+j];
         s3[j] = spins[3*neigh[4*i+2]+j];
         s4[j] = spins[3*neigh[4*i+3]+j];
      }
      for (int j=0; j<3; j++) beff[j] = Jvec[j]*(s1[j]+s2[j]+s3[j]+s4[j]);
      return beff;
   }

   //-------------------------------------------------------------------------
   std::vector<double>  gen_vrtx(int i, vector<double> spins)
   {
      std::vector<double> vrtx(3);
      int plusX = neigh[4*i+3];
      int plusY = neigh[4*i+2];
      // vrtx_x =   (  Sz[i+Y] -  Sz[i]  )
      // vrtx_y = - (  Sz[i+X] -  Sz[i]  )
      // vrtx_z =   (  Sy[i+X] -  Sy[i]  ) - (  Sx[i+Y] -  Sx[i]  )
      vrtx[0] = + ( spins[3*plusY+2]-spins[3*i+2] );
      vrtx[1] = - ( spins[3*plusX+2]-spins[3*i+2] );
      vrtx[2] = + ( spins[3*plusX+1]-spins[3*i+1] ) - ( spins[3*plusY+0]-spins[3*i+0] );
      return vrtx;
   }

   //-------------------------------------------------------------------------

   std::vector<double> deriv(double t, vector<double> spins, vector<double> Jvec)
   {
      std::vector<double> ds(3*nx*ny);
      std::vector<double> beff_i(3),spin_i(3);
      //
      for(int i=0; i < nx*ny; i++)
      {
         //
         beff_i = gen_beff(i,spins,Jvec);
         spin_i = get_spin(i,spins);
         //
         std::vector<double> BxS = cross(beff_i, spin_i);
         std::vector<double> gilb    = cross(spin_i, BxS);
         //
         for (int j=0; j<3; j++) ds[3*i+j] = BxS[j] - damping * gilb[j];
         //
      }
      return ds;
   }

   //-------------------------------------------------------------------------

   void timestep(int tstp, double h, vector<double> &spins, vector<double> Jvec, double Thermostat=0.0 , int renorm_M=1)
   {
     double t1 = tstp*h;
     if(Thermostat==0.0)
     {
        step_rungekutta4<spin2D_t>(t1, h, spins, *this, Jvec);
     }
     else
     {
        step_Ito_alphaInvariant(t1, h, spins, *this, Jvec, Thermostat, renorm_M);
     }
   }

   //-------------------------------------------------------------------------

   double flat_distribution(double m2)
   {
      double a = sqrt(3.*m2);
      std::uniform_real_distribution<double> distribution_flat(-a,+a);
      double r = distribution_flat(generator);
      return r;
   }

   double flat_distribution_bis(double m2)
   {
      std::uniform_real_distribution<double> distribution_flat(0.0,m2);
      double r = distribution_flat(generator);
      return r;
   }

   double gauss_distribution(double m2)
   {
      //std::uniform_real_distribution<double> distribution_flat(0,m2);
      //double p = distribution_flat(generator);
      std::normal_distribution<double> distribution_norm(0.0,m2);
      double r = distribution_norm(generator);
      /*
      FILE *fdist = fopen("gauss_dist.dat", "a");
      fprintf(fdist,"%g\n",r);
      fclose(fdist);
      */
      return r;
   }

   //-------------------------------------------------------------------------

   int nx;
   int ny;
   double damping;
   vector<int> coords;
   vector<int> neigh;
};


//============================================================================
//                    ------  main program ------
//============================================================================

int main(int argc, char *argv[])
{
   //..................................................
   //               input variables
   //..................................................
   int Lside,Ntpts;
   int nsnap,thetaM,phiM;
   int put_noise,tstp_ini,renorm_M;
   double dt,speed,damping;
   double Jex_ini, Jex_fin;
   double Jey_ini, Jey_fin;
   double Jez_ini, Jez_fin;
   double Temp_ini,Temp_fin;
   double J_o, Temp_o;
   double tstp_0=0.0;
   char* flout;
   //..................................................
   //             internal variables
   //..................................................
   int nx,ny,nsize;
   double Tmax,Jex,Jey,Jez,Thermostat;
   double Temp_Numerator;
   std::vector<double> spins;
   std::vector<double> vrtx;
   std::vector<double> Tlat;
   double Boltzmann=8.617333262e-5;
   //..................................................
   //             output variables
   //..................................................
   std::vector<double> Svec(3);
   std::vector<double> Svec2(3);
   std::vector<double> Polar(3);
   double Etot, Temp;
   //
   //..................................................
   //                timer
   //..................................................
   std::chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;
   print_line_star(60);
   cout << "     Test progam: 2D spin system (MF)" << endl;
   print_line_star(60);
   start_tot = std::chrono::system_clock::now();

   try
   {
      //..................................................
      //           READ GENERAL INPUT
      //..................................................
      {
         if(argc<3) throw("COMMAND LINE ARGUMENT MISSING");
         if (argc < 4)
         {
            std::cerr << " Please provide a prefix for the output files. Exiting ..." << std::endl;
            return 1;
         }
         flout = argv[3];
         //
         find_param(argv[1], "__Lside="      , Lside      );
         find_param(argv[1], "__Ntpts="      , Ntpts      );
         find_param(argv[1], "__nsnap="      , nsnap      );
         find_param(argv[1], "__dt="         , dt         );
         find_param(argv[1], "__Jex_ini="    , Jex_ini    );
         find_param(argv[1], "__Jex_fin="    , Jex_fin    );
         find_param(argv[1], "__Jey_ini="    , Jey_ini    );
         find_param(argv[1], "__Jey_fin="    , Jey_fin    );
         find_param(argv[1], "__Jez_ini="    , Jez_ini    );
         find_param(argv[1], "__Jez_fin="    , Jez_fin    );
         find_param(argv[1], "__J_o="        , J_o        );
         find_param(argv[1], "__speed="      , speed      );
         find_param(argv[1], "__damping="    , damping    );
         find_param(argv[1], "__thetaM="     , thetaM     );
         find_param(argv[1], "__phiM="       , phiM       );
         find_param(argv[1], "__Temp_ini="   , Temp_ini   );
         find_param(argv[1], "__Temp_fin="   , Temp_fin   );
         find_param(argv[1], "__Temp_o="     , Temp_o     );
         find_param(argv[1], "__put_noise="  , put_noise  );
         find_param(argv[1], "__tstp_ini="   , tstp_ini   );
         find_param(argv[1], "__renorm_M="   , renorm_M   );
      }
      Tmax = dt*Ntpts;
      nx = Lside;
      ny = Lside;
      nsize = 3*nx*ny;

      //..................................................
      //           SETUP THE RANDOM GENERATOR SEED
      //..................................................
      unsigned long seed = std::chrono::duration_cast<std::chrono::milliseconds> (std::chrono::system_clock::now().time_since_epoch()).count()-1566563000000;
      generator.seed(abs(seed));

      //..................................................
      cout << "nx = " << nx << " ny = " << ny << endl;
      cout << "Initial Jex = " << Jex_ini << " Initial Jez = " << Jez_ini << endl;
      cout << "Final Jex = " << Jex_fin << " Final Jez = " << Jez_fin << endl;
      cout << "The speed of the quench is: " << speed << endl;
      cout << "The seed is: " << seed << endl;

      //..................................................
      spins.resize(nsize);
      vrtx.resize(nx*ny);
      Tlat.resize(nx*ny);

      //..................................................
      //           Initial state
      //..................................................
      FILE *file = fopen(argv[2], "r");
      int i = 0;
      if (file != NULL)
      {
         while (!feof(file) && i < 3*nx*ny)
         {
            fscanf(file, "%lf", &spins[i]);
            i++;
         }
         fclose(file);
      }
      else
      {
         printf("Unable to open the spin initial state file.");
         return 0;
      }
      //..................................................
      //           Initial state Randomness
      //..................................................
      if(tstp_ini==0)
      {
         if(put_noise==1)
         {
            //for Z ferro or antif
            if (spins[2]==0.5 || spins[2]==-0.5)
            {
               for (int l=0; l<nx*ny; l++)
               {
                  int theta=rand() % (thetaM+1);
                  int phi=rand() % (phiM+1);
                  double tempX=spins[3*l];
                  double tempY=spins[3*l+1];
                  double tempZ=spins[3*l+2];
                  //
                  spins[3*l]=( cos(phi*PI/180)*tempX - cos(theta*PI/180)*sin(phi*PI/180)*tempY +sin(theta*PI/180)*sin(phi*PI/180)*tempZ);
                  spins[3*l+1]=( sin(phi*PI/180)*tempX + cos(theta*PI/180)*cos(phi*PI/180)*tempY +sin(theta*PI/180)*cos(phi*PI/180)*tempZ );
                  spins[3*l+2]=( sin(theta*PI/180)*tempY + cos(theta*PI/180)*tempZ );
               }
            }
            //for X ferro or antif
            else
            {
               for (int l=0; l<nx*ny; l++)
               {
                  int theta=rand() % (thetaM+1);
                  int phi=rand() % (phiM+1);
                  double tempX=spins[3*l];
                  double tempY=spins[3*l+1];
                  double tempZ=spins[3*l+2];
                  //
                  spins[3*l]=( cos(phi*PI/180)*tempY - cos(theta*PI/180)*sin(phi*PI/180)*tempZ +sin(theta*PI/180)*sin(phi*PI/180)*tempX);
                  spins[3*l+1]=( sin(phi*PI/180)*tempY + cos(theta*PI/180)*cos(phi*PI/180)*tempZ +sin(theta*PI/180)*cos(phi*PI/180)*tempX );
                  spins[3*l+2]=( sin(theta*PI/180)*tempZ + cos(theta*PI/180)*tempX );
               }
            }
         }
      }
      else
      {
         // opening old file
         char flspins_ini[255];
         strcpy(flspins_ini,flout);
         strcat(flspins_ini,"_spins.dat");
         //
         fstream fileoldconf(flspins_ini);
         ofstream filelastconf("starting_configuration.dat");
         GotoLine(fileoldconf, int(tstp_ini/nsnap));
         fileoldconf >> tstp_0;
         filelastconf << tstp_0;
         cout << "Restarting old calculation at tstp = " << tstp_0 << endl;
         for(int k=0; k<nsize; k++)
         {
            fileoldconf >> spins[k];
            filelastconf << " " << spins[k] << " ";
         }
         filelastconf << "\n";
         fileoldconf.close();
         filelastconf.close();
         exit(1);
      }


      //..................................................
      //           time propagation
      //..................................................
      int tstp;
      //..................................................
      spin2D_t spinsys(nx, ny, damping);

      spinsys.gen_latt();
      spinsys.gen_neigh();

      //..................................................
      char flspins[255];
      strcpy(flspins,flout);
      strcat(flspins,"_spins.dat");

      char flbeff[255];
      strcpy(flbeff,flout);
      strcat(flbeff,"_beff.dat");

      char flexpec[255];
      strcpy(flexpec,flout);
      strcat(flexpec,"_expec.dat");

      char flj[255];
      strcpy(flj,flout);
      strcat(flj,"_J.dat");

      char flvrtx[255]              ; char flvrty[255]              ; char flvrtz[255];
      strcpy(flvrtx,flout)          ; strcpy(flvrty,flout)          ; strcpy(flvrtz,flout);
      strcat(flvrtx,"_vortex_x.dat"); strcat(flvrty,"_vortex_y.dat"); strcat(flvrtz,"_vortex_z.dat");

      char fltlat[255];
      strcpy(fltlat,flout);
      strcat(fltlat,"_Tlat.dat");


      //..................................................TIME LOOP
      for(tstp=0; tstp < Ntpts; tstp++)
      {
         Thermostat = gen_Temp(tstp*dt, Temp_ini, Temp_fin, Temp_o, speed);
         Jex  = gen_Jex(tstp*dt, Jex_ini, Jex_fin ,J_o ,speed);
         Jey  = gen_Jey(tstp*dt, Jey_ini, Jey_fin ,J_o ,speed);
         Jez  = gen_Jez(tstp*dt, Jez_ini, Jez_fin ,J_o ,speed);
         //
         //
         spinsys.timestep(tstp, dt, spins, {Jex,Jey,Jez} , Thermostat, renorm_M);
         //
         //...............................................
         if(tstp % nsnap == 0)
         {
            FILE *fspins=fopen(flspins,"a");
            FILE *fj=fopen(flj,"a");
            /*
               FILE *fbeff=fopen(flbeff,"a");
               FILE *fvrtx=fopen(flvrtx,"a");
               FILE *fvrty=fopen(flvrty,"a");
               FILE *fvrtz=fopen(flvrtz,"a");
            */
            FILE *ftlat=fopen(fltlat,"a");
            FILE *fexpec=fopen(flexpec,"a");

            //................................................LATTICE
            cout << "print latt" <<endl;
            fprintf(fspins,"%.10g",dt*(tstp+1));
            for(int k=0; k<nsize; k++)fprintf(fspins," %.10g",spins[k]);
            fprintf(fspins,"\n");

            //................................................PARAMETERS
            fprintf(fj,"%.10g\t%lf\t%lf\t%lf\t%lf\n",dt*(tstp+1),Jex,Jey,Jez,Thermostat);

            //................................................OBSERVABLES
            cout << "print obs" <<endl;
            Etot=0.0;
            Temp=0.0;
            Temp_Numerator=0.0;
            Svec ={0.0,0.0,0.0};
            Svec2={0.0,0.0,0.0};
            Polar={0.0,0.0,0.0};
            //
            /*
               fprintf(fbeff,"%.10g",dt*(tstp+1));
               fprintf(fvrtx,"%.10g",dt*(tstp+1));
               fprintf(fvrty,"%.10g",dt*(tstp+1));
               fprintf(fvrtz,"%.10g",dt*(tstp+1));
            */
            fprintf(ftlat,"%.10g",dt*(tstp+1));
            //
            for (int k=0;k<nx*ny;k++)
            {
               // I'm defining the local spin and Beff
               std::vector<double> spin_i = spinsys.get_spin(k,spins);
               std::vector<double> beff_i = spinsys.gen_beff(k,spins,{Jex,Jey,Jez});
               // Beff on the lattice
               /*
                  fprintf(fbeff," %.10g %.10g %.10g",beff_i[0],beff_i[1],beff_i[2]);
               */
               // Total energy
               double E_i = dot(spin_i,beff_i);
               Etot += E_i/(nx*ny);
               // Vortex on the lattice
               /*
                  std::vector<double> vrtx = spinsys.gen_vrtx(k,spins);
                  fprintf(fvrtx," %.10g",vrtx[0]);
                  fprintf(fvrty," %.10g",vrtx[1]);
                  fprintf(fvrtz," %.10g",vrtx[2]);
               */
               // Spin expectations
               for (int j=0;j<3;j++) Svec[j] += spin_i[j]/(nx*ny);
               // Spin**2 expectations
               for (int j=0;j<3;j++) Svec2[j] += pow(spin_i[j],2)/(nx*ny);
               // Spin expectations in polar coordinates
               for (int j=0;j<3;j++) Polar[j] += compute_polar(k,j,spins)/(nx*ny);
               // Effective Temperature
               std::vector<double> BxS = cross(beff_i, spin_i);
               double BxSsq_i = dot(BxS,BxS);
               Temp_Numerator += BxSsq_i/(nx*ny);
               // Temperature on the lattice
               Tlat[k] = BxSsq_i/(-2.0*E_i);
               fprintf(ftlat," %.10g",Tlat[k]);
            }
            Temp = Temp_Numerator / (-2.0*Etot);
            //
            fprintf(fexpec,"%.10g",dt*(tstp+1));
            fprintf(fexpec,"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf"
                          ,Svec[0] ,Svec[1] ,Svec[2]
                          ,Svec2[0],Svec2[1],Svec2[2]
                          ,Etot,Temp/Boltzmann
                          ,Polar[0],Polar[1],Polar[2]);
            //
            /*
               fprintf(fbeff,"\n");
               fprintf(fvrtx,"\n");
               fprintf(fvrty,"\n");
               fprintf(fvrtz,"\n");
            */
            fprintf(ftlat,"\n");
            fprintf(fexpec,"\n");
            //
            fclose(fj);
            fclose(fspins);
            /*
               fclose(fbeff);
               fclose(fvrtx);
               fclose(fvrty);
               fclose(fvrtz);
            */
            fclose(ftlat);
            fclose(fexpec);
            //
            cout << Etot << "  " << Temp/Boltzmann << endl;
            //
         }
         cout << (tstp+1) << " /" << Ntpts << "  T: " << Thermostat << "  Jx: " << Jex << "  Jy: " << Jey << "  Jz: " << Jez <<endl;
      }

      //......................................................STOP CHRONO
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
