#include <vector>
#include <cmath>
#include <math.h>
#include <eigen3/Eigen/Geometry>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;


//-------------------------------------------------------------------------


std::vector<double> cross(vector<double> x, vector<double> y)
{
   std::vector<double> z(3);
   z = {x[1]*y[2] - x[2]*y[1],x[2]*y[0] - x[0]*y[2],x[0]*y[1] - x[1]*y[0]};
   return z;
}


//-------------------------------------------------------------------------


double dot(vector<double> x, vector<double> y)
{
   double z;
   z = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
   return z;
}


//-------------------------------------------------------------------------


std::vector<double> euler_rot(double phi, double theta, double psi, std::vector<double> beff, int diff)
{
   Eigen::Vector3d vec_in; vec_in.setZero(3);
   for(int j=0;j<3;j++) vec_in(j)=beff[j];
   //
   Eigen::AngleAxis< double > AngleAxis_X(phi  ,Eigen::Vector3d::UnitX());
   Eigen::AngleAxis< double > AngleAxis_Y(theta,Eigen::Vector3d::UnitY());
   Eigen::AngleAxis< double > AngleAxis_Z(psi  ,Eigen::Vector3d::UnitZ());
   //
   Eigen::Matrix3d rotX = AngleAxis_X.toRotationMatrix();
   Eigen::Matrix3d rotY = AngleAxis_Y.toRotationMatrix();
   Eigen::Matrix3d rotZ = AngleAxis_Z.toRotationMatrix();
   Eigen::Matrix3d rotation = rotZ*rotY*rotX;
   if(rotation.isUnitary()!=1) cout << " WARNING ROTATION NON UNITARY " << endl;
   //
   Eigen::Vector3d vec_out; vec_out.setZero(3);
   vec_out = rotation * vec_in;
   //
   std::vector<double> beff_rot(3);
   if(diff==1)
   {
      for(int j=0;j<3;j++) beff_rot[j] = vec_out(j) - vec_in(j);
   }
   else
   {
      for(int j=0;j<3;j++) beff_rot[j] = vec_out(j);
   }
   return beff_rot;
}


//-------------------------------------------------------------------------


std::vector<double> polar_fluct(double dr, double dphi, double dtheta, std::vector<double> vec_in, int diff)
{
   std::vector<double> vec_out(3);
   double rho = sqrt(pow(vec_in[0],2) + pow(vec_in[1],2) + pow(vec_in[2],2));
   double theta = acos(vec_in[2]/rho);
   double phi = atan2 (vec_in[1],vec_in[0]);
   //
   vec_out[0] = (rho+dr)*sin(theta+dtheta)*cos(phi+dphi);
   vec_out[1] = (rho+dr)*sin(theta+dtheta)*sin(phi+dphi);
   vec_out[2] = (rho+dr)*cos(theta+dtheta);
   //
   std::vector<double> beff_rot(3);
   if(diff==1)
   {
      for(int j=0;j<3;j++) beff_rot[j] = vec_out[j] - vec_in[j];
   }
   else
   {
      for(int j=0;j<3;j++) beff_rot[j] = vec_out[j];
   }
   return beff_rot;
}


//-------------------------------------------------------------------------


std::vector<double> renorm_spin(std::vector<double> vec_in, double lenght)
{
   std::vector<double> vec_out(3);
   double rho = sqrt(pow(vec_in[0],2) + pow(vec_in[1],2) + pow(vec_in[2],2));
   double theta = acos(vec_in[2]/rho);
   double phi = atan2 (vec_in[1],vec_in[0]);
   //
   vec_out[0] = lenght*sin(theta)*cos(phi);
   vec_out[1] = lenght*sin(theta)*sin(phi);
   vec_out[2] = lenght*cos(theta);
   //
   return vec_out;
}


//-------------------------------------------------------------------------


template <typename T>
int sgn(T val)
{
   return (T(0) < val) - (val < T(0));
}


//-------------------------------------------------------------------------


template<class odesys>
void step_rungekutta4(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec)
{
   //VERSION 1
   /*
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n),k3(n),k4(n);
   //
   // K1 = h * f( tn , Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<n;i++) k1[i] = h*dy[i];
   //
   // K2 = h * f( tn+0.5*h , Sn + 0.5*K1)
   for(int i=0;i<n;i++)y1[i] = y0[i] + 0.5*k1[i];
   dy = ode.deriv(t+0.5*h,y1,Jvec);
   for(int i=0;i<n;i++) k2[i] = h*dy[i];
   //
   //  K3 = h * ( tn+0.5*h , Sn + 0.5*K2)
   for(int i=0;i<n;i++) y1[i] = y0[i] + 0.5*k2[i];
   dy = ode.deriv(t+0.5*h,y1,Jvec);
   for(int i=0;i<n;i++) k3[i] = h*dy[i];
   //
   //  K4 = h * ( tn+h , Sn + K3)
   for(int i=0;i<n;i++) y1[i] = y0[i] + k3[i];
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<n;i++) k4[i] = h*dy[i];
   //
   //
   for(int i=0;i<n;i++) y0[i] = y0[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
   */
   //
   //VERSION 2
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n),k3(n),k4(n);
   //
   // K1 = h * f( tn , Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<n;i++)
   {
      k1[i] = h*dy[i];
      y1[i] = y0[i] + 0.5*k1[i];
   }
   dy = ode.deriv(t+0.5*h,y1,Jvec);
   for(int i=0;i<n;i++)
   {
      k2[i] = h*dy[i];
      y1[i] = y0[i] + 0.5*k2[i];
   }
   dy = ode.deriv(t+0.5*h,y1,Jvec);
   for(int i=0;i<n;i++)
   {
      k3[i] = h*dy[i];
      y1[i] = y0[i] + k3[i];
   }
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<n;i++)
   {
      k4[i] = h*dy[i];
      y0[i] = y0[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
   }
}


//-------------------------------------------------------------------------


template<class odesys>
void step_Ito(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat, int renorm=0)
{
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   //
   std::vector<double> Wk, Sk;
   for(int i=0;i<dim;i++)
   {
      Wk.push_back(ode.gauss_distribution(1.0));
      //Sk.push_back(ode.flat_distribution(1.0));
      Sk.push_back(0.0);
   }

   //
   // K1 = h * f( tn , Sn ) + ( Wk - Sk*sqrt(h)) * -bk( Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<dim;i++)
   {
      std::vector<double> hrand(3);
      for(int j=0;j<3;j++) hrand[j]=ode.gauss_distribution(var);
      std::vector<double> spin_i = ode.get_spin(i,y0);
      std::vector<double> bk = cross(spin_i,hrand);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] - sgn<double>(Sk[i])) * bk[j];
   }
   //
   for(int i=0;i<n;i++)y1[i] = y0[i] + k1[i];
   //
   // K2 = h * f( tn+h , Sn+K1 ) + ( Wk + Sk*sqrt(h)) * -bk( Sn+K1 )
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<dim;i++)
   {
      std::vector<double> hrand(3);
      for(int j=0;j<3;j++) hrand[j]=ode.gauss_distribution(var);
      std::vector<double> spin_i = ode.get_spin(i,y1);
      std::vector<double> bk = cross(spin_i,hrand);
      //
      for(int j=0;j<3;j++) k2[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] + sgn<double>(Sk[i])) * bk[j];
   }
   //
   for(int i=0;i<dim;i++)
   {
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + (k1[3*i+j] + k2[3*i+j])/2.0;
   }
   //
   if(renorm==1)
   {
      for(int i=0;i<dim;i++)
      {
         std::vector<double> Swrong = ode.get_spin(i,y0);
         std::vector<double> Sright = renorm_spin(Swrong, 0.5);
         for(int j=0;j<3;j++) y0[3*i+j] = Sright[j];
      }
   }
}


//-------------------------------------------------------------------------


template<class odesys>
void step_Ito_alphaInvariant(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat, int renorm=0)
{
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n);
   //
   int dim = ode.get_dim();
   double var = 2.*ode.get_damping()*Thermostat;
   double devstd = sqrt(var);
   //
   std::vector<double> Wk, Sk;
   for(int i=0;i<dim;i++)
   {
      Wk.push_back(ode.gauss_distribution(1.0));
      Sk.push_back(sgn<double>(ode.flat_distribution(1.0)));
   }
   //
   // K1 = h * f( tn , Sn ) + ( Wk - Sk*sqrt(h)) * -bk( Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<dim;i++)
   {
      std::vector<double> hrand(3);
      for(int j=0;j<3;j++) hrand[j]=ode.gauss_distribution(devstd);
      std::vector<double> spin_i = ode.get_spin(i,y0);
      std::vector<double> bk = cross(spin_i,hrand);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*(dy[3*i+j]-y0[3*i+j]*var) + sqrt(h)*( Wk[i] - Sk[i] ) * bk[j];//-y0[3*i+j]*var
   }
   //
   for(int i=0;i<n;i++)y1[i] = y0[i] + k1[i];
   //
   // K2 = h * f( tn+h , Sn+K1 ) + ( Wk + Sk*sqrt(h)) * -bk( Sn+K1 )
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<dim;i++)
   {
      std::vector<double> hrand(3);
      for(int j=0;j<3;j++) hrand[j]=ode.gauss_distribution(devstd);
      std::vector<double> spin_i = ode.get_spin(i,y1);
      std::vector<double> bk = cross(spin_i,hrand);
      //
      for(int j=0;j<3;j++) k2[3*i+j] = h*(dy[3*i+j]-y1[3*i+j]*var) + sqrt(h)*( Wk[i] + Sk[i] ) * bk[j];//-y1[3*i+j]*var
   }
   //
   for(int i=0;i<dim;i++)
   {
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + (k1[3*i+j] + k2[3*i+j])/2.0;
   }
   //
   if(renorm==1)
   {
      for(int i=0;i<dim;i++)
      {
         std::vector<double> Swrong = ode.get_spin(i,y0);
         std::vector<double> Sright = renorm_spin(Swrong, 0.5);
         for(int j=0;j<3;j++) y0[3*i+j] = Sright[j];
      }
   }
}


//-------------------------------------------------------------------------



struct rparams
{
   double h;
   double gamma;
   double mu;
   std::vector<double> Sn;
   std::vector<double> Bk;
   std::vector<double> Wn;
};

int
rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double h = ((struct rparams *) params)->h;
  double gamma = ((struct rparams *) params)->gamma;
  double mu = ((struct rparams *) params)->mu;
  std::vector<double> Sn = ((struct rparams *) params)->Sn;
  std::vector<double> Bk = ((struct rparams *) params)->Bk;
  std::vector<double> Wn = ((struct rparams *) params)->Wn;
  int dim = Bk.size()/3;
  //
  std::vector<double> DeltaP_i(3);
  std::vector<double> DeltaM_i(3);
  std::vector<double> beff_i(3);
  std::vector<double> Wn_i(3);
  std::vector<double> Fn_i(3);
  //
  for(int i=0; i < dim; i++)
  {
     for(int j=0;j<3;j++)
     {
        DeltaP_i[j] = (gsl_vector_get(x, 3*i+j) + Sn[3*i+j])/2.0;
        DeltaM_i[j] = (gsl_vector_get(x, 3*i+j) - Sn[3*i+j]);
        beff_i[j] = Bk[3*i+j];
        Wn_i[j]   = Wn[3*i+j];
     }
     //
     const std::vector<double> Fn_1 = cross(DeltaP_i, beff_i);
     const std::vector<double> Fn_2 = cross(Fn_1, DeltaP_i);
     const std::vector<double> Fn_3 = cross(DeltaP_i, Wn_i);
     //
     for(int j=0;j<3;j++)
     {
        const double Fn = DeltaM_i[j] + Fn_1[j]*h + gamma*Fn_2[j]*h + mu*Fn_3[j];
        gsl_vector_set (f, 3*i+j, Fn);
     }
     //
     return GSL_SUCCESS;
  }
}

template<class odesys>
void step_midpoint(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> Beff_yavrg(n),yavrg(n);
   //
   int dim = ode.get_dim();
   double gamma = ode.get_damping();
   double mu = sqrt(2.*gamma*Thermostat);
   //
   std::vector<double> Wn;
   for(int i=0;i<n;i++) Wn.push_back(ode.gauss_distribution(h));
   //
   // GSL stuff
   int status, iter = 0;
   size_t i;
   const size_t ngsl = n;
   //
   // declare GSL stuff
   gsl_vector *x = gsl_vector_alloc (ngsl);
   struct rparams p;
   gsl_multiroot_function f;
   //
   // solver init
   const gsl_multiroot_fsolver_type *T;
   gsl_multiroot_fsolver *s;
   //T = gsl_multiroot_fsolver_hybrids; NON SI MUOVE ultra lento
   T = gsl_multiroot_fsolver_hybrid;// NON SI MUOVE
   //T = gsl_multiroot_fsolver_dnewton; SINGULAR MATRIX
   //T = gsl_multiroot_fsolver_broyden; SINGULAR MATRIX
   s = gsl_multiroot_fsolver_alloc (T, ngsl);
   //
   // Bfield init
   for(int i=0; i < dim; i++)
   {
     std::vector<double> beff_i = ode.gen_beff(i,y0,Jvec);
     for(int j=0;j<3;j++) Beff_yavrg[3*i+j] = beff_i[j];
   }
   //
   // functional init
   for(int i=0;i<n;i++)gsl_vector_set (x, i, y0[i]);
   p = {h, gamma, mu, y0, Beff_yavrg, Wn};
   f = {&rosenbrock_f, ngsl, &p};
   gsl_multiroot_fsolver_set (s, &f, x);
   //
   //
   // Loop over root convergence
   do
   {
      //
      //Solve the root finding step
      printf ("iter = %i\n", iter);
      iter++;
      status = gsl_multiroot_fsolver_iterate(s);
      if (status) break;
      status = gsl_multiroot_test_residual(s->f, 1e-7);
      //
      // getting out yavrg = (yk+Sn)/2
      for(int i=0;i<n;i++) yavrg[i] = ( gsl_vector_get(s->x, i) + y0[i] )/2.0;
      //
      // Loop to recalculate Beff(yavrg)
      for(int i=0; i < dim; i++)
      {
         std::vector<double> beff_i = ode.gen_beff(i,yavrg,Jvec);
         for(int j=0;j<3;j++) Beff_yavrg[3*i+j] = beff_i[j];
      }
      printf ("B = %.10g  %.10g  %.10g\n", Beff_yavrg[234],yavrg[234],y0[234],gsl_vector_get(s->x, 234));
      printf ("B = %.10g  %.10g  %.10g\n", Beff_yavrg[235],yavrg[235],y0[235],gsl_vector_get(s->x, 235));
      printf ("B = %.10g  %.10g  %.10g\n", Beff_yavrg[236],yavrg[236],y0[236],gsl_vector_get(s->x, 236));
      printf("\n\n");
      //
      for(int i=0;i<n;i++)gsl_vector_set(x, i, gsl_vector_get(s->x, i) );
      p = {h, gamma, mu, y0, Beff_yavrg, Wn};
      f = {&rosenbrock_f, ngsl, &p};
      gsl_multiroot_fsolver_set (s, &f, x);
      //
   }
   while (status == GSL_CONTINUE && iter < 50);
   printf ("status = %s  itermax = %i\n", gsl_strerror (status), iter);
   //
   for(int i=0;i<n;i++) y0[n] = gsl_vector_get (s->x, i);
   //
   gsl_multiroot_fsolver_free (s);
   gsl_vector_free (x);
   //return 0;
 }












/*


   do
   {
     iter++;
     for(int i=0;i<dim;i++)
     {

     }




     status = gsl_multiroot_fsolver_iterate(s);
     if (status) break;
     status = gsl_multiroot_test_residual(s->f, 1e-7);
   }





template<class odesys>
void step_rungekutta2_v2(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   double dr,dphi,dtheta;
   //
   std::vector<double> Wk, Sk;
   double sumWk=0.0,sumSk=0.0;
   for(int i=0;i<dim;i++)
   {
      double w = ode.gauss_distribution(1.0); Wk.push_back(w); sumWk+=w/dim;
      double s = ode.flat_distribution(1.0) ; Sk.push_back(s); sumSk+=s/dim;
   }

   //
   // K1 = h * f( tn , Sn ) + ( Wk - Sk*sqrt(h)) * -bk( Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<dim;i++)
   {
      dr     = ode.gauss_distribution(var);
      dphi   = ode.gauss_distribution(var);
      dtheta = ode.gauss_distribution(var);
      //
      std::vector<double> drift_i = ode.get_spin(i,dy);
      std::vector<double> bk_i = polar_fluct(dr,dphi,dtheta,drift_i,1);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] - sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<n;i++)y1[i] = y0[i] + k1[i];
   //
   // K2 = h * f( tn+h , Sn+K1 ) + ( Wk + Sk*sqrt(h)) * -bk( Sn+K1 )
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<dim;i++)
   {
      dr     = ode.gauss_distribution(var);
      dphi   = ode.gauss_distribution(var);
      dtheta = ode.gauss_distribution(var);
      //
      std::vector<double> drift_i = ode.get_spin(i,dy);
      std::vector<double> bk_i = polar_fluct(dr,dphi,dtheta,drift_i,1);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] + sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<dim;i++)
   {
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + (k1[3*i+j] + k2[3*i+j])/2.0;
   }
   //
   FILE *fdist = fopen("dist.dat", "a");
   fprintf(fdist,"%g  %g\n",sumWk,sumSk);
   fclose(fdist);
}



template<class odesys>
void step_rungekutta2_v3(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   double dr,dphi,dtheta;
   //
   std::vector<double> Wk, Sk;
   double sumWk=0.0,sumSk=0.0;
   for(int i=0;i<dim;i++)
   {
      double w = ode.gauss_distribution(1.0); Wk.push_back(w); sumWk+=w/dim;
      double s = ode.flat_distribution(1.0) ; Sk.push_back(s); sumSk+=s/dim;
   }

   //
   // K1 = h * f( tn , Sn ) + ( Wk - Sk*sqrt(h)) * -bk( Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<dim;i++)
   {
      dr     = ode.gauss_distribution(var);
      dphi   = ode.gauss_distribution(var);
      dtheta = ode.gauss_distribution(var);
      //
      std::vector<double> spin_i = ode.get_spin(i,y0);
      std::vector<double> beff_i = ode.gen_beff(i,y0,Jvec);
      std::vector<double> drift_i = cross(beff_i, spin_i);
      std::vector<double> bk_i = polar_fluct(dr,dphi,dtheta,drift_i,1);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] - sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<n;i++)y1[i] = y0[i] + k1[i];
   //
   // K2 = h * f( tn+h , Sn+K1 ) + ( Wk + Sk*sqrt(h)) * -bk( Sn+K1 )
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<dim;i++)
   {
      dr     = ode.gauss_distribution(var);
      dphi   = ode.gauss_distribution(var);
      dtheta = ode.gauss_distribution(var);
      //
      std::vector<double> spin_i = ode.get_spin(i,y1);
      std::vector<double> beff_i = ode.gen_beff(i,y1,Jvec);
      std::vector<double> drift_i = cross(beff_i, spin_i);
      std::vector<double> bk_i = polar_fluct(dr,dphi,dtheta,drift_i,1);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] + sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<dim;i++)
   {
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + (k1[3*i+j] + k2[3*i+j])/2.0;
   }
   //
   FILE *fdist = fopen("dist.dat", "a");
   fprintf(fdist,"%g  %g\n",sumWk,sumSk);
   fclose(fdist);
}


template<class odesys>
void step_rungekutta2_v4(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   double dr,dphi,dtheta;
   //
   std::vector<double> Wk, Sk;
   double sumWk=0.0,sumSk=0.0;
   for(int i=0;i<dim;i++)
   {
      double w = ode.gauss_distribution(1.0); Wk.push_back(w); sumWk+=w/dim;
      double s = ode.flat_distribution(1.0) ; Sk.push_back(s); sumSk+=s/dim;
   }

   //
   // K1 = h * f( tn , Sn ) + ( Wk - Sk*sqrt(h)) * -bk( Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<dim;i++)
   {
      dr     = ode.gauss_distribution(var);
      dphi   = ode.gauss_distribution(var);
      dtheta = ode.gauss_distribution(var);
      //
      std::vector<double> spin_i = ode.get_spin(i,y0);
      std::vector<double> beff_i = ode.gen_beff(i,y0,Jvec);
      std::vector<double> hk_i = polar_fluct(dr,dphi,dtheta,beff_i,1);
      std::vector<double> bk_i = cross(hk_i, spin_i);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] - sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<n;i++)y1[i] = y0[i] + k1[i];
   //
   // K2 = h * f( tn+h , Sn+K1 ) + ( Wk + Sk*sqrt(h)) * -bk( Sn+K1 )
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<dim;i++)
   {
      dr     = ode.gauss_distribution(var);
      dphi   = ode.gauss_distribution(var);
      dtheta = ode.gauss_distribution(var);
      //
      std::vector<double> spin_i = ode.get_spin(i,y1);
      std::vector<double> beff_i = ode.gen_beff(i,y1,Jvec);
      std::vector<double> hk_i = polar_fluct(dr,dphi,dtheta,beff_i,1);
      std::vector<double> bk_i = cross(hk_i, spin_i);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] + sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<dim;i++)
   {
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + (k1[3*i+j] + k2[3*i+j])/2.0;
   }
   //
   FILE *fdist = fopen("dist.dat", "a");
   fprintf(fdist,"%g  %g\n",sumWk,sumSk);
   fclose(fdist);
}


template<class odesys>
void step_rungekutta2_v4(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   double alpha,beta,gamma;
   //
   std::vector<double> Wk, Sk;
   double sumWk=0.0,sumSk=0.0;
   for(int i=0;i<dim;i++)
   {
      double w = ode.gauss_distribution(1.0); Wk.push_back(w); sumWk+=w/dim;
      double s = ode.flat_distribution(1.0) ; Sk.push_back(s); sumSk+=s/dim;
   }

   //
   // K1 = h * f( tn , Sn ) + ( Wk - Sk*sqrt(h)) * -bk( Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<dim;i++)
   {
      std::vector<double> beff_i = ode.gen_beff(i,y0,Jvec);
      std::vector<double> spin_i = ode.get_spin(i,y0);
      //
      alpha = ode.gauss_distribution(var);
      beta  = ode.gauss_distribution(var);
      gamma = ode.gauss_distribution(var);
      //
      std::vector<double> hrand_i = euler_rot(alpha,beta,gamma,beff_i,1);
      std::vector<double> bk_i = cross(hrand_i,spin_i);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] - sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<n;i++)y1[i] = y0[i] + k1[i];
   //
   // K2 = h * f( tn+h , Sn+K1 ) + ( Wk + Sk*sqrt(h)) * -bk( Sn+K1 )
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<dim;i++)
   {
      std::vector<double> beff_i = ode.gen_beff(i,y1,Jvec);
      std::vector<double> spin_i = ode.get_spin(i,y1);
      //
      alpha = ode.gauss_distribution(var);
      beta  = ode.gauss_distribution(var);
      gamma = ode.gauss_distribution(var);
      //
      std::vector<double> hrand_i = euler_rot(alpha,beta,gamma,beff_i,1);
      //
      std::vector<double> bk_i = cross(hrand_i,spin_i);
      //
      for(int j=0;j<3;j++) k2[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] + sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<dim;i++)
   {
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + (k1[3*i+j] + k2[3*i+j])/2.0;
   }
   //
   FILE *fdist = fopen("dist.dat", "a");
   fprintf(fdist,"%g  %g\n",sumWk,sumSk);
   fclose(fdist);
}



template<class odesys>
void step_Euler(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> dy(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   double alpha,beta,gamma;
   //
   dy = ode.deriv(t,y0,Jvec);
   //
   for(int i=0;i<dim;i++)
   {
      std::vector<double> beff_i = ode.gen_beff(i,y0,Jvec);
      std::vector<double> spin_i = ode.get_spin(i,y0);
      //
      alpha = ode.gauss_distribution(var);
      beta  = ode.gauss_distribution(var);
      gamma = ode.gauss_distribution(var);
      //
      std::vector<double> hrand_i = euler_rot(alpha,beta,gamma,beff_i,1);
      //
      std::vector<double> bk_i = cross(hrand_i,spin_i);
      //
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + h*dy[3*i+j] + bk_i[j];
   }
}



template<class odesys>
void step_rungekutta2_v3(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> dy(n),y1(n),k1(n),k2(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   double alpha,beta,gamma;
   //
   std::vector<double> Wk, Sk;
   double sumWk=0.0,sumSk=0.0;
   for(int i=0;i<dim;i++)
   {
      double w = ode.gauss_distribution(1.0); Wk.push_back(w); sumWk+=w/dim;
      double s = ode.flat_distribution(1.0) ; Sk.push_back(s); sumSk+=s/dim;
   }

   //
   // K1 = h * f( tn , Sn ) + ( Wk - Sk*sqrt(h)) * -bk( Sn )
   dy = ode.deriv(t,y0,Jvec);
   for(int i=0;i<dim;i++)
   {
      alpha = ode.gauss_distribution(var);
      beta  = ode.gauss_distribution(var);
      gamma = ode.gauss_distribution(var);
      //
      std::vector<double> drift_i = ode.get_spin(i,dy);
      std::vector<double> bk_i = euler_rot(alpha,beta,gamma,drift_i,1);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] - sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<n;i++)y1[i] = y0[i] + k1[i];
   //
   // K2 = h * f( tn+h , Sn+K1 ) + ( Wk + Sk*sqrt(h)) * -bk( Sn+K1 )
   dy = ode.deriv(t+h,y1,Jvec);
   for(int i=0;i<dim;i++)
   {
      alpha = ode.gauss_distribution(var);
      beta  = ode.gauss_distribution(var);
      gamma = ode.gauss_distribution(var);
      //
      std::vector<double> drift_i = ode.get_spin(i,dy);
      std::vector<double> bk_i = euler_rot(alpha,beta,gamma,drift_i,1);
      //
      for(int j=0;j<3;j++) k1[3*i+j] = h*dy[3*i+j] + sqrt(h)*( Wk[i] - sgn<double>(Sk[i])) * bk_i[j];
   }
   //
   for(int i=0;i<dim;i++)
   {
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + (k1[3*i+j] + k2[3*i+j])/2.0;
   }
   //
   FILE *fdist = fopen("dist.dat", "a");
   fprintf(fdist,"%g  %g\n",sumWk,sumSk);
   fclose(fdist);
}



template<class odesys>
void step_Euler2(double t, double h, std::vector<double> &y0, odesys &ode, vector<double> Jvec, double Thermostat)
{
   int n = y0.size();
   std::vector<double> dy(n);
   //
   int dim = ode.get_dim();
   double var = sqrt(2.*ode.get_damping()*Thermostat);
   double alpha,beta,gamma;
   //
   dy = ode.deriv(t,y0,Jvec);
   //
   for(int i=0;i<dim;i++)
   {
      alpha = ode.gauss_distribution(var);
      beta  = ode.gauss_distribution(var);
      gamma = ode.gauss_distribution(var);
      //
      std::vector<double> drift_i = ode.get_spin(i,dy);
      std::vector<double> bk_i = euler_rot(alpha,beta,gamma,drift_i,1);
      //
      for(int j=0;j<3;j++) y0[3*i+j] = y0[3*i+j] + h*dy[3*i+j] + bk_i[j];
   }
}

*/
