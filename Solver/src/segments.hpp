#ifndef ___SEGMENTS___
#define ___SEGMENTS___

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
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>


//============================================================================//


std::random_device random_seed;


std::default_random_engine generator;
double rndm_dre()
{
   std::uniform_real_distribution<double> distribution_flat(0.0,1.0);
   double r = distribution_flat(generator);
   return r;
}


double rndm() {return rand()/(double(RAND_MAX)+1);}


//------------------------------------------------------------------------------


bool tau_uniform_D = true ;
bool tau_uniform_K = true ;
bool tau_uniform_Kp = true ;


//------------------------------------------------------------------------------


class times
{
 public:
  times() {t_start_=0; t_end_=0;};
  times(double t_start, double t_end) {t_start_=t_start; t_end_=t_end;};
  double t_start() const {return t_start_;} // begin of segment
  double t_end() const {return t_end_;}     // end of segment
  void set_t_start(double t_start) {t_start_ = t_start;}
  void set_t_end(double t_end) {t_end_ = t_end;}

 private:
  double t_start_, t_end_;
};

inline bool operator<(const times& t1, const times& t2) {
  return t1.t_start() < t2.t_start();
}
inline bool operator<(const times& t1, const double t2) {
  return t1.t_start() < t2;
}

inline bool operator>(times t1, times t2) {
  return t1.t_start() > t2.t_start();
}

inline bool operator==(times t1, times t2) {
  return t1.t_start() == t2.t_start();
}


//------------------------------------------------------------------------------


typedef std::string path;
typedef std::set<times> segment_container_t;
typedef Eigen::MatrixXd Mat;
typedef std::vector<double> Vec;
typedef std::vector<std::vector<double>> VecVec;
typedef std::vector<std::vector<std::vector<double>>> VecVecVec;
typedef std::vector<Eigen::MatrixXd> VecMat;
typedef std::chrono::time_point<std::chrono::system_clock> duration;


//------------------------------------------------------------------------------


inline int cycle(int i, int size) {return (i>0 ? i-1 : size-1);}


//------------------------------------------------------------------------------


template <class G> inline double interpolate_F(double tau, double Beta, G &F)
{
   double sign=1;
   if (tau<0)
   {
      tau += Beta;
      sign = -1;
   }

   //
   int Ntau_m1 = F.size()-1;
   double dt = Beta / Ntau_m1;

   // linear mesh by default
   double n = tau / dt;

   // dense mesh optional
   if(!tau_uniform_D)
   {
      // generic point on the exponential function
      double b = (Beta/2) / ( exp(dt*Ntau_m1/2) - 1 );
      if( tau<=Beta/2 )
      {
         n = log(tau/b+1) / dt;
      }
      else
      {
         n = Ntau_m1 - log((Beta-tau)/b+1) / dt;
      }

      //// boundary even indexes containing an odd one
      //int ndx_1 = (int)n;
      //if(ndx_1%2==1) ndx_1--;
      //int ndx_2 = ndx_1+2;

      //// points on the exponential function corresponding to the indexes
      //double t1 = b * ( exp(dt*ndx_1) - 1 );
      //double t2 = b * ( exp(dt*ndx_2) - 1 );

      //// new n interpolated on the line connecting the two points on the exponential function
      //n = (tau-t1)/(t2-t1)*(ndx_2-ndx_1) + ndx_1;
   }

   // this is the same between the two schemes
   int n_lower = (int)n;

   // interpolate linearly between n_lower and n_lower+1
   double Finterp = sign * ( F[n_lower] + (n-n_lower)*(F[n_lower+1]-F[n_lower]) );

   //
   return Finterp;
}


//------------------------------------------------------------------------------


template <class S, class G> void construct_matrix( Mat &M, S &segments, double Beta,  G &F)
{
   //
   int N = segments.size();
   M.resize(N,N);
   int row=-1;
   int col=-1;

   //
   for (typename S::iterator it1=segments.begin(); it1!=segments.end(); it1++)
   {
      row++;
      for (typename S::iterator it2=segments.begin(); it2!=segments.end(); it2++)
      {
         col++;

         //
         double argument = it1->t_end()-it2->t_start();
         double sign = 1;

         //
         if(argument<0)
         {
            argument += Beta;
            sign = -1;
         }

         //
         M(row,col) = interpolate_F(argument, Beta, F)*sign;
      }
      //
      col = -1;
   }
}


template <class S, class G> double det_inverse_small( Mat &M, S &segments, double Beta,  G&F)
{
    construct_matrix(M, segments, Beta, F);
    Mat invM = M.inverse();
    swap(M,invM);
    double det = invM.determinant();
    return det;
}


template <class S, class G> double det_inverse( Mat &M, S &segments, double Beta,  G&F)
{
    construct_matrix(M, segments, Beta, F);
    Eigen::FullPivLU<Mat> lu(M);
    if(lu.isInvertible())
    {
       Mat invM = lu.inverse();
       swap(M,invM);
       double det = lu.determinant();
       return det;
    }
    else
    {
       std::cerr<<"Hybridization matrix is not invertible."<<std::endl;
       exit(1);
    }
}


//------------------------------------------------------------------------------


template <class S> inline double segment_overlap(times segment, S &other_segments, int other_full_line, double Beta)
{
   // compute overlap between a segment and a list of segments
   // requires segment with 0<=t_begin<t_end<=Beta
   double length = (segment.t_start()<segment.t_end() ? segment.t_end()-segment.t_start() : segment.t_end()-segment.t_start()+Beta);
   double t_final = segment.t_start()+length;
   double t = segment.t_start();
   double t_final_segment;
   double other_length=0;

   //
   if (other_full_line==1)
   {
      other_length=length;
   }
   else if (other_segments.size()>0)
   {
      typename S::iterator it;
      it = lower_bound(other_segments.begin(), other_segments.end(), t);

      //
      if (it!=other_segments.begin())
      {
         it--;
         t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+Beta);
         if (t<t_final_segment) other_length += (t_final_segment<t_final ? t_final_segment-t : t_final-t);
         it++;
      }

      //
      while(it!=other_segments.end() && it->t_start()<t_final)
      {
         t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+Beta);
         other_length += (t_final_segment<t_final ? t_final_segment-it->t_start() : t_final-it->t_start());
         it++;
      }

      // check if last segment overlaps
      it=other_segments.end();
      it--;
      if (it->t_end()<it->t_start() && t<it->t_end()) other_length += (t_final<it->t_end() ? t_final-t : it->t_end()-t);
   }
   //
   return other_length;
}


inline double compute_length(double r, double l_max, double mu)
{
   if (mu == 0) {return r*l_max;}
   else {return 1/mu*log(r*(exp(mu*l_max)-1)+1);}
}


template <class S> inline double compute_overlap(times segment, S &other_segments, int other_full_line, double Beta)
{
   if (segment.t_start()<segment.t_end())
   {
      return segment_overlap(segment, other_segments, other_full_line, Beta);
   }
   else
   {
      double other_length=0;
      times segment1(0,segment.t_end());
      times segment2(segment.t_start(), Beta);
      other_length += segment_overlap(segment1, other_segments, other_full_line, Beta);
      other_length += segment_overlap(segment2, other_segments, other_full_line, Beta);
      return other_length;
   }
}


template <class S> void compute_intervals(double t, double Beta, double &t_up,
   double &t_down, S &segments, typename S::iterator &s_up, typename S::iterator &s_down)
{
   // compute distances up/down to the next segment and iterators of these segments
   // note: s_down always points to a physical segment, while s_up may point to segments.end()
   if (segments.size() == 0)
   {
      t_up = Beta;
      t_down = Beta;
      s_up = segments.end();
      s_down = segments.end();
   }
   else
   {
      s_up = lower_bound(segments.begin(), segments.end(), t);

      //
      if (s_up == segments.begin())
      {
         s_down = segments.end(); s_down--;
         if (s_down->t_end() < s_down->t_start())
         {
            t_down = t - s_down->t_end();
         }
         else
         {
            t_down = t + Beta - s_down->t_end();
         }
      }
      else
      {
         s_down = s_up; s_down--;
         if (s_down->t_end()>s_down->t_start())
         {
            t_down = t - s_down->t_end();
         }
         else
         {
            t_down = t - (Beta+s_down->t_end());
         }
      }

      //
      if(s_up == segments.end())
      {
         t_up = Beta - t + segments.begin()->t_start();
      }
      else
      {
         t_up = s_up->t_start() - t;
      }
   }
}


//------------------------------------------------------------------------------


void compute_M_up(int k, Mat &M, Vec &Fs, Vec &Fe, double det_rat)
{
   Mat M_new(M.rows()+1,M.rows()+1);
   int i_new, j_new;

   // element (k,k)
   M_new(k,k) = 1./det_rat;

   // row k and column k
   for (int i=0; i<M.rows(); i++)
   {
     i_new = (i<k ? i : i+1);
     M_new(i_new,k) = 0;
     M_new(k,i_new) = 0;

     //
     for (int n=0; n<M.rows(); n++)
     {
        M_new(i_new,k) -= M(i,n)*Fs[n];
        M_new(k,i_new) -= M(n,i)*Fe[n];
     }
     M_new(i_new,k) /= det_rat;
     M_new(k,i_new) /= det_rat;
  }

  // remaining elements
  for (int i=0; i<M.rows(); i++)
  {
     i_new = (i<k ? i : i+1);
     for (int j=0; j<M.rows(); j++)
     {
      j_new = (j<k ? j : j+1);
      M_new(i_new,j_new) = M(i,j) + det_rat*M_new(i_new,k)*M_new(k,j_new);
     }
  }
   swap(M_new, M);
   return;
}


void compute_M_down(int k, Mat &M)
{
   Mat M_new(M.rows()-1, M.rows()-1);
   int i_old, j_old;

   //
   for (int i=0; i<M_new.rows(); i++)
   {
      i_old = (i<k ? i : i+1);
      for (int j=0; j<M_new.rows(); j++)
      {
         j_old = (j<k ? j : j+1);
         M_new(i,j) = M(i_old, j_old)-M(i_old,k)*M(k,j_old)/M(k,k);
      }
   }
   swap(M, M_new);
}


template <class G, class S, class V> void compute_M_insert_anti(times &anti_segment,
   int s, int r, Mat &M, S& segments_old, G &F, double Beta, double det_rat, V &R)
{
   //
   Mat M_new(M.rows()+1,M.rows()+1);
   std::vector<double> F_kp1(R.size()), L(R.size());

   //
   typename S::iterator it=segments_old.begin();
   for (int i=0; i<F_kp1.size(); i++)
   {
      F_kp1[i]=interpolate_F(it->t_end()-anti_segment.t_end(), Beta, F);
      it++;
   }

   //
   for (int i=0; i<L.size(); i++)
   {
      L[i]=0;
      for (int l=0; l<L.size(); l++) L[i] += M(i,l)*F_kp1[l];
   }

   //
   int i_new, j_new;
   int size=M.rows();

   // element (k+1,k)
   M_new(r,s) = -1./det_rat;

   //
   if (r!=0) // segments remain in the usual order
   {
      // row k+1 and column k
      for (int i=0; i<size; i++)
      {
        i_new = (i<r ? i : i+1);
        j_new = (i<s ? i : i+1);
        M_new(i_new,s) = L[i]/det_rat;
        M_new(r,j_new) = R[i]/det_rat;
      }
      // remaining elements
      for (int i=0; i<size; i++)
      {
        i_new = (i<r ? i : i+1);
        for (int j=0; j<size; j++)
        {
            j_new = (j<s ? j : j+1);
            M_new(i_new,j_new) = M(i,j) - L[i]*R[j]/det_rat;
        }
      }
   }
   else // need to permute indices of R, L, M
   {
      // row k+1 and column k
      for (int i=0; i<size; i++)
      {
         i_new = (i<r ? i : i+1);
         j_new = (i<s ? i : i+1);
         M_new(i_new,s) = L[i]/det_rat;
         M_new(r,j_new) = R[cycle(i,size)]/det_rat;
      }

      // remaining elements
      for (int i=0; i<size; i++)
      {

         i_new = (i<r ? i : i+1);
         for (int j=0; j<size; j++)
         {
            j_new = (j<s ? j : j+1);
            M_new(i_new,j_new) = M(i,cycle(j,size)) - L[i]*R[cycle(j,size)]/det_rat;
         }
      }
   }

   //
   swap(M_new, M);
   return;
}


void compute_M_remove_anti(Mat &M, int s, int r)
{
   //
   Mat M_new(M.rows()-1,M.rows()-1);
   int i_old, j_old;
   int size=M_new.rows();

   if(r!=0) // order of segments remains unchanged
   {
      for (int i=0; i<size; i++)
      {
         i_old = (i<r ? i : i+1);
         for (int j=0; j<size; j++)
         {
            j_old = (j<s ? j : j+1);
            M_new(i,j) = M(i_old,j_old) - M(i_old, s)*M(r, j_old)/M(r, s);
         }
      }
   }
   else // need to permute indices of M
   {
      for (int i=0; i<size; i++)
      {
         for (int j=0; j<size; j++) M_new(i,cycle(j,size)) = M(i+1,j) - M(i+1, s)*M(r, j)/M(r, s);
      }
   }

   //
   swap(M_new, M);
   return;
}


template <class G, class S> void compute_M_shift( times &new_segment, int k, Mat &M,
   S &segments_old, G &F, double Beta, double det_rat)
{
   //
   std::vector<double> R(M.rows(),0), M_k(M.rows(),0), Fe(M.rows(),0);
   typename S::iterator it=segments_old.begin();

   //
   for (int i=0; i<M_k.size(); i++)
   {
      M_k[i] = M(i,k);
      Fe[i] = interpolate_F(new_segment.t_end()-it->t_start(), Beta, F);
      it++;
   }

   //
   for (int i=0; i<R.size(); i++)
   {
      if (i!=k)
      {
         for (int j=0; j<R.size(); j++) R[i] += Fe[j]*M(j,i);
      }
   }

   //
   for (int m=0; m<M.rows(); m++)
   {
      if (m!=k)
      {
         for (int n=0; n<M.rows(); n++) M(n,m) -= M_k[n]*R[m]/det_rat;
      }
      else
      {
         for (int n=0; n<M.rows(); n++) M(n,m) = M_k[n]/det_rat;
      }
   }

   //
   return;
}


//------------------------------------------------------------------------------


template <class G, class S, class V> double det_rat_up(times &new_segment, Mat  M,
   S &segments_old, G &F, V &Fs, V &Fe, double Beta, double &det_rat_sign, double &overlap)
{
   //
   typename S::iterator it=segments_old.begin();
   for (int i=0; i<segments_old.size(); i++)
   {
      Fe[i] = interpolate_F(new_segment.t_end()-it->t_start(), Beta, F);
      Fs[i] = interpolate_F(it->t_end()-new_segment.t_start(), Beta, F);
      it++;
   }
   double det_rat = interpolate_F(new_segment.t_end()-new_segment.t_start(), Beta, F);

   //
   for (int i=0; i<M.rows(); i++)
   {
      for (int j=0; j<M.rows(); j++) det_rat -= Fe[i]*M(i,j)*Fs[j];
   }

   // take care of sign changes produced by segments which "wind around"
   if (new_segment.t_end() < new_segment.t_start())
   {
      det_rat *= -1;
      overlap = -1;
   }
   else
   {
      overlap = 1;
   }

   //
   if (det_rat < 0)
   {
      det_rat_sign = -1;
      det_rat *= -1;
   }
   else
   {
      det_rat_sign = 1;
   }

   //
   return det_rat;
}


template <class S> double det_rat_down(int k, Mat &M, S &segments_old, double &det_rat_sign)
{
   //
   double det_rat = M(k,k);

   // take care of sign changes produced by segments which "wind around"
   if (k==segments_old.size()-1)
   {
      typename S::iterator it=segments_old.end(); it--;
      if (it->t_end() < it->t_start()) det_rat *= -1;
   }

   if (det_rat < 0)
   {
      det_rat_sign = -1;
      det_rat *= -1;
   }
   else
   {
      det_rat_sign = 1;
   }

   //
   return det_rat;
}


template <class G, class S, class V> double det_rat_insert_anti(times &anti_segment,
   Mat &M, S &segments_old, G &F, double Beta, double &det_rat_sign, double &overlap, V &R)
{
   //
   std::vector<double> F_k(R.size());
   typename S::iterator it=segments_old.begin();

   //
   for (int i=0; i<F_k.size(); i++)
   {
      F_k[i]=interpolate_F(anti_segment.t_start()-it->t_start(), Beta, F);
      it++;
   }
   double det_rat = -interpolate_F(anti_segment.t_start()-anti_segment.t_end(), Beta, F);

   //
   it=segments_old.begin();
   for (int i=0; i<R.size(); i++)
   {
      R[i]=0;
      for (int l=0; l<R.size(); l++) R[i] += F_k[l]*M(l,i);
      det_rat += interpolate_F(it->t_end()-anti_segment.t_end(), Beta, F)*R[i];
      it++;
   }

   //
   // take care of sign changes produced by segments which "wind around"
   // check if anti-segment winds around
   overlap = 1;
   if (anti_segment.t_end()<anti_segment.t_start())
   {
      det_rat *= -1;
      overlap = -1;
   }

   //
   if (det_rat < 0)
   {
      det_rat_sign = -1;
      det_rat *= -1;
   }
   else
   {
      det_rat_sign = 1;
   }

   //
   return det_rat;
}


template <class G, class S> double det_rat_remove_anti(times anti_segment, int r,
   int s, Mat &M, S &segments_old, G &F, double Beta, double &det_rat_sign)
{
   // r is the index of the segment which is removed
   // s is the index of the segment which is shifted

   typename S::iterator it=segments_old.begin();
   typename S::iterator its(it), itr(it);
   advance(its, s);
   advance(itr, r);
   double inv_det_rat = -interpolate_F(its->t_end()-itr->t_start(), Beta, F);

   //
   for (int i=0; i<segments_old.size(); i++)
   {
      if (i!=s) inv_det_rat -= interpolate_F(it->t_end()-itr->t_start(), Beta, F)*M(r,i)/M(r,s);
      it++;
   }

   // take care of sign changes produced by segments which "wind around"
   if (anti_segment.t_end() < anti_segment.t_start())inv_det_rat *= -1;

   //
   if (inv_det_rat < 0)
   {
      det_rat_sign = -1;
      inv_det_rat *= -1;
   }
   else
   {
      det_rat_sign = 1;
   }

   //
   return 1/inv_det_rat;
}


template <class G, class S> double det_rat_shift(times &new_segment, int k, Mat &M,
   S &segments_old, G &F, double Beta, double &det_rat_sign, double &overlap)
{
   //
   typename S::iterator it;
   double det_rat = 0;

   //
   it=segments_old.begin();
   for (int i=0; i<M.rows(); i++)
   {
      det_rat += interpolate_F(new_segment.t_end()-it->t_start(), Beta, F)*M(i,k);
      it++;
   }

   // take care of sign changes produced by segments which "wind around"
   overlap = 1;
   if (k==segments_old.size()-1)
   {
      it--;
      // check if last segment has been shifted across Beta
      if ((new_segment.t_end()-new_segment.t_start())*(it->t_end()-it->t_start())<0)
      {
         det_rat *= -1;
         overlap = -1;
      }
   }

   //
   if (det_rat < 0)
   {
      det_rat_sign = -1;
      det_rat *= -1;
   }
   else
   {
      det_rat_sign = 1;
   }

   //
   return det_rat;
}


//------------------------------------------------------------------------------


inline double H(double tau, double Beta, std::vector<double>& Kab )
{
   if (tau<0) tau += Beta;

   //
   int Ntau_m1 = Kab.size()-1;
   double dt = Beta / Ntau_m1;

   // linear mesh by default
   double i = tau / dt;

   // dense mesh optional
   if( (!tau_uniform_K) || (!tau_uniform_Kp) )
   {
      // generic point on the exponential function
      double b = (Beta/2) / ( exp(dt*Ntau_m1/2) - 1 );
      if( tau<=Beta/2 )
      {
         i = log(tau/b+1) / dt;
      }
      else
      {
         i = Ntau_m1 - log((Beta-tau)/b+1) / dt;
      }

      //// boundary even indexes containing an odd one
      //int ndx_1 = (int)i;
      //if(ndx_1%2==1) ndx_1--;
      //int ndx_2 = ndx_1+2;

      //// points on the exponential function corresponding to the indexes
      //double t1 = b * ( exp(dt*ndx_1) - 1 );
      //double t2 = b * ( exp(dt*ndx_2) - 1 );

      //// new i interpolated on the line connecting the two points on the exponential function
      //i = (tau-t1)/(t2-t1)*(ndx_2-ndx_1) + ndx_1;
   }

   // this is the same
   int i_lower=(int) i;

   // interpolate linearly between i_lower and i_lower+1
   double Kinterp = Kab[i_lower] + (i-i_lower)*(Kab[i_lower+1]-Kab[i_lower]);

   //
   return Kinterp;
}


template <class S>  double nonlocal(double ts, double te, S &other_segments,
   double Beta, int this_flavor, VecVecVec &K_table, int insert_remove)
{
   double nonloc=0.0;
   for (int flavor=0; flavor<other_segments.size(); flavor++)
   {
      for(std::set<times>::iterator it=other_segments[flavor].begin(); it!=other_segments[flavor].end(); it++)
      {
         /*
         nonloc += -H(it->t_end()-te  , Beta, K_table[flavor][this_flavor])
                   +H(it->t_end()-ts  , Beta, K_table[flavor][this_flavor])
                   +H(it->t_start()-te, Beta, K_table[flavor][this_flavor])
                   -H(it->t_start()-ts, Beta, K_table[flavor][this_flavor]);
         */
         int diag_flavor  = std::max(flavor,this_flavor);
         int offdiag_flavor = std::min(flavor,this_flavor);
         nonloc += -H(it->t_end()-te  , Beta, K_table[diag_flavor][offdiag_flavor])
                   +H(it->t_end()-ts  , Beta, K_table[diag_flavor][offdiag_flavor])
                   +H(it->t_start()-te, Beta, K_table[diag_flavor][offdiag_flavor])
                   -H(it->t_start()-ts, Beta, K_table[diag_flavor][offdiag_flavor]);
      }
   }
   //
   if (insert_remove==0) // insert
   {
      nonloc += H(te-ts, Beta, K_table[this_flavor][this_flavor]);
   }
   else // remove
   {
      nonloc -= -2*H(0, Beta, K_table[this_flavor][this_flavor])+H(te-ts, Beta, K_table[this_flavor][this_flavor]);
   }
   //
   return nonloc;
}


template <class S> double nonlocal_shift(double te_ins, double te_rem, S &other_segments,
   double Beta, int this_flavor, VecVecVec &K_table)
{
   double nonloc=0.;
   //
   for (int flavor=0; flavor<other_segments.size(); flavor++)
   {
      for(std::set<times>::iterator it=other_segments[flavor].begin(); it!=other_segments[flavor].end(); it++)
      {
         /*
         nonloc += -H(it->t_end()-te_ins, Beta, K_table[flavor][this_flavor])+H(it->t_start()-te_ins, Beta, K_table[flavor][this_flavor]);
         nonloc -= -H(it->t_end()-te_rem, Beta, K_table[flavor][this_flavor])+H(it->t_start()-te_rem, Beta, K_table[flavor][this_flavor]);
         */
         int diag_flavor  = std::max(flavor,this_flavor);
         int offdiag_flavor = std::min(flavor,this_flavor);
         nonloc += -H(it->t_end()-te_ins, Beta, K_table[diag_flavor][offdiag_flavor])+H(it->t_start()-te_ins, Beta, K_table[diag_flavor][offdiag_flavor]);
         nonloc -= -H(it->t_end()-te_rem, Beta, K_table[diag_flavor][offdiag_flavor])+H(it->t_start()-te_rem, Beta, K_table[diag_flavor][offdiag_flavor]);
      }
   }
   //
   nonloc -= -H(te_rem-te_ins, Beta, K_table[this_flavor][this_flavor]); // inexistent bond
   nonloc += -H(te_rem-te_rem, Beta, K_table[this_flavor][this_flavor]); // local contribution at te doesn't change
   //
   return nonloc;
}


//------------------------------------------------------------------------------


double get_nj( segment_container_t &segment, int &full_line, double &tau )
{
   //
   double nj = 0.0;

   if(segment.size()==0)
   {
      nj = ( full_line==1 ? 1.0 : 0.0 );
   }
   else
   {
      //this is the segment that starts after tau
      std::set<times>::const_iterator it_after = segment.upper_bound( times(tau,0.0) );
      if( it_after==segment.end() ) it_after = segment.begin();

      //this is the segment that has the start time before tau
      std::set<times>::const_iterator it_before = it_after;
      if( it_before==segment.begin() ) it_before = segment.end();
      it_before--;

      //standard overlap
      bool overlap_std = (it_before->t_end() > it_before->t_start()) && ( tau > it_before->t_start() && tau < it_before->t_end() );
      //wrapped overlap
      bool overlap_wrp = (it_before->t_end() < it_before->t_start()) && ( tau > it_before->t_start() || tau < it_before->t_end() );

      //
      nj = ( overlap_std||overlap_wrp ? 1.0 : 0.0 );
   }

   //
   return nj;
}


//============================================================================//


#endif
