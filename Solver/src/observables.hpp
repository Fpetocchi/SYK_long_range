#ifndef ___OBSERVABLES___
#define ___OBSERVABLES___

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
#include <algorithm>
#include <set>
#include "segments.hpp"


//============================================================================//


template<typename T> Vec normalize_Vec( Vec &Vec_in, T &Norm)
{
   int Nrows=Vec_in.size();
   Vec Vec_out(Nrows,0.0);
   for (int irow=0; irow<Nrows; irow++) Vec_out[irow]=Vec_in[irow]/(double) Norm;
   return Vec_out;
}
template<typename T> VecVec normalize_VecVec( VecVec &VecVec_in, T &Norm)
{
   int Ncols=VecVec_in.size();
   int Nrows=VecVec_in[0].size();
   VecVec VecVec_out(Ncols,Vec(Nrows,0.0));
   for (int icol=0; icol<Ncols; icol++)
   {
      for (int irow=0; irow<Nrows; irow++) VecVec_out[icol][irow]=VecVec_in[icol][irow]/(double) Norm;
   }
   return VecVec_out;
}
template<typename T> VecVec normalize_VecVec( VecVec &VecVec_in, std::vector<T> &NormCol)
{
   int Ncols=VecVec_in.size();
   int Nrows=VecVec_in[0].size();
   VecVec VecVec_out(Ncols,Vec(Nrows,0.0));
   for (int icol=0; icol<Ncols; icol++)
   {
      for (int irow=0; irow<Nrows; irow++) VecVec_out[icol][irow]=VecVec_in[icol][irow]/(double) NormCol[icol];
   }
   return VecVec_out;
}


//------------------------------------------------------------------------------


void measure_G( Vec &G, segment_container_t &segment, Mat &M, int &Ntau, double &Beta, double Norm=1.0)
{
   //
   std::set<times>::iterator it1, it2;

   //
   for (int i=0; i<M.rows(); i++)
   {
      //
      (i==0 ? it1 = segment.begin() : it1++);
      //
      for (int k=0; k<M.rows(); k++)
      {
         //
         (k==0 ? it2 = segment.begin() : it2++);
         //
         if (M(k,i)!=0)
         {
            double argument = it1->t_end()-it2->t_start();
            double bubble_sign=1;
            //
            if (argument > 0)
            {
               bubble_sign = 1;
            }
            else
            {
               bubble_sign = -1;
               argument += Beta;
            }
            //
            int index = argument/Beta*(Ntau-1)+0.5;
            G[index] -= Norm*M(k,i)*bubble_sign/(Beta*Beta);
         }
      }
   }
}


void accumulate_G( VecVec &G, VecVec &Gtmp)
{
   //
   int Nflavor = Gtmp.size();
   int Ntau = Gtmp[0].size();

   //
   for (int ifl=0; ifl<Nflavor; ifl++)
   {
      for (int itau=0; itau<Ntau; itau++)
      {
         G[ifl][itau] += Gtmp[ifl][itau];
      }
   }
}


void binAverageVec( std::vector<int> &bins, VecVec &G, VecVec &Gerr)
{
   //
   int Nflavor = G.size();
   int Ntau = G[0].size();
   int binlength = bins[0];
   int binstart = bins[1];

   // bin average - v1 - less efficient memory wise but safer
   if(binlength>0)
   {
      double Nsample = (double)(2*binlength +1);
      for (int ifl=0; ifl<Nflavor; ifl++)
      {
         //
         Vec Gsym(Ntau,0.0);
         for (int itau=binlength+binstart; itau<Ntau-binlength-binstart; itau++)
         {
            //
            double avrg=0.0;
            for(int ibin=itau-binlength; ibin<=itau+binlength; ibin++) avrg += G[ifl][ibin] / Nsample;
            //
            Gsym[itau]=avrg;
            //
            double stderr=0.0;
            for(int ibin=itau-binlength; ibin<=itau+binlength; ibin++) stderr += pow((G[ifl][ibin]-avrg),2) / (Nsample-1);
            Gerr[ifl][itau] = sqrt(stderr);
         }
         //
         for (int itau=binlength+binstart; itau<Ntau-binlength-binstart; itau++) G[ifl][itau] = Gsym[itau];
      }
   }

   // bin average - v2
   /*
   if(binlength>0)
   {
      double Nsample = (double)(2*binlength +1);
      for (int ifl=0; ifl<Nflavor; ifl++)
      {
         Vec Gsym(binlength+1,0.0);
         for (int itau=binlength; itau<Ntau-binlength; itau++)
         {
            //
            double avrg=0.0;
            for(int ibin=itau-binlength; ibin<=itau+binlength; ibin++) avrg += G[ifl][ibin] / Nsample;
            if(itau>2*binlength)
            {
               G[ifl][itau-(binlength+1)] = Gsym[0];
               Gsym.erase(Gsym.begin());
               Gsym.push_back(avrg);
               // put the reamaining
               if(itau==Ntau-binlength-1)
               {
                  for (int ibin=binlength; ibin<binlength+1; ibin++) G[ifl][itau-(binlength+1)+ibin] = Gsym[ibin];
               }
            }
            else
            {
               Gsym.push_back(avrg);
            }
            //
            double stderr=0.0;
            for(int ibin=itau-binlength; ibin<=itau+binlength; ibin++) stderr += pow((G[ifl][ibin]-avrg),2) / (Nsample-1);
            Gerr[ifl][itau] = sqrt(stderr);
         }
      }
   }
   */
}


//------------------------------------------------------------------------------


void spin_symm( Vec &Vec )
{
   int Nflavor = Vec.size();
   for (int ifl=0; ifl<(Nflavor/2); ++ifl)
   {
      double val = ( Vec[2*ifl] + Vec[2*ifl+1] ) /2.0;
      Vec[2*ifl] = val;
      Vec[2*ifl+1] = val;
   }
}

void spin_symm( VecVec &VecVec )
{
   int Nflavor = VecVec.size();
   int Ntau = VecVec[0].size();
   for (int i=0; i<Ntau; ++i)
   {
      for (int ifl=0; ifl<(Nflavor/2); ++ifl)
      {
         double val = ( VecVec[2*ifl][i] + VecVec[2*ifl+1][i] ) /2.0;
         VecVec[2*ifl][i] = val;
         VecVec[2*ifl+1][i] = val;
      }
   }
}


//------------------------------------------------------------------------------


VecVec measure_nt( std::vector<segment_container_t> &segments, std::vector<int> &full_line, int &Ntau, double &Beta)
{
   //
   int Nflavor = segments.size();
   VecVec n_tau(Nflavor,Vec(Ntau,1));
   std::set<times>::iterator it;

   //
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      if (segments[ifl].size()==0)
      {
         if (full_line[ifl]==0)
            for (int i=0; i<(int)n_tau[ifl].size(); ++i) n_tau[ifl][i]=0;
      }
      else
      {
         it=segments[ifl].end();
         it--;
         //
         if (it->t_end()<it->t_start())
            n_tau[ifl][0]=1;
         else
            n_tau[ifl][0]=0;
         //
         // mark segment start and end points
         int index;
         for (it=segments[ifl].begin(); it!=segments[ifl].end(); it++)
         {
             index = it->t_start()/Beta*Ntau;
             n_tau[ifl][index] *= -1;
             index = it->t_end()/Beta*Ntau;
             n_tau[ifl][index] *= -1;
         }
         //
         // fill vector with occupation number
         for (int i=1; i<(int)n_tau[ifl].size(); i++)
         {
            if (n_tau[ifl][i]==-1)
               n_tau[ifl][i]=1-n_tau[ifl][i-1];
            else
               n_tau[ifl][i]=n_tau[ifl][i-1];
         }
      }
   }
   //
   return n_tau;
}


void accumulate_nt( VecVec &n_meas, VecVec &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();
   //
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      for (int i=0; i<Ntau; ++i)
      {
         n_meas[ifl][i]+=n_tau[ifl][i];
      }
   }
}

//------------------------------------------------------------------------------


VecVec measure_nnt_standalone( std::vector<segment_container_t> &segments, std::vector<int> &full_line, int &Ntau, double &Beta)
{
   //
   int Nflavor = segments.size();
   VecVec n_vectors(Nflavor,Vec(Ntau,1));
   std::set<times>::iterator it;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      if (segments[ifl].size()==0)
      {
         if (full_line[ifl]==0)
            for (int i=0; i<(int)n_vectors[ifl].size(); ++i) n_vectors[ifl][i]=0;
      }
      else
      {
         it=segments[ifl].end();
         it--;
         //
         if (it->t_end()<it->t_start())
            n_vectors[ifl][0]=1;
         else
            n_vectors[ifl][0]=0;
         //
         // mark segment start and end points
         int index;
         for (it=segments[ifl].begin(); it!=segments[ifl].end(); it++)
         {
             index = it->t_start()/Beta*Ntau;
             n_vectors[ifl][index] *= -1;
             index = it->t_end()/Beta*Ntau;
             n_vectors[ifl][index] *= -1;
         }
         //
         // fill vector with occupation number
         for (int i=1; i<(int)n_vectors[ifl].size(); i++)
         {
            if (n_vectors[ifl][i]==-1)
               n_vectors[ifl][i]=1-n_vectors[ifl][i-1];
            else
               n_vectors[ifl][i]=n_vectors[ifl][i-1];
         }
      }
   }
   //
   //std::valarray<double> nn_corr_meas(Nflavor*(Nflavor+1)/2*(Ntau));
   VecVec nn_corr_meas(Nflavor*(Nflavor+1)/2,Vec(Ntau,0.0));
   int position=0;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      for (int jfl=0; jfl<=ifl; ++jfl)
      {
         position++;
         for (int i=0; i<Ntau; ++i)
         {
            for (int index=0; index<Ntau; ++index)
            {
               int j=i+index;
               if (j>Ntau) j -= Ntau;
               nn_corr_meas[position][index] += n_vectors[ifl][i]*n_vectors[jfl][j];
            }
         }
         for (int i=0; i<Ntau; ++i)nn_corr_meas[position][i]/=(Ntau);
      }
   }
   //
   return nn_corr_meas;
}


VecVec measure_nnt( VecVec &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();
   VecVec nn_corr_meas(Nflavor*(Nflavor+1)/2,Vec(Ntau,0.0));

   //
   int position=0;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      for (int jfl=0; jfl<=ifl; ++jfl)
      {
         //
         for (int i=0; i<Ntau; ++i)
         {
            for (int index=0; index<Ntau; ++index)
            {
               int j=i+index;
               if ( j>(Ntau-1) ) j -= (Ntau-1);
               nn_corr_meas[position][index] += n_tau[ifl][i]*n_tau[jfl][j];
            }
         }
         for (int i=0; i<Ntau; ++i)nn_corr_meas[position][i]/=(double)Ntau;

         //
         position++;
         //
      }
   }
   //
   return nn_corr_meas;
}


void accumulate_nnt( VecVec &nn_corr_meas, VecVec &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();

   //
   int position=0;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      for (int jfl=0; jfl<=ifl; ++jfl)
      {
         for (int i=0; i<Ntau; ++i)
         {
            for (int index=0; index<Ntau; ++index)
            {
               int j=i+index;
               if (j>Ntau) j -= Ntau;
               nn_corr_meas[position][index] += n_tau[ifl][i]*n_tau[jfl][j] / (double)Ntau;
            }
         }
         position++;
      }
   }
}


//------------------------------------------------------------------------------


Vec measure_Nhist( VecVec &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();
   Vec nhist(Nflavor+1,0.0);
   //
   for (int itau=1; itau<Ntau; itau++)
   {
      int ntmp=0;
      for (int ifl=0; ifl<Nflavor; ifl++) ntmp+=n_tau[ifl][itau];
      nhist[ntmp]+=1./Ntau;
   }
   return nhist;
}


Vec measure_Szhist( VecVec &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();
   Vec szhist(Nflavor/2+1,0.0);
   //
   for (int itau=1; itau<Ntau; itau++)
   {
      int ntmp=0;
      for (int ifl=0; ifl<(Nflavor-1); ifl+=2) ntmp+=n_tau[ifl][itau]-n_tau[ifl+1][itau];
      if (ntmp<0) ntmp*=-1;
      szhist[ntmp]+=1./Ntau;
   }
   return szhist;
}


void accumulate_Nhist( Vec &nhist, VecVec &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();
   //
   for (int itau=1; itau<Ntau; itau++)
   {
      int ntmp=0;
      for (int ifl=0; ifl<Nflavor; ifl++) ntmp+=abs(n_tau[ifl][itau]);
      if(ntmp<nhist.size())nhist[(int)ntmp]+=1./Ntau;
   }
}


void accumulate_Szhist( Vec &szhist, VecVec &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();
   //
   for (int itau=1; itau<Ntau; itau++)
   {
      int ntmp=0;
      for (int ifl=0; ifl<Nflavor; ifl+=2) ntmp+=n_tau[ifl][itau]-n_tau[ifl+1][itau];
      if (ntmp<0) ntmp*=-1;
      if(ntmp<szhist.size())szhist[(int)ntmp]+=1./Ntau;
   }
}


//============================================================================//


#endif
