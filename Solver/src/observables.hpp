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
//
#include <set>
#include "times.hpp"


//============================================================================//


std::vector<double> measure_G( segment_container_t &segment, Mat &M, int &Ntau_p1, double &Beta)
{
   //
   std::vector<double> G(Ntau_p1,0.0);
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
            int index = argument/Beta*(Ntau_p1-1)+0.5;
            //G_meas[ifl*(Ntau+1)+index] += M(k,i)*bubble_sign/(Beta*Beta);
            G[index] += M(k,i)*bubble_sign/(Beta*Beta);
         }
      }
   }

   //
   return G;
}


//------------------------------------------------------------------------------


std::vector<std::vector<double>> measure_nt( std::vector<segment_container_t> &segments, std::vector<int> &full_line, int &Ntau_p1, double &Beta)
{
   //
   int Nflavor = segments.size();
   std::vector<std::vector<double>> n_tau(Nflavor,std::vector<double>(Ntau_p1,1));
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
             index = it->t_start()/Beta*Ntau_p1;
             n_tau[ifl][index] *= -1;
             index = it->t_end()/Beta*Ntau_p1;
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


//------------------------------------------------------------------------------


std::vector<std::vector<double>> measure_nnt( std::vector<std::vector<double>> &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau_p1 = n_tau[0].size();
   std::vector<std::vector<double>> nn_corr_meas(Nflavor*(Nflavor+1)/2,std::vector<double>(Ntau_p1,0.0));

   //
   int position=0;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      for (int jfl=0; jfl<=ifl; ++jfl)
      {
         position++;
         for (int i=0; i<Ntau_p1; ++i)
         {
            for (int index=0; index<Ntau_p1; ++index)
            {
               int j=i+index;
               if (j>Ntau_p1) j -= Ntau_p1;
               nn_corr_meas[position][index] += n_tau[ifl][i]*n_tau[jfl][j];
            }
         }
         for (int i=0; i<Ntau_p1; ++i)nn_corr_meas[position][i]/=(Ntau_p1);
      }
   }
   //
   return nn_corr_meas;
}


//------------------------------------------------------------------------------


std::vector<std::vector<double>> measure_nnt_standalone( std::vector<segment_container_t> &segments, std::vector<int> &full_line, int &Ntau_p1, double &Beta)
{
   //
   int Nflavor = segments.size();
   std::vector<std::vector<double>> n_vectors(Nflavor,std::vector<double>(Ntau_p1,1));
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
             index = it->t_start()/Beta*Ntau_p1;
             n_vectors[ifl][index] *= -1;
             index = it->t_end()/Beta*Ntau_p1;
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
   //std::valarray<double> nn_corr_meas(Nflavor*(Nflavor+1)/2*(Ntau_p1));
   std::vector<std::vector<double>> nn_corr_meas(Nflavor*(Nflavor+1)/2,std::vector<double>(Ntau_p1,0.0));
   int position=0;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      for (int jfl=0; jfl<=ifl; ++jfl)
      {
         position++;
         for (int i=0; i<Ntau_p1; ++i)
         {
            for (int index=0; index<Ntau_p1; ++index)
            {
               int j=i+index;
               if (j>Ntau_p1) j -= Ntau_p1;
               nn_corr_meas[position][index] += n_vectors[ifl][i]*n_vectors[jfl][j];
            }
         }
         for (int i=0; i<Ntau_p1; ++i)nn_corr_meas[position][i]/=(Ntau_p1);
      }
   }

   //
   return nn_corr_meas;
}


//------------------------------------------------------------------------------


std::vector<double> measure_Nhist( std::vector<std::vector<double>> &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau_p1 = n_tau[0].size();
   std::vector<double> nhist(Nflavor+1,0.0);

   //
   for (int itau=1; itau<Ntau_p1; itau++)
   {
      int ntmp=0;
      for (int ifl=0; ifl<Nflavor; ifl++) ntmp+=n_tau[ifl][itau];
      nhist[ntmp]+=1./Ntau_p1;
   }
   return nhist;
}


std::vector<double> measure_Szhist( std::vector<std::vector<double>> &n_tau)
{
   //
   int Nflavor = n_tau.size();
   int Ntau_p1 = n_tau[0].size();
   std::vector<double> szhist(Nflavor/2+1,0.0);

   //
   for (int itau=1; itau<Ntau_p1; itau++)
   {
      int ntmp=0;
      for (int ifl=0; ifl<Nflavor; ifl+=2) ntmp+=n_tau[ifl][itau]-n_tau[ifl+1][itau];
      if (ntmp<0) ntmp*=-1;
      szhist[ntmp]+=1./Ntau_p1;
   }
   return szhist;
}




#endif
