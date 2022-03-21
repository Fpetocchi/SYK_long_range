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



//----------------------------------------------------------------------------//
//                               SYMMETRIZATIONS                              //
//----------------------------------------------------------------------------//
void spin_symm( Vec &Vec_in )
{
   int Nflavor = Vec_in.size();
   for (int ifl=0; ifl<(Nflavor/2); ++ifl)
   {
      //
      int up = 2*ifl;
      int dw = 2*ifl+1;
      //
      //if(para_mode == 2)
      //{
      //   std::vector<double> shift = { abs(Vec_in[up]-0.5), abs(Vec_in[dw]-0.5) };
      //   int ndx = std::min_element( shift.begin(), shift.end() ) - shift.begin();
      //   if(ndx==0) dw = up;
      //   if(ndx==1) up = dw;
      //}
      //
      double val = ( Vec_in[up] + Vec_in[dw] ) /2.0;
      Vec_in[up] = val;
      Vec_in[dw] = val;
   }
}

void spin_symm( VecVec &VecVec_in )
{
   int Nflavor = VecVec_in.size();
   int Ntau = VecVec_in[0].size();
   for (int ifl=0; ifl<(Nflavor/2); ++ifl)
   {
      //
      int up = 2*ifl;
      int dw = 2*ifl+1;
      //
      for (int i=0; i<Ntau; ++i)
      {
         double val = ( VecVec_in[up][i] + VecVec_in[dw][i] ) /2.0;
         VecVec_in[up][i] = val;
         VecVec_in[dw][i] = val;
      }
   }
}

void orb_symm( Vec &Vec_in, std::vector<std::vector<int>> &Lists )
{
   for (int ilist=0; ilist<Lists.size(); ++ilist)
   {
      double val_up=0.0;
      double val_dw=0.0;
      for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
      {
         int iorb = Lists[ilist][iobj];
         val_up += Vec_in[2*iorb]/Lists[ilist].size();
         val_dw += Vec_in[2*iorb+1]/Lists[ilist].size();
      }
      for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
      {
         int iorb = Lists[ilist][iobj];
         Vec_in[2*iorb]=val_up;
         Vec_in[2*iorb+1]=val_dw;
      }
   }
}

void orb_symm( VecVec &VecVec_in, std::vector<std::vector<int>> &Lists, int Norb=0 )
{
   //
   int Ntau = VecVec_in[0].size();

   //
   if(Norb==0) //Green's function like
   {
      for (int i=0; i<Ntau; ++i)
      {
         for (int ilist=0; ilist<Lists.size(); ++ilist)
         {
            double val_up=0.0;
            double val_dw=0.0;
            for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
            {
               int iorb = Lists[ilist][iobj];
               val_up += VecVec_in[2*iorb][i]/Lists[ilist].size();
               val_dw += VecVec_in[2*iorb+1][i]/Lists[ilist].size();
            }
            for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
            {
               int iorb = Lists[ilist][iobj];
               VecVec_in[2*iorb][i]=val_up;
               VecVec_in[2*iorb+1][i]=val_dw;
            }
         }
      }
   }
   else //n(tau)n(0) correlator like
   {
      //enlarge list with missing orbitals
      for (int iorb=0; iorb<Norb; ++iorb)
      {
         bool found=false;
         for (int ilist=0; ilist<Lists.size(); ++ilist)
         {
            for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
            {
               if(iorb == Lists[ilist][iobj]) found=true;
            }
            if(found)break;
         }
         if(!found)Lists.push_back(std::vector<int>(1,iorb));
      }

      //symmetrization at each itau
      int Nflavor=2*Norb;
      for (int i=0; i<Ntau; ++i)
      {
         //reshuffle vector in [Norb]*Nspin block shape
         VecVec Wmat(Nflavor,Vec(Nflavor,0.0));
         int position=0;
         for (int ifl=0; ifl<Nflavor; ++ifl)
         {
            for (int jfl=0; jfl<=ifl; ++jfl)
            {
               // iorb = ifl/2; ispin = ifl%2
               int io = (ifl/2) + Norb*(ifl%2);
               int jo = (jfl/2) + Norb*(jfl%2);
               Wmat[io][jo] =  VecVec_in[position][i];
               Wmat[jo][io] =  Wmat[io][jo];
               position++;
            }
         }

         //loop over indexes within diagonal blocks
         for (int ilist=0; ilist<Lists.size(); ++ilist)
         {
            double D_uu=0.0,D_ud=0.0;
            double D_du=0.0,D_dd=0.0;
            double OD_uu=0.0,OD_ud=0.0;
            double OD_du=0.0,OD_dd=0.0;
            int O_dim = Lists[ilist].size();
            int OD_dim = O_dim*(O_dim-1);
            //compute average
            for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
            {
               for (int jobj=0; jobj<Lists[ilist].size(); ++jobj)
               {
                  int iorb = Lists[ilist][iobj];
                  int jorb = Lists[ilist][jobj];
                  if(iorb==jorb)
                  {
                     D_uu += Wmat[iorb][iorb]/O_dim;
                     D_ud += Wmat[iorb][iorb+Norb]/O_dim;
                     D_du += Wmat[iorb+Norb][iorb]/O_dim;
                     D_dd += Wmat[iorb+Norb][iorb+Norb]/O_dim;
                  }
                  else
                  {
                     OD_uu += Wmat[iorb][jorb]/OD_dim;
                     OD_ud += Wmat[iorb][jorb+Norb]/OD_dim;
                     OD_du += Wmat[iorb+Norb][jorb]/OD_dim;
                     OD_dd += Wmat[iorb+Norb][jorb+Norb]/OD_dim;
                  }
               }
            }
            //reinsert average
            for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
            {
               for (int jobj=0; jobj<Lists[ilist].size(); ++jobj)
               {
                  int iorb = Lists[ilist][iobj];
                  int jorb = Lists[ilist][jobj];
                  if(iorb==jorb)
                  {
                     Wmat[iorb][iorb] = D_uu;
                     Wmat[iorb][iorb+Norb] = D_ud;
                     Wmat[iorb+Norb][iorb] = D_du;
                     Wmat[iorb+Norb][iorb+Norb] = D_dd;
                  }
                  else
                  {
                      Wmat[iorb][jorb] = OD_uu;
                      Wmat[iorb][jorb+Norb] = OD_ud;
                      Wmat[iorb+Norb][jorb] = OD_du;
                      Wmat[iorb+Norb][jorb+Norb] = OD_dd;
                  }
               }
            }
         }

         //loop over indexes within off-diagonal blocks
         for (int ilist=0; ilist<Lists.size(); ++ilist)
         {
            for (int jlist=0; jlist<Lists.size(); ++jlist)
            {
               if(ilist!=jlist)
               {
                  double OD_uu=0.0,OD_ud=0.0;
                  double OD_du=0.0,OD_dd=0.0;
                  int OD_dim = Lists[ilist].size()*Lists[jlist].size();
                  //compute average
                  for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
                  {
                     for (int jobj=0; jobj<Lists[jlist].size(); ++jobj)
                     {
                        int iorb = Lists[ilist][iobj];
                        int jorb = Lists[jlist][jobj];
                        OD_uu += Wmat[iorb][jorb]/OD_dim;
                        OD_ud += Wmat[iorb][jorb+Norb]/OD_dim;
                        OD_du += Wmat[iorb+Norb][jorb]/OD_dim;
                        OD_dd += Wmat[iorb+Norb][jorb+Norb]/OD_dim;
                     }
                  }
                  //reinsert average
                  for (int iobj=0; iobj<Lists[ilist].size(); ++iobj)
                  {
                     for (int jobj=0; jobj<Lists[jlist].size(); ++jobj)
                     {
                        int iorb = Lists[ilist][iobj];
                        int jorb = Lists[jlist][jobj];
                        Wmat[iorb][jorb] = OD_uu;
                        Wmat[iorb][jorb+Norb] = OD_ud;
                        Wmat[iorb+Norb][jorb] = OD_du;
                        Wmat[iorb+Norb][jorb+Norb] = OD_dd;
                     }
                  }
               }
            }
         }

         //put back [Norb]*Nspin block into vector
         position=0;
         for (int ifl=0; ifl<Nflavor; ++ifl)
         {
            for (int jfl=0; jfl<=ifl; ++jfl)
            {
               // iorb = ifl/2; ispin = ifl%2
               int io = (ifl/2) + Norb*(ifl%2);
               int jo = (jfl/2) + Norb*(jfl%2);
               VecVec_in[position][i] = Wmat[io][jo];
               position++;
            }
         }

      } // itau loop
   }
}



//----------------------------------------------------------------------------//
//                                ACCUMULATIONS                               //
//----------------------------------------------------------------------------//
void accumulate_Vec( Vec &Vec_out, Vec &Vec_in, double Norm=1.0)
{
   int d1 = Vec_out.size();
   for (int i=0; i<d1; i++) Vec_out[i] += Vec_in[i]*Norm;
}

void accumulate_VecVec( VecVec &VecVec_out, VecVec &VecVec_in, double Norm=1.0)
{
   int d1 = VecVec_out.size();
   int d2 = VecVec_out[0].size();
   for (int i=0; i<d1; i++)
   {
      for (int j=0; j<d2; j++)
      {
         VecVec_out[i][j] += VecVec_in[i][j]*Norm;
      }
   }
}



//----------------------------------------------------------------------------//
//                              GREEN'S FUNCTION                              //
//----------------------------------------------------------------------------//
void measure_G( Vec &G, segment_container_t &segment, Mat &M, double &Beta, double Norm=1.0 )
{
   //
   std::set<times>::iterator it1, it2;
   int Ntau = G.size();

   //
   for (int i=0; i<M.cols(); i++)
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
            G[index] -= M(k,i)*bubble_sign/(Beta*Beta);
            //safety checks
            bool nanCond = std::isnan(G[index]);
            bool infCond = std::isinf(G[index]);
            if( nanCond || infCond )
            {
               if( nanCond ) printf(" NaN G[%d] \n",index);
               if( infCond ) printf(" Inf G[%d] \n",index);
               printf(" M[%d,%d]= %f \n",k,i,M(k,i));
               printf(" bubble_sign= %f \n",bubble_sign);
               printf(" Beta= %f \n\n\n",Beta);
            }
            //
         }
      }
   }
   //
   if(Norm!=1.0)
   {
      for (int i=0; i<Ntau; i++) G[i] *= Norm;
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



//----------------------------------------------------------------------------//
//                                   N(tau)                                   //
//----------------------------------------------------------------------------//
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
             index = it->t_start()/Beta*(Ntau-1)+1;
             n_tau[ifl][index] *= -1;
             index = it->t_end()/Beta*(Ntau-1)+1;
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



//----------------------------------------------------------------------------//
//                                 N(tau)N(0)                                 //
//----------------------------------------------------------------------------//
VecVec measure_nnt( std::vector<segment_container_t> &segments, std::vector<int> &full_line, int &Ntau, double &Beta)
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
               if (j>Ntau) j -= (Ntau-1);
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
         for (int i=0; i<Ntau; ++i)
         {
            for (int index=0; index<Ntau; ++index)
            {
               int j=i+index;
               if (j>Ntau-1) j -= (Ntau-1);
               nn_corr_meas[position][index] += n_tau[ifl][i]*n_tau[jfl][j] / (double)Ntau;
            }
         }
         position++;
      }
   }
   //
   return nn_corr_meas;
}

void accumulate_nnt( VecVec &nn_corr_meas, VecVec &n_tau, int nnt_shift=0 )
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();

   //
   int position=0;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      double ni = std::accumulate(n_tau[ifl].begin(), n_tau[ifl].end(), 0.0) * nnt_shift/(Ntau-1) ;
      for (int jfl=0; jfl<=ifl; ++jfl)
      {
         double nj = std::accumulate(n_tau[jfl].begin(), n_tau[jfl].end(), 0.0) * nnt_shift/(Ntau-1) ;
         for (int i=0; i<Ntau; ++i)
         {
            for (int index=0; index<Ntau; ++index)
            {
               int j=i+index;
               if (j>Ntau-1) j -= (Ntau-1);
               nn_corr_meas[position][index] += (n_tau[ifl][i]*n_tau[jfl][j] - ni*nj ) / (double)Ntau;
            }
         }
         position++;
      }
   }
}

void accumulate_nnt( VecVec &nn_corr_meas, VecVec &n_tau, double s, int Nmeas, int nnt_shift=0 )
{
   //
   int Nflavor = n_tau.size();
   int Ntau = n_tau[0].size();
   double coeff=s/Nmeas;

   //
   int position=0;
   for (int ifl=0; ifl<Nflavor; ++ifl)
   {
      for (int jfl=0; jfl<=ifl; ++jfl)
      {
         for (int imeas=0; imeas<Nmeas; ++imeas)
         {
            int m=(int)(rndm()*(Ntau-1));
            if(n_tau[jfl][m]>0)
            {
         	  for(int n=m; n<Ntau-1; ++n) nn_corr_meas[position][n-m]+=n_tau[ifl][n]*coeff;
         	  for(int n=0; n<m+1; ++n)  nn_corr_meas[position][n-m+Ntau-1]+=n_tau[ifl][n]*coeff;
         	}
         }
         position++;
      }
   }
}



//----------------------------------------------------------------------------//
//                                 HISTOGRAMS                                 //
//----------------------------------------------------------------------------//
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
