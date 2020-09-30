#ifndef ___MOVES___
#define ___MOVES___

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
#include "segments.hpp"


//============================================================================//


template <class S> void insert_remove_full_line(double Eo, Mat &Umat, double Beta,
   int &full_line, std::vector<S> &other_segments,
   std::vector<int> &other_full_line, int this_flavor)
{
   //
   int insert = (rndm() < 0.5);

   // insert=1(0) means we want to insert(remove) a full line
   if ((insert==1 && full_line==1) || (insert==0 && full_line==0)) return;

   //
   int FLAVORS = other_full_line.size();
   double otherlength_u=0;

   //
   for (int i=0; i<FLAVORS; i++)
   {
      if (i==this_flavor) continue;
      //
      double other_length=0;
      for (typename S::iterator it=other_segments[i].begin(); it!=other_segments[i].end(); it++)
         other_length += (it->t_end()-it->t_start()>0 ? it->t_end()-it->t_start() : it->t_end()-it->t_start()+Beta);
      //
      if (other_full_line[i]==1) other_length = Beta;
      otherlength_u += other_length*Umat(i, this_flavor);
   }

   if (insert) // try to insert full line
   {
      if (log(rndm()) < Beta*Eo-otherlength_u) full_line = 1;
   }
   else // try to remove full line
   {
      if (log(rndm()) < -Beta*Eo+otherlength_u) full_line = 0;
   }

}


//------------------------------------------------------------------------------


template <class S, class G> void insert_remove_antisegment( double t, double Beta,
   double Eo, Mat &Umat, G &F, int &full_line, S &segments, Mat &M, double sign,
   std::vector<S> &other_segments, std::vector<int> &other_full_line,
   int this_flavor, VecVecVec *K_table)
{
   //
   double t_up;                 // distance to next segment up (t_start)
   double t_down;               // distance to next segment down (t_end)
   typename S::iterator s_up;   // iterator of the segment up
   typename S::iterator s_down; // iterator of the segment down

   //
   if(rndm()<0.5) // try to insert an anti-segment
   {
      if(full_line==1)
      {
         //
         t_down = -Beta;
         double length = compute_length(rndm(), Beta, 0);
         double t_end = (t+length < Beta ? t+length : t+length-Beta);
         times segment_insert(t_end, t);
         times segment_remove(t,t_end);

         //
         double log_prob, overlap, det_rat, det_rat_sign;
         std::vector<double> Fs(segments.size()), Fe(segments.size());
         det_rat = det_rat_up(segment_insert, M, segments, F, Fs, Fe, Beta, det_rat_sign, overlap);

         //
         double otherlength_u=0;
         int FLAVORS=other_full_line.size();
         for(int i=0; i<FLAVORS; i++)
         {
            if (i==this_flavor) continue;
            double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], Beta);
            otherlength_u += other_length*Umat(i, this_flavor);
         }

         //
         double nonloc=0.0;
         if(K_table!=NULL)
         {
            VecVecVec K = *K_table;
            nonloc=nonlocal(segment_remove.t_end(), segment_remove.t_start(), other_segments, Beta, this_flavor, K, 0);
         }


         //
         log_prob = log(Beta*Beta*det_rat)-length*Eo+otherlength_u-nonloc;

         //
         if (log(rndm()) < log_prob)
         {
            compute_M_up(0, M, Fs, Fe, det_rat*overlap);
            sign *= det_rat_sign;
            typename S::iterator sit;
            sit=segments.insert(segment_insert).first;
            if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
            full_line = 0;
         }
      }
      else
      {
         //
         compute_intervals(t, Beta, t_up, t_down, segments, s_up, s_down);

         //
         if(t_down<0) // t does lie on a segment -> it's possible to insert an anti-segment starting from t
         {
            //
            double length = compute_length(rndm(), -t_down, 0);
            times segment_shrink(s_down->t_start(),t);

            //
            double t_start = t + length;
            if (t_start > Beta) t_start-=Beta;

            //
            times segment_insert(t_start, s_down->t_end());
            times anti_segment(t,t_start);

            //
            double otherlength_u=0;
            int FLAVORS=other_full_line.size();
            for (int i=0; i<FLAVORS; i++)
            {
                if (i==this_flavor) continue;
                double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], Beta);
                otherlength_u += other_length*Umat(i, this_flavor);
            }

            //
            double nonloc=0.0;
            if(K_table!=NULL)
            {
               VecVecVec K = *K_table;
               nonloc=nonlocal(anti_segment.t_end(), anti_segment.t_start(), other_segments, Beta, this_flavor, K, 0);
            }

            //
            double log_prob, overlap, det_rat, det_rat_sign;
            std::vector<double> R(segments.size());
            det_rat = det_rat_insert_anti(anti_segment, M, segments, F, Beta, det_rat_sign, overlap, R);
            log_prob = log(Beta*(-t_down)/(segments.size()+1)*det_rat)-length*Eo+otherlength_u-nonloc;

            //
            if (log(rndm()) < log_prob)
            {
               // s is the segment which is shifted, r the segment which is inserted
               int s, r;
               s = 0;
               for (typename S::iterator it=segments.begin(); it!=s_down; it++) s++;
               if (anti_segment.t_end() > segment_shrink.t_start())
               {
                  r = s+1;
               }
               else
               {
                  r = 0;
                  s++;
               }

               //
               compute_M_insert_anti(anti_segment, s, r, M, segments, F, Beta, det_rat*overlap, R);

               //
               times segment_new_endpoint(*s_down);
               typename S::iterator prev_segment=s_down;
               if(s_down !=segments.begin())
               {
                  prev_segment--;
               }
               else
               {
                  prev_segment=segments.begin();
               }

               //
               segment_new_endpoint.set_t_end(t);
               segments.erase(s_down); //erase old segment (without shifted end
               s_down=segments.insert(segment_new_endpoint).first; //in

               //
               if(s_down==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
               typename S::iterator sit=segments.insert(segment_insert).first; //insert  new segment
               if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
               if (segment_insert.t_start()>segments.begin()->t_start())
               {
                  s_down++;
                  typename S::iterator sit=segments.insert(s_down, segment_insert);
                  if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
               }
               else
               {
                  typename S::iterator sit=segments.insert(segments.begin(), segment_insert);
                  if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
               }
            } // if (log(rndm()) < log_prob)
         } // if(t_down<0)
      } // if(full_line==1)
   } // if(rndm()<0.5)
   else if (segments.size()>1) // try to remove an anti-segment
   {
      //
      int r = (int)(rndm()*segments.size());
      s_up = segments.begin();
      for (int i=0; i<r; i++) s_up++;

      //
      int s = r-1;
      if (s<0)
      {
         s=segments.size()-1;
         s_down=segments.end();
         s_down--;
      }
      else
      {
         s_down=s_up;
         s_down--;
      }

      //
      double length = s_up->t_start() - s_down->t_end();
      if (length < 0) length += Beta;

      //
      double t_total = s_up->t_end() - s_down->t_end();
      if (t_total < 0) t_total += Beta;
      times anti_segment(s_down->t_end(),s_up->t_start());

      //
      double otherlength_u=0;
      int FLAVORS=other_full_line.size();
      for(int i=0; i<FLAVORS; i++)
      {
         if (i==this_flavor) continue;
         double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], Beta);
         otherlength_u += other_length*Umat(i, this_flavor);
      }

      //
      double nonloc=0.0;
      if(K_table!=NULL)
      {
         VecVecVec K = *K_table;
         nonloc=nonlocal(anti_segment.t_end(), anti_segment.t_start(), other_segments, Beta, this_flavor, K, 1);
      }

      //
      double log_prob, det_rat, det_rat_sign;
      det_rat = det_rat_remove_anti(anti_segment, r, s, M, segments, F, Beta, det_rat_sign);
      log_prob = log(Beta*t_total/segments.size()/det_rat)-length*Eo+otherlength_u-nonloc;

      //
      if(log(rndm()) < -log_prob)
      {
         compute_M_remove_anti(M, s, r);
         double t_end = s_up->t_end();
         segments.erase(s_up);

         //
         if (r>0)
         {
            s_up=segments.begin();
            for (int k=0; k<s; k++) s_up++;
         }
         else
         {
            s=segments.size()-1;
            s_up = segments.end();
            s_up--;
         }

         //
         times s_up_new(*s_up);
         s_up_new.set_t_end(t_end);
         segments.erase(s_up);
         typename S::iterator sit=segments.insert(s_up_new).first;
         if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
      } // if(log(rndm()) < -log_prob)
   } // else if (segments.size()>1)
   else if (segments.size()==1)
   {
      s_down = segments.begin();

      //
      double det_rat = fabs(M(0,0));
      double length = s_down->t_start()-s_down->t_end();
      if (length<0) length += Beta;
      times anti_segment(s_down->t_end(),s_down->t_start());

      //
      double otherlength_u=0;
      int FLAVORS=other_full_line.size();
      for(int i=0; i<FLAVORS; i++)
      {
         if (i==this_flavor) continue;
         double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], Beta);
         otherlength_u += other_length*Umat(i, this_flavor);
      }

      //
      double nonloc=0.0;
      if(K_table!=NULL)
      {
         VecVecVec K = *K_table;
         nonloc=nonlocal(anti_segment.t_end(), anti_segment.t_start(), other_segments, Beta, this_flavor, K, 1);
      }

      //
      double log_prob = log(Beta*Beta/det_rat)-length*Eo+otherlength_u-nonloc;
      if(log(rndm()) < -log_prob)
      {
         full_line=1;
         segments.erase(s_down);
         compute_M_down(0,M); // attention: M.clear() sets elements to zero
      }
   } // else if (segments.size()==1)
}


//------------------------------------------------------------------------------


template <class S, class G> void insert_remove_segment( double t, double Beta,
   double Eo, Mat &Umat, G &F, S &segments, Mat &M, double &sign,
   std::vector<S>& other_segments, std::vector<int> other_full_line,
   int this_flavor, VecVecVec *K_table)
{
   //
   double t_up;                 // distance to next segment up
   double t_down;               // distance to next segment down
   typename S::iterator s_up;   // iterator of the segment up
   typename S::iterator s_down; // iterator of the segment down

   //
   if (rndm()<0.5) // try to insert a segment
   {
      //
      compute_intervals(t, Beta, t_up, t_down, segments,s_up, s_down);

      //
      if (t_down>0) // t does not lie on a segment -> it's possible to insert a new one starting from t
      {
         //
         double length = compute_length(rndm(), t_up, 0);
         double t_final = t + length;
         times segment_insert;
         segment_insert.set_t_start(t);

         //
         if (t_final > Beta)
             segment_insert.set_t_end(t_final-Beta);
         else
             segment_insert.set_t_end(t_final);

         //
         double otherlength_u=0;
         int FLAVORS=other_full_line.size();
         for (int i=0; i<FLAVORS; i++)
         {
             if (i==this_flavor) continue;
             double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], Beta);
             otherlength_u += other_length*Umat(i, this_flavor);
         }

         //
         double nonloc=0.0;
         if(K_table!=NULL)
         {
            VecVecVec K = *K_table;
            nonloc=nonlocal(segment_insert.t_start(), segment_insert.t_end(), other_segments, Beta, this_flavor, K, 0);
         }

         //
         double log_prob, overlap, det_rat, det_rat_sign;
         std::vector<double> Fs(segments.size()), Fe(segments.size());
         det_rat = det_rat_up(segment_insert, M, segments, F, Fs, Fe, Beta, det_rat_sign, overlap);
         log_prob = log(Beta*t_up/(segments.size()+1)*det_rat)+Eo*length-otherlength_u-nonloc;

         //
         if (log(rndm()) < log_prob)
         {
            int position=0;
            for (segment_container_t::iterator it=segments.begin(); it!=s_up; it++) position++;
            compute_M_up(position, M, Fs, Fe, det_rat*overlap);
            sign *= det_rat_sign;
            segment_container_t::iterator sit=segments.insert(s_up, segment_insert);
            if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
         }
      } // if (t_down>0)
   } // if (rndm()<0.5)
   else if (segments.size()>0) // try to remove a segment
   {
      //
      int position = (int)(rndm()*segments.size());
      s_down = segments.begin();
      for (int i=0; i<position; i++) s_down++;
      s_up=s_down;
      s_up++;
      if (s_up==segments.end())s_up = segments.begin();

      //
      double length = s_down->t_end()-s_down->t_start();
      if (length < 0) length += Beta;

      //
      double t_total = s_up->t_start()-s_down->t_start();
      if (t_total <= 0) t_total += Beta;

      //
      times segment_remove = *s_down;

      //
      double otherlength_u=0;
      int FLAVORS=other_full_line.size();
      for (int i=0; i<FLAVORS; i++)
      {
         if (i==this_flavor) continue;
         double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], Beta);
         otherlength_u += other_length*Umat(i, this_flavor);
      }

      //
      double nonloc=0.0;
      if(K_table!=NULL)
      {
         VecVecVec K = *K_table;
         nonloc=nonlocal(segment_remove.t_start(), segment_remove.t_end(), other_segments, Beta, this_flavor, K, 1);
      }

      //
      double log_prob, det_rat, det_rat_sign;
      det_rat = det_rat_down(position, M, segments, det_rat_sign);
      log_prob = log(Beta*t_total/segments.size()/det_rat)+length*Eo-otherlength_u-nonloc;

      //
      if (log(rndm()) < -log_prob)
      {
         compute_M_down(position, M);
         sign *= det_rat_sign;
         segments.erase(s_down);
      }
   }
}


//------------------------------------------------------------------------------


template <class S, class G> void shift_segment( S &segments, double Beta, double Eo,
   Mat &Umat, G &F, Mat &M, double &sign, std::vector<S> &other_segments,
   std::vector<int> &other_full_line, int this_flavor, VecVecVec *K_table)
{
   //
   int size = segments.size();
   if (size < 1) return;
   int n = (int)(size*rndm());

   //
   typename S::iterator s, s_up;
   s=segments.begin();
   for (int i=0; i<n; i++) s++;
   s_up = s; s_up++;
   if (s_up == segments.end()) s_up = segments.begin();

   //
   double interval = s_up->t_start() - s->t_start();
   if (interval <= 0) interval += Beta;

   //
   double length = compute_length(rndm(), interval, 0);
   double length_old = s->t_end()-s->t_start();
   if (length_old<0) length_old += Beta;

   //
   double new_t_end = s->t_start() + length;
   if (new_t_end > Beta) new_t_end -= Beta;
   times segment_insert(s->t_start(), new_t_end);
   times segment_remove=*s;

   //
   double otherlength_u=0;
   int FLAVORS=other_full_line.size();
   for (int i=0; i<FLAVORS; i++)
   {
      if (i==this_flavor) continue;
      double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], Beta)-compute_overlap(segment_remove, other_segments[i], other_full_line[i], Beta);
      otherlength_u += other_length*Umat(i, this_flavor);
   }

   //
   double nonloc=0.0;
   if(K_table!=NULL)
   {
      VecVecVec K = *K_table;
      nonloc=nonlocal_shift(new_t_end, s->t_end(), other_segments, Beta, this_flavor, K);
   }

   //
   double det_rat, det_rat_sign, overlap;
   det_rat = det_rat_shift(segment_insert, n, M, segments, F, Beta, det_rat_sign, overlap);

   //
   if (log(rndm()) < log(det_rat)+(length-length_old)*Eo-otherlength_u-nonloc)
   {
      compute_M_shift(segment_insert, n, M, segments, F, Beta, det_rat*overlap);
      sign *= det_rat_sign;
      times s_new(*s);
      s_new.set_t_end(new_t_end);
      segments.erase(s);
      segment_container_t::iterator sit;
      sit=segments.insert(s_new).first;
      if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
   }
}

#endif
