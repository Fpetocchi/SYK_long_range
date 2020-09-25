#ifndef ___TIMES___
#define ___TIMES___

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
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>


//============================================================================//


std::default_random_engine generator;

double rndm()
{
   std::uniform_real_distribution<double> distribution_flat(0.0,1.0);
   double r = distribution_flat(generator);
   return r;
}


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


#endif
