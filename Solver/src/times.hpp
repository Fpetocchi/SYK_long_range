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

double flat_distribution()
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
typedef Eigen::MatrixXd dense_matrix;
typedef Eigen::VectorXd vector_t;                                               //typedef std::vector<double> vector_t;
typedef std::set<times> segment_container_t;
typedef Eigen::VectorXd hybridization_t;                                        //all the orbitals(diagonal) at given itau
typedef std::vector<Eigen::VectorXd> hybridization_container_t;                 //all the orbitals(diagonal) for a list of Ntau
typedef Eigen::MatrixXd Kfunct_t;                                               //all the orbitals(Norb**2 x Norb**2) at given itau
typedef std::vector<Eigen::MatrixXd> Kfunct_container_t;                        //all the orbitals(Norb**2 x Norb**2) for a list of Ntau


//------------------------------------------------------------------------------


#endif
