/***************************************************************************
* Impurity solver based on a Z-expansion in the impurity-bath hybridization
* (c) 2005 Philipp Werner
***************************************************************************/

#ifndef ___IMP___
#define ___IMP___

#include <alps/scheduler/montecarlo.h>

#include <stack>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <alps/parameter.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>
#include <cmath>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <alps/osiris/std/list.h>

typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
typedef boost::variate_generator<alps::buffered_rng_base&, boost::uniform_real<double> > prng_t;
typedef dense_matrix blas_matrix;
typedef std::vector<double> vector_t;

class times
{
 public:
  times() {t_start_=0; t_end_=0;};
  times(double t_start, double t_end) {t_start_=t_start; t_end_=t_end;};
  double t_start() const {return t_start_;} // begin of segment
  double t_end() const {return t_end_;} // end of segment
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


typedef std::set<times> segment_container_t;
typedef std::vector<double> hybridization_t;
typedef std::vector<std::vector<double> > hybridization_container_t;
typedef std::vector<std::vector<std::vector<double> > > k_table_t;

class IsingSimulation : public alps::scheduler::MCRun
{
public:
  IsingSimulation(const alps::ProcessList&,const alps::Parameters&,int);
  IsingSimulation(int,const alps::ProcessList&,alps::IDump&,int);
  void dostep();
  bool is_thermalized() const;
  double work_done() const;
  bool change_parameter(const std::string& name, const alps::StringValue& value);
  
private:
  int                               sweeps;					// sweeps done
  int                               thermalization_sweeps;	// sweeps to be done for equilibration
  int                               total_sweeps;			// sweeps to be done after equilibration
  double                            mu;                     // chemical potential
  std::vector<double>               mu_e;                   // mu-<\epsilon>
  dense_matrix                      u;                      // U matrix
  double                            t;                      // bandwidth=4t (for semi-circle)
  hybridization_container_t         F;                      // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu)
  std::vector<segment_container_t >	segments;               // stores configurations with 0,1,... segments (but not full line)
  std::vector<int>					full_line;              // if 1 means that particle occupies full time-line
  std::vector<double>				sign;                   // sign of Z_n_up [actually not needed for density-density]
  std::vector<dense_matrix>			M;                      // inverse hybridization matrix
  std::valarray<double> G_meas;                             // measured GF
  k_table_t K_table;                                     // K function matrix for retarded interactions

};

typedef alps::scheduler::SimpleMCFactory<IsingSimulation> IsingFactory;

#endif









