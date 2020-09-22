// -----------------------------------------------------------------------------

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

// -----------------------------------------------------------------------------

#ifndef CNTR_USE_OMP
#define CNTR_USE_OMP
#endif

#ifndef CNTR_USE_MPI
#define CNTR_USE_MPI
#endif

#define NCA_SOLVER_ASSERT_0 0
#define NCA_SOLVER_ASSERT_1 0

#include <cntr/cntr.hpp>
#include <cntr/utils/read_inputfile.hpp>

// -----------------------------------------------------------------------------

#include "./ppsc/ppsc.hpp"
#include "./ppsc/solver.hpp"

#include "./ppsc/hilbert_spaces/three_band_fermi_diag.hpp"
#include "./ppsc/hilbert_spaces/three_band_fermi_densdens.hpp"
#include "./ppsc/hamiltonians/three_band_hubbard_solar.hpp"

#include "./ppsc/baths/non_int_boson_propagator.hpp"

// -----------------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------------

typedef ppsc::operator_type operator_type;
typedef ppsc::mam::dynamic_matrix_type matrix_type;


////////////////////////////////////////////////////////////////////////////////


class single_particle_greensfunction_type
{
  public:
  single_particle_greensfunction_type(int nt, int ntau) : nt(nt), ntau(ntau),
    o1up(nt, ntau, 1, -1), o1do(nt, ntau, 1, -1),
    o2up(nt, ntau, 1, -1), o2do(nt, ntau, 1, -1),
    o3up(nt, ntau, 1, -1), o3do(nt, ntau, 1, -1)
    {}

  // ---------------------------------------------------------------------------

  void update(int tstp, ppsc::gf_tstps_type & gf_tstps, double linear_mixing=0.0)
  {
    if(tstp == -1) // MATSUBARA SYMMETRIES
    {
      ppsc::gf_tstp_type gloc_old(-1, ntau, 1);
      ppsc::gf_tstp_type gloc_mix(-1, ntau, 1);

      gloc_mix.clear();
      this->o1up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[0], 1.0 - linear_mixing);
      this->o1up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o1do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[1], 1.0 - linear_mixing);
      this->o1do.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o2up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[2], 1.0 - linear_mixing);
      this->o2up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o2do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[3], 1.0 - linear_mixing);
      this->o2do.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o3up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[4], 1.0 - linear_mixing);
      this->o3up.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o3do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, linear_mixing);
      gloc_mix.incr(gf_tstps[5], 1.0 - linear_mixing);
      this->o3do.set_timestep(-1, gloc_mix);

      // SPIN SYMMETRY

      gloc_mix.clear();
      this->o1up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o1do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o1up.set_timestep(-1, gloc_mix);
      this->o1do.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o2up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o2do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o2up.set_timestep(-1, gloc_mix);
      this->o2do.set_timestep(-1, gloc_mix);

      gloc_mix.clear();
      this->o3up.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o3do.get_timestep(-1, gloc_old);
      gloc_mix.incr(gloc_old, 0.5);
      this->o3up.set_timestep(-1, gloc_mix);
      this->o3do.set_timestep(-1, gloc_mix);
    }
    else // TSTP UPDATE WITH NO IMPOSED SYMMETRIES
    {
      this->o1up.set_timestep(tstp, gf_tstps[0]);
      this->o1do.set_timestep(tstp, gf_tstps[1]);
      this->o2up.set_timestep(tstp, gf_tstps[2]);
      this->o2do.set_timestep(tstp, gf_tstps[3]);
      this->o3up.set_timestep(tstp, gf_tstps[4]);
      this->o3do.set_timestep(tstp, gf_tstps[5]);
    }
  }

  // ---------------------------------------------------------------------------

  void store(hid_t file_id)
  {
    hid_t group_id;
    group_id = create_group(file_id, "g1u");
    store_herm_greens_function(group_id, this->o1up);
    close_group(group_id);
    group_id = create_group(file_id, "g1d");
    store_herm_greens_function(group_id, this->o1do);
    close_group(group_id);
    group_id = create_group(file_id, "g2u");
    store_herm_greens_function(group_id, this->o2up);
    close_group(group_id);
    group_id = create_group(file_id, "g2d");
    store_herm_greens_function(group_id, this->o2do);
    close_group(group_id);
    group_id = create_group(file_id, "g3u");
    store_herm_greens_function(group_id, this->o3up);
    close_group(group_id);
    group_id = create_group(file_id, "g3d");
    store_herm_greens_function(group_id, this->o3do);
    close_group(group_id);
  }

  // ---------------------------------------------------------------------------

  void load(std::string filename)
  {
    this->o1up.read_from_hdf5(filename.c_str(), "g1u");
    this->o1do.read_from_hdf5(filename.c_str(), "g1d");
    this->o2up.read_from_hdf5(filename.c_str(), "g2u");
    this->o2do.read_from_hdf5(filename.c_str(), "g2d");
    this->o3up.read_from_hdf5(filename.c_str(), "g3u");
    this->o3do.read_from_hdf5(filename.c_str(), "g3d");
  }

  // ---------------------------------------------------------------------------

  int nt, ntau;
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do, o3up, o3do;

};


////////////////////////////////////////////////////////////////////////////////


class hybridization_function_type
{
  public:
  hybridization_function_type(int nt, int ntau) :
    nt(nt), ntau(ntau),
    o1up(nt, ntau, 1, -1),
    o1do(nt, ntau, 1, -1),
    o2up(nt, ntau, 1, -1),
    o2do(nt, ntau, 1, -1),
    o3up(nt, ntau, 1, -1),
    o3do(nt, ntau, 1, -1),
    // -- bwd hybridizations
    o1up_cc(nt, ntau, 1, -1),
    o1do_cc(nt, ntau, 1, -1),
    o2up_cc(nt, ntau, 1, -1),
    o2do_cc(nt, ntau, 1, -1),
    o3up_cc(nt, ntau, 1, -1),
    o3do_cc(nt, ntau, 1, -1)
  {}

  // ---------------------------------------------------------------------------

  void _component_update(int tstp,
                         cntr::herm_matrix<double> & Delta_comp,
                         cntr::herm_matrix<double> & Delta_comp_cc,
                         cntr::herm_matrix<double> & Gaa,
                         cntr::herm_matrix<double> & Gab1,
                         cntr::herm_matrix<double> & Gab2,
                         // Bethe hopping terms
                         cntr::function<double> & taa,  // same orbital
                         cntr::function<double> & tab1, // different orbital
                         cntr::function<double> & tab2, // different orbital
                         // boson bath
                         cntr::herm_matrix<double> & D0, int bath_flag=0)
  {
    // first interband term
    cntr::herm_matrix_timestep<double> S1, tmp1;
    Gab1.get_timestep(tstp, tmp1);
    Delta_comp.set_timestep(tstp, tmp1);
    Delta_comp.left_multiply(tstp, tab1);
    Delta_comp.right_multiply(tstp, tab1);
    Delta_comp.get_timestep(tstp, S1);

    // second interband term
    cntr::herm_matrix_timestep<double> S2, tmp2;
    Gab2.get_timestep(tstp, tmp2);
    Delta_comp.set_timestep(tstp, tmp2);
    Delta_comp.left_multiply(tstp, tab2);
    Delta_comp.right_multiply(tstp, tab2);
    Delta_comp.get_timestep(tstp, S2);

    // intraband term
    cntr::herm_matrix_timestep<double> S, tmp;
    Gaa.get_timestep(tstp, tmp);
    Delta_comp.set_timestep(tstp, tmp);
    Delta_comp.left_multiply(tstp, taa);
    Delta_comp.right_multiply(tstp, taa);
    Delta_comp.get_timestep(tstp, S);

    // add intraband
    Delta_comp.incr_timestep(tstp, S1);
    Delta_comp.incr_timestep(tstp, S2);

    // add boson bath
    if(bath_flag)
    {
      Bubble2(tstp, tmp, Gaa, D0);
      Delta_comp.incr_timestep(tstp, tmp);
    }

    ppsc::set_bwd_from_fwd(tstp, Delta_comp_cc, Delta_comp);

  }

  // ---------------------------------------------------------------------------

  void update(int tstp, single_particle_greensfunction_type & Gloc,
                        // Bethe hopping matrix
                        cntr::function<double> & tmatrix,
                        // boson bath
                        cntr::herm_matrix<double> & D0,
                        int af_bethe_flag=0, int bath_flag=0)
  {
    cntr::function<double> taa(nt);
    cntr::function<double> tab1(nt);
    cntr::function<double> tab2(nt);

    // -- Update Delta1
    tmatrix.get_matrixelement(0, 0, taa);
    tmatrix.get_matrixelement(0, 1, tab1);
    tmatrix.get_matrixelement(0, 2, tab2);
    //up
    if(af_bethe_flag == 1) _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1do, Gloc.o2do, Gloc.o3do, taa, tab1, tab2, D0, bath_flag);
    else                   _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1up, Gloc.o2up, Gloc.o3up, taa, tab1, tab2, D0, bath_flag);
    //do
    if(af_bethe_flag == 1) _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1up, Gloc.o2up, Gloc.o3up, taa, tab1, tab2, D0, bath_flag);
    else                   _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1do, Gloc.o2do, Gloc.o3do, taa, tab1, tab2, D0, bath_flag);

    // -- Update Delta2
    tmatrix.get_matrixelement(1, 1, taa);
    tmatrix.get_matrixelement(1, 0, tab1);
    tmatrix.get_matrixelement(1, 2, tab2);
    //up
    if(af_bethe_flag == 1) _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o2do, Gloc.o1do, Gloc.o3do, taa, tab1, tab2, D0, bath_flag);
    else                   _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o2up, Gloc.o1up, Gloc.o3up, taa, tab1, tab2, D0, bath_flag);
    //do
    if(af_bethe_flag == 1) _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o2up, Gloc.o1up, Gloc.o3up, taa, tab1, tab2, D0, bath_flag);
    else                   _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o2do, Gloc.o1do, Gloc.o3do, taa, tab1, tab2, D0, bath_flag);

    // -- Update Delta3
    tmatrix.get_matrixelement(2, 2, taa);
    tmatrix.get_matrixelement(2, 0, tab1);
    tmatrix.get_matrixelement(2, 1, tab2);
    //up
    if(af_bethe_flag == 1) _component_update(tstp, this->o3up, this->o3up_cc, Gloc.o3do, Gloc.o1do, Gloc.o2do, taa, tab1, tab2, D0, bath_flag);
    else                   _component_update(tstp, this->o3up, this->o3up_cc, Gloc.o3up, Gloc.o1up, Gloc.o2up, taa, tab1, tab2, D0, bath_flag);
    //do
    if(af_bethe_flag == 1) _component_update(tstp, this->o3do, this->o3do_cc, Gloc.o3up, Gloc.o1up, Gloc.o2up, taa, tab1, tab2, D0, bath_flag);
    else                   _component_update(tstp, this->o3do, this->o3do_cc, Gloc.o3do, Gloc.o1do, Gloc.o2do, taa, tab1, tab2, D0, bath_flag);
  }

  // ---------------------------------------------------------------------------

  void store(hid_t file_id)
  {
    hid_t group_id;
    group_id = create_group(file_id, "d1u");
    store_herm_greens_function(group_id, this->o1up);
    close_group(group_id);
    group_id = create_group(file_id, "d1d");
    store_herm_greens_function(group_id, this->o1do);
    close_group(group_id);
    group_id = create_group(file_id, "d2u");
    store_herm_greens_function(group_id, this->o2up);
    close_group(group_id);
    group_id = create_group(file_id, "d2d");
    store_herm_greens_function(group_id, this->o2do);
    close_group(group_id);
    group_id = create_group(file_id, "d3u");
    store_herm_greens_function(group_id, this->o3up);
    close_group(group_id);
    group_id = create_group(file_id, "d3d");
    store_herm_greens_function(group_id, this->o3do);
    close_group(group_id);
    group_id = create_group(file_id, "d1ucc");
    store_herm_greens_function(group_id, this->o1up_cc);
    close_group(group_id);
    group_id = create_group(file_id, "d1dcc");
    store_herm_greens_function(group_id, this->o1do_cc);
    close_group(group_id);
    group_id = create_group(file_id, "d2ucc");
    store_herm_greens_function(group_id, this->o2up_cc);
    close_group(group_id);
    group_id = create_group(file_id, "d2dcc");
    store_herm_greens_function(group_id, this->o2do_cc);
    close_group(group_id);
    group_id = create_group(file_id, "d3ucc");
    store_herm_greens_function(group_id, this->o3up_cc);
    close_group(group_id);
    group_id = create_group(file_id, "d3dcc");
    store_herm_greens_function(group_id, this->o3do_cc);
    close_group(group_id);
  }

  // ---------------------------------------------------------------------------

  void load(std::string filename)
  {
    this->o1up.read_from_hdf5(filename.c_str(), "d1u");
    this->o1do.read_from_hdf5(filename.c_str(), "d1d");
    this->o2up.read_from_hdf5(filename.c_str(), "d2u");
    this->o2do.read_from_hdf5(filename.c_str(), "d2d");
    this->o3up.read_from_hdf5(filename.c_str(), "d3u");
    this->o3do.read_from_hdf5(filename.c_str(), "d3d");

    this->o1up_cc.read_from_hdf5(filename.c_str(), "d1ucc");
    this->o1do_cc.read_from_hdf5(filename.c_str(), "d1dcc");
    this->o2up_cc.read_from_hdf5(filename.c_str(), "d2ucc");
    this->o2do_cc.read_from_hdf5(filename.c_str(), "d2dcc");
    this->o3up_cc.read_from_hdf5(filename.c_str(), "d3ucc");
    this->o3do_cc.read_from_hdf5(filename.c_str(), "d3dcc");
  }

  // ---------------------------------------------------------------------------

  int nt, ntau;
  cntr::herm_matrix<double> o1up, o1do, o2up, o2do, o3up, o3do;
  cntr::herm_matrix<double> o1up_cc, o1do_cc, o2up_cc, o2do_cc, o3up_cc, o3do_cc;

};


////////////////////////////////////////////////////////////////////////////////


template<class HAM>
ppsc::pp_ints_type get_pp_ints(hybridization_function_type & Delta, HAM & h)
{
  // spin (u)p/(d)own and (c)reation/(a)nihilation operators
  int boson=+1, fermion=-1, fwd=+1, bwd=-1;
  ppsc::pp_ints_type pp_ints;
  // orbital no. 1
  pp_ints.push_back(ppsc::pp_int_type(Delta.o1up,    h.c1uc,  h.c1ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o1up_cc, h.c1ua,  h.c1uc,  fermion, bwd)); // spin up bwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o1do,    h.c1dc,  h.c1da,  fermion, fwd)); // spin do fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o1do_cc, h.c1da,  h.c1dc,  fermion, bwd)); // spin do bwd
  // orbital no. 2
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2up,    h.c2uc,  h.c2ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2up_cc, h.c2ua,  h.c2uc,  fermion, bwd)); // spin up bwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2do,    h.c2dc,  h.c2da,  fermion, fwd)); // spin do fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o2do_cc, h.c2da,  h.c2dc,  fermion, bwd)); // spin do bwd
  // orbital no. 3
  pp_ints.push_back(ppsc::pp_int_type(Delta.o3up,    h.c3uc,  h.c3ua,  fermion, fwd)); // spin up fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o3up_cc, h.c3ua,  h.c3uc,  fermion, bwd)); // spin up bwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o3do,    h.c3dc,  h.c3da,  fermion, fwd)); // spin do fwd
  pp_ints.push_back(ppsc::pp_int_type(Delta.o3do_cc, h.c3da,  h.c3dc,  fermion, bwd)); // spin do bwd
  return pp_ints;
}

template<class HAM>
ppsc::gf_verts_type get_gf_verts(HAM & h)
{
  ppsc::gf_verts_type gf_verts;
  // orbital no. 1
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, h.c1ua, h.c1uc)); // spin up orb 1
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, h.c1da, h.c1dc)); // spin do orb 1
  // orbital no. 2
  gf_verts.push_back(ppsc::gf_vert_type(1, 1, h.c2ua, h.c2uc)); // spin up orb 2
  gf_verts.push_back(ppsc::gf_vert_type(1, 1, h.c2da, h.c2dc)); // spin do orb 2
  // orbital no. 3
  gf_verts.push_back(ppsc::gf_vert_type(2, 2, h.c3ua, h.c3uc)); // spin up orb 3
  gf_verts.push_back(ppsc::gf_vert_type(2, 2, h.c3da, h.c3dc)); // spin do orb 3
  return gf_verts;
}


////////////////////////////////////////////////////////////////////////////////


int main(int argc, char *argv[])
{
  int iter_equil, iter_warmup, itermax,
      iter_rt, iter_rtime,
      nt, ntau, kt, tstp, order, nomp, store_pp;
  int read_eq_sym, read_rt_sym, read_state_from_file, af_bethe_flag;
  int bath_flag, sym_flag;
  double beta, h, errmax, linear_mixing, dmfterr, dmfterr_equil, sym_tol;
  bool matsubara_converged = false;
  dmatrix thop(3,3);

  try
  {

    // Read general input - no microscopic parameters---------------------------
    {
      if (argc < 2) throw("COMMAND LINE ARGUMENT MISSING");
      find_param(argv[1], "__nt=", nt);
      find_param(argv[1], "__ntau=", ntau);
      find_param(argv[1], "__beta=", beta);
      find_param(argv[1], "__h=", h);
      find_param(argv[1], "__itermax=", itermax);
      find_param(argv[1], "__errmax=", errmax);
      find_param(argv[1], "__iter_rtime=", iter_rtime);
      find_param(argv[1], "__kt=", kt);
      find_param(argv[1], "__order=", order);
      find_param(argv[1], "__store_pp=", store_pp);
      find_param(argv[1], "__linear_mixing=", linear_mixing);
      find_param(argv[1], "__sym_flag=", sym_flag);
      find_param(argv[1], "__sym_tol=", sym_tol);
      find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__read_state_from_file=", read_state_from_file);
      find_param(argv[1], "__af_bethe_flag=", af_bethe_flag);
      find_param(argv[1], "__bath_flag=", bath_flag);
      find_param(argv[1], "__nomp=", nomp);
    }


    // Setup boson bath---------------------------------------------------------
    cntr::herm_matrix<double> D0(nt, ntau, 1, +1);
    cntr::herm_matrix<double> D0_cc(nt, ntau, 1, +1);
    double omega;
    std::vector<double> g;

    find_param_tvector(argv[1], "__g=", g, nt);
    find_param(argv[1], "__omega=", omega);

    ppsc::boson_utils::green_from_eps_phonon(beta, D0, omega, h);
    cntr::function<double> gfunc(nt);

    for(int tstp=-1;tstp<=nt;tstp++)
    {
      cdmatrix tmp(1,1);
      tmp(0,0)=g[tstp+1];
      gfunc.set_value(tstp,tmp);
    }
    for(int tstp=-1; tstp <= nt; tstp++)
    {
      D0.left_multiply(tstp, gfunc);
      D0.right_multiply(tstp, gfunc);
    }
    for(int tstp = -1; tstp <= nt; tstp++)ppsc::set_bwd_from_fwd(tstp, D0_cc, D0);


    //Setup pp calculator-------------------------------------------------------
    typedef ppsc::hilbert_spaces::three_band_fermi_diag hilbert_space_type;
    typedef ppsc::hamiltonians::three_band_hubbard_solar<
      hilbert_space_type,
      ppsc::hamiltonians::interaction_type::kanamori> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;

    hilbert_space_type hilbert_space;
    hilbert_space.init();

    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);


    // local Hamiltonian microscopic parameters---------------------------------
    {
      find_param(argv[1]        , "__mu="       , imp.hamiltonian.mu        );
      find_param(argv[1]        , "__U0="       , imp.hamiltonian.U[0]      );
      find_param_tvector(argv[1], "__U="        , imp.hamiltonian.U    , nt );
      find_param_tvector(argv[1], "__dU1="      , imp.hamiltonian.dU1  , nt );
      find_param_tvector(argv[1], "__dU2="      , imp.hamiltonian.dU2  , nt );
      find_param_tvector(argv[1], "__dU3="      , imp.hamiltonian.dU3  , nt );
      find_param(argv[1]        , "__J0="       , imp.hamiltonian.J[0]      );
      find_param_tvector(argv[1], "__J="        , imp.hamiltonian.J    , nt );
      find_param(argv[1]        , "__delta0="   , imp.hamiltonian.delta[0]  );
      find_param_tvector(argv[1], "__delta="    , imp.hamiltonian.delta, nt );
      find_param(argv[1]        , "__E10="      , imp.hamiltonian.E1[0]     );
      find_param(argv[1]        , "__E20="      , imp.hamiltonian.E2[0]     );
      find_param(argv[1]        , "__E30="      , imp.hamiltonian.E3[0]     );
      find_param_tvector(argv[1], "__E1="       , imp.hamiltonian.E1   , nt );
      find_param_tvector(argv[1], "__E2="       , imp.hamiltonian.E2   , nt );
      find_param_tvector(argv[1], "__E3="       , imp.hamiltonian.E3   , nt );
      find_param(argv[1]        , "__Bz0="      , imp.hamiltonian.Bz[0]     );
      find_param_tvector(argv[1], "__Bz="       , imp.hamiltonian.Bz   , nt );
      //
      find_param(argv[1]        , "__t11="      , thop(0,0) );
      find_param(argv[1]        , "__t12="      , thop(0,1) );
      find_param(argv[1]        , "__t13="      , thop(0,2) );
      find_param(argv[1]        , "__t21="      , thop(1,0) );
      find_param(argv[1]        , "__t22="      , thop(1,1) );
      find_param(argv[1]        , "__t23="      , thop(1,2) );
      find_param(argv[1]        , "__t31="      , thop(2,0) );
      find_param(argv[1]        , "__t32="      , thop(2,1) );
      find_param(argv[1]        , "__t33="      , thop(2,2) );
      /*
      // find_param_tvector(argv[1], "__phi1="     , phi_vec1, nt);
      // find_param_tvector(argv[1], "__phi2="     , phi_vec2, nt);
      // find_param_tvector(argv[1], "__phi3="     , phi_vec3, nt);
      */
      imp.update_hamiltonian();

      // Hamiltonian Hermicity check
      ppsc::operator_type Htemp;
      imp.hamiltonian.get_hamiltonian(-1, Htemp);
      ppsc::mam::dynamic_matrix_type H = Htemp.to_dense();
      float herm_diff = (H - H.transpose()).cwiseAbs().maxCoeff();
      std::cout << "check hermiticity " << herm_diff << "\n";
    }


    // init Hoppings------------------------------------------------------------
    cntr::function<double> tmatrix(nt,3);
    tmatrix.set_constant(thop);

    // init DMFT fields---------------------------------------------------------
    single_particle_greensfunction_type Gloc(nt, ntau);
    hybridization_function_type Delta(nt, ntau);

    if(read_state_from_file)
    {
      std::string filename = "data_ppsc.h5";
      std::cout << "--> Reading state from file: " << filename << std::endl;

      Gloc.load(filename);
      Delta.load(filename);
      imp.load(filename);

      if(sym_flag && read_eq_sym)
      {
        std::string filename = "sym_eq.txt";
        imp.read_symmetries(filename);
      }

      imp.update_density_matrix(-1);
      imp.hamiltonian.update_exp_vals(-1, imp.rho);
      std::cout << "pp_mu = " << imp.pp_mu << std::endl;
    }

    // MATSUBARA PART - EQUILIBRIUM INITIAL STATE-------------------------------
    {
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);
      if(!read_state_from_file)
      {
        imp.solve_atomic();
      } else
      {
        gtmp.set_timestep(-1, Gloc.o1up);
      }

      // DMFT loops
      for (iter_equil = 1; iter_equil <= itermax; iter_equil++)
      {
        // -- Construct interactions and verticies
        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);
        imp.update_diagrams(pp_ints, gf_verts);
        if(iter_equil == 1 && read_eq_sym && !imp.has_symmetries())
        {
          std::string filename = "sym_eq.txt";
          imp.read_symmetries(filename);
        }

        // -- Solve pseudo particle problem
        imp.pp_step(-1);

        // -- Get spgf and mix
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);
        Gloc.update(-1, gf_tstps, linear_mixing);

        // -- Check error
        dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc.o1up);
        gtmp.set_timestep(-1, Gloc.o1up);

        // -- Update Hybridization
        Delta.update( -1, Gloc, tmatrix, D0, af_bethe_flag, bath_flag);

        cout << "Eq. iter:  " << iter_equil<< " err: " << dmfterr_equil << endl;

        // -- Convergence check
        if(dmfterr_equil < errmax)
        {
          matsubara_converged = true;
          break;
        }
      } // END DMFT loops
      if (iter_equil > itermax)
      {
        cerr << "WARNING: Matsubara not converged  after " << itermax
             << "steps ... abort" << endl;
        cerr << "skip real-time calculation " << endl;
      }
    } // END MATSUBARA


    // save symmetries----------------------------------------------------------
    if(sym_flag && !imp.has_symmetries())
    {
      imp.symmetry_reduction(-1);
      std::string filename = "sym_eq.txt";
      imp.write_symmetries(filename);
    }
    if(sym_flag && nt > 0) imp.clear_symmetries();
    if(sym_flag && read_rt_sym)
    {
      std::string filename = "sym_rt.txt";
      imp.read_symmetries(filename);
    }


    // REAL TIME PART - BOOTSTRAP-----------------------------------------------
    if (nt > 0 && matsubara_converged == true)
    {
      matsubara_converged = false;
      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();

      // DMFT loops
      for (iter_warmup = 1; iter_warmup <= itermax; iter_warmup++)
      {
        // -- Construct interactions and verticies
        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);
        imp.update_diagrams(pp_ints, gf_verts);

        // -- Solve pseudo particle problem
        imp.pp_step(kt);

        // -- Get spgf & Update Hybridization and mix
        for (int n = 0; n <= kt; n++)
        {
          ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
	        Gloc.update(n, gf_tstps);
          Delta.update( n, Gloc, tmatrix, D0, af_bethe_flag, bath_flag);
        }

        // -- Check error
        dmfterr = cntr::distance_norm2(kt, gtmp, Gloc.o1up);
        gtmp.set_timestep(kt, Gloc.o1up);

        cout << "WARMUP: iter:  " << iter_warmup << " err: " << dmfterr << endl;
        if (dmfterr < errmax)
        {
          matsubara_converged = true;
          break;
        }
      } // END DMFT loops
    } // END BOOTSTRAP


    // Reduce number of diagrams using symmetries-------------------------------
    if(sym_flag && matsubara_converged && nt > 0 && !imp.has_symmetries())
    {
      imp.symmetry_reduction(kt);
      std::string filename = "sym_rt.txt";
      imp.write_symmetries(filename);
    }


    // REAL TIME PART - NON EQUILIBRIUM STATE-----------------------------------
    if(matsubara_converged)

    // Real time steps
    for (tstp = kt + 1; tstp <= nt; tstp++)
    {
      imp.extrapolate_timestep(tstp - 1);

      // DMFT loops
      for (iter_rt = 1; iter_rt <= iter_rtime; iter_rt++)
      {
        // -- Construct interactions and verticies
        ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian);
        imp.update_diagrams(pp_ints, gf_verts);

        // -- Solve pseudo particle problem
        imp.pp_step(tstp);

        // -- Get spgf and mix
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
        Gloc.update(tstp, gf_tstps);

        // -- Update Hybridization
        Delta.update( tstp, Gloc, tmatrix, D0, af_bethe_flag, bath_flag);
        //Delta.update( tstp, Gloc, tfunc1sin, tfunc2sin, tfunc3sin
        //                  , tfunc1cos, tfunc2cos, tfunc3cos
        //                  , D0, af_bethe_flag, bath_flag);
      } // END DMFT loops
    } // END time steps


    // KINETIK ENERGY CALCULATION-----------------------------------------------
    std::vector<double> Ekin_o1up = ppsc::get_kinetic_energy(Gloc.o1up, Delta.o1up, beta, h, kt);
    std::vector<double> Ekin_o1do = ppsc::get_kinetic_energy(Gloc.o1do, Delta.o1do, beta, h, kt);
    std::vector<double> Ekin_o2up = ppsc::get_kinetic_energy(Gloc.o2up, Delta.o2up, beta, h, kt);
    std::vector<double> Ekin_o2do = ppsc::get_kinetic_energy(Gloc.o2do, Delta.o2do, beta, h, kt);
    std::vector<double> Ekin_o3up = ppsc::get_kinetic_energy(Gloc.o3up, Delta.o3up, beta, h, kt);
    std::vector<double> Ekin_o3do = ppsc::get_kinetic_energy(Gloc.o3do, Delta.o3do, beta, h, kt);

    std::vector<double> Ekin;
    Ekin.resize(Ekin_o1up.size());

    for( auto t : ppsc::range(0, Ekin.size()) )
    {
      Ekin[t] = Ekin_o1up[t] + Ekin_o1do[t]
              + Ekin_o2up[t] + Ekin_o2do[t]
              + Ekin_o3up[t] + Ekin_o3do[t];
      std::cout << "ekin, epot " << t << " " << Ekin[t] << " "
      << imp.hamiltonian.Eint_exp[t] << "\n";
    }


    // OBSERVABLES STORE--------------------------------------------------------
    {
      Gloc.o1up.print_to_file("Gloc.o1up.out");
      Gloc.o2up.print_to_file("Gloc.o2up.out");
      Gloc.o3up.print_to_file("Gloc.o3up.out");
      Gloc.o1up.print_to_file("Gloc.o1do.out");
      Gloc.o2do.print_to_file("Gloc.o2do.out");
      Gloc.o3do.print_to_file("Gloc.o3do.out");

      std:string filename = "data_ppsc.h5";
      hid_t file_id = open_hdf5_file(filename);
      hid_t group_id;

      imp.store(file_id, store_pp);
      Gloc.store(file_id);
      Delta.store(file_id);

      // -- Bethe general properties
      group_id = create_group(file_id, "bethe");

      // store_real_data_to_hid(group_id, "t_vec", t_vec.data(), t_vec.size());
      store_real_data_to_hid(group_id, "Ekin", Ekin.data(), Ekin.size());
      store_real_data_to_hid(group_id, "Ekin_o1up", Ekin_o1up.data(), Ekin_o1up.size());
      store_real_data_to_hid(group_id, "Ekin_o1do", Ekin_o1do.data(), Ekin_o1do.size());
      store_real_data_to_hid(group_id, "Ekin_o2up", Ekin_o2up.data(), Ekin_o2up.size());
      store_real_data_to_hid(group_id, "Ekin_o2do", Ekin_o2do.data(), Ekin_o2do.size());
      store_real_data_to_hid(group_id, "Ekin_o3up", Ekin_o3up.data(), Ekin_o3up.size());
      store_real_data_to_hid(group_id, "Ekin_o3do", Ekin_o3do.data(), Ekin_o3do.size());
      store_double_attribute_to_hid(group_id, "dmfterr_equil", dmfterr_equil);
      store_int_attribute_to_hid(group_id, "iter_equil", iter_equil);
      close_group(group_id);
      close_hdf5_file(file_id);
    }

  } // END of global try
  catch (char *message)
  {
    cerr << "exception\n**** " << message << " ****" << endl;
    cerr << "CDMFT input_file [ --test ]\n" << endl;
  } catch (...)
  {
    cerr << "unspecified exception " << endl;
    cerr << "\nCDMFT input_file [ --test ]\n" << endl;
  }
  return 0;
}





/*
void _component_update(int tstp,
                       cntr::herm_matrix<double> & Delta_comp,
                       cntr::herm_matrix<double> & Delta_comp_cc,
                       cntr::herm_matrix<double> & Gloc_comp,
                       cntr::function<double> & tsin,
                       cntr::function<double> & tcos,
                       // boson bath
                       cntr::herm_matrix<double> & D0, int bath_flag=0)
  {
    // Get Gloc-->tmpGloc
    cntr::herm_matrix_timestep<double> tmpGloc;
    Gloc_comp.get_timestep(tstp, tmpGloc);
    // Compute Delta1 = sin(t)G(t,t')sin(t')
    Delta_comp.set_timestep(tstp, tmp);
    Delta_comp.left_multiply(tstp, tsin);
    Delta_comp.right_multiply(tstp, tsin);
    // Delta1-->tmpsin
    cntr::herm_matrix_timestep<double> tmpsin;
    Delta_comp.get_timestep(tstp, tmpsin);
    // Compute Delta2 = cos(t)G(t,t')cos(t')
    Delta_comp.set_timestep(tstp, tmpGloc);
    Delta_comp.left_multiply(tstp, tcos);
    Delta_comp.right_multiply(tstp, tcos);
    // Delta = Delta2 + tmpsin
    Delta_comp.incr_timestep(tstp, tmpsin);

    // add boson bath
    if(bath_flag)
    {
      Bubble2(tstp, tmp, Gloc_comp, D0);
      Delta_comp.incr_timestep(tstp, tmp);
    }
    ppsc::set_bwd_from_fwd(tstp, Delta_comp_cc, Delta_comp);
  }

void update(int tstp,
            single_particle_greensfunction_type & Gloc,
            // Bethe hoppings
            cntr::function<double> & tfunc1sin,
            cntr::function<double> & tfunc2sin,
            cntr::function<double> & tfunc3sin,
            cntr::function<double> & tfunc1cos,
            cntr::function<double> & tfunc2cos,
            cntr::function<double> & tfunc3cos,
            // boson bath
            cntr::herm_matrix<double> & D0,
            int af_bethe_flag=0, int bath_flag=0)
{

  // -- Update Delta.o1up
  if(af_bethe_flag == 1) _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);
  else                   _component_update(tstp, this->o1up, this->o1up_cc, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);

  // -- Update Delta.o1do
  if(af_bethe_flag == 1) _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1up, tfunc1sin, tfunc1cos, D0, bath_flag);
  else                   _component_update(tstp, this->o1do, this->o1do_cc, Gloc.o1do, tfunc1sin, tfunc1cos, D0, bath_flag);

  // for solar cell: modify the hopping for orbitals 2,3 to off-diagonal 2 <-> 3 hopping. This implies \Delta_2=v*G_3*v and \Delta_3=v*G_2*v

  // -- Update Delta.o2up
  if(af_bethe_flag == 1) _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o3do, tfunc2sin, tfunc2cos, D0, bath_flag); // assume that tfunc2 is 2-3 hopping
  else                   _component_update(tstp, this->o2up, this->o2up_cc, Gloc.o3up, tfunc2sin, tfunc2cos, D0, bath_flag);

  // -- Update Delta.o2do
  if(af_bethe_flag == 1) _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o3up, tfunc2sin, tfunc2cos, D0, bath_flag);
  else                   _component_update(tstp, this->o2do, this->o2do_cc, Gloc.o3do, tfunc2sin, tfunc2cos, D0, bath_flag);

  // -- Update Delta.o3up
  if(af_bethe_flag == 1) _component_update(tstp, this->o3up, this->o3up_cc, Gloc.o2do, tfunc3sin, tfunc3cos, D0, bath_flag); // assume that tfunc3 is 3-2 hopping
  else                   _component_update(tstp, this->o3up, this->o3up_cc, Gloc.o2up, tfunc3sin, tfunc3cos, D0, bath_flag);

  // -- Update Delta.o3do
  if(af_bethe_flag == 1) _component_update(tstp, this->o3do, this->o3do_cc, Gloc.o2up, tfunc3sin, tfunc3cos, D0, bath_flag);
  else                   _component_update(tstp, this->o3do, this->o3do_cc, Gloc.o2do, tfunc3sin, tfunc3cos, D0, bath_flag);

}

// ---------------------------------------------------------------------------

*/





/*
phi_vec1[0]=0.;
for (int i=1;i<nt+1;i++)
{
  phi_vec1[i]*=h;
  phi_vec1[i]+=phi_vec1[i-1];
}
phi_vec2[0]=0.;
for (int i=1;i<nt+1;i++)
{
  phi_vec2[i]*=h;
  phi_vec2[i]+=phi_vec2[i-1];
}
phi_vec3[0]=0.;
for (int i=1;i<nt+1;i++)
{
  phi_vec3[i]*=h;
  phi_vec3[i]+=phi_vec3[i-1];
}

cntr::function<double> tfunc1sin(nt), tfunc2sin(nt), tfunc3sin(nt);
for(int tstp=-1; tstp <= nt; tstp++) tfunc1sin[tstp] = t_vec1[tstp+1]*sin(phi_vec1[tstp+1]);
for(int tstp=-1; tstp <= nt; tstp++) tfunc2sin[tstp] = t_vec2[tstp+1]*sin(phi_vec2[tstp+1]);
for(int tstp=-1; tstp <= nt; tstp++) tfunc3sin[tstp] = t_vec3[tstp+1]*sin(phi_vec3[tstp+1]);

cntr::function<double> tfunc1cos(nt), tfunc2cos(nt), tfunc3cos(nt);
for(int tstp=-1; tstp <= nt; tstp++) tfunc1cos[tstp] = t_vec1[tstp+1]*cos(phi_vec1[tstp+1]);
for(int tstp=-1; tstp <= nt; tstp++) tfunc2cos[tstp] = t_vec2[tstp+1]*cos(phi_vec2[tstp+1]);
for(int tstp=-1; tstp <= nt; tstp++) tfunc3cos[tstp] = t_vec3[tstp+1]*cos(phi_vec3[tstp+1]);
*/
