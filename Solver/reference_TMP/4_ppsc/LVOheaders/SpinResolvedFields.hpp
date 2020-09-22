
#ifndef _SPINRESOLVEDFIELDS_HPP
#define _SPINRESOLVEDFIELDS_HPP

// -----------------------------------------------------------------------
//
// Spin resolved fields class for ppsc
//
// Author: Francesco Petocchi, francesco.petocchi@gmail.com (2018)
//
// -----------------------------------------------------------------------

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

  void update(int tstp, ppsc::gf_tstps_type & gf_tstps,
    double linear_mixing=0.0, int af_bethe_flag=0)
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
      if(af_bethe_flag==0)
      {
        std::cout << "--> spin symm " << std::endl;

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

  void spinflip(int tstp)
  {
    this->o1up.set_timestep(tstp, this->o1do);
    this->o1do.set_timestep(tstp, this->o1up);
    this->o2up.set_timestep(tstp, this->o2do);
    this->o2do.set_timestep(tstp, this->o2up);
    this->o3up.set_timestep(tstp, this->o3do);
    this->o3do.set_timestep(tstp, this->o3up);
  }

  // ---------------------------------------------------------------------------

  void copy(int tstp, single_particle_greensfunction_type & G)
  {
    this->o1up.set_timestep(tstp, G.o1up);
    this->o1do.set_timestep(tstp, G.o1do);
    this->o2up.set_timestep(tstp, G.o2up);
    this->o2do.set_timestep(tstp, G.o2do);
    this->o3up.set_timestep(tstp, G.o3up);
    this->o3do.set_timestep(tstp, G.o3do);
  }

  // ---------------------------------------------------------------------------

  void null()
  {
    this->o1up.clear();
    this->o1do.clear();
    this->o2up.clear();
    this->o2do.clear();
    this->o3up.clear();
    this->o3do.clear();
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


ppsc::gf_tstps_type vectorize_gf(int tstp,
  single_particle_greensfunction_type & G, int af_bethe_flag=0)
{
  ppsc::gf_tstps_type Gvec;
  ppsc::gf_tstp_type Gtmp(tstp, G.o1up.ntau(), 1, -1);
  if(af_bethe_flag==0)
  {
    Gtmp.clear();G.o1up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o1do.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o2up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o2do.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o3up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o3do.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
  }
  else
  {
    Gtmp.clear();G.o1do.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o1up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o2do.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o2up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o3do.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
    Gtmp.clear();G.o3up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
  }
  return Gvec;
};


////////////////////////////////////////////////////////////////////////////////


ppsc::gf_tstps_type add_gf(int tstp,
                           ppsc::gf_tstps_type & Sloc,
                           ppsc::gf_tstps_type & G,
                           cntr::function<double> & tmatrix, int periodic)
{
  //----------------------------------------------------------------------------
  int Norb=3;
  int ntau=G[0].ntau();
  cntr::herm_matrix_timestep<double> S(tstp, ntau, 3, -1);
  cntr::herm_matrix_timestep<double> tGt(tstp, ntau, 1, -1);
  cntr::function<double> T(tstp,Norb), Tdag(tstp,Norb);
  Eigen::MatrixXcd tmpT(Norb,Norb), tmpTdag(Norb,Norb);
  T=tmatrix;
  for(int it=-1; it <= tstp; it++)
  {
    tmatrix.get_value(it,tmpT);
    tmpTdag = tmpT.adjoint().eval();
    Tdag.set_value(it,tmpTdag);
  }
  //----------------------------------------------------------------------------

  // UP BLOCK
  {
    //FORWARD
    {
      S.clear();
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i;
        S.set_matrixelement( i, i, G[Deltandx], 0, 0 );
      }
      //
      S.left_multiply(T);
      S.right_multiply(Tdag);
      //
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i;
        for(int j=0; j < Norb; j++)
        {
          tGt.clear();
          tGt.set_matrixelement(0,0,S,i,j); Sloc[Deltandx].incr(tGt,1.0);
        }
      }
    }
    //
    //BACKWARD
    if(periodic==1)
    {
      S.clear();
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i;
        S.set_matrixelement( i, i, G[Deltandx], 0, 0 );
      }
      //
      S.left_multiply(Tdag);
      S.right_multiply(T);
      //
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i;
        for(int j=0; j < Norb; j++)
        {
          tGt.clear();
          tGt.set_matrixelement(0,0,S,i,j); Sloc[Deltandx].incr(tGt,1.0);
        }
      }
    }
  }
  // DO BLOCK
  {
    //FORWARD
    {
      S.clear();
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i+1;
        S.set_matrixelement( i, i, G[Deltandx], 0, 0 );
      }
      //
      S.left_multiply(T);
      S.right_multiply(Tdag);
      //
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i+1;
        for(int j=0; j < Norb; j++)
        {
          tGt.clear();
          tGt.set_matrixelement(0,0,S,i,j); Sloc[Deltandx].incr(tGt,1.0);
        }
      }
    }
    //
    //BACKWARD
    if(periodic==1)
    {
      S.clear();
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i+1;
        S.set_matrixelement( i, i, G[Deltandx], 0, 0 );
      }
      //
      S.left_multiply(Tdag);
      S.right_multiply(T);
      //
      for(int i=0; i < Norb; i++)
      {
        int Deltandx = 2*i+1;
        for(int j=0; j < Norb; j++)
        {
          tGt.clear();
          tGt.set_matrixelement(0,0,S,i,j); Sloc[Deltandx].incr(tGt,1.0);
        }
      }
    }
  }
  return Sloc;
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

  void update(int tstp, int ndx,
                        single_particle_greensfunction_type & Gloc,
                        single_particle_greensfunction_type & Gplane,
                        single_particle_greensfunction_type & Gvert_abv,
                        single_particle_greensfunction_type & Gvert_blw,
                        single_particle_greensfunction_type & Gdiag_abv,
                        single_particle_greensfunction_type & Gdiag_blw,
                        // hopping matrix in the plane
                        cntr::function<double> & hnn,
                        cntr::function<double> & tplane,
                        cntr::function<double> & hxnnn,
                        cntr::function<double> & hynnn,
                        // hopping matrix vertical
                        cntr::function<double> & tvert_abv,
                        cntr::function<double> & tvert_blw,
                        cntr::function<double> & tdiag_abv,
                        cntr::function<double> & tdiag_blw,
                        // boson bath
                        cntr::herm_matrix<double> & D0,
                        int af_bethe_flag=0, int phonon_flag=0, int j_flag=0,
                        double beta=40, double h=0.015, int kt=5)
  {
    //std::cout << "delta update site:  "<< ndx << " timestep: " << tstp << endl;
    //
    std::string ilayer = std::to_string(ndx+1);
    std::string hdrp = "F_plane_L_"; hdrp.append(ilayer).append(".out");
    std::string hdrr = "F_right_L_"; hdrr.append(ilayer).append(".out");
    std::string hdrl = "F_left_L_";  hdrl.append(ilayer).append(".out");
    const char * Fleftfile  = hdrl.c_str();// std::cout << Fleftfile << endl;
    const char * Frightfile = hdrr.c_str();// std::cout << Frightfile << endl;
    const char * Fplanefile = hdrp.c_str();// std::cout << Fplanefile << endl;

    // Gfs initialization
    ppsc::gf_tstps_type Gnn,Gnnn,Gv_a,Gv_b,Gd_a,Gd_b;
    // plane nearest-neighbors - the Gf is switched/spin-flipped
    Gnn = vectorize_gf(tstp, Gplane, af_bethe_flag );
    // plane next-nearest-neighbors - the Gf is the same and spin is conserved
    Gnnn= vectorize_gf(tstp, Gloc, 0 );
    // vertical above - the spin is conserved
    Gv_a= vectorize_gf(tstp, Gvert_abv, 0 );
    // vertical below - the spin is conserved
    Gv_b= vectorize_gf(tstp, Gvert_blw, 0 );
    // diagonal above - the spin is flipped
    Gd_a= vectorize_gf(tstp, Gdiag_abv, af_bethe_flag );
    // diagonal below - the spin is flipped
    Gd_b= vectorize_gf(tstp, Gdiag_blw, af_bethe_flag );

    // Potentials iitialization
    ppsc::gf_tstps_type Splane,SR,SL;
    {
      cntr::herm_matrix_timestep<double> Stmp(tstp, ntau, 1, -1);
      for(int i=0; i<6; i++)
      {
        Splane.push_back(Stmp);
        SR.push_back(Stmp);
        SL.push_back(Stmp);
      }
    }

    //periodic part
    // nearest-neighbors in the two directions
    Splane = add_gf(tstp,Splane,Gnn,hnn,1);
    Splane = add_gf(tstp,Splane,Gnn,tplane,1);
    // next-nearest-neighbors in the two directions
    Splane = add_gf(tstp,Splane,Gnnn,hxnnn,1);
    Splane = add_gf(tstp,Splane,Gnnn,hynnn,1);
    //
    if( (tstp>=(kt+2)) && (j_flag==1) )
    {
      std::vector<std::complex<double>> Fvect;
      Fvect = get_S_star_Gles_timestep(tstp,Splane,Gnnn,beta,h,kt);
      FILE * outFile;
      outFile = fopen (Fplanefile,"a");
      fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
      for(int i=0; i<6; i++)
      {
        fprintf (outFile, "\t%.20e",Fvect[i].imag());
      }
      for(int i=0; i<6; i++)
      {
        fprintf (outFile, "\t%.20e",Fvect[i].real());
      }
      fprintf (outFile, "\n");
      fclose(outFile);
    }
    //
    // local part-1
    // corrents to above or left interface
    SL = add_gf(tstp,SL,Gv_a,tvert_abv, 0);
    SL = add_gf(tstp,SL,Gd_a,tdiag_abv, 0);
    //
    if( (tstp>=(kt+2)) && (j_flag==1) )
    {
      std::vector<std::complex<double>> Fvect;
      Fvect = get_S_star_Gles_timestep(tstp,SL,Gnnn,beta,h,kt);
      FILE * outFile;
      outFile = fopen (Fleftfile,"a");
      fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
      for(int i=0; i<6; i++)
      {
        fprintf (outFile, "\t%.20e",Fvect[i].imag());
      }
      for(int i=0; i<6; i++)
      {
        fprintf (outFile, "\t%.20e",Fvect[i].real());
      }
      fprintf (outFile, "\n");
      fclose(outFile);
    }
    // local part-2
    // corrents to below or right interface
    SR = add_gf(tstp,SR,Gv_b,tvert_blw, 0);
    SR = add_gf(tstp,SR,Gd_b,tdiag_blw, 0);
    //
    if( (tstp>=(kt+2)) && (j_flag==1) )
    {
      std::vector<std::complex<double>> Fvect;
      Fvect = get_S_star_Gles_timestep(tstp,SR,Gnnn,beta,h,kt);
      FILE * outFile;
      outFile = fopen (Frightfile,"a");
      fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
      for(int i=0; i<6; i++)
      {
        fprintf (outFile, "\t%.20e",Fvect[i].imag());
      }
      for(int i=0; i<6; i++)
      {
        fprintf (outFile, "\t%.20e",Fvect[i].real());
      }
      fprintf (outFile, "\n");
      fclose(outFile);
    }

    // Filling components without/with magnatic ordering
    // up block
    this->o1up.set_timestep(tstp, Splane[0]);this->o1up.incr_timestep(tstp, SL[0]);this->o1up.incr_timestep(tstp, SR[0]);
    this->o2up.set_timestep(tstp, Splane[2]);this->o2up.incr_timestep(tstp, SL[2]);this->o2up.incr_timestep(tstp, SR[2]);
    this->o3up.set_timestep(tstp, Splane[4]);this->o3up.incr_timestep(tstp, SL[4]);this->o3up.incr_timestep(tstp, SR[4]);
    // do block
    this->o1do.set_timestep(tstp, Splane[1]);this->o1do.incr_timestep(tstp, SL[1]);this->o1do.incr_timestep(tstp, SR[1]);
    this->o2do.set_timestep(tstp, Splane[3]);this->o2do.incr_timestep(tstp, SL[3]);this->o2do.incr_timestep(tstp, SR[3]);
    this->o3do.set_timestep(tstp, Splane[5]);this->o3do.incr_timestep(tstp, SL[5]);this->o3do.incr_timestep(tstp, SR[5]);

    // adding boson fot the specific site
    if(phonon_flag)
    {
      std::cout << "Phonon flag: " << phonon_flag << endl;
      //up-block
		{
        cntr::herm_matrix_timestep<double> tmp_phonon(tstp, ntau, 1, +1);
        Bubble2(tstp, tmp_phonon, Gloc.o1up, D0);
        this->o1up.incr_timestep(tstp, tmp_phonon);
      }
		{
        cntr::herm_matrix_timestep<double> tmp_phonon(tstp, ntau, 1, +1);
        Bubble2(tstp, tmp_phonon, Gloc.o2up, D0);
        this->o1up.incr_timestep(tstp, tmp_phonon);
      }
		{
        cntr::herm_matrix_timestep<double> tmp_phonon(tstp, ntau, 1, +1);
        Bubble2(tstp, tmp_phonon, Gloc.o3up, D0);
        this->o1up.incr_timestep(tstp, tmp_phonon);
      }
      //do-block
		{
        cntr::herm_matrix_timestep<double> tmp_phonon(tstp, ntau, 1, +1);
        Bubble2(tstp, tmp_phonon, Gloc.o1do, D0);
        this->o1do.incr_timestep(tstp, tmp_phonon);
      }
		{
        cntr::herm_matrix_timestep<double> tmp_phonon(tstp, ntau, 1, +1);
        Bubble2(tstp, tmp_phonon, Gloc.o2do, D0);
        this->o1do.incr_timestep(tstp, tmp_phonon);
      }
		{
        cntr::herm_matrix_timestep<double> tmp_phonon(tstp, ntau, 1, +1);
        Bubble2(tstp, tmp_phonon, Gloc.o3do, D0);
        this->o1do.incr_timestep(tstp, tmp_phonon);
      }
    }


    // setting other component
    ppsc::set_bwd_from_fwd(tstp, this->o1up_cc, this->o1up);
    ppsc::set_bwd_from_fwd(tstp, this->o2up_cc, this->o2up);
    ppsc::set_bwd_from_fwd(tstp, this->o3up_cc, this->o3up);
    ppsc::set_bwd_from_fwd(tstp, this->o1do_cc, this->o1do);
    ppsc::set_bwd_from_fwd(tstp, this->o2do_cc, this->o2do);
    ppsc::set_bwd_from_fwd(tstp, this->o3do_cc, this->o3do);
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


void buildBath(single_particle_greensfunction_type & Gbath,
               solver_type & S, double lead_W, double Vbias, double shift)
{
	int nt=S.nt();
	int ntau=S.ntau();
  int kt=S.kt();
  double beta=S.beta();
  double h=S.h();
	//
	//GBATH creation
	//up Block
	{
		//orb 1
		double cl=-lead_W+S.imp.hamiltonian.E1[0]-S.imp.hamiltonian.mu+shift;
		double cr=+lead_W+S.imp.hamiltonian.E1[0]-S.imp.hamiltonian.mu+shift;
		cntr::smooth_box dos(cl,cr,beta);
		std::cout << "computing flat_dos: orb 1 up" << std::endl;
    //std::cout << "beta: " << beta << " W: " << lead_W << std::endl;
    //std::cout << "bot: " << cl << " top: " << cr << std::endl;
		cntr::green_equilibrium(Gbath.o1up,dos,beta,h);
    std::cout << "copying flat_dos: 1up-->1do" << std::endl;
		Gbath.o1do = Gbath.o1up;
	}
	{
		//orb 2
		double cl=-lead_W+S.imp.hamiltonian.E2[0]-S.imp.hamiltonian.mu+shift;
		double cr=+lead_W+S.imp.hamiltonian.E2[0]-S.imp.hamiltonian.mu+shift;
		cntr::smooth_box dos(cl,cr,beta);
		std::cout << "computing flat_dos: orb 2 up" << std::endl;
    //std::cout << "beta: " << beta << " W: " << lead_W << std::endl;
    //std::cout << "bot: " << cl << " top: " << cr << std::endl;
		cntr::green_equilibrium(Gbath.o2up,dos,beta,h);
    std::cout << "copying flat_dos: 2up-->2do" << std::endl;
		Gbath.o2do = Gbath.o2up;
	}
	{
		//orb 3
		double cl=-lead_W+S.imp.hamiltonian.E3[0]-S.imp.hamiltonian.mu+shift;
		double cr=+lead_W+S.imp.hamiltonian.E3[0]-S.imp.hamiltonian.mu+shift;
		cntr::smooth_box dos(cl,cr,beta);
		std::cout << "computing flat_dos: orb 3 up" << std::endl;
    //std::cout << "beta: " << beta << " W: " << lead_W << std::endl;
    //std::cout << "bot: " << cl << " top: " << cr << std::endl;
		cntr::green_equilibrium(Gbath.o3up,dos,beta,h);
    std::cout << "copying flat_dos: 3up-->3do" << std::endl;
		Gbath.o3do = Gbath.o3up;
	}
	//
	//GBATH external bias
	if(Vbias != 0.0)
	{
		std::vector<double> intdE(nt+2,0),dE(nt+2,0);
		cntr::function<double> psi(nt),psidag(nt);
		//
		for(int tstp=-1;tstp<=nt;tstp++) dE[tstp+1]=(
			S.imp.hamiltonian.E1[tstp+1]-S.imp.hamiltonian.E1[0]+
			S.imp.hamiltonian.E2[tstp+1]-S.imp.hamiltonian.E2[0]+
			S.imp.hamiltonian.E3[tstp+1]-S.imp.hamiltonian.E3[0])/3.0;
		integrate_time(intdE,dE,h,kt);
		for(int tstp=-1;tstp<=nt;tstp++)
		{
			double sde=sin(intdE[tstp+1]);
			double cde=cos(intdE[tstp+1]);
			psi[tstp]= std::complex<double>(cde,-sde);
			psidag[tstp]=conj(psi[tstp]);
		}
		//
		for(int tstp=-1;tstp<=nt;tstp++)
		{
			// up block
			Gbath.o1up.left_multiply(tstp,psi); Gbath.o1up.right_multiply(tstp,psidag);
			Gbath.o2up.left_multiply(tstp,psi); Gbath.o2up.right_multiply(tstp,psidag);
			Gbath.o3up.left_multiply(tstp,psi); Gbath.o3up.right_multiply(tstp,psidag);
			//do Block
      Gbath.o1do.left_multiply(tstp,psi); Gbath.o1do.right_multiply(tstp,psidag);
      Gbath.o2do.left_multiply(tstp,psi); Gbath.o2do.right_multiply(tstp,psidag);
      Gbath.o3do.left_multiply(tstp,psi); Gbath.o3do.right_multiply(tstp,psidag);
		}
	}
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
};

template<class HAM>
ppsc::gf_verts_type get_gf_verts(HAM & h)
{
  ppsc::gf_verts_type gf_verts;
  // up Block
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, h.c1ua, h.c1uc)); // spin up orb 1
	gf_verts.push_back(ppsc::gf_vert_type(1, 1, h.c2ua, h.c2uc)); // spin up orb 2
	gf_verts.push_back(ppsc::gf_vert_type(2, 2, h.c3ua, h.c3uc)); // spin up orb 3
  // do Block
	gf_verts.push_back(ppsc::gf_vert_type(0, 0, h.c1da, h.c1dc)); // spin do orb 1
  gf_verts.push_back(ppsc::gf_vert_type(1, 1, h.c2da, h.c2dc)); // spin do orb 2
  gf_verts.push_back(ppsc::gf_vert_type(2, 2, h.c3da, h.c3dc)); // spin do orb 3
  return gf_verts;
};


// -----------------------------------------------------------------------------

#endif // _SPINRESOLVEDFIELDS_HPP
