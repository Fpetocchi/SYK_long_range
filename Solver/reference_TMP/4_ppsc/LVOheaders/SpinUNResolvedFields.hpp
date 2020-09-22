
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
      o1up(nt, ntau, 1, -1),
      o2up(nt, ntau, 1, -1),
      o3up(nt, ntau, 1, -1)
      {}

   // --------------------------------------------------------------------------

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
         this->o2up.get_timestep(-1, gloc_old);
         gloc_mix.incr(gloc_old, linear_mixing);
         gloc_mix.incr(gf_tstps[1], 1.0 - linear_mixing);
         this->o2up.set_timestep(-1, gloc_mix);

         gloc_mix.clear();
         this->o3up.get_timestep(-1, gloc_old);
         gloc_mix.incr(gloc_old, linear_mixing);
         gloc_mix.incr(gf_tstps[2], 1.0 - linear_mixing);
         this->o3up.set_timestep(-1, gloc_mix);
      }
      else // TSTP UPDATE WITH NO IMPOSED SYMMETRIES
      {
         this->o1up.set_timestep(tstp, gf_tstps[0]);
         this->o2up.set_timestep(tstp, gf_tstps[1]);
         this->o3up.set_timestep(tstp, gf_tstps[2]);
      }
   }

   // --------------------------------------------------------------------------

   void copy(int tstp, single_particle_greensfunction_type & G)
   {
      this->o1up.set_timestep(tstp, G.o1up);
      this->o2up.set_timestep(tstp, G.o2up);
      this->o3up.set_timestep(tstp, G.o3up);
   }

   // --------------------------------------------------------------------------

   void null()
   {
      this->o1up.clear();
      this->o2up.clear();
      this->o3up.clear();
   }

   // --------------------------------------------------------------------------

   void store(hid_t file_id)
   {
      hid_t group_id;
      group_id = create_group(file_id, "g1u");
      store_herm_greens_function(group_id, this->o1up);
      close_group(group_id);
      group_id = create_group(file_id, "g2u");
      store_herm_greens_function(group_id, this->o2up);
      close_group(group_id);
      group_id = create_group(file_id, "g3u");
      store_herm_greens_function(group_id, this->o3up);
      close_group(group_id);
   }

   // --------------------------------------------------------------------------

   void load(std::string filename)
   {
      this->o1up.read_from_hdf5(filename.c_str(), "g1u");
      this->o2up.read_from_hdf5(filename.c_str(), "g2u");
      this->o3up.read_from_hdf5(filename.c_str(), "g3u");
   }

   // --------------------------------------------------------------------------

   int nt, ntau;
   cntr::herm_matrix<double> o1up, o2up, o3up;

};


////////////////////////////////////////////////////////////////////////////////


ppsc::gf_tstps_type vectorize_gf(int tstp,
   single_particle_greensfunction_type & G, int af_bethe_flag=0)
{
   ppsc::gf_tstps_type Gvec;
   //
   {
      ppsc::gf_tstp_type Gtmp(tstp, G.o1up.ntau(), 1, -1);
      G.o1up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
   }
   {
      ppsc::gf_tstp_type Gtmp(tstp, G.o2up.ntau(), 1, -1);
      G.o2up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
   }
   {
      ppsc::gf_tstp_type Gtmp(tstp, G.o3up.ntau(), 1, -1);
      G.o3up.get_timestep(tstp,Gtmp);Gvec.push_back(Gtmp);
   }
   //
   return Gvec;
};


////////////////////////////////////////////////////////////////////////////////


ppsc::gf_tstps_type add_gf(int tstp, ppsc::gf_tstps_type & Sloc,
                                     ppsc::gf_tstps_type & G,
                                     cntr::function<double> & tmatrix)

{
   //---------------------------------------------------------------------------
   int Norb=3;
   int ntau=G[0].ntau();
   cntr::function<double> T(tstp,Norb), Tdag(tstp,Norb);
   Eigen::MatrixXcd tmpT(Norb,Norb), tmpTdag(Norb,Norb);
   for(int it=-1; it <= tstp; it++)
   {
      tmatrix.get_value(it,tmpT);
      T.set_value(it,tmpT);
      tmpTdag = tmpT.adjoint().eval();
      Tdag.set_value(it,tmpTdag);
   }
   //---------------------------------------------------------------------------
   // UP BLOCK
   {
      cntr::herm_matrix_timestep<double> tGt(tstp, ntau, 3, -1);
      for(int i=0; i < Norb; i++)
      {
         int Deltandx = i;
         tGt.set_matrixelement( i, i, G[Deltandx], 0, 0 );
      }
      //
      tGt.left_multiply(T);
      tGt.right_multiply(Tdag);
      //
      for(int i=0; i < Norb; i++)
      {
         cntr::herm_matrix_timestep<double> tGt_element(tstp, ntau, 1, -1);
         int Deltandx = i;
         tGt_element.set_matrixelement(0,0,tGt,i,i);
         Sloc[Deltandx].incr(tGt_element,1.0);
      }
      /* OLD WRONG
      for(int i=0; i < Norb; i++)
      {
         int Deltandx = i;
         for(int j=0; j < Norb; j++)
         {
            cntr::herm_matrix_timestep<double> tGt_element(tstp, ntau, 1, -1);
            tGt_element.set_matrixelement(0,0,tGt,i,j);
            Sloc[Deltandx].incr(tGt_element,1.0);
         }
      }
      */
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
      o2up(nt, ntau, 1, -1),
      o3up(nt, ntau, 1, -1),
      // -- bwd hybridizations
      o1up_cc(nt, ntau, 1, -1),
      o2up_cc(nt, ntau, 1, -1),
      o3up_cc(nt, ntau, 1, -1)
      {}

   // --------------------------------------------------------------------------

   void update(int tstp, int ndx,
                        single_particle_greensfunction_type & Gloc,
                        single_particle_greensfunction_type & Gplane,
                        single_particle_greensfunction_type & Gvert_abv,
                        single_particle_greensfunction_type & Gvert_blw,
                        // hopping matrix in the plane
                        cntr::function<double> & t_nn_Uc  ,
                        cntr::function<double> & t_nn_Rd  ,
                        cntr::function<double> & t_nnRyp  ,
                        cntr::function<double> & t_nnRxm  ,
                        cntr::function<double> & tnnnRxp  ,
                        cntr::function<double> & tnnnRyp  ,
                        cntr::function<double> & tnnnRxm  ,
                        cntr::function<double> & tnnnRym  ,
                        // hopping matrix vertical
                        cntr::function<double> & tvert_abv,
                        cntr::function<double> & tvert_blw,
                        //
                        cntr::function<double> & psi,
                        cntr::function<double> & psidag,
                        // boson bath
                        cntr::herm_matrix<double> & D0,
                        int af_bethe_flag=0, int phonon_flag=0, int j_flag=0,
                        double beta=40, double h=0.015, int kt=5,
                        double gradientBulk=0.0, int Nimp=1)
   {
      //std::cout << "delta update site:  "<< ndx << " timestep: " << tstp << endl;
      //
      int Norb=3;
      std::string ilayer = std::to_string(ndx+1);
      //
      std::string hdrp = "F_plane_L_"; hdrp.append(ilayer).append(".out");
      const char * Fplanefile = hdrp.c_str();
      //
      std::string hdrbot = "F_bot_L_"; hdrbot.append(ilayer).append(".out");
      std::string hdrtop = "F_top_L_"; hdrtop.append(ilayer).append(".out");
      const char * Fbotfile = hdrbot.c_str();
      const char * Ftopfile = hdrtop.c_str();

      // ----------------------------------------------------------------------
      // adding the optional phase factor for a probe bias in the bulk
      cntr::function<double> tva(Gloc.o1up.nt(),Norb);tva=tvert_abv;
      cntr::function<double> tvb(Gloc.o1up.nt(),Norb);tvb=tvert_blw;
      if(tstp==1&&ndx==0)std::cout<<"gradientBulk=  "<<gradientBulk<<std::endl;
      if(gradientBulk!=0.0&&tstp>=0)
      {
         cntr::function<double> psiV(Gloc.o1up.nt(),Norb);
         cntr::function<double> psidagV(Gloc.o1up.nt(),Norb);
         for(int i=0; i < Norb; i++)
         {
            psiV.set_matrixelement(i,i,psi,0,0);
            psidagV.set_matrixelement(i,i,psidag,0,0);
         }
         tva.left_multiply(psiV);
         tvb.left_multiply(psidagV);
      }

      // ----------------------------------------------------------------------
      // Gfs initialization
      ppsc::gf_tstps_type Gnn,Gnnn,Gv_a,Gv_b;
      // plane nearest-neighbors - the Gf is switched/spin-flipped
      Gnn = vectorize_gf(tstp, Gplane, af_bethe_flag );
      // plane next-nearest-neighbors - the Gf is the same and spin is conserved
      Gnnn= vectorize_gf(tstp, Gloc, 0 );
      // vertical above - the spin is conserved
      Gv_a= vectorize_gf(tstp, Gvert_abv, 0 );
      // vertical below - the spin is conserved
      Gv_b= vectorize_gf(tstp, Gvert_blw, 0 );

      // ----------------------------------------------------------------------
      // Potentials iitialization
      ppsc::gf_tstps_type Splane,Sbot,Stop;
      {
         for(int i=0; i<Norb; i++)
         {
            Splane.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            Sbot.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            Stop.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
         }
      }
      //in-plane nearest-neighbors
      Splane = add_gf(tstp,Splane,Gnn ,t_nn_Uc);
      Splane = add_gf(tstp,Splane,Gnn ,t_nn_Rd);
      Splane = add_gf(tstp,Splane,Gnn ,t_nnRyp);
      Splane = add_gf(tstp,Splane,Gnn ,t_nnRxm);
      //in-plane next-nearest-neighbors
      Splane = add_gf(tstp,Splane,Gnnn,tnnnRxp);
      Splane = add_gf(tstp,Splane,Gnnn,tnnnRyp);
      Splane = add_gf(tstp,Splane,Gnnn,tnnnRxm);
      Splane = add_gf(tstp,Splane,Gnnn,tnnnRym);
      //
      if( (tstp>=(kt+2)) && (j_flag==1) )
      {
         std::vector<std::complex<double>> Fvect;
         Fvect = get_S_star_Gles_timestep(tstp,Splane,Gnnn,beta,h,kt);
         FILE * outFile;
         outFile = fopen (Fplanefile,"a");
         fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
         for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
         for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
         fprintf (outFile, "\n");
         fclose(outFile);
      }
      //
      //vertical above
      Stop = add_gf(tstp,Stop,Gv_a,tva);
      //
      if( (tstp>=(kt+2)) && (j_flag==1) )
      {
         std::vector<std::complex<double>> Fvect;
         Fvect = get_S_star_Gles_timestep(tstp,Stop,Gnnn,beta,h,kt);
         FILE * outFile;
         outFile = fopen (Ftopfile,"a");
         fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
         for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
         for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
         fprintf (outFile, "\n");
         fclose(outFile);
      }
      //vertical below or right interface
      Sbot = add_gf(tstp,Sbot,Gv_b,tvb);
      //
      if( (tstp>=(kt+2)) && (j_flag==1) )
      {
         std::vector<std::complex<double>> Fvect;
         Fvect = get_S_star_Gles_timestep(tstp,Sbot,Gnnn,beta,h,kt);
         FILE * outFile;
         outFile = fopen (Fbotfile,"a");
         fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
         for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
         for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
         fprintf (outFile, "\n");
         fclose(outFile);
      }

      // ----------------------------------------------------------------------
      // Aditional convolutions for absorbed energy
      if( (tstp>=(kt+2)) && (j_flag==1) && (Nimp==1))
      {
         // NEXT NEIGHBORS ON THE DIAGONAL
         // - diagonal 1 - forward
         {
            std::string conv_Exy_Uc = "conv_Exy_Uc_L_";  conv_Exy_Uc.append(ilayer);
            conv_Exy_Uc.append("_fwd.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnn,t_nn_Uc);
            //
            const char * conv = conv_Exy_Uc.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
         // - diagonal 1 - backward
         {
            std::string conv_Exy_Uc = "conv_Exy_Uc_L_";  conv_Exy_Uc.append(ilayer);
            conv_Exy_Uc.append("_bkw.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnn,t_nn_Rd);
            //
            const char * conv = conv_Exy_Uc.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
         // - diagonal 2 - forward
         {
            std::string conv_Exy_Ry = "conv_Exy_Ry_L_";  conv_Exy_Ry.append(ilayer);
            conv_Exy_Ry.append("_fwd.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnn,t_nnRyp);
            //
            const char * conv = conv_Exy_Ry.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
         // - diagonal 2 - backward
         {
            std::string conv_Exy_Ry = "conv_Exy_Ry_L_";  conv_Exy_Ry.append(ilayer);
            conv_Exy_Ry.append("_bkw.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnn,t_nnRxm);
            //
            const char * conv = conv_Exy_Ry.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
         //
         // NEXT-NEXT NEIGHBORS ON THE HORIZONTAL
         // - horizontal 1  - forward
         {
            std::string conv_Ex_nnn = "conv_Ex_nnn_L_";  conv_Ex_nnn.append(ilayer);
            conv_Ex_nnn.append("_fwd.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnnn,tnnnRxp);
            //
            const char * conv = conv_Ex_nnn.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
         // - horizontal 1  - backward
         {
            std::string conv_Ex_nnn = "conv_Ex_nnn_L_";  conv_Ex_nnn.append(ilayer);
            conv_Ex_nnn.append("_bkw.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnnn,tnnnRxm);
            //
            const char * conv = conv_Ex_nnn.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
         // - horizontal 2  - forward
         {
            std::string conv_Ey_nnn = "conv_Ey_nnn_L_";  conv_Ey_nnn.append(ilayer);
            conv_Ey_nnn.append("_fwd.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnnn,tnnnRyp);
            //
            const char * conv = conv_Ey_nnn.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
         // - horizontal 2  - backward
         {
            std::string conv_Ey_nnn = "conv_Ey_nnn_L_";  conv_Ey_nnn.append(ilayer);
            conv_Ey_nnn.append("_bkw.out");
            //
            ppsc::gf_tstps_type Sconv;
            for(int i=0; i<Norb; i++)Sconv.push_back(cntr::herm_matrix_timestep<double>(tstp, ntau, 1, -1));
            //
            Sconv = add_gf(tstp,Sconv,Gnnn,tnnnRym);
            //
            const char * conv = conv_Ey_nnn.c_str();
            std::vector<std::complex<double>> Fvect;
            Fvect = get_S_star_Gles_timestep(tstp,Sconv,Gnnn,beta,h,kt);
            FILE * outFile;
            outFile = fopen (conv,"a");
            fprintf (outFile, "%4i"  ,ndx);fprintf (outFile, "%8i"  ,tstp);
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].imag());
            for(int i=0; i<3; i++)fprintf (outFile, "\t%.20e",Fvect[i].real());
            fprintf (outFile, "\n");
            fclose(outFile);
         }
      }

      // Filling components without/with magnatic ordering
      // up block
      this->o1up.set_timestep(tstp, Splane[0]);this->o1up.incr_timestep(tstp, Stop[0]);this->o1up.incr_timestep(tstp, Sbot[0]);
      this->o2up.set_timestep(tstp, Splane[1]);this->o2up.incr_timestep(tstp, Stop[1]);this->o2up.incr_timestep(tstp, Sbot[1]);
      this->o3up.set_timestep(tstp, Splane[2]);this->o3up.incr_timestep(tstp, Stop[2]);this->o3up.incr_timestep(tstp, Sbot[2]);

      // adding boson fot the specific site
      if(phonon_flag)
      {
         std::cout << "Phonon flag: " << phonon_flag << endl;
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
      }

      // setting other component
      ppsc::set_bwd_from_fwd(tstp, this->o1up_cc, this->o1up);
      ppsc::set_bwd_from_fwd(tstp, this->o2up_cc, this->o2up);
      ppsc::set_bwd_from_fwd(tstp, this->o3up_cc, this->o3up);
   }

   // --------------------------------------------------------------------------

   void store(hid_t file_id)
   {
      hid_t group_id;
      group_id = create_group(file_id, "d1u");
      store_herm_greens_function(group_id, this->o1up);
      close_group(group_id);
      group_id = create_group(file_id, "d2u");
      store_herm_greens_function(group_id, this->o2up);
      close_group(group_id);
      group_id = create_group(file_id, "d3u");
      store_herm_greens_function(group_id, this->o3up);
      close_group(group_id);
      group_id = create_group(file_id, "d1ucc");
      store_herm_greens_function(group_id, this->o1up_cc);
      close_group(group_id);
      group_id = create_group(file_id, "d2ucc");
      store_herm_greens_function(group_id, this->o2up_cc);
      close_group(group_id);
      group_id = create_group(file_id, "d3ucc");
      store_herm_greens_function(group_id, this->o3up_cc);
      close_group(group_id);
   }

   // --------------------------------------------------------------------------

   void load(std::string filename)
   {
      this->o1up.read_from_hdf5(filename.c_str(), "d1u");
      this->o2up.read_from_hdf5(filename.c_str(), "d2u");
      this->o3up.read_from_hdf5(filename.c_str(), "d3u");
      this->o1up_cc.read_from_hdf5(filename.c_str(), "d1ucc");
      this->o2up_cc.read_from_hdf5(filename.c_str(), "d2ucc");
      this->o3up_cc.read_from_hdf5(filename.c_str(), "d3ucc");

   }

   // --------------------------------------------------------------------------

   void null()
   {
      this->o1up.clear();
      this->o2up.clear();
      this->o3up.clear();
   }

   // --------------------------------------------------------------------------

   int nt, ntau;
   cntr::herm_matrix<double> o1up, o2up, o3up;
   cntr::herm_matrix<double> o1up_cc, o2up_cc, o3up_cc;

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
   pp_ints.push_back(ppsc::pp_int_type(Delta.o1up,    h.c1dc,  h.c1da,  fermion, fwd)); // spin do fwd
   pp_ints.push_back(ppsc::pp_int_type(Delta.o1up_cc, h.c1da,  h.c1dc,  fermion, bwd)); // spin do bwd
   // orbital no. 2
   pp_ints.push_back(ppsc::pp_int_type(Delta.o2up,    h.c2uc,  h.c2ua,  fermion, fwd)); // spin up fwd
   pp_ints.push_back(ppsc::pp_int_type(Delta.o2up_cc, h.c2ua,  h.c2uc,  fermion, bwd)); // spin up bwd
   pp_ints.push_back(ppsc::pp_int_type(Delta.o2up,    h.c2dc,  h.c2da,  fermion, fwd)); // spin do fwd
   pp_ints.push_back(ppsc::pp_int_type(Delta.o2up_cc, h.c2da,  h.c2dc,  fermion, bwd)); // spin do bwd
   // orbital no. 3
   pp_ints.push_back(ppsc::pp_int_type(Delta.o3up,    h.c3uc,  h.c3ua,  fermion, fwd)); // spin up fwd
   pp_ints.push_back(ppsc::pp_int_type(Delta.o3up_cc, h.c3ua,  h.c3uc,  fermion, bwd)); // spin up bwd
   pp_ints.push_back(ppsc::pp_int_type(Delta.o3up,    h.c3dc,  h.c3da,  fermion, fwd)); // spin do fwd
   pp_ints.push_back(ppsc::pp_int_type(Delta.o3up_cc, h.c3da,  h.c3dc,  fermion, bwd)); // spin do bwd
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
