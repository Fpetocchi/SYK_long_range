#include <ctime>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <complex>
#include <string>
#include <algorithm>


#include <typeinfo>


class QuantumIsing_2D
{
public:
   QuantumIsing_2D (int nx, int ny, bool verbose, bool testing):
   nx(nx), ny(ny), verbose(verbose), testing(testing){}
   int get_nx(void)const {return nx;}
   int get_ny(void)const {return ny;}
   int get_dim(void)const {return nx*ny;}
   int get_Hil(void)const {return pow(2,nx*ny);}
   int get_directions(void)const {return directions;}


   //-------------------------------------------------------------------------OK

   std::vector<int> get_neigh(int ilatt, std::vector<int> &lattice)
   {
      int row = (int)(ilatt/nx);
      int col = ilatt - row * nx;
      //
      int top, bottom, left, right;
      if (col==0)   {left   = lattice[ilatt+nx-1];     } else{left   = lattice[ilatt-1]; }
      if (col==nx-1){right  = lattice[ilatt-nx+1];     } else{right  = lattice[ilatt+1]; }
      if (row==0)   {top    = lattice[ilatt+(ny-1)*nx];} else{top    = lattice[ilatt-nx];}
      if (row==ny-1){bottom = lattice[ilatt-(ny-1)*nx];} else{bottom = lattice[ilatt+nx];}
      //
      return {top,bottom,left,right};
   }

   //-------------------------------------------------------------------------OK

   std::vector<int> dec2bin(int n, int maxdim, bool tospin)
   {
       //
       int initval;
       if(tospin==true){initval=1;}else{initval=0;}
       std::vector<int> binaryNum(maxdim,initval);
       //
       int i = maxdim-1;
       int bin = 0;
       while (n > 0)
       {
          bin = n % 2;
          if(tospin==true)
          {
             if(bin==0) binaryNum[i] = +1 ; // up
             if(bin==1) binaryNum[i] = -1 ; // dw
          }
          else
          {
             binaryNum[i] = bin;
          }
          //
          n = n / 2;
          i--;
       }
       return binaryNum;
   }


   //-------------------------------------------------------------------------OK

   std::vector<std::vector<int>> gen_Hilbert_bin(bool repSz)
   {
      // The Hilbert space is the enseble of all the possible lattice vectors
      std::vector<std::vector<int>> Hilbert( dimHilb , std::vector<int> (dimLatt, 0));

      // Loop over all the possible values of 1-->up and -->dw
      for (int istate = 0; istate < dimHilb; istate++)
      {
         std::vector<int> lattice = dec2bin(istate,dimLatt,repSz);
         Hilbert[istate] = lattice;
      }

      //
      if( verbose ==true ) printf("gen_Hilbert_bin. DONE.\n");
      return Hilbert;
   }

   //-------------------------------------------------------------------------OK

   std::vector<double> gen_Sz()
   {
      // The Hilbert space is the enseble of all the possible lattice vectors
      std::vector<double> Sz( dimHilb , 0);

      // Loop over all the possible values of 1-->up and -->dw
      for (int istate = 0; istate < dimHilb; istate++)
      {
         std::vector<int> state = dec2bin(istate,dimLatt,true);
         for (int ilatt = 0; ilatt < dimLatt; ilatt++)
         {
            Sz[istate] += 0.5 * state[ilatt] / (double)dimLatt;
         }
      }

      //
      if( verbose ==true ) printf("gen_Sz. DONE.\n");
      return Sz;
   }

   //-------------------------------------------------------------------------OK

   std::vector<std::vector<int>> gen_Flips()
   {
      //
      std::vector<std::vector<int> > Flip( dimHilb , std::vector<int> (dimLatt, 0));

      //
      for (int istate = 0; istate < dimHilb; istate++)
      {
         std::vector<int> state = dec2bin(istate,dimLatt,false);
         for (int ilatt = 0; ilatt < dimLatt; ilatt++)
         {
            Flip[istate][ilatt] = istate + pow(2,ilatt) * pow(-1,state[dimLatt-1-ilatt]);
         }
      }

      //
      if( verbose ==true ) printf("gen_Flips. DONE.\n");
      return Flip;
   }

   //-------------------------------------------------------------------------OK

   std::vector<double> init_state_coeff(double a, double f)
   {
      // a c_x coefficient for each configuration
      std::vector<double> state_coeff(dimHilb, 0.0);
      double a_prod;
      double f_prod;

      // count the number of up and the number of pairs with opposite spin
      for (int istate = 0; istate < dimHilb; istate++)
      {
         int Nup=0;
         int Nupdw=0;
         std::vector<int> lattice = dec2bin(istate,dimLatt,true);
         for (int ilatt = 0; ilatt < dimLatt; ilatt++)
         {
            // up count
            if(lattice[ilatt]==1) Nup++;
            // check the neighbors for flips
            std::vector<int> neigh = get_neigh(ilatt, lattice);
            for (int ineig = 0; ineig < 4; ineig++)
            {
               if(lattice[ilatt]*neigh[ineig]<0) Nupdw++;
            }
         }
         //
         a_prod=pow((1.0+a),Nup);
         f_prod=pow((1.0+f),Nupdw);
         //
         state_coeff[istate] = a_prod * f_prod;
      }

      //
      if( verbose ==true ) printf("init_state_coeff. DONE.\n");
      return state_coeff;
   }

   //-------------------------------------------------------------------------OK

   std::vector<double> gen_Hamiltonian_diag()
   {
      //
      std::vector<double> Diag_H(dimHilb, 0.0);

      // Diagonal Ising-like contribution
      for (int istate = 0; istate < dimHilb; istate++)
      {
         std::vector<int> bra = dec2bin(istate,dimLatt,true);
         for (int ilatt = 0; ilatt < dimLatt; ilatt++)
         {
            std::vector<int> neigh = get_neigh(ilatt, bra);
            Diag_H[istate] -= 1.0 * bra[ilatt] * (neigh[0]+neigh[1]+neigh[2]+neigh[3]) * 0.5 ;
         }
      }

      //
      if( verbose ==true ) printf("gen_Hamiltonian_diag. DONE.\n");
      return Diag_H;
   }

   //-------------------------------------------------------------------------OK

   std::vector<std::vector<double>> gen_Hamiltonian_offdiag(std::vector<std::vector<int>> &Flip)
   {
      //
      std::vector<std::vector<double>> OffDiag_H( dimHilb , std::vector<double> (dimLatt, 0.0));

      // Spin-flip contribution - each ket has ilatt corresponding bras
      for (int istate = 0; istate < dimHilb; istate++)
      {
         for (int ilatt = 0; ilatt < dimLatt; ilatt++)
         {
            int connected = Flip[istate][ilatt];
            OffDiag_H[istate][connected] -= 1.0 ;
         }
      }

      //
      if( verbose ==true ) printf("gen_Hamiltonian_offdiag. DONE.\n");
      return OffDiag_H;
   }

   //-------------------------------------------------------------------------OK

   std::vector<std::vector<double>> get_positions()
   {
      //
      std::vector<double> distances( 2 , 0.0);
      std::vector<double> ratios( 2 , 0.0);

      // Always present
      // linear
      ratios[0] = 0.0; distances[0] = 1.0;
      // diagonal
      ratios[1] = 1.0; distances[1] = sqrt(2.0);


      // Direction Identifier
      for (int col_i = nx-1 ; col_i >= 2 ; col_i--)
      {
         for (int row_i = 1 ; row_i < col_i; row_i++)
         {
            // compute the distance
            double Rx = (double)col_i;
            double Ry = (double)row_i;
            double distance = sqrt( pow(Rx,2) + pow(Ry,2)  );

            // compute the ratio (for some reason abs is not working with intel compiler)
            double ratio = sqrt(pow(Ry,2)) / sqrt(pow(Rx,2)); //ratio = abs(Ry) / abs(Rx);

            // other angle entries: check if the same angle is already present
            if(std::find(ratios.begin(), ratios.end(), ratio) != ratios.end())
            {
               // angle already present --> skip
               continue;
            }
            else
            {
               // once checked that the angle is new check that the distance is new as well
               if(std::find(distances.begin(), distances.end(), distance) != distances.end())
               {
                  // distance already present --> skip
                  continue;
               }
               else
               {
                  if(distance>(nx-2)*sqrt(2.0)) continue;
                  // angle & distance not present --> insert
                  ratios.push_back(ratio);
                  distances.push_back(distance);
               }
            }
         }
      }

      // Workaround for the L=4 operators degeneracies
      //ratios.erase(ratios.begin()+2,ratios.begin()+3); distances.erase(distances.begin()+2,distances.begin()+3);
      //ratios.erase(ratios.begin()+2,ratios.begin()+3); distances.erase(distances.begin()+2,distances.begin()+3);
      //ratios.erase(ratios.begin()+3,ratios.begin()+4); distances.erase(distances.begin()+3,distances.begin()+4);

      //
      std::vector<std::vector<double>> positions( ratios.size() , std::vector<double> ( 2, 0.0 ) );

      //
      for (int idir = 0; idir < ratios.size(); idir++)
      {
         positions[idir][0] = ratios[idir];
         positions[idir][1] = distances[idir];
      }

      //
      return positions;
   }


   //-------------------------------------------------------------------------OK


   void set_directions()
   {
      std::vector<std::vector<double>> positions = get_positions();
      directions = positions.size();
      if( verbose ==true )
      {
         for (int idir = 0; idir < directions; idir++)
         {
            printf(" direction: %3i - angle: %2.2f - distance: %2.2f  \n",idir,positions[idir][0],positions[idir][1]);
         }
      }
   }


   //------------------------------------------------------------------------- DOUBLE COUNTINGS NOT CLEAR

   std::vector<std::vector<double>> gen_operators()
   {
      //
      std::vector<std::vector<double>> positions = get_positions();
      std::vector<std::vector<double>> operators( dimHilb , std::vector<double> ( positions.size(), 0.0 ) );
      std::vector<int> Oneigh( positions.size(), 0);

      //
      int tiling = 3;

      // Loop over inner nx*ny sites
      for (int ilatt = 0; ilatt < dimLatt; ilatt++)
      {
         // Loop over outer (nx*tiling)*(ny*tiling) sites
         for (int jlatt = 0; jlatt < dimLatt*tiling*tiling; jlatt++)
         {
            //
            int row_i = ((int)(tiling/2))*ny + (int)(ilatt/nx);
            int col_i = ((int)(tiling/2))*nx + ilatt - (int)(ilatt/nx) * nx;
            //
            int row_j = (int)(jlatt / (nx*tiling));
            int col_j = jlatt - row_j * (nx*tiling);

            // compute the distance
            double Rx = (double) ( col_i - col_j );
            double Ry = (double) ( row_i - row_j );
            double distance = sqrt( pow(Rx,2) + pow(Ry,2)  );

            // compute the ratio
            double ratio;
            double abs_Rx =sqrt(pow(Rx,2));
            double abs_Ry =sqrt(pow(Ry,2));
            if( abs_Ry < abs_Rx )
            {
               ratio = abs_Ry / abs_Rx;
            }
            else
            {
               ratio = abs_Rx / abs_Ry;
            }

            // look if this ratio is present in the list
            for (int idir = 0; idir < positions.size(); idir++)
            {
               if(( ratio == positions[idir][0] )&&( distance == positions[idir][1] ))
               {
                  // Increase the neighbor number
                  Oneigh[idir] ++ ;

                  // find the index of the orifinal lattice corresponding to jsite
                  int row_j2i = row_j % ny ;
                  int col_j2i = col_j % nx ;
                  int j2ilatt = col_j2i + row_j2i*nx;

                  // Update this particular direction in all the spin configuration
                  for (int istate = 0; istate < dimHilb; istate++)
                  {
                     std::vector<int> spinconfig = dec2bin(istate,dimLatt,true);
                     operators[istate][idir] += (spinconfig[ilatt]) * (spinconfig[j2ilatt]);
                  }

                  // testing
                  if( testing ==true )
                  {
                     printf("\n");
                     printf("ilat: %i, row: %i, col: %i\n",ilatt,row_i,col_i);
                     printf("jlat: %i, row: %i, col: %i\n",jlatt,row_j,col_j);
                     printf("ratio: %f  \n",positions[idir][0]);
                     printf("distance: %f  \n",positions[idir][1]);
                     printf("j2ilatt: %i, row: %i, col: %i\n",j2ilatt,row_j2i,col_j2i);
                  }

               }
            }
         }
         if( verbose ==true )printf("ilat: %i done\n",ilatt);
      }

      // final normalization
      for (int istate = 0; istate < dimHilb; istate++)
      {
         for (int idir = 0; idir < positions.size(); idir++)
         {
            operators[istate][idir] /= Oneigh[idir];
         }
      }

      //
      if( verbose ==true )
      {
         printf("\n");
         for (int idir = 0; idir < positions.size(); idir++)printf("idir: %i - Neigh: %i\n",idir,Oneigh[idir]);
         printf("gen_operators. DONE.\n");
      }
      return operators;
   }


   //-------------------------------------------------------------------------OK

   std::vector<std::complex<double>> init_alpha_coeff(std::vector<std::complex<double>> &alpha,
                                                      std::vector<std::vector <double>> &Operators)
   {
      // c_x coefficient for each configuration
      std::vector<std::complex<double>> state_coeff(dimHilb, 0.0);

      // just perform the multiplication
      for (int istate = 0; istate < dimHilb; istate++)
      {
         std::complex<double> exponent=0.0;
         for (int idir = 0; idir < directions; idir++) exponent += alpha[idir]*Operators[istate][idir];
         state_coeff[istate] = exp(exponent) / pow(2,dimLatt/2.);
      }

      //
      if( verbose ==true ) printf("init_alpha_coeff. DONE.\n");
      return state_coeff;
   }

   //-------------------------------------------------------------------------

   void findinvec(std::vector<double> &vec, double &elem, bool &found)
   {
      int n =vec.size();
      //bool
      found = false;
      //
      std::cout << "found " << found  << " siz "  << n << std::endl;
      //
      for (int i = 0; i < n; n++)
      {
         if(vec[i]==elem)
         {
            found = true;
            std::cout << found << std::endl;
            break;
         }
      }
      //
      //return found;
   }

   //-------------------------------------------------------------------------

   // system vars
   int nx,ny;
   int dimHilb=pow(2,nx*ny);
   int dimLatt=nx*ny;
   int directions;
   bool verbose;
   bool testing;

};





////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                               NOT USED CODE                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/*

*/



/*
std::vector<double> get_directional_ratios()
{
   //
   std::vector<double> directional_ratios;

   // Direction Identifier
   for (int col_i = 0 ; col_i < nx; col_i++)
   {
      for (int row_i = 0 ; row_i <= col_i; row_i++)
      {
         // remove self
         if((row_i==0)&&(col_i==0)) continue;

         // compute the distance (not really necessary here)
         double Rx = (double)col_i;
         double Ry = (double)row_i;
         double distance = sqrt( pow(Rx,2) + pow(Ry,2)  );
         if( distance > sqrt(2.0)*(nx-1) ) continue; // this is redundant

         // compute the ratio
         double ratio = abs(Ry) / abs(Rx);

         //
         if(directional_ratios.empty())
         {
            // first angle entry
            //directional_ratios.push_back(ratio);
            directional_ratios.resize(1);
            directional_ratios[0] = ratio;
         }
         else
         {
            // other angle entries: check if the same angle is already present
            if(std::find(directional_ratios.begin(), directional_ratios.end(), ratio) != directional_ratios.end())
            {
               // angle already present --> skip
               continue;
            }
            else
            {
               // angle not present --> insert
               directional_ratios.push_back(ratio);
            }
         }
         //
      }
   }

   //
   return directional_ratios;
}
*/



/*
//-------------------------------------------------------------------------OK

std::vector<std::vector<double>> gen_Hamiltonian_diag_OLD()
{
   //
   std::vector<std::vector<double>> Diag_H( dimHilb , std::vector<double> (dimHilb, 0.0));

   // Diagonal Ising-like contribution
   for (int istate = 0; istate < dimHilb; istate++)
   {
      std::vector<int> bra = dec2bin(istate,dimLatt,true);
      for (int ilatt = 0; ilatt < dimLatt; ilatt++)
      {
         std::vector<int> neigh = get_neigh(ilatt, bra);
         Diag_H[istate][istate] -= 1.0 * bra[ilatt] * (neigh[0]+neigh[1]+neigh[2]+neigh[3]) * 0.5 ;
      }
   }

   //
   if( verbose ==true ) printf("gen_Hamiltonian_diag. DONE.\n");
   return Diag_H;
}
std::vector<std::vector<double>> gen_Hamiltonian(double J, double H ,std::vector<std::vector<int>> &Flip)
{
   // This routine builds the Hamiltonian given the parameter.
   // Its not very efficient as for new J and H it has to be rebuilt. Not used in the code.

   //
   std::vector<std::vector<double>> Hamiltonian( dimHilb , std::vector<double> (dimHilb, 0.0));

   // Diagonal Ising-like contribution
   for (int istate = 0; istate < dimHilb; istate++)
   {
      std::vector<int> bra = dec2bin(istate,dimLatt,true);
      for (int ilatt = 0; ilatt < dimLatt; ilatt++)
      {
         std::vector<int> neigh = get_neigh(ilatt, bra);
         Hamiltonian[istate][istate] -= J * bra[ilatt] * (neigh[0]+neigh[1]+neigh[2]+neigh[3]) * 0.5 ;
      }
   }

   // Spin-flip contribution - each ket has ilatt corresponding bras
   for (int jstate = 0; jstate < dimHilb; jstate++)
   {
      for (int ilatt = 0; ilatt < dimLatt; ilatt++)
      {
         int connected = Flip[jstate][ilatt];
         Hamiltonian[connected][jstate] -= H ;
      }
   }

   //
   return Hamiltonian;
}
*/



/*
std::vector<double> ode_integral( int tnow , double h, std::vector<std::vector<double>> alphadot)
{
   //
   int directions = alphadot[0].size();
   std::vector<double> alpha(directions,0.0);

   //
   for (int idir = 0; idir < directions; idir++)
   {
      for (int it = 0; it <= tnow; it++)
      {
         alpha[idir] -= alphadot[it][idir]*h;
      }
   }

   //
   return alpha;
}
*/



/*
std::vector<double> ode_increment( double h, std::vector<double> alphadot_old
                                           , std::vector<double> alpha_old)
{
   //
   int directions = alphadot_old.size();
   std::vector<double> alpha_new(directions,0.0);

   //
   for (int idir = 0; idir < directions; idir++)
   {
      alpha_new[idir] = -alphadot_old[idir]*h + alpha_old[idir];
   }

   //
   return alpha_new;
}
*/
