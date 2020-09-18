/***************************************************************************
* PALM++/scheduler library
*
* scheduler/convert2xml.C   convert old scheduler files to XML
*
* $Id: evaluate.C,v 1.1 2003/05/08 16:10:28 troyer Exp $
*
* Copyright (C) 2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
***************************************************************************/

#include <alps/scheduler.h>
//#include <scheduler/montecarlo.h>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <alps/alea/detailedbinning.h>
#include <cmath>

#include <iostream>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <alps/osiris/std/list.h>

typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major> complex_matrix;

// invert matrix A and calculate its determinant
template<class M> inline void invert(M & A, double & det) {
  
  M B(A.size1(), A.size1());
	
  B = boost::numeric::ublas::identity_matrix<double>(A.size1());

  boost::numeric::bindings::lapack::gesv(A, B);
  swap(A,B);
  det = 1;
  //for (int i=0; i<A.size1(); i++) {
  //  det *= B(i,i);
  //}  
  //det = std::fabs(det);

}

void evaluate(const boost::filesystem::path& p, double z, double J, int therm, int max_meas, std::string plot_type, std::ostream& out, std::string task) {

  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);

  //---------------------------------------------------------------------------
  
  if (plot_type == "ce") { // Greens function
  
	double dummy, dummyre, dummyim;
    std::vector<complex_matrix> h(4096);
  
    std::ifstream infile_h("/Users/werner/projects/diagmultiorbital/test/hamilt_small");
    for (int k=0; k<4096; k++) {
	  infile_h >> dummy;
	  h[k].resize(7,7);
	  for (int i=0; i<7; i++) {
        for (int j=0; j<7; j++) {
	      infile_h >> dummyre >> dummyim;
		  h[k](i,j) = std::complex<double>(dummyre, dummyim);
 	    }
      }  
	}
  
    double mu=0;
	double beta=5;

	std::vector<complex_matrix> gkinv(4096);
    complex_matrix g(7,7);
	
	std::vector<std::complex<double> > g0(0), g1(0), g2(0), g3(0), g4(0), g5(0), g6(0);

	for (int n=1; n<1000; n+=2) {

	for (int k=0; k<4096; k++) {
	  gkinv[k] = -h[k];
	  for (int i=0; i<7; i++) {
        gkinv[k](i,i) += std::complex<double>(mu,n*M_PI/beta);
	  }	  
	}  
  
	
	g *= 0;
	for (int k=0; k<4096; k++) {
	  complex_matrix M(gkinv[k]);
	  invert(M, dummy);
	  g += M;
	}  
	g /= 4096;
	
    for (int i=0; i<7; i++) {
      for (int j=0; j<7; j++) {
		std::cout << n << "     " << i << " " << j << "   " << g(i,j) << "\n";
	  }	          
	}	  
	std::cout << "\n";
	
	std::cout << "datag  " << n*M_PI/beta << "   " << g(0,0).real() << " " << g(0,0).imag() << "   " << g(1,1).real() << " " << g(1,1).imag() << "   " << g(2,2).real() << " " << g(2,2).imag() << "   " << g(3,3).real() << " " << g(3,3).imag() << "   " << g(4,4).real() << " " << g(4,4).imag() << "   " << g(5,5).real() << " " << g(5,5).imag() << "   " << g(6,6).real() << " " << g(6,6).imag() << "\n";
	std::cout << "\n\n";
	    
	g0.push_back(g(0,0));	
	g1.push_back(g(1,1));
	g2.push_back(g(2,2));
	g3.push_back(g(3,3));
	g4.push_back(g(4,4));
	g5.push_back(g(5,5));
	g6.push_back(g(6,6));
		
	}

	std::vector<double> e1(7), e2(7), e3(7);
	
	for (int k=0; k<4096; k++) {
	  for (int i=0; i<7; i++) {
	    e1[i] += h[k](i,i).real();
		for (int j=0; j<7; j++) {
		  e2[i] += h[k](i,j).real()*h[k](j,i).real();
		}
		for (int j=0; j<7; j++) {
		for (int l=0; l<7; l++) {
		  e3[i] += h[k](i,j).real()*h[k](j,l).real()*h[k](l,i).real();
        }
		}	  
	  }
    }
	for (int i=0; i<7; i++) {
	  e1[i] /= 4096;
	  e2[i] /= 4096;
	  e3[i] /= 4096;
	  std::cout << "epsilon " << i << "   " << e1[i] << " " << e2[i] << " " << e3[i] << "   " << e2[i]-e1[i]*e1[i] << "   " << (e3[i]-2*e1[i]*e2[i]+std::pow(e1[i],3))/(e2[i]-e1[i]*e1[i]) << "\n";
	}
	
	std::vector<std::complex<double> > ftrunc0p(g0), ftrunc1p(g1), ftrunc2p(g2), ftrunc3p(g3), ftrunc4p(g4), ftrunc5p(g5), ftrunc6p(g6);
	std::vector<std::complex<double> > gtrunc0p(g0), gtrunc1p(g1), gtrunc2p(g2), gtrunc3p(g3), gtrunc4p(g4), gtrunc5p(g5), gtrunc6p(g6);	
	
	double s0=0;
	
	for (int i=0; i<ftrunc0p.size(); i++) {
	  ftrunc0p[i] = -1./ftrunc0p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[0]-(e2[0]-e1[0]*e1[0])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[0]-e1[0]*e1[0])*(s0-(mu-(e3[0]-2*e1[0]*e2[0]+std::pow(e1[0],3))/(e2[0]-e1[0]*e1[0])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc1p[i] = -1./ftrunc1p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[1]-(e2[1]-e1[1]*e1[1])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[1]-e1[1]*e1[1])*(s0-(mu-(e3[1]-2*e1[1]*e2[1]+std::pow(e1[1],3))/(e2[1]-e1[1]*e1[1])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc2p[i] = -1./ftrunc2p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[2]-(e2[2]-e1[2]*e1[2])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[2]-e1[2]*e1[2])*(s0-(mu-(e3[2]-2*e1[2]*e2[2]+std::pow(e1[2],3))/(e2[2]-e1[2]*e1[2])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc3p[i] = -1./ftrunc3p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[3]-(e2[3]-e1[3]*e1[3])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[3]-e1[3]*e1[3])*(s0-(mu-(e3[3]-2*e1[3]*e2[3]+std::pow(e1[3],3))/(e2[3]-e1[3]*e1[3])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc4p[i] = -1./ftrunc4p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[4]-(e2[4]-e1[4]*e1[4])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[4]-e1[4]*e1[4])*(s0-(mu-(e3[4]-2*e1[4]*e2[4]+std::pow(e1[4],3))/(e2[4]-e1[4]*e1[4])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc5p[i] = -1./ftrunc5p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[5]-(e2[5]-e1[5]*e1[5])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[5]-e1[5]*e1[5])*(s0-(mu-(e3[5]-2*e1[5]*e2[5]+std::pow(e1[5],3))/(e2[5]-e1[5]*e1[5])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc6p[i] = -1./ftrunc6p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[6]-(e2[6]-e1[6]*e1[6])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[6]-e1[6]*e1[6])*(s0-(mu-(e3[6]-2*e1[6]*e2[6]+std::pow(e1[6],3))/(e2[6]-e1[6]*e1[6])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	
	  gtrunc0p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[0]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  gtrunc1p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[1]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
      gtrunc2p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[2]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
      gtrunc3p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[3]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
      gtrunc4p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[4]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  gtrunc5p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[5]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  gtrunc6p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[6]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);

	  //std::cout << "dataf  " << (2*i+1)*M_PI/beta << "   " << ftrunc0[i].real() << " " << ftrunc0[i].imag() << "   " << ftrunc1[i].real() << " " << ftrunc1[i].imag() << "   " << ftrunc2[i].real() << " " << ftrunc2[i].imag() << "   " << ftrunc3[i].real() << " " << ftrunc3[i].imag() << "   " << ftrunc4[i].real() << " " << ftrunc4[i].imag() << "   " << ftrunc5[i].real() << " " << ftrunc5[i].imag() << "   " << ftrunc6[i].real() << " " << ftrunc6[i].imag() << "\n";
	  //std::cout << "\n\n";	
	
	}
	
	std::vector<double> f0t(101), f1t(101), f2t(101), f3t(101), f4t(101), f5t(101), f6t(101);
	std::vector<double> g0t(101), g1t(101), g2t(101), g3t(101), g4t(101), g5t(101), g6t(101);	
	double t_res = 100;
	
	for (int k=0; k<=t_res; k++) {
	  for (int n=0; n<ftrunc0p.size(); n++) {
		f0t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc0p[n]).real();
		f1t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc1p[n]).real();
		f2t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc2p[n]).real();
		f3t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc3p[n]).real();
		f4t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc4p[n]).real();
		f5t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc5p[n]).real();
		f6t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc6p[n]).real();
		
		g0t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc0p[n]).real();
		g1t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc1p[n]).real();
		g2t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc2p[n]).real();
		g3t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc3p[n]).real();
		g4t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc4p[n]).real();
		g5t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc5p[n]).real();
		g6t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc6p[n]).real();		
	  }
	  f0t[k] += -0.5*(e2[0]-e1[0]*e1[0])-0.25*(beta-2*beta*(k/t_res))*(e2[0]-e1[0]*e1[0])*(s0-(mu-(e3[0]-2*e1[0]*e2[0]+std::pow(e1[0],3))/(e2[0]-e1[0]*e1[0])));
	  f1t[k] += -0.5*(e2[1]-e1[1]*e1[1])-0.25*(beta-2*beta*(k/t_res))*(e2[1]-e1[1]*e1[1])*(s0-(mu-(e3[1]-2*e1[1]*e2[1]+std::pow(e1[1],3))/(e2[1]-e1[1]*e1[1])));
	  f2t[k] += -0.5*(e2[2]-e1[2]*e1[2])-0.25*(beta-2*beta*(k/t_res))*(e2[2]-e1[2]*e1[2])*(s0-(mu-(e3[2]-2*e1[2]*e2[2]+std::pow(e1[2],3))/(e2[2]-e1[2]*e1[2])));
	  f3t[k] += -0.5*(e2[3]-e1[3]*e1[3])-0.25*(beta-2*beta*(k/t_res))*(e2[3]-e1[3]*e1[3])*(s0-(mu-(e3[3]-2*e1[3]*e2[3]+std::pow(e1[3],3))/(e2[3]-e1[3]*e1[3])));
	  f4t[k] += -0.5*(e2[4]-e1[4]*e1[4])-0.25*(beta-2*beta*(k/t_res))*(e2[4]-e1[4]*e1[4])*(s0-(mu-(e3[4]-2*e1[4]*e2[4]+std::pow(e1[4],3))/(e2[4]-e1[4]*e1[4])));
      f5t[k] += -0.5*(e2[5]-e1[5]*e1[5])-0.25*(beta-2*beta*(k/t_res))*(e2[5]-e1[5]*e1[5])*(s0-(mu-(e3[5]-2*e1[5]*e2[5]+std::pow(e1[5],3))/(e2[5]-e1[5]*e1[5])));
	  f6t[k] += -0.5*(e2[6]-e1[6]*e1[6])-0.25*(beta-2*beta*(k/t_res))*(e2[6]-e1[6]*e1[6])*(s0-(mu-(e3[6]-2*e1[6]*e2[6]+std::pow(e1[6],3))/(e2[6]-e1[6]*e1[6])));
	
	  g0t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[0]+s0-mu);
	  g1t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[1]+s0-mu);
	  g2t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[2]+s0-mu);
	  g3t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[3]+s0-mu);
	  g4t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[4]+s0-mu);
	  g5t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[5]+s0-mu);
	  g6t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[6]+s0-mu);
	}
	
	for (int k=0; k<=t_res; k++) {	
	  std::cout << "f(t) " << k/t_res << " " << beta*k/t_res << " " << -f0t[t_res-k] << " " << -f1t[t_res-k] << " " << -f2t[t_res-k] << " " << -f3t[t_res-k] << " " << -f4t[t_res-k] << " " << -f5t[t_res-k] << " " << -f6t[t_res-k] << "\n";  
	  std::cout << "-g(t) " << k/t_res << " " << beta*k/t_res << " " << -g0t[k] << " " << -g1t[k] << " " << -g2t[k] << " " << -g3t[k] << " " << -g4t[k] << " " << -g5t[k] << " " << -g6t[k] << "\n";  	
	}
	    
    /*
	
	alps::RealVectorObsevaluator G=sim.get_measurements()["Greens"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);
	
	alps::RealVectorObsevaluator n=sim.get_measurements()["n"];
    
	for (int j=0; j<FLAVORS/2; j++) {
	  std::cout << "# orbital " << j << ", spins up and down\n\n";
	  
	  std::valarray<double> Gsym(N+1);
	  std::valarray<double> Gsym_error(N+1);
	  for (int i=0; i<N+1; i++) {
	    Gsym[i] = (G.mean()[2*j*(N+1)+i]+G.mean()[(2*j+1)*(N+1)+i])/2;
	    Gsym_error[i] = (G.error()[2*j*(N+1)+i]+G.error()[(2*j+1)*(N+1)+i])/(2*sqrt(2.));
	  }
	
	  printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n", 0, (BETA/N)*0., 1-n.mean()[2*j], n.error()[2*j], 1-n.mean()[2*j+1], n.error()[2*j+1], 0.5, 0.);
	
      for (int i=1; i<N; i++)
        printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n", i, (BETA/N)*i, G.mean()[2*j*(N+1)+i], G.error()[2*j*(N+1)+i], G.mean()[(2*j+1)*(N+1)+i], G.error()[(2*j+1)*(N+1)+i], Gsym[i], Gsym_error[i]);

      printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n\n\n", N, (BETA/N)*N, n.mean()[2*j], n.error()[2*j], n.mean()[2*j+1], n.error()[2*j+1], 0.5, 0.);
    
	}
	
	*/
  }
  
  //---------------------------------------------------------------------------
  
  if (plot_type == "ceself") { // Greens function

    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
    int Np1 = N+1;
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);

    double mu=0;
	double beta=5;

    std::vector<std::vector<double> > F;
    F.resize(FLAVORS);
	for (int i=0; i<FLAVORS; i++) {
	  F[i].resize(Np1);
    }
  
    // read F from file
    std::ifstream infile_f(boost::lexical_cast<std::string>(sim.get_parameters()["G"]).c_str());
    for (int i=0; i<Np1; i++) {
      double dummy;
      infile_f >> dummy; 
	  for (int j=0; j<FLAVORS; j++)
	    infile_f >> F[j][i];
    }
	
    std::vector<std::vector<double> > G;
    G.resize(FLAVORS);
	for (int i=0; i<FLAVORS; i++) {
	  G[i].resize(Np1);
    }
	  
    // read F from file
    //ifstream infile_f(boost::lexical_cast<std::string>(parms["G"]).c_str());
	std::ifstream infile_g("G");
    for (int i=0; i<Np1; i++) {
      double dummy;
      infile_g >> dummy; 
	  for (int j=0; j<FLAVORS; j++)
	    infile_g >> G[j][i];
    }
	
    std::vector<std::vector<std::complex<double> > > Fw(FLAVORS);
	for (int i=0; i<FLAVORS; i++) {
	  Fw[i].resize(500);
	  for (int n=0; n<500; n++) {
	    Fw[i][n] = 0.5*beta/N*(F[i][0]-F[i][N]);
		for (int k=1; k<N; k++)
		  Fw[i][n] += beta/N*std::complex<double>(cos((2*n+1)*M_PI*k/N), sin((2*n+1)*M_PI*k/N))*F[i][k]; 
	  }
    }	
	
    std::vector<std::vector<std::complex<double> > > Gw(FLAVORS);
	for (int i=0; i<FLAVORS; i++) {
	  Gw[i].resize(500);
	  for (int n=0; n<500; n++) {
	    Gw[i][n] = 0.5*beta/N*(G[i][0]-G[i][N]);
		for (int k=1; k<N; k++)
		  Gw[i][n] += beta/N*std::complex<double>(cos((2*n+1)*M_PI*k/N), sin((2*n+1)*M_PI*k/N))*G[i][k]; 
	  }
    }		
  
	double dummy, dummyre, dummyim;
    std::vector<complex_matrix> h(4096);
  
    std::ifstream infile_h("/Users/werner/projects/diagmultiorbital/test/hamilt_small");
    for (int k=0; k<4096; k++) {
	  infile_h >> dummy;
	  h[k].resize(7,7);
	  for (int i=0; i<7; i++) {
        for (int j=0; j<7; j++) {
	      infile_h >> dummyre >> dummyim;
		  h[k](i,j) = std::complex<double>(dummyre, dummyim);
 	    }
      }  
	}
  
	std::vector<complex_matrix> gkinv(4096);
    complex_matrix g(7,7);
	
	std::vector<std::complex<double> > g0(0), g1(0), g2(0), g3(0), g4(0), g5(0), g6(0);

	for (int n=1; n<1000; n+=2) {

	for (int k=0; k<4096; k++) {
	  gkinv[k] = -h[k];
	  for (int i=0; i<7; i++) {
        gkinv[k](i,i) += std::complex<double>(mu,n*M_PI/beta);
	  }	  
	}  
  
	
	g *= 0;
	for (int k=0; k<4096; k++) {
	  complex_matrix M(gkinv[k]);
	  invert(M, dummy);
	  g += M;
	}  
	g /= 4096;
	
    for (int i=0; i<7; i++) {
      for (int j=0; j<7; j++) {
		std::cout << n << "     " << i << " " << j << "   " << g(i,j) << "\n";
	  }	          
	}	  
	std::cout << "\n";
	
	std::cout << "datag  " << n*M_PI/beta << "   " << g(0,0).real() << " " << g(0,0).imag() << "   " << g(1,1).real() << " " << g(1,1).imag() << "   " << g(2,2).real() << " " << g(2,2).imag() << "   " << g(3,3).real() << " " << g(3,3).imag() << "   " << g(4,4).real() << " " << g(4,4).imag() << "   " << g(5,5).real() << " " << g(5,5).imag() << "   " << g(6,6).real() << " " << g(6,6).imag() << "\n";
	std::cout << "\n\n";
	    
	g0.push_back(g(0,0));	
	g1.push_back(g(1,1));
	g2.push_back(g(2,2));
	g3.push_back(g(3,3));
	g4.push_back(g(4,4));
	g5.push_back(g(5,5));
	g6.push_back(g(6,6));
		
	}

	std::vector<double> e1(7), e2(7), e3(7);
	
	for (int k=0; k<4096; k++) {
	  for (int i=0; i<7; i++) {
	    e1[i] += h[k](i,i).real();
		for (int j=0; j<7; j++) {
		  e2[i] += h[k](i,j).real()*h[k](j,i).real();
		}
		for (int j=0; j<7; j++) {
		for (int l=0; l<7; l++) {
		  e3[i] += h[k](i,j).real()*h[k](j,l).real()*h[k](l,i).real();
        }
		}	  
	  }
    }
	for (int i=0; i<7; i++) {
	  e1[i] /= 4096;
	  e2[i] /= 4096;
	  e3[i] /= 4096;
	  std::cout << "epsilon " << i << "   " << e1[i] << " " << e2[i] << " " << e3[i] << "   " << e2[i]-e1[i]*e1[i] << "   " << (e3[i]-2*e1[i]*e2[i]+std::pow(e1[i],3))/(e2[i]-e1[i]*e1[i]) << "\n";
	}
	
	std::vector<std::vector<std::complex<double> > > Sw(FLAVORS);
	for (int i=0; i<FLAVORS; i++) {
	  Sw[i].resize(500);
	  for (int n=0; n<500; n++) {
	    //Sw[i][n] = -conj(Fw[i][n])+std::complex<double>(mu-e1[i], (2*n+1)*M_PI/beta)-1./Gw[i][n];
		Sw[i][n] = 1./Gw[i][n];
	  }
    }	
	
    for (int i=0; i<500; i++) {
	  std::cout << "datas  " << (2*i+1)*M_PI/beta << "   " << Sw[0][i].real() << " " << Sw[0][i].imag() << "   " << Sw[1][i].real() << " " << Sw[1][i].imag() << "   " << Sw[2][i].real() << " " << Sw[2][i].imag() << "   " << Sw[3][i].real() << " " << Sw[3][i].imag() << "   " << Sw[4][i].real() << " " << Sw[4][i].imag() << "   " << Sw[5][i].real() << " " << Sw[5][i].imag() << "   " << Sw[6][i].real() << " " << Sw[6][i].imag() << "\n";
	  std::cout << "\n\n";	
	}
	

	
	std::vector<std::complex<double> > ftrunc0p(g0), ftrunc1p(g1), ftrunc2p(g2), ftrunc3p(g3), ftrunc4p(g4), ftrunc5p(g5), ftrunc6p(g6);
	std::vector<std::complex<double> > gtrunc0p(g0), gtrunc1p(g1), gtrunc2p(g2), gtrunc3p(g3), gtrunc4p(g4), gtrunc5p(g5), gtrunc6p(g6);	
	
	double s0=0;
	
	for (int i=0; i<ftrunc0p.size(); i++) {
	  ftrunc0p[i] = -1./ftrunc0p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[0]-(e2[0]-e1[0]*e1[0])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[0]-e1[0]*e1[0])*(s0-(mu-(e3[0]-2*e1[0]*e2[0]+std::pow(e1[0],3))/(e2[0]-e1[0]*e1[0])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc1p[i] = -1./ftrunc1p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[1]-(e2[1]-e1[1]*e1[1])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[1]-e1[1]*e1[1])*(s0-(mu-(e3[1]-2*e1[1]*e2[1]+std::pow(e1[1],3))/(e2[1]-e1[1]*e1[1])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc2p[i] = -1./ftrunc2p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[2]-(e2[2]-e1[2]*e1[2])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[2]-e1[2]*e1[2])*(s0-(mu-(e3[2]-2*e1[2]*e2[2]+std::pow(e1[2],3))/(e2[2]-e1[2]*e1[2])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc3p[i] = -1./ftrunc3p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[3]-(e2[3]-e1[3]*e1[3])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[3]-e1[3]*e1[3])*(s0-(mu-(e3[3]-2*e1[3]*e2[3]+std::pow(e1[3],3))/(e2[3]-e1[3]*e1[3])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc4p[i] = -1./ftrunc4p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[4]-(e2[4]-e1[4]*e1[4])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[4]-e1[4]*e1[4])*(s0-(mu-(e3[4]-2*e1[4]*e2[4]+std::pow(e1[4],3))/(e2[4]-e1[4]*e1[4])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc5p[i] = -1./ftrunc5p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[5]-(e2[5]-e1[5]*e1[5])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[5]-e1[5]*e1[5])*(s0-(mu-(e3[5]-2*e1[5]*e2[5]+std::pow(e1[5],3))/(e2[5]-e1[5]*e1[5])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  ftrunc6p[i] = -1./ftrunc6p[i]+std::complex<double>(mu, (2*i+1)*M_PI/beta)-e1[6]-(e2[6]-e1[6]*e1[6])/std::complex<double>(0, (2*i+1)*M_PI/beta)-(e2[6]-e1[6]*e1[6])*(s0-(mu-(e3[6]-2*e1[6]*e2[6]+std::pow(e1[6],3))/(e2[6]-e1[6]*e1[6])))/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	
	  gtrunc0p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[0]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  gtrunc1p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[1]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
      gtrunc2p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[2]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
      gtrunc3p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[3]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
      gtrunc4p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[4]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  gtrunc5p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[5]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);
	  gtrunc6p[i] -= 1./std::complex<double>(0, (2*i+1)*M_PI/beta)+(e1[6]+s0-mu)/std::pow(std::complex<double>(0, (2*i+1)*M_PI/beta),2);

	  //std::cout << "dataf  " << (2*i+1)*M_PI/beta << "   " << ftrunc0[i].real() << " " << ftrunc0[i].imag() << "   " << ftrunc1[i].real() << " " << ftrunc1[i].imag() << "   " << ftrunc2[i].real() << " " << ftrunc2[i].imag() << "   " << ftrunc3[i].real() << " " << ftrunc3[i].imag() << "   " << ftrunc4[i].real() << " " << ftrunc4[i].imag() << "   " << ftrunc5[i].real() << " " << ftrunc5[i].imag() << "   " << ftrunc6[i].real() << " " << ftrunc6[i].imag() << "\n";
	  //std::cout << "\n\n";	
	
	}
	
	std::vector<double> f0t(101), f1t(101), f2t(101), f3t(101), f4t(101), f5t(101), f6t(101);
	std::vector<double> g0t(101), g1t(101), g2t(101), g3t(101), g4t(101), g5t(101), g6t(101);	
	double t_res = 100;
	
	for (int k=0; k<=t_res; k++) {
	  for (int n=0; n<ftrunc0p.size(); n++) {
		f0t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc0p[n]).real();
		f1t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc1p[n]).real();
		f2t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc2p[n]).real();
		f3t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc3p[n]).real();
		f4t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc4p[n]).real();
		f5t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc5p[n]).real();
		f6t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*ftrunc6p[n]).real();
		
		g0t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc0p[n]).real();
		g1t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc1p[n]).real();
		g2t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc2p[n]).real();
		g3t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc3p[n]).real();
		g4t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc4p[n]).real();
		g5t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc5p[n]).real();
		g6t[k] += 2*(1/beta*std::complex<double>(cos((2*n+1)*M_PI*k/t_res), -sin((2*n+1)*M_PI*k/t_res))*gtrunc6p[n]).real();		
	  }
	  f0t[k] += -0.5*(e2[0]-e1[0]*e1[0])-0.25*(beta-2*beta*(k/t_res))*(e2[0]-e1[0]*e1[0])*(s0-(mu-(e3[0]-2*e1[0]*e2[0]+std::pow(e1[0],3))/(e2[0]-e1[0]*e1[0])));
	  f1t[k] += -0.5*(e2[1]-e1[1]*e1[1])-0.25*(beta-2*beta*(k/t_res))*(e2[1]-e1[1]*e1[1])*(s0-(mu-(e3[1]-2*e1[1]*e2[1]+std::pow(e1[1],3))/(e2[1]-e1[1]*e1[1])));
	  f2t[k] += -0.5*(e2[2]-e1[2]*e1[2])-0.25*(beta-2*beta*(k/t_res))*(e2[2]-e1[2]*e1[2])*(s0-(mu-(e3[2]-2*e1[2]*e2[2]+std::pow(e1[2],3))/(e2[2]-e1[2]*e1[2])));
	  f3t[k] += -0.5*(e2[3]-e1[3]*e1[3])-0.25*(beta-2*beta*(k/t_res))*(e2[3]-e1[3]*e1[3])*(s0-(mu-(e3[3]-2*e1[3]*e2[3]+std::pow(e1[3],3))/(e2[3]-e1[3]*e1[3])));
	  f4t[k] += -0.5*(e2[4]-e1[4]*e1[4])-0.25*(beta-2*beta*(k/t_res))*(e2[4]-e1[4]*e1[4])*(s0-(mu-(e3[4]-2*e1[4]*e2[4]+std::pow(e1[4],3))/(e2[4]-e1[4]*e1[4])));
      f5t[k] += -0.5*(e2[5]-e1[5]*e1[5])-0.25*(beta-2*beta*(k/t_res))*(e2[5]-e1[5]*e1[5])*(s0-(mu-(e3[5]-2*e1[5]*e2[5]+std::pow(e1[5],3))/(e2[5]-e1[5]*e1[5])));
	  f6t[k] += -0.5*(e2[6]-e1[6]*e1[6])-0.25*(beta-2*beta*(k/t_res))*(e2[6]-e1[6]*e1[6])*(s0-(mu-(e3[6]-2*e1[6]*e2[6]+std::pow(e1[6],3))/(e2[6]-e1[6]*e1[6])));
	
	  g0t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[0]+s0-mu);
	  g1t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[1]+s0-mu);
	  g2t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[2]+s0-mu);
	  g3t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[3]+s0-mu);
	  g4t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[4]+s0-mu);
	  g5t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[5]+s0-mu);
	  g6t[k] += -0.5-0.25*(beta-2*beta*(k/t_res))*(e1[6]+s0-mu);
	}
	
	for (int k=0; k<=t_res; k++) {	
	  std::cout << "f(t) " << k/t_res << " " << beta*k/t_res << " " << -f0t[t_res-k] << " " << -f1t[t_res-k] << " " << -f2t[t_res-k] << " " << -f3t[t_res-k] << " " << -f4t[t_res-k] << " " << -f5t[t_res-k] << " " << -f6t[t_res-k] << "\n";  
	  std::cout << "-g(t) " << k/t_res << " " << beta*k/t_res << " " << -g0t[k] << " " << -g1t[k] << " " << -g2t[k] << " " << -g3t[k] << " " << -g4t[k] << " " << -g5t[k] << " " << -g6t[k] << "\n";  	
	}
	    
    /*
	
	alps::RealVectorObsevaluator G=sim.get_measurements()["Greens"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);
	
	alps::RealVectorObsevaluator n=sim.get_measurements()["n"];
    
	for (int j=0; j<FLAVORS/2; j++) {
	  std::cout << "# orbital " << j << ", spins up and down\n\n";
	  
	  std::valarray<double> Gsym(N+1);
	  std::valarray<double> Gsym_error(N+1);
	  for (int i=0; i<N+1; i++) {
	    Gsym[i] = (G.mean()[2*j*(N+1)+i]+G.mean()[(2*j+1)*(N+1)+i])/2;
	    Gsym_error[i] = (G.error()[2*j*(N+1)+i]+G.error()[(2*j+1)*(N+1)+i])/(2*sqrt(2.));
	  }
	
	  printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n", 0, (BETA/N)*0., 1-n.mean()[2*j], n.error()[2*j], 1-n.mean()[2*j+1], n.error()[2*j+1], 0.5, 0.);
	
      for (int i=1; i<N; i++)
        printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n", i, (BETA/N)*i, G.mean()[2*j*(N+1)+i], G.error()[2*j*(N+1)+i], G.mean()[(2*j+1)*(N+1)+i], G.error()[(2*j+1)*(N+1)+i], Gsym[i], Gsym_error[i]);

      printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n\n\n", N, (BETA/N)*N, n.mean()[2*j], n.error()[2*j], n.mean()[2*j+1], n.error()[2*j+1], 0.5, 0.);
    
	}
	
	*/
  }
   
  //---------------------------------------------------------------------------
  
  if (plot_type == "Gce") { // Greens function
  
	alps::RealVectorObsevaluator G=sim.get_measurements()["Greens"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);
	
	alps::RealVectorObsevaluator n=sim.get_measurements()["n"];
    
	for (int j=0; j<FLAVORS; j++) {
	  std::cout << "# orbital " << j << "\n\n";
	
	  printf("index: %4d  tau: %9.5f  G: %9.7e +- %9.7e\n", 0, (BETA/N)*0., 1-n.mean()[j], n.error()[j]);
	
      for (int i=1; i<N; i++)
        printf("index: %4d  tau: %9.5f  G: %9.7e +- %9.7e\n", i, (BETA/N)*i, G.mean()[j*(N+1)+i], G.error()[j*(N+1)+i]);

      printf("index: %4d  tau: %9.5f  G: %9.7e +- %9.7e\n\n\n", N, (BETA/N)*N, n.mean()[j], n.error()[j]);
    
	}
  }
  
  //---------------------------------------------------------------------------
  
  if (plot_type == "Gnoerror") { // Greens function
  
	alps::RealVectorObsevaluator G=sim.get_measurements()["Greens"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);
	
	alps::RealVectorObsevaluator n=sim.get_measurements()["n"];
    
	printf("%9.5f  ", 0.);
	for (int j=0; j<FLAVORS; j++) 
	  printf("%9.7e  ", -(1-n.mean()[j]));
	std::cout << "\n";
	
	for (int i=1; i<N; i++) {
  	  printf("%9.5f  ", (BETA/N)*i);		
	  for (int j=0; j<FLAVORS; j++) 	  
	    printf("%9.7e  ", -G.mean()[j*(N+1)+i]);
	  std::cout << "\n";		
	}

	printf("%9.5f  ", 1.);
	for (int j=0; j<FLAVORS; j++) 
	  printf("%9.7e  ", -n.mean()[j]);
	std::cout << "\n";
  }

  //---------------------------------------------------------------------------
  
  if (plot_type == "G") { // Greens function
  
	alps::RealVectorObsevaluator G=sim.get_measurements()["Greens"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);
	
	alps::RealVectorObsevaluator n=sim.get_measurements()["n"];
    
	for (int j=0; j<FLAVORS/2; j++) {
	  std::cout << "# orbital " << j << ", spins up and down\n\n";
	  
	  std::valarray<double> Gsym(N+1);
	  std::valarray<double> Gsym_error(N+1);
	  for (int i=0; i<N+1; i++) {
	    Gsym[i] = (G.mean()[2*j*(N+1)+i]+G.mean()[(2*j+1)*(N+1)+i])/2;
	    Gsym_error[i] = (G.error()[2*j*(N+1)+i]+G.error()[(2*j+1)*(N+1)+i])/(2*sqrt(2.));
	  }
	
	  printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n", 0, (BETA/N)*0., 1-n.mean()[2*j], n.error()[2*j], 1-n.mean()[2*j+1], n.error()[2*j+1], 0.5, 0.);
	
      for (int i=1; i<N; i++)
        printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n", i, (BETA/N)*i, G.mean()[2*j*(N+1)+i], G.error()[2*j*(N+1)+i], G.mean()[(2*j+1)*(N+1)+i], G.error()[(2*j+1)*(N+1)+i], Gsym[i], Gsym_error[i]);

      printf("index: %4d  tau: %9.5f  G_up: %9.7e +- %9.7e  G_down: %9.7e +- %9.7e  Gsym: %9.7e +- %9.7e\n\n\n", N, (BETA/N)*N, n.mean()[2*j], n.error()[2*j], n.mean()[2*j+1], n.error()[2*j+1], 0.5, 0.);
    
	}
  }

  //---------------------------------------------------------------------------
  
  if (plot_type == "nn") { // Greens function
  
	alps::RealVectorObsevaluator nn_corr=sim.get_measurements()["nn_corr"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N_corr = static_cast<int>(sim.get_parameters()["N_CORR"]);
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);
	
	std::valarray<double> nn_corr_meas(FLAVORS*(FLAVORS+1)/2*(N_corr+1));
    int position=0;
    for (int flavor1=0; flavor1<FLAVORS; ++flavor1) {
      for (int flavor2=0; flavor2<=flavor1; ++flavor2) {
  
        std::cout << "#flavors " << flavor1 << " and " << flavor2 << "\n\n";
  
		for (int index=0; index<N_corr+1; ++index) 
		  printf("tau: %9.5f  nn: %9.7e +- %9.7e\n", index*BETA/N_corr, nn_corr.mean()[position+index], nn_corr.error()[position+index]);
	
	    std::cout << "\n\n";
	    position += (N_corr+1);
      }
    }
 	
  }
  
  //-------------------------------------------------------------------------------------------
 
  if (plot_type == "selfsym2a") { // Greens function

	alps::RealVectorObsevaluator G=sim.get_measurements()["Greens"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	int FLAVORS = static_cast<int>(sim.get_parameters()["FLAVORS"]);
	
	//alps::RealVectorObsevaluator n=sim.get_measurements()["n"];
    
	std::valarray<double> Ga(N+1), Gb(N+1);
	std::valarray<double> Ga_error(N+1), Gb_error(N+1);
	for (int i=0; i<N+1; i++) {
	  Ga[i] = (G.mean()[i]+G.mean()[(N+1)+i])/2;
	  Ga_error[i] = (G.error()[i]+G.error()[(N+1)+i])/(2*sqrt(2.));
	  Gb[i] = (G.mean()[2*(N+1)+i]+G.mean()[3*(N+1)+i])/2;
	  Gb_error[i] = (G.error()[2*(N+1)+i]+G.error()[3*(N+1)+i])/(2*sqrt(2.));
	}
	
	std::valarray<double> Ga_sym(N+1), Gb_sym(N+1);
	for (int i=0; i<N+1; i++) {
	  Ga_sym[i] = (Ga[i]/Ga_error[i]+Ga[N-i]/Ga_error[N-i])/(1/Ga_error[i]+1/Ga_error[N-i]);
	  Gb_sym[i] = (Gb[i]/Gb_error[i]+Gb[N-i]/Gb_error[N-i])/(1/Gb_error[i]+1/Gb_error[N-i]);
	}
	//(G_up.mean()[i]/G_up.error()[i]+G_up.mean()[N-i]/G_up.error()[N-i]+G_down.mean()[i]/G_down.error()[i]+G_down.mean()[N-i]/G_down.error()[N-i])/(1/G_up.error()[i]+1/G_up.error()[N-i]+1/G_down.error()[i]+1/G_down.error()[N-i]);
	
	std::valarray<double> Gsym(Ga_sym), Gcoarse(Ga_sym), Ga_coarse(Ga_sym), Gb_coarse(Gb_sym);
	
	for (int orbital=0; orbital<2; orbital++) {
 
	  if (orbital==0) {
	    Gcoarse=Ga_coarse;
		Gsym=Ga_sym;
	  }
	  else {
	    Gcoarse=Gb_coarse;
	  	Gsym=Gb_sym;
	  }
		
	int k_max=J;
	
	for (int k=1; k<=k_max; k++) {
	for (int i=N/BETA*z*(1+(k-1)*0.2); i<N/BETA*z*(1.+k*0.2); i++) {
	  double coarse=0;
	  for (int n=-k; n<=k; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k+1);
	}

	for (int i=N/BETA*(1-z*(1+k*0.2)); i<N/BETA*(1-z*(1.+(k-1)*0.2)); i++) {
	  double coarse=0;
	  for (int n=-k; n<=k; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k+1);
	}
	}
	
	for (int i=N/BETA*z*(1+k_max*0.2); i<N/BETA*(1-z*(1.+k_max*0.2)); i++) {
	  double coarse=0;
	  for (int n=-k_max; n<=k_max; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k_max+1);
	}
	
	  if (orbital==0) Ga_coarse=Gcoarse;
	  else Gb_coarse=Gcoarse;

	}
	
	printf("%9.5f  %10.8e\n", (BETA/N)*0, -0.5);
	
    for (int i=1; i<N; i++)
      printf("%9.5f  %10.8e\n", (BETA/N)*i, 0.5*(-Ga_coarse[i]-Gb_coarse[i]));
	
	printf("%9.5f  %10.8e\n", (BETA/N)*N, -0.5);
  }


/*
//---------------------------------------------------------------------------
  
  if (plot_type == "total") {
  
	alps::RealVectorObsevaluator t_up=sim.get_measurements()["total_up"];
	alps::RealVectorObsevaluator t_down=sim.get_measurements()["total_down"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);

    for (int i=0; i<N+1; i++)
      printf("index: %4d  tau: %9.5f  t_up: %9.7e +- %9.7e  t_down: %9.7e +- %9.7e  tsym: %9.7e +- %9.7e\n", i, (BETA/N)*i, t_up.mean()[i], t_up.error()[i], t_down.mean()[i], t_down.error()[i], .5*(t_up.mean()[i]+t_down.mean()[i]), .5*(t_up.error()[i]+t_down.error()[i]));
  }
 
//-------------------------------------------------------------------------------------------
 
  if (plot_type == "selfsym") { // Greens function
  
	alps::RealVectorObsevaluator G_up=sim.get_measurements()["Greens_up"];
	alps::RealVectorObsevaluator G_down=sim.get_measurements()["Greens_down"];
	//alps::RealObsevaluator sign=sim.get_measurements()["sign"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
	//double const2 = static_cast<double>(sim.get_parameters()["CONST2"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
    
	//alps::RealVectorObsevaluator Gnorm = G; // /sign;
	
	std::valarray<double> Gsym(N+1);
	for (int i=0; i<N+1; i++)
	  Gsym[i] = (G_up.mean()[i]+G_up.mean()[N-i]+G_down.mean()[i]+G_down.mean()[N-i])/4;
	
	std::valarray<double> Gcoarse(Gsym);
	
	for (int i=N/BETA*z*1; i<N/BETA*z*1.2; i++) {
	  double coarse = (Gsym[i-1]+Gsym[i]+Gsym[i+1])/3.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.2; i<N/BETA*z*1.4; i++) {
	  double coarse = (Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2])/5.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.4; i<N/BETA*z*1.6; i++) {
	  double coarse = (Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3])/7.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.6; i<N/BETA*z*1.8; i++) {
	  double coarse = (Gsym[i-4]+Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3]+Gsym[i+4])/9.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.8; i<N/BETA*(1-z*1.8); i++) {
	  double coarse = (Gsym[i-5]+Gsym[i-4]+Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3]+Gsym[i+4]+Gsym[i+5])/11.;
	  Gcoarse[i] = coarse;
	}
	

	for (int i=N/BETA*(1-z*1.2); i<N/BETA*(1-z*1); i++) {
	  double coarse = (Gsym[i-1]+Gsym[i]+Gsym[i+1])/3.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*(1-z*1.4); i<N/BETA*(1-z*1.2); i++) {
	  double coarse = (Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2])/5.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*(1-z*1.6); i<N/BETA*(1-z*1.4); i++) {
	  double coarse = (Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3])/7.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*(1-z*1.8); i<N/BETA*(1-z*1.6); i++) {
	  double coarse = (Gsym[i-4]+Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3]+Gsym[i+4])/9.;
	  Gcoarse[i] = coarse;
	}
		  
	printf("%9.5f  %10.8e\n", (BETA/N)*0, -0.5);
	
    for (int i=1; i<N; i++)
      printf("%9.5f  %10.8e\n", (BETA/N)*i, -Gcoarse[i]);
	
	printf("%9.5f  %10.8e\n", (BETA/N)*N, -0.5);
  }

//-------------------------------------------------------------------------------------------
 
  if (plot_type == "selfsym2a") { // Greens function

	alps::RealVectorObsevaluator G_up=sim.get_measurements()["Greens_up"];
	alps::RealVectorObsevaluator G_down=sim.get_measurements()["Greens_down"];
	//alps::RealObsevaluator sign=sim.get_measurements()["sign"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
	//double const2 = static_cast<double>(sim.get_parameters()["CONST2"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	
	std::valarray<double> Gsym(N+1);
	for (int i=0; i<N+1; i++)
	  Gsym[i] = (G_up.mean()[i]/G_up.error()[i]+G_up.mean()[N-i]/G_up.error()[N-i]+G_down.mean()[i]/G_down.error()[i]+G_down.mean()[N-i]/G_down.error()[N-i])/(1/G_up.error()[i]+1/G_up.error()[N-i]+1/G_down.error()[i]+1/G_down.error()[N-i]);
	
	std::valarray<double> Gcoarse(Gsym);
	
	int k_max=J;
	
	for (int k=1; k<=k_max; k++) {
	for (int i=N/BETA*z*(1+(k-1)*0.2); i<N/BETA*z*(1.+k*0.2); i++) {
	  double coarse=0;
	  for (int n=-k; n<=k; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k+1);
	}

	for (int i=N/BETA*(1-z*(1+k*0.2)); i<N/BETA*(1-z*(1.+(k-1)*0.2)); i++) {
	  double coarse=0;
	  for (int n=-k; n<=k; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k+1);
	}
	}
	
	for (int i=N/BETA*z*(1+k_max*0.2); i<N/BETA*(1-z*(1.+k_max*0.2)); i++) {
	  double coarse=0;
	  for (int n=-k_max; n<=k_max; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k_max+1);
	}
			  
	printf("%9.5f  %10.8e\n", (BETA/N)*0, -0.5);
	
    for (int i=1; i<N; i++)
      printf("%9.5f  %10.8e\n", (BETA/N)*i, -Gcoarse[i]);
	
	printf("%9.5f  %10.8e\n", (BETA/N)*N, -0.5);
  }

//--------------------------------------------------------------------------------------------

  if (plot_type == "kin") { // Greens function
  
	alps::RealVectorObsevaluator G_up=sim.get_measurements()["Greens_up"];
	alps::RealVectorObsevaluator G_down=sim.get_measurements()["Greens_down"];
	//alps::RealObsevaluator sign=sim.get_measurements()["sign"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	double t = static_cast<double>(sim.get_parameters()["t"]);
	
	alps::RealObsevaluator n_up=sim.get_measurements()["n_up"];
	alps::RealObsevaluator n_down=sim.get_measurements()["n_down"];
    
	std::valarray<double> Gsym(N+1);
	std::valarray<double> Gsym_error(N+1);
	for (int i=0; i<N+1; i++) {
	  Gsym[i] = (G_up.mean()[i]+G_down.mean()[i])/2;
	  Gsym_error[i] = (G_up.error()[i]+G_down.error()[i])/(2*sqrt(2.));
	}
	
	std::valarray<double> G2(N+1);
	for (int i=0; i<N+1; i++) {
	  G2[i] = Gsym[i]*Gsym[N-i];
	}

	double kin_triangle=(0.25+0.25)/2.;
	double kin_simpson=0.25+0.25;
	for (int i=1; i<N; i++)
	  kin_triangle += G2[i];
	
	kin_triangle *= -2*BETA/N*t/BETA;  
	
	for (int i=1; i<N-1; i+=2) {
	  kin_simpson += 4*G2[i];
	  kin_simpson += 2*G2[i+1];
    }	  
	kin_simpson += 4*G2[N-1];
	  
	kin_simpson *= -2*(1./3.)*BETA/N*t/BETA;  
	
	printf("kin_triangle/t: %9.5f   kin_simpson/t: %9.5f\n", kin_triangle, kin_simpson);
	
	if(sim.get_parameters().defined("OVERLAP")) {
	  double u = static_cast<double>(sim.get_parameters()["U"]);
	  alps::RealObsevaluator overlap=sim.get_measurements()["overlap"];
	  printf("E_pot/t: %9.5f +/- %9.5f\n", u/t*overlap.mean(), u/t*overlap.error());
	  printf("E_tot/t: %9.5f\n", kin_simpson+u/t*overlap.mean());

	}
  }

//-------------------------------------------------------------------------------------------
 
  if (plot_type == "self") { // Greens function
  
    alps::RealObsevaluator n_up=sim.get_measurements()["n_up"];
	alps::RealObsevaluator n_down=sim.get_measurements()["n_down"];
	alps::RealVectorObsevaluator G_up=sim.get_measurements()["Greens_up"];
	alps::RealVectorObsevaluator G_down=sim.get_measurements()["Greens_down"];
	//alps::RealObsevaluator sign=sim.get_measurements()["sign"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
	//double const2 = static_cast<double>(sim.get_parameters()["CONST2"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	
	std::valarray<double> Gsym(N+1);
	for (int i=0; i<N+1; i++) {
	  Gsym[i] = (G_up.mean()[i]+G_down.mean()[i])/2.;
	}
	
	std::valarray<double> Gcoarse(Gsym);

	for (int i=N/BETA*z*1; i<N/BETA*z*1.2; i++) {
	  double coarse = (Gsym[i-1]+Gsym[i]+Gsym[i+1])/3.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.2; i<N/BETA*z*1.4; i++) {
	  double coarse = (Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2])/5.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.4; i<N/BETA*z*1.6; i++) {
	  double coarse = (Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3])/7.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.6; i<N/BETA*z*1.8; i++) {
	  double coarse = (Gsym[i-4]+Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3]+Gsym[i+4])/9.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*z*1.8; i<N/BETA*(1-z*1.8); i++) {
	  double coarse = (Gsym[i-5]+Gsym[i-4]+Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3]+Gsym[i+4]+Gsym[i+5])/11.;
	  Gcoarse[i] = coarse;
	}
	

	for (int i=N/BETA*(1-z*1.2); i<N/BETA*(1-z*1); i++) {
	  double coarse = (Gsym[i-1]+Gsym[i]+Gsym[i+1])/3.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*(1-z*1.4); i<N/BETA*(1-z*1.2); i++) {
	  double coarse = (Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2])/5.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*(1-z*1.6); i<N/BETA*(1-z*1.4); i++) {
	  double coarse = (Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3])/7.;
	  Gcoarse[i] = coarse;
	}
	for (int i=N/BETA*(1-z*1.8); i<N/BETA*(1-z*1.6); i++) {
	  double coarse = (Gsym[i-4]+Gsym[i-3]+Gsym[i-2]+Gsym[i-1]+Gsym[i]+Gsym[i+1]+Gsym[i+2]+Gsym[i+3]+Gsym[i+4])/9.;
	  Gcoarse[i] = coarse;
	}
	
	printf("%9.5f  %10.8e\n", (BETA/N)*0, -(1-n_up.mean()+1-n_down.mean())/2);	  
	
    for (int i=1; i<N; i++)
      printf("%9.5f  %10.8e\n", (BETA/N)*i, -Gcoarse[i]);
	
	printf("%9.5f  %10.8e\n", (BETA/N)*N, -(n_up.mean()+n_down.mean())/2);

  }
  

  //-------------------------------------------------------------------------------------------------------------------------------------------
  
  if (plot_type == "self2a") { // Greens function

    alps::RealObsevaluator n_up=sim.get_measurements()["n_up"];
	alps::RealObsevaluator n_down=sim.get_measurements()["n_down"];
	alps::RealVectorObsevaluator G_up=sim.get_measurements()["Greens_up"];
	alps::RealVectorObsevaluator G_down=sim.get_measurements()["Greens_down"];
	//alps::RealObsevaluator sign=sim.get_measurements()["sign"];
    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
	//double const2 = static_cast<double>(sim.get_parameters()["CONST2"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	
	std::valarray<double> Gsym(N+1);
	for (int i=0; i<N+1; i++)
	  Gsym[i] = (G_up.mean()[i]/G_up.error()[i]+G_down.mean()[i]/G_down.error()[i])/(1/G_up.error()[i]+1/G_down.error()[i]);
	
	std::valarray<double> Gcoarse(Gsym);
	
	int k_max=J;
	
	for (int k=1; k<=k_max; k++) {
	for (int i=N/BETA*z*(1+(k-1)*0.2); i<N/BETA*z*(1.+k*0.2); i++) {
	  double coarse=0;
	  for (int n=-k; n<=k; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k+1);
	}

	for (int i=N/BETA*(1-z*(1+k*0.2)); i<N/BETA*(1-z*(1.+(k-1)*0.2)); i++) {
	  double coarse=0;
	  for (int n=-k; n<=k; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k+1);
	}
	}
	
	for (int i=N/BETA*z*(1+k_max*0.2); i<N/BETA*(1-z*(1.+k_max*0.2)); i++) {
	  double coarse=0;
	  for (int n=-k_max; n<=k_max; n++)
	    coarse += Gsym[i+n];
	  
	  Gcoarse[i] = coarse/(2*k_max+1);
	}
			  
	printf("%9.5f  %10.8e\n", (BETA/N)*0, -(1-n_up.mean()+1-n_down.mean())/2);
    
	for (int i=1; i<N; i++)
      printf("%9.5f  %10.8e\n", (BETA/N)*i, -Gcoarse[i]);
	
	printf("%9.5f  %10.8e\n", (BETA/N)*N, -(n_up.mean()+n_down.mean())/2);
  }

//------------------------------------------------------------------------------------
    

  if (plot_type == "o") { // Greens function
    alps::RealVectorObsevaluator order_up=sim.get_measurements()["order_up"];
    alps::RealVectorObsevaluator order_down=sim.get_measurements()["order_down"];

    double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    int N = static_cast<int>(sim.get_parameters()["N"]);
	int N_order = static_cast<int>(sim.get_parameters()["N_ORDER"]);
    
    for (int i=0; i<N_order; i++)
      printf("index: %4d  order_down: %9.7e +- %9.7e   order_down: %9.7e +- %9.7e\n", i, order_up.mean()[i], order_up.error()[i], order_down.mean()[i], order_down.error()[i]);

  }
    

  if (plot_type == "s") { // Greens function
    //alps::RealObsevaluator f=sim.get_measurements()["full"];
	alps::RealObsevaluator sign=sim.get_measurements().get<alps::AbstractSimpleObservable<double> >("sign");
	alps::RealObsevaluator sign_up=sim.get_measurements().get<alps::AbstractSimpleObservable<double> >("sign_up");
    alps::RealObsevaluator sign_down=sim.get_measurements().get<alps::AbstractSimpleObservable<double> >("sign_down");

	//alps::RealObsevaluator s3=sim.get_measurements().get<alps::AbstractSimpleObservable<double> >("s3");
	//alps::RealObsevaluator s4=sim.get_measurements().get<alps::AbstractSimpleObservable<double> >("s4");
	//alps::RealObsevaluator s5=sim.get_measurements().get<alps::AbstractSimpleObservable<double> >("s5");

    //double BETA = static_cast<double>(sim.get_parameters()["BETA"]);
    //int N = static_cast<int>(sim.get_parameters()["N"]);
    
	printf("sign: %9.7e +- %9.7e  sign_up: %9.7e +- %9.7e  sign_down: %9.7e +- %9.7e\n", sign.mean(), sign.error(), sign_up.mean(), sign_up.error(), sign_down.mean(), sign_down.error());
	//printf("s3: %9.7e +- %9.7e\n", s3.mean(), s3.error());
	//printf("s4: %9.7e +- %9.7e\n", s4.mean(), s4.error());
	//printf("s5: %9.7e +- %9.7e\n", s5.mean(), s5.error());

  }
  

*/

}

int main(int argc, char** argv)
{
  try {
    alps::scheduler::SimpleMCFactory<alps::scheduler::DummyMCRun> factory;
    alps::scheduler::init(factory);

    if (argc<3) {
      std::cerr << "Usage: " << argv[0] << " selfsym tau inputfile (will average over 3, then 5, 7, 9, 11 bins in the interval[tau, BETA-tau])\n";	
	  std::cerr << "Usage: " << argv[0] << " selfsym2 tau n inputfile (will average over 3, then 5, 7, 9, ... (n increases from tau to tau*(1+n/5)) bins in the interval[tau, BETA-tau])\n";	
      std::cerr << "Usage: " << argv[0] << " new_J z J [therm [max_meas]] inputfile\n";
      std::cerr << "Usage: " << argv[0] << " Binder_L z inputfile ({BETA,L}, alpha fixed)\n";
      std::cerr << "or:    " << argv[0] << " MagSqr_alpha inputfile ({BETA,alpha}, L fixed)\n";
      std::cerr << "or:    " << argv[0] << " MagSqr_beta inputfile ({alpha,BETA}, L fixed)\n";
      std::cerr << "or:    " << argv[0] << " MagSqr_G inputfile ({BETA,G}, L fixed)\n";
      std::cerr << "or:    " << argv[0] << " nu_J [ thermalisation [J] ] inputfile\n";
      std::exit(-1);
    }    

    std::ifstream infile;
    std::string plot_type = argv[1];
    std::string name_base = argv[argc-1];
    double z(0), J(0);
    int therm(0), max_meas(10000000);

    if (plot_type == "Binder_L")
      z = boost::lexical_cast<double>(argv[2]);
    else if (plot_type == "selfsym" || plot_type == "self")
      z = boost::lexical_cast<double>(argv[2]);	  
	else if (plot_type == "selfsym2" || plot_type == "self2" || plot_type == "selfsym2a" || plot_type == "self2a") {
      z = boost::lexical_cast<double>(argv[2]);
      J = boost::lexical_cast<double>(argv[3]);
    }
	else if (plot_type == "new_J") {
      z = boost::lexical_cast<double>(argv[2]);
      J = boost::lexical_cast<double>(argv[3]);
      if (argc-1>4) therm = boost::lexical_cast<int>(argv[4]);
      if (argc-1>5) max_meas = boost::lexical_cast<int>(argv[5]);
    }
    else if (plot_type == "nu_J") {
      if (argc-1>2) therm = boost::lexical_cast<int>(argv[2]);
      if (argc-1>3) J = boost::lexical_cast<double>(argv[3]);
    }

    // prepare filenames and start evaluation
    for(int i=1;;i++) {

      std::string name = name_base;
      std::string task = ".task";
      task.append(boost::lexical_cast<std::string>(i)).append(".out.xml");
      name.append(task);

      infile.open(name.c_str());

      if (!infile.good()) break;

      boost::filesystem::path p(name.c_str(),boost::filesystem::native);

      evaluate(p, z, J, therm, max_meas, plot_type, std::cout, task);

      infile.close();

    }

  }
  catch (std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << "\n";
    std::exit(-5);
  }
}
