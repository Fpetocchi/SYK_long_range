/***************************************************************************
 * Impurity solver based on a Z-expansion in the impurity-bath hybridization
 * (c) 2005 Philipp Werner
 ***************************************************************************/

#include "impurity.h"
#include "moves.h"
#include <alps/alea.h>

using namespace std;
using namespace alps;


IsingSimulation::IsingSimulation(const alps::ProcessList& where,const alps::Parameters& p,int node)
: alps::scheduler::MCRun(where,p,node),
sweeps(0),
thermalization_sweeps(static_cast<int>(parms["THERMALIZATION"])),
total_sweeps(static_cast<int>(parms["SWEEPS"])),
mu(static_cast<double>(parms["MU"])),
mu_e(static_cast<int>(parms["FLAVORS"])),
t(static_cast<double>(parms["t"])),
full_line(static_cast<int>(parms["FLAVORS"]),0),
sign(static_cast<int>(parms["FLAVORS"])),
G_meas(static_cast<int>(parms["FLAVORS"])*(static_cast<int>(parms["N"])+1))

{

    int N = static_cast<int>(parms["N"]);
    int FLAVORS = static_cast<int>(parms["FLAVORS"]);
    int Np1 = N+1;
    double BETA = static_cast<double>(parms["BETA"]);

    //fill K_table
    ifstream input_K(boost::lexical_cast<std::string>(parms["K_FILE"]).c_str());
    int N_K_TABLE = static_cast<int>(parms["N_K_TABLE"]);
    K_table.resize(FLAVORS);
    for (int i=0; i<FLAVORS; i++) {
        K_table[i].resize(FLAVORS);
        for (int j=0; j<FLAVORS; j++) {
            K_table[i][j].resize(N_K_TABLE+1);
            for (int k=0; k<N_K_TABLE+1; k++) {
                double dummy;
                input_K >> dummy;
                input_K >> K_table[i][j][k]; // need to choose reasonable format for K_FILE
            }
        }
    }


    /*
     // IN CASE OF A SINGLE BOSON MODEL
     double lambda = static_cast<int>(parms["LAMBDA"]);
     int omega0 = static_cast<int>(parms["W0"]);
     int N_K_TABLE=10000;
     K_table.resize(N_K_TABLE+1);
     for (int i=0; i<N_K_TABLE+1; i++) {
     double tau=i*BETA/N_K_TABLE;
     K_table[i] = -(lambda/omega0)*(lambda/omega0)*(cosh((BETA/2-tau)*omega0)-cosh((BETA/2)*omega0))/sinh((BETA/2)*omega0);
     }
    */


    // define interactions between flavors
    u.resize(FLAVORS, FLAVORS);

    ifstream infile_u(boost::lexical_cast<std::string>(parms["U_MATRIX"]).c_str());
    for (int i=0; i<FLAVORS; i++) {
        for (int j=0; j<FLAVORS; j++) {
            infile_u >> u(i,j);
            // HERE COMES THE CORRECTION DUE TO DISCONTINUITY (OR SIMPLY READ IN SCEENED U)
            //if (i!=j) u(i,j)-=2*lambda*lambda/omega0;
        }
    }

    // for semi-circular density of states we can directly compute F_{\sigma}(\tau) = t^2*G_{-\sigma}(-\tau)
    std::vector<std::vector<double> > G;

    G.resize(FLAVORS);
    F.resize(FLAVORS);
    for (int i=0; i<FLAVORS; i++) {
        G[i].resize(Np1);
        F[i].resize(Np1);
    }

    // read G from file
    ifstream infile_G(boost::lexical_cast<std::string>(parms["G"]).c_str());
    for (int i=0; i<Np1; i++) {
        double dummy;
        infile_G >> dummy;
        for (int j=0; j<FLAVORS; j++) {
            infile_G >> G[j][i];
        }
    }

    // THIS IS THE DMFT SELFCONSISTENCY FOR SEMI-CIRCLE
    // WE SHOULD READ IN DIRECTLY DELTA
    for (int i=0; i<Np1; i++) {
        for (int j=0; j<FLAVORS; j++) {
            F[j][i] = -t*t*G[j][N-i];
        }
    }

    for (int i=0; i<FLAVORS; i++) {
        //infile_e >> mu_e[i];
        mu_e[i] = mu;
        // HERE COMES THE CORRECTION DUE TO DISCONTINUITY (OR SIMPLY READ IN SHIFTED MU)
        //mu_e[i] -= (FLAVORS-1)*lambda*lambda/omega0;
    }


    /*
     // E_VEC IS AVERAGE BAND ENERGY, NEED TO CHECK IF WE NEED THIS
     ifstream infile_e(boost::lexical_cast<std::string>(parms["E_VEC"]).c_str());
     for (int i=0; i<FLAVORS; i++) {
     infile_e >> mu_e[i];
     mu_e[i] *= -t;
     mu_e[i] += mu;
     }

     F.resize(FLAVORS);
     for (int i=0; i<FLAVORS; i++) {
     F[i].resize(Np1);
     }

     // read F from file
     ifstream infile(boost::lexical_cast<std::string>(parms["G"]).c_str());
     for (int i=0; i<Np1; i++) {
     double dummy;
     infile >> dummy;
     for (int j=0; j<FLAVORS; j++)
     infile >> F[j][i];
     }
     for (int i=0; i<Np1; i++) {
     for (int j=0; j<FLAVORS; j++) {
     F[j][i] *= t*t;
     }
     }
    */

    segments.resize(FLAVORS);
    M.resize(FLAVORS);
    // initialize list of segments
    for (int i=0; i<FLAVORS; i++) {
        segments[i].clear();
        M[i].clear();
    }

    // create measurement objects
    measurements << RealVectorObservable("n");
    measurements << RealVectorObservable("order");
    measurements << RealVectorObservable("Greens");
    measurements << RealVectorObservable("nn_corr");
    measurements << RealObservable("sign");

    measurements << RealVectorObservable("nhist");
    measurements << RealVectorObservable("szabshist");


    //if (static_cast<int>(parms["OVERLAP"]))
    //measurements << RealObservable("overlap");

}

/*
 void IsingSimulation::load(alps::IDump& dump) {
 dump >> sweeps;
 if(!where.empty()) // skip reading the spins if we are just evaluating
 dump >> segments_up >> full_line_up >> segments_down >> full_line_down;
 else
 measurements.compact(); // and throw away time series
 }

 void IsingSimulation::save(alps::ODump& dump) const
 {
 dump << sweeps << segments_up << full_line_up << segments_down << full_line_down;
 }
*/

bool IsingSimulation::change_parameter(const std::string& name, const alps::StringValue& value)
{
    if(name=="SWEEPS")
        total_sweeps=static_cast<uint32_t>(value);
    else if (name=="THERMALIZATION" && !is_thermalized())
        thermalization_sweeps=static_cast<uint32_t>(value);
    else
        return false; // cannot do it
    return true; // could do it
}

bool IsingSimulation::is_thermalized() const
{
    return (sweeps >= thermalization_sweeps);
}

double IsingSimulation::work_done() const
{
    return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) : 0.);

}

void IsingSimulation::dostep()
{

    // increment sweep count
    sweeps++;

    double BETA = static_cast<double>(parms["BETA"]);
    int N = static_cast<int>(parms["N"]);
    int FLAVORS = static_cast<int>(parms["FLAVORS"]);
    int N_order = static_cast<int>(parms["N_ORDER"]);
    int N_corr = static_cast<int>(parms["N_CORR"]);
    int N_meas = static_cast<int>(parms["N_MEAS"]);
    int N_shift = static_cast<int>(parms["N_SHIFT"]);
    int N_flip = static_cast<int>(parms["N_FLIP"]);
    int N_move = static_cast<int>(parms["N_MOVE"]);
    int N_swap = static_cast<int>(parms["N_SWAP"]);
    int overlap = static_cast<int>(parms["OVERLAP"]);

    times full_segment(0,BETA);

    std::set<times>::iterator it1, it2;

    std::valarray<double> order_meas(N_order*FLAVORS);
    std::valarray<double> n_meas(FLAVORS);
    double sign_meas=0, s=1;
    G_meas = 0;
    n_meas = 0;

    for (int i=0; i<N_meas; i++) {

        for (int j=0; j<FLAVORS; j++) {

            if (segments[j].size() == 0) {
                // insert or remove full line
                insert_remove_full_line(random_01, mu_e[j], u, BETA, full_line[j], segments, full_line, j);
            }

            insert_remove_antisegment(random_01, BETA*random_01(), N, BETA, mu_e[j], u, F[j], full_line[j], segments[j], M[j], sign[j], segments, full_line, j, K_table);

            if (!full_line[j]) {
                // local update
                insert_remove_segment(random_01, BETA*random_01(), N, BETA, mu_e[j], u, F[j], segments[j], M[j], sign[j], segments, full_line, j, K_table);

                // shift segments
                for (int k=0; k<N_shift; k++)
                    shift_segment(random_01, segments[j], N, BETA, mu_e[j], u, F[j], M[j], sign[j], segments, full_line, j, K_table);

                /*
                 // flip segment
                 for (int i=0; i<N_flip; i++)
                 flip_segment(random_01, segments_up, N, BETA, M_up, sign_up, sign_down, F_down, M_down, segments_down, full_line_down);
                */
            }

            // measure perturbation order

            if (segments[j].size()<N_order)
                order_meas[j*N_order+segments[j].size()] += 1;


            // measure Green functions

            if (segments[j].size()>0) {
                for (int i=0; i<M[j].size1(); i++) {
                    (i==0 ? it1 = segments[j].begin() : it1++);
                    for (int k=0; k<M[j].size1(); k++) {
                        (k==0 ? it2 = segments[j].begin() : it2++);
                        if (M[j](k,i)!=0) {
                            double argument = it1->t_end()-it2->t_start();
                            double bubble_sign=1;
                            if (argument > 0) {
                                bubble_sign = 1;
                            }
                            else {
                                bubble_sign = -1;
                                argument += BETA;
                            }

                            int index = argument/BETA*N+0.5;
                            G_meas[j*(N+1)+index] += M[j](k,i)*bubble_sign/(BETA*BETA);
                        }
                    }
                }
            }

            s *= sign[j];

            n_meas[j] += compute_overlap(full_segment, segments[j], full_line[j], BETA)/BETA;

        }

        sign_meas += s;

        // swap updates may be needed in case of trapping in a symmetry-broken state
        /*
         if (i%N_swap==1) {
         int orbital=random_01()*(FLAVORS/2);
         if (M[2*orbital].size1()!=0 && M[2*orbital+1].size1()!=0)
         swap_segments(random_01, BETA, F[orbital], F[orbital+1], segments[orbital], segments[orbital+1], full_line[orbital], full_line[orbital+1], sign[orbital], sign[orbital+1], M[orbital], M[orbital+1]);
         }
        */

        /*
         // global spin flip (does not require overlap calculation)
         if (i%N_swap==1){
         //check if doable
         bool ok=true;
         for(int j=0; j<FLAVORS; j++) if(M[j].size1()==0) ok=false;
         if(ok) swap_spins(random_01, BETA, FLAVORS, F, segments, full_line, sign, M);
         }
        */


    }

    // measure density correlation functions

    std::vector<std::vector<double> > n_vectors(FLAVORS);
    std::set<times>::iterator it;

    for (int flavor=0; flavor<FLAVORS; ++flavor) {

        n_vectors[flavor].resize(N_corr+1, 1);

        if (segments[flavor].size()==0) {
            if (full_line[flavor]==0) {
                for (int i=0; i<n_vectors[flavor].size(); ++i)
                    n_vectors[flavor][i]=0;
            }
        }
        else {
            it=segments[flavor].end(); it--;
            if (it->t_end()<it->t_start())
                n_vectors[flavor][0]=1;
            else
                n_vectors[flavor][0]=0;

            // mark segment start and end points
            int index;
            for (it=segments[flavor].begin(); it!=segments[flavor].end(); it++) {
                index = it->t_start()/BETA*N_corr+1;
                n_vectors[flavor][index] *= -1;
                index = it->t_end()/BETA*N_corr+1;
                n_vectors[flavor][index] *= -1;
            }
            // fill vector with occupation number
            for (int i=1; i<n_vectors[flavor].size(); i++) {
                if (n_vectors[flavor][i]==-1)
                    n_vectors[flavor][i]=1-n_vectors[flavor][i-1];
                else
                    n_vectors[flavor][i]=n_vectors[flavor][i-1];
            }
        }

    }

    // compute n(\tau)n(0)
    std::valarray<double> nn_corr_meas(FLAVORS*(FLAVORS+1)/2*(N_corr+1));
    int position=0;
    for (int flavor1=0; flavor1<FLAVORS; ++flavor1) {
        for (int flavor2=0; flavor2<=flavor1; ++flavor2) {

            for (int i=0; i<N_corr+1; ++i) {
                for (int index=0; index<N_corr+1; ++index) {

                    int j=i+index;
                    if (j>N_corr) j -= N_corr;

                    nn_corr_meas[position+index] += n_vectors[flavor1][i]*n_vectors[flavor2][j];

                }
            }

            position += (N_corr+1);
        }
    }

    nn_corr_meas /= (N_corr+1);
    measurements.get<RealVectorObservable>("nn_corr") << nn_corr_meas;


    std::valarray<double> nhist_meas(FLAVORS+1);
    nhist_meas *= 0.;
    for (int i=1; i<n_vectors[0].size(); i++) {
        int ntmp=0;
        for (int j=0; j<FLAVORS; j++) {
            ntmp+=n_vectors[j][i];
        }
        nhist_meas[ntmp]+=1./n_vectors[0].size();
    }

    measurements.get<RealVectorObservable>("nhist") << nhist_meas;


    std::valarray<double> szabshist_meas(FLAVORS/2+1);
    szabshist_meas *= 0.;
    for (int i=1; i<n_vectors[0].size(); i++) {
        int ntmp=0;
        for (int j=0; j<FLAVORS; j+=2 ) {
            ntmp+=n_vectors[j][i]-n_vectors[j+1][i];
        }
        if (ntmp<0) ntmp*=-1;
        szabshist_meas[ntmp]+=1./n_vectors[0].size();
    }

    measurements.get<RealVectorObservable>("szabshist") << szabshist_meas;


    order_meas /= N_meas;
    measurements.get<RealVectorObservable>("order") << order_meas;

    G_meas *= (1.*N)/N_meas;
    measurements.get<RealVectorObservable>("Greens") << G_meas;

    sign_meas /= N_meas;
    measurements.get<RealObservable>("sign") << sign_meas;

    n_meas /= N_meas;
    measurements.get<RealVectorObservable>("n") << n_meas;

    /*
     // for measurement of double occupancy (needs to be adapted)
     if (overlap) {
     if (segments[j].size()>0) {
     for (it1=segments[j].begin(); it1 != segments[j].end(); it1++) {
     overlap_meas += compute_overlap(*it1, segments_down, full_line_down, BETA)/BETA;
     }
     }
     else if (full_line_up){
     for (it1=segments_down.begin(); it1 != segments_down.end(); it1++) {
     overlap_meas += compute_overlap(full_segment, segments_down, full_line_down, BETA)/BETA;
     }
     }
     measurements.get<RealObservable>("overlap") << overlap_meas;
     }
    */

}
