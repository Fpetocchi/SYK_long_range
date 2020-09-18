#ifndef ___UPDATE___
#define ___UPDATE___

#include "impurity.h"


// invert matrix A and calculate its determinant
inline void invert(dense_matrix & A, double & det) {
    
    dense_matrix B(A.size1(), A.size1());
    if (A.size1()>0) {
        B = boost::numeric::ublas::identity_matrix<double>(A.size1());
        
        boost::numeric::bindings::lapack::gesv(A, B);
        swap(A,B);
        det = 1;
        for (int i=0; i<A.size1(); i++) {
            det *= B(i,i);
        }
        det = std::fabs(det);
    }
}

template <class G> void construct_matrix(dense_matrix & M, std::list<times> & segments, double BETA,  G& F) {
    int N = segments.size();
    M.resize(N,N);
    int row=-1;
    int column=-1;
    for (std::list<times>::iterator it1=segments.begin(); it1!=segments.end(); it1++) {
        row++;
        for (std::list<times>::iterator it2=segments.begin(); it2!=segments.end(); it2++) {
            column++;
            
            double argument = it1->t_end()-it2->t_start();
            double sign = 1;
            if (argument<0) {
                argument += BETA;
                sign = -1;
            }
            M(row,column) = interpolate_F(argument, BETA, F)*sign;
        }
        column = -1;
    }
    
}

template <class G> double construct_inverse(dense_matrix & M, std::list<times> & segments, double BETA,  G& F) {
    construct_matrix(M, segments, BETA, F);
    double dummy;
    invert(M, dummy);
    return dummy;
}



// determine F(\tau)
template <class G> inline double interpolate_F(double t, double BETA, G& F) {
    
    double sign=1;
    if (t<0) {
        t += BETA;
        sign=-1;
    }
    
    int N = F.size()-1;
    double n = t/BETA*N;
    int n_lower = n; // interpolate linearly between n_lower and n_lower+1
    
    return sign*(F[n_lower] + (n-n_lower)*(F[n_lower+1]-F[n_lower]));
    
}



// nonlocal stuff

inline double H(double tau, double beta, std::vector<double>& K_table) {
    
    if (tau<0) tau+=beta;
    
    double i=tau/beta*(K_table.size()-1);
    int i_lower=(int) i;
    
    return K_table[i_lower]+(i-i_lower)*(K_table[i_lower+1]-K_table[i_lower]);
}

// this is the version without self-interaction
double nonlocal(double ts, double te, std::vector<segment_container_t>& other_segments, double beta, int this_flavor, k_table_t& K_table, int insert_remove) {
    
    double nonloc=0.;
    
    for (int flavor=0; flavor<other_segments.size(); flavor++) {
        for(std::set<times>::iterator it=other_segments[flavor].begin(); it!=other_segments[flavor].end(); it++) {
            nonloc += -H(it->t_end()-te, beta, K_table[flavor][this_flavor])+H(it->t_end()-ts, beta, K_table[flavor][this_flavor])+H(it->t_start()-te, beta, K_table[flavor][this_flavor])-H(it->t_start()-ts, beta, K_table[flavor][this_flavor]);
        }
    }
    
    if (insert_remove==0) { // insert
        nonloc += H(te-ts, beta, K_table[this_flavor][this_flavor]);
    }
    else{ // remove
        nonloc -= -2*H(0, beta, K_table[this_flavor][this_flavor])+H(te-ts, beta, K_table[this_flavor][this_flavor]); // note: H(0)=0 in the model without self-interaction, so the first term is zero
    }
    
    return nonloc;
}


double nonlocal_shift(double te_ins, double te_rem, std::vector<segment_container_t>& other_segments, double beta, int this_flavor, k_table_t& K_table) {
    
    double nonloc=0.;
    
    for (int flavor=0; flavor<other_segments.size(); flavor++) {
        for(std::set<times>::iterator it=other_segments[flavor].begin(); it!=other_segments[flavor].end(); it++) {
            nonloc += -H(it->t_end()-te_ins, beta, K_table[flavor][this_flavor])+H(it->t_start()-te_ins, beta, K_table[flavor][this_flavor]);
            nonloc -= -H(it->t_end()-te_rem, beta, K_table[flavor][this_flavor])+H(it->t_start()-te_rem, beta, K_table[flavor][this_flavor]);
        }
    }
    
    nonloc -= -H(te_rem-te_ins, beta, K_table[this_flavor][this_flavor]); // inexistent bond
    nonloc += -H(te_rem-te_rem, beta, K_table[this_flavor][this_flavor]); // local contribution at te doesn't change
    
    return nonloc;
}


// compute distances up/down to the next segment and iterators of these segments
// note: s_down always points to a physical segment, while s_up may point to segments.end()
template <class S> void compute_intervals(double t, double BETA, double& t_up, double& t_down, S& segments, typename S::iterator& s_up, typename S::iterator& s_down) {
    
    if (segments.size() == 0) {
        t_up = BETA;
        t_down = BETA;
        s_up = segments.end();
        s_down = segments.end();
    }
    else {
        
        s_up = lower_bound(segments.begin(), segments.end(), t);
        
        if (s_up == segments.begin()) {
            s_down = segments.end(); s_down--;
            if (s_down->t_end() < s_down->t_start())
                t_down = t - s_down->t_end();
            else
                t_down = t + BETA - s_down->t_end();
        }
        else {
            s_down = s_up; s_down--;
            if (s_down->t_end()>s_down->t_start())
                t_down = t - s_down->t_end();
            else
                t_down = t - (BETA+s_down->t_end());
        }
        
        if(s_up == segments.end()) {
            t_up = BETA - t + segments.begin()->t_start();
        }
        else {
            t_up = s_up->t_start() - t;
        }
        
    }
    
}

// compute overlap between a segment and a list of segments
// requires segment with 0<=t_begin<t_end<=BETA
template <class S> inline double segment_overlap(times segment, S& other_segments, int other_full_line, double BETA) {
    
    double length = (segment.t_start()<segment.t_end() ? segment.t_end()-segment.t_start() : segment.t_end()-segment.t_start()+BETA);
    double t_final = segment.t_start()+length;
    double t = segment.t_start();
    double t_final_segment;
    double other_length=0;
    if (other_full_line==1)
        other_length=length;
    else if (other_segments.size()>0){
        typename S::iterator it;
        it = lower_bound(other_segments.begin(), other_segments.end(), t);
        
        if (it!=other_segments.begin()) {
            it--;
            t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
            if (t<t_final_segment) {
                other_length += (t_final_segment<t_final ? t_final_segment-t : t_final-t);
            }
            it++;
            
        }
        while(it!=other_segments.end() && it->t_start()<t_final) {
            t_final_segment = (it->t_start()<it->t_end() ? it->t_end() : it->t_end()+BETA);
            other_length += (t_final_segment<t_final ? t_final_segment-it->t_start() : t_final-it->t_start());
            it++;
        }
        // check if last segment overlaps
        it=other_segments.end();
        it--;
        if (it->t_end()<it->t_start() && t<it->t_end()) {
            other_length += (t_final<it->t_end() ? t_final-t : it->t_end()-t);
        }
    }
    return other_length;
}


template <class S> inline double compute_overlap(times segment, S& other_segments, int other_full_line, double BETA) {
    if (segment.t_start()<segment.t_end())
        return segment_overlap(segment, other_segments, other_full_line, BETA);
    else {
        double other_length=0;
        times segment1(0,segment.t_end());
        times segment2(segment.t_start(), BETA);
        other_length += segment_overlap(segment1, other_segments, other_full_line, BETA);
        other_length += segment_overlap(segment2, other_segments, other_full_line, BETA);
        return other_length;
    }
}


// functions required to compute determinant ratios and perform fast matrix updates

template <class G, class S, class V> double det_rat_up(times & new_segment, dense_matrix & M, S& segments_old, G& F, V& Fs, V& Fe, double BETA, double & det_rat_sign, double & overlap) {
    
    typename S::iterator it=segments_old.begin();
    for (int i=0; i<segments_old.size(); i++) {
        Fe[i] = interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);
        Fs[i] = interpolate_F(it->t_end()-new_segment.t_start(), BETA, F);
        it++;
    }
    
    double det_rat = interpolate_F(new_segment.t_end()-new_segment.t_start(), BETA, F);
    
    for (int i=0; i<M.size1(); i++) {
        for (int j=0; j<M.size1(); j++) {
            det_rat -= Fe[i]*M(i,j)*Fs[j];
        }
    }
    
    // take care of sign changes produced by segments which "wind around"
    if (new_segment.t_end() < new_segment.t_start()) {
        det_rat *= -1;
        overlap = -1;
    }
    else {
        overlap = 1;
    }
    
    if (det_rat < 0) {
        det_rat_sign = -1;
        det_rat *= -1;
    }
    else {
        det_rat_sign = 1;
    }
    
    return det_rat;
}


void compute_M_up(int k, dense_matrix & M, vector_t& Fs, vector_t &Fe, double det_rat) {
    
    dense_matrix M_new(M.size1()+1,M.size1()+1);
    int i_new, j_new;
    
    // element (k,k)
    M_new(k,k) = 1./det_rat;
    
    // row k and column k
    for (int i=0; i<M.size1(); i++) {
        i_new = (i<k ? i : i+1);
        M_new(i_new,k) = 0;
        M_new(k,i_new) = 0;
        
        for (int n=0; n<M.size1(); n++) {
            M_new(i_new,k) -= M(i,n)*Fs[n];
            M_new(k,i_new) -= M(n,i)*Fe[n];
        }
        M_new(i_new,k) /= det_rat;
        M_new(k,i_new) /= det_rat;
    }
    
    // remaining elements
    for (int i=0; i<M.size1(); i++) {
        i_new = (i<k ? i : i+1);
        for (int j=0; j<M.size1(); j++) {
            j_new = (j<k ? j : j+1);
            M_new(i_new,j_new) = M(i,j) + det_rat*M_new(i_new,k)*M_new(k,j_new);
        }
    }
    
    swap(M_new, M);
    return;
}


template <class S> double det_rat_down(int k, dense_matrix & M, S& segments_old, double & det_rat_sign) {
    
    double det_rat = M(k,k);
    
    // take care of sign changes produced by segments which "wind around"
    if (k==segments_old.size()-1) {
        typename S::iterator it=segments_old.end(); it--;
        if (it->t_end() < it->t_start())
            det_rat *= -1;
    }
    
    if (det_rat < 0) {
        det_rat_sign = -1;
        det_rat *= -1;
    }
    else {
        det_rat_sign = 1;
    }
    
    return det_rat;
}


void compute_M_down(int k, dense_matrix & M) {
    
    dense_matrix M_new(M.size1()-1, M.size1()-1);
    int i_old, j_old;
    
    for (int i=0; i<M_new.size1(); i++) {
        i_old = (i<k ? i : i+1);
        for (int j=0; j<M_new.size1(); j++) {
            j_old = (j<k ? j : j+1);
            M_new(i,j) = M(i_old, j_old)-M(i_old,k)*M(k,j_old)/M(k,k);
        }
    }
    
    swap(M, M_new);
    
}

// move segment without changin its length
template <class G, class S> double det_rat_move(times & new_segment, int k, dense_matrix & M, S& segments_old, G& F, double BETA, double & det_rat_sign, double & overlap) {
    
    double F_i, F_j;
    typename S::iterator it1, it2;
    
    double det_rat = M(k,k)*interpolate_F(new_segment.t_end()-new_segment.t_start(), BETA, F);
    
    it1=segments_old.begin();
    for (int i=0; i<M.size1(); i++) {
        if (i != k) {
            F_i = interpolate_F(new_segment.t_end()-it1->t_start(), BETA, F);
            
            it2=segments_old.begin();
            for (int j=0; j<M.size1(); j++) {
                if (j != k) {
                    F_j = interpolate_F(it2->t_end()-new_segment.t_start(), BETA, F);
                    det_rat -= F_i*(M(k,k)*M(i,j)-M(i,k)*M(k,j))*F_j;
                }
                it2++;
            }
        }
        it1++;
    }
    
    overlap = 1;
    // take care of sign changes produced by segments which "wind around"
    if (k==segments_old.size()-1) {
        it1--;
        // check if last segment has been shifted across beta
        if ((new_segment.t_end()-new_segment.t_start())*(it1->t_end()-it1->t_start())<0) {
            det_rat *= -1;
            overlap = -1;
        }
    }
    
    if (det_rat < 0) {
        det_rat_sign = -1;
        det_rat *= -1;
    }
    else {
        det_rat_sign = 1;
    }
    
    return det_rat;
}


template <class G, class S> void compute_M_move(times & new_segment, int k, dense_matrix & M, S& segments_old, G& F, double BETA, double det_rat) {
    
    dense_matrix M_new(M.size1(),M.size1());
    //double argument;
    
    // row k and column k
    for (int i=0; i<M.size1(); i++) {
        if (i!=k) {
            M_new(i,k) = 0;
            M_new(k,i) = 0;
            
            typename S::iterator it=segments_old.begin();
            for (int n=0; n<M.size1(); n++) {
                if (n!=k) {
                    M_new(i,k) -= 1/det_rat*(M(k,k)*M(i,n)-M(i,k)*M(k,n))*interpolate_F(it->t_end()-new_segment.t_start(), BETA, F);
                    M_new(k,i) -= 1/det_rat*(M(k,k)*M(n,i)-M(n,k)*M(k,i))*interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);
                }
                it++;
            }
        }
        else {
            M_new(k,k) = M(k,k)/det_rat;
        }
    }
    
    // remaining elements
    for (int i=0; i<M.size1(); i++) {
        if (i!=k) {
            for (int j=0; j<M.size1(); j++) {
                if (j!=k)
                    M_new(i,j) = M(i,j) + (-M(i,k)*M(k,j)+det_rat*M_new(i,k)*M_new(k,j))/M(k,k);
            }
        }
    }
    
    swap(M_new, M);
    return;
}

// shift end point of segment
template <class G, class S> double det_rat_shift(times & new_segment, int k, dense_matrix & M, S& segments_old, G& F, double BETA, double & det_rat_sign, double & overlap) {
    
    typename S::iterator it;
    double det_rat = 0;
    
    it=segments_old.begin();
    for (int i=0; i<M.size1(); i++) {
        det_rat += interpolate_F(new_segment.t_end()-it->t_start(), BETA, F)*M(i,k);
        it++;
    }
    
    overlap = 1;
    // take care of sign changes produced by segments which "wind around"
    if (k==segments_old.size()-1) {
        it--;
        // check if last segment has been shifted across beta
        if ((new_segment.t_end()-new_segment.t_start())*(it->t_end()-it->t_start())<0) {
            det_rat *= -1;
            overlap = -1;
        }
    }
    
    if (det_rat < 0) {
        det_rat_sign = -1;
        det_rat *= -1;
    }
    else {
        det_rat_sign = 1;
    }
    
    return det_rat;
}


template <class G, class S> void compute_M_shift(times & new_segment, int k, dense_matrix & M, S& segments_old, G& F, double BETA, double det_rat) {
    
    std::vector<double> R(M.size1(),0), M_k(M.size1(),0), Fe(M.size1(),0);
    
    typename S::iterator it=segments_old.begin();
    for (int i=0; i<M_k.size(); i++) {
        M_k[i] = M(i,k);
        Fe[i] = interpolate_F(new_segment.t_end()-it->t_start(), BETA, F);
        it++;
    }
    
    for (int i=0; i<R.size(); i++) {
        if (i!=k) {
            for (int j=0; j<R.size(); j++)
                R[i] += Fe[j]*M(j,i);
        }
    }
    
    for (int m=0; m<M.size1(); m++) {
        if (m!=k) {
            for (int n=0; n<M.size1(); n++) {
                M(n,m) -= M_k[n]*R[m]/det_rat;
            }
        }
        else {
            for (int n=0; n<M.size1(); n++) {
                M(n,m) = M_k[n]/det_rat;
            }
        }
    }
    
    return;
}


template <class G, class S, class V> double det_rat_insert_anti(times & anti_segment, dense_matrix & M, S& segments_old, G& F, double BETA, double & det_rat_sign, double & overlap, V& R) {
    
    std::vector<double> F_k(R.size());
    
    typename S::iterator it=segments_old.begin();
    for (int i=0; i<F_k.size(); i++) {
        F_k[i]=interpolate_F(anti_segment.t_start()-it->t_start(), BETA, F);
        it++;
    }
    
    double det_rat = -interpolate_F(anti_segment.t_start()-anti_segment.t_end(), BETA, F);
    
    it=segments_old.begin();
    for (int i=0; i<R.size(); i++) {
        R[i]=0;
        for (int l=0; l<R.size(); l++) {
            R[i] += F_k[l]*M(l,i);
        }
        det_rat += interpolate_F(it->t_end()-anti_segment.t_end(), BETA, F)*R[i];
        it++;
    }
    
    overlap = 1;
    // take care of sign changes produced by segments which "wind around"
    // check if anti-segment winds around
    if (anti_segment.t_end()<anti_segment.t_start()) {
        det_rat *= -1;
        overlap = -1;
    }
    
    if (det_rat < 0) {
        det_rat_sign = -1;
        det_rat *= -1;
    }
    else {
        det_rat_sign = 1;
    }
    
    return det_rat;
    
}


inline int cycle(int i, int size) {
    return (i>0 ? i-1 : size-1);
}

template <class G, class S, class V> void compute_M_insert_anti(times & anti_segment, int s, int r, dense_matrix & M, S& segments_old, G& F, double BETA, double det_rat, V& R) {
    
    dense_matrix M_new(M.size1()+1,M.size1()+1);
    std::vector<double> F_kp1(R.size()), L(R.size());
    
    typename S::iterator it=segments_old.begin();
    for (int i=0; i<F_kp1.size(); i++) {
        F_kp1[i]=interpolate_F(it->t_end()-anti_segment.t_end(), BETA, F);
        it++;
    }
    
    for (int i=0; i<L.size(); i++) {
        L[i]=0;
        for (int l=0; l<L.size(); l++) {
            L[i] += M(i,l)*F_kp1[l];
        }
    }
    
    int i_new, j_new;
    int size=M.size1();
    
    // element (k+1,k)
    M_new(r,s) = -1./det_rat;
    
    if (r!=0) { // segments remain in the usual order
        
        // row k+1 and column k
        for (int i=0; i<size; i++) {
            i_new = (i<r ? i : i+1);
            j_new = (i<s ? i : i+1);
            
            M_new(i_new,s) = L[i]/det_rat;
            M_new(r,j_new) = R[i]/det_rat;
        }
        
        // remaining elements
        for (int i=0; i<size; i++) {
            i_new = (i<r ? i : i+1);
            for (int j=0; j<size; j++) {
                j_new = (j<s ? j : j+1);
                M_new(i_new,j_new) = M(i,j) - L[i]*R[j]/det_rat;
            }
        }
    }
    else { // need to permute indices of R, L, M
        
        // row k+1 and column k
        for (int i=0; i<size; i++) {
            i_new = (i<r ? i : i+1);
            j_new = (i<s ? i : i+1);
            
            M_new(i_new,s) = L[i]/det_rat;
            M_new(r,j_new) = R[cycle(i,size)]/det_rat;
        }
        
        // remaining elements
        for (int i=0; i<size; i++) {
            i_new = (i<r ? i : i+1);
            for (int j=0; j<size; j++) {
                j_new = (j<s ? j : j+1);
                M_new(i_new,j_new) = M(i,cycle(j,size)) - L[i]*R[cycle(j,size)]/det_rat;
            }
        }  
    }
    
    swap(M_new, M);
    return;
}

template <class G, class S> double det_rat_remove_anti(times anti_segment, int r, int s, dense_matrix & M, S& segments_old, G& F, double BETA, double & det_rat_sign) {
    
    // r is the index of the segment which is removed
    // s is the index of the segment which is shifted
    
    typename S::iterator it=segments_old.begin();
    typename S::iterator its(it), itr(it);
    advance(its, s); 
    advance(itr, r);
    
    double inv_det_rat = -interpolate_F(its->t_end()-itr->t_start(), BETA, F);
    
    for (int i=0; i<segments_old.size(); i++) {
        if (i!=s) {
            inv_det_rat -= interpolate_F(it->t_end()-itr->t_start(), BETA, F)*M(r,i)/M(r,s);
        }
        it++;
    }
    
    // take care of sign changes produced by segments which "wind around"
    if (anti_segment.t_end() < anti_segment.t_start()) {
        inv_det_rat *= -1;
    }
    
    if (inv_det_rat < 0) {
        det_rat_sign = -1;
        inv_det_rat *= -1;
    }
    else {
        det_rat_sign = 1;
    }
    
    return 1/inv_det_rat;
    
}


void compute_M_remove_anti(dense_matrix & M, int s, int r) {
    
    dense_matrix M_new(M.size1()-1,M.size1()-1);
    
    int i_old, j_old;
    int size=M_new.size1();
    
    if(r!=0) { // order of segments remains unchanged
        for (int i=0; i<size; i++) {
            i_old = (i<r ? i : i+1);
            for (int j=0; j<size; j++) {
                j_old = (j<s ? j : j+1);
                M_new(i,j) = M(i_old,j_old) - M(i_old, s)*M(r, j_old)/M(r, s);
            }
        }
    }
    else { // need to permute indices of M
        for (int i=0; i<size; i++) {
            for (int j=0; j<size; j++) {
                M_new(i,cycle(j,size)) = M(i+1,j) - M(i+1, s)*M(r, j)/M(r, s);
            }
        }  
    }
    
    swap(M_new, M);
    return;
}



#endif
