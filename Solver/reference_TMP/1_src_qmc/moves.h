#ifndef ___MOVES___
#define ___MOVES___

#include "impurity.h"
#include "update.h"

// length of inserted segment (mu>0) or anti-segment (mu<0)
inline double compute_length(double r, double l_max, double mu) {

    if (mu == 0)
        return r*l_max;
    else
        return 1/mu*log(r*(exp(mu*l_max)-1)+1);
}


template <class RNG, class S> void insert_remove_full_line(RNG& rng, double mu, dense_matrix& u, double BETA, int& full_line, std::vector<S>& other_segments, std::vector<int>& other_full_line, int this_flavor) {

    int insert = (rng() < 0.5);

    if ((insert==1 && full_line==1) || (insert==0 && full_line==0)) return; // insert=1(0) means we want to insert(remove) a full line

    int FLAVOR = other_full_line.size();

    double otherlength_u=0;
    for (int i=0; i<FLAVOR; i++) {
        if (i==this_flavor) continue;

        double other_length=0;
        for (typename S::iterator it=other_segments[i].begin(); it!=other_segments[i].end(); it++)
            other_length += (it->t_end()-it->t_start()>0 ? it->t_end()-it->t_start() : it->t_end()-it->t_start()+BETA);

        if (other_full_line[i]==1)
            other_length = BETA;

        otherlength_u += other_length*u(i, this_flavor);

    }

    if (insert) { // try to insert full line
        if (log(rng()) < BETA*mu-otherlength_u)
            full_line = 1;
    }
    else { // try to remove full line
        if (log(rng()) < -BETA*mu+otherlength_u)
            full_line = 0;
    }

}


template <class RNG, class S, class G> void insert_remove_segment(RNG& rng, double t, int N, double BETA, double mu, dense_matrix& u, G& F, S& segments, dense_matrix& M, double & sign, std::vector<S>& other_segments, std::vector<int> other_full_line, int this_flavor, k_table_t& K_table) {

    double t_up; // distance to next segment up
    double t_down; // distance to next segment down
    segment_container_t::iterator s_up; // iterator of the segment up
    segment_container_t::iterator s_down; // iterator of the segment down

    if (rng()<0.5) { // try to insert a segment
        compute_intervals(t, BETA, t_up, t_down, segments,s_up, s_down);

        if (t_down>0) { // t does not lie on a segment -> it's possible to insert a new one starting from t

            double length = compute_length(rng(), t_up, 0);

            times segment_insert;
            segment_insert.set_t_start(t);
            double t_final = t + length;
            if (t_final > BETA)
                segment_insert.set_t_end(t_final-BETA);
            else
                segment_insert.set_t_end(t_final);

            double otherlength_u=0;
            int FLAVORS=other_full_line.size();
            for (int i=0; i<FLAVORS; i++) {
                if (i==this_flavor) continue;
                double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA);
                otherlength_u += other_length*u(i, this_flavor);
            }

            double nonloc=nonlocal(segment_insert.t_start(), segment_insert.t_end(), other_segments, BETA, this_flavor, K_table, 0);

            double log_prob, overlap, det_rat, det_rat_sign;
            std::vector<double> Fs(segments.size()), Fe(segments.size());

            det_rat = det_rat_up(segment_insert, M, segments, F, Fs, Fe, BETA, det_rat_sign, overlap);

            log_prob = log(BETA*t_up/(segments.size()+1)*det_rat)+mu*length-otherlength_u-nonloc;

            if (log(rng()) < log_prob) {
                int position=0;
                for (segment_container_t::iterator it=segments.begin(); it!=s_up; it++)
                    position++;
                compute_M_up(position, M, Fs, Fe, det_rat*overlap);
                sign *= det_rat_sign;
                segment_container_t::iterator sit=segments.insert(s_up, segment_insert);
                if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
            }
        }
    }

    else if (segments.size()>0) { // try to remove a segment
        int position = (int)(rng()*segments.size());
        s_down = segments.begin();
        for (int i=0; i<position; i++)
            s_down++;
        s_up=s_down;
        s_up++;
        if (s_up==segments.end())
            s_up = segments.begin();

        double length = s_down->t_end()-s_down->t_start();
        if (length < 0) length += BETA;

        double t_total = s_up->t_start()-s_down->t_start();
        if (t_total <= 0) t_total += BETA;

        times segment_remove = *s_down;

        double otherlength_u=0;
        int FLAVORS=other_full_line.size();
        for (int i=0; i<FLAVORS; i++) {
            if (i==this_flavor) continue;
            double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
            otherlength_u += other_length*u(i, this_flavor);
        }

        double nonloc=nonlocal(segment_remove.t_start(), segment_remove.t_end(), other_segments, BETA, this_flavor, K_table, 1);

        double log_prob, det_rat, det_rat_sign;

        det_rat = det_rat_down(position, M, segments, det_rat_sign);

        log_prob = log(BETA*t_total/segments.size()/det_rat)+length*mu-otherlength_u-nonloc;

        if (log(rng()) < -log_prob) {
            compute_M_down(position, M);
            sign *= det_rat_sign;
            segments.erase(s_down);
        }
    }

}


template <class RNG, class S, class G> void insert_remove_antisegment(RNG& rng, double t, int N, double BETA, double mu, dense_matrix& u, G& F, int& full_line, S& segments, dense_matrix& M, double & sign, std::vector<S>& other_segments, std::vector<int> other_full_line, int this_flavor, k_table_t& K_table) {

    double t_up; // distance to next segment up (t_start)
    double t_down; // distance to next segment down (t_end)
    segment_container_t::iterator s_up; // iterator of the segment up
    segment_container_t::iterator s_down; // iterator of the segment down

    if (rng()<0.5) { // try to insert an anti-segment

        if (full_line==1) {
            t_down = -BETA;
            double length = compute_length(rng(), BETA, 0);
            double t_end = (t+length < BETA ? t+length : t+length-BETA);
            times segment_insert(t_end, t);
            times segment_remove(t,t_end);

            double log_prob, overlap, det_rat, det_rat_sign;
            std::vector<double> Fs(segments.size()), Fe(segments.size());
            det_rat = det_rat_up(segment_insert, M, segments, F, Fs, Fe, BETA, det_rat_sign, overlap);

            double otherlength_u=0;
            int FLAVORS=other_full_line.size();
            for (int i=0; i<FLAVORS; i++) {
                if (i==this_flavor) continue;
                double other_length = compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
                otherlength_u += other_length*u(i, this_flavor);
            }

            double nonloc=nonlocal(segment_remove.t_end(), segment_remove.t_start(), other_segments, BETA, this_flavor, K_table, 0);

            log_prob = log(BETA*BETA*det_rat)-length*mu+otherlength_u-nonloc;

            if (log(rng()) < log_prob) {
                compute_M_up(0, M, Fs, Fe, det_rat*overlap);
                sign *= det_rat_sign;
                segment_container_t::iterator sit;
                sit=segments.insert(segment_insert).first;
                if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
                full_line = 0;
            }
        }

        else {
            compute_intervals(t, BETA, t_up, t_down, segments, s_up, s_down);

            if (t_down<0) { // t does lie on a segment -> it's possible to insert an anti-segment starting from t

                double length = compute_length(rng(), -t_down, 0);

                times segment_shrink(s_down->t_start(),t);

                double t_start = t + length;
                if (t_start > BETA)
                    t_start-=BETA;

                times segment_insert(t_start, s_down->t_end());
                times anti_segment(t,t_start);

                double otherlength_u=0;
                int FLAVORS=other_full_line.size();
                for (int i=0; i<FLAVORS; i++) {
                    if (i==this_flavor) continue;
                    double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
                    otherlength_u += other_length*u(i, this_flavor);
                }

                double nonloc=nonlocal(anti_segment.t_end(), anti_segment.t_start(), other_segments, BETA, this_flavor, K_table, 0);

                double log_prob, overlap, det_rat, det_rat_sign;
                std::vector<double> R(segments.size());
                det_rat = det_rat_insert_anti(anti_segment, M, segments, F, BETA, det_rat_sign, overlap, R);

                log_prob = log(BETA*(-t_down)/(segments.size()+1)*det_rat)-length*mu+otherlength_u-nonloc;

                if (log(rng()) < log_prob) {

                    int s, r; // s is the segment which is shifted, r the segment which is inserted
                    s = 0;
                    for (segment_container_t::iterator it=segments.begin(); it!=s_down; it++)
                        s++;
                    if (anti_segment.t_end() > segment_shrink.t_start())
                        r = s+1;
                    else {
                        r = 0;
                        s++;
                    }

                    compute_M_insert_anti(anti_segment, s, r, M, segments, F, BETA, det_rat*overlap, R);
                    times segment_new_endpoint(*s_down);
                    segment_container_t::iterator prev_segment=s_down;
                    if(s_down !=segments.begin()) prev_segment--;
                    else prev_segment=segments.begin();
                    segment_new_endpoint.set_t_end(t);
                    segments.erase(s_down); //erase old segment (without shifted end
                    s_down=segments.insert(segment_new_endpoint).first; //in
                    if(s_down==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
                    segment_container_t::iterator sit=segments.insert(segment_insert).first; //insert  new segment
                    if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
                    if (segment_insert.t_start()>segments.begin()->t_start()) {
                        s_down++;
                        segment_container_t::iterator sit=segments.insert(s_down, segment_insert);
                        if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
                    }
                    else {
                        segment_container_t::iterator sit=segments.insert(segments.begin(), segment_insert);
                        if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
                    }
                }
            }
        }

    }
    else if (segments.size()>1) { // try to remove an anti-segment
        int r = (int)(rng()*segments.size());
        s_up = segments.begin();
        for (int i=0; i<r; i++) s_up++;

        int s = r-1;
        if (s<0) {
            s=segments.size()-1;
            s_down=segments.end();
            s_down--;
        }
        else {
            s_down = s_up;
            s_down--;
        }

        double length = s_up->t_start() - s_down->t_end();
        if (length < 0) length += BETA;

        double t_total = s_up->t_end() - s_down->t_end();
        if (t_total < 0) t_total += BETA;

        times anti_segment(s_down->t_end(),s_up->t_start());

        double otherlength_u=0;
        int FLAVORS=other_full_line.size();
        for (int i=0; i<FLAVORS; i++) {
            if (i==this_flavor) continue;
            double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
            otherlength_u += other_length*u(i, this_flavor);
        }

        double nonloc=nonlocal(anti_segment.t_end(), anti_segment.t_start(), other_segments, BETA, this_flavor, K_table, 1);

        double log_prob, det_rat, det_rat_sign;

        det_rat = det_rat_remove_anti(anti_segment, r, s, M, segments, F, BETA, det_rat_sign);

        log_prob = log(BETA*t_total/segments.size()/det_rat)-length*mu+otherlength_u-nonloc;

        if (log(rng()) < -log_prob) {

            compute_M_remove_anti(M, s, r);

            double t_end = s_up->t_end();
            segments.erase(s_up);

            if (r>0) {
                s_up=segments.begin();
                for (int k=0; k<s; k++)
                    s_up++;
            }
            else {
                s=segments.size()-1;
                s_up = segments.end();
                s_up--;
            }
            times s_up_new(*s_up);
            s_up_new.set_t_end(t_end);
            segments.erase(s_up);
            segment_container_t::iterator sit=segments.insert(s_up_new).first;
            if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
        }
    }

    else if (segments.size()==1) {

        s_down = segments.begin();

        double det_rat = fabs(M(0,0));
        double length = s_down->t_start()-s_down->t_end();
        if (length<0) length += BETA;
        times anti_segment(s_down->t_end(),s_down->t_start());

        double otherlength_u=0;
        int FLAVORS=other_full_line.size();
        for (int i=0; i<FLAVORS; i++) {
            if (i==this_flavor) continue;
            double other_length = compute_overlap(anti_segment, other_segments[i], other_full_line[i], BETA);
            otherlength_u += other_length*u(i, this_flavor);
        }

        double nonloc=nonlocal(anti_segment.t_end(), anti_segment.t_start(), other_segments, BETA, this_flavor, K_table, 1);

        double log_prob = log(BETA*BETA/det_rat)-length*mu+otherlength_u-nonloc;

        if (log(rng()) < -log_prob) {
            full_line=1;
            segments.erase(s_down);
            compute_M_down(0,M); // attention: M.clear() sets elements to zero
        }
    }

}



// shift segment
template <class RNG, class S, class G> void shift_segment(RNG& rng, S& segments, int N, double BETA, double mu, dense_matrix& u, G& F, dense_matrix& M, double & sign, std::vector<S>& other_segments, std::vector<int>& other_full_line, int this_flavor, k_table_t& K_table) {

    int size = segments.size();

    if (size < 1) return;

    int n = (int)(size*rng());

    segment_container_t::iterator s, s_up;
    s=segments.begin();
    for (int i=0; i<n; i++) s++;
    s_up = s; s_up++;
    if (s_up == segments.end()) s_up = segments.begin();

    double interval = s_up->t_start() - s->t_start();
    if (interval <= 0) interval += BETA;

    double length = compute_length(rng(), interval, 0);
    double length_old = s->t_end()-s->t_start();
    if (length_old<0)
        length_old += BETA;

    double new_t_end = s->t_start() + length;
    if (new_t_end > BETA)
        new_t_end -= BETA;

    times segment_insert(s->t_start(), new_t_end);
    times segment_remove=*s;

    double otherlength_u=0;
    int FLAVORS=other_full_line.size();
    for (int i=0; i<FLAVORS; i++) {
        if (i==this_flavor) continue;
        double other_length = compute_overlap(segment_insert, other_segments[i], other_full_line[i], BETA)-compute_overlap(segment_remove, other_segments[i], other_full_line[i], BETA);
        otherlength_u += other_length*u(i, this_flavor);
    }

    double nonloc=nonlocal_shift(new_t_end, s->t_end(), other_segments, BETA, this_flavor, K_table);

    double det_rat, det_rat_sign, overlap;

    det_rat = det_rat_shift(segment_insert, n, M, segments, F, BETA, det_rat_sign, overlap);

    if (log(rng()) < log(det_rat)+(length-length_old)*mu-otherlength_u-nonloc) {

        compute_M_shift(segment_insert, n, M, segments, F, BETA, det_rat*overlap);
        sign *= det_rat_sign;
        times s_new(*s);
        s_new.set_t_end(new_t_end);
        segments.erase(s);
        segment_container_t::iterator sit;
        sit=segments.insert(s_new).first;
        if(sit==segments.end()){std::cerr<<"segment could not be inserted! exiting."<<std::endl;}
    }

}

// swap segment configurations
template <class RNG, class S, class G> void swap_segments(RNG& rng, double BETA, G& F_up, G& F_down, S& segments_up, S& segments_down, int& full_line_up, int& full_line_down, double & sign_up, double & sign_down, dense_matrix& M_up, dense_matrix& M_down) {

    // ATTENTION!!!! IN THE GENERAL CASE THIS NEEDS IMPLEMENTATION OF OVERLAP PART!!!!!!

    dense_matrix M_new_up, M_new_down;

    // before swap
    double det_old_up = construct_inverse(M_new_up, segments_up, BETA,  F_up); // here M_new_up is just a dummy
    double det_old_down = construct_inverse(M_new_down, segments_down, BETA,  F_down); // here M_new_down is just a dummy
    // before swap
    double det_new_up = construct_inverse(M_new_up, segments_down, BETA,  F_up);
    double det_new_down = construct_inverse(M_new_down, segments_up, BETA,  F_down);

    double det_rat = (det_new_up/det_old_up)*(det_new_down/det_old_down);

    // length of segments, overlap and phonon part are not changed
    if (rng() < fabs(det_rat)) {
        /*
         swap(M_new_up, M_up);
         swap(M_new_down, M_down);
         swap(segments_up, segments_down);
         double dummy=full_line_up;
         full_line_down=full_line_up;
         full_line_up=dummy;
         dummy=sign_up;
         sign_down=sign_up;
         sign_up=dummy;
         */

        swap(M_new_up, M_up);
        swap(M_new_down, M_down);
        swap(segments_up, segments_down);
        std::swap(full_line_up, full_line_down);
        std::swap(sign_up, sign_down);

    }

}

/*
 // HH swap ALL SPINS -> need not recompute overlap
 template <class RNG, class S, class G> void swap_spins(RNG& rng, double BETA, int FLAVORS, G& F_vec,
 S& segments_vec, std::vector<int> &full_line_vec,
 std::vector<double>& sign_vec, std::vector<dense_matrix> & M_vec){

	std::vector<dense_matrix> M_new_vec(FLAVORS);
	double det_rat=1.0;

	for(int j=0; j<FLAVORS/2; j++){//for all orbitals swap spins //2j=orb j up, 2j+1=orb j down
 // before swap
 double det_old_up = construct_inverse(M_new_vec[2*j], segments_vec[2*j], BETA, F_vec[2*j]);
 double det_old_down = construct_inverse(M_new_vec[2*j+1], segments_vec[2*j+1], BETA,  F_vec[2*j+1]); // here M_new is just a dummy
 // after swap
 double det_new_up = construct_inverse(M_new_vec[2*j], segments_vec[2*j+1], BETA,  F_vec[2*j]);
 double det_new_down = construct_inverse(M_new_vec[2*j+1], segments_vec[2*j], BETA,  F_vec[2*j+1]);

 det_rat *= (det_new_up/det_old_up)*(det_new_down/det_old_down);
	}
	// length of segments, overlap and phonon part are not changed
	//    std::cout << "det_rat=" << det_rat << std::endl << std::flush;
	if (rng() < fabs(det_rat)) {
 //    std::cout << "success\n";
 for(int j=0; j<FLAVORS/2; j++){
 swap(M_new_vec[2*j], M_vec[2*j]);
 swap(M_new_vec[2*j+1], M_vec[2*j+1]);
 swap(segments_vec[2*j], segments_vec[2*j+1]);
 std::swap(full_line_vec[2*j], full_line_vec[2*j+1]);
 std::swap(sign_vec[2*j], sign_vec[2*j+1]);
 }

	}

 }
 */

/*
 // flip segment between up/down
 template <class RNG, class S, class G> void flip_segment(RNG& rng, S& segments, int N, double BETA, dense_matrix& M, double & sign, double & other_sign, G& other_F, dense_matrix& other_M, S& other_segments, int other_full_line) {

 int size = segments.size();

 if (size < 1) return;

 int position = size*rng();

 typename S::iterator s;
 s=segments.begin();
 for (int i=0; i<position; i++) s++;

 if (compute_overlap(*s, other_segments, other_full_line, BETA)>0)
 return;


 // no overlap -> can flip segment (mu and u contributions don't change)

 double det_rat_remove, det_rat_remove_sign, det_rat_insert, det_rat_insert_sign, overlap;

 det_rat_remove = det_rat_down(position, M, segments, det_rat_remove_sign);

 std::vector<double> other_Fs(other_segments.size()), other_Fe(other_segments.size());
 det_rat_insert = det_rat_up(*s, other_M, other_segments, other_F, other_Fs, other_Fe, BETA, det_rat_insert_sign, overlap);

 if (rng() < det_rat_remove*det_rat_insert*segments.size()/(other_segments.size()+1)) {

	typename S::iterator s_up, s_down;
	double t_up, t_down;
	compute_intervals(s->t_start(), BETA, t_up, t_down, other_segments, s_up, s_down);

	int n=0;
	for (typename S::iterator it=other_segments.begin(); it!=s_up; it++)
 n++;
	compute_M_up(*s, n, other_M, other_segments, other_F, other_Fs, other_Fe, BETA, det_rat_insert*overlap);
	other_segments.insert(s_up, *s);

	compute_M_down(position, M);
	segments.erase(s);

	sign *= det_rat_remove_sign;
	other_sign *= det_rat_insert_sign;

 }

 }
 */

#endif
