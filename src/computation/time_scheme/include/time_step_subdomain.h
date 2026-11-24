#ifndef TIME_STEP_SUBDOMAIN_
#define TIME_STEP_SUBDOMAIN_

#include "common.h"
#include "config_structures_ts.h"
#include "time_scheme.h"

class time_step_subdomain
{
    const time_scheme_config * tsc;

    subdomain_system next_time;
    subdomain_system current_time;

    shared_ptr<SparseLU> solver;
    vector<subdomain_values> current_values;

    public:
    
    time_step_subdomain(const time_scheme_config & tsc_) : tsc(&tsc_) {};
    time_step_subdomain(const time_step_subdomain & tss_) : tsc(tss_.tsc) {};
    
    void setup_time_step_solver(const subdomain_system &);
    void setup_first_solutions(const vector<subdomain_values> & first_solutions);

    void build_schur_complement(vector<triplet> &);
    void build_schur_rhs(vector<vect> &);
    void solve(vector<subdomain_values> &);
};

#endif