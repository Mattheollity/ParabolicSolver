#ifndef TIME_STEP_SUBDOMAIN_
#define TIME_STEP_SUBDOMAIN_

#include "common.h"
#include "time_scheme.h"


class time_step_subdomain
{
    const time_scheme & scheme;

    subdomain_system space_system;
    subdomain_system current_time;
    subdomain_system next_time;

    unique_ptr<SparseLU> solver;

    vect build_inner_rhs(
        const subdomain_values & current
    ) const;



    public:

    time_step_subdomain(const time_scheme &, subdomain_system);

    double compute_local_time_step(const subdomain_values &) const;

    void setup_time_step_solver(double T);

    void build_schur_complement(vector<triplet> &) const;

    vect build_schur_rhs(
        const subdomain_values & current,
        const subdomain_values & next
    ) const;

    void solve(
        const subdomain_values & current,
        subdomain_values & next
    ) const;
};

#endif