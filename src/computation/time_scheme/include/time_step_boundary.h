#ifndef TIME_STEP_BOUNDARY_
#define TIME_STEP_BOUNDARY_

#include "common.h"
#include "time_scheme.h"


class time_step_boundary
{
    const time_scheme & scheme;

    spmat boundary_space_system;
    spmat boundary_current_time;

    vector<int> csum;

    SparseLU solver;



    public:

    time_step_boundary(
        const time_scheme & scheme_,
        const list<triplet> & trp,
        const vector<int> & boundary_sizes
    );



    void build_boundary_system(
        const vector<vector<triplet>> & schur,
        double T
    );



    void solve(
        const vector<subdomain_values> & current,
        vector<subdomain_values> & next,
        const vector<vect> & rhs
    ) const;
};

#endif