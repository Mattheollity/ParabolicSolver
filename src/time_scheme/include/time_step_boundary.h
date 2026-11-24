#ifndef TIME_STEP_BOUNDARY_
#define TIME_STEP_BOUNDARY_

#include "config_structures_ts.h"
#include "space_matrix_assembler.h"

class time_step_boundary
{
    SparseLU solver;

    public:

    time_step_boundary();
    void build_boundary_system(const vector<triplet> &, int amount);
    void solve(vector<subdomain_values> &, const vector<vect> &);
};

#endif