#include "time_step_boundary.h"

void time_step_boundary::build_boundary_system(const vector<triplet> & trp, int amount)
{
    spmat B(amount, amount);
    B.setFromTriplets(trp.begin(), trp.end());
    solver.compute(B);
}

void time_step_boundary::solve(vector<subdomain_values> & vsv, const vector<vect> & rhs)
{
    vect concate_rhs;
    concate_vector(rhs, concate_rhs);

    vect concate_result = solver.solve(concate_rhs);
    vector<vect> result(rhs);
    split_vector(result, concate_result);
    
    for(int i=0;i<result.size();i++)
    vsv[i].local_boundary=result[i];
}