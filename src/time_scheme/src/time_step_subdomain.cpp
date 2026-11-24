#include "time_step_subdomain.h"

void time_step_subdomain::setup_time_step_solver(const subdomain_system & ss)
{
    tsc->scheme->build_next_time_step_system(ss, next_time);
    tsc->scheme->build_current_time_step_system(ss, current_time);

    solver->compute(next_time.ii);
}

void time_step_subdomain::setup_first_solutions(const vector<subdomain_values> & first_solutions)
{
    current_values = first_solutions;
}

void time_step_subdomain::build_schur_complement(vector<triplet> & trp)
{
    spmat schur_complement = next_time.bb - next_time.ib.transpose()*solver->solve(next_time.ib);

    trp.reserve(schur_complement.nonZeros());

    for (int k = 0; k < schur_complement.outerSize(); ++k)
    for (spmat::InnerIterator it(schur_complement, k); it; ++it)
    trp.emplace_back(it.row(), it.col(), it.value());
}

void time_step_subdomain::build_schur_rhs(vector<vect> & boundaries)
{
    for(int i=0;i<boundaries.size();i++)
    boundaries[i] = next_time.bb - next_time.ib.transpose()*solver->solve(next_time.ib);
}

void time_step_subdomain::solve(vector<subdomain_values> & next_values)
{
    for(int i=0;i<next_values.size();i++)
    {
        vect & y0 = current_values[i].inner;
        vect & ylb0 = current_values[i].local_boundary;
        vect & ygb0 = current_values[i].global_boundary;

        vect & y1 = next_values[i].inner;
        vect & ylb1 = next_values[i].local_boundary;
        vect & ygb1 = next_values[i].global_boundary;

        vect rhs = current_time.ii*y0
                   - current_time.ib*ylb0 - next_time.ib*ylb1 
                   - current_time.bb*ygb0 - next_time.bb*ygb1;

        y1 = solver->solve(rhs);
        current_values[i] = next_values[i];
    }
}