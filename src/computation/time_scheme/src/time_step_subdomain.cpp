#include "time_step_subdomain.h"

#include <stdexcept>


time_step_subdomain::time_step_subdomain(
    const time_scheme & scheme_,
    subdomain_system ss
) :
scheme(scheme_),
space_system(std::move(ss)),
solver(std::make_unique<SparseLU>())
{}



double time_step_subdomain::compute_local_time_step(
    const subdomain_values & current
) const
{
    return scheme.compute_local_time_step(space_system, current);
}



void time_step_subdomain::setup_time_step_solver(double T)
{
    scheme.build_current_time_step_system(space_system, current_time, T);
    scheme.build_next_time_step_system(space_system, next_time, T);

    if(next_time.ii.rows()==0) return;

    next_time.ii.makeCompressed();

    solver->compute(next_time.ii);

    if(solver->info()!=Eigen::Success)
    throw std::runtime_error("Subdomain factorization failed");
}



void time_step_subdomain::build_schur_complement(vector<triplet> & trp) const
{
    spmat schur = next_time.bb;

    if(next_time.ii.rows()!=0 && next_time.bb.rows()!=0)
    {
        spmat X = solver->solve(next_time.bi);

        schur -= next_time.ib*X;
    }

    trp.clear();
    trp.reserve(schur.nonZeros());

    for(int k=0; k<schur.outerSize(); k++)
    for(spmat::InnerIterator it(schur, k); it; ++it)
    trp.emplace_back(it.row(), it.col(), it.value());
}



vect time_step_subdomain::build_inner_rhs(
    const subdomain_values & current
) const
{
    return current_time.ii*current.inner
         + current_time.bi*current.local_boundary;
}



vect time_step_subdomain::build_schur_rhs(
    const subdomain_values & current,
    const subdomain_values & next
) const
{
    vect rhs = current_time.ib*current.inner
             + current_time.bb*current.local_boundary
             + current_time.gb_subd*current.global_boundary
             - next_time.gb_subd*next.global_boundary;

    if(next_time.ii.rows()!=0 && next_time.bb.rows()!=0)
    {
        vect inner_rhs = build_inner_rhs(current);
        vect x = solver->solve(inner_rhs);

        rhs -= next_time.ib*x;
    }

    return rhs;
}



void time_step_subdomain::solve(
    const subdomain_values & current,
    subdomain_values & next
) const
{
    if(next_time.ii.rows()==0)
    {
        next.inner.resize(0);
        return;
    }

    vect rhs = build_inner_rhs(current)
             - next_time.bi*next.local_boundary;

    next.inner = solver->solve(rhs);
}