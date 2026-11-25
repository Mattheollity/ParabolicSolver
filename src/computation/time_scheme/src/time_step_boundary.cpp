#include "time_step_boundary.h"

#include <stdexcept>


time_step_boundary::time_step_boundary(
    const time_scheme & scheme_,
    const list<triplet> & trp,
    const vector<int> & boundary_sizes
) :
scheme(scheme_)
{
    csum = {0};

    for(int size : boundary_sizes)
    csum.push_back(csum.back()+size);

    boundary_space_system.resize(csum.back(), csum.back());

    boundary_space_system.setFromTriplets(trp.begin(), trp.end());
}



void time_step_boundary::build_boundary_system(
    const vector<vector<triplet>> & schur,
    double T
)
{
    spmat boundary_next_time;

    scheme.build_boundary_time_step_systems(
        boundary_space_system,
        boundary_current_time,
        boundary_next_time,
        T
    );

    if(csum.back()==0) return;

    size_t amount = 0;

    for(const auto & s : schur)
    amount += s.size();

    vector<triplet> trp;
    trp.reserve(amount);

    for(int i=0; i<schur.size(); i++)
    for(const triplet & t : schur[i])
    trp.emplace_back(
        t.row()+csum[i],
        t.col()+csum[i],
        t.value()
    );

    spmat S(csum.back(), csum.back());

    S.setFromTriplets(trp.begin(), trp.end());

    S += boundary_next_time;
    S.makeCompressed();

    solver.compute(S);

    if(solver.info()!=Eigen::Success)
    throw std::runtime_error("Boundary factorization failed");
}



void time_step_boundary::solve(
    const vector<subdomain_values> & current,
    vector<subdomain_values> & next,
    const vector<vect> & rhs
) const
{
    vect b(csum.back());
    vect u(csum.back());

    for(int i=0; i<rhs.size(); i++)
    {
        int size = csum[i+1]-csum[i];

        b.segment(csum[i], size) = rhs[i];
        u.segment(csum[i], size) = current[i].local_boundary;
    }

    if(csum.back()!=0)
    {
        b += boundary_current_time*u;
        u = solver.solve(b);
    }

    for(int i=0; i<next.size(); i++)
    next[i].local_boundary = u.segment(
        csum[i],
        csum[i+1]-csum[i]
    );
}