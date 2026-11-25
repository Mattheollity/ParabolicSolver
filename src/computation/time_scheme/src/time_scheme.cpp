#include "time_scheme.h"


double weighted_time_scheme::compute_local_time_step(
    const subdomain_system &,
    const subdomain_values &
) const
{
    return T;
}



void weighted_time_scheme::build_next_time_step_system(
    const subdomain_system & ss,
    subdomain_system & tss,
    double T
) const
{
    spmat iA(ss.ii.rows(), ss.ii.cols());
    spmat iB(ss.bb.rows(), ss.bb.cols());

    iA.setIdentity();
    iB.setIdentity();

    double a = theta*T;

    tss.ii = iA + a*ss.ii;
    tss.ib = a*ss.ib;
    tss.bi = a*ss.bi;
    tss.bb = iB + a*ss.bb;
    tss.gb_subd = a*ss.gb_subd;
}



void weighted_time_scheme::build_current_time_step_system(
    const subdomain_system & ss,
    subdomain_system & tss,
    double T
) const
{
    spmat iA(ss.ii.rows(), ss.ii.cols());
    spmat iB(ss.bb.rows(), ss.bb.cols());

    iA.setIdentity();
    iB.setIdentity();

    double a = -(1.0-theta)*T;

    tss.ii = iA + a*ss.ii;
    tss.ib = a*ss.ib;
    tss.bi = a*ss.bi;
    tss.bb = iB + a*ss.bb;
    tss.gb_subd = a*ss.gb_subd;
}

void weighted_time_scheme::build_boundary_time_step_systems(
    const spmat & ss,
    spmat & current_time,
    spmat & next_time,
    double T
) const
{
    current_time = -(1.0-theta)*T*ss;
    next_time = theta*T*ss;
}