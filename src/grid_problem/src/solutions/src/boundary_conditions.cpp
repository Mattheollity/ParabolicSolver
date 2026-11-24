#include "boundary_conditions.h"

const_dirichlet::const_dirichlet(double C) : C(C) {}

double const_dirichlet::operator()(const point &, double) const
{
    return C;
}