#ifndef BOUNDARY_CONDITIONS_
#define BOUNDARY_CONDITIONS_

#include "common.h"

struct dirichlet
{
    virtual ~dirichlet() = default;

    virtual double operator()(const point &, double t) const = 0;
};

struct const_dirichlet : dirichlet
{
    double C;

    const_dirichlet(double C);

    double operator()(const point &, double t) const override;
};

#endif