#ifndef SOLUTIONS_
#define SOLUTIONS_

#include "common.h"

double zero_value(const point &);
double dimension_sine(const point &);

struct dirichlet
{
    virtual double operator()(const point &, double t) const = 0;
};

struct const_dirichlet : dirichlet
{
    double C;
    double operator()(const point &, double t) {return C;}
};


#endif