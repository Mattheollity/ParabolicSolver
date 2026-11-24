#include "initial_conditions.h"

double zero_value::operator()(const point &) const
{
    return 0.0;
}

double dimension_sine::operator()(const point & p) const
{
    double result = 1.0;

    for(double d : p)
    result *= std::sin(M_PI * d);

    return result;
}