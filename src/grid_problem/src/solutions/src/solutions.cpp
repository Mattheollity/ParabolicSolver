#include "solutions.h"

double zero_value(const point &)
{
    return 0;
}

double dimension_sine(const point & p)
{
    double result = 1.0;
    for(double d : p)
    result*=sin(M_PI*d);

    return result;
}