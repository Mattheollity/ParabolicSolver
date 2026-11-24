#ifndef INITIAL_CONDITIONS_
#define INITIAL_CONDITIONS_

#include "common.h"

struct initial_condition
{
    virtual ~initial_condition() = default;

    virtual double operator()(const point &) const = 0;
};

struct zero_value : initial_condition
{
    double operator()(const point &) const override;
};

struct dimension_sine : initial_condition
{
    double operator()(const point &) const override;
};

#endif