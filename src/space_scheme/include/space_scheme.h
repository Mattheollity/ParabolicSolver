#ifndef SPACE_SCHEME_
#define SPACE_SCHEME_

#include "common.h"
#include "interactions.h"

struct space_scheme
{
    virtual double operator()(const space_scheme_request &) const = 0;
};

struct space_scheme_request
{
    const voxel_point & base;
    const voxel_point & sub;
    point h;
    point compute_abs_h(const point & a, const point & b);
    space_scheme_request(int i, int j, const interaction &);
};

struct anisotrop_laplace : space_scheme
{
    vector<double> k;
    virtual double operator()(const space_scheme_request &) const override;
};

#endif