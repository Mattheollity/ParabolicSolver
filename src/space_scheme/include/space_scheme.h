#ifndef SPACE_SCHEME_
#define SPACE_SCHEME_

#include "common.h"
#include "config_structures_ss.h"
#include "interactions.h"

struct space_scheme_request
{
    const voxel_point & base;
    const voxel_point & s_point;
    point h;
    point compute_abs_h(const point & a, const point & b);
};

space_scheme_request ss_request(int i, int j, const interaction &);

struct anisotrop_laplace : space_scheme
{
    vector<double> k;
    virtual double operator()(const space_scheme_request &) const override;

    anisotrop_laplace(string CASE_DIR);
};



struct const_dirichlet
{
    double C;
    double operator()(const voxel_point &) {return C;}

    const_dirichlet(string SPACE_APPROX_DIR);
};
#endif