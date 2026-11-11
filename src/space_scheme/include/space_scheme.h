#ifndef SPACE_SCHEME_
#define SPACE_SCHEME_

#include "common.h"

#include "problem_conditions.h"
#include "file_io.h"

struct space_scheme
{
    virtual double operator()(const voxel_point & base, const voxel_point & s_point);
};

space_scheme load_space_scheme(string SCHEME_DIR);

struct anisotrop_laplace : space_scheme
{
    vector<double> k;
    vector<double> h;
    virtual double operator()(const vector<voxel_point &>) override;

    anisotrop_laplace(string SPACE_APPROX_DIR);
};



struct const_dirichlet
{
    double C;
    double operator()(const voxel_point &) {return C;}

    const_dirichlet(string SPACE_APPROX_DIR);
};
#endif