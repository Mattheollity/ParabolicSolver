#ifndef SPACE_SCHEME_
#define SPACE_SCHEME_

#include "common.h"

#include "grid.h"
#include "problem_conditions.h"
#include "file_io.h"

struct space_scheme_request
{
    const voxel_point & base;
    const voxel_point & s_point;
    const point h;
};

space_scheme_request generate_space_scheme_request(int i, int j, const interaction &);

struct space_scheme
{
    virtual double operator()(const space_scheme_request &);
};

space_scheme load_space_scheme(string SCHEME_DIR);

struct anisotrop_laplace : space_scheme
{
    vector<double> k;
    virtual double operator()(const space_scheme_request &) override;

    anisotrop_laplace(string SPACE_APPROX_DIR);
};



struct const_dirichlet
{
    double C;
    double operator()(const voxel_point &) {return C;}

    const_dirichlet(string SPACE_APPROX_DIR);
};
#endif