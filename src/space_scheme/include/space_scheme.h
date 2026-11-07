#ifndef SPACE_SCHEME_
#define SPACE_SCHEME_

#include "common.h"

#include "problem_conditions.h"
#include "file_io.h"

struct anisotrop_laplace
{
    vector<double> k;
    vector<double> h;
    double operator()(const voxel_point & a, const voxel_point & b);

    anisotrop_laplace(string SPACE_APPROX_DIR);
};



struct const_dirichlet
{
    double C;
    double operator()(const voxel_point &) {return C;}

    const_dirichlet(string SPACE_APPROX_DIR);
};



double real_solution(int dim, const point p, double t);
vect y0(const grid & gr, vect & v);

double approx_func(const approx &, const vertex & main, const vertex & neighbor);
double dirich_boundary_conds(point, double t);
vect get_solution_at_time_step(const grid &, const time_integr_conds & tic, int N);
vect get_time_deriv_at_time_step(const grid &, const time_integr_conds & tic, int N);
vect get_dummy_time_deriv_at_time_step(const grid &, const time_integr_conds & tic, int N);
vect get_laplace_at_time_step(const grid &, const time_integr_conds & tic, const sparse_system & p_system, int N);
double get_laplace_eqlid_at_time_step(const grid &, const time_integr_conds & tic, const sparse_system & p_system, int N);

#endif