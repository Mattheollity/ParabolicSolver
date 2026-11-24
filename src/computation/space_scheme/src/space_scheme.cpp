#include "space_scheme.h"

space_scheme_request::space_scheme_request(int i, int j, const interaction & interact) :
    base(interact.vertexes_io[i].output_coord),
    sub(interact.vertexes_io[i].input_coords[j])
{
    h = compute_abs_h(
        interact.vertexes_io[i].output_real_coord,
        interact.vertexes_io[i].input_real_coords[j]
    );
}

double anisotrop_laplace::operator()(const space_scheme_request & ssr) const
{
    if(k.size()!=ssr.h.size()) throw simple_error("dimensions are uncompatible");

    double coeff = 0;

    int dim = k.size();
    auto & h = ssr.h;

    if(ssr.base==ssr.sub)
    {
        for(int d=0;d<dim;d++)
        coeff -= 2*k[d]/(h[d]*h[d]);
    }
    else
    {
        for(int d=0;d<dim;d++)
        if(ssr.base[d]!=ssr.sub[d])
        coeff=k[dim]/(std::pow(((ssr.base[d]-ssr.sub[d])*h[d]),2));
    }

    return -coeff;
}

// vect y0(const grid & gr, vect & v)
// {
//     const vector<vertex> & vtx = gr.get_vertexes();
//     auto inner = gr.get_inner();
//     v = vect(inner.size());
//     for(int i=0;i<inner.size(); i++)
//     v.coeffRef(i)=trip_sin(
//         gr.to_real_coord(vtx[inner[i]].voxel_coord)
//     );
//     return v;
// }

// double real_solution(int dim, const point p, double t)
// {
//     double k = 10101;

//     return exp(-M_PI*M_PI*k*t)*trip_sin(p);
// }

// vect get_solution_at_time_step(const grid & g, const time_integr_conds & tic, int N)
// {
//     auto points = g.get_inner_points();
//     double time = tic.time.first+tic.T*N;

//     vect sol(points.size());

//     for(int i=0; i<points.size(); i++)
//     sol.coeffRef(i) = real_solution(g.get_dim(), points[i], time);

//     return sol;
// }



// vect get_time_deriv_at_time_step(const grid & g, const time_integr_conds & tic, int N)
// {
//     vect td = (-g.get_dim())*M_PI*M_PI*get_solution_at_time_step(g, tic, N);

//     return td;
// }

// vect get_dummy_time_deriv_at_time_step(const grid & g, const time_integr_conds & tic, int N)
// {
//     vect a = get_solution_at_time_step(g, tic, N);
//     vect b = get_solution_at_time_step(g, tic, N+1);
//     vect td = (1.0/tic.T)*(b-a);

//     return td;
// }

// vect get_laplace_at_time_step(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N)
// {
//     auto points = g.get_inner_points();
//     vect yn = get_solution_at_time_step(g, tic, N);
//     vect ynn = yn;
//     ynn.setZero();

//     const mat P = p_system.permutation;
//     yn = P*yn;

//     vector<int> csum = {0};
//     for(int i : p_system.sizes)
//     csum.push_back(csum.back()+i);

//     for(int i=0; i<p_system.amount_of_parts; i++)
//     for(int j=0; j<p_system.amount_of_parts; j++)
//     ynn.block(csum[i], 0, csum[i+1]-csum[i], 1) += p_system.system[i][j]*yn.block(csum[j], 0, csum[j+1]-csum[j], 1);

//     return P.transpose()*ynn;
// }

// double get_laplace_eqlid_at_time_step(const grid & g, const time_integr_conds & tic, const sparse_system & p_system, int N)
// {
//     vect lpv = get_laplace_at_time_step(g, tic, p_system, N);
//     double norm = 0;
    
//     for(int i=0;i<lpv.rows(); i++)
//     norm+=lpv.coeff(i)*lpv.coeff(i);

//     norm/=lpv.rows();

//     return sqrt(norm);

// }
