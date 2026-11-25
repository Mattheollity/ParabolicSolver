#include "space_scheme.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>

point space_scheme::compute_h(const point & base, const point & sub) const
{
    point h(base.size());

    for(int d=0; d<base.size(); d++)
    h[d] = sub[d]-base[d];

    return h;
}



vector<double> anisotrop_laplace::operator()(const vertex_interaction & vi) const
{
    auto base = std::find(vi.interacting_coords.begin(),
                          vi.interacting_coords.end(),
                          vi.output_coord);

    if(base==vi.interacting_coords.end())
    throw std::invalid_argument("Base vertex is missing");

    int base_index = base-vi.interacting_coords.begin();

    vector<double> H(k.size(), 0.0);
    vector<double> q(k.size(), 0.0);

    vector<double> coefficients(vi.interacting_coords.size(), 0.0);

    for(const point & p : vi.interacting_coords)
    {
        point h = compute_h(vi.output_coord, p);

        for(int d=0; d<k.size(); d++)
        H[d] += std::abs(h[d]);
    }

    for(int i=0; i<vi.interacting_coords.size(); i++)
    {
        if(i==base_index) continue;

        point h = compute_h(vi.output_coord, vi.interacting_coords[i]);

        for(int d=0; d<k.size(); d++)
        q[d] = h[d]!=0.0 ? 2.0/(std::abs(h[d])*H[d]) : 0.0;

        double a = -std::inner_product(k.begin(), k.end(), q.begin(), 0.0);

        coefficients[i] += a;
        coefficients[base_index] -= a;
    }

    vector<double> result;
    result.reserve(vi.input_coord_indexes.size());

    for(int i : vi.input_coord_indexes)
    result.push_back(coefficients[i]);

    return result;
}