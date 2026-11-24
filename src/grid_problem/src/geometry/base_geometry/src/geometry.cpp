#include "geometry.h"

#include <algorithm>
#include <stdexcept>

void geometry::build_permutation_vector()
{
    vector<int> result;
    result.reserve(inner.size());

    vector<bool> used(inner.size(), false);

    for(const subdomain & s : subdomains)
    for(int i : s.inner)
    result.push_back(find_local_index(i, inner));

    for(const subdomain & s : subdomains)
    for(int i : s.boundary)
    result.push_back(find_local_index(i, inner));

    permutation = std::move(result);
}



int geometry::find_local_index(int i, const vector<int> & local) const
{
    auto it = std::lower_bound
    (
        local.begin(),
        local.end(),
        i
    );

    if(it==local.end() || *it!=i)
        throw std::runtime_error("Vertex index was not found");

    return it - local.begin();
}



point geometry::compute_abs_h(const point & a, const point & b) const
{
    if(a.size()!=b.size())
        throw std::runtime_error("Points have different dimensions");

    point h(a.size());

    for(std::size_t d=0; d<a.size(); d++)
        h[d] = std::abs(a[d]-b[d]);

    return h;
}