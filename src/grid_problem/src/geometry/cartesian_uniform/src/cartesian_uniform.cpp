#include "cartesian_uniform.h"

#include <limits>
#include <stdexcept>

cartesian_uniform::cartesian_uniform(const grid_config & gc_) :
    gc(gc_)
{
    setup_geometry();
}



void cartesian_uniform::setup_geometry()
{
    if(gc.dim<=0 || gc.N.size()!=gc.dim ||
        gc.bbox.first.size()!=gc.dim || gc.bbox.second.size()!=gc.dim)
    throw std::runtime_error("Invalid grid dimensions");

    h.resize(gc.dim);
    csum = vector<int>(gc.dim, 1);

    int amount = 1;

    for(int d=0; d<gc.dim; d++)
    {
        if(gc.N[d]<2 || !(gc.bbox.second[d]>gc.bbox.first[d]))
        throw std::runtime_error("Invalid grid parameters");

        if(amount>std::numeric_limits<int>::max()/gc.N[d])
        throw std::runtime_error("Grid is too large");

        csum[d] = amount;
        amount *= gc.N[d];
        h[d] = (gc.bbox.second[d]-gc.bbox.first[d])/(gc.N[d]-1);

        if(!std::isfinite(h[d]) || h[d]<=0)
        throw std::runtime_error("Invalid grid step");
    }

    vertexes = vector<vertex>(amount);
    subdomains = vector<subdomain>(gc.subdomains.size());

    inner.clear();
    boundary.clear();
    permutation.clear();

    if(gc.try_to_import) return;

    fill_domain();
    setup_subdomains();
    build_permutation_vector();
}