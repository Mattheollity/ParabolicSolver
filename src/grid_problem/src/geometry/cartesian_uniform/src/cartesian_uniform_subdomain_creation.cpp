#include "cartesian_uniform.h"

#include <algorithm>
#include <stdexcept>

void cartesian_uniform::setup_subdomains()
{
    for(int i=0; i<gc.subdomains.size(); i++)
    create_subdomain(i, gc.subdomains[i]);

    assert_all_vertexes_has_subdomain();
    find_inner_and_boundary_at_subdomains();
}

void cartesian_uniform::create_subdomain(int i, const pair<point, point> & bbox_)
{
    subdomain & s = subdomains[i];
    s.index = i;

    mark_subdomains_vertexes(s, setup_subdomain_N_bbox(bbox_));
}

pair<voxel_point, voxel_point> cartesian_uniform::setup_subdomain_N_bbox
(const pair<point, point> & bbox_) const
{
    const point & a = bbox_.first;
    const point & b = bbox_.second;

    if(a.size()!=gc.dim || b.size()!=gc.dim)
    throw std::runtime_error("Invalid subdomain dimensions");

    voxel_point va(gc.dim);
    voxel_point vb(gc.dim);

    for(int d=0; d<gc.dim; d++)
    {
        if(!std::isfinite(a[d]) || !std::isfinite(b[d]) || a[d]>b[d])
        throw std::runtime_error("Invalid subdomain bounds");

        double A = std::clamp(a[d], gc.bbox.first[d], gc.bbox.second[d]);
        double B = std::clamp(b[d], gc.bbox.first[d], gc.bbox.second[d]);

        va[d] = std::clamp(int(std::round((A-gc.bbox.first[d])/h[d])), 1, gc.N[d]-1);
        vb[d] = std::clamp(int(std::round((B-gc.bbox.first[d])/h[d])), 1, gc.N[d]-1);
    }

    return {va, vb};
}

void cartesian_uniform::mark_subdomains_vertexes
(
    subdomain & s,
    const pair<voxel_point, voxel_point> & N_bbox
)
{
    for(dim_iterator di(N_bbox); !di.is_overfull(); ++di)
    {
        vertex & v = vertexes[global_index(di.voxel_coord)];

        if(v.subdomain!=-1)
        throw std::runtime_error("Subdomains have common vertices");

        v.subdomain = s.index;
    }
}

void cartesian_uniform::find_inner_and_boundary_at_subdomains()
{
    for(int i : inner)
    {
        const vertex & v = vertexes[i];
        subdomain & s = subdomains[v.subdomain];

        if(amount_neighbors_at_subdomain(v)==v.neighbors_indexes.size())
        s.inner.push_back(i);
        else s.boundary.push_back(i);
    }

    for(int i : boundary)
    {
        set<int> adjacent_subdomains;

        for(int ni : vertexes[i].neighbors_indexes)
        if(vertexes[ni].subdomain!=-1)
        adjacent_subdomains.insert(vertexes[ni].subdomain);

        for(int s : adjacent_subdomains)
        subdomains[s].global_boundary.push_back(i);
    }
}

int cartesian_uniform::amount_neighbors_at_subdomain(const vertex & v) const
{
    int amount = 0;

    for(int ni : v.neighbors_indexes)
    if(vertexes[ni].subdomain==v.subdomain)
    amount++;

    return amount;
}

void cartesian_uniform::assert_all_vertexes_has_subdomain() const
{
    for(int i : inner)
    if(vertexes[i].subdomain==-1)
    throw std::runtime_error("Internal vertex has no subdomain");
}