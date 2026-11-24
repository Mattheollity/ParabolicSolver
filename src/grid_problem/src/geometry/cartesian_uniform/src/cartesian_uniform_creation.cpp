#include "cartesian_uniform.h"

void cartesian_uniform::fill_domain()
{
    for(dim_iterator di({voxel_point(gc.dim, 0), gc.N}); !di.is_overfull(); ++di)
    {
        int i = di.local_index();
        point p(gc.dim);

        for(int d=0; d<gc.dim; d++)
        p[d] = gc.bbox.first[d]+di.voxel_coord[d]*h[d];

        vertexes[i] = vertex(i, std::move(p));
        find_neighbors(vertexes[i], di.voxel_coord);
    }
}

void cartesian_uniform::find_neighbors(vertex & v, const voxel_point & vc)
{
    for(int d=0; d<gc.dim; d++)
    if(vc[d]+1<gc.N[d])
    v.neighbors_indexes.push_back(v.index+csum[d]);

    for(int d=0; d<gc.dim; d++)
    if(vc[d]>0)
    v.neighbors_indexes.push_back(v.index-csum[d]);

    if(v.neighbors_indexes.size()==2*gc.dim) inner.push_back(v.index);
    else boundary.push_back(v.index);
}



int cartesian_uniform::global_index(const voxel_point & vc) const
{
    int index = 0;

    for(int d=0; d<gc.dim; d++)
    index += csum[d]*vc[d];

    return index;
}