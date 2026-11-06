#include "dim_iterator.h"

dim_iterator::dim_iterator(const pair<voxel_point, voxel_point> & bbox_) : 
    vbbox(bbox_),
    dim(bbox_.first.size()),
    voxel_coord(vbbox.first),
    conv_sizes(dim, 1)
{
    for(unsigned d=1;d<dim;d++)
        conv_sizes[d]=conv_sizes[d-1]*(vbbox.second[d-1]-vbbox.first[d-1]);
}

dim_iterator & dim_iterator::operator++() {
    voxel_coord[0]++;
    for(unsigned d=0; d<dim-1;d++)
    {
        if(voxel_coord[d]==vbbox.second[d])
        {
            voxel_coord[d]=vbbox.first[d];
            voxel_coord[d+1]++;
        }
    }
    return *this;
}

bool dim_iterator::is_overfull()
{
    return voxel_coord[dim-1]>=vbbox.second[dim-1];
}

int dim_iterator::local_index() const
{
    unsigned s=0;

    for(unsigned d=0;d<dim;d++)
    s+=conv_sizes[d]*voxel_coord[d];

    return s;
}

pair<vector<int>, vector<int>> dim_iterator::find_neighbors_and_dims() const
{
    vector<int> neighbors;
    vector<int> dims;

    int s = local_index();

    for(unsigned d=0;d<dim;d++)
    if(voxel_coord[d]+1<vbbox.second[d])
    {
        neighbors.push_back(s+conv_sizes[d]);
        dims.push_back(d);
    }

    for(unsigned d=0;d<dim;d++)
    if(voxel_coord[d]-1>=vbbox.first[d])
    {
        neighbors.push_back(s-conv_sizes[d]);
        dims.push_back(-d);
    }

    return {neighbors,dims};
}