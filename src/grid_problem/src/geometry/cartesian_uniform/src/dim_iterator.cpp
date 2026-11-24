#include "dim_iterator.h"

dim_iterator::dim_iterator(const pair<voxel_point, voxel_point> & bbox_) :
    dim(bbox_.first.size()),
    vbbox(bbox_),
    voxel_coord(vbbox.first),
    conv_sizes(dim, 1)
{
    for(int d=1; d<dim; d++)
    conv_sizes[d] = conv_sizes[d-1]*(vbbox.second[d-1]-vbbox.first[d-1]);
}



int dim_iterator::local_index() const
{
    int index = 0;

    for(int d=0; d<dim; d++)
    index += conv_sizes[d]*(voxel_coord[d]-vbbox.first[d]);

    return index;
}



dim_iterator & dim_iterator::operator++()
{
    voxel_coord[0]++;

    for(int d=0; d<dim-1; d++)
    if(voxel_coord[d]==vbbox.second[d])
    {
        voxel_coord[d] = vbbox.first[d];
        voxel_coord[d+1]++;
    }

    return *this;
}



bool dim_iterator::is_overfull() const
{
    if(dim==0) return true;

    for(int d=0; d<dim; d++)
    if(vbbox.first[d]>=vbbox.second[d]) return true;

    return voxel_coord[dim-1]>=vbbox.second[dim-1];
}