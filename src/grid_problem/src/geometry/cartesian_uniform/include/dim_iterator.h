#ifndef DIM_ITERATOR_
#define DIM_ITERATOR_

#include "common.h"

struct dim_iterator
{
    int dim;
    pair<voxel_point, voxel_point> vbbox;
    voxel_point voxel_coord;
    voxel_point conv_sizes;



    dim_iterator(const pair<voxel_point, voxel_point> &);

    int local_index() const;
    dim_iterator & operator++();
    bool is_overfull() const;
};

#endif