#ifndef DIM_ITERATOR_
#define DIM_ITERATOR_

#include "common.h"
#include "problem_conditions.h"

struct dim_iterator
{
    int dim;
    pair<voxel_point, voxel_point> vbbox;
    voxel_point voxel_coord;
    voxel_point conv_sizes;
    


    dim_iterator(const pair<voxel_point, voxel_point> & vbbox);

    pair<vector<int>, vector<int>> find_neighbors_and_dims();
    dim_iterator & operator++();
    bool is_overfull();
};

#endif