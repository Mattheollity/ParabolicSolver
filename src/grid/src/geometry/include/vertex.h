#ifndef VERTEX_
#define VERTEX_

#include "common.h"
#include "config_structures.h"

struct vertex
{
    int dim;
    int index;
    int subdomain;
    voxel_point voxel_coord;
    point real_coord;

    vector<int> neighbors_indexes;
    vector<int> neighbors_dims;

    vertex() {};
    vertex(int index_, voxel_point vp, point p);

    string string_export() const;
    vertex & string_import(string);
};



#endif