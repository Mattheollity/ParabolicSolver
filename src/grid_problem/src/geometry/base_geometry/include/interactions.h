#ifndef INTERACTIONS_
#define INTERACTIONS_

#include "common.h"

struct vertex_interaction
{
    int output_index;
    point output_coord;

    vector<point> interacting_coords;

    vector<int> input_indexes;
    vector<int> input_coord_indexes;
};

struct interaction
{
    int input_size;
    int output_size;

    vector<vertex_interaction> vertexes_io;
};

struct subdomain_interaction
{
    interaction ii;
    interaction ib;
    interaction bb;
    interaction global_gb;
};

#endif