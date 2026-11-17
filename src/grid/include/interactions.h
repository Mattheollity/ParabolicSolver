#include "common.h"
#include "problem_conditions.h"

struct vertex_interaction
{
    int output_index;
    voxel_point output_coord;

    vector<int> input_indexes;
    vector<voxel_point> input_coords;
};

struct vertex_interaction_request
{
    int output_index;
    int vertex_index;
    int input_subdomain;

    const vector<int> & input;
};



struct interaction
{
    int input_subdomain;
    int output_subdomain;

    int input_size;
    int output_size;
    
    vector<vertex_interaction> vertexes_io;
};

struct interaction_request
{
    int input_subdomain;
    int output_subdomain;

    const vector<int> & input;
    const vector<int> & output;
};



struct subdomain_interaction
{
    interaction ii;
    interaction ib;
    interaction bb;
    interaction global_gb;
};