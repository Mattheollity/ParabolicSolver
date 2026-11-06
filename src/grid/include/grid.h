#ifndef GRID_
#define GRID_

#include "common.h"

#include "vertex.h"
#include "dim_iterator.h"
#include "problem_conditions.h"
#include "file_io.h"
#include <vector>

struct vertex;
struct subdomain;
struct interaction;
struct vertex_interactions;
struct interaction_request;
struct vertex_interactions_request;
class grid;

class grid
{
    int dim;
    pair<point, point> bbox;
    voxel_point N_bbox;
    vector<int> csum;
    vector<double> h;
    
    vector<vertex> vertexes;
    vector<int> inner;
    vector<int> boundary;

    vector<subdomain> subdomains;

    vector<int> permutation;

    string GRID_DIR;


    void create_grid(const grid_conds &);
    void fill_domain();
    void find_neighbors_and_axis(vertex &);
    int global_index(const voxel_point &);

    void create_subdomain(int i, pair<point, point> bbox_);
    void fill_subdomain(subdomain &);
    void setup_subdomain_N_bbox(subdomain &, pair<point, point> &);
    void find_inner_boundary_at_subdomains();
    int amount_neighbors_at_subdomain(const vertex &) const;
    
    void assert_all_vertexes_has_subdomain() const;
    
    subdomain import_subdomain(string filename);
    subdomain export_subdomain(string filename);
    vector<int> build_permutation_vector();

    bool is_boundary(const vertex &) const;
    bool is_boundary(const subdomain &, const vertex &) const;
    point to_real_coord(voxel_point) const;
    voxel_point to_voxel_coord(point) const;

    bool import_grid();
    bool import_domain();
    bool import_subdomains();

    void export_grid() const;
    void export_domain() const;
    void export_subdomains() const;



    public:

    grid(const grid_conds &);
    grid(const string GRID_FILENAME);

    vector<interaction> get_interaction_between(int i, int j) const;
    vector<interaction> get_self_interaction(int i) const;
    interaction create_interaction(const interaction_request &) const;
    vertex_interactions create_vertex_interactions(const vertex_interactions_request &) const;

    interaction get_inner_boundary_interaction() const;



    int get_dim() const {return dim;}
    const pair<point, point> & get_bbox() const {return bbox;}
    const voxel_point & get_vbbox() const {return N_bbox;}
    const vector<double> & get_h() const {return h;}
    const vector<point> get_boundary_points() const;
    const vector<point> get_inner_points() const;
    const vector<int> get_permutation_vector() const {return permutation;}
    
};


struct subdomain
{
    int index;
    pair<voxel_point, voxel_point> N_bbox;
    vector<int> inner;
    vector<int> boundary;

    set<int> neighboring_subdomains;
};



struct vertex_interactions
{
    int output_index;
    voxel_point output_coord;

    vector<int> input_indexes;
    vector<voxel_point> input_coords;
};

struct interaction
{
    int input_subdomain;
    int output_subdomain;

    int input_size;
    int output_size;
    
    vector<vertex_interactions> vertexes_io;
};

struct vertex_interactions_request
{
    int output_index;
    int vertex_index;
    int input_subdomain;

    const vector<int> & input;
};

struct interaction_request
{
    int input_subdomain;
    int output_subdomain;

    const vector<int> & input;
    const vector<int> & output;
};






#endif
