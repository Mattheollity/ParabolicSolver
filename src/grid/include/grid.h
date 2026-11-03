#ifndef GRID_
#define GRID_

#include "common.h"
#include "vertex.h"
#include "dim_iterator.h"
#include "problem_conditions.h"
#include <ostream>
#include <sstream>

struct vertex;
struct subdomain;
struct interaction;
class grid;

class grid
{
    int dim;
    pair<point, point> bbox;
    pair<voxel_point, voxel_point> N_bbox;
    vector<double> h;
    
    vector<vertex> vertexes;
    vector<subdomain> subdomains;

    vector<int> permutation;



    void fill_domain();

    subdomain & create_subdomain(int i, pair<point, point> bbox_);
    void fill_subdomain(subdomain &);
    void truncate_subdomain_bbox(pair<point, point> &);
    
    void connect_neighboring_subdomains();
    void assert_all_vertexes_has_subdomain() const;
    
    subdomain import_subdomain(string filename);
    subdomain export_subdomain(string filename);
    vector<int> build_permutation_vector();

    bool is_boundary(const vertex &) const;
    bool is_boundary(const subdomain &, const vertex &) const;
    point to_real_coord(voxel_point) const;
    voxel_point to_voxel_coord(point) const;



    public:

    grid(const grid_conds &);

    void fill_interaction(interaction &, const vector<int> & input, const vector<int> & output) const;
    vector<interaction> get_interaction_between(int i, int j) const;
    vector<interaction> get_self_interaction(int i) const;
    interaction get_inner_boundary_transpose_interaction() const;

    void import_grid(string filename) const;
    void export_grid(string filename);

    int get_dim() const {return dim;}
    const pair<point, point> & get_bbox() const {return bbox;}
    const pair<voxel_point, voxel_point> & get_vbbox() const {return N_bbox;}
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

struct interaction
{
    vector<vector<int>> interaction_db;
    int input_size;
    int output_size;
};






#endif
