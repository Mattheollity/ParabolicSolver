#ifndef GRID_
#define GRID_

#include "common.h"

#include "problem_conditions.h"
#include "file_io.h"

#include "vertex.h"
#include "dim_iterator.h"
#include "interactions.h"

struct subdomain_interaction;

class grid;
struct grid_conds;
struct subdomain;



class grid
{
    grid_conds gc;

    vector<int> csum;
    vector<double> h;
    
    vector<vertex> vertexes;
    vector<int> inner;
    vector<int> boundary;

    vector<subdomain> subdomains;

    vector<subdomain_interaction> subd_intr;
    interaction boundaries_intr;

    vector<int> permutation;

    case_config cc;
    string GRID_DIR;



    void setup_grid_conds();



    void create_grid();
    void fill_domain();
    void find_neighbors_and_axis(vertex &);
    int global_index(const voxel_point &);

    void setup_subdomains();
    void create_subdomain(int i, pair<point, point> bbox_);
    void mark_subdomains_vertexes(subdomain &);
    void setup_subdomain_N_bbox(subdomain &, pair<point, point> &);
    void find_inner_boundary_at_subdomains();
    void assert_all_vertexes_has_subdomain() const;
    int amount_neighbors_at_subdomain(const vertex &) const;
    


    void build_permutation_vector();

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



    void build_self_subdomain_interactions();
    void build_subdomain_interaction(int i);
    vector<interaction_request> build_subdomain_interaction_requests(int i);
    void build_part_of_subdomain_interaction(interaction &, const interaction_request &);
    vertex_interaction create_vertex_interactions(const vertex_interaction_request &) const;



    void build_boundaries_interaction();

    vector<int> bvio_csum() const;
    interaction setup_boundaries_interaction() const;

    void fill_bvio(vector<vertex_interaction> &);
    vector<int> local_indexes_of_neighbors_for_bi(const vector<int> & gni, const vector<int> & csum) const;
    vector<int> global_indexes_of_neighbors_for_bi(int i) const;
    vector<voxel_point> global_coords_of_neighbors_for_bi(const vector<int> &) const;

    int find_local_index(int i, const vector<int> &) const;

    
    interaction create_global_interaction();
    



    public:

    grid(case_config cc_);

    vector<subdomain_interaction> get_subdomain_interactions() const;
    interaction get_boundaries_interaction() const;

    const grid_conds get_gc() const {return gc;}
    const vector<double> & get_h() const {return h;}
    const vector<point> & get_boundary_points() const;
    const vector<point> & get_inner_points() const;
    const vector<int> & get_permutation_vector() const {return permutation;}
    int get_amount_subdomains() const {return subdomains.size();}
    
};



struct subdomain
{
    int index;
    pair<voxel_point, voxel_point> N_bbox;
    vector<int> inner;
    vector<int> boundary;

    set<int> neighboring_subdomains;
};






#endif