#ifndef GEOMETRY_
#define GEOMETRY_

#include "common.h"
#include "vertex.h"
#include "dim_iterator.h"
#include "problem_conditions.h"
#include "file_io.h"
#include <vector>

struct subdomain;
class geometry
{
    case_config cc;
    string GRID_DIR;

    grid_conds gc;
    vector<int> csum;
    vector<double> h;

    vector<vertex> vertexes;
    vector<int> inner;
    vector<int> boundary;
    vector<int> permutation;

    vector<subdomain> subdomains;



    void setup_grid_conds();

    void create_grid();
    void fill_domain();
    void find_neighbors_and_axis(vertex &);
    int global_index(const voxel_point &);

    void setup_subdomains();
    void create_subdomain(int i, pair<point, point> bbox_);
    void mark_subdomains_vertexes(subdomain &);
    void setup_subdomain_N_bbox(subdomain &, pair<point, point> &);
    void find_inner_and_boundary_at_subdomains();
    void assert_all_vertexes_has_subdomain() const;
    int amount_neighbors_at_subdomain(const vertex &) const;
    
    void build_permutation_vector();



    bool import_grid();
    bool import_domain();
    bool import_subdomains();

    void export_grid() const;
    void export_domain() const;
    void export_subdomains() const;



    public:

    geometry(case_config cc_);

    int find_local_index(int i, const vector<int> &) const;
    point to_real_coord(voxel_point) const;
    vector<point> to_real_coord(vector<voxel_point>) const;
    voxel_point to_voxel_coord(point) const;
    point compute_abs_h(const point & a, const point & b) const;

    const vector<vertex> & get_vertexes() const {return vertexes;}
    const vector<int> & get_inner() const {return inner;}
    const vector<int> & get_boundary() const {return boundary;}
    const vector<int> & get_permutation() const {return permutation;}

    const vector<subdomain> & get_subdomains() const {return subdomains;}
};

struct subdomain
{
    int index;
    pair<voxel_point, voxel_point> N_bbox;
    vector<int> inner;
    vector<int> boundary;
    vector<int> global_boundary;

    set<int> neighboring_subdomains;
};

#endif