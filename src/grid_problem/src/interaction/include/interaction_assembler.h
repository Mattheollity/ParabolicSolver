#ifndef GRID_
#define GRID_

#include "common.h"
#include "config_structures.h"

#include "geometry.h"
#include "interactions.h"

class grid
{
    unique_ptr<geometry> grid_geometry;
    
    vector<subdomain_interaction> subd_intr;
    interaction boundaries_intr;

    

    void build_self_subdomain_interactions();
    void build_subdomain_interaction(int i);
    vector<interaction_request> build_subdomain_interaction_requests(int i);
    void build_part_of_subdomain_interaction(interaction &, const interaction_request &);
    vertex_interaction create_vertex_interactions(const vertex_interaction_request &) const;



    void build_boundaries_interaction();

    interaction setup_boundaries_interaction() const;
    vector<int> bvio_csum() const;

    void fill_bvio(vector<vertex_interaction> &);
    vector<int> local_indexes_of_neighbors_for_bi(const vector<int> & gni, const vector<int> & csum) const;
    vector<int> global_indexes_of_neighbors_for_bi(int i) const;
    vector<voxel_point> global_coords_of_neighbors_for_bi(const vector<int> &) const;



    public:

    const grid_config & gc;

    grid(const grid_config &);

    vector<subdomain_interaction> get_subdomain_interactions() const;
    interaction get_boundaries_interaction() const;
    vector<vector<point>> get_global_boundaries_points() const;

    const vector<int> & get_permutation_vector() const {return grid_geometry->get_permutation();}
};






#endif