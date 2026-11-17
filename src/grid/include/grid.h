#ifndef GRID_
#define GRID_

#include "common.h"

#include "geometry.h"
#include "interactions.h"

struct subdomain_interaction;

class grid;
struct grid_conds;
struct subdomain;



class grid
{
    case_config cc;
    string GRID_DIR;

    geometry * g;
    
    vector<subdomain_interaction> subd_intr;
    interaction boundaries_intr;

    

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

    
    interaction create_global_interaction();
    



    public:

    grid(case_config cc_);
    ~grid();

    vector<subdomain_interaction> get_subdomain_interactions() const;
    interaction get_boundaries_interaction() const;

    const vector<int> & get_permutation_vector() const {return g->get_permutation();}
};






#endif