#ifndef CARTESIAN_UNIFORM_
#define CARTESIAN_UNIFORM_

#include "geometry.h"
#include "config_structures.h"
#include "dim_iterator.h"


class cartesian_uniform : public geometry
{
    grid_config gc;
    point h;
    vector<int> csum;

    void fill_domain();
    void find_neighbors(vertex &, const voxel_point &);
    int global_index(const voxel_point &) const;

    void setup_subdomains();
    void create_subdomain(int, const pair<point, point> &);
    pair<voxel_point, voxel_point> setup_subdomain_N_bbox(const pair<point, point> &) const;
    void mark_subdomains_vertexes(subdomain &, const pair<voxel_point, voxel_point> &);
    void find_inner_and_boundary_at_subdomains();
    int amount_neighbors_at_subdomain(const vertex &) const;
    void assert_all_vertexes_has_subdomain() const;



    public:

    cartesian_uniform(const grid_config &);

    void setup_geometry() override;
};

#endif