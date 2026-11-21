#include "grid.h"

void grid::build_boundaries_interaction()
{
    boundaries_intr = setup_boundaries_interaction();
    auto bvio = boundaries_intr.vertexes_io;
    fill_bvio(bvio);
}

vector<int> grid::bvio_csum() const
{
    auto & subdomains = g->get_subdomains();

    vector<int> subdb_csum = {0};
    for(auto s : subdomains)
    subdb_csum.push_back(s.boundary.size()+subdb_csum.back());

    return subdb_csum;
}


interaction grid::setup_boundaries_interaction() const
{
    vector<int> subdomain_boundaries_csum = bvio_csum();

    interaction boundary_io;
    boundary_io.input_size = subdomain_boundaries_csum.back();
    boundary_io.output_size = subdomain_boundaries_csum.back();

    auto & bvio = boundary_io.vertexes_io;
    bvio.resize(boundary_io.output_size);

    return boundary_io;
}

void grid::fill_bvio(vector<vertex_interaction> & bvio)
{
    auto & subdomains = g->get_subdomains();
    auto & vertexes = g->get_vertexes();

    vector<int> subdb_csum = bvio_csum();

    for(auto s : subdomains)
    for(int i=0;i<s.boundary.size();i++)
    {
        const vertex & v = vertexes[s.boundary[i]];
        auto & bvi = bvio[i+subdb_csum[s.index]];

        bvi.output_index = i+subdb_csum[s.index];
        bvi.output_coord = v.voxel_coord;

        auto gin = global_indexes_of_neighbors_for_bi(v.index);
        auto gcn = global_coords_of_neighbors_for_bi(gin);
        auto grcn = g->to_real_coord(gcn);
        auto lin = local_indexes_of_neighbors_for_bi(gin, subdb_csum);        

        bvi.input_indexes.insert(bvi.input_indexes.end(), lin.begin(), lin.end());
        bvi.input_coords.insert(bvi.input_coords.end(), gcn.begin(), gcn.end());
        bvi.input_real_coords.insert(bvi.input_real_coords.end(), grcn.begin(), grcn.end());
    }
}

vector<int> grid::global_indexes_of_neighbors_for_bi(int i) const
{
    auto & vertexes = g->get_vertexes();

    const vertex & v = vertexes[i];
    const vector<int> & neighbors = v.neighbors_indexes;
    vector<int> result;

    for(int ni : neighbors)
    if(vertexes[ni].subdomain>v.subdomain)
    result.push_back(ni);
    
    return result;
}

vector<voxel_point> grid::global_coords_of_neighbors_for_bi(const vector<int> & gin) const
{
    auto & vertexes = g->get_vertexes();

    vector<voxel_point> gcn(gin.size());

    for(int i=0;i<gin.size();i++)
    gcn[i]=vertexes[gin[i]].voxel_coord;

    return gcn;
}

vector<int> grid::local_indexes_of_neighbors_for_bi(const vector<int> & gin, const vector<int> & csum) const
{
    auto & subdomains = g->get_subdomains();
    auto & vertexes = g->get_vertexes();

    vector<int> lin(gin.size());

    for(int i=0;i<gin.size();i++)
    {
        const vertex & v = vertexes[gin[i]];
        const subdomain & s = subdomains[v.subdomain];
        lin[i] = g->find_local_index(v.index, s.boundary) + csum[i];
    }

    return lin;
}