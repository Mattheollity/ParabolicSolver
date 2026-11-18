#include "grid.h"
#include "file_io.h"
#include "problem_conditions.h"
#include <cstdio>
#include <omp.h>
#include <vector>

grid::grid(case_config cc_) :
    GRID_DIR(cc.GRID_DIR)
{
    g = new geometry(cc);
    build_self_subdomain_interactions();
    build_boundaries_interaction();
}

void grid::build_self_subdomain_interactions()
{
    auto & subdomains = g->get_subdomains();

    subd_intr.resize(subdomains.size());

    for(int i=0;i<subd_intr.size();i++)
    build_subdomain_interaction(i);
}

void grid::build_subdomain_interaction(int i)
{
    auto si = subd_intr[i];
    auto requests = build_subdomain_interaction_requests(i);

    build_part_of_subdomain_interaction(si.ii, requests[0]);
    build_part_of_subdomain_interaction(si.ib, requests[1]);
    build_part_of_subdomain_interaction(si.bb, requests[2]);
    build_part_of_subdomain_interaction(si.global_gb, requests[3]);
}

vector<interaction_request> grid::build_subdomain_interaction_requests(int i)
{
    auto & subdomains = g->get_subdomains();
    auto & boundary = g->get_boundary();

    const vector<int> & inner_ = subdomains[i].inner;
    const vector<int> & boundary_ = subdomains[i].boundary;
    const vector<int> & global_boundary = boundary;

    interaction_request ii = {
        i,
        i,
        inner_,
        inner_
    };

    interaction_request ib = {
        i,
        i,
        inner_,
        boundary_
    };

    interaction_request bb = {
        i,
        i,
        boundary_,
        boundary_
    };

    interaction_request gbb = {
        i,
        i,
        global_boundary,
        boundary_
    };

    return {ii, ib, bb, gbb};
}

void grid::build_part_of_subdomain_interaction(interaction & io, const interaction_request & ior)
{
    io.input_size = ior.input.size();
    io.output_size = ior.output.size();
    io.vertexes_io = vector<vertex_interaction>(io.output_size);

    #pragma omp parallel shared(io, ior)
    #pragma omp for
    for(int i=0; i<io.output_size; i++)
    {
        vertex_interaction_request vir = {
            i,
            ior.output[i],
            ior.input_subdomain,
            ior.input
        };
        io.vertexes_io[i] = create_vertex_interactions(vir);
    }
}

vertex_interaction grid::create_vertex_interactions(const vertex_interaction_request & vir) const
{
    auto & vertexes = g->get_vertexes();

    const vertex & v = vertexes[vir.vertex_index];
    vertex_interaction vio;
    vio.output_index = vir.output_index;
    vio.output_coord = v.voxel_coord;

    vector<int> input_vertex_indexes;
    for(auto nv : v.neighbors_indexes)
    if(vertexes[nv].subdomain==vir.input_subdomain)
    input_vertex_indexes.push_back(nv);

    for(int i : input_vertex_indexes)
    vio.input_indexes.push_back(g->find_local_index(i, vir.input));
    
    for(int i : input_vertex_indexes)
    vio.input_coords.push_back(vertexes[i].voxel_coord);

    for(int i : input_vertex_indexes)
    vio.input_real_coords.push_back(vertexes[i].real_coord);

    return vio;
}



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


grid::~grid()
{
    delete g;
}