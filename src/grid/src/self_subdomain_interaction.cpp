#include "grid.h"

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