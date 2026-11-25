#include "geometry.h"

#include <algorithm>

void geometry::build_interactions()
{
    build_self_subdomain_interactions();
    build_boundaries_interaction();
}

void geometry::build_self_subdomain_interactions()
{
    subdomain_intr.resize(subdomains.size());

    for(int i=0; i<subdomains.size(); i++)
    build_subdomain_interaction(i);
}

void geometry::build_subdomain_interaction(int i)
{
    const subdomain & s = subdomains[i];
    subdomain_interaction & si = subdomain_intr[i];

    build_part_of_subdomain_interaction(si.ii, s.inner, s.inner);
    build_part_of_subdomain_interaction(si.ib, s.inner, s.boundary);
    build_part_of_subdomain_interaction(si.bb, s.boundary, s.boundary);
    build_part_of_subdomain_interaction(si.global_gb, s.global_boundary, s.boundary);
}

void geometry::build_part_of_subdomain_interaction(
    interaction & io,
    const vector<int> & input,
    const vector<int> & output
)
{
    io.input_size = input.size();
    io.output_size = output.size();

    io.vertexes_io.resize(io.output_size);

    for(int i=0; i<io.output_size; i++)
    io.vertexes_io[i] = create_vertex_interactions(i, output[i], input);
}

vertex_interaction geometry::create_vertex_interactions(
    int output_index,
    int vertex_index,
    const vector<int> & input
) const
{
    const vertex & v = vertexes[vertex_index];

    vertex_interaction vio;

    vio.output_index = output_index;
    vio.output_coord = v.real_coord;

    vector<int> indexes = v.neighbors_indexes;
    indexes.push_back(v.index);

    vio.interacting_coords.reserve(indexes.size());

    for(int i=0; i<indexes.size(); i++)
    {
        int ni = indexes[i];

        vio.interacting_coords.push_back(vertexes[ni].real_coord);

        auto it = std::lower_bound(input.begin(), input.end(), ni);

        if(it==input.end() || *it!=ni) continue;

        vio.input_indexes.push_back(it-input.begin());
        vio.input_coord_indexes.push_back(i);
    }

    return vio;
}

void geometry::build_boundaries_interaction()
{
    vector<int> csum = bvio_csum();

    boundaries_intr.input_size = csum.back();
    boundaries_intr.output_size = csum.back();

    boundaries_intr.vertexes_io = vector<vertex_interaction>(csum.back());

    fill_bvio(boundaries_intr.vertexes_io);
}

vector<int> geometry::bvio_csum() const
{
    vector<int> csum = {0};

    for(const subdomain & s : subdomains)
    csum.push_back(csum.back()+s.boundary.size());

    return csum;
}

void geometry::fill_bvio(vector<vertex_interaction> & bvio)
{
    vector<int> csum = bvio_csum();

    for(const subdomain & s : subdomains)
    for(int i=0; i<s.boundary.size(); i++)
    {
        const vertex & v = vertexes[s.boundary[i]];

        int index = i+csum[s.index];

        vertex_interaction & vi = bvio[index];

        vi.output_index = index;
        vi.output_coord = v.real_coord;

        vector<int> indexes = v.neighbors_indexes;
        indexes.push_back(v.index);

        vi.interacting_coords.reserve(indexes.size());

        for(int j=0; j<indexes.size(); j++)
        {
            int ni = indexes[j];

            const vertex & nv = vertexes[ni];

            vi.interacting_coords.push_back(nv.real_coord);

            if(nv.subdomain<=v.subdomain) continue;

            int local_index = find_local_index(
                ni,
                subdomains[nv.subdomain].boundary
            );

            vi.input_indexes.push_back(local_index+csum[nv.subdomain]);
            vi.input_coord_indexes.push_back(j);
        }
    }
}