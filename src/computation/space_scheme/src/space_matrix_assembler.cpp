#include "space_matrix_assembler.h"


space_matrix_assembler::space_matrix_assembler(
    const grid_problem & gp_,
    const space_scheme_config & ssc_
) :
gp(gp_), ssc(ssc_)
{
    if(!ssc.try_to_import)
    {
        build_subdomain_systems();
        build_base_boundary_system();
    }
}



void space_matrix_assembler::build_subdomain_systems()
{
    const auto & subdomains_interaction = gp.get_subdomain_interactions();

    subdomain_systems.resize(subdomains_interaction.size());

    for(int i=0; i<subdomain_systems.size(); i++)
    build_subdomain_system(i, subdomains_interaction[i]);
}



void space_matrix_assembler::build_subdomain_system(int i, const subdomain_interaction & si)
{
    subdomain_system & ss = subdomain_systems[i];

    build_spmat_from_interaction(si.ii, ss.ii);
    build_spmat_from_interaction(si.ib, ss.ib);
    build_spmat_from_interaction(si.bi, ss.bi);
    build_spmat_from_interaction(si.bb, ss.bb);
    build_spmat_from_interaction(si.global_gb, ss.gb_subd);
}



void space_matrix_assembler::build_spmat_from_interaction(const interaction & interact, spmat & target)
{
    target = spmat(interact.output_size, interact.input_size);

    vector<triplet> trp;

    for(const vertex_interaction & vi : interact.vertexes_io)
    {
        if(vi.input_indexes.empty()) continue;

        vector<double> coefficients = (*ssc.scheme)(vi);

        for(int i=0; i<vi.input_indexes.size(); i++)
        trp.emplace_back(vi.output_index, vi.input_indexes[i], coefficients[i]);
    }

    target.setFromTriplets(trp.begin(), trp.end());
}



void space_matrix_assembler::build_base_boundary_system()
{
    const interaction & interact = gp.get_base_boundary_interaction();

    base_boundary_system.clear();

    for(const vertex_interaction & vi : interact.vertexes_io)
    {
        if(vi.input_indexes.empty()) continue;

        vector<double> coefficients = (*ssc.scheme)(vi);

        for(int i=0; i<vi.input_indexes.size(); i++)
        base_boundary_system.emplace_back(vi.output_index, vi.input_indexes[i], coefficients[i]);
    }
}