#include "space_matrix_assembler.h"

space_matrix_assembler::space_matrix_assembler(const space_scheme_config & ssc_) :
    ssc(ssc_)
{
    if(!ssc.try_to_import)
    {
        build_subdomain_systems();
        build_base_boundary_system();
    }
}

void space_matrix_assembler::build_subdomain_systems()
{
    const auto & subdomains_interaction = ssc.g->get_subdomain_interactions();
    subdomain_systems.resize(subdomains_interaction.size());
    
    for(int i=0;i<subdomain_systems.size();i++)
    build_subdomain_system(i, subdomains_interaction[i]);
}

void space_matrix_assembler::build_subdomain_system(int i, const subdomain_interaction & si)
{
    subdomain_system & ss = subdomain_systems[i];
    build_spmat_from_interaction(si.ii, ss.ii);
    build_spmat_from_interaction(si.ib, ss.ib);
    build_spmat_from_interaction(si.bb, ss.bb);
    build_spmat_from_interaction(si.global_gb, ss.gb_subd);
}

void space_matrix_assembler::build_spmat_from_interaction(const interaction & interact, spmat & target)
{
    target = spmat(interact.output_size, interact.input_size);
    vector<triplet> trp;

    for(int i=0;i<interact.output_size;i++)
    for(int j=0;j<interact.vertexes_io[i].input_indexes.size();j++)
    trp.push_back
    (triplet(
        interact.vertexes_io[i].output_index,
        interact.vertexes_io[i].input_indexes[j],
        ssc.scheme->operator()(space_scheme_request(i, j, interact))
    ));

    target.setFromTriplets(trp.begin(), trp.end());
}

void space_matrix_assembler::build_base_boundary_system()
{
    const interaction & interact = ssc.g->get_boundaries_interaction();

    for(int i=0;i<interact.output_size;i++)
    for(int j=0;j<interact.vertexes_io[i].input_indexes.size();j++)
    base_boundary_system.push_back
    (triplet(
        interact.vertexes_io[i].output_index,
        interact.vertexes_io[i].input_indexes[j],
        ssc.scheme->operator()(space_scheme_request(i, j, interact))
    ));
}