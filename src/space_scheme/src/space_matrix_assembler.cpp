#include "space_matrix_assembler.h"
#include "file_io.h"
#include "grid.h"
#include "problem_conditions.h"
#include "space_scheme.h"
#include <string>

space_matrix_assembler::space_matrix_assembler(case_config cc_) :
    cc(cc_)
{
    if(ability_to_import()) import_matrixes();
    else create_from_grid();
}

void space_matrix_assembler::create_from_grid()
{
    g = new grid(cc);
    set_scheme();

    build_permutation();
    build_subdomain_systems();
    build_base_boundary_system();

    export_matrixes();

    delete g;
}

void space_matrix_assembler::set_scheme()
{
    scheme = load_space_scheme(cc.SPACE_DIR);
}

void space_matrix_assembler::build_permutation()
{
    auto permutation_vector = g->get_permutation_vector();
    vector<triplet> trp (permutation_vector.size());

    #pragma omp parallel for shared(permutation_vector, trp)
    for(int i=0; i<trp.size(); i++)
    trp[i] = triplet(i, permutation_vector[i], 1);

    permutation.setFromTriplets(trp.begin(), trp.end());
}

void space_matrix_assembler::build_subdomain_systems()
{
    const auto & subdomains_interaction = g->get_subdomain_interactions();
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
        scheme(interact.vertexes_io[i].output_coord, interact.vertexes_io[i].input_coords[j])
    ));

    target.setFromTriplets(trp.begin(), trp.end());
}

void space_matrix_assembler::build_base_boundary_system()
{
    const interaction & interact = g->get_boundaries_interaction();

    for(int i=0;i<interact.output_size;i++)
    for(int j=0;j<interact.vertexes_io[i].input_indexes.size();j++)
    base_boundary_system.push_back
    (triplet(
        interact.vertexes_io[i].output_index,
        interact.vertexes_io[i].input_indexes[j],
        scheme(interact.vertexes_io[i].output_coord, interact.vertexes_io[i].input_coords[j])
    ));
}

bool space_matrix_assembler::ability_to_import()
{
    if
    (
        fs::exists(cc.SPACE_DIR+"/used_grid_conds") &&
        fs::exists(cc.GRID_DIR+"/exp_grid_conds") &&
        fs::exists(cc.GRID_DIR+"grid_conds")
    )
    if
    (
        file_equality(cc.GRID_DIR+"/grid_conds", cc.GRID_DIR+"/exp_grid_conds") &&
        file_equality(cc.SPACE_DIR+"/used_grid_conds", cc.GRID_DIR+"/grid_conds")
    )
    return true;

    return false;
}

void space_matrix_assembler::import_matrixes()
{
    grid_conds gc = load_grid_conds(cc.GRID_DIR+"/used_grid_conds");

    subdomain_systems.resize(gc.subdomains.size());

    #pragma omp parallel for shared(subdomain_systems)
    for(int i=0;i<subdomain_systems.size();i++)
    {
        import_spmat(to_string(i)+"_ii", subdomain_systems[i].ii);
        import_spmat(to_string(i)+"_ib", subdomain_systems[i].ib);
        import_spmat(to_string(i)+"_bb", subdomain_systems[i].bb);
        import_spmat(to_string(i)+"_gb", subdomain_systems[i].gb_subd);
    }

    import_spmat("base_boundary", base_boundary_system);
}

void space_matrix_assembler::export_matrixes()
{
    fs::copy_file(cc.GRID_DIR+"/grid_conds",cc.SPACE_DIR+"/used_grid_conds");

    export_spmat(cc.SPACE_DIR+"/permutation", permutation);

    #pragma omp parallel for shared(subdomain_systems)
    for(int i=0;i<subdomain_systems.size();i++)
    {
        export_spmat(to_string(i)+"_ii", subdomain_systems[i].ii);
        export_spmat(to_string(i)+"_ib", subdomain_systems[i].ib);
        export_spmat(to_string(i)+"_bb", subdomain_systems[i].bb);
        export_spmat(to_string(i)+"_gb", subdomain_systems[i].gb_subd);
    }

    export_spmat("base_boundary_system", base_boundary_system);
}

