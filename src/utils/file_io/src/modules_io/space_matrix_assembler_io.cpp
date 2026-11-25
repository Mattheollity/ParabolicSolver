#include "space_matrix_assembler_io.h"


namespace
{
    void check_matrix_size(
        const spmat & m,
        const interaction & intr
    )
    {
        if(m.rows()!=intr.output_size ||
           m.cols()!=intr.input_size)
        throw simple_error("Invalid spatial matrix size");
    }
}



void space_matrix_assembler_io::import_matrixes(
    space_matrix_assembler & sma,
    string path
)
{
    const auto & interactions =
        sma.gp.get_subdomain_interactions();

    auto & systems = sma.subdomain_systems;

    systems.resize(interactions.size());

    for(int i=0; i<systems.size(); i++)
    {
        subdomain_system & ss = systems[i];

        string prefix = path+"/"+to_string(i);

        import_spmat(prefix+"_ii", ss.ii);
        import_spmat(prefix+"_ib", ss.ib);
        import_spmat(prefix+"_bi", ss.bi);
        import_spmat(prefix+"_bb", ss.bb);
        import_spmat(prefix+"_gb", ss.gb_subd);

        const subdomain_interaction & si = interactions[i];

        check_matrix_size(ss.ii, si.ii);
        check_matrix_size(ss.ib, si.ib);
        check_matrix_size(ss.bi, si.bi);
        check_matrix_size(ss.bb, si.bb);
        check_matrix_size(ss.gb_subd, si.global_gb);
    }

    import_spmat(
        path+"/base_boundary_system",
        sma.base_boundary_system
    );

    const interaction & bi =
        sma.gp.get_base_boundary_interaction();

    for(const triplet & t : sma.base_boundary_system)
    if(t.row()<0 || t.row()>=bi.output_size ||
       t.col()<0 || t.col()>=bi.input_size)
    throw simple_error("Invalid boundary triplet");
}



void space_matrix_assembler_io::export_matrixes(
    const space_matrix_assembler & sma,
    string path
) const
{
    fs::create_directories(path);

    const auto & systems = sma.subdomain_systems;

    for(int i=0; i<systems.size(); i++)
    {
        const subdomain_system & ss = systems[i];

        string prefix = path+"/"+to_string(i);

        export_spmat(prefix+"_ii", ss.ii);
        export_spmat(prefix+"_ib", ss.ib);
        export_spmat(prefix+"_bi", ss.bi);
        export_spmat(prefix+"_bb", ss.bb);
        export_spmat(prefix+"_gb", ss.gb_subd);
    }

    export_spmat(
        path+"/base_boundary_system",
        sma.base_boundary_system
    );
}