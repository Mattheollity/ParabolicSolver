#include "space_matrix_assembler_io.h"

void space_matrix_assembler_io::import_matrixes(space_matrix_assembler & sma, string path)
{
    auto & gc = sma.ssc.g->gc;
    auto & subdomain_systems = sma.subdomain_systems;

    subdomain_systems.resize(gc.subdomains.size());

    #pragma omp parallel for shared(subdomain_systems)
    for(int i=0;i<subdomain_systems.size();i++)
    {
        import_spmat(to_string(i)+"_ii", subdomain_systems[i].ii);
        import_spmat(to_string(i)+"_ib", subdomain_systems[i].ib);
        import_spmat(to_string(i)+"_bb", subdomain_systems[i].bb);
        import_spmat(to_string(i)+"_gb", subdomain_systems[i].gb_subd);
    }

    import_spmat("base_boundary", sma.base_boundary_system);
}

void space_matrix_assembler_io::export_matrixes(const space_matrix_assembler & sma, string path) const
{
    string sma_export_path = path+"/space_matrix";
    fs::create_directory(sma_export_path);
    auto & subdomain_systems = sma.subdomain_systems;

    #pragma omp parallel for shared(subdomain_systems)
    for(int i=0;i<subdomain_systems.size();i++)
    {
        export_spmat(sma_export_path+to_string(i)+"_ii", subdomain_systems[i].ii);
        export_spmat(sma_export_path+to_string(i)+"_ib", subdomain_systems[i].ib);
        export_spmat(sma_export_path+to_string(i)+"_bb", subdomain_systems[i].bb);
        export_spmat(sma_export_path+to_string(i)+"_gb", subdomain_systems[i].gb_subd);
    }

    export_spmat("base_boundary_system", sma.base_boundary_system);
}