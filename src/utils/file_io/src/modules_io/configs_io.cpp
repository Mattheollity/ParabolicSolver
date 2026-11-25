#include "configs_io.h"

#include <iomanip>
#include <limits>


namespace
{
    std::ofstream open_config(string CASE_DIR, string NAME)
    {
        fs::create_directories(CASE_DIR);

        string path = (fs::path(CASE_DIR)/NAME).string();

        std::ofstream file(path);

        if(!file)
        throw simple_error("Cannot open "+path);

        file<<std::setprecision(
            std::numeric_limits<double>::max_digits10
        );

        return file;
    }



    template<typename DATATYPE>
    void write_vector_values(
        const vector<DATATYPE> & v,
        std::ofstream & file
    )
    {
        for(const DATATYPE & value : v)
        file<<value<<" ";

        file<<"\n";
    }
}



case_config load_case(string CASE_DIR)
{
    return
    {
        load_space_scheme(CASE_DIR),
        load_time_scheme(CASE_DIR),
        load_initial_boundary(CASE_DIR),
        load_grid(CASE_DIR),
        load_computation(CASE_DIR),
        load_writter(CASE_DIR)
    };
}



space_scheme_config load_space_scheme(string CASE_DIR)
{
    auto file = load_file(
        (fs::path(CASE_DIR)/"space_scheme").string()
    );

    space_scheme_config ssc;

    string scheme_type = load_string("scheme_type", file);

    if(scheme_type=="anisotrop_laplace")
    {
        unique_ptr<anisotrop_laplace> scheme;

        build_anisotrop_laplace(scheme, file);

        ssc.scheme = std::move(scheme);
    }
    else
    throw simple_error("Undefined scheme_type "+scheme_type);

    file>>std::ws;

    if(!file.eof())
    ssc.try_to_import = load_bool("try_to_import", file);

    return ssc;
}



time_scheme_config load_time_scheme(string CASE_DIR)
{
    auto file = load_file(
        (fs::path(CASE_DIR)/"time_scheme").string()
    );

    time_scheme_config tsc;

    tsc.start_time = load_field<double>("start_time", file);
    tsc.end_time = load_field<double>("end_time", file);

    double T = load_field<double>("T", file);
    double theta = load_field<double>("theta", file);

    tsc.try_to_import = load_bool("try_to_import", file);

    tsc.scheme = std::make_unique<weighted_time_scheme>(
        theta, T
    );

    return tsc;
}



grid_config load_grid(string CASE_DIR)
{
    auto file = load_file(
        (fs::path(CASE_DIR)/"grid").string()
    );

    grid_config gc;

    gc.try_to_import = load_bool("try_to_import", file);
    gc.pass_creationg = load_bool("pass_creation", file);

    gc.dim = load_int("dim", file);

    gc.bbox = load_pair_of_vectors<double>(
        "bbox", gc.dim, file
    );

    gc.N = load_vector<int>("N", gc.dim, file);

    expected_section("subdomains", file);

    int amount = load_int("amount_of_subdomains", file);

    if(amount<0)
    throw simple_error("Invalid amount of subdomains");

    gc.subdomains.resize(amount);

    for(int i=0; i<amount; i++)
    gc.subdomains[i] = load_pair_of_vectors<double>(
        "subdomain_"+to_string(i),
        gc.dim,
        file
    );

    return gc;
}



initial_boundary_config load_initial_boundary(
    string CASE_DIR
)
{
    auto file = load_file(
        (fs::path(CASE_DIR)/"zero_boundary").string()
    );

    initial_boundary_config ibc;

    string y0 = load_string("y0", file);

    if(y0=="dimension_sine")
    ibc.y0 = std::make_unique<dimension_sine>();

    else if(y0=="zero_value")
    ibc.y0 = std::make_unique<zero_value>();

    else
    throw simple_error("Undefined initial condition "+y0);



    string boundary_type = load_string(
        "boundary_type", file
    );

    load_string("boundary", file);

    if(boundary_type=="const_dirichlet")
    {
        unique_ptr<const_dirichlet> cd;

        build_const_dirichlet(cd, file);

        ibc.dirichlet_condition = std::move(cd);
    }
    else
    throw simple_error("Undefined boundary type "+boundary_type);

    return ibc;
}



computation_config load_computation(string CASE_DIR)
{
    auto file = load_file(
        (fs::path(CASE_DIR)/"computation").string()
    );

    int total_memory = load_int("total_memory(MB)", file);
    int total_time = load_int("total_time(s)", file);

    int OpenMP_procs = load_int("OpenMP_procs", file);
    int MPI_procs = load_int("MPI_procs", file);

    return computation_config(
        total_memory,
        total_time,
        OpenMP_procs,
        MPI_procs
    );
}



writter_config load_writter(string CASE_DIR)
{
    auto file = load_file(
        (fs::path(CASE_DIR)/"writter").string()
    );

    int write_interval = load_int(
        "write_interval(steps)", file
    );

    string WRITE_DIR = load_string("write_dir", file);

    bool paraview_output = load_bool(
        "paraview_output", file
    );

    return writter_config(
        write_interval,
        WRITE_DIR,
        paraview_output
    );
}


void write_case(string CASE_DIR, const case_config & cfg)
{
    write_space_scheme(CASE_DIR, cfg.ssc);
    write_time_scheme(CASE_DIR, cfg.tsc);
    write_initial_boundary(CASE_DIR, cfg.ibc);
    write_grid(CASE_DIR, cfg.gc);
    write_computation(CASE_DIR, cfg.cc);
    write_writter(CASE_DIR, cfg.wc);
}



void write_space_scheme(
    string CASE_DIR,
    const space_scheme_config & ssc
)
{
    const auto * scheme =
        dynamic_cast<const anisotrop_laplace *>(
            ssc.scheme.get()
        );

    if(!scheme)
    throw simple_error("Undefined space scheme");

    auto file = open_config(CASE_DIR, "space_scheme");

    file<<"scheme_type anisotrop_laplace\n";

    write_anisotrop_laplace(*scheme, file);

    file<<"try_to_import "<<ssc.try_to_import<<"\n";
}



void write_time_scheme(
    string CASE_DIR,
    const time_scheme_config & tsc
)
{
    const auto * scheme =
        dynamic_cast<const weighted_time_scheme *>(
            tsc.scheme.get()
        );

    if(!scheme)
    throw simple_error("Undefined time scheme");

    auto file = open_config(CASE_DIR, "time_scheme");

    file<<"start_time "<<tsc.start_time<<"\n";
    file<<"end_time "<<tsc.end_time<<"\n";

    file<<"T "<<scheme->T<<"\n";
    file<<"theta "<<scheme->theta<<"\n";

    file<<"try_to_import "<<tsc.try_to_import<<"\n";
}



void write_grid(
    string CASE_DIR,
    const grid_config & gc
)
{
    auto file = open_config(CASE_DIR, "grid");

    file<<"try_to_import "<<gc.try_to_import<<"\n";
    file<<"pass_creation "<<gc.pass_creationg<<"\n";

    file<<"dim "<<gc.dim<<"\n";

    file<<"bbox\n";

    write_vector_values(gc.bbox.first, file);
    write_vector_values(gc.bbox.second, file);

    file<<"N ";

    write_vector_values(gc.N, file);



    file<<"subdomains\n";
    file<<"amount_of_subdomains "<<gc.subdomains.size()<<"\n";

    for(int i=0; i<gc.subdomains.size(); i++)
    {
        file<<"subdomain_"<<i<<"\n";

        write_vector_values(gc.subdomains[i].first, file);
        write_vector_values(gc.subdomains[i].second, file);
    }
}



void write_initial_boundary(
    string CASE_DIR,
    const initial_boundary_config & ibc
)
{
    string y0;

    if(dynamic_cast<const dimension_sine *>(ibc.y0.get()))
    y0 = "dimension_sine";

    else if(dynamic_cast<const zero_value *>(ibc.y0.get()))
    y0 = "zero_value";

    else
    throw simple_error("Undefined initial condition");



    const auto * cd =
        dynamic_cast<const const_dirichlet *>(
            ibc.dirichlet_condition.get()
        );

    if(!cd)
    throw simple_error("Undefined boundary condition");



    auto file = open_config(CASE_DIR, "zero_boundary");

    file<<"y0 "<<y0<<"\n";

    file<<"boundary_type const_dirichlet\n";
    file<<"boundary const_dirichlet\n";

    write_const_dirichlet(*cd, file);
}



void write_computation(
    string CASE_DIR,
    const computation_config & cc
)
{
    auto file = open_config(CASE_DIR, "computation");

    file<<"total_memory(MB) "<<cc.total_memory<<"\n";
    file<<"total_time(s) "<<cc.total_time<<"\n";

    file<<"OpenMP_procs "<<cc.OpenMP_procs<<"\n";
    file<<"MPI_procs "<<cc.MPI_procs<<"\n";
}



void write_writter(
    string CASE_DIR,
    const writter_config & wc
)
{
    auto file = open_config(CASE_DIR, "writter");

    file<<"write_interval(steps) "<<wc.write_interval<<"\n";
    file<<"write_dir "<<wc.WRITE_DIR<<"\n";

    file<<"paraview_output "<<wc.paraview_output<<"\n";
}