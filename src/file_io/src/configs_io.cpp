#include "configs_io.h"

case_config load_case(string CASE_DIR)
{
    case_config case_obj
    (
        load_space_scheme(CASE_DIR),
        load_time_scheme(CASE_DIR),
        load_zero_boundary(CASE_DIR),
        load_grid(CASE_DIR),
        load_computation(CASE_DIR),
        load_writter(CASE_DIR)
    );
}

space_scheme_config load_space_scheme(string CASE_DIR)
{
    string scheme_path = CASE_DIR+"/space_scheme";
    auto file = load_file(scheme_path);

    space_scheme_config ssc;
    string scheme_type = load_string("scheme_type", file);

    if(scheme_type=="anisotrop_laplace")
    ssc.scheme = std::make_unique<anisotrop_laplace>(scheme_path);
    else throw simple_error(string("undefined scheme_type ")+scheme_type);

    return ssc;
}


time_scheme_config load_time_scheme(string CASE_DIR)
{
    string time_scheme_path = CASE_DIR+"/time_scheme";
    auto file = load_file(time_scheme_path);

    double start_time = load_field<double>("start_time", file);
    double end_time = load_field<double>("end_time", file);
    double T = load_field<double>("T", file);
    double theta = load_field<double>("theta", file);
    bool try_to_import = load_bool("try_to_import", file);

    return time_scheme_config(start_time, end_time, T, theta, try_to_import);
}

grid_config load_grid(string CASE_DIR)
{
    string grid_path = CASE_DIR+"grid";
    auto file= load_file(grid_path);

    grid_config gc = {
        0,
        {{},{}},
        {},
        {},
        false
    };

    gc.try_to_import = load_bool("try_to_import", file);
    gc.dim = load_int("dim", file);
    gc.bbox = load_pair_of_vectors<double>("bbox", gc.dim, file);
    gc.N = load_vector<int>("N", gc.dim, file);

    expected_section("subdomains", file);
    string s = "subdomain_";
    gc.subdomains.resize(load_int("amount_of_subdomains", file));
    for(int i=0;i<gc.subdomains.size();i++)
    gc.subdomains[i]=load_pair_of_vectors<double>(s+to_string(i), gc.dim, file);

    return gc;
}

zero_boundary_config load_zero_boundary(string CASE_DIR)
{
    string zero_boundary_path = CASE_DIR+"zero_boundary";
    auto file(load_file(zero_boundary_path));
    string y0 = load_string("y0", file);
    string boundary_type = load_string("boundary_type", file);
    string boundary = load_string("boundary", file);

    zero_boundary_config zbc;
    
    if(y0=="trip_sin") zbc.y0=dimension_sine;
    else throw simple_error("undefined start_time function");

    if(boundary_type=="dirichlet")
    {
        if(boundary=="zero") zbc.dirichlet=zero_value;
        else throw simple_error("undefined dirichlet boundary function");
    }
    else throw simple_error("undefined boundary type");

    return zbc;
}

computation_config load_computation(string CASE_DIR)
{
    string computation_path = CASE_DIR+"computation";
    auto file = load_file(computation_path);

    int total_memory = load_int("total_memory(MB)", file);
    int total_time = load_int("total_time(s)", file);
    int OpenMP_procs = load_int("OpenMP_procs", file);
    int MPI_procs = load_int("MPI_procs", file);

    return computation_config
    (
        total_memory, total_time, 
        OpenMP_procs, MPI_procs
    );
}

writter_config load_writter(string CASE_DIR)
{
    string writter_path = CASE_DIR+"writter";
    auto file = load_file(writter_path);

    int write_interval = load_int("write_interval(steps)",file);
    string WRITE_DIR = load_string("write_dir", file);
    bool paraview_output = load_bool("paraview_output", file);

    return writter_config
    (
        write_interval,
        WRITE_DIR,
        paraview_output
    );
}