#include "parabolic_solver.h"
#include "cartesian_uniform.h"

#include <mpi.h>
#include <omp.h>

#include <filesystem>
#include <charconv>
#include <algorithm>
#include <stdexcept>

namespace fs = std::filesystem;



parabolic_solver::parabolic_solver(string case_dir_) :
case_dir(fs::absolute(case_dir_).lexically_normal().string()),
cfg(load_case(case_dir))
{
    int initialized = 0;

    MPI_Initialized(&initialized);

    if(initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    }

    if(cfg.cc.OpenMP_procs<=0)
    throw simple_error("Invalid OpenMP threads amount");

    if(cfg.wc.write_interval<=0)
    throw simple_error("Invalid write interval");

    if(cfg.tsc.end_time<cfg.tsc.start_time)
    throw simple_error("Invalid time interval");

    omp_set_num_threads(cfg.cc.OpenMP_procs);



    fs::path write_dir(cfg.wc.WRITE_DIR);

    if(write_dir.is_relative())
    write_dir = fs::path(case_dir)/write_dir;

    cfg.wc.WRITE_DIR = write_dir.lexically_normal().string();
}



bool parabolic_solver::check_geometry() const
{
    string path = case_dir+"/geometry";

    if(!fs::exists(path)) return false;

    if(!fs::is_directory(path))
    throw simple_error("Invalid geometry directory");

    if(!fs::is_regular_file(path+"/vertexes") ||
       !fs::is_regular_file(path+"/permutation") ||
       !fs::is_regular_file(path+"/subdomains"))
    throw simple_error("Incomplete geometry");



    auto file = load_file(path+"/subdomains");

    int amount = load_int("amount", file);

    if(amount<=0 || amount!=cfg.gc.subdomains.size())
    throw simple_error("Invalid amount of saved subdomains");

    for(int i=0; i<amount; i++)
    if(!fs::is_regular_file(path+"/"+to_string(i)))
    throw simple_error("Incomplete subdomain geometry");

    return true;
}



bool parabolic_solver::check_space_matrixes() const
{
    string path = case_dir+"/space_matrix";

    if(!fs::exists(path)) return false;

    if(!fs::is_directory(path))
    throw simple_error("Invalid spatial matrix directory");

    if(!fs::is_regular_file(path+"/base_boundary_system"))
    throw simple_error("Missing base boundary system");



    vector<string> names =
    {
        "_ii",
        "_ib",
        "_bi",
        "_bb",
        "_gb"
    };

    for(int i=0; i<cfg.gc.subdomains.size(); i++)
    for(const string & name : names)
    if(!fs::is_regular_file(path+"/"+to_string(i)+name))
    throw simple_error("Incomplete spatial approximation");

    return true;
}



int parabolic_solver::find_last_time_step() const
{
    if(!fs::exists(cfg.wc.WRITE_DIR))
    return -1;

    int last_time_step = -1;

    for(const auto & entry :
        fs::directory_iterator(cfg.wc.WRITE_DIR))
    {
        if(!entry.is_directory()) continue;

        string name = entry.path().filename().string();

        int step;

        auto result = std::from_chars(
            name.data(),
            name.data()+name.size(),
            step
        );

        if(result.ec!=std::errc() ||
           result.ptr!=name.data()+name.size())
        continue;

        if(step<0 || to_string(step)!=name)
        continue;

        if(!fs::is_regular_file(entry.path()/"solution"))
        continue;

        last_time_step = std::max(last_time_step, step);
    }

    return last_time_step;
}



void parabolic_solver::setup_solver()
{
    if(gp)
    throw simple_error("Solver is already initialized");

    int import_plan[3] = {0, 0, -1};

    if(rank==0)
    {
        import_plan[0] =
            cfg.gc.try_to_import &&
            check_geometry();

        import_plan[1] =
            import_plan[0] &&
            cfg.ssc.try_to_import &&
            check_space_matrixes();

        import_plan[2] = find_last_time_step();



        if(import_plan[2]>=0)
        {
            if(!cfg.tsc.try_to_import)
            throw simple_error("Previous calculation exists");

            if(!import_plan[0] || !import_plan[1])
            throw simple_error(
                "Cannot restore solution without saved approximation"
            );
        }

        if(cfg.gc.pass_creationg && !import_plan[0])
        throw simple_error("Geometry creation is disabled");
    }



    if(process_count>1)
    MPI_Bcast(
        import_plan,
        3,
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );



    setup_grid_problem(import_plan[0]);

    setup_space_scheme(import_plan[1]);

    setup_time_scheme();

    restore_solution(import_plan[2]);
}



void parabolic_solver::setup_grid_problem(
    bool import_geometry
)
{
    grid_config gc = cfg.gc;

    gc.try_to_import = import_geometry || rank!=0;

    unique_ptr<geometry> g =
        std::make_unique<cartesian_uniform>(gc);



    if(rank==0)
    {
        geometry_io gio;

        string path = case_dir+"/geometry";

        if(import_geometry)
        gio.import_geometry(*g, path);

        else
        gio.export_geometry(*g, path);
    }



    gp = std::make_unique<grid_problem>(
        std::move(g),
        std::move(cfg.ibc)
    );
}



void parabolic_solver::setup_space_scheme(
    bool import_matrixes
)
{
    cfg.ssc.try_to_import = import_matrixes || rank!=0;

    sma = std::make_unique<space_matrix_assembler>(
        *gp,
        cfg.ssc
    );

    if(rank!=0) return;



    space_matrix_assembler_io smaio;

    string path = case_dir+"/space_matrix";

    if(import_matrixes)
    smaio.import_matrixes(*sma, path);

    else
    smaio.export_matrixes(*sma, path);
}



void parabolic_solver::setup_time_scheme()
{
    tsa = std::make_unique<time_scheme_assembler>(
        *gp,
        cfg.tsc
    );

    tsa->setup_solver(*sma);

    tsio = std::make_unique<time_scheme_io>(
        cfg.wc
    );
}



void parabolic_solver::restore_solution(
    int last_time_step
)
{
    if(last_time_step<0) return;

    string path = cfg.wc.WRITE_DIR+"/"+
        to_string(last_time_step);

    tsio->import_solution(*tsa, path);
}



void parabolic_solver::import_solution(string path)
{
    if(!tsa)
    throw simple_error("Solver is not initialized");

    fs::path solution_path(path);

    if(solution_path.is_relative())
    solution_path = fs::path(case_dir)/solution_path;

    tsio->import_solution(
        *tsa,
        solution_path.string()
    );
}



void parabolic_solver::export_solution(string path) const
{
    if(!tsa)
    throw simple_error("Solver is not initialized");

    if(rank!=0) return;

    fs::path solution_path(path);

    if(solution_path.is_relative())
    solution_path = fs::path(case_dir)/solution_path;

    tsio->export_solution(
        *tsa,
        solution_path.string()
    );
}



void parabolic_solver::restart_solver()
{
    if(!tsa)
    throw simple_error("Solver is not initialized");

    if(rank==0 && find_last_time_step()>=0)
    throw simple_error(
        "Remove previous results before restarting"
    );

    tsa->restart_solver();
}



void parabolic_solver::time_step()
{
    if(!tsa)
    throw simple_error("Solver is not initialized");

    if(tsa->get_current_time()>=cfg.tsc.end_time)
    return;

    double previous_time = tsa->get_current_time();

    tsa->time_step();

    if(tsa->get_current_time()<=previous_time)
    throw simple_error("Time integration has stopped");
}



void parabolic_solver::solve()
{
    if(!tsa)
    throw simple_error("Solver is not initialized");

    if(rank==0)
    {
        int last_time_step = find_last_time_step();

        if(last_time_step>tsa->get_current_time_step())
        throw simple_error("Newer saved solution exists");

        if(last_time_step<tsa->get_current_time_step())
        tsio->write_solution(*tsa);
    }



    while(tsa->get_current_time()<cfg.tsc.end_time)
    {
        time_step();

        if(rank!=0) continue;

        int step = tsa->get_current_time_step();

        double time = tsa->get_current_time();

        if(step%cfg.wc.write_interval==0 ||
           time>=cfg.tsc.end_time)
        {
            tsio->write_solution(*tsa);

            cout<<"Time step: "<<step
                <<", time: "<<time<<endl;
        }
    }
}