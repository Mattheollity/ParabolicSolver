#include "time_scheme_assembler.h"

#include <mpi.h>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cmath>


time_scheme_assembler::time_scheme_assembler(
    const grid_problem & gp_,
    const time_scheme_config & tsc_
) :
gp(gp_),
tsc(tsc_),
current_time(tsc.start_time)
{
    int initialized = 0;

    MPI_Initialized(&initialized);

    if(initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    }
}



void time_scheme_assembler::setup_solver(
    const space_matrix_assembler & sma
)
{
    int amount = 0;

    if(rank==0)
    amount = sma.get_subdomain_systems().size();

    if(process_count>1)
    MPI_Bcast(&amount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(amount==0)
    throw std::runtime_error("No subdomains");

    if(process_count>1 && process_count<amount+1)
    throw std::runtime_error("Not enough MPI processes");



    subdomains_solvers.clear();
    boundary_solver.reset();

    T = 0.0;

    if(rank==0)
    {
        const auto & systems = sma.get_subdomain_systems();

        vector<int> boundary_sizes;
        boundary_sizes.reserve(systems.size());

        for(const subdomain_system & ss : systems)
        boundary_sizes.push_back(ss.bb.rows());

        boundary_solver = std::make_unique<time_step_boundary>(
            *tsc.scheme,
            sma.get_base_boundary_system(),
            boundary_sizes
        );
    }

    distribute_subdomain_systems(sma);

    restart_solver();
}



void time_scheme_assembler::restart_solver()
{
    current_time = tsc.start_time;
    current_time_step = 0;

    if(rank==0)
    {
        current_values = gp.get_initial_values();

        if(current_time!=0.0)
        {
            vector<subdomain_values> boundary_values =
                gp.get_boundary_values(current_time);

            for(int i=0; i<current_values.size(); i++)
            current_values[i].global_boundary =
                std::move(boundary_values[i].global_boundary);
        }
    }

    distribute_initial_values();

    next_values.clear();
    next_values.resize(current_values.size());

    if(current_time>=tsc.end_time) return;

    double new_T = compute_global_time_step();

    new_T = std::min(new_T, tsc.end_time-current_time);

    if(new_T!=T)
    rebuild_time_step_systems(new_T);
}



double time_scheme_assembler::compute_global_time_step() const
{
    double local_T = std::numeric_limits<double>::infinity();

    for(int i=0; i<subdomains_solvers.size(); i++)
    {
        double subdomain_T =
            subdomains_solvers[i].compute_local_time_step(
                current_values[i]
            );

        if(!std::isfinite(subdomain_T) || subdomain_T<=0.0)
        subdomain_T = 0.0;

        local_T = std::min(local_T, subdomain_T);
    }

    double global_T = local_T;

    if(process_count>1)
    MPI_Allreduce(
        &local_T,
        &global_T,
        1,
        MPI_DOUBLE,
        MPI_MIN,
        MPI_COMM_WORLD
    );

    if(!std::isfinite(global_T) || global_T<=0.0)
    throw std::runtime_error("Invalid time step");

    return global_T;
}



void time_scheme_assembler::rebuild_time_step_systems(
    double new_T
)
{
    for(time_step_subdomain & solver : subdomains_solvers)
    solver.setup_time_step_solver(new_T);

    vector<vector<triplet>> schur =
        collect_schur_complements();

    if(rank==0)
    boundary_solver->build_boundary_system(schur, new_T);

    T = new_T;
}



void time_scheme_assembler::time_step()
{
    double remaining_time = tsc.end_time-current_time;

    if(remaining_time<=0.0) return;

    double new_T = T;

    if(!tsc.scheme->static_time_step)
    new_T = compute_global_time_step();

    new_T = std::min(new_T, remaining_time);

    if(new_T!=T)
    rebuild_time_step_systems(new_T);



    double next_time = current_time+T;

    if(T==remaining_time)
    next_time = tsc.end_time;

    if(rank==0)
    next_values = gp.get_boundary_values(next_time);

    distribute_global_boundary_values();



    vector<vect> rhs = collect_schur_rhs();

    if(rank==0)
    boundary_solver->solve(
        current_values,
        next_values,
        rhs
    );

    distribute_local_boundary_values();



    for(int i=0; i<subdomains_solvers.size(); i++)
    subdomains_solvers[i].solve(
        current_values[i],
        next_values[i]
    );

    collect_inner_values();



    current_values.swap(next_values);

    current_time = next_time;
    current_time_step++;
}