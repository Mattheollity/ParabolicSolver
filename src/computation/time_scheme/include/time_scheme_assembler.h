#ifndef TIME_SCHEME_ASSEMBLER_
#define TIME_SCHEME_ASSEMBLER_

#include "common.h"
#include "config_structures_ts.h"
#include "space_matrix_assembler.h"
#include "time_step_subdomain.h"
#include "time_step_boundary.h"


class time_scheme_assembler
{
    const grid_problem & gp;
    const time_scheme_config & tsc;

    int rank = 0;
    int process_count = 1;

    int current_time_step = 0;

    double current_time;
    double T = 0.0;

    vector<time_step_subdomain> subdomains_solvers;
    unique_ptr<time_step_boundary> boundary_solver;

    vector<subdomain_values> current_values;
    vector<subdomain_values> next_values;



    void distribute_subdomain_systems(
        const space_matrix_assembler & sma
    );

    void distribute_initial_values();

    void rebuild_time_step_systems(double new_T);

    vector<vector<triplet>> collect_schur_complements() const;



    double compute_global_time_step() const;

    void distribute_global_boundary_values();

    vector<vect> collect_schur_rhs() const;

    void distribute_local_boundary_values();

    void collect_inner_values();



    public:

    time_scheme_assembler(
        const grid_problem & gp_,
        const time_scheme_config & tsc_
    );

    void setup_solver(const space_matrix_assembler & sma);

    void restart_solver();

    void time_step();
};

#endif