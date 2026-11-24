#ifndef TIME_SCHEME_ASSEMBLER_
#define TIME_SCHEME_ASSEMBLER_

#include "mpi_utils.h"
#include "space_matrix_assembler.h"
#include "time_step_boundary.h"
#include "time_step_subdomain.h"

class time_scheme_assembler
{   
    int MPI_is_active;
    int thread_num;
    vector<vector<int>> subdomains_per_MPI_threads;

    vector<time_step_subdomain> subdomains_solvers;
    time_step_boundary boundary_solver;
    int current_time_step = 0;



    void add_subdomain(const subdomain_system &);

    void reset_MPI(const vector<subdomain_system> & systems);

    void send_subdomain_systems(const vector<subdomain_system> & systems);
    void recv_subdomain_systems();

    void send_subdomain_values(const vector<subdomain_values> & values);
    void recv_subdomain_values();

    



    public:

    const time_scheme_config & tsc;

    time_scheme_assembler(const time_scheme_config & tsc_);
    void setup_solver(const vector<subdomain_system> & systems);
    void restart_solver(const vector<subdomain_values> & start_value);
    void time_step(vector<subdomain_values> &);
};

#endif