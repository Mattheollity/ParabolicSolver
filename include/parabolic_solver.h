#ifndef PARABOLIC_SOLVER_
#define PARABOLIC_SOLVER_

#include "configs_io.h"
#include "mpi_utils.h"

class parabolic_solver
{
    case_config cc;
    unique_ptr<grid> g;
    unique_ptr<space_matrix_assembler> ss;
    unique_ptr<time_scheme_assembler> ts;

    int MPI_is_active;
    int thread;


    public:

    parabolic_solver(string CASE_DIR);
};

#endif