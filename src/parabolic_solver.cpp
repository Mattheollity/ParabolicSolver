#include "parabolic_solver.h"

parabolic_solver::parabolic_solver(string CASE_DIR) :
    cc(load_case(CASE_DIR)),
    thread(0)
{
    MPI_Initialized(&MPI_is_active);
    if(MPI_is_active) MPI_Comm_rank(MPI_COMM_WORLD, &thread);
    if(thread==0)
    {
        g=std::make_unique<grid>(cc.gc);
        ss=std::make_unique<space_matrix_assembler>(cc.ssc);
    }
    ts = std::make_unique<time_scheme_assembler>(cc.tsc);


    for(int i=0;i<cc.tsc.)
}