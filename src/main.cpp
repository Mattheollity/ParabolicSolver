#include "parabolic_solver.h"

#include <mpi.h>
#include <iostream>


int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);

    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc!=2)
    {
        if(rank==0)
        std::cerr<<"Usage: "<<argv[0]<<" <case_directory>\n";

        MPI_Finalize();

        return 1;
    }

    {
        parabolic_solver solver(argv[1]);

        solver.setup_solver();

        solver.solve();
    }

    MPI_Finalize();

    return 0;
}