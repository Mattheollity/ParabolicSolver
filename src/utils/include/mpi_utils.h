#ifndef MPI_UTILS_
#define MPI_UTILS_

#include "common.h"
#include <mpi.h>

vector<double> serializeSparseMatrixFlat(spmat & target);
spmat deserializeSparseMatrixFlat(const std::vector<double>& buffer);
void sendMatrix(const spmat & target, int dest, int tag, MPI_Comm comm);
void recvMatrix(spmat & target, int src, int tag, MPI_Comm comm);



#endif