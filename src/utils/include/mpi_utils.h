#ifndef MPI_UTILS_
#define MPI_UTILS_

#include "common.h"
#include <mpi.h>

vector<double> serializeSparseMatrixFlat(spmat & target);
spmat deserializeSparseMatrixFlat(const std::vector<double>& buffer);




#endif