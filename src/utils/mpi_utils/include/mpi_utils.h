#ifndef MPI_UTILS_
#define MPI_UTILS_

#include "common.h"
#include "eigen_parts.h"

#include <mpi.h>


vector<double> serializeSparseMatrixFlat(const spmat & target);

spmat deserializeSparseMatrixFlat(const vector<double> & buffer);



void sendMatrix(
    const spmat & target,
    int dest,
    int tag,
    MPI_Comm comm
);

void recvMatrix(
    spmat & target,
    int src,
    int tag,
    MPI_Comm comm
);



void sendVector(
    const vect & target,
    int dest,
    int tag,
    MPI_Comm comm
);

void recvVector(
    vect & target,
    int src,
    int tag,
    MPI_Comm comm
);



void sendTriplets(
    const vector<triplet> & target,
    int dest,
    int tag,
    MPI_Comm comm
);

void recvTriplets(
    vector<triplet> & target,
    int src,
    int tag,
    MPI_Comm comm
);

#endif