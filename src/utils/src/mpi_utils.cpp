#include "mpi_utils.h"
#include <mpi.h>

vector<double> serializeSparseMatrixFlat(const spmat & target) {
    vector<double> buffer;
    buffer.resize(2 + target.nonZeros() * 2 + target.cols());

    buffer[0] = target.rows();
    buffer[1] = target.cols();

    int position=2;
    for (int col = 0; col < target.outerSize(); ++col) {
        int nnz_col = 0;

        int col_pos = position;
        position++;

        for (Eigen::SparseMatrix<double>::InnerIterator it(target, col); it; ++it) {
            buffer[position] = it.row();
            buffer[position+1] = it.value();
            position+=2;
            ++nnz_col;
        }

        buffer[col_pos]=nnz_col;
    }
    return buffer;
}

spmat deserializeSparseMatrixFlat(const std::vector<double>& buffer) {
    int rows = buffer[0];
    int cols = buffer[1];
    spmat target(rows, cols);

    size_t idx = 2;
    int trp_idx = 0;
    std::vector<Eigen::Triplet<double>> triplets ((buffer.size()-2-cols)/2);

    for (int col = 0; col < cols; ++col) {
        int nnz_col = buffer[idx];
        idx++;
        for (int k = 0; k < nnz_col; ++k) {
            int row = buffer[idx];
            double val = buffer[idx+1];
            idx+=2;
            triplets[trp_idx]={row, col, val};
            trp_idx++;
        }
    }
    target.setFromTriplets(triplets.begin(), triplets.end());
    return target;
}

void sendMatrix(const spmat & target, int dest, int tag, MPI_Comm comm) {
    auto buf = serializeSparseMatrixFlat(target);
    int n = buf.size();
    MPI_Send(&n, 1, MPI_INT, dest, tag, comm);
    MPI_Send(buf.data(), n, MPI_DOUBLE, dest, tag + 1, comm);
}

void recvMatrix(spmat & target, int src, int tag, MPI_Comm comm) {
    int n = 0;
    MPI_Recv(&n, 1, MPI_INT, src, tag, comm, MPI_STATUS_IGNORE);
    std::vector<double> buf(n);
    MPI_Recv(buf.data(), n, MPI_DOUBLE, src, tag + 1, comm, MPI_STATUS_IGNORE);
    target = deserializeSparseMatrixFlat(buf);
}