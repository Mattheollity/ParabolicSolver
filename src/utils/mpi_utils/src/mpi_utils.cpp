#include "mpi_utils.h"


vector<double> serializeSparseMatrixFlat(const spmat & target)
{
    vector<double> buffer(
        2 + target.nonZeros()*2 + target.cols()
    );

    buffer[0] = target.rows();
    buffer[1] = target.cols();

    int position = 2;

    for(int col=0; col<target.outerSize(); col++)
    {
        int nnz_col = 0;
        int col_pos = position++;

        for(spmat::InnerIterator it(target, col); it; ++it)
        {
            buffer[position++] = it.row();
            buffer[position++] = it.value();

            nnz_col++;
        }

        buffer[col_pos] = nnz_col;
    }

    return buffer;
}



spmat deserializeSparseMatrixFlat(const vector<double> & buffer)
{
    int rows = buffer[0];
    int cols = buffer[1];

    spmat target(rows, cols);

    vector<triplet> trp;
    trp.reserve((buffer.size()-2-cols)/2);

    int position = 2;

    for(int col=0; col<cols; col++)
    {
        int nnz_col = buffer[position++];

        for(int k=0; k<nnz_col; k++)
        {
            int row = buffer[position++];
            double value = buffer[position++];

            trp.emplace_back(row, col, value);
        }
    }

    target.setFromTriplets(trp.begin(), trp.end());

    return target;
}



void sendMatrix(
    const spmat & target,
    int dest,
    int tag,
    MPI_Comm comm
)
{
    vector<double> buffer = serializeSparseMatrixFlat(target);

    int n = buffer.size();

    MPI_Send(&n, 1, MPI_INT, dest, tag, comm);

    MPI_Send(
        buffer.data(), n, MPI_DOUBLE,
        dest, tag+1, comm
    );
}



void recvMatrix(
    spmat & target,
    int src,
    int tag,
    MPI_Comm comm
)
{
    int n = 0;

    MPI_Recv(
        &n, 1, MPI_INT,
        src, tag, comm, MPI_STATUS_IGNORE
    );

    vector<double> buffer(n);

    MPI_Recv(
        buffer.data(), n, MPI_DOUBLE,
        src, tag+1, comm, MPI_STATUS_IGNORE
    );

    target = deserializeSparseMatrixFlat(buffer);
}



void sendVector(
    const vect & target,
    int dest,
    int tag,
    MPI_Comm comm
)
{
    int n = target.size();

    MPI_Send(&n, 1, MPI_INT, dest, tag, comm);

    MPI_Send(
        target.data(), n, MPI_DOUBLE,
        dest, tag+1, comm
    );
}



void recvVector(
    vect & target,
    int src,
    int tag,
    MPI_Comm comm
)
{
    int n = 0;

    MPI_Recv(
        &n, 1, MPI_INT,
        src, tag, comm, MPI_STATUS_IGNORE
    );

    target.resize(n);

    MPI_Recv(
        target.data(), n, MPI_DOUBLE,
        src, tag+1, comm, MPI_STATUS_IGNORE
    );
}



void sendTriplets(
    const vector<triplet> & target,
    int dest,
    int tag,
    MPI_Comm comm
)
{
    int n = target.size();

    vector<double> buffer(3*n);

    for(int i=0; i<n; i++)
    {
        buffer[3*i] = target[i].row();
        buffer[3*i+1] = target[i].col();
        buffer[3*i+2] = target[i].value();
    }

    MPI_Send(&n, 1, MPI_INT, dest, tag, comm);

    MPI_Send(
        buffer.data(), 3*n, MPI_DOUBLE,
        dest, tag+1, comm
    );
}



void recvTriplets(
    vector<triplet> & target,
    int src,
    int tag,
    MPI_Comm comm
)
{
    int n = 0;

    MPI_Recv(
        &n, 1, MPI_INT,
        src, tag, comm, MPI_STATUS_IGNORE
    );

    vector<double> buffer(3*n);

    MPI_Recv(
        buffer.data(), 3*n, MPI_DOUBLE,
        src, tag+1, comm, MPI_STATUS_IGNORE
    );

    target.clear();
    target.reserve(n);

    for(int i=0; i<n; i++)
    target.emplace_back(
        static_cast<int>(buffer[3*i]),
        static_cast<int>(buffer[3*i+1]),
        buffer[3*i+2]
    );
}