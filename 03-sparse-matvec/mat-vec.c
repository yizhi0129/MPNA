#include <sys/time.h>
#include <mpi.h>

#include "files.h"

#define MatrixFile1 "./bcsstk03/bcsstk03.mtx"
#define MatrixFile2 "./cfd1/cfd1.mtx"


void matvec(int n, int *index, int *col_id, double *val, double *x, double *y)
{
    for (int i = 0; i < n; i ++)
    {
        y[i] = 0;
        for (int j = index[i]; j < index[i + 1]; j ++)
        {
            y[i] += val[col_id[j]] * x[j];
        }
    }
}

int main(int argc, char** argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n1, nnz1, n2, nnz2;
    int *index1, *col_id1, *index2, *col_id2;
    double *val1, *val2, *x1, *y1, *x2, *y2;

    read_matrix(MatrixFile1, n1, nnz1, index1, col_id1, val1);
    read_matrix(MatrixFile2, n2, nnz2, index2, col_id2, val2);
    
    x1 = (double *) malloc(n1 * sizeof(double));
    x2 = (double *) malloc(n2 * sizeof(double));
    y1 = (double *) malloc(n1 * sizeof(double));
    y2 = (double *) malloc(n2 * sizeof(double));

    gen_vector("vector1.txt", n1, x1);
    gen_vector("vector2.txt", n2, x2);

    matvec(n1, index1, col_id1, val1, x1, y1);

    
    matvec(n2, index2, col_id2, val2, x2, y2);

    return 0;
}