#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "files.h"

#define MatrixFile "./data/bcsstk03/bcsstk03.mtx"
//#define MatrixFile "./data/cfd1/cfd1.mtx"

#define VectorFile "./vector1.txt"
//#define VectorFile "./vector2.txt"

#define ResultFile "./parallel_result1"
//#define ResultFile "./parallel_result2"

// timer
double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

void write_vector_p(int size, int n, double *y)
{
    char filename[MAX_SIZE];
    sprintf(filename, "%s_%d.txt", ResultFile, size);
    write_vector(filename, n, y);
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 0, nnz = 0;
    int *index = NULL, *col_id = NULL;
    double *val = NULL, *x = NULL, *y = NULL;

    if (rank == 0)
    {
        read_size(MatrixFile, &n, &nnz);
        nnz = nnz * 2 - n;
    }
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int n_local = n / size;
    int remainder = n % size;
    int block_start = rank * n_local + (rank < remainder ? rank : remainder);
    if (rank < remainder)
        n_local++;

    index = (int *)malloc((n + 1) * sizeof(int));
    col_id = (int *)malloc(nnz * sizeof(int));
    val = (double *)malloc(nnz * sizeof(double));
    x = (double *)malloc(n * sizeof(double));
    y = (double *)malloc(n * sizeof(double));
    double *y_local = (double *)malloc(n_local * sizeof(double));

    if (rank == 0)
    {
        read_matrix(MatrixFile, n, nnz, index, col_id, val);
        read_vector(VectorFile, n, x);
    }
    
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    double start = get_time();
    matvec_block(n, n_local, block_start, index, col_id, val, x, y_local);
    MPI_Gather(y_local, n_local, MPI_DOUBLE, y, n_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double end = get_time();
    
    if (rank == 0)
{
    printf("%.10f\t%d\t%d\t%d\n", end - start, size, n, nnz);
    write_vector_p(size, n, y);
}
    
    free(index);
    free(col_id);
    free(val);
    free(x);
    free(y);
    free(y_local);
    
    MPI_Finalize();
    return 0;
}
