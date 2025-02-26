#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "files.h"

#define MatrixFile1 "./data/bcsstk03/bcsstk03.mtx"
#define MatrixFile2 "./data/cfd1/cfd1.mtx"

#define VectorFile1 "./vector1.txt"
#define VectorFile2 "./vector2.txt"

#define ResultFile1 "./parallel_result1.txt"
#define ResultFile2 "./parallel_result2.txt"

// timer
double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

// y_seg = A_block * x ( A in CSR format)
void matvec_block(int n_local, int block_begin, int *index_block, int *col_id_block, double *val_block, double *x_local, double *y)
{
    for (int i = block_begin; i < block_begin + n_local; i ++)
    {
        y[i] = 0;
        for (int j = index_block[i]; j < index_block[i + 1]; j ++)
        {
            y[i] += val_block[col_id_block[j]] * x_local[j];
        }
    }
}

// compute local size
int local_size(int n, int rank, int n_para)
{
    int n_local = n / n_para;
    int n_extra = n % n_para;
    return n_local + (rank < n_extra ? 1 : 0);
}

// compute local begin
int local_begin(int n, int rank, int n_para)
{
    int n_local = n / n_para;
    int n_extra = n % n_para;
    return rank * n_local + (rank < n_extra ? rank : n_extra);
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        printf("Please enter the number of processes\n");
        return 1;
    }   

    int n = 0, nnz = 0;
    int *index = NULL, *col_id = NULL;
    double *val = NULL, *x = NULL, *y = NULL;

    int n_para = atoi(argv[1]);

    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_para);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // read the data
    if (rank == 0)
    {
        read_size(MatrixFile1, &n, &nnz);
        
        nnz = nnz * 2 - n;
    
        x = (double *) malloc(n * sizeof(double));
        y = (double *) malloc(n * sizeof(double));
        col_id = (int *) malloc(nnz * sizeof(int));
        index = (int *) malloc((n + 1) * sizeof(int));
        val = (double *) malloc(nnz * sizeof(double));    
   
        read_matrix(MatrixFile1, n, nnz, index, col_id, val);

        read_vector(VectorFile1, n, x);
    }

    // distribute the data
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) 
    {
        x = malloc(n * sizeof(double));
    }
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // local calculation
    int n_local = local_size(n, rank, n_para);
    int begin = local_begin(n, rank, n_para);

    double *y_local = malloc(n_local * sizeof(double));

    double global_start = get_time();
    matvec_block(n_local, begin, index + begin, col_id, val, x, y_local);
    MPI_Barrier(MPI_COMM_WORLD);
    double global_end = get_time();

    // collect the results and print in a txt file
    int *recvcounts = NULL, *displs = NULL;
    if (rank == 0) 
    {
        recvcounts = malloc(n_para * sizeof(int));
        displs = malloc(n_para * sizeof(int));
        int offset = 0;
        for (int i = 0; i < n_para; i++) 
        {
            recvcounts[i] = local_size(n, i, n_para);
            displs[i] = offset;
            offset += recvcounts[i];
        }
    }
    MPI_Gatherv(y_local, n_local, MPI_DOUBLE, y, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // local free
    free(x);
    free(y_local); 

    // print and global free
    if (rank == 0)
    {
        printf("Time;\tNumber of proscesses;\tMatrix size;\tNon-zeros\n");
        printf("---------------------------------------------------------\n");
        printf("%.6f\t%d\t%d\t%d\n", global_end - global_start, n_para, n, nnz);
        write_vector(ResultFile1, n, y);

        free(y); 
        free(col_id); 
        free(index); 
        free(val);

        free(recvcounts); 
        free(displs); 

    }

    MPI_Finalize();
    return 0;
}