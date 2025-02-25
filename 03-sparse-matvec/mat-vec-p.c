#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "files.h"

#define MatrixFile1 "./data/bcsstk03/bcsstk03.mtx"
#define MatrixFile2 "./data/cfd1/cfd1.mtx"

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

    int n1, nnz1, n2, nnz2;
    int *index1 = NULL, *col_id1 = NULL, *index2 = NULL, *col_id2 = NULL;
    double *val1 = NULL, *val2 = NULL, *x1 = NULL, *y1 = NULL, *x2 = NULL, *y2 = NULL;

    int n_para = atoi(argv[1]);

    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_para);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // read the data
    if (rank == 0)
    {
        read_size(MatrixFile1, &n1, &nnz1);
        read_size(MatrixFile2, &n2, &nnz2);
        
        nnz1 = nnz1 * 2 - n1;
        nnz2 = nnz2 * 2 - n2;
    
        x1 = (double *) malloc(n1 * sizeof(double));
        y1 = (double *) malloc(n1 * sizeof(double));
        col_id1 = (int *) malloc(nnz1 * sizeof(int));
        index1 = (int *) malloc((n1 + 1) * sizeof(int));
        val1 = (double *) malloc(nnz1 * sizeof(double));    

        x2 = (double *) malloc(n2 * sizeof(double));
        y2 = (double *) malloc(n2 * sizeof(double));
        col_id2 = (int *) malloc(nnz2 * sizeof(int));
        index2 = (int *) malloc((n2 + 1) * sizeof(int));
        val2 = (double *) malloc(nnz2 * sizeof(double));
   
        read_matrix(MatrixFile1, n1, nnz1, index1, col_id1, val1);
        read_matrix(MatrixFile2, n2, nnz2, index2, col_id2, val2);

        gen_vector("vector1.txt", n1, x1);
        gen_vector("vector2.txt", n2, x2);
    }

    // distribute the data
    MPI_Bcast(&n1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) 
    {
        x1 = malloc(n1 * sizeof(double));
    }
    MPI_Bcast(x1, n1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&n2, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) 
    {
        x2 = malloc(n2 * sizeof(double));
    }
    MPI_Bcast(x2, n2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // local calculation
    int n_local1 = local_size(n1, rank, n_para);
    int begin1 = local_begin(n1, rank, n_para);
    int n_local2 = local_size(n2, rank, n_para);
    int begin2 = local_begin(n2, rank, n_para);

    double *y_local1 = malloc(n_local1 * sizeof(double));
    double *y_local2 = malloc(n_local2 * sizeof(double));

    double global_start1 = get_time();
    matvec_block(n_local1, begin1, index1 + begin1, col_id1, val1, x1, y_local1);
    MPI_Barrier(MPI_COMM_WORLD);
    double global_end1 = get_time();
    
    double global_start2 = get_time();
    matvec_block(n_local2, begin2, index2 + begin2, col_id2, val2, x2, y_local2);
    MPI_Barrier(MPI_COMM_WORLD);
    double global_end2 = get_time();

    // collect the results and print in a txt file
    int *recvcounts1 = NULL, *displs1 = NULL;
    if (rank == 0) 
    {
        recvcounts1 = malloc(n_para * sizeof(int));
        displs1 = malloc(n_para * sizeof(int));
        int offset1 = 0;
        for (int i = 0; i < n_para; i++) 
        {
            recvcounts1[i] = local_size(n1, i, n_para);
            displs1[i] = offset1;
            offset1 += recvcounts1[i];
        }
    }
    MPI_Gatherv(y_local1, n_local1, MPI_DOUBLE, y1, recvcounts1, displs1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int *recvcounts2 = NULL, *displs2 = NULL;
    if (rank == 0) 
    {
        recvcounts2 = malloc(n_para * sizeof(int));
        displs2 = malloc(n_para * sizeof(int));
        int offset2 = 0;
        for (int i = 0; i < n_para; i ++) 
        {
            recvcounts2[i] = local_size(n2, i, n_para);
            displs2[i] = offset2;
            offset2 += recvcounts2[i];
        }
    }
    MPI_Gatherv(y_local2, n_local2, MPI_DOUBLE, y2, recvcounts2, displs2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // local free
    free(x1);
    free(y_local1); 
    free(x2); 
    free(y_local2);

    // print and global free
    if (rank == 0)
    {
        printf("Time;\tNumber of proscesses;\tMatrix size;\tNon-zeros\n");
        printf("---------------------------------------------------------\n");
        printf("%.6f\t%d\t%d\t%d\n", global_end1 - global_start1, n_para, n1, nnz1);
        printf("%.6f\t%d\t%d\t%d\n", global_end2 - global_start2, n_para, n2, nnz2);
        write_vector("parallel_result1.txt", n1, y1);
        write_vector("parallel_result2.txt", n2, y2);

        free(y1); 
        free(col_id1); 
        free(index1); 
        free(val1);

        free(y2); 
        free(col_id2); 
        free(index2); 
        free(val2);

        free(recvcounts1); 
        free(displs1); 

        free(recvcounts2); 
        free(displs2);
    }

    MPI_Finalize();
    return 0;
}