#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "files.h"

#define MatrixFile "./data/bcsstk03/bcsstk03.mtx"
//#define MatrixFile "./data/cfd1/cfd1.mtx"

#define VectorFile "./vector1.txt"
//#define VectorFile "./vector2.txt"

#define PerfFile "./n1.txt"
//#define PerfFile "./n2.txt"

// timer
double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
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

    int *index_block = NULL, *col_id_block = NULL;
    double *val_block = NULL, *y_local = NULL;
    int n_local = 0, nnz_local = 0;

    int *sendc_ind = NULL, *sendc_v = NULL, *disp_ind = NULL, *disp_v = NULL;
    int *recvc_y = NULL, *disp_y = NULL;

    if (rank == 0)
    {
        read_size(MatrixFile, &n, &nnz);
        nnz = nnz * 2 - n;
    }
        
    // broadcast n
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    x = (double *)malloc(n * sizeof(double));

    if (rank == 0)
    {
        int remainder = n % size;
        n_local = n / size + (rank < remainder);
        printf("compute rank0: n = %d, n_local = %d\n", n, n_local);

        sendc_ind = (int *)malloc(size * sizeof(int));
        sendc_v = (int *)malloc(size * sizeof(int));
        recvc_y = (int *)malloc(size * sizeof(int));
        disp_y = (int *)malloc(size * sizeof(int));
        disp_ind = (int *)malloc(size * sizeof(int));
        disp_v = (int *)malloc(size * sizeof(int));

        index = (int *)malloc((n + 1) * sizeof(int));
        col_id = (int *)malloc(nnz * sizeof(int));
        val = (double *)malloc(nnz * sizeof(double));     
        y = (double *)malloc(n * sizeof(double));

        read_matrix(MatrixFile, n, nnz, index, col_id, val);

        // compute n_loacl for index_block and nnz_local for col_id_block and val_block
        // send n_local and nnz_local to other processes
        for (int r = 0; r < size; r ++)
        {
            n_local = n / size + (r < remainder);
            int block_start = r * n_local + (r < remainder ? r : remainder);
            sendc_ind[r] = n_local;
            disp_ind[r] = block_start;           
            MPI_Send(&n_local, 1, MPI_INT, r, 10, MPI_COMM_WORLD);

            recvc_y[r] = n_local;
            disp_y[r] = block_start;

            nnz_local = index[block_start + n_local] - index[block_start];
            sendc_v[r] = nnz_local;
            disp_v[r] = index[block_start];
            MPI_Send(&nnz_local, 1, MPI_INT, r, 20, MPI_COMM_WORLD);
            printf("send rank %d: block_start: %d, index[f] - index[i]: (%d - %d), n_local = %d, nnz_local = %d\n", r, block_start, index[block_start + n_local], index[block_start], n_local, nnz_local);
        }
    }
    else
    {
        MPI_Recv(&n_local, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&nnz_local, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("receive rank %d: n_local = %d, nnz_local = %d\n", rank, n_local, nnz_local);
    }
   
    index_block = (int *)malloc((n_local) * sizeof(int));
    col_id_block = (int *)malloc(nnz_local * sizeof(int));
    val_block = (double *)malloc(nnz_local * sizeof(double));
    y_local = (double *)malloc(n_local * sizeof(double));

    // scatter index_block
    MPI_Scatterv(index, sendc_ind, disp_ind, MPI_INT, index_block, n_local + 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("receive rank %d: index[f] = %d, index[i] = %d, n_local = %d, nnz_local = %d\n", rank, index_block[n_local-1], index_block[0], n_local, nnz_local);

    // read x and copy loacl block of matrix A   
    if (rank == 0)
    {
        read_vector(VectorFile, n, x);

        for (int i = 0; i < sendc_ind[0]; i ++)
        {
            index_block[i] = index[i];
        }
        for (int i = 0; i < index[disp_ind[1]]; i ++)
        {
            col_id_block[i] = col_id[i];
            val_block[i] = val[i];
        }
    }

    // broadcast x
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // normalize index_block
    if (index_block != NULL) 
    {
        int offset = index_block[0];
        for (int i = 0; i < n_local + 1; i++) 
        {
            index_block[i] -= offset;
        }
    }

    // scatter local block of matrix A
    MPI_Scatterv(col_id, sendc_v, disp_v, MPI_INT, col_id_block, nnz_local, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(val, sendc_v, disp_v, MPI_DOUBLE, val_block, nnz_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // compute y_local and measure time
    double start = get_time();
    matvec(n_local, index_block, col_id_block, val_block, x, y_local);
    MPI_Barrier(MPI_COMM_WORLD);
    double end = get_time();

    // gather y_local to y
    MPI_Gatherv(y_local, n_local, MPI_DOUBLE, y, recvc_y, disp_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
    if (rank == 0)
    {
        FILE *fp = fopen(PerfFile, "a");
        if (!fp)
        {
            perror("Error opening file: main");
            exit(1);
        }
        fprintf(fp, "%.10f\t%d\t%d\t%d\n", end - start, size, n, nnz);
        fclose(fp);

        for (int i = 0; i < n; i ++)
        {
            printf("%.10f\n", y[i]);
        }
    }

    if (rank == 0)
    {
        free(index);
        free(col_id);
        free(val);
        free(y);
        free(sendc_ind);
        free(sendc_v);
        free(recvc_y);
        free(disp_y);
        free(disp_ind);
        free(disp_v);
    }
    
    free(index_block);
    free(col_id_block);
    free(val_block);
    free(x);
    free(y_local);
    
    MPI_Finalize();
    return 0;
}
