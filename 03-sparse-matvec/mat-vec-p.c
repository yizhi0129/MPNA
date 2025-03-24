#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "files.h"

//#define MatrixFile "./data/bcsstk03/bcsstk03.mtx"
//#define VectorFile "./vector1.txt"
//#define PerfFile "./mat_vec_time_n1.txt"

#define MatrixFile "./data/cfd1/cfd1.mtx"
#define VectorFile "./vector2.txt"
#define PerfFile "./mat_vec_time_n2.txt"

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

    int *sendc_v = NULL, *disp_v = NULL;
    int *recvc_y = NULL, *disp_y = NULL;

    double start = 0.0, end = 0.0;

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
        index = (int *)malloc((n + 1) * sizeof(int));
        col_id = (int *)malloc(nnz * sizeof(int));
        val = (double *)malloc(nnz * sizeof(double));     
        y = (double *)malloc(n * sizeof(double));

        sendc_v = (int *)malloc(size * sizeof(int));
        disp_v = (int *)malloc(size * sizeof(int));
    
        recvc_y = (int *)malloc(size * sizeof(int));
        disp_y = (int *)malloc(size * sizeof(int));

        read_matrix(MatrixFile, n, nnz, index, col_id, val);
       
        for (int r = 0; r < size; r ++)
        {
            // compute n_loacl for index_block and nnz_local for col_id_block and val_block
            int n_r = n / size + (r < (n % size));
            int block_start = r * (n / size) + (r < (n % size) ? r : (n % size)); 

            int nnz_r = index[block_start + n_r] - index[block_start];

            sendc_v[r] = nnz_r;
            disp_v[r] = index[block_start];

            recvc_y[r] = n_r;
            disp_y[r] = block_start;

            // send n_local, nnz_local and index_block to other processes
            MPI_Send(&n_r, 1, MPI_INT, r, 10, MPI_COMM_WORLD);
            MPI_Send(&nnz_r, 1, MPI_INT, r, 20, MPI_COMM_WORLD);

            int *index_tmp = (int *)malloc((n_r + 1) * sizeof(int));
            int offset = index[block_start];
            for (int i = 0; i <= n_r; i ++) 
            {
                index_tmp[i] = index[block_start + i] - offset;
            }
            MPI_Send(index_tmp, n_r + 1, MPI_INT, r, 30, MPI_COMM_WORLD);
            free(index_tmp);
        }
        
        // local index for rank 0
        n_local = n / size + (rank < (n % size));
        nnz_local = index[n_local] - index[0];
        index_block = (int *)malloc((n_local + 1) * sizeof(int));
        for (int i = 0; i < n_local + 1; i ++)
        {
            index_block[i] = index[i];
        }

        // local col_id and local val for rank 0
        col_id_block = (int *)malloc(nnz_local * sizeof(int));
        val_block = (double *)malloc(nnz_local * sizeof(double));
        for (int i = 0; i < nnz_local; i ++)
        {
            col_id_block[i] = col_id[i];
            val_block[i] = val[i];
        }

        read_vector(VectorFile, n, x);
    }
    else
    {
        // receive n_local and nnz_local from rank 0
        MPI_Recv(&n_local, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&nnz_local, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        

        // receive local index from rank 0 
        index_block = (int *)malloc((n_local + 1) * sizeof(int));
        MPI_Recv(index_block, n_local + 1, MPI_INT, 0, 30, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // malloc local col_id and local val from rank 0
        if (nnz_local > 0) 
        {
            col_id_block = (int *)malloc(nnz_local * sizeof(int));
            val_block = (double *)malloc(nnz_local * sizeof(double));
        }
    }
   
    y_local = (double *)malloc(n_local * sizeof(double));

    // broadcast x
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // scatter col_id_block and val_block
    MPI_Scatterv(col_id, sendc_v, disp_v, MPI_INT, col_id_block, nnz_local, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(val, sendc_v, disp_v, MPI_DOUBLE, val_block, nnz_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // compute y_local and measure time
    if (rank == 0) start = get_time();
    matvec(n_local, index_block, col_id_block, val_block, x, y_local);
    if (rank == 0) end = get_time();


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
        fprintf(fp, "%.10e\t%d\t%d\t%d\n", end - start, size, n, nnz);
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
        free(sendc_v);
        free(disp_v);
        free(recvc_y);
        free(disp_y);
    }
    
    free(index_block);
    free(col_id_block);
    free(val_block);
    free(x);
    free(y_local);
    
    MPI_Finalize();
    return 0;
}
