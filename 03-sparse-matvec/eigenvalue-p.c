#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "files.h"

#define MatrixFile "./data/bcsstk03/bcsstk03.mtx"
#define VectorFile "./vector1.txt"
#define IterFile "./parallel_iter1"
#define EigFile "./parallel_eig1"

//#define MatrixFile "./data/cfd1/cfd1.mtx"
//#define VectorFile "./vector2.txt"
//#define IterFile "./parallel_iter2"
//#define EigFile "./parallel_eig2"

#define EPSILON 1e-6

// timer
double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

double distributed_eigenvalue(
    int n_local, int global_n,
    int *index, int *col_id, double *val,
    double *x, double *x_new, MPI_Comm comm,
    const char* iter_filename) 
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double lambda_max = 0.0, r = 1.0;
    int count = 0;
    double *diff = (double *)malloc(global_n * sizeof(double));
    double *x_local = (double *)malloc(n_local * sizeof(double));
    double *recvbuf = (rank == 0) ? malloc(global_n * sizeof(double)) : NULL;

    int *recvcounts = (int *)malloc(size * sizeof(int));
    int *displs = (int *)malloc(size * sizeof(int));

    MPI_Allgather(&n_local, 1, MPI_INT, recvcounts, 1, MPI_INT, comm);
    displs[0] = 0;
    for (int i = 1; i < size; i++) displs[i] = displs[i-1] + recvcounts[i-1];

    normalize(global_n, x);

    FILE *fp = NULL;
    if (rank == 0) 
    {
        fp = fopen(iter_filename, "w");
        if (!fp) 
        {
            perror("Error opening iteration file");
            MPI_Abort(comm, 1);
        }
    }

    while (r > EPSILON) 
    {
        matvec(n_local, index, col_id, val, x, x_local);
        MPI_Allgatherv(x_local, n_local, MPI_DOUBLE, x_new, recvcounts, displs, MPI_DOUBLE, comm);

        normalize(global_n, x_new);

        for (int i = 0; i < global_n; i ++) 
        {
            diff[i] = x_new[i] - x[i];
        }

        double local_r = dot_vec(global_n, diff, diff);
        MPI_Allreduce(&local_r, &r, 1, MPI_DOUBLE, MPI_SUM, comm);
        r = sqrt(r);

        if (rank == 0) 
        {
            fprintf(fp, "%d\t%.10f\n", ++count, r);
        }

        for (int i = 0; i < global_n; i ++) 
        {
            x[i] = x_new[i];
        }
    }

    matvec(n_local, index, col_id, val, x_new, x_local);
    MPI_Allgatherv(x_local, n_local, MPI_DOUBLE, x, recvcounts, displs, MPI_DOUBLE, comm);

    double local_num = dot_vec(n_local, x_local, x_new + displs[rank]);
    double global_num = 0.0;
    MPI_Allreduce(&local_num, &global_num, 1, MPI_DOUBLE, MPI_SUM, comm);

    double denom = dot_vec(global_n, x_new, x_new);
    MPI_Allreduce(MPI_IN_PLACE, &denom, 1, MPI_DOUBLE, MPI_SUM, comm);

    lambda_max = global_num / denom;

    if (rank == 0 && fp) 
    {
        fclose(fp);
    }

    free(diff); 
    free(x_local); 
    free(recvcounts); 
    free(displs);

    if (recvbuf) 
    {
        free(recvbuf);
    }

    return lambda_max;
}

void distribute_csr(
    int *global_index, int *global_col, double *global_val,
    int n,
    int **local_index, int **local_col, double **local_val,
    int *n_local, int *nnz_local,
    MPI_Comm comm) 
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int *sendcounts = NULL, *displs = NULL;
    if (rank == 0) 
    {
        sendcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
    }

    int rows_per_proc = n / size;
    int rem = n % size;
    int row_start = rank * rows_per_proc + (rank < rem ? rank : rem);
    int row_end = row_start + rows_per_proc + (rank < rem);
    *n_local = row_end - row_start;

    int start_nnz = 0, end_nnz = 0;
    if (rank == 0) 
    {
        for (int r = 0; r < size; r ++) 
        {
            int rs = r * rows_per_proc + (r < rem ? r : rem);
            int re = rs + rows_per_proc + (r < rem);
            sendcounts[r] = global_index[re] - global_index[rs];
            displs[r] = global_index[rs];
        }
    }

    MPI_Scatter(sendcounts, 1, MPI_INT, nnz_local, 1, MPI_INT, 0, comm);

    *local_index = (int *)malloc((*n_local + 1) * sizeof(int));
    *local_col = (int *)malloc(*nnz_local * sizeof(int));
    *local_val = (double *)malloc(*nnz_local * sizeof(double));

    if (rank == 0) 
    {
        for (int r = 0; r < size; r ++) 
        {
            int rs = r * rows_per_proc + (r < rem ? r : rem);
            int offset = global_index[rs];
            for (int i = 0; i <= (r < rem ? rows_per_proc + 1 : rows_per_proc); i ++) 
            {
                global_index[rs + i] -= offset;
            }
        }
        for (int i = 0; i <= *n_local; i ++) 
        {
            (*local_index)[i] = global_index[row_start + i];
        }
    }
    MPI_Bcast(*local_index, *n_local + 1, MPI_INT, 0, comm);
    MPI_Scatterv(global_col, sendcounts, displs, MPI_INT, *local_col, *nnz_local, MPI_INT, 0, comm);
    MPI_Scatterv(global_val, sendcounts, displs, MPI_DOUBLE, *local_val, *nnz_local, MPI_DOUBLE, 0, comm);

    if (rank == 0) 
    {
        free(sendcounts);
        free(displs);
    }
}

int main(int argc, char **argv) 
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char iter_filename[256], eig_filename[256];
    if (rank == 0) 
    {
        sprintf(iter_filename, "%s_%d.txt", IterFile, size);
        sprintf(eig_filename, "%s_%d.txt", EigFile, size);
    }

    int n = 0, nnz = 0;
    int *index = NULL, *col_id = NULL;
    double *val = NULL, *x = NULL, *x_new = NULL;

    int *index_block = NULL, *col_id_block = NULL;
    double *val_block = NULL;
    int n_local = 0, nnz_local = 0;

    if (rank == 0) 
    {
        read_size(MatrixFile, &n, &nnz);
        nnz = nnz * 2 - n;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    x = (double *)malloc(n * sizeof(double));
    x_new = (double *)malloc(n * sizeof(double));
    read_vector(VectorFile, n, x);

    if (rank == 0) 
    {
        index = (int *)malloc((n + 1) * sizeof(int));
        col_id = (int *)malloc(nnz * sizeof(int));
        val = (double *)malloc(nnz * sizeof(double));
        read_matrix(MatrixFile, n, nnz, index, col_id, val);
    }

    distribute_csr(index, col_id, val, n, &index_block, &col_id_block, &val_block, &n_local, &nnz_local, MPI_COMM_WORLD);

    double start = get_time();
    double lambda = distributed_eigenvalue(n_local, n, index_block, col_id_block, val_block, x, x_new, MPI_COMM_WORLD, iter_filename);
    double end = get_time();

    if (rank == 0) 
    {
        FILE *fp = fopen(eig_filename, "w");
        if (!fp) 
        {
            perror("Error opening EigFile");
            exit(1);
        }
        fprintf(fp, "%.10f\n", lambda);
        fclose(fp);
        printf("%.10f\t%d\t%d\t%d\n", end - start, size, n, nnz);

        free(index); 
        free(col_id); 
        free(val);
    }

    free(x); 
    free(x_new);
    free(index_block); 
    free(col_id_block); 
    free(val_block);

    MPI_Finalize();
    return 0;
}
