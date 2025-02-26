#include <sys/time.h>
#include <time.h>

#include "files.h"

#define MatrixFile1 "./data/bcsstk03/bcsstk03.mtx"
#define MatrixFile2 "./data/cfd1/cfd1.mtx"

#define VectorFile1 "./vector1.txt"
#define VectorFile2 "./vector2.txt"

#define ResultFile1 "./serial_result1.txt"
#define ResultFile2 "./serial_result2.txt"

// timer
double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

// y = A * x ( A in CSR format)
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
    int n = 0, nnz = 0;
    int *index = NULL, *col_id = NULL;
    double *val = NULL, *x = NULL, *y = NULL;

    read_size(MatrixFile1, &n, &nnz);
    nnz = nnz * 2 - n;
    
    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    col_id = (int *) malloc(nnz * sizeof(int));
    index = (int *) malloc((n + 1) * sizeof(int));
    val = (double *) malloc(nnz * sizeof(double));    

    read_matrix(MatrixFile1, n, nnz, index, col_id, val);

    gen_vector(VectorFile1, n, x);

    double start = get_time();
    matvec(n, index, col_id, val, x, y);  
    double end = get_time();

    printf("Time;\tNumber of proscesses;\tMatrix size;\tNon-zeros\n");
    printf("---------------------------------------------------------\n");
    printf("%.6f\t%d\t%d\t%d\n", end - start, 1, n, nnz);
    write_vector(ResultFile1, n, y);

    free(x);
    free(y);
    free(col_id);
    free(index);
    free(val);

    return 0;
}