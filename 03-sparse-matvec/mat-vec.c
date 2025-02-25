#include <sys/time.h>
#include <time.h>

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
    int n1, nnz1, n2, nnz2;
    int *index1, *col_id1, *index2, *col_id2;
    double *val1, *val2, *x1, *y1, *x2, *y2;

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

    printf("matrix1:");
    for (int i = 0; i < n1 + 1; i ++)
    {
        printf("%d ", index1[i]);
    }
    printf("\n");
    for (int i = 0; i < nnz1; i ++)
    {  
        printf("%d %d %lf\n", i, col_id1[i], val1[i]);        
    }
    printf("------------------------------\n");

    printf("matrix2:");
    for (int i = 0; i < n2 + 1; i ++)
    {
        printf("%d ", index2[i]);
    }
    printf("\n");
    for (int i = 0; i < nnz2; i ++)
    {  
        printf("%d %lf\n", col_id2[i], val2[i]);        
    }
    printf("------------------------------\n");

    gen_vector("vector1.txt", n1, x1);
    gen_vector("vector2.txt", n2, x2);
    
    double start1 = get_time();
    matvec(n1, index1, col_id1, val1, x1, y1);  
    double end1 = get_time();
    printf("Time;\tNumber of proscesses;\tMatrix size;\tNon-zeros\n");
    printf("---------------------------------------------------------\n");
    printf("%.6f\t%d\t%d\t%d\n", end1 - start1, 1, n1, nnz1);

    double start2 = get_time();
    matvec(n2, index2, col_id2, val2, x2, y2);
    double end2 = get_time();
    printf("%.6f\t%d\t%d\t%d\n", end2 - start2, 1, n2, nnz2);

    write_vector("serial_result1.txt", n1, y1);
    write_vector("serial_result2.txt", n2, y2);

    free(x1);
    free(y1);
    free(col_id1);
    free(index1);
    free(val1);

    free(x2);
    free(y2);
    free(col_id2);
    free(index2);
    free(val2);

    return 0;
}