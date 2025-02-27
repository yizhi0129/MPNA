#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "files.h"

#define MatrixFile "./data/bcsstk03/bcsstk03.mtx"
//#define MatrixFile "./data/cfd1/cfd1.mtx"

#define VectorFile "./vector1.txt"
//#define VectorFile "./vector2.txt"

#define IterFile "./iter1"
//#define IterFile "./iter2"

#define EigFile "./eig1"
//#define EigFile "./eig2"

#define EPSILON 1e-6

// timer
double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

double eigenvalue(int n, int nnz, int *index, int *col_id, double *val, double *x)
{
    double lambda_max = 0.0;
    double r = 1.0;
    int count = 0;
    double *x_new = (double *) malloc(n * sizeof(double));
    double *diff = (double *) malloc(n * sizeof(double));

    // x = x / ||x||
    normalize(n, x); 

    while (r > EPSILON)
    {
        // x_new = A * x
        matvec(n, index, col_id, val, x, x_new); 

        // x_new = x_new / ||x_new||
        normalize(n, x_new); 
        diff_vec(n, x_new, x, diff);

        // r = ||x_new - x||
        r = dot_vec(n, diff, diff);
        r = sqrt(r);      
        
        // x = x_new
        eq_vec(n, x_new, x);                                            
        count ++;
        printf("Iteration %d, residual = %.10f\n", count, r);
    }

    //  x = A * x_new
    matvec(n, index, col_id, val, x_new, x);                     
    
    // lambda_max = x^T * A * x / x^T * x
    lambda_max = dot_vec(n, x, x_new) / dot_vec(n, x_new, x_new); 

    free(x_new);
    free(diff);
    return lambda_max;
}

int main(int argc, char** argv)
{
    int n = 0, nnz = 0;
    int *index = NULL, *col_id = NULL;
    double *val = NULL, *x = NULL;
    double lambda_max = 0.0;

    read_size(MatrixFile, &n, &nnz);
    nnz = nnz * 2 - n;
    
    x = (double *) malloc(n * sizeof(double));
    col_id = (int *) malloc(nnz * sizeof(int));
    index = (int *) malloc((n + 1) * sizeof(int));
    val = (double *) malloc(nnz * sizeof(double));    

    read_matrix(MatrixFile, n, nnz, index, col_id, val);

    read_vector(VectorFile, n, x);

    double start = get_time();
    lambda_max = eigenvalue(n, nnz, index, col_id, val, x);  
    double end = get_time();

    printf("Max eigenvalue = %.10f\n", lambda_max);

    printf("Time;\tNumber of proscesses;\tMatrix size;\tNon-zeros\n");
    printf("%.10f\t%d\t%d\t%d\n", end - start, 1, n, nnz);

    free(x);
    free(col_id);
    free(index);
    free(val);

    return 0;
}