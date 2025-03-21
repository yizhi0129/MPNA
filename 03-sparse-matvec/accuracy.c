#include "files.h"

#define MatrixFile "./data/bcsstk03/bcsstk03.mtx"
#define VectorFileRef "./serial_result1.txt"
#define ResultFile "./parallel_result1"

//#define MatrixFile "./data/cfd1/cfd1.mtx"
//#define VectorFileRef "./serial_result2.txt"
//#define ResultFile "./parallel_result2"

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <number of processes>\n", argv[0]);
        exit(1);
    }

    int n = 0, nnz = 0;
    int n_proc = atoi(argv[1]);

    double *xs = NULL, *xp = NULL, *diff = NULL;
    read_size(MatrixFile, &n, &nnz);
    
    xs = (double *) malloc(n * sizeof(double));
    xp = (double *) malloc(n * sizeof(double));
    diff = (double *) malloc(n * sizeof(double));
    
    read_vector(VectorFileRef, n, xs);
    double ref = dot_vec(n, xs, xs);

    char filename[50];
    sprintf(filename, "%s_%d.txt", ResultFile, n_proc);    
    read_vector(filename, n, xp);

    diff_vec(n, xs, xp, diff);
    double error = dot_vec(n, diff, diff);
    error /= ref;
    error = sqrt(error);
    printf("%d\t%.10f\n", n_proc, error);

    free(xs);
    free(xp);
    free(diff);
    
    return 0;
}
