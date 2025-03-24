#include <math.h>

#include "files.h"

//#define MatrixFile "./data/bcsstk03/bcsstk03.mtx"
#define MatrixFile "./data/cfd1/cfd1.mtx"

void circle(int n, int *index, int *col_id, double *val)
{
    double *center = (double *) malloc(n * sizeof(double));
    double *radius = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i ++)
    {
        double r = 0;
        for (int j = index[i]; j < index[i + 1]; j ++)
        {
            if (col_id[j] % n == i) // is diagonal
            {   
                center[i] = val[j];
                continue;
            }
            r += fabs(val[j]);
        }
        radius[i] = r;
        printf("%d\t%lf\t%lf\n", i, center[i], radius[i]);
    }
    free(center);
    free(radius);
}

int main(int argc, char** argv)
{
    int n = 0, nnz = 0;
    int *index = NULL, *col_id = NULL;
    double *val = NULL;

    read_size(MatrixFile, &n, &nnz);
    nnz = nnz * 2 - n;

    index = (int *) malloc((n + 1) * sizeof(int));
    col_id = (int *) malloc(nnz * sizeof(int));
    val = (double *) malloc(nnz * sizeof(double));
    read_matrix(MatrixFile, n, nnz, index, col_id, val);

    circle(n, index, col_id, val);

    free(index);
    free(col_id);
    free(val);

    return 0;
}