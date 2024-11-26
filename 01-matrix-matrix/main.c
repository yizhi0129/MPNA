#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void matrix_matrix_multiplication(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    // Write the matrix-matrix multiplication here
}


void fill_matrix(int n, int m, double A[n][m])
{
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < m; k++)
        {
            A[i][k] = 1.0;
        }
    }
}


int main(int argc, char *argv[])
{
    int n = 100;

    // Parse dimensions of matrices
    if (argc > 1)
    {
        n = atoi(argv[1]);
    }

    int m = n;
    int p = n;

    if (argc > 2)
    {
        m = atoi(argv[2]);
    }

    if (argc > 3)
    {
        p = atoi(argv[3]);
    }

    printf("Working with matrices of size %d x %d and %d x %d\n", n, m, m, p);

    // Allocate memory for matrices
    // Note how we allocate memory for a 2D array in a single block
    // This is done to ensure that the memory is contiguous
    // Can you comment on how data are stored in memory for this 2D array?
    double **A = (double **)malloc(n * m * sizeof(double));
    double **B = (double **)malloc(m * p * sizeof(double));
    double **C = (double **)malloc(n * p * sizeof(double));

    // Initialize matrices
    fill_matrix(n, m, A);
    fill_matrix(m, p, B);
    fill_matrix(n, p, C);

    // Call matrix_matrix_multiplication
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    matrix_matrix_multiplication(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

    double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken: %e s\n", time_taken);

    // Free memory
    free(A);
    free(B);
    free(C);

    return EXIT_SUCCESS;
}