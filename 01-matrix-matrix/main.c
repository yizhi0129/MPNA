#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void matrix_matrix_multiplication(int n, int m, int, int p, double A[n][m], double B[m][p], double C[n][p])
{
    // Write the matrix-matrix multiplication here
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < p; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < m; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


int main(int argc, char *argv[])
{
    int n = 100;
    int m = 100;
    int p = 100;

    // Parse dimensions of matrices
    if (argc > 1)
    {
        n = atoi(argv[1]);
    }
    if (argc > 2)
    {
        m = atoi(argv[2]);
    }
    if (argc > 3)
    {
        p = atoi(argv[3]);
    }

    // Allocate memory for matrices
    double **A;
    double **B;
    double **C;

    // Initialize matrices


    // Call matrix_matrix_multiplication
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    // matrix_matrix_multiplication(n, m, p, A, B, C);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

    double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
    printf("Time taken: %e s\n", time_taken);

    return EXIT_SUCCESS;
}