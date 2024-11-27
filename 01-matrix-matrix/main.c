#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//#include </usr/local/Cellar/openblas/0.3.26/include/cblas.h>

// naive version ijk
void matrix_matrix_multiplication(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for (int i = 0; i < n; i ++)
    {
        for (int j = 0; j < p; j ++)
        {
            for (int k = 0; k < m; k ++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// switched version kij
void matrix_matrix_multiplication2(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for (int k = 0; k < m; k ++)
    {
        for (int i = 0; i < n; i ++)
        {
            for (int j = 0; j < p; j ++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// switched version jki
void matrix_matrix_multiplication3(int n, int m, int p, double A[n][m], double B[m][p], double C[n][p])
{
    for (int j = 0; j < p; j ++)
    {
        for (int k = 0; k < m; k ++)
        {
            for (int i = 0; i < n; i ++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// square matrix version (B transposed)
void matrix_matrix_multiplication4(int n, double A[n][n], double B[n][n], double C[n][n])
{
    for (int j = 0; j < n; j ++)
    {
        for (int i = 0; i < n; i ++)
        {
            for (int k = 0; k < n; k ++)
            {
                C[i][j] += A[i][k] * B[j][k];
            }
        }
    }
}

// block version

// BLAS version
//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, p, m, 1.0, A, m, B, p, 0.0, C, p);

void fill_matrix(int n, int m, double A[n][m])
{
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < m; k++)
        {
            A[i][k] = (double)rand() / RAND_MAX;
        }
    }
}

void init_matrix(int n, int m, double A[n][m])
{
    for (int i = 0; i < n; i ++)
    {
        for (int k = 0; k < m; k++)
        {
            A[i][k] = 0.0;
        }
    }
}

int main(int argc, char *argv[])
{
    char *mode = argv[1];
    int n = 100;

    // Parse dimensions of matrices
    if (argc > 2)
    {
        n = atoi(argv[2]);
    }

    int m = n;
    int p = n;

    if (argc > 3)
    {
        m = atoi(argv[3]);
    }

    if (argc > 4)
    {
        p = atoi(argv[4]);
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
    init_matrix(n, p, C);

    if (strcmp(mode, "naive") == 0)
    {
        printf("Running naive version\n");
        
        // Call matrix_matrix_multiplication
        struct timespec start, end;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        matrix_matrix_multiplication(n, m, p, A, B, C);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        
        double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
        printf("Time taken: %e s\n", time_taken);
    }
    else if (strcmp(mode, "switched") == 0)
    {
        printf("Running switched version\n");
        // Call matrix_matrix_multiplication
        struct timespec start, end;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        matrix_matrix_multiplication2(n, m, p, A, B, C);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        
        double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
        printf("Time taken: %e s\n", time_taken);
    }
    else if (strcmp(mode, "another_switched") == 0)
    {
        printf("Running another switched version\n");
        // Call matrix_matrix_multiplication
        struct timespec start, end;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        matrix_matrix_multiplication3(n, m, p, A, B, C);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        
        double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
        printf("Time taken: %e s\n", time_taken);
    }
    else if (strcmp(mode, "square") == 0)
    {
        printf("Running square version\n");
        // Call matrix_matrix_multiplication
        struct timespec start, end;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        matrix_matrix_multiplication4(n, A, B, C);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        
        double time_taken = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1.e9;
        printf("Time taken: %e s\n", time_taken);
    }
    else
    {
        printf("Please indicate calculate version: naive, switched, another_switched or square\n");
    }

    // Free memory
    free(A);
    free(B);
    free(C);

    return EXIT_SUCCESS;
}