#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <HYPRE.h>
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"

#define MAX_ITER 1000
#define TOL 1e-6
#define KAPPA_0 0.01

double kappa1(double u)
{
    return KAPPA_0 * sqrt(u);
}

double kappa2(double u)
{
    return KAPPA_0 * u * u;
}

double time_step(double gamma, double sigma, int N, double u_max, double kappa_max)
{
    return  gamma / (2 * (sigma * u_max * u_max * u_max + kappa_max * N * N));
}

double maximun(double *u, int n)
{
    double max = 0;
    for (int i = 0; i < n; i ++)
    {
        if (u[i] > max)
        {
            max = u[i];
        }
    }
    return max;
}

double heaviside(int i, int N, double beta)
{
    return (i <= 0.2 * N) * beta;
}


void local_diff(double *u, double *u_new, int local_size, double *diff, double *sum)
{
    *diff = 0;
    *sum = 0;
    for (int i = 0; i < local_size + 1; i ++)
    {
        *sum += u_new[i] * u_new[i];
        *diff += (u_new[i] - u[i]) * (u_new[i] - u[i]);
    }
}


int main(int argc, char **argv)
{
    if (argc != 6)
    {
        printf("Please enter parameters: <kappa_type: 1 or 2> <sigma> <beta> <gamma> <N>\n");
        return 1;
    }
    double sigma = atof(argv[2]);
    double beta = atof(argv[3]);
    double gamma = atof(argv[4]);
    int N = atoi(argv[5]);
    double h = 1.0 / N;

    double *u = (double *)malloc((N + 1) * sizeof(double));
    double *u_new = (double *)malloc((N + 1) * sizeof(double));

    MPI_init();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_size = (N + 1) / size + (rank < (N + 1) % size);
    int ilower = rank * (N + 1) / size + (rank < (N + 1) % size);
    int iupper = ilower + local_size - 1;
    
    HYPRE_IJMatrix ij_A;
    HYPRE_IJVector ij_b, ij_u;
    HYPRE_ParCSRMatrix A;
    HYPRE_ParVector b, u;
}


