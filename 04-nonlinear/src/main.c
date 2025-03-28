#include "../include/solver.h"
#include "../include/matrix_utils.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char **argv) 
{
    MPI_Init(&argc, &argv);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int N = 100;
    int steps = 20000;
    double tol = 1e-6;
    int max_iter = 20000;
    int use_newton = 0;

    for (int i = 1; i < argc; i ++) 
    {
        if (strcmp(argv[i], "--N") == 0) N = atoi(argv[++ i]);
        else if (strcmp(argv[i], "--steps") == 0) steps = atoi(argv[++ i]);
        else if (strcmp(argv[i], "--tol") == 0) tol = atof(argv[++ i]);
        else if (strcmp(argv[i], "--max_iter") == 0) max_iter = atoi(argv[++ i]);
        else if (strcmp(argv[i], "--newton") == 0) use_newton = 1;
    }

    init_hypre();

    if (use_newton) 
    {
        solve_newton_method(N, tol, max_iter);
    } 
    else 
    {
        solve_linearized_implicit(N, steps);
    }

    finalize_hypre();

    MPI_Finalize();
    return 0;
}
