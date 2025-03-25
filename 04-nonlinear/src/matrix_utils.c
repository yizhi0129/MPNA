#include "../include/matrix_utils.h"
#include <stdio.h>


void init_hypre() 
{
    HYPRE_Init();
}

void finalize_hypre() 
{
    HYPRE_Finalize();
}

void solve_with_hypre(HYPRE_IJMatrix A_ij, HYPRE_IJVector b_ij, HYPRE_IJVector x_ij) 
{
    HYPRE_ParCSRMatrix A;
    HYPRE_ParVector b, x;

    // Get the underlying ParCSR objects
    HYPRE_IJMatrixGetObject(A_ij, (void **) &A);
    HYPRE_IJVectorGetObject(b_ij, (void **) &b);
    HYPRE_IJVectorGetObject(x_ij, (void **) &x);

    // Create solver and preconditioner
    HYPRE_Solver solver, precond;

    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_ParCSRGMRESSetMaxIter(solver, 1000);
    HYPRE_ParCSRGMRESSetTol(solver, 1e-8);
    HYPRE_ParCSRGMRESSetLogging(solver, 1);

    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, 1);      // Print AMG info
    HYPRE_BoomerAMGSetNumSweeps(precond, 1);       // Sweeps per level
    HYPRE_BoomerAMGSetMaxLevels(precond, 25);
    HYPRE_BoomerAMGSetTol(precond, 0.0);           // AMG as preconditioner

    HYPRE_ParCSRGMRESSetPrecond(solver,
                                HYPRE_BoomerAMGSolve,
                                HYPRE_BoomerAMGSetup,
                                precond);

    // Setup and solve
    HYPRE_ParCSRGMRESSetup(solver, A, b, x);
    HYPRE_ParCSRGMRESSolve(solver, A, b, x);

    // Optional: print iterations info
    int num_iterations;
    double final_res_norm;
    HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iterations);
    HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0) 
    {
        printf("GMRES iterations = %d\n", num_iterations);
        printf("Final relative residual norm = %e\n", final_res_norm);
    }

    // Clean up
    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
}

