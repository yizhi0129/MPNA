#include "../include/solver.h"
#include "../include/matrix_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "HYPRE_utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"

void solve_newton_method(int N, double tol, int max_iter) 
{
    int myid, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    double L = 1.0;
    double dx = L / (N - 1);
    double K0 = KAPPA0, sigma = SIGMA, delta = 0.2;

    int ilower = myid * N / num_procs;
    int iupper = (myid + 1) * N / num_procs - 1;
    int local_n = iupper - ilower + 1;

    double *X = malloc(local_n * sizeof(double));
    double *Q = malloc(local_n * sizeof(double));
    double *u = calloc(local_n + 2, sizeof(double));  // include ghost cells
    double *Fn = calloc(local_n, sizeof(double));
    double *u_new = calloc(local_n, sizeof(double));

    for (int i = 0; i < local_n; i ++) 
    {
        int gi = ilower + i;
        double x = gi * dx;
        X[i] = x;
        Q[i] = heaviside(x);
        u[i+1] = 1.0;
    }

    for (int iter = 0; iter < max_iter; iter ++) 
    {
        // ghost cell communication
        double send_left = u[1], recv_left = 0.0;
        double send_right = u[local_n], recv_right = 0.0;
        MPI_Request reqs[4];
        if (myid > 0)
            MPI_Isend(&send_left, 1, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, &reqs[0]);
        if (myid < num_procs - 1)
            MPI_Isend(&send_right, 1, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, &reqs[1]);
        if (myid > 0)
            MPI_Irecv(&recv_left, 1, MPI_DOUBLE, myid - 1, 1, MPI_COMM_WORLD, &reqs[2]);
        if (myid < num_procs - 1)
            MPI_Irecv(&recv_right, 1, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, &reqs[3]);
        MPI_Waitall((myid == 0 || myid == num_procs - 1) ? 2 : 4, reqs, MPI_STATUSES_IGNORE);
        u[0] = (myid == 0) ? u[2] : recv_left;
        u[local_n + 1] = (myid == num_procs - 1) ? 1.0 : recv_right;

        double *Kn12 = malloc((local_n + 1) * sizeof(double));
        for (int i = 0; i < local_n + 1; i ++) 
        {
            Kn12[i] = 0.5 * (kappa(u[i]) + kappa(u[i + 1]));
        }

        HYPRE_IJMatrix J;
        HYPRE_IJVector F, delta;
        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &J);
        HYPRE_IJMatrixSetObjectType(J, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(J);

        HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &F);
        HYPRE_IJVectorSetObjectType(F, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(F);

        HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &delta);
        HYPRE_IJVectorSetObjectType(delta, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(delta);

        for (int i = 0; i < local_n; i ++) 
        {
            int gi = ilower + i;
            double Fi;
            int cols[3];
            double vals[3];
            int ncols = 0;

            if (gi == 0) 
            {
                Fi = K0 / (dx * dx) * (Kn12[0] * (u[2] - u[1]) + Kn12[0] * (u[2] - u[1]))
                     - sigma * (pow(u[1], 4) - 1.0) + Q[i];
                cols[0] = gi;
                cols[1] = gi + 1;
                vals[0] = K0 / (dx * dx) * (2 * Kn12[0]) + 4 * sigma * pow(u[1], 3);
                vals[1] = -2 * K0 * Kn12[0] / (dx * dx);
                ncols = 2;
            } 
            else if (gi == N - 1) 
            {
                Fi = u[local_n] - 1.0;
                cols[0] = gi;
                vals[0] = 1.0;
                ncols = 1;
            } 
            else 
            {
                Fi = K0 / (dx * dx) * (Kn12[i + 1] * (u[i + 2] - u[i + 1]) + Kn12[i] * (u[i] - u[i + 1]))
                     - sigma * (pow(u[i + 1], 4) - 1.0) + Q[i];
                cols[0] = gi - 1;
                cols[1] = gi;
                cols[2] = gi + 1;
                vals[0] = -K0 * Kn12[i] / (dx * dx);
                vals[1] = K0 / (dx * dx) * (Kn12[i + 1] + Kn12[i]) + 4 * sigma * pow(u[i + 1], 3);
                vals[2] = -K0 * Kn12[i + 1] / (dx * dx);
                ncols = 3;
            }

            HYPRE_IJMatrixSetValues(J, 1, &ncols, &gi, cols, vals);
            HYPRE_IJVectorSetValues(F, 1, &gi, &Fi);
        }

        HYPRE_IJMatrixAssemble(J);
        HYPRE_IJVectorAssemble(F);
        HYPRE_IJVectorAssemble(delta);

        solve_with_hypre(J, F, delta);

        double err2_local = 0.0;
        for (int i = 0; i < local_n; i ++) 
        {
            int gi = ilower + i;
            double d;
            HYPRE_IJVectorGetValues(delta, 1, &gi, &d);
            u[i+1] -= d;
            err2_local += d * d;
        }

        double err2_global;
        MPI_Allreduce(&err2_local, &err2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double err = sqrt(err2_global / N);

        if (myid == 0)
            printf("Newton iter %d: residual = %.2e\n", iter, err);

        HYPRE_IJMatrixDestroy(J);
        HYPRE_IJVectorDestroy(F);
        HYPRE_IJVectorDestroy(delta);
        free(Kn12);

        if (err < tol) break;
    }

    free(X); free(Q); free(u); free(Fn); free(u_new);
}
