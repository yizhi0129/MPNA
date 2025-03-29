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

    int ilower = myid * N / num_procs;
    int iupper = (myid + 1) * N / num_procs - 1;
    int local_n = iupper - ilower + 1;

    double *X = malloc(local_n * sizeof(double));
    double *Q = malloc(local_n * sizeof(double));
    double *u = calloc(local_n + 2, sizeof(double));  // include ghost cells
    double *Fn = calloc(local_n, sizeof(double));

    for (int i = 0; i < local_n; i ++) 
    {
        int gi = ilower + i;
        double x = gi * dx;
        X[i] = x;
        Q[i] = heaviside(x);
        u[i + 1] = 1.0 + 4.0 * SIGMA * Q[i];
    }

    for (int iter = 0; iter < max_iter; iter ++) 
    {
        // ghost cell communication
        double send_l = u[1], recv_l = 0.0;
        double send_r = u[local_n], recv_r = 0.0;
        MPI_Request reqs[4];
        int req_idx = 0;
        if (myid > 0)
            MPI_Isend(&send_l, 1, MPI_DOUBLE, myid - 1, 0, MPI_COMM_WORLD, &reqs[req_idx ++]);
        if (myid < num_procs - 1)
            MPI_Isend(&send_r, 1, MPI_DOUBLE, myid + 1, 1, MPI_COMM_WORLD, &reqs[req_idx ++]);
        if (myid > 0)
            MPI_Irecv(&recv_l, 1, MPI_DOUBLE, myid - 1, 1, MPI_COMM_WORLD, &reqs[req_idx ++]);
        if (myid < num_procs - 1)
            MPI_Irecv(&recv_r, 1, MPI_DOUBLE, myid + 1, 0, MPI_COMM_WORLD, &reqs[req_idx++]);
        MPI_Waitall(req_idx, reqs, MPI_STATUSES_IGNORE);

        u[0] = (myid != 0) ? recv_l : u[1];
        u[local_n + 1] = (myid != num_procs - 1) ? recv_r : 1.0;

        double *Kn12 = malloc((local_n + 1) * sizeof(double));
        for (int i = 0; i < local_n + 1; i ++) 
        {
            double ku1 = kappa(u[i]);
            double ku2 = kappa(u[i + 1]);
            if (isnan(ku1) || isnan(ku2) || ku1 < 0 || ku2 < 0) 
            {
                ku1 = fmax(ku1, 0.0);
                ku2 = fmax(ku2, 0.0);
            }
            Kn12[i] = 0.5 * (ku1 + ku2);
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

        for (int i = 0; i < local_n; i++) 
        {
            int gi = ilower + i;
            double Fi = 0.0;
            int cols[3];
            double vals[3];
            int ncols = 0;

            double ui = u[i + 1];
            double uim1 = u[i];
            double uip1 = u[i + 2];

            if (isnan(ui)) ui = 1e-6;

            if (gi == N - 1) 
            {
                Fi = ui - 1.0;
                cols[0] = gi;
                vals[0] = 1.0;
                ncols = 1;
            } 
            else if (gi == 0) 
            {
                Fi = (2 * Kn12[0] * (uip1 - ui)) / (dx * dx) - SIGMA * (pow(ui, 4.0) - 1.0) + Q[i];
                cols[0] = gi;
                cols[1] = gi + 1;
                vals[0] = 2 * Kn12[0] / (dx * dx) + 4 * SIGMA * pow(ui, 3);
                vals[1] = - 2 * Kn12[0] / (dx * dx);
                ncols = 2;
            } 
            else 
            {
                Fi = (Kn12[i - 1] * (uim1 - ui) + Kn12[i] * (uip1 - ui)) / (dx * dx) - SIGMA * (pow(ui, 4.0) - 1.0) + Q[i];

                cols[0] = gi - 1;
                cols[1] = gi;
                cols[2] = gi + 1;
                vals[0] = - Kn12[i - 1] / (dx * dx);
                vals[1] = (Kn12[i - 1] + Kn12[i]) / (dx * dx) + 4 * SIGMA * pow(ui, 3);
                vals[2] = - Kn12[i] / (dx * dx);
                ncols = 3;
            }

            if (isnan(Fi) || isinf(Fi)) 
            {
                printf("Rank %d: bad Fi at gi = %d, u = %.3e, Kn12 = %.3e\n", myid, gi, ui, Kn12[i]);
                MPI_Abort(MPI_COMM_WORLD, 1);
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
            u[i + 1] += d;
            err2_local += d * d;
        }

        double err2_global;
        MPI_Allreduce(&err2_local, &err2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double err = sqrt(err2_global / N);

        if (myid == 0)
            printf("%d\t%.2e\n", iter, err);

        HYPRE_IJMatrixDestroy(J);
        HYPRE_IJVectorDestroy(F);
        HYPRE_IJVectorDestroy(delta);
        free(Kn12);

        if (err < tol) break;
    }

    free(X); 
    free(Q); 
    free(u); 
    free(Fn); 
}
