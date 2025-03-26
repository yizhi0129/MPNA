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

void solve_linearized_implicit(int N, int max_steps) 
{
    int myid, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    double gamma = GAMMA;
    double L = 1.0;
    double dx = L / (N - 1);
    double sigma = SIGMA;

    int ilower = myid * N / num_procs;
    int iupper = (myid + 1) * N / num_procs - 1;
    int local_n = iupper - ilower + 1;

    double *X = calloc(local_n, sizeof(double));
    double *Q = calloc(local_n, sizeof(double));
    double *u = calloc(local_n + 2, sizeof(double));  // ghost cells

    for (int i = 0; i < local_n; i ++) 
    {
        int gi = ilower + i;
        X[i] = gi * dx;
        Q[i] = heaviside(X[i]);
        u[i + 1] = 1.0;
    }

    for (int step = 0; step < max_steps; step ++) 
    {
        // ghost communication
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

        u[0] = (myid != 0) * recv_l;
        u[local_n + 1] = (myid != num_procs - 1) * recv_r;

        // compute dt
        double umax_local = 0.0;
        for (int i = 1; i <= local_n; i ++)
            if (u[i] > umax_local) umax_local = u[i];
        double umax;
        MPI_Allreduce(&umax_local, &umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double dt = gamma * 2.0 / (4.0 * sigma * pow(umax, 3) + 4.0 * kappa(umax) / (dx * dx));

        // Kn+1/2
        double *Kn12 = calloc(local_n + 1, sizeof(double));
        for (int i = 0; i < local_n + 1; i ++)
            Kn12[i] = 0.5 * (kappa(u[i + 1]) + kappa(u[i + 2]));

        // HYPRE setup
        HYPRE_IJMatrix A;
        HYPRE_IJVector b, x;
        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
        HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A);

        HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
        HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b);

        HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
        HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(x);

        for (int i = 0; i < local_n; i ++) 
        {
            int gi = ilower + i;
            double rhs;
            int cols[3];
            double vals[3];
            int ncols = 0;

            if (gi == 0) 
            {
                rhs = u[1] + dt * (Q[0] + sigma);
                cols[0] = 0; 
                cols[1] = 1;
                vals[0] = 1 + dt * (2 * Kn12[0] / (dx * dx) + 4 * sigma * pow(u[1], 3));
                vals[1] = -2 * dt * Kn12[0] / (dx * dx);
                ncols = 2;
            } 
            else if (gi == N - 1) 
            {
                rhs = 1.0;
                cols[0] = gi; 
                vals[0] = 1.0;
                ncols = 1;
            } 
            else 
            {
                rhs = u[i + 1] + dt * (Q[i] + sigma);
                cols[0] = gi - 1; 
                cols[1] = gi; 
                cols[2] = gi + 1;
                vals[0] = - dt * Kn12[i - 1] / (dx * dx);
                vals[1] = 1 + dt * ((Kn12[i - 1] + Kn12[i]) / (dx * dx) + 4 * sigma * pow(u[i + 1], 3));
                vals[2] = - dt * Kn12[i] / (dx * dx);
                ncols = 3;
            }

            HYPRE_IJMatrixSetValues(A, 1, &ncols, &gi, cols, vals);
            HYPRE_IJVectorSetValues(b, 1, &gi, &rhs);
        }

        HYPRE_IJMatrixAssemble(A);
        HYPRE_IJVectorAssemble(b);
        HYPRE_IJVectorAssemble(x);

        solve_with_hypre(A, b, x);

        double err2_local = 0.0;
        for (int i = 0; i < local_n; i ++) 
        {
            int gi = ilower + i;
            double val;
            HYPRE_IJVectorGetValues(x, 1, &gi, &val);
            double diff = (val - u[i + 1]) / dt;
            err2_local += diff * diff;
            u[i + 1] = val;
        }

        double err2_global;
        MPI_Allreduce(&err2_local, &err2_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double err = sqrt(err2_global / N);
        if (myid == 0 && step % 10 == 0)
            printf("Step %d: dt = %.2e, error = %.2e\n", step, dt, err);

        free(Kn12);
        HYPRE_IJMatrixDestroy(A);
        HYPRE_IJVectorDestroy(b);
        HYPRE_IJVectorDestroy(x);

        if (err < 1e-5) break;
    }

    free(X); 
    free(Q); 
    free(u); 
}
