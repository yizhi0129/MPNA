#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void CG(double * A, double * b, double * x, int n, int max_iter, double tol)
{
    double * r = (double *) malloc(n * sizeof(double));
    double * p = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i ++)
    {        
        p[i] = x[i];
        double sum = 0.0;
        for (int j = 0; j < n; j ++)
        {
            sum += A[i * n + j] * x[j];
        }
        r[i] = b[i] - sum;
    }
    double alpha = 0.0;
    double beta = 0.0;
    for (int k = 0; k < max_iter; k ++)
    {
        double r_dot_r = 0.0;
        for (int i = 0; i < n; i ++)
        {
            r_dot_r += r[i] * r[i];
        }
        double r_dot_r_new = 0.0;
        for (int i = 0; i < n; i ++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j ++)
            {
                sum += A[i * n + j] * p[j];
            }
            r_dot_r_new += r[i] * sum;
        }
        alpha = r_dot_r / r_dot_r_new;
        for (int i = 0; i < n; i ++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * A[i * n + i];
        }
        double r_dot_r_new = 0.0;
        for (int i = 0; i < n; i ++)
        {
            r_dot_r_new += r[i] * r[i];
        }
        if (r_dot_r_new < tol)
        {
            break;
        }
        beta = r_dot_r_new / r_dot_r;
        for (int i = 0; i < n; i ++)
        {
            p[i] = r[i] + beta * p[i];
        }
    }
}