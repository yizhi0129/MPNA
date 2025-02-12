#ifndef TOOLS_H
#define TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define EPSILON 1e-7 // convergence criterion
#define MAX_ITER 1000 // maximum number of iterations

// store the matrix in CSR format
void generate_grid_CSR(int n, int * indices, int * col_id, double * values)
{
    int count = 0; // total number of non zeros to be stored
    indices[0] = 0;

    for (int j = 0; j < n; j ++)
    {
        for (int i = 0; i < n; i ++)
        {
            int index = i + j * n; // global index
            if (j > 0)
            {
                col_id[count] = index - n; // neighbour below
                values[count] = -1.0;
                count ++;
            }
            if (i > 0)
            {
                col_id[count] = index - 1; // neighbour to the left
                values[count] = -1.0;
                count ++;
            }
            col_id[count] = i + j * n; // diagonal
            values[count] = 4.0;
            count ++;
            if (i < n - 1)
            {
                col_id[count] = index + 1; // neighbour to the right
                values[count] = -1.0;
                count ++;
            }
            if (j < n - 1)
            {
                col_id[count] = index + n; // neighbour above
                values[count] = -1.0;
                count ++;
            }
            indices[index + 1] = count;         
        }               
    }
}

// Jacobi method
void jacobi(int N, int * indices, int * col_id, double * values, double * b, double * x, double * x_new, double b_abs, double r_rel)
{
    // initialization
    for(int i = 0; i < N; i ++)
    {
        x[i] = 0.0;
        x_new[i] = 0.0;
    }

    int iter = 0;
    while (r_rel > EPSILON && iter < MAX_ITER)
    { 
        for (int i = 0; i < N; i ++)
        {
            int self = -1;
            double sum = 0.0;
            for (int j = indices[i]; j < indices[i + 1]; j ++)
            {
                if (col_id[j] != i)
                {
                    sum += values[j] * x[col_id[j]];
                }
                else
                {
                    self = j; // find the diagonal element
                }
            }  
            if (self != -1)
            {
                x_new[i] = (b[i] - sum) / values[self]; 
                x[i] = x_new[i]; // update x 
            } 
        }

        // calculate the residual ||Ax - b||
        double res_sum2 = 0.0;
        for (int i = 0; i < N; i ++)
        {
            double ax_i = 0.0;
            for (int j = indices[i]; j < indices[i + 1]; j ++)
            {
                ax_i += x_new[col_id[j]] * values[j];
            }
            res_sum2 += (ax_i - b[i]) * (ax_i - b[i]);
        }
        r_rel = sqrt(res_sum2) / b_abs; // relative residual
        iter ++; // iteration count
        printf("%d %.9f\n", iter, r_rel); // write the interation count and relative residual
    }
}

// Gauss-Seidel method
void gauss_seidel(int N, int * indices, int * col_id, double * values, double * b, double * x, double b_abs, double r_rel)
{
    // initialization
    for(int i = 0; i < N; i ++)
    {
        x[i] = 0.0;
    }

    int iter = 0;
    while (r_rel > EPSILON && iter < MAX_ITER)
    {    
        for (int i = 0; i < N; i ++)
        {
            double sum = 0.0;
            int self = -1;
            for (int j = indices[i]; j < indices[i + 1]; j ++)
            {
                if (col_id[j] == i)
                {
                    self = j; // find the diagonal element
                }
            }
            // split the sum into two parts: do not invert the order
            for (int j = indices[i]; j < self; j ++)
            {
                sum += values[j] * x[col_id[j]]; 
            }
            for (int j = self + 1; j < indices[i + 1]; j ++)
            {
                sum += values[j] * x[col_id[j]];
            }
            if (self != -1)
            {
                x[i] = (b[i] - sum) / values[self]; 
            }    
        }

        // calculate the residual ||Ax - b||
        double res_sum2 = 0.0;
        for (int i = 0; i < N; i ++)
        {
            double ax_i = 0.0;
            for (int j = indices[i]; j < indices[i + 1]; j ++)
            {
                ax_i += x[col_id[j]] * values[j];
            }
            res_sum2 += (ax_i - b[i]) * (ax_i - b[i]);     
        }
        r_rel = sqrt(res_sum2) / b_abs;
        iter ++;
        printf("%d %.9f\n", iter, r_rel); // write the interation count and relative residual
    }
}

// Conjugate Gradient method
void CG(int n, int * indices, int * col_id, double * values, double * b, double * x, int max_iter, double tol)
{
    
    
}

// Generalized Minimal Residual method
void GMRES(int n, int * indices, int * col_id, double * values, double * b, double * x, int max_iter, double tol)
{
    
}

#endif