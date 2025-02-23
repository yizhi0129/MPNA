#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>
#include <time.h>

#define EPSILON 1e-7 // convergence criterion
#define MAX_ITER 10000 // maximum number of iterations

// timer
double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}


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

// initialize the vector
void init(int N, double * x)
{
    for (int i = 0; i < N; i ++)
    {
        x[i] = 0.0;
    }
}

// Jacobi method
void jacobi(FILE *Jac, int N, int * indices, int * col_id, double * values, double * b, double * x, double * x_new, double b_abs, double r_rel)
{
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
        fprintf(Jac, "%d %.9f\n", iter, r_rel); // write the interation count and relative residual
    }
}

// Gauss-Seidel method
void gauss_seidel(FILE *GS, int N, int * indices, int * col_id, double * values, double * b, double * x, double b_abs, double r_rel)
{
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
        fprintf(GS, "%d %.9f\n", iter, r_rel); // write the interation count and relative residual
    }
}

// matrix-vector multiplication
void mat_vec_SCR(int N, int * indices, int * col_id, double * values, double * x, double * y)
{
    for (int i = 0; i < N; i ++)
    {
        y[i] = 0.0;
        for (int j = indices[i]; j < indices[i + 1]; j ++)
        {
            y[i] += values[j] * x[col_id[j]];
        }
    }
}

// Conjugate Gradient method
void CG_CSR(FILE *CG, int N, int * indices, int * col_id, double * values, double * b, double * x, double b_abs, double r_rel)
{
    double d[N], r[N], Ad[N];
    for (int i = 0; i < N; i ++)
    {
        r[i] = b[i];
        d[i] = r[i];
    }
    int count = 0;
    double alpha = 0.0, beta = 0.0, r_d = 0.0, r_Ad = 0.0, d_Ad = 0.0;
    while (count < 1000 && r_rel > EPSILON)
    {
        mat_vec_SCR(N, indices, col_id, values, d, Ad);
        for (int i = 0; i < N; i ++)
        {
            r_d += r[i] * d[i];
            r_Ad += r[i] * Ad[i];
            d_Ad += d[i] * Ad[i];
        }
        alpha = r_d / d_Ad;
        for (int i = 0; i < N; i ++)
        {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ad[i];
        }
        r_rel = sqrt(r_d) / b_abs;
        beta = r_d / r_Ad;
        for (int i = 0; i < N; i ++)
        {
            d[i] = r[i] + beta * d[i];
        }
        count ++;
        fprintf(CG, "%d %.9f\n", count, r_rel); // write the interation count and relative residual
    }
}

// norm of a vector
double norm(int N, double * x)
{
    double sum = 0.0;
    for (int i = 0; i < N; i ++)
    {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}

// dot product
double dot_prod(int N, double *x, double *y)
{
    double result = 0.0;
    for (int i = 0; i < N; i ++)
    {
        result += x[i] * y[i];
    }
    return result;
}

// CSR storage
void CSR(int N, double *V, int *indices, int *col_id, double *values)
{
    indices[0] = 0;
    int count = 0;
    for (int i = 0; i < N; i ++)
    {
        for (int j = 0; j < N; j ++)
        {
            if (V[N * i + j] != 0)
            {
                col_id[count] = N * i + j;
                values[count] = V[N * i + j];
                count ++;
            }
        }
        indices[i + 1] = count;
    }
}

// initialization r0
void init_r0(int N, double *b, double b_abs, double *r)
{
    for (int i = 0; i < N; i ++)
    {
        r[i] = b[i] / b_abs;
    }
}

// GMRES method
void GMRES_CSR(FILE *GMRES, int N, int * indices, int * col_id, double * values, double * b, double * x, double b_abs, double * r)
{
    double V[MAX_ITER + 1][N], H[MAX_ITER][MAX_ITER] = {0};
    double g[MAX_ITER + 1] = {0}, y[MAX_ITER] = {0};
    
    g[0] = b_abs;

    // Arnoldi
    for (int k = 0; k < MAX_ITER; k ++)
    {
        double w[N];

        // w = A * V[k]
        mat_vec_SCR(N, indices, col_id, values, V[k], w);

        // orthogonalization
        for (int j = 0; j <= k; j++) 
        {
            H[j][k] = dot_prod(N, w, V[j]);
            for (int i = 0; i < N; i++)
            {
                w[i] -= H[j][k] * V[j][i];
            }
        }
        H[k + 1][k] = norm(N, w);
        
        if (H[k + 1][k] < EPSILON)
        {
            break;
        }

        // normalization
        for (int i = 0; i < N; i++) 
        {
            V[k+1][i] = w[i] / H[k+1][k];
        }

        // least square
        for (int i = 0; i <= k; i++) 
        {
            double temp = H[i][k];
            H[i][k] = H[k][k] * temp - H[i][k] * H[k + 1][k];
            H[k + 1][k] = H[k + 1][k] * temp + H[i][k] * H[k + 1][k];
        }

        g[k + 1] = -g[k] * H[k + 1][k] / H[k][k];
        g[k] = g[k] / H[k][k];

        if (fabs(g[k + 1]) < EPSILON) 
        {
            break;
        } 
    }

    // compute y
    for (int i = MAX_ITER - 1; i >= 0; i --) 
    {
        y[i] = g[i];
        for (int j = i + 1; j < MAX_ITER; j ++)
        {
            y[i] -= H[i][j] * y[j];
        }         
        y[i] /= H[i][i];
    }

    // compute x = x0 + V * y
    for (int i = 0; i < N; i ++) 
    {
        for (int j = 0; j < MAX_ITER; j ++) 
        {
            x[i] += V[j][i] * y[j];
        }
    }
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        printf("Please enter a number n: grid size");
        return 1;
    }

    // grid size
    int n = atoi(argv[1]);
    int N = n * n;

    // 3 neighbours per corner, 4 neighbours for edges, 5 neighbours for interiors
    int nonzeros = 3 * 4 + 4 * 4 * (n - 2) + 5 * (n - 2) * (n - 2); 
    
    // allocate memory
    int * indices = (int *)malloc((N + 1) * sizeof(int));
    int * col_id = (int *)malloc(nonzeros * sizeof(int));
    double * values = (double *)malloc(nonzeros * sizeof(double));

    double * b = (double *)malloc(N * sizeof(double));
    double b_abs2 = 0.0, r_rel = 1.0; // relative residual initialized as 1.0 (could be any value greater than EPSILON)
    double * r = (double *)malloc(N * sizeof(double));
    double * x = (double *)malloc(N * sizeof(double));
    double * x_new = (double *)malloc(N * sizeof(double));

    // initialization 
    // printf("b: ");
    for(int i = 0; i < N; i ++)
    {
        b[i] = (double)rand() / RAND_MAX + 1;  // random values [1, 2] 
        // printf("%.3f ", b[i]);
        b_abs2 += b[i] * b[i];
    }
    // printf("\n");
    double b_abs = sqrt(b_abs2); // pre-calculate ||b||
    printf("b_abs: %.3f\n", b_abs);

    generate_grid_CSR(n, indices, col_id, values);
    //printf("Indices: ");
/*    for (int i = 0; i <= N; i ++) 
    {
        printf("%d ", indices[i]);
    }*/
    //printf("\n");
    //printf("Col_id: ");
/*    for (int i = 0; i < indices[N]; i ++) 
    {
        printf("%d ", col_id[i]);
    }*/
    //printf("\n");
    //printf("Values: ");
/*    for (int i = 0; i < indices[N]; i ++) 
    {
        printf("%.3f ", values[i]);
    } */
    //printf("\n");

    char file1[50], file2[50], file3[50], file4[50];
    sprintf(file1, "jacobi_%d.txt", n);
    sprintf(file2, "gauss_seidel_%d.txt", n);
    sprintf(file3, "conjugate_gradient_%d.txt", n);
    sprintf(file4, "gmres_%d.txt", n);

    FILE *Jac = fopen(file1, "w");
    if (!Jac) 
    {
        perror("Error opening file: jacobi_n.txt");
        return 1;
    }

    init(N, x);
    init(N, x_new);

    // compute Jacobi
    double start_jac = get_time();
    jacobi(Jac, N, indices, col_id, values, b, x, x_new, b_abs, r_rel);
    double end_jac = get_time();

    //printf("Jacobi: ");
/*    for (int i = 0; i < N; i ++)
    {
        printf("%.3f ", x[i]);
    } */
    //printf("\n");
    printf("Jacobi runtime: %.6f seconds\n", end_jac - start_jac);
    fclose(Jac);

    FILE *GS = fopen(file2, "w");
    if (!GS) 
    {
        perror("Error opening file: gauss_seidel_n.txt");
        return 1;
    }

    init(N, x);

    // compute Gauss-Seidel
    double start_GS = get_time();
    gauss_seidel(GS, N, indices, col_id, values, b, x, b_abs, r_rel);
    double end_GS = get_time();

    //printf("Gauss-Seidel: ");
/*    for (int i = 0; i < N; i ++)
    {
        printf("%.3f ", x[i]);
    }*/
    //printf("\n");
    printf("Gauss-Seidel runtime: %.6f seconds\n", end_GS - start_GS);
    fclose(GS);

    FILE *CG = fopen(file3, "w");
    if (!CG) 
    {
        perror("Error opening file: conjugate_gradient_n.txt");
        return 1;
    }

    init(N, x);

    // compute conjugate gradient
    double start_CG = get_time();
    CG_CSR(CG, N, indices, col_id, values, b, x, b_abs, r_rel);
    double end_CG = get_time();

    //printf("Conjugate Gradient: ");
/*    for (int i = 0; i < N; i ++)
    {
        printf("%.3f ", x[i]);
    }*/
    //printf("\n");
    printf("Conjugate Gradient runtime: %.6f seconds\n", end_CG - start_CG);
    fclose(CG);

    FILE *GMRES = fopen(file4, "w");
    if (!GMRES) 
    {
        perror("Error opening file: gmres_n.txt");
        return 1;
    }

    init(N, x);
    init_r0(N, b, b_abs, r);
    
    // compute generaized minimal residual
    double start_GMRES = get_time();
    GMRES_CSR(GMRES, N, indices, col_id, values, b, x, b_abs, r);
    double end_GMRES = get_time();

    //printf("GMRES: ");
/*    for (int i = 0; i < N; i ++)
    {
        printf("%.3f ", x[i]);
    } */
    //printf("\n");
    printf("GMRES runtime: %.6f seconds\n", end_GMRES - start_GMRES);
    fclose(GMRES);

    free(indices);
    free(col_id);
    free(values);
    free(b);
    free(x);
    free(x_new);
    free(r);

    return 0;
}
