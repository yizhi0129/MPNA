#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
void jacobi(FILE *Jac, int N, int *indices, int *col_id, double *values, double *b, double *x, double *x_new, double b_abs, double r_rel)
{
    int iter = 0;
    while (r_rel > EPSILON && iter < MAX_ITER) 
    { 
        for (int i = 0; i < N; i ++) 
        {
            int self = -1;
            double sum = 0.0;
            for (int j = indices[i]; j < indices[i + 1]; j++) 
            {
                if (col_id[j] != i) 
                {
                    sum += values[j] * x[col_id[j]];
                } 
                else 
                {
                    self = j; // diagonal element
                }
            }  
            if (self != -1) 
            {
                x_new[i] = (b[i] - sum) / values[self];  
            }
        }
        
        for (int i = 0; i < N; i ++) 
        {
            x[i] = x_new[i];
        }

        double res_sum2 = 0.0;
        for (int i = 0; i < N; i++) 
        {
            double ax_i = 0.0;
            for (int j = indices[i]; j < indices[i + 1]; j ++) 
            {
                ax_i += values[j] * x[col_id[j]];
            }
            res_sum2 += (ax_i - b[i]) * (ax_i - b[i]);
        }
        r_rel = sqrt(res_sum2) / b_abs;
        iter++;
        fprintf(Jac, "%d %.9f\n", iter, r_rel);
    }
}


// Gauss-Seidel method
void gauss_seidel(FILE *GS, int N, int *indices, int *col_id, double *values, double *b, double *x, double b_abs, double r_rel)
{
    int iter = 0;
    while (r_rel > EPSILON && iter < MAX_ITER) 
    {    
        for (int i = 0; i < N; i ++) 
        {
            double sum = 0.0;
            int self = -1;
            for (int j = indices[i]; j < indices[i + 1]; j++) 
            {
                if (col_id[j] == i) 
                {
                    self = j; // diagonal element
                }
            }

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

        double res_sum2 = 0.0;
        for (int i = 0; i < N; i ++) 
        {
            double ax_i = 0.0;
            for (int j = indices[i]; j < indices[i + 1]; j ++) 
            {
                ax_i += values[j] * x[col_id[j]];
            }
            res_sum2 += (ax_i - b[i]) * (ax_i - b[i]);     
        }
        r_rel = sqrt(res_sum2) / b_abs;
        iter++;
        fprintf(GS, "%d %.9f\n", iter, r_rel);
    }
}



// matrix-vector multiplication
void mat_vec_CSR(int N, int * indices, int * col_id, double * values, double * x, double * y)
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

// Conjugate Gradient Method
void CG_CSR(FILE *CG, int N, int *indices, int *col_id, double *values, double *b, double *x, double b_abs) 
{
    double *r = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));
    double *Ad = malloc(N * sizeof(double));

    mat_vec_CSR(N, indices, col_id, values, x, r);
    for (int i = 0; i < N; i ++) 
    {
        r[i] = b[i] - r[i];
        d[i] = r[i];
    }

    double r_d_old = dot_prod(N, r, r);
    double r_d_new;
    int count = 0;

    while (count < MAX_ITER && sqrt(r_d_old) / b_abs > EPSILON) 
    {
        mat_vec_CSR(N, indices, col_id, values, d, Ad);
        double d_Ad = dot_prod(N, d, Ad);
        double alpha = r_d_old / d_Ad;
        
        for (int i = 0; i < N; i++) 
        {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ad[i];
        }

        r_d_new = dot_prod(N, r, r);
        double beta = r_d_new / r_d_old;
        for (int i = 0; i < N; i++) 
        {
            d[i] = r[i] + beta * d[i];
        }

        r_d_old = r_d_new;
        count++;
        fprintf(CG, "%d %.9f\n", count, sqrt(r_d_new) / b_abs);
    }

    free(r);
    free(d);
    free(Ad);
}



// GMRES Method
void GMRES_CSR(FILE *GMRES, int N, int *indices, int *col_id, double *values, double *b, double *x, double b_abs) 
{
    double *r = malloc(N * sizeof(double));
    double *v[MAX_ITER + 1];
    double *h = malloc((MAX_ITER + 1) * MAX_ITER * sizeof(double));
    double *g = malloc((MAX_ITER + 1) * sizeof(double));
    
    for (int i = 0; i < MAX_ITER + 1; i ++) 
    {
        v[i] = malloc(N * sizeof(double));
    }
    
    mat_vec_CSR(N, indices, col_id, values, x, r);
    for (int i = 0; i < N; i ++) 
    {
        r[i] = b[i] - r[i];
    }
    
    g[0] = norm(N, r);
    for (int i = 0; i < N; i++) 
    {
        v[0][i] = r[i] / g[0];
    }
    
    for (int k = 0; k < MAX_ITER; k++) 
    {
        double *w = malloc(N * sizeof(double));
        mat_vec_CSR(N, indices, col_id, values, v[k], w);
        
        for (int j = 0; j <= k; j++) 
        {
            h[j * MAX_ITER + k] = dot_prod(N, w, v[j]);
            for (int i = 0; i < N; i++) 
            {
                w[i] -= h[j * MAX_ITER + k] * v[j][i];
            }
        }
        
        h[(k + 1) * MAX_ITER + k] = norm(N, w);
        if (h[(k + 1) * MAX_ITER + k] < EPSILON) 
        {
            break;
        }
        
        for (int i = 0; i < N; i++) 
        {
            v[k + 1][i] = w[i] / h[(k + 1) * MAX_ITER + k];
        }
        
        g[k + 1] = -g[k] * h[(k + 1) * MAX_ITER + k] / h[k * MAX_ITER + k];
        g[k] /= h[k * MAX_ITER + k];

        double r_rel = fabs(g[k + 1]) / b_abs;
        fprintf(GMRES, "%d %.9f\n", k + 1, r_rel);
        
        if (fabs(g[k + 1]) < EPSILON) 
        {
            break;
        }
        free(w);
    }
    
    for (int i = 0; i < MAX_ITER + 1; i ++) 
    {
        free(v[i]);
    }
    free(r);
    free(h);
    free(g);
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
    for(int i = 0; i < N; i ++)
    {
        b[i] = (double)rand() / RAND_MAX + 1;  // random values [1, 2] 
        b_abs2 += b[i] * b[i];
    }
    double b_abs = sqrt(b_abs2); // pre-calculate ||b||
    printf("b_abs: %.3f\n", b_abs);

    generate_grid_CSR(n, indices, col_id, values);

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
    CG_CSR(CG, N, indices, col_id, values, b, x, b_abs);
    double end_CG = get_time();

    printf("Conjugate Gradient runtime: %.6f seconds\n", end_CG - start_CG);
    fclose(CG);

    FILE *GMRES = fopen(file4, "w");
    if (!GMRES) 
    {
        perror("Error opening file: gmres_n.txt");
        return 1;
    }

    init(N, x);
    
    // compute generaized minimal residual
    double start_GMRES = get_time();
    GMRES_CSR(GMRES, N, indices, col_id, values, b, x, b_abs);
    double end_GMRES = get_time();

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
