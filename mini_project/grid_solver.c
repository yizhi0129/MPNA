#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSILON 1e-6 // convergence criterion
#define MAX_ITER 1000 // maximum number of iterations

// store the matrix in CSR format
void generate_grid_CSR(int n, int * indices, int * col_id, float * values)
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
void jacobi(FILE *Jac, int N, int * indices, int * col_id, float * values, float * b, float * x, float * x_new, float b_abs, float r_rel)
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
            float sum = 0.0;
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
        float res_sum2 = 0.0;
        for (int i = 0; i < N; i ++)
        {
            float ax_i = 0.0;
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
void gauss_seidel(FILE *GS, int N, int * indices, int * col_id, float * values, float * b, float * x, float b_abs, float r_rel)
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
            float sum = 0.0;
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
        float res_sum2 = 0.0;
        for (int i = 0; i < N; i ++)
        {
            float ax_i = 0.0;
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
    float * values = (float *)malloc(nonzeros * sizeof(float));

    float * b = (float *)malloc(N * sizeof(float));
    float b_abs2 = 0.0, r_rel = 1.0; // relative residual initialized as 1.0 (could be any value greater than EPSILON)
    float * x = (float *)malloc(N * sizeof(float));
    float * x_new = (float *)malloc(N * sizeof(float));

    printf("b: ");
    for(int i = 0; i < N; i ++)
    {
        b[i] = (float)rand() / RAND_MAX + 1;  // random values [1, 2] 
        printf("%.2f ", b[i]);
        b_abs2 += b[i] * b[i];
    }
    printf("\n");
    float b_abs = sqrt(b_abs2); // pre-calculate ||b||
    printf("b_abs: %.2f\n", b_abs);

    generate_grid_CSR(n, indices, col_id, values);
    printf("Indices: ");
    for (int i = 0; i <= N; i ++) 
    {
        printf("%d ", indices[i]);
    }
    printf("\n");
    printf("Col_id: ");
    for (int i = 0; i < indices[N]; i ++) 
    {
        printf("%d ", col_id[i]);
    }
    printf("\n");
    printf("Values: ");
    for (int i = 0; i < indices[N]; i ++) 
    {
        printf("%.2f ", values[i]);
    }
    printf("\n");

    char file1[50], file2[50];
    sprintf(file1, "jacobi_%d.txt", n);
    sprintf(file2, "gauss_seidel_%d.txt", n);

    FILE *Jac = fopen(file1, "w");
    if (!Jac) {
        perror("Error opening file: jacobi_n.txt");
        return 1;
    }

    jacobi(Jac, N, indices, col_id, values, b, x, x_new, b_abs, r_rel);
    printf("Jacobi: ");
    for (int i = 0; i < N; i ++)
    {
        printf("%.2f ", x[i]);
    }
    printf("\n");
    fclose(Jac);

    FILE *GS = fopen(file2, "w");
    if (!GS) {
        perror("Error opening file: gauss_seidel_n.txt");
        return 1;
    }

    gauss_seidel(GS, N, indices, col_id, values, b, x, b_abs, r_rel);
    printf("Gauss-Seidel: ");
    for (int i = 0; i < N; i ++)
    {
        printf("%.2f ", x[i]);
    }
    printf("\n");
    fclose(GS);

    free(indices);
    free(col_id);
    free(values);
    free(b);
    free(x);
    free(x_new);

    return 0;
}
