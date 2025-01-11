#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSILON 1e-7 // convergence criterion

void generate_grid_CSR(int n, int * indices, int * col_id, float * values)
{
    int count = 0; // total number of non zeros to be stored
    indices[0] = 0;

    for (int i = 0; i < n; i ++)
    {
        for (int j = 0; j < n; j ++)
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
            else if (j < n - 1)
            {
                col_id[count] = index + n; // neighbour above
                values[count] = -1.0;
                count ++;
            }
            indices[index + 1] = count;         
        }               
    }
    // allocate memory for col_id and values according to non zeros' count
    col_id = (int *)malloc(count * sizeof(int)); 
    values = (float *)malloc(count * sizeof(float));
}

void jacobi(int N, int * indices, int * col_id, float * values, float * b, float * x, float * x_new, float b_abs)
{
    float r_relative = 1.0;
    for(int i = 0; i < N; i ++)
    {
        x[i] = b[i];
    }
    while (r_relative > EPSILON)
    {
        float temp = 0.0;
        for (int i = 0; i < N; i ++)
        {
            x_new[i] = x[i];
            for (int j = 0; j < i; j ++)
            {
                x_new[i] -= values[indices[i]] * x[i]; // check index 
            }
            for (int j = i + 1; j < N; j ++)
            {
                x_new[i] -= values[indices[i]] * x[i];
            }
            x_new[i] /= values[indices[i]]; 
            x[i] = x_new[i];  
            temp += (x[i] * values[i] - b[i]) * (x[i] * values[i] - b[i]);
        }
        r_relative = sqrt(temp) / b_abs;
    }
}

void gauss_seidel(int N, int * indices, int * col_id, float * values, float * b, float * x, float b_abs)
{
    float r_relative = 1.0;
    for(int i = 0; i < N; i ++)
    {
        x[i] = b[i];
    }
    while (r_relative > EPSILON)
    {
        float temp = 0.0;
        for (int i = 0; i < N; i ++)
        {
            for (int j = 0; j < i; j ++)
            {
                x[i] -= values[indices[i]] * x[i]; // check index
            }
            for (int j = i + 1; j < N; j ++)
            {
                x[i] -= values[indices[i]] * x[i];
            }
            x[i] /= values[indices[i]]; 
            temp += (x[i] * values[i] - b[i]) * (x[i] * values[i] - b[i]);     
        }
        r_relative = sqrt(temp) / b_abs;
    }
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        printf("Please enter a number n: grid size");
        return 1;
    }
    int n = atoi(argv[1]);
    int N = n * n;
    int * indices = (int *)malloc((N + 1) * sizeof(int));
    int * col_id;
    float * values;

    float * b = (float *)malloc(N * sizeof(float));
    float b_abs2 = 0.0;
    float * x = (float *)malloc(N * sizeof(float));
    float * x_new = (float *)malloc(N * sizeof(float));
    for(int i = 0; i < N; i ++)
    {
        b[i] = rand() % 100;
        b_abs2 += b[i] * b[i];
    }
    float b_abs = sqrt(b_abs2);

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
        printf("%.1f ", values[i]);
    }
    printf("\n");

    //jacobi(N, indices, col_id, values, b, x, x_new, b_abs);
    //gauss_seidel(N, indices, col_id, values, b, x, b_abs);

    free(indices);
    free(col_id);
    free(values);
    free(b);
    free(x);
    free(x_new);

    return 0;
}
