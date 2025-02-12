#include "tools.h"


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
    double * x = (double *)malloc(N * sizeof(double));
    double * x_new = (double *)malloc(N * sizeof(double));

    printf("b: ");
    for(int i = 0; i < N; i ++)
    {
        b[i] = (double)rand() / RAND_MAX + 1;  // random values [1, 2] 
        printf("%.2f ", b[i]);
        b_abs2 += b[i] * b[i];
    }
    printf("\n");
    double b_abs = sqrt(b_abs2); // pre-calculate ||b||
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
