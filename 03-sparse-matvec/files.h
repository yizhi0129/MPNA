#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE 100


void read_size(char *filename, int *n, int *nnz)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        perror("Error opening file");
        exit(1);
    }

    char line[MAX_SIZE];
    int count = 0;

    while (fgets(line, sizeof(line), fp))
    {
        if (line[0] == '%')
        {
            continue;
        }
        else
        {
            if (count == 0)
            {
                sscanf(line, "%d %*d %d", n, nnz);
            }
            count ++;
        }
    }
    fclose(fp);   
}


// naïve sort
void bubble_sort(int *id, double *value, int n) 
{
    for (int i = 0; i < n - 1; i ++) 
    {
        for (int j = 0; j < n - 1 - i; j ++) 
        {
            if (id[j] > id[j + 1]) 
            {
                int temp_id = id[j];
                id[j] = id[j + 1];
                id[j + 1] = temp_id;

                double temp_val = value[j];
                value[j] = value[j + 1];
                value[j + 1] = temp_val;
            }
        }
    }
}

void read_matrix(char *filename, int n, int nnz, int *index, int *col_id, double *val)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        perror("Error opening file");
        exit(1);
    }

    char line[MAX_SIZE];
    int count = 0;
    int *row = malloc(nnz * sizeof(int));
    int *col = malloc(nnz * sizeof(int));
    double *A = malloc(nnz * sizeof(double));

    while (fgets(line, sizeof(line), fp))
    {
        if (line[0] == '%')
        {
            continue;
        }
        else
        {
            if (count == 0)
            {
                int n_check = 0, nnz_check = 0;
                sscanf(line, "%d %*d %d", &n_check, &nnz_check);
                int nnz_check1 = nnz_check * 2 - n_check;
                if (n_check != n || nnz_check1 != nnz)
                {
                    fprintf(stderr, "Matrix size does not match\n");
                    exit(1);
                }
            }
            else
            { 
                for (int k = 0; k < nnz; k ++)
                {
                    sscanf(line, "%d %d %lf", &row[k], &col[k], &A[k]);
                    if (row[k] != col[k])
                    {
                        row[k + 1] = col[k];
                        col[k + 1] = row[k];
                        A[k + 1] = A[k];
                        k ++;
                    }
                }
            }
            count ++;
        }
    }
    fclose(fp);

    // convert to CSR format (not sorted)
    for (int i = 0; i < nnz; i ++) 
    {
        col_id[i] = col[i] + row[i] * n;
        val[i] = A[i];
    }    

    free(A);
    free(row);
    free(col);

    // sort the CSR format and complete index array
    bubble_sort(col_id, val, nnz);
    int count_ind = 0;    
    index[0] = 0;
    for (int k = 0; k < n; k ++)
    {
        for (int i = n * k; i < n * (k + 1); i ++)
        {
    
        }
    }   
}

void gen_vector(char *filename, int n, double *x)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < n; i ++)
    {
        x[i] = rand() / (double)RAND_MAX;
        fprintf(fp, "%lf\n", x[i]);
    }
    fclose(fp);
}


void read_vector(char *filename, int n, double *x)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        perror("Error opening file");
        exit(1);
    }

    char line[MAX_SIZE];
    int count = 0;

    while (fgets(line, sizeof(line), fp))
    {
        sscanf(line, "%lf", &x[count]);
        count ++;
    }
    fclose(fp);
}

void write_vector(char *filename, int n, double *y)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < n; i ++)
    {
        fprintf(fp, "%lf\n", y[i]);
    }
    fclose(fp);
}