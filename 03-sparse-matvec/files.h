#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE 1024


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


// naive bubble sort
void bubble_sort(int *id, double *value, int n) 
{
    for (int i = 0; i < n; i ++) 
    {
        for (int j = 0; j < n - i; j ++) 
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
        int n_check = 0, nnz_check = 0;
        if (line[0] == '%')
        {
            continue;
        }
        else
        {
            if (count == 0)
            {
                sscanf(line, "%d %*d %d", &n_check, &nnz_check);
                int nnz_check1 = nnz_check * 2 - n_check;
                if (n_check != n || nnz_check1 != nnz)
                {
                    fprintf(stderr, "Matrix size does not match\n");
                    exit(1);
                }
            }
            count ++;
            
            int sym = nnz_check; 
            for (int k = 0; k < nnz_check; k ++)
            {
                if (fgets(line, sizeof(line), fp) == NULL) 
                {
                    fprintf(stderr, "Fail to read line %d\n", k + 1);
                    break;
                }
            
                int temp_r = 0, temp_c = 0;
                double temp_A = 0.0;
                if (sscanf(line, "%d %d %lf", &temp_r, &temp_c, &temp_A) != 3) 
                {
                    fprintf(stderr, "Fail to load data line %d %s\n", k + 1, line);
                    continue;
                }

                row[k] = temp_r - 1;
                col[k] = temp_c - 1;
                A[k] = temp_A;    
                printf("Read: k = %d, row = %d, col = %d, A = %.lf\n", k, row[k], col[k], A[k]);
    
                if (temp_r != temp_c)
                {
                    row[sym] = temp_c - 1;
                    col[sym] = temp_r - 1;
                    A[sym] = temp_A;
                    printf("Complete symmetry: k = %d, row = %d, col = %d, A = %.lf\n", sym, row[sym], col[sym], A[sym]);
                    sym ++;
                }
            }           
        }

        // convert to CSR format (not sorted)
        printf("Matrix: CSR not sorted\n");
        for (int i = 0; i < nnz; i ++) 
        {
            col_id[i] = col[i] * n + row[i];
            val[i] = A[i];
            printf("%d %d %lf\n", i, col_id[i], val[i]);
        }    
        printf("------------------------------\n");

        // sort the CSR format and complete index array
        bubble_sort(col_id, val, nnz);
        printf("Matrix: CSR sorted\n");
        for (int i = 0; i < nnz; i ++)
        {
            printf("%d %d %lf\n", i, col_id[i], val[i]);
        }
        printf("------------------------------\n");

        // complete the index array
        int count_ind = 0;    
        index[0] = 0;
        printf("Index: ");
        for (int k = 1; k < n; k ++)
        {
            for (int i = 0; i < nnz; i ++)
            {
                if (col_id[i] < n * k)
                {
                    count_ind ++;
                }
                else
                {
                    break;
                }
                
            }
            index[k] = col_id[count_ind + 1];
            printf("%d %d\n", k + 1, index[k]);
        }
        index[n] = nnz;   
        printf("------------------------------\n");
    }
    
    free(A);
    free(row);
    free(col);

    fclose(fp);
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