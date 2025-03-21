#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SIZE 1024


void read_size(char *filename, int *n, int *nnz)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        perror("Error opening file: read_size");
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


// bubble_sort row, column
void bubble_sort(int *row, int *col, double *value, int nnz) 
{
    for (int i = 0; i < nnz - 1; i++) 
    {
        for (int j = 0; j < nnz - i - 1; j++) 
        {
            if ((row[j] > row[j + 1]) || (row[j] == row[j + 1] && col[j] > col[j + 1])) 
            {
                int temp_r = row[j];
                row[j] = row[j + 1];
                row[j + 1] = temp_r;

                int temp_c = col[j];
                col[j] = col[j + 1];
                col[j + 1] = temp_c;

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
        perror("Error opening file: read_matrix");
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
        count++;

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

            if (temp_r != temp_c)
            {
                row[sym] = temp_c - 1;
                col[sym] = temp_r - 1;
                A[sym] = temp_A;
                sym ++;
            }
        }           

        // sort
        bubble_sort(row, col, A, nnz);

        // CSR
        for (int i = 0; i < nnz; i++) 
        {
            col_id[i] = col[i]; // change to column index instead of global index
            val[i] = A[i];
        }    

        // complete index
        index[0] = 0;
        int current_row = 0;
        for (int i = 0; i < nnz; i ++)
        {
            while (row[i] > current_row)
            {
                current_row ++;
                index[current_row] = i;
            }
        }
        for (; current_row < n; current_row ++)
        {
            index[current_row + 1] = nnz;
        }
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
        perror("Error opening file: gen_vector");
        exit(1);
    }

    for (int i = 0; i < n; i ++)
    {
        x[i] = rand() / (double)RAND_MAX;
        fprintf(fp, "%.10f\n", x[i]);
    }
    fclose(fp);
}


void read_vector(char *filename, int n, double *x)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        perror("Error opening file: read_vector");
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
        perror("Error opening file: write_vector");
        exit(1);
    }

    for (int i = 0; i < n; i ++)
    {
        fprintf(fp, "%.10f\n", y[i]);
    }
    fclose(fp);
}

// y = A * x ( A in CSR format)
void matvec(int n, int *index, int *col_id, double *val, double *x, double *y)
{
    for (int i = 0; i < n; i ++)
    {
        y[i] = 0;
        for (int j = index[i]; j < index[i + 1]; j ++)
        {
            y[i] += val[j] * x[col_id[j]];
        }
    }
}

void normalize(int n, double *x)
{
    double sum = 0;
    for (int i = 0; i < n; i ++)
    {
        sum += x[i] * x[i];
    }
    sum = sqrt(sum);
    for (int i = 0; i < n; i ++)
    {
        x[i] /= sum;
    }
}

double dot_vec(int n, double *x, double *y)
{
    double result = 0.0;
    for (int i = 0; i < n; i ++)
    {
        result += x[i] * y[i];
    }
    return result;
}

void diff_vec(int n, double *x, double *y, double *result)
{
    for (int i = 0; i < n; i ++)
    {
        result[i] = x[i] - y[i];
    }
}

void eq_vec(int n, double *x, double *y)
{
    for (int i = 0; i < n; i ++)
    {
        y[i] = x[i];
    }
}