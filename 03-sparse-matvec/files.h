#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE 100

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
    int count2 = 0;
    int count3 = 0;
    index[0] = 0;
    int i = 0;
    int j = 0;    
    int temp = 1;

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
            else
            { 
                sscanf(line, "%d %d %lf", i, j, val[count - 1]);
                if (temp == j)
                {
                    count2 ++;
                }
                else
                {
                    index[count3 + 1] = index[count3] + count2;
                    count3 ++;
                    count2 = 1;
                    temp = j;
                }
                col_id[count - 1] = (j - 1) * n + (i - 1);
            }
            count ++;
        }
    }
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
        sscanf(line, "%lf", x[count]);
        count ++;
    }
    fclose(fp);
}