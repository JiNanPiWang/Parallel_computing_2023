/******************************************************************************
* FILE: matmul_0.c
* DESCRIPTION:  
*   A simple program for Matrix Multiply C = AB (ijk version) for 
*   students to modify
* AUTHOR: Bing Bing Zhou
* Last revised: 07/07/2022
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

void print_matrix(double **T, int rows, int cols);

int main(int argc, char *argv[])
{
    double *a0; //auxiliary 1D for 2D matrix a
    double *b0; //auxiliary 1D for 2D matrix b
    double *c0; //auxiliary 1D for 2D matrix c
    double **a; //the two-dimensional input matrix
    double **b; //the two-dimensional input matrix
    double **c; //the resulting matrix (ijk version)

    int NRA, NCA, NCB; //matrices lengths

    int i, j, k;
    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    if (argc == 4)
    {
        NRA = atoi(argv[1]);
        NCA = atoi(argv[2]);
        NCB = atoi(argv[3]);

        printf("NRA = %d, NCA = %d, NCB = %d\n\n", NRA, NCA, NCB);
    }
    else
    {
        printf("Usage: %s NRA NCA NCB \n\n"
               " NRA: matrix a row length\n"
               " NCA: matrix a column (or b row) length\n"
               " NCB:  matrix b column length\n\n", argv[0]);
        return 1;
    }

    // Allocate contiguous memory for 2D matrices
    a0 = (double *) malloc(NRA * NCA * sizeof(double));
    a = (double **) malloc(NRA * sizeof(double *));
    for (int i = 0; i < NRA; i++)
    {
        a[i] = a0 + i * NCA;
    }

    b0 = (double *) malloc(NCA * NCB * sizeof(double));
    b = (double **) malloc(NCA * sizeof(double *));
    for (int i = 0; i < NCA; i++)
    {
        b[i] = &(b0[i * NCB]);
    }

    c0 = (double *) malloc(NRA * NCB * sizeof(double));
    c = (double **) malloc(NRA * sizeof(double *));
    for (int i = 0; i < NRA; i++)
    {
        c[i] = &(c0[i * NCB]);
    }

    printf("Initializing matrices\n\n");
    srand(time(0)); // Seed the random number generator
    for (i = 0; i < NRA; i++)
    {
        for (j = 0; j < NCA; j++)
        {
            a[i][j] = (double) rand() / RAND_MAX;
        }
    }

    for (i = 0; i < NCA; i++)
    {
        for (j = 0; j < NCB; j++)
        {
            b[i][j] = (double) rand() / RAND_MAX;
        }
    }

    for (i = 0; i < NRA; i++)
    {
        for (j = 0; j < NCB; j++)
        {
            c[i][j] = 0.0;
        }
    }

//  printf ("matrix a:\n");
//  print_matrix(a, NRA, NCA);
//  printf ("matrix b:\n");
//  print_matrix(b, NCA, NCB);
//  printf ("matrix C:\n");
//  print_matrix(c, NRA, NCB);

    printf("Starting matrix multiplication - ijk version\n\n");
    gettimeofday(&start_time, 0);
    double ct;
    for (i = 0; i < NRA; i++)
    {
        for (j = 0; j < NCB; j++)
        {
            ct = c[i][j];
            for (k = 0; k < NCA; k++)
            {
                ct += a[i][k] * b[k][j];
            }
            c[i][j] = ct;
        }
    }
    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("ijk version takes %f seconds to finish the computation.\n\n", elapsed);

// printf("******************************************************\n");
// printf("Result Matrix:\n");
// print_matrix(c, NRA, NCB);
// printf("******************************************************\n");

}

void print_matrix(double **T, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%.2f  ", T[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}
