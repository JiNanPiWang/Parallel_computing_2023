/******************************************************************************
* FILE: matvec_0.c
* DESCRIPTION:  
*   A simple program for Matrix-vector Multiply b = Ax for
*   students to modify
* AUTHOR: Bing Bing Zhou
* Last revised: 05/07/2022
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

void print_matrix(double **T, int rows, int cols);

void print_vector(double *T, int cols);

int main(int argc, char *argv[])
{
    double *a0; //auxiliary 1D array to make a contiguously allocated
    double **a; //the two-dimensional input matrix
    double *x; //input vector
    double *b; //the resulting vector

    int NRA, NCA; //matrix size

    int i, k;
    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    if (argc == 3)
    {
        NRA = atoi(argv[1]);
        NCA = atoi(argv[2]);

        printf("NRA = %d, NCA = %d\n", NRA, NCA);
    }
    else
    {
        printf("Usage: %s NRA NCA\n\n"
               " NRA: matrix a row length\n"
               " NCA: matrix a column (or x) length\n\n", argv[0]);
        return 1;
    }

    // Allocate contiguous memory for 2D matrices
    a0 = (double *) malloc(NRA * NCA * sizeof(double));
    a = (double **) malloc(NRA * sizeof(double *));
    for (int i = 0; i < NRA; i++)
    {
        a[i] = a0 + i * NCA;
    }

    //Allocate memory for vectors
    x = (double *) malloc(NCA * sizeof(double));
    b = (double *) malloc(NRA * sizeof(double));

    printf("Initializing matrix and vectors\n\n");
    srand(time(0)); // Seed the random number generator
    /*** Initialize matrix and vectors ***/
    for (i = 0; i < NRA; i++)
    {
        for (k = 0; k < NCA; k++)
        {
            a[i][k] = (double) rand() / RAND_MAX;
        }
    }

    for (i = 0; i < NCA; i++)
    {
        x[i] = (double) rand() / RAND_MAX;
    }

    for (i = 0; i < NRA; i++)
    {
        b[i] = 0.0;
    }

/* 
  printf ("matrix a:\n");
  print_matrix(a, NRA, NCA);
  printf ("vector x:\n");
  print_vector(x, NCA);
  printf ("vector b:\n");
  print_vector(b, NRA);
*/

    printf("Starting matrix-vector multiplication\n\n");
    gettimeofday(&start_time, 0);
    for (i = 0; i < NRA; i++)
    {
        for (k = 0; k < NCA; k++)
        {
            b[i] += a[i][k] * x[k];
        }
    }
    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("The computation takes %f seconds to complete.\n\n", elapsed);


/*** Print results ***/
// printf("******************************************************\n");
// printf("Resulting vector:\n");
// print_vector(b, NRA);
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

void print_vector(double *T, int cols)
{
    for (int i = 0; i < cols; i++)
    {
        printf("%.2f  ", T[i]);
    }
    printf("\n\n");
    return;
}


