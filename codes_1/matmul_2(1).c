/*************************************************************************************
* FILE: matmul_2.c
* DESCRIPTION:  
*   This serial program compute Matrix Multiply C = C + AB
*   It dynamically allocates contiguous memory space for 2D matrices,
*   uses ikj version, and
*   unroll i and j loops with unrolling factor = 4
* AUTHOR: Bing Bing Zhou
* Last revised: 12/03/2023
*************************************************************************************/
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
    double *c10; //auxiliary 1D for 2D matrix c1
    double **a; //the two-dimensional input matrix
    double **b; //the two-dimensional input matrix
    double **c; //the two-dimensional input/output matrix
    double **c1; //the two-dimensional input/output matrix

    int NRA, NCA, NCB; //matrices lengths
    int NRA0, NCB0; //divisible by 4 for loop unrolling

    double ai0, ai1, ai2, ai3, bj0, bj1, bj2, bj3;

    int i, j, k;

    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    if (argc == 4)
    {
        NRA = atoi(argv[1]);
        NCA = atoi(argv[2]);
        NCB = atoi(argv[3]);

        printf("NRA = %d, NCA = %d, NCB = %d\n", NRA, NCA, NCB);
    }
    else
    {
        printf("Usage: %s NRA NCA NCB \n\n"
               " NRA: matrix a row length\n"
               " NCA: matrix a column (or b row) length\n"
               " NCB:  matrix b column length\n\n", argv[0]);
        return 1;
    }

    NRA0 = NRA / 4 * 4;
    NCB0 = NCB / 4 * 4;

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

    c10 = (double *) malloc(NRA * NCB * sizeof(double));
    c1 = (double **) malloc(NRA * sizeof(double *));
    for (int i = 0; i < NRA; i++)
    {
        c1[i] = &(c10[i * NCB]);
    }

    printf("Initializing matrices...\n\n");
    srand(time(0)); // Seed the random number generator
//   srand(0); 
    for (i = 0; i < NRA; i++)
    {
        for (j = 0; j < NCA; j++)
        {
            a[i][j] = (double) rand() / RAND_MAX;
        }
    }
//         a[i][j]= (double) i + 1.0;

    for (i = 0; i < NCA; i++)
    {
        for (j = 0; j < NCB; j++)
        {
            b[i][j] = (double) rand() / RAND_MAX;
        }
    }
//         b[i][j]= (double) j;

    for (i = 0; i < NRA; i++)
    {
        for (j = 0; j < NCB; j++)
        {
            c[i][j] = 0.0;
        }
    }

    for (i = 0; i < NRA; i++)
    {
        for (j = 0; j < NCB; j++)
        {
            c1[i][j] = 0.0;
        }
    }

//  printf ("matrix a:\n");
//  print_matrix(a, NRA, NCA);
//  printf ("matrix b:\n");
//  print_matrix(b, NCA, NCB);
//  printf ("matrix C:\n");
//  print_matrix(c, NRA, NCB);

    /*** start matrix multiplication (ikj version) ***/
    printf("Starting matrix multiplication - ikj version\n\n");
    gettimeofday(&start_time, 0);
    for (i = 0; i < NRA; i++)
    {
        for (k = 0; k < NCA; k++)
        {
            ai0 = a[i][k];
            for (j = 0; j < NCB; j++)
            {
                c1[i][j] += ai0 * b[k][j];
            }
        }
    }
    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("ikj version takes %f seconds to finish the computation.\n\n", elapsed);


    /*** start matrix multiplication (ikj version with loop unrolling) ***/
    printf("Starting matrix multiplication - ikj version with loop unrolling\n\n");
    gettimeofday(&start_time, 0);
    //unroll i and j loops with unrolling factor = 4
    for (i = 0; i < NRA0; i += 4)
    {
        for (k = 0; k < NCA; k++)
        {
            ai0 = a[i][k];
            ai1 = a[i + 1][k];
            ai2 = a[i + 2][k];
            ai3 = a[i + 3][k];

            for (j = 0; j < NCB0; j += 4)
            {
                bj0 = b[k][j];
                bj1 = b[k][j + 1];
                bj2 = b[k][j + 2];
                bj3 = b[k][j + 3];

                c[i][j] += ai0 * bj0;
                c[i][j + 1] += ai0 * bj1;
                c[i][j + 2] += ai0 * bj2;
                c[i][j + 3] += ai0 * bj3;

                c[i + 1][j] += ai1 * bj0;
                c[i + 1][j + 1] += ai1 * bj1;
                c[i + 1][j + 2] += ai1 * bj2;
                c[i + 1][j + 3] += ai1 * bj3;

                c[i + 2][j] += ai2 * bj0;
                c[i + 2][j + 1] += ai2 * bj1;
                c[i + 2][j + 2] += ai2 * bj2;
                c[i + 2][j + 3] += ai2 * bj3;

                c[i + 3][j] += ai3 * bj0;
                c[i + 3][j + 1] += ai3 * bj1;
                c[i + 3][j + 2] += ai3 * bj2;
                c[i + 3][j + 3] += ai3 * bj3;
            }

            //For elelments in remaining j columns
            for (j = NCB0; j < NCB; j++)
            {
                bj0 = b[k][j];
                c[i][j] += ai0 * bj0;
                c[i + 1][j] += ai1 * bj0;
                c[i + 2][j] += ai2 * bj0;
                c[i + 3][j] += ai3 * bj0;
            }
        }
    }

    //For elements in remaining i rows
    for (i = NRA0; i < NRA; i++)
    {
        for (k = 0; k < NCA; k++)
        {
            ai0 = a[i][k];
            for (j = 0; j < NCB; j++)
            {
                c[i][j] += ai0 * b[k][j];
            }
        }
    }

    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("sequential calculation with loop unrolling time: %f\n\n", elapsed);

    //print the result matrix c
//   printf ("matrix C:\n");
//   print_matrix(c, NRA, NCB);


/*** check the correctness ***/
    printf("Starting comparison\n\n");
    int cnt = 0;
    for (i = 0; i < NRA; i++)
    {
        for (j = 0; j < NCB; j++)
        {
            if ((c[i][j] - c1[i][j]) * (c[i][j] - c1[i][j]) > 1.0E-10)
            {
                cnt++;
            }
        }
    }
    if (cnt == 0)
    {
        printf("Done. There are no differences!\n");
    }
    else
    {
        printf("Results are incorrect! The number of different elements is %d\n", cnt);
    }
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
