/******************************************************************************
* FILE: matmul_1.c
* DESCRIPTION:  
*   Sequential Matrix Multiply C = AB
*   Compare the performance of ijk and ikj versions
* AUTHOR: Tianyi Zhang
* Last revised: 13/07/2023
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
    double *c10; //auxiliary 1D for 2D matrix c1
    double **a; //the two-dimensional input matrix
    double **b; //the two-dimensional input matrix
    double **c; //the resulting matrix (ijk version)
    double **c1; //the resulting matrix (ikj version)

    int NRA, NCA, NCB; //matrices lengths

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

    c10 = (double *) malloc(NRA * NCB * sizeof(double));
    c1 = (double **) malloc(NRA * sizeof(double *));
    for (int i = 0; i < NRA; i++)
    {
        c1[i] = &(c10[i * NCB]);
    }

    printf("Initializing matrices\n\n");
    srand(time(0)); // Seed the random number generator
    for (int i = 0; i < NRA; i++)
    {
        for (int j = 0; j < NCA; j++)
        {
            a[i][j] = (double) rand() / RAND_MAX;
        }
    }

    for (int i = 0; i < NCA; i++)
    {
        for (int j = 0; j < NCB; j++)
        {
            b[i][j] = (double) rand() / RAND_MAX;
        }
    }

    for (int i = 0; i < NRA; i++)
    {
        for (int j = 0; j < NCB; j++)
        {
            c[i][j] = 0.0;
        }
    }

    for (int i = 0; i < NRA; i++)
    {
        for (int j = 0; j < NCB; j++)
        {
            c1[i][j] = 0.0;
        }
    }

    printf("Starting matrix multiplication - ijk version\n\n");
    double ct;
    gettimeofday(&start_time, 0);
    for (int i = 0; i < NRA; i++)
    {
        for (int j = 0; j < NCB; j++)
        {
            ct = c[i][j];
            for (int k = 0; k < NCA; k++)
            {
                ct += a[i][k] * b[k][j];
            }
//            c[i][j] +=  a[i][k] * b[k][j];
//            由于数组的访问模式不同，可能无法有效地利用缓存。
            c[i][j] = ct;
        }
    }
    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("ijk version takes %f seconds to finish the computation.\n\n", elapsed);






// printf("******************************************************\n");
// DOWN HERE IS ANOTHER PROGRAM FOR COMPARISON，下面的代码用于对比
// DOWN HERE IS ANOTHER PROGRAM FOR COMPARISON，下面的代码用于对比
// printf("******************************************************\n");

    printf("Starting matrix multiplication - ikj version\n\n");
    double at;
    gettimeofday(&start_time, 0);
    double ct1, ct2, ct3, ct4;
    for (int i = 0; i < NRA; i++)
    {
        for (int k = 0; k < NCA; k++)
        {
            // 内层的无法弄出来的就不用变量代换
            at = a[i][k];
            for (int j = 0; j < NCB; j++)
            {
                c1[i][j] += at * b[k][j];
            }
        }
    }

    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("ikj version takes %f seconds to finish the computation.\n\n", elapsed);






// printf("******************************************************\n");
// 下面的代码用于检查
// 下面的代码用于检查
// printf("******************************************************\n");

/*** check the correctness ***/
    printf("Starting comparison\n\n");
    int cnt = 0;
    for (int i = 0; i < NRA; i++)
    {
        for (int j = 0; j < NCB; j++)
        {
            if ((c[i][j] - c1[i][j]) * (c[i][j] - c1[i][j]) > 1.0E-16)
            {
                cnt++;
            }
        }
    }

    if (cnt)
    {
        printf("results are not the same, the number of different elements is %d\n", cnt);
    }
    else
    {
        printf("Done. There are no differences!\n");
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
