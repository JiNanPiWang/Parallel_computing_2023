/******************************************************************************
* FILE: matvec_0.c
* DESCRIPTION:  
*   A simple program for Matrix-vector Multiply b = Ax for
*   students to modify
* AUTHOR: Bing Bing Zhou
* Last revised: 12/07/2023 by Tianyi Zhang
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

    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    // 接收参数，行列大小
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

    // 分配矩阵空间 Allocate contiguous memory for 2D matrices
    a0 = (double *) malloc(NRA * NCA * sizeof(double)); // 指向一堆内容的指针，普通指针
    a = (double **) malloc(NRA * sizeof(double *)); // 指向一堆指针的指针，指针数组。数组大小为矩阵行数
    for (int i = 0; i < NRA; i++)
    {
        a[i] = a0 + i * NCA; // 初始化指针数组
    }

    // 分配向量空间 Allocate memory for vectors
    x = (double *) malloc(NCA * sizeof(double));
    b = (double *) malloc(NRA * sizeof(double));

    printf("Initializing matrix and vectors\n\n");
    srand(time(0)); // Seed the random number generator
    /*** Initialize matrix and vectors ***/
    // 初始化矩阵
    for (int i = 0; i < NRA; i++)
    {
        for (int k = 0; k < NCA; k++)
        {
            a[i][k] = (double) rand() / RAND_MAX;
        }
    }
    // 初始化向量
    for (int i = 0; i < NCA; i++)
    {
        x[i] = (double) rand() / RAND_MAX;
    }

    for (int i = 0; i < NRA; i++)
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
    gettimeofday(&start_time, 0); // gettimeofday：获取当前的系统时间。
    // 矩阵向量相乘，一行一行的算
    for (int i = 0; i < NRA; i++)
    {
        for (int k = 0; k < NCA; k++)
        {
            b[i] += a[i][k] * x[k];
        }
    }
    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("The computation takes %f seconds to complete.\n\n", elapsed);








    printf("Initializing matrix and vectors\n\n");
    srand(time(0)); // Seed the random number generator
    /*** Initialize matrix and vectors ***/
    // 初始化矩阵
    for (int i = 0; i < NRA; i++)
    {
        for (int k = 0; k < NCA; k++)
        {
            a[i][k] = (double) rand() / RAND_MAX;
        }
    }
    // 初始化向量
    for (int i = 0; i < NCA; i++)
    {
        x[i] = (double) rand() / RAND_MAX;
    }

    for (int i = 0; i < NRA; i++)
    {
        b[i] = 0.0;
    }

    // 使用循环展开的方法
    printf("Starting loop unrolling matrix-vector multiplication\n\n");
    gettimeofday(&start_time, 0); // gettimeofday：获取当前的系统时间。
    // 矩阵向量相乘，一行一行的算
    double tmp0, tmp1, tmp2, tmp3;
    for (int i = 0; i < NRA; i+=4)
    {
        tmp0 = b[i];
        tmp1 = b[i + 1];
        tmp2 = b[i + 2];
        tmp3 = b[i + 3];
        for (int k = 0; k < NCA; k++)
        {
            // 四列一起算
            tmp0 += a[i][k] * x[k];
            tmp1 += a[i + 1][k] * x[k];
            tmp2 += a[i + 2][k] * x[k];
            tmp3 += a[i + 3][k] * x[k];
        }
        b[i] = tmp0;
        b[i + 1] = tmp1;
        b[i + 2] = tmp2;
        b[i + 3] = tmp3;
    }
    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("The loop unrolling computation takes %f seconds to complete.\n\n", elapsed);

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


