/************************************************************************************
* FILE: gepp_0.c
* DESCRIPTION:
* sequential program for Gaussian elimination with partial pivoting 
* for student to modify
* AUTHOR: Bing Bing Zhou
* LAST REVISED: 01/06/2023
*************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

void print_matrix(double **T, int rows, int cols);

void init(double **T1, double **T2, int rows, int cols);

void calculate_matrix(double **T, int rows, int cols);

int compare_matrix(double **T1, double **T2, int rows, int cols);

int main(int agrc, char *agrv[])
{
    double *a0, *a1; //auxiliary 1D for 2D matrix a
    double **a, **a_unroll; //2D matrix for sequential computation

    int n; //input size
    int i, j, k;
    int indk;
    double c, amax;
    const int unrol_fac = 4;

    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    if (agrc == 2)
    {
        n = atoi(agrv[1]);
        printf("The matrix size:  %d * %d \n", n, n);
    }
    else
    {
        printf("Usage: %s n\n\n"
               " n: the matrix size\n\n", agrv[0]);
        return 1;
    }

    printf("Creating and initializing matrices...\n\n");
    /*** Allocate contiguous memory for 2D matrices ***/
    a0 = (double *) malloc(n * n * sizeof(double));
    a1 = (double *) malloc(n * n * sizeof(double));
    a = (double **) malloc(n * sizeof(double *));
    a_unroll = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        a[i] = a0 + i * n;
        a_unroll[i] = a1 + i * n;
    }

    // 使用随机数填充矩阵a的元素。
    srand(time(0));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            a[i][j] = (double) rand() / RAND_MAX;
            a_unroll[i][j] = a[i][j];
        }
    }

    printf("Starting sequential computation...\n\n");
    /**** 无展开的高斯消元 *****/
    /**** 无展开的高斯消元 *****/
    /**** 无展开的高斯消元 *****/
    /**** 无展开的高斯消元 *****/
    /**** 无展开的高斯消元 *****/
    gettimeofday(&start_time, 0);
    for (i = 0; i < n - 1; i++)
    {
        // 制造一个从左上到右下对角线值为1的上三角矩阵

        // find and record k where |a(k,i)|=𝑚ax|a(j,i)|
        // 对于每一列i，找到绝对值最大的元素a(k,i)，并记录其行号indk。
        amax = a[i][i];
        indk = i;
        for (k = i + 1; k < n; k++)
        {
            if (fabs(a[k][i]) > fabs(amax))
            {
                amax = a[k][i];
                indk = k;
            }
        }

        // exit with a warning that a is singular
        // 如果最大元素为0，则说明矩阵是奇异的（不可逆），程序退出。
        if (amax == 0)
        {
            printf("matrix is singular!\n");
            exit(1);
        }
        else if (indk != i) //swap row i and row k
        {
            // 如果最大元素不在当前行i上，则交换行i和行k。
            for (j = 0; j < n; j++)
            {
                c = a[i][j];
                a[i][j] = a[indk][j];
                a[indk][j] = c;
            }
        }

        // store multiplier in place of A(k,i)
        // 将除第i行外的第i列元素都除以a(i,i)，得到乘数。
        for (k = i + 1; k < n; k++)
        {
            a[k][i] = a[k][i] / a[i][i];
        }

        // subtract multiple of row a(i,:) to zero out a(j,i)
        // 将除第i行外的第i列元素都除以a(i,i)，得到乘数。
        for (k = i + 1; k < n; k++)
        {
            c = a[k][i];
            for (j = i + 1; j < n; j++)
            {
                a[k][j] -= c * a[i][j];
            }
        }
    }
    gettimeofday(&end_time, 0);

    //print the running time
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("sequential calculation time: %f\n\n", elapsed);





    printf("Starting sequential computation...\n\n");
    /**** 有展开的高斯消元 *****/
    /**** 有展开的高斯消元 *****/
    /**** 有展开的高斯消元 *****/
    /**** 有展开的高斯消元 *****/
    /**** 有展开的高斯消元 *****/
    gettimeofday(&start_time, 0);
    for (i = 0; i < n - 1; i++)
    {
        // 制造一个从左上到右下对角线值为1的上三角矩阵

        // find and record k where |a(k,i)|=𝑚ax|a(j,i)|
        // 对于每一列i，找到绝对值最大的元素a(k,i)，并记录其行号indk。
        amax = a_unroll[i][i];
        indk = i;
        for (k = i + 1; k < n; k++)
        {
            if (fabs(a_unroll[k][i]) > fabs(amax))
            {
                amax = a_unroll[k][i];
                indk = k;
            }
        }

        // exit with a warning that a is singular
        // 如果最大元素为0，则说明矩阵是奇异的（不可逆），程序退出。
        if (amax == 0)
        {
            printf("matrix is singular!\n");
            exit(1);
        }
        else if (indk != i) //swap row i and row k
        {
            // 如果最大元素不在当前行i上，则交换行i和行indk。
            for (j = 0; j < n; j++)
            {
                c = a_unroll[i][j];
                a_unroll[i][j] = a_unroll[indk][j];
                a_unroll[indk][j] = c;
            }
        }

        // store multiplier in place of A(k,i)
        // 将除第i行外的第i列元素都除以a(i,i)，得到乘数，第i行i列即为1。
        for (k = i + 1; k < n / unrol_fac * (unrol_fac - 1); k += unrol_fac)
        {
            a_unroll[k][i] = a_unroll[k][i] / a_unroll[i][i];
            a_unroll[k + 1][i] = a_unroll[k + 1][i] / a_unroll[i][i];
            a_unroll[k + 2][i] = a_unroll[k + 2][i] / a_unroll[i][i];
            a_unroll[k + 3][i] = a_unroll[k + 3][i] / a_unroll[i][i];
        }
        for (; k < n; k++)
        {
            a_unroll[k][i] = a_unroll[k][i] / a_unroll[i][i];
        }

        // subtract multiple of row a(i,:) to zero out a(j,i)
        // 将除第i行外的第i列元素都除以a(i,i)，得到乘数。
        double c0, c1, c2, c3;
        for (k = i + 1; k < n / unrol_fac * (unrol_fac - 1); k += unrol_fac)
        {
            c0 = a_unroll[k][i];
            c1 = a_unroll[k + 1][i];
            c2 = a_unroll[k + 2][i];
            c3 = a_unroll[k + 3][i];
            for (j = i + 1; j < n; j++)
            {
                a_unroll[k][j] -= c0 * a_unroll[i][j];
                a_unroll[k + 1][j] -= c1 * a_unroll[i][j];
                a_unroll[k + 2][j] -= c2 * a_unroll[i][j];
                a_unroll[k + 3][j] -= c3 * a_unroll[i][j];
            }
        }
        for (; k < n; k++)
        {
            c = a[k][i];
            for (j = i + 1; j < n; j++)
            {
                a_unroll[k][j] -= c * a_unroll[i][j];
            }
        }
    }
    gettimeofday(&end_time, 0);

    //print the running time
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("sequential calculation(with unrolling) time: %f\n\n", elapsed);

    if (compare_matrix(a, a_unroll, n, n) == 1)
        printf("Results are the same\n");
    else
        printf("Results are not the same\n");

}

void print_matrix(double **T, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%.2f   ", T[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

int compare_matrix(double **T1, double **T2, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (T1[i][j] != T2[i][j])
                return 0;
    return 1;
}

