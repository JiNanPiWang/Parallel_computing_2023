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

void print_matrix(double** T, int rows, int cols);

int main(int agrc, char* agrv[])
{
    double* a0; //auxiliary 1D for 2D matrix a
    double** a; //2D matrix for sequential computation
 
    int n; //input size
    int i, j, k;
    int indk;
    double c, amax;

    struct timeval start_time, end_time;
    long seconds, microseconds;
    double elapsed;

    if (agrc == 2)
    {
        n = atoi(agrv[1]);
        printf("The matrix size:  %d * %d \n", n, n);
    }
    else{
        printf("Usage: %s n\n\n"
               " n: the matrix size\n\n", agrv[0]);
        return 1;
    }

    printf("Creating and initializing matrices...\n\n"); 
    /*** Allocate contiguous memory for 2D matrices ***/
    a0 = (double*)malloc(n*n*sizeof(double));
    a = (double**)malloc(n*sizeof(double*));
    for (i=0; i<n; i++)
    {
        a[i] = a0 +  i*n;
    }

    srand(time(0));
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            a[i][j] = (double)rand()/RAND_MAX;

//    printf("matrix a: \n");
//    print_matrix(a, n, n);

    printf("Starting sequential computation...\n\n"); 
    /**** Sequential computation *****/
    gettimeofday(&start_time, 0);
    for (i=0; i<n-1; i++)
    {
        //find and record k where |a(k,i)|=𝑚ax|a(j,i)|
        amax = a[i][i];
        indk = i;
        for (k=i+1; k<n; k++)
        {
            if (fabs(a[k][i]) > fabs(amax))
            {
                amax = a[k][i];
                indk = k;
            }
        }

        //exit with a warning that a is singular
        if (amax == 0)
        {
            printf("matrix is singular!\n");
            exit(1);
        }  
	else if (indk != i) //swap row i and row k 
        {
            for (j=0; j<n; j++)
            {
                c = a[i][j];
                a[i][j] = a[indk][j];
                a[indk][j] = c;
            }
        } 

        //store multiplier in place of A(k,i)
        for (k=i+1; k<n; k++)
        {
            a[k][i] = a[k][i]/a[i][i];
        }

        //subtract multiple of row a(i,:) to zero out a(j,i)
        for (k=i+1; k<n; k++)
        { 
            c = a[k][i]; 
            for (j=i+1; j<n; j++)
            {
                a[k][j] -= c*a[i][j];
            }
        }
    }
    gettimeofday(&end_time, 0);
 
    //print the running time
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("sequential calculation time: %f\n\n",elapsed); 

}

void print_matrix(double** T, int rows, int cols)
{
	for (int i=0; i < rows; i++)
	{
		for (int j=0; j < cols; j++)
		{
			printf("%.2f   ", T[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
	return;
}

