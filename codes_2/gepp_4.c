/********************************************************************************
* FILE: gepp_4.c
* DESCRIPTION:
* The C program for Gaussian elimination with partial pivoting
* Try to use loop unrolling to improve the performance - fourth attempt
* The matrix is divided into blocks, each being of 4X4 (unrolling factor = 4) 
* The algorithm takes several steps:
* 1. eliminate the main block column
* 2. update the main block row using the multipliers stored in lower triangle of
*    the main block submatrix (4X4)
* 3. rank 4 updating for the trailing submatrix
* 4. repeat 1, 2, and 3 and finally, eliminate a small trailing submatrix < 4X4
* Performance is significantly improved - over twice as fast
* AUTHOR: Bing Bing Zhou
* LAST REVISED: 06/06/2023
*********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

void print_matrix(double** T, int rows, int cols);
int test(double** t1, double** t2, int rows);

int main(int agrc, char* agrv[])
{
    double* a0; //auxiliary 1D for 2D matrix a
    double** a; //2D matrix for sequential computation
    double* d0; //auxiliary 1D for 2D matrix d
    double** d; //2D matrix, same initial data as a for the computation with unrolling
    int n; //input size
    int nb;
    int i, ib, j, k;
    int indk;
    double amax;

    double di00, di01, di02, di03;
    double di10, di11, di12, di13;
    double di20, di21, di22, di23;
    double di30, di31, di32, di33;

    double dj00, dj01, dj02, dj03;
    double dj10, dj11, dj12, dj13;
    double dj20, dj21, dj22, dj23;
    double dj30, dj31, dj32, dj33;

    double c;
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
               " n: the matrix size\n", agrv[0]);
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
    d0 = (double*)malloc(n*n*sizeof(double));
    d = (double**)malloc(n*sizeof(double*));
    for (i=0; i<n; i++)
    {
        d[i] = d0 +  i*n;
    }

    srand(time(0));
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            a[i][j] = (double)rand()/RAND_MAX;
            d[i][j] = a[i][j];
        }
    }
//    printf("matrix a: \n");
//    print_matrix(a, n, n);
//    printf("matrix d: \n");
//    print_matrix(d, n, n);

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

       //store multiplier in place of A(j,i)
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


    printf("Starting sequential computation with loop unrolling...\n\n"); 

    /***sequential computation with loop unrolling with unrolling factor = 4***/
    gettimeofday(&start_time, 0);
    nb = n/4 * 4;
    for (i=0; i<nb; i+=4)
    {
       /***elimination for the main block of 4 columns***/
       for (ib=i; ib<i+4; ib++)
       {
          amax = d[ib][ib];
          indk = ib;
          for (k=ib+1; k<n; k++)
          if (fabs(d[k][ib]) > fabs(amax))
          {
             amax = d[k][ib];
             indk = k;
          }

          if (amax == 0.0)
          {
             printf("the matrix is singular\n");
             exit(1);
          }
          else if (indk != ib) //swap row ii and row k 
             {
             for (j=0; j<n; j++)
             {
                c = d[ib][j];
                d[ib][j] = d[indk][j];
                d[indk][j] = c;
             }
          } 
 
          for (k=ib+1; k<n; k++)
             d[k][ib] = d[k][ib]/d[ib][ib];

          //subtract multiple of row d(ib,:) to zero out d(k, ib)
          for (k=ib+1; k<n; k++)
          { 
             c = d[k][ib]; 
             for (j=ib+1; j<i+4; j++)
             {
                d[k][j] -= c*d[ib][j];
             }
          }
       }

       /***update the main block of 4 rows***/       di10 = d[i+1][i]; 
       di20 = d[i+2][i]; di21 = d[i+2][i+1]; 
       di30 = d[i+3][i]; di31 = d[i+3][i+1]; di32 = d[i+3][i+2]; 

       for (j=i+4; j<nb; j+=4)
       {
          dj00 = d[i][j];   dj01 = d[i][j+1];   dj02 = d[i][j+2];   dj03 = d[i][j+3];

          d[i+1][j]   = d[i+1][j]   - di10*dj00;
          d[i+1][j+1] = d[i+1][j+1] - di10*dj01;
          d[i+1][j+2] = d[i+1][j+2] - di10*dj02;
          d[i+1][j+3] = d[i+1][j+3] - di10*dj03;

          dj10 = d[i+1][j]; dj11 = d[i+1][j+1]; dj12 = d[i+1][j+2]; dj13 = d[i+1][j+3];

          d[i+2][j]   = d[i+2][j]   - di20*dj00 - di21*dj10;
          d[i+2][j+1] = d[i+2][j+1] - di20*dj01 - di21*dj11;
          d[i+2][j+2] = d[i+2][j+2] - di20*dj02 - di21*dj12;
          d[i+2][j+3] = d[i+2][j+3] - di20*dj03 - di21*dj13;

          dj20 = d[i+2][j]; dj21 = d[i+2][j+1]; dj22 = d[i+2][j+2]; dj23 = d[i+2][j+3];

          d[i+3][j]   = d[i+3][j]   - di30*dj00 - di31*dj10 - di32*dj20;
          d[i+3][j+1] = d[i+3][j+1] - di30*dj01 - di31*dj11 - di32*dj21;
          d[i+3][j+2] = d[i+3][j+2] - di30*dj02 - di31*dj12 - di32*dj22;
          d[i+3][j+3] = d[i+3][j+3] - di30*dj03 - di31*dj13 - di32*dj23;
       }

       for (j=nb; j<n; j++)
       {
          d[i+1][j] = d[i+1][j] - di10*d[i][j];
          d[i+2][j] = d[i+2][j] - di20*d[i][j] - di21*d[i+1][j];
          d[i+3][j] = d[i+3][j] - di30*d[i][j] - di31*d[i+1][j] - di32*d[i+2][j];
       }

      /***update the trailing submatrix (rank-4 updating)***/
      for (k=i+4; k<nb; k+=4)
       {
          di00 = d[k][i];   di01 = d[k][i+1];   di02 = d[k][i+2];   di03 = d[k][i+3];
          di10 = d[k+1][i]; di11 = d[k+1][i+1]; di12 = d[k+1][i+2]; di13 = d[k+1][i+3];
          di20 = d[k+2][i]; di21 = d[k+2][i+1]; di22 = d[k+2][i+2]; di23 = d[k+2][i+3];
          di30 = d[k+3][i]; di31 = d[k+3][i+1]; di32 = d[k+3][i+2]; di33 = d[k+3][i+3];

          for (j=i+4; j<nb; j+=4)
          {
             dj00 = d[i][j];   dj01 = d[i][j+1];   dj02 = d[i][j+2];   dj03 = d[i][j+3];
             dj10 = d[i+1][j]; dj11 = d[i+1][j+1]; dj12 = d[i+1][j+2]; dj13 = d[i+1][j+3];
             dj20 = d[i+2][j]; dj21 = d[i+2][j+1]; dj22 = d[i+2][j+2]; dj23 = d[i+2][j+3];
             dj30 = d[i+3][j]; dj31 = d[i+3][j+1]; dj32 = d[i+3][j+2]; dj33 = d[i+3][j+3];

             d[k][j]   = d[k][j]   - di00*dj00 - di01*dj10 - di02*dj20 - di03*dj30;  
             d[k][j+1] = d[k][j+1] - di00*dj01 - di01*dj11 - di02*dj21 - di03*dj31;
             d[k][j+2] = d[k][j+2] - di00*dj02 - di01*dj12 - di02*dj22 - di03*dj32;
             d[k][j+3] = d[k][j+3] - di00*dj03 - di01*dj13 - di02*dj23 - di03*dj33;

             d[k+1][j]   = d[k+1][j]   - di10*dj00 - di11*dj10 - di12*dj20 - di13*dj30;  
             d[k+1][j+1] = d[k+1][j+1] - di10*dj01 - di11*dj11 - di12*dj21 - di13*dj31;
             d[k+1][j+2] = d[k+1][j+2] - di10*dj02 - di11*dj12 - di12*dj22 - di13*dj32;
             d[k+1][j+3] = d[k+1][j+3] - di10*dj03 - di11*dj13 - di12*dj23 - di13*dj33;

             d[k+2][j]   = d[k+2][j]   - di20*dj00 - di21*dj10 - di22*dj20 - di23*dj30;  
             d[k+2][j+1] = d[k+2][j+1] - di20*dj01 - di21*dj11 - di22*dj21 - di23*dj31;
             d[k+2][j+2] = d[k+2][j+2] - di20*dj02 - di21*dj12 - di22*dj22 - di23*dj32;
             d[k+2][j+3] = d[k+2][j+3] - di20*dj03 - di21*dj13 - di22*dj23 - di23*dj33;

             d[k+3][j]   = d[k+3][j]   - di30*dj00 - di31*dj10 - di32*dj20 - di33*dj30;  
             d[k+3][j+1] = d[k+3][j+1] - di30*dj01 - di31*dj11 - di32*dj21 - di33*dj31;
             d[k+3][j+2] = d[k+3][j+2] - di30*dj02 - di31*dj12 - di32*dj22 - di33*dj32;
             d[k+3][j+3] = d[k+3][j+3] - di30*dj03 - di31*dj13 - di32*dj23 - di33*dj33;
          }

          //remaining columns
          for (j=nb; j<n; j++)
          {
            dj00 = d[i][j]; dj10 = d[i+1][j]; dj20 = d[i+2][j]; dj30 = d[i+3][j]; 

             d[k][j]   = d[k][j]   - di00*dj00 - di01*dj10 - di02*dj20 - di03*dj30; 
             d[k+1][j] = d[k+1][j] - di10*dj00 - di11*dj10 - di12*dj20 - di13*dj30;
             d[k+2][j] = d[k+2][j] - di20*dj00 - di21*dj10 - di22*dj20 - di23*dj30;
             d[k+3][j] = d[k+3][j] - di30*dj00 - di31*dj10 - di32*dj20 - di33*dj30; 
          }
       }

       //remaining rows
       for (k=nb; k<n; k++)
       {
         di00 = d[k][i];   di01 = d[k][i+1];   di02 = d[k][i+2];   di03 = d[k][i+3];

          for (j=i+4; j<nb; j+=4)
          {
             dj00 = d[i][j];   dj01 = d[i][j+1];   dj02 = d[i][j+2];   dj03 = d[i][j+3];
             dj10 = d[i+1][j]; dj11 = d[i+1][j+1]; dj12 = d[i+1][j+2]; dj13 = d[i+1][j+3];
             dj20 = d[i+2][j]; dj21 = d[i+2][j+1]; dj22 = d[i+2][j+2]; dj23 = d[i+2][j+3];
             dj30 = d[i+3][j]; dj31 = d[i+3][j+1]; dj32 = d[i+3][j+2]; dj33 = d[i+3][j+3];

             d[k][j]   = d[k][j]   - di00*dj00 - di01*dj10 - di02*dj20 - di03*dj30;  
             d[k][j+1] = d[k][j+1] - di00*dj01 - di01*dj11 - di02*dj21 - di03*dj31;
             d[k][j+2] = d[k][j+2] - di00*dj02 - di01*dj12 - di02*dj22 - di03*dj32;
             d[k][j+3] = d[k][j+3] - di00*dj03 - di01*dj13 - di02*dj23 - di03*dj33;
          }

          //remaining columns
          for (j=nb; j<n; j++)
          {
             dj00 = d[i][j]; dj10 = d[i+1][j]; dj20 = d[i+2][j]; dj30 = d[i+3][j]; 

             d[k][j] = d[k][j] - di00*dj00 - di01*dj10 - di02*dj20 - di03*dj30; 
          }
       } 
    }      

   /***elimination for a small trailing submatrix < 4X4***/
   for (i=nb; i<n-1; i++)
    {
       //find and record k where |d(k,i)|=𝑚ax|d(j,i)|
        amax = d[i][i];
        indk = i;
        for (k=i+1; k<n; k++)
        {
            if (fabs(d[k][i]) > fabs(amax))
            {
                amax = d[k][i];
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
                c = d[i][j];
                d[i][j] = d[indk][j];
                d[indk][j] = c;
            }
        } 

        //store multiplier in place of d(k,i)
        for (k=i+1; k<n; k++)
        {
            d[k][i] = d[k][i]/d[i][i];
        }

        //subtract multiple of row a(i,:) to zero out d(k,i)
        for (k=i+1; k<n; k++)
        { 
            c = d[k][i]; 
            for (j=i+1; j<n; j++)
            {
                d[k][j] -= c*d[i][j];
            }
        }
    }
    gettimeofday(&end_time, 0);
 
    //print the running time
    seconds = end_time.tv_sec - start_time.tv_sec;
    microseconds = end_time.tv_usec - start_time.tv_usec;
    elapsed = seconds + 1e-6 * microseconds;
    printf("sequential calculation with loop unrolling time: %f\n\n",elapsed);

    printf("Starting comparison...\n\n"); 
    int cnt;
    cnt = test(a,d,n);
    if (cnt == 0)
        printf ("Done. There are no differences!\n");
    else
        printf ("Results are incorrect! The number of different elements is %d\n", cnt);

    return 0;
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

int test(double** t1, double** t2, int n)
{
    int i, j;
    int cnt;
    cnt = 0;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            if ((t1[i][j] - t2[i][j])*(t1[i][j] - t2[i][j]) > 1.0e-10)
            {
                cnt += 1;
            }
        }
    }

    return cnt;
}
