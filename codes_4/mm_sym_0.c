/******************************************************************************
* FILE: mm_sym_0.c
* DESCRIPTION:
* A simple sequential program for Computing C = A*A^T
* for students to modify
* to compile: gcc -fopenmp -O3 -o mm_sym_0 mm_sym_0.c
* AUTHOR: Bing Bing Zhou
* LAST REVISED: 18/06/2022
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_matrix(double **T, int rows, int cols);

int main(int argc, char *argv[])
{
	double *a0; //auxiliary 1D for 2D matrix a
	double *c0; //auxiliary 1D for 2D matrix c
	double **a; //the two-dimensional input matrix
	double **c; //the resulting matrix

	int NRA, NCA; //matrix size

	int nthreads, i, j, k;
	double ct;

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
		       " NCA: matrix a column length\n\n", argv[0]);
		return 1;
	}

	printf("Creating and initializing matrices...\n\n");
	// Allocate contiguous memory for 2D matrices
	a0 = (double *) malloc(NRA * NCA * sizeof(double));
	a = (double **) malloc(NRA * sizeof(double *));
	for (int i = 0; i < NRA; i++)
	{
		a[i] = a0 + i * NCA;
	}

	c0 = (double *) malloc(NRA * NRA * sizeof(double));
	c = (double **) malloc(NRA * sizeof(double *));
	for (int i = 0; i < NRA; i++)
	{
		c[i] = &(c0[i * NRA]);
	}

	/*** Initialize matrices ***/
	srand(time(0)); // Seed the random number generator
	for (i = 0; i < NRA; i++)
	{
		for (j = 0; j < NCA; j++)
		{
			a[i][j] = (double) rand() / RAND_MAX;
		}
	}

	for (i = 0; i < NRA; i++)
	{
		for (j = 0; j < NRA; j++)
		{
			c[i][j] = 0.0;
		}
	}
/* 
   printf ("matrix a:\n");
   print_matrix(a, NRA, NCA);
   printf ("matrix C:\n");
   print_matrix(c, NRA, NRA);
*/
	printf("Starting sequential matrix multiplication C = C + A*A^T...\n\n");
	double tt = omp_get_wtime();
	for (i = 0; i < NRA; i++)
	{
		for (j = i; j < NRA; j++)
		{
			ct = c[i][j];
			for (k = 0; k < NCA; k++)
			{
				ct += a[i][k] * a[j][k];
			}
			c[i][j] = ct;
		}
	}
	tt = omp_get_wtime() - tt;
	printf("It takes %f seconds to finish the computation\n\n", tt);


/*** Print results ***/
/*  printf("Result Matrix:\n");
  print_matrix(c, NRA, NRA);
*/

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

