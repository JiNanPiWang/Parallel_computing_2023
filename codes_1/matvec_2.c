/******************************************************************************
* FILE: matvec_2.c
* DESCRIPTION:  
*   A sequential program for Matrix-vector Multiply b = Ax with loop unrolling
*   unroll both i and k loops with unrolling factor = 4
* AUTHOR: Bing Bing Zhou
* Last revised: 05/07/2022
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

void print_matrix(double** T, int rows, int cols);
void print_vector(double* T, int cols);

int main (int argc, char *argv[]) 
{
   double* a0; //auxiliary 1D array to make a contiguously allocated
   double** a; //the two-dimensional input matrix
   double* x; //input vector
   double* b; //the resulting vector

   int NRA, NCA, NRA0, NCA0; //matrix size and NRA0/NCA0 divisible by unrolling factor

   int i, k;
   struct timeval start_time, end_time;
   long seconds, microseconds;
   double elapsed;

   if (argc == 3){
      NRA = atoi(argv[1]); 
      NCA = atoi(argv[2]);

      printf("NRA = %d, NCA = %d\n", NRA, NCA);
    }  
    else{
            printf("Usage: %s NRA NCA\n\n"
                   " NRA: matrix a row length\n"
                   " NCA: matrix a column (or x) length\n\n", argv[0]);
        return 1;
    }

   // Allocate contiguous memory for 2D matrices
   a0 = (double*)malloc(NRA*NCA*sizeof(double));
   a = (double**)malloc(NRA*sizeof(double*));
   for (int i=0; i<NRA; i++){
      a[i] = a0 + i*NCA;
   }

   //Allocate memory for vectors      
   x = (double*)malloc(NCA*sizeof(double));
   b = (double*)malloc(NRA*sizeof(double));

  printf("Initializing matrix and vectors\n\n");
  srand(time(0)); // Seed the random number generator
  /*** Initialize matrix and vectors ***/
  for (i=0; i<NRA; i++)
    for (k=0; k<NCA; k++)
      a[i][k]= (double) rand() / RAND_MAX;

  for (i=0; i<NCA; i++)
      x[i] = (double) rand() / RAND_MAX;

  for (i=0; i<NRA; i++)
      b[i]= 0.0;

/* 
  printf ("matrix a:\n");
  print_matrix(a, NRA, NCA);
  printf ("vector x:\n");
  print_vector(x, NCA);
  printf ("vector b:\n");
  print_vector(b, NRA);
*/

  NRA0 = NRA/4 * 4; //unrolling factor is set to 4
  NCA0 = NCA/4 * 4; //unrolling factor is set to 4


  printf("Starting matrix-vector multiplication with loop unrolling\n\n");
  gettimeofday(&start_time, 0);
  double b0, b1, b2, b3;
  double x0, x1, x2, x3;
  for (i=0; i<NRA0; i+=4)
  {
    b0 = b[i]; b1 = b[i+1]; b2 = b[i+2]; b3 = b[i+3];
    for (k=0; k<NCA0; k+=4)
    {
      x0 = x[k]; x1 = x[k+1]; x2 = x[k+2]; x3 = x[k+3];
      b0 = b0 + a[i][k]*x0   + a[i][k+1]*x1   + a[i][k+2]*x2   + a[i][k+3]*x3;
      b1 = b1 + a[i+1][k]*x0 + a[i+1][k+1]*x1 + a[i+1][k+2]*x2 + a[i+1][k+3]*x3;
      b2 = b2 + a[i+2][k]*x0 + a[i+2][k+1]*x1 + a[i+2][k+2]*x2 + a[i+2][k+3]*x3;
      b3 = b3 + a[i+3][k]*x0 + a[i+3][k+1]*x1 + a[i+3][k+2]*x2 + a[i+3][k+3]*x3;
    }

    /*** for remaining k ***/
    for (k=NCA0; k<NCA; k++)
    {
      x0 = x[k];   
      b0 +=  a[i][k]   * x0;
      b1 +=  a[i+1][k] * x0;
      b2 +=  a[i+2][k] * x0; 
      b3 +=  a[i+3][k] * x0;
    }

    b[i] = b0; b[i+1] = b1; b[i+2] = b2; b[i+3] = b3;
  }

  /*** Compute remaining b[i] from NRA0 to NRA-1 ***/
  for (i=NRA0; i<NRA; i++){
    b0 = b[i];
    for (k=0; k<NCA; k++)   
      b0 +=  a[i][k] * x[k];
    b[i] = b0;
  }
  gettimeofday(&end_time, 0);
  seconds = end_time.tv_sec - start_time.tv_sec;
  microseconds = end_time.tv_usec - start_time.tv_usec;
  elapsed = seconds + 1e-6 * microseconds;
  printf("The computation takes %f seconds to complete.\n\n", elapsed);  

  /*** Check the correctness ***/
  /*** compare the result with one without unrolling ***/
  printf("Checking correctness\n\n");
  double cc, er;
  int cnt = 0; //number of elements which are different
  for (i=0; i<NRA; i++)
  {
    cc = 0.0;
    for (k=0; k<NCA; k++)
      cc += a[i][k] * x[k];
    er = (b[i] - cc) * (b[i] - cc);
    if (er > 1.0E-10)
    {
      cnt += 1;
    }
  }
  
/*** Print results ***/

// printf("******************************************************\n");
// printf("Resulting vector:\n");
// print_vector(b, NRA);
// printf("******************************************************\n");

  if (cnt == 0)
    printf ("Done.\n");
  else
    printf ("Results are incorrect! The number of different elements is %d\n", cnt);
  
/*** Print results ***/
// printf("******************************************************\n");
// printf("Resulting vector:\n");
// print_vector(b, NRA);
// printf("******************************************************\n");

}

void print_matrix(double** T, int rows, int cols){
    for (int i=0; i < rows; i++){
        for (int j=0; j < cols; j++)
            printf("%.2f  ", T[i][j]);
        printf("\n");
    }
    printf("\n\n");
    return;
}

void print_vector(double* T, int cols){
    for (int i=0; i < cols; i++)
       printf("%.2f  ", T[i]);
    printf("\n\n");
    return;
}


