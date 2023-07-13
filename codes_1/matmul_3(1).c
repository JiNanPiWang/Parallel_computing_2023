/*************************************************************************************
* FILE: matmul_3.c
* DESCRIPTION:  
*   This serial program compute Matrix Multiply C = C + AB
*   It dynamically allocates contiguous memory space for 2D matrices,
*   uses ikj version to contiguously access matrices to reduce cache miss rate, and
*   unroll all i, j and k loops with unrolling factor = 4
* AUTHOR: Bing Bing Zhou
* Last revised: 12/03/2023
*************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

void print_matrix(double** T, int rows, int cols);

int main (int argc, char *argv[])
{
   double* a0; //auxiliary 1D for 2D matrix a
   double* b0; //auxiliary 1D for 2D matrix b
   double* c0; //auxiliary 1D for 2D matrix c
   double* c10; //auxiliary 1D for 2D matrix c1
   double** a; //the two-dimensional input matrix
   double** b; //the two-dimensional input matrix
   double** c; //the two-dimensional input/output matrix
   double** c1; //the two-dimensional input/output matrix

   int NRA, NCA, NCB; //matrices lengths
   int NRA0, NCA0, NCB0; //divisible by 4 for loop unrolling

   double ai00, ai01, ai02, ai03;
   double ai10, ai11, ai12, ai13;
   double ai20, ai21, ai22, ai23;
   double ai30, ai31, ai32, ai33;
   double bj00, bj01, bj02, bj03;
   double bj10, bj11, bj12, bj13;
   double bj20, bj21, bj22, bj23;
   double bj30, bj31, bj32, bj33;

   int i, j, k;

   struct timeval start_time, end_time;
   long seconds, microseconds;
   double elapsed;

   if (argc == 4){
      NRA = atoi(argv[1]);
      NCA = atoi(argv[2]);
      NCB = atoi(argv[3]);

      printf("NRA = %d, NCA = %d, NCB = %d\n", NRA, NCA, NCB);
    }
    else{
            printf("Usage: %s NRA NCA NCB \n\n"
                   " NRA: matrix a row length\n"
                   " NCA: matrix a column (or b row) length\n"
                   " NCB:  matrix b column length\n\n", argv[0]);
        return 1;
    }

    NRA0 = NRA/4 * 4;
    NCA0 = NCA/4 * 4;
    NCB0 = NCB/4 * 4;

    // Allocate contiguous memory for 2D matrices
    a0 = (double*)malloc(NRA*NCA*sizeof(double));
    a = (double**)malloc(NRA*sizeof(double*));
    for (int i=0; i<NRA; i++){
        a[i] = a0 + i*NCA;
    }

    b0 = (double*)malloc(NCA*NCB*sizeof(double));
    b = (double**)malloc(NCA*sizeof(double*));
    for (int i=0; i<NCA; i++){
        b[i] = &(b0[i*NCB]);
    }

    c0 = (double*)malloc(NRA*NCB*sizeof(double));
    c = (double**)malloc(NRA*sizeof(double*));
    for (int i=0; i<NRA; i++){
        c[i] = &(c0[i*NCB]);
    }

    c10 = (double*)malloc(NRA*NCB*sizeof(double));
    c1 = (double**)malloc(NRA*sizeof(double*));
    for (int i=0; i<NRA; i++){
        c1[i] = &(c10[i*NCB]);
    }

   printf("Initializing matrices...\n\n");
   srand(time(0)); // Seed the random number generator
//   srand(0); 
   for (i=0; i<NRA; i++)
      for (j=0; j<NCA; j++)
         a[i][j]= (double) rand() / RAND_MAX;
//         a[i][j]= (double) i + 1.0;

   for (i=0; i<NCA; i++)
      for (j=0; j<NCB; j++)
         b[i][j]= (double) rand() / RAND_MAX;
//         b[i][j]= (double) j;

   for (i=0; i<NRA; i++)
      for (j=0; j<NCB; j++)
         c[i][j]= 0.0;

   for (i=0; i<NRA; i++)
      for (j=0; j<NCB; j++)
         c1[i][j]= 0.0;

//  printf ("matrix a:\n");
//  print_matrix(a, NRA, NCA);
//  printf ("matrix b:\n");
//  print_matrix(b, NCA, NCB);
//  printf ("matrix C:\n");
//  print_matrix(c, NRA, NCB);

  /*** start matrix multiplication (ikj version) ***/
   printf("Starting matrix multiplication - ikj version\n\n");
   gettimeofday(&start_time, 0);
   for (i=0; i<NRA; i++)
      for (k=0; k<NCA; k++)
      {
         ai00 = a[i][k];
         for (j=0; j<NCB; j++)
            c1[i][j] +=  ai00 * b[k][j];
//            c1[i][j] +=  a[i][k] * b[k][j];
      }
   gettimeofday(&end_time, 0);
   seconds = end_time.tv_sec - start_time.tv_sec;
   microseconds = end_time.tv_usec - start_time.tv_usec;
   elapsed = seconds + 1e-6 * microseconds;
   printf("ikj version takes %f seconds to finish the computation.\n\n", elapsed);


  /*** start matrix multiplication (ikj version with loop unrolling) ***/
  printf("Starting matrix multiplication - ikj version with loop unrolling\n\n");
  gettimeofday(&start_time, 0);
  //unroll all i, j and k loops with unrolling factor = 4
  for (i=0; i<NRA0; i+=4)
  {
    for (k=0; k<NCA0; k+=4)
    {
      ai00 = a[i][k];   ai01 = a[i][k+1];   ai02 = a[i][k+2];   ai03 = a[i][k+3];
      ai10 = a[i+1][k]; ai11 = a[i+1][k+1]; ai12 = a[i+1][k+2]; ai13 = a[i+1][k+3];
      ai20 = a[i+2][k]; ai21 = a[i+2][k+1]; ai22 = a[i+2][k+2]; ai23 = a[i+2][k+3];
      ai30 = a[i+3][k]; ai31 = a[i+3][k+1]; ai32 = a[i+3][k+2]; ai33 = a[i+3][k+3];

      for(j=0; j<NCB0; j+=4)
      {
        bj00 = b[k][j];   bj01 = b[k][j+1];   bj02 = b[k][j+2];   bj03 = b[k][j+3];
        bj10 = b[k+1][j]; bj11 = b[k+1][j+1]; bj12 = b[k+1][j+2]; bj13 = b[k+1][j+3];
        bj20 = b[k+2][j]; bj21 = b[k+2][j+1]; bj22 = b[k+2][j+2]; bj23 = b[k+2][j+3];
        bj30 = b[k+3][j]; bj31 = b[k+3][j+1]; bj32 = b[k+3][j+2]; bj33 = b[k+3][j+3];

        c[i][j]   = c[i][j]   + ai00*bj00 + ai01*bj10 + ai02*bj20 + ai03*bj30;
        c[i][j+1] = c[i][j+1] + ai00*bj01 + ai01*bj11 + ai02*bj21 + ai03*bj31;
        c[i][j+2] = c[i][j+2] + ai00*bj02 + ai01*bj12 + ai02*bj22 + ai03*bj32;
        c[i][j+3] = c[i][j+3] + ai00*bj03 + ai01*bj13 + ai02*bj23 + ai03*bj33;

        c[i+1][j]   = c[i+1][j]   + ai10*bj00 + ai11*bj10 + ai12*bj20 + ai13*bj30;
        c[i+1][j+1] = c[i+1][j+1] + ai10*bj01 + ai11*bj11 + ai12*bj21 + ai13*bj31;
        c[i+1][j+2] = c[i+1][j+2] + ai10*bj02 + ai11*bj12 + ai12*bj22 + ai13*bj32;
        c[i+1][j+3] = c[i+1][j+3] + ai10*bj03 + ai11*bj13 + ai12*bj23 + ai13*bj33;

        c[i+2][j]   = c[i+2][j]   + ai20*bj00 + ai21*bj10 + ai22*bj20 + ai23*bj30;
        c[i+2][j+1] = c[i+2][j+1] + ai20*bj01 + ai21*bj11 + ai22*bj21 + ai23*bj31;
        c[i+2][j+2] = c[i+2][j+2] + ai20*bj02 + ai21*bj12 + ai22*bj22 + ai23*bj32;
        c[i+2][j+3] = c[i+2][j+3] + ai20*bj03 + ai21*bj13 + ai22*bj23 + ai23*bj33;

        c[i+3][j]   = c[i+3][j]   + ai30*bj00 + ai31*bj10 + ai32*bj20 + ai33*bj30;
        c[i+3][j+1] = c[i+3][j+1] + ai30*bj01 + ai31*bj11 + ai32*bj21 + ai33*bj31;
        c[i+3][j+2] = c[i+3][j+2] + ai30*bj02 + ai31*bj12 + ai32*bj22 + ai33*bj32;
        c[i+3][j+3] = c[i+3][j+3] + ai30*bj03 + ai31*bj13 + ai32*bj23 + ai33*bj33;

      }

     //For elelments in remaining j columns
      for (j=NCB0; j<NCB; j++)
      {
         bj00 = b[k][j]; bj10 = b[k+1][j]; bj20 = b[k+2][j]; bj30 = b[k+3][j];
         c[i][j]   = c[i][j]   + ai00*bj00 + ai01*bj10 + ai02*bj20 + ai03*bj30;
         c[i+1][j] = c[i+1][j] + ai10*bj00 + ai11*bj10 + ai12*bj20 + ai13*bj30;
         c[i+2][j] = c[i+2][j] + ai20*bj00 + ai21*bj10 + ai22*bj20 + ai23*bj30;
         c[i+3][j] = c[i+3][j] + ai30*bj00 + ai31*bj10 + ai32*bj20 + ai33*bj30;
      }
    }

    // for the remaining k
    for (k=NCA0; k<NCA; k++)
    {
      ai00 = a[i][k];
      ai10 = a[i+1][k];
      ai20 = a[i+2][k];
      ai30 = a[i+3][k];

      for(j=0; j<NCB0; j+=4)
      {
        bj00 = b[k][j];
        bj01 = b[k][j+1];
        bj02 = b[k][j+2];
        bj03 = b[k][j+3];

        c[i][j]   += ai00*bj00;
        c[i][j+1] += ai00*bj01;
        c[i][j+2] += ai00*bj02;
        c[i][j+3] += ai00*bj03;

        c[i+1][j]   += ai10*bj00;
        c[i+1][j+1] += ai10*bj01;
        c[i+1][j+2] += ai10*bj02;
        c[i+1][j+3] += ai10*bj03;

        c[i+2][j]   += ai20*bj00;
        c[i+2][j+1] += ai20*bj01;
        c[i+2][j+2] += ai20*bj02;
        c[i+2][j+3] += ai20*bj03;

        c[i+3][j]   += ai30*bj00;
        c[i+3][j+1] += ai30*bj01;
        c[i+3][j+2] += ai30*bj02;
        c[i+3][j+3] += ai30*bj03;
      }

      //For elelments in remaining j columns
      for (j=NCB0; j<NCB; j++)
      {
        bj00 = b[k][j];
        c[i][j]   += ai00 * bj00;
        c[i+1][j] += ai10 * bj00;
        c[i+2][j] += ai20 * bj00;
        c[i+3][j] += ai30 * bj00;
      }
    }
  }

  //For elements in remaining i rows
  for (i=NRA0; i<NRA; i++)
  {
    for (k=0; k<NCA0; k+=4)
    {
      ai00 = a[i][k]; ai01 = a[i][k+1]; ai02 = a[i][k+2]; ai03 = a[i][k+3];
      for (j=0; j<NCB0; j+=4)
      {
        bj00 = b[k][j];   bj01 = b[k][j+1];   bj02 = b[k][j+2];   bj03 = b[k][j+3];
        bj10 = b[k+1][j]; bj11 = b[k+1][j+1]; bj12 = b[k+1][j+2]; bj13 = b[k+1][j+3];
        bj20 = b[k+2][j]; bj21 = b[k+2][j+1]; bj22 = b[k+2][j+2]; bj23 = b[k+2][j+3];
        bj30 = b[k+3][j]; bj31 = b[k+3][j+1]; bj32 = b[k+3][j+2]; bj33 = b[k+3][j+3];

        c[i][j]   = c[i][j]   + ai00*bj00 + ai01*bj10 + ai02*bj20 + ai03*bj30;
        c[i][j+1] = c[i][j+1] + ai00*bj01 + ai01*bj11 + ai02*bj21 + ai03*bj31;
        c[i][j+2] = c[i][j+2] + ai00*bj02 + ai01*bj12 + ai02*bj22 + ai03*bj32;
        c[i][j+3] = c[i][j+3] + ai00*bj03 + ai01*bj13 + ai02*bj23 + ai03*bj33;
      }

     //For elelments in remaining j columns
      for (j=NCB0; j<NCB; j++)
      {
         bj00 = b[k][j]; bj10 = b[k+1][j]; bj20 = b[k+2][j]; bj30 = b[k+3][j];
         c[i][j]   = c[i][j]   + ai00*bj00 + ai01*bj10 + ai02*bj20 + ai03*bj30;
      }
   }

   // for the remaining k
    for (k=NCA0; k<NCA; k++)
    {
      ai00 = a[i][k];
     for (j=0; j<NCB; j++)
        c[i][j]   += ai00 * b[k][j];
    }
  }

   gettimeofday(&end_time, 0);
   seconds = end_time.tv_sec - start_time.tv_sec;
   microseconds = end_time.tv_usec - start_time.tv_usec;
   elapsed = seconds + 1e-6 * microseconds;
   printf("sequential calculation with loop unrolling time: %f\n\n",elapsed);

   //print the result matrix c
//   printf ("matrix C:\n");
//   print_matrix(c, NRA, NCB);

/*** check the correctness ***/
   printf("Starting comparison\n\n");
   int cnt= 0;
   for (i=0; i<NRA; i++)
      for (j=0; j<NCB; j++)
         if((c[i][j] - c1[i][j])*(c[i][j] - c1[i][j])>1.0E-10) cnt++;
   if (cnt == 0)
      printf ("Done. There are no differences!\n");
   else
      printf ("Results are incorrect! The number of different elements is %d\n", cnt);
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