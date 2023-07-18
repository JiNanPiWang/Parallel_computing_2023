/******************************************************************************
* FILE: omp_scan.c
* DESCRIPTION:  
*   OpenMp Program - prefix scan: s[i+1] += s[i] recursion starting from i=0
*   in place computation, i.e., output stored in the same input vector 
* AUTHOR: Bing Bing Zhou
* Last revised: 09/07/2022
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_vector(double* T, int cols);

int main (int argc, char *argv[]) 
{
   double* s; //input and output vector, i.e, inplace computation
   double* c; //store input data for correctness checking
   double* sw; //working vector of size nthreads
   int N; //vector size

   int	nthreads, i;
   double tt;

   if (argc == 3){
      N = atoi(argv[1]); 
      nthreads = atoi(argv[2]);

      printf("N = %d\n", N);
      printf("nthreads = %d\n\n", nthreads);
    }  
    else{
            printf("Usage: %s N nthreads\n\n"
                   " N: vector length\n"
                   " nthreads: the number of threads\n\n", argv[0]);
        return 1;
    }

   omp_set_num_threads(nthreads); //set the number of threads

  printf("initializing the vector...\n\n");
  // Allocate memory for the vector 
  s = (double*)malloc(N*sizeof(double));
  c = (double*)malloc(N*sizeof(double));
  sw = (double*)malloc(nthreads*sizeof(double)); 

  srand(time(0)); // Seed the random number generator
  /*** Initialize vector ***/
  for (i=0; i<N; i++)
  {
     s[i] = (double) rand() / RAND_MAX;
//     s[i] = (double) i + 1.0;
     c[i] = s[i];
  }

  for (i=0; i<nthreads; i++)
    sw[i] = 0.0;

//  printf ("vector s:\n");
//  print_vector(s, N);

  printf("Parallel computation starting...\n\n");
  tt = omp_get_wtime();
  #pragma omp parallel shared(s,nthreads) private(i) 
  {
    int tid = omp_get_thread_num();
    double sl = 0.0; //local scan

    //local scan
    #pragma omp for schedule(static)
    for (i=0; i<N; i++)
    {
      sl += s[i];
      s[i] = sl;
    }
    
    if (tid < nthreads-1)
      sw[tid+1] = sl; //store the last value sl

    //synchronize the threads
    // - very important to avoid data race condition!
    #pragma omp barrier 

    //seq scan on sw
    #pragma omp single 
    {
      for (i=1; i<nthreads; i++)
        sw[i] += sw[i-1];
    } //implicit barrier here!

    //scale elements by sw[tid] 
    sl = sw[tid];
    #pragma omp for schedule(static)
    for (i=0; i<N; i++)
      s[i] += sl;

  }
  tt = omp_get_wtime() - tt;
  printf("It takes %f seconds to finish the computation.\n\n", tt); 


//  printf ("vector s:\n");
//  print_vector(s, N);

  printf("Sequential computation starting...\n\n");  
  tt = omp_get_wtime();
  /*** sequential computation***/
  double st = 0.0;  
  for (i=0; i<N; i++)
  {
    st += c[i];
    c[i] = st;
  }
  tt = omp_get_wtime() - tt;
  printf("It takes %f seconds to finish the sequential computation.\n\n", tt); 


  printf("Checking the result...\n\n");
  /*** check the correctness ***/
  int cnt = 0;
  for (i=0; i<N; i++)
    if ((c[i]-s[i])*(c[i]-s[i])>1.0E-10) cnt++;
 
  if (cnt==0)
    printf ("Done.\n");
  else
    printf("results are incorrect! the number of different elements is %d\n", cnt);

}

void print_vector(double* T, int cols){
    for (int i=0; i < cols; i++)
       printf("%.2f  ", T[i]);
    printf("\n\n");
    return;
}


