/******************************************************************************
* FILE: omp_hello.c
* DESCRIPTION:
*   A simple "Hello World" openmp program for students to modify
* AUTHOR: Bing Bing Zhou
* LAST REVISED: 02/07/2022
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
   /* Fork threads, each print out a "Hello World" statement */
   #pragma omp parallel
   {
      int tid; // local variable

      tid = omp_get_thread_num(); // get thread ID
      printf("Hello World from thread = %d\n", tid);

   }  /* end of parallel section */

}
