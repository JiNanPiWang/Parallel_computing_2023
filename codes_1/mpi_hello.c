/******************************************************************************
* FILE: mpi_hello.c
* DESCRIPTION:
*   A simple "Hello World" mpi program 
*      also print out processor name on which each process runs
* AUTHOR: Bing Bing Zhou
* LAST REVISED: 18/07/2022
******************************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
void main(int argc, char *argv[]){
   int nprocs, myid;
   int name_len;
   char processor_name[MPI_MAX_PROCESSOR_NAME];

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);

   MPI_Get_processor_name(processor_name, &name_len); //get processor name

   printf("Hello world from %s, rank %d/%d\n", processor_name, myid, nprocs);

   MPI_Finalize();
}
