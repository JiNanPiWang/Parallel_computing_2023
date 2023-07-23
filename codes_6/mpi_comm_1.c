/******************************************************************************
* FILE: mpi_comm_o.c
* DESCRIPTION:
*   A mpi program for passing a token around the processes:
*     Processes are organized as a 1D ring. At the beginning process 0 sets
*     the token (integer) to 1, and then passes the token to process numprocs-1.
*     When process numprocs-1 receives the token, it passes it to its immediate
*     left neighbour, and so on as follows (numprocs = 5):
*       process 0 sends the token to process 4
*       process 4 sends the token to process 3
*       process 3 sends the token to process 2
*       process 2 sends the token to process 1
*       process 1 sends the token bake to process 0
* AUTHOR: Bing Bing Zhou
* LAST REVISED: 18/07/2022
******************************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	MPI_Status status;

	int myid, numprocs;
	int token = 0; //token initialized to 0

	int *l;
	l = (int *) malloc(100 * sizeof(int));

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (numprocs == 1)
	{
		//only one process, no communication
		printf("I'm the lonely process with id = %d :(\n", myid);
		MPI_Finalize();
		return 0;
	}
	else
	{
		if (myid == 0)
			printf("Processes num: %d\n", numprocs);
	}

	int recv_id;
	MPI_Send(&myid, 1, MPI_INT, (myid + numprocs - 1) % numprocs, 0, MPI_COMM_WORLD);
	MPI_Recv(&recv_id, 1, MPI_INT, (myid + 1) % numprocs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(&recv_id, 1, MPI_INT, (myid + numprocs - 1) % numprocs, 0, MPI_COMM_WORLD);

	MPI_Gather(&recv_id, 1, MPI_INT, l + myid, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
//	printf("NO.%d: left: %d right: %d | l[i]: %d\n", myid, (myid + numprocs - 1) % numprocs, (myid + 1) % numprocs, l[myid]);

	if (myid == 0)
	{
		for (int i = 0; i < numprocs; ++i)
		{
			printf("No.%d process receives ID: %d\n", i, l[i]);
		}
	}

	MPI_Finalize();

	return 0;
}


