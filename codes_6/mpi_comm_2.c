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

typedef struct
{
	int val;
	int nums;
} node;

int main(int argc, char **argv)
{
	int myid, numprocs;

	node *l;
	l = (node *) malloc(100 * sizeof(node));

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

	node recv_data;
	MPI_Send(&((node) {myid, myid + 1}), sizeof(node), MPI_BYTE, (myid + numprocs - 1) % numprocs, 0, MPI_COMM_WORLD);
	MPI_Recv(&recv_data, sizeof(node), MPI_BYTE, (myid + 1) % numprocs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	MPI_Gather(&recv_data, sizeof(node), MPI_BYTE, l + myid, sizeof(node), MPI_BYTE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
//	printf("NO.%d: left: %d right: %d | l[i]val: %d, l[i]num: %d\n",
//		   myid, (myid + numprocs - 1) % numprocs, (myid + 1) % numprocs, l[myid].val, l[myid].nums);

	if (myid == 0)
	{
		for (int i = 0; i < numprocs; ++i)
		{
			printf("No.%d process receives ID:", i);
			for (int j = 0; j < l[i].nums; ++j)
			{
				printf(" %d", l[i].val);
			}
			printf("\n");
		}
	}

	MPI_Finalize();

	return 0;
}


