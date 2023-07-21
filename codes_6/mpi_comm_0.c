#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	int rank, size;
	int my_number, other_number, sum;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (size != 2)
	{
		fprintf(stderr, "This program requires exactly 2 processes.\n");
		MPI_Finalize();
		return 1;
	}

	// 假设第一个进程(rank=0)持有整数1，第二个进程(rank=1)持有整数2
	if (rank == 0)
	{
		my_number = 1;
	}
	else
	{
		my_number = 2;
	}

	// 使用MPI_Send和MPI_Recv进行通信
	if (rank == 0)
	{
//		MPI_Send函数：
//		函数原型：int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
//		功能：将数据从一个进程发送到另一个进程。
//		参数：
//		buf：指向要发送数据的缓冲区的指针。
//		count：发送的数据元素数量。
//		datatype：发送数据的类型。在这里，我们使用MPI_INT，表示发送整数类型的数据。
//		dest：目标进程的rank（进程号），表示数据将发送给哪个进程。
//		tag：标签，用于标识发送的消息，可以为任意整数。
//		comm：通信域（communicator），表示这个通信的范围。在这里，我们使用MPI_COMM_WORLD，表示全局通信域，包含了所有的进程。

//		MPI_Recv函数：
//		函数原型：int
//		MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
//		功能：接收来自另一个进程的数据。
//		参数：
//		buf：指向接收数据的缓冲区的指针。
//		count：接收的数据元素数量。
//		datatype：接收数据的类型。在这里，我们使用MPI_INT，表示接收整数类型的数据。
//		source：源进程的rank（进程号），表示从哪个进程接收数据。
//		tag：标签，用于标识接收的消息，必须和发送时的tag匹配，可以为任意整数。
//		comm：通信域（communicator），表示这个通信的范围。在这里，我们使用MPI_COMM_WORLD，表示全局通信域，包含了所有的进程。
//		status：指向MPI_Status结构体的指针，用于获取接收操作的状态信息。在这里，我们使用MPI_STATUS_IGNORE，表示忽略状态信息。

		MPI_Send(&my_number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&other_number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else
	{
		MPI_Recv(&other_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(&my_number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	// 对接收到的数据进行相加
	sum = my_number + other_number;

	// 输出结果
	if (rank == 0)
	{
		printf("Process 0: Received %d from Process 1. Sum = %d\n", other_number, sum);
	}
	else
	{
		printf("Process 1: Received %d from Process 0. Sum = %d\n", other_number, sum);
	}

	MPI_Finalize();
	return 0;
}
