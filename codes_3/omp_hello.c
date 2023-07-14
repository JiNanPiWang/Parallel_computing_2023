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

int main(int argc, char *argv[])
{
	int j = 0;
	#pragma omp parallel for
	for (int i = 0; i < 10; ++i)
	{
		j += 2;
		printf("    %d: ", omp_get_thread_num());
//		int j = 2 * (i + 1);
		printf("%d ", j);
	}
//  上方代码输出异常，输出的数都很大，如12 12 16 16 20等
//  原因是变量j是在并行区域的外部定义的，它被所有线程共享。这可能导致多个线程同时读取和更新j的值，从而产生竞争条件。
//  如多个线程都先执行了+2操作，最后再输出

	printf("\n\n\n");
	#pragma omp parallel for
	for (int i = 0; i < 10; ++i)
	{
		printf("    %d: ", omp_get_thread_num());
		int j = 2 * (i + 1);
		printf("%d ", j);
	}
}
