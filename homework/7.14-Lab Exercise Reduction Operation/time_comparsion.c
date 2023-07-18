#include <stdio.h>
#include <omp.h>

int arr[100000000];

int main()
{
	int n = 100000000; // 数字的数量
	int sum1 = 0; // 使用reduction语句计算的和
	int sum2 = 0; // 手动累加部分和计算的和

	// 初始化数组的值
	for (int i = 0; i < n; i++)
	{
		arr[i] = i + 1;
	}

	double start_time1 = omp_get_wtime(); // 使用reduction语句方法的开始时间
	// 使用OpenMP reduction语句计算和
	#pragma omp parallel for reduction(+: sum1)
	for (int i = 0; i < n; i++)
	{
		sum1 += arr[i];
	}
	double end_time1 = omp_get_wtime(); // 使用reduction语句方法的结束时间

	double start_time2 = omp_get_wtime(); // 手动部分和方法的开始时间
	// 通过手动累加部分和计算和
	for (int i = 0; i < n; i++)
	{
		sum2 += arr[i];
	}
	double end_time2 = omp_get_wtime(); // 手动部分和方法的结束时间

	printf("Use openmp reduction: %f s\n", end_time1 - start_time1);
	printf("Use normal way: %f s\n", end_time2 - start_time2);

	return 0;
}
