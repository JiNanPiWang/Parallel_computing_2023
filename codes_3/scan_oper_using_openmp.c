//
// Created by Tianyi Zhang on 2023/7/15.
//
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

// parseArgument用于处理输入参数，并判断是否有效
int parseArgument(const char *arg, int *value);

// 随机初始化数组
void init_random_array(double *array, int array_len);

// 输出数组
void print_array(double *array, int array_len);

// scan_normal用于计算普通扫描，并保存运算用时
void scan_normal(double *array, int array_len, double *result, double *time);

// 用于比较两个数组是否相等
int compare_two_array(double *array_a, double *array_b, int array_len);


int main(int argc, char *argv[])
{
	double *A;
	double *result_normal; //the result (normal version)
	double *result_openmp; //the result (openmp version)
	int A_len, result_len;
	double elapsed_time_normal, elapsed_time_openmp;

	if (argc == 2)
	{
		if (parseArgument(argv[1], &A_len) == -1)
		{
			printf("Invalid value for A_len, should be an integer\n");
			return 2;
		}
		result_len = A_len;
		printf("A_len = %d\n\n", A_len);
	}
	else
	{
		printf("Incorrect number of arguments, valid input: A_len\n");
		return 2;
	}

	A = (double *) malloc(A_len * sizeof(int));
	result_normal = (double *) malloc(A_len * sizeof(int));
	result_openmp = (double *) malloc(A_len * sizeof(int));
	init_random_array(A, A_len);
	scan_normal(A, A_len, result_normal, &elapsed_time_normal);


	print_array(A, A_len);
	print_array(result_normal, result_len);

	return 0;
}


// parseArgument用于处理输入参数，并判断是否有效
int parseArgument(const char *arg, int *value) // 处理参数
{
	char *endptr; // 指向字符指针的参数endptr，用于指示转换过程中停止的位置。
	long result = strtol(arg, &endptr, 10);
	if (endptr == arg || errno == ERANGE || errno == EINVAL)
	{
//		如果endptr等于arg，表示没有成功转换任何字符，或者字符串为空。这意味着转换失败，返回-1作为错误码。
//      如果errno的值为ERANGE，表示转换结果超出了可表示的范围。这意味着转换失败，返回-1作为错误码。
//      如果errno的值为EINVAL，表示传递给strtol的参数无效。这意味着转换失败，返回-1作为错误码。
		return -1; // 转换错误，返回错误码
	}
	*value = (int) result;
	return 0; // 转换成功，返回成功码
}


// 输出数组
void print_array(double *array, int array_len)
{
	for (int i = 0; i < array_len; ++i)
	{
		printf("%.2f ", array[i]);
	}
	printf("\n");
}


// 随机初始化数组
void init_random_array(double *array, int array_len)
{
	unsigned int seed = (unsigned int) (time(NULL) ^ getpid());
	srand(seed);

	for (int i = 0; i < array_len; ++i)
	{
		array[i] = (double) rand() / RAND_MAX;;
	}
}


// scan_normal用于计算普通扫描，并保存运算用时
void scan_normal(double *array, int array_len, double *result, double *time)
{
	double start_time = omp_get_wtime();
	result[0] = array[0];
	for (int i = 1; i < array_len; ++i)
	{
		result[i] = result[i - 1] + array[i];
	}
	double end_time = omp_get_wtime();
	*time = end_time - start_time;
}


// 用于比较两个数组是否相等
int compare_two_array(double *array_a, double *array_b, int array_len)
{
	for (int i = 0; i < array_len; ++i)
	{
		if (fabs(array_a[i] - array_b[i]) > 1e-6)
		{
			return 0;
		}
	}
	return 1;
}
