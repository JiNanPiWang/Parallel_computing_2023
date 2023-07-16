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

void scan_openmp(double *a, int a_len, double *result, double *time)
{
	double start_time = omp_get_wtime();
	double *c, *W;
	c = (double *) malloc(a_len * sizeof(double));
	W = (double *) malloc(a_len * sizeof(double));

	#pragma omp parallel
	{
		// STAGE 1 存储分段前缀和
		int is_first = 1; // 判断当前循环是否是开始
		#pragma omp for
		for (int i = 0; i < a_len; ++i)
		{
			if (is_first) // 判断是否为分段的第一个，如果是，则为第一个 c0(c[0]) = a0
				c[i] = a[i];
			else // 如果不是，则为c0−1(c[1]) = c0 ⊕ a1, c0−2(c[2]) = c0−1 ⊕ a2
				c[i] = c[i - 1] + a[i];

			// 还需要判断是不是最后一个

			is_first = 0;
			printf("%d %.2f\n", omp_get_thread_num(), c[i]);
		}

		// STAGE 2 得到分段前缀和的前缀和，存储于W中
		#pragma omp barrier
		#pragma omp single
		{
			// 首先获取W数组 𝑊:[0, 𝑐0−2, 𝑐3−5]，需要使用c数组
			//     然后    𝑊:[0, 𝑐0−2, 𝑐0−5]，到边界时，W[i] = W[i - 1] + c[i]
			double last_Wi = 0; // 用于记录W[i - 1]
			W[0] = 0.0; // 初始化W
			for (int i = 1; i < a_len; ++i)
			{
				W[i] = 0.0; // 初始化W
				if (i == (a_len - 1) || c[i + 1] == a[i + 1]) // 到达最末端或中间的边界
				{
//					printf("%d ", i);
					W[i] = last_Wi + c[i];
					last_Wi = W[i];
				}
			}
		}

		// STAGE 3


//		for (int i = 0; i < 100; ++i)
//		{
////          检测效果
//			printf("%d %d %d\n", omp_get_thread_num(), i, is_first);
//			is_first = 1;
//		}

	}

			printf("W: ");
			print_array(W, a_len);
	printf("c: ");
	print_array(c, a_len);
	printf("a: ");
	print_array(a, a_len);

	double end_time = omp_get_wtime();
	*time = end_time - start_time;
}

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

	A = (double *) malloc(A_len * sizeof(double));
	result_normal = (double *) malloc(A_len * sizeof(double));
	result_openmp = (double *) malloc(A_len * sizeof(double));

	init_random_array(A, A_len);

	scan_normal(A, A_len, result_normal, &elapsed_time_normal);
	scan_openmp(A, A_len, result_openmp, &elapsed_time_openmp);

//	print_array(A, A_len);
//	print_array(result_normal, result_len);

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
