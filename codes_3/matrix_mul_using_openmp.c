//
// Created by Tianyi Zhang on 2023/7/14.
//
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// parseArgument用于处理输入参数，并判断是否有效
int parseArgument(const char *arg, int *value);

// init_random_matrix用于生成随机矩阵
void init_random_matrix(double **matrix, int rows, int cols);

// init_zero_matrix用于生成0矩阵
void init_zero_matrix(double **matrix, int rows, int cols);

// print_matrix用于输出矩阵
void print_matrix(double **matrix, int rows, int cols);

// matrix_mul_normal用于计算普通矩阵乘法，并保存运算用时
void matrix_mul_normal(double **ma_a, int rows_a, int cols_a, double **ma_b, int cols_b, double **result, double *time);

// matrix_mul_normal用于计算普通矩阵乘法，并保存运算用时
void matrix_mul_openmp(double **ma_a, int rows_a, int cols_a, double **ma_b, int cols_b, double **result, double *time);

// 用于比较两个矩阵是否相等
int compare_two_matrix(double **ma_a, double **ma_b, int rows, int cols);

int main(int argc, char *argv[])
{
	// 计算矩阵A×B，使用openmp和普通方法，并进行对比
	int A_rows, A_cols, B_rows, B_cols, result_rows, result_cols; // B_rows = A_cols，不用输入
	double **A; //the two-dimensional input matrix
	double **B; //the two-dimensional input matrix
	double **result_normal; //the resulting matrix (normal version)
	double **result_openmp; //the resulting matrix (openmp version)
	double elapsed_time_normal, elapsed_time_openmp;

	// 处理输入参数
	if (argc == 4)
	{
		if (parseArgument(argv[1], &A_rows) == -1)
		{
			printf("Invalid value for A_rows, should be an integer\n");
			return 2;
		}
		if (parseArgument(argv[2], &A_cols) == -1)
		{
			printf("Invalid value for A_cols, should be an integer\n");
			return 2;
		}
		if (parseArgument(argv[3], &B_cols) == -1)
		{
			printf("Invalid value for B_cols, should be an integer\n");
			return 2;
		}
		B_rows = A_cols;
		result_rows = A_rows;
		result_cols = B_cols;
		printf("A_rows = %d, A_cols = %d, B_cols = %d\n\n", A_rows, A_cols, B_cols);
	}
	else
	{
		printf("Incorrect number of arguments, valid input: A_rows, A_cols, B_cols\n");
		return 2;
	}

	// 初始化
	A = (double **) malloc(A_rows * sizeof(double *));
	B = (double **) malloc(B_rows * sizeof(double *));
	init_random_matrix(A, A_rows, A_cols);
	init_random_matrix(B, B_rows, B_cols);

	result_normal = (double **) malloc(A_rows * sizeof(double *));
	result_openmp = (double **) malloc(A_rows * sizeof(double *));
	init_zero_matrix(result_normal, A_rows, B_cols);
	init_zero_matrix(result_openmp, A_rows, B_cols);

	// 测试十次普通运算时间取平均值
	double sum_time_normal = 0;
	for (int i = 0; i < 10; ++i)
	{
		matrix_mul_normal(A, A_rows, A_cols, B, B_cols,
		                  result_normal, &elapsed_time_normal);
		sum_time_normal += elapsed_time_normal;
		printf("Normal matrix multiplication %d times use : %.2f seconds to calculate\n", i + 1, elapsed_time_normal);

		init_zero_matrix(result_normal, A_rows, B_cols);
		elapsed_time_normal = 0;
	}
	printf("\n");

	// 测试十次openmp运算时间取平均值
	double sum_time_openmp = 0;
	for (int i = 0; i < 10; ++i)
	{
		matrix_mul_openmp(A, A_rows, A_cols, B, B_cols,
		                  result_openmp, &elapsed_time_openmp);
		sum_time_openmp += elapsed_time_openmp;
		printf("Normal matrix multiplication %d times use : %.2f seconds to calculate\n", i + 1, elapsed_time_openmp);

		init_zero_matrix(result_openmp, A_rows, B_cols);
		elapsed_time_openmp = 0;
	}
	printf("\n");

	printf("Normal matrix multiplication used average %.2f seconds to calculate\n", sum_time_normal / 10);
	printf("matrix multiplication using openmp used average  %.2f seconds to calculate\n", sum_time_openmp / 10);
	printf("openmp is %.2f times faster than normal\n\n", sum_time_normal / sum_time_openmp);

	if (compare_two_matrix(result_openmp, result_normal, result_rows, result_cols) == 0)
		printf("matrix multiplication using openmp's answer is wrong.\n");
	else
		printf("matrix multiplication using openmp's answer is right.\n");

//	print_matrix(A, A_rows, A_cols);
//	print_matrix(B, B_rows, B_cols);
//	print_matrix(result_normal, result_rows, result_cols);

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


// init_random_matrix用于生成随机矩阵
void init_random_matrix(double **matrix, int rows, int cols)
{
	double *aux_matrix = (double *) malloc(rows * cols * sizeof(double));

	for (int i = 0; i < rows; ++i)
	{
		matrix[i] = aux_matrix + i * cols;
//		matrix[i] = &(aux_matrix[i * cols]);
	}

	unsigned int seed = (unsigned int) (time(NULL) ^ getpid());
	srand(seed);

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			matrix[i][j] = (double) rand() / RAND_MAX;
		}
	}
}


// init_zero_matrix用于生成0矩阵
void init_zero_matrix(double **matrix, int rows, int cols)
{
	double *aux_matrix = (double *) malloc(rows * cols * sizeof(double));

	for (int i = 0; i < rows; ++i)
	{
		matrix[i] = aux_matrix + i * cols;
//		matrix[i] = &(aux_matrix[i * cols]);
	}

	unsigned int seed = (unsigned int) (time(NULL) ^ getpid());
	srand(seed);

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			matrix[i][j] = 0.0;
		}
	}
}


// print_matrix用于输出矩阵
void print_matrix(double **matrix, int rows, int cols)
{
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			printf("%.2f ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


// matrix_mul_normal用于计算普通矩阵乘法，并保存运算用时
void matrix_mul_normal(double **ma_a, int rows_a, int cols_a, double **ma_b, int cols_b, double **result, double *time)
{
	double start_time = omp_get_wtime();
	for (int i = 0; i < rows_a; ++i)
	{
		for (int j = 0; j < cols_b; ++j)
		{
			for (int k = 0; k < cols_a; k++)
			{
				result[i][j] += ma_a[i][k] * ma_b[k][j];
			}
		}
	}
	double end_time = omp_get_wtime();
	*time = end_time - start_time;
}

void matrix_mul_openmp(double **ma_a, int rows_a, int cols_a, double **ma_b, int cols_b, double **result, double *time)
{
	double start_time = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < rows_a; ++i)
	{
//		#pragma omp parallel for
//      慢了一点，不到四倍
		for (int j = 0; j < cols_b; ++j)
		{
//			#pragma omp parallel for
//          很慢，不到两倍
			for (int k = 0; k < cols_a; k++)
			{
				result[i][j] += ma_a[i][k] * ma_b[k][j];
			}
		}
	}
	double end_time = omp_get_wtime();
	*time = end_time - start_time;
}


// 用于比较两个矩阵是否相等
int compare_two_matrix(double **ma_a, double **ma_b, int rows, int cols)
{
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			if (ma_a[i][j] != ma_b[i][j])
			{
				return 0;
			}
		}
	}
	return 1;
}
