//
// Created by Administrator on 2023/7/14.
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

int matrix_mul_normal(double **A, double **B);

int main(int argc, char *argv[])
{
	// 计算矩阵A×B，使用openmp和普通方法，并进行对比
	int A_rows, A_cols, B_rows, B_cols; // B_rows = A_cols，不用输入
	double **A; //the two-dimensional input matrix
	double **B; //the two-dimensional input matrix
	double **result_normal; //the resulting matrix (normal version)
	double **result_openmp; //the resulting matrix (openmp version)

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

//	print_matrix(result_normal, A_rows, B_cols);

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