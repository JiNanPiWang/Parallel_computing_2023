//
// Created by Tianyi Zhang on 2023/7/23.
// an MPI program to compute matrix-vector multiplication y = Ax
//
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>


// parseArgument用于处理输入参数，并判断是否有效
int parseArgument(const char *arg, int *value);

// init_random_matrix用于生成随机矩阵
void init_random_matrix(double **matrix, int rows, int cols);

// init_zero_matrix用于生成0矩阵
void init_zero_matrix(double **matrix, int rows, int cols);

// 随机初始化数组
void init_random_vec(double *vec, int vec_len);

// 随机初始化数组
void init_zero_vec(double *vec, int vec_len);

// print_matrix用于输出矩阵
void print_matrix(double **matrix, int rows, int cols);

// 输出数组
void print_array(double *array, int array_len);

// matrix_mul_normal用于计算普通矩阵乘法，并保存运算用时
void matrix_mul_normal(double **A, int rows_A, int cols_A, double *x, int rows_x, double *result, double *time);

// matrix_mul_openmp用于计算使用了openmp的矩阵乘法，并保存运算用时
void matrix_mul_openmp(double **ma_a, int rows_a, int cols_a, double **ma_b, int cols_b, double **result, double *time);

// 用于比较两个矩阵是否相等
int compare_two_matrix(double **ma_a, double **ma_b, int rows, int cols);

int main(int argc, char *argv[])
{
	// 计算矩阵A×B，使用openmp和普通方法，并进行对比
	int A_rows, A_cols, x_len, result_len; // B_rows = A_cols，不用输入
	double **A; //the two-dimensional input matrix
	double *x;
	double *result_normal; //the resulting matrix (normal version)
	double *result_openmp; //the resulting matrix (openmp version)
	double elapsed_time_normal, elapsed_time_openmp;

	// 处理输入参数
	if (argc == 3)
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
		x_len = A_cols;
		result_len = A_cols;
		printf("A_rows = %d, A_cols = %d\n\n", A_rows, A_cols);
	}
	else
	{
		printf("Incorrect number of arguments, valid input: A_rows, A_cols\n");
		return 2;
	}

	// 初始化
	A = (double **) malloc(A_rows * sizeof(double *));
	x = (double *) malloc(x_len * sizeof(double *));
	init_random_matrix(A, A_rows, A_cols);
	init_random_vec(x, x_len);

	result_normal = (double *) malloc(result_len * sizeof(double *));
	result_openmp = (double *) malloc(result_len * sizeof(double *));
	init_zero_vec(result_normal, result_len);
	init_zero_vec(result_openmp, result_len);

	matrix_mul_normal(A, A_rows, A_cols, x, x_len, result_normal, &elapsed_time_normal);

	print_matrix(A, A_rows, A_cols);
	print_array(x, x_len);
	print_array(result_normal, result_len);

//	int numprocs, myid;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


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


// 随机初始化数组
void init_random_vec(double *vec, int vec_len)
{
	unsigned int seed = (unsigned int) (time(NULL) ^ getpid());
	srand(seed);

	for (int i = 0; i < vec_len; ++i)
	{
		vec[i] = (double) rand() / RAND_MAX;
	}
}


// 随机初始化数组
void init_zero_vec(double *vec, int vec_len)
{
	memset(vec, 0, vec_len);
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


// 输出数组
void print_array(double *array, int array_len)
{
	for (int i = 0; i < array_len; ++i)
	{
		printf("%.2f ", array[i]);
	}
	printf("\n\n");
}


// matrix_mul_normal用于计算普通矩阵乘法，并保存运算用时
void matrix_mul_normal(double **A, int rows_A, int cols_A, double *x, int rows_x, double *result, double *time)
{
	double start_time = omp_get_wtime();
	for (int i = 0; i < rows_A; ++i)
	{
		for (int k = 0; k < cols_A; k++)
		{
			result[i] += A[i][k] * x[k];
		}
	}
	double end_time = omp_get_wtime();
	*time = end_time - start_time;
}


// matrix_mul_openmp用于计算使用了openmp的矩阵乘法，并保存运算用时
void matrix_mul_openmp(double **ma_a, int rows_a, int cols_a, double **ma_b, int cols_b, double **result, double *time)
{
	double start_time = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < rows_a; ++i)
	{
//		#pragma omp parallel for
		for (int j = 0; j < cols_b; ++j)
		{
//			#pragma omp parallel for
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
			if (fabs(ma_a[i][j] - ma_b[i][j]) > 1e-6)
			{
				return 0;
			}
		}
	}
	return 1;
}