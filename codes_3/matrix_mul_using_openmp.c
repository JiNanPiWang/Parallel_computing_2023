//
// Created by Administrator on 2023/7/14.
//
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

// 用于处理输入参数，并判断是否有效
int parseArgument(const char *arg, int *value);

int main(int argc, char *argv[])
{
	// 计算矩阵A×B，使用openmp和普通方法，并进行对比
	int A_rows, A_cols, B_cols;
	double *A_aux; //auxiliary 1D for 2D matrix a
	double *B_aux; //auxiliary 1D for 2D matrix b
	double *result_normal_aux; //auxiliary 1D for 2D matrix c
	double *result_openmp_aux; //auxiliary 1D for 2D matrix c
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
		printf("A_rows = %d, A_cols = %d, B_cols = %d\n\n", A_rows, A_cols, B_cols);
	}
	else
	{
		printf("Incorrect number of arguments\n");
		return 2;
		// 执行适当的错误处理措施
	}

	return 0;
}

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