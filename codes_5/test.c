#include <stdio.h>
#include <omp.h>
#include <math.h>

#define N 4

void printMatrix(double A[N][N + 1])
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			printf("%.2f\t", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void gaussianElimination(double A[N][N + 1])
{
	// 主元的遍历
	for (int k = 0; k < N; k++)
	{
		int max_row = k;
		double max_val = A[k][k];

		// 查找主元列中绝对值最大的行
		#pragma omp parallel for
		for (int i = k + 1; i < N; i++)
		{
			if (fabs(A[i][k]) > max_val)
			{
				// 通过临界区，确保对max_val和max_row的更新是线程安全的，避免竞态条件。
				#pragma omp critical
				{
					if (fabs(A[i][k]) > max_val)
					{
						max_val = fabs(A[i][k]);
						max_row = i;
					}
				}
			}
		}

		// 交换行，将绝对值最大的元素放到主元位置
		if (max_row != k)
		{
			#pragma omp parallel for
			for (int j = k; j <= N; j++)
			{
				double temp = A[k][j];
				A[k][j] = A[max_row][j];
				A[max_row][j] = temp;
			}
		}

		// 消元操作
		#pragma omp parallel for
		for (int i = k + 1; i < N; i++)
		{
			double factor = A[i][k] / A[k][k];
			for (int j = k; j <= N; j++)
			{
				A[i][j] -= factor * A[k][j];
			}
		}
	}
}

void backSubstitution(double A[N][N + 1], double x[N])
{
	for (int i = N - 1; i >= 0; i--)
	{
		x[i] = A[i][N];
		for (int j = i + 1; j < N; j++)
		{
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
}

int main()
{
	double A[N][N + 1] = {
			{2,  -1, 1, 3,  10},
			{1,  3,  2, 8,  13},
			{-2, -2, 2, -3, -2},
			{3,  1,  3, -2, 15}
	};

	double x[N];

	printf("Original Augmented Matrix:\n");
	printMatrix(A);

	gaussianElimination(A);

	printf("Matrix after Gaussian Elimination:\n");
	printMatrix(A);

	backSubstitution(A, x);

	printf("Solution:\n");
	for (int i = 0; i < N; i++)
	{
		printf("x[%d] = %.2f\n", i, x[i]);
	}

	return 0;
}
