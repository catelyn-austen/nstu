#include "stdafx.h"

const int n = 2500;
int main()
{
	int i = 0, j = 0, k = 0;
	double norm = 0;
	// создание матриц А и В
	double* A = NULL;
	double* B = NULL;
	double* C = NULL;
	if (A == NULL) A = new double[n * n];
	if (B == NULL) B = new double[n * n];
	if (C == NULL) C = new double[n * n];

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			A[i * n + j] = i <= j ? 1 : 0;
			B[i * n + j] = i >= j ? 1 : 0;
		}
	}

	// перемножение
	clock_t start = clock();
#pragma omp parallel for private(j,k) num_threads(8)
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			C[i * n + j] = 0;
			for (k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[k * n + j];
		}
	}
	clock_t end = clock();


	// норма матрицы С
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			norm = norm + pow(C[i * n + j], 2);
	norm = sqrt(norm);

	double s = (double)(end - start) / CLOCKS_PER_SEC;
	cout << s << endl;
	printf("norm = %f", norm);

	// удаление памяти
	if (A) { delete[] A; A = NULL; }
	if (B) { delete[] B; B = NULL; }
	if (C) { delete[] C; C = NULL; }
}