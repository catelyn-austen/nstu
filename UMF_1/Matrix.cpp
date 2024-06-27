#include "matrix.h"

void SLAE::IterativeMethod(int NumOfW)
{
	double* curX, * prevX;
	VecotorCopy(x0, local_x0);
	if (method == 0) //Якоби
	{
		prevX = local_x0;
		curX = x;
	}
	else { //Зейдель
		curX = local_x0;
		prevX = local_x0;
	}
	double normB = VectorNorm(b);
	double RelativeDiscrepancy = CalculateRelativeDiscrepancy(x0);
	int curIteration = 0;
	double DiscrepancyF_Ax = 0;
	for (; curIteration < MaxNumOfIterations and RelativeDiscrepancy > AccuracyOfTheSolution; curIteration++)
	{
		DiscrepancyF_Ax = 0;
		for (int i = 0; i < n; i++)
		{
			int indX = 0;
			double sum = 0;
			for (int j = 0; j < 2; j++)
			{
				indX = i + TableOfShifts[j];
				if (indX >= 0)
				{
					sum += prevX[indX] * matrix[i][j];
				}
			}
			sum += prevX[i] * matrix[i][2];
			for (int j = 3; j < 5; j++)
			{
				indX = i + TableOfShifts[j];
				if (indX < n)
				{
					sum += prevX[indX] * matrix[i][j];
				}
			}
			curX[i] = prevX[i] + (b[i] - sum) / matrix[i][2];
			DiscrepancyF_Ax += (b[i] - sum) * (b[i] - sum);
		}
		std::swap(curX, prevX);
		RelativeDiscrepancy = sqrt(DiscrepancyF_Ax) / normB;
		//printf_s("%.15lf\n", RelativeDiscrepancy);
		if (isinf(RelativeDiscrepancy) or isnan(RelativeDiscrepancy))
			break;
	}

	VecotorCopy(prevX, x);

	TableOfNumOfConditionality[NumOfW - 1] = CalculateNumOfConditionality(RelativeDiscrepancy);
	if (curIteration >= MaxNumOfIterations)
	{
		//printf("Выход по причиние выхода за макс. число итераций\n");
		NumOfIterationsDependingOnW[NumOfW - 1] = -1;
	}
	else if (RelativeDiscrepancy < AccuracyOfTheSolution) {
		//printf("Выход по причиние достигнутой желаемой погрешности\n");
		printf("%.2lf\n", w);
		printf("%E\n", RelativeDiscrepancy);
		NumOfIterationsDependingOnW[NumOfW - 1] = curIteration;
		printf("%d\n", curIteration);
		VectorOutput(x);
	}
	else if (isnan(RelativeDiscrepancy)) {
		//printf("Невязка nan\n");
		NumOfIterationsDependingOnW[NumOfW - 1] = -2;
	}
	else if (isinf(RelativeDiscrepancy)) {
		//printf("Невязка inf\n");
		NumOfIterationsDependingOnW[NumOfW - 1] = -2;
	}
}

void SLAE::Input(FILE* matrixFile, FILE* vectorFile, FILE* paramFile)
{
	fscanf_s(matrixFile, "%d", &n);
	fscanf_s(matrixFile, "%d", &m);

	fscanf_s(paramFile, "%d", &method);
	fscanf_s(paramFile, "%lf", &AccuracyOfTheSolution);
	fscanf_s(paramFile, "%d", &MaxNumOfIterations);
	AllocateMemory();
	InitializeShiftsTable();
	for (int Icount = 0; Icount < 9; Icount++)
	{
		int curDiag = TableOfShifts[Icount];
		if (curDiag <= 0)
		{
			for (int i = abs(TableOfShifts[Icount]); i < n; i++)
			{
				fscanf_s(matrixFile, "%lf", &matrix[i][Icount]);
			}
		}
		else
		{
			for (int i = 0; i < (n - abs(curDiag)); i++)
			{
				fscanf_s(matrixFile, "%lf", &matrix[i][Icount]);
			}
		}

	}

	for (int i = 0; i < n; i++)
		fscanf_s(vectorFile, "%lf", &b[i]);
	for (int i = 0; i < n; i++)
		fscanf_s(vectorFile, "%lf", &x0[i]);

	for (int i = 0; i < n; i++)
	{
		xtrue[i] = (double)(i + 1);
	}
}

void SLAE::Input(double** inMatrix, double* inX0, double* inB, int inN, int inM, double eps, int maxIter)
{
	n = inN;
	m = inM;

	method = 1;
	AccuracyOfTheSolution = eps;
	MaxNumOfIterations = maxIter;
	AllocateMemory();
	InitializeShiftsTable();
	
	matrix = inMatrix;

	b = inB;
	x0 = inX0;

	for (int i = 0; i < n; i++)
	{
		xtrue[i] = (double)(i + 1);
	}
}

void SLAE::AllocateMemory()
{
	matrix = new double* [n];
	for (int i = 0; i < n; ++i) {
		matrix[i] = new double[n];
	}
	b = new double[n];
	x = new double[n];
	x0 = new double[n];
	xtrue = new double[n];
	local_x0 = new double[n];
	vectorForDiscrepancy = new double[n];
	NumOfIterationsDependingOnW = new int[200];
	TableOfNumOfConditionality = new double[200];
}



void SLAE::MatrixVectorMultiplicationForDiscrepancy(double* first)
{
	for (int i = 0; i < n; i++)
	{
		int indX = 0;
		double sum = 0;
		for (int j = 2; j < 5; j++)
		{
			indX = i + TableOfShifts[j];
			if (indX < n)
			{
				sum += first[indX] * matrix[i][j];
			}
		}
		sum += first[i] * matrix[i][4];
		for (int j = 0; j < 2; j++)
		{
			indX = i + TableOfShifts[j];
			if (indX >= 0)
			{
				sum += first[indX] * matrix[i][j];
			}
		}
		vectorForDiscrepancy[i] = b[i] - sum;
		//cout << vectorOut[i] << endl;
	}

}

double SLAE::CalculateRelativeDiscrepancy(double* first)
{
	MatrixVectorMultiplicationForDiscrepancy(first);
	return VectorNorm(vectorForDiscrepancy) / VectorNorm(b);
}

double SLAE::CalculateNumOfConditionality(double RelativeDiscrepancy)
{
	VecotorSubtract(x, xtrue);
	double VectorXRelDiscrepancy = VectorNorm(vectorForDiscrepancy) / VectorNorm(xtrue);
	return VectorXRelDiscrepancy / RelativeDiscrepancy;
}

void SLAE::OutputDense()
{
	double** matrixDense;
	matrixDense = new double* [n];
	for (int i = 0; i < n; ++i) {
		matrixDense[i] = new double[n];
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			matrixDense[i][j] = 0.0;
		}
	}

	for (int Icount = 0; Icount < 5; Icount++)
	{
		int curDiagonal = TableOfShifts[Icount];
		if (curDiagonal <= 0)
		{
			for (int i = abs(curDiagonal), j = 0; i < n and j < n; i++, j++)
			{
				matrixDense[i][j] = matrix[i][Icount];
			}
		}
		else {
			for (int i = 0, j = abs(curDiagonal); i < n and j < n; i++, j++)
			{
				matrixDense[i][j] = matrix[i][Icount];
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf(REALOUTD, matrixDense[i][j]);
		}
		printf("\n");
	}
}

double SLAE::VectorNorm(double* first)
{
	double norm = 0;
	for (int i = 0; i < n; i++)
	{
		norm += first[i] * first[i];
	}
	return sqrt(norm);
}

void SLAE::VecotorSubtract(double* first, double* second)
{
	for (int i = 0; i < n; i++)
	{
		vectorForDiscrepancy[i] = first[i] - second[i];
	}
}

void SLAE::VecotorCopy(double* first, double* second)
{
	for (int i = 0; i < n; i++)
	{
		second[i] = first[i];
	}
}

void SLAE::InitializeShiftsTable()
{
	TableOfShifts[0] -= m;
	TableOfShifts[4] += m;
}

void SLAE::OutputSolutionVector(FILE* out) {
	for (int i = 0; i < n; i++)
		fprintf_s(out, REALOUT, x[i]);
	fprintf_s(out, "\n");
}

void SLAE::VectorOutput(double* curX)
{
	for (int i = 0; i < n; i++)
		printf_s(REALOUT, curX[i]);
	printf_s("\n");
}

void SLAE::OutputResultParametrs()
{
	for (int i = 1; i <= (200); i++)
	{
		printf("%.2lf ", 0.01 * i);
		printf("%d ", NumOfIterationsDependingOnW[i - 1]);
		printf("%lf\n", TableOfNumOfConditionality[i - 1]);
	}
}