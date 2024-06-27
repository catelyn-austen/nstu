#pragma once
#include "stdafx.h"

#define REALOUT "%.15lf\n"
#define REALOUTD "%.3lf\t"
#define EPS 1e-13
using namespace std;

class SLAE {
public:
	int n, m, MaxNumOfIterations, method = 1;
	int* NumOfIterationsDependingOnW;
	double AccuracyOfTheSolution, w = 0.50;
	double* x, * x0, * b, * vectorForDiscrepancy, * local_x0, * TableOfNumOfConditionality, *xtrue;
	double** matrix;
	int TableOfShifts[5] = {-2, -1, 0, 1, 2};

	void Input(FILE* matrixFile, FILE* vectorFile, FILE* paramFile);
	void Input(double** inMatrix, double* inX0, double* inB, int inN, int inM, double eps, int maxIter);
	void OutputDense();
	void IterativeMethod(int NumOfW);
	void VectorOutput(double* curX);
	void OutputSolutionVector(FILE* out);
	void InitializeShiftsTable();
	double VectorNorm(double* first);
	double CalculateRelativeDiscrepancy(double* first);
	double CalculateNumOfConditionality(double RelativeDiscrepancy);
	void MatrixVectorMultiplicationForDiscrepancy(double* vectorMult);
	void VecotorCopy(double* first, double* second);
	void VecotorSubtract(double* first, double* second);
	void OutputResultParametrs();
protected:
	void AllocateMemory();
};