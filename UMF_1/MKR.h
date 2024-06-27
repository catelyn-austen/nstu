#pragma once
#include "stdafx.h"

using namespace std;

struct BoundaryCondition {
	int priority;
	int type;
	char orientation;

	bool operator<(BoundaryCondition const& b2) const
	{
		return priority > b2.priority; //Поменял знак с < на >
	}
};

class MKR
{
public:
	int variant;

	void Input(FILE* inParams, FILE* inPoints, FILE* inBoundary);
	void Solve();
	void OutputSolutionVector();
	void PrintTrue();
	void PrintDes();
private:
	int k1, k2, total_n, total_m;
	double lambda, gamma, hx, hy;
	double** A;
	double *q, *f;
	vector<double> x;
	vector<double> y;
	priority_queue<BoundaryCondition> boundaryConditions;

	double Func(int i, int j);
	double U(int i, int j);
	double Theta(int i, int j, int nx, int ny);
	void AllocateMemory(int totalPoints);
	void MainPass(); // Главный проход по матрице, т.е без краевых
	void BoundaryPass(); // Второй проход по матрице для краевых
	void AddDerivativeX(double** &A, int i, int j, int max);
	void AddDerivativeY(double** &A, int i, int j, int max);
	int GetPointNum(int i, int j) const;
};	