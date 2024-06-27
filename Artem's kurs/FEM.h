#pragma once
#include "stdafx.h"

using namespace std;

struct Vertex
{
	double x;
	double y;
};

struct Triangle
{
	int vert1, vert2, vert3;
	int region;
};

struct FirstBoundaryCondition
{
	int vert;
	int equationNum;
};

struct SecondBoundaryCondition
{
	int vert1, vert2;
	int equationNum;
};

struct ThirdBoundaryCondition
{
	int vert1, vert2;
	int betaEquationNum, UbetaEquationNum;
};

void GenerateGrid(int factor);
double AverageRate(vector<double> q0, vector<double> q1);

class FEM
{
public:
	vector<Vertex> vertices; // ���� �������
	vector<Triangle> tris; // �������� ��������
	vector<FirstBoundaryCondition> firstBoundary; // ������ ������� �������
	vector<SecondBoundaryCondition> secondBoundary; // ������ ������� �������
	vector<ThirdBoundaryCondition> thirdBoundary; // ������ ������� �������
	vector<double> timeStamps; // ���� �� �������
	vector<double> q, errVec; // ������ �������
	vector<double> q_1, q_2, q_3; // ���������� ������� q(j-1) and q(j-2)
	vector<vector<double>> qs;

	void Input(); // ���� ������ �� �����
	void Solve(); // ����� ������� ������� �������
	void PrintSolution();
	double CalculateError(int r_ind);
private:
	int regionsNum, globalN; // ���������� �������� � �����
	vector<int> ig, jg; // ���������� �������
	vector<double> gglA, gguA, diA, gglM, gguM, diM, gglG, gguG, diG, b, pureB;
	double G[3][3]{}; // ������ ������� G
	double M[3][3]{}; // ������ ������� M
	const double pureM[3][3] = { {2, 1, 1}, {1, 2, 1}, {1, 1, 2} }; // ������ ������� M ��� ����������� �� � �������� ��������� �� ������ ��������
	double localB[3]{}; // ��������� ������ b

	double Lamda(int vert, int region); // ���������� �������� ������
	double Sigma(int vert, int region); // ���������� �������� �����
	double Function(int vert, int region, int tInd); // ���������� �������� ������� f
	double Beta(int vert, int eqNum); // ���������� �������� ����
	double Ubeta(int vert, int eqNum); // ���������� �������� U ���� ��� 3 �������� �������
	double Theta(int vert, int eqNum); // ���������� �������� ���� ��� 2 �������� �������
	double Ug(int vert, int tInd); // ���������� �������� � ���� ��� 1 �������� �������
	double Ug(double x, double y, double t);
	double Uq(double x, double y, Triangle tri, vector<double> resQ);
	double GetAverageLamda(Triangle tri); // ���������� �������� ������ �� ��������
	double GetAverageSigma(Triangle tri); // ���������� �������� ����� �� ��������
	double DetD(Triangle tri); // ���������� ������������ D (��������� �������) ��������
	double Alpha(Triangle tri, int k, int i);  // ���������� �������� ����� ��� ���������� ������� G
	double EdgeLength(int vert1, int vert2); // ���������� ����� �����
	double DivGrad(int vert, int tInd, double h); // ���������� div(gradu) ��� �������������� ��������� f
	double CalcNorm(vector<double> resQ, int tInd); // ���������� ����� ����������� ������� � ������� ���������
	
	double Uq(double t, int r_ind, int tj);
	int IndexOfUnknown(Triangle tri, int i); // ��������� ����������� ������ ���� �� ���������� � ��������
	void FormLocalM(Triangle tri); // ������������ ������� G
	void FormLocalG(Triangle tri); // ������������ ������� M
	void FormGlobalMatrices(int tInd); // ������������ ���������� M � G
	void FormGlobalB(double deltaT0, double deltaT1, double deltaT2, double deltaT3, double deltaT4, double deltaT5);
	void FormPortrait(); // ������������ �������� ���������� �������
	void ResolveBoundaries(int tInd); // ���� ���� ������� �������
	void AddToGlobalA(int i, int j, double add); // ���������� �������� � ���������� ������� A
	void AddToGlobalG(int i, int j, double add); // ���������� �������� � ���������� ������� G
	void AddToGlobalM(int i, int j, double add); // ���������� �������� � ���������� ������� M
	void AllocateGlobalMatrices(); // ��������� ������ ��� ���������� �������
	void FormLocalB(Triangle tri, int tInd); // ������������ ���������� ������� b
	void TwoLayerScheme(int j); // ����������� ����� ��� ������� ������
	void ClearMatrices(); // ������� ���������� ������

};