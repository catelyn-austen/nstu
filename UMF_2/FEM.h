#pragma once
#include <stdio.h>
#include <vector>
#include <queue>
#include <math.h>
#include <iostream>
#define REALOUT "%.15lf\t"

using namespace std;

/*------ ��������� ��� ������ � �������� ��������� ------*/

struct BoundCond {
	int priority;
	int type;
	char orientation;
	bool operator < (BoundCond const& b2) const
	{
		return priority > b2.priority;
	}
};

/*------ ��������� ��� ������ � ��� ------*/

struct FEM
{
	unsigned int total_n; // ���������� ����� �� �����������
	unsigned int total_n_t; // ���������� ����� �� �������
	unsigned int time_layers; // ��������� ����
	unsigned int test; // ����� �����

	double eps = 1e-13;

	// ������� ��� ������ � ��� LU
	vector<double> x_grid;
	vector<double> y;
	vector<double> t;
	vector<double> p_i_prev;
	vector<double> p_i;
	vector<double> b;
	vector<double> al;
	vector<double> au;
	vector<double> di;
	vector<double> alLU;
	vector<double> auLU;
	vector<double> diLU;
	vector<double> tmp;
	vector<int> ia;
	vector<int> iter;
	vector<vector<double>> q_layers;
	priority_queue<BoundCond> boundcond;

	/*------ ������� � ��������� ------*/

	double Uf(int number, double x, double t);
	double f(int number, double x, double t);
	double lambda(int number, double x);
	double gamma(int number, double x);

	/*------ ������ ������ � ���������� �������� ------*/

	void Input(FILE* Params, FILE* Points, FILE* Bounds, FILE* Time);
	void Allocate_Memory(int n);

	/*------ ������� ������ ------*/

	void Solve();
	void Main_Part(int curTime);
	void Condition_account(int curTime);
	void Clearing();

	/*------ LU-���������� ------*/

	void LU_decompose();
	void Copy_vector(vector<double>& v1, vector<double>& v2);
	double Vector_norm(vector<double> v);
	double Rel_discrepancy();
	void Vector_sub(vector<double> first, vector<double> second, vector<double>& result);
	void Matrix_mult_vector(vector<double> v1, vector<double>& res);
	void CalculateX();
	void CalculateY();

	/*------ ����� �� �2 ------*/

	double Lebeg2(int timelayer);
	double Gauss3(int timelayer);
	double Lebeg_function(double x, int timelayer);
	double Find_inner_points(double point, double t_layer);
	
	/*----- ����� ������� ------*/

	void Results();
};