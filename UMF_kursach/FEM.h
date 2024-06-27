#pragma once
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "vector"

using namespace std;

// ��������� � ������
struct Vertex
{
	double x;
	double y;
};
struct Triangle
{
	int v1, v2, v3;
	int area;
};
struct cond_1
{
	// ������ � �������
	int v;
	int eq_n;
};
struct cond_2
{
	int v1, v2;
	int eq_n;
};
struct cond_3
{
	int v1, v2;
	int beta_eq_n, ubeta_eq_n;
};

class FEM
{
public:

	/* -------- ������� � ��������� -------- */
	vector<Vertex> points; // �����
	vector<Triangle> finits; // �������� ��������
	vector<cond_1> border_cond1; // ������ � ����������� � 1 �������
	vector<cond_2> border_cond2; // � 2 �������
	vector<cond_3> border_cond3; // � 3 �������
	vector<double> time_points; // ���� �� �������
	vector<double> q; // ������� ������
	vector<double> errVec; // ������ ������
	vector<double> q_1, q_2; // q(j-1), q(j-2)
	vector<double> gglA, gguA, diA, gglM, gguM, diM, gglG, gguG, diG, b, F; // ��� ����
	vector<int> ig, jg; // ���������� �������
	vector<vector<double>> qs; // ������ ��������
	int areas_count, points_count; // ����� �������� � �����
	
	/* -------- ������� ��������������� -------- */
	double G[3][3]{}; // ��������� G
	double M[3][3]{}; // ��������� M
	const double coef_M[3][3] = { {2, 1, 1}, {1, 2, 1}, {1, 1, 2} }; // ������� ����-��� ��� M
	double local_F[3]{}; // ��������� ������ b

	/* -------- ���� ������ �� ������ -------- */
	void Input();

	/* -------- ������� -------- */
	void MainPass(); // �������� �������, ���������� � main
	int global_index(Triangle triangle, int i); // ���������� ����� ���� �� ����������
	void M_matrix(Triangle triangle); // ��������� �
	void G_matrix(Triangle triangle); // ��������� G
	void F_local(Triangle triangle, int ind_t); // ��������� F
	void Global_matrix(int ind_t); // ���������� M � G
	void Global_F(double T, double T0, double T1); // ���������� F
	void Portrait(); // �������
	void cond_123(int ind_t); // �������
	void add_to_A(int i, int j, double to_add); // ���������� � ���������� �
	void add_to_G(int i, int j, double to_add); // ���������� � ���������� G
	void add_to_M(int i, int j, double to_add); // ���������� � ���������� �
	void allocate_memory(); // ��������� ������ (��� ���������� ������)
	void clearing(); // �������� ���������� ������

	/* -------- ������� � ����-�� -------- */
	double coef_lambda(int p, int area); // ������
	double coef_sigma(int p, int area); // �����
	double f(int p, int area, int ind_t); // �������� � �����
	double coef_beta(int p, int number); // ����
	double coef_ubeta(int p, int number); // u����
	double coef_tetta(int p, int number); // �����
	double Ug(int p, int ind_t); // Ug
	double Ug(double x, double y, double t);
	double Uq(double x, double y, Triangle triangle, vector<double> res); // ������� ��� ����� ����������� �� �����������
	double aver_coef_sigma(Triangle triangle); // ����������� �����
	double Det(Triangle triangle); // ������������ ������� D
	double coef_alpha(Triangle triangle, int k, int i);  // ����-�� ������� D
	double Hm(int p1, int p2); // ����� ����� hm
	double Int_norm_error(vector<double> res, int ind_t); // ����� ����������� ������� �� �����������

	/* -------- ����������� ������� �� ������� -------- */
	double Int_norm_error_time(int p); // ����� ����������� ������� �� �������
	double Q_resh(double t, int p, int ind_t); // ������� �� ind_t ��������� �������
};