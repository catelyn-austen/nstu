#pragma once
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "vector"

using namespace std;

// структуры в помощь
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
	// только в вершине
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

	/* -------- Вектора и константы -------- */
	vector<Vertex> points; // точки
	vector<Triangle> finits; // конечные элементы
	vector<cond_1> border_cond1; // вектор с информацией о 1 краевых
	vector<cond_2> border_cond2; // о 2 краевых
	vector<cond_3> border_cond3; // о 3 краевых
	vector<double> time_points; // узлы по времени
	vector<double> q; // искомый вектор
	vector<double> errVec; // вектор ошибок
	vector<double> q_1, q_2; // q(j-1), q(j-2)
	vector<double> gglA, gguA, diA, gglM, gguM, diM, gglG, gguG, diG, b, F; // для СЛАУ
	vector<int> ig, jg; // глобальная матрица
	vector<vector<double>> qs; // вектор векторов
	int areas_count, points_count; // число областей и точек
	
	/* -------- Матрицы вспомогательные -------- */
	double G[3][3]{}; // локальная G
	double M[3][3]{}; // локальная M
	const double coef_M[3][3] = { {2, 1, 1}, {1, 2, 1}, {1, 1, 2} }; // матрица коэф-тов при M
	double local_F[3]{}; // Локальный вектор b

	/* -------- Ввод данных из файлов -------- */
	void Input();

	/* -------- Решение -------- */
	void MainPass(); // основная функция, вызываемая в main
	int global_index(Triangle triangle, int i); // глобальный номер узла из локального
	void M_matrix(Triangle triangle); // локальная М
	void G_matrix(Triangle triangle); // локальная G
	void F_local(Triangle triangle, int ind_t); // локальный F
	void Global_matrix(int ind_t); // глобальные M и G
	void Global_F(double T, double T0, double T1); // глобальный F
	void Portrait(); // Портрет
	void cond_123(int ind_t); // краевые
	void add_to_A(int i, int j, double to_add); // добавление в глобальную А
	void add_to_G(int i, int j, double to_add); // добавление в глобальную G
	void add_to_M(int i, int j, double to_add); // добавление в глобальную М
	void allocate_memory(); // выделение памяти (для глобальных матриц)
	void clearing(); // очищение глобальных матриц

	/* -------- Функции и коэф-ты -------- */
	double coef_lambda(int p, int area); // лямбда
	double coef_sigma(int p, int area); // сигма
	double f(int p, int area, int ind_t); // значения в точке
	double coef_beta(int p, int number); // бета
	double coef_ubeta(int p, int number); // uбета
	double coef_tetta(int p, int number); // тетта
	double Ug(int p, int ind_t); // Ug
	double Ug(double x, double y, double t);
	double Uq(double x, double y, Triangle triangle, vector<double> res); // функция для нормы погрешности по координатам
	double aver_coef_sigma(Triangle triangle); // усредненная сигма
	double Det(Triangle triangle); // определитель матрицы D
	double coef_alpha(Triangle triangle, int k, int i);  // коэф-ты матрицы D
	double Hm(int p1, int p2); // длина ребра hm
	double Int_norm_error(vector<double> res, int ind_t); // норма погрешности решения по координатам

	/* -------- Погрешность решения по времени -------- */
	double Int_norm_error_time(int p); // норма погрешности решения по времени
	double Q_resh(double t, int p, int ind_t); // решение на ind_t временном отрезке
};