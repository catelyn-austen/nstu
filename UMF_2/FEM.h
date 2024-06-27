#pragma once
#include <stdio.h>
#include <vector>
#include <queue>
#include <math.h>
#include <iostream>
#define REALOUT "%.15lf\t"

using namespace std;

/*------ Структура для работы с краевыми условиями ------*/

struct BoundCond {
	int priority;
	int type;
	char orientation;
	bool operator < (BoundCond const& b2) const
	{
		return priority > b2.priority;
	}
};

/*------ Структура для работы с МКЭ ------*/

struct FEM
{
	unsigned int total_n; // количество узлов по координатам
	unsigned int total_n_t; // количество узлов по времени
	unsigned int time_layers; // временные слои
	unsigned int test; // номер теста

	double eps = 1e-13;

	// вектора для метода и для LU
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

	/*------ Функции и константы ------*/

	double Uf(int number, double x, double t);
	double f(int number, double x, double t);
	double lambda(int number, double x);
	double gamma(int number, double x);

	/*------ Чтение файлов и заполнение массивов ------*/

	void Input(FILE* Params, FILE* Points, FILE* Bounds, FILE* Time);
	void Allocate_Memory(int n);

	/*------ Решение задачи ------*/

	void Solve();
	void Main_Part(int curTime);
	void Condition_account(int curTime);
	void Clearing();

	/*------ LU-разложение ------*/

	void LU_decompose();
	void Copy_vector(vector<double>& v1, vector<double>& v2);
	double Vector_norm(vector<double> v);
	double Rel_discrepancy();
	void Vector_sub(vector<double> first, vector<double> second, vector<double>& result);
	void Matrix_mult_vector(vector<double> v1, vector<double>& res);
	void CalculateX();
	void CalculateY();

	/*------ Норма по Л2 ------*/

	double Lebeg2(int timelayer);
	double Gauss3(int timelayer);
	double Lebeg_function(double x, int timelayer);
	double Find_inner_points(double point, double t_layer);
	
	/*----- Вывод решения ------*/

	void Results();
};