#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>

using namespace std;

#pragma warning(disable:4996)

struct rectangle
{
    //Номер элемента
    int number;
    //Координаты для х
    int x_left, x_right;
    //Координаты для y
    int y_down, y_up;
};

class MKR
{
public:
    /*-------- Массивы, вектора переменные --------*/
    const int max_iter = 100000;
    const double eps = 1e-13;
    const double w = 1;
    int x_points_n, y_points_n, areas_n; //Длины массивов Xw и Yw, кол-во подобластей
    int x_part_n, y_part_n; //Размеры массивов X и Y
    int bound_n; //Количество краевых условий
    int N; //Размерность глобальной матрицы
    vector<double> x_points; //Массив X без разбиений (просто исходные точки по Х)
    vector<double> y_points; //Массив Y без разбиений (просто исзодные точки по У)
    vector<double> x_part; //Массив X c разбиениями (точки после применения q)
    vector<double> y_part; //Массив Y с разбиениями (точки после применения q)
    vector<double> di1;
    vector<double> di2;
    vector<double> di3;
    vector<double> di4;
    vector<double> di5;
    vector<double> b;
    vector<double> u;
    vector<double> uk;
    vector<int> x_i; //Массив индекса элементов в Х, которые есть в X_w (буквально 0, 1, 2...)
    vector<int> y_i; //Массив индекса элементов в Y, которые есть в Y_w
    vector<int> ig;
    vector<vector<int>> bounds; //Матрица краевых условий
    vector<rectangle> sub_areas; //Массив подобластей
    
    /*-------- Функции --------*/
    double Uf(double x, double y, int i);
    double lambda(int i);
    double gamma(int i);
    double f(double x, double y, int i);
    double Ug(double x, double y, int i);
    double Thetta(double x, double y, int i);

    /*-------- Считывание данных из файлов --------*/
    void input_grid();
    void input_part();
    void input_bounds();

    /*-------- Функции для составления СЛАУ --------*/
    void AllocateMemory();
    void MainPass();
    void BoundaruPass();
    void FictionNodes();
    bool IsFict(int s, int p);
    int GetAreaNum(int s, int p);
    //int GetPointNum(double x, double y); 
    void AddDerivative(int i, int j, double elem);
    void ClearString(int i, int diValue);
    //bool compareRows(const vector<int>& row1, const vector<int>& row2);
    void kraevoe_1(vector<int>& kr_1);
    void kraevoe_2(vector<int>& kr_2);

    /*-------- Функции для решения СЛАУ Гауссом-Зейделем --------*/
    double NormVector(vector <double>& a);
    double RelDiscrepancy(vector <double>& u);
    void Solve();

    /*-------- Вывод --------*/
    void Output();
};
