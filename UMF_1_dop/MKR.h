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
    //����� ��������
    int number;
    //���������� ��� �
    int x_left, x_right;
    //���������� ��� y
    int y_down, y_up;
};

class MKR
{
public:
    /*-------- �������, ������� ���������� --------*/
    const int max_iter = 100000;
    const double eps = 1e-13;
    const double w = 1;
    int x_points_n, y_points_n, areas_n; //����� �������� Xw � Yw, ���-�� �����������
    int x_part_n, y_part_n; //������� �������� X � Y
    int bound_n; //���������� ������� �������
    int N; //����������� ���������� �������
    vector<double> x_points; //������ X ��� ��������� (������ �������� ����� �� �)
    vector<double> y_points; //������ Y ��� ��������� (������ �������� ����� �� �)
    vector<double> x_part; //������ X c ����������� (����� ����� ���������� q)
    vector<double> y_part; //������ Y � ����������� (����� ����� ���������� q)
    vector<double> di1;
    vector<double> di2;
    vector<double> di3;
    vector<double> di4;
    vector<double> di5;
    vector<double> b;
    vector<double> u;
    vector<double> uk;
    vector<int> x_i; //������ ������� ��������� � �, ������� ���� � X_w (��������� 0, 1, 2...)
    vector<int> y_i; //������ ������� ��������� � Y, ������� ���� � Y_w
    vector<int> ig;
    vector<vector<int>> bounds; //������� ������� �������
    vector<rectangle> sub_areas; //������ �����������
    
    /*-------- ������� --------*/
    double Uf(double x, double y, int i);
    double lambda(int i);
    double gamma(int i);
    double f(double x, double y, int i);
    double Ug(double x, double y, int i);
    double Thetta(double x, double y, int i);

    /*-------- ���������� ������ �� ������ --------*/
    void input_grid();
    void input_part();
    void input_bounds();

    /*-------- ������� ��� ����������� ���� --------*/
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

    /*-------- ������� ��� ������� ���� �������-�������� --------*/
    double NormVector(vector <double>& a);
    double RelDiscrepancy(vector <double>& u);
    void Solve();

    /*-------- ����� --------*/
    void Output();
};
