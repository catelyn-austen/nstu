#pragma once
#include <vector>;

using namespace std;

class SLAE {
public:
	int n, aulS;
	int max_iter;
	double eps;
	vector<double> di;
	vector<double> au;
	vector<double> al;
	vector<double> b;
	vector<double> d1;
	vector<int> ja;
	vector<int> ia;


	vector<double> r_prev;
	vector<double> r_new;
	vector<double> p_prev;
	vector<double> p_new;
	vector<double> z_prev;
	vector<double> z_new;
	vector<double> s_prev;
	vector<double> s_new;
	vector<double> x_prev;
	vector<double> x_new;

	void MatrixMultVector(vector<double>& v, vector<double>& res);
	double VectorMultVector(vector<double>& v1, vector<double>& v2);
	double VectorNorm(vector<double>& v);
	//void CopyVector(vector<double>& v1, vector<double>& v2, int size);
	void VectSub(vector<double>& v1, vector<double>& v2, vector<double>& res);
	void TransponMultVector(vector<double>& v, vector<double>& res);

	int Init();
	//void Allocate_memory(int size);
	void BCG();

	void PrintX();
};