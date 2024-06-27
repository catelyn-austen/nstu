#include <iostream>
#include <fstream>
#include <math.h>
#include "SLAE.h"

using namespace std;

int SLAE::Init()
{
	int err = 0;

	ifstream in;
	in.open("test1/kuslau.txt");
	if (!in.is_open())
		err++;
	in >> n >> max_iter >> eps;
	in.close();

	di.resize(n);
	ia.resize(n + 1);
	d1.resize(n);
	b.resize(n);
	x_prev.resize(n);
	x_new.resize(n);
	r_prev.resize(n);
	r_new.resize(n);
	z_prev.resize(n);
	z_new.resize(n);
	p_prev.resize(n);
	p_new.resize(n);
	s_prev.resize(n);
	s_new.resize(n);

	for (int i = 0; i < n; i++)
	{
		x_prev[i] = 0;
		x_new[i] = 0;
		r_prev[i] = 0;
		r_new[i] = 0;
		z_prev[i] = 0;
		z_new[i] = 0;
		p_prev[i] = 0;
		p_new[i] = 0;
		s_prev[i] = 0;
		s_new[i] = 0;
	}

	in.open("test1/ia.txt");
	if (!in.is_open())
		err++;
	for (int i = 0; i < n + 1; i++)
		in >> ia[i];
	aulS = ia[n] - ia[0];
	in.close();

	au.resize(aulS);
	al.resize(aulS);
	ja.resize(aulS);

	in.open("test1/di.txt");
	if (!in.is_open())
		err++;
	for (int i = 0; i < n; i++)
		in >> di[i];
	in.close();

	in.open("test1/au.txt");
	if (!in.is_open())
		err++;
	for (int i = 0; i < aulS; i++)
		in >> au[i];
	in.close();

	in.open("test1/al.txt");
	if (!in.is_open())
		err++;
	for (int i = 0; i < aulS; i++)
		in >> al[i];
	in.close();

	in.open("test1/ja.txt");
	if (!in.is_open())
		err++;
	for (int i = 0; i < aulS; i++)
		in >> ja[i];
	in.close();

	in.open("test1/f.txt");
	if (!in.is_open())
		err++;
	for (int i = 0; i < n; i++)
		in >> b[i];
	in.close();

	/*for (int i = 0; i < aulS; i++)
		ja[i]--;
	for (int i = 0; i < n + 1; i++)
		ia[i]--;*/

	for (int i = 0; i < n; i++)
		d1[i] = 1. / di[i];

	if (err)
	{
		cout << "Error while reading files" << endl;
		return 1;
	}


	return 0;
}

void SLAE::BCG()
{
	vector<double> Ax0;
	Ax0.resize(n);
	vector<double> res;
	res.resize(n);
	vector<double> ATs;
	ATs.resize(n);
	vector<double> Az;
	Az.resize(n);
	vector<double> Ax_b;
	Ax_b.resize(n);
	vector<double> forNorm;
	forNorm.resize(n);

	//x_prev = { 1,2,3,4,5 };
	//TransponMultVector(x_prev, Ax0);
	MatrixMultVector(x_prev, Ax0);
	double ak = 0, bk = 0;
	for (int i = 0; i < r_prev.size(); i++)
	{
		r_prev[i] = b[i] - Ax0[i];
		p_prev[i] = r_prev[i];
		z_prev[i] = r_prev[i];
		s_prev[i] = r_prev[i];
	}
	double residual = 0;
	MatrixMultVector(x_prev, Ax_b);
	VectSub(Ax_b, b, forNorm);
	residual = VectorNorm(forNorm) / VectorNorm(b);

	int iter;
	for (iter = 0; iter < max_iter && residual > eps; iter++)
	{


		MatrixMultVector(z_prev, Az);
		ak = VectorMultVector(p_prev, r_prev) / VectorMultVector(s_prev, Az);

		TransponMultVector(s_prev, ATs);
		for (int j = 0; j < n; j++)
		{
			x_new[j] = x_prev[j] + ak * z_prev[j];
			r_new[j] = r_prev[j] - ak * Az[j];
			p_new[j] = p_prev[j] - ak * ATs[j];
		}
		bk = VectorMultVector(p_new, r_new) / VectorMultVector(p_prev, r_prev);
		for (int j = 0; j < n; j++)
		{
			z_new[j] = r_new[j] + bk * z_prev[j];
			s_new[j] = p_new[j] + bk * s_prev[j];
		}



		MatrixMultVector(x_new, Ax_b);
		VectSub(Ax_b, b, forNorm);
		residual = VectorNorm(forNorm) / VectorNorm(b);
		//if (iter % 1000 == 0)
			cout << iter << " | " << residual << endl;

		z_prev = z_new;
		p_prev = p_new;
		s_prev = s_new;
		r_prev = r_new;
		x_prev = x_new;
	}

	cout << "Iterations: " << iter << endl;
	cout << "Residual: " << residual << endl;
}

void SLAE::MatrixMultVector(vector<double>& v, vector<double>& res)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = di[i] * v[i];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			res[i] += al[j] * v[ja[j]];
			res[ja[j]] += au[j] * v[i];
		}
	}
}

void SLAE::TransponMultVector(vector<double>& v, vector<double>& res)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = di[i] * v[i];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			res[i] += au[j] * v[ja[j]];
			res[ja[j]] += al[j] * v[i];
		}
	}
}

double SLAE::VectorNorm(vector<double>& v)
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += v[i] * v[i];
	}
	return sqrt(res);
}

double SLAE::VectorMultVector(vector<double>& v1, vector<double>& v2)
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += v1[i] * v2[i];
	}
	return res;
}

void SLAE::VectSub(vector<double>& v1, vector<double>& v2, vector<double>& res)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = v1[i] - v2[i];
	}
}

void SLAE::PrintX()
{
	ofstream out;
	out.open("test1/x.txt");

	vector<double> x_real;
	x_real.resize(4394);
	double j = -1;
	for (int i = 0; i < 4394; i++)
	{
		if (i % 25 == 0) j++;
		x_real[i] = j;
	}

	vector<double> tmp;
	tmp.resize(n);
	VectSub(x_real, x_new, tmp);
	double ist_nev = VectorNorm(tmp) / VectorNorm(x_real);

	for (int i = 0; i < n; i++)
	{
		//printf_s("% .15lf\n", x_new[i]);
		out << x_new[i] << endl;
	}
	cout << ist_nev << endl;
	out.close();

}
