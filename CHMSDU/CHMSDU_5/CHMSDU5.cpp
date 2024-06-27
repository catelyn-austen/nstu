#include <iostream>
#include <math.h>
#include <stdio.h>
using namespace std;

double q0 = 0.001, l = 1.0, d = 0.01, ro = 1000.0, Csnd = 1260, patm = 1e+5;
const double pi = 3.14159265358979323846;
double S = pi * d * d / 4.0;
double z = 1 - cos(pi / 4);

double P(double t, double q)
{
	double qn, C;
	C = (l * S) / (ro * (Csnd * Csnd));

	if (t <= 1)
		qn = q0 * t;
	else
		qn = q0;
	return (qn - q) / C;
}

double Q(double t, double p, double q)
{
	double sign;
	if ((p - patm) > 0)
		sign = 1.0;
	else if ((p - patm) == 0)
		sign = 0.0;
	else
		sign = -1.0;
	return sqrt(z * abs(p - patm) / (2 * ro)) * (S * sqrt(2 * abs(p - patm) / (ro * z)) * sign - q);
}


void Direct(double h, FILE* out)
{
	//fprintf_s(out, "                                           Direct                                     \n");
	int size = 20 / h;
	double* arr_e1 = new double[size + 2];
	double* arr_e2 = new double[size + 2];
	double pn = 0, tn = 0, pkn = 0, pk1 = 0, pk2 = 0, pk3 = 0, pk4 = 0;
	double qn = 0, qkn = 0, qk1 = 0, qk2 = 0, qk3 = 0, qk4 = 0;
	arr_e1[0] = patm;
	arr_e2[0] = 0;

	for (int i = 0; i < size + 1; i++)
	{
		tn = h * i;
		pn = arr_e1[i];
		qn = arr_e2[i];
		pk1 = P(tn, qn);
		qk1 = Q(tn, pn, qn);

		pk2 = P(tn + h / 2, qn + h / 2 * qk1);
		qk2 = Q(tn + h / 2, pn + h / 2 * pk1, qn + h / 2 * qk1);

		pk3 = P(tn + h / 2, qn + h / 2 * qk2);
		qk3 = Q(tn + h / 2, pn + h / 2 * pk2, qn + h / 2 * qk2);

		pk4 = P(tn + h, qn + h * qk3);
		qk4 = Q(tn + h, pn + h * pk3, qn + h * qk3);

		pkn = (pk1 + 2 * pk2 + 2 * pk3 + pk4) / 6.;
		qkn = (qk1 + 2 * qk2 + 2 * qk3 + qk4) / 6.;
		arr_e1[i + 1] = pn + h * pkn;
		arr_e2[i + 1] = qn + h * qkn;
	}

	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i+=10000)
	{
		//yiE = arr_e1[i];
		yiA = arr_e2[i];
		fprintf_s(out, "%.15f\n", yiA);
	}
	//fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}

void Sequantial(double h, FILE* out)
{
	//fprintf_s(out, "                                           Sequantial                                     \n");
	int size = 20 / h;
	double* arr_e1 = new double[size + 2];
	double* arr_e2 = new double[size + 2];
	double pn = 0, tn = 0, pkn = 0, pk1 = 0, pk2 = 0, pk3 = 0, pk4 = 0;
	double qn = 0, qkn = 0, qk1 = 0, qk2 = 0, qk3 = 0, qk4 = 0;
	arr_e1[0] = patm;
	arr_e2[0] = 0;

	for (int i = 0; i < size + 1; i++)
	{
		tn = h * i;
		pn = arr_e1[i];
		qn = arr_e2[i];

		pk1 = P(tn, qn);
		pk2 = P(tn + h / 2, qn + h / 2 * qn);
		pk3 = P(tn + h / 2, qn + h / 2 * qn);
		pk4 = P(tn + h, qn + h * qn);
		pkn = (pk1 + 2 * pk2 + 2 * pk3 + pk4) / 6.;
		arr_e1[i + 1] = pn + h * pkn;

		pn = arr_e1[i + 1];
		qk1 = Q(tn, pn, qn);
		qk2 = Q(tn + h / 2, pn, qn + h / 2 * qk1);
		qk3 = Q(tn + h / 2, pn, qn + h / 2 * qk2);
		qk4 = Q(tn + h, pn, qn + h * qk3);
		qkn = (qk1 + 2 * qk2 + 2 * qk3 + qk4) / 6.;
		arr_e2[i + 1] = qn + h * qkn;
	}

	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i+=10000)
	{
		
		//yiE = arr_e1[i];
		yiA = arr_e2[i];
		fprintf_s(out, "%.17f\n", yiA);
	}
	//fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}

int main()
{
	FILE* out;
	if (fopen_s(&out, "out.txt", "a"))
		return 1;
	double h1 = 1e-4, h2 = 1e-5;
	//Direct(h1, out);
	//Direct(h2, out);
	//Sequantial(h1, out);
	Sequantial(h2, out);
	fclose(out);
}
