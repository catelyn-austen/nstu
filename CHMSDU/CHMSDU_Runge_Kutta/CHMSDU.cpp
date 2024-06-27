#include <iostream>
#include <math.h>
#include <stdio.h>
using namespace std;

double fun1(double t, double y)
{
	return 2 * t * y;
}
double fun2(double t, double y)
{
	return -25 * y + cos(t) + 25 * sin(t);
}
void explEuler(double h, FILE* out)
{
	fprintf_s(out, "                                           Explicit Euler                                    \n");
	int size = 2 / h;
	double* arr_e1 = new double[size + 1];
	double yn = 0, tn = 0;
	arr_e1[0] = 1;
	for (int i = 0; i < size; i++)
	{
		tn = h * i;
		yn = arr_e1[i];
		arr_e1[i + 1] = yn + h * fun2(tn, yn);
	}
	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr_e1[i];
		yiA = exp(-25 * i * h) + sin(i * h);
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void modEuler(double h, FILE* out)
{
	fprintf_s(out, "                                           Modified Euler                                    \n");
	int size = 2 / h;
	double* arr_e1 = new double[size + 1];
	double yn = 0, tn = 0, tn1 = 0;
	arr_e1[0] = 1;
	for (int i = 0; i < size; i++)
	{
		tn = h * i;
		tn1 = h * (i + 1);
		yn = arr_e1[i];
		arr_e1[i + 1] = yn + h / 2. * (fun2(tn, yn) + fun2(tn1, yn + h * fun2(tn, yn)));
	}
	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr_e1[i];
		yiA = exp(-25 * i * h) + sin(i * h);
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void implEuler(double h, FILE* out)
{
	fprintf_s(out, "                                           Implicit Euler                                    \n");
	double eps = 1e-14;
	double yn = 1, tn = 0;
	double an0 = 0, an1 = 0;
	double yiE = 0, yiA = 0;
	an0 = yn + h * fun2(tn, yn);
	an1 = an0 - (yn + h * fun2(tn, an0) - an0) / ((-25 * h) - 1);
	for (tn = 0; tn <= 2 + h/2; tn += h)
	{
		yiE = yn;
		yiA = exp(-25 * tn) + sin(tn);
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n",
			tn, yiE, yiA, abs(yiA - yiE));
		do {
			an0 = an1;
			an1 = an0 - (yn + h * fun2(tn, an0) - an0) / ((-25 * h) - 1);
		} while (abs(an1 - an0) > eps);
		yn = an1;
	}
	fprintf_s(out,
		"--------------------------------------------------------------------------------------------\n");
}
void Trapez(double h, FILE* out) 
{
	fprintf_s(out, "                                           Trapezoid                              \n");
	double eps = 1e-14;
	double yn = 1, tn = 0;
	double an0 = 0, an1 = 0;
	double yiE = 0, yiA = 0;
	an0 = yn + h * fun2(tn, yn);
	an1 = an0 - (yn + (h / 2) * (fun2(tn, yn) + fun2(tn + h, an0)) - an0) / ((-25 * h) - 1);
	for (tn = 0; tn <= 2 + eps; tn += h)
	{
		yiE = yn;
		yiA = exp(-25 * tn) + sin(tn);
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n",
			tn, yiE, yiA, abs(yiA - yiE));
		do {
			an0 = an1;
			an1 = an0 - (yn + (h / 2) * (fun2(tn, yn) + fun2(tn + h, an0)) - an0) / ((-25 * h) - 1);
		} while (abs(an1 - an0) > eps);
		yn = an1;
	}
	fprintf_s(out,
		"--------------------------------------------------------------------------------------------\n");
}
void RungeKutta1(double h, FILE* out)
{
	fprintf_s(out, "                                           Runge-Kutta                                     \n");
	int size = 1 / h;
	double* arr_e1 = new double[size + 1];
	double yn = 0, tn = 0, kn = 0, k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	arr_e1[0] = 1;
	for (int i = 0; i < size + 1; i++)
	{
		tn = h * i;
		yn = arr_e1[i];
		k1 = fun1(tn, yn);
		k2 = fun1(tn + h / 2, yn + h / 2 * k1);
		k3 = fun1(tn + h / 2, yn + h / 2 * k2);
		k4 = fun1(tn + h, yn + h * k3);
		kn = (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		arr_e1[i + 1] = yn + h * kn;
	}
	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr_e1[i];
		yiA = exp(h * h * i * i);
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void RungeKutta2(double h, FILE* out)
{
	fprintf_s(out, "                                           Runge-Kutta                                     \n");
	int size = 2 / h;
	double* arr_e1 = new double[size + 1];
	double yn = 0, tn = 0, kn = 0, k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	arr_e1[0] = 1;
	for (int i = 0; i < size + 1; i++)
	{
		tn = h * i;
		yn = arr_e1[i];
		k1 = fun2(tn, yn);
		k2 = fun2(tn + h / 2, yn + h / 2 * k1);
		k3 = fun2(tn + h / 2, yn + h / 2 * k2);
		k4 = fun2(tn + h, yn + h * k3);
		kn = (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
		arr_e1[i + 1] = yn + h * kn;
	}
	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr_e1[i];
		yiA = exp(-25 * i * h) + sin(i * h);
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
int main()
{
	FILE* out;
	if (fopen_s(&out, "out.txt", "a"))
		return 1;
	double h1 = 0.2, h2 = 0.1, h3 = 0.05, h4 = 0.025, h5 = 0.5, h6 = 0.25, h7 = 0.125;
	RungeKutta1(h2, out);
	RungeKutta1(h3, out);
	RungeKutta1(h4, out);
	explEuler(h2, out);
	explEuler(h3, out);
	explEuler(h4, out);
	implEuler(h2, out);
	implEuler(h3, out);
	implEuler(h4, out);
	modEuler(h1, out);
	modEuler(h2, out);
	modEuler(h3, out);
	Trapez(h1, out);
	Trapez(h2, out);
	Trapez(h3, out);
	RungeKutta2(h5, out);
	RungeKutta2(h6, out);
	RungeKutta2(h7, out);
	RungeKutta2(h2, out);
	fclose(out);
}