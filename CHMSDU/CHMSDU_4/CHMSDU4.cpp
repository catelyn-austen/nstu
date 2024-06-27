#include <iostream>
#include <math.h>
#include <stdio.h>
using namespace std;

double fun(double t, double y)
{
	return 2 * t * y;
}
double kn(double h, double t, double y)
{
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0;
	k1 = fun(t, y);
	k2 = fun(t + h / 2, y + h * k1 / 2);
	k3 = fun(t + h / 2, y + h * k2 / 2);
	k4 = fun(t + h, y + h * k3);
	return 1.0 / 6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}
void explAdams3(double h, FILE* out)
{
	fprintf_s(out, "                                           Explicit Adams 3                                   \n");
	int size = 1 / h;
	double yn = 1, tn = 0;
	if (size < 2)
		return;
	double* fn = new double[3];
	double* arr = new double[size + 1];

	fn[0] = 1;
	for (int i = 1; i < 3; i++)
	{
		arr[i - 1] = fn[i - 1];
		fn[i] = fn[i - 1] + h * kn(h, tn, fn[i - 1]);
		fn[i - 1] = fun(tn, fn[i - 1]);
		tn = i * h;
	}
	arr[2] = fn[2];
	yn = fn[2];
	fn[2] = fun(tn, yn);

	for (int i = 3; i <= size; i++)
	{
		tn = i * h;
		yn += h / 12.0 * (23.0 * fn[2] - 16.0 * fn[1] + 5.0 * fn[0]);
		arr[i] = yn;
		fn[0] = fn[1];
		fn[1] = fn[2];
		fn[2] = fun(tn, yn);
	}

	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr[i];
		yiA = exp(pow(i * h, 2));
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void explAdams4(double h, FILE* out)
{
	fprintf_s(out, "                                           Explicit Adams 4                                   \n");
	int size = 1 / h;
	double yn = 1, tn = 0;
	if (size < 3)
		return;
	double* fn = new double[4];
	double* arr = new double[size + 1];

	fn[0] = 1;
	for (int i = 1; i < 4; i++)
	{
		arr[i - 1] = fn[i - 1];
		fn[i] = fn[i - 1] + h * kn(h, tn, fn[i - 1]);
		fn[i - 1] = fun(tn, fn[i - 1]);
		tn = i * h;
	}
	arr[3] = fn[3];
	yn = fn[3];
	fn[3] = fun(tn, yn);

	for (int i = 4; i <= size; i++)
	{
		tn = i * h;
		yn += h / 24.0 * (55.0 * fn[3] - 59.0 * fn[2] + 37.0 * fn[1] - 9.0 * fn[0]);
		arr[i] = yn;
		fn[0] = fn[1];
		fn[1] = fn[2];
		fn[2] = fn[3];
		fn[3] = fun(tn, yn);
	}

	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr[i];
		yiA = exp(pow(i * h, 2));
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void implAdams3(double h, FILE* out)
{
	fprintf_s(out, "                                           Implicit Adams 3                                   \n");
	int size = 1 / h;
	double yn = 1, tn = 0, eps = 1e-14, y0, y1;
	if (size < 2)
		return;
	double* fn = new double[2];
	double* arr = new double[size + 1];

	arr[0] = fn[0] = 1;
	fn[1] = fn[0] + h * kn(h, tn, fn[0]);
	fn[0] = fun(tn, fn[0]);
	arr[1] = fn[1];
	yn = fn[1];
	tn = h;
	fn[1] = fun(tn, yn);


	for (int i = 2; i <= size; i++)
	{
		tn = i * h;
		y0 = yn + fn[1];
		y1 = y0 + fun(tn, y0);
		while (abs(y1 - y0) > eps)
		{
			y0 = y1;
			y1 = yn + h / 12.0 * (5.0 * fun(tn, y0) + 8.0 * fn[1] - fn[0]);
		}
		yn = y1;
		arr[i] = yn;
		fn[0] = fn[1];
		fn[1] = fun(tn, yn);
	}

	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr[i];
		yiA = exp(pow(i * h, 2));
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void implAdams4(double h, FILE* out)
{
	fprintf_s(out, "                                           Implicit Adams 4                                   \n");
	int size = 1 / h;
	double yn = 1, tn = 0, eps = 1e-14, y0, y1;
	if (size < 3)
		return;
	double* fn = new double[3];
	double* arr = new double[size + 1];

	arr[0] = fn[0] = 1;
	for (int i = 1; i < 3; i++)
	{
		arr[i - 1] = fn[i - 1];
		fn[i] = fn[i - 1] + h * kn(h, tn, fn[i - 1]);
		fn[i - 1] = fun(tn, fn[i - 1]);
		tn = i * h;
	}
	arr[2] = fn[2];
	yn = fn[2];
	fn[2] = fun(tn, yn);


	for (int i = 3; i <= size; i++)
	{
		tn = i * h;
		y0 = yn + fn[2];
		y1 = y0 + fun(tn, y0);
		while (abs(y1 - y0) > eps)
		{
			y0 = y1;
			y1 = yn + h / 24.0 * (9.0 * fun(tn, y0) + 19.0 * fn[2] - 5 * fn[1] + fn[0]);
		}
		yn = y1;
		arr[i] = yn;
		fn[0] = fn[1];
		fn[1] = fn[2];
		fn[2] = fun(tn, yn);
	}

	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr[i];
		yiA = exp(pow(i * h, 2));
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void forecast_corr3(double h, FILE* out)
{
	fprintf_s(out, "                                           Forecast & Correction 3                                   \n");
	int size = 1 / h;
	double yn = 1, tn = 0, yn1;
	if (size < 2)
		return;
	double* fn = new double[3];
	double* arr = new double[size + 1];

	fn[0] = 1;
	for (int i = 1; i < 3; i++)
	{
		arr[i - 1] = fn[i - 1];
		fn[i] = fn[i - 1] + h * kn(h, tn, fn[i - 1]);
		fn[i - 1] = fun(tn, fn[i - 1]);
		tn = i * h;
	}
	arr[2] = fn[2];
	yn = fn[2];
	fn[2] = fun(tn, yn);

	for (int i = 3; i <= size; i++)
	{
		tn = i * h;
		yn1 = yn + h / 12.0 * (23.0 * fn[2] - 16.0 * fn[1] + 5.0 * fn[0]);
		yn += h / 12.0 * (5.0 * fun(tn, yn1) + 8.0 * fn[2] - fn[1]);
		arr[i] = yn;

		fn[0] = fn[1];
		fn[1] = fn[2];
		fn[2] = fun(tn, yn);
	}
	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr[i];
		yiA = exp(pow(i * h, 2));
		fprintf_s(out, "%.3f: %.17f Analyt: %.17f Pogreshnost': %.2e\n", h * i, yiE,
			yiA, abs(yiA - yiE));
	}
	fprintf_s(out, "---------------------------------------------------------------------------------------------\n");
}
void forecast_corr4(double h, FILE* out)
{
	fprintf_s(out, "                                           Forecast & Correction 4                                   \n");
	int size = 1 / h;
	double yn = 1, tn = 0, yn1;
	if (size < 3)
		return;
	double* fn = new double[4];
	double* arr = new double[size + 1];

	fn[0] = 1;
	for (int i = 1; i < 4; i++)
	{
		arr[i - 1] = fn[i - 1];
		fn[i] = fn[i - 1] + h * kn(h, tn, fn[i - 1]);
		fn[i - 1] = fun(tn, fn[i - 1]);
		tn = i * h;
	}
	arr[3] = fn[3];
	yn = fn[3];
	fn[3] = fun(tn, yn);

	for (int i = 4; i <= size; i++)
	{
		tn = i * h;
		yn1 = yn + h / 24.0 * (55.0 * fn[3] - 59.0 * fn[2] + 37.0 * fn[1] - 9.0 * fn[0]);
		yn += h / 24.0 * (9.0 * fun(tn, yn1) + 19.0 * fn[3] - 5.0 * fn[2] + fn[1]);
		arr[i] = yn;

		fn[0] = fn[1];
		fn[1] = fn[2];
		fn[2] = fn[3];
		fn[3] = fun(tn, yn);
	}
	double yiE = 0, yiA = 0;
	for (int i = 0; i < size + 1; i++)
	{
		yiE = arr[i];
		yiA = exp(pow(i * h, 2));
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
	double h1 = 0.1, h2 = 0.05, h3 = 0.025;
	explAdams3(h1, out);
	explAdams3(h2, out);
	explAdams3(h3, out);

	explAdams4(h1, out);
	explAdams4(h2, out);
	explAdams4(h3, out);

	implAdams3(h1, out);
	implAdams3(h2, out);
	implAdams3(h3, out);

	implAdams4(h1, out);
	implAdams4(h2, out);
	implAdams4(h3, out);

	forecast_corr3(h1, out);
	forecast_corr3(h2, out);
	forecast_corr3(h3, out);

	forecast_corr4(h1, out);
	forecast_corr4(h2, out);
	forecast_corr4(h3, out);
	fclose(out);
}