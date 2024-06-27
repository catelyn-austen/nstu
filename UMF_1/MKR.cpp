#include "stdafx.h"
#include "MKR.h"
#include "matrix.h"

double MKR::Func(int i, int j)
{
	double curX = x[j];
	double curY = y[i];
	switch (variant)
	{
	case 1:
		return curX + curY;
	case 2:
		return -4 + curX * curX + curY * curY;
	case 3:
		return -6 * (curX + curY) + pow(curX, 3) + pow(curY, 3);
	case 4:
		return -12 * (curX * curX + curY * curY) + pow(curX, 4) + pow(curY, 4);
	case 5:
		return 9 * cos(2 * curX + 2 * curY);
	default:
		break;
	}
	return NAN;
}

double MKR::U(int i, int j)
{
	double curX = x[j];
	double curY = y[i];
	switch (variant)
	{
	case 1:
		return curX + curY;
	case 2:
		return curX * curX + curY * curY;
	case 3:
		return pow(curX, 3) + pow(curY, 3);
	case 4:
		return pow(curX, 4) + pow(curY, 4);
	case 5:
		return cos(2 * curX + 2 * curY);
	default:
		break;
	}
	return NAN;
}

double MKR::Theta(int i, int j, int nx, int ny)
{
	double curX = x[j];
	double curY = y[i];
	switch (variant)
	{
	case 1:
		return 1 * nx + 1 * ny;
	case 2:
		return 2 * curX * nx + 2 * curY * ny;
	case 3:
		return 3 * curX * curX * nx + 3 * curY * curY * ny;
	case 4:
		return 4 * pow(curX, 3) * nx + 4 * pow(curY, 3) * ny;
	case 5:
		return -2 * sin(2 * curX + 2 * curY) * nx + -2 * sin(2 * curX + 2 * curY) * ny;
	default:
		break;
	}
	return NAN;
}

void MKR::Input(FILE* inParams, FILE* inPoints, FILE* inBoundary)
{
	fscanf_s(inParams, "%lf %lf %d", &lambda, &gamma, &variant);

	double q, ax, bx, ay, by;
	int curXi, curYj;
	int lvlN, numOfX, numOfY, n, m;
	total_n = 0, total_m = 0;
	if (inPoints)
	{
		fscanf_s(inPoints, "%d", &lvlN); // Вложенность сетки
		fscanf_s(inPoints, "%d", &numOfX);
		fscanf_s(inPoints, "%d", &numOfY);
		for (int i = 0; i < numOfX; i++)
		{
			fscanf_s(inPoints, "%lf", &q);
			if (q == 1)
			{
				fscanf_s(inPoints, "%lf %d %lf", &ax, &n, &bx);
				total_n += n + 1; // так как n - интервалы
				int NinNested = n;
				for (int i = 1; i < lvlN; i++)
					NinNested *= 2;

				double hx = (bx - ax) / NinNested;
				if (x.empty()) // мне здесь не нравится
					x.push_back(ax);
				else if (ax != x.back())
					x.push_back(ax);
				else if (ax == x.back())
					total_n -= 1; /// ахахахах ну и намудрил я. Это если даны интервалы, (2, 4) (4, 8), то надо убрать один узел, если (2, 4) (5, 8) то не надо

				for (int j = 1; j < NinNested; j++)
					x.push_back(ax + hx * j);
				x.push_back(bx);
			}
			else
			{

				fscanf_s(inPoints, "%lf %d %lf", &ax, &n, &bx);
				total_n += n + 1; // так как n - интервалы
				int NinNested = n;
				for (int i = 1; i < lvlN; i++)
					NinNested *= 2;

				for (int j = 1; j < lvlN; j++)
					q = sqrt(q);

				double h0 = (bx - ax) / (powl(q, NinNested) - 1) * (q - 1);

				if (x.empty()) // мне здесь не нравится
					x.push_back(ax);
				else if (ax != x.back())
					x.push_back(ax);
				else if (ax == x.back())
					total_n -= 1;

				for (int j = 1; j < NinNested; j++) {
					x.push_back(x.back() + h0 * pow(q, j - 1)); // Вот тут не уверен, что правильно
				}

				x.push_back(bx);
			}
		}

		for (int i = 0; i < numOfY; i++)
		{
			fscanf_s(inPoints, "%lf", &q);
			if (q == 1)
			{
				fscanf_s(inPoints, "%lf %d %lf", &ay, &m, &by);
				total_m += m + 1;  // так как m - интервалы
				int MinNested = m;
				for (int i = 1; i < lvlN; i++)
					MinNested *= 2;

				double hy = (by - ay) / MinNested;
				if (y.empty()) // мне здесь не нравится
					y.push_back(ay);
				else if (ay != y.back())
					y.push_back(ay);
				else if (ay == y.back())
					total_m -= 1;

				for (int j = 1; j < MinNested; j++)
				{
					y.push_back(ay + hy * j);
				}
				y.push_back(by);
			}
			else
			{
				fscanf_s(inPoints, "%lf %d %lf", &ay, &m, &by);
				total_m += m + 1;  // так как m - интервалы
				int MinNested = m;
				for (int i = 1; i < lvlN; i++)
					MinNested *= 2;

				for (int j = 1; j < lvlN; j++)
					q = sqrtl(q);

				double h0 = (by - ay) / (powl(q, MinNested) - 1) * (q - 1);

				if (y.empty()) // мне здесь не нравится
					y.push_back(ay);
				else if (ay != y.back())
					y.push_back(ay);
				else if (ay == y.back())
					total_m -= 1;

				for (int j = 1; j < MinNested; j++) {
					y.push_back(y.back() + h0 * pow(q, j - 1)); // Вот тут не уверен, что правильно
				}
				y.push_back(by);

			}
		}

	}
	else
		throw new exception();

	for (int i = 1; i < lvlN; i++) // Смотря какая вложенность сетки
	{
		total_n += total_n - 1;
		total_m += total_m - 1;
	}
	int num = total_n * total_m; 
	AllocateMemory(num);

	fscanf_s(inPoints, "%d %d", &curXi, &curYj); // считывание углового узла L
	k1 = curYj;
	k2 = curXi;
	for (int i = 1; i < lvlN; i++)
	{
		k1 += k1 - 1;
		k2 += k2 - 1;
	}


	BoundaryCondition bound;
	for (int i = 0; i < 4; i++)
	{
		fscanf_s(inBoundary, "%d %c", &bound.type, &bound.orientation, 1);
		if (bound.type > 1)
			bound.priority = 1;
		else
			bound.priority = 2;
		boundaryConditions.push(bound);
	}
}

void MKR::Solve()
{
	SLAE slae;
	MainPass();
	BoundaryPass();
	slae.Input(A, q, f, total_n * total_m, total_n - 2, 1e-15, 10000);
	slae.OutputDense();
	slae.IterativeMethod(1);
	q = slae.x;
}

void MKR::AllocateMemory(int totalPoints)
{
	A = new double* [totalPoints];
	for (int i = 0; i < totalPoints; i++)
	{
		A[i] = new double[5]();
	}
	q = new double[totalPoints]();
	f = new double[totalPoints]();
}

void MKR::MainPass()
{
	for (int i = 0; i < k1; i++)
	{
		for (int j = 0; j < total_n; j++)
		{
			// по x
			AddDerivativeX(A, i, j, total_n - 1);
			// по y
			AddDerivativeY(A, i, j, j < k2 ? total_m - 1 : k1 - 1);
			// учет гаммы
			A[GetPointNum(i, j)][2] += gamma;
			// Правая часть
			f[GetPointNum(i, j)] = Func(i, j);
		}
	}
	for (int i = k1; i < total_m; i++)
	{
		for (int j = 0; j < k2; j++)
		{
			// по x
			AddDerivativeX(A, i, j, k2 - 1);
			// по y
			AddDerivativeY(A, i, j, total_m - 1);
			// учет гаммы
			A[GetPointNum(i, j)][2] += gamma;
			// Правая часть
			f[GetPointNum(i, j)] = Func(i, j);
		}
		for (int k = k2; k < total_n; k++) {
			A[GetPointNum(i, k)][2] = 1;
			f[GetPointNum(i, k)] = 0;
		}
	}
}

void MKR::BoundaryPass()
{
	for (int i = 0; i < 4; i++)
	{
		BoundaryCondition current = boundaryConditions.top();
		boundaryConditions.pop();
		if (current.type == 1) {
			switch (current.orientation)
			{
			case 'l':
				for (int i = 0; i < total_m; i++)
				{
					int pointNum = GetPointNum(i, 0);
					for (int k = 0; k < 5; k++)
						A[pointNum][k] = 0;
					A[pointNum][2] = 1;
					f[pointNum] = U(i, 0);
				}
				break;
			case 'u':
				for (int j = 0; j < k2; j++)
				{
					int pointNum = GetPointNum(total_m - 1, j);
					for (int k = 0; k < 5; k++)
						A[pointNum][k] = 0;
					A[pointNum][2] = 1;
					f[pointNum] = U(total_m - 1, j);
				}
				for (int j = k2 - 1; j < total_n; j++)
				{
					int pointNum = GetPointNum(k1 - 1, j);
					for (int k = 0; k < 5; k++)
						A[pointNum][k] = 0;
					A[pointNum][2] = 1;
					f[pointNum] = U(k1 - 1, j);
				}
				break;
			case 'r':
				for (int i = 0; i < k1; i++)
				{
					int pointNum = GetPointNum(i, total_n - 1);
					for (int k = 0; k < 5; k++)
						A[pointNum][k] = 0;
					A[pointNum][2] = 1;
					f[pointNum] = U(i, total_n - 1);
				}
				for (int i = k1; i < total_m; i++)
				{
					int pointNum = GetPointNum(i, k2 - 1);
					for (int k = 0; k < 5; k++)
						A[pointNum][k] = 0;
					A[pointNum][2] = 1;
					f[pointNum] = U(i, k2 - 1);
				}
				break;
			case 'd':
				for (int j = 0; j < total_n; j++)
				{
					int pointNum = GetPointNum(0, j);
					for (int k = 0; k < 5; k++)
						A[pointNum][k] = 0;
					A[pointNum][2] = 1;
					f[pointNum] = U(0, j);
				}
				break;
			default:
				printf_s("Invalid boundary condition orientation");
				break;
			}
		}
		else if (current.type == 2) {
			//так же сделать так, чтобы в тетту x и y передавался
			switch (current.orientation)
			{
			case 'l':
				hx = x[1] - x[0];
				for (int i = 0; i < total_m; i++)
				{
				    int pointNum = GetPointNum(i, 0);
					A[pointNum][2] = -lambda / hx;
					A[pointNum][3] = lambda / hx;
					f[pointNum] = Theta(i, 0, -1, 0);
				}
				break;
			case 'u':
				hy = y[total_m - 1] - y[total_m - 2];
				for (int j = 0; j < k2; j++)
				{
					int pointNum = GetPointNum(total_m - 1, j);
					A[pointNum][2] = lambda / hy;
					A[pointNum][0] = -lambda / hy;
					f[pointNum] = Theta(total_m - 1, j, 0, 1);
				}
				hy = y[k1 - 1] - y[k1 - 2];
				for (int j = k2; j < total_n; j++)
				{
					int pointNum = GetPointNum(k1 - 1, j);
					A[pointNum][2] = lambda / hy;
					A[pointNum][0] = -lambda / hy;
					f[pointNum] = Theta(k1 - 1, j, 0, 1);
				}
				break;
			case 'r':
				hx = x[total_n - 1] - x[total_n - 2];
				for (int i = 0; i < k1; i++)
				{
					int pointNum = GetPointNum(i, total_n - 1);
					A[pointNum][2] = lambda / hx;
					A[pointNum][1] = -lambda / hx;
					f[pointNum] = Theta(i, total_n - 1, 1, 0);
				}
				hx = x[k2 - 1] - x[k2 - 2];
				for (int i = k1; i < total_m; i++)
				{
					int pointNum = GetPointNum(i, k2 - 1);
					A[pointNum][2] = lambda / hx;
					A[pointNum][1] = -lambda / hx;
					f[pointNum] = Theta(i, k2 - 1, 1, 0);
				}
				break;
			case 'd':
				hy = y[1] - y[0];
				for (int j = 0; j < total_n; j++)
				{
					int pointNum = GetPointNum(0, j);
					A[pointNum][2] = -lambda / hy;
					A[pointNum][4] = lambda / hy;
					f[pointNum] = Theta(0, j, 0, -1);
				}
				break;
			default:
				printf_s("Invalid boundary condition orientation");
				break;
			}
		}
		else
			printf_s("Invalid boundary condition type");
	}
}

void MKR::AddDerivativeX(double** &A, int i, int j, int max)
{
	int num = GetPointNum(i, j);
	if (j == 0) {
		hx = x[1] - x[0];
		A[num][2] += lambda / hx;
		A[num][3] += -lambda / hx;
	}
	else if (j == max) {
		hx = x[total_n - 1] - x[total_n - 2];
		A[num][2] += -lambda / hx;
		A[num][1] += lambda / hx;
	}
	else {
		double hx_i_1 = x[j] - x[j - 1]; // По формулам i, а дефактом j у нас отвечает за x
		double hx_i = x[j + 1] - x[j];
		A[num][1] += -2 * lambda / (hx_i_1 * (hx_i + hx_i_1));
		A[num][3] += -2 * lambda / (hx_i * (hx_i + hx_i_1));
		A[num][2] += 2 * lambda / (hx_i_1 * hx_i);
	}
}

void MKR::AddDerivativeY(double** &A, int i, int j, int max)
{
	int num = GetPointNum(i, j);
	if (i == 0) {
		hy = y[1] - y[0];
		A[num][2] += lambda / hy;
		A[num][4] += -lambda / hy;
	}
	else if (i == max) {
		hy = y[total_m - 1] - y[total_m - 2];
		A[num][2] += -lambda / hy;
		A[num][0] += lambda / hy;
	}
	else {
		double hy_j = y[i + 1] - y[i];
		double hy_j_1 = y[i] - y[i - 1];
		A[num][0] += -2 * lambda / (hy_j_1 * (hy_j + hy_j_1));
		A[num][4] += -2 * lambda / (hy_j * (hy_j + hy_j_1));
		A[num][2] += 2 * lambda / (hy_j * hy_j_1);
	}
}

int MKR::GetPointNum(int i, int j) const
{
	return total_n * i + j;
}


void MKR::OutputSolutionVector()
{
	for (int i = 0; i < total_n*total_m; i++)
	{
		printf("%.7lf\n", q[i]);
	}
}

void MKR::PrintTrue()
{
	for (int i = 0; i < total_m; i++)
	{
		for (int j = 0; j < total_n; j++)
		{
			printf_s("%.15lf\n", U(i, j));
		}
	}
}

void MKR::PrintDes()
{
	for (int i = 0; i < total_m; i++)
	{
		for (int j = 0; j < total_n; j++)
		{
			printf_s("(%.15lf, %.15lf)\n", x[j], y[i]);
		}
	}
}
