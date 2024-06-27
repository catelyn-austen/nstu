#include "FEM.h"

using namespace std;

/*----- Функции и константы ------*/

double FEM::Uf(int n, double x, double t)
{
	switch (n)
	{
	case 1:
		return x;
		break;
	case 2:
		return x * x;
		break;
	case 3:
		return cos(x);
		break;

		// U = x
	case 4:
		return x;
		break;
	case 5:
		return x;
		break;
	case 6:
		return x;
		break;
	case 7:
		return x;
		break;

		// U = x*x
	case 8:
		return x * x;
		break;
	case 9:
		return x * x;
		break;
	case 10:
		return x * x;
		break;
	case 11:
		return x * x;
		break;

		// U = sin(x)
	case 12:
		return sin(x);
		break;
	case 13:
		return sin(x);
		break;
	case 14:
		return sin(x);
		break;
	case 15:
		return sin(x);
		break;

		// U = exp(x)
	case 16:
		return exp(x);
		break;
	case 17:
		return exp(x);
		break;
	case 18:
		return exp(x);
		break;
	case 19:
		return exp(x);
		break;
		// U = cos(x)
	case 20:
		return cos(x);
		break;
	case 21:
		return cos(x);
		break;
	case 22:
		return cos(x);
		break;

	default:
		fprintf(stderr, "UfunStat err");
		break;
	}
}
double FEM::f(int n, double x, double t)
{
	double tmp1;
	double tmp2;
	double tmp3;
	switch (n)
	{
	case 1:
		return x;
		break;
	case 2:
		return -2 + x * x;
		break;
	case 3:
		return 2 * cos(x);
		break;

		// U = x
	case 4:
		return -1 + x;
		break;
	case 5:
		return -2 * x + x;
		break;
	case 6:
		return -exp(x) + x;
		break;
	case 7:
		return sin(x) + x;
		break;

		// U = x*x
	case 8:
		return -3 * 2 * x * x + x * x;
		break;
	case 9:
		return -2 * 5 * pow(x, 4) + x * x;
		break;
	case 10:
		return -(4 * x * x + 2) * exp(x * x) + x * x;
		break;
	case 11:
		return -(2 * cos(x * x) - 4 * x * x * sin(x * x)) + x * x;
		break;

		// U = sin(x)
	case 12:
		return -(cos(x) * cos(x) - sin(x) * sin(x)) + sin(x);
		break;
	case 13:
		return -(2 * cos(x) * cos(x) * sin(x) - sin(x) * sin(x) * sin(x)) + sin(x);
		break;
	case 14:
		return -(exp(sin(x)) * cos(x) * cos(x) - exp(sin(x)) * sin(x)) + sin(x);
		break;
	case 15:
		return -(-cos(x) * cos(x) * sin(sin(x)) - sin(x) * cos(sin(x))) + sin(x);
		break;

		// U = exp(x)
	case 16:
		return -2 * exp(2 * x) + exp(x);
		break;
	case 17:
		return  -3 * exp(3 * x) + exp(x);
		break;
	case 18:
		return -(exp(2 * x) + exp(x)) * exp(exp(x)) + exp(x);
		break;
	case 19:
		return -(exp(x) * cos(exp(x)) - exp(2 * x) * sin(exp(x))) + exp(x);
		break;
		// U = cos(x)
	case 20:
		return cos(x) + cos(2 * x);
		break;
	case 21:
		return cos(x) + pow(cos(x), 3) - 2 * cos(x) * sin(x) * sin(x);
		break;
	case 22:
		return cos(x) + cos(x) * cos(cos(x)) + sin(cos(x)) * sin(x) * sin(x);
		break;
	default:
		throw "Incorrect f";
		break;
	}
}
double FEM::lambda(int n, double qh)
{
	switch (n)
	{
	case 1:
		return 1;
		break;
	case 2:
		return 1;
		break;
	case 3:
		return 1;
		break;
		// U = x
	case 4:
		return qh;
		break;
	case 5:
		return qh * qh;
		break;
	case 6:
		return exp(qh);
		break;
	case 7:
		return cos(qh);
		break;

		// U = x*x
	case 8:
		return qh;
		break;
	case 9:
		return qh * qh;
		break;
	case 10:
		return exp(qh);
		break;
	case 11:
		return cos(qh);
		break;

		// U = sin(x)
	case 12:
		return qh;
		break;
	case 13:
		return qh * qh;
		break;
	case 14:
		return exp(qh);
		break;
	case 15:
		return cos(qh);
		break;

		// U = exp(x)
	case 16:
		return qh;
		break;
	case 17:
		return qh * qh;
		break;
	case 18:
		return exp(qh);
		break;
	case 19:
		return cos(qh);
		break;
		// U = cos(x)
	case 20:
		return qh;
		break;
	case 21:
		return qh * qh;
		break;
	case 22:
		return cos(qh);
		break;
	default:
		throw "Incorrect lambda";
		break;
	}
}
double FEM::gamma(int n, double x)
{
	switch (n)
	{
	case 1:
		return 1;
		break;
	case 2:
		return 1;
		break;
	case 3:
		return 1;
		break;
	case 4:
		return 1;
		break;
	case 5:
		return 1;
		break;
	case 6:
		return 1;
		break;
	case 7:
		return 1;
		break;
	case 8:
		return 1;
		break;
	case 9:
		return 1;
		break;
	case 10:
		return 1;
		break;
	case 11:
		return 1;
		break;
	case 12:
		return 1;
		break;
	case 13:
		return 1;
		break;
	case 14:
		return 1;
		break;
	case 15:
		return 1;
		break;
	case 16:
		return 1;
		break;
	case 17:
		return 1;
		break;
	case 18:
		return 1;
		break;
	case 19:
		return 1;
		break;
	case 20:
		return 1;
		break;
	case 21:
		return 1;
		break;
	case 22:
		return 1;
		break;
	default:
		throw "Incorrect gamma";
		break;
	}
}

/*----- Чтение файлов и заполнение массивов ------*/

void FEM::Input(FILE* Params, FILE* Points, FILE* Bounds)
{
	double q, ax, bx; // коэф. растяжения, начало и конец отрезка
	int split, areas, intervals;  // кол-во дроблений отрезка, его области и кол-во интервалов на отрезке
	total_n = 0;
	fscanf_s(Points, "%d", &split);
	fscanf_s(Points, "%d", &areas);

	// Просматриваем каждую область по Х
	for (int i = 0; i < areas; i++)
	{
		fscanf_s(Points, "%lf", &q);
		if (q == 1) // равномерная сетка
		{
			fscanf_s(Points, "%lf %d %lf", &ax, &intervals, &bx);
			total_n += intervals + 1;
			int nest = intervals;
			for (int i = 1; i < split; i++)
				nest *= 2;
			double hx = (bx - ax) / nest;
			if (x_grid.empty())
				x_grid.push_back(ax);
			else if (ax != x_grid.back())
				x_grid.push_back(ax);
			else if (ax == x_grid.back())
				total_n -= 1;
			for (int j = 1; j < nest; j++) // заполняем сетку
				x_grid.push_back(ax + hx * j);
			x_grid.push_back(bx);
		}
		else // неравномерная сетка
		{
			fscanf_s(Points, "%lf %d %lf", &ax, &intervals, &bx);
			total_n += intervals + 1;
			int nest = intervals;
			for (int i = 1; i < split; i++)
				nest *= 2;
			for (int j = 1; j < split; j++)
				q = sqrtl(q);
			double h0 = (bx - ax) / (1 - powl(q, nest)) * (1 - q);
			if (x_grid.empty())
				x_grid.push_back(ax);
			else if (ax != x_grid.back())
				x_grid.push_back(ax);
			else if (ax == x_grid.back())
				total_n -= 1;
			for (int j = 1; j < nest; j++)
			{
				double tmp1 = h0 * pow(q, j - 1);
				double tmp2 = x_grid.back() + tmp1;
				x_grid.push_back(tmp2);
			}
			x_grid.push_back(bx);
		}
	}
	for (int i = 1; i < split; i++)
		total_n += total_n - 1;

	fscanf_s(Params, "%d", &test);
	/*double t_q, t0, tn;
	int t_split, t_areas, t_intervals;
	total_n_t = 0;
	fscanf_s(Time, "%d", &t_split);
	fscanf_s(Time, "%d", &t_areas);

	// Просматриваем каждую область по t
	for (int i = 0; i < t_areas; i++)
	{
		fscanf_s(Time, "%lf", &t_q);
		if (t_q == 1) // равномерная сетка
		{
			fscanf_s(Time, "%lf %d %lf", &t0, &t_intervals, &tn);
			total_n_t += t_intervals + 1;
			int nest = t_intervals;
			for (int i = 1; i < t_split; i++)
				nest *= 2;
			double ht = (tn - t0) / nest;
			if (t.empty())
				t.push_back(t0);
			else if (t0 != t.back())
				t.push_back(t0);
			else if (tn == t.back())
				total_n_t -= 1;
			for (int j = 1; j < nest; j++)
				t.push_back(t0 + ht * j);
			t.push_back(tn);
		}
		else // неравномерная сетка
		{
			fscanf_s(Time, "%lf %d %lf", &t0, &t_intervals, &tn);
			total_n_t += t_intervals + 1;
			int nest = t_intervals;
			for (int i = 1; i < t_split; i++)
				nest *= 2;
			for (int j = 1; j < t_split; j++)
				t_q = sqrtl(t_q);
			double h0 = (tn - t0) / (powl(t_q, nest) - 1) * (t_q - 1);
			if (t.empty())
				t.push_back(t0);
			else if (t0 != t.back())
				t.push_back(t0);
			else if (t0 == t.back())
				total_n_t -= 1;
			for (int j = 1; j < nest; j++) {
				t.push_back(t.back() + h0 * pow(t_q, j - 1));
			}
			t.push_back(tn);
		}
	}
	for (int i = 1; i < t_split; i++)
		total_n_t += total_n_t - 1;

	time_layers = t.size() - 1;
	*/
	// выделяем память для всех векторов
	Allocate_Memory(total_n);

	// заполняем pi и 1ый временной слой
	for (int i = 0; i < total_n; i++)
	{
		p_i_prev[i] = Uf(test, x_grid[0], 0);
		//q_layers[0][i] = p_i_prev[i];
	}

	// заполнение очереди с краевыми условиями
	BoundCond bounds;
	for (int i = 0; i < 2; i++)
	{
		fscanf_s(Bounds, "%d %c", &bounds.type, &bounds.orientation, 1);
		if (bounds.type > 1)
			bounds.priority = 1;
		else
			bounds.priority = 2;
		boundcond.push(bounds);
	}
}
void FEM::Allocate_Memory(int n)
{
	p_i_prev.resize(n);
	b.resize(n);
	y.resize(n);
	tmp.resize(n);
	p_i.resize(n);
	di.resize(n);
	diLU.resize(n);
	al.resize(n - 1);
	au.resize(n - 1);
	alLU.resize(n - 1);
	auLU.resize(n - 1);
	//iter.resize(time_layers);
	alDeriv.resize(n - 1);
	auDeriv.resize(n - 1);
	diDeriv.resize(n);
	bDeriv.resize(n);
	alLine.resize(n - 1);
	auLine.resize(n - 1);
	diLine.resize(n);
	bLine.resize(n);
	//q_layers.resize(time_layers + 1);
	/*for (int i = 0; i <= time_layers; i++)
	{
		q_layers[i].resize(total_n);
	}*/
	ia.resize(n + 1);
	ia[0] = 0;
	ia[1] = 0;
	for (int k = 1, i = 2; i <= n; i++, k++)
	{
		ia[i] = k;
	}
}

/*----- Решение задачи ------*/

void FEM::Relaxation(vector<double>& v1, vector<double>& v2, vector<double>& v3, double w)
{
	for (int i = 0; i < v1.size(); i++)
	{
		v3[i] = w * v1[i] + (1. - w) * v2[i];
	}
}
double FEM::NumDer(int test, double x) // численная производная
{
	double h = 1e-4;
	double doubleplus = lambda(test, x + 2. * h);
	double doubleminus = lambda(test, x - 2. * h);
	double plus = lambda(test, x + h);
	double minus = lambda(test, x - h);
	return (-doubleplus + 8. * plus - 8. * minus + doubleminus) / (12. * h);
}
void FEM::Solve()
{
	double NewRelD, PrevRelD, w;
	int iter = 0;
	Main_Part();
	Condition_account();

	NewRelD = Rel_discrepancy();
	PrevRelD = NewRelD;
	vector<double> temporary;
	vector<double> p_i_prev_iter; // q_prev
	temporary.resize(total_n);
	p_i_prev_iter.resize(total_n);

	// основная часть
	for (; iter <= max_iter && PrevRelD > eps; iter++)
	{
		LU_decompose();
		CalculateY();
		CalculateX();
		Copy_vector(p_i_prev, p_i_prev_iter);

		int i = 10; // пробуем каждый коэф. релаксации
		do {
			w = 0.1 * i;
			Relaxation(p_i, p_i_prev_iter, temporary, w);
			Copy_vector(temporary, p_i_prev);
			Clearing(); // на всякий случай
			Main_Part();
			Condition_account();
			NewRelD = Rel_discrepancy();
			i--;
		} while (NewRelD > PrevRelD && i != 0);
		PrevRelD = NewRelD;
	}
	printf("Number of iterations: %d", iter);
}
void FEM::Main_Part()
{
	// Создание локальных матриц массы, жестковсти
	vector< vector<double>> M_coef = { {2, 1}, {1, 2} };
	vector< vector<double>> G_coef = { {1, -1}, {-1, 1} };
	vector< vector<double>> local_M = { {0, 0}, {0, 0} };
	vector< vector<double>> local_G = { {0, 0}, {0, 0} };
	vector< vector<double>> local_A = { {0, 0}, {0, 0} };

	// Создание локальных векторов правой части
	vector<double> function = { 0, 0 };
	vector<double> local_b = { 0, 0 };
	//vector<double> right_part = { 0, 0 };

	vector<double> lambda_uh;
	lambda_uh.resize(total_n);
	for (int i = 0; i < total_n; i++)
	{
		lambda_uh[i] = lambda(test, p_i_prev[i]);
	}

	for (int finit_elem = 0; finit_elem < total_n - 1; finit_elem++)
	{
		double x1 = x_grid[finit_elem];
		double x2 = x_grid[finit_elem + 1];
		double hk = x_grid[finit_elem + 1] - x_grid[finit_elem];

		// вычисление матриц
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				local_M[i][j] = gamma(test, 0) * hk / 6 * M_coef[i][j];

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				local_G[i][j] = (lambda_uh[finit_elem] + lambda_uh[finit_elem + 1]) / (2 * hk) * G_coef[i][j];

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				local_A[i][j] = local_M[i][j] + local_G[i][j];

		// вычисление правой части
		function[0] = f(test, x_grid[finit_elem], 0);
		function[1] = f(test, x_grid[finit_elem + 1], 0);
		local_b[0] = hk / 6 * (2 * function[0] + function[1]);
		local_b[1] = hk / 6 * (function[0] + 2 * function[1]);
		

		// заполнение векторов для решения СЛАУ
		di[finit_elem] += local_A[0][0];
		di[finit_elem + 1] += local_A[1][1];
		al[finit_elem] += local_A[1][0];
		au[finit_elem] += local_A[0][1];
		b[finit_elem] += local_b[0];
		b[finit_elem + 1] += local_b[1];

		// находим численные производные лямбды
		double lambda_der1 = NumDer(test, p_i_prev[finit_elem]);
		double lambda_der2 = NumDer(test, p_i_prev[finit_elem + 1]);
		// производные правой части
		double f1 = f(test, x_grid[finit_elem], 0);
		double f2 = f(test, x_grid[finit_elem + 1], 0);
		// производные компонент локальной матрицы
		double A11q1, A12q1, A21q1, A22q1, A11q2, A12q2, A21q2, A22q2;
		double hk2 = 2 * hk;
		A11q1 = lambda_der1 / hk2;
		A12q1 = -lambda_der1 / hk2;
		A21q1 = -lambda_der1 / hk2;
		A22q1 = lambda_der1 / hk2;
		A11q2 = lambda_der2 / hk2;
		A12q2 = -lambda_der2 / hk2;
		A21q2 = -lambda_der2 / hk2;
		A22q2 = lambda_der2 / hk2;

		double q1 = p_i_prev[finit_elem];
		double q2 = p_i_prev[finit_elem + 1];
		// Считаем компоненты локальной матрицы
		double A11, A12, A21, A22, b1, b2;
		A11 = q1 * A11q1 + q2 * A12q1;
		A12 = q1 * A11q2 + q2 * A12q2;
		A21 = q1 * A21q1 + q2 * A22q1;
		A22 = q1 * A21q2 + q2 * A22q2;
		b1 = q1 * (q1 * A11q1 + q2 * A11q2) + q2 * (q1 * A12q1 + q2 * A12q2);
		b2 = q1 * (q1 * A21q1 + q2 * A21q2) + q2 * (q1 * A22q1 + q2 * A22q2);

		// вставляем их в вектора
		diDeriv[finit_elem] += A11;
		diDeriv[finit_elem + 1] += A22;
		alDeriv[finit_elem] += A21;
		auDeriv[finit_elem] += A12;
		bDeriv[finit_elem] += b1;
		bDeriv[finit_elem + 1] += b2;
	}

	// склеивание матриц и векторов
	for (int i = 0; i < total_n; i++)
	{
		diLine[i] = di[i] + diDeriv[i];
		bLine[i] = b[i] + bDeriv[i];
	}
	for (int i = 0; i < total_n - 1; i++)
	{
		alLine[i] = al[i] + alDeriv[i];
		auLine[i] = au[i] + auDeriv[i];
	}
}
void FEM::Condition_account()
{
	priority_queue<BoundCond> bounds{ boundcond };

	for (int cond_number = 0; cond_number < 2; cond_number++)
	{
		BoundCond cond = bounds.top();
		bounds.pop();
		if (cond.type == 1) // условие
		{
			switch (cond.orientation)
			{
			case 'l':
				diLine[0] = 1;
				auLine[0] = 0;
				bLine[0] = Uf(test, x_grid[0], 0);
				di[0] = 1;
				au[0] = 0;
				b[0] = Uf(test, x_grid[0], 0);
				break;
			case 'r':
				diLine[total_n - 1] = 1;
				alLine[total_n - 2] = 0;
				bLine[total_n - 1] = Uf(test, x_grid[total_n - 1], 0);
				di[total_n - 1] = 1;
				al[total_n - 2] = 0;
				b[total_n - 1] = Uf(test, x_grid[total_n - 1], 0);
				break;
			default:
				printf_s("Incorrect bound");
				break;
			}
		}
		else
			printf_s("Invalid type");
	}
}
void FEM::Clearing()
{
	for (int i = 0; i < total_n; i++)
	{
		di[i] = 0;
		diLU[i] = 0;
		diDeriv[i] = 0;
		diLine[i] = 0;

		b[i] = 0;
		bDeriv[i] = 0;
		bLine[i] = 0;
	}
	for (int i = 0; i < total_n - 1; i++)
	{
		al[i] = 0;
		alLU[i] = 0;
		alDeriv[i] = 0;
		alLine[i] = 0;

		au[i] = 0;
		auLU[i] = 0;
		auDeriv[i] = 0;
		auLine[i] = 0;
	}
}

/*----- LU-разложение ------*/

void FEM::LU_decompose()
{

	for (int i = 0; i < total_n; i++)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		double sumD = 0;
		int j = i - (i1 - i0);
		for (int k = i0; k < i1; j++, k++)
		{
			double sumL = 0, sumU = 0;

			int j0 = ia[j], j1 = ia[j + 1];

			int difI = i1 - i0 - (i - j);
			int difJ = j1 - j0;
			int difIJ = difI - difJ;

			int ki = i0, kj = j0;

			if (difIJ < 0)
				kj -= difIJ;
			else
				ki += difIJ;

			for (; ki < k; ki++, kj++)
			{
				sumL += alLU[ki] * auLU[kj];
				sumU += alLU[kj] * auLU[ki];
			}

			alLU[k] = alLine[k] - sumL;
			auLU[k] = (auLine[k] - sumU) / diLU[j];
			sumD += alLU[k] * auLU[k];
		}
		diLU[i] = diLine[i] - sumD;
	}
}
void FEM::Copy_vector(vector<double>& v1, vector<double>& v2)
{
	for (int i = 0; i < v1.size(); i++)
	{
		v2[i] = v1[i];
	}
}
double FEM::Vector_norm(vector<double>v)
{
	double norm = 0;
	for (int i = 0; i < total_n; i++)
	{
		norm += v[i] * v[i];
	}
	return sqrt(norm);
}
double FEM::Rel_discrepancy()
{
	Matrix_mult_vector(p_i_prev, tmp);
	Vector_sub(tmp, b, tmp);
	double Aq = Vector_norm(tmp);
	double b1 = Vector_norm(b);
	return Aq / b1;
}
void FEM::Vector_sub(vector<double> v1, vector<double> v2, vector<double>& res)
{
	for (int i = 0; i < total_n; i++)
	{
		res[i] = v1[i] - v2[i];
	}
}
void FEM::Matrix_mult_vector(vector<double> v1, vector<double>& res)
{
	for (int i = 0; i < total_n; i++)
	{
		res[i] = di[i] * v1[i];
		int j = i - (ia[i + 1] - ia[i]);
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			res[i] += al[k] * v1[j];
			res[j] += au[k] * v1[i];
			j++;
		}
	}
}
void FEM::CalculateY() {
	for (int i = 0; i < total_n; i++)
	{
		double sum = 0;
		int i0 = ia[i], i1 = ia[i + 1];
		int j = i - (i1 - i0);
		for (int k = i0; k < i1; k++, j++)
			sum += y[j] * alLU[k];
		y[i] = (bLine[i] - sum) / diLU[i];
	}
}
void FEM::CalculateX() {
	for (int i = total_n - 1; i >= 0; i--)
	{
		int i0 = ia[i], i1 = ia[i + 1];
		int j = i - (i1 - i0);
		double xi = y[i];
		for (int ij = i0; ij < i1; ij++, j++)
			y[j] -= auLU[ij] * xi;
		p_i[i] = xi;
	}
}

/*----- Норма по Л2 ------*/

double FEM::Lebeg2()
{
	double res = Gauss3();
	return sqrt(res);
}
double FEM::Gauss3()
{
	double sum = 0;
	double p = sqrt(3. / 5.);
	vector<double> q_weight = { 5. / 9., 8. / 9., 5. / 9. };
	for (int i = 0; i < total_n - 1; i++)
	{
		double C = (x_grid[i + 1] - x_grid[i]) / 2;
		double d = (x_grid[i + 1] + x_grid[i]) / 2;

		double x1 = C * p + d;
		double x2 = -C * p + d;
		double x3 = d;

		double y1 = Lebeg_function(x1);
		double y2 = Lebeg_function(x2);
		double y3 = Lebeg_function(x3);

		sum += (q_weight[0] * (y1 + y2) + q_weight[1] * y3) * C;
	}

	return sum;
}
double FEM::Lebeg_function(double x)
{
	//double curTime = t[timelayer];
	return (Find_inner_points(x) - Uf(test, x, 0)) * (Find_inner_points(x) - Uf(test, x, 0));
}
double FEM::Find_inner_points(double point)
{
	int index = 0;
	if (point < x_grid[0] || point > x_grid[total_n - 1])
		throw new exception;

	for (int flag = 0; index < total_n && flag != 1; )
	{
		if (x_grid[index] > point)
			flag = 1;
		else
			index++;
	}
	int i_prev = index - 1;
	int i = index;
	double h = x_grid[i] - x_grid[i_prev];
	double xi_1 = (x_grid[i] - point) / h;
	double xi_2 = (point - x_grid[i_prev]) / h;
	double q1 = p_i[i_prev];
	double q2 = p_i[i];
	return q1 * xi_1 + q2 * xi_2;
}

/*----- Вывод решения ------*/

void FEM::Results()
{
	vector<double> Rel_discrepancy;
	//Rel_discrepancy.resize(time_layers);
	
	double tempU = Uf(test, x_grid[0], 0);
	//printf("%.12lf\t\n", x_grid[0]);
	printf("\n%.12lf\t\n", abs(tempU - p_i[0]));
	for (int i = 1; i < total_n; i++)
	{
		double tempU = Uf(test, x_grid[i], 0);
		double diff = abs(tempU - p_i[i]);
		//printf("%.12lf\t\n", x_grid[i]);
		printf("%.12lf\t\n", diff);
		/*printf("%.14lf\n", abs(uf - q1));
		for (int i = 1; i < total_n; i++)
		{
			double tmp_q = q_layers[layer][i];
			double tmp_Uf = Uf(test, x_grid[i], t[layer]);
			double differ = abs(tmp_Uf - tmp_q);
			Rel_discrepancy[layer - 1] += differ * differ;
			//printf("\t|\t|   %.14lf\t|   %.14lf\t|   %.14lf\t|   %.14lf\t\n", x_grid[i], tmp_q, tmp_Uf, differ);
			//printf("%.14lf\n", x_grid[i]);
			printf("%.14lf\n", differ);
		}
		Rel_discrepancy[layer - 1] = sqrtl(Rel_discrepancy[layer - 1]);*/
	}

	printf("Lebeg2: %.15lf \n", Lebeg2());

}