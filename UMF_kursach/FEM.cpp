#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vector"
#include "FEM.h"
#include "slae.h"

/* -------- Коэф-ты и функции -------- */
double FEM::coef_lambda(int p, int area)
{
	return 1.0;
}
double FEM::coef_sigma(int p, int area)
{
	return 1.0;
}
double FEM::f(int p, int area, int ind_t)
{
	double x = points[p].x;
	double y = points[p].y;
	double t = time_points[ind_t];

	return 1;
}
double FEM::coef_beta(int p, int number)
{
	return NAN;
}
double FEM::coef_ubeta(int p, int number)
{
	return NAN;
}
double FEM::coef_tetta(int p, int number)
{
	switch (number) {
	case 1:
		return -1;
	case 2:
		return 1;
	}
}
double FEM::Ug(int p, int ind_t)
{
	double x = points[p].x;
	double y = points[p].y;
	double t = time_points[ind_t];

	return x + y + t;
}
double FEM::Ug(double x, double y, double t)
{
	return x + y + t;
}
double FEM::Uq(double x, double y, Triangle triangles, vector<double> result)
{
	double res = 0.0;

	double x0 = points[triangles.v1].x;
	double x1 = points[triangles.v2].x;
	double x2 = points[triangles.v3].x;

	double y0 = points[triangles.v1].y;
	double y1 = points[triangles.v2].y;
	double y2 = points[triangles.v3].y;

	auto S = [&](double x, double y, Vertex v1, Vertex v2) -> double {
		return abs((v2.x - v1.x) * (y - v1.y) - (x - v1.x) * (v2.y - v1.y));
		};

	double detD = abs(Det(triangles));
	double L1 = (S(x, y, points[triangles.v1], points[triangles.v3])) / detD;
	double L2 = (S(x, y, points[triangles.v3], points[triangles.v1])) / detD;
	double L3 = (S(x, y, points[triangles.v1], points[triangles.v2])) / detD;

	res += result[triangles.v1] * L1;
	res += result[triangles.v2] * L2;
	res += result[triangles.v3] * L3;

	return res;
}

/* -------- Ввод -------- */

void FEM::Input()
{
	ifstream file1("Vertices.txt");
	double st, en, step;
	if (file1.is_open()) 
	{
		file1 >> st >> step >> en;
	}
	else
		std::cout << "Error: Unable to open file" << std::endl;

	int intervals = (en - st) / step;
	time_points.reserve(intervals + 1);
	for (int i = 0; i < intervals; i++)
		time_points.push_back(st + step * i);
	time_points.push_back(en);

	// Узлы
	file1 >> points_count;
	Vertex temporary_verts{};
	for (int i = 0; i < points_count; i++)
	{
		file1 >> temporary_verts.x >> temporary_verts.y;
		points.push_back(temporary_verts);
	}
	file1.close();

	// Треугольники
	ifstream file2 ("Triangles.txt");
	file2 >> areas_count;
	Triangle temporary_triangles{};
	for (int i = 0; i < areas_count; i++)
	{
		file2 >> temporary_triangles.v1 >> temporary_triangles.v2 >> temporary_triangles.v3 >> temporary_triangles.area;
		temporary_triangles.v1--;
		temporary_triangles.v2--;
		temporary_triangles.v3--;
		finits.push_back(temporary_triangles);
	}
	file2.close();

	// Краевые
	ifstream file3("BoundaryConditions.txt");
	int num;
	file3 >> num;
	cond_1 first_bc{};
	cond_2 second_bc{};
	cond_3 third_bc{};
	for (int i = 0; i < num; i++)
	{
		int type;
		file3 >> type;
		switch (type)
		{
		case 1:
			file3 >> first_bc.v >> first_bc.eq_n;
			first_bc.v--;
			border_cond1.push_back(first_bc);
			break;
		case 2:
			file3 >> second_bc.v1 >> second_bc.v2 >> second_bc.eq_n;
			second_bc.v1--;
			second_bc.v2--;
			border_cond2.push_back(second_bc);
			break;
		case 3:
			file3 >> third_bc.v1 >> third_bc.v2 >> third_bc.beta_eq_n >> third_bc.ubeta_eq_n;
			third_bc.v1--;
			third_bc.v2--;
			border_cond3.push_back(third_bc);
			break;
		default:
			break;
		}
	}
	file3.close();
}

/* -------- Решение -------- */
void FEM::MainPass()
{
	SLAE slae;
	q_2.resize(points_count);
	// заполняем начальный вектор q(j-2)
	for (int i = 0; i < points_count; i++)
		q_2[i] = Ug(i, 0);
	Portrait();
	allocate_memory();
	q_1.resize(points_count);
	// заполняем вектор q(j-1)
	for (int i = 0; i < points_count; i++)
		q_1[i] = Ug(i, 1);
	// вектор векторов решений
	qs.push_back(q_2);
	qs.push_back(q_1);

	// считаем q на каждом временном слое
	for (int j = 2; j < time_points.size(); j++)
	{
		// согласно формулам вычисляем шаги по времени
		double deltaT, deltaT0, deltaT1;
		deltaT = time_points[j] - time_points[j - 2];
		deltaT1 = time_points[j - 1] - time_points[j - 2];
		deltaT0 = time_points[j] - time_points[j - 1];

		// очищаем и заполняем матрицу глабальную для нового слоя
		clearing();
		Global_matrix(j);

		// (T + T0)/ TT0 * M - первое слагаемое в А
		double coef = (deltaT + deltaT0) / (deltaT * deltaT0);
		for (int i = 0; i < ig[points_count]; i++) {
			gglA[i] = coef * gglM[i] + gglG[i];
			gguA[i] = coef * gguM[i] + gguG[i];
		}
		for (int i = 0; i < points_count; i++)
			diA[i] = coef * diM[i] + diG[i];

		// правая часть
		Global_F(deltaT, deltaT0, deltaT1);
		
		// краевые
		cond_123(j);

		// размерность, макс_итер, епс, вектора. считаем через МСГ
		slae.Input(points_count, 10000, 1e-32, ig.data(), jg.data(), gglA.data(), gguA.data(), diA.data(), b.data());
		slae.MethodOfConjugateGradientsForSymMatrix();

		//q.insert(q.end(), &slae.x[0], &slae.x[slae.n]);
		for (int i = points_count - 1; i > -1; i--)
			q.push_back(slae.x[i]);
		// заносим получившийся q в вектор векторов
		qs.push_back(q);

		cout << "T = " << time_points[j] << " error: " << scientific << Int_norm_error(q, j) << endl;
		//printf_s("%.15e\n", timeStamps[j]);
		for (double el : q)
			printf_s("%.14e\n", el);
		printf_s("\n");
		for (int i = 0; i < points.size(); i++)
			printf_s("%.14e\n", Ug(i, j));
		printf_s("\n\n");
		//errVec.push_back(CalcNorm(q, j));

		q_2 = q_1;
		q_1 = q;
		q.clear();
	}
}
void FEM::Portrait()
{
	ig.resize(points_count + 1);
	int* list[2]{};
	list[0] = new int[2 * points_count * (points_count - 2)];
	list[1] = new int[2 * points_count * (points_count - 2)];
	int* listbeg = new int[points_count];

	int listsize = -1;
	for (int i = 0; i < points_count; i++)
		listbeg[i] = -1;
	for (int ielem = 0; ielem < areas_count; ielem++)
	{
		for (int i = 0; i < 3; i++)
		{
			int k = global_index(finits[ielem], i);
			for (int j = i + 1; j < 3; j++)
			{
				int ind1 = k;
				int ind2 = global_index(finits[ielem], j);
				if (ind2 < ind1) {
					ind1 = ind2;
					ind2 = k;
				}
				int iaddr = listbeg[ind2];
				if (iaddr == -1) {
					listsize++;
					listbeg[ind2] = listsize;
					list[0][listsize] = ind1;
					list[1][listsize] = -1;
				}
				else {
					while (list[0][iaddr] < ind1 && list[1][iaddr] > 0)
						iaddr = list[1][iaddr];

					if (list[0][iaddr] > ind1) {
						listsize++;
						list[0][listsize] = list[0][iaddr];
						list[1][listsize] = list[1][iaddr];
						list[0][iaddr] = ind1;
						list[1][iaddr] = listsize;
					}
					else {
						if (list[0][iaddr] < ind1) {
							listsize++;
							list[1][iaddr] = listsize;
							list[0][listsize] = ind1;
							list[1][listsize] = -1;
						}
					}
				}
			}
		}
	}

	jg.resize(listsize + 1);
	ig[0] = 0;
	for (int i = 0; i < points_count; i++)
	{
		ig[i + 1] = ig[i];
		int iaddr = listbeg[i];
		while (iaddr != -1) {
			jg[ig[i + 1]] = list[0][iaddr];
			ig[i + 1]++;
			iaddr = list[1][iaddr];
		}
	}
	delete[] list[0];
	delete[] list[1];
	delete[] listbeg;
}
int FEM::global_index(Triangle triangles, int i)
{
	switch (i)
	{
	case 0:
		return triangles.v1;
	case 1:
		return triangles.v2;
	case 2:
		return triangles.v3;
	default:
		return NAN;
	}
}
void FEM::allocate_memory()
{
	diA.resize(points_count);
	diM.resize(points_count);
	diG.resize(points_count);
	b.resize(points_count);
	F.resize(points_count);
	gglA.resize(ig[points_count] - ig[0]);
	gguA.resize(ig[points_count] - ig[0]);
	gglM.resize(ig[points_count] - ig[0]);
	gguM.resize(ig[points_count] - ig[0]);
	gglG.resize(ig[points_count] - ig[0]);
	gguG.resize(ig[points_count] - ig[0]);
}
void FEM::clearing()
{
	fill(diA.begin(), diA.end(), 0);
	fill(diM.begin(), diM.end(), 0);
	fill(diG.begin(), diG.end(), 0);
	fill(F.begin(), F.end(), 0);
	fill(b.begin(), b.end(), 0);
	fill(gglA.begin(), gglA.end(), 0);
	fill(gguA.begin(), gguA.end(), 0);
	fill(gglM.begin(), gglM.end(), 0);
	fill(gguM.begin(), gguM.end(), 0);
	fill(gglG.begin(), gglG.end(), 0);
	fill(gguG.begin(), gguG.end(), 0);
}
void FEM::Global_matrix(int ind_t)
{
	// заполняем матрицы на текущем временном слое и разных треугольниках
	for (int r = 0; r < areas_count; r++)
	{
		Triangle temporary_tri = finits[r];
		// формируем матрицы G M и вектор b
		G_matrix(temporary_tri);
		M_matrix(temporary_tri);
		F_local(temporary_tri, ind_t);
		int globalBasis[3] = { temporary_tri.v1, temporary_tri.v2, temporary_tri.v3 };
		// заполняем их
		for (int i = 0; i < 3; i++)
		{
			F[globalBasis[i]] += local_F[i];
			for (int j = 0; j < 3; j++) {
				add_to_G(globalBasis[i], globalBasis[j], G[i][j]);
				add_to_M(globalBasis[i], globalBasis[j], M[i][j]);
			}
		}
	}
}
void FEM::G_matrix(Triangle triangles)
{
	double coef = 1.0 / (fabs(Det(triangles)) * 6.0);
	//double coef = fabs(DetD(tri)) / 6.0;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++) {
			G[i][j] = 0;
			int a[3] = { triangles.v1, triangles.v2, triangles.v3 };
			for (int k = 0; k < 3; k++)
			{
				// согласно формуле: сумм(lambda[k] * (alfa1i * alfa1j + alfa2i * alfa2j) * detD / 6)
				G[i][j] += coef_lambda(a[k], triangles.area) * coef * (coef_alpha(triangles, 1, i) * coef_alpha(triangles, 1, j) + coef_alpha(triangles, 2, i) * coef_alpha(triangles, 2, j));
			}
		}
	}
}
void FEM::M_matrix(Triangle triangles)
{
	std::copy(&coef_M[0][0], &coef_M[0][0] + 9, &M[0][0]);
	// согласно формуле: sigma / 24 * (m)
	double factor = (aver_coef_sigma(triangles) * fabs(Det(triangles))) / 24.0;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			M[i][j] *= factor;
	}
}
void FEM::F_local(Triangle triangles, int ind_t)
{
	double coef = fabs(Det(triangles)) / 24.0;
	// согласно формулам: detD * (f1/12 + f2/24 + f3/24) и тд
	local_F[0] = coef * (2.0 * f(triangles.v1, triangles.area, ind_t) + f(triangles.v2, triangles.area, ind_t) + f(triangles.v3, triangles.area, ind_t));
	local_F[1] = coef * (f(triangles.v1, triangles.area, ind_t) + 2.0 * f(triangles.v2, triangles.area, ind_t) + f(triangles.v3, triangles.area, ind_t));
	local_F[2] = coef * (f(triangles.v1, triangles.area, ind_t) + f(triangles.v2, triangles.area, ind_t) + 2.0 * f(triangles.v3, triangles.area, ind_t));
}
void FEM::add_to_G(int i, int j, double to_add)
{
	if (i == j)
		diG[i] += to_add;
	else if (i < j) {
		int ind;
		for (ind = ig[j]; ind < ig[j + 1]; ind++)
		{
			if (jg[ind] == i)
				break;
		}
		gguG[ind] += to_add;
	}
	else {
		int ind;
		for (ind = ig[i]; ind < ig[i + 1]; ind++)
		{
			if (jg[ind] == j)
				break;
		}
		gglG[ind] += to_add;
	}
}
void FEM::add_to_M(int i, int j, double to_add)
{
	if (i == j)
		diM[i] += to_add;
	else if (i < j) {
		int ind;
		for (ind = ig[j]; ind < ig[j + 1]; ind++)
		{
			if (jg[ind] == i)
				break;
		}
		gguM[ind] += to_add;
	}
	else {
		int ind;
		for (ind = ig[i]; ind < ig[i + 1]; ind++)
		{
			if (jg[ind] == j)
				break;
		}
		gglM[ind] += to_add;
	}
}
void FEM::Global_F(double T, double T0, double T1)
{
	vector<double> temp1, temp2;
	temp1.resize(points_count);
	temp2.resize(points_count);
	double coef1 = T0 / (T1 * T);
	double coef2 = T / (T1 * T0);
	
	// coef1 * M * q(j-2)
	for (int i = 0; i < points_count; i++)
	{
		temp1[i] = coef1 * diM[i] * q_2[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			temp1[i] += coef1 * gglM[k] * q_2[j];
			temp1[j] += coef1 * gguM[k] * q_2[i];
		}
	}

	// coef2 * M * q(j-1)
	for (int i = 0; i < points_count; i++)
	{
		temp2[i] = coef2 * diM[i] * q_1[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			temp2[i] += coef2 * gglM[k] * q_1[j];
			temp2[j] += coef2 * gguM[k] * q_1[i];
		}
	}

	// весь вектор правой части
	for (int i = 0; i < points_count; i++)
	{
		b[i] = F[i] - temp1[i] + temp2[i];
	}
}
void FEM::cond_123(int ind_t)
{
	// сначала 3 краевые
	const double localA[2][2] = { {2.0, 1.0}, {1.0, 2.0} };
	for (int i = 0; i < border_cond3.size(); i++)
	{
		cond_3 temp = border_cond3[i];
		// усредненная beta * hm / 6 * (m)
		double coef1 = (((coef_beta(temp.v1, temp.beta_eq_n) + coef_beta(temp.v2, temp.beta_eq_n)) / 2.0) * Hm(temp.v1, temp.v2)) / 6.0;
		b[temp.v1] += coef1 * (2 * coef_ubeta(temp.v1, temp.ubeta_eq_n) + coef_ubeta(temp.v2, temp.ubeta_eq_n));
		b[temp.v2] += coef1 * (coef_ubeta(temp.v1, temp.ubeta_eq_n) + 2 * coef_ubeta(temp.v2, temp.ubeta_eq_n));

		// меняем теперь А
		int globalBasis[2] = { temp.v1, temp.v2 };
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				add_to_A(globalBasis[i], globalBasis[j], localA[i][j] * coef1);
		}
	}

	// затем 2 краевые
	for (int i = 0; i < border_cond2.size(); i++)
	{
		cond_2 temp = border_cond2[i];
		// hm / 6 * (m)
		double coef = Hm(temp.v1, temp.v2) / 6.0;
		b[temp.v1] += coef * (2 * coef_tetta(temp.v1, temp.eq_n) + coef_tetta(temp.v2, temp.eq_n));
		b[temp.v2] += coef * (coef_tetta(temp.v1, temp.eq_n) + 2 * coef_tetta(temp.v2, temp.eq_n));
	}

	// последними идут 1 краевые
	for (int i = 0; i < border_cond1.size(); i++)
	{
		cond_1 temp = border_cond1[i];
		// b = Ug
		b[temp.v] = Ug(temp.v, ind_t) * 1e+20;
		diA[temp.v] += 1e+20;
	}
}
void FEM::add_to_A(int i, int j, double to_add)
{
	if (i == j)
		diA[i] += to_add;
	else if (i < j) {
		int ind;
		for (ind = ig[j]; ind < ig[j + 1]; ind++)
		{
			if (jg[ind] == i)
				break;
		}
		gguA[ind] += to_add;
	}
	else {
		int ind;
		for (ind = ig[i]; ind < ig[i + 1]; ind++)
		{
			if (jg[ind] == j)
				break;
		}
		gglA[ind] += to_add;
	}
}
double FEM::aver_coef_sigma(Triangle triangles)
{
	return (coef_sigma(triangles.v1, triangles.area) + coef_sigma(triangles.v2, triangles.area) + coef_sigma(triangles.v3, triangles.area)) / 3.0;
}
double FEM::Det(Triangle triangles)
{
	return (points[triangles.v2].x - points[triangles.v1].x) * (points[triangles.v3].y - points[triangles.v1].y) - (points[triangles.v3].x - points[triangles.v1].x) * (points[triangles.v2].y - points[triangles.v1].y);
}
double FEM::coef_alpha(Triangle triangle, int k, int i)
{
	// D^-1
	if (k == 1) {
		switch (i)
		{
		case 0:
			return points[triangle.v2].y - points[triangle.v3].y;
		case 1:
			return points[triangle.v3].y - points[triangle.v1].y;
		case 2:
			return points[triangle.v1].y - points[triangle.v2].y;
		default:
			return NAN;
		}
	}
	else if (k == 2) {
		switch (i)
		{
		case 0:
			return points[triangle.v3].x - points[triangle.v2].x;
		case 1:
			return points[triangle.v1].x - points[triangle.v3].x;
		case 2:
			return points[triangle.v2].x - points[triangle.v1].x;
		default:
			return NAN;
		}
	}
	else
		return NAN;
}
double FEM::Hm(int v1, int v2)
{
	double x = points[v1].x - points[v2].x;
	double y = points[v1].y - points[v2].y;
	return sqrt(x * x + y * y);
}
double FEM::Int_norm_error(vector<double> result, int ind_t)
{
	double res = 0;
	for (int r = 0; r < finits.size(); r++)
	{
		Triangle curTri = finits[r];
		double x0 = points[curTri.v1].x;
		double x1 = points[curTri.v2].x;
		double x2 = points[curTri.v3].x;

		double y0 = points[curTri.v1].y;
		double y1 = points[curTri.v2].y;
		double y2 = points[curTri.v3].y;

		double x = (x0 + x1 + x2) / 3.0;
		double y = (y0 + y1 + y2) / 3.0;
		double J = abs(Det(curTri));
		double f = std::pow(Uq(x, y, curTri, result) - Ug(x, y, time_points[ind_t]), 2.0);
		res += 0.5 * J * f;
	}
	return sqrt(res);
}

/* -------- Вывод погрешности в виде интеграла от ||u - u*|| -------- */
double FEM::Int_norm_error_time(int p)
{
	double res = 0;
	// считаем норму на каждом временном слое
	for (int i = 2; i < time_points.size() - 1; i++)
	{
		double t0 = time_points[i - 1];
		double t1 = time_points[i];

		double t = (t1 - t0) / 2.0 + t0;
		double J = t1 - t0;
		// подинтегральная функция в середине отрезка времени
		double f = pow(Q_resh(t, p, i) - Ug(p, t), 2.0);
		res += 0.5 * J * f;
	}
	return sqrt(res);
}
double FEM::Q_resh(double t, int p, int ind_t)
{
	double deltaT = time_points[ind_t] - time_points[ind_t - 2];
	double deltaT1 = time_points[ind_t - 1] - time_points[ind_t - 2];
	double deltaT0 = time_points[ind_t] - time_points[ind_t - 1];

	auto n2 = [&](double t) -> double {
		return (1.0 / (deltaT * deltaT1)) * (t - time_points[ind_t - 1]) * (t - time_points[ind_t]);
		};
	auto n1 = [&](double t) -> double {
		return (-1.0 / (deltaT1 * deltaT0)) * (t - time_points[ind_t - 2]) * (t - time_points[ind_t]);
		};
	auto n0 = [&](double t) -> double {
		return (1.0 / (deltaT * deltaT0)) * (t - time_points[ind_t - 2]) * (t - time_points[ind_t - 1]);
		};
	return qs[ind_t - 2][p] * n2(t) + qs[ind_t - 1][p] * n1(t) + qs[ind_t][p] * n0(t);
}
