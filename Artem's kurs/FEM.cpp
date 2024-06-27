#include "stdafx.h"
#include "FEM.h"
#include "slae.h"

double FEM::Lamda(int vert, int region)
{
	// ������� �� ���������� ������
	return 1.0;
}

double FEM::Sigma(int vert, int region)
{
	// ������� �� ���������� ������
	return 1.0;
}

double FEM::Function(int vert, int region, int tInd)
{
	double x = vertices[vert].x;
	double y = vertices[vert].y;
	double t = timeStamps[tInd];
	double h = 1e-4;

	return 4 * t * t * t;


	/*double deltaT = timeStamps[tInd] - timeStamps[tInd - 2];
	double deltaT1 = timeStamps[tInd - 1] - timeStamps[tInd - 2];
	double deltaT0 = timeStamps[tInd] - timeStamps[tInd - 1];

	double factor1 = deltaT0 / (deltaT * deltaT1);
	double factor2 = -deltaT / (deltaT1 * deltaT0);
	double factor3 = (deltaT + deltaT0) / (deltaT * deltaT0);*/

	//return factor1 * Ug(x, y, timeStamps[tInd - 2]) + factor2 * Ug(x, y, timeStamps[tInd - 1]) + factor3 * Ug(x, y, t) - Lamda(vert, region) * DivGrad(vert, tInd, h);
}

double FEM::Beta(int vert, int eqNum)
{
	// ������� �� ���������� ������
	return NAN;
}

double FEM::Ubeta(int vert, int eqNum)
{
	// ������� �� ���������� ������
	return NAN;
}

double FEM::Theta(int vert, int eqNum)
{
	// ������� �� ���������� ������
	switch (eqNum) {
	case 1:
		return -1;
	case 2:
		return 1;
	}
}

double FEM::Ug(int vert, int tInd)
{
	return Ug(vertices[vert].x, vertices[vert].y, timeStamps[tInd]);
}

double FEM::Ug(double x, double y, double t)
{
	return x + y + t * t * t * t;
}

double FEM::Uq(double x, double y, Triangle tri, vector<double> resQ)
{
	double res = 0.0;

	double x0 = vertices[tri.vert1].x;
	double x1 = vertices[tri.vert2].x;
	double x2 = vertices[tri.vert3].x;

	double y0 = vertices[tri.vert1].y;
	double y1 = vertices[tri.vert2].y;
	double y2 = vertices[tri.vert3].y;

	auto S = [&](double x, double y, Vertex vert1, Vertex vert2) -> double {
		return abs((vert2.x - vert1.x) * (y - vert1.y) - (x - vert1.x) * (vert2.y - vert1.y));
		};

	double detD = abs(DetD(tri));
	double L1 = (S(x, y, vertices[tri.vert2], vertices[tri.vert3])) / detD;
	double L2 = (S(x, y, vertices[tri.vert3], vertices[tri.vert1])) / detD;
	double L3 = (S(x, y, vertices[tri.vert1], vertices[tri.vert2])) / detD;

	res += resQ[tri.vert1] * L1;
	res += resQ[tri.vert2] * L2;
	res += resQ[tri.vert3] * L3;

	return res;
}



void FEM::Input()
{
	FILE* file;
	fopen_s(&file, "Vertices.txt", "r");
	// ���������� ��������� �����
	double tStart, tEnd, tStep;

	fscanf_s(file, "%lf %lf %lf", &tStart, &tStep, &tEnd);
	int tIntervalNum = (tEnd - tStart) / tStep;
	timeStamps.reserve(tIntervalNum + 1);
	for (int i = 0; i < tIntervalNum; i++)
		timeStamps.push_back(tStart + tStep * i);
	timeStamps.push_back(tEnd);
	// ���������� �����
	int num;
	fscanf_s(file, "%d", &num);
	Vertex tempVert{};

	globalN = num;

	for (int i = 0; i < num; i++)
	{
		fscanf_s(file, "%lf %lf", &tempVert.x, &tempVert.y);
		vertices.push_back(tempVert);
	}
	fclose(file);
	// ���������� �������������
	fopen_s(&file, "Triangles.txt", "r");

	fscanf_s(file, "%d", &num);
	Triangle tempTri{};

	regionsNum = num;

	for (int i = 0; i < num; i++)
	{
		fscanf_s(file, "%d %d %d %d", &tempTri.vert1, &tempTri.vert2, &tempTri.vert3, &tempTri.region);
		tempTri.vert1--;
		tempTri.vert2--;
		tempTri.vert3--;
		tris.push_back(tempTri);
	}
	fclose(file);

	// ���������� �������
	fopen_s(&file, "BoundaryConditions.txt", "r");

	fscanf_s(file, "%d", &num);
	FirstBoundaryCondition firstBoundTemp{};
	SecondBoundaryCondition secondBoundTemp{};
	ThirdBoundaryCondition thirdBoundTemp{};

	for (int i = 0; i < num; i++)
	{
		int type;
		fscanf_s(file, "%d", &type);
		switch (type)
		{
		case 1:
			fscanf_s(file, "%d %d", &firstBoundTemp.vert, &firstBoundTemp.equationNum);
			firstBoundTemp.vert--;
			firstBoundary.push_back(firstBoundTemp);
			break;
		case 2:
			fscanf_s(file, "%d %d %d", &secondBoundTemp.vert1, &secondBoundTemp.vert2, &secondBoundTemp.equationNum);
			secondBoundTemp.vert1--;
			secondBoundTemp.vert2--;
			secondBoundary.push_back(secondBoundTemp);
			break;
		case 3:
			fscanf_s(file, "%d %d %d %d", &thirdBoundTemp.vert1, &thirdBoundTemp.vert2, &thirdBoundTemp.betaEquationNum, &thirdBoundTemp.UbetaEquationNum);
			thirdBoundTemp.vert1--;
			thirdBoundTemp.vert2--;
			thirdBoundary.push_back(thirdBoundTemp);
			break;
		default:
			break;
		}
	}
	fclose(file);
}

void FEM::Solve()
{
	SLAE slae;
	q_3.resize(globalN);
	for (int i = 0; i < globalN; i++)
		q_3[i] = Ug(vertices[i].x, vertices[i].y, timeStamps[0]);

	q_2.resize(globalN);
	for (int i = 0; i < globalN; i++)
		q_2[i] = Ug(vertices[i].x, vertices[i].y, timeStamps[1]);

	FormPortrait();
	AllocateGlobalMatrices();
	q_1.resize(globalN);
	for (int i = 0; i < globalN; i++)
		q_1[i] = Ug(vertices[i].x, vertices[i].y, timeStamps[2]);

	qs.push_back(q_3);
	qs.push_back(q_2);
	qs.push_back(q_1);

	double sss = 0;

	for (int j = 3; j < timeStamps.size(); j++)
	{
		double deltaT0, deltaT1, deltaT2, deltaT3, deltaT4, deltaT5;

		deltaT0 = timeStamps[j] - timeStamps[j - 1];
		deltaT1 = timeStamps[j] - timeStamps[j - 2];
		deltaT2 = timeStamps[j] - timeStamps[j - 3];
		deltaT3 = timeStamps[j - 1] - timeStamps[j - 2];
		deltaT4 = timeStamps[j - 1] - timeStamps[j - 3];
		deltaT5 = timeStamps[j - 2] - timeStamps[j - 3];

		/*deltaT = timeStamps[j] - timeStamps[j - 2];
		deltaT1 = timeStamps[j - 1] - timeStamps[j - 2];
		deltaT0 = timeStamps[j] - timeStamps[j - 1];*/

		ClearMatrices();
		FormGlobalMatrices(j);
		// �������� ���������� A
		//double factor = (deltaT + deltaT0) / (deltaT * deltaT0);
		double factor = (deltaT0 * deltaT1 + deltaT0 * deltaT2 + deltaT1 * deltaT2) / (deltaT0 * deltaT1 * deltaT2);

		for (int i = 0; i < ig[globalN]; i++) {
			gglA[i] = factor * gglM[i] + gglG[i];
			gguA[i] = factor * gguM[i] + gguG[i];
		}
		for (int i = 0; i < globalN; i++)
			diA[i] = factor * diM[i] + diG[i];

		// �������� ������ �����
		FormGlobalB(deltaT0, deltaT1, deltaT2, deltaT3, deltaT4, deltaT5);


		ResolveBoundaries(j);
		slae.Input(globalN, 10000, 1e-32, ig.data(), jg.data(), gglA.data(), gguA.data(), diA.data(), b.data());
		slae.MethodOfConjugateGradientsForSymMatrix();
		//qs.push_back(q);
		//slae.MethodOfConjugateGradientsForSymMatrix();
		q.insert(q.end(), &slae.x[0], &slae.x[slae.n]);
		qs.push_back(q);
		cout << "T = " << timeStamps[j] << " error: " << scientific << CalcNorm(q, j) << endl;
		//printf_s("%.15e\n", timeStamps[j]);
		//printf_s("%.15e\n", q[globalN - 1]);
		//double kaka = pow(abs(Ug(vertices.size() - 1, j) - q[globalN - 1]), 2);
		//sss += kaka;
		//printf_s("%.15e\n", Ug(vertices.size() - 1, j));
		//printf_s("%.15e\n", kaka);

		for (double el : q)
			printf_s("%.14e\n", el);
		printf_s("\n");

		for (int i = 0; i < vertices.size(); i++)
			printf_s("%.14e\n", Ug(i, j));
		printf_s("\n\n");
		errVec.push_back(CalcNorm(q, j));
		q_3 = q_2;
		q_2 = q_1;
		q_1 = q;
		q.clear();
	}
	//printf_s("%.15e\n", sqrt(sss));
	//cout << qs.size();
}

void FEM::PrintSolution()
{
	for (int i = 0; i < globalN; i++)
		printf_s("%.14lf\n", q[i]);
}

double FEM::GetAverageLamda(Triangle tri)
{
	return (Lamda(tri.vert1, tri.region) + Lamda(tri.vert2, tri.region) + Lamda(tri.vert3, tri.region)) / 3.0;
}

double FEM::GetAverageSigma(Triangle tri)
{
	return (Sigma(tri.vert1, tri.region) + Sigma(tri.vert2, tri.region) + Sigma(tri.vert3, tri.region)) / 3.0;
}

double FEM::DetD(Triangle tri)
{
	return (vertices[tri.vert2].x - vertices[tri.vert1].x) * (vertices[tri.vert3].y - vertices[tri.vert1].y) - (vertices[tri.vert3].x - vertices[tri.vert1].x) * (vertices[tri.vert2].y - vertices[tri.vert1].y);
}

double FEM::Alpha(Triangle tri, int k, int i)
{

	if (k == 1) {
		switch (i)
		{
		case 0:
			return vertices[tri.vert2].y - vertices[tri.vert3].y;
		case 1:
			return vertices[tri.vert3].y - vertices[tri.vert1].y;
		case 2:
			return vertices[tri.vert1].y - vertices[tri.vert2].y;
		default:
			return NAN;
		}
	}
	else if (k == 2) {
		switch (i)
		{
		case 0:
			return vertices[tri.vert3].x - vertices[tri.vert2].x;
		case 1:
			return vertices[tri.vert1].x - vertices[tri.vert3].x;
		case 2:
			return vertices[tri.vert2].x - vertices[tri.vert1].x;
		default:
			return NAN;
		}
	}
	else
		return NAN;
}

double FEM::EdgeLength(int vert1, int vert2)
{
	double x = vertices[vert1].x - vertices[vert2].x;
	double y = vertices[vert1].y - vertices[vert2].y;
	return sqrt(x * x + y * y);
}

double FEM::DivGrad(int vert, int tInd, double h)
{
	double x = vertices[vert].x;
	double y = vertices[vert].y;
	double t = timeStamps[tInd];
	double d2udx2 = (Ug(x - h, y, t) - 2 * Ug(x, y, t) + Ug(x + h, y, t)) / (h * h);
	double d2udy2 = (Ug(x, y - h, t) - 2 * Ug(x, y, t) + Ug(x, y + h, t)) / (h * h);
	//double d2udt2 = (Ug(x, y, t - h) - 2 * Ug(x, y, t) + Ug(x, y, t + h)) / (h * h);

	return d2udx2 + d2udy2;
}

double FEM::CalcNorm(vector<double> resQ, int tInd)
{
	double res = 0;
	for (int r = 0; r < tris.size(); r++)
	{
		Triangle curTri = tris[r];
		double x0 = vertices[curTri.vert1].x;
		double x1 = vertices[curTri.vert2].x;
		double x2 = vertices[curTri.vert3].x;

		double y0 = vertices[curTri.vert1].y;
		double y1 = vertices[curTri.vert2].y;
		double y2 = vertices[curTri.vert3].y;

		double x = (x0 + x1 + x2) / 3.0;
		double y = (y0 + y1 + y2) / 3.0;
		double J = abs(DetD(curTri));
		double f = std::pow(Uq(x, y, curTri, resQ) - Ug(x, y, timeStamps[tInd]), 2.0);
		res += 0.5 * J * f;
	}
	return sqrt(res);
}

int FEM::IndexOfUnknown(Triangle tri, int i)
{
	switch (i)
	{
	case 0:
		return tri.vert1;
	case 1:
		return tri.vert2;
	case 2:
		return tri.vert3;
	default:
		return NAN;
	}
}

void FEM::FormLocalM(Triangle tri)
{
	std::copy(&pureM[0][0], &pureM[0][0] + 9, &M[0][0]);
	double factor = (GetAverageSigma(tri) * fabs(DetD(tri))) / 24.0;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			M[i][j] *= factor;
	}
}

void FEM::FormLocalG(Triangle tri)
{
	double factor = 1.0 / (fabs(DetD(tri)) * 6.0);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++) {
			G[i][j] = 0;
			int helper[3] = { tri.vert1, tri.vert2, tri.vert3 };
			for (int k = 0; k < 3; k++)
			{
				G[i][j] += Lamda(helper[k], tri.region) * factor * (Alpha(tri, 1, i) * Alpha(tri, 1, j) + Alpha(tri, 2, i) * Alpha(tri, 2, j));
			}
		}
	}
}

void FEM::FormGlobalMatrices(int tInd)
{
	for (int r = 0; r < regionsNum; r++)
	{
		Triangle currTri = tris[r];
		FormLocalG(currTri);
		FormLocalM(currTri);
		FormLocalB(currTri, tInd);
		int globalBasis[3] = { currTri.vert1, currTri.vert2, currTri.vert3 };
		for (int i = 0; i < 3; i++)
		{
			pureB[globalBasis[i]] += localB[i];
			for (int j = 0; j < 3; j++) {
				AddToGlobalG(globalBasis[i], globalBasis[j], G[i][j]);
				AddToGlobalM(globalBasis[i], globalBasis[j], M[i][j]);
			}
		}
	}
}

void FEM::FormGlobalB(double deltaT0, double deltaT1, double deltaT2, double deltaT3, double deltaT4, double deltaT5)
{
	vector<double> temp1, temp2, temp3;
	temp1.resize(globalN);
	temp2.resize(globalN);
	temp3.resize(globalN);

	double factor1 = (deltaT0 * deltaT1) / (deltaT2 * deltaT4 * deltaT5);
	double factor2 = (deltaT0 * deltaT2) / (deltaT1 * deltaT3 * deltaT5);
	double factor3 = (deltaT1 * deltaT2) / (deltaT0 * deltaT3 * deltaT4);

	// ��� ��������� factor1*M+G �� q(j-1)
	// factor1 * M * q(j-2)
	/*for (int i = 0; i < globalN; i++)
	{
		temp1[i] = factor1 * diM[i] * q_2[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			temp1[i] += factor1 * gglM[k] * q_2[j];
			temp1[j] += factor1 * gguM[k] * q_2[i];
		}
	}*/

	// factor1 * M * q(j-3)
	for (int i = 0; i < globalN; i++)
	{
		temp1[i] = factor1 * diM[i] * q_3[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			temp1[i] += factor1 * gglM[k] * q_3[j];
			temp1[j] += factor1 * gguM[k] * q_3[i];
		}
	}

	// factor2 * M * q(j-1)
	/*for (int i = 0; i < globalN; i++)
	{
		temp2[i] = factor2 * diM[i] * q_1[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			temp2[i] += factor2 * gglM[k] * q_1[j];
			temp2[j] += factor2 * gguM[k] * q_1[i];
		}
	}*/

	// factor2 * M * q(j-2)
	for (int i = 0; i < globalN; i++)
	{
		temp2[i] = factor2 * diM[i] * q_2[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			temp2[i] += factor2 * gglM[k] * q_2[j];
			temp2[j] += factor2 * gguM[k] * q_2[i];
		}
	}

	// factor3 * M * q(j-1)
	for (int i = 0; i < globalN; i++)
	{
		temp3[i] = factor3 * diM[i] * q_1[i];
		for (int k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			temp3[i] += factor3 * gglM[k] * q_1[j];
			temp3[j] += factor3 * gguM[k] * q_1[i];
		}
	}

	// ����������� ������ �����

	for (int i = 0; i < globalN; i++)
	{
		b[i] = pureB[i] + temp1[i] - temp2[i] + temp3[i];
	}
}

void FEM::FormPortrait()
{
	ig.resize(globalN + 1);
	int* list[2]{};
	list[0] = new int[2 * globalN * (globalN - 2)];
	list[1] = new int[2 * globalN * (globalN - 2)];
	int* listbeg = new int[globalN];

	int listsize = -1;
	for (int i = 0; i < globalN; i++)
		listbeg[i] = -1;
	for (int ielem = 0; ielem < regionsNum; ielem++)
	{
		for (int i = 0; i < 3; i++)
		{
			int k = IndexOfUnknown(tris[ielem], i);
			for (int j = i + 1; j < 3; j++)
			{
				int ind1 = k;
				int ind2 = IndexOfUnknown(tris[ielem], j);
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
	for (int i = 0; i < globalN; i++)
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

void FEM::ResolveBoundaries(int tInd)
{
	// ���� 3 �������
	const double localA[2][2] = { {2.0, 1.0}, {1.0, 2.0} };
	for (int i = 0; i < thirdBoundary.size(); i++)
	{
		ThirdBoundaryCondition temp = thirdBoundary[i];
		double factor = (((Beta(temp.vert1, temp.betaEquationNum) + Beta(temp.vert2, temp.betaEquationNum)) / 2.0) * EdgeLength(temp.vert1, temp.vert2)) / 6.0;
		b[temp.vert1] += factor * (2 * Ubeta(temp.vert1, temp.UbetaEquationNum) + Ubeta(temp.vert2, temp.UbetaEquationNum));
		b[temp.vert2] += factor * (Ubeta(temp.vert1, temp.UbetaEquationNum) + 2 * Ubeta(temp.vert2, temp.UbetaEquationNum));

		int globalBasis[2] = { temp.vert1, temp.vert2 };
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				AddToGlobalA(globalBasis[i], globalBasis[j], localA[i][j] * factor);
		}
	}
	// ���� 2 �������
	for (int i = 0; i < secondBoundary.size(); i++)
	{
		SecondBoundaryCondition temp = secondBoundary[i];
		double factor = EdgeLength(temp.vert1, temp.vert2) / 6.0;
		double t1 = factor * (2 * Theta(temp.vert1, temp.equationNum) + Theta(temp.vert2, temp.equationNum));
		b[temp.vert1] += factor * (2 * Theta(temp.vert1, temp.equationNum) + Theta(temp.vert2, temp.equationNum));
		b[temp.vert2] += factor * (Theta(temp.vert1, temp.equationNum) + 2 * Theta(temp.vert2, temp.equationNum));
	}
	// ���� 1 �������
	for (int i = 0; i < firstBoundary.size(); i++)
	{
		FirstBoundaryCondition temp = firstBoundary[i];
		b[temp.vert] = Ug(temp.vert, tInd) * 1e+20;

		diA[temp.vert] += 1e+20;
	}
}

void FEM::AddToGlobalA(int i, int j, double add)
{
	if (i == j)
		diA[i] += add;
	else if (i < j) {
		int ind;
		for (ind = ig[j]; ind < ig[j + 1]; ind++)
		{
			if (jg[ind] == i)
				break;
		}
		gguA[ind] += add;
	}
	else {
		int ind;
		for (ind = ig[i]; ind < ig[i + 1]; ind++)
		{
			if (jg[ind] == j)
				break;
		}
		gglA[ind] += add;
	}
}

void FEM::AddToGlobalG(int i, int j, double add)
{
	if (i == j)
		diG[i] += add;
	else if (i < j) {
		int ind;
		for (ind = ig[j]; ind < ig[j + 1]; ind++)
		{
			if (jg[ind] == i)
				break;
		}
		gguG[ind] += add;
	}
	else {
		int ind;
		for (ind = ig[i]; ind < ig[i + 1]; ind++)
		{
			if (jg[ind] == j)
				break;
		}
		gglG[ind] += add;
	}
}

void FEM::AddToGlobalM(int i, int j, double add)
{
	if (i == j)
		diM[i] += add;
	else if (i < j) {
		int ind;
		for (ind = ig[j]; ind < ig[j + 1]; ind++)
		{
			if (jg[ind] == i)
				break;
		}
		gguM[ind] += add;
	}
	else {
		int ind;
		for (ind = ig[i]; ind < ig[i + 1]; ind++)
		{
			if (jg[ind] == j)
				break;
		}
		gglM[ind] += add;
	}
}

void FEM::AllocateGlobalMatrices()
{
	diA.resize(globalN);
	diM.resize(globalN);
	diG.resize(globalN);
	b.resize(globalN);
	pureB.resize(globalN);

	gglA.resize(ig[globalN] - ig[0]);
	gguA.resize(ig[globalN] - ig[0]);
	gglM.resize(ig[globalN] - ig[0]);
	gguM.resize(ig[globalN] - ig[0]);
	gglG.resize(ig[globalN] - ig[0]);
	gguG.resize(ig[globalN] - ig[0]);
}

void FEM::FormLocalB(Triangle tri, int tInd)
{
	double factor = fabs(DetD(tri)) / 24.0;

	localB[0] = factor * (2.0 * Function(tri.vert1, tri.region, tInd) + Function(tri.vert2, tri.region, tInd) + Function(tri.vert3, tri.region, tInd));
	localB[1] = factor * (Function(tri.vert1, tri.region, tInd) + 2.0 * Function(tri.vert2, tri.region, tInd) + Function(tri.vert3, tri.region, tInd));
	localB[2] = factor * (Function(tri.vert1, tri.region, tInd) + Function(tri.vert2, tri.region, tInd) + 2.0 * Function(tri.vert3, tri.region, tInd));
}

void FEM::ClearMatrices()
{
	fill(diA.begin(), diA.end(), 0);
	fill(diM.begin(), diM.end(), 0);
	fill(diG.begin(), diG.end(), 0);
	fill(pureB.begin(), pureB.end(), 0);
	fill(b.begin(), b.end(), 0);

	fill(gglA.begin(), gglA.end(), 0);
	fill(gguA.begin(), gguA.end(), 0);
	fill(gglM.begin(), gglM.end(), 0);
	fill(gguM.begin(), gguM.end(), 0);
	fill(gglG.begin(), gglG.end(), 0);
	fill(gguG.begin(), gguG.end(), 0);
}

double FEM::CalculateError(int r_ind)
{
	double res = 0;
	for (int i = 3; i < timeStamps.size() - 1; i++)
	{
		double t0 = timeStamps[i - 1];
		double t1 = timeStamps[i];

		double t = (t1 - t0) / 2.0 + t0;
		double J = t1 - t0;
		double f = pow(Uq(t, r_ind, i) - Ug(r_ind, t), 2.0);
		res += J * f;
	}
	return sqrt(res);
}

double FEM::Uq(double t, int r_ind, int tInd)
{
	/*double deltaT = timeStamps[tInd] - timeStamps[tInd - 2];
	double deltaT1 = timeStamps[tInd - 1] - timeStamps[tInd - 2];
	double deltaT0 = timeStamps[tInd] - timeStamps[tInd - 1];*/

	double deltaT0, deltaT1, deltaT2, deltaT3, deltaT4, deltaT5;

	deltaT0 = timeStamps[tInd] - timeStamps[tInd - 1];
	deltaT1 = timeStamps[tInd] - timeStamps[tInd - 2];
	deltaT2 = timeStamps[tInd] - timeStamps[tInd - 3];
	deltaT3 = timeStamps[tInd - 1] - timeStamps[tInd - 2];
	deltaT4 = timeStamps[tInd - 1] - timeStamps[tInd - 3];
	deltaT5 = timeStamps[tInd - 2] - timeStamps[tInd - 3];

	auto n3 = [&](double t) -> double {
		return (-1.0 / (deltaT5 * deltaT4 * deltaT2)) * (t - timeStamps[tInd - 2]) * (t - timeStamps[tInd - 1]) * (t - timeStamps[tInd]);
		};
	auto n2 = [&](double t) -> double {
		return (1.0 / (deltaT5 * deltaT3 * deltaT1)) * (t - timeStamps[tInd - 3]) * (t - timeStamps[tInd - 1]) * (t - timeStamps[tInd]);
		};
	auto n1 = [&](double t) -> double {
		return (-1.0 / (deltaT4 * deltaT3 * deltaT0)) * (t - timeStamps[tInd - 3]) * (t - timeStamps[tInd - 2]) * (t - timeStamps[tInd]);
		};
	auto n0 = [&](double t) -> double {
		return (1.0 / (deltaT0 * deltaT1 * deltaT2)) * (t - timeStamps[tInd - 3]) * (t - timeStamps[tInd - 2]) * (t - timeStamps[tInd - 1]);
		};

	return qs[tInd - 3][r_ind] * n3(t) + qs[tInd - 2][r_ind] * n2(t) + qs[tInd - 1][r_ind] * n1(t) + qs[tInd][r_ind] * n0(t);

	/*auto n2 = [&](double t) -> double {
		return (1.0 / (deltaT * deltaT1)) * (t - timeStamps[tInd - 1]) * (t - timeStamps[tInd]);
		};
	auto n1 = [&](double t) -> double {
		return (-1.0 / (deltaT1 * deltaT0)) * (t - timeStamps[tInd - 2]) * (t - timeStamps[tInd]);
		};
	auto n0 = [&](double t) -> double {
		return (1.0 / (deltaT * deltaT0)) * (t - timeStamps[tInd - 2]) * (t - timeStamps[tInd - 1]);
		};*/
	//return qs[tInd - 2][r_ind] * n2(t) + qs[tInd - 1][r_ind] * n1(t) + qs[tInd][r_ind] * n0(t);
}
