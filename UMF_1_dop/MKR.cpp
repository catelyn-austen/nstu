#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include "MKR.h"

using namespace std;

/*-------- Функции --------*/
double MKR::Uf(double x, double y, int i)
{
    switch (i)
    {
        case 1: return x + y;
        default: return x + y;
        //case 1: return y * y + x * x;
        //default: return x * x + y * y;
        //case 1: return pow(y, 3) + pow(x, 3);
        // //default: return pow(y, 3) + pow(x, 3);
        //case 1: return cos(2 * x + 2 * y);
        //default: return cos(2 * x + 2 * y);
        
        
    }
}
double MKR::lambda(int i)
{
    switch (i)
    {
    case 1: return 1;
    default: return 1;
    }
}
double MKR::gamma(int i)
{
    switch (i)
    {
    case 1: return 1;
    default: return 1;
    }
}
double MKR::f(double x, double y, int i)
{
    switch (i)
    {
        case 1: return x + y;
        default: return x + y;
        //case 1: return y * y + x * x - 4;
        //default: return x * x + y * y - 4;
        //case 1: return pow(y, 3) + pow(x, 3) - 6 * x * y;
        //default: return pow(y, 3) + pow(x, 3) - 6 * x * y;
        //case 1: return 9 * cos(2 * x + 2 * y);
        //default: return 9 * cos(2 * x + 2 * y);
    
    }
}
double MKR::Ug(double x, double y, int i)
{
    switch (i)
    {
        case 1: return y + x;
        default: return x + y;
        //case 1: return y * y + x * x;
        //default: return x * x + y * y;
        //case 1: return pow(y, 3) + pow(x, 3);
        //default: return pow(y, 3) + pow(x, 3);
        //case 1: return cos(2 * x + 2 * y);
        //default: return cos(2 * x + 2 * y);
    }
}
double MKR::Thetta(double x, double y, int i)
{
    switch (i)
    {
        case 1: return -1;
        case 2: return 1;
        case 3: return -1;
        case 4: return 1;
        case 5: return 1;
        case 6: return -1;
        case 7: return 1;
        case 8: return -1;
        default: return 1;

        /*case 1: return -1 * 2 * y;
        case 2: return 2 * x;
        case 3: return - 1 * 2 * y;
        case 4: return 2 * x;
        case 5: return 2 * y;
        case 6: return -1 * 2 * x;
        case 7: return 2 * y;
        case 8: return - 1 * 2 * x;*/

        /*case 1: return -1 * 3 * y * y;
        case 2: return 3 * x * x;
        case 3: return -1 * 3 * y * y;
        case 4: return 3 * x * x;
        case 5: return 3 * y * y;
        case 6: return -1 * 3 * x * x;
        case 7: return 3 * y * y;
        case 8: return -1 * 3 * x * x;
        default: return 1;*/

        /*case 1: return 2 * sin(2 * x + 2 * y);
        case 2: return -2 * sin(2 * x + 2 * y);
        case 3: return 2 * sin(2 * x + 2 * y);
        case 4: return -2 * sin(2 * x + 2 * y);
        case 5: return -2 * sin(2 * x + 2 * y);
        case 6: return 2 * sin(2 * x + 2 * y);
        case 7: return -2 * sin(2 * x + 2 * y);
        case 8: return -2 * sin(2 * x + 2 * y);
        default: return 1;*/
    }
}

/*-------- Чтение из файлов --------*/
void MKR::input_grid()
{
    ifstream in("grid.txt");
    if (!in) exit(-1);

    // Считываем Х
    in >> x_points_n;
    x_points.resize(x_points_n);
    for (int i = 0; i < x_points_n; i++)
    {
        in >> x_points[i];
    }

    // Считываем У
    in >> y_points_n;
    y_points.resize(y_points_n);
    for (int i = 0; i < y_points_n; i++)
    {
        in >> y_points[i];
    }

    // Подобласти (прям-ники): сюда заносятся именно индексы х-в и у-в, соответствующие массивам
    in >> areas_n;
    sub_areas.resize(areas_n);
    for (int i = 0; i < areas_n; i++)
    {
        in >> sub_areas[i].number
            >> sub_areas[i].x_left
            >> sub_areas[i].x_right
            >> sub_areas[i].y_down
            >> sub_areas[i].y_up;
    }
    in.close();
}
void MKR::input_part()
{
    ifstream in("partitions.txt");
    if (!in) exit(-1);

    x_i.resize(x_points_n);
    y_i.resize(y_points_n);
    int cur = 0;
    // Подобласти отрезков по х
    for (int i = 0; i < x_points_n - 1; i++)
    {
        int k; // Количество подынтервалов на каждом отрезке
        double q; // Коэффициент q
        double tmp = 1, fin_coef = 1, step;
        in >> k >> q;
        // Находим количество отрезков на глобальном отрезке
        for (int j = 1; j < k; j++)
        {
            tmp *= q;
            fin_coef += tmp;
        }
        step = (x_points[i + 1] - x_points[i]) / fin_coef;
        x_part.resize(cur + k);
        x_part[cur] = x_points[i];
        x_i[i] = cur;
        tmp = 1;
        // заполняем массив с с разбитым на отрезки глобальным отрезком
        for (int j = cur + 1; j < cur + k; j++)
        {
            x_part[j] = x_part[j - 1] + step * tmp;
            tmp *= q;
        }
        cur += k;
    }
    x_part.resize(cur + 1);
    x_part[cur] = x_points[x_points_n - 1];
    x_i[x_points_n - 1] = cur;
    x_part_n = cur + 1;

    //Построение разбиения по y
    cur = 0;
    for (int i = 0; i < y_points_n - 1; i++)
    {
        int k; // Количество подынтервалов на каждом отрезке
        double q; // Коэффициент
        double tmp = 1, fin_coef = 1, step;
        in >> k >> q;
        for (int j = 1; j < k; j++)
        {
            tmp *= q;
            fin_coef += tmp;
        }
        step = (y_points[i + 1] - y_points[i]) / fin_coef;
        y_part.resize(cur + k);
        y_part[cur] = y_points[i];
        y_i[i] = cur;
        tmp = 1;
        for (int j = cur + 1; j < cur + k; j++)
        {
            y_part[j] = y_part[j - 1] + step * tmp;
            tmp *= q;
        }
        cur += k;
    }
    y_part.resize(cur + 1);
    y_part[cur] = y_points[y_points_n - 1];
    y_i[y_points_n - 1] = cur;
    y_part_n = cur + 1;
}
void MKR::input_bounds()
{
    ifstream in("kraevie.txt");
    in >> bound_n;
    bounds.resize(bound_n);
    for (int i = 0; i < bound_n; i++)
    {
        bounds[i].resize(6);
        // Точки, краевое (1 или 2), номер границы
        int point, bound, border;
        in >> bound >> border;
        
        bounds[i][0] = bound;
        bounds[i][1] = border;
        for (int j = 2; j < 6; j++)
        {
            in >> point;
            bounds[i][j] = point - 1;
        }
    }
}

/*-------- Формирование матрицы --------*/
void MKR::AllocateMemory()
{
    N = x_part_n * y_part_n;
    ig = { -x_part_n, -1, 1, x_part_n }; // Так как неглавные диаг-ли сдвинуты таким вот образом в матрице в нормальном её виде
    di1.resize(N - x_part_n);
    di2.resize(N - 1);
    di3.resize(N);
    di4.resize(N - 1);
    di5.resize(N - x_part_n);
    b.resize(N);
}
void MKR::MainPass()
{
    // нумерация с 1
    for (int y = 1; y < y_part_n - 1; y++)
    {
        for (int x = 1; x < x_part_n - 1; x++)
        {
            int numArea = GetAreaNum(y, x);
            // Если область найдена
            if (numArea != -1)
            {
                vector <int> gn = {
                    x_part_n * (y - 1) + x, // номер нижней по у точки
                    x_part_n * y + x - 1, // номер левой по х точки
                    x_part_n * y + x, // номер тек. точки
                    x_part_n * y + x + 1, // номер правой по х точки
                    x_part_n * (y + 1) + x // номер верхней по у точки
                };

                // Шаги по x
                double hx_l = x_part[x] - x_part[x - 1];
                double hx_r = x_part[x + 1] - x_part[x];

                // Шаги по y
                double hy_l = y_part[y] - y_part[y - 1];
                double hy_r = y_part[y + 1] - y_part[y];

                // Находим коэф-ы
                double lam = lambda(numArea);
                double gam = gamma(numArea);

                // Находим коэф-ы перед uшками
                double kl_x = 2.0 * lam / (hx_l * (hx_l + hx_r)); // 2 диагональ
                double kx = 2.0 * lam / (hx_l * hx_r); // главная диагональ 
                double kr_x = 2.0 * lam / (hx_r * (hx_l + hx_r)); // 4 диагональ

                double kl_y = 2.0 * lam / (hy_l * (hy_l + hy_r)); // 1 диагональ
                double ky = 2.0 * lam / (hy_l * hy_r); // главная диагональ
                double kr_y = 2.0 * lam / (hy_r * (hy_l + hy_r)); // 5 диагональ

                // Добавляем в элемент текцщей точки производные
                AddDerivative(gn[2], gn[0], -kl_y);
                AddDerivative(gn[2], gn[1], -kl_x);
                AddDerivative(gn[2], gn[2], kx + ky + gam);
                AddDerivative(gn[2], gn[3], -kr_x);
                AddDerivative(gn[2], gn[4], -kr_y);

                b[gn[2]] += f(x_part[x], y_part[y], numArea);
            }
        }
    }
}
int MKR::GetAreaNum(int y, int x)
{
    for (int i = 0; i < areas_n; i++)
    {
        // Вычисляем координаты прямоугольника
        int x0 = x_i[sub_areas[i].x_left - 1], x1 = x_i[sub_areas[i].x_right - 1];
        int y0 = y_i[sub_areas[i].y_down - 1], y1 = y_i[sub_areas[i].y_up - 1];
        // Проверяем, содержится ли точка в данной области
        if ((x >= x0 && x <= x1) && (y >= y0 && y <= y1))
            if ((x + 1 >= x0 && x + 1 <= x1) && (y >= y0 && y <= y1))
                if ((x >= x0 && x <= x1) && (y + 1 >= y0 && y + 1 <= y1))
                    if ((x + 1 >= x0 && x + 1 <= x1) && (y + 1 >= y0 && y + 1 <= y1))
                        return sub_areas[i].number;
    }
    return -1;
}
void MKR::AddDerivative(int i, int j, double elem)
{
    if (i == j)
    {
        di3[i] += elem;
    }
    else
    {
        int tmp = j - i;
        if (tmp == ig[0])
            di1[j] += elem;
        if (tmp == ig[1])
            di2[j] += elem;
        if (tmp == ig[2])
            di4[i] += elem;
        if (tmp == ig[3])
            di5[i] += elem;
    }
}
void MKR::BoundaruPass()
{
    //sort(matrixKraev.begin(), matrixKraev.end());
    // сначала 2 краевое, потом 1
    for (int i = 0; i < bound_n; i++)
    {
        if (bounds[i][0] == 2)
        {
            kraevoe_2(bounds[i]);
        }
    }
    for (int i = 0; i < bound_n; i++)
    {
        if (bounds[i][0] == 1)
        {
            kraevoe_1(bounds[i]);
        }
    }
}

void MKR::kraevoe_2(vector<int>& kr_2)
{
    int numForm = kr_2[1];
    // Если границы вертикальная, т. е. x0 = x1
    if (kr_2[2] == kr_2[3])
    {
        // Проходим по y снизу вверх
        for (int k = y_i[kr_2[4]]; k < y_i[kr_2[5]]; k++)
        {
            double coef = 0;
            // Если х не на самой крайней границе
            if ((kr_2[2] + 1) < x_i.size() && x_i[kr_2[2] + 1] < x_part.size())
            {
                /*if ((k + 1) != y_i[kr_2[5]])
                {
                    coef = lambda(numForm) / (x_part[x_i[kr_2[2]]] - x_part[x_i[kr_2[2] - 1]]);
                    vector<int> globalNum = {
                       x_part_n * k + x_i[kr_2[2]] + x_part_n,
                       x_part_n * k + x_i[kr_2[2]] + x_part_n - 1
                    };
                    //ClearString(globalNum[0], 0);
                    AddDerivative(globalNum[0], globalNum[0], coef);
                    AddDerivative(globalNum[0], globalNum[1], -coef);
                    b[globalNum[0]] = Thetta(x_part[x_i[kr_2[3]]], y_part[k], numForm);
                }*/
                coef = lambda(numForm) / (x_part[x_i[kr_2[2] + 1]] - x_part[x_i[kr_2[2]]]);
                vector<int> gn = {
                 x_part_n * k + x_i[kr_2[2]],
                 x_part_n * k + x_i[kr_2[2]] + x_part_n
                };
                ClearString(gn[0], 0);
                AddDerivative(gn[0], gn[0], coef);
                AddDerivative(gn[0], gn[1], -coef);
                b[gn[0]] = Thetta(x_part[x_i[kr_2[2]]], y_part[k], numForm);
            }
            else {
                /*if ((k + 1) != y_i[kr_2[5]])
                {
                    coef = lambda(numForm) / (x_part[x_i[kr_2[2]]] - x_part[x_i[kr_2[2] - 1]]);
                    vector<int> globalNum = {
                       x_part_n * k + x_i[kr_2[2]] + x_part_n,
                       x_part_n * k + x_i[kr_2[2]] + x_part_n - 1
                    };
                    //ClearString(globalNum[0], 0);
                    AddDerivative(globalNum[0], globalNum[0], coef);
                    AddDerivative(globalNum[0], globalNum[1], -coef);
                    b[globalNum[0]] = Thetta(x_part[x_i[kr_2[3]]], y_part[k], numForm);
                }*/
                coef = lambda(numForm) / (x_part[x_i[kr_2[2]]] - x_part[x_i[kr_2[2] - 1]]);
                vector<int> gn = {
                   x_part_n * k + x_i[kr_2[2]] + x_part_n,
                   x_part_n * k + x_i[kr_2[2]] + x_part_n - 1
                };
                ClearString(gn[0], 0);
                AddDerivative(gn[0], gn[0], coef);
                AddDerivative(gn[0], gn[1], -coef);
                b[gn[0]] = Thetta(x_part[x_i[kr_2[3]]], y_part[k], numForm);
            }
        }
    }
    // Горизонталь
    else
    {
        for (int k = x_i[kr_2[2]]; k < x_i[kr_2[3]]; k++)
        {
            double coef = 0;
            if (kr_2[4] + 1 < y_i.size() && y_i[kr_2[4] + 1] < y_part.size()) {

                coef = lambda(numForm) / (y_part[y_i[kr_2[4] + 1]] - y_part[y_i[kr_2[4]]]);
                vector<int> gn = {
                 x_part_n * y_i[kr_2[4]] + k + 1,
                 x_part_n * (y_i[kr_2[4]]) + x_part_n
                };
                ClearString(gn[0], 0);
                AddDerivative(gn[0], gn[0], coef);
                AddDerivative(gn[0], gn[1], -coef);

                b[gn[0]] = Thetta(x_part[x_i[kr_2[2]]], y_part[k], numForm);
            }
            else {

                coef = lambda(numForm) / (y_part[y_i[kr_2[4]]] - y_part[y_i[kr_2[4]] - 1]);
                vector<int> gn = {
                   x_part_n * y_i[kr_2[4]] + k + 1,
                   x_part_n * (y_i[kr_2[4]] - 1) + x_part_n
                };
                ClearString(gn[0], 0);
                AddDerivative(gn[0], gn[0], coef);
                AddDerivative(gn[0], gn[1], -coef);
                b[gn[0]] = Thetta(x_part[x_i[kr_2[3]]], y_part[k], numForm);
            }

        }
    }
}

void MKR::kraevoe_1(vector<int>& kr_1)
{
    int numForm = kr_1[1];
    // Если граница вертикальная
    if (kr_1[2] == kr_1[3])
    {
        for (int k = y_i[kr_1[4]]; k <= y_i[kr_1[5]]; k++)
        {
            int x = x_i[kr_1[2]];
            int gn = k * x_part_n + x;
            ClearString(gn, 1);
            b[gn] = Ug(x_part[x], y_part[k], numForm);
        }
    }
    // Если граница горизонтальная
    else
    {
        for (int k = x_i[kr_1[2]]; k <= x_i[kr_1[3]]; k++)
        {
            int y = y_i[kr_1[4]];
            int gn = x_part_n * y + k;
            ClearString(gn, 1);
            b[gn] = Ug(x_part[k], y_part[y], numForm);
        }
    }
}

void MKR::ClearString(int i, int diValue)
{
    if (i + ig[0] >= 0)
        di1[i + ig[0]] = 0;
    if (i + ig[1] >= 0)
        di2[i + ig[1]] = 0;
    di3[i] = diValue;
    if (i + ig[2] < N)
        di4[i] = 0;
    if (i + ig[3] < N)
        di5[i] = 0;
    b[i] = 0;
}

void MKR::FictionNodes()
{
    for (int y = 0; y < y_part_n; y++)
    {
        for (int x = 0; x < x_part_n; x++)
        {
            int gn = y * x_part_n + x;
            if (IsFict(y, x)) {
                di3[gn] = 1;
                b[gn] = 0;
            }
        }
    }
}
bool MKR::IsFict(int y, int x)
{
    // Проверяем вхождение точки в каждую область
    for (int i = 0; i < areas_n; i++)
    {
        int x0 = x_i[sub_areas[i].x_left - 1];
        int x1 = x_i[sub_areas[i].x_right - 1];
        int y0 = y_i[sub_areas[i].y_down - 1];
        int y1 = y_i[sub_areas[i].y_up - 1];
        if (((x >= x0) && (x <= x1) && (y >= y0) && (y <= y1)))
            return false;
    }
    return true;
}

/*-------- СЛАУ Г_З --------*/
void MKR::Solve()
{
    u.resize(N);
    double normb = NormVector(b);
    double relDiscrepancy = sqrt(RelDiscrepancy(u) / normb);
    uk = u;

    for (int k = 0; k < max_iter && relDiscrepancy > eps; k++)
    {
        relDiscrepancy = 0;
        for (int i = 0; i < N; i++)
        {
            double sum = di3[i] * u[i];
            if (i + ig[0] >= 0) sum += di1[i + ig[0]] * u[i + ig[0]];
            if (i + ig[1] >= 0) sum += di2[i + ig[1]] * u[i + ig[1]];
            if (i + ig[2] < N) sum += di4[i] * u[i + ig[2]];
            if (i + ig[3] < N) sum += di5[i] * u[i + ig[3]];
            u[i] = u[i] + w * (b[i] - sum) / di3[i];

            relDiscrepancy += (b[i] - sum) * (b[i] - sum);
        }
        relDiscrepancy = sqrt(relDiscrepancy / normb);
        swap(u, u);
    }

    /*FILE* f = fopen("output.txt", "w");
    if (f != NULL)
    {
        for (int i = 0; i < N; i++)
            fprintf(f, "%.15lf\n", u[i]);
        fclose(f);
    }*/
}
double MKR::NormVector(vector <double>& a)
{
    double norm = 0;

    for (int i = 0; i < a.size(); i++)
        norm += a[i] * a[i];
    return norm;
}
double MKR::RelDiscrepancy(vector <double>& u)
{
    double discrepancy = 0;
    for (int i = 0; i < N; i++)
    {
        double summ = di3[i] * u[i];
        if (i + ig[0] >= 0) 
            summ += di1[i + ig[0]] * u[i + ig[0]];
        if (i + ig[1] >= 0) 
            summ += di2[i + ig[1]] * u[i + ig[1]];
        if (i + ig[2] < N) 
            summ += di4[i] * u[i + ig[2]];
        if (i + ig[3] < N) 
            summ += di5[i] * u[i + ig[3]];
        discrepancy += (b[i] - summ) * (b[i] - summ);
    }
    return discrepancy;
}
void MKR::Output()
{
    vector <double> uR(N);
    for (int s = 0; s < y_part_n; s++)
    {
        for (int p = 0; p < x_part_n; p++)
        {
            int globalNum = s * x_part_n + p;
            int numArea = GetAreaNum(s, p);
            if (IsFict(s, p))
            {
                uR[globalNum] = 0;
                u[globalNum] = 0;
            }
            else
            {
                uR[globalNum] = Uf(x_part[p], y_part[s], numArea);
            }
            //printf_s("% .15lf % .15lf\n", x_part[p], y_part[s]);
            //printf_s("% .15lf\n",u[globalNum]);
            //printf_s("% .15lf\n",uR[globalNum]);
            printf_s("% .15lf % .15lf % .15lf\n", x_part[p], y_part[s], abs(u[globalNum] - uR[globalNum]));
        }
    }
}