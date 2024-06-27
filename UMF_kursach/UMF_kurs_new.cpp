#include <iostream>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "vector"
#include "FEM.h"

int main()
{
    FEM A;
    //vector<double> u_resh;

    A.Input();
    A.MainPass();
    //u_resh = A.errVec;

    double summ = 0;
    for (int p = 0; p < A.points.size(); p++)
    {
        summ += A.Int_norm_error_time(p);
    }
    printf_s("Error: %.15e\n", summ / A.time_points.size());
}