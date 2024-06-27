﻿#include <iostream>
#include "stdafx.h"
#include "FEM.h"

int main()
{
    FEM fem, fem2;
    vector<double> u0, u1;
    int factor = 1;
    // На сетке N
    //GenerateGrid(factor);
    fem.Input();
    fem.Solve();
    u0 = fem.errVec;

    double sum1 = 0;
    /*for (int i = 0; i < u0.size(); i++)
        sum1 += u0[i];*/
    for (int p = 0; p < fem.vertices.size(); p++)
    {
        sum1 += fem.CalculateError(p);
    }
    printf_s("N: %.15e\n", sum1 / fem.timeStamps.size());

    // На сетке 2N
   /* GenerateGrid(factor + 1);
    fem2.Input();
    fem2.Solve();
    u1 = fem2.errVec;
    double sum2 = 0;*/
    /*for (int i = 0; i < u1.size(); i++)
        sum2 += u1[i];*/
    /*for (int p = 0; p < fem2.vertices.size(); p++)
    {
        sum2 += fem2.CalculateError(p);
    }
    printf_s("2N: %.15e\n", sum2 / fem2.timeStamps.size());
    printf_s("Log2: %lf", log2((sum1 / fem.timeStamps.size()) / (sum2 / fem2.timeStamps.size())));*/
}