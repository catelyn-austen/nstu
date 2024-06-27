#include <iostream>
#include "FEM.h"
int main()
{
    FEM object;
    FILE* Params, * Points, * Boundary, * Time;

    fopen_s(&Params, "Params.txt", "r");
    fopen_s(&Points, "Points.txt", "r");
    fopen_s(&Boundary, "Bounds.txt", "r");
    fopen_s(&Time, "Time.txt", "r");

    object.Input(Params, Points, Boundary, Time);
    object.Solve();
    object.Results();
}
