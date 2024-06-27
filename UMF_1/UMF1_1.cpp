#include "MKR.h"

int main()
{
    MKR solver;
    FILE* params, * points, * boundaries;

    fopen_s(&params, "Params.txt", "r");
    fopen_s(&points, "Points.txt", "r");
    fopen_s(&boundaries, "Bounds.txt", "r");

    solver.Input(params, points, boundaries);

    fclose(params);
    fclose(points);
    fclose(boundaries);

    solver.Solve();
    //solver.PrintTrue();
    //solver.PrintDes();
}