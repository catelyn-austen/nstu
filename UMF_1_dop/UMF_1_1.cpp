#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include "MKR.h"

using namespace std;

#pragma warning(disable:4996)

int main()
{
    MKR A;
    A.input_grid();
    A.input_part();
    A.input_bounds();
    A.AllocateMemory();
    A.MainPass();
    A.BoundaruPass();
    A.FictionNodes();
    A.Solve();
    A.Output();
}
