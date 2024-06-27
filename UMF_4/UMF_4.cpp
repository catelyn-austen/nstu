#include <iostream>
#include "SLAE.h"
#include <ctime>

using namespace std;

int main()
{
	SLAE A;
	if (A.Init())
	{
		cout << "Stopped." << endl;
		return 1;
	}
	time_t t1 = 0, t2 = 0;
	t1 = clock();
	A.BCG();
	t2 = clock();
	A.PrintX();
	cout << "Ended in " << t2 - t1 << " ms" << endl;
}
