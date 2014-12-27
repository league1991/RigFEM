// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

using namespace RigFEM;

void testLineSearch()
{
	SimpleFunction sf;
	LineSearch ls(&sf);

	EigVec x(1), dx(1);
	x[0] = 0;
	dx[0] = 1;
	double a;
	int res = ls.lineSearch(x, dx, a);
	PRINT_F("%lf\n");
}

int main(int argc, char* argv[])
{
	testLineSearch();
	return 0;
}

