// test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

using namespace RigFEM;

void testLineSearch()
{
// 	SimpleFunction sf;
// 	LineSearch ls(&sf);
// 
// 	EigVec x(1), dx(1);
// 	x[0] = 0;
// 	dx[0] = 1;
// 	double a;
// 	int res = ls.lineSearch(x, dx, a);
// 	PRINT_F("%lf\n");
}

void testMat()
{
	string s = "a = [1, 2, 3;4, 5, 6]; b = [1, \n 2]";
	int endPos;
	EigDense mat;
	string   matName;
	Utilities::stringToDense(s, 0, endPos, mat, matName);
	int begPos = endPos;
	Utilities::stringToDense(s, begPos, endPos, mat, matName);

}

void testMatFile()
{
	map<string, EigDense> matMap;
	Utilities::fileToDense("I:/Programs/VegaFEM-v2.1/myProject/RigPlugin/result.m", matMap);
}
int main(int argc, char* argv[])
{
	//testLineSearch();
	//testMat();
	testMatFile();
	return 0;
}

