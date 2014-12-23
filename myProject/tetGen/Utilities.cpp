#include "StdAfx.h"
#include "Utilities.h"

Utilities::Utilities(void)
{
}

Utilities::~Utilities(void)
{
}

void Utilities::vegaSparse2Eigen( const SparseMatrix& src, EigSparse& tar , int nCols)
{
	int nRows = src.GetNumRows();
	if (nCols == -1)
		nCols = src.GetNumColumns();
	typedef Eigen::Triplet<double> Tri;
	vector<Tri> triList;
	for (int r = 0; r < nRows; ++r)
	{
		for (int ithEntry = 0; ithEntry < src.GetRowLength(r); ++ithEntry)
		{
			int c = src.GetColumnIndex(r, ithEntry);
			double v = src.GetEntry(r, ithEntry);
			triList.push_back(Tri(r,c,v));
		}
	}
	tar = EigSparse(nRows, nCols);
	tar.setFromTriplets(triList.begin(), triList.end());
}

void Utilities::mergeVec( const EigVec& v1, const EigVec& v2, EigVec& v )
{
	int l1 = v1.size();
	int l2 = v2.size();
	v.resize(l1 + l2);
	for (int i = 0; i < l1; ++i)
	{
		v[i] = v1[i];
	}
	for (int i = 0; i < l2; ++i)
	{
		v[i+l1] = v2[i];
	}
}

void MathUtilities::testMath()
{
	double a0 = 2.291, fa0 = 2.027, dfa0 = 1.934;
	double a1 = 0.657, fa1 = 1.062, dfa1 = 1.638;

	double xMin, fxMin;
	bool res = findCubicMin(a0, fa0, dfa0, 
							a1, fa1, dfa1, 
							xMin, fxMin);

	printf("xmin = %lf  fxmin = %lf\n", xMin, fxMin);
}
