#include "StdAfx.h"
#include "Utilities.h"
using namespace std;
using namespace RigFEM;

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

void Utilities::vegaSparse2Eigen( const SparseMatrix& src, 
								  const vector<int>& rowID, const vector<int>& colID,
								  EigSparse& tar)
{
	typedef Eigen::Triplet<double> Tri;
	vector<Tri> triList;
	for (int ithTarRow = 0; ithTarRow < rowID.size(); ++ithTarRow)
	{
		// 找到原矩阵对应的行
		int ithSrcRow = rowID[ithTarRow];	
		for (int ithSrcEntry = 0, ithTarCol = 0; ithSrcEntry < src.GetRowLength(ithSrcRow); ++ithSrcEntry)
		{
			int ithSrcCol = src.GetColumnIndex(ithSrcRow, ithSrcEntry);		// 原矩阵元素的列号
			while (ithTarCol < colID.size()-1 && colID[ithTarCol] < ithSrcCol)// 找到大于等于这个列号的最小目标列号
				ithTarCol++;

			if (ithSrcCol == colID[ithTarCol])								// 若当前列在目标列号中，加入目标矩阵
			{
				double v = src.GetEntry(ithSrcRow, ithSrcEntry);
				triList.push_back(Tri(ithTarRow,ithTarCol,v));
			}
		}
	}
	tar = EigSparse(rowID.size(), colID.size());
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

	PRINT_F("xmin = %lf  fxmin = %lf\n", xMin, fxMin);
}



void RigFEM::Utilities::transformBBox( const double srcMin[3], const double srcMax[3], const double mat[4][4], double dstMin[3], double dstMax[3] )
{
	double x[2] = {srcMin[0], srcMax[0]};
	double y[2] = {srcMin[1], srcMax[1]};
	double z[2] = {srcMin[2], srcMax[2]};

	dstMin[0] = dstMin[1] = dstMin[2] = DBL_MAX;
	dstMax[0] = dstMax[1] = dstMax[2] = -DBL_MAX;
	for(int ithPnt = 0; ithPnt < 8; ++ithPnt)
	{
		int xIdx =  ithPnt & 0x1;
		int yIdx = (ithPnt >> 1) & 0x1;
		int zIdx = (ithPnt >> 2) & 0x1;

		for(int j = 0; j < 3; ++j)
		{
			double rj = x[xIdx]*mat[0][j] + y[yIdx]*mat[1][j] + z[zIdx]*mat[2][j] + mat[3][j];
			dstMin[j] = dstMin[j] < rj ? dstMin[j] : rj;
			dstMax[j] = dstMax[j] > rj ? dstMax[j] : rj;
		}
	}
}
