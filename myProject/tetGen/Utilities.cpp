#include "StdAfx.h"
#include "Utilities.h"
using namespace std;
using namespace RigFEM;

std::string RigFEM::Utilities::s_whiteSpace = " \t\n\r\f";

void Utilities::vegaSparse2Eigen( const SparseMatrix& src, EigSparse& tar , int nCols)
{
	int nRows = src.GetNumRows();
	if (nCols == -1)
		nCols = src.GetNumColumns();
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

bool RigFEM::Utilities::stringToDense( const std::string& str, int begPos, int& endPos, EigDense& mat, std::string& matName )
{
	// get position of each part
	int trueBeg     = str.find_first_not_of(s_whiteSpace, begPos);
	if (trueBeg     == string::npos)
		return false;
	int equalPos    = str.find('=', trueBeg);
	if (equalPos    == string::npos)
		return false;
	int lBracketPos = str.find('[', equalPos+1);
	if (lBracketPos == string::npos)
		return false;
	int rBracketPos = str.find(']', lBracketPos+1);
	if (rBracketPos == string::npos)
		return false;

	// extract content
	string nameStr = str.substr(trueBeg, equalPos-trueBeg);
	nameStr = nameStr.substr(0, nameStr.find_first_of(s_whiteSpace, 0));
	string numberStr = str.substr(lBracketPos+1, rBracketPos-lBracketPos-1);

	// parse lines
	int lineBeg = 0;
	int nRows, nCols;
	typedef vector<double> MatrixRow;
	vector<MatrixRow> matrixBuf;
	bool hasComma = numberStr.find(',') != string::npos;
	for (int lineEnd = 0;lineBeg < numberStr.size();lineBeg = lineEnd+1)
	{
		// extract line
		lineEnd = numberStr.find_first_of("\n;", lineBeg);
		string lineStr = numberStr.substr(lineBeg, lineEnd-lineBeg);
		if (lineStr == "")
			continue;

		// parse a line
		stringstream lineStream(lineStr);
		MatrixRow matRow;
		while(true)
		{
			double v;
			lineStream >> v;
			if (!lineStream)
				break;
			matRow.push_back(v);

			// eat up commas
			if (hasComma)
			{
				char comma;
				lineStream >> comma;
			}
		}

		if (matRow.size() != 0) 
		{
			if (matrixBuf.size() == 0)
			{
				matrixBuf.push_back(matRow);
				nCols = matRow.size();
			}
			else if(nCols == matRow.size())
			{
				matrixBuf.push_back(matRow);
			}
			else
				return false;
		}

		if (lineEnd == string::npos)
			break;
	}

	// fill matrix
	nRows = matrixBuf.size();
	mat.resize(nRows, nCols);
	for (int r = 0; r < nRows; ++r)
	{
		for (int c = 0; c < nCols; ++c)
		{
			mat(r,c) = matrixBuf[r][c];
		}
	}

	// find end positon
	int semicolonPos= str.find_first_not_of(s_whiteSpace, rBracketPos+1);
	if (semicolonPos != string::npos && str[semicolonPos] == ';')
	{
		endPos = str.find_first_not_of(s_whiteSpace, semicolonPos+1);
	}
	else
	{
		endPos = semicolonPos;
	}
	matName = nameStr;
	return true;
}

int RigFEM::Utilities::skipMatlabComment( const string& str, int begPos )
{
	while (begPos < str.length())
	{
		begPos = str.find_first_not_of(s_whiteSpace, begPos);
		if (begPos == string::npos)
			return str.length();
		if (str[begPos] != '%')
			return begPos;
		begPos = str.find_first_of('\n', begPos);
	}
	return str.length();
}

bool RigFEM::Utilities::fileToDense( const char* fileName, map<string, EigDense>& denseMap )
{
	ifstream file(fileName);
	if (!file)
		return false;
	std::string str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
	file.close();

	int curPos = 0;
	denseMap.clear();
	while (curPos < str.length())
	{
		curPos = skipMatlabComment(str, curPos);
		if (curPos >= str.length())
			break;
		int newPos = curPos;
		EigDense mat;
		string   matName;
		if(!stringToDense(str, curPos, newPos, mat, matName))
			return false;

		denseMap[matName] = mat;
		curPos = newPos;
	}
	return true;
}

bool RigFEM::Utilities::eigDense2Sparse( const EigDense& denseMat, EigSparse& sparseMat )
{
	vector<Tri> triList;

	for (int i = 0; i < denseMat.rows(); ++i)
	{
		for (int j = 0; j < denseMat.cols(); ++j)
		{
			double v = denseMat(i,j);
			if (v == 0.0)
				continue;

			triList.push_back(Tri(i,j,v));
		}
	}
	sparseMat = EigSparse(denseMat.rows(), denseMat.cols());
	sparseMat.setFromTriplets(triList.begin(), triList.end());
	return true;
}

bool RigFEM::Utilities::kronecker3X( EigSparse& src, EigSparse& dst )
{
	int nRow = src.rows();
	int nCol = src.cols();

	vector<Tri> triList;
	for (int c = 0; c < nCol; ++c)
	{
		for (EigSparse::InnerIterator it(src, c); it; ++it)  
		{
			double v = it.value();
			int rId = it.row() * 3;   // row index
			int cId = it.col() * 3;   // col index (here it is equal to k)
			triList.push_back(Tri(rId, cId, v));	rId++; cId++;
			triList.push_back(Tri(rId, cId, v));	rId++; cId++;
			triList.push_back(Tri(rId, cId, v));
		}
	}

	dst = EigSparse(nRow*3, nCol*3);
	dst.setFromTriplets(triList.begin(), triList.end());
	return true;
}

void RigFEM::Utilities::double2EigenDiag( const double* m, int n, EigSparse& diag )
{
	vector<Tri> triList;
	for (int r = 0; r < n; ++r)
	{
		triList.push_back(Tri(r,r,m[r]));
	}
	diag = EigSparse(n,n);
	diag.setFromTriplets(triList.begin(), triList.end());
}


void RigFEM::MathUtilities::clampVector( EigVec& v, double maxElement )
{
	double maxAbsVal = v.cwiseAbs().maxCoeff();
	if (maxAbsVal > maxElement)
		v *= (maxElement / maxAbsVal);
}
