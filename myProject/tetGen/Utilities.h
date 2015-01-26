#pragma once


namespace RigFEM
{

class Utilities
{
public:
	static void vegaSparse2Eigen(const SparseMatrix& src, EigSparse& tar, int nCols = -1);
	
	// 把vega矩阵的某个子矩阵转成Eigen矩阵，rowID，colID指定选中的行和列
	static void vegaSparse2Eigen( const SparseMatrix& src, const std::vector<int>& rowID, const std::vector<int>& colID, EigSparse& tar);

	template<class MatrixType>
	static double maxError(const MatrixType& m0, const MatrixType& m1)
	{
		MatrixType dm = m0 - m1;
		double error = 0;
		for (int r = 0; r < dm.rows(); ++r)
		{
			for (int c = 0; c < dm.cols(); ++c)
			{
				double absVal = abs(dm.coeff(r,c));
				error = error > absVal ? error : absVal;
			}
		}
		return error;
	}

	// 把Eigen矩阵写到字符串
	template<class MatrixType>
	static std::string matToString(const MatrixType& mat, const char* name)
	{
		std::string str;
		str += (string(name) + "=[\n");
		char digitStr[50];
		for (int ithRow = 0; ithRow < mat.rows(); ++ithRow)
		{
			for (int ithCol = 0; ithCol < mat.cols(); ++ithCol)
			{
				sprintf(digitStr, "%lf ", mat.coeff(ithRow, ithCol));
				str += digitStr;
			}
			str += ";\n";
		}
		str += "\n];\n";
		return str;
	}

	template<class VectorType>
	static std::string vecToString( const VectorType&vec, const char*name )
	{
		std::string str;
		str += (string(name) + "=[\n");
		char digitStr[50];
		for (int i = 0; i < vec.size(); ++i)
		{
			double val = vec[i];
			sprintf(digitStr, "%lf ", val);
			str += digitStr;
		}
		str += "\n];\n";
		return str;
	}

	template<class VectorType>
	static bool	saveVector(const VectorType& vec, const char* fileName)
	{
		std::ofstream file(fileName);
		if (file)
		{
			int length = vec.size();
			file << length;
			for (int i = 0; i < length; ++i)
			{
				file << ' ' << vec[i];
			}
			file.close();
			return true;
		}
		return false;
	}
	template<class VectorType>
	static bool	loadVector(VectorType& vec, const char* fileName)
	{
		std::ifstream file(fileName);
		if (file)
		{
			int length = 0;
			file >> length;
			if (length <= 0)
				return false;

			vec.resize(length);
			for (int i = 0; i < length; ++i)
				file >> vec[i];
			return true;
		}
		return false;
	}
	// v = [v1 v2]
	static void mergeVec(const EigVec& v1, const EigVec& v2, EigVec& v);

	static void transformBBox(	const double srcMin[3], 
								const double srcMax[3], 
								const double mat[4][4], 
								double dstMin[3], 
								double dstMax[3]);
};

class MathUtilities
{
public:
	static double absSin(double t)
	{
		return 3 * abs(sin(t));
	}
	static double zero(double t)
	{
		return 0;
	}
	static double linear(double t)
	{
		return t;
	}

	// 找出插值二次函数f的最值点x,x未必在区间中
	// a   > 0
	// f0  = f(0)
	// df0 = f'(0) < 0 
	// fa  = f(a)
	static inline bool findQuadricMin(double f0, double df0, double a, double fa, double& xMin, double& fxMin)
	{
		double coefA = (fa - f0 - a * df0) / (a * a);
		if (coefA < 0)
			return false;
		double coefB = df0;
		double coefC = f0;
		double x = -coefB / (2 * coefA);
		xMin = x;
		fxMin= coefA * x * x + coefB * x + coefC;
		return true;
	}

	// 找出插值三次函数f的极小值点
	// a0,a1 > 0
	// f0  = f(0)
	// df0 = f'(0)  < 0
	// fa0 = f(a0)
	// fa1 = f(a1)
	// 若存在这样的点，返回true，否则返回false
	static inline bool findCubicMin0(	double f0, double df0, double a0, double fa0, double a1, double fa1, 
										double& xMin, double& fxMin)
	{
		double a02 = a0 * a0;
		double a12 = a1 * a1;
		double a03 = a02 * a0;
		double a13 = a12 * a1;

		double devisor = a02*a12*(a1 - a0);
		double v0 = fa1 - f0 - df0 * a1;
		double v1 = fa0 - f0 - df0 * a0;

		double coefA = (a02 * v0 - a12 * v1) / devisor;
		double coefB = (-a03* v0 + a13 * v1) / devisor;
		double coefC = df0;
		double coefD = f0;

		double r    = coefB * coefB - 3 * coefA * df0;
		if (r < 0)
			return false;
		xMin  = (-coefB + sqrt(r)) / (3 * coefA);
		fxMin = coefD + xMin * (coefC + xMin * (coefB + xMin * coefA));
		return true;
	}

	// 找出插值三次函数f的极小值点
	// a0,a1 > 0
	// fa0  = f(a0)
	// fa1  = f(a1)
	// dfa0 = f'(a0)
	// dfa1 = f'(a1)
	static inline bool findCubicMin( double a0, double fa0, double dfa0, 
									 double a1, double fa1, double dfa1,
									 double& xMin, double& fxMin)
	{
		// 映射到另一个函数g(x), 此函数有以下特点
		// f(x) = g((x-a1)/(a0 - a1))
		// f(a1)   = g(0)
		// f(a0)   = g(1)
		// f'(a1)  = g'(0) / (a0 - a1)
		// f'(a0)  = g'(1) / (a0 - a1)
		double scale = a0 - a1;
		double dg0   = dfa1 * scale;
		double dg1   = dfa0 * scale;
		double g0    = fa1;
		double g1    = fa0;

		// 计算g(x) = coefA * x^3 + coefB * x^2 + coefC * x + coefD
		double coefA = dg1 - 2*g1 + dg0 + 2 * g0;
		double coefB = 3*g1 - 3*g0 - 2*dg0 - dg1;
		double coefC = dg0;
		double coefD = g0;

		// 计算g(x) 的极小值点
		double delta = coefB * coefB - 3 * coefA * coefC;
		if (delta < 0)
			return false;
		double gxMin = (-coefB + sqrt(delta)) / (3 * coefA);
		double gMin  = coefD + gxMin*(coefC + gxMin*(coefB + gxMin* coefA));

		// 映射回原来的值
		xMin = a1 + gxMin * scale;
		fxMin= gMin;
		return true;
	}


	static void testMath();
};

}