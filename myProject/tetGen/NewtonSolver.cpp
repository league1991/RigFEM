#include "StdAfx.h"
#include "newtonSolver.h"

using namespace RigFEM;

LineSearch::~LineSearch(void)
{
}

int RigFEM::LineSearch::guessMinPnt( double a0, double fa0, double *dfa0, double a1, double fa1, double *dfa1, double& aMin )
{
	double a, fa;
	// 利用 f(a0) f(a1) f'(a0) f'(a1) 拟合出一个三次函数
	if(dfa0 && dfa1 && MathUtilities::findCubicMin(a0, fa0, *dfa0, a1, fa1, *dfa1, a, fa))
	{
		if ((a0 - a) * (a1 - a) < 0)
		{
			aMin = a;			
			return 0;
		}
	}
	// 利用 f(0) f(a0) f(a1) f'(0) 拟合出一个三次函数
	if (a0 != 0.0 && a1 != 0.0 && MathUtilities::findCubicMin0(m_f0, m_df0, a, fa0, a1, fa1, a, fa))
	{
		if ((a0 - a) * (a1 - a) < 0)
		{
			aMin = a;
			return 0;
		}
	}
	// 利用 f(0)  f(a1) f'(0) 拟合出一个二次函数
	if (MathUtilities::findQuadricMin(m_f0, m_df0, a1, fa1, a, fa))
	{
		if ((a0 - a) * (a1 - a) < 0)
		{
			aMin = a;
			return 0;
		}
	}
	// 各种拟合方法都失败,取两个端点
	if (fa0 < fa1)
	{
		aMin = a0;
		return 1;
	}

	aMin = a1;
	return 2;
}

void RigFEM::LineSearch::setOriginAndDir( const EigVec& x0, const EigVec& dx )
{
	m_x0 = x0;
	m_dx = dx;
}

int RigFEM::LineSearch::zoom( double al, double fal, double* dfal, double ah, double fah, double* dfah, double& a )
{
	int maxIter = 20;
	while(maxIter--)
	{
		double aj;
		int res = guessMinPnt(al, fal, dfal, ah, fah, dfah, aj);
		double fj;
		EigVec xj = m_x0 + m_dx * aj;
		computeValueAndDeri(xj, &fj, NULL);

		if (fj > m_f0 + m_c1 * aj * m_df0 || fj >= fal)
		{
			// 插值后得出的新点不满足wolfe充分下降条件，或函数值比fal大
			ah = aj;
		}
		else
		{
			// 满足wolfe充分下降条件，且函数值比fal小
			double dfj;
			computeValueAndDeri(xj, NULL, &dfj);

			// 满足强wolfe条件,终止搜索
			if (abs(dfj) <= m_c2 * abs(m_df0))
			{
				a = aj;
				return 0;
			}
			if (dfj * (ah - al) >= 0)
			{
				ah = al;
			}
			al = aj;
		}

		if (res)
		{
			break;
		}
	}

	// 没有找到符合wolfe条件的新点
	a = al;
	return 1;
}

int RigFEM::LineSearch::lineSearch( const EigVec& x0, const EigVec& dx, const EigVec& param, double& aFinal,
								    double* f0, double *df0)
{
	setOriginAndDir(x0, dx);
	m_param = param;
	double amax = m_maxStep;

	// 必要时计算起始点函数
	double* pF0 = f0 ? NULL : &m_f0;
	double* pDf0= df0 ? NULL : &m_df0;
	computeValueAndDeri(m_x0, pF0, pDf0);
	m_f0 = f0 ? *f0 : m_f0;
	m_df0= df0? *df0: m_df0;

	double ai_1 = 0;
	double fi_1 = m_f0;
	double dfi_1= m_df0;

	double ai   = m_initStep;

	for (int i = 1; ai < amax; ++i)
	{
		PRINT_F("%dth line search iter, a = %lf", i-1,ai);
		EigVec x = m_x0 + m_dx * ai;

		// 只计算函数值，尽量提高效率
		double fi;
		double dfi;
		int t0 = clock();
		computeValueAndDeri(x, &fi, NULL);
		int t1 = clock();
		PRINT_F("line search compute value and deri #1: %f", (t1-t0)/1000.f);

		// 若wolfe充分下降条件被违反或函数值较上次增长,终止迭代
		if (fi > m_f0 + m_c1 * ai * m_df0 ||
			(i > 1 && fi >= fi_1))
		{
			int res = zoom(ai_1, fi_1, &dfi_1, ai, fi, &dfi, aFinal);
			return res;
		}

		computeValueAndDeri(x, NULL, &dfi);
		int t2 = clock();
		PRINT_F("line search compute value and deri #2: %f", (t2-t1)/1000.f);

		// wolfe充分下降条件被满足，函数值较上次下降				
		if (abs(dfi) <= m_c2 * abs(m_df0))
		{	
			// wolfe曲率条件被满足,直接退出
			aFinal = ai;
			return 0;
		}

		if (dfi >= 0)
		{
			// 函数已经开始上升，因此在这之间必能找到符合曲率条件的点
			int res = zoom(ai, fi, &dfi, ai_1, fi_1, &dfi_1, aFinal);
			return res;
		}

		// 更新步长
		ai_1 = ai;
		fi_1 = fi;
		dfi_1= dfi;

		ai *= 1.5;
	}
	return 1;
}

bool RigFEM::LineSearch::computeValueAndDeri( const EigVec& x, double* f, double* df )
{
	if (!f && !df)
		return false;

	EigVec grad;
	EigVec* pG = df ? &grad : NULL;
	m_objFunc->computeValueAndGrad(x, m_param, f, pG);

	if (df)
	{
		*df = grad.dot(m_dx); 
	}
	return true;
}

RigFEM::LineSearch::LineSearch( ObjectFunction* objFun /*= NULL*/, double initStep /*= 1.0*/, double maxStep /*= 10*/ ) :m_objFunc(objFun), m_initStep(initStep), m_maxStep(maxStep), m_c1(1e-4), m_c2(0.9)
{

}

bool RigFEM::NewtonSolver::step()
{
	PRINT_F("########### newton iteration begin. #############");
	EigVec p,n;
	m_fem->getDof(n,p);
	double t = m_fem->getCurTime();
	int nIntDof = n.size();
	int nParam  = p.size();
	EigVec tVec(1);
	tVec[0] = t;

	bool isSucceed = true;
	EigDense  Hpp, Hnp, Hpn;
	EigSparse Hnn;
	EigVec	 G,G0,x,dx,dN,dP,dGN,dGP;
	int maxIter = m_maxIter;
	for (int ithIter = 0; ithIter < maxIter; ++ithIter)
	{
		// 计算当前q p值的函数值、梯度值和Hessian
		PRINT_F("%dth Iteration", ithIter);
		int t0 = clock();
		double   f;
		Utilities::mergeVec(n,p,x);							// x = [n p]
		isSucceed &= m_fem->computeValueAndGrad(x, tVec, &f, &G);
		int t1 = clock();
		PRINT_F("compute value and grad:%f", (t1-t0)/1000.f);
		if (!isSucceed)
			return false;

		Eigen::Map<EigVec> Gn(&G[0], nIntDof);				// G = [Gn Gp]
		Eigen::Map<EigVec> Gp(&G[0]+nIntDof, nParam);

		// 计算Hessian
		EigDense* pHpp = &Hpp;//ithIter == 0 ? &Hpp : NULL;
		isSucceed &= m_fem->computeHessian(n,p, t, Hnn, Hnp, Hpn, pHpp);
		int t2 = clock();
		PRINT_F("compute hessian:%f",(t2-t1)/1000.f);
		if (!isSucceed)
			return false;

		if (ithIter > 0 && 0)
		{
			// 利用BFGS方法计算Hpp
			// yk   = [dGn dGp]
			// sk   = [dn  dp]
			// Hk+1 = Hk - Hk*sk*sk^T*Hk / (sk^T*Hk*sk) + yk*yk^T / (yk^T*sk)
			EigVec y = G - G0;
			EigVec&s = dx;

			double sHs= dN.dot(dGN) + dP.dot(dGP);
			double yTs = 1.0 / (y.dot(s));

			for (int i = 0; i < nParam; ++i)
			{
				for (int j = 0; j < nParam; ++j)
				{
					Hpp(i,j) += y(i+nIntDof) * y(j+nIntDof) / yTs;
					Hpp(i,j) -= Gp(i) * Gp(j) / sHs;
				}
			}
		}

		Eigen::SuperLU<EigSparse> solver;
		solver.compute(Hnn);
		if(solver.info()!=Eigen::Success) 
		{
			PRINT_F("LU factorization FAILED!\n");
			return false;
		}
		EigDense invHnnHnp(nIntDof, nParam);
		for (int ithParam = 0; ithParam < nParam; ++ithParam)
		{
			EigVec c = solver.solve(Hnp.col(ithParam));
			invHnnHnp.col(ithParam) = c;
		}
		EigVec   invHnnN = solver.solve(Gn);

		EigDense A = Hpp - Hpn * invHnnHnp;
		EigDense b = -Gp + Hpn * invHnnN;
		dP= A.colPivHouseholderQr().solve(b);

		EigVec   b2= -Gn - Hnp * dP;
		dN= solver.solve(b2);
		int t3 = clock();
		PRINT_F("compute dN, dP:%f",(t3-t2)/1000.f);

		// 进行一维搜索
		Utilities::mergeVec(dN, dP, dx);
		double df = G.dot(dx);
		double a;
		//m_lineSearch.setInitStep(3.05);
		int res = m_lineSearch.lineSearch(x,dx,tVec,a, &f, &df);
		if (res == 0 && a > 1e-4)
		{
			PRINT_F("a = %lf", a);
			if (a != 1.0)
			{
				dN *= a;
				dP *= a;
			}
		}
		else
		{
			PRINT_F("line search failed: a = %lf", a);
		}
		int t4 = clock();
		PRINT_F("line search:%f", (t4-t3)/1000.f);

		// 计算残差
		dGN = Hnn * dN + Hnp * dP;
		dGP = Hpn * dN + Hpp * dP;
		EigVec resiN = dGN + Gn;
		EigVec resiP = dGP + Gp;
		resiN = resiN.cwiseAbs();
		resiP = resiP.cwiseAbs();

		// 更新自变量
		n += dN;
		p += dP;

		PRINT_F("|dN|=%le |dP|=%le |resiN|∞=%le |resiP|∞=%le", dN.norm(), dP.norm(), resiN.maxCoeff(), resiP.maxCoeff());

		
		if (resiN.maxCoeff() < 1e-5 && resiP.maxCoeff() < 1e-5)
			break;
		if (dN.norm() < 1e-12 && dP.norm() < 1e-15)
			break;
			

		G0 = G;
	}
	// 更新为最终迭代状态
	m_fem->setDof(n, p);

	// 记录新状态
	int nTotParam = m_fem->m_transRig->getNParam();
	EigVec totParam(nTotParam);
	m_fem->m_transRig->getParam(&totParam[0]);
	m_paramResult.push_back(totParam);

	// 输出状态
	PRINT_F("time: %.2lf", t);
	PRINT_F("new params:");
	for (int i = 0; i < nParam; ++i)
	{
		PRINT_F("%.2lf ", p[i]);
	}
	PRINT_F("newton iteration end.");

	return true;
}

RigFEM::NewtonSolver::NewtonSolver( RiggedMesh* fem ) :m_fem(fem), m_maxIter(10), m_lineSearch(fem)
{

}

void RigFEM::NewtonSolver::saveResult( const char* fileName, const char* paramName /*= "param"*/ )
{
	ofstream file(fileName);
	file << paramName << "=[\n";
	for (int i = 0; i < m_paramResult.size(); ++i)
	{
		EigVec& vec = m_paramResult[i];
		for (int ithEle = 0; ithEle < vec.size(); ++ithEle)
		{
			file << vec[ithEle] << ' ';
		}
		file << ";\n";
	}
	file << "];";
	file.close();
}

void RigFEM::NewtonSolver::setMesh( RiggedMesh* fem )
{
	m_fem = fem;
	m_lineSearch = LineSearch(fem);
}

bool RigFEM::SimpleFunction::computeValueAndGrad( const EigVec& x, double* v /*= NULL*/, EigVec* grad /*= NULL*/ )
{
	double xVal = x[0];

	if (v)
	{
		*v = computeFuncVal(xVal);
	}
	if (grad)
	{
		// 中心差商计算梯度（导数）
		double e = 1e-4;
		grad->resize(1);
		double f0 = computeFuncVal(xVal-e);
		double f1 = computeFuncVal(xVal+e);
		(*grad)[0] = (f1 - f0) / (2*e);
	}
	return true;
}

double RigFEM::SimpleFunction::computeFuncVal( double x )
{
	double coefA = 0.35;
	double coefB = -1.63;
	double coefC = 1.82;
	double coefD = 0.25;

	double x2 = x*x;
	double x3 = x*x2;
	return coefA * x3 + coefB * x2 + coefC * x + coefD;
}
