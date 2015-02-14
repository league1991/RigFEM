#include "StdAfx.h"
#include "newtonSolver.h"
#include "RiggedSkinMesh.h"
using namespace RigFEM;

LineSearcher::~LineSearcher(void)
{
}

int RigFEM::LineSearcher::guessMinPnt( double a0, double fa0, double *dfa0, double a1, double fa1, double *dfa1, double& aMin )
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
	// 导数插值
	if (dfa0 && dfa1 && *dfa0 * *dfa1 < 0)
	{
		a = (a0 * *dfa1 - a1 * *dfa0) / (*dfa1 - *dfa0);
		if ((a0 - a) * (a1 - a) < 0)
		{
			aMin = a;
			return 0;
		}
	}
	// 各种拟合方法都失败,取中点
	aMin = 0.5 * (a0 + a1);
	return 1;
}

void RigFEM::LineSearcher::setOriginAndDir( const EigVec& x0, const EigVec& dx )
{
	m_x0 = x0;
	m_dx = dx;
}

int RigFEM::LineSearcher::zoom( double al, double fal, double* dfal, double ah, double fah, double* dfah, double& a )
{
	int maxIter = m_maxZoomIter;
	int problem = LS_NONE;
	while(maxIter--)
	{
		double aj;
		guessMinPnt(al, fal, dfal, ah, fah, dfah, aj);

		// 保证aj在al和ah之间

		double fj;
		EigVec xj = m_x0 + m_dx * aj;
		computeValueAndDeri(xj, &fj, NULL);

		problem = LS_NONE;
		if (fj > m_f0 + m_c1 * aj * m_df0 || fj >= fal)
		{
			// 插值后得出的新点不满足wolfe充分下降条件，或函数值比fal大
			ah = aj;
			problem |= LS_WOLFE_DESCENT;	
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
				return LS_NONE;
			}

			problem |= LS_WOLFE_CURVATURE;
			if (dfj * (ah - al) >= 0)
			{
				ah = al;
			}
			al = aj;
		}
	}

	// 没有找到符合wolfe条件的新点
	a = 0.5 * (al + ah);
	return problem;
}

int RigFEM::LineSearcher::lineSearch( const EigVec& x0, const EigVec& dx, const EigVec& param, double& aFinal,
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
	int problem = LS_NONE;
	for (int i = 1; ai < amax; ++i)
	{
		//PRINT_F("%dth line search iter, a = %lf", i-1,ai);
		EigVec x = m_x0 + m_dx * ai;

		// 只计算函数值，尽量提高效率
		double fi;
		double dfi;
		int t0 = clock();
		computeValueAndDeri(x, &fi, NULL);
		int t1 = clock();
		//PRINT_F("line search compute value and deri #1: %f", (t1-t0)/1000.f);
		problem = LS_NONE;
		if (fi > m_f0 + m_c1 * ai * m_df0)
			problem |= LS_WOLFE_DESCENT;
		// 若wolfe充分下降条件被违反或函数值较上次增长,终止迭代
		if (problem & LS_WOLFE_DESCENT ||
			(i > 1 && fi - fi_1 > 1e-8))
		{
			int res = zoom(ai_1, fi_1, &dfi_1, ai, fi, NULL, aFinal);
			return res;
		}

		computeValueAndDeri(x, NULL, &dfi);
		int t2 = clock();
		//PRINT_F("line search compute value and deri #2: %f", (t2-t1)/1000.f);

		// wolfe充分下降条件被满足，函数值较上次下降				
		if (abs(dfi) <= m_c2 * abs(m_df0))
		{	
			// wolfe曲率条件被满足,直接退出
			aFinal = ai;
			return 0;
		}
		problem |= LS_WOLFE_CURVATURE;

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
	return problem;
}

bool RigFEM::LineSearcher::computeValueAndDeri( const EigVec& x, double* f, double* df )
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

RigFEM::LineSearcher::LineSearcher( ObjectFunction* objFun /*= NULL*/, double initStep /*= 1.0*/, double maxStep /*= 10*/ ) :
m_objFunc(objFun), 
m_initStep(initStep), m_maxStep(maxStep), m_maxZoomIter(15),
m_c1(1e-4), m_c2(0.9)
{

}

bool RigFEM::LineSearcher::setC( double c1, double c2 )
{
	if (0 < c1 && c1 < c2 && c2 < 1)
	{
		m_c1 = c1;
		m_c2 = c2;
		return true;
	}
	return false;
}

bool RigFEM::PointParamSolver::step()
{
	PRINT_F("############################################ newton iteration begin. ############################################");
	m_fem->setStatus(m_initStatus);
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

	char fileName[100];
	sprintf(fileName, "I:/Programs/VegaFEM-v2.1/myProject/RigPlugin/%dFrame.m", m_paramResult.size());
	//ofstream file(fileName);
	for (int ithIter = 0; ithIter < maxIter; ++ithIter)
	{
		// 计算当前q p值的函数值、梯度值和Hessian
		PRINT_F("######################### %dth Iteration #########################", ithIter);
		int t0 = clock();
		double   f;
		Utilities::mergeVec(n,p,x);							// x = [n p]
		isSucceed &= m_fem->computeValueAndGrad(x, tVec, &f, &G);
		int t1 = clock();
		//PRINT_F("compute value and grad:%f", (t1-t0)/1000.f);
		if (!isSucceed)
			return false;

		Eigen::Map<EigVec> Gn(&G[0], nIntDof);				// G = [Gn Gp]
		Eigen::Map<EigVec> Gp(&G[0]+nIntDof, nParam);

		// 如果梯度接近0，终止迭代
		double GnNorm = Gn.norm();
		double GpNorm = Gp.norm();
		PRINT_F("f = %lf, |Gn| = %lf, |GP| = %lf, minGradSize = %e", f, GnNorm, GpNorm, m_minGradSize);
		if(GnNorm < m_minGradSize && GpNorm < m_minGradSize)
			break;
		// 如果步长太小，也终止迭代
		if (ithIter > 0 && dN.norm() < m_minStepSize && dP.norm() < m_minStepSize)
			break;

		// 计算Hessian
		EigDense* pHpp = &Hpp;//ithIter == 0 ? &Hpp : NULL;
		isSucceed &= m_fem->computeHessian(n,p, t, Hnn, Hnp, Hpn, pHpp);
		int t2 = clock();
		//PRINT_F("compute hessian:%f",(t2-t1)/1000.f);
		if (!isSucceed)
			return false;
		// 若某参数梯度为0，此时为预防方程退化，强制修改方程系数，令该参数增量为0
// 		for (int ithParam = 0; ithParam < nParam; ++ithParam)
// 		{ 
// 			if (abs(Gp[ithParam]) < 1e-4)
// 			{
// 				const double v = 1e5;
// 				for (int ithIntDof = 0; ithIntDof < nIntDof; ++ithIntDof)
// 				{
// 					Hnp(ithIntDof, ithParam) = v;
// 				}
// 				for (int jthParam = 0; jthParam < nParam; ++jthParam)
// 				{
// 					Hpp(jthParam, ithParam) = v;
// 				}
// 			}
// 		}
// 		string HnnStr = Utilities::matToString(Hnn, "Hnn");
// 		string HnpStr = Utilities::matToString(Hnp, "Hnp");
// 		string HpnStr = Utilities::matToString(Hpn, "Hpn");
// 		string HppStr = Utilities::matToString(Hpp, "Hpp");
// 		file << HnnStr;
// 		file << HnpStr;
// 		file << HpnStr;
// 		file << HppStr;
		

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
		string Astr = Utilities::matToString(A, "A");
		string bStr = Utilities::matToString(b, "b");
		//file << Astr << bStr;

		EigVec   b2= -Gn - Hnp * dP;
		dN= solver.solve(b2);
		//Global::showVector(Gp, "Gp");
		//Global::showVector(dP, "dP");
		int t3 = clock();
		//PRINT_F("compute dN, dP:%f",(t3-t2)/1000.f);

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
		//PRINT_F("line search:%f", (t4-t3)/1000.f);

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
		PRINT_F("minStepSize = %e", m_minStepSize);
		//Global::showVector(p, "p");
		
		//Utilities::saveVector(n, "I:/Programs/VegaFEM-v2.1/myProject/RigPlugin/ExperimentData/1.16/N.txt");
		//Utilities::saveVector(p, "I:/Programs/VegaFEM-v2.1/myProject/RigPlugin/ExperimentData/1.16/P.txt");
		//Utilities::loadVector(n, "I:/Programs/VegaFEM-v2.1/myProject/RigPlugin/ExperimentData/1.16/N_nodir.txt");
		//Utilities::loadVector(p, "I:/Programs/VegaFEM-v2.1/myProject/RigPlugin/ExperimentData/1.16/P_nodir.txt");

		G0 = G;
		PRINT_F("************ load N P ***************");
	}
	//file.close();
	// 更新为最终迭代状态
	m_fem->setDof(n, p);

	// 记录新状态
	int nTotParam = m_fem->getRigObj()->getNParam();
	EigVec totParam(nTotParam);
	m_fem->getRigObj()->getParam(&totParam[0]);
	m_paramResult.push_back(totParam);
	m_finalStatus = m_fem->getStatus();

	// 输出状态
	PRINT_F("time: %.2lf", t);
// 	PRINT_F("new params:");
// 	for (int i = 0; i < nParam; ++i)
// 	{
// 		PRINT_F("%.2lf ", p[i]);
// 	}
	PRINT_F("############################################# newton iteration end. #############################################");

	return true;
}

RigFEM::PointParamSolver::PointParamSolver( RiggedMesh* fem ) :m_fem(fem), m_lineSearch(fem)
{

}

void RigFEM::PointParamSolver::saveResult( const char* fileName, const char* paramName /*= "param"*/ )
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

void RigFEM::PointParamSolver::setMesh( RiggedMesh* fem )
{
	m_fem = fem;
	m_lineSearch = LineSearcher(fem);
}

void RigFEM::NewtonSolver::setTerminateCond( 
	int maxIter, double minStepSize, double minGradSize,
	int maxCGIter, double minCGStepSize)
{
	m_maxIter = maxIter;
	m_minStepSize = minStepSize;
	m_minGradSize = minGradSize;
	m_maxCGIter = maxCGIter;
	m_minCGStepSize = minCGStepSize;
}

RigFEM::NewtonSolver::NewtonSolver() :
m_minStepSize(1e-3), m_minGradSize(1e-2), m_maxIter(10),
m_maxCGIter(10), m_minCGStepSize(1e-2),
m_iterMaxStepSize(1e0)
{

}

bool RigFEM::NewtonSolver::setInitStatus( const RigStatus&s )
{
	m_initStatus = s;
	return true;
}

void RigFEM::NewtonSolver::setIterationMaxStepSize( double maxStep )
{
	m_iterMaxStepSize = maxStep;
}



RigFEM::PointParamSolver::~PointParamSolver()
{

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

bool RigFEM::ParamSolver::step()
{
	PRINT_F("############################################ newton iteration begin. ############################################");
	m_fem->setStatus(m_initStatus);
	EigVec P = m_initStatus.getP();
	EigVec PkPk_1;
	if (m_initStatus.getCustom("Pi-Pi-1", PkPk_1))
	{
		// 根据上一帧的变化趋势推测参数初始值
		//Global::showVector(PkPk_1,"dP0");
		P += PkPk_1;
	}
	double t = m_fem->getCurTime();
	int nParam  = P.size();
	EigVec tVec(1);
	tVec[0] = t;

	bool isSucceed = true;
	RigStatus	resultStatus;
	EigDense  H;
	EigVec	 G,  dP, dG;
	EigVec	 G0, dP0;				// 上一次迭代的值
	double    initStepSize;
	if (!m_initStatus.getCustom("initStep", initStepSize))
		initStepSize = 1e-6;
	double	  stepSize = initStepSize;

	m_lineSearch.setC(1e-4, 0.45);
	m_lineSearch.setMaxZoomIter(15);
	const double minPStep = 1e-2, maxPStep = 10;		// 限制每步参数变化范围
	double iterMaxStep = m_iterMaxStepSize;
	for (int ithIter = 0; ithIter < m_maxCGIter; ++ithIter)
	{
		PRINT_F("##################### %dth conjugate gradient ####################", ithIter);
		double f;
		isSucceed &= m_fem->computeValueAndGrad(P, tVec, &f, &G);
		double GNorm = G.norm();
		PRINT_F("f = %lf, |GP| = %lf, minGradSize = %e", f, GNorm, m_minGradSize);
		//Global::showVector(G, "G");
		//Global::showVector(P, "p");

		if(G.norm() < m_minGradSize * 20) 
			break;
		if (ithIter != 0 && dP.norm() < m_minCGStepSize)
			break;

		if (ithIter == 0)
			dP = -G;
		else
		{
			dP = -G + dP0 * G.dot(G) / G0.dot(G0);
			if (G.dot(dP) > 0)
				dP = -G;
		}

		// 选择合适的初始步长
		if (ithIter != 0)
			initStepSize = stepSize * G0.dot(dP0) / G.dot(dP);
		double dPMax = dP.cwiseAbs().maxCoeff();
		initStepSize = CLAMP_DOUBLE(minPStep/dPMax, maxPStep/dPMax, initStepSize);
		m_lineSearch.setInitStep(initStepSize);

		// 一维搜索
		double df = G.dot(dP);					// 此为自变量变化G时函数增量
		int res = m_lineSearch.lineSearch(P, dP, tVec, stepSize, &f, &df);
		if (ithIter == 0 && res == LineSearcher::LS_NONE)
		{
			resultStatus.addOrSetCustom("initStep", CLAMP_DOUBLE(1e-9, 1e-1, stepSize));
		}
		if (res == 0)
		{
			PRINT_F("line search succeed, step = %le", stepSize);
		}
		else
		{
			PRINT_F("line search FAILED");
		}

		dP *= stepSize;
		MathUtilities::clampVector(dP, iterMaxStep);
		P += dP;
		PRINT_F("|dP|=%le", dP.norm());

		G0 = G;
		dP0 = dP;
	}

	m_lineSearch.setC(1e-4, 0.9);
	for (int ithIter = 0; ithIter < m_maxIter; ++ithIter)
	{
		// 计算当前q p值的函数值、梯度值和Hessian
		PRINT_F("#########################  %dth Newton   #########################", ithIter);
		int t0 = clock();
		double   f;
		isSucceed &= m_fem->computeValueAndGrad(P, tVec, &f, &G);
		int t1 = clock();
		//PRINT_F("compute value and grad:%f", (t1-t0)/1000.f);
		if (!isSucceed)
			return false;

		// 如果梯度接近0，终止迭代
		double GNorm = G.norm();
		//Global::showVector(G, "G");
		PRINT_F("f = %lf, |GP| = %lf, minGradSize = %e", f, GNorm, m_minGradSize);
		if(GNorm < m_minGradSize)
			break;
		// 如果步长太小，也终止迭代
		if (ithIter > 0 && dP.norm() < m_minStepSize)
			break;

		// 计算Hessian
		isSucceed &= m_fem->computeHessian(P, tVec, H);
		int t2 = clock();
		//Global::showMatrix(H, "H");
		//PRINT_F("compute hessian:%f",(t2-t1)/1000.f);
		if (!isSucceed)
			return false;

		dP = H.colPivHouseholderQr().solve(-G);
		//Global::showVector(dP, "dP");
		int t3 = clock();
		//PRINT_F("compute dN, dP:%f",(t3-t2)/1000.f);

		// 进行一维搜索
		double df = G.dot(dP);		// 自变量变化dP时函数增量
		if (df > 0)
		{
			PRINT_F("failed to find newton direction, use gradient directly");
			dP = -G;
			m_lineSearch.setInitStep(initStepSize);
			df = G.dot(dP);
		}
		else
		{
			m_lineSearch.setInitStep(1.0);
		}
		double a = 1.0;
		//m_lineSearch.setInitStep(3.05);
		int res = m_lineSearch.lineSearch(P,dP,tVec,a, &f, &df);
		if (res == 0)
		{
			PRINT_F("a = %lf", a);
			if (a != 1.0)
			{
				dP *= a;
			}
		}
		else
		{
			PRINT_F("line search failed: a = %lf", a);
		}
		MathUtilities::clampVector(dP, iterMaxStep);
		int t4 = clock();
		//PRINT_F("line search:%f", (t4-t3)/1000.f);

		// 计算残差
		dG = H * dP;
		EigVec resiP = dG + G;
		resiP = resiP.cwiseAbs();

		// 更新自变量
		P += dP;

		PRINT_F("|dP|=%le |resiP|∞=%le", dP.norm(), resiP.maxCoeff());
		PRINT_F("minStepSize = %e", m_minStepSize);
		//Global::showVector(P, "p");

		G0 = G;
	}
	//file.close();
	// 更新为最终迭代状态
	m_fem->setDof(P);

	// 记录新状态
	m_paramResult.push_back(P);
	m_finalStatus = m_fem->getStatus();
	m_finalStatus.mergeCustom(resultStatus);
	m_finalStatus.addOrSetCustom("Pi-Pi-1", m_finalStatus.getP() - m_initStatus.getP());

	// 输出状态
	PRINT_F("time: %.2lf", t);
	//Global::showVector(P, "p");
	PRINT_F("############################################# newton iteration end. #############################################");

	return true;
}

RigFEM::ParamSolver::ParamSolver( RiggedSkinMesh* fem /*= NULL*/ ) :m_fem(fem), m_lineSearch(fem)
{

}
