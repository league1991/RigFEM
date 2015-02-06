#include "stdafx.h"

using namespace RigFEM;
void RigFEM::RiggedSkinMesh::setWeight( EigSparse& sparse )
{
	m_weightMat = sparse;
	m_weightMatTran = sparse.transpose();
}

bool RigFEM::RiggedSkinMesh::setWeight( const char* weightFile )
{
	std::map<string, EigDense> matMap;
	Utilities::fileToDense(weightFile, matMap);
	if (matMap.find("weight") == matMap.end())
		return false;

	EigSparse weightMat;
	Utilities::eigDense2Sparse(matMap["weight"], weightMat);
	if (weightMat.rows() == m_intDofIdx.size() &&
		weightMat.cols() == m_surfDofIdx.size())
	{
		m_weightMat = weightMat;
	}
	else if (weightMat.rows() * 3 == m_intDofIdx.size() &&
		weightMat.cols() * 3 == m_surfDofIdx.size())
	{
		Utilities::kronecker3X(weightMat, m_weightMat);
	}
	else
		return false;
	m_weightMatTran = m_weightMat.transpose();
	return true;
}

bool RigFEM::RiggedSkinMesh::computeSkinQ( const double* p, double t, double* q )
{
	EigVec s;
	bool res = computeSurfOffset(p, t, s);
	if (!res)
		return res;

	// 求出内部点坐标
	EigVec n = m_weightMat * s;
	for (int i = 0; i < m_nIntPnt; ++i)
	{
		int idx = m_intPntIdx[i];
		q[idx*3]   = n[i*3];
		q[idx*3+1] = n[i*3+1];
		q[idx*3+2] = n[i*3+2];
	}
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx = m_surfPntIdx[i];
		q[idx*3]   = s[i*3];
		q[idx*3+1] = s[i*3+1];
		q[idx*3+2] = s[i*3+2];
	}
	return true;
}

bool RigFEM::RiggedSkinMesh::computeValueAndGrad( const EigVec& x, const EigVec& param, double* v, EigVec* grad )
{
	bool res = true;
	EigVec q(m_nTotPnt*3);
	res &= computeSkinQ(&x[0], param[0], &q[0]);
	if (!res)
		return false;

	// 计算函数值
	if (v)
	{
		*v = computeValue(q);
	}
	if (grad)
	{
		EigVec ma = (q - m_q) / (m_h * m_h) - m_v / m_h;
		for (int i = 0; i < m_nTotPnt; ++i)
		{
			ma[i*3]   *= m_mass[i];
			ma[i*3+1] *= m_mass[i];
			ma[i*3+2] *= m_mass[i];
		}

		// 计算此状态下的力
		EigVec f(m_nTotPnt*3);
		m_forceModel->GetInternalForce(&q[0], &f[0]);
		f *= -1.0;
		EigVec residual = ma - f;

		EigVec gn(m_nIntPnt*3);
		for (int ithN = 0; ithN < m_nIntPnt; ++ithN)
		{
			int idx = m_intPntIdx[ithN];
			gn[ithN*3]   = residual[idx*3];
			gn[ithN*3+1] = residual[idx*3+1];
			gn[ithN*3+2] = residual[idx*3+2];
		}

		EigVec gs(m_nSurfPnt*3);
		for (int ithS = 0; ithS < m_nSurfPnt; ++ithS)
		{
			int idx = m_surfPntIdx[ithS];
			gs[ithS*3]   = residual[idx*3];
			gs[ithS*3+1] = residual[idx*3+1];
			gs[ithS*3+2] = residual[idx*3+2];
		}

		EigDense J;
		res &= m_transRig->computeJacobian(J);
		*grad = J.transpose() * (gs + m_weightMatTran * gn);
	}
	return res;
}

bool RigFEM::RiggedSkinMesh::computeHessian( const EigVec&x, const EigVec& param, EigDense& H )
{
	// 计算当前各个自由度的值q
	bool res = true;
	EigVec q(m_nTotPnt*3);
	res &= computeSkinQ(&x[0], param[0], &q[0]);

	// 提取内力、tangent Stiffness matrix、质量矩阵
	// 注意tangent stiffness matrix为实际的负值，因为系统计算出的弹力为实际的负值,因此需要先反转
	EigVec force(m_nTotPnt*3);
	if(!m_tangentStiffnessMatrix)
		m_forceModel->GetTangentStiffnessMatrixTopology(&m_tangentStiffnessMatrix);
	m_forceModel->GetForceAndMatrix(&q[0], &force[0], m_tangentStiffnessMatrix);
	//*m_tangentStiffnessMatrix *= -1;
	force *= -1;

	// 计算 b = ms*as - fs + W^T * (mq*aq - fq)
	EigVec ma = (q - m_q) / (m_h * m_h) - m_v / m_h;
	for (int i = 0; i < m_nTotPnt; ++i)
	{
		ma[i*3]   *= m_mass[i];
		ma[i*3+1] *= m_mass[i];
		ma[i*3+2] *= m_mass[i];
	}
	EigVec residual = ma - force;
	EigVec gn(m_nIntPnt*3);
	EigVec gs(m_nSurfPnt*3);
	for (int ithN = 0; ithN < m_nIntPnt*3; ++ithN)
	{
		int idx = m_intDofIdx[ithN];
		gn[ithN]   = residual[idx];
	}
	for (int ithS = 0; ithS < m_nSurfPnt*3; ++ithS)
	{
		int idx = m_surfDofIdx[ithS];
		gs[ithS]   = residual[idx];
	}
	EigVec b = gs + m_weightMatTran * gn;

	// 计算雅可比矩阵
	EigDense J;
	res &= m_transRig->computeJacobian(J);

	// 计算Hessian第一项
	// dJki/dp * bk
	H.resize(m_nParam, m_nParam);
	EigVec jacobianDeri(m_nSurfPnt*3);	
	for (int i = 0; i < m_nParam; ++i)
	{
		for (int l = 0; l < m_nParam; ++l)
		{
			double& Hval = H(i,l);
			res &= m_transRig->computeJacobianDerivative(i,l,&jacobianDeri[0]);
			Hval = jacobianDeri.dot(b);
		}
	}

	// 计算Hessian 第二项
	// Jki * dbk/dp
	// 提取 dFn/dn dFn/ds dFs/dn dFs/ds,
	// 其中Fn为内部节点受到的力，Fs为表面节点受到的力
	// n为内部节点位置，s为表面节点位置
	EigSparse dFss, dFsn, dFns, dFnn;
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_surfDofIdx, dFss);
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_intDofIdx,  dFsn);
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_intDofIdx,  m_surfDofIdx, dFns);
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_intDofIdx,  m_intDofIdx,  dFnn);

	// 计算 db/dp 
	// 大小为(表面点数*3，参数数)
	double h2 = m_h * m_h;
	for (int i = 0; i < m_surfPntIdx.size(); ++i)
	{
		int idx = m_surfPntIdx[i];
		double m = m_mass[idx] / h2;
		int i3 = i*3;
		dFss.coeffRef(i3,i3) += m;	++i3;
		dFss.coeffRef(i3,i3) += m;	++i3;
		dFss.coeffRef(i3,i3) += m;
	}
	for (int i = 0; i < m_intPntIdx.size(); ++i)
	{
		int idx = m_intPntIdx[i];
		double m = m_mass[idx] / h2;
		int i3 = i*3;
		dFnn.coeffRef(i3,i3) += m;	++i3;
		dFnn.coeffRef(i3,i3) += m;	++i3;
		dFnn.coeffRef(i3,i3) += m;
	}
	EigDense JT     = J.transpose();
	EigDense JTdFss = JT * dFss;
	EigDense JTdFsn = JT * dFsn;
	EigDense dFnsJ  = dFns * J;
	EigDense JTWT   = JT * m_weightMatTran;
	EigDense WJ     = m_weightMat * J;
	EigDense term2  = JTdFss * J + JTWT * (dFnsJ +  dFnn * WJ) + JTdFsn * WJ;
	//EigDense term20 = J.transpose() * (dFss + m_weightMatTran * (dFns + dFnn*m_weightMat) + dFsn*m_weightMat) * J;

	// 计算Hessian第二项Jki * dbk/dp
	H += term2;
	return true;
}

bool RigFEM::RiggedSkinMesh::testCurFrameGrad( RigStatus& lastFrame, RigStatus& curFrame, double noiseP /*= 1.0*/ )
{
	PRINT_F("Test Gradient and Function Value");
	PRINT_F("noiseP = %lf", noiseP);
	setStatus(lastFrame);
	EigVec p = curFrame.getP();

	EigVec dP = EigVec::Random(m_nParam)* (noiseP * 2) - EigVec::Constant(m_nParam, noiseP);
	EigVec gp;
	EigVec tParam(1);
	tParam[0] = m_t;
	double f0;
	computeValueAndGrad(p, tParam, &f0, &gp);
	PRINT_F("funcVal:%lf gradient: |gp| = %lf", f0, gp.norm());

	for (double step = 1; step > 1e-15; step *= 0.1)
	{
		EigVec dPi = dP * step;

		// 计算近似的函数值
		double fi_app = f0 + gp.dot(dPi);

		// 计算准确的函数值
		EigVec pi = p + dPi;
		double fi;
		computeValueAndGrad(pi, tParam, &fi, NULL);

		double error = abs(fi-fi_app);
		PRINT_F("step = %le funVal = %le approxVal = %le error = %le error/dx^2 = %le", step, fi, fi_app, error, error/(step*step));
	}
	PRINT_F("\n");
	return true;
}

bool RigFEM::RiggedSkinMesh::testCurFrameHessian( RigStatus& lastFrame, RigStatus& curFrame, double noiseP/*=1.0*/ )
{
	// 当前帧的配置
	PRINT_F("Test Gradient and Hessian");
	PRINT_F("noiseP = %lf", noiseP);
	setStatus(lastFrame);
	EigVec p;
	p = curFrame.getP();

	EigVec dP = EigVec::Random(m_nParam) * (noiseP * 2.0) - EigVec::Constant(m_nParam, noiseP);
	EigVec gp;
	EigVec tParam(1);
	tParam[0] = m_t;
	EigDense H;

	computeValueAndGrad(p, tParam, NULL, &gp);
	computeHessian(p, tParam, H);
	for (double step = 1; step > 1e-15; step *= 0.1)
	{
		EigVec dPi = dP * step;

		// 计算近似的梯度值
		EigVec gpi_app = gp + H * dPi;

		// 计算准确的梯度值
		EigVec pi = p + dPi;
		EigVec gpi;
		computeValueAndGrad(pi, tParam, NULL, &gpi);

		EigVec resP = gpi - gpi_app;

		double resPNorm = resP.norm();
		resP = resP.cwiseAbs();
		double resPMax  = resP.maxCoeff();

		double invE = 1.0 / (step * step);
		PRINT_F("ε=%le  maxResP=%le  |maxP|/ε^2=%le",
			step, resPMax, resPMax* invE);
	}
	return true;
}

void RigFEM::RiggedSkinMesh::setDof( EigVec& p, bool proceedTime /*= true */ )
{
	EigVec q(m_nTotPnt*3);
	computeSkinQ(&p[0], m_t, &q[0]);

	EigVec v = (q - m_q) / m_h;
	EigVec a = (v - m_v) / m_h; 
	m_param = p;
	m_v = v;
	m_a = a;
	m_q = q;
	if (proceedTime)
	{
		m_t += m_h;
	}
}

void RigFEM::RiggedSkinMesh::getDof( EigVec& p )
{
	p = m_param;
}

