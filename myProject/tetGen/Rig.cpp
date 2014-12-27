#include "StdAfx.h"
#include "Rig.h"
using namespace RigFEM;

RigBase::RigBase(int nParam):m_t(0), m_keyFrameFunc(NULL), m_nParam(nParam)
{
	m_param = new double[nParam];
	m_keyFrameFunc = new KeyFrameFunc[nParam];
	for (int i = 0; i < nParam; ++i)
	{
		m_keyFrameFunc[i] = NULL;
		m_param[i] = 0;
	}
}

RigBase::~RigBase(void)
{
	delete[] m_param;
	delete[] m_keyFrameFunc;
}



void RigFEM::TransformRig::setInitPnts( double* pnts, int nPnts )
{
	m_initPntList.clear();
	Vec3d totPnt(0.0);
	for (int i = 0; i < nPnts; ++i, pnts += 3)
	{
		m_initPntList.push_back(Vec3d(pnts));
		totPnt += Vec3d(pnts);
	}
	m_localCenter = totPnt / nPnts;
}

void RigFEM::TransformRig::transform()
{
	m_translation = Vec3d(m_param);
	m_rotate      = Vec3d(m_param+3);
	m_scale       = Vec3d(m_param+6);
	m_curPntList.resize(m_initPntList.size());
	for (int ithPnt = 0; ithPnt < m_initPntList.size(); ++ithPnt)
	{
		Vec3d p = m_initPntList[ithPnt] - m_localCenter;
		p[0] *= m_scale[0];
		p[1] *= m_scale[1];
		p[2] *= m_scale[2];

		p = rotateX(p, m_rotate[0]);
		p = rotateY(p, m_rotate[1]);
		p = rotateZ(p, m_rotate[2]);

		p += m_localCenter;
		m_curPntList[ithPnt] = p + m_translation;
	}
}



void RigFEM::TransformRig::compute( const Vec3d& trans, const Vec3d& scl )
{
	m_translation = trans;
	m_scale = scl;
	transform();
}

void RigFEM::TransformRig::compute()
{
	transform();
}

void RigFEM::TransformRig::setAllParam( const Vec3d& trans, const Vec3d& rotate, const Vec3d& scale )
{
	m_translation = trans;
	m_rotate = rotate;
	m_scale = scale;
	m_translation.convertToArray(m_param);
	m_rotate.convertToArray(m_param+3);
	m_scale.convertToArray(m_param+6);
}

/*
void RigFEM::TransformRig::computeJacobian( EigDense& J )
{
	EigDense jacobian = EigDense::Zero(m_initPntList.size()*3, 6);
	for (int ithPnt = 0; ithPnt < m_initPntList.size(); ++ithPnt)
	{
		int idx = ithPnt * 3;
		jacobian(idx,0) = jacobian(idx+1,1) = jacobian(idx+2,2) = 1;

		Vec3d del = m_initPntList[ithPnt] - m_localCenter;
		jacobian(idx,3) = del[0];	idx++;
		jacobian(idx,4) = del[1];	idx++;
		jacobian(idx,5) = del[2];
	}

	int nFreeParam = getNFreeParam();
	J = EigDense(m_initPntList.size() * 3, nFreeParam);
	for (int i = 0, ithFree = 0; i < m_nParam; ++i)
	{
		if (!m_keyFrameFunc[i])
		{
			J.col(ithFree) = jacobian.col(i);
			ithFree++;
		}
	}
}

void RigFEM::TransformRig::computeJacobianDerivative( int i, int j, double* res )
{
	int nDoF   = getNPoints() * 3;
	for (int i = 0; i < nDoF; ++i)
	{
		res[i] = 0.0;
	}
}*/

bool RigFEM::TransformRig::computeValue( double* result, const double* params /*= 0*/ )
{
	if (params)
	{
		setFreeParam(params);
	}
	transform();
	for (int i = 0; i < m_curPntList.size(); ++i)
	{
		result[i*3] = m_curPntList[i][0];
		result[i*3+1] = m_curPntList[i][1];
		result[i*3+2] = m_curPntList[i][2];
	}
	return true;
}

RigFEM::TransformRig::TransformRig() :RigBase(9), m_translation(0.0), m_localCenter(0.0),m_scale(1.0), m_rotate(0.0)
{
	m_translation.convertToArray(m_param);
	m_rotate.convertToArray(m_param+3);
	m_scale.convertToArray(m_param+6);
}


bool RigFEM::RigBase::computeJacobianDerivative( int i, int j, double* res )
{
	double e = 1e-4;

	int nParam = getNFreeParam();
	int nDoF   = getNPoints() * 3;

	vector<double> p(nParam);
	getFreeParam(&p[0]);
	double oldPi = p[i];
	double oldPj = p[j];
	vector<double> s00(nDoF), s01(nDoF), s10(nDoF), s11(nDoF);
	bool isSucceed = true;

	// p00
	p[i] = oldPi - e;
	p[j] = oldPj - e;
	isSucceed &= computeValue(&s00[0], &p[0]);

	// p01
	p[i] = oldPi - e;
	p[j] = oldPj + e;
	isSucceed &= computeValue(&s01[0], &p[0]);

	// p10
	p[i] = oldPi + e;
	p[j] = oldPj - e;
	isSucceed &= computeValue(&s10[0], &p[0]);

	// p11
	p[i] = oldPi + e;
	p[j] = oldPj + e;
	isSucceed &= computeValue(&s11[0], &p[0]);

	// 中心差商法计算
	double invE2 = 1 / (4.0 * e * e);
	for (int d = 0; d < nDoF; ++d)
	{
		res[d] = ((s00[d] + s11[d]) - (s10[d] + s01[d])) * invE2;
	}

	// 恢复原有参数
	p[i] = oldPi;
	p[j] = oldPj;
	setFreeParam(&p[0]);
	return isSucceed;
}

void RigFEM::RigBase::keyParam( int ithParam, KeyFrameFunc func )
{
	m_keyFrameFunc[ithParam] = func;
}

void RigFEM::RigBase::unKeyParam( int ithParam )
{
	m_keyFrameFunc[ithParam] = NULL;
}

int RigFEM::RigBase::getNFreeParam() const
{
	int nFree = 0;
	for (int i = 0; i < m_nParam; ++i)
		nFree += (m_keyFrameFunc[i] == NULL);
	return nFree;
}

void RigFEM::RigBase::setFreeParam( const double* params )
{
	for (int ithParam = 0, ithFree = 0; ithParam < m_nParam; ++ithParam)
	{
		if (!m_keyFrameFunc[ithParam])
		{
			m_param[ithParam] = params[ithFree++];
		}
	}
}

void RigFEM::RigBase::getFreeParam( double* params ) const
{
	for (int ithParam = 0, ithFree = 0; ithParam < m_nParam; ++ithParam)
	{
		if (!m_keyFrameFunc[ithParam])
		{
			params[ithFree++] = m_param[ithParam];
		}
	}
}

void RigFEM::RigBase::setTime( double t )
{
	m_t = t;
	// 根据给定的关键帧函数计算参数值
	for (int i = 0; i < m_nParam; ++i)
	{
		if (m_keyFrameFunc[i])
		{
			m_param[i] = m_keyFrameFunc[i](t);
		}
	}
}

bool RigFEM::RigBase::computeJacobian( Eigen::MatrixXd& jacobian )
{
	double e = 1e-4;
	double inv2E = 0.5 / e;

	int nParam = getNFreeParam();
	int nDof   = getNPoints() * 3;

	jacobian = EigDense(nDof, nParam);
	EigVec p(nParam);
	getFreeParam(&p[0]);

	EigVec s0(nDof), s1(nDof);
	double* pJ = jacobian.data();
	bool res = true;
	for (int ithParam = 0; res && ithParam < nParam; ++ithParam)
	{
		double oldVal = p[ithParam];

		p[ithParam] = oldVal - e;
		res &= computeValue(&s0[0], &p[0]);

		p[ithParam] = oldVal + e;
		res &= computeValue(&s1[0], &p[0]);

		for (int ithDof = 0; ithDof < nDof; ++ithDof)
		{
			*pJ = (s1[ithDof] - s0[ithDof]) * inv2E;
			pJ++;
		}

		p[ithParam] = oldVal;
	}

	setFreeParam(&p[0]);
	return res;
}

void RigFEM::RigBase::getParam( double* params )
{
	for (int i = 0; i < m_nParam; ++i)
	{
		params[i] = m_param[i];
	}
}
