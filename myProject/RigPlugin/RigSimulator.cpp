#include "StdAfx.h"
#include "RigSimulator.h"

using namespace RigFEM;


RigSimulator::~RigSimulator(void)
{
	freeSimObj();
	delete m_rig;
	m_rig = NULL;
}

void RigFEM::RigSimulator::freeSimObj()
{
	delete m_rigMesh;
	delete m_solver;
	delete m_recorder;
	m_rigMesh = NULL;
	m_solver = NULL;
	m_recorder = NULL;
}

void RigFEM::RigSimulator::allocateSimObj()
{
	freeSimObj();
	m_rigMesh = new RiggedMesh();
	m_solver  = new PointParamSolver(m_rigMesh);
	m_recorder= new StatusRecorder();
}

bool RigFEM::RigSimulator::testHessian( int curFrame )
{
	if (!m_recorder || !m_rigMesh || !m_solver)
		return false;

	// 设置当前帧状态
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame-1, lastStatus);
	if (!res)
		res = getInitStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame, curStatus);
	if (!res)
		res = getInitStatus(curStatus);
	if (!res)
		return false;
	double noiseN = 1.0, noiseP = 1.0;
	m_rigMesh->testCurFrameHessian(lastStatus, curStatus, noiseN, noiseP);
	return true;
}

bool RigFEM::RigSimulator::testGradient( int curFrame )
{
	if (!m_recorder || !m_rigMesh || !m_solver)
		return false;

	// 设置当前帧状态
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame-1, lastStatus);
	if (!res)
		res = getInitStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame, curStatus);
	if (!res)
		res = getInitStatus(curStatus);
	if (!res)
		return false;

	double noiseN = 1.0, noiseP = 1.0;
	m_rigMesh->testCurFrameGrad(lastStatus, curStatus, noiseN, noiseP);
	return true;
}

bool RigFEM::RigSimulator::saveResult( const char* fileName )
{
	if (!m_recorder)
		return false;
	return m_recorder->saveToFile(fileName);
}

bool RigFEM::RigSimulator::stepRig( int curFrame )
{
	if (!m_recorder || !m_rigMesh || !m_solver)
		return false;

	// 设置上一帧状态
	RigStatus s;
	bool res = m_recorder->getStatus(curFrame-1, s);
	if (!res)
		res = getInitStatus(s);
	if (res)
	{
		m_solver->setInitStatus(s);
	}

	if(!m_rig)
		return false; 

	// 模拟
	if(!m_solver->step())
	{
		PRINT_F("simulation failed.");
		return false;
	}

	// 记录结果
	s = m_solver->getFinalStatus();
	return m_recorder->setStatus(curFrame, s);
}

bool RigFEM::RigSimulator::getInitStatus( RigStatus& status )
{
	int paramLength = m_recorder->getParamVecLength();
	int pntLength = m_recorder->getPointVecLength();

	if (paramLength <= 0 || pntLength <= 0)
		return false;

	EigVec q,v,a;
	q.setZero(pntLength);
	v.setZero(pntLength);
	a.setZero(pntLength);
	status = RigStatus(q,v,a,m_initParam);
	return true;
}

bool RigSimulator::init( tetgenio& surfMesh, RigSimulationNode* node, EigVec initParam, int curFrame /*= 0*/, double maxVolume /*= 1*/, double edgeRatio /*= 2 */, double youngModulus /*= 1e6*/, double nu /*= 0.45*/, double density /*= 1000*/, double timeStep /*= 1/24.0*/ , int maxStaticSolveIter)
{
	m_maxStaticSolveIter = maxStaticSolveIter;

	int nValidParam = initParam.size();
	m_initParam = initParam;
	int nSurfVtx = surfMesh.numberofpoints;
	GeneralRig* rig = new GeneralRig(node, nValidParam, nSurfVtx);
	m_rig = rig;
	allocateSimObj();

	bool res = m_rigMesh->init(surfMesh, rig, maxVolume, edgeRatio, youngModulus, nu, density);
	m_rigMesh->setStepTime(timeStep);

	// 初始化结果记录器
	int totPnt = m_rigMesh->getNTotPnt();
	vector<double> pnts;
	m_rigMesh->getMeshPntPos(pnts);
	m_recorder->init(
		curFrame, totPnt*3, nValidParam, 
		m_rigMesh->getInternalPntIdx(), m_rigMesh->getSurfacePntIdx(),
		pnts);
	return res;
}

bool RigFEM::RigSimulator::setStaticSolveMaxIter( int maxIter )
{
	m_maxStaticSolveIter = maxIter;
	return true;
}

bool RigFEM::RigSimulator::staticStepRig( int curFrame, const EigVec& curParam )
{
	if (!m_rig || !m_rigMesh || !m_solver || !m_recorder)
		return false;
	if (curParam.size() != m_recorder->getParamVecLength())
		return false;
	// 获得上一帧状态
	RigStatus s;
	bool res = m_recorder->getStatus(curFrame-1, s);
	if (!res)
		res =getInitStatus(s);
	if (!res)
		return false;

	// 求解平衡位置
	EigVec p = curParam;
	EigVec q;
	const EigVec* initQ = &s.getQ();		// 以上一帧的位置开始迭代，加快速度
	if(!m_rigMesh->computeStaticPos(p, 0, q, m_maxStaticSolveIter, initQ))
		return false;

	// 记录结果
	double dt = m_rigMesh->getStepTime();
	EigVec v = (q - s.getQ()) / dt;
	EigVec a = (v - s.getV()) / dt;
	return m_recorder->setStatus(curFrame, RigStatus(q,v,a,p));
}

bool RigFEM::RigSimulator::setDeriStepSize( double step )
{
	m_rig->setDelta(step);
	return true;
}

bool RigFEM::RigSimulator::showStatus( int curFrame, double* bbox /*= 0*/ )
{
	RigStatus status;
	bool res = m_recorder->getStatus(curFrame, status);
	if (!res)
		res = getInitStatus(status);
	if (!res)
		return false;
	m_rigMesh->showStatus(status, bbox);
	return true;
}

bool RigFEM::RigSimulator::getParam( int curFrame, EigVec& curParam )
{
	RigStatus status;
	bool res = m_recorder->getStatus(curFrame, status);
	if (!res)
		res = getInitStatus(status);
	if (!res)
		return false;
	curParam = status.getP();
	return true;
}

bool RigFEM::RigSimulator::isReady() const
{
	return m_recorder != NULL && m_rig != NULL && m_rigMesh != NULL && m_solver != NULL;
}

void RigFEM::RigSkinSimulator::allocateSimObj()
{
	RiggedSkinMesh* skinMesh = new RiggedSkinMesh();
	m_solver  = new ParamSolver(skinMesh);
	m_recorder= new StatusRecorder();
	m_rigMesh = skinMesh;
}

bool RigFEM::RigSkinSimulator::testHessian( int curFrame )
{
	if (!isReady())
		return false;

	// 设置当前帧状态
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame-1, lastStatus);
	if (!res)
		res = getInitStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame, curStatus);
	if (!res)
		res = getInitStatus(curStatus);
	if (!res)
		return false;
	double noiseP = 1.0;
	RiggedSkinMesh* mesh = (RiggedSkinMesh*) m_rigMesh;
	mesh->testCurFrameHessian(lastStatus, curStatus, noiseP);
	return true;
}

bool RigFEM::RigSkinSimulator::testGradient( int curFrame )
{
	if (!isReady())
		return false;

	// 设置当前帧状态
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame-1, lastStatus);
	if (!res)
		res = getInitStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame, curStatus);
	if (!res)
		res = getInitStatus(curStatus);
	if (!res)
		return false;
	double noiseP = 1.0;
	RiggedSkinMesh* mesh = (RiggedSkinMesh*)m_rigMesh;
	mesh->testCurFrameGrad(lastStatus, curStatus, noiseP);
	return true;
}

bool RigFEM::RigSkinSimulator::staticStepRig( int curFrame, const EigVec& curParam )
{
	PRINT_F("Static solve is not supported currently.");
	return false;
}

bool RigFEM::RigSkinSimulator::init( tetgenio& surfMesh, RigSimulationNode* node, EigVec initParam, const char* weightFile, int curFrame /*= 0*/, double maxVolume /*= 1*/, double edgeRatio /*= 2 */, double youngModulus /*= 1e6*/, double nu /*= 0.45*/, double density /*= 1000*/, double timeStep /*= 1/24.0*/, int maxStaticSolveIter/*=20*/ )
{
	bool res = RigSimulator::init(surfMesh, node, initParam, curFrame, maxVolume, edgeRatio, youngModulus, nu, density, timeStep, maxStaticSolveIter);
	if (!res)
		return false;
	RiggedSkinMesh* mesh = (RiggedSkinMesh*)m_rigMesh;
	res = mesh->setWeight(weightFile);
	if (!res)
		return false;
	return true;
}
