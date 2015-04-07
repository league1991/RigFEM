#include "StdAfx.h"
#include "RigSimulator.h"

using namespace RigFEM;

const char* RigFEM::SimulatorBase::s_materialName = "material";


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
	m_rigMesh = NULL;
	m_solver = NULL;
	m_recorder = NULL;
}

void RigFEM::RigSimulator::allocateSimObj()
{
	freeSimObj();
	m_rigMesh = new RiggedMesh();
	m_solver  = new PointParamSolver(m_rigMesh);
	m_recorder= &m_solver->getRecorder();
}

bool RigFEM::RigSimulator::testHessian( int curFrame, double noiseN, double noiseP )
{
	if (!m_recorder || !m_rigMesh || !m_solver)
		return false;

	// 设置当前帧状态
	m_solver->setCurrentFrame(curFrame);
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame, lastStatus);
	if (!res)
		res = m_solver->getRestStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame+1, curStatus);
	if (!res)
		res = m_solver->getRestStatus(curStatus);
	if (!res)
		return false;
	m_rigMesh->testCurFrameHessian(lastStatus, curStatus, noiseN, noiseP);
	return true;
}

bool RigFEM::RigSimulator::testGradient( int curFrame, double noiseN, double noiseP )
{
	if (!m_recorder || !m_rigMesh || !m_solver)
		return false;

	// 设置当前帧状态
	m_solver->setCurrentFrame(curFrame);
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame, lastStatus);
	if (!res)
		res = m_solver->getRestStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame+1, curStatus);
	if (!res)
		res = m_solver->getRestStatus(curStatus);
	if (!res)
		return false;

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
	m_solver->setCurrentFrame(curFrame);
	RigStatus s;
	bool res = m_recorder->getStatus(curFrame, s);
	if (!res)
		res = m_solver->getRestStatus(s);
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
	return m_recorder->setStatus(curFrame+1, s);
}


bool RigSimulator::init( tetgenio& surfMesh, RigSimulationNode* node, EigVec initParam, int curFrame /*= 0*/, double maxVolume /*= 1*/, double edgeRatio /*= 2 */, double youngModulus /*= 1e6*/, double nu /*= 0.45*/, double density /*= 1000*/, double timeStep /*= 1/24.0*/ , int maxStaticSolveIter)
{

	int nValidParam = initParam.size();
	int nSurfVtx = surfMesh.numberofpoints;
	GeneralRig* rig = new GeneralRig(node, nValidParam, nSurfVtx);
	m_rig = rig;
	allocateSimObj();

	bool res = m_rigMesh->init(surfMesh, rig, maxVolume, edgeRatio, youngModulus, nu, density);
	m_rigMesh->setStepTime(timeStep);

	m_solver->setRestParam(initParam);
	m_solver->setStaticSolveMaxIter(maxStaticSolveIter);

	// 初始化结果记录器
	int totPnt = m_rigMesh->getNTotPnt();
	vector<double> pnts;
	m_rigMesh->getMeshPntPos(pnts);
	m_recorder->init(
		curFrame, totPnt*3, nValidParam, 
		m_rigMesh->getInternalPntIdx(), m_rigMesh->getSurfacePntIdx(),
		pnts);

	RigStatus initStatus;
	m_solver->getRestStatus(initStatus);
	m_recorder->setStatus(curFrame, initStatus);
	return res;
}

bool RigFEM::RigSimulator::setStaticSolveMaxIter( int maxIter )
{
	if (!m_solver)
	{
		return false;
	}
	m_solver->setStaticSolveMaxIter(  maxIter );
	return true;
}

bool RigFEM::RigSimulator::staticStepRig( int curFrame, const EigVec& curParam )
{
	if (!m_rig || !m_rigMesh || !m_solver || !m_recorder)
		return false;
	if (curParam.size() != m_recorder->getParamVecLength())
		return false;
	m_solver->setCurrentFrame(curFrame);

	// 获得上一帧状态
	return m_solver->staticSolve(curParam);
}

bool RigFEM::RigSimulator::setDeriStepSize( double step )
{
	m_rig->setDelta(step);
	return true;
}

bool RigFEM::RigSimulator::showStatus( int curFrame, double* bbox)
{
	RigStatus status;
	bool res = m_recorder->getStatus(curFrame, status);
	if (!res)
		res = m_solver->getRestStatus(status);
	if (!res)
		return false;
	m_rigMesh->showStatus(status, m_dispConfig, bbox);
	return true;
}

bool RigFEM::RigSimulator::getParam( int curFrame, EigVec& curParam )
{
	RigStatus status;
	bool res = m_recorder->getStatus(curFrame, status);
	if (!res)
		res = m_solver->getRestStatus(status);
	if (!res)
		return false;
	curParam = status.getP();
	return true;
}

bool RigFEM::RigSimulator::isReady() const
{
	return m_recorder != NULL && m_rig != NULL && m_rigMesh != NULL && m_solver != NULL;
}

bool RigFEM::RigSimulator::setControlType( RigControlType type )
{
	if (!m_solver)
	{
		return false;
	}
	m_solver->setControlType(type);
	return true;
}

bool RigFEM::RigSimulator::stepAndSaveEleGFRig( int curFrame, const EigVec& curParam )
{
	if (!m_rig || !m_rigMesh || !m_solver || !m_recorder)
		return false;
	if (curParam.size() != m_recorder->getParamVecLength())
		return false;
	m_solver->setCurrentFrame(curFrame);

	// 获得上一帧状态
	return m_solver->staticSolveWithEleGF(curParam);
}

bool RigFEM::RigSimulator::saveGFResult( const char* fileName )
{
	if (!isReady())
		return false;

	ofstream file(fileName);
	if (!file)
		return false;

	m_recorder->addCustomToFile(m_solver->s_reducedElementGF, file);

	char buf[50];
	int nParam = m_rigMesh->getNParam();
	sprintf(buf, "nParam = %d;\n", nParam);
	file << buf;
	file.close();

	return true;
}

void RigFEM::RigSkinSimulator::allocateSimObj()
{
	RiggedSkinMesh* skinMesh = new RiggedSkinMesh();
	m_solver  = new ParamSolver(skinMesh);
	m_recorder= &m_solver->getRecorder();
	m_rigMesh = skinMesh;
}

bool RigFEM::RigSkinSimulator::testHessian( int curFrame )
{
	if (!isReady())
		return false;

	// 设置当前帧状态
	m_solver->setCurrentFrame(curFrame);
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame, lastStatus);
	if (!res)
		res = m_solver->getRestStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame+1, curStatus);
	if (!res)
		res = m_solver->getRestStatus(curStatus);
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
	m_solver->setCurrentFrame(curFrame);
	RigStatus lastStatus, curStatus;
	bool res = m_recorder->getStatus(curFrame, lastStatus);
	if (!res)
		res = m_solver->getRestStatus(lastStatus);
	res&= m_recorder->getStatus(curFrame+1, curStatus);
	if (!res)
		res = m_solver->getRestStatus(curStatus);
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

void RigFEM::SimulatorBase::setExternalForceDispFactor( double factor )
{
	m_dispConfig.m_extForceDispFactor = factor;
}

bool RigFEM::SimulatorBase::loadElementMaterialFactor( const char* fileName )
{
	if (!isReady())
	{
		return false;
	}
	std::map<string, EigDense> matMap;
	Utilities::fileToDense(fileName, matMap);
	if (matMap.find(s_materialName) == matMap.end())
		return false;

	EigDense& mat = matMap[s_materialName];
	EigVec vec;
	Utilities::eigDense2Vec(mat, vec);
	m_rigMesh->loadElementMaterialFactor(vec);
	return true;
}

bool RigFEM::SimulatorBase::resetElementMaterialFactor()
{
	if (!isReady())
	{
		return false;
	}
	m_rigMesh->clearElementMaterialFactor();
	return true;
}
