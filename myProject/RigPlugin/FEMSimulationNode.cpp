#include "StdAfx.h"
#include "FEMSimulationNode.h"

using namespace RigFEM; 
 
MTypeId RigSimulationNode::m_id(NODE_FEM_SIMULATION_ID); 
const char* RigSimulationNode::m_nodeName = NODE_FEM_SIMULATION_NAME;

const char* RigSimulationNode::m_initParamName[2] = {"rigInitParameter","rInitParam"};
const char* RigSimulationNode::m_paramName[2] = {"rigParameter","rParam"};
const char* RigSimulationNode::m_meshName[2] = {"riggedMesh", "rMesh"};
const char* RigSimulationNode::m_transformMatrixName[2] = {"meshTransform", "trans"};
const char*	RigSimulationNode::m_tetEdgeRatioName[2]={"tetEdgeRatio", "edgeRatio"};			// 体网格化参数，四面体边比例
const char*	RigSimulationNode::m_tetMaxVolumeName[2]={"tetMaxVolume", "maxVolume"};			// 体网格化参数，四面体最大体积
const char*	RigSimulationNode::m_youngModulusName[2]={"youngModulus", "yModulus"};			// 杨氏模量
const char*	RigSimulationNode::m_nuName[2] = {"nu","nu"};					// 控制不同轴向变形影响程度的参数
const char*	RigSimulationNode::m_densityName[2]={"density","den"};				// 密度
const char*	RigSimulationNode::m_timeStepName[2]={"stepTime","dt"};				// 时间步长
const char* RigSimulationNode::m_deriStepName[2]={"derivativeStep", "dStep"};   // 导数步长
const char* RigSimulationNode::m_minGradSizeName[2] = {"inverseGradTolerance", "invGradSize"};
const char* RigSimulationNode::m_minStepSizeName[2] = {"inverseStepTolerance", "invStepSize"};
const char* RigSimulationNode::m_maxIterName[2] = {"maxIteration", "maxIter"};
const char* RigSimulationNode::m_simTypeName[2] = {"simulationType", "simType"};
const char* RigSimulationNode::m_weightPathName[2]={"weightPath","wPath"};

MObject RigSimulationNode::m_weightPath;
MObject RigSimulationNode::m_simType;
MObject RigSimulationNode::m_minGradSize;
MObject RigSimulationNode::m_minStepSize;
MObject RigSimulationNode::m_maxIter; 
MObject RigSimulationNode::m_deriStep;
MObject RigSimulationNode::m_transformMatrix;
MObject RigSimulationNode::m_param;   
MObject	RigSimulationNode::m_initParam;
MObject RigSimulationNode::m_mesh;
MObject	RigSimulationNode::m_tetEdgeRatio;			// 体网格化参数，四面体边比例
MObject	RigSimulationNode::m_tetMaxVolume;			// 体网格化参数，四面体最大体积
MObject	RigSimulationNode::m_youngModulus;			// 杨氏模量
MObject	RigSimulationNode::m_nu;					// 控制不同轴向变形相互影响的参数
MObject	RigSimulationNode::m_density;				// 密度
MObject	RigSimulationNode::m_timeStep;				// 时间步长

RigSimulationNode::RigSimulationNode(void):
m_box(MPoint(-1.1,-0.5,-1.1), MPoint(4.1,0.5,1.1)),
m_simTypeFlag(RigSimulationNode::SIM_STANDARD),
m_simulator(NULL)
{
}

RigSimulationNode::~RigSimulationNode(void)
{
	clearRig();
}

void RigSimulationNode::postConstructor()
{
	MStatus s;
	MFnDependencyNode nodeFn(thisMObject(), &s);
	if (s)
	{
		nodeFn.setName( "rigSimulationShape#", &s);
	}
}

void RigSimulationNode::draw( M3dView & view, const MDagPath & path, M3dView::DisplayStyle style, M3dView:: DisplayStatus )
{
	view.beginGL();
	glPushAttrib(GL_CURRENT_BIT);

	drawIcon();

	m_box = MBoundingBox(MPoint(-1.1,-0.5,-1.1), MPoint(4.1,0.5,1.1));
	if (m_simulator)
	{
		int curFrame = getCurFrame();

		MMatrix mat  = path.inclusiveMatrixInverse();
		double  matBuf[4][4];
		mat.get(matBuf);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glMultMatrixd(&matBuf[0][0]);

		// 更新参数值
		EigVec p;
		if(m_simulator->getParam(curFrame, p))
			setParamPlug(&p[0], p.size());
		else
			setParamToInit();
		double bbox[6];
		bool res = m_simulator->showStatus(curFrame, bbox);
		if (res)
		{
			// 更新包围盒
			double xformBox[6];
			Utilities::transformBBox(bbox, bbox+3, matBuf, xformBox, xformBox+3);
			m_box.expand(MPoint(xformBox[0], xformBox[1], xformBox[2]));
			m_box.expand(MPoint(xformBox[3], xformBox[4], xformBox[5]));
 
		}

		glPopMatrix();
	}
	glPopAttrib();
	view.endGL();
}

MBoundingBox RigSimulationNode::boundingBox() const
{
	return m_box;
}

MStatus RigSimulationNode::initialize()
{
	MStatus s;
	MFnNumericAttribute nAttr;
	MFnMatrixAttribute  mAttr;
	MFnTypedAttribute   tAttr;
	MFnEnumAttribute	eAttr;

	m_param = nAttr.create(m_paramName[0], m_paramName[1], MFnNumericData::kDouble,0, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(false); 
	nAttr.setReadable(true);
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true); 
	nAttr.setUsesArrayDataBuilder(true);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_initParam = nAttr.create(m_initParamName[0], m_initParamName[1], MFnNumericData::kDouble,0, &s);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true); 
	nAttr.setUsesArrayDataBuilder(true);
	nAttr.setAffectsAppearance(true);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_tetEdgeRatio = nAttr.create(m_tetEdgeRatioName[0], m_tetEdgeRatioName[1], MFnNumericData::kDouble,2.0, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(1.01);
	nAttr.setMax(10);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_tetMaxVolume = nAttr.create(m_tetMaxVolumeName[0], m_tetMaxVolumeName[1], MFnNumericData::kDouble,5.0, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(0.01);
	nAttr.setSoftMax(10);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_youngModulus = nAttr.create(m_youngModulusName[0], m_youngModulusName[1], MFnNumericData::kDouble,1e6, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(0.01);
	nAttr.setSoftMin(1e5);
	nAttr.setSoftMax(1e6);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_nu = nAttr.create(m_nuName[0], m_nuName[1], MFnNumericData::kDouble,0.45, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(0.01);
	nAttr.setMax(0.5);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_density = nAttr.create(m_densityName[0], m_densityName[1], MFnNumericData::kDouble,1000, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(0.01);
	nAttr.setSoftMin(500);
	nAttr.setSoftMax(5000);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_timeStep = nAttr.create(m_timeStepName[0], m_timeStepName[1], MFnNumericData::kDouble,1/24.0, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(0.001);
	nAttr.setSoftMin(0.01);
	nAttr.setSoftMax(1);
	CHECK_MSTATUS_AND_RETURN_IT(s);


	m_deriStep = nAttr.create(m_deriStepName[0], m_deriStepName[1], MFnNumericData::kDouble,1e-3, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(1e-5);
	nAttr.setMax(1);
	nAttr.setSoftMin(1e-4);
	nAttr.setSoftMax(1e-2);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_minGradSize = nAttr.create(m_minGradSizeName[0], m_minGradSizeName[1], MFnNumericData::kDouble,100, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(0.001);
	nAttr.setMax(1e15);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_minStepSize = nAttr.create(m_minStepSizeName[0], m_minStepSizeName[1], MFnNumericData::kDouble,100000, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(0.001);
	nAttr.setMax(1e7);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_maxIter = nAttr.create(m_maxIterName[0], m_maxIterName[1], MFnNumericData::kInt,10, &s);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setMin(1);
	nAttr.setMax(50);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_mesh = tAttr.create(m_meshName[0], m_meshName[1], MFnData::kMesh, &s);
	tAttr.setStorable(false);
	tAttr.setHidden(false);
	tAttr.setWritable(true);
	tAttr.setReadable(false);
	tAttr.setArray(false);
	tAttr.setKeyable(true);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_transformMatrix = mAttr.create(m_transformMatrixName[0], m_transformMatrixName[1]);
	mAttr.setHidden(false);
	mAttr.setReadable(false);
	mAttr.setWritable(true);
	mAttr.setKeyable(true);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_simType = eAttr.create(m_simTypeName[0], m_simTypeName[1]);
	eAttr.addField("Standard", RigSimulationNode::SIM_STANDARD);
	eAttr.addField("Skin",   RigSimulationNode::SIM_SKIN);
	eAttr.setHidden(false);
	eAttr.setReadable(true);
	eAttr.setWritable(true);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	m_weightPath = tAttr.create(m_weightPathName[0], m_weightPathName[1], MFnData::kString);
	tAttr.setKeyable(true);
	tAttr.setStorable(true);
	tAttr.setHidden(false);
	tAttr.setWritable(true);
	tAttr.setUsedAsFilename(true);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	s = addAttribute(m_param);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_initParam);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_mesh);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_transformMatrix);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_tetEdgeRatio);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_tetMaxVolume);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_youngModulus);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_nu);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_density);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_timeStep);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_deriStep);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_maxIter);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_minGradSize);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_minStepSize);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_simType);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_weightPath);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	s = attributeAffects(m_initParam, m_param);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	return s;
}

MStatus RigSimulationNode::setParamToInit()
{
	MStatus s;
	MPlug paramArrayPlug = Global::getPlug(this, m_paramName[0]);
	MPlug inParamArrayPlug= Global::getPlug(this, m_initParamName[0]);
	for (int phyIdx = 0; phyIdx < paramArrayPlug.numElements(&s); ++phyIdx)
	{
		MPlug paramPlug = paramArrayPlug.elementByPhysicalIndex(phyIdx,&s);
		CHECK_MSTATUS_AND_RETURN_IT(s);

		int logIdx = paramPlug.logicalIndex();
		MPlug inParamPlug = inParamArrayPlug.elementByLogicalIndex(logIdx, &s);
		CHECK_MSTATUS_AND_RETURN_IT(s);

		double val = 0;
		inParamPlug.getValue(val);
		paramPlug.setValue(val);
	}
	return s;
}
MStatus RigSimulationNode::compute( const MPlug& plug, MDataBlock& data )
{
	MStatus s;
	if (plug == m_param || plug.parent() == m_param)
	{
		// 设置当前状态的参数值
		int curFrame = getCurFrame();
		EigVec p;
		if(m_simulator && m_simulator->isReady() && m_simulator->getParam(curFrame, p))
			setParamPlug(&p[0], p.size());
		else
			setParamToInit();
	
	}
	data.setClean(plug);
	return s;
}

void* RigSimulationNode::creator()
{
	return new RigSimulationNode;
}

void RigSimulationNode::testRig()
{
	MStatus s;
	MPlug paramArrayPlug = Global::getPlug(this, m_paramName[0]);
	for (int ithPlug = 0; ithPlug < 10; ++ithPlug)
	{
		MPlug paramPlug = paramArrayPlug.elementByLogicalIndex(ithPlug, &s);
		int v = rand();
		paramPlug.setValue((double)v / 65535.0);
	}

	MPlug meshPlug = Global::getPlug(this, m_meshName[0]);
	MObject meshObj = meshPlug.asMObject();

	MFnMesh meshFn(meshObj);
	const int nVtx = meshFn.numVertices(&s);

	MItMeshVertex it(meshObj, &s);
	for (int ithVtx = 0; !it.isDone(&s) && ithVtx < 1; it.next(), ++ithVtx)
	{
		MPoint p = it.position(MSpace::kWorld);
		PRINT_F("%lf \t %lf \t %lf\n", p.x, p.y, p.z);
	}
} 

int RigSimulationNode::getNumParam()
{
	MStatus s;
	MPlug paramArrayPlug = Global::getPlug(this, m_paramName[0]);
	int nPlugs = paramArrayPlug.numElements(&s);
	for (int logIdx = 0; logIdx < nPlugs; ++logIdx)
	{
		MPlug paramPlug = paramArrayPlug.elementByLogicalIndex(logIdx,&s);
		if (!s || !paramPlug.isConnected(&s)) 
			return logIdx;
	}
	return nPlugs;
}

void RigSimulationNode::clearRig()
{
	delete m_simulator;
	m_simulator = NULL;
}

bool RigSimulationNode::resetRig()
{
	clearRig();

	int nValidParam = getNumParam();
	if (nValidParam <= 0)
	{
		PRINT_F("no param");
		return false;
	}

	// 获取网格对象
	MStatus s;
	setParamToInit();			// 设置参数值为初始参数指定的值
	MPlug meshPlug = Global::getPlug(this, m_meshName[0]);
	MObject meshObj = meshPlug.asMObject(MDGContext::fsNormal, &s);
	if (s != MS::kSuccess)
	{
		PRINT_F("failed to get mesh obj.");
		return false;
	}
	MFnMesh meshFn(meshObj);
	int nSurfVtx = meshFn.numVertices();
	if (nSurfVtx <= 0 || meshFn.numPolygons() <= 0)
	{
		PRINT_F("mesh is empty.");
		return false;
	}

	// 初始化RiggedMesh对象
	MPlug tetEdgeRatioPlug = Global::getPlug(this, m_tetEdgeRatioName[0]);
	MPlug tetMaxVolumePlug = Global::getPlug(this, m_tetMaxVolumeName[0]);
	MPlug youngModulusPlug = Global::getPlug(this, m_youngModulusName[0]);
	MPlug nuPlug = Global::getPlug(this, m_nuName[0]);
	MPlug densityPlug = Global::getPlug(this, m_densityName[0]);
	MPlug stepTimePlug = Global::getPlug(this, m_timeStepName[0]);
	MPlug weightPathPlug = Global::getPlug(this, m_weightPathName[0]);

	tetgenio tetMesh;
	MMatrix mat = getMatrix();
	bool res = MS::kSuccess == Global::maya2TetgenMesh(meshObj, tetMesh, mat);
	m_simTypeFlag = (SimulationType)Global::getPlug(this, m_simTypeName[0]).asShort();
	MString weightStr = weightPathPlug.asString();
	EigVec initParam;
	int curFrame = getCurFrame();
	getInitParam(initParam);
	if (res)
	{
		switch(m_simTypeFlag)
		{
		case RigSimulationNode::SIM_STANDARD:
			{
				RigFEM::RigSimulator* sim = new RigFEM::RigSimulator();
				res &= sim->init(
					tetMesh, this,
					initParam,
					curFrame,
					tetMaxVolumePlug.asDouble(), 
					tetEdgeRatioPlug.asDouble(), 
					youngModulusPlug.asDouble(), 
					nuPlug.asDouble(), 
					densityPlug.asDouble(),
					stepTimePlug.asDouble()
					);
				m_simulator = sim;
			}
			break;
		case RigSimulationNode::SIM_SKIN:
			{
				RigFEM::RigSkinSimulator* sim = new RigFEM::RigSkinSimulator();
				res &= sim->init(
					tetMesh, this,
					initParam,
					weightStr.asChar(),
					curFrame,
					tetMaxVolumePlug.asDouble(), 
					tetEdgeRatioPlug.asDouble(), 
					youngModulusPlug.asDouble(), 
					nuPlug.asDouble(), 
					densityPlug.asDouble(),
					stepTimePlug.asDouble()
					);
				m_simulator = sim;
			}
			break;
		}
	}

	if (!res)
	{
		PRINT_F("failed to get rig mesh obj.");
		clearRig();
		return false;
	}

	// 初始化rig对象
	MPlug deriStepPlug = Global::getPlug(this, m_deriStepName[0]);
	double deriStepVal = 1e-3;
	deriStepPlug.getValue(deriStepVal);
	m_simulator->setDeriStepSize(deriStepVal);
	return true;
}

MMatrix RigSimulationNode::getMatrix()
{
	// 获取变换矩阵
	MPlug matrixPlug  = Global::getPlug(this, RigSimulationNode::transformLongName());
	MObject matrixObj = matrixPlug.asMObject();
	MFnMatrixData matrixData(matrixObj);
	return matrixData.matrix();
}

bool RigSimulationNode::testHessian(double noiseN, double noiseP)
{
	if (!m_simulator || !m_simulator->isReady())
		return false;

	updateDeriStepSize();

	// 设置当前帧状态
	int curFrame = getCurFrame();
	return m_simulator->testHessian(curFrame);
}
bool RigSimulationNode::updateDeriStepSize()
{
	if(!m_simulator || !m_simulator->isReady())
		return false;
	MPlug deriStepPlug = Global::getPlug(this, m_deriStepName[0]);
	double deriStepVal = 1e-3;
	deriStepPlug.getValue(deriStepVal);
	m_simulator->setDeriStepSize(deriStepVal);
	return true;
}
bool RigSimulationNode::stepRig()
{
	if (!m_simulator || !m_simulator->isReady())
		return false;

	int curFrame = getCurFrame();

	// 设置上一帧状态
	updateDeriStepSize();		// 设置导数步长
	updateTerminationCond();	// 设置终止条件

	// 模拟
	return m_simulator->stepRig(curFrame);
}


int RigSimulationNode::getCurFrame()
{
	int curFrame = 1;
	MGlobal::executeCommand("currentTime -q", curFrame);
	return curFrame;
// 	return 1;
// 	MTime t = MAnimControl::currentTime();
// 	return t.value();
}

MStatus RigSimulationNode::setParamPlug( const double* param, int nParam )
{
	MStatus s;
	MPlug paramArrayPlug = Global::getPlug(this, m_paramName[0]);
	for (int ithParam = 0; ithParam < nParam; ++ithParam)
	{
		MPlug paramPlug = paramArrayPlug.elementByLogicalIndex(ithParam,&s);
		if (!s)
			break;
		
		// 参数发生更改时才更新，避免死循环
		double v;
		paramPlug.getValue(v);
		if (v != param[ithParam])
			paramPlug.setValue(param[ithParam]);
	}
	return s;
}

void RigSimulationNode::drawIcon()
{
	float p[][3] = {{0.5,0,-1},{0.8,0,1},{-1,0,1},{0.2,0,0.4}};

	glBegin(GL_LINE_LOOP);
	glVertex3fv(p[0]);
	glVertex3fv(p[1]);
	glVertex3fv(p[2]);
	glEnd();

	glBegin(GL_LINES);
	glVertex3fv(p[0]);
	glVertex3fv(p[3]);

	glVertex3fv(p[1]);
	glVertex3fv(p[3]);

	glVertex3fv(p[2]);
	glVertex3fv(p[3]); 
	glEnd();

	// draw "FEM"
	glBegin(GL_LINE_STRIP);
	glVertex3f(1.8,0,0);
	glVertex3f(1,0,0);
	glVertex3f(1,0,1);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3f(2.8,0,0);
	glVertex3f(2,0,0);
	glVertex3f(2,0,1);
	glVertex3f(2.8,0,1);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex3f(3,0,1);
	glVertex3f(3,0,0);
	glVertex3f(3.4,0,0.7);
	glVertex3f(3.8,0,0);
	glVertex3f(3.8,0,1);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(1,0,0.5);
	glVertex3f(1.6,0,0.5);
	glVertex3f(2,0,0.5); 
	glVertex3f(2.6,0,0.5);
	glEnd();
}

bool RigSimulationNode::updateTerminationCond()
{
	if (!m_simulator || !m_simulator->isReady())
		return false;
	if (RigSimulator* sim = dynamic_cast<RigSimulator*>(m_simulator))
	{
		MPlug maxIterPlug = Global::getPlug(this, m_maxIterName[0]);
		MPlug minGradPlug = Global::getPlug(this, m_minGradSizeName[0]);
		MPlug minStepPlug = Global::getPlug(this, m_minStepSizeName[0]);

		int maxIter = maxIterPlug.asInt();
		double gradSize= 1.0 / minGradPlug.asDouble();
		double stepSize= 1.0 / minStepPlug.asDouble();
		NewtonSolver* solver = sim->getSolver();
		if (solver)
			solver->setTerminateCond(maxIter, stepSize, gradSize);
	}
	return true;
}

bool RigSimulationNode::testGrad( double noiseN, double noiseP )
{
	if (!m_simulator || !m_simulator->isReady())
		return false;
	updateDeriStepSize();

	// 设置当前帧状态
	int curFrame = getCurFrame();
	return m_simulator->testGradient(curFrame);
}

bool RigSimulationNode::saveSimulationData( const char* fileName )
{
	if (!m_simulator || !m_simulator->isReady())
		return false;
	return m_simulator->saveResult(fileName);
}

bool RigSimulationNode::getInitParam( EigVec& p )
{
	MStatus s;
	int paramLength = getNumParam();

	if (paramLength <= 0)
		return false;

	p.resize(paramLength);
	MPlug paramArrayPlug = Global::getPlug(this, m_initParamName[0]);
	for (int ithParam = 0; ithParam < paramLength; ++ithParam)
	{
		MPlug paramPlug = paramArrayPlug.elementByLogicalIndex(ithParam,&s);
		if (!s)
			return false;
		paramPlug.getValue(p[ithParam]);
	}
	return true;
}

bool RigSimulationNode::staticStepRig()
{
	if (!m_simulator || !m_simulator->isReady())
		return false;

	updateDeriStepSize();

	int curFrame = getCurFrame();
	EigVec initParam;
	getInitParam(initParam);
	m_simulator->staticStepRig(curFrame, initParam);
	return true;
}












