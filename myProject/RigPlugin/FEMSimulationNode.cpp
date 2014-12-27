#include "StdAfx.h"
#include "FEMSimulationNode.h"

using namespace RigFEM;

MTypeId RigSimulationNode::m_id(NODE_FEM_SIMULATION_ID);
const char* RigSimulationNode::m_nodeName = NODE_FEM_SIMULATION_NAME;

const char* RigSimulationNode::m_paramName[2] = {"rigParameter","rParam"};
const char* RigSimulationNode::m_meshName[2] = {"riggedMesh", "rMesh"};
const char* RigSimulationNode::m_transformMatrixName[2] = {"meshTransform", "trans"};

MObject RigSimulationNode::m_transformMatrix;
MObject RigSimulationNode::m_param; 
MObject RigSimulationNode::m_mesh;

RigSimulationNode::RigSimulationNode(void):m_box(MPoint(-0.1,-0.5,-0.1), MPoint(1.1,0.5,1.1)),m_rigMesh(NULL), m_rig(NULL), m_solver(NULL)
{
}

RigSimulationNode::~RigSimulationNode(void)
{
	clearRig();
}

void RigSimulationNode::postConstructor()
{
	MFnDependencyNode nodeFn(thisMObject());
	nodeFn.setName( "rigSimulationShape#");
}

void RigSimulationNode::draw( M3dView & view, const MDagPath & path, M3dView::DisplayStyle style, M3dView:: DisplayStatus )
{
	view.beginGL();
	glPushAttrib(GL_CURRENT_BIT);

	glBegin(GL_LINE_LOOP);
	glVertex3f(0,0,0);
	glVertex3f(1,0,0);
	glVertex3f(1,0,1);
	glVertex3f(0,0,1);
	glEnd();

	if (m_rigMesh)
	{
		m_rigMesh->show();
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

	m_param = nAttr.create(m_paramName[0], m_paramName[1], MFnNumericData::kDouble,0, &s);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true); 
	nAttr.setReadable(true);
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true); 
	nAttr.setUsesArrayDataBuilder(true);


	m_mesh = tAttr.create(m_meshName[0], m_meshName[1], MFnData::kMesh, &s);
	tAttr.setStorable(false);
	tAttr.setHidden(false);
	tAttr.setWritable(true);
	tAttr.setReadable(false);
	tAttr.setArray(false);
	tAttr.setKeyable(true);


	m_transformMatrix = mAttr.create(m_transformMatrixName[0], m_transformMatrixName[1]);
	mAttr.setHidden(false);
	mAttr.setReadable(false);
	mAttr.setWritable(true);
	mAttr.setKeyable(true);

	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_param);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_mesh);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_transformMatrix);

	return s;
}

MStatus RigSimulationNode::compute( const MPlug& plug, MDataBlock& data )
{
	MStatus s;
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
		char str[50];
		sprintf(str, "%lf \t %lf \t %lf\n", p.x, p.y, p.z);
		MGlobal::displayInfo(str);
	}
}

int RigSimulationNode::getNumParam()
{
	MStatus s;
	MPlug paramArrayPlug = Global::getPlug(this, m_paramName[0]);
	int maxLogicalIdx = -1;
	for (int phyIdx = 0; phyIdx < paramArrayPlug.numElements(&s); ++phyIdx)
	{
		MPlug paramPlug = paramArrayPlug.elementByPhysicalIndex(phyIdx,&s);
		if (!s)
			return 0;
		int logIdx = paramPlug.logicalIndex();
		maxLogicalIdx = max(maxLogicalIdx, logIdx);
	}
	return maxLogicalIdx+1;
}

void RigSimulationNode::clearRig()
{
	delete m_rig;
	m_rig = NULL;
	delete m_rigMesh;
	m_rigMesh = NULL;
	delete m_solver;
	m_solver = NULL;
}

bool RigSimulationNode::resetRig()
{
	clearRig();

	int nValidParam = getNumParam();
	if (nValidParam <= 0)
		return false;

	// 获取网格对象
	MStatus s;
	MPlug meshPlug = Global::getPlug(this, m_meshName[0]);
	MObject meshObj = meshPlug.asMObject(MDGContext::fsNormal, &s);
	if (s != MS::kSuccess)
		return false;
	MFnMesh meshFn(meshObj);
	int nVtx = meshFn.numVertices();
	if (nVtx <= 0 || meshFn.numPolygons() <= 0)
		return false;

	// 初始化rig对象
	m_rig = new RigFEM::GeneralRig(this, nValidParam, nVtx);
	m_rig->fetchParamFromNode();

	// 初始化RiggedMesh对象
	tetgenio tetMesh;
	MBoundingBox box = meshFn.boundingBox();
	MMatrix mat = getMatrix();
	Global::maya2TetgenMesh(meshObj, tetMesh, mat);
	m_rigMesh = new RigFEM::RiggedMesh();
	bool res = m_rigMesh->init(tetMesh, m_rig);

	if (!res)
	{
		clearRig();
		return false;
	}

	m_solver = new RigFEM::NewtonSolver(m_rigMesh);
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

bool RigSimulationNode::stepRig()
{
	if (!m_rig || !m_rigMesh || !m_solver)
		return false;

	bool res = m_solver->step();
	return res;
}





