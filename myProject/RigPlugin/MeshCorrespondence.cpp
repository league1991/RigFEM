#include "StdAfx.h"
#include "MeshCorrespondence.h"

const char* CageDeformerNode::m_outMeshName[2] = {"outMesh", "outMesh"};
const char* CageDeformerNode::m_inControlMeshName[2] = {"inControlMesh", "inCtrlMesh"};
const char* CageDeformerNode::m_refMeshName[2] = {"refMesh", "refMesh"};
const char* CageDeformerNode::m_refControlMeshName[2]={"refControlMesh","refCtrlMesh"};
const char* CageDeformerNode::m_maxNeighName[2]= {"maxNeighbourNum", "maxN"};
const char* CageDeformerNode::m_powerName[2] = {"weightPower","pow"};

MObject CageDeformerNode::m_power;
MObject CageDeformerNode::m_maxNeigh; 
MObject CageDeformerNode::m_outMesh;
MObject CageDeformerNode::m_inControlMesh;
MObject CageDeformerNode::m_refMesh;
MObject CageDeformerNode::m_refControlMesh;

const char* CageDeformerNode::m_nodeName = NODE_CAGE_DEFORMER_NAME;

const MTypeId CageDeformerNode::m_id(NODE_CAGE_DEFORMER_ID);

CageDeformerNode::CageDeformerNode(void)
{
}

CageDeformerNode::~CageDeformerNode(void)
{
}

MStatus CageDeformerNode::compute( const MPlug& plug, MDataBlock& data )
{
	if (plug == m_outMesh)
	{
		computeNewPnt();
		data.setClean(plug);
	}
	return MS::kSuccess;
}

MStatus CageDeformerNode::initialize()
{
	MStatus s; 
	MFnTypedAttribute tAttr;
	MFnNumericAttribute nAttr;

	m_refControlMesh = tAttr.create(m_refControlMeshName[0], m_refControlMeshName[1], MFnData::kMesh, &s);
	tAttr.setStorable(true);
	tAttr.setHidden(false);
	tAttr.setWritable(true);
	tAttr.setReadable(false);
	tAttr.setArray(false);

	m_refMesh = tAttr.create(m_refMeshName[0], m_refMeshName[1], MFnData::kMesh, &s);
	tAttr.setStorable(true);
	tAttr.setHidden(false);
	tAttr.setWritable(true);
	tAttr.setReadable(false);
	tAttr.setArray(false);

	m_inControlMesh = tAttr.create(m_inControlMeshName[0], m_inControlMeshName[1], MFnData::kMesh, &s);
	tAttr.setStorable(false);
	tAttr.setHidden(false);
	tAttr.setWritable(true);
	tAttr.setReadable(false);
	tAttr.setArray(false);

	m_outMesh = tAttr.create(m_outMeshName[0], m_outMeshName[1], MFnData::kMesh, &s);
	tAttr.setStorable(false);
	tAttr.setHidden(false);
	tAttr.setWritable(false);
	tAttr.setReadable(true);
	tAttr.setArray(false);
	tAttr.setAffectsAppearance(true);

	m_maxNeigh = nAttr.create(m_maxNeighName[0], m_maxNeighName[1], MFnNumericData::kInt,3);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true);
	nAttr.setMin(1);
	nAttr.setMax(10000);
	nAttr.setSoftMin(3);
	nAttr.setSoftMax(10);

	m_power = nAttr.create(m_powerName[0], m_powerName[1], MFnNumericData::kDouble,1);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true);
	nAttr.setMin(1e-3);
	nAttr.setMax(1e6);
	nAttr.setSoftMin(1);
	nAttr.setSoftMax(5);

	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute( m_refControlMesh );
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute( m_refMesh );
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute( m_inControlMesh );
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute( m_outMesh );
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute( m_maxNeigh );
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute( m_power );

	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_inControlMesh, m_outMesh);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_refControlMesh, m_outMesh);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_refMesh, m_outMesh);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_power, m_outMesh);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_maxNeigh, m_outMesh);

	return MS::kSuccess;
}

MStatus CageDeformerNode::getVertices( MObject& meshObj, EigDense& vertMat ) const
{
	MStatus s;
	MFnMesh meshFn(meshObj, &s);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	const int nVertices = meshFn.numVertices(&s);
	if (nVertices <= 0)
		return MS::kFailure;

	// fill vertex array	
	vertMat.resize(nVertices, 3);
	MItMeshVertex it(meshObj, &s);
	for (int ithVert=0; !it.isDone(&s); it.next(),ithVert++)
	{
		MPoint p;
		meshFn.getPoint(ithVert, p);
		vertMat(ithVert, 0) = p.x;
		vertMat(ithVert, 1) = p.y;
		vertMat(ithVert, 2) = p.z;
	}
	return s;
}

MStatus CageDeformerNode::computeWeight()
{
	MStatus s;
	MPlug refControlMeshPlug = Global::getPlug(this, m_refControlMeshName[0]);
	MPlug refMeshPlug		 = Global::getPlug(this, m_refMeshName[0]);
	MPlug nNeighPlug         = Global::getPlug(this, m_maxNeighName[0]);
	MPlug powerPlug			 = Global::getPlug(this, m_powerName[0]);

	MObject refControlMeshObj = refControlMeshPlug.asMObject();
	MObject refMeshObj        = refMeshPlug.asMObject();

	EigDense controlPnt, meshPnt;
	s = getVertices(refControlMeshObj, controlPnt);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = getVertices(refMeshObj, meshPnt);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	int nCtrlMeshPnts = controlPnt.rows();
	int nMeshPnts = meshPnt.rows();

	nanoflann::KDTreeEigenMatrixAdaptor<EigDense> kdtree(3, controlPnt);

	const int nNearest = nNeighPlug.asInt();
	double exp = powerPlug.asDouble();
	typedef Eigen::Triplet<double> Tri;
	vector<Tri> triList;
	vector<double> weight(nNearest);
	vector<bool>   isAccepted(nNearest);
	vector<double> dist(nNearest);
	vector<EigDense::Index> resultIdx(nNearest);
	for (int ithMeshPnt = 0; ithMeshPnt < nMeshPnts; ++ithMeshPnt)
	{
		// 找近邻点
		double pnt[3] = {meshPnt(ithMeshPnt,0), meshPnt(ithMeshPnt,1), meshPnt(ithMeshPnt,2)};
		int nFound = kdtree.query(pnt, nNearest, &resultIdx[0], &dist[0]);

		if (!nFound)
			continue;

		// 计算权重
		double maxWeight = 0;
		//PRINT_F("neighbour of %lf %lf %lf", pnt[0],pnt[1],pnt[2]);
		for (int ithFound = 0; ithFound < nFound; ++ithFound)
		{
			double d = dist[ithFound];
			//PRINT_F("%d :%lf", resultIdx[ithFound], d);
			weight[ithFound] = pow(1.0 / (d + 1e-10), exp);
			maxWeight = max(maxWeight, weight[ithFound]);
		}
		// 筛选权重
		double sumWeight = 0;
		for (int ithFound = 0; ithFound < nFound; ++ithFound)
		{
			isAccepted[ithFound] = weight[ithFound] > maxWeight  * 0.01;
			if (isAccepted[ithFound])
			{
				sumWeight += weight[ithFound];
			}
		}

		// 添加最终权重
		for (int ithFound = 0; ithFound < nFound; ++ithFound)
		{
			if (isAccepted[ithFound])
			{
				double w = weight[ithFound] / sumWeight;
				triList.push_back(Tri(ithMeshPnt, resultIdx[ithFound], w));
			}
		}
	}

	m_weightMat = EigSparse(nMeshPnts, nCtrlMeshPnts);
	m_weightMat.setFromTriplets(triList.begin(), triList.end());
	PRINT_F("weight matrix: %d x %d, %d weights", nMeshPnts, nCtrlMeshPnts, triList.size());
	return MS::kSuccess;
}

MStatus CageDeformerNode::computeNewPnt()
{
	MStatus s;
	MPlug refControlMeshPlug = Global::getPlug(this, m_refControlMeshName[0]);
	MPlug refMeshPlug		 = Global::getPlug(this, m_refMeshName[0]);
	MPlug inControlMeshPlug    = Global::getPlug(this, m_inControlMeshName[0]);
	MPlug outMeshPlug           = Global::getPlug(this, m_outMeshName[0]);

	MObject refControlMeshObj = refControlMeshPlug.asMObject();
	MObject refMeshObj        = refMeshPlug.asMObject();
	MObject inControlMeshObj    = inControlMeshPlug.asMObject();

	EigDense refControlMeshPnt, refMeshPnt, inControlMeshPnt, outMeshPnt;
	s = getVertices(refControlMeshObj, refControlMeshPnt);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = getVertices(refMeshObj, refMeshPnt);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = getVertices(inControlMeshObj, inControlMeshPnt);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	MFnMeshData meshData;
	MObject     outMeshObj = meshData.create();
	MFnMesh     outMeshFn;
	outMeshFn.copy(refMeshObj, outMeshObj,&s);
	outMeshFn.setObject(outMeshObj);

	if (refControlMeshPnt.rows() != inControlMeshPnt.rows() ||			
		m_weightMat.rows() != refMeshPnt.rows() ||
		m_weightMat.cols() != refControlMeshPnt.rows())
	{
		MGlobal::displayError("Mesh resolution doesn't match.");
	}
	else
	{
		EigDense outPnt = refMeshPnt + m_weightMat * (inControlMeshPnt - refControlMeshPnt);
		MItMeshVertex it(outMeshObj, &s);
		for (int ithVtx = 0; !it.isDone(&s); it.next(), ++ithVtx)
		{
			MVector pos(outPnt(ithVtx,0), outPnt(ithVtx,1), outPnt(ithVtx, 2));
			outMeshFn.setPoint(it.index(),pos);
		}
		//MGlobal::displayInfo("Computation Complete.");
	}

	outMeshFn.updateSurface();
	outMeshPlug.setMObject(outMeshObj);
	return MS::kSuccess;
}



