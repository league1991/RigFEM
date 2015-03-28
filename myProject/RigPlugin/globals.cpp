#include "stdafx.h"
#include "globals.h"

MPlug Global::getPlug( MPxNode* node, const char* longName )
{
	MObject nodeObj = node->thisMObject();
	MFnDependencyNode nodeFn(nodeObj);
	return nodeFn.findPlug(longName, true);
}

MColor Global::getColor( const MPlug& plug, MStatus* s /*= NULL*/ )
{
	MColor color;
	int nC = plug.numChildren(s);
	if (nC == 3)
	{
		float r = plug.child(0).asFloat();
		float g = plug.child(1).asFloat();
		float b = plug.child(2).asFloat();
		color = MColor(r,g,b);
	}
	else if (s)
		*s = MS::kFailure;
	return color;
}

void Global::displayError( const MPxNode* node, const MString& errorMsg )
{
	MGlobal::displayError(node->name() + ":\t" + errorMsg);
}

MStatus Global::maya2TetgenMesh( const MObject& meshObj, tetgenio& tetMesh, MMatrix transform /*= MMatrix()*/ )
{
	MStatus s;
	MFnMesh meshFn(meshObj,&s);
	tetMesh.deinitialize();
	tetMesh.initialize();
	int nVtx = meshFn.numVertices(&s);
	if (nVtx <= 0 || s != MS::kSuccess)
		return MS::kFailure;

	tetMesh.numberofpoints = nVtx;
	tetMesh.pointlist = new REAL[nVtx*3];
	for (int ithVtx = 0; ithVtx < nVtx; ++ithVtx)
	{
		MPoint p;
		meshFn.getPoint(ithVtx, p, MSpace::kWorld);
		p = p * transform;
		if (!_finite(p.x) || !_finite(p.y) || !_finite(p.z))
			return MS::kFailure;
		tetMesh.pointlist[ithVtx*3+0] = p.x;
		tetMesh.pointlist[ithVtx*3+1] = p.y;
		tetMesh.pointlist[ithVtx*3+2] = p.z;
	}

	int nPolys = meshFn.numPolygons(&s);
	if (nPolys <= 0 || s != MS::kSuccess)
		return MS::kFailure;

	tetMesh.numberoffacets = nPolys;
	tetMesh.facetlist = new tetgenio::facet[nPolys];
	MItMeshPolygon it(meshObj, &s);
	int ithFace = 0;
	for (; !it.isDone(&s); it.next())
	{
		MIntArray vIDs;
		it.getVertices(vIDs);
		int nPolyVtx = vIDs.length();

		tetgenio::facet& f = tetMesh.facetlist[ithFace];
		f.holelist = NULL;
		f.numberofholes = 0;

		if (nPolyVtx <= 0)
		{
			f.numberofpolygons = 0;
			f.polygonlist = NULL;
			continue;
		}

		f.numberofpolygons = 1;
		f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
		tetgenio::polygon& p = f.polygonlist[0];
		p.numberofvertices = vIDs.length();
		p.vertexlist = new int[p.numberofvertices];
		for (unsigned ithV = 0; ithV < vIDs.length(); ++ithV)
		{
			p.vertexlist[ithV] = vIDs[ithV];
		}
		ithFace++;
	}
	tetMesh.numberoffacets = ithFace;
	return s;
}

MStatus Global::getVectorArrayData( const MPlug& plug, EigVec& data )
{
	MStatus s;
	MObject vecObj = plug.asMObject();
	MFnVectorArrayData vecFn(vecObj, &s);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	int vecLength = vecFn.length();
	if (vecLength <= 0)
		return MS::kFailure;

	data.setZero(vecLength*3);
	for (int ithVec = 0, ithDof = 0; ithVec < vecLength; ++ithVec, ithDof += 3)
	{
		MVector v = vecFn[ithVec];
		data[ithDof]   = v.x;
		data[ithDof+1] = v.y;
		data[ithDof+2] = v.z;
	}
	return s;
}

MStatus Global::setVectorArrayData( MPlug& plug, const EigVec& data )
{
	MStatus s;
	MVectorArray vectorArray;
	int nVec = data.size()/3;
	vectorArray.setLength(nVec);

	for (int ithVec = 0, ithDof = 0; ithVec < nVec; ++ithVec, ithDof += 3)
	{
		vectorArray[ithVec] = MVector(data[ithDof], data[ithDof+1], data[ithDof+2]);
	}

	MFnVectorArrayData dataFn;
	MObject vectorObj = dataFn.create(vectorArray, &s);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = plug.setMObject(vectorObj);
	return s;
}

MStatus Global::getDoubleArrayData( const MPlug& plug, EigVec& data )
{
	MStatus s;
	MObject doubleObj = plug.asMObject();
	MFnDoubleArrayData doubleFn(doubleObj, &s);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	int length = doubleFn.length();
	if (length <= 0)
		return MS::kFailure;

	data.setZero(length);
	for (int i = 0; i < length; ++i)
	{
		data[i] = doubleFn[i];
	}
	return s;
}

MStatus Global::setDoubleArrayData( MPlug& plug, const EigVec& data )
{
	MStatus s;
	MDoubleArray doubleArray(data.data(), data.size());

	MFnDoubleArrayData dataFn;
	MObject doubleObj = dataFn.create(doubleArray, &s);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = plug.setMObject(doubleObj);
	return s;
}
