#pragma once

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define CLAMP_INT(minValue,maxValue,value) (MAX(int(minValue),MIN(int(value),int(maxValue))))
#define CLAMP_FLOAT(minValue,maxValue,value) (MAX(float(minValue),MIN(float(value),float(maxValue))))

#define CMD_EXEC_SIMULATION			"rigSimulate"

#define NODE_FEM_SIMULATION_NAME	"rigSimulator"
#define NODE_FEM_SIMULATION_ID		4500001

class Global
{
public:
	static MPlug getPlug(MPxNode* node, const char* longName)
	{
		MObject nodeObj = node->thisMObject();
		MFnDependencyNode nodeFn(nodeObj);
		return nodeFn.findPlug(longName, true);
	}

	static MColor getColor(const MPlug& plug, MStatus* s = NULL)
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

	static void displayError(const MPxNode* node, const MString& errorMsg)
	{
		MGlobal::displayError(node->name() + ":\t" + errorMsg);
	}

	static MStatus maya2TetgenMesh(const MObject& meshObj, tetgenio& tetMesh, MMatrix transform = MMatrix())
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
		for (int ithFace=0; !it.isDone(&s); it.next(),ithFace++)
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
		}
		return s;
	}

	template<class MatrixType>
	static void showMatrix(const MatrixType& mat, const char* name)
	{
		std::string str = Utilities::matToString(mat, name);
		MGlobal::displayInfo(str.data());
	}

	template<class VectorType>
	static void showVector(const VectorType& vec, const char* name)
	{
		std::string str = Utilities::vecToString(vec, name);
		MGlobal::displayInfo(str.data());
	}
};
