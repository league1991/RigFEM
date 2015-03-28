#pragma once


#define CMD_EXEC_SIMULATION			"rigSimulate"
#define CMD_MESH_CONTROL			"cageDeform"

#define NODE_FEM_SIMULATION_NAME	"rigSimulator"
#define NODE_FEM_SIMULATION_ID		4500001

#define NODE_CAGE_DEFORMER_NAME		"cageDeformer"
#define NODE_CAGE_DEFORMER_ID		4500002

#define NODE_DATA_SWITCH_NAME		"dataSwitch"
#define NODE_DATA_SWITCH_ID			4500003

class Global
{
public:
	static MPlug getPlug(MPxNode* node, const char* longName);

	static MColor getColor(const MPlug& plug, MStatus* s = NULL);

	static void displayError(const MPxNode* node, const MString& errorMsg);

	static MStatus maya2TetgenMesh(const MObject& meshObj, tetgenio& tetMesh, MMatrix transform = MMatrix());

	static MStatus	getVectorArrayData(const MPlug& plug, EigVec& data);
	static MStatus	setVectorArrayData(MPlug& plug, const EigVec& data);

	static MStatus  getDoubleArrayData(const MPlug& plug, EigVec& data);
	static MStatus	setDoubleArrayData(MPlug& plug, const EigVec& data);

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
