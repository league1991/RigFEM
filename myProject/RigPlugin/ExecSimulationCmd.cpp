#include "StdAfx.h"
#include "ExecSimulationCmd.h"

const char* RigSimulateCmd::m_initFlag[2] = {"-init", "-i"};
const char* RigSimulateCmd::m_stepFlag[2] = {"-step", "-s"};
const char* RigSimulateCmd::m_stepStaticFlag[2] = {"-stepStatic", "-ss"};
const char* RigSimulateCmd::m_nameFlag[2] = {"-name", "-n"};
const char* RigSimulateCmd::m_hessianFlag[2] = {"-hessian", "-h"};
const char* RigSimulateCmd::m_gradFlag[2] = {"-gradient", "-g"};
const char* RigSimulateCmd::m_saveFlag[2] = {"-save", "-sa"};
const char* RigSimulateCmd::m_surfPntFlag[2] = {"-numSurfacePnt",	"-nsp"};
const char* RigSimulateCmd::m_intPntFlag[2] = {"-numInternalPnt",   "-nip"};
const char* RigSimulateCmd::m_createFlag[2] = {"-create", "-c"};
const char* RigSimulateCmd::m_recordEleGFFlag[2] = {"-stepEleGF", "-rgf"};
const char* RigSimulateCmd::m_saveEleGFFlag[2] = {"-saveGF", "-sgf"};
const char* RigSimulateCmd::m_loadEleMatFlag[2] = {"-loadEleMat", "-lm"};
const char* RigSimulateCmd::m_resetEleMatFlag[2] = {"-resetEleMat","-rsm"};



RigSimulateCmd::RigSimulateCmd(void)
{
} 

RigSimulateCmd::~RigSimulateCmd(void)
{
}

MStatus RigSimulateCmd::doIt( const MArgList& args )
{
	MStatus s; 

	MSelectionList selection;
	MSyntax syn = syntax();
	MArgDatabase argData(syn, args);
	setResult(0);
	if (argData.isFlagSet(m_nameFlag[1], &s))
	{
		MString nodeName;											// 提取某个名字的节点
		s = argData.getFlagArgument(m_nameFlag[1],0, nodeName);
		s = MGlobal::getSelectionListByName(nodeName, selection);
	}
	else
		s = MGlobal::getActiveSelectionList(selection);				// 提取选中节点
	if (!s)return MS::kFailure;

	bool isInitFlagSet = argData.isFlagSet(m_initFlag[1], &s);
	bool isStepFlagSet = argData.isFlagSet(m_stepFlag[1], &s); 
	bool isHessianFlagSet = argData.isFlagSet(m_hessianFlag[1], &s);
	bool isGradFlagSet = argData.isFlagSet(m_gradFlag[1], &s);
	bool isSaveFlagSet = argData.isFlagSet(m_saveFlag[1], &s);
	bool isStepStaticSet = argData.isFlagSet(m_stepStaticFlag[1], &s);
	bool isIntPFlagSet = argData.isFlagSet(m_intPntFlag[1], &s);
	bool isSurfPFlagSet= argData.isFlagSet(m_surfPntFlag[1], &s);
	bool isCreateFlagSet=  argData.isFlagSet(m_createFlag[1], &s);
	bool isRecordGFFlagSet=argData.isFlagSet(m_recordEleGFFlag[1], &s);
	bool isSaveEleGFFlagSet=argData.isFlagSet(m_saveEleGFFlag[1], &s);
	bool isLoadEleMatFlagSet= argData.isFlagSet(m_loadEleMatFlag[1], &s);
	bool isResetEleMatFlagSet= argData.isFlagSet(m_resetEleMatFlag[1], &s);

	if (isCreateFlagSet)
	{
		const char* createCmd = 
			"{													\
			string $nodeName = `createNode rigSimulator`;		\
			connectAttr time1.outTime ($nodeName + \".time\");	\
			}";
		PRINT_F("command: %s", createCmd);
		MGlobal::executeCommand(createCmd);
	}

	MItSelectionList pSel(selection, MFn::kDependencyNode , &s);
	MObject obj;
	for (; !pSel.isDone(); pSel.next())
	{
		MObject depNode;
		pSel.getDependNode(depNode);
		MFnDependencyNode nodeFn(depNode, &s);
		MString typeName = nodeFn.typeName(&s);

		// 找到模拟节点，执行之
		if (typeName == NODE_FEM_SIMULATION_NAME)
		{
			RigSimulationNode* node = (RigSimulationNode*)nodeFn.userNode(&s);
			if (!s)
				continue;

			bool res = true;
			MString name = nodeFn.name(&s);
			if (isInitFlagSet)
			{
				res &= node->resetRig();
			}
			else if (isStepFlagSet)
			{
				res &= node->stepRig();
			}
			else if (isStepStaticSet)
			{
				res &= node->staticStepRig();
			}
			else if (isRecordGFFlagSet)
			{
				res &= node->stepWithEleHessianOrGF();
			}
			else if (isHessianFlagSet)
			{
				double noiseN = 1, noiseP = 1;
				argData.getFlagArgument(m_hessianFlag[1], 0, noiseN);
				argData.getFlagArgument(m_hessianFlag[1], 1, noiseP);
				res &= node->testHessian(noiseN, noiseP);
			}
			else if (isGradFlagSet)
			{
				double noiseN = 1, noiseP = 1;
				argData.getFlagArgument(m_gradFlag[1], 0, noiseN);
				argData.getFlagArgument(m_gradFlag[1], 1, noiseP);
				res &= node->testGrad(noiseN, noiseP);
			}
			else if (isSaveFlagSet)
			{
				MString filePath;
				argData.getFlagArgument(m_saveFlag[1], 0, filePath);
				res &= node->saveSimulationData(filePath.asChar());
			}
			else if (isSaveEleGFFlagSet)
			{
				MString filePath;
				argData.getFlagArgument(m_saveEleGFFlag[1], 0, filePath);
				res &= node->saveGFResult(filePath.asChar());
			}
			else if (isLoadEleMatFlagSet)
			{
				res &= node->loadElementMaterial();
			}
			else if (isResetEleMatFlagSet)
			{
				res &= node->resetElementMaterial();
			}
			else if (isIntPFlagSet)
			{
				int nPnt = node->getNumInternalPnt();
				setResult(nPnt);
				return MS::kSuccess;
			}
			else if (isSurfPFlagSet)
			{
				int nPnt = node->getNumSurfacePnt();
				setResult(nPnt);
				return MS::kSuccess;
			}

			if (res)
			{
				MString info = "rigSimulate command succeed in " + name;
				MGlobal::displayInfo(info);
			}
			else
			{
				MString info = "rigSimulate command failed in " + name;
				MGlobal::displayError(info);
				return MS::kSuccess;
			}
		}
	}
	setResult(1);
	return MS::kSuccess;
}

void* RigSimulateCmd::creator()
{
	return new RigSimulateCmd;
}

MSyntax RigSimulateCmd::newSyntax()
{
	MSyntax syntax;
	MStatus s;
	s = syntax.addFlag(m_initFlag[1], m_initFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_stepFlag[1], m_stepFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_recordEleGFFlag[1], m_recordEleGFFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_saveEleGFFlag[1], m_saveEleGFFlag[0], MSyntax::kString);
	s = syntax.addFlag(m_intPntFlag[1], m_intPntFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_surfPntFlag[1], m_surfPntFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_createFlag[1], m_createFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_stepStaticFlag[1], m_stepStaticFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_nameFlag[1], m_nameFlag[0], MSyntax::kString);
	s = syntax.addFlag(m_hessianFlag[1], m_hessianFlag[0], MSyntax::kDouble, MSyntax::kDouble);
	s = syntax.addFlag(m_gradFlag[1], m_gradFlag[0], MSyntax::kDouble, MSyntax::kDouble);
	s = syntax.addFlag(m_saveFlag[1], m_saveFlag[0], MSyntax::kString);
	s = syntax.addFlag(m_loadEleMatFlag[1], m_loadEleMatFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_resetEleMatFlag[1], m_resetEleMatFlag[0], MSyntax::kNoArg);
	return syntax;
}


