#include "StdAfx.h"
#include "ExecSimulationCmd.h"

const char* RigSimulateCmd::m_initFlag[2] = {"-init", "-i"};
const char* RigSimulateCmd::m_stepFlag[2] = {"-step", "-s"};
const char* RigSimulateCmd::m_nameFlag[2] = {"-name", "-n"};

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
	if (argData.isFlagSet(m_nameFlag[1], &s))
	{
		MString nodeName;											// 提取某个名字的节点
		s = argData.getFlagArgument(m_nameFlag[1],0, nodeName);
		s = MGlobal::getSelectionListByName(nodeName, selection);
	}
	else
		s = MGlobal::getActiveSelectionList(selection);				// 提取选中节点
	CHECK_MSTATUS_AND_RETURN_IT(s);

	bool isInitFlagSet = argData.isFlagSet(m_initFlag[1], &s);
	bool isStepFlagSet = argData.isFlagSet(m_stepFlag[1], &s);

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
			CHECK_MSTATUS_AND_RETURN_IT(s);

			bool res = true;
			MString name = nodeFn.name(&s);
			if (isInitFlagSet)
			{
				res &= node->resetRig();
			}
			if (isStepFlagSet)
			{
				res &= node->stepRig();
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
			}
		}
	}
	return s;
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
	s = syntax.addFlag(m_nameFlag[1], m_nameFlag[0], MSyntax::kString);
	return syntax;
}

