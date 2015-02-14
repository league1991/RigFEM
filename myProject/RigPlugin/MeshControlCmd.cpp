#include "StdAfx.h"
#include "MeshControlCmd.h"

const char* MeshControlCmd::m_initFlag[2] = {"-init", "-i"};
const char* MeshControlCmd::m_nameFlag[2] = {"-name", "-n"};

MeshControlCmd::MeshControlCmd(void)
{
}

MeshControlCmd::~MeshControlCmd(void)
{
}

MSyntax MeshControlCmd::newSyntax()
{
	MSyntax syntax;
	MStatus s;
	s = syntax.addFlag(m_initFlag[1], m_initFlag[0], MSyntax::kNoArg);
	s = syntax.addFlag(m_nameFlag[1], m_nameFlag[0], MSyntax::kString);
	return syntax;
}

MStatus MeshControlCmd::doIt( const MArgList& args )
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
	if (!s)return MS::kSuccess;

	bool isInitFlagSet = argData.isFlagSet(m_initFlag[1], &s);

	MItSelectionList pSel(selection, MFn::kDependencyNode , &s);
	MObject obj;
	for (; !pSel.isDone(); pSel.next())
	{
		MObject depNode;
		pSel.getDependNode(depNode);
		MFnDependencyNode nodeFn(depNode, &s);
		MString typeName = nodeFn.typeName(&s);

		if (typeName != NODE_CAGE_DEFORMER_NAME)
			continue;

		CageDeformerNode* node = (CageDeformerNode*)nodeFn.userNode(&s);
		if (!s)
			continue;

		if (isInitFlagSet)
		{
			s = node->computeWeight();
		}
		if (!s)
		{
			MGlobal::displayError("faild to compute weight");
		}
		else
			MGlobal::displayInfo("succeed to compute weight");
	}
	return MS::kSuccess;
}

