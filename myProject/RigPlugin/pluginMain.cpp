#include "stdafx.h"
#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is loaded into Maya.  It 
//		registers all of the services that this plug-in provides with 
//		Maya.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{ 
	MStatus   status;
	MFnPlugin plugin( obj, "Ouyang", "1.0", "Any");

	status = plugin.registerCommand(CMD_EXEC_SIMULATION, RigSimulateCmd::creator, RigSimulateCmd::newSyntax);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	status = plugin.registerNode(	RigSimulationNode::m_nodeName, RigSimulationNode::m_id, 
									RigSimulationNode::creator,	
									RigSimulationNode::initialize, 
									MPxNode::kLocatorNode);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	MGlobal::displayInfo("Simulation plugin is loaded successfully!");
	return status;
}

MStatus uninitializePlugin( MObject obj)
//
//	Description:
//		this method is called when the plug-in is unloaded from Maya. It 
//		deregisters all of the services that it was providing.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{
	MStatus   status;
	MFnPlugin plugin( obj );

	status = plugin.deregisterCommand(CMD_EXEC_SIMULATION);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	status = plugin.deregisterNode(RigSimulationNode::m_id);
	CHECK_MSTATUS_AND_RETURN_IT(status);

	MGlobal::displayInfo("Simulation plugin is unloaded successfully!");

	return status;
}
