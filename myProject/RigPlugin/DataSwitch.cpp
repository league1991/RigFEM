#include "StdAfx.h"
#include "DataSwitch.h"

const char* DataSwitchNode::m_sourceName[2] = {"source", "src"};
const char* DataSwitchNode::m_outName[2]    = {"out","out"};
const char* DataSwitchNode::m_secondInName[2]  = {"inSecond", "inSec"};
const char* DataSwitchNode::m_firstInName[2]   = {"inFirst", "inFir"};

MObject DataSwitchNode::m_source;
MObject DataSwitchNode::m_out;
MObject DataSwitchNode::m_secondIn;
MObject DataSwitchNode::m_firstIn;

const char* DataSwitchNode::m_nodeName = NODE_DATA_SWITCH_NAME;
MTypeId DataSwitchNode::m_id(NODE_DATA_SWITCH_ID);

DataSwitchNode::DataSwitchNode(void)
{
}

DataSwitchNode::~DataSwitchNode(void)
{
}

MStatus DataSwitchNode::compute( const MPlug& plug, MDataBlock& data )
{
	MStatus s;
	if (plug == m_out || plug.parent() == m_out)
	{
		s = transferData();
	}
	data.setClean(plug);
	return s;
}

MStatus DataSwitchNode::initialize()
{
	MStatus s;
	MFnNumericAttribute nAttr;
	MFnEnumAttribute  eAttr;

	m_firstIn = nAttr.create(m_firstInName[0], m_firstInName[1], MFnNumericData::kDouble,0);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);

	m_secondIn = nAttr.create(m_secondInName[0], m_secondInName[1], MFnNumericData::kDouble,0);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);

	m_out = nAttr.create(m_outName[0], m_outName[1], MFnNumericData::kDouble,0);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(false);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);


	m_source = eAttr.create(m_sourceName[0], m_sourceName[1], SRC_1);
	eAttr.addField("1st In", SRC_1);
	eAttr.addField("2nd In", SRC_2);
	eAttr.setHidden(false);
	eAttr.setReadable(true);
	eAttr.setWritable(true);

	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_firstIn);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_secondIn);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_out);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_source);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	s = attributeAffects(m_firstIn, m_out);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_secondIn, m_out);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_source, m_out);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	return MStatus::kSuccess;
}

MStatus DataSwitchNode::transferData()
{
	MStatus s;
	MPlug sourcePlug = Global::getPlug(this, m_sourceName[0]);
	const char* inName;
	switch (sourcePlug.asShort())
	{
	case SRC_1:
		inName = m_firstInName[0];	break;
	case SRC_2:
		inName = m_secondInName[0]; break;
	default:
		inName = m_firstInName[0];	break;
	}
	MPlug inPlug     = Global::getPlug(this, inName);
	MPlug outPlug    = Global::getPlug(this, m_outName[0]);

	for (int ithOut = 0; ithOut < outPlug.numElements(); ++ithOut)
	{
		MPlug outElePlug = outPlug.elementByPhysicalIndex(ithOut, &s);
		CHECK_MSTATUS_AND_RETURN_IT(s);
		if (!outElePlug.isConnected())
			continue;
		int logIdx = outElePlug.logicalIndex(&s);
		MPlug inElePlug  = inPlug.elementByLogicalIndex(logIdx, &s);
		CHECK_MSTATUS_AND_RETURN_IT(s);

		outElePlug.setValue(inElePlug.asDouble());
	}
	return s;
}

void DataSwitchNode::postConstructor()
{
	this->setExistWithoutInConnections(true);
	this->setExistWithoutOutConnections(true);
}
