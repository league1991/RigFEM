#include "StdAfx.h"
#include "Simple1DNoiseNode.h"


MObject Simple1DNoiseNode::m_phase;
MObject Simple1DNoiseNode::m_input;
MObject Simple1DNoiseNode::m_outputVal;
MObject Simple1DNoiseNode::m_period;
MObject Simple1DNoiseNode::m_maxVal;
MObject Simple1DNoiseNode::m_minVal;

const char* Simple1DNoiseNode::m_nodeName = NODE_SIMPLE_1D_NOISE_NAME;
MTypeId Simple1DNoiseNode::m_id(NODE_SIMPLE_1D_NOISE_ID);

const char* Simple1DNoiseNode::m_phaseName[2] ={"phase", "phase"};
const char* Simple1DNoiseNode::m_inputName[2] = {"input", "in"};
const char* Simple1DNoiseNode::m_outputValName[2] = {"output", "out"};
const char* Simple1DNoiseNode::m_periodName[2] = {"period", "period"};
const char* Simple1DNoiseNode::m_maxValName[2] = {"maxValue", "maxVal"};
const char* Simple1DNoiseNode::m_minValName[2] = {"minValue", "minVal"};

Simple1DNoiseNode::Simple1DNoiseNode(void)
{
}

Simple1DNoiseNode::~Simple1DNoiseNode(void)
{
}


MStatus Simple1DNoiseNode::initialize()
{
	MStatus s; 
	MFnTypedAttribute tAttr;
	MFnNumericAttribute nAttr;
	MFnMatrixAttribute  mAttr;

	m_minVal = nAttr.create(m_minValName[0], m_minValName[1], MFnNumericData::kDouble,0);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);


	m_maxVal = nAttr.create(m_maxValName[0], m_maxValName[1], MFnNumericData::kDouble,1);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);


	m_period = nAttr.create(m_periodName[0], m_periodName[1], MFnNumericData::kDouble,1);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setMin(0);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);



	m_phase = nAttr.create(m_phaseName[0], m_phaseName[1], MFnNumericData::kDouble,0);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);

	m_outputVal = nAttr.create(m_outputValName[0], m_outputValName[1], MFnNumericData::kDouble,0);
	nAttr.setKeyable(false);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(false);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);

	m_input = nAttr.create(m_inputName[0], m_inputName[1], MFnNumericData::kDouble,0);
	nAttr.setKeyable(true);
	nAttr.setStorable(true);
	nAttr.setHidden(false);
	nAttr.setWritable(true);
	nAttr.setReadable(true); 
	nAttr.setSoftMin(0);
	nAttr.setSoftMax(1);

	s = addAttribute(m_minVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_maxVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_phase);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_period);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_outputVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = addAttribute(m_input);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	s = attributeAffects(m_input, m_outputVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_minVal, m_outputVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_maxVal, m_outputVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_phase, m_outputVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);
	s = attributeAffects(m_period, m_outputVal);
	CHECK_MSTATUS_AND_RETURN_IT(s);

	return s;
}

double Simple1DNoiseNode::computeNoise( double x, double minVal, double maxVal, double period, double phase )
{
	if (period == 0)
		return (minVal + maxVal) * 0.5;
	double periodCoord = (x+phase) / period;
	int periodIdx = (int)floor(periodCoord);
	int nextPeriodIdx = periodIdx+1;
	double periodFrag  = periodCoord - periodIdx;

	double v0 = randDouble(periodIdx);
	double v1 = randDouble(nextPeriodIdx);

	double w  = periodFrag * (3 * periodFrag - 2 * periodFrag * periodFrag);
	double v  = v0 * (1-w) + v1 * w;

	double norV= (v + 1) / 2.0;
	return minVal + (maxVal - minVal) * norV;
}

double Simple1DNoiseNode::randDouble( unsigned x )
{
	x = (x<<13) ^ x;
	return ( 1.0 - ( (x * (x * x * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
}

MStatus Simple1DNoiseNode::compute( const MPlug& plug, MDataBlock& data )
{
	MStatus s;
	if (plug == m_outputVal)
	{
		MPlug minArrayPlug = Global::getPlug(this, m_minValName[0]);
		MPlug maxArrayPlug = Global::getPlug(this, m_maxValName[0]);
		MPlug periodArrayPlug = Global::getPlug(this, m_periodName[0]);
		MPlug phaseArrayPlug = Global::getPlug(this, m_phaseName[0]);
		MPlug outputArrayPlug = Global::getPlug(this, m_outputValName[0]);
		MPlug inputPlug    = Global::getPlug(this, m_inputName[0]);

		double x = inputPlug.asDouble();

		for (int i = 0; i < outputArrayPlug.numElements(); ++i)
		{
			MPlug outPlug = outputArrayPlug.elementByPhysicalIndex(i);
			int logIdx = outPlug.logicalIndex();

			MPlug minPlug = minArrayPlug.elementByLogicalIndex(logIdx);
			MPlug maxPlug = maxArrayPlug.elementByLogicalIndex(logIdx);
			MPlug periodPlug= periodArrayPlug.elementByLogicalIndex(logIdx);
			MPlug phasePlug = phaseArrayPlug.elementByLogicalIndex(logIdx);

			double res = computeNoise(x, minPlug.asDouble(), maxPlug.asDouble(), periodPlug.asDouble(), phasePlug.asDouble());
			outPlug.setDouble(res);
		}
	}
	data.setClean(plug);
	return s;
}



