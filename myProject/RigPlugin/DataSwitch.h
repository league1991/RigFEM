#pragma once

class DataSwitchNode :
	public MPxNode
{
public:
	enum DataSource
	{
		SRC_1 = 0,
		SRC_2 = 1
	};
	DataSwitchNode(void);
	~DataSwitchNode(void);
	void postConstructor();
	static void*		creator(){return new DataSwitchNode;}
	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );
	static  MStatus		initialize();

	static MTypeId		m_id;
	static const char*  m_nodeName;
private:
	MStatus				transferData();
	static MObject		m_firstIn;
	static MObject		m_secondIn;
	static MObject		m_out;
	static MObject		m_source;

	static const char*  m_firstInName[2];
	static const char*  m_secondInName[2];
	static const char*  m_outName[2];
	static const char*  m_sourceName[2];
};
