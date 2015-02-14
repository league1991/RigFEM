#pragma once

class CageDeformerNode :
	public MPxNode
{
public:
	CageDeformerNode(void);
	~CageDeformerNode(void);

	static void*		creator(){return new CageDeformerNode;}
	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );
	static  MStatus		initialize();

	MStatus				computeWeight();
	static const MTypeId m_id;
	static const char*  m_nodeName;
private:
	MStatus				computeNewPnt();

	MStatus				getVertices(MObject& meshObj, EigDense& vertMat)const;
	static MObject		m_refControlMesh;
	static MObject		m_refMesh;
	static MObject		m_inControlMesh;
	static MObject		m_outMesh;
	static MObject		m_maxNeigh;
	static MObject		m_power;

	static const char*  m_refControlMeshName[2];
	static const char*  m_refMeshName[2];
	static const char*  m_inControlMeshName[2];
	static const char*  m_outMeshName[2];
	static const char*  m_maxNeighName[2];
	static const char*  m_powerName[2];

	EigSparse			m_weightMat;
};
