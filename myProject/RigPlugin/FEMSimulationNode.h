#pragma once


class RigSimulationNode: public MPxLocatorNode
{
public:
	RigSimulationNode(void);
	virtual ~RigSimulationNode(void);
	void postConstructor();

	virtual void draw( M3dView & view, const MDagPath & path, M3dView::DisplayStyle style, M3dView:: DisplayStatus );
	virtual bool isBounded() const{return true;}
	virtual MBoundingBox boundingBox() const;

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );
	static  void*		creator();
	static  MStatus		initialize();

	static MTypeId		m_id;
	static const char*  m_nodeName;

	// 各种Plug的名字
	static const char*	paramLongName(){return m_paramName[0];}
	static const char*	meshLongName(){return m_meshName[0];}
	static const char*  transformLongName(){return m_transformMatrixName[0];}

	// 模拟相关函数
	void				clearRig();
	bool				resetRig();
	bool				stepRig();
	MStatus				setParamPlug( const double* param, int nParam);

	// 用于测试
	void				testRig();
private:
	// 获取当前有效的参数个数
	int					getNumParam();
	MMatrix				getMatrix();
	bool				getInitStatus(RigFEM::RigStatus& s);
	int					getCurFrame();

	MBoundingBox		m_box;

	static MObject		m_initParam;
	static MObject		m_param;
	static MObject		m_mesh;
	static MObject      m_transformMatrix;		// mesh transform matrix

	static const char*  m_initParamName[2];
	static const char*  m_paramName[2];
	static const char*  m_meshName[2];
	static const char*  m_transformMatrixName[2];

				
	RigFEM::GeneralRig*	m_rig;					// 封装了节点求值机制	
	RigFEM::RiggedMesh*	m_rigMesh;				// fem 数据
	RigFEM::NewtonSolver* m_solver;				// 求解器
	RigFEM::StatusRecorder m_recorder;			// 记录求解状态		
};
