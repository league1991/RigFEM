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
	static const char*  initParamLongName(){return m_initParamName[0];}

	// 模拟相关函数
	void				clearRig();
	bool				resetRig();
	bool				stepRig();
	MStatus				setParamPlug( const double* param, int nParam);

	// 用于测试
	void				testRig();

private:
	
	int					getNumParam();			// 获取当前有效的参数个数
	MMatrix				getMatrix();
	bool				getInitStatus(RigFEM::RigStatus& s);
	int					getCurFrame();
	MStatus				setParamToInit();		// 设置参数为初始参数
	void				drawIcon();
	MBoundingBox		m_box;

	static MObject		m_initParam;
	static MObject		m_param;
	static MObject		m_mesh;
	static MObject      m_transformMatrix;		// mesh transform matrix
	static MObject		m_tetEdgeRatio;			// 体网格化参数，四面体边比例
	static MObject		m_tetMaxVolume;			// 体网格化参数，四面体最大体积
	static MObject		m_youngModulus;			// 杨氏模量
	static MObject		m_nu;					// 控制不同轴向变形相互影响的参数
	static MObject		m_density;				// 密度
	static MObject		m_timeStep;				// 时间步长
	static MObject		m_deriStep;				// 参数求导时的有限差商大小

	static const char*  m_initParamName[2];
	static const char*  m_paramName[2];
	static const char*  m_meshName[2];
	static const char*  m_transformMatrixName[2];
	static const char*	m_tetEdgeRatioName[2];			// 体网格化参数，四面体边比例
	static const char*	m_tetMaxVolumeName[2];			// 体网格化参数，四面体最大体积
	static const char*	m_youngModulusName[2];			// 杨氏模量
	static const char*	m_nuName[2];					// 控制不同轴向变形影响程度的参数
	static const char*	m_densityName[2];				// 密度
	static const char*	m_timeStepName[2];				// 时间步长
	static const char*  m_deriStepName[2];					// 参数求导时的有限差商大小
				
	RigFEM::GeneralRig*	m_rig;					// 封装了节点求值机制	
	RigFEM::RiggedMesh*	m_rigMesh;				// fem 数据
	RigFEM::NewtonSolver* m_solver;				// 求解器
	RigFEM::StatusRecorder m_recorder;			// 记录求解状态		
};
