#pragma once


class RigSimulationNode: public MPxLocatorNode
{
public:
	enum SimulationType
	{
		SIM_STANDARD = 0,
		SIM_SKIN     = 1
	};
	enum DisplayType
	{
		DISP_SIM	 = 0,
		DISP_INIT    = 1
	};
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
	bool				stepRig();					// 动态模拟
	bool				staticStepRig();			// 静态模拟
	bool				stepWithEleGF();			// 静态模拟，并保持每个元素的广义力
	bool				saveGFResult(const char* fileName);
	bool				loadElementMaterial();
	bool				resetElementMaterial();

	MStatus				setParamPlug( const double* param, int nParam);

	// 用于测试
	void				testRig();	
	bool				testHessian(double noiseN, double noiseP);					// 测试当前帧的Hessian
	bool				testGrad(double noiseN, double noiseP);

	// 保存数据
	bool				saveSimulationData(const char* fileName);
	bool				computeExternalForce( const EigVec& pos, const EigVec& vel, const EigVec& m, double time, EigVec& extForce , EigVec& surfForce);
	bool				getControlParams(EigVec& targetParam, EigVec& propGain, EigVec& deriGain);
	bool				getControlGain(EigVec& propGain, EigVec& deriGain);
	bool				getControlTarget(EigVec& targetParam);

	int					getNumInternalPnt();
	int					getNumSurfacePnt();
private:
	
	int					getNumParam();			// 获取当前有效的参数个数
	MMatrix				getMatrix();
	bool				getInitStatus(RigFEM::RigStatus& s);
	bool				getInitParam(EigVec& p);
	int					getCurFrame();
	MStatus				setParamToInit();		// 设置参数为初始参数
	void				drawIcon();
	DisplayType			getDispType();
	
	bool				updateDeriStepSize();	// 更新有限差商导数步长
	bool				updateTerminationCond();// 更新终止条件

	// 节点属性
	void				printFieldData();
	MStatus				getInputForce( EigVec& fieldForce, EigVec& surfForce );
	MStatus				testField();

	// 显示配置
	RigFEM::MeshDispConfig getDisplayConfig();

	MBoundingBox		m_box;

	static MObject		m_initParam;			// 参数初始值
	static MObject		m_param;				// 参数值
	static MObject		m_mesh;					// 表面网格
	static MObject      m_transformMatrix;		// mesh transform matrix
	static MObject		m_tetEdgeRatio;			// 体网格化参数，四面体边比例
	static MObject		m_tetMaxVolume;			// 体网格化参数，四面体最大体积
	static MObject		m_youngModulus;			// 杨氏模量
	static MObject		m_nu;					// 控制不同轴向变形相互影响的参数
	static MObject		m_density;				// 密度
	static MObject		m_timeStep;				// 时间步长
	static MObject		m_dispType;				// 显示方式,直接显示模拟结果还是初始参数
	static MObject		m_dispFemMesh;			// 是否显示fem网格
	static MObject		m_dispVertex;			// 是否显示顶点
	static MObject		m_dispEdge;				// 是否显示边
	static MObject		m_dispBBox;				// 是否显示包围盒
	static MObject      m_fieldForceFactor;		// 力线长度因子
	static MObject		m_time;					// 输入的时间，用于驱动节点求值

	// 牛顿法参数
	static MObject		m_deriStep;				// 参数求导时的有限差商大小
	static MObject		m_maxIter;				// 最大迭代次数
	static MObject		m_minStepSize;			// 迭代的最小步长，若步长小于此值，终止迭代
	static MObject		m_minGradSize;			// 最小的梯度值，如梯度小于此值，终止迭代
	static MObject		m_simType;				// 模拟算法选择
	static MObject		m_weightPath;			// 权重路径
	static MObject		m_resultPath;			// 结果路径
	static MObject		m_maxParamStep;			// 每次迭代参数的最大增量

	// 自适应材料硬度
	static MObject		m_GFResultPath;			// 广义力的保存路径
	static MObject		m_materialPath;			// 材料相对硬度文件
	static MObject		m_displayMaterial;		// 是否显示材料硬度
	static MObject		m_minFactor;			// 颜色范围
	static MObject		m_maxFactor;			// 
	static MObject		m_cutAxis;				// 切割轴向
	static MObject		m_cutRatio;				// 切割比例

	// 共轭梯度法参数
	static MObject		m_cgMinStepSize;		// 迭代的最小步长，若步长小于此值，终止迭代
	static MObject		m_maxCGIter;			// 共轭梯度法最大迭代次数

	// maya场参数
	static MObject		m_inputForce;			// 场的输入力（含内部与表面节点）
	static MObject		m_fieldData;			// 输出的有限元节点状态
	static MObject		m_fieldDataPosition;	// 位置子属性
	static MObject		m_fieldDataVelocity;	// 速度子属性
	static MObject		m_fieldDataMass;		// 质量子属性
	static MObject		m_fieldDataDeltaTime;	// 时间子属性

	// 碰撞
	static MObject		m_surfaceForce;			// 表面受到的力
	static MObject		m_surfaceForceBitmap;	// 记录各个点表面力是否有效的位图
	static MObject      m_surfaceForceFactor;	// 表面力强度因子
	static MObject		m_collisionParticle;	// 用于代理碰撞检测的n粒子形状节点
	static MObject		m_collisionMesh;		// 粒子追踪的碰撞物体

	// PD自动控制参数
	static MObject		m_controlType;			// 自动控制开关
	static MObject		m_targetParam;			// 期望的状态	
	static MObject		m_proportionalGain;		// 比例控制强度
	static MObject		m_derivativeGainRatio;	// 微分控制强度比率

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
	static const char*  m_deriStepName[2];				// 参数求导时的有限差商大小
	static const char*  m_maxIterName[2];				// 最大迭代次数
	static const char*  m_minStepSizeName[2];			// 迭代的最小步长，若步长小于此值，终止迭代
	static const char*  m_minGradSizeName[2];			// 最小的梯度值，如梯度小于此值，终止迭代
	static const char*  m_simTypeName[2];				// 模拟算法选择
	static const char*	m_weightPathName[2];			// 权重路径
	static const char*  m_resultPathName[2];			// 结果路径
	static const char*  m_GFResultPathName[2];			// 广义力保存路径
	static const char*	m_maxParamStepName[2];			// 每次迭代参数的最大增量
	static const char*  m_dispTypeName[2];
	static const char*  m_dispFemMeshName[2];			// 是否显示fem网格
	static const char*	m_dispVertexName[2];			// 是否显示顶点
	static const char*	m_dispEdgeName[2];				// 是否显示边
	static const char*	m_dispBBoxName[2];				// 是否显示包围盒
	static const char*  m_cgMinStepSizeName[2];
	static const char*  m_maxCGIterName[2];				// 共轭梯度法最大迭代次数
	static const char*  m_inputForceName[2];			// 场的输入力
	static const char*  m_fieldDataName[2];				// 输出的有限元节点状态
	static const char*  m_fieldDataPositionName[2];		// 子属性
	static const char*  m_fieldDataVelocityName[2];
	static const char*  m_fieldDataMassName[2];
	static const char*  m_fieldDataDeltaTimeName[2];
	static const char*  m_fieldForceFactorName[2];
	static const char*	m_controlTypeName[2];			// 自动控制开关
	static const char*	m_targetParamName[2];			// 期望的状态	
	static const char*	m_proportionalGainName[2];		// 比例控制强度
	static const char*	m_derivativeGainRatioName[2];	// 微分控制强度
	static const char*  m_surfaceForceName[2];
	static const char*  m_collisionParticleName[2];
	static const char*  m_collisionMeshName[2];
	static const char*	m_surfaceForceBitmapName[2];	// 记录各个点表面力是否有效的位图
	static const char*  m_surfaceForceFactorName[2];	// 表面力强度因子
	static const char*  m_timeName[2];
	static const char*	m_displayMaterialName[2];		// 是否显示材料硬度
	static const char*	m_minFactorName[2];				// 颜色范围
	static const char*	m_maxFactorName[2];				// 
	static const char*	m_cutAxisName[2];				// 切割轴向
	static const char*	m_cutRatioName[2];				// 切割比例
	static const char*  m_materialPathName[2];			// 材料相对硬度文件


	SimulationType		m_simTypeFlag;

	RigFEM::SimulatorBase*	m_simulator;				// 模拟器
};
