#pragma once

using namespace std;
//#define SHOW_TETGEN_RESULT		// 是否显示tetgen生成的out数据
#define SHOW_TETMESH_RESULT			// 是否显示Vega的体网格数据
namespace RigFEM
{

enum RigControlType
{
	CONTROL_NONE			= (0x1),			// 不加入任何控制力
	CONTROL_EXPLICIT_FORCE	= (0x1) << 1,		// 在每一帧求解的过程中用相同的控制力
	CONTROL_IMPLICIT_FORCE	= (0x1) << 2		// 在每一帧求解的过程中用可变的控制力
};
struct MeshDispConfig
{
	MeshDispConfig();
	void setDefault();

	double m_extForceDispFactor;
	double m_ctrlForceDispFactor;
};
class RiggedMesh :public ObjectFunction
{
public:
	RiggedMesh(void);
	~RiggedMesh(void);
 
	void init();
	// 用表面网格和rig对象初始化，rig由外部负责释放
	bool init( tetgenio& surfMesh, RigBase* rig, double maxVolume = 1, double edgeRatio = 2 , double youngModulus = 1e6, double nu = 0.45, double density = 1000);
	void clear();
	void show();
	bool showStatus( RigStatus& s, const MeshDispConfig& config , double* bbox = NULL);
	void computeRig();

	// 各种计算函数
	// 计算给定状态x = [n，p] 以及时间参数param下的函数值，以及梯度
	bool computeValueAndGrad(const EigVec& x, const EigVec& param, double* v, EigVec* grad);
	// 根据给定的内部点的偏移量n以及参数值p计算函数值、梯度以及Hessian矩阵
	double computeValue(const EigVec& n, const EigVec& p, double t);
	bool computeGradient(const EigVec& n, const EigVec& p, double t, EigVec& gn, EigVec& gp);
	bool computeHessian(const EigVec& n, const EigVec& p, double t, EigSparse& Hnn, EigDense& Hnp, EigDense& Hpn, EigDense* pHpp = NULL);		// pHpp 若为NULL，不计算
	// 给定参数以确定表面点位置，计算内部点受力平衡的位置
	// initN可以指定初始迭代位置
	bool computeStaticPos(const EigVec& p, double t, EigVec& q, int maxIter = 20, const EigVec* initQ = NULL);

	// 各种状态变量
	// 获得当前各个自由度的状态
	void getDof(EigVec& n, EigVec& p);
	// 更新各个自由度的状态，同时更新的有顶点位置、速度、加速度
	void setDof(EigVec&n, EigVec&p, bool proceedTime = true);
	// 外力，在每一时刻，每个节点的外力视为恒定
	bool setExternalForce(const EigVec&f);
	const EigVec& getExternalForce()const;
	// 质量
	const EigVec& getMass()const;
	// 时间
	double getCurTime(){return m_t;}
	void setStepTime(double dt){m_h = dt;}
	double getStepTime()const{return m_h;}
	// 顶点坐标
	void getVertexPosition(EigVec& pos);
	// 更新外力，包括场力，此外力在求解过程中保持不变
	bool updateExternalAndControlForce();
	// 控制
	void setRigControlType(RigControlType type);
	const EigVec& getControlTargetParam()const;
	const EigVec& getControlTargetVelocity()const;
	bool getControlTargetFromRigNode(EigVec& target);
	void setControlTarget(const EigVec& targetParam, const EigVec& targetParamVelocity);

	RigBase* getRigObj(){return m_rigObj;}

	// 返回封装的状态
	RigStatus getStatus()const;
	bool setStatus(const RigStatus& s);

	const vector<int>& getInternalPntIdx()const{return m_intPntIdx;}
	const vector<int>& getSurfacePntIdx()const{return m_surfPntIdx;}
	void getMeshPntPos(vector<double>& pnts)const;
	int getNTotPnt()const{return m_nTotPnt;}
	int getNSurfPnt()const{return m_nSurfPnt;}
	int getNIntPnt()const{return m_nIntPnt;}

	// 各种测试函数，调试专用
	// 给定当前的配置n,p,检查Hessian是否正确逼近
	bool testHessian();
	// 测试当前能量函数对位移的梯度是否等于内力
	bool testElasticEnergy();
	// 测试函数值计算是否正确
	bool testValue();
	void testStep(double dt);

	// 以下函数被maya mel命令调用
	// 测试当前帧的Hessian
	bool testCurFrameHessian( RigStatus& lastFrame, RigStatus& curFrame, double noiseN = 1.0, double noiseP=1.0);
	bool testCurFrameGrad( RigStatus& lastFrame, RigStatus& curFrame, double noiseN = 1.0, double noiseP = 1.0);


protected:
	struct PntPair
	{
		Vec3d m_pnt;
		int   m_idx;
		bool  operator<(const PntPair& other)const;
	};

	bool findOriPoints(tetgenio& in, tetgenio& out);
	bool buildVegaData(double E = 1e6, double nu = 0.45, double density = 1000);
	// 由内部点偏移n和参数p更新所有点的偏移q
	bool computeQ(const double* n, const double* p, double t, double* q);
	// 只更新被参数p影响的点的位置
	bool computeQ(const double* p, double t, double* q);
	// 计算表面点的偏移量
	bool computeSurfOffset(const double* p, double t, EigVec& s);

	double computeValue(const EigVec& q, const EigVec& p);

	bool getN(RigStatus& s, EigVec& n);

	// 计算施加到每个节点的控制力（参数数个分量）
	// 此函数用于显式控制，其根据给定参数值和参数速度直接计算出控制力
	bool computeControlForce( EigVec& controlForce, const EigVec& param, const EigVec& paramVelocity);
	// 计算施加到每个节点的控制力（参数数个分量），
	// 此函数用于隐式控制,其以当前的参数值为上一帧的参数值，以此计算参数速度
	bool computeControlForce( EigVec& generalForce, const EigVec& param );

	// Tengen生成的四面体网格数据,这些数据一旦生成，在模拟过程中除坐标值外不再修改
	tetgenio		m_out;
	vector<int>		m_surfPntIdx;						// 来自m_in的顶点（视为表面点）在m_out的索引
	vector<int>		m_intPntIdx;						// 新加入的顶点(视为内部点)在m_out的索引
	vector<int>		m_surfDofIdx;						// 表面点各个自由度（xyzxyzxyz）索引
	vector<int>		m_intDofIdx;						// 内部点各个自由度（xyzxyzxyz）索引
	int				m_nTotPnt;							// 总点数
	int				m_nIntPnt;							// 内部点（也就是自由运动的点）个数
	int				m_nSurfPnt;							// 表面点（也就是被参数控制的点）个数
	int				m_nParam;							// 自由参数个数
	TetMesh*		m_tetMesh;							// 四面体网格

	// 各种Vega数据结构
	ModelWrapper*			m_modelwrapper;
	ForceModel*				m_forceModel;

	// rig数据
	RigBase*				m_rigObj;					// rig 数据

	// 以下状态变量长度为点数*3, 
	// 在迭代过程中，这些值保持为上一帧的值，直到新一帧的值计算出来后才更新
	EigVec					m_q,m_v,m_a;				// 每个顶点的偏移量,速度，加速度
	EigVec					m_force;					// 受的弹力
	EigVec					m_extForce;					// 受的外力
	EigVec					m_param;					// 所有rig 参数排成的向量
	EigVec					m_paramVelocity;			// 所有rig 参数速度
	SparseMatrix*			m_tangentStiffnessMatrix;

	EigVec					m_mass;						// 质量
	double					m_h;						// 时间步长
	double					m_t;						// 当前时刻

	// 控制参数
	RigControlType			m_controlType;				// 是否引入控制力
	EigVec					m_controlForce;				// 广义控制力，长度为参数数*3,
														// 在显式控制的情况下，控制力在求解进行前由当前状态给定
														// 在隐式控制的情况下，控制力依赖于解，在求解过程中不断变化
	// 在迭代过程中，这些值保持为上一帧的值，直到新一帧的值计算出来后才更新
	EigVec					m_targetParam;				// 目标参数值
	EigVec					m_targetParamVelocity;		// 目标参数速度
	EigVec					m_propGain;					// 比例控制增益
	EigVec					m_deriGain;					// 微分控制增益
};



class FEMSystem
{
public:
	void init();
	void show(){m_mesh.show();}
	void saveResult(const char* fileName);
	void clearResult();
	void step();
	RiggedMesh				m_mesh;						// 体网格对象
	PointParamSolver			m_solver;					// 牛顿法求解器
};

}