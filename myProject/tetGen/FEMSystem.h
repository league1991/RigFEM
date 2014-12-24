#pragma once

using namespace std;
//#define SHOW_TETGEN_RESULT		// 是否显示tetgen生成的out数据
#define SHOW_TETMESH_RESULT			// 是否显示Vega的体网格数据
namespace RigFEM
{

class RiggedMesh :public ObjectFunction
{
public:
	RiggedMesh(void);
	~RiggedMesh(void);

	void init();
	void show();
	void computeRig();

	// 计算给定状态x = [n，p] 以及时间参数param下的函数值，以及梯度
	bool computeValueAndGrad(const EigVec& x, const EigVec& param, double* v, EigVec* grad);

	// 根据给定的内部点的偏移量n以及参数值p计算函数值、梯度以及Hessian矩阵
	double computeValue(const EigVec& n, const EigVec& p, double t);
	bool computeGradient(const EigVec& n, const EigVec& p, double t, EigVec& gn, EigVec& gp);
	bool computeHessian(const EigVec& n, const EigVec& p, double t, EigSparse& Hnn, EigDense& Hnp, EigDense& Hpn, EigDense* pHpp = NULL);		// pHpp 若为NULL，不计算

	// 各种状态变量
	// 获得当前各个自由度的状态
	void getDof(EigVec& n, EigVec& p);
	// 更新各个自由度的状态，同时更新的有顶点位置、速度、加速度
	void setDof(EigVec&n, EigVec&p, bool proceedTime = true);
	double getCurTime(){return m_t;}

	// 各种测试函数，调试专用
	// 给定当前的配置n,p,检查Hessian是否正确逼近
	bool testHessian();
	// 测试当前能量函数对位移的梯度是否等于内力
	bool testElasticEnergy();
	// 测试函数值计算是否正确
	bool testValue();
	void testStep(double dt);
	
	// rig数据
	TransformRig	m_transRig;

private:
	struct PntPair
	{
		Vec3d m_pnt;
		int   m_idx;
		bool  operator<(const PntPair& other)const;
	};

	void findOriPoints();
	bool buildTetMesh();

	// 由内部点偏移n和参数p更新所有点的偏移q
	void computeQ(const double* n, const double* p, double t, double* q);
	// 只更新被参数p影响的点的位置
	void computeQ(const double* p, double t, double* q);
	// 计算表面点的偏移量
	void computeSurfOffset(const double* p, double t, EigVec& s);

	double computeValue(const EigVec& q);

	bool computeHnn(vector<int>& nIdx, vector<int>& sIdx)
	{
		if (!m_tangentStiffnessMatrix)
			return false;
		/*
		if (!m_HnnCache)
		{
			m_HnnCache = new SparseMatrix();
			SparseMatrix& Hnn = *m_HnnCache;

			SparseMatrix dFnn(*m_tangentStiffnessMatrix);
			dFnn.RemoveRowsColumns(sIdx.size(), &sIdx[0]);
			Utilities::vegaSparse2Eigen(dFnn, Hnn, m_nIntPnt*3);
			Hnn *= -1.0;

			for (int ithInt = 0; ithInt < m_nIntPnt; ++ithInt)
			{
				int begIdx = ithInt*3;
				double massTerm = m_mass[m_intPntIdx[ithInt]] / (m_h * m_h);
				Hnn.coeffRef(begIdx, begIdx) += massTerm;	begIdx++;
				Hnn.coeffRef(begIdx, begIdx) += massTerm;	begIdx++;
				Hnn.coeffRef(begIdx, begIdx) += massTerm;
			}
		}
		else
		{

		}*/
		return true;
	}


	// Tengen生成的四面体网格数据,这些数据一旦生成，在模拟过程中除坐标值外不再修改
	tetgenio		m_in, m_addin, m_bgmin, m_out;
	vector<int>		m_surfPntIdx;						// 来自m_in的顶点（视为表面点）在m_out的索引
	vector<int>		m_intPntIdx;						// 新加入的顶点(视为内部点)在m_out的索引
	int				m_nTotPnt;							// 总点数
	int				m_nIntPnt;							// 内部点（也就是自由运动的点）个数
	int				m_nSurfPnt;							// 表面点（也就是被参数控制的点）个数
	int				m_nParam;							// 自由参数个数
	TetMesh*		m_tetMesh;							// 四面体网格

	// 各种Vega数据结构
	ModelWrapper*			m_modelwrapper;
	ForceModel*				m_forceModel;

	// 以下状态变量长度为点数*3, 
	// 在迭代过程中，这些值保持为上一帧的值，直到新一帧的值计算出来后才更新
	EigVec					m_q,m_v,m_a;				// 每个顶点的偏移量,速度，加速度
	EigVec					m_force;					// 受力
	EigVec					m_param;					// 所有rig 参数排成的向量
	SparseMatrix*			m_tangentStiffnessMatrix;

	// 以下是暂存的Hessian各个分块
	SparseMatrix*			m_HnnCache;					// 暂存的Hnn
	

	EigVec					m_mass;						// 质量
	double					m_h;						// 时间步长
	double					m_t;						// 当前时刻

	// 记录模拟出来的参数变化
	vector<EigVec>			m_paramResult;				// 参数向量，包括关键帧驱动的参数，和模拟出来的参数
};

class FEMSystem
{
public:
	void init(){m_mesh.init();m_solver.setMesh(&m_mesh);}
	void show(){m_mesh.show();}
	void saveResult(const char* fileName);
	void clearResult();
	void step();
	RiggedMesh				m_mesh;						// 体网格对象
	NewtonSolver			m_solver;					// 牛顿法求解器
};

}