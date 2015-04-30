#pragma once
#include "StatusRecorder.h"

// 定义不同的共轭梯度算法
#define FR_CG		0			// fletcher-reeves  法
#define PR_CG		1			// polak-ribiere	法
#define CG_METHOD	FR_CG
namespace RigFEM
{
	// 进行一维搜索的函数对象
	class ObjectFunction
	{
	public:
		virtual ~ObjectFunction(){}
		virtual bool computeValueAndGrad(const EigVec& x, const EigVec& param, double* v = NULL, EigVec* grad = NULL) = 0;

	protected:
	};

	// 一维函数，用于测试
	class SimpleFunction:public ObjectFunction
	{
	public:
		virtual bool computeValueAndGrad(const EigVec& x, double* v = NULL, EigVec* grad = NULL);
	private:
		double computeFuncVal(double x);
	};
	
	class LineSearcher
	{
	public:
		LineSearcher(ObjectFunction* objFun = NULL, double initStep = 1.0, double maxStep = 10);
		~LineSearcher(void);
		enum Status
		{
			LS_NONE = 0,
			LS_WOLFE_DESCENT	= (0x1) << 1,
			LS_WOLFE_CURVATURE	= (0x1) << 2
		};
	
		bool setC(double c1, double c2);
		
		// 一维搜索函数
		// param 是目标函数的一些参数，例如时间。在搜索过程中参数保持不变
		// 若返回值等于0， 成功找到符合wolfe条件的函数
		// 若返回值不等于0， 没有找到符合wolfe条件的函数，可根据Status各个枚举变量获得问题所在
		// f0 df0为起始点的函数值和导数，若不提供，则另行计算
		int lineSearch( const EigVec& x0, const EigVec& dx, const EigVec& param, double& aFinal,
						double* f0 = NULL, double *df0 = NULL);
		void setInitStep(double initStep){m_initStep = initStep;}
		void setMaxStep(double maxStep){m_maxStep = maxStep;}
		void setMaxZoomIter(int maxZoomIter){m_maxZoomIter = maxZoomIter;}
	private:

		void setOriginAndDir(const EigVec& x0, const EigVec& dx);

		// 计算函数值和方向导数值，
		// 若f或df为NULL，不进行对应计算
		bool computeValueAndDeri( const EigVec& x, double* f, double* df );
		// 较准确的放大函数
		// 在al，ah确定的区间之间（al不一定小于ah）找一个符合wolfe条件的a,若找到，返回true
		// al 为一维搜索过程中，满足充分下降条件的点中函数值最小的一个
		// ah 使得 f'(al)(ah - al) < 0
		// dfal,dfah 为al ah点的导数，若当前没有，可以设为NULL
		// 若返回值等于0， 成功找到符合wolfe条件的函数
		// 若返回值等于1， 没有找到符合wolfe条件的函数
		int zoom(	double al, double fal, double* dfal, 
					double ah, double fah, double* dfah, 
					double& a );
		// 在区间[a0,a1]或[a1,a0]中估计一个极小值点
		// 若返回0，表示成功在区间找到一个极小值点
		// 若返回1，表示失败，此时返回中点
		int guessMinPnt( double a0, double fa0, double *dfa0, double a1, double fa1, double *dfa1, double& aMin );

		EigVec			m_param;			// 目标函数参数
		EigVec			m_dx;				// 搜索方向
		EigVec			m_x0;				// 搜索起始点
		double			m_f0;				// 搜索起始点的函数值
		double			m_df0;				// 搜索起始点函数的方向导数

		double			m_c1;				// wolfe充分下降条件系数
		double			m_c2;				// wolfe曲率条件系数

		double			m_initStep;			// 初始步长
		double			m_maxStep;			// 最大步长
		int				m_maxZoomIter;

		ObjectFunction* m_objFunc;			// 目标函数对象
	};

	class RiggedMesh;
	class RiggedSkinMesh;
	enum  RigControlType;

	class NewtonSolver
	{
	public:
		NewtonSolver(RiggedMesh* fem);
		virtual ~NewtonSolver(){}

		// 设置不发生变形时的参数值
		void setRestParam(EigVec& restParam);
		// 终止条件
		void setTerminateCond( int maxIter, double minStepSize, double minGradSize, int maxCGIter, double minCGStepSize);

		bool setInitStatus(const RigStatus&s);
		const RigStatus& getFinalStatus()const{return m_finalStatus;}
		virtual bool step()=0;
		virtual bool staticSolve(const EigVec& curParam){return false;}
		virtual bool staticSolveWithEleGF(const EigVec& curParam){return false;}

		void setIterationMaxStepSize(double maxStep);
		
		virtual void setControlType(RigControlType type){m_controlType = type;}

		StatusRecorder& getRecorder(){return m_recorder;}

		void setCurrentFrame(int curFrame){m_curFrame = curFrame;}
		int  getCurrentFrame(){return m_curFrame;}

		void setStaticSolveMaxIter(int maxIter);

		bool getRestStatus( RigStatus& status );

		// 额外传递的状态记录名
		static const char*  const	s_dPName;					// 参数前后帧的值增量
		static const char*	const	s_initStepName;				// 迭代初始步长
		static const char*	const	s_reducedElementGFName;			// 每个元素对参数空间产生的广义力
	protected:
		// 获得当前帧的目标控制参数和目标控制参数的速度
		bool getCurCtrlParam(EigVec& tarParam, EigVec& tarParamVelocity);

		EigVec			m_restParam;

		RiggedMesh*		m_femBase;

		RigStatus	m_initStatus;		// 模拟初始状态
		RigStatus	m_finalStatus;		// 模拟结束状态
		StatusRecorder	m_recorder;		// 状态记录器

		int			m_curFrame;			// 当前帧号

		RigControlType m_controlType;	// 是否加入控制

		// 以下是迭代的终止条件
		int			m_maxIter;			// 最大迭代次数			
		double		m_minStepSize;		// 最小步长，步长长度小于此值时迭代终止
		double		m_minGradSize;		// 梯度长度小于此值时迭代终止
		int			m_maxCGIter;		// 共轭梯度法最大迭代次数
		double      m_minCGStepSize;	// 共轭梯度法中，梯度长度小于此值时迭代终止

		double		m_iterMaxStepSize;	// 迭代过程中最大步长,每次参数各个分量的增量不得超过此值

		int			m_maxStaticSolveIter;

	};
	class PointParamSolver:public NewtonSolver
	{
	public:
		PointParamSolver(RiggedMesh* fem = NULL);
		~PointParamSolver();

		void setMesh(RiggedMesh* fem);
		// 计算函数
		bool step();
		bool staticSolve(const EigVec& curParam);
		// 静态求解,并算出各个有限元的广义力记录下来
		bool staticSolveWithEleGF(const EigVec& curParam);

		void clearResult(){m_paramResult.clear();}
		void saveResult(const char* fileName, const char* paramName = "param");

		virtual void setControlType(RigControlType type);
	private:
		// 给定参数值，利用力平衡条件计算出内部点位置，
		// 参数值发生变化时，把表面点，内部点都看成参数的函数
		// 计算出这个向量值函数的雅可比矩阵
		bool computeStaticJacobian(const EigVec& curParam, EigDense& J);

		RiggedMesh*		m_fem;
		LineSearcher	m_lineSearch;

		// 记录模拟出来的参数变化
		vector<EigVec>			m_paramResult;				// 参数向量，包括关键帧驱动的参数，和模拟出来的参数
	};

	class ParamSolver:public NewtonSolver
	{
	public:
		enum HessianType
		{
			HESSIAN_TRUE,
			HESSIAN_CONST_JACOBIAN
		};
		ParamSolver(RiggedSkinMesh* fem = NULL, HessianType type = HESSIAN_CONST_JACOBIAN);

		bool step();
		bool staticSolveWithEleGF(const EigVec& curParam);

		virtual void setControlType(RigControlType type);
	private:
		RiggedSkinMesh*	m_fem;
		LineSearcher	m_lineSearch;


		// 记录模拟出来的参数变化
		vector<EigVec>			m_paramResult;				// 参数向量，包括关键帧驱动的参数，和模拟出来的参数
		HessianType				m_hessianType;
	};
}
