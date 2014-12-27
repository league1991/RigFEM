#pragma once

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
	
	class LineSearch
	{
	public:
		LineSearch(ObjectFunction* objFun = NULL, double initStep = 1.0, double maxStep = 10);
		~LineSearch(void);

		
		// 一维搜索函数
		// param 是目标函数的一些参数，例如时间。在搜索过程中参数保持不变
		// 若返回值等于0， 成功找到符合wolfe条件的函数
		// 若返回值等于1， 没有找到符合wolfe条件的函数
		// f0 df0为起始点的函数值和导数，若不提供，则另行计算
		int lineSearch( const EigVec& x0, const EigVec& dx, const EigVec& param, double& aFinal,
						double* f0 = NULL, double *df0 = NULL);
		void setInitStep(double initStep){m_initStep = initStep;}
		void setMaxStep(double maxStep){m_maxStep = maxStep;}
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
		// 若返回1，表示取a0为最小点
		// 若返回2，表示取a1为最小点
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

		ObjectFunction* m_objFunc;			// 目标函数对象
	};

	class RiggedMesh;
	class NewtonSolver
	{
	public:
		NewtonSolver(RiggedMesh* fem = NULL);
		void setMesh(RiggedMesh* fem);
		// 计算函数
		bool step();

		void clearResult(){m_paramResult.clear();}
		void saveResult(const char* fileName, const char* paramName = "param");

	private:
		RiggedMesh*	m_fem;
		LineSearch	m_lineSearch;
		int			m_maxIter;

		// 记录模拟出来的参数变化
		vector<EigVec>			m_paramResult;				// 参数向量，包括关键帧驱动的参数，和模拟出来的参数
	};
}
