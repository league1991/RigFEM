#pragma once

using namespace std;

namespace RigFEM
{
	class RiggedSkinMesh:public RiggedMesh
	{
	public:
		void setWeight(EigSparse& sparse);
		// 调用了init之后再被调用
		bool setWeight(const char* weightFile);

		void getDof( EigVec& p );
		void setDof( EigVec& p, bool proceedTime = true );

		// 各种计算函数
		// 计算给定状态x = p 以及时间参数param下的函数值，以及梯度
		bool computeValueAndGrad(const EigVec& x, const EigVec& param, double* v = NULL, EigVec* grad = NULL);
		bool computeHessian(const EigVec&x, const EigVec& param, EigDense& H);
		// 计算近似的Hessian矩阵
		// 假定雅可比矩阵为常数
		bool computeApproxHessian(const EigVec&x, const EigVec& param, EigDense& H);

		bool testCurFrameGrad(RigStatus& lastFrame, RigStatus& curFrame, double noiseP = 1.0);
		bool testCurFrameHessian( RigStatus& lastFrame, RigStatus& curFrame, double noiseP=1.0);
	protected:
		bool computeSkinQ(const double* p, double t, double* q);

		EigSparse				m_weightMat;				// 权重矩阵，大小为（内部点数*3, 表面点数*3）
		EigSparse				m_weightMatTran;			// 权重矩阵的转置
	};
}