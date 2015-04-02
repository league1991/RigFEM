#pragma once
using namespace std;
namespace RigFEM
{
	// rig 黑盒
	class RigBase
	{
	public:
		typedef double (*KeyFrameFunc)(double);
		RigBase(int nParam);
		virtual ~RigBase(void);

		// 顶点,
		virtual int  getNPoints()const=0;							// 返回顶点个数

		int  getNParam()const{return m_nParam;}						// 返回总的参数个数
		void getParam(double* params);

		// 关键帧驱动的参数
		void keyParam(int ithParam, KeyFrameFunc func);
		void unKeyParam(int ithParam);

		// 自由参数
		virtual int  getNFreeParam()const;
		void setFreeParam(const double* params);
		void getFreeParam(double* params)const;;

		// 时间
		void setTime(double t);
		double		 getTime()const{return m_t;}

		// 设置有限差商求导时的步长
		void setDelta(double delta){m_delta = delta;}

		// 以下是各个计算函数, 根据当前的参数和时间计算各种结果
		// 计算当前参数下点的新坐标
		virtual bool computeValue(double* result, const double* params = 0) = 0;

		// 有限差商计算当前参数下的雅可比矩阵 dS/dP 
		virtual bool computeJacobian(Eigen::MatrixXd& jacobian);

		// 有限差商计算当前参数下的 d2S / (dPi * dPj)  其中S为点向量，Pi、Pj为第i、j个参数
		virtual bool computeJacobianDerivative(int i, int j, double* res);

		// 计算外力
		virtual bool computeExternalForce(const EigVec& pos, const EigVec& vel, const EigVec& m, 
			double time, EigVec& extForce, EigVec& surfForce){return false;}
		
		virtual bool getControlGain(EigVec& propGain, EigVec& deriGain){return false;}
		virtual bool getControlParams(EigVec& targetParam, EigVec& propGain, EigVec& deriGain){return false;}
		virtual bool getControlTarget(EigVec& targetParam){return false;}
	protected:
		double	m_delta;											// 有限差商求导时的微小位移
		int     m_nParam;											// 参数个数
		double  m_t;												// 当前时间
		double* m_param;											// 参数向量
		KeyFrameFunc* m_keyFrameFunc;								// 关键帧向量
	};

	class TransformRig:public RigBase
	{
	public:
		TransformRig();
		void setInitPnts(double* pnts, int nPnts);

		int  getNPoints()const{return m_initPntList.size();}

		Vec3d getTranslation()const{return m_translation;}
		Vec3d getScale()const{return m_scale;}
		Vec3d getRotation()const{return m_rotate;}
		void setAllParam( const Vec3d& trans, const Vec3d& rotate, const Vec3d& scale );
		const vector<Vec3d>& getCurPnt()const{return m_curPntList;}
		const vector<Vec3d>& getInitPnt()const{return m_initPntList;}

		virtual bool computeValue(double* result, const double* params = 0);
		//void computeJacobian(Eigen::MatrixXd& jacobian);
		//void computeJacobianDerivative(int i, int j, double* res);

		static inline Vec3d rotateX(const Vec3d& v, double deg)
		{
			Vec3d res;
			res[0] = v[0];
			double c = cos(deg), s = sin(deg);
			res[1] = c * v[1] - s * v[2];
			res[2] = s * v[1] + c * v[2];
			return res;
		}
		static inline Vec3d rotateY(const Vec3d& v, double deg)
		{
			Vec3d res;
			res[1] = v[1];
			double c = cos(deg), s = sin(deg);
			res[2] = c * v[2] - s * v[0];
			res[0] = s * v[2] + c * v[0];
			return res;
		}
		static inline Vec3d rotateZ(const Vec3d& v, double deg)
		{
			Vec3d res;
			res[2] = v[2];
			double c = cos(deg), s = sin(deg);
			res[0] = c * v[0] - s * v[1];
			res[1] = s * v[0] + c * v[1];
			return res;
		}
	private:
		// 根据当前的平移缩放参数计算变换后的点
		void transform();

		void compute();
		void compute(const Vec3d& trans, const Vec3d& scl);

		vector<Vec3d>	m_initPntList;		// 初始点位置（无平移，原大小）
		vector<Vec3d>   m_curPntList;		// 当前参数下点的位置

		// 当前rig 参数
		Vec3d			m_translation;
		Vec3d			m_scale;
		Vec3d			m_rotate;
		Vec3d			m_localCenter;		// 局部坐标系原点位置
	};
}

