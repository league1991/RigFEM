#pragma once

using namespace std;
namespace RigFEM
{	
	class RigStatus
	{
	public:
		enum Status
		{
			// 各种状态的枚举类型 
			STATUS_Q,
			STATUS_V,
			STATUS_A,
			STATUS_P,
			STATUS_F,
			STATUS_PV,
			STATUS_TAR_P
		};
		RigStatus(){}
		RigStatus(	const EigVec& q, 
					const EigVec& v,
					const EigVec& a,
					const EigVec& param,
					const EigVec& paramV,
					const EigVec& extF,
					const EigVec& targetParam);
		~RigStatus();

		// 返回点数向量长度，若q v a长度不一，返回-1
		int getPointVecLength()const;
		int getParamVecLength()const;

		const EigVec& getQ()const{return m_q;}
		const EigVec& getV()const{return m_v;}
		const EigVec& getA()const{return m_a;}
		const EigVec& getF()const{return m_f;}
		const EigVec& getP()const{return m_param;}
		const EigVec& getPV()const{return m_paramVelocity;}
		const EigVec& getTarP()const{return m_targetParam;}
		const EigVec* getByName(Status s)const;

		bool matchLength(int pntVecLength, int paramVecLength)const;

		// 自定义属性
		void addOrSetCustom(const string& name, double v);
		void addOrSetCustom(const string& name, const EigVec& v);
		void addOrSetCustom(const string& name, const EigDense& v);
		bool getCustom(const string& name, double& v)const;
		bool getCustom(const string& name, EigVec& v)const;
		bool getCustom(const string& name, EigDense& v)const;

		void mergeCustom(const RigStatus& s);

	private:
		EigVec	m_q, m_v, m_a, m_f;
		EigVec	m_param, m_paramVelocity;
		EigVec	m_targetParam;
		map<string, double> m_customScalar;			// 一些自定义参数
		map<string, EigVec> m_customVector;
		map<string, EigDense>m_customDenseMat;
	};
	class StatusRecorder
	{
	public:
		StatusRecorder(void): m_startFrame(0), m_pntLength(0), m_paramLength(0){}
		~StatusRecorder(void);

		void clear(){m_statusList.clear();}
		void init( int startFrame, int pntLength, int paramLength , 
			const vector<int>& intPntIdx, const vector<int>& surfPntIdx,
			const vector<double>& initPntPos);
		bool getStatus(int ithFrame, RigStatus& s)const;
		const RigStatus* getStatus(int ithFrame);
		bool setStatus(int ithFrame, const RigStatus& s);
		void setPntIdx(const vector<int>& intPntIdx, const vector<int>& surfPntIdx);

		int  getPointVecLength()const{return m_pntLength;}
		int  getParamVecLength()const{return m_paramLength;}
		int  getRecordFrameLength()const{return m_statusList.size();}

		bool saveToFile(const char* fileName)const;
		bool saveCustomToFile(const char* customParamName, const char* fileName)const;
		bool addCustomToFile( const char* customParamName, ofstream& file ) const;
	private:
		void state2Str(RigStatus::Status s, const char* name, string& str)const;
		void customMat2Str(const char* stateName, const char* matlabVarName, string& str)const;

		int				  m_pntLength;
		int				  m_paramLength;
		int				  m_startFrame;
		vector<RigStatus> m_statusList;
		vector<int>		  m_intPntIdx, m_surfPntIdx;				// 索引从1开始
		vector<double>	  m_initPntPos;
	};
}
