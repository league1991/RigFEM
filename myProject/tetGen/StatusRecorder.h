#pragma once

namespace RigFEM
{
	
	class RigStatus
	{
	public:
		RigStatus(){}
		RigStatus(	const EigVec& q, 
					const EigVec& v,
					const EigVec& a,
					const EigVec& param):m_q(q), m_v(v), m_a(a), m_param(param)
		{}

		// 返回点数向量长度，若q v a长度不一，返回-1
		int getPointVecLength()const;
		int getParamVecLength()const;

		const EigVec& getQ()const{return m_q;}
		const EigVec& getV()const{return m_v;}
		const EigVec& getA()const{return m_a;}
		const EigVec& getParam()const{return m_param;}

		bool matchLength(int pntVecLength, int paramVecLength)const;
	private:
		EigVec	m_q, m_v, m_a;
		EigVec	m_param;
	};
	class StatusRecorder
	{
	public:
		StatusRecorder(void): m_startFrame(0), m_pntLength(0), m_paramLength(0){}
		~StatusRecorder(void);

		void clear(){m_statusList.clear();}
		void init(int startFrame, int pntLength, int paramLength);
		bool appendStatus(const RigStatus& s);
		bool getStatus(int ithFrame, RigStatus& s)const;
		bool setStatus(int ithFrame, const RigStatus& s);
		int  getPointVecLength()const{return m_pntLength;}
		int  getParamVecLength()const{return m_paramLength;}
	private:
		int				  m_pntLength;
		int				  m_paramLength;
		int				  m_startFrame;
		vector<RigStatus> m_statusList;
	};
}
