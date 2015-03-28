#include "StdAfx.h"
#include "StatusRecorder.h"

using namespace RigFEM;



StatusRecorder::~StatusRecorder(void)
{
}

int RigFEM::RigStatus::getPointVecLength() const
{
	if (m_q.size() != m_v.size() || m_v.size() != m_a.size() || m_a.size() != m_f.size())
	{
		return -1;
	}
	return m_q.size();
}

int RigFEM::RigStatus::getParamVecLength() const
{
	return m_param.size();
}

bool RigFEM::RigStatus::matchLength( int pntVecLength, int paramVecLength )const
{
	return (m_q.size() == pntVecLength) &&
		(m_v.size() == pntVecLength) &&
		(m_a.size() == pntVecLength) &&
		(m_param.size() == paramVecLength) &&
		(m_paramVelocity.size() == paramVecLength);
}

const EigVec* RigFEM::RigStatus::getByName( Status s ) const
{
	switch(s)
	{
	case STATUS_Q:
		return &m_q;
	case STATUS_V:
		return &m_v;
	case STATUS_A:
		return &m_a;
	case STATUS_P:
		return &m_param;
	case STATUS_F:
		return &m_f;
	case STATUS_PV:
		return &m_paramVelocity;
	case STATUS_TAR_P:
		return &m_targetParam;
	default:
		return NULL;
	}
}

RigFEM::RigStatus::~RigStatus()
{
}

RigFEM::RigStatus::RigStatus( 
	const EigVec& q, const EigVec& v, const EigVec& a, 
	const EigVec& param, const EigVec& paramVelocity, const EigVec& f,
	const EigVec& targetParam) :
m_q(q), m_v(v), m_a(a), m_param(param), m_f(f), m_paramVelocity(paramVelocity), m_targetParam(targetParam)
{
}


bool RigFEM::RigStatus::getCustom( const string& name, double& v ) const
{
	map<string, double>::const_iterator p = m_customScalar.find(name);
	if (p == m_customScalar.end())
	{
		return false;
	}
	v = p->second;
	return true;
}

bool RigFEM::RigStatus::getCustom( const string& name, EigVec& v ) const
{
	map<string, EigVec>::const_iterator p = m_customVector.find(name);
	if (p == m_customVector.end())
	{
		return false;
	}
	v = p->second;
	return true;
}

void RigFEM::RigStatus::addOrSetCustom( const string& name, double v )
{
	m_customScalar[name] = v;
}

void RigFEM::RigStatus::addOrSetCustom( const string& name, const EigVec& v )
{
	m_customVector[name] = v;
}

void RigFEM::RigStatus::mergeCustom( const RigStatus& s )
{
	for (map<string, double>::const_iterator it = s.m_customScalar.begin();
		it != s.m_customScalar.end(); ++it)
	{
		m_customScalar[it->first] = it->second;
	}
	for (map<string, EigVec>::const_iterator it = s.m_customVector.begin();
		it != s.m_customVector.end(); ++it)
	{
		m_customVector[it->first] = it->second;
	}
}

bool RigFEM::StatusRecorder::getStatus( int ithFrame, RigStatus& s ) const
{
	if (ithFrame < m_startFrame || 
		ithFrame >= m_startFrame + m_statusList.size())
	{
		return false;
	}
	s = m_statusList[ithFrame - m_startFrame];
	return true;
}

const RigStatus* RigFEM::StatusRecorder::getStatus( int ithFrame )
{
	if (ithFrame < m_startFrame || 
		ithFrame >= m_startFrame + m_statusList.size())
	{
		return NULL;
	}
	return &m_statusList[ithFrame - m_startFrame];
}



void RigFEM::StatusRecorder::init( int startFrame, int pntLength, int paramLength , 
								  const vector<int>& intPntIdx,
								  const vector<int>& surfPntIdx,
								  const vector<double>& initPntPos)
{
	m_pntLength = pntLength;
	m_startFrame = startFrame;
	m_paramLength = paramLength;
	setPntIdx(intPntIdx, surfPntIdx);
	m_initPntPos = initPntPos;
	m_statusList.clear();
}

bool RigFEM::StatusRecorder::setStatus( int ithFrame, const RigStatus& s )
{
	if (!s.matchLength(m_pntLength, m_paramLength))
		return false;

	int endFrame = m_startFrame + m_statusList.size() - 1;
	if (ithFrame >= m_startFrame && ithFrame <= endFrame)
	{
		int idx = ithFrame - m_startFrame;
		m_statusList[ithFrame - m_startFrame] = s;		
		return true;
	}
	if (ithFrame == endFrame+1)
	{
		int idx = ithFrame - m_startFrame;
		m_statusList.push_back(s);
		return true;
	}
	return false;
}

bool RigFEM::StatusRecorder::saveToFile( const char* fileName ) const
{
	ofstream file(fileName);
	if (!file)
	{
		return false;
	}

	string buffer;
	buffer = "\
% simulation result\n\
% q,v,a    (nFrame, nVtx*3)     offset,velocity,acceleration of every DOF, listed like xyzxyz...\n\
% param    (nFrame, nParam)     parameter value\n\
% intPntIdx(1, nIntVtx)			index of internal vertices(1-indexed)\n\
% surfPntIdx(1, nSurfVtx)       index of surface vertices(1-indexed)\n\
% initPos  (1, nVtx*3)			initial point position\n\
";

	file << buffer;
	state2Str(RigStatus::STATUS_Q, "q", buffer);
	file << buffer;
	state2Str(RigStatus::STATUS_V, "v", buffer);
	file << buffer;
	state2Str(RigStatus::STATUS_A, "a", buffer);
	file << buffer;
	state2Str(RigStatus::STATUS_P, "param", buffer);
	file << buffer;
	buffer = Utilities::vecToString(m_intPntIdx, "intPntIdx");
	file << buffer;
	buffer = Utilities::vecToString(m_surfPntIdx, "surfPntIdx");
	file << buffer;
	buffer = Utilities::vecToString(m_initPntPos, "initPos");
	file << buffer;

	file.close();
	return true;
}

void RigFEM::StatusRecorder::state2Str( RigStatus::Status s, const char* name, string& str ) const
{
	str = "";
	str += name;
	str += "=[\n";
	for (int ithFrame = 0; ithFrame < m_statusList.size(); ++ithFrame)
	{
		const RigStatus& status = m_statusList[ithFrame];
		if (const EigVec* vec = status.getByName(s))
		{
			char numberBuf[100];
			for (int i = 0; i < vec->size(); ++i)
			{
				double val = (*vec)[i];
				sprintf_s(numberBuf, 100, "%lf ", val);
				str += numberBuf;
			}
			str += "\n";
		}
	}
	str += "\n];\n";
}

void RigFEM::StatusRecorder::setPntIdx( const vector<int>& intPntIdx, const vector<int>& surfPntIdx )
{
	m_intPntIdx = intPntIdx;
	m_surfPntIdx = surfPntIdx;
	for (int i = 0; i < m_intPntIdx.size(); ++i)
	{
		m_intPntIdx[i]++;
	}
	for (int i = 0; i < m_surfPntIdx.size(); ++i)
	{
		m_surfPntIdx[i]++;
	}
}
