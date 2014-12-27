#include "StdAfx.h"
#include "StatusRecorder.h"

using namespace RigFEM;



StatusRecorder::~StatusRecorder(void)
{
}

int RigFEM::RigStatus::getPointVecLength() const
{
	if (m_q.size() != m_v.size() || m_v.size() != m_a.size())
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
		(m_param.size() == paramVecLength);
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

bool RigFEM::StatusRecorder::appendStatus( const RigStatus& s )
{
	if (s.matchLength(m_pntLength, m_paramLength))
	{
		m_statusList.push_back(s);
		return true;
	}
	return false;
}

void RigFEM::StatusRecorder::init( int startFrame, int pntLength, int paramLength )
{
	m_pntLength = pntLength;
	m_startFrame = startFrame;
	m_paramLength = paramLength;
	m_statusList.clear();
}

bool RigFEM::StatusRecorder::setStatus( int ithFrame, const RigStatus& s )
{
	if (!s.matchLength(m_pntLength, m_paramLength))
		return false;

	int endFrame = m_startFrame + m_statusList.size() - 1;
	if (ithFrame >= m_startFrame && ithFrame <= endFrame)
	{
		m_statusList[ithFrame - m_startFrame] = s;
		return true;
	}
	if (ithFrame == endFrame+1)
	{
		m_statusList.push_back(s);
		return true;
	}
	return false;
}
