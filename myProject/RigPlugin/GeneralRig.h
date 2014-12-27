#pragma once

class RigSimulationNode;
namespace RigFEM
{

	class RigException
	{
	public:
		enum Type
		{
			MESH_TOPOLOGY_CHANGED = 0x1,	// 网格拓扑发生改变
		};
		RigException(Type type, const MString& info = ""): m_type(type), m_info(info){}
		Type getType()const{return m_type;}
		const MString& getInfo()const{return m_info;}
	private:
		Type m_type;
		MString m_info;
	};

	class GeneralRig : public RigBase
	{
	public:
		GeneralRig(RigSimulationNode* node, int nParam, int nPnts);
		~GeneralRig(void);

		virtual int getNPoints()const{return m_nPnts;}
		virtual int getNFreeParam()const{return m_nParam;}

		// 计算函数
		virtual bool computeValue(double* result, const double* params = 0);
		
		// 从节点获取参数值，一般用于初始化节点
		MStatus		fetchParamFromNode();

		MStatus		getMesh(MObject& meshObj);
	private:
		// 禁止关键帧驱动
		void keyParam(int ithParam, KeyFrameFunc func);
		void unKeyParam(int ithParam);

		// 获取网格顶点,数组预先分配好空间
		MStatus getMeshPoints(double* points);
		// 设置参数值
		MStatus setParamPlug();

		RigSimulationNode*			m_node;
		int							m_nPnts;
	};
}
