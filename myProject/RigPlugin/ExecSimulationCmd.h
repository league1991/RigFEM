#pragma once
/***********************************************************************
   MEL 命令: rigSimulate
   -init			-i	初始化模拟数据
   -step			-s	进行一步模拟，返回1代表正常，0代表有错误
   -stepStatic		-ss 进行静态模拟，返回1代表正常，0代表有错误
   -stepEleGF		-rgf进行模拟，并记录每个元素的广义力
   -saveGF			-sgf 保持广义力
   -loadEleMat		-lm 加载元素硬度
   -resetEleMat	    -rsm重设默认硬度
   -name			-n  指定节点名称
   -hessian			-h  测试Hessian矩阵，指定参数noiseN, noiseP
   -grad			-g  测试梯度是否与函数值增量吻合
   -save			-sa  保存模拟结果
   -numInternalPnt	-nip	获得内部点数目
   -numSurfacePnt	-nsp    获得表面点数目
   -create			-c		创建模拟节点
/************************************************************************/

class RigSimulateCmd: public MPxCommand
{
public:
	RigSimulateCmd(void);
	~RigSimulateCmd(void);

	MStatus doIt(const MArgList& args);

	static void* creator();

	static MSyntax newSyntax();
private:
	static const char*				m_initFlag[2];
	static const char*				m_stepFlag[2];
	static const char*				m_stepStaticFlag[2];
	static const char*				m_recordEleGFFlag[2];
	static const char*				m_saveFlag[2];
	static const char*				m_saveEleGFFlag[2];
	static const char*				m_loadEleMatFlag[2];
	static const char*				m_resetEleMatFlag[2];

	static const char*				m_nameFlag[2];
	static const char*				m_hessianFlag[2];
	static const char*				m_gradFlag[2];
	static const char*				m_intPntFlag[2];
	static const char*				m_surfPntFlag[2];
	static const char*				m_createFlag[2];
};
