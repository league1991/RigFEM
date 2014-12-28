#pragma once
/***********************************************************************
   MEL 命令: rigSimulate
   -init			-i	初始化模拟数据
   -step			-s	进行一步模拟，返回1代表正常，0代表有错误
   -name			-n  指定节点名称
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
	static const char*				m_nameFlag[2];
};
