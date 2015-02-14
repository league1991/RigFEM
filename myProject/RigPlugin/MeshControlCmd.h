#pragma once

class MeshControlCmd :
	public MPxCommand
{
public:
	MeshControlCmd(void);
	~MeshControlCmd(void);

	MStatus doIt(const MArgList& args);

	static void* creator(){return new MeshControlCmd;}

	static MSyntax newSyntax();

private:
	static const char*		m_initFlag[2];
	static const char*      m_nameFlag[2];
};
