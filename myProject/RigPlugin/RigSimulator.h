#pragma once

class RigSimulationNode;
namespace RigFEM
{
	class SimulatorBase
	{
	public:
		SimulatorBase(): m_rigMesh(NULL){}
		virtual ~SimulatorBase(){}

		virtual bool testGradient(int curFrame, double noiseN, double noiseP )=0;
		virtual bool testHessian(int curFrame, double noiseN, double noiseP )=0;

		// 正常模拟各项功能
		virtual bool stepRig(int curFrame)=0;
		virtual bool staticStepRig(int curFrame, const EigVec& curParam)=0;
		virtual bool saveResult(const char* fileName)=0;

		// 修改四面体硬度各项功能
		virtual bool stepAndSaveEleGFRig(int curFrame, const EigVec& curParam){return false;}
		virtual bool saveGFResult(const char* fileName) = 0;
		virtual bool loadElementMaterialFactor(const char* fileName);
		virtual bool resetElementMaterialFactor();

		virtual bool showStatus(int curFrame, double* bbox = 0)=0;
		virtual bool getParam(int curFrame, EigVec& curParam)=0;
		virtual bool setDeriStepSize(double step)=0;
		virtual bool setStaticSolveMaxIter(int maxIter)=0;
		virtual bool isReady()const=0;

		RiggedMesh*		getRiggedMesh(){return m_rigMesh;}

		// 一些显示状态
		MeshDispConfig& getMeshConfig(){return m_dispConfig;}
		void			setMeshConfig(const MeshDispConfig& config){m_dispConfig = config;}
		void			setExternalForceDispFactor(double factor);
	protected:
		// 一些显示状态
		MeshDispConfig	m_dispConfig;
		
		// rig网格
		RiggedMesh*		m_rigMesh;

		static const char*s_materialName;
	};
	class GeneralRig;
	class RigSimulator:public SimulatorBase
	{
	public:
		RigSimulator(): m_rig(NULL), m_solver(NULL), m_recorder(NULL){};
		~RigSimulator(void);

		bool			init(	tetgenio& surfMesh, RigSimulationNode* node, 
								EigVec initParam,
								int    curFrame = 0,
								double maxVolume = 1, double edgeRatio = 2 , 
								double youngModulus = 1e6, double nu = 0.45, 
								double density = 1000,
								double timeStep = 1/24.0,
								int	   maxStaticSolveIter=20);

		bool			isReady()const;
		NewtonSolver*   getSolver(){return m_solver;}
		bool			getParam(int curFrame, EigVec& curParam);
		bool			testHessian(int curFrame, double noiseN, double noiseP );
		bool			testGradient(int curFrame, double noiseN, double noiseP );

		bool			stepRig(int curFrame);
		bool			staticStepRig(int curFrame, const EigVec& curParam);
		bool			stepAndSaveEleGFRig(int curFrame, const EigVec& curParam);

		bool			setStaticSolveMaxIter(int maxIter);
		bool			setDeriStepSize(double step);
		bool			showStatus( int curFrame, double* bbox = NULL);
		bool			saveResult(const char* fileName);
		bool			saveGFResult(const char* fileName);

		bool			setControlType(RigControlType type);
	protected:
		virtual void	allocateSimObj();
		virtual void	freeSimObj();
		
		GeneralRig*		m_rig;
		NewtonSolver*	m_solver;
		StatusRecorder*	m_recorder;
	};

	class RigSkinSimulator:public RigSimulator
	{
	public:

		bool			init(	tetgenio& surfMesh, RigSimulationNode* node, 
								EigVec initParam,
								const char*  weightFile,
								int    curFrame = 0,
								double maxVolume = 1, double edgeRatio = 2 , 
								double youngModulus = 1e6, double nu = 0.45, 
								double density = 1000,
								double timeStep = 1/24.0,
								int	   maxStaticSolveIter=20);
		bool			testHessian(int curFrame);
		bool			testGradient(int curFrame);

		bool			staticStepRig(int curFrame, const EigVec& curParam);
	protected:
		void			allocateSimObj();
	};
}
