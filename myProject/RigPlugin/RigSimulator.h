#pragma once

class RigSimulationNode;
namespace RigFEM
{
	class SimulatorBase
	{
	public:
		SimulatorBase(){}
		virtual ~SimulatorBase(){}

		virtual bool testGradient(int curFrame)=0;
		virtual bool testHessian(int curFrame)=0;
		virtual bool saveResult(const char* fileName)=0;
		virtual bool stepRig(int curFrame)=0;
		virtual bool staticStepRig(int curFrame, const EigVec& curParam)=0;
		virtual bool showStatus(int curFrame, double* bbox = 0)=0;
		virtual bool getParam(int curFrame, EigVec& curParam)=0;
		virtual bool setDeriStepSize(double step)=0;
		virtual bool setStaticSolveMaxIter(int maxIter)=0;
		virtual bool isReady()const=0;
	};
	class GeneralRig;
	class RigSimulator:public SimulatorBase
	{
	public:
		RigSimulator():m_rigMesh(NULL), m_rig(NULL), m_solver(NULL), m_recorder(NULL){};
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
		RiggedMesh*		getRiggedMesh(){return m_rigMesh;}
		NewtonSolver*   getSolver(){return m_solver;}
		bool			testHessian(int curFrame);
		bool			testGradient(int curFrame);
		bool			saveResult(const char* fileName);
		bool			stepRig(int curFrame);
		bool			staticStepRig(int curFrame, const EigVec& curParam);
		bool			setStaticSolveMaxIter(int maxIter);
		bool			setDeriStepSize(double step);
		bool			showStatus(int curFrame, double* bbox = 0);
		bool			getParam(int curFrame, EigVec& curParam);
	protected:
		bool			getInitStatus(RigStatus& status);
		virtual void	allocateSimObj();
		virtual void	freeSimObj();
		
		int				m_maxStaticSolveIter;
		EigVec			m_initParam;
		GeneralRig*		m_rig;
		RiggedMesh*		m_rigMesh;
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
