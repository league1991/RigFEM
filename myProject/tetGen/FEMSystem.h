#pragma once

using namespace std;
//#define SHOW_TETGEN_RESULT		// �Ƿ���ʾtetgen���ɵ�out����
#define SHOW_TETMESH_RESULT			// �Ƿ���ʾVega������������
namespace RigFEM
{

class RiggedMesh :public ObjectFunction
{
public:
	RiggedMesh(void);
	~RiggedMesh(void);
 
	void init();
	// �ñ��������rig�����ʼ����rig���ⲿ�����ͷ�
	bool init(tetgenio& surfMesh, RigBase* rig, double maxVolume = 1, double edgeRatio = 2);
	void clear();
	void show();
	bool showStatus(RigStatus& s);
	void computeRig();

	// �������״̬x = [n��p] �Լ�ʱ�����param�µĺ���ֵ���Լ��ݶ�
	bool computeValueAndGrad(const EigVec& x, const EigVec& param, double* v, EigVec* grad);

	// ���ݸ������ڲ����ƫ����n�Լ�����ֵp���㺯��ֵ���ݶ��Լ�Hessian����
	double computeValue(const EigVec& n, const EigVec& p, double t);
	bool computeGradient(const EigVec& n, const EigVec& p, double t, EigVec& gn, EigVec& gp);
	bool computeHessian(const EigVec& n, const EigVec& p, double t, EigSparse& Hnn, EigDense& Hnp, EigDense& Hpn, EigDense* pHpp = NULL);		// pHpp ��ΪNULL��������

	// ����״̬����
	// ��õ�ǰ�������ɶȵ�״̬
	void getDof(EigVec& n, EigVec& p);
	// ���¸������ɶȵ�״̬��ͬʱ���µ��ж���λ�á��ٶȡ����ٶ�
	void setDof(EigVec&n, EigVec&p, bool proceedTime = true);
	double getCurTime(){return m_t;}

	// ���ط�װ��״̬
	RigStatus getStatus()const;
	bool setStatus(const RigStatus& s);

	int getNTotPnt()const{return m_nTotPnt;}

	// ���ֲ��Ժ���������ר��
	// ������ǰ������n,p,���Hessian�Ƿ���ȷ�ƽ�
	bool testHessian();
	// ���Ե�ǰ����������λ�Ƶ��ݶ��Ƿ��������
	bool testElasticEnergy();
	// ���Ժ���ֵ�����Ƿ���ȷ
	bool testValue();
	void testStep(double dt);
	
	// rig����
	RigBase*	m_transRig;

private:
	struct PntPair
	{
		Vec3d m_pnt;
		int   m_idx;
		bool  operator<(const PntPair& other)const;
	};

	bool findOriPoints(tetgenio& in, tetgenio& out);
	bool buildVegaData();

	// ���ڲ���ƫ��n�Ͳ���p�������е��ƫ��q
	bool computeQ(const double* n, const double* p, double t, double* q);
	// ֻ���±�����pӰ��ĵ��λ��
	bool computeQ(const double* p, double t, double* q);
	// ���������ƫ����
	bool computeSurfOffset(const double* p, double t, EigVec& s);

	double computeValue(const EigVec& q);

	bool computeHnn(vector<int>& nIdx, vector<int>& sIdx)
	{
		if (!m_tangentStiffnessMatrix)
			return false;
		/*
		if (!m_HnnCache)
		{
			m_HnnCache = new SparseMatrix();
			SparseMatrix& Hnn = *m_HnnCache;

			SparseMatrix dFnn(*m_tangentStiffnessMatrix);
			dFnn.RemoveRowsColumns(sIdx.size(), &sIdx[0]);
			Utilities::vegaSparse2Eigen(dFnn, Hnn, m_nIntPnt*3);
			Hnn *= -1.0;

			for (int ithInt = 0; ithInt < m_nIntPnt; ++ithInt)
			{
				int begIdx = ithInt*3;
				double massTerm = m_mass[m_intPntIdx[ithInt]] / (m_h * m_h);
				Hnn.coeffRef(begIdx, begIdx) += massTerm;	begIdx++;
				Hnn.coeffRef(begIdx, begIdx) += massTerm;	begIdx++;
				Hnn.coeffRef(begIdx, begIdx) += massTerm;
			}
		}
		else
		{

		}*/
		return true;
	}


	// Tengen���ɵ���������������,��Щ����һ�����ɣ���ģ������г�����ֵ�ⲻ���޸�
	tetgenio		m_out;
	vector<int>		m_surfPntIdx;						// ����m_in�Ķ��㣨��Ϊ����㣩��m_out������
	vector<int>		m_intPntIdx;						// �¼���Ķ���(��Ϊ�ڲ���)��m_out������
	int				m_nTotPnt;							// �ܵ���
	int				m_nIntPnt;							// �ڲ��㣨Ҳ���������˶��ĵ㣩����
	int				m_nSurfPnt;							// ����㣨Ҳ���Ǳ��������Ƶĵ㣩����
	int				m_nParam;							// ���ɲ�������
	TetMesh*		m_tetMesh;							// ����������

	// ����Vega���ݽṹ
	ModelWrapper*			m_modelwrapper;
	ForceModel*				m_forceModel;

	// ����״̬��������Ϊ����*3, 
	// �ڵ��������У���Щֵ����Ϊ��һ֡��ֵ��ֱ����һ֡��ֵ���������Ÿ���
	EigVec					m_q,m_v,m_a;				// ÿ�������ƫ����,�ٶȣ����ٶ�
	EigVec					m_force;					// ����
	EigVec					m_param;					// ����rig �����ųɵ�����
	SparseMatrix*			m_tangentStiffnessMatrix;

	EigVec					m_mass;						// ����
	double					m_h;						// ʱ�䲽��
	double					m_t;						// ��ǰʱ��
};

class FEMSystem
{
public:
	void init(){m_mesh.init();m_solver.setMesh(&m_mesh);}
	void show(){m_mesh.show();}
	void saveResult(const char* fileName);
	void clearResult();
	void step();
	RiggedMesh				m_mesh;						// ���������
	NewtonSolver			m_solver;					// ţ�ٷ������
};

}