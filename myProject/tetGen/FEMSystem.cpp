#include "StdAfx.h"
#include "FEMSystem.h"

using namespace RigFEM;

unsigned char MeshDispConfig::s_colorRuler[] = 
{
	0,   64,  143,
	2,   66,  144,
	6,   68,  145,
	9,   71,  147,
	13,   74,  149,
	17,   77,  151,
	20,   80,  154,
	25,   83,  156,
	29,   88,  158,
	33,   92,  162,
	39,   95,  164,
	44,   99,  167,
	50,  104,  171,
	55,  108,  173,
	60,  113,  177,
	66,  118,  180,
	72,  123,  183,
	77,  128,  187,
	83,  133,  191,
	90,  137,  194,
	96,  143,  197,
	102,  148,  201,
	109,  153,  204,
	115,  159,  208,
	121,  163,  211,
	128,  169,  215,
	134,  173,  219,
	140,  179,  222,
	147,  184,  225,
	153,  189,  228,
	159,  194,  231,
	165,  198,  234,
	172,  204,  237,
	178,  208,  240,
	183,  212,  242,
	189,  216,  245,
	195,  221,  247,
	201,  225,  249,
	205,  228,  252,
	211,  232,  253,
	216,  235,  254,
	221,  239,  255,
	225,  241,  255,
	230,  244,  255,
	235,  246,  255,
	238,  249,  255,
	243,  251,  255,
	246,  252,  255,
	249,  253,  255,
	253,  255,  255,
	255,  255,  255,
	255,  255,  253,
	255,  255,  251,
	255,  255,  249,
	255,  255,  246,
	255,  253,  243,
	255,  252,  240,
	255,  251,  236,
	255,  249,  232,
	255,  247,  228,
	255,  245,  224,
	255,  243,  219,
	255,  240,  213,
	255,  238,  208,
	255,  234,  203,
	255,  231,  198,
	255,  228,  191,
	255,  224,  185,
	255,  221,  180,
	255,  218,  174,
	255,  213,  167,
	255,  209,  160,
	255,  206,  154,
	255,  202,  148,
	255,  197,  141,
	255,  194,  135,
	255,  189,  128,
	255,  185,  122,
	255,  181,  115,
	255,  177,  108,
	255,  173,  101,
	255,  168,   94,
	255,  164,   89,
	255,  160,   82,
	255,  155,   75,
	255,  152,   69,
	255,  148,   63,
	255,  144,   57,
	255,  141,   52,
	255,  137,   46,
	255,  133,   40,
	255,  130,   36,
	255,  127,   30,
	255,  124,   26,
	255,  121,   21,
	255,  118,   17,
	255,  116,   13,
	255,  114,   10,
	255,  111,   6,
	255,  110,   3
};

RiggedMesh::RiggedMesh(void):
m_tetMesh(NULL),  m_forceModel(NULL), 
m_h(0.03), m_t(0), m_tangentStiffnessMatrix(NULL),  
m_modelwrapper(NULL), 
m_nTotPnt(0), m_nIntPnt(0), m_nSurfPnt(0), m_nIntDof(0), m_nSurfDof(0), m_nParam(0),
m_rigObj(NULL), m_controlType(CONTROL_NONE),
m_surfPntIdx(NULL), m_intPntIdx(NULL), m_surfDofIdx(NULL), m_intDofIdx(NULL)
{
}

RiggedMesh::~RiggedMesh(void)
{
	clear();
}

void RiggedMesh::init()
{
	// 先释放内存
	clear();

	char*  myArgv[] = {
		"i:/Programs/VegaFEM-v2.1/myProject/tetGen/Debug/tetGen.exe",
		"-pq1.3a0.5mR",
		"model/torus.off"
	};
	int myArgc = sizeof(myArgv)/4;

	tetgenbehavior b;
	tetgenio in, addin, bgmin;

	if (!b.parse_commandline(myArgc, myArgv)) {
		terminatetetgen(NULL, 10);
	}

	// Read input files.
	if (b.refine) { // -r
		if (!in.load_tetmesh(b.infilename, (int) b.object)) {
			terminatetetgen(NULL, 10);
		}
	} else { // -p
		if (!in.load_plc(b.infilename, (int) b.object)) {
			terminatetetgen(NULL, 10);
		}
	}
	if (b.insertaddpoints) { // -i
		// Try to read a .a.node file.
		addin.load_node(b.addinfilename);
	}
	if (b.metric) { // -m
		// Try to read a background mesh in files .b.node, .b.ele.
		bgmin.load_tetmesh(b.bgmeshfilename, (int) b.object);
	}

	tetrahedralize(&b, &in, &m_out, &addin, &bgmin);

	// 找出原来的点
	findOriPoints(in, m_out);

	// 构建四面体网格
	buildVegaData();

	// 设置rig初始点
	m_rigObj = new TransformRig();
	if (TransformRig* transRig = dynamic_cast<TransformRig*>(m_rigObj))
	{
		double* oriPnts = new double[m_nSurfPnt*3];
		double* p = oriPnts;
		for (int i = 0; i < m_nSurfPnt; ++i, p+=3)
		{
			double* pnt = m_out.pointlist + m_surfPntIdx[i]*3;
			p[0] = pnt[0];
			p[1] = pnt[1];
			p[2] = pnt[2];
		}
		transRig->setInitPnts(oriPnts, m_nSurfPnt);
		delete[] oriPnts;
	}

	// 计算当前参数下的配置
	m_rigObj->keyParam(0, MathUtilities::zero);
	m_rigObj->keyParam(1, MathUtilities::zero);
	m_rigObj->keyParam(2, MathUtilities::absSin);
	m_rigObj->keyParam(3, MathUtilities::zero);
	m_rigObj->keyParam(4, MathUtilities::zero);
	m_rigObj->keyParam(5, MathUtilities::zero);

	m_nParam = m_rigObj->getNFreeParam();
	m_param.setZero(m_nParam);
	m_controlForce.setZero(m_nParam);
	m_paramVelocity.setZero(m_nParam);
	m_param[0] = m_param[1] = m_param[2] = 1.0;
	m_rigObj->setFreeParam(&m_param[0]);

	computeRig();
}

void RiggedMesh::show()
{
	if (!m_tetMesh)
		return;
	glColor3f(0,0,1);
	glPointSize(3.f);

#ifdef SHOW_TETGEN_RESULT
	// 画出所有顶点
	glBegin(GL_POINTS);
	for (int i = 0; i < m_nTotPnt; ++i)
	{
		glVertex3dv(m_out.pointlist + i * 3);
	}
	glEnd();

	// 画出原来面网格的顶点
	glColor3f(1,0,0);
	glPointSize(5.f);
	glBegin(GL_POINTS);
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		glVertex3dv(m_out.pointlist + m_surfPntIdx[i] * 3);
	}
	glEnd();

	// 画出原来面网格的边
	glColor3f(1,0,0);
	glBegin(GL_LINES);
	for (int i = 0; i < m_out.numberofedges; ++i)
	{
		int eSrc = m_out.edgelist[i*2];
		int eTar = m_out.edgelist[i*2+1];
		glVertex3dv(m_out.pointlist + eSrc*3);
		glVertex3dv(m_out.pointlist + eTar*3);
	}
	glEnd();
#endif

#ifdef SHOW_TETMESH_RESULT
	glBegin(GL_POINTS);
	for (int ithVtx = 0; ithVtx < m_tetMesh->getNumVertices(); ++ithVtx)
	{
		Vec3d pV = *m_tetMesh->getVertex(ithVtx);
		//glVertex3dv(&pV[0]);
		double* dV = &m_q[ithVtx*3];
		glVertex3d(pV[0] + dV[0], pV[1] + dV[1], pV[2] + dV[2]);
	}
	glEnd();

	glColor3f(0.5, 0.5,0.5);
	glBegin(GL_LINES);
	for(int ithTet = 0; ithTet < m_tetMesh->getNumElements(); ++ithTet)
	{
		Vec3d  pnt[4];
		for (int i = 0; i < 4; ++i)
		{
			pnt[i] = *m_tetMesh->getVertex(ithTet, i);
			int idx = m_tetMesh->getVertexIndex(ithTet, i);
			pnt[i] += Vec3d(m_q[idx*3],m_q[idx*3+1],m_q[idx*3+2]);
		}

		int edge[6][2] = {{0,1},{0,2},{0,3},{1,2},{2,3},{3,1}};
		for (int e = 0; e < 6; ++e)
		{
			glVertex3dv(&(pnt[edge[e][0]])[0]);
			glVertex3dv(&(pnt[edge[e][1]])[0]);
		}
	}
	glEnd();
#endif


	// 画出受到rig控制的顶点
	glColor3f(0,1,0);
	glPointSize(3.f);
	glBegin(GL_POINTS);
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx = m_surfPntIdx[i];
		Vec3d p = *m_tetMesh->getVertex(idx);
		p += Vec3d(m_q[idx*3],m_q[idx*3+1],m_q[idx*3+2]);
		glVertex3dv(&p[0]);
	}
	glEnd();

	// 画出受力情况
	glColor3f(1,0,0);
	glBegin(GL_LINES);
	for (int i = 0; i < m_force.size()/3; ++i)
	{
		Vec3d* v = m_tetMesh->getVertex(i);
		double endPnt[3];
		v->convertToArray(endPnt);
		endPnt[0] += m_q[i*3];
		endPnt[1] += m_q[i*3+1];
		endPnt[2] += m_q[i*3+2];
		glVertex3dv(endPnt);
		endPnt[0] -= m_force[i*3]   * 1e-5;
		endPnt[1] -= m_force[i*3+1] * 1e-5;
		endPnt[2] -= m_force[i*3+2] * 1e-5;
		glVertex3dv(endPnt);
	}
	glEnd();
}

bool RiggedMesh::findOriPoints(tetgenio& in, tetgenio& out)
{
	if (out.numberofpoints <= 0)
		return false;
	std::set<PntPair> pntSet;

	double* p = in.pointlist;
	for (int i = 0; i < in.numberofpoints; ++i, p+=3)
	{
		PntPair pair = {Vec3d (p[0], p[1], p[2]), i};
		pntSet.insert(pair);
	}

	p = out.pointlist;
	vector<int> surfPntIdx, intPntIdx;
	for (int i = 0; i < out.numberofpoints; ++i, p+=3)
	{ 
		PntPair pair = {Vec3d (p[0], p[1], p[2]), i};
		std::set<PntPair>::iterator pPair = pntSet.find(pair);
		if (pPair != pntSet.end())
		{
			surfPntIdx.push_back(i);
			pntSet.erase(pPair);
		}
		else
			intPntIdx.push_back(i);
	}

	m_nIntPnt = intPntIdx.size();
	m_nSurfPnt = surfPntIdx.size();
	m_nIntDof = m_nIntPnt*3;
	m_nSurfDof= m_nSurfPnt*3;
	m_nTotPnt = out.numberofpoints;

	m_surfPntIdx = new int[m_nSurfPnt];
	m_intPntIdx  = new int[m_nIntPnt];
	m_surfDofIdx = new int[m_nSurfDof];
	m_intDofIdx  = new int[m_nIntDof];

	std::copy(surfPntIdx.begin(), surfPntIdx.end(), m_surfPntIdx);
	std::copy(intPntIdx.begin(), intPntIdx.end(), m_intPntIdx);

	for (int i = 0; i < m_nIntPnt; ++i)
	{
		m_intDofIdx[i*3]   = 3 * m_intPntIdx[i];
		m_intDofIdx[i*3+1] = 3 * m_intPntIdx[i] + 1;
		m_intDofIdx[i*3+2] = 3 * m_intPntIdx[i] + 2;
	}
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		m_surfDofIdx[i*3]   = 3 * m_surfPntIdx[i];
		m_surfDofIdx[i*3+1] = 3 * m_surfPntIdx[i] + 1;
		m_surfDofIdx[i*3+2] = 3 * m_surfPntIdx[i] + 2;
	}

	if (m_nIntPnt <= 0)
		PRINT_F("no internal points, please increase tetrahedral resolution");
	if (m_nSurfPnt <= 0)
		PRINT_F("no rigged points");
	return m_nIntPnt > 0 && m_nSurfPnt > 0 && m_nTotPnt > 0;
}

bool RiggedMesh::PntPair::operator<( const PntPair& other ) const
{
	if (m_pnt[0] != other.m_pnt[0])
		return m_pnt[0] < other.m_pnt[0];
	if (m_pnt[1] != other.m_pnt[1])
		return m_pnt[1] < other.m_pnt[1];
	return m_pnt[2] < other.m_pnt[2];
}

bool RiggedMesh::buildVegaData(double E, double nu, double density)
{
	if (!m_out.numberofpoints || !m_out.numberoftetrahedra || !m_out.numberofcorners)
		return false;
	
	delete m_tetMesh;
	delete m_forceModel;

	int  numTet = m_out.numberoftetrahedra;
	int* tetIDs = new int[numTet * 4];
	int* srcTet = m_out.tetrahedronlist;
	int* tarTet = tetIDs;
	for (int ithTet = 0; ithTet < numTet; ++ithTet, tarTet+=4, srcTet += m_out.numberofcorners)
	{
		tarTet[0] = srcTet[0];
		tarTet[1] = srcTet[1];
		tarTet[2] = srcTet[2];
		tarTet[3] = srcTet[3];
	}
 
	int nPnts = m_out.numberofpoints;
	m_tetMesh = new TetMesh(nPnts, m_out.pointlist, numTet, tetIDs, E, nu, density);
	delete[] tetIDs;


	int wrap = 2;
	CorotationalLinearFEMWrapper* wrapper = new CorotationalLinearFEMWrapper(m_tetMesh, wrap);
	m_modelwrapper    = wrapper;
	m_forceModel      = new CorotationalLinearFEMForceModel(wrapper, wrap);
	
/*
	NeoHookeanIsotropicMaterial* hookeanMat = new NeoHookeanIsotropicMaterial(m_tetMesh, 1, 0.5);
	IsotropicHyperelasticFEM* defModel = new IsotropicHyperelasticFEM(m_tetMesh, hookeanMat);
	m_forceModel = new IsotropicHyperelasticFEMForceModel(defModel);*/

	m_q = m_v = m_a = m_force = m_extForce = EigVec::Constant(nPnts*3, 0.0);
	m_mass = EigVec(nPnts);
	GenerateMassMatrix::computeVertexMasses(m_tetMesh, &m_mass[0]);


	if (m_tangentStiffnessMatrix)
		delete m_tangentStiffnessMatrix;
	allocateTempData();

	computeElementLaplacian(m_eleAdjList);
	return true;
}

void RiggedMesh::computeRig()
{
	EigVec res(m_nSurfPnt * 3);
	computeSurfOffset(&m_param[0], m_t, res);

	double* pQ = &m_q[0];
	for (int i = 0; i < m_nTotPnt * 3; ++i)
		pQ[i] = 0;

	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		Vec3d d(res[i*3], res[i*3+1], res[i*3+2]);
		int idx = m_surfPntIdx[i];
		pQ[idx*3] = d[0];
		pQ[idx*3+1] = d[1];
		pQ[idx*3+2] = d[2];
	}

	m_forceModel->GetInternalForce(&m_q[0], &m_force[0]);
}

bool RiggedMesh::computeHessian(const EigVec& n, const EigVec& p, double t, EigSparse& Hnn, EigDense& Hnp, EigDense& Hpn, EigDense* pHpp)
{ 
	bool res = true;
	// 计算指定参数下各个点的偏移量q
	EigVec q(m_nTotPnt*3);
	res &= computeQ(&n[0], &p[0], t, &q[0]);

	// 提取内力、tangent Stiffness matrix、质量矩阵
	// 注意tangent stiffness matrix为实际的负值，因为系统计算出的弹力为实际的负值,因此需要先反转
	EigVec force(m_nTotPnt*3);
	m_forceModel->GetForceAndMatrix(&q[0], &force[0], m_tangentStiffnessMatrix);
	*m_tangentStiffnessMatrix *= -1;
	force *= -1;

	// 提取 dFn/dn dFn/ds dFs/dn dFs/ds,
	// 其中Fn为内部节点受到的力，Fs为表面节点受到的力
	// n为内部节点位置，s为表面节点位置

	// 构建矩阵Hnn = Mn / h^2 - dFn/dn
// 	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_intDofIdx, m_intDofIdx, m_nIntDof, m_nIntDof, Hnn);
// 	Hnn *= -1.0;
	Utilities::vegaSparse2AllocatedEigen(*m_tangentStiffnessMatrix, m_intDofIdx, m_intDofIdx, m_nIntDof, m_nIntDof, m_dFnn);
	Hnn  = -m_dFnn;
	for (int ithInt = 0; ithInt < m_nIntPnt; ++ithInt)
	{
		int begIdx = ithInt*3;
		double massTerm = m_mass[m_intPntIdx[ithInt]] / (m_h * m_h);
		Hnn.coeffRef(begIdx, begIdx) += massTerm;	begIdx++;
		Hnn.coeffRef(begIdx, begIdx) += massTerm;	begIdx++;
		Hnn.coeffRef(begIdx, begIdx) += massTerm;
	}

	// 构建矩阵Hnp = - dFn/ds * J
	EigDense J;
	res &= m_rigObj->computeJacobian(J);
	EigSparse dFns;
// 	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_intDofIdx, m_surfDofIdx, m_nIntDof, m_nSurfDof, dFns);
// 	Hnp = -1.f * dFns * J;
	Utilities::vegaSparse2AllocatedEigen(*m_tangentStiffnessMatrix, m_intDofIdx, m_surfDofIdx, m_nIntDof, m_nSurfDof, m_dFns);
	Hnp = -1.f * m_dFns * J;

	// 构建矩阵Hpn = - J' * dFs/dn
	EigSparse dFsn;
// 	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_intDofIdx, m_nSurfDof, m_nIntDof,dFsn);
// 	Hpn = - J.transpose() * dFsn;
	Utilities::vegaSparse2AllocatedEigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_intDofIdx, m_nSurfDof, m_nIntDof, m_dFsn);
	Hpn = -J.transpose() * m_dFsn;

	// 构建矩阵Hpp,此矩阵计算量最大
	if (pHpp)
	{
		EigVec jacobianDeri(m_nSurfPnt*3);
		EigVec sResi(m_nSurfPnt*3);
		double invH  = 1.0 / m_h;
		double invH2 = 1.0 / (m_h * m_h);
		for (int ithSurf = 0; ithSurf < m_nSurfPnt; ++ithSurf)
		{
			int idx = m_surfPntIdx[ithSurf];
			int idx3 = idx * 3;
			int ithSurf3 = ithSurf * 3;
			double m = m_mass[idx];
			// 某个表面点的xyz
			sResi[ithSurf3] = m * ((q[idx3] - m_q[idx3]) * invH2 - m_v[idx3] * invH) - force[idx3]; idx3++; ithSurf3++;
			sResi[ithSurf3] = m * ((q[idx3] - m_q[idx3]) * invH2 - m_v[idx3] * invH) - force[idx3]; idx3++; ithSurf3++;
			sResi[ithSurf3] = m * ((q[idx3] - m_q[idx3]) * invH2 - m_v[idx3] * invH) - force[idx3];
		}
		EigDense& Hpp = *pHpp;
		Hpp = EigDense::Zero(m_nParam, m_nParam);
// 		EigSparse dFss;
// 		Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_surfDofIdx, m_nSurfDof, m_nSurfDof,dFss);
		Utilities::vegaSparse2AllocatedEigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_surfDofIdx, m_nSurfDof, m_nSurfDof,m_dFss);
		EigSparse& dFss = m_dFss;
		for (int i = 0; i < m_nParam; ++i)
		{
			for (int l = 0; l < m_nParam; ++l)
			{
				double& Hppil = Hpp(i,l);
				res &= m_rigObj->computeJacobianDerivative(i,l, &jacobianDeri[0]);
				for (int k = 0; k < m_nSurfPnt*3; ++k)
				{
					Hppil += jacobianDeri[k] * sResi[k];

					double mass = m_mass[m_surfPntIdx[k/3]];
					double v = 0;
#ifdef ROW_MAJOR
					EigSparse::InnerIterator it(dFss, k);
					for (; it; ++it)  
					{
						int m = it.col();
						v += it.value() * J(m,l);
					}
#else
					for (int m = 0; m < m_nSurfPnt*3; ++m)
					{
						v += dFss.coeff(k,m) * J(m,l);
					}
#endif
					Hppil += J(k,i) *(mass * J(k,l) * invH2 - v);
				}
			}
		}

		// 隐式控制的情况下，控制力是位移差和速度差的线性函数
		if (m_controlType == CONTROL_IMPLICIT_FORCE)
		{
			for (int ithP = 0; ithP < m_nParam; ++ithP)
			{
				Hpp(ithP, ithP) += m_propGain[ithP] + m_deriGain[ithP] / m_h;
			}
		}
	}

	return res;
}

void RiggedMesh::testStep( double dt )
{
	int nDeg = m_nTotPnt*3;
	if (m_q.size() != nDeg || m_force.size() != nDeg)
		return;
	if (m_v.size() != nDeg || m_a.size() != nDeg)
	{
		m_v = m_a = EigVec::Constant(nDeg,0.0);
	}

	vector<double> mass(m_nTotPnt);
	int nIter = 1;
	dt /= nIter;
	while(nIter--)
	{
		m_forceModel->GetInternalForce(&m_q[0], &m_force[0]);
		GenerateMassMatrix::computeVertexMasses(m_tetMesh, &mass[0]);
		for (int i = 0; i < m_nTotPnt; ++i)
		{
			int idx = i*3;
			m_a[idx] = -m_force[idx] / mass[i] + m_v[idx] * -10.1;  idx++;
			m_a[idx] = -m_force[idx] / mass[i] + m_v[idx] * -10.1;  idx++;
			m_a[idx] = -m_force[idx] / mass[i] + m_v[idx] * -10.1;
		}

		for (int i = 0; i < nDeg; ++i)
		{
			m_q[i] += m_v[i] * dt;
			m_v[i] += m_a[i] * dt;
		}
	}
}

bool RiggedMesh::computeQ( const double* n, const double* p, double t, double* q )
{
	EigVec& s = m_surfDofBuf;
	bool res = computeSurfOffset(p, t, s);
	if (!res)
		return res;
	for (int i = 0; i < m_nIntPnt; ++i)
	{
		int idx = m_intPntIdx[i];
		q[idx*3]   = n[i*3];
		q[idx*3+1] = n[i*3+1];
		q[idx*3+2] = n[i*3+2];
	}
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx = m_surfPntIdx[i];
		q[idx*3]   = s[i*3];
		q[idx*3+1] = s[i*3+1];
		q[idx*3+2] = s[i*3+2];
	}
	return true;
}

bool RiggedMesh::computeQ( const double* p, double t, double* q )
{
	EigVec& s = m_surfDofBuf;
	bool res = computeSurfOffset(p, t, s);
	if (!res)
		return res;
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx = m_surfPntIdx[i];
		q[idx*3]   = s[i*3];
		q[idx*3+1] = s[i*3+1];
		q[idx*3+2] = s[i*3+2];
	}
	return true;
}

bool RiggedMesh::computeSurfOffset( const double* p, double t, EigVec& s )
{
	s.resize(m_nSurfPnt*3);
	m_rigObj->setTime(t);
	bool res = m_rigObj->computeValue(&s[0], p);
	if (!res)
		return res;
	double *ps = &s[0];
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx3 = m_surfPntIdx[i] * 3;
		ps[i*3]   -= m_out.pointlist[idx3++];
		ps[i*3+1] -= m_out.pointlist[idx3++];
		ps[i*3+2] -= m_out.pointlist[idx3];
	}
	return true;
}

bool RiggedMesh::computeGradient( const EigVec& n, const EigVec& p, double t, EigVec& gn, EigVec& gs )
{
	// 计算给定参数下状态
	EigVec q(m_nTotPnt*3);
	bool res = computeQ(&n[0], &p[0], t, &q[0]);

	EigVec ma = (q - m_q) / (m_h * m_h) - m_v / m_h;
	for (int i = 0; i < m_nTotPnt; ++i)
	{
		ma[i*3]   *= m_mass[i];
		ma[i*3+1] *= m_mass[i];
		ma[i*3+2] *= m_mass[i];
	}

	// 计算此状态下的力
	EigVec f(m_nTotPnt*3);
	m_forceModel->GetInternalForce(&q[0], &f[0]);
	f *= -1.0;
	EigVec residual = ma - f - m_extForce;

	gn.setZero(m_nIntPnt*3);
	for (int ithN = 0; ithN < m_nIntPnt; ++ithN)
	{
		int idx = m_intPntIdx[ithN];
		gn[ithN*3]   = residual[idx*3];
		gn[ithN*3+1] = residual[idx*3+1];
		gn[ithN*3+2] = residual[idx*3+2];
	}

	gs.setZero(m_nSurfPnt*3);
	for (int ithS = 0; ithS < m_nSurfPnt; ++ithS)
	{
		int idx = m_surfPntIdx[ithS];
		gs[ithS*3]   = residual[idx*3];
		gs[ithS*3+1] = residual[idx*3+1];
		gs[ithS*3+2] = residual[idx*3+2];
	}
	EigDense J;
	res &= m_rigObj->computeJacobian(J);
	gs = J.transpose() * gs;
	if (m_controlType == CONTROL_EXPLICIT_FORCE)
	{
		gs -= m_controlForce;
	}
	else if (m_controlType == CONTROL_IMPLICIT_FORCE)
	{
		computeControlForce(m_controlForce, p);
		gs -= m_controlForce;
	}


	return res;
}

bool RiggedMesh::testHessian()
{
	// 扰动上一帧状态
	double lastT = 0;
	double dir = (rand() % 2) ? 1.0 : -1.0;
	double noise = 0.02 * dir;
	m_q += EigVec::Random(m_nTotPnt*3) * noise;
	m_v += EigVec::Random(m_nTotPnt*3) * noise;
	m_a += EigVec::Random(m_nTotPnt*3) * noise ;
	m_param += EigVec::Random(m_nParam)* noise;
	computeQ(&m_param[0], lastT, &m_q[0]);

	// 在上一帧的基础上计算新的假想位置
	double noise1 = ((rand() % 2) ? 1.0 : -1.0) * 0.02;
	EigVec n(m_nIntPnt*3);
	for (int ithN = 0; ithN < m_nIntPnt; ++ithN)
	{
		int idx = m_intPntIdx[ithN];
		n[ithN*3] = m_q[idx*3];
		n[ithN*3+1] = m_q[idx*3+1];
		n[ithN*3+2] = m_q[idx*3+2];
	}
	n += EigVec::Random(m_nIntPnt*3)  * noise1;

	EigVec p = m_param;
	p += EigVec::Random(m_nParam) * noise1;
	EigVec dN = EigVec::Random(m_nIntPnt*3) * noise1;
	EigVec dP = EigVec::Random(m_nParam)  * noise1;

	EigVec gn, gp;
	EigDense Hnp, Hpn, Hpp;
	EigSparse Hnn;
	computeGradient(n,p, lastT, gn,gp);
	computeHessian(n,p, lastT, Hnn, Hnp, Hpn, &Hpp);
	PRINT_F("test Hessian\n");
	for (double step = 100; step > 1e-15; step *= 0.1)
	{
		EigVec dNi = dN * step;
		EigVec dPi = dP * step;

		// 计算近似的梯度值
		EigVec gni_app = gn + Hnn * dNi + Hnp * dPi;
		EigVec gpi_app = gp + Hpn * dNi + Hpp * dPi;

		// 计算准确的梯度值
		EigVec ni = n + dNi;
		EigVec pi = p + dPi;
		EigVec gni, gpi;
		computeGradient(ni, pi, lastT, gni, gpi);

		EigVec resN = gni - gni_app;
		EigVec resP = gpi - gpi_app;

		double resNNorm = resN.norm();
		double resPNorm = resP.norm();
		resN = resN.cwiseAbs();
		resP = resP.cwiseAbs();
		double resNMax  = resN.maxCoeff();
		double resPMax  = resP.maxCoeff();

		double invE = 1.0 / (step * step);
		PRINT_F("ε=%le  maxResN=%le  maxResP=%le  |maxN|/ε^2=%le  |maxP|/ε^2=%le\n",
			step, resNMax, resPMax, resNMax * invE, resPMax* invE);
	}
	return true;
}






bool RiggedMesh::testElasticEnergy()
{
	int nDof = m_nTotPnt * 3;

	// 扰动当前位置
	double noise = 0.1;
	EigVec q = m_q;
	q += (EigVec::Random(nDof) - EigVec::Constant(nDof, 0.5))* noise;

	// 中心差商计算梯度，也就是力
	double dq = 1e-5;
	EigVec  f(nDof);
	for (int ithDof = 0; ithDof < nDof; ++ithDof)
	{
		double oriPos = q[ithDof];

		q[ithDof] = oriPos - dq;
		double e0 = m_modelwrapper->computeElasticEnergy(&q[0]);

		q[ithDof] = oriPos + dq;
		double e1 = m_modelwrapper->computeElasticEnergy(&q[0]);

		f[ithDof] = (e1 - e0) / (2 * dq);

		q[ithDof] = oriPos;
	}

	EigVec fRef(nDof);
	m_forceModel->GetInternalForce(&q[0], &fRef[0]);

	EigVec error = fRef - f;
	PRINT_F("test result\n");
	for (int ithDof = 0; ithDof < nDof; ++ithDof)
	{
		PRINT_F("f=%lf \t fRef=%lf \t error = %lf \n", 
			f[ithDof], fRef[ithDof], error[ithDof]);
	}
	PRINT_F("\n");
	PRINT_F("fNorm = %lf \t fRefNorm = %lf errorNorm = %lf\n", f.norm(), fRef.norm(), error.norm());
	PRINT_F("errorMax = %lf\n", error.cwiseAbs().maxCoeff());

	
	m_q = q;
	return true;
}

double RiggedMesh::computeValue( const EigVec& n, const EigVec& p, double t)
{

	// 计算给定参数下状态
	EigVec q(m_nTotPnt*3);
	computeQ(&n[0], &p[0], t, &q[0]);

	return computeValue(q,p);
}

double RiggedMesh::computeValue( const EigVec& q, const EigVec& p )
{
	// 计算动能增量
	double kinetic = 0;
	EigVec acc = (q - m_q) / (m_h * m_h) - m_v / m_h;
	for (int i = 0, idx = 0; i < m_nTotPnt; ++i)
	{
		kinetic += acc[idx] * m_mass[i] * acc[idx];	idx++;
		kinetic += acc[idx] * m_mass[i] * acc[idx];	idx++;
		kinetic += acc[idx] * m_mass[i] * acc[idx];	idx++;
	}
	kinetic = kinetic * m_h * m_h / 2;

	// 计算弹性能量
	double elasticEnergy = m_modelwrapper->computeElasticEnergy(&q[0]);

	// 计算外力能量 = -外力做功 = -f * (q - m_q)
	double externalEnergy = m_extForce.dot(m_q-q);
	double totEnergy = elasticEnergy + kinetic + externalEnergy;
	if (m_controlType == CONTROL_EXPLICIT_FORCE)
	{
		double ctrlEnergy     = m_controlForce.dot(m_param);
		totEnergy -= ctrlEnergy;
	}
	else if (m_controlType == CONTROL_IMPLICIT_FORCE)
	{
		double energy = 0;
		if(computeControlEnergy(p, energy))
			totEnergy += energy;
	}
	return totEnergy;
}

bool RiggedMesh::testValue()
{
	// 扰动上一帧状态
	double lastT = 0;
	double dir = (rand() % 2) ? 1.0 : -1.0;
	double noise = 0.05 * dir;
	m_q += EigVec::Random(m_nTotPnt*3) * noise;
	m_v += EigVec::Random(m_nTotPnt*3) * noise;
	m_a += EigVec::Random(m_nTotPnt*3) * noise ;
	m_param += EigVec::Random(m_nParam)* noise;
	computeQ(&m_param[0], lastT, &m_q[0]);

	// 在上一帧的基础上计算新的假想位置 n,p
	double noise1 = ((rand() % 2) ? 1.0 : -1.0) * 0.05;
	EigVec n(m_nIntPnt*3);
	for (int ithN = 0; ithN < m_nIntPnt; ++ithN)
	{
		int idx = m_intPntIdx[ithN];
		n[ithN*3] = m_q[idx*3];
		n[ithN*3+1] = m_q[idx*3+1];
		n[ithN*3+2] = m_q[idx*3+2];
	}
	n += EigVec::Random(m_nIntPnt*3)  * noise1;

	EigVec p = m_param;
	p += EigVec::Random(m_nParam) * noise1;
	EigVec dN = EigVec::Random(m_nIntPnt*3) * noise1;
	EigVec dP = EigVec::Random(m_nParam)  * noise1;

	// 计算函数值、梯度值
	EigVec gn, gp;
	computeGradient(n,p, lastT, gn,gp);
	double f0 = computeValue(n, p, lastT);
	PRINT_F("test func val\n");
	for (double step = 100; step > 1e-15; step *= 0.1)
	{
		EigVec dNi = dN * step;
		EigVec dPi = dP * step;

		// 计算近似的函数值
		double fi_app = f0 + gn.dot(dNi) + gp.dot(dPi);

		// 计算准确的函数值
		EigVec ni = n + dNi;
		EigVec pi = p + dPi;
		double fi = computeValue(ni, pi, lastT);

		double error = abs(fi-fi_app);
		PRINT_F("step = %le funVal = %le approxVal = %le error = %le error/dx^2 = %le\n", step, fi, fi_app, error, error/(step*step));
	}
	PRINT_F("\n");
	return true;
}

bool RiggedMesh::computeValueAndGrad(const EigVec& x, const EigVec& param,  double* v, EigVec* grad)
{
	if (!v && !grad)
		return false;

	double t = param[0];
	bool res = true;
	// 计算给定参数下状态
	EigVec n(m_nIntPnt*3), p(m_nParam);
	for (int i = 0; i < m_nIntPnt*3; ++i)
	{
		n[i] = x[i];
	}
	for (int i = 0; i < m_nParam; ++i)
	{
		p[i] = x[i+m_nIntPnt*3];
	}

	EigVec q(m_nTotPnt*3);
	res &= computeQ(&n[0], &p[0], t, &q[0]);

	// 计算函数值
	if (v)
	{
		*v = computeValue(q, p);
	}

	// 计算梯度值
	if (grad)
	{
		EigVec ma = (q - m_q) / (m_h * m_h) - m_v / m_h;
		for (int i = 0; i < m_nTotPnt; ++i)
		{
			ma[i*3]   *= m_mass[i];
			ma[i*3+1] *= m_mass[i];
			ma[i*3+2] *= m_mass[i];
		}
		//PRINT_F("|m_q| = %lf, |m_v| = %lf, |m_a| = %lf, |m_param| = %lf", m_q.norm(), m_v.norm(), m_a.norm(), m_param.norm());

		// 计算此状态下的力
		EigVec f(m_nTotPnt*3);
		m_forceModel->GetInternalForce(&q[0], &f[0]);
		f *= -1.0;
		EigVec residual = ma - f - m_extForce;

		EigVec gn(m_nIntPnt*3);
		for (int ithN = 0; ithN < m_nIntPnt; ++ithN)
		{
			int idx = m_intPntIdx[ithN];
			gn[ithN*3]   = residual[idx*3];
			gn[ithN*3+1] = residual[idx*3+1];
			gn[ithN*3+2] = residual[idx*3+2];
		}

		EigVec gs(m_nSurfPnt*3);
		for (int ithS = 0; ithS < m_nSurfPnt; ++ithS)
		{
			int idx = m_surfPntIdx[ithS];
			gs[ithS*3]   = residual[idx*3];
			gs[ithS*3+1] = residual[idx*3+1];
			gs[ithS*3+2] = residual[idx*3+2];
		}
		EigDense J;
		res &= m_rigObj->computeJacobian(J);
		gs = J.transpose() * gs;
		if (m_controlType == CONTROL_EXPLICIT_FORCE)
		{
			gs -= m_controlForce;
		}
		else if (m_controlType == CONTROL_IMPLICIT_FORCE)
		{
			computeControlForce(m_controlForce, p);
			gs -= m_controlForce;
		}

		grad->resizeLike(x);
		EigVec& g = *grad;
		for (int i = 0; i < m_nIntPnt*3; ++i)
		{
			g[i] = gn[i];
		}
		for (int i = 0; i < m_nParam; ++i)
		{
			g[i+m_nIntPnt*3] = gs[i];
		}

	}
	return res;
}

void RiggedMesh::setDof( EigVec&n, EigVec&p, bool proceedTime /*= true*/ )
{
	EigVec q(m_nTotPnt*3);
	computeQ(&n[0], &p[0], m_t, &q[0]);
	EigVec v = (q - m_q) / m_h;
	EigVec a = (v - m_v) / m_h; 
	EigVec pV= (p - m_param) / m_h;
	m_param = p;
	m_paramVelocity = pV;
	m_v = v;
	m_a = a;
	m_q = q;
	if (proceedTime)
	{
		m_t += m_h;
	}
}

void RiggedMesh::getDof( EigVec& n, EigVec& p )
{
	p = m_param;
	n.resize(m_nIntPnt*3);
	for (int ithInt = 0; ithInt < m_nIntPnt; ++ithInt)
	{
		int idx = m_intPntIdx[ithInt];
		n(ithInt*3)   = m_q(idx*3);
		n(ithInt*3+1) = m_q(idx*3+1);
		n(ithInt*3+2) = m_q(idx*3+2);
	}
}

bool RiggedMesh::getN(RigStatus& s, EigVec& n)
{
	const EigVec& q = s.getQ();
	if (q.size() != m_nTotPnt * 3)
		return false;

	n.resize(m_nIntPnt*3);
	for (int ithInt = 0; ithInt < m_nIntPnt; ++ithInt)
	{
		int idx = m_intPntIdx[ithInt];
		n(ithInt*3)   = q(idx*3);
		n(ithInt*3+1) = q(idx*3+1);
		n(ithInt*3+2) = q(idx*3+2);
	}
	return true;
}


bool RiggedMesh::init( tetgenio& surfMesh, RigBase* rig, double maxVolume /*= 1*/, double edgeRatio /*= 2*/ , double youngModulus, double nu, double density)
{
	if (surfMesh.numberofpoints <= 0 || surfMesh.numberoffacets <= 0 ||
		!surfMesh.pointlist || !surfMesh.facetlist ||
		!rig || !rig->getNFreeParam())
	{
		PRINT_F("mesh data is null or no rig param");
		return false;
	}
	clear();

	// 估计四面体个数，拒绝个数太多的参数
	double minPnt[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, maxPnt[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
	for (double* p = surfMesh.pointlist; 
		 p < surfMesh.pointlist + surfMesh.numberofpoints*3; p+=3)
	{
		for (int i = 0; i < 3; ++i)
		{
			minPnt[i] = p[i] < minPnt[i] ? p[i] : minPnt[i];
			maxPnt[i] = p[i] > maxPnt[i] ? p[i] : maxPnt[i];
		}
	}
	double bboxVol =	(maxPnt[0] - minPnt[0]) *
						(maxPnt[1] - minPnt[1]) * 
						(maxPnt[2] - minPnt[2]);
	double approxTetNum = bboxVol / maxVolume;
	if (approxTetNum > 50000)
	{
		maxVolume = bboxVol / 50000;
		PRINT_F("max Volume is too small, set to %lf", maxVolume);
	}

	char cmd[100];
	sprintf(cmd, "pq%lfa%lf", edgeRatio, maxVolume);

	tetrahedralize(cmd, &surfMesh, &m_out);

	// 找出原来的点
	bool res = findOriPoints(surfMesh, m_out);

	if (!res)
	{
		PRINT_F("failed to find surface points");
		clear();
		return false;
	}

	// 构建四面体网格,为模拟的位置、速度、加速度等参数分配空间
	res &= buildVegaData(youngModulus, nu, density);

	if (!res)
	{
		PRINT_F("failed to build fem data");
		clear();
		return false;
	}

	// 设置rig参数
	m_rigObj = rig;
	m_nParam = rig->getNFreeParam();
	m_param.resize(m_nParam);
	rig->getFreeParam(&m_param[0]);
	m_paramVelocity.setZero(m_nParam);
	PRINT_F("Initialized completed.\n internal Pnts:%d\n surface Pnts:%d\n params:%d\n", 
			m_nIntPnt, m_nSurfPnt, m_nParam);
	return true;
}

void RiggedMesh::clear()
{
	m_rigObj = NULL;

	m_out.deinitialize();
	m_out.initialize();
	delete[] m_surfPntIdx;
	delete[] m_intPntIdx;
	delete[] m_surfDofIdx;
	delete[] m_intDofIdx;
	m_surfPntIdx = m_intPntIdx = m_surfDofIdx = m_intDofIdx = NULL;


	m_nTotPnt = 0;
	m_nIntPnt = 0;
	m_nSurfPnt = 0;
	m_nIntDof = 0;
	m_nSurfDof = 0;
	m_nParam = 0;

	delete m_forceModel;
	m_forceModel = NULL;

	delete m_modelwrapper;
	m_modelwrapper = NULL;

	delete m_tetMesh;
	m_tetMesh = NULL;

	m_q = m_v = m_a = m_force = m_param = EigVec();

	m_mass = EigVec();
	m_t = 0;

	freeTempData();
}

RigStatus RiggedMesh::getStatus() const
{
	return RigStatus(m_q, m_v, m_a, m_param, m_paramVelocity, m_extForce, m_targetParam);
}

bool RiggedMesh::setStatus( const RigStatus& s )
{
	int pntLength = s.getPointVecLength();
	int paramLength = s.getParamVecLength();
	if (pntLength == m_nTotPnt * 3 &&
		paramLength == m_nParam &&
		pntLength == m_q.size() &&
		pntLength == m_v.size() &&
		pntLength == m_a.size() &&
		paramLength == m_param.size() &&
		paramLength == m_paramVelocity.size())
	{
		m_q = s.getQ();
		m_v = s.getV();
		m_a = s.getA();
		m_param = s.getP();
		m_paramVelocity = s.getPV();
		m_extForce = s.getF();
		m_targetParam = s.getTarP();
		return true;
	}
	return false;
}

bool RiggedMesh::showStatus( RigStatus& s, const MeshDispConfig& config , double* boundingbox)
{
	if (!s.matchLength(m_nTotPnt*3, m_nParam))
		return false;

	glPushAttrib(GL_CURRENT_BIT);
	const EigVec& q = s.getQ();
	const EigVec& f = s.getF();

	// 画点
	glPointSize(3.f);
	glBegin(GL_POINTS);
	glColor3f(0,0,1);
	if (config.m_showPoint)
	{
		for (int i = 0; i < m_nIntPnt; ++i)
		{
			int idx = m_intPntIdx[i];
			Vec3d v = *m_tetMesh->getVertex(idx);
			Vec3d dV(q[idx*3], q[idx*3+1], q[idx*3+2]);
			v += dV;
			glVertex3d(v[0], v[1], v[2]);
		}

		glColor3f(0,1,0);
		for (int i = 0; i < m_nSurfPnt; ++i)
		{
			int idx = m_surfPntIdx[i];
			Vec3d v = *m_tetMesh->getVertex(idx);
			Vec3d dV(q[idx*3], q[idx*3+1], q[idx*3+2]);
			v += dV;
			glVertex3d(v[0], v[1], v[2]);
		}
	}
	glEnd();

	// 计算包围盒
	double bbox[6];
	bbox[0] = bbox[1] = bbox[2] = DBL_MAX;
	bbox[3] = bbox[4] = bbox[5] = -DBL_MAX;
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx = m_surfPntIdx[i];
		Vec3d v = *m_tetMesh->getVertex(idx);
		Vec3d dV(q[idx*3], q[idx*3+1], q[idx*3+2]);
		v += dV;
		for (int j = 0; j < 3; ++j)
		{
			bbox[j] = min(bbox[j], v[j]);
			bbox[j+3]=max(bbox[j+3], v[j]);
		}
	}
	if (boundingbox)
	{
		for(int i = 0; i < 6; ++i)
			boundingbox[i] = bbox[i];
	}

	// 画力
	double extForceFactor = config.m_extForceDispFactor;
	if (extForceFactor > 0)
	{
		glColor3f(1,1,0);
		glBegin(GL_LINES);
		for (int i = 0; i < m_nTotPnt; ++i)
		{
			int i3 = i*3;
			Vec3d v = *m_tetMesh->getVertex(i);
			v[0] += q[i3];
			v[1] += q[i3+1];
			v[2] += q[i3+2];
			glVertex3d(v[0], v[1], v[2]);
			v[0] += f[i3] * extForceFactor;
			v[1] += f[i3+1] * extForceFactor;
			v[2] += f[i3+2] * extForceFactor;
			glVertex3d(v[0], v[1], v[2]);
		}
		glEnd();
	}

	// 画边
	if (config.m_showEdge)
	{
		glColor3f(0.5, 0.5,0.5);
		glBegin(GL_LINES);
		for(int ithTet = 0; ithTet < m_tetMesh->getNumElements(); ++ithTet)
		{
			Vec3d  pnt[4];
			for (int i = 0; i < 4; ++i)
			{
				pnt[i] = *m_tetMesh->getVertex(ithTet, i);
				int idx = m_tetMesh->getVertexIndex(ithTet, i);
				pnt[i] += Vec3d(q[idx*3],q[idx*3+1],q[idx*3+2]);
			}
			int edge[6][2] = {{0,1},{0,2},{0,3},{1,2},{2,3},{3,1}};
			for (int e = 0; e < 6; ++e)
			{
				glVertex3dv(&(pnt[edge[e][0]])[0]);
				glVertex3dv(&(pnt[edge[e][1]])[0]);
			}
		}
		glEnd();
	}

	// 画面
	if (config.m_showMaterial)
	{
		int axis = config.m_axis;
		double minPlane = bbox[axis] - 0.01;
		double maxPlane = bbox[axis+3] + 0.01;
		double threshold = minPlane + (maxPlane - minPlane) * config.m_axisFactor;

		config.beginColorMaterial();
		glBegin(GL_TRIANGLES);
		for(int ithTet = 0; ithTet < m_tetMesh->getNumElements(); ++ithTet)
		{
			Vec3d  pnt[4];
			for (int i = 0; i < 4; ++i)
			{
				pnt[i] = *m_tetMesh->getVertex(ithTet, i);
				int idx = m_tetMesh->getVertexIndex(ithTet, i);
				pnt[i] += Vec3d(q[idx*3],q[idx*3+1],q[idx*3+2]);
			}

			int face[4][3] = {{0,1,3},{1,2,3},{0,2,3},{0,2,1}};
			if (pnt[0][axis] <= threshold ||
				pnt[1][axis] <= threshold ||
				pnt[2][axis] <= threshold ||
				pnt[3][axis] <= threshold)
			{
				double factor = m_modelwrapper->getElementMaterialFactor(ithTet);
				double norFactor = 
					(factor - config.m_minMaterialFactor) / 
					(config.m_maxMaterialFactor - config.m_minMaterialFactor);

				int nColor = 100;
				int idx = CLAMP_INT(0,nColor-1, int(norFactor * (nColor-1)));

				unsigned char* color = &config.s_colorRuler[idx*3];
				glColor3ubv(color);
				for (int f = 0; f < 4; ++f)
				{
					Vec3d normal = MathUtilities::computeNormal(pnt[face[f][0]],pnt[face[f][1]],pnt[face[f][2]]);
					glNormal3dv(&normal[0]);
					glVertex3dv(&(pnt[face[f][0]])[0]);
					glNormal3dv(&normal[0]);
					glVertex3dv(&(pnt[face[f][1]])[0]);
					glNormal3dv(&normal[0]);
					glVertex3dv(&(pnt[face[f][2]])[0]);
				}
			}
		}
		glEnd();
		config.endColorMaterial();
	}

	if (config.m_showBBox)
	{
		int pntIdx[][3] = {
			{0,1,2},{3,1,2},{3,4,2},{0,4,2},
			{0,1,5},{3,1,5},{3,4,5},{0,4,5}
		};
		
		int edgeIdx[][2] = {
			{0,1},{1,2},{2,3},{3,0},
			{0,4},{1,5},{2,6},{3,7},
			{4,5},{5,6},{6,7},{7,4}
		};
		glColor3f(0.8,0.8,0.8);
		glBegin(GL_LINES);
		for(int e = 0; e < 12; ++e)
		{
			int i0 = edgeIdx[e][0];
			int i1 = edgeIdx[e][1];

			int* p0i = pntIdx[i0];
			int* p1i = pntIdx[i1];

			double p0[] = {bbox[p0i[0]], bbox[p0i[1]], bbox[p0i[2]]};
			double p1[] = {bbox[p1i[0]], bbox[p1i[1]], bbox[p1i[2]]};

			glVertex3dv(p0);
			glVertex3dv(p1);
		}
		glEnd();

	}

	glPopAttrib();
	return true;
}




void FEMSystem::clearResult()
{
	m_solver.clearResult();
}

void FEMSystem::saveResult( const char* fileName )
{
	m_solver.saveResult(fileName);
}

void FEMSystem::step()
{
	m_solver.step();
}

bool RiggedMesh::testCurFrameHessian( RigStatus& lastFrame, RigStatus& curFrame, double noiseN /*= 1.0*/, double noiseP /*= 1.0*/ )
{
	// 当前帧的配置
	PRINT_F("Test Gradient and Hessian");
	PRINT_F("noiseN = %lf, noiseP = %lf", noiseN, noiseP);
	setStatus(lastFrame);
	EigVec n, p;
	if(!getN(curFrame, n))
		return false;
	p = curFrame.getP();

	int nIntDof = m_nIntPnt * 3;
	EigVec dN = EigVec::Random(nIntDof) * (noiseN * 2.0) - EigVec::Constant(nIntDof, noiseN);
	EigVec dP = EigVec::Random(m_nParam) * (noiseP * 2.0) - EigVec::Constant(m_nParam, noiseP);

	EigVec gn, gp;
	EigDense Hnp, Hpn, Hpp;
	EigSparse Hnn;

	computeGradient(n, p, m_t, gn, gp);
	computeHessian(n, p, m_t, Hnn, Hnp, Hpn, &Hpp);
	for (double step = 1; step > 1e-15; step *= 0.1)
	{
		EigVec dNi = dN * step;
		EigVec dPi = dP * step;

		// 计算近似的梯度值
		EigVec gni_app = gn + Hnn * dNi + Hnp * dPi;
		EigVec gpi_app = gp + Hpn * dNi + Hpp * dPi;

		// 计算准确的梯度值
		EigVec ni = n + dNi;
		EigVec pi = p + dPi;
		EigVec gni, gpi;
		computeGradient(ni, pi, m_t, gni, gpi);

		EigVec resN = gni - gni_app;
		EigVec resP = gpi - gpi_app;

		double resNNorm = resN.norm();
		double resPNorm = resP.norm();
		resN = resN.cwiseAbs();
		resP = resP.cwiseAbs();
		double resNMax  = resN.maxCoeff();
		double resPMax  = resP.maxCoeff();

		double invE = 1.0 / (step * step);
		PRINT_F("ε=%le  maxResN=%le  maxResP=%le  |maxN|/ε^2=%le  |maxP|/ε^2=%le",
			step, resNMax, resPMax, resNNorm * invE, resPNorm* invE);
	}
	return true;
}

bool RiggedMesh::testCurFrameGrad( RigStatus& lastFrame, RigStatus& curFrame, double noiseN, double noiseP)
{
	// 当前帧的配置
	PRINT_F("Test Gradient and Function Value");
	PRINT_F("noiseN = %lf, noiseP = %lf", noiseN, noiseP);
	setStatus(lastFrame);
	EigVec n, p;
	if(!getN(curFrame, n))
		return false;
	p = curFrame.getP();

	int nIntDof = m_nIntPnt * 3;
	EigVec dN = EigVec::Random(nIntDof) * (noiseN * 2) - EigVec::Constant(nIntDof, noiseN);
	EigVec dP = EigVec::Random(m_nParam)* (noiseP * 2) - EigVec::Constant(m_nParam, noiseP);
	EigVec gn, gp;
	computeGradient(n,p, m_t, gn, gp);
	double f0 = computeValue(n, p, m_t);
	PRINT_F("funcVal:%lf gradient: |gn| = %lf, |gp| = %lf", f0, gn.norm(), gp.norm());

	for (double step = 1; step > 1e-15; step *= 0.1)
	{
		EigVec dNi = dN * step;
		EigVec dPi = dP * step;

		// 计算近似的函数值
		double fi_app = f0 + gn.dot(dNi) + gp.dot(dPi);

		// 计算准确的函数值
		EigVec ni = n + dNi;
		EigVec pi = p + dPi;
		double fi = computeValue(ni, pi, m_t);

		double error = abs(fi-fi_app);
		PRINT_F("step = %le funVal = %le approxVal = %le error = %le error/dx^2 = %le", step, fi, fi_app, error, error/(step*step));
	}
	PRINT_F("\n");
	return true;
}

void RigFEM::RiggedMesh::getMeshPntPos( vector<double>& pnts ) const
{
	pnts.clear();
	if (m_nTotPnt == 0)
		return;
	pnts.resize(m_nTotPnt*3);

	for (int i = 0; i < m_nTotPnt; ++i)
	{
		Vec3d v = *m_tetMesh->getVertex(i);
		pnts[i*3+0] = v[0];
		pnts[i*3+1] = v[1];
		pnts[i*3+2] = v[2];
	}
}

bool RigFEM::RiggedMesh::computeStaticPos( const EigVec& p, double t, EigVec& q, int maxIter /*= 20*/, const EigVec* initQ /*= NULL*/ )
{
	if (p.size() != m_nParam || (initQ && initQ->size() != m_nTotPnt*3))
		return false;

	PRINT_F("compute static position");

	EigVec force(m_nTotPnt*3);
	q = initQ ? *initQ : EigVec::Zero(m_nTotPnt*3);
	computeQ(&p[0], t, &q[0]);

	for (int iter = 0; iter < maxIter; ++iter)
	{
		// 计算受力和切向刚度矩阵
		m_forceModel->GetForceAndMatrix(&q[0], &force[0], m_tangentStiffnessMatrix);
		//EigSparse tangentStiffnessMat;
		//Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_intDofIdx, m_intDofIdx, m_nIntDof, m_nIntDof, tangentStiffnessMat);
		EigSparse& tangentStiffnessMat = m_dFnn;
		Utilities::vegaSparse2AllocatedEigen(*m_tangentStiffnessMatrix, m_intDofIdx, m_intDofIdx, m_nIntDof, m_nIntDof, tangentStiffnessMat);

		// 求解新位置
		EigVec intForce(3 * m_nIntPnt);
		for (int ithIntPnt = 0; ithIntPnt < m_nIntPnt; ++ithIntPnt)
		{
			int idx = m_intPntIdx[ithIntPnt] * 3;
			intForce[ithIntPnt*3+0] = -force[idx+0];
			intForce[ithIntPnt*3+1] = -force[idx+1];
			intForce[ithIntPnt*3+2] = -force[idx+2];
		}
		Eigen::SuperLU<EigSparse> solver;
		solver.compute(tangentStiffnessMat);
		if(solver.info()!=Eigen::Success) 
		{
			PRINT_F("LU factorization FAILED!\n");
			return false;
		}
		EigVec dN = solver.solve(intForce);

		double intPntLengthSq = 0;				// 内部点向量长度平方
		for (int ithN = 0; ithN < m_nIntPnt; ++ithN)
		{
			int intIdx = m_intPntIdx[ithN];
			double* pDq = &dN[ithN*3];
			q[intIdx*3 + 0] += pDq[0];
			q[intIdx*3 + 1] += pDq[1];
			q[intIdx*3 + 2] += pDq[2];
			intPntLengthSq += pDq[0]*pDq[0] + pDq[1]*pDq[1] + pDq[2]*pDq[2];
		}
		intPntLengthSq = sqrt(intPntLengthSq);

		PRINT_F("%d th iteration, |dq| = %lf", iter, intPntLengthSq);
		if (intPntLengthSq < 1e-3)
		{
			break;
		}

	}
	return true;
}

const EigVec& RigFEM::RiggedMesh::getExternalForce() const
{
	return m_extForce;
}

bool RigFEM::RiggedMesh::setExternalForce( const EigVec&f )
{
	if (f.size() == m_nTotPnt*3)
	{
		m_extForce = f;
		return true;
	}
	return false;
}

const EigVec& RigFEM::RiggedMesh::getMass() const
{
	return m_mass;
}

void RigFEM::RiggedMesh::getVertexPosition( EigVec& pos )
{
	pos = m_q;

	for (int i = 0; i < m_nIntPnt; ++i)
	{
		int idx = m_intPntIdx[i];
		int idx3 = idx * 3;
		Vec3d v = *m_tetMesh->getVertex(idx);
		pos[idx3]   += v[0];
		pos[idx3+1] += v[1];
		pos[idx3+2] += v[2];
	}


	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx = m_surfPntIdx[i];
		int idx3 = idx * 3;
		Vec3d v = *m_tetMesh->getVertex(idx);
		pos[idx3]   += v[0];
		pos[idx3+1] += v[1];
		pos[idx3+2] += v[2];
	}
}

bool RigFEM::RiggedMesh::updateExternalAndControlForce()
{
	EigVec pos, extForce, surfForce;
	getVertexPosition(pos);
	bool res = m_rigObj->computeExternalForce(pos, m_v, m_mass, m_h, extForce, surfForce);

	if (!res || extForce.size() != m_nTotPnt * 3)
	{
		return false;
	}

	m_extForce = extForce;
	for (int i = 0; i < m_nSurfDof; ++i)
	{
		int idx = m_surfDofIdx[i];
		m_extForce[idx] += surfForce[i];
	}
	if (m_controlType == CONTROL_EXPLICIT_FORCE)
	{
		if(!computeControlForce(m_controlForce, m_param, m_paramVelocity))
		{
			m_controlForce.setZero(m_nParam);
		}
		Global::showVector(m_controlForce, "ctrlForce");
	}
	if (m_controlType != CONTROL_NONE)
	{
		m_rigObj->getControlGain(m_propGain, m_deriGain);
	}

	return true;
}

bool RigFEM::RiggedMesh::computeControlForce( 
	EigVec& generalForce, 
	const EigVec& param, const EigVec& paramVelocity)
{
	EigVec pGain, dGain;
	if(!m_rigObj->getControlGain(pGain, dGain))
		return false;

	int nParam = m_targetParam.size();
	if (param.size() != nParam || 
		paramVelocity.size() != nParam)
		return false;

	generalForce = pGain.cwiseProduct(m_targetParam - param);
	generalForce += dGain.cwiseProduct(m_targetParamVelocity - paramVelocity);
	return true;
}

bool RigFEM::RiggedMesh::computeControlForce( EigVec& generalForce, const EigVec& param )
{
	EigVec paramVel = (param - m_param) / m_h;
	return computeControlForce(generalForce, param, paramVel);
}

void RigFEM::RiggedMesh::setRigControlType( RigControlType type )
{
	m_controlType = type;
	m_controlForce.setZero(m_nTotPnt * 3);
	m_targetParam.setZero(m_nTotPnt*3);
	m_targetParamVelocity.setZero(m_nTotPnt*3);
}

void RigFEM::RiggedMesh::setControlTarget( const EigVec& targetParam, const EigVec& targetParamVelocity )
{
	m_targetParam = targetParam;
	m_targetParamVelocity = targetParamVelocity;
}

const EigVec& RigFEM::RiggedMesh::getControlTargetParam() const
{
	return m_targetParam;
}

const EigVec& RigFEM::RiggedMesh::getControlTargetVelocity() const
{
	return m_targetParamVelocity;
}

bool RigFEM::RiggedMesh::getControlTargetFromRigNode( EigVec& target )
{
	return m_rigObj->getControlTarget(target);
}

void RigFEM::RiggedMesh::computeReducedForceMatrix( const EigVec& q, const EigDense& T, EigDense& A )
{
	m_modelwrapper->computeReducedForceMatrix(q.data(), T, A);
}

bool RigFEM::RiggedMesh::loadElementMaterialFactor( EigVec& factor )
{
	return m_modelwrapper->setElementMaterialFactor(factor);
}

bool RigFEM::RiggedMesh::clearElementMaterialFactor()
{
	m_modelwrapper->clearElementMaterialFactor();
	return true;
}

bool RigFEM::RiggedMesh::computeControlEnergy( const EigVec& param, double& energy )
{
	EigVec pGain, dGain;
	if(!m_rigObj->getControlGain(pGain, dGain))
		return false;

	int nParam = m_targetParam.size();
	if (param.size() != nParam)
		return false;

	EigVec quadCoef(nParam);
	EigVec k0(nParam);
	for (int i = 0; i < nParam; ++i)
	{
		k0[i] = pGain[i] * m_targetParam[i] + 
			dGain[i] * m_targetParamVelocity[i] +
			dGain[i] * m_param[i] / m_h;

		quadCoef[i] = pGain[i] + dGain[i] / m_h;
	}

	energy = 0;
	for (int i = 0; i < nParam; ++i)
	{
		energy += 0.5 * param[i] * param[i] * quadCoef[i] - k0[i] * param[i];
	}
	return true;
}

bool RigFEM::RiggedMesh::computeElementLaplacian( vector<NeighIdx>& adjList )
{
	if (!m_tetMesh)
		return false;
	int nEle = m_tetMesh->getNumElements();
	if (nEle <= 0)
		return false;

	map<OrderedPair, vector<int>> neighMap;
	for (int ithEle = 0; ithEle < nEle; ++ithEle)
	{
		int vtxID[4] = {
			m_tetMesh->getVertexIndex(ithEle,0),
			m_tetMesh->getVertexIndex(ithEle,1),
			m_tetMesh->getVertexIndex(ithEle,2),
			m_tetMesh->getVertexIndex(ithEle,3)
		};

		std::sort(vtxID,&vtxID[0]+4);
		neighMap[OrderedPair(vtxID[0],vtxID[1],vtxID[2])].push_back(ithEle);
		neighMap[OrderedPair(vtxID[1],vtxID[2],vtxID[3])].push_back(ithEle);
		neighMap[OrderedPair(vtxID[0],vtxID[2],vtxID[3])].push_back(ithEle);
		neighMap[OrderedPair(vtxID[0],vtxID[1],vtxID[3])].push_back(ithEle);
	}

	adjList.clear();
	adjList.resize(nEle);
	for (map<OrderedPair, vector<int>>::iterator pFace = neighMap.begin();
		pFace != neighMap.end(); ++pFace)
	{
		vector<int>& eleID = pFace->second;
		if (eleID.size() < 2)
			continue;
		if (eleID.size() != 2)
		{
			PRINT_F("wrong face has %d neighbour", eleID.size());
			continue;
		}

		int ele1 = eleID[0];
		int ele2 = eleID[1];
		adjList[ele1].push_back(ele2);
		adjList[ele2].push_back(ele1);
	}
	return false;
}

const vector<int> RigFEM::RiggedMesh::getInternalPntIdx() const
{
	vector<int> intPntIdx(m_nIntPnt);
	std::copy(m_intPntIdx, m_intPntIdx+m_nIntPnt, intPntIdx.begin());
	return intPntIdx;
}

const vector<int> RigFEM::RiggedMesh::getSurfacePntIdx() const
{
	vector<int> surfPntIdx(m_nSurfPnt);
	std::copy(m_surfPntIdx, m_surfPntIdx+m_nSurfPnt, surfPntIdx.begin());
	return surfPntIdx;
}

bool RigFEM::RiggedMesh::allocateTempData()
{
	m_forceModel->GetTangentStiffnessMatrixTopology(&m_tangentStiffnessMatrix);
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_surfDofIdx, m_nSurfDof, m_nSurfDof, m_dFss);
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_surfDofIdx, m_intDofIdx,  m_nSurfDof,  m_nIntDof, m_dFsn);
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_intDofIdx,  m_surfDofIdx, m_nIntDof,  m_nSurfDof, m_dFns);
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, m_intDofIdx,  m_intDofIdx,  m_nIntDof,  m_nIntDof,  m_dFnn);

	m_intDofBuf.resize(m_nIntDof);
	m_surfDofBuf.resize(m_nSurfDof);
	m_totDofBuf.resize(m_nTotPnt*3);
	return true;
}

bool RigFEM::RiggedMesh::freeTempData()
{
	delete m_tangentStiffnessMatrix;
	m_tangentStiffnessMatrix = NULL;
	m_dFnn = m_dFns = m_dFsn = m_dFss = EigSparse();

	m_intDofBuf = m_surfDofBuf = m_totDofBuf = EigVec();
	return true;
}


void RigFEM::FEMSystem::init()
{
	m_mesh.init();
	m_solver.setMesh(&m_mesh);
}

void RigFEM::MeshDispConfig::setDefault()
{
	m_extForceDispFactor = 0;
	
	m_showEdge = true;
	m_showPoint = true;
	m_showBBox = true;
	m_showMaterial = false;

	m_minMaterialFactor = 0;
	m_maxMaterialFactor = 2;
	m_axis = AXIS_X;
	m_axisFactor = 1;
}

RigFEM::MeshDispConfig::MeshDispConfig()
{
	setDefault();
}

void RigFEM::MeshDispConfig::endColorMaterial()
{
	glPopAttrib();
}

void RigFEM::MeshDispConfig::beginColorMaterial()
{
	float zeros[] = {0,0,0,0};

	glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT);

 	glEnable(GL_LIGHTING); // 使用灯光

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);  

	GLfloat specular[] = {0,0,0,0};
	GLfloat shininess[]  = {0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

bool RigFEM::RiggedMesh::OrderedPair::operator<( const OrderedPair& o ) const
{
	if (m_k1 != o.m_k1)
		return m_k1 < o.m_k1;
	if (m_k2 != o.m_k2)
		return m_k2 < o.m_k2;
	return m_k3 < o.m_k3;
}
