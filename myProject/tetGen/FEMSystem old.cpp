#include "StdAfx.h"
#include "FEMSystem.h"

using namespace RigFEM;

RiggedMesh::RiggedMesh(void):m_tetMesh(NULL),  m_forceModel(NULL), m_h(0.1), m_t(0), m_tangentStiffnessMatrix(NULL),  m_modelwrapper(NULL), m_nTotPnt(0), m_nIntPnt(0), m_nSurfPnt(0), m_nParam(0)
{
}

RiggedMesh::~RiggedMesh(void)
{
}

void RiggedMesh::init()
{
	char*  myArgv[] = {
		"i:/Programs/VegaFEM-v2.1/myProject/tetGen/Debug/tetGen.exe",
		"-pq1.3a0.1m",
		"model/torus.off"
	};
	int myArgc = sizeof(myArgv)/4;

	tetgenbehavior b;

	if (!b.parse_commandline(myArgc, myArgv)) {
		terminatetetgen(NULL, 10);
	}

	// Read input files.
	if (b.refine) { // -r
		if (!m_in.load_tetmesh(b.infilename, (int) b.object)) {
			terminatetetgen(NULL, 10);
		}
	} else { // -p
		if (!m_in.load_plc(b.infilename, (int) b.object)) {
			terminatetetgen(NULL, 10);
		}
	}
	if (b.insertaddpoints) { // -i
		// Try to read a .a.node file.
		m_addin.load_node(b.addinfilename);
	}
	if (b.metric) { // -m
		// Try to read a background mesh in files .b.node, .b.ele.
		m_bgmin.load_tetmesh(b.bgmeshfilename, (int) b.object);
	}

	tetrahedralize(&b, &m_in, &m_out, &m_addin, &m_bgmin);

	// 找出原来的点
	findOriPoints();

	// 设置rig初始点
	double* oriPnts = new double[m_nSurfPnt*3];
	double* p = oriPnts;
	for (int i = 0; i < m_nSurfPnt; ++i, p+=3)
	{
		double* pnt = m_out.pointlist + m_surfPntIdx[i]*3;
		p[0] = pnt[0];
		p[1] = pnt[1];
		p[2] = pnt[2];
	}
	m_transRig.setInitPnts(oriPnts, m_nSurfPnt);
	delete[] oriPnts;

	// 构建四面体网格
	buildTetMesh();

	// 计算当前参数下的配置
	double param[9] = {0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0};
	m_transRig.setFreeParam(param);

	m_transRig.keyParam(0, MathUtilities::zero);
	m_transRig.keyParam(1, MathUtilities::zero);
	m_transRig.keyParam(2, MathUtilities::absSin);
	m_transRig.keyParam(3, MathUtilities::zero);
	m_transRig.keyParam(4, MathUtilities::zero);
	m_transRig.keyParam(5, MathUtilities::zero);
	m_nParam = m_transRig.getNFreeParam();
	computeRig();
}

void RiggedMesh::show()
{
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

	glColor3f(0.8, 0.8,0.8);
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

void RiggedMesh::findOriPoints()
{
	std::set<PntPair> pntSet;

	double* p = m_in.pointlist;
	for (int i = 0; i < m_in.numberofpoints; ++i, p+=3)
	{
		PntPair pair = {Vec3d (p[0], p[1], p[2]), i};
		pntSet.insert(pair);
	}

	p = m_out.pointlist;
	for (int i = 0; i < m_out.numberofpoints; ++i, p+=3)
	{ 
		PntPair pair = {Vec3d (p[0], p[1], p[2]), i};
		std::set<PntPair>::iterator pPair = pntSet.find(pair);
		if (pPair != pntSet.end())
		{
			m_surfPntIdx.push_back(i);
			pntSet.erase(pPair);
		}
		else
			m_intPntIdx.push_back(i);
	}
	m_nIntPnt = m_intPntIdx.size();
	m_nSurfPnt = m_surfPntIdx.size();
	m_nTotPnt = m_out.numberofpoints;
}

bool RiggedMesh::PntPair::operator<( const PntPair& other ) const
{
	if (m_pnt[0] != other.m_pnt[0])
		return m_pnt[0] < other.m_pnt[0];
	if (m_pnt[1] != other.m_pnt[1])
		return m_pnt[1] < other.m_pnt[1];
	return m_pnt[2] < other.m_pnt[2];
}

bool RigFEM::RiggedMesh::buildTetMesh()
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
	m_tetMesh = new TetMesh(nPnts, m_out.pointlist, numTet, tetIDs, 10000, 0.45);
	delete[] tetIDs;


	int wrap = 2;
	CorotationalLinearFEMWrapper* wrapper = new CorotationalLinearFEMWrapper(m_tetMesh, wrap);
	m_modelwrapper    = wrapper;
	m_forceModel      = new CorotationalLinearFEMForceModel(wrapper, wrap);
	
/*
	NeoHookeanIsotropicMaterial* hookeanMat = new NeoHookeanIsotropicMaterial(m_tetMesh, 1, 0.5);
	IsotropicHyperelasticFEM* defModel = new IsotropicHyperelasticFEM(m_tetMesh, hookeanMat);
	m_forceModel = new IsotropicHyperelasticFEMForceModel(defModel);*/

	m_q = m_v = m_a = EigVec::Constant(nPnts*3, 0.0);
	m_mass = EigVec(nPnts);
	GenerateMassMatrix::computeVertexMasses(m_tetMesh, &m_mass[0]);
	if (!m_tangentStiffnessMatrix)
		m_forceModel->GetTangentStiffnessMatrixTopology(&m_tangentStiffnessMatrix);
	return true;
}

void RigFEM::RiggedMesh::computeRig()
{
	EigVec res(m_nSurfPnt * 3);
	m_param.setZero(m_transRig.getNFreeParam());
	m_transRig.getFreeParam(&m_param[0]);
	m_transRig.computeValue(&res[0]);

	m_q.resize(m_nTotPnt*3);
	double* pQ = &m_q[0];
	for (int i = 0; i < m_nTotPnt * 3; ++i)
		pQ[i] = 0;

	const vector<Vec3d>& pntArray = m_transRig.getCurPnt();
	const vector<Vec3d>& initPntArray = m_transRig.getInitPnt();

	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		Vec3d d = pntArray[i] - initPntArray[i];
		int idx = m_surfPntIdx[i];
		pQ[idx*3] = d[0];
		pQ[idx*3+1] = d[1];
		pQ[idx*3+2] = d[2];
	}

	m_force.resize(m_nTotPnt*3);
	m_forceModel->GetInternalForce(&m_q[0], &m_force[0]);

	/*
	// 测试tangent stiffness matrix	
	vector<double> dq(nTotPnt*3);
	vector<double> q(nTotPnt*3);
	for (int i = 0; i < nTotPnt*3; ++i)
	{
		q[i]  = rand() / 65535.0 * 1;
		dq[i] = rand() / 65535.0 * 2;
	}
	m_forceModel->TestStiffnessMatrix(&q[0], &dq[0]);*/
	

	//buildHessians(&m_q[0]);
	/*
	vector<double> d2sdpidpj(nSurfPnt*3);
	for (int i = 0; i < m_transRig.getNParam(); ++i)
	{
		for (int j = 0; j < m_transRig.getNParam(); ++j)
		{
			m_transRig.computeJacobianDerivative(i,j, &d2sdpidpj[0]);
			printf("(%d %d)\n", i,j);
			for (int k = 0; k < d2sdpidpj.size(); ++k)
			{
				printf("%lf\t", d2sdpidpj[k]);
			}
			printf("\n");
		}
	}*/
}

bool RigFEM::RiggedMesh::computeHessian(const EigVec& n, const EigVec& p, double t, EigSparse& Hnn, EigDense& Hnp, EigDense& Hpn, EigDense* pHpp)
{
	// 计算指定参数下各个点的偏移量q
	EigVec q(m_nTotPnt*3);
	computeQ(&n[0], &p[0], t, &q[0]);

	// 提取内力、tangent Stiffness matrix、质量矩阵
	// 注意tangent stiffness matrix为实际的负值，因为系统计算出的弹力为实际的负值,因此需要先反转
	EigVec force(m_nTotPnt*3);
	if(!m_tangentStiffnessMatrix)
		m_forceModel->GetTangentStiffnessMatrixTopology(&m_tangentStiffnessMatrix);
	m_forceModel->GetForceAndMatrix(&q[0], &force[0], m_tangentStiffnessMatrix);
	*m_tangentStiffnessMatrix *= -1;
	force *= -1;

	// 提取 dFn/dn dFn/ds dFs/dn dFs/ds,
	// 其中Fn为内部节点受到的力，Fs为表面节点受到的力
	// n为内部节点位置，s为表面节点位置
	vector<int> nIdx(m_nIntPnt * 3), sIdx(m_nSurfPnt * 3);
	for (int i = 0; i < m_nIntPnt; ++i)
	{
		nIdx[i*3]   = 3 * m_intPntIdx[i];
		nIdx[i*3+1] = 3 * m_intPntIdx[i] + 1;
		nIdx[i*3+2] = 3 * m_intPntIdx[i] + 2;
	}
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		sIdx[i*3]   = 3 * m_surfPntIdx[i];
		sIdx[i*3+1] = 3 * m_surfPntIdx[i] + 1;
		sIdx[i*3+2] = 3 * m_surfPntIdx[i] + 2;
	}
	SparseMatrix dFnn(*m_tangentStiffnessMatrix);
	SparseMatrix dFns(*m_tangentStiffnessMatrix);
	SparseMatrix dFsn(*m_tangentStiffnessMatrix);
	SparseMatrix dFss(*m_tangentStiffnessMatrix);
	dFnn.RemoveRowsColumns(sIdx.size(), &sIdx[0]);
	dFss.RemoveRowsColumns(nIdx.size(), &nIdx[0]);
	dFns.RemoveRows(sIdx.size(), &sIdx[0]);
	dFns.RemoveColumns(nIdx.size(), &nIdx[0]);
	dFsn.RemoveRows(nIdx.size(), &nIdx[0]);
	dFsn.RemoveColumns(sIdx.size(), &sIdx[0]);

	// 构建矩阵Hnn = Mn / h^2 - dFn/dn
	Utilities::vegaSparse2Eigen(dFnn, Hnn, m_nIntPnt*3);
	EigSparse Hnn_e;
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, nIdx, nIdx, Hnn_e);
	double maxError = Utilities::maxError(Hnn, Hnn_e);
	printf("Hnn maxError = %lf\n", maxError);
	Hnn *= -1.0;
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
	m_transRig.computeJacobian(J);
	EigSparse dFns_e;
	Utilities::vegaSparse2Eigen(dFns, dFns_e, m_nSurfPnt*3);
	EigSparse dFnsee;
	Utilities::vegaSparse2Eigen(*m_tangentStiffnessMatrix, nIdx, sIdx, dFnsee);
	maxError = Utilities::maxError(dFns_e, dFnsee);
	printf("dfns maxError = %lf\n", maxError);
	Hnp = -1.f * dFns_e * J;

	// 构建矩阵Hpn = - J' * dFs/dn
	EigSparse dFsn_e;
	Utilities::vegaSparse2Eigen(dFsn, dFsn_e, m_nIntPnt*3);
	Hpn = - J.transpose() * dFsn_e;

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
		EigSparse dFss_e;
		Utilities::vegaSparse2Eigen(dFss, dFss_e, m_nSurfPnt*3);
		for (int i = 0; i < m_nParam; ++i)
		{
			for (int l = 0; l < m_nParam; ++l)
			{
				double& Hppil = Hpp(i,l);
				m_transRig.computeJacobianDerivative(i,l, &jacobianDeri[0]);
				for (int k = 0; k < m_nSurfPnt*3; ++k)
				{
					Hppil += jacobianDeri[k] * sResi[k];

					double m = m_mass[m_surfPntIdx[k/3]];
					double v = 0;
					for (int m = 0; m < m_nSurfPnt*3; ++m)
					{
						v += dFss_e.coeff(k,m) * J(m,l);
					}
					Hppil += J(k,i) *(m * J(k,l) * invH2 - v);
				}
			}
		}
	}

	return true;
}

void RigFEM::RiggedMesh::testStep( double dt )
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

void RigFEM::RiggedMesh::computeQ( const double* n, const double* p, double t, double* q )
{
	EigVec s;
	computeSurfOffset(p, t, s);
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
}

void RigFEM::RiggedMesh::computeQ( const double* p, double t, double* q )
{
	EigVec s;
	computeSurfOffset(p, t, s);
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx = m_surfPntIdx[i];
		q[idx*3]   = s[i*3];
		q[idx*3+1] = s[i*3+1];
		q[idx*3+2] = s[i*3+2];
	}
}

void RigFEM::RiggedMesh::computeSurfOffset( const double* p, double t, EigVec& s )
{
	s = EigVec(m_nSurfPnt*3);
	m_transRig.setTime(t);
	m_transRig.computeValue(&s[0], p);
	for (int i = 0; i < m_nSurfPnt; ++i)
	{
		int idx3 = m_surfPntIdx[i] * 3;
		s[i*3]   -= m_out.pointlist[idx3++];
		s[i*3+1] -= m_out.pointlist[idx3++];
		s[i*3+2] -= m_out.pointlist[idx3];
	}
}

bool RigFEM::RiggedMesh::computeGradient( const EigVec& n, const EigVec& p, double t, EigVec& gn, EigVec& gs )
{
	// 计算给定参数下状态
	EigVec q(m_nTotPnt*3);
	computeQ(&n[0], &p[0], t, &q[0]);

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
	EigVec residual = ma - f;

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
	m_transRig.computeJacobian(J);
	gs = J.transpose() * gs;
	return true;
}

bool RigFEM::RiggedMesh::testHessian()
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
	printf("test Hessian\n");
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
		printf("ε=%le  maxResN=%le  maxResP=%le  |maxN|/ε^2=%le  |maxP|/ε^2=%le\n",
			step, resNMax, resPMax, resNMax * invE, resPMax* invE);
	}
	return true;
}






bool RigFEM::RiggedMesh::testElasticEnergy()
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
	printf("test result\n");
	for (int ithDof = 0; ithDof < nDof; ++ithDof)
	{
		printf("f=%lf \t fRef=%lf \t error = %lf \n", 
			f[ithDof], fRef[ithDof], error[ithDof]);
	}
	printf("\n");
	printf("fNorm = %lf \t fRefNorm = %lf errorNorm = %lf\n", f.norm(), fRef.norm(), error.norm());
	printf("errorMax = %lf\n", error.cwiseAbs().maxCoeff());

	
	m_q = q;
	return true;
}

double RigFEM::RiggedMesh::computeValue( const EigVec& n, const EigVec& p, double t)
{

	// 计算给定参数下状态
	EigVec q(m_nTotPnt*3);
	computeQ(&n[0], &p[0], t, &q[0]);

	return computeValue(q);
}

double RigFEM::RiggedMesh::computeValue( const EigVec& q )
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
	return elasticEnergy + kinetic;
}

bool RigFEM::RiggedMesh::testValue()
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
	printf("test func val\n");
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
		printf("step = %le funVal = %le approxVal = %le error = %le error/dx^2 = %le\n", step, fi, fi_app, error, error/(step*step));
	}
	printf("\n");
	return true;
}

bool RigFEM::RiggedMesh::computeValueAndGrad(const EigVec& x, const EigVec& param,  double* v, EigVec* grad)
{
	if (!v && !grad)
		return false;

	double t = param[0];
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
	computeQ(&n[0], &p[0], t, &q[0]);

	// 计算函数值
	if (v)
	{
		*v = computeValue(q);
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

		// 计算此状态下的力
		EigVec f(m_nTotPnt*3);
		m_forceModel->GetInternalForce(&q[0], &f[0]);
		f *= -1.0;
		EigVec residual = ma - f;

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
		m_transRig.computeJacobian(J);
		gs = J.transpose() * gs;

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
	return true;
}

void RigFEM::RiggedMesh::setDof( EigVec&n, EigVec&p, bool proceedTime /*= true*/ )
{
	EigVec q(m_nTotPnt*3);
	computeQ(&n[0], &p[0], m_t, &q[0]);
	EigVec v = (q - m_q) / m_h;
	EigVec a = (v - m_v) / m_h; 
	m_param = p;
	m_v = v;
	m_a = a;
	m_q = q;
	if (proceedTime)
	{
		m_t += m_h;
	}
}

void RigFEM::RiggedMesh::getDof( EigVec& n, EigVec& p )
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




void RigFEM::FEMSystem::clearResult()
{
	m_solver.clearResult();
}

void RigFEM::FEMSystem::saveResult( const char* fileName )
{
	m_solver.saveResult(fileName);
}

void RigFEM::FEMSystem::step()
{
	m_solver.step();
}
