#include "StdAfx.h"
#include "CorotationalLinearFEMWrapper.h"


CorotationalLinearFEMWrapper::~CorotationalLinearFEMWrapper(void)
{
}

double CorotationalLinearFEMWrapper::computeElasticEnergy(const double* u)
{
	int numElements = tetMesh->getNumElements();
	double energy = 0;
	for (int el=0; el < numElements; el++)
	{
		double materialFactor = 1;
		if (el < m_eleMatFactor.size())
		{
			materialFactor = m_eleMatFactor[el];
		}

		int vtxIndex[4];
		for (int vtx=0; vtx<4; vtx++)
			vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);

		double KElement[144]; // element stiffness matrix, to be computed below; row-major

		// 线性弹力模型
		if (m_wrap == 0)
		{
			// 找到单元刚度矩阵
			memcpy(KElement, KElementUndeformed[el], sizeof(double) * 144);

			// f = K u 计算受力
			double fElement[12];
			for(int i=0; i<12; i++)
			{
				fElement[i] = 0;
				for(int j=0; j<4; j++)
				{
					// 乘以缩放材料硬度的因子
					fElement[i] += (
						KElement[12 * i + 3 * j + 0] * u[3 * vtxIndex[j] + 0] +
						KElement[12 * i + 3 * j + 1] * u[3 * vtxIndex[j] + 1] +
						KElement[12 * i + 3 * j + 2] * u[3 * vtxIndex[j] + 2]) * materialFactor;

				}
			}

			// 单元弹性能量E = 1/2 * u^T * K * u = 1/2 * u^T * f 
			double eleEnergy = 0.0;
			for (int ithVtx = 0; ithVtx < 4; ++ithVtx)
			{
				eleEnergy += fElement[ithVtx*3 + 0] * u[3*vtxIndex[ithVtx] + 0];
				eleEnergy += fElement[ithVtx*3 + 1] * u[3*vtxIndex[ithVtx] + 1];
				eleEnergy += fElement[ithVtx*3 + 2] * u[3*vtxIndex[ithVtx] + 2];
			}
			eleEnergy /= 2;
			energy += eleEnergy;
		}
		else if (m_wrap > 0)
		{
			double P[16]; // the current world-coordinate positions (row-major)
			/*
			P = [ v0   v1   v2   v3 ]
			[  1    1    1    1 ]
			*/
			// rows 1,2,3
			for(int i=0; i<3; i++)
				for(int j=0; j<4; j++)
					P[4 * i + j] = undeformedPositions[3 * vtxIndex[j] + i] + u[3 * vtxIndex[j] + i];
			// row 4
			for(int j=0; j<4; j++)
				P[12 + j] = 1;

			// F = P * Inverse(M)
			double F[9]; // upper-left 3x3 block
			for(int i=0; i<3; i++) 
				for(int j=0; j<3; j++) 
				{
					F[3 * i + j] = 0;
					for(int k=0; k<4; k++)
						F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
				}

				double R[9]; // rotation (row-major)
				double S[9]; // symmetric (row-major)
				double tolerance = 1E-15;
				int forceRotation = 1;
				PolarDecomposition::Compute(F, R, S, tolerance, forceRotation);

#if (ENERGY_METHOD == ENERGY_FROM_STIFFNESSMAT)
				// 方法1:抵消旋转后用刚度矩阵计算

				// 计算 dx = R^T * x - x0
				double rotOffset[12];
				for (int ithPnt = 0; ithPnt < 4; ++ithPnt)
				{
					int pntIdx = ithPnt * 3;
					double pnt[] = {P[0 + ithPnt], P[4 + ithPnt], P[8 + ithPnt]};

					rotOffset[pntIdx+0] = R[0] * pnt[0] + R[3] * pnt[1] + R[6] * pnt[2] - undeformedPositions[3 * vtxIndex[ithPnt] + 0];
					rotOffset[pntIdx+1] = R[1] * pnt[0] + R[4] * pnt[1] + R[7] * pnt[2] - undeformedPositions[3 * vtxIndex[ithPnt] + 1];
					rotOffset[pntIdx+2] = R[2] * pnt[0] + R[5] * pnt[1] + R[8] * pnt[2] - undeformedPositions[3 * vtxIndex[ithPnt] + 2];
				}

				// 计算f = K * dx
				double f[12];
				memcpy(KElement, KElementUndeformed[el], sizeof(double) * 144);
				for(int ithDof = 0; ithDof < 12; ++ithDof)
				{
					f[ithDof] = 0;
					for (int jthDof = 0; jthDof < 12; ++jthDof)
			  {
				  // 乘以缩放因子
				  f[ithDof] += KElement[ithDof*12 + jthDof] * rotOffset[jthDof] * materialFactor;
			  }
				}

				// 计算弹性能量 E= 1/2 * dx^T * f
				double eleEnergy = 0;
				for(int ithDof = 0; ithDof < 12; ++ithDof)
				{
					eleEnergy += f[ithDof] * rotOffset[ithDof];
				}
				eleEnergy /= 2.0;
				energy += eleEnergy;

#elif ENERGY_METHOD == ENERGY_FROM_MODEL

				// 方法2:直接代入corotational模型的能量密度公式计算
				VolumetricMesh::Material * material = tetMesh->getElementMaterial(el);
				VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
				if (eNuMaterial == NULL)
					continue;

				double lambda = eNuMaterial->getLambda() * materialFactor;
				double mu     = eNuMaterial->getMu()     * materialFactor;
				double volume = tetMesh->getElementVolume(el);
				/*
				Vec3d v0(P[0],P[4],P[8]);
				Vec3d v1(P[1],P[5],P[9]);
				Vec3d v2(P[2],P[6],P[10]);
				Vec3d v3(P[3],P[7],P[11]);
				TetMesh::getTetVolume(&v0, &v1, &v2, &v3);
				*/

				S[0] -= 1.0;
				S[4] -= 1.0;
				S[8] -= 1.0;

				double ss = 0;
				for (int i = 0; i < 9; ++i)
					ss += S[i] * S[i];

				double trace = S[0] + S[4] + S[8];
				double eleEnergy2 = (mu * ss + lambda * 0.5 * trace * trace) * volume;
				//printf("error %lf%%\n", (eleEnergy - eleEnergy2) / eleEnergy);

				energy += eleEnergy2;
				if (abs(eleEnergy2 - eleEnergy) / eleEnergy > 0.01)
				{
					PRINT_F("%dth element: energy1 = %lf energy2 = %lf error = %lf\n", el, eleEnergy, eleEnergy2, eleEnergy2 - eleEnergy);
					for (int i = 0; i < 9; ++i)
			  {
				  PRINT_F("%lf ", S[i]);
			  }
					PRINT_F("\n");
				}
#endif
		}
	}
	return energy;
}

void CorotationalLinearFEMWrapper::computeReducedForceMatrix(const double * vertexDisplacements, 
															 const EigDense& reduceMat,
															 EigDense& reducedEleForceMat)
{
	// 初始化参数
	const double* u = vertexDisplacements;
	int elementLo = 0;
	int elementHi = tetMesh->getNumElements();
	int nReducedDim = reduceMat.rows();
	int warp = m_wrap;

	// 按照vega原来的方法计算内力
	EigVec force;
	reducedEleForceMat.resize(nReducedDim, elementHi - elementLo);
	double* reducedEleForceData = reducedEleForceMat.data();
	for (int el=elementLo; el < elementHi; el++)
	{
		force.setZero(3*numVertices);

		int vtxIndex[4];
		for (int vtx=0; vtx<4; vtx++)
			vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);

		double KElement[144]; // element stiffness matrix, to be computed below; row-major

		if (warp > 0)
		{
			double P[16]; // the current world-coordinate positions (row-major)
			/*
			P = [ v0   v1   v2   v3 ]
			[  1    1    1    1 ]
			*/
			// rows 1,2,3
			for(int i=0; i<3; i++)
				for(int j=0; j<4; j++)
					P[4 * i + j] = undeformedPositions[3 * vtxIndex[j] + i] + u[3 * vtxIndex[j] + i];
			// row 4
			for(int j=0; j<4; j++)
				P[12 + j] = 1;

			// F = P * Inverse(M)
			double F[9]; // upper-left 3x3 block
			for(int i=0; i<3; i++) 
			{
				for(int j=0; j<3; j++) 
				{
					F[3 * i + j] = 0;
					for(int k=0; k<4; k++)
						F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
				}
			}

			double R[9]; // rotation (row-major)
			double S[9]; // symmetric (row-major)
			double tolerance = 1E-6;
			int forceRotation = 1;
			PolarDecomposition::Compute(F, R, S, tolerance, forceRotation);

			// RK = R * K
			// KElement = R * K * R^T
			double RK[144]; // row-major
			WarpMatrix(KElementUndeformed[el], R, RK, KElement);

			// f = RK (RT x - x0)
			double fElement[12];
			for(int i=0; i<12; i++)
			{
				fElement[i] = 0;
				for(int j=0; j<4; j++)
					for(int l=0; l<3; l++)
						fElement[i] += KElement[12 * i + 3 * j + l] * P[4 * l + j] - RK[12 * i + 3 * j + l] * undeformedPositions[3 * vtxIndex[j] + l];
			}

			// add fElement into the global f
			for(int j=0; j<4; j++)
				for(int l=0; l<3; l++)
					force[3 * vtxIndex[j] + l] += fElement[3 * j + l];
		}
		else
		{
			// no warp
			memcpy(KElement, KElementUndeformed[el], sizeof(double) * 144);
			// f = K u
			double fElement[12];
			for(int i=0; i<12; i++)
			{
				fElement[i] = 0;
				for(int j=0; j<4; j++)
				{
					fElement[i] += 
						KElement[12 * i + 3 * j + 0] * u[3 * vtxIndex[j] + 0] +
						KElement[12 * i + 3 * j + 1] * u[3 * vtxIndex[j] + 1] +
						KElement[12 * i + 3 * j + 2] * u[3 * vtxIndex[j] + 2];
				}
			}

			// add fElement into the global f
			for(int j=0; j<4; j++)
			{
				force[3 * vtxIndex[j] + 0] += fElement[3 * j + 0];
				force[3 * vtxIndex[j] + 1] += fElement[3 * j + 1];
				force[3 * vtxIndex[j] + 2] += fElement[3 * j + 2];
			}
		}

		EigVec reducedElementForce = reduceMat * force;
		for (int i = 0; i < nReducedDim; ++i)
		{
			*reducedEleForceData = reducedElementForce[i];
			reducedEleForceData++;
		}
	}
}

bool CorotationalLinearFEMWrapper::setElementMaterialFactor( EigVec& factor )
{
	if (factor.size() == tetMesh->getNumElements())
	{
		m_eleMatFactor = factor;
		return true;
	}
	return false;
}

double CorotationalLinearFEMWrapper::getElementMaterialFactor( int ithElement )
{
	if (ithElement >= 0 && ithElement < m_eleMatFactor.size())
	{
		return m_eleMatFactor[ithElement];
	}
	return 1;
}

void CorotationalLinearFEMWrapper::ComputeForceAndStiffnessMatrixOfSubmesh( double * u, double * f, SparseMatrix * stiffnessMatrix, int warp, int elementLo, int elementHi )
{
	// clear f to zero
	if (f != NULL)
		memset(f, 0, sizeof(double) * 3 * numVertices);

	// clear stiffness matrix to zero
	if (stiffnessMatrix != NULL)
		stiffnessMatrix->ResetToZero();

	for (int el=elementLo; el < elementHi; el++)
	{

		double matFactor = 1;
		if (m_eleMatFactor.size() == tetMesh->getNumElements())
		{
			matFactor = m_eleMatFactor[el];
		}

		int vtxIndex[4];
		for (int vtx=0; vtx<4; vtx++)
			vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);

		double KElement[144]; // element stiffness matrix, to be computed below; row-major

		if (warp > 0)
		{
			double P[16]; // the current world-coordinate positions (row-major)
			/*
			P = [ v0   v1   v2   v3 ]
			[  1    1    1    1 ]
			*/
			// rows 1,2,3
			for(int i=0; i<3; i++)
				for(int j=0; j<4; j++)
					P[4 * i + j] = undeformedPositions[3 * vtxIndex[j] + i] + u[3 * vtxIndex[j] + i];
			// row 4
			for(int j=0; j<4; j++)
				P[12 + j] = 1;

			// F = P * Inverse(M)
			double F[9]; // upper-left 3x3 block
			for(int i=0; i<3; i++) 
				for(int j=0; j<3; j++) 
				{
					F[3 * i + j] = 0;
					for(int k=0; k<4; k++)
						F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
				}

				double R[9]; // rotation (row-major)
				double S[9]; // symmetric (row-major)
				double tolerance = 1E-6;
				int forceRotation = 1;
				PolarDecomposition::Compute(F, R, S, tolerance, forceRotation);

				// RK = R * K
				// KElement = R * K * R^T
				double RK[144]; // row-major
				WarpMatrix(KElementUndeformed[el], R, RK, KElement);

				// f = RK (RT x - x0)
				double fElement[12];
				for(int i=0; i<12; i++)
				{
					fElement[i] = 0;
					for(int j=0; j<4; j++)
						for(int l=0; l<3; l++)
							fElement[i] += KElement[12 * i + 3 * j + l] * P[4 * l + j] - RK[12 * i + 3 * j + l] * undeformedPositions[3 * vtxIndex[j] + l];
				}

				// add fElement into the global f
				if (f != NULL)
				{
					for(int j=0; j<4; j++)
						for(int l=0; l<3; l++)
							f[3 * vtxIndex[j] + l] += fElement[3 * j + l] * matFactor;
				}

				// compute exact stiffness matrix
				if (warp == 2)
				{
					// compute G = (tr(S) I - S) R^T
					double G[9]; 
					double tr = S[0] + S[4] + S[8];
					double temp[9];
					for(int i=0; i<9; i++)
						temp[i] = -S[i];
					temp[0] += tr;
					temp[4] += tr;
					temp[8] += tr;
					// G = temp * R^T
					MATRIX_MULTIPLY3X3ABT(temp, R, G);

					double invG[9]; // invG = G^{-1}
					inverse3x3(G, invG);

					double rhs[27]; // 3 x 9 matrix (column-major)
					for(int i=0; i<3; i++)
						for(int j=0; j<3; j++)
						{
							double temp[9];
							for(int k=0; k<9; k++)
								temp[k] = 0.0;
							// copy i-th row of R into column j of temp      
							for(int k=0; k<3; k++)
								temp[3 * k + j] = R[3 * i + k];
							// extract the skew-symmetric part
							SKEW_PART(temp, &rhs[3 * (3 * i + j)]);
						}
						// must undo division by 2 from inside the SKEW_PART macro
						for(int i=0; i<27; i++)
							rhs[i] *= 2.0;

						// solve G * omega = rhs
						double omega[27]; // column-major
						for(int i=0; i<9; i++)
						{
							MATRIX_VECTOR_MULTIPLY3X3(invG, &rhs[3 * i], &omega[3 * i]);
						}

						double dRdF[81]; // each column is skew(omega) * R ; column-major
						for(int i=0; i<9; i++)
						{
							double skew[9];
							SKEW_MATRIX(&omega[3 * i], skew);
							MATRIX_MULTIPLY3X3(skew, R, &dRdF[9 * i]);
						}

						double B[3][3][9];
						// re-arrange dRdF into B, for easier dRdF * dFdx multiplication (to exploit sparsity of dFdx)
						for(int i=0; i<3; i++)
							for(int j=0; j<3; j++)
								for(int k=0; k<3; k++)
									for(int l=0; l<3; l++)
									{
										int row = 3 * i + k;
										int column = 3 * j + l;
										B[i][j][3 * k + l] = dRdF[9 * column + row];
									}

									// four pointers to a 3-vector
									double * minv[4] = { &MInverse[el][0], &MInverse[el][4], &MInverse[el][8], &MInverse[el][12] }; // the four rows of MInverse (last column ignored)

									double dRdx[108]; // derivative of the element rotation matrix with respect to the positions of the tet vertices; column-major
									for(int k=0; k<4; k++)
										for(int i=0; i<3; i++)
											for(int j=0; j<3; j++)
											{
												double temp[3];
												MATRIX_VECTOR_MULTIPLY3X3(B[i][j], minv[k], temp);
												int row = 3 * i;
												int column = 3 * k + j;
												VECTOR_SET3(&dRdx[9 * column + row], temp);
											}

											// add contribution of dRdx to KElement

											// term 1: \hat{dR/dxl} K (R^T x - m)

											// compute K (R^T x - m)
											double tempVec[12]; // R^T x - m
											for(int vtx=0; vtx<4; vtx++)
											{
												double pos[3];
												for(int i=0; i<3; i++)
													pos[i] = P[4 * i + vtx];
												MATRIX_VECTOR_MULTIPLY3X3T(R, pos, &tempVec[3*vtx]);
												// subtract m
												for(int i=0; i<3; i++)
													tempVec[3*vtx+i] -= undeformedPositions[3 * vtxIndex[vtx] + i];
											}
											double a[12]; // a = K * tempVec
											for (int i=0; i<12; i++)
											{
												a[i] = 0.0;
												for (int j=0; j<12; j++)
													a[i] += KElementUndeformed[el][12 * i + j] * tempVec[j];
											}

											// add [\hat{dR/dxl} K R^T x]_l, l=1 to 12
											for(int column=0; column<12; column++)
											{
												double b[12]; // b = \hat{dR/dxl} * a
												for(int j=0; j<4; j++)
												{
													MATRIX_VECTOR_MULTIPLY3X3(&dRdx[9 * column], &a[3*j], &b[3*j]);
												}
												// write b into KElement (add b to i-th column)
												for(int row=0; row<12; row++)
													KElement[12 * row + column] += b[row]; // KElement is row-major
											}

											// term 2: (R K \hat{dRdxl}^T)x

											// re-write positions into a
											for(int vtx=0; vtx<4; vtx++)
											{
												for(int i=0; i<3; i++)
													a[3 * vtx + i] = P[4 * i + vtx];
											}

											// compute [\hat{dRdxl}^T x)]_l, l=1 to 12
											for(int column=0; column<12; column++)
											{
												double b[12]; // b = \hat{dRdxl}^T * a
												for(int j=0; j<4; j++)
												{
													MATRIX_VECTOR_MULTIPLY3X3T(&dRdx[9 * column], &a[3*j], &b[3*j]);
												}

												// add RK * b to column of KElement
												int rowStart = 0;
												double* pRK = RK;
												for (int row=0; row<12; row++)
												{
													double contrib = 0.0;
													// 													for (int j=0; j<12; j++)
													// 														contrib += RK[rowStart + j] * b[j];

													// 进行优化
													double* pb  = b;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;

													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;

													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;

													KElement[rowStart + column] += contrib;
													rowStart += 12;
												}
											}
				}
		}
		else
		{
			// no warp
			memcpy(KElement, KElementUndeformed[el], sizeof(double) * 144);
			// f = K u
			double fElement[12];
			for(int i=0; i<12; i++)
			{
				fElement[i] = 0;
				for(int j=0; j<4; j++)
				{
					fElement[i] += 
						KElement[12 * i + 3 * j + 0] * u[3 * vtxIndex[j] + 0] +
						KElement[12 * i + 3 * j + 1] * u[3 * vtxIndex[j] + 1] +
						KElement[12 * i + 3 * j + 2] * u[3 * vtxIndex[j] + 2];
				}
			}

			// add fElement into the global f
			if (f != NULL)
			{
				for(int j=0; j<4; j++)
				{
					f[3 * vtxIndex[j] + 0] += fElement[3 * j + 0] * matFactor;
					f[3 * vtxIndex[j] + 1] += fElement[3 * j + 1] * matFactor;
					f[3 * vtxIndex[j] + 2] += fElement[3 * j + 2] * matFactor;
				}
			}
		}

		if (stiffnessMatrix != NULL)
		{
			int * rowIndex = rowIndices[el];
			int * columnIndex = columnIndices[el];

			// add KElement to the global stiffness matrix
			for (int i=0; i<4; i++)
				for (int j=0; j<4; j++)
					for(int k=0; k<3; k++)
						for(int l=0; l<3; l++)
						{
							double value = KElement[12 * (3 * i + k) + 3 * j + l] * matFactor;
							stiffnessMatrix->AddEntry(3 * rowIndex[i] + k, 3 * columnIndex[4 * i + j] + l, value);
						}
		}
	}
}

void CorotationalLinearFEMWrapper::ComputeForceAndStiffnessMatrix( double * u, double * f, SparseMatrix * stiffnessMatrix, int warp )
{
	CorotationalLinearFEMWrapper::ComputeForceAndStiffnessMatrixOfSubmesh(u, f, stiffnessMatrix, warp, 0, tetMesh->getNumElements());
}

void CorotationalLinearFEMWrapper::computeReducedHessianMatrix( const double * vertexDisplacements, const EigDense& reduceMat, EigDense& reducedEleHessianMat, SparseMatrix * stiffnessMatrixBuffer )
{
	// 初始化参数
	const double* u = vertexDisplacements;
	int elementLo = 0;
	int elementHi = tetMesh->getNumElements();
	int nReducedDim = reduceMat.rows();
	int nOriginalDim= reduceMat.cols();
	int warp = m_wrap;
	SparseMatrix * stiffnessMatrix = stiffnessMatrixBuffer;
	EigDense tempMat(nOriginalDim, nReducedDim);
	reducedEleHessianMat.resize(nReducedDim*nReducedDim, elementHi - elementLo);


	for (int el=elementLo; el < elementHi; el++)
	{
		double matFactor = 1;

		int vtxIndex[4];
		for (int vtx=0; vtx<4; vtx++)
			vtxIndex[vtx] = tetMesh->getVertexIndex(el, vtx);

		double KElement[144]; // element stiffness matrix, to be computed below; row-major

		if (warp > 0)
		{
			double P[16]; // the current world-coordinate positions (row-major)
			/*
			P = [ v0   v1   v2   v3 ]
			[  1    1    1    1 ]
			*/
			// rows 1,2,3
			for(int i=0; i<3; i++)
				for(int j=0; j<4; j++)
					P[4 * i + j] = undeformedPositions[3 * vtxIndex[j] + i] + u[3 * vtxIndex[j] + i];
			// row 4
			for(int j=0; j<4; j++)
				P[12 + j] = 1;

			// F = P * Inverse(M)
			double F[9]; // upper-left 3x3 block
			for(int i=0; i<3; i++) 
				for(int j=0; j<3; j++) 
				{
					F[3 * i + j] = 0;
					for(int k=0; k<4; k++)
						F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
				}

				double R[9]; // rotation (row-major)
				double S[9]; // symmetric (row-major)
				double tolerance = 1E-6;
				int forceRotation = 1;
				PolarDecomposition::Compute(F, R, S, tolerance, forceRotation);

				// RK = R * K
				// KElement = R * K * R^T
				double RK[144]; // row-major
				WarpMatrix(KElementUndeformed[el], R, RK, KElement);

				// compute exact stiffness matrix
				if (warp == 2)
				{
					// compute G = (tr(S) I - S) R^T
					double G[9]; 
					double tr = S[0] + S[4] + S[8];
					double temp[9];
					for(int i=0; i<9; i++)
						temp[i] = -S[i];
					temp[0] += tr;
					temp[4] += tr;
					temp[8] += tr;
					// G = temp * R^T
					MATRIX_MULTIPLY3X3ABT(temp, R, G);

					double invG[9]; // invG = G^{-1}
					inverse3x3(G, invG);

					double rhs[27]; // 3 x 9 matrix (column-major)
					for(int i=0; i<3; i++)
						for(int j=0; j<3; j++)
						{
							double temp[9];
							for(int k=0; k<9; k++)
								temp[k] = 0.0;
							// copy i-th row of R into column j of temp      
							for(int k=0; k<3; k++)
								temp[3 * k + j] = R[3 * i + k];
							// extract the skew-symmetric part
							SKEW_PART(temp, &rhs[3 * (3 * i + j)]);
						}
						// must undo division by 2 from inside the SKEW_PART macro
						for(int i=0; i<27; i++)
							rhs[i] *= 2.0;

						// solve G * omega = rhs
						double omega[27]; // column-major
						for(int i=0; i<9; i++)
						{
							MATRIX_VECTOR_MULTIPLY3X3(invG, &rhs[3 * i], &omega[3 * i]);
						}

						double dRdF[81]; // each column is skew(omega) * R ; column-major
						for(int i=0; i<9; i++)
						{
							double skew[9];
							SKEW_MATRIX(&omega[3 * i], skew);
							MATRIX_MULTIPLY3X3(skew, R, &dRdF[9 * i]);
						}

						double B[3][3][9];
						// re-arrange dRdF into B, for easier dRdF * dFdx multiplication (to exploit sparsity of dFdx)
						for(int i=0; i<3; i++)
							for(int j=0; j<3; j++)
								for(int k=0; k<3; k++)
									for(int l=0; l<3; l++)
									{
										int row = 3 * i + k;
										int column = 3 * j + l;
										B[i][j][3 * k + l] = dRdF[9 * column + row];
									}

									// four pointers to a 3-vector
									double * minv[4] = { &MInverse[el][0], &MInverse[el][4], &MInverse[el][8], &MInverse[el][12] }; // the four rows of MInverse (last column ignored)

									double dRdx[108]; // derivative of the element rotation matrix with respect to the positions of the tet vertices; column-major
									for(int k=0; k<4; k++)
										for(int i=0; i<3; i++)
											for(int j=0; j<3; j++)
											{
												double temp[3];
												MATRIX_VECTOR_MULTIPLY3X3(B[i][j], minv[k], temp);
												int row = 3 * i;
												int column = 3 * k + j;
												VECTOR_SET3(&dRdx[9 * column + row], temp);
											}

											// add contribution of dRdx to KElement

											// term 1: \hat{dR/dxl} K (R^T x - m)

											// compute K (R^T x - m)
											double tempVec[12]; // R^T x - m
											for(int vtx=0; vtx<4; vtx++)
											{
												double pos[3];
												for(int i=0; i<3; i++)
													pos[i] = P[4 * i + vtx];
												MATRIX_VECTOR_MULTIPLY3X3T(R, pos, &tempVec[3*vtx]);
												// subtract m
												for(int i=0; i<3; i++)
													tempVec[3*vtx+i] -= undeformedPositions[3 * vtxIndex[vtx] + i];
											}
											double a[12]; // a = K * tempVec
											for (int i=0; i<12; i++)
											{
												a[i] = 0.0;
												for (int j=0; j<12; j++)
													a[i] += KElementUndeformed[el][12 * i + j] * tempVec[j];
											}

											// add [\hat{dR/dxl} K R^T x]_l, l=1 to 12
											for(int column=0; column<12; column++)
											{
												double b[12]; // b = \hat{dR/dxl} * a
												for(int j=0; j<4; j++)
												{
													MATRIX_VECTOR_MULTIPLY3X3(&dRdx[9 * column], &a[3*j], &b[3*j]);
												}
												// write b into KElement (add b to i-th column)
												for(int row=0; row<12; row++)
													KElement[12 * row + column] += b[row]; // KElement is row-major
											}

											// term 2: (R K \hat{dRdxl}^T)x

											// re-write positions into a
											for(int vtx=0; vtx<4; vtx++)
											{
												for(int i=0; i<3; i++)
													a[3 * vtx + i] = P[4 * i + vtx];
											}

											// compute [\hat{dRdxl}^T x)]_l, l=1 to 12
											for(int column=0; column<12; column++)
											{
												double b[12]; // b = \hat{dRdxl}^T * a
												for(int j=0; j<4; j++)
												{
													MATRIX_VECTOR_MULTIPLY3X3T(&dRdx[9 * column], &a[3*j], &b[3*j]);
												}

												// add RK * b to column of KElement
												int rowStart = 0;
												double* pRK = RK;
												for (int row=0; row<12; row++)
												{
													double contrib = 0.0;
													// 													for (int j=0; j<12; j++)
													// 														contrib += RK[rowStart + j] * b[j];

													// 进行优化
													double* pb  = b;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;

													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;

													contrib += *pRK * *pb;	pRK++; pb++;
													contrib += *pRK * *pb;	pRK++; pb++;

													KElement[rowStart + column] += contrib;
													rowStart += 12;
												}
											}
				}
		}
		else
		{
			// no warp
			memcpy(KElement, KElementUndeformed[el], sizeof(double) * 144);
		}

		int * rowIndex = rowIndices[el];
		int * columnIndex = columnIndices[el];

		// add KElement to the global stiffness matrix
		// clear stiffness matrix to zero
		stiffnessMatrix->ResetToZero();
		for (int i=0; i<4; i++)
		{
			for (int j=0; j<4; j++)
			{
				for(int k=0; k<3; k++)
					for(int l=0; l<3; l++)
					{
						double value = KElement[12 * (3 * i + k) + 3 * j + l];
						stiffnessMatrix->AddEntry(3 * rowIndex[i] + k, 3 * columnIndex[4 * i + j] + l, value);
					}
			}
		}

		// 计算 R * K * R^T
		tempMat.setZero(nOriginalDim, nReducedDim);
		for (int i=0; i<4; i++)
		{
			for(int k=0; k<3; k++)
			{
				int rowID = 3 * rowIndex[i] + k;
				for (int ithParam = 0; ithParam < nReducedDim; ++ithParam)
				{
					int rowLength = stiffnessMatrix->GetRowLength(rowID);
					double res = 0;
					for (int ithRowEntry = 0; ithRowEntry < rowLength; ++ithRowEntry)
					{
						int colID = stiffnessMatrix->GetColumnIndex(rowID, ithRowEntry);
						res += stiffnessMatrix->GetEntry(rowID, ithRowEntry) * reduceMat(ithParam, colID);
					}

					tempMat(rowID, ithParam) = res;
				}
			}
		}

		EigDense reducedStiffnessMat = reduceMat * tempMat;
		for (int ithParam = 0, ithTar = 0; ithParam < nReducedDim; ++ithParam)
		{
			for (int jthParam = 0; jthParam < nReducedDim; ++jthParam, ++ithTar)
			{
				reducedEleHessianMat(ithTar, el) = reducedStiffnessMat(jthParam, ithParam);
			}
		}
	}
}
