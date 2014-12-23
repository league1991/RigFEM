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
					fElement[i] += 
						KElement[12 * i + 3 * j + 0] * u[3 * vtxIndex[j] + 0] +
						KElement[12 * i + 3 * j + 1] * u[3 * vtxIndex[j] + 1] +
						KElement[12 * i + 3 * j + 2] * u[3 * vtxIndex[j] + 2];
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
				  f[ithDof] += KElement[ithDof*12 + jthDof] * rotOffset[jthDof];
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

		  double lambda = eNuMaterial->getLambda();
		  double mu     = eNuMaterial->getMu();
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
			  printf("%dth element: energy1 = %lf energy2 = %lf error = %lf\n", el, eleEnergy, eleEnergy2, eleEnergy2 - eleEnergy);
			  for (int i = 0; i < 9; ++i)
			  {
				  printf("%lf ", S[i]);
			  }
			  printf("\n");
		  }
#endif
		}
	}
	return energy;
}
