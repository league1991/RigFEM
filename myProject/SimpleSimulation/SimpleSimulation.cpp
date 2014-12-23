// SimpleSimulation.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


int main(int argc, char* argv[])
{
	char* inputFileName = "I:/Programs/VegaFEM-v2.1/models/turtle/turtle-volumetric-homogeneous.veg";
	VolumetricMesh* volumetricMesh = VolumetricMeshLoader::load(inputFileName);
	if (volumetricMesh == NULL)
	{
		printf("load failed!");
	}
	else 
	{
		printf("%d vertices, %d elements", volumetricMesh->getNumVertices(), volumetricMesh->getNumElements());
	}

	TetMesh* tetMesh;
	if (volumetricMesh->getElementType() == VolumetricMesh::TET)
	{
		tetMesh = (TetMesh*) volumetricMesh;
	}
	else 
		printf("not a tet mesh\n");

	CorotationalLinearFEM* deformableModel = new CorotationalLinearFEM(tetMesh);
	ForceModel* forceModel = new CorotationalLinearFEMForceModel(deformableModel);

	int nVtx = tetMesh->getNumVertices();
	int r = 3 * nVtx;
	double timestep = 0.0333;

	SparseMatrix* massMatrix;
	GenerateMassMatrix::computeMassMatrix(tetMesh, &massMatrix, true);
	massMatrix->SaveToMatlabFormat("massMatrix.m");

	printf("%d rows, %d cols\n", massMatrix->GetNumRows(), massMatrix->GetNumColumns());
	int positiveDefiniteSolver = 0;

	int numConstrainedDOFs = 9;
	int constrainedDOFs[9]= {12,13,14,30,31,32,42,43,44};

	double dampingMassCoef = 0.0;
	double dampingStiffnessCoef = 0.01;


	ImplicitBackwardEulerSparse* integrator = new ImplicitBackwardEulerSparse(r, timestep, massMatrix, forceModel, positiveDefiniteSolver, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef);

	//CentralDifferencesSparse* integrator = new CentralDifferencesSparse(r, timestep, massMatrix, forceModel, numConstrainedDOFs, constrainedDOFs, dampingMassCoef, dampingStiffnessCoef);

	double * f = new double[r];
	int numTimesteps = 10;
	double*u = new double[r];
	for (int i = 0; i < numTimesteps; ++i)
	{
		integrator->SetExternalForcesToZero();
		if (i==0)
		{
			for (int j = 0; j < r; j++)
				f[j] = 0;
			f[37] = -500;
			integrator->SetExternalForces(f);
		}
		integrator->GetqState(u);
		printf("v = [", i);
		for (int ithVtx = 0; ithVtx < nVtx; ++ithVtx)
			printf("%lf, %lf, %lf\n", u[ithVtx*3],u[ithVtx*3+1],u[ithVtx*3+2]);
		printf("];");

		integrator->DoTimestep();
	}

	return 0;
}

 