// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // 从 Windows 头中排除极少使用的资料
// Windows 头文件:
#include <windows.h>
#include <stdlib.h>


#include <maya/MGlobal.h>
#include <maya/MQtUtil.h>

#include <maya/MPxNode.h>
#include <maya/MPxCommand.h>
#include <maya/MPxData.h>
#include <maya/MPxLocatorNode.h>

#include <maya/MIOStream.h>
#include <maya/MDGModifier.h>
#include <maya/MDagModifier.h>
#include <maya/MTypeId.h> 
#include <maya/MPlug.h>
#include <maya/MPlugArray.h> 
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MPointArray.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h> 
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MMatrix.h>
#include <maya/MArgList.h>
#include <maya/MAnimControl.h>
#include <maya/MVectorArray.h>

#include <maya/MFnNurbsSurfaceData.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnPluginData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnDoubleArrayData.h>


#include <maya/MItMeshVertex.h>
#include <maya/MItMeshEdge.h>
#include <maya/MItSelectionList.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MItMeshPolygon.h>


#include <QDataStream>
#include <QPolygon>
#include <QHostAddress>
#include <QSharedPointer>
#include <QTcpServer>
#include <QTcpSocket>
#include <QThread>
#include <QPainterPath>
#include <QNetworkInterface>
#include <QFont>

// 系统
#include <stdio.h>
#include <tchar.h>
#include <stdlib.h>
#include <math.h>
#include <GL/GLUT.H>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <float.h>

// Vega
#include "vec3d.h"
#include "tetMesh.h"
#include "polarDecomposition.h"
#include "volumetricMeshLoader.h"
#include "volumetricMeshENuMaterial.h"
#include "corotationalLinearFEM.h"
#include "corotationalLinearFEMForceModel.h"
#include "generateMassMatrix.h"
#include "implicitBackwardEulerSparse.h"
#include "centraldifferencessparse.h"
#include "neoHookeanIsotropicMaterial.h"
#include "isotropicHyperelasticFEMForceModel.h"
#include "matrixMultiplyMacros.h"

// Eigen
#include "Eigen/sparse"
#include "Eigen/LU"
#include "Eigen/QR"
#include "Eigen/SuperLUSupport"

typedef Eigen::SparseMatrix<double> EigSparse;
typedef Eigen::MatrixXd				EigDense;
typedef Eigen::VectorXd				EigVec;

// 体网格生成
#include "tetgen.h"

#define PRINT_F(format,...)			{char buf[300];sprintf_s(buf, 299, format, ##__VA_ARGS__);MGlobal::displayInfo(buf);}

// 本地文件
#include "Utilities.h"
#include "CorotationalLinearFEMWrapper.h"
#include "rig.h"
#include "StatusRecorder.h"
#include "NewtonSolver.h"
#include "FEMSystem.h"
#include "RiggedSkinMesh.h"
#include "RigSimulator.h"
#include "nanoflann.hpp"

#include "globals.h"
#include "Simple1DNoiseNode.h"
#include "GeneralRig.h"
#include "FEMSimulationNode.h"
#include "ExecSimulationCmd.h"
#include "MeshCorrespondence.h"
#include "MeshControlCmd.h"
#include "DataSwitch.h"