/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "Large Modal Deformation Factory",                                    *
 * a pre-processing utility for model reduction of                       *
 * deformable objects undergoing large deformations.                     *
 *                                                                       *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2014 USC                           *
 *                                                                       *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef _STATES_H_
#define _STATES_H_

// the current state of the user interface
typedef struct
{
  // viewing 
  bool showAxes;
  enum viewModeType { VIEW_NONE, VIEW_RENDERING_MESH, 
    VIEW_SIMULATION_MESH, VIEW_LINEAR_MODES, 
    VIEW_NONLINEAR_MODES, VIEW_RUNTIME_SIMULATION } viewMode;

  // rendering mesh
  wxString renderingMeshFilename;

  // sq mesh generation
  bool fillInteriorChambers;
  int cubicMeshResolution;

  // vertex selection
  bool vertexSelectionActivated;

  // linear modes
  int numComputedLinearModes;
  bool firstModalComputation;

  // frequencies
  double lastFrequencyScaleFactor;
  int eraseRangeLo;
  int eraseRangeHi;

  // nonlinear modes
  int numComputedNonLinearModes;

  // num threads for the cubic polynomial computation
  int numComputationThreads;

  wxString projectName;
  bool scaleRealTimeTo1HzCheckBox;

  wxString manySimulationsscriptFilename;
  wxString currentWorkingDirectory;

  bool enableBatchOutput;
  bool showMaterialGroups;
  bool renderMesh;

} UIState;

// the current state of the precomputation pipeline
typedef struct
{
  // rendering mesh
  bool renderingMeshAvailable;
  ObjMesh * renderingMesh;

  // fixed vertices
  bool fixedVerticesAvailable;
  set<int> fixedVertices;
  
  // simulation mesh
  bool simulationMeshAvailable;
  VolumetricMesh * simulationMesh;
  int * FFDinterpolationVertices;
  double * FFDinterpolationWeights;
  ObjMesh * simulationSurfaceMesh;
  struct {
    double centerX, centerY, centerZ;
    double radius;
  } simulationMeshGeometricParameters;

  // linear modes
  bool linearModesAvailable;
  ModalMatrix * linearModalMatrix;
  int rLin;

  // num rigid modes
  int numRigidModes;

  // frequencies
  bool frequenciesAvailable;
  double * frequencies;

  // nonlinear modes
  bool nonLinearModesAvailable;
  ModalMatrix * nonLinearModalMatrix;
  int rNonLin;

  // modal derivatives
  bool modalDerivativesAvailable;
  int numDeriv;
  ModalMatrix * modalDerivativesMatrix;

  // sketch
  bool sketchDataAvailable;
  ModalMatrix * sketchDataMatrix;

  // cubic polynomials
  bool cubicPolynomialsAvailable;
  StVKReducedInternalForces * cubicPolynomials;

  // interpolation
  bool interpolationDataAvailable;
  int * interpolationData_vertices;
  double * interpolationData_weights;

  bool runtimeSimReady;

} PrecomputationState;

#endif

