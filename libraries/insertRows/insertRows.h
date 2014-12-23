/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "insertRows" library , Copyright (C) 2007 CMU, 2009 MIT               *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 * Research: Jernej Barbic, Doug L. James, Jovan Popovic                 *
 * Funding: NSF, Link Foundation, Singapore-MIT GAMBIT Game Lab          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef _INSERT_ROWS_H_
#define _INSERT_ROWS_H_

/*
  Insert or remove given components from a linear array.
  The dynamic solver uses these routines to fix the specified vertices and remove
  rows from mass and stiffness matrices, as necessary.
*/

// Note: these routines use the terminology "Rows" because they were mainly used to condense mass and stiffness matrices. They, however, operate on 1D arrays, not directly on 2D matrices

// inserts zero entries into an array, at the specified locations
// the locations must be given with respect to the full array 
// input: xConstrained
// output: x
void InsertRows(int mFull, double * xConstrained, double * x, int numFixedRows, int * fixedRows, int oneIndexed=0); 

// removes entries at the specified locations from an array
// the locations must be given with respect to the full array 
// input: x
// output: xConstrained
void RemoveRows(int mFull, double * xConstrained, double * x, int numFixedRows, int * fixedRows, int oneIndexed=0); 

// translates the array indices from original indices to indices after removal of the specified entries
// input: DOFs (must be sorted) (0-indexed)
// output: DOFsConstrained (0-indexed)
// oneIndexed applies only to fixedRows array, NOT to DOFsConstrained or DOFs
void FullDOFsToConstrainedDOFs(int mFull, int numDOFs, int * DOFsConstrained, int * DOFs, int numFixedRows, int * fixedRows, int oneIndexed=0); 

#endif

