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

#include "insertRows.h"
#include "stdio.h"
#include "stdlib.h"

void InsertRows(int mFull, double * xConstrained, double * x, int numFixedRows, int * fixedRows, int oneIndexed)
{
  int destRow = 0; // current row in the dest matrix
  int sourceRow = 0; // in source

  for(int i=0; i<numFixedRows; i++)
  {
    int index = fixedRows[i] + 1 - oneIndexed;
    if ((index > mFull) || (index < 1))
    {
      printf("Error: invalid index %d specified.\n",index);
      exit(1);
    }
    index--;

    while (destRow < index)
    {
      // while row index smaller than index, keep on copying from source
      x[destRow] = xConstrained[sourceRow];
      destRow++;
      sourceRow++;
    }

    // insert zero row
    x[destRow] = 0.0;
    destRow++;
  }
  
  while (destRow < mFull)
  { 
    x[destRow] =  xConstrained[sourceRow];
    destRow++;
    sourceRow++;
  }
}

void RemoveRows(int mFull, double * xConstrained, double * x, int numFixedRows, int * fixedRows, int oneIndexed)
{
  int numrows = 0;
  int row = 0;

  for(int i=0; i<numFixedRows; i++)
  {
    int index = fixedRows[i] + 1 - oneIndexed;
    if ((index > mFull) || (index < 1))
    {
      printf("Error: invalid index %d specified.\n",index);
      exit(1);
    }
    index--;

    while (row<index)
    {
      xConstrained[numrows] = x[row];
      numrows++;
      row++;
    }

    row++; // skip the deselected row

    if (numrows > mFull)
    {
      printf("Error: too many rows specified.\n");
      exit(1);
    }
  }

  while (row < mFull)
  {
    xConstrained[numrows] = x[row];

    numrows++;
    row++;

    if (numrows > mFull)
    {
      printf("Error: too many rows specified.\n");
      exit(1);
    }
  }
}

void FullDOFsToConstrainedDOFs(int mFull, int numDOFs, int * DOFsConstrained, int * DOFs, int numFixedRows, int * fixedRows, int oneIndexed)
{
  if (numDOFs == 0)
    return;

  int dof = 0;

  for(int i=0; i<numFixedRows; i++)
  {
    // each iteration processes one bucket of fixed vertices

    // correct for (optional) 1-indexing
    int index = fixedRows[i] + 1 - oneIndexed;
    if ((index > mFull) || (index < 1))
    {
      printf("Error: invalid index %d specified.\n",index);
      exit(1);
    }
    index--;

    while (DOFs[dof] < index)
    {
      DOFsConstrained[dof] = DOFs[dof] - i; 
      dof++;
      if (dof >= numDOFs)
        return;
    }

    // assign -1 to the fixed DOF
    if (DOFs[dof] == index)
    {
      DOFsConstrained[dof] = -1;
      dof++;
    }
  }

  while (dof < numDOFs)
  {
    DOFsConstrained[dof] = DOFs[dof] - numFixedRows; 
    dof++;
  }
}

