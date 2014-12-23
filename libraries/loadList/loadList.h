/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 2.1                               *
 *                                                                       *
 * "loadList" library , Copyright (C) 2007 CMU, 2009 MIT                 *
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

#ifndef _LOADLIST_H_
#define _LOADLIST_H_

/*
  A class to load an integer list from a text file into memory. 
*/

#include <stdio.h>

class LoadList
{
public:

  // returns 0 on success and non-zero otherwise
  static int load(const char * filename, int * numListEntries, int ** listEntries, int offset=0);
  static int save(const char * filename, int numListEntries, int * listEntries, int offset=0);

  static int loadBinary(const char * filename, int * numListEntries, int ** listEntries, int offset=0);
  static int loadBinary(FILE * fin, int * numListEntries, int ** listEntries, int offset=0);
  static int saveBinary(const char * filename, int numListEntries, int * listEntries, int offset=0);
  static int saveBinary(FILE * fout, int numListEntries, int * listEntries, int offset=0);

  // loads/saves multiple lists to one binary file
  static int loadBinaryMulti(const char * filename, int * numLists, int ** numListEntries, int *** listEntries, int offset=0);
  static int saveBinaryMulti(const char * filename, int numLists, int * numListEntries, int ** listEntries, int offset=0);

  static void sort(int numListEntries, int * listEntries);
  static void print(int numListEntries, int * listEntries);
  static void stripBlanks(char * s);

protected:
};

#endif

