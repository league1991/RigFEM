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

#include "loadList.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

// removes all whitespace characters from string s
void LoadList::stripBlanks(char * s)
{
  char * w = s;
  while (*w != '\0')
  {
    while (*w == ' ') // erase blank
    {
      char * u = w;
      while (*u != '\0') // shift everything left one char
      {
        *u = *(u+1);
        u++;
      }
    }
    w++;
  }
}

int compareLoadList(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int LoadList::load(const char * filename, int * numListEntries, int ** listEntries, int offset)
{
   // comma-separated text file of fixed vertices
  FILE * fin;
  fin = fopen(filename,"r");
  if (!fin)
  {
    printf("Error: could not open file %s.\n",filename);
    return 1;
  }

  *numListEntries = 0;

  char s[4096];
  while (fgets(s,4096,fin) != NULL)
  { 
    stripBlanks(s);

    char * pch;
    pch = strtok (s,",");
    while ((pch != NULL) && (isdigit(*pch)))
    {
      (*numListEntries)++;
      pch = strtok (NULL, ",");
    } 
  }

  *listEntries = (int*) malloc (sizeof(int) * (*numListEntries));

  rewind(fin);

  (*numListEntries) = 0;

  while (fgets(s,4096,fin) != NULL)
  {
    stripBlanks(s);
    char * pch;
    pch = strtok (s,",");
    while ((pch != NULL) && (isdigit(*pch)))
    {
      (*listEntries)[*numListEntries] = atoi(pch) - offset;
      (*numListEntries)++;
      pch = strtok (NULL, ",");
    }
  }

  // sort the list entries
  qsort ((*listEntries), *numListEntries, sizeof(int), compareLoadList);

  fclose(fin);

  return 0;
}

int loadListComparator(const void * a, const void * b)
{
  if (*(int*)a < *(int*)b)
    return -1;

  if (*(int*)a == *(int*)b)
    return 0;

  return 1;
}

void LoadList::sort(int numListEntries, int * listEntries)
{
  qsort(listEntries,numListEntries,sizeof(int),loadListComparator);
}

int LoadList::save(const char * filename, int numListEntries, int * listEntries, int offset)
{
  // comma-separated text file of fixed vertices
  FILE * fout;
  fout = fopen(filename, "w");
  if (!fout)
  {
    printf("Error: could not open list file %s.\n",filename);
    return 1;
  }

  for(int nv=0; nv < numListEntries; nv++)
  {
    fprintf(fout, "%d,", listEntries[nv] + offset);
    if (nv % 8 == 7)
      fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  fclose(fout);

  return 0;
}

void LoadList::print(int numListEntries, int * listEntries)
{
  for(int nv=0; nv < numListEntries; nv++)
  {
    printf("%d,", listEntries[nv]);
    if (nv % 8 == 7)
      printf("\n");
  }
  printf("\n"); fflush(NULL);
}

int LoadList::loadBinary(const char * filename, int * numListEntries, int ** listEntries, int offset)
{
  FILE * fin;
  fin = fopen(filename, "rb");
  if (!fin)
  {
    printf("Error: could not open list file %s.\n",filename);
    return 1;
  }

  int code = loadBinary(fin, numListEntries, listEntries, offset);
  fclose(fin);
  return code;
}

int LoadList::loadBinary(FILE * fin, int * numListEntries, int ** listEntries, int offset)
{
  int item = fread(numListEntries, sizeof(int), 1, fin);
  if (item != 1)
  {
    printf("Error: could not read the number of list entries.\n");
    return 1;
  }
  *listEntries = (int*) malloc (sizeof(int) * (*numListEntries));

  for(int nv=0; nv < *numListEntries; nv++)
  {
    item = fread(&(*listEntries)[nv], sizeof(int), 1, fin);
    if (item != 1)
    {
      printf("Error: could not read the list entry %d.\n", nv);
      return 1;
    }
  }

  for(int nv=0; nv < *numListEntries; nv++)
    (*listEntries)[nv] -= offset;

  return 0;
}

int LoadList::saveBinary(const char * filename, int numListEntries, int * listEntries, int offset)
{
  FILE * fout;
  fout = fopen(filename, "wb");
  if (!fout)
  {
    printf("Error: could not open list file %s.\n",filename);
    return 1;
  }

  int code = saveBinary(fout, numListEntries, listEntries, offset);
  fclose(fout);
  return code;
}

int LoadList::saveBinary(FILE * fout, int numListEntries, int * listEntries, int offset)
{
  fwrite(&numListEntries, sizeof(int), 1, fout);

  for(int nv=0; nv < numListEntries; nv++)
  {
    int value = listEntries[nv] + offset;
    fwrite(&value, sizeof(int), 1, fout);
  }

  return 0;
}

// loads/saves multiple lists to one binary file
int LoadList::loadBinaryMulti(const char * filename, int * numLists, int ** numListEntries, int *** listEntries, int offset)
{
  FILE * fin;
  fin = fopen(filename, "rb");
  if (!fin)
  {
    printf("Error: could not open list file %s.\n",filename);
    return 1;
  }

  int item = fread(numLists, sizeof(int), 1, fin);
  if (item != 1)
  {
    printf("Error: could not read the number of list.\n");
    return 1;
  }

  *numListEntries = (int*) malloc (sizeof(int) * *numLists);
  *listEntries = (int**) malloc (sizeof(int*) * *numLists);

  int code = 0;
  for(int i=0; i<*numLists; i++)
    code = code || loadBinary(fin, &((*numListEntries)[i]), &((*listEntries)[i]), offset);

  fclose(fin);

  return code;
}

int LoadList::saveBinaryMulti(const char * filename, int numLists, int * numListEntries, int ** listEntries, int offset)
{
  FILE * fout;
  fout = fopen(filename, "wb");
  if (!fout)
  {
    printf("Error: could not open list file %s.\n",filename);
    return 1;
  }

  fwrite(&numLists, sizeof(int), 1, fout);

  int code = 0;
  for(int i=0; i<numLists; i++)
    code = code || saveBinary(fout, numListEntries[i], listEntries[i], offset);

  fclose(fout);

  return code;
}

