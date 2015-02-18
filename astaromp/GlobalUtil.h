/*
Copyright 2011 Nazim Burak Karahanoglu

This source code is provided as a part of AStarOMP project. 

Using, altering and redistributing this software is permitted to anyone for academical purposes, 
with to the following restrictions:

1 - Original code shall not be misrepresented.

2 - Modifications made to the code should be clearly indicated.

3 - You must not claim that this is your own code.

4 - This note may not be removed or modified. 

In case you use this code in a product, an acknowledgment in documentation would be appreciated.

The author cannot be held responsible for any damages that arise from using this software.

Nazim Burak Karahanoglu 
karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
*/


/// GlobalUtil is a collection of general purpose functions that are globally used
/// by other classes in this implementation. It includes some file input/output, command line output
/// and memory deallocation functions.
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com

#pragma once
#include <iostream>
#include <stdlib.h>
#include <map>
#include <fstream>
#include <vector>
#include <algorithm>
#include <deque>
#include <cstring>
#include <cmath>

using namespace std;

//************************************************************************//
// Function declarations

/// Function to read a float array from a binary file. 
int freadbin(const char * pFilename, float* pFloatReadArray, int pN);
/// Function to read a float matrix from a binary file. 
int freadbin(const char * pFilename, float** pFloatReadMatrix, int pRows, int pCols);
/// Function to read a float matrix from a text file. 
int freadtxt(char* pFilename, int pNoRows, int pNoCols, float** pDstMatrix);
/// Function to read a float matrix from a text ifstream. 
int freadtxt(ifstream *pFile, float **pDstMatrix ,int pNoCols, int pNoRows);
/// Function to display contents of a float matrix on command line. 
void displayMatrix(float** pMatrix, int pNoRows, int pNoCols);
/// Function to display contents of a integer matrix on command line. 
void displayMatrix(int** pMatrix, int pNoRows, int pNoCols);
/// Function to delete a vector of pointers to float arrays. 
void deletevector(vector<float*> &pVector);
/// Function to delete a vector of pointers to int arrays.
void deletevector(vector<int*> &pVector);
/// Function to delete a vector of pointers to unsigned int arrays.
void deletevector(vector<unsigned int*> &pVector);
/// Function to delete a vector of pointers of type void*.
void deletevector( vector<void*> &pVector );
/// Function to write a float array to a text file. 
int fwritetxt(const char * pFilename, float* pWriteFloatArray,int pSize);
/// Function to write a float array to an ofstream in text format.
int fwritetxt(ofstream* pFile, float* pWriteFloatArray,int pSize);
/// Function to write a float array to a binary file.
int fwritebin(const char * pFilename, float* pWriteFloatArray,int pSize);
/// Function to write a float array to an ofstream in binary format.
int fwritebin(ofstream* pFile, float* pWriteFloatArray,int pSize);