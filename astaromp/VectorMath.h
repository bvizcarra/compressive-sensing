/*
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

/// VectorMath is a collection of functions that perform vector and matrix operations.
/// Note that vector is used as a mathematical term instead of the C++ data structure.
/// A vector is represented by a 1D array. A matrix represented by an array
/// of pointers to its columns, which are also stored as arrays.
/// The suffix _I at the end of some function names denote functions that perform in-place
/// operations, i.e. overwrite the result on the source vector.
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
#pragma once

#include <vector>
#include <map>
#include <math.h>

#include <string.h>

using namespace std;

/// Function to allocate a float matrix
float** allocateFloatMatrix(int pCols, int pRows);

/// Function to allocate an integer matrix
int** allocateIntMatrix(int pCols, int pRows);

/// Function to multiply a matrix with a vector.
void multMatrixWithVector(float* pSrcMatrix, float* pSrcVec, int pRows, int pCols, float* pDestVec);

/// Function to normalize a vector
float normalizeVector(float* pSrc, float* pDst, int pSize);

/// Function to normalize a vector (in place)
float normalizeVector_I(float* pSrcDst, int pSize);

/// Function to delete a float matrix
void deleteFloatMatrix( float** pMatrix, int pCols );

/// Function to compute \f$l_2\f$ norm of a vector
float l2Norm(float* pVector, int pSize);

/// Function to subtract one vector from another
void subtractVectorfromVector(float* pSrc1, float* pSrc2, float* pDst,int pSize);

/// Function to subtract one vector from another (in-place operation)
void subtractVectorfromVector_I(float* pSrcDst, float* pSrc2, int pSize);

/// Function to subtract product of a scalar and a vector from a vector (in-place operation)
void subtractProductScalarfromVector_I(float* pSrcDst, float* pSrc2, float pScalar, int pSize);

/// Function to subtract product of a scalar and a vector from a vector
void subtractProductScalarfromVector(float* pSrc1, float* pSrc2, float *pDst, float pScalar, int pSize);

/// Function to subtract element-wise product of two vectors from a vector (in-place operation)
void subtractProductfromVector_I(float* pSrcDst, float* pSrc1, float* pSrc2, int pSize);

/// Function to multiply a vector with a scalar
void multVectorwithScalar(float *pSrc, float* pDst, float pScalar, int pSize);

/// Function to multiply a vector with a scalar (in-place operation)
void multVectorwithScalar_I( float *pSrcDst, float pScalar, int pSize );

/// Function to sum two vectors
void sumVectors( float *pSrc1, float *pSrc2, float* pDst, int pSize );

/// Function to sum two vectors (in-place operation)
void sumVectors_I( float *pSrcDst, float *pSrc2, int pSize );

/// Function to divide a vector by a scalar (in-place operation)
void divideVectorByScalar_I(float *pSrcDst, float pScalar, int pSize);

/// Function to divide a vector by a scalar
void divideVectorByScalar( float *pSrc, float* pDst, float pScalar, int pSize );

/// Function to copy a vector into another
void copyVector(float* pSrc, float* pDst, int pSize);

/// Function to delete an integer matrix
void deleteIntMatrix( int** pMatrix, int pCols );

/// Function to compute inner-product of two vectors
float computeInnerProd(float* pSrc1, float* pSrc2, int pSize);

/// Function to solve a linear system with upper triangular (UT) system matrix
void linsolve_UT(float **pA, float *pb, float *px, int pNoCols);

/// Function to compute mean of a vector
float computeMean(float *pSrc, int pSize);

/// Function to sum elements of a vector
float sumVector(float *pSrc, int pSize);

/// Function to copy a matrix into another
void copyMatrix( float** pSrc, float** pDst, int pCols, int pRows );
