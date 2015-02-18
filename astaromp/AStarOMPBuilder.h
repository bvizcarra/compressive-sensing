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


#pragma once

#include "BaseAStar.h"
#include "BaseOMP.h"
#include <time.h>
#include <time.h>
#include "GlobalUtil.h"
#include <set>
#include <math.h>
#include "ConfigFile.h"
#include <cstring>

#define myIntend "   "

/// This is the builder class for A*OMP. 
/// It initializes, runs and evaluates A*OMP.
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
class AStarOMPBuilder
{
public:
	/// Default constructor
	AStarOMPBuilder();

	/// Function for initialization with config file
	int init( std::string const pConfigFile );

	/// Destructor
	~AStarOMPBuilder();

	/// Function to run AStarOMP algorithm
	int run();

	/// Function to print evaluation results
	void printEvaluation();

	/// Function to get mTargetVectorsProvided
	bool getTargetVectorsProvided();

private:
	/// Function to read pCols vectors of size pRows from the file pFileName)
	int InitMatrixFromFile(float** pMatrix,  string &pFileName, int pRows, int pCols );

	/// Function to evaluate reconstruction of a single vector
	void evaluateSingleVector(int VectorInd);

	/// Function to check if a vector is K-sparse
	bool isKSparse(float *pVector, int pSize, int pK);

	int countNonzeroElements(float *pVector, int pSize);

	int mK;		///< desired sparsity (max. number of nonzero elements)
	int mN;		///< vector dimension
	int mM;			///< number of observations
	int mInitPL;	///< initial path length
	float mEps;		///< error toleration for terminating the search
	int mNoVectors;  ///< number of test vectors
	float mAlpha;	///< alpha for adaptive-additive cost model
	float mBeta;	///< beta for multiplicative and adaptive-multiplicative cost models
	AuxiliaryFunctionMode mAuxiliaryFunctionMode;	///< parameter for cost function choice
	int mI;		///< I: number of initial A*OMP paths
	int mB;		///< B: number of expanded A*OMP branches per iteration 
	int mP;		///< P: number of maximum search paths in the A* tree

	float **mX;		///< pointer to the matrix of target sparse vectors
	float **mY;		///< pointer to the matrix of observed vectors
	float **mDict;	///< pointer to the matrix containing the dictionary (holographic basis)

	bool* mExRec;	///< pointer to the vector for evaluation of exact reconstruction of test vectors
	int mNoExRecVec;	///< number of exactly reconstructed vectors
	float* mErr;		///< pointer to the reconstruction error of an individual vector
	double mTime;		///< processing time for A*OMP
	float* mNMSE;		///< pointer to the vector of normalized mean squared errors of test vectors

	BaseOMP *mBaseOMP;		///< pointer to the instance of BaseOMP class
	BaseAStar *mBaseAStar;	///< pointer to the instance of BaseAStar class 

	ofstream mRecVectOfstream;		///< ofstream for reconstructed vectors
	bool mBinOutput;				///< parameter for selecting binary or text output

	ofstream mResultOfstream;		///< ofstream for results file

	int mNoIterations;		///< total number of iterations during the search
	int mNoEqBranch;		///< total number of branches that are found to be equivalent to some other path during the search 
	int mNoBranchAdded;		///< total number of branches that are added to the tree during the search 
	int mNoBranchIgnored;	///< total number of branches that are ignored (via Stack size pruning) during the search
	int mNoBranchReplaced;	///< total number of branches that are replaced by their first extensions during the search

	bool mTargetVectorsProvided; ///< states if target vectors are provided
	string mRecVectorsFileName;	///< filename for writing reconstructed vectors
	bool mCompExactRec;			///< decision for computing exact reconstruction (vectors should be K-sparse).
	bool mMultiDict;			///< flag for reading a seperate dictionary for each test vector (multi-dictionary mode)
	ifstream mDictIfstream;		///< ifstream for reading dictionaries in multi-dictionary mode
	bool mBinDictIfstream;		///< parameter for selecting binary multi-dictionary ifstream (true:binary, false: txt)
};


