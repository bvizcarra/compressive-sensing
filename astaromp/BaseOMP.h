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
#include "GlobalUtil.h"
#include "VectorMath.h"
#include "AStarDefinitions.h"

using namespace std;

/// Struct that defines the side info (necessary for the BaseOMP class) assigned individually to each path
/// This struct holds all info for a path that is necessary for BaseOMP. The info in this struct is not known by BaseAStar, 
/// BaseAStar keeps only void* type pointers to the side info of each path.
/// SideInfo keeps the QR decomposition of the selected parts of the dictionary, indices of selected atoms,
/// and the residue.
struct SideInfo
{
	float** mQ; ///< pointer to the Q (orthogonal) matrix from QR decomposition
	float** mR; ///< pointer to the R (upper triangular) matrix from QR decomposition
	float* mZ;  ///< pointer to the Z vector (\f$z = Q^Ty\f$)
	float* mRes; ///< pointer to the residue of a path 
	vector<unsigned int> mIndList; ///< vector holding ID's of the nodes in a path
};

/// BaseOMP is the base class for the OMP part of the A*OMP. In this A*OMP implementation, this class is the problem class,
/// that defines the problem A* search tries to solve.
/// BaseOMP involves operations based on OMP, 
/// such as finding initial paths, finding the best candidates, computing priorities and costs,
/// performing post operations to extract the solution from the side info and etc.
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
class BaseOMP
{
public:
	/// Constructor
	BaseOMP(int pK,int pM,int pN, float pEps, int pNodesPerInitPath);

	/// Destructor
	~BaseOMP(void);
	
	/// Function to allocate an instance of SideInfo struct 
	SideInfo* allocateSideInfo();

	/// Function to set mDict
	void setDict(float** pDict);

	/// Function to set y
	void sety(float* py);

	/// Function to compute priorities of dictionary members
	void computePriorities(priority* pPriority);

	/// Function to compute the pre-cost of a path from SideInfo
	float computeCost(SideInfo* pSideInfo, int pNewElementID);

	/// Function to find the best candidates for expansion of a path
	void findBestCandidates(  int pNoCand, SideInfo* pSideInfo, elementID* pCandList  );

	/// Function to perform necessary operations after A* search is terminated
	void performPostOperations(  SideInfo* pSideInfo  );

	/// Function to create a duplicate of a SideInfo struct
	SideInfo* duplicateSideInfo(SideInfo* pSrc);

	/// Function to copy the SideInfo struct pSrc into pDst
	void copySideInfo(SideInfo* pSrc, SideInfo* pDst);

	/// Function to return solution of A*OMP
	float* getSolution();

	/// Function to reset a SideInfo struct
	void resetSideInfo(SideInfo* pSideInfo);

	/// Function to return \f$l_2\f$ norm of y
	float getNorm_y();
	
	/// Function to check if the termination criterion of A*OMP is satisfied
	bool isSearchComplete(int pPL, float pErr);

	/// Function to initialize paths for A*OMP
	int findInitialPaths( int pNoInitialPaths, vector<unsigned int*> *pNodeList, int &pNodesPerPath);

private:
	/// Function to solve for the coefficients from a QR decomposition
	void solveCoefs(vector<unsigned int> *pIndList, float* pZ, float** pR);

	/// Function to add a new element to the QR decomposition of a path
	void addElementToRepresentation( float* pNewElement, int i, float** pQ, float* pR, float * pZ, float *pRes);

	/// Function to find a sorted list of vectors which lie closest to a vector among an array of vectors 
	map<float,int,greater<float> > * findClosestVectorsIndList( float** pVectorArray, float* pVectorNorm, float* pVector, int pNoVectors, int pReturnSize, int pSize );

	/// Function to find a sorted list of vectors which lie closest to a vector among an array of vectors 
	void findClosestVectorsIndList( float** pVectorArray, float* pVectorNorm, float* pVector,unsigned int* pReturnList ,int pNoVectors, int pReturnSize, int pSize );
	
	/// Function to find the initial paths for A*OMP
	int findInitialPaths1( int pNoInitialPaths, vector<unsigned int*> *pNodeList );

	/// Function to find the initial paths of length 2 for image reconstruction
	int findInitialPaths2( int pNoInitialPaths, vector<unsigned int*> *pNodeList );


	int mK;				///< desired sparsity (max. number of nonzero elements)
	int mN;				///< dimension of the desired sparse vectors
	int mM;				///< number of observations
	float **mDict;		///< pointer to the matrix containing the dictionary (holographic basis)
	float* mDictNorm;	///< pointer to the vector holding norms of dictionary elements
	float* my;			///< pointer to the observed vector
	float mNorm_y;		///< norm of the observation vector
	float* mSolution;	///< pointer to the vector holding the solution
	float mEps;			///< error toleration for terminating the search
	int mNodesPerInitPath; ///< number of nodes in each initial path (1 or 2)
};

