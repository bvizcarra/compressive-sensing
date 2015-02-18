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

#include "Trie.h"
#include "AlgorithmInterface.h"
#include "VectorMath.h"

/// This class is the base A* implementation and the main interface of the algorithm.
/// It provides functions to initialize, reset and run A*. It holds pointers to necessary 
/// data structures such as the search stack, search tree and interface to the search problem (AlgorithmInterface).
/// Search tree includes all nodes that are opened during the search. It is mainly for keeping all explored paths and
/// finding equivalencies in the tree. Search stack keeps pointers to only "active" paths in the tree, ordered wrt. ascending cost. 
/// In this structure, pruned paths are not removed from the tree, but from the search stack. 
///
/// This class also handles assignment of SideInfo's to paths in the stack, where contents of these are provided by the problem class.
/// SideInfo instances are allocated when necessary by the problem class via mAlgInterface. When a path is removed from the stack,
/// its SideInfo is not deleted, but stored in the vector mFreeSideInfoList for later use. A nes SideInfo is allocated only when 
/// mFreeSideInfoList contains no free SideInfo.
///
/// Note that this class calls no functions from the problem class directly, but runs these via
/// the AlgorithmInterface class . This provides flexibilty to change the search problem 
/// only by modifying the function calls in the AlgorithmInterface, without the necessity of modifying this BaseAStar implementation.
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
class BaseAStar
{
public:
	/// Constructor
	BaseAStar(int pB, int pP, int pI, int pK, int pN, int pM,float pAlpha, float pBeta, AuxiliaryFunctionMode pAuxiliaryFunctionMode);
	/// Destructor
	~BaseAStar(void);

	/// Function to reset and initialize BaseAStar for a new search
	int initialize();

	/// Function to run a new search
	int run();
	
	/// Function to set mB
	void setB(int pB);

	/// Function to set mP
	void setP(int pP);

	/// Function to set mI
	void setI(int pI);

	/// Function to set mK
	void setK(int pK);

	/// Function to set mAlpha
	void setAlpha(float pAlpha);

	/// Function to set mPriority
	void setPriority(float* pPriority);

	/// Function to get mSearchStack
	searchStack* getSearchStack();

	/// Function to get mAlgorithmInterface
	AlgorithmInterface* getAlgorithmInterface();

	/// Function to get mNoEqBranch
	int getNoEqBranch();

	/// Function to get mNoBranchAdded
	int getNoBranchAdded();
	/// Function to get mNoBranchIgnored
	int getNoBranchIgnored();
	/// Function to get mNoBranchReplaced
	int getNoBranchReplaced();
	/// Function to get mNoIterations
	int getNoIterations();

	/// Function to get the best path in search stack
	path* getBestPath();

	/// Function to clear the paths in search stack
	void clearSearchStack();

	/// Function to return search result
	void* getSolution();

private:

	/// Function to add a new path to the search stack when the stack is not full
	void addPath();

	/// Function to add a new path to the search stack when the stack is full
	void addPath_StackFull();

	/// Function to compute Multiplicate cost function
	cost compansatePathLengthMult(cost pPreCost, int pPathLength);

	/// Function to compute Adaptive-Additive cost function
	cost compansatePathLengthAdap( cost pPreCost, cost pOldPreCost,int pPathLength );

	/// Function to compute Adaptive*Multiplicate cost function
	cost compansatePathLengthAdapMul( cost pPreCost, cost pOldPreCost,int pPathLength );

	/// Function to perform one A* iteration
	int iterate();
	
	/// Function to perform necessary operations when a candidate branch is expanded
	void processBranch( int pBranchNo );

	/// Function to get a new SideInfo
	void* getNewSideInfo();

	/// Function to compute cost of a path
	cost ComputeCost(path* pPath, elementID pNewElementID);

	Trie mSearchTrie;			///< Search tree
	searchStack mSearchStack;	///< Search stack
	int mB;				///< Number of extensions per path
	int mP;				///< Number of maximum paths in search stack
	int mI;				///< Number of initial paths
	int mK;				///< Desired path length
	int mN;				///< Number of elements in the dictionary
	int mM;				///< Size of elements in the dictionary
	float mAlpha;		///< alpha for additive cost functions
	float mBeta;		///< beta for multiplicative cost functions
	cost* mCostCorrVec;		///< Cost correction vector for multiplicative cost function
	priority* mPriority;	///<  Pointer to the priority vector indexed by elementID's of dictionary vectors
	AlgorithmInterface mAlgInterface;	///< Algorithm interface to the search problem
	vector<void*> mFreeSideInfoList;	///< Vector holding unused SideInfos (for later use) 
	int mNoIterations;		///< Number of search iterations
	int mNoEqBranch;		///< Number of equivalent branches found during the search
	int mNoBranchAdded;		///< Number of branches added to the stack during the search
	int mNoBranchIgnored;	///< Number of branches ignored during the search
	int mNoBranchReplaced;	///< Number of branches replaced during the search
	//temporary for run...
	elementID* mCandList;	///< Temporary storage for best candidates
	path mTempBestPath;		///< Temporary storage for best path
	path mTempPath;			///< Temporary storage for a path
	elementID mActualCand;	///< elementID of the actual candidate
	AuxiliaryFunctionMode mAuxiliaryFunctionMode;	///<Auxiliary function mode

};																																
																												  