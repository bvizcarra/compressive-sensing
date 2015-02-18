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
#include "AStarDefinitions.h"
#include "BaseOMP.h"

/// This class is the interface of A* search to the search problem. 
/// BaseAStar class has no direct access to the search problem, but via this interface to provide flexibility for chancing the test problem.
/// Functions of this interface may be modified to call the functions of the desired problem.
/// Except this class, anything at the problem side is unknown to BaseAStar. Any info about paths are passed to the problem class via SideInfo 
/// struct that should be modified accordingly.
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
class AlgorithmInterface
{
public:
	/// Default constructor
	AlgorithmInterface(void);

	/// Default destructor
	~AlgorithmInterface(void);

	/// Function interface for initialization of the A* search tree
	int getInitialPaths( int pNoInitialPaths, vector<unsigned int*> *pNodeList, int &pNodesPerPath );

	/// Function interface for computing best B candidates for expansion of the path pPath
	void getBestCandidates(int pB, path* pPath, elementID* pCandList);

	/// Function interface for computing pre-cost of pPath after expansion with pNewElementID
	cost getPreCost( path* pPath, int pNewElementID);

	/// Function interface for computing priorities of dictionary elements for sorting nodes in a path
	void getPriorities(priority* pPriority);

	/// Function to set mProblem
	void setProblem(BaseOMP* pProblem);

	/// Function to allocate a new SideInfo
	void* getNewSideInfo();

	/// Function interface to perform any necessary operations to extract the solution after A* search is terminated
	void performPostOperations(path *pPath);

	/// Function interface to copy a SideInfo structure into another
	void copySideInfo(void *pSrc, void* pDst);

	/// Function interface to return the final solution
	void* getSolution();

	/// Function interface to reset a SideInfo structure
	void resetSideInfo(void* pSideInfo);

	/// Function interface to get the pre-cost of an empty path
	cost getInitialCost();

	bool isSearchComplete(path *pPath);

private:
	BaseOMP* mProblem; ///< pointer to the problem class


};
