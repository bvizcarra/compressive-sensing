#include "AlgorithmInterface.h"

AlgorithmInterface::AlgorithmInterface(void)
{
}

AlgorithmInterface::~AlgorithmInterface(void)
{
}

// Function interface for initialization
///
/// It should provide pNoInitialPaths initial paths to the BaseAStar class in pNodeList and number of nodes per initial path in pNodesPerPath.
/// @param pNoInitialPaths number of initial paths in the search tree.
/// @param pNodeList Each element in pNodeList is a pointer to an array holding the ID's of nodes in an initial path.
/// @param pNodesPerPath number of nodes per initial path (set by the problem class)
/// @return 1 if initialization successful, 0 otherwise
int AlgorithmInterface::getInitialPaths( int pNoInitialPaths, vector<unsigned int*> *pNodeList, int &pNodesPerPath )
{
	return mProblem->findInitialPaths(pNoInitialPaths, pNodeList, pNodesPerPath);
}

void AlgorithmInterface::getBestCandidates( int pB, path* pPath, elementID* pCandList )
{
	mProblem->findBestCandidates(pB, (SideInfo*)pPath->mSideInfo, pCandList);
}

// Function interface for computing pre-cost of pPath after expansion with pNewElementID.
/// This function should provide BaseAStar the pre-cost (i.e.  without path length compensation via auxiliary function)
/// of pPath after expansion by the node pNewElementID.
/// If side info of the path (pPath->mSideInfo) is used, it should be updated in order to cover the new element.
/// The function may use any info in pPath, but should not alter anything except mSideInfo, which in fact contains all 
/// info that should be passed to the problem class.
/// @param pPath pointer to the path whose pre-score is to be computed
/// @param pNewElementID ID of the new element to be added to the path
/// @return computed pre-cost
cost AlgorithmInterface::getPreCost( path* pPath, int pNewElementID )
{
	return mProblem->computeCost((SideInfo*)pPath->mSideInfo, pNewElementID);
}

// Function interface for computing priorities of dictionary elements for sorting nodes in a path
/// This function should provide BaseAStar class with the priorities of dictionary elements. These priorities should be stored 
/// in pPriority, which is indexed by the order of elements in the dictionary. (i.e Dictionary element mDict[i] has priority pPriority[i].)
/// @param pPriority pointer to the array holding priorities for dictionary elements
void AlgorithmInterface::getPriorities( priority* pPriority )
{
	mProblem->computePriorities(pPriority);
}

void AlgorithmInterface::setProblem( BaseOMP* pProblem )
{
	mProblem = pProblem;
}

// Function interface to allocate a new SideInfo
/// This function should allocate a new SideInfo structure and return a pointer to it as a void*.
/// @return void* pointer to the allocated SideInfo
void* AlgorithmInterface::getNewSideInfo()
{
	return  (void*)(mProblem->allocateSideInfo());
}

// Function interface to perform any necessary operations after A* search is terminated
/// This function should perform the necessary operations on the solution after A* search is terminated. 
/// It should extract the solution from pPath and store it in the problem class.
/// @param pPath pointer to the optimum path found by A* search
void AlgorithmInterface::performPostOperations( path *pPath )
{
	mProblem->performPostOperations((SideInfo*)pPath->mSideInfo);
}


// Function interface to copy a SideInfo structure into another
/// This function should copy the contents of SideInfo pSrc into SideInfo pDst.
/// @param pSrc pointer to the source SideInfo to be copied
/// @param pDst pointer to the destination SideInfo
void AlgorithmInterface::copySideInfo( void *pSrc, void* pDst )
{
	mProblem->copySideInfo((SideInfo*)pSrc,(SideInfo*)pDst);
}

// Function interface to return the final solution
/// this function should return a pointer to the solution. Its type is unknown to BaseAStar, hence it is returned as void*.
/// @return pointer to the solution
void* AlgorithmInterface::getSolution()
{
	return (void*)(mProblem->getSolution());
}

// Function interface to reset a SideInfo structure
/// This function should reset pSideInfo for later use.
/// @param pSideInfo pointer to the SideInfo structure
void AlgorithmInterface::resetSideInfo( void* pSideInfo )
{
	mProblem->resetSideInfo((SideInfo*)pSideInfo);
}

// Function interface to get the pre-cost of an empty path
/// This function should return the pre-cost (i.e. without any path length compensation) of an empty path.
/// @return pre-cost of an empty path
cost AlgorithmInterface::getInitialCost()
{
	return mProblem->getNorm_y();
}

bool AlgorithmInterface::isSearchComplete( path *pPath )
{
	return mProblem->isSearchComplete(pPath->mPathLength, pPath->mPreCost);
}