#include "BaseAStar.h"

// Constructor
/// @param pB Number of extensions per path (mB)
/// @param pP Number of maximum paths in search stack (mP)
/// @param pI Number of initial paths (mI)
/// @param pK Desired path length (mK)
/// @param pN Number of elements in the dictionary (mN)
/// @param pM Size of elements in the dictionary (mM)
/// @param pAlpha alpha for additive cost functions (mAlpha)
/// @param pBeta beta for multiplicative cost functions (mBeta)
/// @param pAuxiliaryFunctionMode Auxiliary function mode (mAuxiliaryFunctionMode)
BaseAStar::BaseAStar(int pB, int pP, int pI, int pK, int pN, int pM,float pAlpha, float pBeta, 
					 AuxiliaryFunctionMode pAuxiliaryFunctionMode)
{
	mB = pB;
	mP = pP;
	mI = pI;
	mK = pK;
	mN = pN;
	mM = pM;
	
	mNoBranchIgnored = 0;
	mNoEqBranch = 0;
	mNoBranchAdded = 0;
	mNoBranchReplaced = 0;
	mNoIterations = 0;
	
	mPriority = new priority[mN];
	mSearchTrie.setPriority(mPriority);

	//for run
	mCandList = new elementID[mB];
	mFreeSideInfoList.reserve(mP+mI);

	mAuxiliaryFunctionMode = pAuxiliaryFunctionMode;

	mAlpha = pAlpha;
	mBeta = pBeta;
	mCostCorrVec = NULL;

	if (mAuxiliaryFunctionMode == MUL)
	{
		mCostCorrVec = new cost[mK];
		mCostCorrVec[0] = 1;
		for(int i = 1; i<mK; i++)
			mCostCorrVec[i]= mCostCorrVec[i-1]*mAlpha;
	}
}

// Destructor
BaseAStar::~BaseAStar(void)
{
	if(mCostCorrVec)
		delete mCostCorrVec;
	if(mPriority)
		delete mPriority;
	if(mCandList)
		delete mCandList;
	if(mSearchStack.size())
	{
		searchStackIter myIter;	
		for(myIter = mSearchStack.begin(); myIter!=mSearchStack.end(); myIter++)
		{
			if(myIter->second.mSideInfo)
				delete myIter->second.mSideInfo;
		}
	}

	deletevector(mFreeSideInfoList);
}

// Function to set mB
/// @param pB new value of mB
void BaseAStar::setB( int pB )
{
	mB = pB;
}

// Function to set mP
/// @param pP new value of mP
void BaseAStar::setP( int pP )
{
	mP = pP;
}

// Function to set mI
/// @param pI new value of mI
void BaseAStar::setI( int pI )
{
	mI = pI;
}
// Function to set mAlpha
/// @param pAlpha new value of mAlpha
void BaseAStar::setAlpha( float pAlpha )
{
	if(mAlpha > 0.0f && mAlpha<1.0f )
		mAlpha = pAlpha;
	else
		cout<<"Alpha should satisfy 0 < Alpha < 1 "<<endl;
}

// Function to compute multiplicate cost function
/// This function computes the multiplicative cost function for a path with path length pPathLenght and precost (i.e cost
/// without any path length compensation) pPreCost.
/// @param pPreCost	precost (cst without path length compensation) of the path
/// @param pPathLength length of the path
/// @return cost of the path after path length compensation with multiplicative cost model
cost BaseAStar::compansatePathLengthMult( cost pPreCost, int pPathLength )
{
	return pPreCost* mCostCorrVec[mK-pPathLength];
}

// Function to compute adaptive-additive cost function
/// This function computes the adaptive-additive cost function for a path with path length pPathLenght and precost (i.e cost
/// without any path length compensation) pNewPreCost, where pOldPreCost is the precost of the path before 
/// previos expansion of the path.
/// @param pPreCost	precost (cst without path length compensation) of the path
/// @param pOldPreCost	precost of the path before the previous expansion of the path
/// @param pPathLength length of the path
/// @return cost of the path after path length compensation with adaptive-additive cost model
cost BaseAStar::compansatePathLengthAdap( cost pPreCost, cost pOldPreCost, int pPathLength )
{
	return pPreCost - mBeta*(pOldPreCost - pPreCost)*(mK-pPathLength);
}

// Function to compute adaptive-multiplicative cost function
/// This function computes the adaptive-multiplicative cost function for a path with path length pPathLenght and precost (i.e cost
/// without any path length compensation) pNewPreCost, where pOldPreCost is the precost of the path before 
/// previos expansion of the path.
/// @param pPreCost	precost (cst without path length compensation) of the path
/// @param pOldPreCost	precost of the path before the previous expansion of the path
/// @param pPathLength length of the path
/// @return cost of the path after path length compensation with adaptive-multiplicative cost model
cost BaseAStar::compansatePathLengthAdapMul( cost pPreCost, cost pOldPreCost,int pPathLength )
{
	float decRate = pPreCost / pOldPreCost;
	float expDec = pow(mAlpha*decRate,(float)(mK-pPathLength));
	return pPreCost*expDec;

}

// Function to reset and initialize BaseAStar for a new search
/// This function resets and initializes BaseAStar instance for a new search.
/// It clears the search tree and the search stack. Then these are initialized to
/// contain mI initial paths, which are obtained from the search problem class via
/// mAlgorithmicInterface.
/// @return 1 if initialization is successful, 0 otherwise
int BaseAStar::initialize()
{
	//first clear all from last search
	mSearchTrie.clearTrie();
	clearSearchStack();
	mAlgInterface.getPriorities(mPriority);
	
	mNoBranchIgnored = 0;
	mNoEqBranch = 0;
	mNoBranchAdded = 0;
	mNoBranchReplaced = 0;
	mNoIterations = 0;

	//initialize the algorithm
	vector<elementID *> nodeList;
	nodeList.reserve(mI);
	
	int nodesPerPath;
	cost tempCost;
	TrieNode* tempNode;

	//get initial paths and place into stack

	if(mAlgInterface.getInitialPaths(mI, &nodeList, nodesPerPath))
	{
		for(int i=0;i<mI;i++)
		{
			//for(int j=0; j<nodesPerPath;j++)
			//	tempPath2.mMap.insert(pair<priority,elementID>(mPriority[nodeList[i][j]],nodeList[i][j]));
			//mTempPath.mSideInfo = getNewSideInfo();
			//mAlgInterface.resetSideInfo(mTempPath.mSideInfo);
			//mTempPath.mLeaf = mSearchTrie.addPath(mSearchTrie.getRootNode(), nodeList[i], nodesPerPath);
			//mTempPath.mPathLength = nodesPerPath;
			//mTempPath.mPreCost = mAlgInterface.getInitialCost();
			//cost tempCost;
			//for(int j =0; j<nodesPerPath; j++)
			//	tempCost = ComputeCost(&mTempPath, nodeList[i][j]);
			//mSearchStack.insert(pair<cost, path>(tempCost,mTempPath));
			//mNoBranchAdded++;

			mTempPath.mSideInfo = getNewSideInfo();
			mAlgInterface.resetSideInfo(mTempPath.mSideInfo);
			mTempPath.mPreCost = mAlgInterface.getInitialCost();
			mTempPath.mLeaf = mSearchTrie.getRootNode();
			for(int j =0; j<nodesPerPath; j++)
			{
				tempNode = mSearchTrie.addPath(mTempPath.mLeaf, nodeList[i][j]);
				if(tempNode  )
					mTempPath.mLeaf = tempNode;
				else
					mTempPath.mLeaf = mTempPath.mLeaf->getChild(nodeList[i][j]);
				mTempPath.mPathLength = j+1;
				tempCost = ComputeCost(&mTempPath, nodeList[i][j]);
			}
			mSearchStack.insert(pair<cost, path>(tempCost,mTempPath));
			mNoBranchAdded++;
		}
	}
	else
	{
		deletevector(nodeList);
		cerr<<"Initialization Failed: No initial paths returned!"<<endl;
		return 0;
	}

	deletevector(nodeList);
	return 1;
}

// Function to run a new search
/// This function runs a new search. Search is not run, and 0 is returned if search stack contains no initial paths.
/// After the search terminates, AlgorithmInterface::performPostOperations is called for the search problem to extract 
/// the solution from the SideInfo of the returned path. (This solution should be stored in the search problem class.)
/// @return 0 if search is not performed as there are no initial paths in search stack, 1 otherwise
int BaseAStar::run()
{
	if(mSearchStack.empty())
	{
		cerr<<"BaseAStar.run(): Search Stack contains no initial paths. Call BaseAStar.initialize first!"<<endl;
		return 0;
	}
	//	while((int)(mSearchStack.begin()->second.mPathLength) < mK)
		while (!mAlgInterface.isSearchComplete(&(mSearchStack.begin()->second)))
		{
			iterate();
			mNoIterations++;	
		}

	//deletevector(mFreeSideInfoList);
	mAlgInterface.performPostOperations(&(mSearchStack.begin()->second));
	return 1;
}

// Function to perform one A* iteration
/// This function performs one iteration of A* search. It gets the best path from the search stack, best candidates
/// from the search problem class via AlgorithmInterface class and expands these.
int BaseAStar::iterate()
{

	mTempBestPath = mSearchStack.begin()->second; //copy bestPath 
	mSearchStack.erase(mSearchStack.begin());  //and remove it from stack!
	mAlgInterface.getBestCandidates(mB,&mTempBestPath,mCandList);

	for(int branchNo=0;branchNo<mB;branchNo++)
	{
		mTempPath = mTempBestPath;
		mActualCand = mCandList[branchNo];
		TrieNode *newNode = mSearchTrie.addPath(mTempPath.mLeaf, mActualCand);
		if( newNode )
		{
			mTempPath.mPathLength++;
			mTempPath.mLeaf = newNode;
			processBranch(branchNo);
		}
		else
		{
			mNoEqBranch++;
		}
	}

	mFreeSideInfoList.push_back(mTempBestPath.mSideInfo);
	return 1;
}
// Function to add a new path to the search stack when the stack is not full
/// This function adds a new path to the search stack if the stack is not full.
/// The new path is passed to the function via class member mTempPath.
void BaseAStar::addPath()
{
	mTempPath.mSideInfo = getNewSideInfo();
	mAlgInterface.copySideInfo(mTempBestPath.mSideInfo, mTempPath.mSideInfo);
	mSearchStack.insert(pair<cost,path>(ComputeCost(&mTempPath,mActualCand),mTempPath));
}

// Function to compute cost of a path
/// This function computes the cost of pPath with the selected cost model mAuxiliaryFunctionMode.
/// It gets the precost (cost without any path length compensation) from the search problem class via AlgorithmInterface
/// and computes final cost wrt. mAuxiliaryFunctionMode.
/// @param pPath pointer to the path whose cost is inquired
/// @param pNewElementID elementID of the last node added to pPath
/// @return cost of pPath
cost BaseAStar::ComputeCost(path* pPath, elementID pNewElementID )
{
	switch(mAuxiliaryFunctionMode)
	{
	case MUL :
		{
			pPath->mPreCost = mAlgInterface.getPreCost(&mTempPath,pNewElementID);
			return compansatePathLengthMult(pPath->mPreCost, pPath->mPathLength);	
			break;
		}

	case ADAP :
		{
			cost oldPreCost =  pPath->mPreCost;
			pPath->mPreCost = mAlgInterface.getPreCost(pPath,pNewElementID);
			return compansatePathLengthAdap(pPath->mPreCost, oldPreCost, pPath->mPathLength);
			break;
		}
	case ADAPMUL :
		{
			cost oldPreCost =  pPath->mPreCost;
			pPath->mPreCost = mAlgInterface.getPreCost(pPath,pNewElementID);
			return compansatePathLengthAdapMul(pPath->mPreCost, oldPreCost, pPath->mPathLength);
			break;
		}

	default:
		{
			cerr<<"Invalid AuxiliaryFunctionMode in BaseAStar::ComputeCost()"<<endl;
			return -1;
			break;
		}
	}
	return -1;
}

// Function to get a new SideInfo
/// This function returns a pointer to a SideInfo instance that can be assigned to a new path in the search stack.
/// If there is a free SideInfo in mFreeSideInfoList, a pointer to it is returned. Otherwise, a new SideInfo is created
/// via AlgorithmInterface by the search problem and a pointer to it is returned. (The contents of the returned SideInfo
/// are not cleared.)
/// @return pointer to the new SideInfo struct
void* BaseAStar::getNewSideInfo()
{
	if((int)mFreeSideInfoList.size()>0)
	{
		void *mySideInfo;
		mySideInfo = mFreeSideInfoList.back();
		mFreeSideInfoList.pop_back();
		return mySideInfo;
	}
	else
	{
		return mAlgInterface.getNewSideInfo();	//new residue
	}
}

// Function to add a new path to the search stack when the stack is full
/// This function adds a new path to the search stack if the stack is full (i.e. has mP paths).
/// The new path is passed to the function via class member mTempPath. This path is added to the search stack iff 
/// its cost is lower than the worst path in the tree, which forces removal of the worst path from the stack.
void BaseAStar::addPath_StackFull()
{
	mTempPath.mSideInfo = getNewSideInfo();
	mAlgInterface.copySideInfo(mTempBestPath.mSideInfo, mTempPath.mSideInfo);
	cost pathScore = ComputeCost(&mTempPath,mActualCand);
	searchStackIter myIter = mSearchStack.end();
	myIter--;
	if(myIter->first >= pathScore)
	{	//if we are here, residue was surely used... otherwise, there cannot be mP paths.
		mFreeSideInfoList.push_back(myIter->second.mSideInfo);	//we will use this space later (avoid reallocation)
		mSearchStack.erase(myIter);
		mSearchStack.insert(pair<cost,path>(pathScore,mTempPath));
		mNoBranchAdded++;
	}
	else
	{
		mFreeSideInfoList.push_back(mTempPath.mSideInfo);
		mNoBranchIgnored++;
	}
}

// Function to get mSearchStack
/// @return pointer to mSearchStack
searchStack* BaseAStar::getSearchStack()
{
	return &mSearchStack;
}

// Function to perform necessary operations when a candidate branch is expanded
/// This function performs necessary operations to add an expanded branch to the search stack.
/// This new branch is either added to the search stack, or neglected wrt. pruning rules.
/// @param pBranchNo order of the branch in the ordered best candidate list
void BaseAStar::processBranch(int pBranchNo)
{
	if(pBranchNo==0)		//replace best path
	{	
		addPath();
		mNoBranchReplaced++;
	}
	else
	{		
		if((int)mSearchStack.size() <= mP)
		{
			addPath();
			mNoBranchAdded++;
		}
		else
		{
			addPath_StackFull();
		}
	}
}

// Function to get mAlgorithmInterface
/// @return pointer to mAlgInterface
AlgorithmInterface* BaseAStar::getAlgorithmInterface()
{
	return &mAlgInterface;
}

// Function to get mNoEqBranch
/// @return mNoEqBranch
int BaseAStar::getNoEqBranch()
{
	return mNoEqBranch;
}

// Function to get mNoBranchAdded
/// @return mNoBranchAdded
int BaseAStar::getNoBranchAdded()
{
	return mNoBranchAdded;
}
// Function to get mNoBranchIgnored
/// @return mNoBranchIgnored
int BaseAStar::getNoBranchIgnored()
{
	return mNoBranchIgnored;
}
// Function to get mNoBranchReplaced
/// @return mNoBranchReplaced
int BaseAStar::getNoBranchReplaced()
{
	return mNoBranchReplaced;
}
// Function to get mNoIterations
/// @return mNoIterations
int BaseAStar::getNoIterations()
{
	return mNoIterations;
}

// Function to get the best path in search stack
/// This function returns the best path in the search stack. As paths in mSearchStack are ordered wrt. ascending cost,
/// the first path in the stack is the best path.
/// @return pointer to the best(first) path in the search stack
path* BaseAStar::getBestPath()
{
	return &(mSearchStack.begin()->second);
}

// Function to clear the paths in search stack
/// This function clears the search stack.
void BaseAStar::clearSearchStack()
{
	searchStackIter myIter;
	for(myIter = mSearchStack.begin(); myIter!= mSearchStack.end(); myIter++)
	{		
		//delete myIter->second.mSideInfo;
		mFreeSideInfoList.push_back( myIter->second.mSideInfo);
	}
	mSearchStack.clear();
}

// Function to return search result
/// This function gets the final solution from the problem class via AlgorithmInterface as a void*, and returns it.
/// @return void pointer to the solution
void* BaseAStar::getSolution()
{
	return mAlgInterface.getSolution();
}

