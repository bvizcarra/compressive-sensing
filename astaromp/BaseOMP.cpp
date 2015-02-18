#include "BaseOMP.h"

// Constructor
/// This is the constructor function for BaseOMP class.
/// @param pK value for member variable mK (desired sparsity)
/// @param pM value for member variable mM (number of observations)
/// @param pN value for member variable mN (vector dimension)
/// @param pEps value for member variable mEps (error toleration)
/// @param pNodesPerInitPath value for member variable mNodesPerInitPath (initial path length)
BaseOMP::BaseOMP( int pK,int pM,int pN, float pEps, int pNodesPerInitPath )
{
	mK = pK;
	mM = pM;
	mN = pN;
	mNodesPerInitPath = pNodesPerInitPath;
	mEps = pEps;
	mSolution = new float[mN];
	mDictNorm = new float[mN];
}

// Destructor
/// This is the destructor for BaseOMP class.
BaseOMP::~BaseOMP(void)
{
	delete mSolution;
	delete mDictNorm;
}

// Function to solve for the sparse target vector from the QR decomposition.
/// This function solves for the sparse target vector from the QR decomposition provided by pR and pZ. 
/// Nonzero coefs (c) are solved from pR*c = pZ, where the upper triangular matrix pR and vector pZ have been
/// obtained via QR decomposition.  The nonzero coefs are placed in their positions in mSolution, where these 
/// locations are indexed by pIndList.
/// @param pIndList pointer to the vector containing locations of nonzero coefs in the target sparse vector.
/// @param pZ pointer to the z vector obtained from QR decomposition
/// @param pR pointer to the upper triangular matrix R from QR decomposition
void BaseOMP::solveCoefs( vector<unsigned int> *pIndList, float* pZ, float** pR )
{
	memset(mSolution,0,sizeof(float)*mN);
	int mSparsity = pIndList->size();
	for(int i = mSparsity-1;i>0;i--)
	{
		mSolution[pIndList->at(i)] = pZ[i]/pR[i][i];
		subtractProductScalarfromVector_I(pZ, pR[i],mSolution[pIndList->at(i)],i);
	}
	mSolution[pIndList->at(0)] = pZ[0]/pR[0][0];
}

// Function to add a new element to the QR decomposition of a path
/// This function updates the QR decomposition and residue of a path by adding a new element.
/// R is an upper triangular matrix. Q is an orthogonal matrix. z is the vector that holds inner products of the 
/// columns of Q with the observation vector y.
/// @param pNewElement pointer to the dictionary vector that is to be added to the representation
/// @param i number of elements in the representation (incuding the new element)
/// @param pQ pointer to the matrix Q (orthogonal matrix)
/// @param pR pointer to the matrix R (upper triangular matrix)
/// @param pZ pointer to the vector Z (\f$Z = Q^Ty\f$)
/// @param pRes pointer to the residue vector
void BaseOMP::addElementToRepresentation( float* pNewElement, int i, float** pQ, float* pR, float * pZ, float *pRes )
{
	copyVector(pNewElement,pQ[i],mM);

	//orthogonalize new element and update QR
	//R(1:i-1,i)=Q(:,1:i-1)'*new_element;
	for(int j=0;j<i;j++)
		pR[j] = computeInnerProd(pQ[j],pQ[i],mM);

	//Q(:,i)=new_element-Q(:,1:i-1)*R(1:i-1,i);
	for(int j=0;j<i;j++)
		subtractProductScalarfromVector_I(pQ[i],pQ[j],pR[j],mM);

	pR[i] = normalizeVector_I(pQ[i],mM);
	pZ[i] = computeInnerProd(pQ[i],pRes,mM);

	//update residue
	subtractProductScalarfromVector_I(pRes,pQ[i],pZ[i],mM);
}

// Function to set mDict
/// This function sets mDict. It computes length of the vectors in the dictionary
/// and stores them in mDictNorm.
void BaseOMP::setDict( float** pDict )
{
	mDict = pDict;
	for(int i = 0; i<mN; i++)
		mDictNorm[i] = l2Norm(mDict[i],mM);
}

// Function to set y
void BaseOMP::sety( float* py )
{
	my = py;
	mNorm_y = l2Norm(my,mM);
}

// Function to initialize paths for A*OMP
/// This function finds pNoInitialPaths initial paths for A*OMP. Length of these initial paths is returned in pNodesPerPath.
/// Pointers to arrays containing initial paths are stored in vector pNodeList.
/// @param pNoInitialPaths number of initial paths
/// @param pNodeList pointer to the vector of initial paths.
/// @param pNodesPerPath number of nodes in each path (returned by address)
/// @return 1 if initialization successfull, 0 otherwise
int BaseOMP::findInitialPaths( int pNoInitialPaths, vector<unsigned int*> *pNodeList, int &pNodesPerPath)
{
	pNodesPerPath = mNodesPerInitPath;
	if(mNodesPerInitPath == 1)
		return findInitialPaths1( pNoInitialPaths, pNodeList );
	else if (mNodesPerInitPath == 2)
		return findInitialPaths2( pNoInitialPaths, pNodeList );
	else 
	{
		cerr<<"Unsupported inital path length..."<<endl;
		return 0;
	}
}

// Function to find the initial paths (of length 1) for A*OMP
/// This function finds pNoInitialPaths initial paths for A*OMP wrt. correlation of dictionary elements to y.
/// Each initial path has length 1 and is stored as unsigned int * (elementID). All paths are stored in vector pNodeList. 
/// @param pNoInitialPaths number of initial paths
/// @param pNodeList pointer to the vector of initial paths.
/// @return 1 if initialization successfull, 0 otherwise
int BaseOMP::findInitialPaths1( int pNoInitialPaths, vector<unsigned int*> *pNodeList )
{	//safe to run?

	if((int)pNodeList->capacity() < pNoInitialPaths)
	{
		cerr<<"Not enough space reserved in pNodeList..."<<endl;
		return 0;
	}

	//find best candidates
	map<float,int,greater<float> >* tempCost = findClosestVectorsIndList( mDict, mDictNorm, my, mN, pNoInitialPaths, mM);

	//initialize pNodeList
	map<float,int,greater<float> >::iterator myIter = tempCost->begin();
	for(int i = 0;i<pNoInitialPaths;i++,myIter++ )
	{
		unsigned int* myInt = new unsigned int;
		*myInt = myIter->second;
		pNodeList->push_back(myInt);
	}
	delete tempCost;
	return 1;
}

/// Function to find the initial paths (of length 2) for A*OMP
/// This function finds pNoInitialPaths initial paths for A*OMP wrt. correlation of dictionary elements to y.
/// Each initial path has length 2. The first of these is fixed as the 0'th vector in the dictionary,
/// which is usually the DC term. ElementID's of the nodes in each path are stored in an unsigned int array. 
/// Pointers to arrays containing initial paths are stored in vector pNodeList. 
/// @param pNoInitialPaths number of initial paths
/// @param pNodeList pointer to the vector of initial paths.
/// @return 1 if initialization successfull, 0 otherwise
int BaseOMP::findInitialPaths2( int pNoInitialPaths, vector<unsigned int*> *pNodeList )
{	//safe to run?

	if((int)pNodeList->capacity() < pNoInitialPaths)
	{
		cerr<<"Not enough space reserved in pNodeList..."<<endl;
		return 0;
	}

	// add DC component
	float DCCorr = computeInnerProd(mDict[0], my, mM)/mDictNorm[0];
	float *res = new float[mM];
	subtractProductScalarfromVector(my,mDict[0],res,DCCorr/mDictNorm[0],mM);

	//find best candidates
	map<float,int,greater<float> >* tempCost = findClosestVectorsIndList( mDict, mDictNorm, res, mN, pNoInitialPaths, mM);
	
	//initialize pNodeList
	map<float,int,greater<float> >::iterator myIter = tempCost->begin();
	for(int i = 0;i<pNoInitialPaths;i++,myIter++ )
	{
		unsigned int* myInt = new unsigned int[2];
		myInt[0] = 0;
		myInt[1] = myIter->second;
		pNodeList->push_back(myInt);
	}
	delete tempCost;
	delete res;
	return 1;
}
// Function to find the best candidates for expansion of a path
/// This function finds and returns the pNoCand best candidates in the dictionary for the expansion of the path in pSideInfo.
/// Candidates are selected wrt. their correlations to the residue in pSideInfo and are returned in pCandList.
/// @param pNoCand number of requested candidates
/// @param pSideInfo pointer to the SideInfo struct containing info about the path to be expanded
/// @param pCandList pointer to the list of selected dictionary atoms
void BaseOMP::findBestCandidates( int pNoCand, SideInfo* pSideInfo, elementID* pCandList )
{
	findClosestVectorsIndList( mDict, mDictNorm, pSideInfo->mRes, pCandList, mN, pNoCand, mM );
}

// Function to compute the pre-cost of a path from SideInfo
/// This function computes the pre-cost (i.e. without path length compensation / auxiliary function) of the path
/// in SideInfo after addition of the pNewElementID. It adds pNewElementID to the representation SideInfo and 
/// performs orthogonal projection of the residue over the path via QR decomposition. Pre-cost of the new path is
/// the \f$l_2\f$ norm of the residue after orthogonal projection.
/// @param pSideInfo pointer to the SideInfo struct to which new element is added
/// @param pNewElementID ElementID of the new node in the path
/// @return pre-cost of the new path
float BaseOMP::computeCost( SideInfo* pSideInfo, int pNewElementID )
{
	int stepNo = (int)pSideInfo->mIndList.size();
	pSideInfo->mIndList.push_back(pNewElementID);
	addElementToRepresentation( mDict[pNewElementID], stepNo, pSideInfo->mQ, pSideInfo->mR[stepNo], pSideInfo->mZ, pSideInfo->mRes );
	return l2Norm(pSideInfo->mRes,mM);
}
// Function to compute priorities of dictionary members
/// This function computes the priorities of dictionary members wrt. their correlation to y.
/// Members with high correlation to y get higher priorities in the A* trie, which stores nodes 
/// in a path with descending priority.
/// @param pPriority pointer to the array holding priorities (in the same order as the vectors in mDict)
void BaseOMP::computePriorities( priority* pPriority )
{
	multimap<float,int,greater<float> > tempCost;
	
	for(int i=0; i<mN;i++)
		tempCost.insert(pair<float,int>(abs( computeInnerProd(mDict[i],my,mM)/mDictNorm[i]),i));

	multimap<float,int,greater<float> >::iterator tempCostIter = tempCost.begin();
	for(int i=0;i<mN;i++,tempCostIter++)
	{
		pPriority[tempCostIter->second]=mN-i;
	}
}
// Function to perform necessary operations after A* search is terminated
/// This function solves the reconstructed sparse vector from the path found by
/// A* search. Sparse coefficients of the vector are computed from the QR decomposition 
/// into in pSideInfo.
/// @param pSideInfo pointer to the SideInfo struct belonging to the solution found by BaseAStar
void BaseOMP::performPostOperations( SideInfo* pSideInfo )
{
	solveCoefs(&(pSideInfo->mIndList), pSideInfo->mZ, pSideInfo->mR);
}

// Function to allocate an instance of SideInfo struct 
/// This function allocates a new instance of the SideInfo struct.
/// @return pointer to the allocated instance
SideInfo* BaseOMP::allocateSideInfo()
{
	SideInfo* newSideInfo = new SideInfo;
	newSideInfo->mQ = allocateFloatMatrix(mK,mM);  //K column vectors of size M each
	newSideInfo->mR = new float*[mK];		//K column vectors with varying size (UT matrix)
	for(int i = 0; i<mK; i++)
	{
		newSideInfo->mR[i] = new float[i+1];		//keep nonzero entries only...
		memset(newSideInfo->mR[i],0,(i+1)*sizeof(float));
	}
	newSideInfo->mIndList.reserve(mK);
	newSideInfo->mZ = new float[mK];
	memset(newSideInfo->mZ,0,mK*sizeof(float));
	newSideInfo->mRes = new float[mM];
	memset(newSideInfo->mRes,0,mM*sizeof(float));
	return newSideInfo;
}

// Function to create a duplicate of a SideInfo struct
/// This function creates a duplicate of the SideInfo struct pSrc and 
/// returns a pointer to the duplicate.
/// @return pointer to the duplicate of pSrc
SideInfo* BaseOMP::duplicateSideInfo( SideInfo* pSrc )
{
	SideInfo* dst = allocateSideInfo();
	dst->mIndList=pSrc->mIndList;
	copyVector(pSrc->mRes, dst->mRes, mM);
	copyVector(pSrc->mZ,dst->mZ,mK);
	copyMatrix(pSrc->mQ,dst->mQ,mK,mM);
	for(int i = 0; i<mK; i++)
	{
		copyVector(pSrc->mR[i], dst->mR[i],i+1);
	}
	return dst;
}

// Function to copy the SideInfo struct pSrc into pDst
/// This function copies the contents of the SideInfo struct pSrc 
/// into pDst. This operation does not simply copy the pointers in
/// pSrc, but the arrays themselves. Space for all pointers in
/// pDst should have been allocated prior to calling this function.
void BaseOMP::copySideInfo( SideInfo* pSrc, SideInfo* pDst )
{
	pDst->mIndList = pSrc->mIndList;
	copyVector(pSrc->mRes, pDst->mRes, mM);
	copyVector(pSrc->mZ,pDst->mZ,mK);
	copyMatrix(pSrc->mQ,pDst->mQ,mK,mM);
	for(int i = 0; i<mK; i++)
	{
		copyVector(pSrc->mR[i], pDst->mR[i],i+1);
	}
}

// Function to return solution of A*OMP
/// This function returns the solution of A*OMP.
/// @return pointer to the reconstructed vector
float* BaseOMP::getSolution()
{
	return mSolution;
}

// Function to reset a SideInfo struct
/// This function resets the contents of pSideInfo. mIndList in pSideInfo
/// is cleared and residue (mRes) is set equal to y.
/// @param pSideInfo pointer to the SideInfo struct to be reseted
void BaseOMP::resetSideInfo( SideInfo* pSideInfo )
{
	pSideInfo->mIndList.clear();
	copyVector(my,pSideInfo->mRes,mM);
}

// Function to return \f$l_2\f$ norm of y
/// @return \f$l_2\f$ norm of y 
float BaseOMP::getNorm_y()
{
	return mNorm_y;
}

// Function to find a sorted list of vectors which lie closest to a vector among an array of vectors 
/// This function computes the inner products of pVector with the vectors stored in pVectorArray, and divides these by the norms of the
/// vectors in pVectorArray to find the set that matches pVector the most. It returns the indices of pReturnSize vectors which are 
/// the closest to pVector. These indices are stored in a map that is sorted wrt. decreasing
/// inner-product.
/// @param pVectorArray pointer to the vector array
/// @param pVectorNorm pointer to the array containing norms of the vectors in pVectorArray
/// @param pVector	pointer to the vector
/// @param pNoVectors number of vectors in pVectorArray
/// @param pReturnSize number of vectors to be returned
/// @param pSize length of vectors
/// @return map that stores indices of vectors having maximum inner-product with pVector.
map<float,int,greater<float> >* BaseOMP::findClosestVectorsIndList( float** pVectorArray, float* pVectorNorm, 
																 float* pVector, int pNoVectors, int pReturnSize, int pSize )
{
	map<float,int,greater<float> >* tempCost = new map<float,int,greater<float> >;
	map<float,int,greater<float> >::iterator tempCostIter;
	float tempCorr;
	for(int i=0; i<pNoVectors;i++)
	{
		tempCorr = computeInnerProd(pVectorArray[i],pVector,pSize)/pVectorNorm[i];
		tempCost->insert(pair<float,int>(abs(tempCorr),i));
		if((int)tempCost->size() > pReturnSize)
		{
			tempCostIter = tempCost->end();
			tempCostIter--;
			tempCost->erase(tempCostIter);
		}
	}
	return tempCost;
}

// Function to find a sorted list of vectors which lie closest to a vector among an array of vectors 
/// This function computes the inner products of pVector with the vectors stored in pVectorArray, and divides these by the norms of the
/// vectors in pVectorArray to find the set that matches pVector the most. It returns in pReturnList the indices of pReturnSize vectors which are 
/// the closest to pVector. Indices in pReturnList are sorted wrt. decreasing inner-product.
/// @param pVectorArray pointer to the vector array
/// @param pVectorNorm pointer to the array containing norms of the vectors in pVectorArray
/// @param pVector pointer to the vector
/// @param pReturnList pointer to the array that stores returned index list
/// @param pNoVectors number of vectors in pVectorArray
/// @param pReturnSize number of vectors to be returned
/// @param pSize length of vectors
void BaseOMP::findClosestVectorsIndList( float** pVectorArray, float* pVectorNorm, float* pVector, 
									   unsigned int* pReturnList ,int pNoVectors, int pReturnSize, int pSize )
{
	map<float,int,greater<float> > tempCost;
	map<float,int,greater<float> >::iterator tempCostIter;
	float tempCorr;
	for(int i=0; i<pNoVectors;i++)
	{
		tempCorr = computeInnerProd(pVectorArray[i],pVector,pSize)/pVectorNorm[i];
		tempCost.insert(pair<float,int>(abs(tempCorr),i));
		if((int)tempCost.size() > pReturnSize)
		{
			tempCostIter = tempCost.end();
			tempCostIter--;
			tempCost.erase(tempCostIter);
		}
	}
	tempCostIter = tempCost.begin();

	for(int i = 0; i<pReturnSize; i++, tempCostIter++)
		pReturnList[i] = tempCostIter->second;
}

// Function to check if the termination criteria of A*OMP is satisfied
/// This function checks if a specific path in search satisfies the termination criteria for A*OMP. The termination criteria
/// are the maximum desired sparsity (mK) or allowable amount of error (mEps) in the reconstruction of the 
/// observation my, whichever earlier satisfied. It returns true if it is satisfied and 
/// false otherwise. 
/// @param pPL		length of the path for which termination criteria are checked
/// @param pErr		observation error norm of the path for which termination criteria are checked
/// @return true if termination criterion is satisfied, false otherwise
bool BaseOMP::isSearchComplete( int pPL, float pErr )
{
	if ((pPL >= mK) || (pErr/mNorm_y <= mEps))
	{
		return true;
	}
	else
		return false;
}
