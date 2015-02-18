#include "Trie.h"

// Default constructor
Trie::Trie(void)
{
	mRoot = new TrieNode;
	mRoot->makeRoot();
}

// Default destructor
Trie::~Trie(void)
{
	delete mRoot;
}

// Function to set the root node
/// This function sets pRootNode as root of the tree.
/// @param pRootNode pointer to the new root node
void Trie::setRootNode( TrieNode* pRootNode )
{
	mRoot = pRootNode;
	mRoot->makeRoot();
}

// Function to get the root node
/// @return pointer to the root
TrieNode* Trie::getRootNode()
{
	return mRoot;
}

// Function to add a new path by addition of a single node
/// This function adds a new path to the tree by adding a new node with elementID pNewNode 
/// to the path specified by the leaf node pLeafNode. If new path does not already exist, it is added to
/// the trie and the function returns a pointer to the new leaf. Otherwise, it returns NULL.
/// Starting with the leaf node, the function searches the nodes aboe pLeafNode for the correct location of the new node 
/// wrt. its priority and adds it as a child of the first node (say n1) which has higher priority than it if it is 
/// not already a child of it. Then it copies the nodes below n1 under the new node if they don't already exist. 
/// @param pLeafNode pointer to the leaf node of the path to which new node is to be added
/// @param pNewNode elementID of the new node
/// @returns pointer to the leaf node of the new path if a new path is added, NULL otherwise.
TrieNode *Trie::addPath( TrieNode* pLeafNode, elementID pNewNode )
{
	TrieNode* currentNode = pLeafNode;
	vector<elementID> subPath;
	while(getNodePriority(currentNode) < mPriority[pNewNode])
	{
		subPath.push_back(currentNode->getElementID());
		currentNode = currentNode->getParent();
	}
	if(getNodePriority(currentNode) == mPriority[pNewNode])		//if node is already on the trie, path is already considered...
	{
		return NULL;
	}
	else
	{
		subPath.push_back(pNewNode);
		return addSubPathToNode(currentNode, &subPath);
	}	
}

// Function that adds a new path below a node
/// This function adds a new path of nodes with elementID's in pSubPath as a child of pCurrentNode if
/// it does not already exist and returns a pointer to the new leaf node. If the nodes already exist under 
/// pCurrentNode, it returns NULL. Node ID's in pSubPath should be sorted in pSubPath wrt. descending 
/// priority.
/// @param pCurrentNode pointer to the leaf node under which the new path is to be added
/// @param pSubPath list of node ID's to be added wrt. descending priority
/// @return pointer to the new leaf node if a new path is added, NULL otherwise
TrieNode *Trie:: addSubPathToNode(TrieNode* pCurrentNode, vector<elementID> *pSubPath)
{
	while(pSubPath->size()>0)
	{
		//get next nodeID in pSubPath
		elementID myID = pSubPath->back();
		pSubPath->pop_back();
		// check if myID is already a child of pCurrentNode
		TrieNode *childNode = pCurrentNode->getChild(myID);
		if(childNode)
		{	
			pCurrentNode = childNode;
		}
		else
		{
			TrieNode *newNode = new TrieNode;
			newNode->setElementID(myID);
			newNode->setParent(pCurrentNode);
			pCurrentNode->setChild(newNode);
			pCurrentNode = newNode;
		}
	}
	if(pCurrentNode->isFinal())
	{
		return NULL;
	}
	else
	{
		pCurrentNode->setFinal(true);
		return pCurrentNode;
	}
}

// Function to clear the trie
/// This function clears the trie by deleting all nodes under the root (excluding root node.)
void Trie::clearTrie()
{
	childrenMapIter myIter;
	while(mRoot->getChildren()->size()!=0)
	{
		myIter = mRoot->getChildren()->begin();
		delete myIter->second;
	}
}

// Function to set the mPriority vector
/// This function sets the pointer to the priority vector mPriority equal to pPriority.
/// @param pPriority pointer to the array holding priorities indexed wrt. ascending Node ID's
void Trie::setPriority( priority* pPriority )
{
	mPriority = pPriority;
}

// Function that returns the priority of a node
/// This function returns the priority of pNode.
/// @param pNode pointer to the node whose priority is requested
/// @return priority of pNode
priority Trie::getNodePriority( TrieNode* pNode )
{	
	if(pNode->isRoot())
		return UINT_MAX;
	else
		return mPriority[pNode->getElementID()];

}



