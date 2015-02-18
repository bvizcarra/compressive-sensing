#include "TrieNode.h"

// Default constructor
/// This constructor creates a node with mElementID=0 that is neither final, nor root, and that has no children nor parent.
TrieNode::TrieNode(void)
{
	mFinalFlag = false;
	mRootFlag = false;
	mElementID = 0;
	mParent = NULL;
}

// Constructor with an elementID
/// This constructor creates a node with mElementID=ID that is neither final, nor root, and that has no children nor parent.
TrieNode::TrieNode(elementID ID)
{
	mFinalFlag = false;
	mRootFlag = false;
	mElementID = ID;
	mParent = NULL;
}

// Default destructor
TrieNode::~TrieNode(void)
{
	childrenMapIter myIter;
	while(mChildren.size()!=0)
	{
		myIter = mChildren.begin();
		delete myIter->second;
	}
	if(mParent && !isRoot())
		mParent->removeChild(mElementID);
}

// Function to set a child to the node
/// This function adds pChild to mChildren of the node.
/// @param pChild pointer to the new child of the node
void TrieNode::setChild( TrieNode* pChild )
{
	mChildren.insert(pair<int,TrieNode*>(pChild->getElementID(),pChild));
}

// Function to check if an element is child of the node
/// @param pElementID elementID of the element questioned to be a child of the node
/// @return 1 if pElementID is a child of the node, 0 otherwise
bool TrieNode::isChild(elementID pElementID )
{
	if(mChildren.size() == 0)
		return 0;

	childrenMapIter myIter;
	myIter = mChildren.find(pElementID);
	if(myIter != mChildren.end())
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Function to return a child of the node
/// This function searchs among the children of the node for a node with pElementID.
/// If found, the function returns the pointer to this child. Otherwise, it returns NULL.
/// @param pElementID elementID of the element questioned to be a child of the node
/// @return pointer to the cihld if it exists, NULL otherwise
TrieNode* TrieNode::getChild(elementID pElementID)
{
	if(mChildren.size() == 0)
		return NULL;

	childrenMapIter myIter;
	myIter = mChildren.find(pElementID);
	
	if(myIter != mChildren.end())
		return myIter->second;
	else 
		return NULL;
}

// Function to set elementID of the node
/// @param pElementID new elementID
void TrieNode::setElementID(elementID pElementID )
{
	mElementID = pElementID;
}

// Function to get elementID of the node
/// @return elementID of the node
elementID TrieNode::getElementID()
{
	return mElementID;
}

// Function to make the node root
void TrieNode::makeRoot()
{
	mRootFlag = true;
}

// Function to get the parent of the node
/// @return pointer to the parent node
TrieNode* TrieNode::getParent()
{
	return mParent;
}

// Function to set the parent of the node
/// @param pParent pointer to the new parent node
void TrieNode::setParent( TrieNode* pParent )
{
	mParent = pParent;
}

// Function to set mFinal
/// @param pFinal new value for mFinal
void TrieNode::setFinal( bool pFinal )
{
	mFinalFlag = pFinal;
}

// Function to check if the node is final
/// @return 1 if the node is the final node in the path, 0 otherwise
bool TrieNode::isFinal()
{
	return mFinalFlag;
}

// Function to get children of the node
/// @return pointer to the mChildrenMap of the node
childrenMap* TrieNode::getChildren()
{
	return &mChildren;
}

// Function to check if the node is root
/// @return 1 if the node is root, 0 otherwise
bool TrieNode::isRoot()
{
	return mRootFlag;
}

// Function to disable the root flag
void TrieNode::disableRootFlag()
{
	mRootFlag = false;
}

// Function that removes a child of the node
/// This function removes the child with pElementID among the children of the node. 
/// It returns 1 if the node is removed, and 0 if not found among the children.
///	@param pElementID elementID of the child to be removed
/// @return 1 if the node is removed, and 0 if not found among the children
bool TrieNode::removeChild( elementID pElementID )
{
	return (bool)(mChildren.erase(pElementID));
}



