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

#include "TrieNode.h"
#include <limits.h>

/// This class implements a trie structure to be used as A* search tree.
/// "Trie" class handles addition of nodes and paths to the trie in addition to ordering of nodes in a path wrt. their priorities. 
/// By sorting nodes wrt. a fixed priority, equivalent paths can be easily found and neglected. In addition, if the priorities are 
/// chosen appropriately, the size of the tree can be reduced by increasing the number of shared nodes. (i.e. 
/// we would like to keep only one track of nodes that are common in many paths by placing them close to the root.)
/// "Trie" only stores a pointer to the root node to access the tree. The root node has the highest priority, UINT_MAX.
// In order to add a node or nodes to a path, one should explicitely 
/// provide a pointer to the leaf node of the path to which the node/nodes should be added. Once a node/nodes are added, Trie class returns the leaf 
/// node to the new path, which should be handled outside this class to add further nodes to the same path.
/// Removing nodes has not been implemented as nodes are not removed at all. We keep the path in the trie eventhough it is removed from the search 
/// as we would like to avoid later equivalent paths. This does not affect active paths, as they are stored in SearchStack.
///
/// This class has been partially implemented by Umut Sen (umutsen@sabanciuniv.edu).
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
class Trie
{
public:
	/// Default constructor
	Trie(void);

	/// Default destructor
	~Trie(void);

	/// Function to set the root node
	void setRootNode(TrieNode* pRootNode);

	/// Function to get the root node
	TrieNode* getRootNode();

	// Check if the given path is already in trie:
	// ElementID's for the given path are passed in pPath (already sorted).
	// pos is the iterator of the path
	// if path exists return the last node
	// if it does not exist, return the last node of the existing part of the tree
	// pos contains the corresponding value (first pPath value that is not included in the trie)
	//TrieNode* checkEquivalentPath( pathMap *pPath,pathMapIter &pos);

	// Add the given path to the tree:
	// ElementID's for the given path are passed in pPath (already sorted).
	// check if path exists, 
	// if path exists, 
	//		if not final, make it final. return pointer to the last node. (path added)
	//		else if final, return NULL (path ignored.)
	// if does not exist, 
	//		add the path, make it final, return pointer to the last node. (path added)
	//TrieNode* addPath(pathMap *pPath);

	/// Function to add a new path by addition of a single node
	TrieNode *addPath( TrieNode* pLeafNode, elementID pNewNode );

	/// Function to clear the trie
	void clearTrie();

	/// Function to set the mPriority vector
	void setPriority(priority* pPriority);
	
private:

	/// Function that adds a new path below a node
	TrieNode *addSubPathToNode(TrieNode* pCurrentNode, vector<elementID> *pSubPath);
	
	/// Function that returns the priority of a node
	priority getNodePriority(TrieNode* pNode);

	TrieNode* mRoot;		///< pointer to the root node
	priority* mPriority;	///< pointer to the array holding priorities of dictionary elements
};
