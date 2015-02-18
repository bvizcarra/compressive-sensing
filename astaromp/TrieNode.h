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

/// TrieNode class defines a node in the A* search trie. 
/// Each node has a parent, some children, an elementID and root/final node flags.
/// This class has been partially implemented by Umut Sen (umutsen@sabanciuniv.edu).
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
class TrieNode
{
public:
	/// Default constructor
	TrieNode(void);
	/// Constructor with an elementID
	TrieNode(elementID ID);

	/// Default destructor
	~TrieNode(void);
	
	/// Function to set a child to the node
	void setChild(TrieNode* pChild);

	/// Function to check if an element is child of the node
	bool isChild(elementID pElementID);	

	/// Function to return a child of the node
	TrieNode* getChild(elementID pElementID);

	/// Function that removes a child of the node
	bool removeChild(elementID pElementID);

	/// Function to set elementID of the node
	void setElementID(elementID pElementID);

	/// Function to get elementID of the node
	elementID getElementID();

	/// Function to make the node root
	void makeRoot();

	/// Function to set mFinal
	void setFinal(bool pFinal);

	/// Function to check if the node is final
	bool isFinal();

	/// Function to get the parent of the node
	TrieNode* getParent();

	/// Function to set the parent of the node
	void setParent(TrieNode* pParent);

	/// Function to get children of the node
	childrenMap* getChildren();
	
	/// Function to check if the node is root
	bool isRoot();

	/// Function to disable the root flag
	void disableRootFlag();

private:
	//Trie node elements:
	bool mFinalFlag;	///< flag indicating the last node of a path
	bool mRootFlag;		///< flag indicating the root node of the trie
	TrieNode* mParent;	///< pointer to the parent node
	childrenMap mChildren;		///< children of the node
	elementID mElementID;  ///< ID of the element the node represents
};