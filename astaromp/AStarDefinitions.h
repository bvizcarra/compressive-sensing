/*
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

/// AStarDefinitions is the collection of declarations of data structures that are 
/// used by the BaseAStar and BaseOMP classes in this A*OMP implementation.
///
///
/// Copyright 2011 Nazim Burak Karahanoglu
///
/// karahanoglu@sabanciuniv.edu,  burak.karahanoglu@gmail.com
#pragma once

#include "GlobalUtil.h"

using namespace std;

class TrieNode;

//************************************************************************//
// new types
typedef unsigned int priority;  ///< priority of an element
typedef unsigned int elementID; ///< ID of an element (in the dictionary)
typedef float cost;	///< cost of a path

/// map that holds children of a node ordered by ascending elementID
typedef map<elementID, TrieNode*, less<elementID> > childrenMap;
typedef childrenMap::iterator childrenMapIter;	///< iterator for childrenMap

/// Path defines a path in the search stack.
/// A path is stored by a pointer to its final node, and includes all ancestors of this node.
/// In addition to this pointer, the path's precost (cost without path length compensation),
/// its length and a pointer to the SideInfo, which contains necessary info for the problem class
/// are also stored.
struct path{
	TrieNode* mLeaf;	///< pointer to the leaf node of the path
	void* mSideInfo;	///< pointer to the SideInfo of the path (contents unknown to BaseAStar, used by BaseOMP)
 	unsigned int mPathLength;		///< length of a path
	float mPreCost;			///< pre-cost of a path: cost without path length compensation (without the auxiliary function)
};

typedef std::multimap<cost, path, less<cost> > searchStack;		///< list of active search paths ordered by ascending cost
typedef searchStack::iterator searchStackIter;		///<iterator for searchStack


/// enum that defines the auxiliary function type.
enum AuxiliaryFunctionMode
{
	ADAP,		///< Adaptive-additive auxiliary function
	MUL,		///< Aultiplicative auxiliary function	
	ADAPMUL		///< Adaptive-multiplicative auxiliary function
};
