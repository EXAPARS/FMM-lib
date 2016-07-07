/*
  Copyright 2015 - UVSQ
  Authors list: Nathalie Möller, Eric Petit

  This file is part of the FMM-lib.

  FMM-lib is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later version.

  FMM-lib is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License along with
  the FMM-lib. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NODE_HPP
#define NODE_HPP

#include "mpi.h"
#include "GASPI.h"
#include "assert.h"
#include <vector>
#include <list>
#include <unistd.h>

#include "Decomposition.hpp"
#include "Gaspi_communicator.hpp"
#include "types.hpp"
#include "vec3D.hpp"

#include "fmm_tools.hpp"
// TEMPORARILY USE a STL STACK
#include <stack>

using namespace std;


template<typename T>
class Node
{
public:	
	static int cptLeaves;
	
protected:
	Node<T> * _parent;
	Node<T> ** _children;
	int _nbChildren;	
	int _depth;
	int64_t _id;	
	T _content;

public:
	// constructors
	Node<T>();
	Node<T>(const Node<T> & n);
	Node<T>(const T & content);
	
	// destructors
	~Node<T>();

	// read octree
	void read_octree(i64 * nbElemPerNode, i64 * firstElemAdress, i64 * nbSonsPerNode, i64 * firstSonId, double * centers);

	// getters & setters
	Node<T> * getParent() const { return _parent; }
	int getDepth() const { return _depth; }
	Node<T> ** getChildren() const { return _children; }
	int getNbChildren() const { return _nbChildren; }
	int64_t getId() const { return _id; }
	T getContent() const { return _content; }
	T & getContent() { return _content; }
	Node<T> * getNodePtr(const int64_t & id);
	Node<T> * getNodePtrF(const int64_t & id);	 // Fortran flat incomplete Octree
	Node<T> * goToTarget(const int64_t & id);
	Node<T> * getFirstLeafDescendant();
	Node<T> * getLastLeafDescendant();
	int getFirstIndex() const { return getContent().getFirstIndex(); }
	ui32 getEdge() const { return getContent().getEdge(); }
	vec3D getOrigin() const { return getContent().getOrigin(); }
	int getNbItems() {return getContent().getNbParticles(); }
	void setChildren(T * content, const int & nbChilds);
	void setContent(T * content) { _content = getContent(); }
	void setAttributes(const int & index, const int & nbParticles,  const int & edge, const vec3D & o);
	double getCoeff() const { return getContent().getCoeff(); }
	double getTranslate() const { return getContent().getTranslate(); }


	// tree functions
	int64_t computeId(const int & height, const int64_t & root, const int & LeafIndex);
	int64_t getAncesterAtLevel(const int & level);
	bool isLeaf(){ return (_nbChildren <= 0); } // spectre = -1 => LEAF
	int countLeaves();
	ui32 getHeight(); 

	// recursive Octree Division
	void divideOctree();
	void divideOctreeNTimes(const int & nTimes);
	void recDivideOctreeNbP(const int & MaxNbItems);
	void recDivideOctreeH(const int & Height);
	
	// tree leaves infos	
	void FillSendBuffer(int * buffer, int targetLevel);
	void FillSendBufferAndIds(int * buffer, i64 * IDs, int targetLevel);
	void traverseAndFillLeaves(int * buffer, int targetLevel);
	void traverseAndFillLeavesAndIDs(int * buffer, i64 * IDs, int targetLevel);
	void findLastParticleIndex(const int64_t & myNodeID, int & particleIndex);
	void recSearchLastParticleBeforeNode(Node<T> * treeHead, const int64_t & myID, const int & myDepth, int & value);
		
	// OpenGL tools
	vector<double> compBB() const;	
	
	// display
	void recDisplayInfoShort();
	void displayInfo();
	void displayInfoShort();
	
	// tests parcours avec pointeur sur fonction
	template < typename Fx >
	void traverse(Fx f, Node<T> * p);
};

// class static variable
template<typename T>
int Node<T>::cptLeaves = -1;

/*----------------------------------------------------------------------
							Constructors
----------------------------------------------------------------------*/

/**
* Default constructor.
*/
template<typename T>
Node<T>::Node()
	: _parent(nullptr)
	, _children(nullptr)
	, _nbChildren(0)
	, _depth(0)	
	, _id(0)
{}

/**
* Copy constructor.
*/
template<typename T>
Node<T>::Node(const Node<T> & n)
	: _parent(n.getParent())
	, _nbChildren(n.getNbChildren())
	, _depth(n.getDepth())
	, _id(n.getId())	
	, _content(n.getContent())
{
	// children
	if(n.getNbChildren()) 
	{
		_children = new Node<T>* [n.getNbChildren()];
		for (int i=0; i<_nbChildren; ++i)
			_children[i] = n.getChildren()[i];
	}
	else
		_children = nullptr;
}

/**
* Constructor with <T> content
*/
template<typename T>
Node<T>::Node(const T & content)
	: _parent(NULL)
	, _children(NULL)
	, _nbChildren(0)
	, _depth(0)
	, _id(0)
	, _content(content)	
{}

/*----------------------------------------------------------------------
							Destructors
----------------------------------------------------------------------*/
template<typename T>
Node<T>::~Node<T>()
{
	// if current node has children
	if(getChildren() != NULL) 
	{
		for (int i=0; i<getNbChildren(); ++i)
		{
			if (getChildren()[i])
				delete getChildren()[i];
		}
		_nbChildren = 0;
	}
	///FIXME : delete *this ?
}

/*----------------------------------------------------------------------
			Construct Octree by reading other octree format
----------------------------------------------------------------------*/
template<typename T>
void Node<T>::read_octree(i64 * nbElemPerNode, i64 * firstElemAdress, i64 * nbSonsPerNode, i64 * firstSonId, double * centers)
{
	i64 nodeID = this->_id;
	i64 nbChildren = nbSonsPerNode[nodeID];
	_nbChildren = nbChildren;
	
	if (nbChildren > 0)
	{
		i64 fSonID = firstSonId[nodeID]-1; // Fortran array indices start with 1 instead of 0
		
		// allocate the children nodes, update them and recursive call
		T * pTab = new Particles[nbChildren];
		_children = new Node<T>* [nbChildren];

		for (int i=0; i<nbChildren; ++i)
		{
			_children[i] = new Node<T>();
			_children[i]->_parent = this;				
			_children[i]->_depth = _depth + 1;
			_children[i]->_id = (fSonID + i);
			
			// Content update
			
			// if it is a Leaf, we have _first and _last indexes
			int first = -1;
			if (nbSonsPerNode[_children[i]->_id]<= 0) // this child is a leaf
			{
				first = firstElemAdress[fSonID+i]-1; // From F to C
			}
			
			pTab[i].setAttributes ( first, 							/*index (_first)*/
									nbElemPerNode[(fSonID+i)],  /*nbparticles*/
									0, 							/*edge*/
									vec3D(centers[(fSonID*3)], centers[(fSonID*3)+1], centers[(fSonID*3)+2]) /*origin*/
								   ); 
			_children[i]->_content = pTab[i];
			
			// recursive call
			_children[i]->read_octree(nbElemPerNode, firstElemAdress, nbSonsPerNode, firstSonId, centers);
		}

		delete [] pTab;
	}
}
/*----------------------------------------------------------------------
							Getters & Setters
----------------------------------------------------------------------*/

template<typename T>
Node<T> * Node<T>::getNodePtr(const int64_t & id)
{
	// Pointer on the root of the tree
	Node<T> * n = this;
	while(n->getId() != 0)
		n = n->getParent();

	// Compute path 
	int64_t * path = new int64_t [32];
	int64_t tmpID = id;
	
	if (id > 8)
	{
		// Compute the path down to the father
		int index = 0;
		while(tmpID > 8)
		{
			tmpID -= 1;
			tmpID = (tmpID >> 3);
			path[index] = tmpID;
			index++;
		}
		
		// Index ends 1 box too far
		index--; 
		
		int64_t sonID;
		// Go down to the father
		for (int i = index; i>=0; i--)
		{
			sonID = path[i] - (8*n->getId()) - 1;
			if (n->getChildren())			
				n = n->getChildren()[sonID];
			else
			{
				cerr << "Has no children !" ;
				exit(4);
			}
		}
	}
	
	// Pick the right child
	int childIndex = (id-1)%8;
	n = n->getChildren()[childIndex];
	
	return n;
}

template<typename T>
Node<T> * Node<T>::getNodePtrF(const int64_t & targetID)
{
	// Pointer on the root of the tree
	Node<T> * n = this;
	while(n->getId() != 0)
		n = n->getParent();

	// Initiate the DFS search
	if (targetID == 0)
		return n;
	else
	{
		n = goToTarget(targetID);
		return n;
	}
}

template<typename T>
Node<T> * Node<T>::goToTarget(const int64_t & targetID)
{
	stack<Node<T>*> st;
	Node<T> * tmp = nullptr;
	Node<T> * res = nullptr;
	bool found = false;
	
	st.push(this);
	while (!st.empty() && !found)
	{
		// pop a node and if it is the target
		tmp = st.top();
		st.pop();
		if (tmp->_id == targetID)
		{
			//printf( "Found the target node. id = %ld, depth = %d\n",tmp->getId(), tmp->getDepth());	
			res = tmp;
			found = true;
		}
		// push the children into the stack
		else
		{
			for (int i=0; i<tmp->getNbChildren(); i++)
				st.push(tmp->getChildren()[i]);
		}
	}
	return res;
}

template<typename T>
Node<T> * Node<T>::getFirstLeafDescendant()
{
	Node<T> * tmp = this;
	
	while (!tmp->isLeaf())
	{
		tmp = tmp->getChildren()[0];
	}
	
	return tmp;
}

template<typename T>
Node<T> * Node<T>::getLastLeafDescendant()
{
	Node<T> * tmp = this;
	
	while (!tmp->isLeaf())
	{
		tmp = tmp->getChildren()[tmp->getNbChildren()-1];
	}
	
	return tmp;
}


template<typename T>
void Node<T>::setChildren(T * content, const int & nbChilds)
{
	_nbChildren = nbChilds;
	_children = new Node<T>* [nbChilds];

	for (int i=0; i<nbChilds; ++i)
	{
		_children[i] = new Node();
		_children[i]->_parent = this;
		_children[i]->_depth = _depth + 1;
		_children[i]->_id = (_id*8) + (i+1); 
		_children[i]->_content = content[i];
	}
}

template<typename T>
void Node<T>::setAttributes(const int & index, const int & nbParticles,  const int & edge, const vec3D & o)
{
	getContent().setAttributes(index, nbParticles, edge, o);
}

/*----------------------------------------------------------------------
							Tree functions
----------------------------------------------------------------------*/
template<typename T>
int Node<T>::countLeaves()
{
	static int ctr = 0;
	if ( _children)
	{
		// recursive call
		for (int i=0; i<_nbChildren; ++i)
			_children[i]->countLeaves();
	}
	else
		ctr++;
	
	return ctr;
}

template<typename T>
int64_t Node<T>::computeId(const int & height, const int64_t & rootNodeID, const int & LeafIndex)
{
	int64_t firstSon = (((8*rootNodeID+1)*8)+1)*8+1;
	int64_t id = firstSon + LeafIndex;
	return id;
}

template<typename T>
int64_t Node<T>::getAncesterAtLevel(const int & level)
{
	int64_t ancester = -1;
	
	if (getDepth() <= level)
		exit(6);

	ancester = getId();
	int diff = getDepth() - level;
	while (diff > 0)
	{
		ancester = (ancester - 1) / 8;
		diff--;
	}
	
	return ancester;
}

template <typename T> 
ui32 Node<T>::getHeight()
{
	stack<Node<T>*> st;
	Node<T> * tmp = nullptr;

	int maxDepth = 0;
	st.push(this);
	
	while (!st.empty())
	{
		// pop a node
		tmp = st.top();
		st.pop();
		if (tmp->_depth > maxDepth)
		{
			maxDepth = tmp->_depth;
		}
		// push the children into the stack
		else
		{
			for (int i=0; i<tmp->getNbChildren(); i++)
				st.push(tmp->getChildren()[i]);
		}
	}
	return maxDepth;
}


/*----------------------------------------------------------------------
							OCTREE
----------------------------------------------------------------------*/

/**
* This method divides the node into 8 smaller nodes. \n
* It calls the divideOctree() method on the content. \n
* \warning : The generated T sub contents are \b copied into the children nodes.
*/
template<typename T>
void Node<T>::divideOctree()
{	
//DEBUG virer Gaspi	
	gaspi_rank_t rank;
	gaspi_proc_rank(&rank);
	
	// get the content of the node
	T p = getContent();
	
	// Declaration of an array of contents
	T * tabContents = nullptr;
	
	// Get the content and call division method on it. The result is an array of 8 sub-contents.
	tabContents = p.divideOctree();
	
	// Fill the children
	if (_children)
		verbose(rank, "There were already children here !!!");
		
	_children = new Node<T>* [8];
	_nbChildren = 8;

	
	for (int i=0; i<8; ++i)
	{
		_children[i] = new Node();
		_children[i]->_parent = this;				
		_children[i]->_depth = _depth + 1;
		_children[i]->_id = (_id*8) + (i+1); 
		_children[i]->_content = tabContents[i];
/*		if(tabContents[i].getNbParticles() > 0)
			dump_tree_add_child(rank, _id, _children[i]->_id, _children[i]->getNbItems());*/
	}	
	delete [] tabContents;
}	


/**
* This method divides the node recursively.
* It uses the max number of elements per code as threshold.
* @param int MaxNbItems : max number of elements per node.
*/
template<typename T>
void Node<T>::recDivideOctreeNbP(const int & MaxNbItems)
{	
	if (_content.getNbParticles() > MaxNbItems)
	{
		divideOctree();
		for (int i=0; i<8; ++i)
			_children[i]->recDivideOctreeNbP(MaxNbItems);
	}
}

/**
* This method divides the node recursively until reaching the targeted depth.
* @param int height : targeted depth.
*/
template<typename T>
void Node<T>::recDivideOctreeH(const int & height)
{	

	if (_depth < height)
	{
		divideOctree();
		for (int i=0; i<8; ++i)
			_children[i]->recDivideOctreeH(height);
	}
}

template<typename T>
void Node<T>::divideOctreeNTimes(const int & nTimes)
{
	if (nTimes)
	{
		divideOctree();
		for (int i=0; i<8; ++i)
			_children[i]->divideOctreeNTimes(nTimes-1);
	}
}

/*----------------------------------------------------------------------
							LEAVES INFO
----------------------------------------------------------------------*/

/**
* Initialize the static counter and calls the traverseAndFill Method. \n
*/
template<typename T>
void Node<T>::FillSendBuffer(int * buffer, int targetLevel)
{
	// init static counter
	cptLeaves = -1;	
	
	// traverse the tree and fill the buffer with the leaves information
	traverseAndFillLeaves(buffer, targetLevel);	
}

template<typename T>
void Node<T>::FillSendBufferAndIds(int * buffer, i64 * IDs, int targetLevel)
{
	// init static counter
	cptLeaves = -1;	
	
	// traverse the tree and fill the buffer with the leaves information
	traverseAndFillLeavesAndIDs(buffer, IDs, targetLevel);
}

/**
* Recursive method that fills a buffer with the number of particles in the leaves. \n
*/
template<typename T>
void Node<T>::traverseAndFillLeavesAndIDs(int * buffer, i64 * IDs, int targetLevel)
{
	if( getDepth() < targetLevel )
	{
		// Recursive calls
		for (int i=0; i<_nbChildren; ++i)
			_children[i]->traverseAndFillLeavesAndIDs(buffer, IDs, targetLevel);
	}
	else if ( getDepth() == targetLevel )
	{
		// Update buffer
		Node<T>::cptLeaves++;
		buffer[cptLeaves] = getContent().getNbParticles();
		IDs[cptLeaves] = getId();
	}
	else if ( getDepth() > targetLevel )
	{
		cerr << "[traverseAndFillLeavesAndIDs] depth > targetLevel" << endl;
		exit(1);
	}
}

template<typename T>
void Node<T>::traverseAndFillLeaves(int * buffer, int targetLevel)
{
	if( getDepth() < targetLevel )
	{
		// Recursive calls
		for (int i=0; i<_nbChildren; ++i)
			_children[i]->traverseAndFillLeaves(buffer, targetLevel);
	}
	else if ( getDepth() == targetLevel )
	{
		// Update buffer
		Node<T>::cptLeaves++;
		buffer[cptLeaves] = getContent().getNbParticles();
	}
	else if ( getDepth() > targetLevel )
	{
		cerr << "[traverseAndFillLeaves] depth > targetLevel" << endl;
		exit(1);
	}
}


template<typename T>
void Node<T>::findLastParticleIndex(const int64_t & myNodeID, int & particleIndex)
{
	// Get depth
	int depth = getNodePtr(myNodeID)->getDepth();
	
	// Get a pointer on the root of the tree
	Node<T> * treeHead = this;
	while(treeHead->getId() != 0)
		treeHead = treeHead->getParent();
		
	recSearchLastParticleBeforeNode(treeHead, myNodeID, depth, particleIndex);			
}

/**
 * Recursive search for the last node containing particles before the myNodeID node.
 * "Before" in the sense of "on the left side"
 * Since Morton order corresponds to a dfs traversal.
 **/
template<typename T>
void Node<T>::recSearchLastParticleBeforeNode(Node<T> * n, const int64_t & myNodeID, const int & myNodeDepth, int & value)
{
	if ( n->getNbChildren() > 0 )
	{
		// recursive call
		for (int i=0; i<n->getNbChildren(); ++i)
			n->getChildren()[i]->recSearchLastParticleBeforeNode(n->getChildren()[i], myNodeID, myNodeDepth, value);
	}
	else
	{
		if (n->getNbItems())
		{
			int depth = n->getDepth();
			int64_t id = n->getId();
			
			// same depth && id < myNodeID
			if ( ( depth == myNodeDepth) && (id < myNodeID) )
				value = n->getContent().getLastIndex();

			// deeper than me  &&  his ancestor < me
			else if ( depth > myNodeDepth )
			{ 
				int64_t hisAncesterAtMyLevel = n->getAncesterAtLevel(myNodeDepth);			
				if ( hisAncesterAtMyLevel < myNodeID )
					value = n->getContent().getLastIndex();
			}
			
			// higher than me  &&  id < my ancestor 
			else if ( depth < myNodeDepth) 
			{
				int64_t myAncesterAtHisLevel = getNodePtr(myNodeID)->getAncesterAtLevel(depth);
				if ( id < myAncesterAtHisLevel )
						value = n->getContent().getLastIndex();				
			}
		}
	}
}

/*----------------------------------------------------------------------
							Display
----------------------------------------------------------------------*/

/**
* Recursive method which traverses the tree until reaching the leaves. \n
* Displays the coordinates of the first nb T contents. \n
* @param : int nb : number of contents to display.
*/
template<typename T>
void Node<T>::recDisplayInfoShort()
{
	static int cptLeaves = -1;
	if ( _children)
	{
		// recursive call
		for (int i=0; i<_nbChildren; ++i)
			_children[i]->recDisplayInfoShort();
	}
	else
	{
		// display the leaves only
		cptLeaves++;
		cout << "L_" << cptLeaves << "  ";
		cout << "(" << _id << ") ";
		getContent().displayInfoShort();
	}
}

/**
* Recursive method which traverses the tree until reaching leaves. \n
* Displays the Origin, Edge, nbElements, and Index of the first element in the global array. \n
*/
template<typename T>
void Node<T>::displayInfo()
{
	getContent().displayInfo();
	cout << "Depth = " << _depth << endl;
}

template<typename T>
void Node<T>::displayInfoShort()
{
	cout << std::dec ;
	getContent().displayInfoShort();
}

/*----------------------------------------------------------------------
							OpenGL tools
----------------------------------------------------------------------*/

/**
* This methods computes the bounding box of the holded content from the node.
* @return : Returns a vector of 24 doubles, corresponding to the 8 3-dimensional coordinates.
*/
template<typename T>
vector<double> Node<T>::compBB() const
{
	vector<double> v;
	v = getContent().compBB();
	return v;
}

#endif
