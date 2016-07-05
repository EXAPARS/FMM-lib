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
  the FMM-lib. If not, see <http://www.gnu.ogit rg/licenses/>.
*/

#ifndef LB_MORTON_SYNC_MPI_HPP
#define LB_MORTON_SYNC_MPI_HPP

#include "Node.hpp"
#include "LBMortonBase.hpp"

using namespace std;

class MortonSyncMPI : public LBMortonBase
{
public:
	template <typename T>
	void loadBalance(Node<T> * n, const decompo & nb1ers, const double & dist, double tolerance, 
		const int & first, const int & last, const double & maxEdge, const vec3D & center, Gaspi_communicator & gComm, 
		double * nodeCenters, i64 * nodeOwners, int nbLeaves) const;
	
	template <typename T>
	void loadBalance2(Node<T> * n, const decompo & nb1ers, const double & dist, double tolerance, 
		const int & first, const int & last, Gaspi_communicator & gComm, i64 * nodeOwners) const;
	
	template <typename T>
	void initializeSeps(Node<T> * n, int *globalBuffer, const int & nbLeaves, const int & nbSeps, 
		int * targets, int * nbUntilNode, int64_t * sepNodes, const int64_t & rootNodeID, 
		const int & divHeight) const;	

	template <typename T>
	void initializeSeps2(Node<T> * n, int *globalBuffer, i64* IDs, const int & nbLeaves, const int & nbSeps, 
		int * targets, int * nbUntilNode, int64_t * sepNodes, const int64_t & rootNodeID, 
		const int & divHeight) const;

	template <typename T>
	void loopRefine(Node<T> * n, const int & nbSeps, int64_t * sepNodes, int * nbUntilNode, 
		double * diff, int * targets, const double & tolerance, const int & meanNbItems) const;

	template <typename T>
	void loopRefine2(Node<T> * n, const int & nbSeps, int64_t * sepNodes, int * nbUntilNode, 
		double * diff, int * targets, const double & tolerance, const int & meanNbItems, int * stop) const;

	template <typename T>		
	void refine(Node<T> * n, int64_t * sepNodes, int * nbUntilNode, double * diff, int * targets, 
		const double & tolerance) const;

	template<typename T>
	void refine2(Node<T> * n, int64_t * sepNodes, int * nbUntilNode, double * diff, int * targets, 
	const double & tolerance) const;
	
	template <typename T>
	void mortonParticlesExchange(Node<T> * n, const int & nbSeps, int64_t * sepNodes, 
		const int & first, const int & last) const;	
	
	template <typename T>		
	void fillOwnersArray(Node<T> *n, const int & nbSeps, int64_t * sepNodes, i64 * nodeOwners) const;
};


template <typename T>
void MortonSyncMPI::loadBalance(Node<T> * n, const decompo & nb1ers, const double & dist, double tolerance, 
		const int & first, const int & last, const double & maxEdge, const vec3D & center, Gaspi_communicator & gComm, 
		double * nodeCenters, i64 * nodeOwners, int nbLeaves) const
{ 
	// TEMPORARILY call lb2 for Fortran octree
	
	loadBalance2(n, nb1ers, dist, tolerance, first, last, gComm, nodeOwners);
	/*cout << "--> Synchronous Morton MPI load balancing" << endl; 

	// MPI init
	int wsize, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	

	// LB init
	int nbParts = wsize;
	int nbSeps = nbParts - 1;
			
	// Octree construction with height = 3	
	int64_t rootNodeID = 0;	
	int divHeight = 3;
	n->recDivideOctreeH(3);
	
	// Fill the send buffer
	int nbLeaves = 1 << 9;						//512
	int *sendBuffer = new int [nbLeaves]();
	n->FillSendBuffer(sendBuffer, 3);
		
	// sepNodes, nbUntilNode
	int64_t * sepNodes = new int64_t [nbSeps];	
	for(int i=0; i<nbSeps; i++)
		sepNodes[i] = -1;
	int * nbUntilNode = new int [nbSeps](); 
	double * diff = new double [nbSeps](); 	
	
	// Rank 0 gathers all the sendbuffers
	int *globalBuffer = nullptr;
	if (rank == 0)
		globalBuffer = new int [nbLeaves]();
	MPI_Reduce(sendBuffer, globalBuffer, nbLeaves, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	delete sendBuffer;
	
	// Rank 0 computes sumNbItems and broadcasts it
	int sumNbItems = 0;
	if (rank == 0)
		for (int i=0; i<nbLeaves; i++)
			sumNbItems += globalBuffer[i] ;
	MPI_Bcast(&sumNbItems, 1, MPI_INT, 0, MPI_COMM_WORLD);		

	// Everybody updates targets and meanNbItems
	int meanNbItems = sumNbItems / wsize;			
	int * targets = new int [nbSeps];
	for(int i=0; i<nbSeps; i++)
		targets[i] = meanNbItems*(i+1); 
		
	// Initialize the separators computation
	initializeSeps(n, globalBuffer, nbLeaves, nbSeps, targets, nbUntilNode, sepNodes, rootNodeID, divHeight);
	
	// Refine until reaching the tolerance constraint
	loopRefine(n, nbSeps, sepNodes, nbUntilNode, diff, targets, tolerance, meanNbItems);	
	
	// Exchange the particles
	mortonParticlesExchange(n, nbSeps, sepNodes, first, last);*/
}

template <typename T>
void MortonSyncMPI::loadBalance2(Node<T> * n, const decompo & nb1ers, const double & dist, double tolerance, 
		const int & first, const int & last, Gaspi_communicator & gComm, i64 * nodeOwners) const
{ 

	// MPI init
	int wsize; MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	

	// LB init
	int nbParts = wsize;
	int nbSeps = nbParts - 1;
			
	// Initialize the first Leaves buffer
	// The algorithm starts at level 2 -> 64 leaves
	int64_t rootNodeID = 0;	
	int divHeight = 2;	
	int nbLeaves = 1 << 6;	
	
	// buffer and corresponding IDs
	int * globalBuffer = nullptr;
	i64 * IDsBuffer = nullptr;
	if (rank == 0)
	{
		// alloc
		globalBuffer = new int [nbLeaves]();
		IDsBuffer = new i64[nbLeaves]();
		
		// fill with information from level divHeight
		n->FillSendBufferAndIds(globalBuffer, IDsBuffer, divHeight);
		
		// SPECTRE : everybody already has the same informations, so quantity is multiplied by wsize.
		for (int i=0; i<nbLeaves; i++)
			globalBuffer[i] *= wsize;
	} 
	
	// Everybody updates targets and meanNbItems, SPECTRE : same input -> meanNbItems=number of items per rank
	int meanNbItems = n->getNbItems();
	int * targets = new int [nbSeps];
	for(int i=0; i<nbSeps; i++)
		targets[i] = meanNbItems*(i+1); 

	// sepNodes, nbUntilNode and stop criteria
	int64_t * sepNodes = new int64_t [nbSeps];	
	for(int i=0; i<nbSeps; i++)
		sepNodes[i] = -1;
	int * nbUntilNode = new int [nbSeps](); 
	double * diff = new double [nbSeps]();
	int * stop = new int [nbSeps]();

	// Initialize the separators computation
	initializeSeps2(n, globalBuffer, IDsBuffer, nbLeaves, nbSeps, targets, nbUntilNode, sepNodes, rootNodeID, divHeight);

	// Refine until reaching the tolerance constraint
	loopRefine2(n, nbSeps, sepNodes, nbUntilNode, diff, targets, tolerance, meanNbItems, stop);	
	
	// complete the array of leave owners with the new computed owners
	fillOwnersArray(n, nbSeps, sepNodes, nodeOwners);
	
	// delete arrays
	if (rank == 0)
	{
		delete globalBuffer;
		delete IDsBuffer;
	}
}

template<typename T>
void MortonSyncMPI::initializeSeps2(Node<T> * n, int *globalBuffer, i64 * IDs, const int & nbLeaves, const int & nbSeps, 
	int * targets, int * nbUntilNode, int64_t * sepNodes, const int64_t & rootNodeID, const int & divHeight) const
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Rank 0 initializes all separators
	if (rank == 0)
		computeMortonSeps2(n, globalBuffer, IDs, nbLeaves, nbSeps, targets, nbUntilNode, sepNodes, rootNodeID, divHeight);
	
	// Broadcasts sepNodes and nbBeforeBox
	MPI_Bcast(sepNodes, nbSeps, MPI_INT64_T, 0, MPI_COMM_WORLD);
	MPI_Bcast(nbUntilNode, nbSeps, MPI_INT, 0, MPI_COMM_WORLD);
}

template<typename T>
void MortonSyncMPI::initializeSeps(Node<T> * n, int *globalBuffer, const int & nbLeaves, const int & nbSeps, 
	int * targets, int * nbUntilNode, int64_t * sepNodes, const int64_t & rootNodeID, const int & divHeight) const
{
	int rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Rank 0 initializes all separators
	if (rank == 0)
	{
		computeMortonSeps(n, globalBuffer, nbLeaves, nbSeps, targets, nbUntilNode, sepNodes, rootNodeID, divHeight);
		delete globalBuffer;
	}
	
	// Broadcasts sepNodes and nbBeforeBox
	MPI_Bcast(sepNodes, nbSeps, MPI_INT64_T, 0, MPI_COMM_WORLD);
	MPI_Bcast(nbUntilNode, nbSeps, MPI_INT, 0, MPI_COMM_WORLD);
}




template<typename T>
void MortonSyncMPI::loopRefine(Node<T> * n, const int & nbSeps, int64_t * sepNodes, int * nbUntilNode, 
	double * diff, int * targets, const double & tolerance, const int & meanNbItems) const
{
	int finished = 0;
	
	while (!finished)
	{		
		// update diff array - (diff without box)
		for (int i=0; i<nbSeps; i++)
			diff[i] = abs((targets[i] - nbUntilNode[i]))/(meanNbItems*1.0);		
		
		// test if there is a sep to refine
		int refinement = 0;
		for (int i=0; i<nbSeps; i++)
			if (diff[i] > tolerance)
			{
				refinement = 1;
				break;
			}
		
		// refine if necessary AND if Possible
		if (refinement)
		{	
			refine(n, sepNodes, nbUntilNode, diff, targets, tolerance);
		}
		else
		{
			finished = 1;
		}
//		cout << "While (! finished) - loop again" << endl;
	}
}

template<typename T>
void MortonSyncMPI::loopRefine2(Node<T> * n, const int & nbSeps, int64_t * sepNodes, int * nbUntilNode, 
	double * diff, int * targets, const double & tolerance,  const int & meanNbItems, int * stop) const
{
	// MPI DEBUG
	int wsize, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int finished = 0;
	
	while (!finished)
	{
		// update diff array - (diff without box)
		for (int i=0; i<nbSeps; i++)
			diff[i] = abs((targets[i] - nbUntilNode[i]))/(meanNbItems*1.0);		
		
		// update stop array - test if last level is reached
		Node<T> * tmp;
		for (int i=0; i<nbSeps; i++)
		{
			if (!stop[i])
			{
				tmp = n->getNodePtrF(sepNodes[i]);
				if (tmp->isLeaf())
					stop[i] = 1;
			}
		}
		// display sepNodes
		/**for (int i=0; i<nbSeps; i++)
		{
			cout << "[" << rank << "] " << "sepNode : " << i << " " << sepNodes[i] << " " << stop[i]<< endl;
		}**/

		// test if there is a sep to refine
		int refinement = 0;
		for (int i=0; i<nbSeps; i++)
		{
			if (!stop[i])
			{
				if (diff[i] > tolerance)
				{
					refinement = 1;
					break;
				}
			}
		}
		
		// refine if necessary AND if Possible
		if (refinement)
		{	
			refine2(n, sepNodes, nbUntilNode, diff, targets, tolerance);
		}
		else
		{
			finished = 1;
		}
		/**cout << "[" << rank << "] " << "In LoopRefine2, finished = " << finished << endl;**/
	}
	
	/** TODO --> Ameliorer la suite pour ne pas essayer de raffiner des noeuds deja finis */
}



template<typename T>
void MortonSyncMPI::refine(Node<T> * n, int64_t * sepNodes, int * nbUntilNode, double * diff, int * targets, 
	const double & tolerance) const
{	
	// MPI init
	int wsize, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Get global information on which nodes are to refine
/* TODO - Simplifier, il n'est pas nécessaire de faire un Allgather, tout le monde a toutes les infos necessaires **/

	int64_t * NodesToRefine = new int64_t[wsize];
	int mySep = rank-1;
	int64_t myNode = -1;
	if ((rank > 0) &&  (diff[mySep] > tolerance))
		myNode = sepNodes[mySep];

	MPI_Allgather(&myNode, 1, MPI_INT64_T, NodesToRefine, 1, MPI_INT64_T, MPI_COMM_WORLD);
	
	// Everybody makes the necessary refinement and sends their information to the handler
	int nbLeaves = 512; // 3 LEVELS
	int * sendBuffer = new int[nbLeaves];
	int * globalBuffer = new int[nbLeaves];
	int divHeight = 3;

	// Send global information to each node to refine handler
	int64_t nodeID;
	Node<T> * node;
	for (int i=0; i<wsize; i++)
		if(NodesToRefine[i] >= 0)
		{
			// get the node
			nodeID = NodesToRefine[i];
			node = n->getNodePtr(nodeID);
//			cout << rank << " From refine function, ID : " << node->getId() << ", level : " << node->getDepth() << endl;
			// refinement
			node->divideOctreeNTimes(3);
			
			if (node->getChildren() > 0)
			{
				// compute leaves buffer
				node->FillSendBuffer(sendBuffer, node->getDepth() + divHeight);
				
				// Reduce buffer on the handler
				MPI_Reduce(sendBuffer, globalBuffer, nbLeaves, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
			}
		}

	// Each handler makes the refinement
	if (myNode >= 0)
		computeMortonOneSep(n, globalBuffer, nbLeaves, targets[rank-1], nbUntilNode[rank-1], sepNodes[rank-1], divHeight);

/// TODO : supprimer le broadcast de nbUNtilNode si ça ne sert à rien
	// Broadcast all new sepNodes and nbBeforeBox
	for (int i=1; i<wsize; i++)
	{
		MPI_Bcast(&sepNodes[i-1], 1, MPI_INT64_T, i, MPI_COMM_WORLD);
		MPI_Bcast(&nbUntilNode[i-1], 1, MPI_INT, i, MPI_COMM_WORLD);		
	}
	delete sendBuffer;
	delete globalBuffer;
}

template<typename T>
void MortonSyncMPI::refine2(Node<T> * n, int64_t * sepNodes, int * nbUntilNode, double * diff, int * targets, 
	const double & tolerance) const
{	
	// MPI init
	int wsize, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Get global information on which nodes are to refine
	int64_t * NodesToRefine = new int64_t[wsize];
	int mySep = rank-1;
	int64_t myNode = -1;
	
	// Test the node
	Node<T> * tmpNode;
	if ((rank > 0) &&  (diff[mySep] > tolerance))
		myNode = sepNodes[mySep];

/* TODO - Simplifier, il n'est pas nécessaire de faire un Allgather, tout le monde a toutes les infos necessaires **/
	MPI_Allgather(&myNode, 1, MPI_INT64_T, NodesToRefine, 1, MPI_INT64_T, MPI_COMM_WORLD);
	
	// Everybody makes the necessary refinement and sends their information to the handler
	int nbLeaves = 8; // 1 LEVEL ONLY
	int * sendBuffer = new int[nbLeaves];
	int * globalBuffer = new int[nbLeaves];
	i64 * IDsBuffer = new i64[nbLeaves]();
	i64 * myIDsBuffer = new i64[nbLeaves]();
	
	int divHeight = 1;

	// Send global information to each node to refine handler
	int64_t nodeID;
	Node<T> * node;
	for (int i=0; i<wsize; i++)
		if(NodesToRefine[i] >= 0)
		{
			// get the node
			nodeID = NodesToRefine[i];
			node = n->getNodePtrF(nodeID);		// TEMPORARILY : NEW FUNCTION TO GET THE POINTER ON NODE, SINCE TREE IS NOT COMPLETE
//			cout << rank << " From refine function, ID : " << node->getId() << ", level : " << node->getDepth() << endl;
			
			// refinement
			if (node->getChildren() > 0)
			{
				// compute leaves buffer
				//node->FillSendBuffer(sendBuffer, node->getDepth() + divHeight);
				node->FillSendBufferAndIds(sendBuffer, IDsBuffer, node->getDepth() + divHeight);
//				cout << rank << " From refine function, fill send buffer -> ok" << endl;
				
				for (int ii=0; ii<node->getNbChildren() ; ii++)
				{
//					cout << "child : " << ii << " id : " << node->getChildren()[ii]->getId() << endl;
//					cout << "nbElements : " << sendBuffer[ii] << ", nodeID : " << IDsBuffer[ii] << endl;
				}
				
				// Reduce buffer on the handler
				/*-- Ce reduce ne sert à rien, virer, multiplier les quantités par wsize FIXME**/
				MPI_Reduce(sendBuffer, globalBuffer, node->getNbChildren(), MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
//				cout << rank << " From refine function, mpi reduce  on " << i << "--> ok" << endl;
				
				// The handler keeps its IDs Buffer
				if (i==rank)
					for (int ii=0; ii<node->getNbChildren(); ii++)
						myIDsBuffer[ii] = IDsBuffer[ii];
				
			}
/*			else
			{
				cout << rank << " Impossible to refine that node, is already a LEAF" << endl;
			}*/
		}

	// Each handler makes the refinement
//	cout << rank << " Before Morton Comp One Sep, nbUntilNode : " << nbUntilNode[rank-1] << endl;
	


	if (myNode >= 0)
	{
		Node<T> * myNodePtr = n->getNodePtrF(myNode);
		computeMortonOneSep2(myNodePtr, globalBuffer, myIDsBuffer, myNodePtr->getNbChildren(), targets[rank-1], nbUntilNode[rank-1], sepNodes[rank-1], divHeight);
	}
//	cout << rank << " After Morton Comp One Sep " << endl;


/// TODO : supprimer le broadcast de nbUNtilNode si ça ne sert à rien
	// Broadcast all new sepNodes and nbBeforeBox
	for (int i=1; i<wsize; i++)
	{
		MPI_Bcast(&sepNodes[i-1], 1, MPI_INT64_T, i, MPI_COMM_WORLD);
		MPI_Bcast(&nbUntilNode[i-1], 1, MPI_INT, i, MPI_COMM_WORLD);		
	}

//	cout << rank << " After BroadCasts " << endl;

	delete sendBuffer;
	delete globalBuffer;
	
//	cout << rank << " After Delete " << endl;
}

template<typename T>
void MortonSyncMPI::mortonParticlesExchange(Node<T> * n, const int & nbSeps, int64_t * sepNodes, 
	const int & first, const int & last) const
{	
	// prepare flatIdxes array
	int flatIdxSize = nbSeps + 2;	
	int * flatIdxes = new int[flatIdxSize];
	
	flatIdxes[0] = first-1;
	for (int i=0; i<nbSeps; i++)
	{
		int lastParticleIndex = n->getNodePtr(sepNodes[i])->getContent().getLastIndex();

		// if there is no particle in this node, find the last particle before
		if (lastParticleIndex == -1)
		{	
			n->findLastParticleIndex(sepNodes[i], lastParticleIndex);
		}

		flatIdxes[i+1] = lastParticleIndex;

	}
	flatIdxes[flatIdxSize-1] = last;
	
	// call MPI exchange function
	n->getContent().exchangeMPI(flatIdxes);
	
	delete [] flatIdxes;		
}


template<typename T>
void MortonSyncMPI::fillOwnersArray(Node<T> * n, const int & nbSeps, int64_t * sepNodes, i64 * nodeOwners) const
{
//	cout << "****************** NEW GAME !!! *******************" << endl;
	
	// Extremities IDs
	Node<T> * firstLeaf = n->getFirstLeafDescendant();
	Node<T> * lastLeaf = n->getLastLeafDescendant();
	i64 firstLeafID = firstLeaf->getId();
	i64 lastLeafID = lastLeaf->getId();
	int wsize; MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	
	// from 0 to first Leaf :
	for (i64 i=0; i<firstLeafID; i++)
		nodeOwners[i]=0;
	
	i64 lastSepID = firstLeafID;
	
	// from firstSep to lastSep
	for (int i=0; i<nbSeps; i++)
	{
		Node<T> * sepNode = n->getNodePtrF(sepNodes[i]);
		i64 nextSepID = sepNode->getId();
		
		if (!sepNode->isLeaf())
			nextSepID = sepNode->getLastLeafDescendant()->getId();

		for (i64 j=lastSepID; j<=nextSepID; j++)
			nodeOwners[j] = i + 1;

		lastSepID=nextSepID + 1;
	}

	// from lastSep to last Leaf
	for (i64 i=lastSepID; i<=lastLeafID; i++)
		nodeOwners[i]=wsize;	
}


#endif
