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

#ifndef LB_HIST_APPROX_HPP
#define LB_HIST_APPROX_HPP

#include <list>
#include <iomanip>
#include "Node.hpp"

using namespace std;

class HistApprox
{
public:
	template <typename T>
	void loadBalance(Node<T> * octree, const decompo & nb1ers, const double & dist, double tol,
		const int & first, const int & last, const double & maxEdge, const vec3D & center, Gaspi_communicator & gComm, 
		double * nodeCenters, i64 * nodeOwners, int nbLeaves) const 
	{ 
		cout << "--> Approx Histogram load balancing" << endl;

		/**
		 * Compute the octree characteristics
		 **/
		
		// compute separators list
		ui32 height = octree->getHeight();
		double coeff = octree->getCoeff();
		double edge = maxEdge / (1<<(height));
		edge *= coeff;
		int nbGridAxis = (1 << height) - 1;
		int firstAxisIdx = nbGridAxis/2*(-1);

		// Allocate and init the grid of octree separators, in 3D, X, Y and Z
		double ** grid = new double* [3]();
		for(int i=0; i<3; i++)
			grid[i] = new double[nbGridAxis]();
		for (int i=0; i<3; i++)
		{
			int index = firstAxisIdx;
			for (int j=0; j<nbGridAxis; j++)
			{
				grid[i][j] = center[i] + (index*edge);
				index++;
			}
		}

		
		/**
		* Traverse the tree and recursively the leaves until reaching the targeted depth.
		* The targeted depth corresponds to the size of the prime numbers decomposition list.
		**/
				
		// KD Tree
		Node<T> * kdTree = nullptr;
		kdTree = new Node<T>();
		kdTree->setAttributes(octree->getFirstIndex(), octree->getNbItems(), octree->getEdge(), octree->getOrigin());

		// BFS list
		list <Node<T>*> bfsList;
		bfsList.push_back(kdTree);
		
		// flatten Indexes initialization
		int flatIdxSize = 2; 
		int * flatIdxes = new int [flatIdxSize];
		flatIdxes[0] = kdTree->getContent().getFirstIndex()-1;
		flatIdxes[1] = kdTree->getContent().getLastIndex();


		// traverse the the tree in BFS way
		while(!bfsList.empty())
		{
			// Get a pointer on the first node
			Node<T> * ptr = bfsList.front();
			
			// if it is a leaf and targeted depth is not reached -> recurse on the complete level
			if ( (ptr->isLeaf()) && (static_cast<unsigned int>(ptr->getDepth()) < nb1ers._list.size() ) ) 
			{			
				// create neighboring vectors
				vector< Node<T>* > levelNodes;
				
				// initialize neighboring vectors
				for (auto it = bfsList.begin(); it != bfsList.end(); it++)
					if ((*it)->getDepth() == ptr->getDepth())
					{
						levelNodes.push_back(*it);
					}			
				// call the load balance function, on a complete level, update an array of T
				T ** p;
				ptr->getContent().compSepHistApprox(ptr->getDepth(), nb1ers, levelNodes.size(), p, flatIdxes, flatIdxSize, edge, height, grid, nbGridAxis, nodeCenters, nodeOwners, nbLeaves); 

				// create all the children nodes, for the complete level
				int nbChilds = nb1ers._list[ptr->getDepth()];
				int nbWorkers = levelNodes.size();
				for (int i=0; i<nbWorkers; i++)
					levelNodes[i]->setChildren(p[i], nbChilds);			
				
				// insert all the children nodes into the bfs list, for the complete level
				for (int i=0; i<nbWorkers; i++)
					for (int j=0; j<nbChilds; j++)
						bfsList.push_back(levelNodes[i]->getChildren()[j]);
				
				// pop the current node
				bfsList.pop_front();
			}
			else
				bfsList.pop_front();
		}
		delete [] flatIdxes;

		// C to F +1 pour les ranks MPI au niveau des feuilles
		for (int i=0; i<nbLeaves; i++)
			nodeOwners[i]+= 1;

	}
};

#endif
