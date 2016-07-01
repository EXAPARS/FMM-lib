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
	void loadBalance(Node<T> * tree, const decompo & nb1ers, const double & dist, double tol,
		const int & first, const int & last, const double & maxEdge, const vec3D & center, Gaspi_communicator & gComm, 
		double * nodeCenters, i64 * nodeOwners, int nbLeaves) const 
	{ 
		cout << "--> Approx Histogram load balancing" << endl; 

		// -- Temporarily call another function
		int under_construction = 1;
		if (under_construction == 1)
		{
			loadBalance2(tree, nb1ers, dist, tol, first, last, maxEdge, center, nodeCenters, nodeOwners, nbLeaves);
		}
		else // original function 
		{

			/**
			 * Compute the octree characteristics
			 **/	
			// edge = square side, has to be the upper power of 2 : 2^log2(dist) upper rounded
			// 31 = height of the octree if d = 1, since coordinates in [0 - 2^31]
			// h = height of the octree if size of a leaf is c
			// h = 31 - log2(c)= 31 - log2d;

			ui32 log2d = ceil(log2(dist));
			ui32 edge = (1 << log2d); 	//pow(2, log2d)	
			ui32 height = 31 - log2d;	

			/**
			* Traverse the tree and recursively the leaves until reaching the targeted depth.
			* The targeted depth corresponds to the size of the prime numbers decomposition list.
			**/
					
			// BFS list
			list <Node<T>*> bfsList;
			bfsList.push_back(tree);
			
			// flatten Indexes initialization
			int flatIdxSize = 2; 
			int * flatIdxes = new int [flatIdxSize];
			flatIdxes[0] = tree->getContent().getFirstIndex()-1;
			flatIdxes[1] = tree->getContent().getLastIndex();

			// traverse the the tree in BFS way
			while(!bfsList.empty())
			{
				// update ptr
				Node<T> * ptr = bfsList.front();
				
				// if it is a leaf and targeted depth is not reached -> recurse on the complete level
				if ( (ptr->isLeaf()) && (static_cast<unsigned int>(ptr->getDepth()) < nb1ers._list.size() ) ) 
				{			
					// create neighboring vectors
					vector< Node<T>* > levelNodes;
					
					// initialize neighboring vectors
					for (auto it = bfsList.begin(); it != bfsList.end(); it++)
						if ((*it)->getDepth() == ptr->getDepth())
							levelNodes.push_back(*it);
									
					// call the load balance function, on a complete level, update an array of T
					T ** p;
					ptr->getContent().compSepHistApprox(ptr->getDepth(), nb1ers, levelNodes.size(), 
						p, flatIdxes, flatIdxSize, edge, height); 

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

			/**
			* Exchange the particles
			**/
			
			// MPI Exchange
			tree->getContent().exchangeMPI(flatIdxes);
			
			// Dealloc
			delete [] flatIdxes;
		}
	}
	
	/// Under construction, Alternative version to branch with spectre
	template <typename T>
	void loadBalance2(Node<T> * tree, const decompo & nb1ers, const double & dist, double tol,
		const int & first, const int & last, const double & maxEdge, const vec3D & center, double * nodeCenters, i64 * nodeOwners, int nbLeaves) const 
	{ 
		//cout << "--> Rerouted to 'Under Construction' " << endl; 

		/**
		 * Compute the octree characteristics
		 **/

//FIXME : fix height !!!
		ui32 height = 8;//tree->getHeight(); 3 for mini-sphere
		cout << setprecision(9);
		
		// compute separators list
		double edge = maxEdge / (1<<(height));
		double coeff = tree->getCoeff();
		//cout << "--- LIBFMM --- edge : " << edge << endl;		

/* TODO : clean after testing*/
		edge *= coeff;
		//cout << "--- LIBFMM --- new edge with coeff : " << edge << endl;		
		int nbGridAxis = (1 << height) - 1;		
		int firstAxisIdx = nbGridAxis/2*-1;

		// in 3D, X, Y and Z
		double ** grid = new double* [3]();
		for(int i=0; i<3; i++)
			grid[i] = new double[nbGridAxis]();
		
		// complete grid array
		for (int i=0; i<3; i++)
		{
			int index = firstAxisIdx;
			for (int j=0; j<nbGridAxis; j++)
			{
				grid[i][j] = center[i] + (index*edge);
				index++;
			}
		}

		//affichage de debug
		int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		//if (rank == 0)
		//if (tree->getDepth()==0)
		/**for (int i=0; i<3; i++)
		{
			int index = 0;
			if (i==0)
				cout << "X axis" << endl;
			else if (i ==1)
				cout << "Y axis" << endl;
			else if (i==2)
				cout << "Z axis" << endl;
			else
				cout << "Unexpected axis" << endl;			

			for (int j=0; j<nbGridAxis; j++)
			{
				cout << grid[i][j] << " " << scaleBackDB(grid[i][j]) << endl;
				index++;
			}
		}**/
		
		/**
		* Traverse the tree and recursively the leaves until reaching the targeted depth.
		* The targeted depth corresponds to the size of the prime numbers decomposition list.
		**/
				
		// BFS list
		list <Node<T>*> bfsList;
		bfsList.push_back(tree);
		
		// flatten Indexes initialization
		//cout << "Allocate flatIdxes" << endl;
		int flatIdxSize = 2; 
		int * flatIdxes = new int [flatIdxSize];
		
//		cout << "flatIdxes has been allocated" << endl;
//		cout << "complete flatidxes" << endl;
		flatIdxes[0] = tree->getContent().getFirstIndex()-1;
		flatIdxes[1] = tree->getContent().getLastIndex();
/* TODO : On n'a peut-être plus besoin de flatIdxes, à voir + tard */		

		//cout << "traverse the tree in BFS way" << endl;/*
		// traverse the the tree in BFS way
		while(!bfsList.empty())
		{
			// update ptr
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
				//ptr->getContent().compSepHistApprox(ptr->getDepth(), nb1ers, levelNodes.size(), p, flatIdxes, flatIdxSize, edge, height); 
				ptr->getContent().compSepHistApprox2(ptr->getDepth(), nb1ers, levelNodes.size(), p, flatIdxes, flatIdxSize, edge, height, grid, nbGridAxis, nodeCenters, nodeOwners, nbLeaves); 

				/*if (rank == 0)
				{
					cout << "********************************" << endl;
					cout << "***** INTERMEDIATE RESULTS *****" << endl;
					cout << "********************************" << endl;
					for (int i=0; i<nbLeaves; i++)
					{
						cout <<"\t(" << i << ")" << nodeOwners[i];
						if (i%10 == 0)
							cout << endl;
					}	
				}*/

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

		//cout << " ---------------------------- Sortie de load balance2" << endl;
		/**
		* Exchange the particles
		**/

		// MPI Exchange
		//n->getContent().exchangeMPI(flatIdxes);
		//n->getContent().exchangeMPI2(flatIdxes);
		// Dealloc
		delete [] flatIdxes;

		// C to F +1 pour les ranks MPI au niveau des feuilles
		for (int i=0; i<nbLeaves; i++)
			nodeOwners[i]+= 1;

	}
};

#endif
