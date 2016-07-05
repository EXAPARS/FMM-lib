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

#ifndef LB_MORTON_BASE_HPP
#define LB_MORTON_BASE_HPP

#include "Node.hpp"

#include <iostream>
using namespace std;

class LBMortonBase{
public:	

	template<typename T>		
	void computeMortonOneSep(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbLeaves, 
		const int & target, int & nbUntilNode, int64_t & separator, const int & divHeight) const;

	template<typename T>
	void computeMortonSeps(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbLeaves, 
		const int & nbSeps, int * targets,  int * nbUntilNode, int64_t * sepNodes, 
		const int64_t & rootNodeID, const int & divHeight) const;
	
	// Keep this version for the Gaspi asynchronous Morton version
	// TODO : to be updated
	template<typename T>
	void computeMortonOneSepG(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
		const int & target, int & nbUntilNode, int64_t & separator, const int & divHeight) const;
	
	template<typename T>
	void computeMortonSepsG(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
		const int & nbSeps, int * targets,  int * nbUntilNode, int64_t * sepNodes, 
		const int64_t & rootNodeID, const int & divHeight) const;		
};


template<typename T>
void LBMortonBase::computeMortonOneSep(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbChildren, 
	const int & target, int & nbUntilNode, int64_t & separator, const int & divHeight) const
{	
	// restart
	int sum = 0;
	for (int i=0; i<nbChildren; i++)
		sum += globalBuffer[i];
	nbUntilNode -= sum;
	
	// for each leave
	for (int i=0; i<nbChildren; i++)
	{
		// increase the counter
		nbUntilNode += globalBuffer[i];	
		
		// if the counter exceeds the max		
		if (nbUntilNode >= target )
		{	

			// update localSep, counters and exit
			separator = IDs[i];
			break;
		}
	}
}

template<typename T>
void LBMortonBase::computeMortonSeps(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbLeaves, 
	const int & nbSeps, int * targets,  int * nbUntilNode, int64_t * sepNodes, 
	const int64_t & rootNodeID, const int & divHeight) const
{
	int * tabDone = new int[nbSeps]();				// Flag indicating the separator's state		

	for (int i=0; i<nbLeaves; i++) 					// for each box
		for (int k=0; k<nbSeps; k++)				// For each separator
			if(!tabDone[k])							// If not already treated
			{
				nbUntilNode[k] += globalBuffer[i];
				
				// If the number of particles exceeds the max
				if (nbUntilNode[k] >= targets[k])
				{	
					// Update the separators tab
					sepNodes[k] = IDs[i];

					// and mark the separator as treated
					tabDone[k]=1;
				}
			}
	// Dealloc
	delete [] tabDone;
}



	/*** Gaspi versions **/

// Keep this version for the Gaspi asynchronous Morton version
// TODO : to be updated
template<typename T>
void LBMortonBase::computeMortonOneSepG(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
	const int & target, int & nbUntilNode, int64_t & separator, const int & divHeight) const
{	
	// restart
	int sum = 0;
	for (int i=0; i<nbLeaves; i++)
		sum += globalBuffer[i];
	nbUntilNode -= sum;
		
	// for each leave
	for (int i=0; i<nbLeaves; i++)
	{
		// increase the counter
		nbUntilNode += globalBuffer[i];			
		
		// if the counter exceeds the max		
		if (nbUntilNode >= target )
		{	
			// update localSep, counters and exit
			separator = n->computeId(divHeight, separator, i);
			break;
		}
	}
}

// Keep this version for the Gaspi asynchronous Morton version
// TODO : to be updated
template<typename T>
void LBMortonBase::computeMortonSepsG(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
	const int & nbSeps, int * targets,  int * nbUntilNode, int64_t * sepNodes, 
	const int64_t & rootNodeID, const int & divHeight) const
{	
	int * tabDone = new int[nbSeps]();				// Flag indicating the separator's state		

	for (int i=0; i<nbLeaves; i++) 					// for each box
		for (int k=0; k<nbSeps; k++)				// For each separator
			if(!tabDone[k])							// If not already treated
			{
				nbUntilNode[k] += globalBuffer[i];
				
				// If the number of particles exceeds the max
				if (nbUntilNode[k] >= targets[k])
				{	
					// Update the separators tab
					sepNodes[k] = n->computeId(divHeight, rootNodeID, i);

					// and mark the separator as treated
					tabDone[k]=1;
				}
			}
	// Dealloc
	delete [] tabDone;
}

#endif
