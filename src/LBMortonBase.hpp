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
	void computeMortonSeps(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
		const int & nbSeps, int * targets,  int * nbUntilNode, int64_t * sepNodes, 
		const int64_t & rootNodeID, const int & divHeight) const;

	template<typename T>
	void computeMortonSeps2(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbLeaves, 
		const int & nbSeps, int * targets,  int * nbUntilNode, int64_t * sepNodes, 
		const int64_t & rootNodeID, const int & divHeight) const;

	template<typename T>		
	void computeMortonOneSep(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
		const int & target, int & nbUntilNode, int64_t & separator, const int & divHeight) const;
	
	template<typename T>		
	void computeMortonOneSep2(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbLeaves, 
		const int & target, int & nbUntilNode, int64_t & separator, const int & divHeight) const;

};



template<typename T>
void LBMortonBase::computeMortonSeps(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
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

template<typename T>
void LBMortonBase::computeMortonSeps2(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbLeaves, 
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


template<typename T>
void LBMortonBase::computeMortonOneSep(Node<T> * n, const int * globalBuffer, const int & nbLeaves, 
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

template<typename T>
void LBMortonBase::computeMortonOneSep2(Node<T> * n, const int * globalBuffer, i64 * IDs, const int & nbChildren, 
	const int & target, int & nbUntilNode, int64_t & separator, const int & divHeight) const
{	
//	printf("---------> ENTERING FromComputeMortonOneSep2, at NODE : %ld, with nbUntilNode = %d, nbChildren = %d\n", n->getId(), nbUntilNode, nbChildren);

	// restart
	int sum = 0;
	for (int i=0; i<nbChildren; i++)
	{
//		printf("---------> BUFFER, child n° %d contains : %d\n", i, globalBuffer[i]);
		sum += globalBuffer[i];
	}
	nbUntilNode -= sum;
	
//	printf("---------> FromComputeMortonOneSep2, sum = %d, nbUntilNode = %d\n", sum, nbUntilNode);
	
	
	// for each leave
	for (int i=0; i<nbChildren; i++)
	{
		// increase the counter
		nbUntilNode += globalBuffer[i];	
//		printf("---------> Just added %d, nbUntilNode = %d\n", globalBuffer[i], nbUntilNode);
		
		// if the counter exceeds the max		
		if (nbUntilNode >= target )
		{	
//			printf("---------> NOW EXCEEDS, new sep = %ld\n", IDs[i]);

			// update localSep, counters and exit
			//separator = n->computeId(divHeight, separator, i);
			separator = IDs[i];
			break;
		}
	}
}




#endif
