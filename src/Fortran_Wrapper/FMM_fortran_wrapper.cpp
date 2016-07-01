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

#include <string>
#include <iostream>

#include "mpi.h"
#include "vec3D.hpp"
#include "Particles.hpp"
#include "Node.hpp"
#include "Decomposition.hpp"
#include "LoadBalancerBase.hpp"
#include "LoadBalancer.hpp"
#include "LBMortonSyncMPI.hpp"
#include "LBHistApprox.hpp"
#include "FMM_fortran_wrapper.hpp"
using namespace std;

// Load Balancing strategy
#define	MORTON_MPI_SYNC 1
#define HIST_APPROX 2
#define MORTON_GASPI_ASYNC 3

// Global variables
static vec3D * elements;


/**
* This function updates the elements array by computing the gravity center of each element.
* @param elemToNode : For each element, number of nodes f indexes of each node
* @param nbElem : Total number of elements
* @param nodesXcoords : X node coordinates list
* @param nodesYcoords : Y node coordinates list
* @param nodesZcoords : Z node coordinates list
*/
void fmm_get_elem_coors_ (int * elemToNode, i64 * nbElem, double * nodesXcoords, double * nodesYcoords, double * nodesZcoords)
{
	// allocate global elements array
	elements = new vec3D[*nbElem];
	
	// for each element
	int index = 0;
	int elemCpt = 0;
	while(elemCpt<*nbElem)
	{
		// get number of current element's nodes
		int nbNodes = elemToNode[index];
		index++;
		
		// compute the element's Gravity Center
		double Gx = 0;
		double Gy = 0;
		double Gz = 0;
		
		for (int i=0; i<nbNodes; i++)
		{
			int nodeIdx = elemToNode[index] - 1; // from F to C
			Gx += nodesXcoords[nodeIdx];
			Gy += nodesYcoords[nodeIdx];
			Gz += nodesZcoords[nodeIdx];
			index++;
		}
		
		Gx /= (nbNodes)*1.0;
		Gy /= (nbNodes)*1.0;
		Gz /= (nbNodes)*1.0;
		
		// update elements array
		elements[elemCpt] = vec3D(Gx, Gy, Gz);
		elemCpt++;
	}	
}	

/**
* This function 
* @param * nbElemPerNode : Array of number of elements, per octree node
* @param * firstElemAdress : /!\ adresse du 1er element de chaque cellule dans le tableau tab
* @param * nbSonsPerNode : Array of number of sons, per octree node
* @param * firstSonId : Array of number of first son Id, per octree node
* @param * nodeOwners : Array of responsible mpi ranks, per octree node. (From 1 to nbRanks, Fortran's style)
* @param * nodeCenters : Array of octree node centers
* @param * endlev : Array of last octree node id, per octree level
* @param * nbLevels : Number of Octree levels
* @param * maxEdge : Edge of the first octree node. (the biggest)
* @param *LBstrategy : integer used to choose the load balancing strategy
*/
void fmm_load_balance_(
	i64 * nbElemPerNode, 
	i64 * firstElemAdress, 
	i64 * nbSonsPerNode, 
	i64 * firstSonId, 
	i64 * nodeOwners, 
	double * nodeCenters, 
	i64 * endlev, 
	i64 * nbLevels,
	double * maxEdge, 
	int * LBstrategy)
{
	// création de l'arbre et lecture des particules
	Node<Particles> * treeHead = nullptr;
	Particles p;
	p.setAttributes(0,nbElemPerNode[0],*maxEdge, vec3D(nodeCenters[0] << " " << nodeCenters[1] << " " << nodeCenters[2]));
	
	treeHead = new Node<Particles>(p);
	int nbElements = nbElemPerNode[0];
	p.setNewCoordinates(elements, nbElements);
	p.scale();
	
	int height = (*nbLevels) - 1;
	int nbLeaves = endlev[height] - endlev[height-1];
	int firstLeave = endlev[height-1]; // endlev[height-1] +1(next) -1(F to C) 
	
	// Apply Morton Load Balancing
	int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	decompo nb1ers(size);
	int first = 0;
	int last = nbElemPerNode[0]-1;
	vec3D center(*Xc, *Yc, *Zc);

	// scale the center too
	scale(center);
	
	int nbNodes2 = endlev[height]*3;
	double * centers = new double[nbNodes2];
	copyAndScaleArray(nodeCenters, centers, nbNodes2);
	
	//scaleArray(nodeCenters, endlev[height]*3);

	// Load Balance
	LB_Base * LBB = nullptr;
	cout << *LBstrategy << endl;
	switch (*LBstrategy)
	{
		case MORTON_MPI_SYNC : 
			cout << "Load Balancing with Morton Space Filling Curve" << endl;
			treeHead -> read_octree(nbElemPerNode, firstElemAdress, nbSonsPerNode, firstSonId, nodeCenters);
			LBB = new LoadBalancer<Particles, MortonSyncMPI>(treeHead, nb1ers, 0, 0, first, last, *maxEdge, center, nullptr, nullptr, nodeOwners, 0);
			break;
		case HIST_APPROX :
			cout << "Load Balancing with Histograms, maxEdge = " << *maxEdge << endl;
			LBB = new LoadBalancer<Particles, HistApprox>(treeHead, nb1ers, 0, 0, first, last, *maxEdge, center, nullptr, &centers[firstLeave*3], &nodeOwners[firstLeave], nbLeaves);
			break;
		default :
			cerr << "No identified Load Balancing strategy" << endl;
			exit(0);
	}
	LBB->run();
	
	delete Particles::_coordinates;
}
