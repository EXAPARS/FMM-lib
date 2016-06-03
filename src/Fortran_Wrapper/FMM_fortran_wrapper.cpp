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
#include "Particles.hpp"
#include "Node.hpp"
#include "Decomposition.hpp"
#include "LoadBalancerBase.hpp"
#include "LoadBalancer.hpp"
#include "LBMortonSyncMPI.hpp"
#include "FMM_fortran_wrapper.hpp"
using namespace std;

// Load Balancing strategy
#define HIST_EXACT 0
#define HIST_APPROX 1
#define	MORTON_MPI_SYNC 2
#define MORTON_GASPI_ASYNC 3

void fmm_load_balance_(i64 * nbElemPerNode, i64 * firstElemAdress, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nbNodes, i64 * nodeOwners, double * centers, i64 * endlev, i64 * nbLevels)
{
	// Dump octree info with a recursive DFS traversal


	// Read octree from fortran application into lib
	Node<Particles> * treeHead = nullptr;
	Particles p;
	p.setAttributes(0,nbElemPerNode[0],0, vec3D(0,0,0));
	treeHead = new Node<Particles>(p);
	treeHead -> read_octree(nbElemPerNode, firstElemAdress, nbSonsPerNode, firstSonId, centers);
	
	// Apply Morton Load Balancing
	int size; MPI_Comm_size(MPI_COMM_WORLD,&size);
	decompo nb1ers(size);
	int first = 0;
	int last = nbElemPerNode[0]-1;
	
	LB_Base * LBB = nullptr;
	LBB = new LoadBalancer<Particles, MortonSyncMPI>(treeHead, nb1ers, 0, 0, first, last, nullptr, nodeOwners);
	LBB->run();
	
	
}

// Possible dumps not used anymore made by Root only
	//dfs_dump_spectre_octree("dump/lib/after_morton_graphViz_", nbElemPerNode, nbSonsPerNode, firstSonId, nbNodes, nodeOwners, 0, centers);
	//dfs_dump_centers("dump/lib/centers_", nbSonsPerNode, firstSonId, 0, centers);
	//bfs_dump_centers_level_by_level("dump/lib/centers_", nbElemPerNode, nbSonsPerNode, firstSonId, nbNodes, nodeOwners, 0, centers, endlev, nbLevels);
	/** dfs_dump_spectre_octree("before_lib_", nbElemPerNode, nbSonsPerNode, firstSonId, nbNodes, nodeOwners, 0);**/
