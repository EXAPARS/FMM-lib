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

#include "FMM.hpp"
#include "mpi.h"


#include <iostream>
using namespace std;


#include "Node.hpp"
#include "Particles.hpp"
#include "Gaspi_communicator.hpp"
#include "Decomposition.hpp"
#include "LoadBalancerBase.hpp"
#include "LoadBalancer.hpp"
#include "LBHistApprox.hpp"
#include "LBHistExact.hpp"
#include "LBMortonSyncMPI.hpp"
#include "LBMortonAsyncGASPI.hpp"

void FMM_load_balance(string file, int nbParticles, double dist, double tolerance, int LBStrategy)
{
	int size; MPI_Comm_size(MPI_COMM_WORLD,&size);
	decompo nb1ers(size);
	int first = 0;
	int last = nbParticles-1;
	
	Node<Particles> * treeHead = nullptr;
	Gaspi_communicator * gComm = nullptr;

	// tree creation
	Particles p;
	if ( LBStrategy == HIST_EXACT ||
		 LBStrategy == HIST_APPROX ||
		 LBStrategy == MORTON_MPI_SYNC)
	{
		p.loadCoordinatesASCII(nbParticles, file);		
	}
	else
	{	
		gComm = new Gaspi_communicator(512, nbParticles);		
		p.loadCoordinatesASCII(nbParticles, file, gComm);		
	}
	
	p.scale();
	treeHead = new Node<Particles>(p);
	
	// Node Owners array
	i64 * nodeOwners = nullptr;
	// Load Balancing
	LB_Base * LBB = nullptr;
	switch (LBStrategy)
	{
		case HIST_EXACT : 
			cout <<"--> Exact Histograms" << endl;
			LBB = new LoadBalancer<Particles, HistExact> (treeHead, nb1ers, dist, tolerance, first, last, gComm, nodeOwners);
			LBB->run();		
			break;
		
		case HIST_APPROX :		
			cout <<"--> Approx Histograms" << endl;			
			LBB = new LoadBalancer<Particles, HistApprox> (treeHead, nb1ers, dist, tolerance, first, last, gComm, nodeOwners);
			LBB->run();			
			break;
		
		case MORTON_MPI_SYNC :
			cout <<"--> Morton DFS, with tolerance : " << tolerance << endl;			
			LBB = new LoadBalancer<Particles, MortonSyncMPI>(treeHead, nb1ers, dist, tolerance, first, last, gComm, nodeOwners);
			LBB->run();
			break;
		
		case MORTON_GASPI_ASYNC :			
			cout <<"--> Morton DFS Async, with tolerance : " << tolerance << endl;
			LBB = new LoadBalancer<Particles, MortonAsyncGASPI>(treeHead, nb1ers, dist, tolerance, first, last, gComm, nodeOwners);
			LBB->run();			
			break;
		
		default :
			cerr << "No load balancing strategy detected !";
			exit(5);
	}	
}
