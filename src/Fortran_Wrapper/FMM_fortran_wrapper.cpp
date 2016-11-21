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


#include "FMM_fortran_wrapper.hpp"


using namespace std;

// Load Balancing strategy
#define	MORTON_MPI_SYNC 1
#define HIST_APPROX 2
#define HIST_APPROX_X 3
#define MORTON_GASPI_ASYNC 4 

// Global variables
static vec3D * elements = nullptr;
static Gaspi_m2l_communicator * gCommM2L = nullptr;


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
* @param * firstElemAdress : Array of Indexes. 
* 							 Gives the index(address) first element's data in the Id's array fmmData%elem, per octree node.
* @param * nbSonsPerNode : Array of number of sons, per octree node
* @param * firstSonId : Array of number of first son Id, per octree node
* @param * nodeOwners : Array of responsible mpi ranks, per octree node. (From 1 to nbRanks, Fortran's style)
* 						The Result of the load balancing is read here.
* @param * nodeCenters : Array of octree node centers
* @param * endlev : Array of last octree node id, per octree level
* @param * nbLevels : Number of Octree levels
* @param * maxEdge : Edge of the first octree node. (the biggest)
* @param *LBstrategy : integer used to choose the load balancing strategy
*/
void fmm_load_balance_(	i64 * nbElemPerNode, i64 * firstElemAdress, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nodeOwners, double * nodeCenters, i64 * endlev, 
	i64 * nbLevels, double * maxEdge, int * LBstrategy)
{

	// MPI parameters
	int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// create Particles
	vec3D center(nodeCenters[0], nodeCenters[1], nodeCenters[2]);
	scale(center);
	Particles p; 
	p.setAttributes(0,nbElemPerNode[0],*maxEdge, center); 
	int nbElements = nbElemPerNode[0];
	p.setNewCoordinates(elements, nbElements);
	p.scale();
	
	// creatre Octree
	Node<Particles> * octree = nullptr;
	octree = new Node<Particles>(p);
	octree -> read_octree(nbElemPerNode, firstElemAdress, nbSonsPerNode, firstSonId, nodeCenters);
	
	// Octree useful characteristics
	int height = (*nbLevels) - 1;
	int nbLeaves = endlev[height] - endlev[height-1];
	int firstLeave = endlev[height-1]; // endlev[height-1] +1(next node) -1(F to C) 
	int firstElem = 0;
	int lastElem = nbElemPerNode[0]-1;
	int nbNodes = endlev[height]*3;
	double * centers = new double[nbNodes];
	
	// Load Balancing useful parameters
	decompo nb1ers(size);
	
	// Load Balance
	LB_Base * LBB = nullptr;
	switch (*LBstrategy)
	{
		case MORTON_MPI_SYNC :
			cout << "*** ----- MORTON -----" << endl;
			LBB = new LoadBalancer<Particles, MortonSyncMPI>(octree, nb1ers, 0, 0, firstElem, lastElem, *maxEdge, center, nullptr, nullptr/*centers*/, nodeOwners, 0);
			break;
		case HIST_APPROX :
			cout << "*** ----- KD TREE -----" << endl;		
			// Get a copy of the octree centers, scale them and apply load balancing strategy
			copyAndScaleArray(nodeCenters, centers, nbNodes);
			LBB = new LoadBalancer<Particles, HistApprox>(octree, nb1ers, 0, 0, firstElem, lastElem, *maxEdge, center, nullptr, &centers[firstLeave*3], &nodeOwners[firstLeave], nbLeaves);
			break;
		case HIST_APPROX_X :
			cout << "*** ----- SAUCISONNAGE EN X -----" << endl;		
			// Get a copy of the octree centers, scale them and apply load balancing strategy
			copyAndScaleArray(nodeCenters, centers, nbNodes);
			// test rapide - saucissonnage en X
			cout << "Decomposition has been modified to : " << endl;
			nb1ers._list.resize(1);
			nb1ers._list[0] = size;
			nb1ers.display();			
			LBB = new LoadBalancer<Particles, HistApprox>(octree, nb1ers, 0, 0, firstElem, lastElem, *maxEdge, center, nullptr, &centers[firstLeave*3], &nodeOwners[firstLeave], nbLeaves);
			break;
		default :
			cerr << "No identified Load Balancing strategy" << endl;
			exit(0);
	}
	LBB->run();

	delete Particles::_coordinates;
	cout << "---- Load Balancing ----> TERMINATED !" <<endl;
	int * counters = new int[size + 1]();
}

void fmm_handle_comms_gaspi_(i64 * recvnode, 	i64 * recvnode_sz, 
							 i64 * sendnode, 	i64 * sendnode_sz,
							 i64 * nb_recv, 	i64 * nb_recv_sz,
							 i64 * nb_send, 	i64 * nb_send_sz,
							 i64 * nivterm, 
							 i64 * levcom, 
							 i64 * fniv, 
							 i64 * nst, 
							 i64 * nsp,
							 complex * ff,
							 i64 * fsend, 		i64 * send,
							 i64 * frecv,		i64 * recv,
							 i64 * endlev,		i64 * codech,
							 complex * bufsave)
{
	// Passage en Gaspi
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	int mpi_rank, mpi_wsize;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_wsize);
    MPI_Barrier(MPI_COMM_WORLD);
    gaspi_rank_t rank, wsize;
    SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
    SUCCESS_OR_DIE(gaspi_proc_num(&wsize));
    assert(rank == mpi_rank);
    assert(wsize == mpi_wsize);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_switch_interop", t_end - t_begin);

	// First call : Class instantiation, allocations and Gaspi segment creation
	if (! gCommM2L)
	{	
		t_begin = MPI_Wtime();
    	construct_m2l_communicator(nb_send, (int)(*nb_send_sz), 
							  nb_recv, (int)(*nb_recv_sz),
							  sendnode, (int)(*sendnode_sz),
							  recvnode, (int)(*recvnode_sz),
							  gCommM2L);
		t_end = MPI_Wtime();
		add_time_sec("GASPI_init_class_and_create_segements", t_end - t_begin);
	}

	// All other calls :
	t_begin = MPI_Wtime();
	gCommM2L->runM2LCommunications(sendnode, (int)(*sendnode_sz),
								nb_send, (int) (*levcom), (int) (*nivterm), endlev,
								frecv, recv, fsend, send,  nst, nsp, fniv, codech, bufsave, ff);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_func_echanges", t_end - t_begin);

	// Rend la main au MPI
	t_begin = MPI_Wtime();	
    SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	t_end = MPI_Wtime();
	add_time_sec("GASPI_switch_interop", t_end - t_begin);    
}

void fmm_gaspi_init_()
{
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD); 
	SUCCESS_OR_DIE (gaspi_proc_init (GASPI_BLOCK));
	SUCCESS_OR_DIE (gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	t_end = MPI_Wtime();
	add_time_sec("GASPI_proc_init", t_end - t_begin);
}

void fmm_gaspi_finalize_()
{
	// TODO - dealloc les segments dans le destructeur, appelé seulement à la fin du programme
	delete gCommM2L;

	double t_begin, t_end;
	t_begin = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));	
	SUCCESS_OR_DIE(gaspi_proc_term(GASPI_BLOCK));
	MPI_Barrier(MPI_COMM_WORLD);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_proc_term", t_end - t_begin);	
}
