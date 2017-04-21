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


#include "FMM_Gaspi_wrapper.hpp"

using namespace std;


// Global variables
static Gaspi_FF_communicator * gCommFF = nullptr;

// start stop switch 
void fmm_gaspi_init_()
{
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD); 
	SUCCESS_OR_DIE (gaspi_proc_init (GASPI_BLOCK));
	SUCCESS_OR_DIE (gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	t_end = MPI_Wtime();
	add_time_sec("GASPI_proc_init", t_end - t_begin);
	add_time_sec("GASPI_FF_sendrecv", t_end - t_begin);	
}

void fmm_gaspi_finalize_()
{
	// TODO - dealloc les segments dans le destructeur, appelé seulement à la fin du programme
	delete gCommFF;
	MPI_Barrier(MPI_COMM_WORLD);
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));	
	SUCCESS_OR_DIE(gaspi_proc_term(GASPI_BLOCK));
	MPI_Barrier(MPI_COMM_WORLD);
}

void fmm_switch_to_mpi_()
{
	// Rend la main au MPI
	double t_begin, t_end;
	t_begin = MPI_Wtime();
    SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	t_end = MPI_Wtime();
	add_time_sec("GASPI_switch_interop", t_end - t_begin);
	add_time_sec("GASPI_FF_sendrecv", t_end - t_begin);	 
}

void fmm_switch_to_gaspi_()
{
	// Passe la main au Gaspi
	double t_begin, t_end;
	t_begin = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_switch_interop", t_end - t_begin);
	add_time_sec("GASPI_FF_sendrecv", t_end - t_begin);
}

// Init
/*
void init_gaspi_ff_communicator_(i64 * recvnode,	i64 * recvnode_sz,
	i64 * sendnode, i64 * sendnode_sz, i64 * nb_recv, i64 * nb_recv_sz,
	i64 * nb_send, i64 * nb_send_sz, i64 * nivterm, i64 * levcom,
	i64 * fniv, i64 * nst, i64 * nsp, i64 * fsend, i64 * send,
	i64 * frecv, i64 * recv, i64 * endlev,	i64 * codech, i64 * includeLevcom)
{
	// switch to Gaspi
	fmm_switch_to_gaspi_();
	
	// Construction du communicateur Gaspi M2L, lors du 1er appel
	if (! gCommFF)
	{
    	construct_m2l_communicator(
			nb_send, (int)(*nb_send_sz), 
			nb_recv, (int)(*nb_recv_sz),
			sendnode, (int)(*sendnode_sz),
			recvnode, (int)(*recvnode_sz),
			(int)(*nivterm), (int)(*levcom),
			fsend, send, frecv, recv, nst, nsp, fniv, endlev, codech, (int)(*includeLevcom),
			gCommFF);
	}
	// switch back to mpi
	fmm_switch_to_mpi_();
}*/

// FF
void fmm_handle_ff_gaspi_bulk_(complex * ff, complex * bufsave, i64 * idom)
{
	//cout << "wrapper, bulk" << endl;
	gCommFF->exchangeFFBulk(bufsave, ff, (int)(*idom)-1);
}

void gaspi_send_ff_(i64 * niv, complex * ff, i64 * idom)
{
	//cout << "toto" << endl;
	//fflush(stdout);	
	//int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	//printf("[%d][ENTER] gaspi_send_ff_ , level : %d, domain : %d\n", rank,((int)(*niv)-1), ((int)(*idom)-1));	
	if(gCommFF)
	{
		gCommFF->send_ff_level((int)(*niv)-1, ff, (int)(*idom)-1);
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
	//printf("[%d][EXIT] gaspi_send_ff_ , level : %d, domain : %d\n", rank, (int)(*niv)-1, (int)(*idom)-1);	
	//printf("[%d] sent level : %d, octree : %d\n", rank, (int)(*niv)-1, (int)(*idom)-1);	
	//fflush(stdout);
}

void gaspi_task_send_ff_(i64 * start, i64 * stop, i64 * niv, complex * ff, i64 * idom)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	/*printf("[%d][ENTER] gaspi_task_send_ff_ , level : %d, domain : %d\n", rank, ((int)(*niv)-1), ((int)(*idom)-1));	
	fflush(stdout);*/
	
	if(gCommFF)
	{
		gCommFF->send_task_ff_level((int)(*niv)-1, ff, (int)(*idom)-1, (int)(*start), (int)(*stop));
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
	
	//~ printf("[%d][EXIT] gaspi_task_send_ff_ , level : %d, domain : %d\n", rank, ((int)(*niv)-1), ((int)(*idom)-1));	
	//~ fflush(stdout);
}

void gaspi_task_recv_ff_(i64 * niv, complex * ff, i64 * idom)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	//printf("[%d][ENTER] gaspi_task_recv_ff_ , level : %d, domain : %d\n", rank, ((int)(*niv)-1), ((int)(*idom)-1));	
	//fflush(stdout);
	
	if(gCommFF)
	{
		gCommFF->recv_task_ff_level((int)(*niv)-1, ff, (int)(*idom)-1);
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
	//printf("[%d][EXIT] gaspi_task_recv_ff_ , level : %d, domain : %d\n", rank, ((int)(*niv)-1), ((int)(*idom)-1));	
	//fflush(stdout);	
}

void gaspi_recv_ff_(i64 * niv, complex * ff, i64 * idom)
{

	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	//printf("[%d][ENTER] gaspi_recv_ff_ , level : %d, octree : %d\n",rank, ((int)(*niv)-1), ((int)(*idom)-1));	
	//fflush(stdout);
	if (gCommFF)
	{
		//cout << "into recv level" << (int)(*niv)-1 << "octree : " << (int)(*idom)-1 << endl;
		gCommFF->recv_ff_level((int)(*niv)-1, ff, (int)(*idom)-1);
	}
	else
	{
		cerr << "[wrapper gaspi_recv_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}

//	printf("[%d][EXIT] gaspi_recv_ff_ , level : %d, domain : %d\n",rank,((int)(*niv)-1), ((int)(*idom)-1));
	//printf("[%d] received level : %d, octree : %d\n",rank,((int)(*niv)-1), ((int)(*idom)-1));
	//fflush(stdout);
}

// multimat version
void gaspi_init_ff_(i64 * max_send_terms, i64 * max_recv_terms, i64 * nbMat, i64 * max_send_nodes, i64 * max_recv_nodes, i64 * incLevcom)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
 
	if (!gCommFF)
	{
		init_gaspi_ff((int)(*max_send_terms), (int)(*max_recv_terms), (int)(*max_send_nodes), (int)(*max_recv_nodes), (int)(*nbMat), (int)(*incLevcom), gCommFF);
	}
}

void gaspi_init_offsets_(i64 * recvnode, i64 * recvnode_sz, i64 * sendnode, i64 * sendnode_sz, i64 * nb_recv, 
	i64 * nb_recv_sz, i64 * nb_send, i64 * nb_send_sz, i64 * idom, i64 * ndom,
	i64 * nivterm, i64 * frecv, i64 * recv, i64 * levcom, i64 * endlev, 
	i64 * fniv, i64 * fsend, i64 * send, i64 * nst, i64 * nsp, 
	i64 * codech)
{
	//cout << "gaspi_init_offsets_ [ENTER]" << endl;
	if(gCommFF)
	{
		gCommFF->init_gaspi_offsets(recvnode, (int)(*recvnode_sz), sendnode, (int)(*sendnode_sz), nb_recv, 
			(int)(*nb_recv_sz), nb_send, (int)(*nb_send_sz), (int)(*idom)-1, (int)(*ndom), 
			(int)(*nivterm), frecv, recv, (int)(*levcom), endlev, 
			fniv, fsend, send, nst, nsp, 
			codech);
	}
	else
	{
		cerr << "[wrapper gaspi_init_offsets]Gaspi M2L Communicator is not initialized !" << endl; exit(-1);
	}

	//cout << "gaspi_init_offsets_ [EXIT]" << endl;
	fflush(stdout);
}

// Unk
/*
void fmm_handle_unknowns_broadcast_(complex * xtmp, complex * xtmp2, i64 * size)
{
	// Passage en Gaspi
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	
	// First call : Class instantiation, allocations and Gaspi segment creation
	if (! gCommUNK)
	{	
		gCommUNK = new Gaspi_UNK_communicator(xtmp, xtmp2, (int) (*size));
	}
	
	// run broadcast
	gCommUNK->runBroadcastUnknowns();
	t_end = MPI_Wtime();
	add_time_sec("GASPI_broadcast", t_end - t_begin);
}

void fmm_handle_unknowns_allreduce_()
{
	
	// run allreduce
	double t_begin, t_end;
	t_begin = MPI_Wtime();	
	gCommUNK->runAllReduceUnknowns();
	t_end = MPI_Wtime();
	add_time_sec("GASPI_allReduce_unk", t_end - t_begin);

	// Rend la main au MPI
	t_begin = MPI_Wtime();
    SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	t_end = MPI_Wtime();
	add_time_sec("GASPI_switch_interop", t_end - t_begin); 
}
*/


// debug tools
void fmm_dump_(complex * tab)
{
	int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	dumpBuffer(mpi_rank, tab, 10, "fortran", "ff");
}
