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
static Gaspi_m2l_communicator * gCommM2L = nullptr;
static Gaspi_unknowns_communicator * gCommUNK = nullptr;


void fmm_handle_unknowns_broadcast_(complex * xtmp, complex * xtmp2, i64 * size)
{
	// Passage en Gaspi
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	
	// First call : Class instantiation, allocations and Gaspi segment creation
	if (! gCommUNK)
	{	
		gCommUNK = new Gaspi_unknowns_communicator(xtmp, xtmp2, (int) (*size));
	}
	
	// run broadcast
	gCommUNK->runBroadcastUnknowns();
	t_end = MPI_Wtime();
	add_time_sec("GASPI_broadcast", t_end - t_begin);
}

void init_gaspi_ff_communicator_(
							i64 * recvnode,	i64 * recvnode_sz,
							i64 * sendnode,	i64 * sendnode_sz,
							i64 * nb_recv, 	i64 * nb_recv_sz,
							i64 * nb_send, 	i64 * nb_send_sz,
							i64 * nivterm,	i64 * levcom,
							i64 * fniv,		i64 * nst,			i64 * nsp,
							complex * ff, 	complex * ne, 		i64 * allreduce_sz,
							i64 * fsend, 	i64 * send,
							i64 * frecv,	i64 * recv,
							i64 * endlev,	i64 * codech, 
							i64 * ff_sz)
{
	// switch to Gaspi
	fmm_switch_to_gaspi_();
	
	// Construction du communicateur Gaspi M2L, si nécessaire
	int indexToC = -1;
	int l = (int)(*levcom);
	
	if (! gCommM2L)
	{
    	construct_m2l_communicator(
			nb_send, (int)(*nb_send_sz), 
			nb_recv, (int)(*nb_recv_sz),
			sendnode, (int)(*sendnode_sz),
			recvnode, (int)(*recvnode_sz),
			(int)(*nivterm), (int)(*levcom),
			/*&ff[fniv[l]],
			ne, (int)(*allreduce_sz),*/
			fsend, send, frecv, recv, nst, nsp, fniv, endlev, codech, /*(int)(*ff_sz),*/
			gCommM2L);	
	}
	// switch back to mpi
	fmm_switch_to_mpi_();
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
	gCommM2L->runM2LCommunications(bufsave, ff);
}

void gaspi_send_ff_(i64 * niv, complex * ff)
{
	gCommM2L->send_ff_level((int)(*niv)-1, ff);
}

void gaspi_recv_ff_(i64 * niv, complex * ff)
{
	gCommM2L->recv_ff_level((int)(*niv)-1, ff);
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
	delete gCommM2L;
	MPI_Barrier(MPI_COMM_WORLD);
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));	
	SUCCESS_OR_DIE(gaspi_proc_term(GASPI_BLOCK));
	MPI_Barrier(MPI_COMM_WORLD);
}

// debug tools
void fmm_dump_(complex * tab)
{
	int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	dumpBuffer(mpi_rank, tab, 10, "fortran", "ne, ff");
}
