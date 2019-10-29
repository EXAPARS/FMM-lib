/*
  Copyright 2015 - UVSQ
  Authors list: Nathalie M�ller, Eric Petit

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
#include <algorithm>
#include <list>
#include <numeric>
#include <random>
#include <vector>
#include <iterator>

#include "FMM_Gaspi_wrapper.hpp"

using namespace std;


// Global variables
static Gaspi_FF_communicator * gCommFF = nullptr;
static pthread_mutex_t mutex;

// start stop switch 
void fmm_gaspi_init_()
{
	MPI_Barrier(MPI_COMM_WORLD); 
	SUCCESS_OR_DIE (gaspi_proc_init (GASPI_BLOCK));
	SUCCESS_OR_DIE (gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
	pthread_mutex_init(&mutex, NULL);
}

void fmm_gaspi_finalize_()
{
	// TODO - dealloc les segments dans le destructeur, appel� seulement � la fin du programme
	delete gCommFF;
	MPI_Barrier(MPI_COMM_WORLD);
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));	
	SUCCESS_OR_DIE(gaspi_proc_term(GASPI_BLOCK));
	MPI_Barrier(MPI_COMM_WORLD);
}

void fmm_switch_to_mpi_()
{
	// Rend la main au MPI
    SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK)); 
}

void fmm_switch_to_gaspi_()
{
	// Passe la main au Gaspi
    MPI_Barrier(MPI_COMM_WORLD);
}

// Init
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
	if(gCommFF)
	{
		gCommFF->init_gaspi_offsets(recvnode, (int)(*recvnode_sz), sendnode, (int)(*sendnode_sz), nb_recv, 
			(int)(*nb_recv_sz), nb_send, (int)(*nb_send_sz), (int)(*idom)-1, (int)(*ndom), 
			(int)(*nivterm), frecv, recv, (int)(*levcom), endlev, 
			fniv, fsend, send, nst, nsp, codech);
	}
	else
	{
		cerr << "[wrapper gaspi_init_offsets]Gaspi M2L Communicator is not initialized !" << endl; exit(-1);
	}

	//cout << "gaspi_init_offsets_ [EXIT]" << endl;
	fflush(stdout);
}

/* *************************
 * 	Gaspi BULK
 * *************************/
void fmm_handle_ff_gaspi_bulk_(complex * ff, complex * bufsave, i64 * idom)
{
	//cout << "wrapper, bulk" << endl;
	gCommFF->exchangeFFBulk(bufsave, ff, (int)(*idom)-1);
}

/* *************************
 * GASPI - complete level
 * *************************/
void gaspi_recv_ff_(i64 * niv, complex * ff, i64 * idom)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	if (gCommFF)
	{
		gCommFF->recv_ff_level((int)(*niv)-1, ff, (int)(*idom)-1);

	}
	else
	{
		cerr << "[wrapper gaspi_recv_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
}

void gaspi_send_ff_(i64 * niv, complex * ff, i64 * idom)
{

	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	if(gCommFF)
	{
		//cout << "send ff call" << endl;
		gCommFF->send_ff_level((int)(*niv)-1, ff, (int)(*idom)-1);
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
}

/* *************************
 * GASPI - complete level - multithreaded comm only
 * *************************/

void gaspi_recv_ff_multi_(i64 * niv, complex * ff, i64 * idom)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	if (gCommFF)
	{
		gCommFF->recv_ff_level_multithreaded_2((int)(*niv)-1, ff, (int)(*idom)-1);
		//gCommFF->recv_ff_level((int)(*niv)-1, ff, (int)(*idom)-1);
	}
	else
	{
		cerr << "[wrapper gaspi_recv_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
}

void gaspi_send_ff_multi_(i64 * niv, complex * ff, i64 * idom)
{

	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	if(gCommFF)
	{
		//cout << "send ff call" << endl;
		gCommFF->send_ff_level_multithreaded((int)(*niv)-1, ff, (int)(*idom)-1);
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
}


/* *************************
 * 	Gaspi TASK ---  Complete Level
 * *************************/
void gaspi_task_send_ff_(i64 * start, i64 * stop, i64 * niv, complex * ff, i64 * idom)
{
	//~ double t0,t1;
	//~ t0 = MPI_Wtime();
	//~ cout << "In --> gaspi send\n" << endl;
	if(gCommFF)
	{
		gCommFF->send_task_ff_level((int)(*niv)-1, ff, (int)(*idom)-1, (int)(*start), (int)(*stop));
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
	//~ cout << "Out --> gaspi send\n" << endl;
	//~ t1 = MPI_Wtime();
	//~ add_time_sec("gaspi_send", t1-t0);
}

void gaspi_task_recv_ff_(i64 * niv, complex * ff, i64 * idom)
{
	//~ double t0,t1;
	//~ t0 = MPI_Wtime();
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	/*printf("[%d][ENTER] gaspi_task_recv_ff_ , level : %d, domain : %d\n", rank, ((int)(*niv)-1), ((int)(*idom)-1));	
	fflush(stdout);*/
	//debug("enter gaspi_task_recv_ff_");

	if(gCommFF)
	{
		//debug("wait_for", "enter level : " + itoa((int)(*niv)-1));
		gCommFF->recv_task_ff_level((int)(*niv)-1, ff, (int)(*idom)-1);		
		//gCommFF->recv_task_ff_level_rest((int)(*niv)-1, ff, (int)(*idom)-1);	
		//debug("wait_for", "exit level : " + itoa((int)(*niv)-1));

	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}

}


/* *************************
 * 	Gaspi TASK ---  Chunk
 * *************************/

void gaspi_task_chunk_send_ff_(i64 * start, i64 * stop, i64 * niv, complex * ff, i64 * idom)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//int tid = omp_get_thread_num();
	
	if(gCommFF)
	{
		gCommFF->send_task_ff_dbg2((int)(*niv)-1, ff, (int)(*idom)-1, (int)(*start), (int)(*stop));
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
	//debug("debug", itoa(tid) + " exit send, level " + itoa ((int)(*niv)-1));
}

void gaspi_task_chunk_recv_ff_(i64 * niv, complex * ff, i64 * idom)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//int tid = omp_get_thread_num();

	if(gCommFF)
	{
		printf("%d call recv level %d \n", rank, ((int)(*niv)-1)); fflush(stdout);
		gCommFF->recv_task_ff(((int)(*niv)-1), ff, ((int)(*idom)-1));		
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}

}


/* *************************		 *
 * 	Gaspi TASK ---  CHUNK --- ANY 	 *
 * *************************		 */

/*
 * Cette fonction envoie les FF.
 * Ecrit dans le buffer d'envoi.
 * Envoi lorsque COMM_SIZE est atteint, ou complete (= dernier node � envoyer)
 */
void gaspi_task_chunk_send_(i64 * start, i64 * stop, i64 * niv, complex * ff, i64 * idom)
{

	if(gCommFF)
	{
		///gCommFF->send_task_ff_dbg2((int)(*niv)-1, ff, (int)(*idom)-1, (int)(*start), (int)(*stop)); 
		// ==> test run complet, UAV, 4 mpi ok *** checkpoint A ***

		///gCommFF->send_task_ff_dbg2_last((int)(*niv)-1, ff, (int)(*idom)-1, (int)(*start), (int)(*stop)); 
		// ===> test run complet, UAV, 4 mpi, ok *** checkpoint B ***
		// ajout de la modification : notifID indique la src, Value : qty ==>  *** checkpoint C : ok *** (elimination du bug � iter 456)

		gCommFF->send_task_ff_chunk((int)(*niv)-1, ff, (int)(*idom)-1, (int)(*start), (int)(*stop)); 
		// ==> checkpoint D, OK
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
}

/*
 * Cette fonction d�pile les comms.
 * Elle re�oit n'importe quoi de n'importe qui.
 * Arr�t lorsque le niveau est complet
 */
void gaspi_task_chunk_recv_(i64 * niv, complex * ff, i64 * idom)
{

	if(gCommFF)
	{
		///gCommFF->recv_task_ff(((int)(*niv)-1), ff, ((int)(*idom)-1));		 
		// Version d'origine
		// ==> test run complet, UAV, 4 mpi ok *** checkpoint A ***
		
		/// gCommFF->recv_task_ff_last(((int)(*niv)-1), ff, ((int)(*idom)-1));	 // En cours de dev ok, fonctionne
		// avec decodage
		// ==> test run complet, UAV, 4 mpi ok *** checkpoint B *** 
		
		gCommFF->recv_task_ff_any_but_complete_level(((int)(*niv)-1), ff, ((int)(*idom)-1)); 
		///*** checkpoint C : ok ***	
		// CHECKPOINT D : OK
	}
	else
	{
		cerr << "[wrapper gaspi_send_ff] Gaspi M2L Communicator is not initialized !" << endl;
		exit(-1);
	}
}




/* *************************
 * 	Debug Tools
 * *************************/


void fmm_finalize_dump_vector_()
{
	dumpMSG("send_task_ff");
}

void fmm_dump_cplx_(complex * tab, i64* size, i64 * fileNum)
{
	int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	dumpBuffer(mpi_rank, tab, int(*size), int (*fileNum), "rank", "");
}

void fmm_dump_i8_(i64 * tab, i64* size, i64 * fileNum)
{
	int mpi_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	dumpBuffer(mpi_rank, tab, int(*size), int (*fileNum), "rank", "");
}

void fmm_dump_2ble_(double * tab, i64 * size, i64 * fileNum)
{
	int mpi_rank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	dumpBuffer(mpi_rank, tab, int(*size), int (*fileNum), "rank", "");
}

void fmm_raz_i8_(i64 * tab, i64* size)
{
	
	for (i64 i=0; i< (*size); i++)
		tab[i] = 0;
}

void randomizempi_(i64 * wsize, i64 * procIDs)
{

	vector<i64> v(*wsize);	
	for (int i=0; i< *wsize; i++)
		v[i] = i+1;
	
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	
	// master
	if (rank == 0)
	{
		random_shuffle (v.begin(), v.end());
	}
	// broadcast to all ranks
	MPI_Bcast(v.data(), *wsize, MPI_LONG, 0, MPI_COMM_WORLD);
	
	// quick and dirty copy
	for (int i=0; i< *wsize; i++)
		procIDs[i] = v[i];
	
	// output check
	std::cout << "OUTPUT myvector contains:";
	for (int i=0; i< *wsize; i++)
		cout << procIDs[i] << " ";
	cout << endl;
}

