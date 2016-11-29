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

#ifndef GASPI_M2L_COMMUNICATOR
#define GASPI_M2L_COMMUNICATOR

#include "GASPI.h"
#include "mpi.h"
#include "types.hpp"
#include "Complex.hpp"
#include "fmm_tools.hpp"
#include <cilk/cilk.h>


#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/measure.hpp"
#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/byteCounter.hpp"



#include <iostream>
using namespace std;
/* ---- RAPPELS DES SEGMENTS UTILISES
 *  global recv buffer
 *  global send buffer
 *  remote recv buffer adresses : where to write
 *  data to send adresses : what to write
 */



class Gaspi_m2l_communicator
{	
private:
	// segments sizes (in bytes for segment creation)
	gaspi_size_t _seg_globalRecvBuffer_size; 
	gaspi_size_t _seg_globalSendBuffer_size; 
    gaspi_size_t _seg_globalRecvBufIdxPerRank_size;    
    gaspi_size_t _seg_remoteBufferIndexes_size;
    gaspi_size_t _seg_reduce_size;

	// gaspi pointers on segments
	gaspi_pointer_t _ptr_seg_globalRecvBuffer = nullptr;
    gaspi_pointer_t _ptr_seg_globalSendBuffer = nullptr;
	gaspi_pointer_t _ptr_seg_globalRecvBufIdxPerRank = nullptr;
    gaspi_pointer_t _ptr_seg_remoteBufferIndexes = nullptr;
	gaspi_pointer_t _ptr_seg_ff_allreduce = nullptr;
	gaspi_pointer_t _ptr_seg_ne_allreduce = nullptr;

public:
	// processes info
	gaspi_rank_t _wsize;
	gaspi_rank_t _rank;

	// segments Ids
    gaspi_segment_id_t _seg_globalRecvBuffer_id; //	7 global Receive Buffer 
    gaspi_segment_id_t _seg_globalSendBuffer_id; //	8
    gaspi_segment_id_t _seg_globalRecvBufIdxPerRank_id; // 9
    gaspi_segment_id_t _seg_remoteBufferIndexes_id; // 10
    gaspi_segment_id_t _seg_ff_allreduce_id;
	gaspi_segment_id_t _seg_ne_allreduce_id;

	// usual pointers on values in segments
	complex * _globalRecvBuffer = nullptr;
	complex * _globalSendBuffer = nullptr;
    int * _globalRecvBufIdxPerRank = nullptr;
    int * _remoteBufferIndexes = nullptr;
    complex * _reduceNE = nullptr;
    complex * _reduceFF = nullptr;
    
    // other arrays
    int * _sendBufferIndexes = nullptr;

public:
	Gaspi_m2l_communicator(
		i64 * nb_send, int nb_send_sz, 
		i64 * nb_recv, int nb_recv_sz,
		i64 * sendnode, int sendnode_sz,
		i64 * recvnode, int recvnode_sz,
		complex * ff, complex * ne, int allreduce_sz);

	void create_allReduceBuffers(complex * ff, complex * ne, int nbEltsToReduce);
	void create_globalRecvBuffer(i64 * nb_recv, int nb_recv_sz);
	void create_remoteBufferIndexes(i64 * recvnode, int recvnode_sz, i64 * nb_recv);
	void create_globalSendBuffer(i64 * nb_send, int nb_send_sz);
	void init_sendBufferIndexes(i64 * sendnode, int sendnode_sz, i64 * nb_send);	


	void runM2LallReduce(complex * ff, complex * ne);
	void runM2LCommunications (i64 * sendnode, int sendnode_sz, i64 * nb_send,int levcom, int nivterm, 
		i64 * endlev, i64 * frecv, i64 * recv, i64 * fsend, i64 * send, i64 * nst, i64 * nsp, i64 * fniv, i64 * codech, complex * bufsave, complex * ff);

	void initGlobalSendSegment(i64 * sendnode, int sendnode_sz, i64 * nb_send, int nivterm, int levcom, i64 * fsend, i64 * send, i64 * endlev,
		i64 * codech, i64 * nst, i64 * nsp, complex * bufsave, i64 * fniv, complex * ff);
	void initAllReduceBuffers(complex * ff, complex * ne);

	
	void updateFarFields(int src, int levcom, int nivterm, i64 * endlev, i64 * frecv, i64 * recv, i64 * nst, i64 * nsp, i64 * fniv, complex * ff);
};



void construct_m2l_communicator(i64 * nb_send, int nb_send_sz, 
							 i64 * nb_recv, int nb_recv_sz,
							 i64 * sendnode, int sendnode_sz,
							 i64 * recvnode, int recvnode_sz,
							 complex * ff, complex * ne, int size,
                             Gaspi_m2l_communicator *& gCommM2L);

/** GASPI TOOLS **/
void print_gaspi_config();
//void print_gaspi_config(char ** out);

/**
 * GASPI TAGS - Notification values
 **/

// seg remoteIndexes
#define REMOTE_ADDRESS 1


// seg global recv bugger
#define SEND_DATA 10
#define NO_DATA 20
#define ALLREDUCE 30

#endif
