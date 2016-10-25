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
	// segments sizes
	gaspi_size_t _seg_globalRecvBuffer_size; 
    gaspi_size_t _seg_globalRecvBufIdxPerRank_size;    
    gaspi_size_t _seg_remoteBufferIndexes_size;

	// gaspi pointers on segments
	gaspi_pointer_t _ptr_seg_globalRecvBuffer = nullptr;
    gaspi_pointer_t _ptr_seg_globalRecvBufIdxPerRank = nullptr;
    gaspi_pointer_t _ptr_seg_remoteBufferIndexes = nullptr;

public:
	// processes info
	gaspi_rank_t _wsize;
	gaspi_rank_t _rank;

	// segments Ids
    gaspi_segment_id_t _seg_globalRecvBuffer_id; //	7 global Receive Buffer 
    gaspi_segment_id_t _seg_globalRecvBufIdxPerRank_id; // 8
    gaspi_segment_id_t _seg_remoteBufferIndexes_id; // 9
	
	// usual pointers on values in segments
	complex * _globalRecvBuffer = nullptr;
    int * _globalRecvBufIdxPerRank = nullptr;
    int * _remoteBufferIndexes = nullptr;

public:
	Gaspi_m2l_communicator(i64 * nb_recv, int nb_recv_sz, i64 * recvnode, int recvnode_sz);
	void init_globalRecvBuffer(i64 * nb_recv, int nb_recv_sz);
	void init_remoteBufferIndexes(i64 * recvnode, int recvnode_sz, i64 * nb_recv);

};


void init_gaspi_m2l_segments(i64 * nb_recv, int nb_recv_sz, i64 * recvnode, int recvnode_sz) ;

void create_gaspi_m2l_segments(i64 * nb_send, int nb_send_sz, 
							 i64 * nb_recv, int nb_recv_sz,
							 i64 * sendnode, int sendnode_sz,
							 i64 * recvnode, int recvnode_sz,
							 int levcom, int nivterm,
							 i64 * fsend, i64 * send, i64 * endlev, i64 * codech,
							 i64 * nst, i64 * nsp, complex * bufsave, i64 * fniv,
							 complex * ff, 
                             Gaspi_m2l_communicator *& gCommM2L);




void init_globalSendBuffer();
void init_dataToSendIndexes();








/**
 * GASPI TAGS - Notification values
 **/

// seg remoteIndexes
#define REMOTE_ADDRESS 1


  
/*// seg recvBuffer
#define INIT_LEAVES_BUFFER 1
#define LEAVES_BUFFER_ANSWER 2

// seg sepNodes
#define INIT_SEP_NODES 1
#define UPDATE_SEP_NODES 2
#define LEAVES_BUFFER_REQUEST 3

// seg nbUntil
#define INIT_NB_UNTIL_NODE 1

// seg newCoords
#define COORDS_COMPLETED 1
#define COORDS_TO_BE_CONTINUED 2
#define COORDS_EMPTY 3
*/




#endif
