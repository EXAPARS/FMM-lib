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
#include "Gaspi_tools.hpp"
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
public:
	// segments sizes (in bytes for segment creation)
	gaspi_size_t _seg_globalRecvBuffer_size; 
	gaspi_size_t _seg_globalSendBuffer_size; 
    gaspi_size_t _seg_globalRecvBufIdxPerRank_size;    
    gaspi_size_t _seg_remoteBufferIndexes_size;

	// gaspi pointers on segments
	gaspi_pointer_t _ptr_seg_globalRecvBuffer = nullptr;
    gaspi_pointer_t _ptr_seg_globalSendBuffer = nullptr;
	gaspi_pointer_t _ptr_seg_globalRecvBufIdxPerRank = nullptr;
    gaspi_pointer_t _ptr_seg_remoteBufferIndexes = nullptr;

	// processes info
	gaspi_rank_t _wsize;
	gaspi_rank_t _rank;

	// segments Ids
    gaspi_segment_id_t _seg_globalRecvBuffer_id; //	7 global Receive Buffer 
    gaspi_segment_id_t _seg_globalSendBuffer_id; //	8
    gaspi_segment_id_t _seg_globalRecvBufIdxPerRank_id; // 9
    gaspi_segment_id_t _seg_remoteBufferIndexes_id; // 10

	// usual pointers on values in segments
	complex * _globalRecvBuffer = nullptr;
	complex * _globalSendBuffer = nullptr;
    int * _globalRecvBufIdxPerRank = nullptr;
    int * _remoteBufferIndexes = nullptr;
    
    // other arrays
    int * _sendBufferIndexes = nullptr;
    int * _offsetKeeper = nullptr;
    int ** _nbToRecvPerSrcAndLevel = nullptr;
    
    /* arrays from Fortran*/
    int _nivterm;					// hauteur de l'octree = dernier niveau de l'arbre
    int _levcom;					// niveau le + haut où il y a des comms M2L (initialement, toutes les comms)
    i64 * _fniv = nullptr;			// pour chaque niveau, @ dans FF du dernier terme de ce niveau
    i64 * _fsend = nullptr;			// pour chaque rank, @ dans send de la 1ere cellule à échanger
    i64 * _send = nullptr;			// liste des cellules à echanger
    i64 * _frecv = nullptr;
    i64 * _recv = nullptr;
    i64 * _nst = nullptr;			// nb d'angles en theta	
    i64 * _nsp = nullptr;			// nb d'angles en phi
    i64 * _endlev = nullptr;
	i64 * _codech = nullptr;		// 	
	i64 * _nb_send = nullptr;
	i64 * _nb_recv = nullptr;
	i64 * _sendnode = nullptr;
	i64 * _recvnode = nullptr;
	int _nb_send_sz;
	int _nb_recv_sz;
	int _sendnode_sz;
	int _recvnode_sz;    

public:
	Gaspi_m2l_communicator(
		i64 * nb_send, int nb_send_sz, 
		i64 * nb_recv, int nb_recv_sz,
		i64 * sendnode, int sendnode_sz,
		i64 * recvnode, int recvnode_sz,
		int nivterm, int levcom,
		i64 * fsend, i64 * send, i64 * frecv, i64 * recv, i64 * nst, i64 * nsp, i64 * fniv, i64 * endelv, i64 * codech);

	void create_allReduceBuffers(complex * ff, complex * ne, int nbEltsToReduce);
	void create_globalRecvBuffer(i64 * nb_recv, int nb_recv_sz);
	void create_remoteBufferIndexes(i64 * recvnode, int recvnode_sz, i64 * nb_recv);
	void create_globalSendBuffer(i64 * nb_send, int nb_send_sz);
	void init_sendBufferIndexes(i64 * sendnode, int sendnode_sz, i64 * nb_send);	
	
	//void runM2LallReduce(complex * ff, complex * ne);
	void runM2LCommunications (complex * bufsave, complex * ff);

	void initGlobalSendSegment(complex * bufsave, complex * ff);
	void initAllReduceBuffers(complex * ff, complex * ne);
	void init_nbExchangeArrays();
	void updateFarFields(int src, complex * ff);
	void updateFarFields(int src, int level, complex * ff);

	
	// gaspi overlap
	void send_ff_level(int level, complex * ff);
	void recv_ff_level(int level, complex * ff);
};

void construct_m2l_communicator(
	i64 * nb_send, int nb_send_sz, 
	i64 * nb_recv, int nb_recv_sz,
	i64 * sendnode, int sendnode_sz,
	i64 * recvnode, int recvnode_sz,
	int nivterm, int levcom,
	i64 * fsend, i64 * send, i64 * frecv, i64 * recv, 
	i64 * nst, i64 * nsp, i64 * fniv, i64 * endlev, i64 * codech,
	Gaspi_m2l_communicator *& gCommM2L);

/**
 * GASPI TAGS - Notification values
 **/

// seg remoteIndexes
#define REMOTE_ADDRESS 1

// seg global recv buffer

/* 
values 0 to 31 reserved for the octree level
*/

#define SEND_DATA 100
#define NO_DATA 101
#define ALLREDUCE 102

#endif
