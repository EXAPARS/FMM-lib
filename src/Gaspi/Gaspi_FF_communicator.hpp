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

#ifndef GASPI_FF_COMMUNICATOR
#define GASPI_FF_COMMUNICATOR

#include "GASPI.h"
#include "Gaspi_tools.hpp"
#include "mpi.h"
#include <omp.h>
#include <pthread.h>
#include "../Tools/types.hpp"
#include "../Tools/Complex.hpp"
#include "../Tools/fmm_tools.hpp"
#include <cilk/cilk.h>

#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/measure.hpp"
#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/byteCounter.hpp"


#include <iostream>
using namespace std;

class Gaspi_FF_communicator
{	
public:
	// processes info
	gaspi_rank_t _wsize;
	gaspi_rank_t _rank;

/******************************
 * GASPI SEGMENTS --> BUFFERS *
 ******************************/
	/* FF Buffers */
	gaspi_segment_id_t	_FF_sendBuf_seg_id;  						// FF send Buffer
	gaspi_pointer_t 	_FF_sendBuf_seg_ptr = nullptr;	
	complex * 			_FF_sendBuffer = nullptr;
	gaspi_segment_id_t 	_FF_recvBuf_seg_id;  						// FF receive Buffer
	gaspi_pointer_t 	_FF_recvBuf_seg_ptr = nullptr;	
	complex * 			_FF_recvBuffer = nullptr;
	
	/* INFOS Buffers */
	gaspi_segment_id_t 	_Infos_sendbuf_seg_id;						// Infos send Buffer
    gaspi_pointer_t 	_Infos_sendbuf_seg_ptr = nullptr;
    int * 				_Infos_sendbuffer = nullptr;
    gaspi_segment_id_t 	_Infos_recvbuf_seg_id;						// Infos receives Buffer
	gaspi_pointer_t 	_Infos_recvbuf_seg_ptr = nullptr;
	int * 				_Infos_recvbuffer = nullptr;

/***********************************
 * GASPI SEGMENTS --> SEND OFFSETS *
 ***********************************/
 
	/* FF remote send OFFSETS*/
    gaspi_segment_id_t 	_FF_sendRemoteOffsets_seg_id; 				// Remote Offsets for FF send buffer
    gaspi_pointer_t 	_FF_sendRemoteOffsets_seg_ptr = nullptr;
	int * 				_FF_sendRemoteOffsets = nullptr;			
	
	/* FF local send OFFSETS*/
	int ** 				_FF_sendLocalOffsets = nullptr;				// Local Offsets for FF send buffer
																	

	/* INFOS remote send OFFSETS*/
    gaspi_segment_id_t 	_Infos_sendRemoteOffsets_seg_id; 			// Remote Offsets for Infos Send Buffer
	gaspi_pointer_t 	_Infos_sendRemoteOffsets_seg_ptr = nullptr;   
    int * 				_Infos_sendRemoteOffsets = nullptr;			
    
    /* INFOS local send OFFSETS*/
    int **				_Infos_sendLocalOffsets = nullptr;			// Local Offsets for Infos Send Buffer

    
/************************************
 * GASPI SEGMENTS --> RECV OFFSETS  *
 ************************************/
	
	/* FF recv OFFSETS*/		
	gaspi_segment_id_t	_FF_recvOffsets_seg_id; 					// Offsets for FF recv buffer
	gaspi_pointer_t 	_FF_recvOffsets_seg_ptr = nullptr;
	int * 				_FF_recvOffsets = nullptr;
	
	/* INFOS RECV OFFSETS */    
    gaspi_segment_id_t 	_Infos_recvOffsets_seg_id; 					// Offsets for Infos Recv Buffer
	gaspi_pointer_t 	_Infos_recvOffsets_seg_ptr = nullptr;
    int * 				_Infos_recvOffsets = nullptr;

/************************************
 * COUNTERS  *
 ************************************/

    int ** _FF_sendLocalOffsets_counter = nullptr;  					// Counter to add to _FF_sendLocalOffsets
    int ** _FF_sendRemoteOffsets_counter = nullptr; 					// Counter to add to _FF_sendRemoteOffsets
	
	int ** _Infos_sendLocalOffsets_counter = nullptr; 				// Counter to add to _Infos_sendLocalOffsets
	int ** _Infos_sendRemoteOffsets_counter = nullptr;				// Counter to add to _Infos_sendRemoteOffsets


/*************************
 * GASPI ASYNC COMM DATA *
 *************************/    

    int *** _Expect = nullptr;		// per domain, per rank, per level
    int *** _start_send = nullptr;  // per domain, per rank, per level
    int *** _stop_send  = nullptr;  // per domain, per rank, per level
    int *** _count_send = nullptr;  // per domain, per rank, per level
    int *** _start_recv = nullptr;
    int *** _stop_recv  = nullptr;
    int _nbQueues;
          
    
/******************
 * FORTRAN DATA   *
 ******************/
    i64 ** _fniv = nullptr;			// pour chaque niveau, @ dans FF du dernier terme de ce niveau
    i64 ** _fsend = nullptr;		// pour chaque rank, @ dans send de la 1ere cellule à échanger
    i64 ** _send = nullptr;			// liste des cellules à envoyer,  par octree
    i64 ** _frecv = nullptr;
    i64 ** _recv = nullptr;			// liste des cellules à recevoir, par octree
    i64 ** _nst = nullptr;			// nb d'angles en theta	
    i64 ** _nsp = nullptr;			// nb d'angles en phi
    i64 ** _endlev = nullptr;
	i64 ** _codech = nullptr;
	i64 ** _nb_send = nullptr;		// per octree, dest
	i64 ** _nb_recv = nullptr;
	i64 ** _sendnode = nullptr;
	i64 ** _recvnode = nullptr;
	int * _nb_send_sz;
	int * _nb_recv_sz;
	int * _sendnode_sz;
	int * _recvnode_sz;
    int * _nivterm;					// hauteur de l'octree = dernier niveau de l'arbre
    int * _levcom;					// niveau le + haut où il y a des comms M2L (initialement, toutes les comms)	
	int _incLevcom;
	int _nbOct; 

public:
	Gaspi_FF_communicator(i64 * nb_send, int nb_send_sz, i64 * nb_recv, int nb_recv_sz, i64 * sendnode, int sendnode_sz,
		i64 * recvnode, int recvnode_sz, int nivterm, int levcom, i64 * fsend, i64 * send, i64 * frecv, i64 * recv, 
		i64 * nst, i64 * nsp, i64 * fniv, i64 * endlev, i64 * codech, int includeLevcom);
		
	Gaspi_FF_communicator(int max_send_terms, int max_recv_terms, int max_send_nodes, int max_recv_nodes, int incLevcom, int nbMat);

	void create_allReduceBuffers(complex * ff, complex * ne, int nbEltsToReduce);
	void create_FF_recvBuffer(i64 * nb_recv, int nb_recv_sz);
	void create_remoteBufferIndexes(i64 * recvnode, int recvnode_sz, i64 * nb_recv);
	void create_globalSendBuffer(i64 * nb_send, int nb_send_sz);
	void init_FF_sendBufferIndexes(i64 * sendnode, int sendnode_sz, i64 * nb_send);	
	
	void exchangeFFBulk (complex * bufsave, complex * ff, int iOct);

	void initGlobalSendSegment(complex * bufsave, complex * ff, int iOct);
	void initAllReduceBuffers(complex * ff, complex * ne);
	void init_expectPerSrcAndLevel();
	void updateFarFields(int src, complex * ff, int iOct);
	void updateFarFields(int src, int level, complex * ff);
	void updateFarFields(int src, int level, complex * ff, int iOct);
	void updateFarFieldsFromInfos(int src, int level, complex * ff, int counter, int iOct);
	
	// gaspi overlap
	void send_ff_level(int level, complex * ff, int iOct);
	void recv_ff_level(int level, complex * ff, int iOct);
	void send_task_ff_level(int level, complex * ff, int iOct, int start, int stop);
	void recv_task_ff_level(int level, complex * ff, int iOct);

	// multimat version
	void create_segments(int max_send_terms, int max_recv_terms, int max_send_nodes, int max_recv_nodes);
	void init_gaspi_offsets(i64 * recvnode, int recvnode_sz, i64 * sendnode, int sendnode_sz, 
		i64 * nb_recv, int nb_recv_sz, i64 * nb_send, int nb_send_sz, int mat, int nbMat,
		int nivterm, i64 * frecv, i64 * recv, int levcom, i64 * endlev, i64 * fniv, i64 * fsend, i64 * send, i64 * nst, i64 * nsp, i64 * codech);
				
	void fill_FF_sendRemoteOffsets(i64 * recvnode, int recvnode_sz, i64 * nb_recv, int nb_recv_sz, int iOct);
	void fill_FF_sendLocalOffsets(i64 * sendnode, int sendnode_sz, i64 * nb_recv, int nb_recv_sz, int iOct);
	void fill_expectations(int iOct);
	void fill_send_start_stop_count(int iOct);
	void fill_recv_start_stop(int iOct);
	void fill_send_info_ptrs(int iOct);
	void fill_Infos_recvOffsets(int iOct);
	void fill_Infos_sendRemoteOffsets(int iOct);
	void fill_Info_sendLocalOffsets(int iOct);
	
	void fill_attributes(int iOct, int nivterm, int levcom, i64 * fniv, i64 * fsend, i64 * send, i64 * frecv, i64 * recv, i64 * nst, i64 * nsp,
		i64 * endlev, i64 * codech, i64 * nb_send, i64 * nb_recv, i64 * sendnode, i64 * recvnode, int nb_send_sz, int nb_recv_sz, int sendnode_sz, int recvnode_sz);
	void alloc_attributes();
	
	//debug
	void fmm_raz_gaspi_segments(int max_send_terms, int max_recv_terms, int max_send_nodes, int max_recv_nodes);
	void send_chunk(int iOct, int dest, int level);
};

void construct_m2l_communicator(
	i64 * nb_send, int nb_send_sz, 
	i64 * nb_recv, int nb_recv_sz,
	i64 * sendnode, int sendnode_sz,
	i64 * recvnode, int recvnode_sz,
	int nivterm, int levcom,
	i64 * fsend, i64 * send, i64 * frecv, i64 * recv, 
	i64 * nst, i64 * nsp, i64 * fniv, i64 * endlev, i64 * codech,
	int includeLevcom, Gaspi_FF_communicator *& gCommFF);

	// multimat version
	void init_gaspi_ff(int max_send_terms, int max_recv_terms, int max_send_nodes, int max_recv_nodes, int nbMat, int incLevcom, Gaspi_FF_communicator *& gCommFF);



/**
 * GASPI TAGS - Notification values
 **/

// seg remoteIndexes
#define REMOTE_ADDRESS 1
#define REMOTE_INFO_ADDRESS 2

// seg global recv buffer
#define SEND_DATA 100
#define NO_DATA 101

#endif
