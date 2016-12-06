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

#ifndef GASPI_UNKNOWNS_COMMUNICATOR
#define GASPI_UNKNOWNS_COMMUNICATOR

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

class Gaspi_unknowns_communicator
{	
public:
	// segments sizes (in bytes for segment creation)
	gaspi_size_t _seg_local_unk_size; 
	gaspi_size_t _seg_global_unk_size; 

	// gaspi pointers on segments
	gaspi_pointer_t _ptr_seg_loc_unk = nullptr;
	gaspi_pointer_t _ptr_seg_loc_unk_tmp = nullptr;
    gaspi_pointer_t _ptr_seg_glob_unk = nullptr;

	// processes info
	gaspi_rank_t _wsize;
	gaspi_rank_t _rank;
	int _nbUnknowns;

	// segments Ids
    gaspi_segment_id_t _seg_loc_unk_id; 
    gaspi_segment_id_t _seg_loc_unk_tmp_id; 
    gaspi_segment_id_t _seg_glob_unk_id; 
    
	// usual pointers on values in segments
	complex * _unknowns = nullptr;
    complex * _unknownsTmp = nullptr;
	complex * _globalUnknowns = nullptr;
	
public:
	Gaspi_unknowns_communicator(complex * xtmp, complex * xtmp2, int nbUnk);
	void runBroadcastUnknowns();
	void runAllReduceUnknowns();
	void receive_and_update_allReduce();
};


/**
 * GASPI TAGS - Notification values
 **/

/* 
// seg remoteIndexes
#define REMOTE_ADDRESS 1
// seg global recv bugger
#define SEND_DATA 10
#define NO_DATA 20
#define ALLREDUCE 30
*/

// seg unknowns
#define BROADCAST_UNKNOWNS 40
#define ALLREDUCE_UNKNOWNS 50

#endif
