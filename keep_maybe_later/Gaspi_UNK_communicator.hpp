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

#ifndef Gaspi_UNK_communicator
#define Gaspi_UNK_communicator

#include "GASPI.h"
#include "Gaspi_tools.hpp"

#include "mpi.h"
#include "../Tools/types.hpp"
/*#include "../Tools/Complex.hpp"
#include "../Tools/fmm_tools.hpp"

#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/measure.hpp"
#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/byteCounter.hpp"
*/

#include <iostream>
using namespace std;
/*
class Gaspi_UNK_communicator
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
	//Gaspi_UNK_communicator(complex * xtmp, complex * xtmp2, int nbUnk);
	Gaspi_UNK_communicator();
	void runBroadcastUnknowns();
	void runAllReduceUnknowns();
	void receive_and_update_allReduce();
};*/

#endif
