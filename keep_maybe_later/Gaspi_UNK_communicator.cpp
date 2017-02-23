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

#include "Gaspi_UNK_communicator.hpp"
using namespace std;
/*
Gaspi_UNK_communicator::Gaspi_UNK_communicator(complex * xtmp, complex * xtmp2, int nbUnk)
{
	// update class attributes
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);
	_nbUnknowns = nbUnk;

	// update segments ids
	_seg_loc_unk_id	= 29;
	_seg_glob_unk_id = 30;
	_seg_loc_unk_tmp_id	= 31;

	// class size attributes
	gaspi_size_t _seg_local_unk_size  = _nbUnknowns * sizeof(complex);
	gaspi_size_t _seg_global_unk_size = _nbUnknowns * _wsize * sizeof(complex);

	// create gaspi segments, by using F buffers
	SUCCESS_OR_DIE(
		gaspi_segment_use(
			_seg_loc_unk_id,
			xtmp,
			_seg_local_unk_size,
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			0
		)
	);

	SUCCESS_OR_DIE(
		gaspi_segment_use(
			_seg_loc_unk_tmp_id,
			xtmp2,
			_seg_local_unk_size,
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			0
		)
	);
	
	// create gaspi segment
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_glob_unk_id,
			_seg_global_unk_size,
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			0
		)
	);
	
	// gaspi pointers
	gaspi_segment_ptr(_seg_loc_unk_id, &_ptr_seg_loc_unk);
	gaspi_segment_ptr(_seg_loc_unk_tmp_id, &_ptr_seg_loc_unk_tmp);
	gaspi_segment_ptr(_seg_glob_unk_id, &_ptr_seg_glob_unk);

	// user pointers
	_unknowns 		= (complex *) _ptr_seg_loc_unk;
	_unknownsTmp 	= (complex *) _ptr_seg_loc_unk_tmp;
	_globalUnknowns = (complex *) _ptr_seg_glob_unk;
}*/
/*
void::Gaspi_UNK_communicator::runAllReduceUnknowns()
{
	int nbQueues = 1;
	int offsetMultiple = 4;
	int localOffset = 0;
	double t_begin, t_end;

	// SEND
	broadcast_to_global_buffer(
		nbQueues, localOffset, offsetMultiple, _nbUnknowns, sizeof(complex), 
		_rank, _wsize, _seg_loc_unk_tmp_id, _seg_glob_unk_id, ALLREDUCE_UNKNOWNS, 
		"GASPI_REDUCE_UNK_write_notify");

	// RECV and reduce on _unknowns array
	copy_local_data<complex>(_unknowns, _unknownsTmp, _nbUnknowns, "GASPI_REDUCE_UNK");
	receive_allReduce(offsetMultiple, "GASPI_REDUCE_UNK", _nbUnknowns, _wsize, _seg_glob_unk_id, ALLREDUCE_UNKNOWNS, _unknowns, _globalUnknowns);
}

void Gaspi_UNK_communicator::runBroadcastUnknowns()
{
	int nbQueues = 1;
	int offsetMultiple = 3;
	broadcast_buffer(nbQueues, offsetMultiple, _nbUnknowns, sizeof(complex), _rank, _wsize, _seg_loc_unk_id, BROADCAST_UNKNOWNS);
}
*/
