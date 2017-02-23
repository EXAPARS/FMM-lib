/*
  Copyright 2016 - UVSQ
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

#ifndef GASPI_TOOLS
#define GASPI_TOOLS

#include "mpi.h"
#include "GASPI.h"
#include <iostream>

#include "../Tools/Complex.hpp"
#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/measure.hpp"
#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/byteCounter.hpp"

void flush_queues(int nbQueues);


// allReduce
void broadcast_to_global_buffer(int nbQueues, int localOffset, int offsetMultiple, int nbElts, int sizeOfElem,
	gaspi_rank_t _rank, gaspi_rank_t _wsize, gaspi_segment_id_t srcSeg, gaspi_segment_id_t destSeg, int notifValue, string timingMsg);

void receive_allReduce(int offsetMultiple, string timingPrefix, int nbElts,
	gaspi_rank_t _wsize, gaspi_segment_id_t destSeg, int notifValue, complex * buffer, complex * globalBuffer);

void gaspi_loop_broadcast();

template <typename T>
void copy_local_data(T * destBuffer, T * srcBuffer, int nbElts, string timingMsg)
{
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	
	int i;
	#pragma omp parallel for default(shared) private(i)
	for (i=0; i<nbElts; i++)
	{   
		destBuffer[i] = srcBuffer[i];
	}
	
	t_end = MPI_Wtime();
	add_time_sec(timingMsg + "_wb", t_end - t_begin);
}


// broadcast
void broadcast_buffer(int nbQueues, int offsetMultiple, int nbElts, int sizeOfElem,
	gaspi_rank_t _rank, gaspi_rank_t _wsize, gaspi_segment_id_t seg, int notifValue);

// clem 
void print_gaspi_config();

#define SUCCESS_OR_DIE(f...)\
do\
{\
	const gaspi_return_t r = f;\
	if (r != GASPI_SUCCESS)\
	{\
        printf("SUCCESS_OR_DIE failed at [%s:%i]\n",__FILE__,__LINE__);\
        gaspi_rank_t rank;\
        gaspi_proc_rank(&rank);\
		printf ("Error rank %d: '%s' [%s:%i]: %i\n",rank, #f, __FILE__, __LINE__, r);\
		exit (EXIT_FAILURE);\
	}\
} while (0)
#endif
