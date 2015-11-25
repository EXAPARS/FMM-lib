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

#ifndef GASPI_COMMUNICATOR
#define GASPI_COMMUNICATOR

#include "GASPI.h"
#include "vec3D.hpp"

#include <iostream>
using namespace std;


class Gaspi_communicator
{	
private:
	// segments sizes
	gaspi_size_t _seg_RecvBuffer_size;	
	gaspi_size_t _seg_LocalBuffer_size;
	gaspi_size_t _seg_SepNodes_size;
	gaspi_size_t _seg_NbUntilNode_size;
	gaspi_size_t _seg_InitCoords_size;
	gaspi_size_t _seg_NewCoords_size;
	gaspi_size_t _seg_CommInfos_size;

	// gaspi pointers on segments
	gaspi_pointer_t _ptr_seg_RecvBuffer = 	nullptr;
	gaspi_pointer_t _ptr_seg_LocalBuffer = 	nullptr;
	gaspi_pointer_t _ptr_seg_SepNodes = 	nullptr;
	gaspi_pointer_t _ptr_seg_NbUntilNode = 	nullptr;
	gaspi_pointer_t _ptr_seg_InitCoords = 	nullptr;
	gaspi_pointer_t _ptr_seg_NewCoords = 	nullptr;
	gaspi_pointer_t _ptr_seg_CommInfos = 	nullptr;

public:
	// processes info
	gaspi_rank_t _wsize;
	gaspi_rank_t _rank;

	// segments Ids
	const gaspi_segment_id_t _seg_RecvBuffer_id = 	0;// Receive Buffer 
	const gaspi_segment_id_t _seg_LocalBuffer_id = 	1;// Local Buffer
	const gaspi_segment_id_t _seg_SepNodes_id = 	2;// SepNodes Buffer
	const gaspi_segment_id_t _seg_NbUntilNode_id = 	3;// NbUntilNode Buffer
	const gaspi_segment_id_t _seg_InitCoords_id = 	4;// Coords
	const gaspi_segment_id_t _seg_NewCoords_id = 	5;// newCoords
	const gaspi_segment_id_t _seg_CommInfos_id = 	6;// info
	
	// usual pointers on values in segments
	int * _recvBuffer = nullptr;
	int * _localBuffer = nullptr;
	int64_t * _sepNodes = nullptr;
	int * _nbUntilNode = nullptr;
	vec3D * _initCoords = nullptr;
	vec3D * _newCoords = nullptr;
	int * _commInfos = nullptr;

public:
	Gaspi_communicator(int nbLeaves, int nbParticles);
	void initNewCoords(int maxNbItems);

};


/**
 * GASPI TAGS - Notification value
 **/ 
// seg recvBuffer
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
	
// seg commInfos
/// nb de particules envoyées /!\ max = unsigned int : 2 milliards and bananas



#endif
