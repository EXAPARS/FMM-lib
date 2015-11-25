#include "Gaspi_communicator.hpp"
using namespace std;


Gaspi_communicator::Gaspi_communicator(int nbLeaves, int nbParticles)	
{		
	// processes info
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);
	
	//nbSeps
	int nbSeps = _wsize - 1;
	
	// segments sizes
	_seg_RecvBuffer_size = _wsize * nbLeaves * sizeof(int);
	_seg_LocalBuffer_size = _wsize * nbLeaves * sizeof(int);
	_seg_SepNodes_size = nbSeps * sizeof(int64_t);
	_seg_NbUntilNode_size = nbSeps * sizeof(int);
	_seg_InitCoords_size = nbParticles * 3 * sizeof(double);
	_seg_CommInfos_size = sizeof(int);

	// segments creations	
	gaspi_segment_create( _seg_RecvBuffer_id, _seg_RecvBuffer_size, 
		GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT);			
	gaspi_segment_create( _seg_LocalBuffer_id, _seg_LocalBuffer_size, 
		GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT);
	gaspi_segment_create( _seg_SepNodes_id, _seg_SepNodes_size, 
		GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT);
	gaspi_segment_create( _seg_NbUntilNode_id, _seg_NbUntilNode_size, 
		GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT);
	gaspi_segment_create( _seg_InitCoords_id, _seg_InitCoords_size, 
		GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT);
	gaspi_segment_create( _seg_CommInfos_id, _seg_CommInfos_size, 
		GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT);		
	
	// update segpointers
	gaspi_segment_ptr( _seg_RecvBuffer_id, 	&_ptr_seg_RecvBuffer 	);
	gaspi_segment_ptr( _seg_LocalBuffer_id,	&_ptr_seg_LocalBuffer 	);
	gaspi_segment_ptr( _seg_SepNodes_id,   	&_ptr_seg_SepNodes 		);
	gaspi_segment_ptr( _seg_NbUntilNode_id,	&_ptr_seg_NbUntilNode 	);
	gaspi_segment_ptr( _seg_InitCoords_id, 	&_ptr_seg_InitCoords 	);
	gaspi_segment_ptr( _seg_CommInfos_id,	&_ptr_seg_CommInfos 	);	
	
	// update usual pointers
	_recvBuffer = (int *)	_ptr_seg_RecvBuffer;
	_localBuffer = (int *)	_ptr_seg_LocalBuffer;
	_sepNodes = (int64_t *)	_ptr_seg_SepNodes;
	_nbUntilNode = (int *)	_ptr_seg_NbUntilNode;
	_initCoords = (vec3D *)	_ptr_seg_InitCoords;
	_commInfos = (int *)	_ptr_seg_CommInfos;
		
	// segments initializations, /!\ : uninitialized segments are already allocated with calloc
	for(int i=0; i<nbSeps; i++)
		_sepNodes[i] = -1;		
}


void Gaspi_communicator::initNewCoords(int newNbParticles)
{
	_seg_NewCoords_size = newNbParticles * 3 * sizeof(double);			// size
	gaspi_segment_create( _seg_NewCoords_id, _seg_NewCoords_size, 		// create
		GASPI_GROUP_ALL, GASPI_BLOCK, GASPI_ALLOC_DEFAULT);	
	gaspi_segment_ptr( _seg_NewCoords_id,	&_ptr_seg_NewCoords	);		// pointers
	_newCoords = (vec3D *)	_ptr_seg_NewCoords;		
}	


/// TODO : pour pouvoir virer la barriere dans Morton ASYNC, faire le segment create en 2 fonctions
/// /!\ part du principe que DIM = 3
