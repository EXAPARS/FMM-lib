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

Gaspi_unknowns_communicator::Gaspi_unknowns_communicator(complex * xtmp, complex * xtmp2, int nbUnk)
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
}

void::Gaspi_unknowns_communicator::runAllReduceUnknowns()
{
	// Gaspi queue
	int nbQueues = 1;
	/*flush_queues(nbQueues);*/

	double t_begin, t_end;

	// Gaspi broadcast to global buffer
	broadcast_to_global_buffer(nbQueues, 0, 4, _nbUnknowns, sizeof(complex), _rank, _wsize, _seg_loc_unk_tmp_id, _seg_glob_unk_id, ALLREDUCE_UNKNOWNS, "GASPI_REDUCE_UNK_write_notify");

/*
	 send data 
	int nbQueues = 1;
	gaspi_queue_id_t queue = 0;
	gaspi_offset_t local_offset = 0;
    gaspi_offset_t remote_offset;
   	gaspi_notification_id_t notif_offset = _wsize * 4;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_size_t qty= _nbUnknowns * sizeof(complex);
    
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank ) // if not current rank
		{
			remote_offset = _rank * qty;
			
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					_seg_loc_unk_tmp_id,  				// local seg ID
					local_offset,						// local offset
					i,									// receiver rank
					_seg_glob_unk_id,					// remote seg ID
					remote_offset,						// remote offset
					qty,								// size of data to write
					notifyID,							// remote notif ID
					ALLREDUCE_UNKNOWNS,					// value of the notif to write
					queue,								// queue
					GASPI_BLOCK							// Gaspi block
				)
			);
		}		
		queue=(queue+1)%nbQueues;
	}
	t_end = MPI_Wtime();
	
	add_time_sec("GASPI_REDUCE_UNK_write_notify", t_end - t_begin);
*/	

   	gaspi_notification_id_t notif_offset = _wsize * 4;   
    
	// receptions
	 // write its own data into the result !
	t_begin = MPI_Wtime();    
	
	int i;
	#pragma omp parallel for default(shared) private(i)
	for (i=0; i<_nbUnknowns; i++)
	{   
		_unknowns[i] = _unknownsTmp[i];
	}
	t_end = MPI_Wtime();
	add_time_sec("GASPI_REDUCE_UNK_write_back_unk", t_end - t_begin);

	// wait to receive all messages from the others
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
    int sender;
    double t_begin_loop, t_end_loop;
    
    for(int i=0; i< (_wsize-1); i++)
	{
		//methode 1 - avec GASPI_BLOCK
		while(1)
		{
			t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_seg_glob_unk_id,
					notif_offset,			// surveille les notifications depuis offset _wsize
					_wsize,					// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_REDUCE_UNK_notify_waitsome", t_end_loop - t_begin_loop);
			t_begin_loop = MPI_Wtime();
			
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_seg_glob_unk_id, 
					new_notif_id, 
					&new_notif_val
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_REDUCE_UNK_notify_reset", t_end_loop - t_begin_loop);
			t_begin_loop = MPI_Wtime();
			
			if (new_notif_val) 
				break;
		}
		
		sender = new_notif_id - notif_offset;
		
		// test the notification value and update counter
		t_begin_loop = MPI_Wtime();
	    if (new_notif_val == ALLREDUCE_UNKNOWNS)
	    {
		    // update the far field array		    		     
            int offset = _nbUnknowns * sender;
            int j;
            #pragma omp parallel for default(shared) private (j)
	        for (j=0; j<_nbUnknowns; j++)
	        {   
				_unknowns[j] = _unknowns[j] + _globalUnknowns[offset + j];
	        }
        }	
		t_end_loop = MPI_Wtime();
		add_time_sec("GASPI_REDUCE_UNK_write_back_unk", t_end_loop - t_begin_loop);
	} 
}

void Gaspi_unknowns_communicator::runBroadcastUnknowns()
{
	double t_begin, t_end;
	
	/** avec gestion de la queue **/
    gaspi_queue_id_t queue=0;
    int nbQueues = 1;        
	t_begin = MPI_Wtime();
	gaspi_number_t queueSizeMax;
	gaspi_number_t queueSize;
	
	SUCCESS_OR_DIE (gaspi_queue_size_max (&queueSizeMax));
	
	for(int i=0; i<nbQueues; i++)
	{
		queue=i;
		SUCCESS_OR_DIE (gaspi_queue_size (queue, &queueSize));
		if (queueSize > queueSizeMax) 
		{
			cerr << "Rank " << _rank << " has exceeded its queue capacity." << endl;
			exit(1);
		}
		else if (queueSize >= queueSizeMax/2) 
		{
			SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
		}
	}
	t_end = MPI_Wtime();
	add_time_sec("GASPI_BROADCAST_flush_queue", t_end - t_begin);
	t_begin = MPI_Wtime();
	queue = 0;
	
    // send infos  
    gaspi_offset_t local_offset = 0;
    gaspi_offset_t remote_offset = 0;
   	gaspi_notification_id_t notif_offset = _wsize * 3;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_size_t qty = _nbUnknowns * sizeof(complex);  
    
    if (_rank == 0)
    {
		for (int i=1; i<_wsize; i++) // broadcast to others only
		{
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					_seg_loc_unk_id,					// local seg ID
					local_offset,						// local offset
					i,									// receiver rank
					_seg_loc_unk_id,					// remote seg ID
					remote_offset,						// remote offset
					qty,								// size of data to write
					notifyID,							// remote notif ID
					BROADCAST_UNKNOWNS,					// value of the notif to write
					queue,								// queue
					GASPI_BLOCK							// Gaspi block
				)
			);
		}
		queue=(queue+1)%nbQueues;
	}

	t_end = MPI_Wtime();
	add_time_sec("GASPI_BROADCAST_write_notify", t_end - t_begin);
	t_begin = MPI_Wtime();

	// wait to receive the broadcasted unknowns array
	if (_rank != 0)
	{
		gaspi_notification_id_t new_notif_id;
		gaspi_notification_t new_notif_val;
		double t_begin_loop, t_end_loop;

		//methode 1 - avec GASPI_BLOCK
		while(1)
		{
			t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_seg_loc_unk_id,
					notif_offset,				// surveille les notifications depuis 0
					1,							// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_BROADCAST_notify_waitsome", t_end_loop - t_begin_loop);
			t_begin_loop = MPI_Wtime();
			
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_seg_loc_unk_id, 
					new_notif_id, 
					&new_notif_val
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_BROADCAST_notify_reset", t_end_loop - t_begin_loop);
			
			if (new_notif_val == BROADCAST_UNKNOWNS) 
				break;
		}
	}
}
