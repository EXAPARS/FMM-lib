#include "Gaspi_tools.hpp"
using namespace std;


void flush_queues(int nbQueues)
{
	
	double t_begin, t_end;
	t_begin = MPI_Wtime();
    
    gaspi_queue_id_t queue=0;       
	gaspi_number_t queueSizeMax;
	gaspi_number_t queueSize;
	
	// get maximum size
	SUCCESS_OR_DIE (gaspi_queue_size_max (&queueSizeMax));	
	
	// for each queue
	for(int i=0; i<nbQueues; i++)
	{
		queue=i;
		
		// get size
		SUCCESS_OR_DIE (gaspi_queue_size (queue, &queueSize));
		
		// if too much, throw error
		if (queueSize > queueSizeMax) 
		{
			gaspi_rank_t _rank; 
			gaspi_proc_rank(&_rank);
			cerr << "Rank " << _rank << " has exceeded its queue capacity." << endl;
			exit(1);
		}
		// if more than half, flush
		else if (queueSize >= queueSizeMax/2) 
		{
			SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
		}
	}
	
	// measure
	t_end = MPI_Wtime();
	add_time_sec("GASPI_REDUCE_UNK_flush_queue", t_end - t_begin);
}


void broadcast_to_global_buffer(int nbQueues, int localOffset, int offsetMultiple, int nbElts, int sizeOfElem,
	gaspi_rank_t _rank, gaspi_rank_t _wsize, gaspi_segment_id_t srcSeg, gaspi_segment_id_t destSeg, int notifValue, string timingMsg)
{
	
	// flush
	flush_queues(nbQueues);
	
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	gaspi_queue_id_t queue = 0;
	gaspi_offset_t local_offset = localOffset;
    gaspi_offset_t remote_offset;
   	gaspi_notification_id_t notif_offset = _wsize * offsetMultiple;
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_size_t qty= nbElts * sizeOfElem;
    
    // for each process, except the current one
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank ) // if not current rank
		{
			remote_offset = _rank * qty;
			
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					srcSeg,  				// local seg ID
					local_offset,			// local offset
					i,						// receiver rank
					destSeg,				// remote seg ID
					remote_offset,			// remote offset
					qty,					// size of data to write
					notifyID,				// remote notif ID
					notification,		// value of the notif to write
					queue,					// queue
					GASPI_BLOCK				// Gaspi block
				)
			);
		}		
		queue=(queue+1)%nbQueues;
	}
	
	//timing
	t_end = MPI_Wtime();
	add_time_sec(timingMsg, t_end - t_begin);
}
