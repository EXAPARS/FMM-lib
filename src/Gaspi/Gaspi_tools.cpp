#include "Gaspi_tools.hpp"
using namespace std;


/* Attention si 512 MPI, queue déjà pleine ... */
void flush_queues(int nbQueues)
{
	   
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
		//cout << "queueSize : " << queueSize << endl;
		
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
	
}


void broadcast_to_global_buffer(int nbQueues, int localOffset, int offsetMultiple, int nbElts, int sizeOfElem,
	gaspi_rank_t _rank, gaspi_rank_t _wsize, gaspi_segment_id_t srcSeg, gaspi_segment_id_t destSeg, int notifValue, string timingMsg)
{
	
	// flush
	flush_queues(nbQueues);
	
	//double t_begin, t_end;
	//t_begin = MPI_Wtime();
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
					notifValue,				// value of the notif to write
					queue,					// queue
					GASPI_BLOCK				// Gaspi block
				)
			);
		}		
		queue=(queue+1)%nbQueues;
	}
	
	//timing
	//~ t_end = MPI_Wtime();
//~ #ifdef TIMING	
	//~ add_time_sec(timingMsg, t_end - t_begin);
//~ #endif
}

void broadcast_buffer(int nbQueues, int offsetMultiple, int nbElts, int sizeOfElem,
	gaspi_rank_t _rank, gaspi_rank_t _wsize, gaspi_segment_id_t seg, unsigned int notifValue)
{
	//double t_begin, t_end;
	
	// flush
	flush_queues(nbQueues);
	
    // send infos  
	gaspi_queue_id_t queue = 0;
    gaspi_offset_t local_offset = 0;
    gaspi_offset_t remote_offset = 0;
   	gaspi_notification_id_t notif_offset = _wsize * offsetMultiple;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_size_t qty = nbElts * sizeOfElem;
    
    
    if (_rank == 0)
    {
		for (int i=1; i<_wsize; i++) // broadcast to others only
		{
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					seg,					// local seg ID
					local_offset,			// local offset
					i,						// receiver rank
					seg,					// remote seg ID
					remote_offset,			// remote offset
					qty,					// size of data to write
					notifyID,				// remote notif ID
					notifValue,				// value of the notif to write
					queue,					// queue
					GASPI_BLOCK				// Gaspi block
				)
			);
		}
		queue=(queue+1)%nbQueues;
	}

	// wait to receive the broadcasted unknowns array
	if (_rank != 0)
	{
		gaspi_notification_id_t new_notif_id;
		gaspi_notification_t new_notif_val;
		//double t_begin_loop, t_end_loop;

		//methode 1 - avec GASPI_BLOCK
		while(1)
		{
			//t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					seg,
					notif_offset,				// surveille les notifications depuis 0
					1,							// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					seg, 
					new_notif_id, 
					&new_notif_val
				)
			);

			if (new_notif_val == notifValue) 
				break;
		}
	}
}

void receive_allReduce(int offsetMultiple, string timingPrefix, int nbElts,
	gaspi_rank_t _wsize, gaspi_segment_id_t destSeg, unsigned int notifValue, complex * buffer, complex * globalBuffer)
{
	//double t_begin, t_end;
	//t_begin = MPI_Wtime();
	
	// wait to receive all messages from the others
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	gaspi_notification_id_t notif_offset = _wsize * offsetMultiple;   

    int sender;
    //double t_begin_loop, t_end_loop;
    
    // receive from all ranks except 1
    for(int i=0; i< (_wsize-1); i++)
	{
		while(1)
		{
			//t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					destSeg,
					notif_offset,			// surveille les notifications depuis offset _wsize
					_wsize,					// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			//t_end_loop = MPI_Wtime();
//~ #ifdef TIMING			
			//~ add_time_sec(timingPrefix + "_notify_waitsome", t_end_loop - t_begin_loop);
//~ #endif
			//t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					destSeg, 
					new_notif_id, 
					&new_notif_val
				)
			);
			//t_end_loop = MPI_Wtime();
//~ #ifdef TIMING			
			//~ add_time_sec(timingPrefix + "_notify_reset", t_end_loop - t_begin_loop);
//~ #endif
			if (new_notif_val) 
				break;
		}
		
		sender = new_notif_id - notif_offset;
		
		// test the notification value and update buffer
		//t_begin_loop = MPI_Wtime();
	    if (new_notif_val == notifValue)
	    {
		    // update the far field array		    		     
            int offset = nbElts * sender;
            int j;
            //#pragma omp parallel for default(shared) private (j)
	        for (j=0; j<nbElts; j++)
	        {   
				buffer[j] = buffer[j] + globalBuffer[offset + j];
	        }
        }	
		//t_end_loop = MPI_Wtime();
//~ #ifdef TIMING		
		//~ add_time_sec(timingPrefix + "_write_back_unk", t_end_loop - t_begin_loop);
//~ #endif
	}
	
	// timing
	//~ t_end = MPI_Wtime();
//~ #ifdef TIMING	
	//~ add_time_sec(timingPrefix + "_write_back_unk", t_end - t_begin); 
//~ #endif
}

/*void gaspi_loop_broadcast(int nbElts)
{
	flush_queues(nbQueues);
	
    // send infos  
	gaspi_queue_id_t queue = 0;
    gaspi_offset_t local_offset = 0;
    gaspi_offset_t remote_offset = 0;
   	gaspi_notification_id_t notif_offset = _wsize * offsetMultiple;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_size_t qty = nbElts * sizeOfElem;
    
    for (int i=0; i<nbElts; i+=255)
    {
		// tail
		if (nbElts-1 <=255)
		{
			gaspi_broadcast()
		}
		else // chunk
		{
			255
		}
	}
    
    
}*/

/** CLEM GASPI TOOLS **/
//~ void print_gaspi_config()
//~ {
	//~ gaspi_config_t config; 
	//~ gaspi_config_get(&config);
	
	//~ cout << "logger=" << config.logger << std::endl;
	//~ cout << "sn_port=" << config.sn_port << std::endl;
	//~ cout << "net_info=" << config.net_info << std::endl;
	//~ cout << "netdev_id=" << config.netdev_id << std::endl;
	//~ cout << "mtu=" << config.mtu << std::endl;
	//~ cout << "port_check=" << config.port_check << std::endl;
	//~ cout << "user_net=" << config.user_net << std::endl;
	//~ cout << "network=" << config.network << std::endl;
	//~ cout << "queue_depth=" << config.queue_depth << std::endl;
	//~ cout << "queue_num=" << config.queue_num << std::endl;
	//~ cout << "group_max=" << config.group_max << std::endl;
	//~ cout << "segment_max=" << config.segment_max << std::endl;
	//~ cout << "transfer_size_max=" << config.transfer_size_max << std::endl;
	//~ cout << "notification_num=" << config.notification_num << std::endl;
	//~ cout << "passive_queue_size_max=" << config.passive_queue_size_max << std::endl;
	//~ cout << "passive_transfer_size_max=" << config.passive_transfer_size_max << std::endl;
	//~ cout << "allreduce_buf_size=" << config.allreduce_buf_size << std::endl;
	//~ cout << "allreduce_elem_max=" << config.allreduce_elem_max << std::endl;
	//~ cout << "build_infrastructure=" << config.build_infrastructure << std::endl;
//~ }
