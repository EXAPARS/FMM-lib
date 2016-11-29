#include "Gaspi_M2L_communicator.hpp"
using namespace std;

Gaspi_m2l_communicator::Gaspi_m2l_communicator(
	i64 * nb_send, int nb_send_sz, 
	i64 * nb_recv, int nb_recv_sz,
	i64 * sendnode, int sendnode_sz,
	i64 * recvnode, int recvnode_sz,
	complex * ff, complex * ne, int nbEltsToReduce)
{
	// update class attributes
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);
		
	// update segments ids
	_seg_ff_allreduce_id			= 15;
	_seg_ne_allreduce_id			= 16;
	_seg_globalRecvBuffer_id 	    = 17; // Receive Buffer 
	_seg_globalSendBuffer_id 	    = 18; // Send Buffer 
    _seg_remoteBufferIndexes_id     = 19; // Index where to write on other ranks	
    _seg_globalRecvBufIdxPerRank_id = 20; // Index in global recv buffer, per RANK
										  // temporary, only for initialization computations
										  
	// allocate array
	_sendBufferIndexes = new int[_wsize]();
	
	// create gaspi segments, and initialize them
	create_allReduceBuffers(ff, ne, nbEltsToReduce);
	create_globalRecvBuffer(nb_recv, nb_recv_sz);
	create_remoteBufferIndexes(recvnode, recvnode_sz, nb_recv);
	create_globalSendBuffer(nb_send, nb_send_sz);
	init_sendBufferIndexes(sendnode, sendnode_sz, nb_send);
}

void Gaspi_m2l_communicator::create_allReduceBuffers(complex * ff, complex * ne, int nbEltsToReduce)
{
	// create Gaspi segments
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_ff_allreduce_id,
			nbEltsToReduce*sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_ne_allreduce_id,
			nbEltsToReduce*sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	
	// gaspi pointers
	gaspi_segment_ptr(_seg_ff_allreduce_id, &_ptr_seg_ff_allreduce);
	gaspi_segment_ptr(_seg_ne_allreduce_id, &_ptr_seg_ne_allreduce);
	
	// user pointers
	_reduceNE= (complex *)_ptr_seg_ne_allreduce;
	_reduceFF = (complex *) _ptr_seg_ff_allreduce;
	
	// class size attributes
	_seg_reduce_size = nbEltsToReduce * sizeof(complex);
	_seg_globalRecvBuffer_size = _seg_reduce_size * _wsize;
}

void Gaspi_m2l_communicator::create_globalRecvBuffer(i64 * nb_recv, int nb_recv_sz)
{
 	// create segment
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_globalRecvBuffer_id,
			_seg_globalRecvBuffer_size,
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);

	// update gaspi segpointer
	gaspi_segment_ptr(_seg_globalRecvBuffer_id, &_ptr_seg_globalRecvBuffer);
	
	// update user pointer
	_globalRecvBuffer = (complex *)_ptr_seg_globalRecvBuffer;
}

void Gaspi_m2l_communicator::create_remoteBufferIndexes(i64 * recvnode, int recvnode_sz, i64 * nb_recv)
{
	// From Fortran to C/C++
	int indexToC = -1;
	
    /* 
     * GlobalRecvBufferIndexes organised by communication ROUND,
     * since recvnode is organised by communication round
     */
     
	int * globalRecvBufferIndexPerROUND = nullptr;
	globalRecvBufferIndexPerROUND = new int[recvnode_sz]();
	globalRecvBufferIndexPerROUND[0]=0;

	for (int i=1; i<recvnode_sz; i++)
	{
		if (recvnode[i-1] == 0)
			globalRecvBufferIndexPerROUND[i] = globalRecvBufferIndexPerROUND[i-1];
	    else
			globalRecvBufferIndexPerROUND[i] = globalRecvBufferIndexPerROUND[i-1] + nb_recv[recvnode[i-1]+indexToC]; 
    }
    
    /* 
     * Prepare GlobalRecvBufferIndexes SEGMENT organised by process RANK 
     */
    
    // update segment size
    _seg_globalRecvBufIdxPerRank_size = _wsize * sizeof(int);
    
    // create segment
    SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_globalRecvBufIdxPerRank_id,
			_seg_globalRecvBufIdxPerRank_size,
			GASPI_GROUP_ALL,
			GASPI_BLOCK,
			GASPI_ALLOC_DEFAULT
		)
	);
    
    // update gaspi segpointers
    gaspi_segment_ptr( _seg_globalRecvBufIdxPerRank_id, &_ptr_seg_globalRecvBufIdxPerRank);
	
	// update user pointers
    _globalRecvBufIdxPerRank = (int *)_ptr_seg_globalRecvBufIdxPerRank;    
    
    // initialize segment with values
    for (int k=0; k<recvnode_sz; k++)
    {
		if (recvnode[k] > 0)
		{
			int from = recvnode[k]-1;
			if (from != _rank)
			{
				_globalRecvBufIdxPerRank[from] = globalRecvBufferIndexPerROUND[k];
			}
			else
			{
				cerr << "[Create_remoteBufferIndexes] Fortran recvnode array should not contain the current rank !"; exit(0);
			}
		}
    }
	_globalRecvBufIdxPerRank[_rank] = -1;
    
    /* 
     * Prepare Remote Buffer Indexes SEGMENT 
     */
     
    // update segment size
    _seg_remoteBufferIndexes_size = _wsize * sizeof(int);

	// create gaspi segment
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_remoteBufferIndexes_id, 
			_seg_remoteBufferIndexes_size,
			GASPI_GROUP_ALL, 
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
    );
    
    // update gaspi segpointer
    gaspi_segment_ptr(_seg_remoteBufferIndexes_id, &_ptr_seg_remoteBufferIndexes); 
	
	// update user pointer
    _remoteBufferIndexes = (int *) _ptr_seg_remoteBufferIndexes;
	
    /* 
     * Data Exchange to Update the Remote Buffer Indexes segment
     */
    
    int local_offset;
    int remote_offset = _rank * sizeof(int);
    
    // send infos
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank)
		{
			local_offset = i * sizeof(int);
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					_seg_globalRecvBufIdxPerRank_id,	// local seg ID
					local_offset,						// local offset
					i,									// receiver rank
					_seg_remoteBufferIndexes_id,		// remote seg ID
					remote_offset,						// remote offset
					sizeof(int),						// size of data to write
					_rank,								// remote notif ID
					REMOTE_ADDRESS,						// value of the notif to write
					0,									// queue
					GASPI_BLOCK							// Gaspi block
				)
			);
		}
	}
	
	// wait to receive all infos
	int cpt = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;

	while (cpt < (_wsize-1))
	{
		// wait for notification
		if(gaspi_notify_waitsome(
			_seg_remoteBufferIndexes_id,
			0,							            // surveille les notifications depuis 0
			_wsize,						            // en surveille wsize
			&new_notif_id,
			GASPI_BLOCK) == GASPI_SUCCESS)
		{

			// get notification value, and reset
			gaspi_notify_reset(_seg_remoteBufferIndexes_id, new_notif_id, &new_notif_val);
		
			// test notification and increase counter
			if (new_notif_val == REMOTE_ADDRESS)
				cpt++;
			else
			{
				cerr << "[Create_remoteBufferIndexes] Unexpected gaspi msg."; 
				exit(0);
			}
		}
	}
}

void Gaspi_m2l_communicator::create_globalSendBuffer(i64 * nb_send, int nb_send_sz)
{
	// compute size
	int globalSendBufferSize = 0;
	for (int i=0; i<nb_send_sz; i++)
	    globalSendBufferSize = globalSendBufferSize + nb_send[i];
   
 	// update segment size
 	_seg_globalSendBuffer_size = globalSendBufferSize * sizeof(complex);
 	//debug("SendSegment size : " + convert(_seg_globalSendBuffer_size));

 	// create segment
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_globalSendBuffer_id,
			_seg_globalSendBuffer_size,
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	
	// update gaspi segpointer
	gaspi_segment_ptr(_seg_globalSendBuffer_id, &_ptr_seg_globalSendBuffer);
	
	// update user pointer
	_globalSendBuffer = (complex *)_ptr_seg_globalSendBuffer;
}

void Gaspi_m2l_communicator::initGlobalSendSegment(
	i64 * sendnode, int sendnode_sz, i64 * nb_send, 
	int nivterm, int levcom, i64 * fsend, i64 * send, i64 * endlev,
	i64 * codech, i64 * nst, i64 * nsp, complex * bufsave, i64 * fniv, complex * ff)
{
	int indexToC = -1;
	
	// Global send buffer indexes
	int * globalSendBufferIndex = nullptr;
	globalSendBufferIndex = new int[sendnode_sz]();
	
	for (int i=1; i<sendnode_sz; i++)
	{
		if (sendnode[i-1] == 0)
			globalSendBufferIndex[i] = globalSendBufferIndex[i-1];
	    else
			globalSendBufferIndex[i] = globalSendBufferIndex[i-1] + nb_send[sendnode[i-1]+indexToC];
	}
	
	// Fill Global send Buffer
	for (int i=0; i<sendnode_sz; i++)
	{
		int iDest = sendnode[i];
		int k = nivterm;
		int q = globalSendBufferIndex[i]-1;
		
		if (iDest > 0)
		{
			for (int j=fsend[iDest+1+indexToC]-1; j>=fsend[iDest+indexToC]; j--)
			{
				int jc = j + indexToC;
				
				bool todo = true;
				int p;
				int l = send[jc];
				if ((levcom<=3 || (levcom == 4 && nivterm>8))  && (l < endlev[levcom + indexToC]))
					todo = false;
				if ((levcom>4  || (levcom == 4 && nivterm<=8)) && (l <= endlev[levcom-1 + indexToC]))
					todo = false;
					
				if (todo)
				{
					// se place au bon niveau
					while(l <= endlev[k-1 + indexToC])
						k--;
					
					// cas 1
					if (codech[l + indexToC] > 1)
					{
						p=codech[l + indexToC]-2;
						for (int st=0; st<nst[k + indexToC]; st++)
						{
							for (int sp=0; sp<nsp[k + indexToC]; sp++)
							{
								p++;
								q++;
								_globalSendBuffer[q]=bufsave[p + indexToC];
							}
						}
					}
					// cas 2
					else
					{
						p=fniv[k+1 + indexToC]+(endlev[k + indexToC]-l)*nst[k + indexToC]*nsp[k + indexToC];
						for (int st=0; st<nst[k + indexToC]; st++)
						{
							for (int sp=0; sp<nsp[k + indexToC];  sp++)
							{
								p++;
								q++;
								_globalSendBuffer[q]=ff[p + indexToC];
							}
						}
					}
				}
			}
		}
	} // end fill global send buffer
	
	delete [] globalSendBufferIndex;
}


void Gaspi_m2l_communicator::initAllReduceBuffers(
	complex * ff, complex * ne)
{
	int nbElts = _seg_reduce_size / sizeof(complex);
	// ff
	for (int i=0; i<nbElts; i++)
		_reduceFF[i] = ff[i];
	
	// ne
	for (int i=0; i<nbElts; i++)
		_reduceNE[i] = 0;
}

void Gaspi_m2l_communicator::init_sendBufferIndexes(i64 * sendnode, int sendnode_sz, i64 * nb_send)
{
	// From Fortran to C/C++
	int indexToC = -1;
	
    /* 
     * GlobalSendBufferIndexes organised by communication ROUND,
     * since sendnode is organised by communication round
     */
     
	int * globalSendBufferIndexPerROUND = nullptr;
	globalSendBufferIndexPerROUND = new int[sendnode_sz]();
	globalSendBufferIndexPerROUND[0]=0;

	for (int i=1; i<sendnode_sz; i++)
	{
		if (sendnode[i-1] == 0)
			globalSendBufferIndexPerROUND[i] = globalSendBufferIndexPerROUND[i-1];
	    else
			globalSendBufferIndexPerROUND[i] = globalSendBufferIndexPerROUND[i-1] + nb_send[sendnode[i-1]+indexToC]; 
    }
    
    /* 
     * Prepare _sendBufferIndexes organised by process RANK 
     */

    for (int k=0; k<sendnode_sz; k++)
    {
		if(sendnode[k] > 0)
		{
			int to = sendnode[k]-1;
			if (to != _rank)
			{
				_sendBufferIndexes[to] = globalSendBufferIndexPerROUND[k];
			}
			else
			{
				cerr << "Error in Gaspi_m2l_communicator::init_dataToSendIndexes.\nFortran sendnode array should not contain the current rank !";
				exit(0);
			}
		}
    }  
	_sendBufferIndexes[_rank] = -1;
	
	// delete
	delete [] globalSendBufferIndexPerROUND;
}

void Gaspi_m2l_communicator::runM2LallReduce(complex * ff, complex * ne)
{
	double t_begin, t_end;

	// init ff and ne segments
	t_begin = MPI_Wtime();
	initAllReduceBuffers(ff, ne);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_REDUCE_write_reduce_sendSegment", t_end - t_begin);
	

	/** gestion de la queue **/
	t_begin = MPI_Wtime();
    int nbQueues = 1;
    gaspi_queue_id_t queue=0;       
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
	add_time_sec("GASPI_REDUCE_flush_queue", t_end - t_begin);
	t_begin = MPI_Wtime();
	queue = 0;

	/** send data **/
	gaspi_offset_t local_offset=0;
    gaspi_offset_t remote_offset;
   	gaspi_notification_id_t notif_offset = _wsize;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_size_t qty= _seg_reduce_size;
	
	// send infos    
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank ) // if not current rank
		{
			remote_offset = _rank * _seg_reduce_size;
			
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					_seg_ff_allreduce_id,  			// local seg ID
					local_offset,						// local offset
					i,									// receiver rank
					_seg_globalRecvBuffer_id,			// remote seg ID
					remote_offset,						// remote offset
					qty,								// size of data to write
					notifyID,							// remote notif ID
					ALLREDUCE,							// value of the notif to write
					queue,								// queue
					GASPI_BLOCK							// Gaspi block
				)
			);
		}	
		queue=(queue+1)%nbQueues;
	}
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_write_notify", t_end - t_begin);

    
    // write its own data into the result !
	t_begin = MPI_Wtime();    
    int nbElts = _seg_reduce_size / sizeof(complex);
	cilk_for (int i=0; i<nbElts; i++)
	{   
		ne[i] = _reduceFF[i];
	}
	t_end = MPI_Wtime();
	add_time_sec("GASPI_REDUCE_write_back_ne", t_end - t_begin);

	 
	// wait to receive all messages from the others
	int recvCpt = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
    int sender;
    double t_begin_loop, t_end_loop;
	
	while (recvCpt < (_wsize-1))
	{
		//methode 1 - avec GASPI_BLOCK
		while(1)
		{
			t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_seg_globalRecvBuffer_id,
					notif_offset,						// surveille les notifications depuis offset _wsize
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_REDUCE_notify_waitsome", t_end_loop - t_begin_loop);
			t_begin_loop = MPI_Wtime();
			
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_seg_globalRecvBuffer_id, 
					new_notif_id, 
					&new_notif_val
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_REDUCE_notify_reset", t_end_loop - t_begin_loop);
			t_begin_loop = MPI_Wtime();
			
			if (new_notif_val) 
				break;
		}
		
		sender = new_notif_id - notif_offset;
		// test the notification value and update counter
		t_begin_loop = MPI_Wtime();
	    if (new_notif_val == ALLREDUCE)
	    {
	        // update counter
		    recvCpt++;
		
		    // update the far field array		    		     
            int offset = nbElts * sender;
	        cilk_for (int i=0; i<nbElts; i++)
	        {   
				ne[i] = ne[i] + _globalRecvBuffer[offset + i];
	        }
        }	
        t_end_loop = MPI_Wtime();
		add_time_sec("GASPI_REDUCE_write_back_ne", t_end_loop - t_begin_loop);
	}        
}	
	
	
void Gaspi_m2l_communicator::runM2LCommunications(i64 * sendnode, int sendnode_sz, i64 * nb_send,int levcom, int nivterm, 
	i64 * endlev, i64 * frecv, i64 * recv, i64 * fsend, i64 * send, i64 * nst, i64 * nsp, i64 * fniv, i64 * codech, complex * bufsave, complex * ff)
{
	double t_begin, t_end;

	// init global send buffer
	t_begin = MPI_Wtime();
	initGlobalSendSegment(sendnode, sendnode_sz, nb_send, nivterm, levcom, fsend, send, endlev, codech, nst, nsp, bufsave, fniv, ff);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_write_global_sendSegment", t_end - t_begin);

	
	gaspi_offset_t local_offset;
    gaspi_offset_t remote_offset;
   	gaspi_notification_id_t notif_offset = _wsize * 2;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_queue_id_t queue=0;
    gaspi_size_t qty;
    
    int nbQueues = 1;
        
	/** avec gestion de la queue **/
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
	add_time_sec("GASPI_SEND_flush_queue", t_end - t_begin);
	t_begin = MPI_Wtime();
	queue = 0;
	
    // send infos    
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank ) // if not current rank
		{
			if (nb_send[i] > 0) // if something to send
			{
				local_offset  = _sendBufferIndexes[i] * sizeof(complex);
				remote_offset = _remoteBufferIndexes[i] * sizeof(complex);
				qty = nb_send[i] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_seg_globalSendBuffer_id,			// local seg ID
						local_offset,						// local offset
						i,									// receiver rank
						_seg_globalRecvBuffer_id,			// remote seg ID
						remote_offset,						// remote offset
						qty,								// size of data to write
						notifyID,							// remote notif ID
						SEND_DATA,							// value of the notif to write
						queue,								// queue
						GASPI_BLOCK							// Gaspi block
					)
				);
				
				register_write(_rank, i, qty);
			}
			else // no data to send, notify anyway
			{
				SUCCESS_OR_DIE(
					gaspi_notify(
						_seg_globalRecvBuffer_id,			// seg
						i,		 							// receiver rank
						notifyID,							// remote notif ID
						NO_DATA,							// value of the notif to write
						queue, 								// queue
						GASPI_BLOCK
					)
				);
			}
		}
		queue=(queue+1)%nbQueues;
	}

	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_write_notify", t_end - t_begin);
	t_begin = MPI_Wtime();


	// wait to receive all infos
	int recvCpt = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int sender;

	double t_begin_loop, t_end_loop;

	while (recvCpt < (_wsize-1))
	{
		
		//methode 1 - avec GASPI_BLOCK
		while(1)
		{
			t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_seg_globalRecvBuffer_id,
					notif_offset,				// surveille les notifications depuis 0
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_SEND_notify_waitsome", t_end_loop - t_begin_loop);
			t_begin_loop = MPI_Wtime();
			
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_seg_globalRecvBuffer_id, 
					new_notif_id, 
					&new_notif_val
				)
			);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_SEND_notify_reset", t_end_loop - t_begin_loop);
			
			if (new_notif_val) 
				break;
		}
		
		// test the notification value and update counter
		if (new_notif_val == SEND_DATA || new_notif_val == NO_DATA)
			recvCpt++;

		// update the far field array
		if (new_notif_val == SEND_DATA)
		{
			sender = new_notif_id - notif_offset;
			t_begin_loop = MPI_Wtime();			
			updateFarFields(sender, levcom, nivterm, endlev, frecv, recv, nst, nsp, fniv, ff);
			t_end_loop = MPI_Wtime();
			add_time_sec("GASPI_SEND_write_back_ff", t_end_loop - t_begin_loop);
		}
	}
}

void Gaspi_m2l_communicator::updateFarFields(int src, int levcom, int nivterm, 
	i64 * endlev, i64 * frecv, i64 * recv, i64 * nst, i64 * nsp, i64 * fniv, complex * ff)
{
	int indexToC = -1;
	int k = nivterm;
	
	int srcF=src+1;
	int qp = _globalRecvBufIdxPerRank[src]-1; // se mettre 1 case avant l'index à lire
	
	for (int j=frecv[srcF+1+indexToC]-1; j>=frecv[srcF+indexToC]; j--)
	{
		int jc = j + indexToC;
		int p = recv[jc];
		bool todo = true;
		
		if ((levcom<=3 || (levcom == 4 && nivterm>8))  && (p < endlev[levcom + indexToC]))
			todo = false;
		if ((levcom>4  || (levcom == 4 && nivterm<=8)) && (p <= endlev[levcom-1 + indexToC]))
			todo = false;
		
		if (todo)
		{
			// go to the right octree level
			while(p <= endlev[k-1 + indexToC])
				k--;
			
			// update far fields
			int pp = fniv[k+1 + indexToC]+(endlev[k + indexToC]-p)*nst[k + indexToC]*nsp[k + indexToC];
			for (int st=0; st<nst[k + indexToC]; st++)
			{
				for (int sp=0; sp<nsp[k + indexToC]; sp++)
				{
					pp++;
					qp++;
					ff[pp + indexToC] = ff[pp + indexToC] + _globalRecvBuffer[qp];
				}
			}
		}
	}
}


void construct_m2l_communicator(i64 * nb_send, int nb_send_sz, 
							 i64 * nb_recv, int nb_recv_sz,
							 i64 * sendnode, int sendnode_sz,
							 i64 * recvnode, int recvnode_sz,
							 complex * ff, complex * ne, int allreduce_sz, 
                             Gaspi_m2l_communicator *& gCommM2L)
{
    // Class Constructor
    gCommM2L = new Gaspi_m2l_communicator(nb_send, nb_send_sz, 	nb_recv, nb_recv_sz,
		sendnode, sendnode_sz, recvnode, recvnode_sz, ff, ne, allreduce_sz);
}

/** CLEM GASPI TOOLS **/
void print_gaspi_config()
{
	gaspi_config_t config; 
	gaspi_config_get(&config);
	
	cout << "logger=" << config.logger << std::endl;
	cout << "sn_port=" << config.sn_port << std::endl;
	cout << "net_info=" << config.net_info << std::endl;
	cout << "netdev_id=" << config.netdev_id << std::endl;
	cout << "mtu=" << config.mtu << std::endl;
	cout << "port_check=" << config.port_check << std::endl;
	cout << "user_net=" << config.user_net << std::endl;
	cout << "network=" << config.network << std::endl;
	cout << "queue_depth=" << config.queue_depth << std::endl;
	cout << "queue_num=" << config.queue_num << std::endl;
	cout << "group_max=" << config.group_max << std::endl;
	cout << "segment_max=" << config.segment_max << std::endl;
	cout << "transfer_size_max=" << config.transfer_size_max << std::endl;
	cout << "notification_num=" << config.notification_num << std::endl;
	cout << "passive_queue_size_max=" << config.passive_queue_size_max << std::endl;
	cout << "passive_transfer_size_max=" << config.passive_transfer_size_max << std::endl;
	cout << "allreduce_buf_size=" << config.allreduce_buf_size << std::endl;
	cout << "allreduce_elem_max=" << config.allreduce_elem_max << std::endl;
	cout << "build_infrastructure=" << config.build_infrastructure << std::endl;
}
