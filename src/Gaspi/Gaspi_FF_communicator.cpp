#include "Gaspi_FF_communicator.hpp"
using namespace std;

Gaspi_FF_communicator::Gaspi_FF_communicator(
	i64 * nb_send, int nb_send_sz, i64 * nb_recv, int nb_recv_sz, i64 * sendnode, int sendnode_sz,
	i64 * recvnode, int recvnode_sz, int nivterm, int levcom, i64 * fsend, i64 * send, i64 * frecv, 
	i64 * recv, i64 * nst, i64 * nsp, i64 * fniv, i64 * endlev, i64 * codech, int includeLevcom)
{
	// update class attributes
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);
	_nivterm = nivterm;
    _levcom = levcom;
	_fniv = fniv;
    _fsend = fsend;
    _send = send;
    _frecv = frecv;
    _recv = recv;
    _nst = nst;
    _nsp = nsp;
	_codech = codech;
	_nb_send = nb_send;
	_nb_send_sz = nb_send_sz;
	_nb_recv = nb_recv;
	_nb_recv_sz = nb_recv_sz;
	_sendnode = sendnode;
	_sendnode_sz = sendnode_sz;
	_recvnode = recvnode;
	_recvnode_sz = recvnode_sz;
	_endlev = endlev;
	_incLevcom = includeLevcom;

	// update segments ids
	_seg_globalRecvBuffer_id 	    = 7; // Receive Buffer 
	_seg_globalSendBuffer_id 	    = 8; // Send Buffer 
    _seg_remoteBufferIndexes_id     = 9; // Index where to write on other ranks	
    _seg_globalRecvBufIdxPerRank_id = 10; // Index in global recv buffer, per RANK
										  // temporary, only for initialization computations 
	// allocate array
	_sendBufferIndexes = new int[_wsize]();
	_offsetKeeper = new int[_wsize]();

	// create gaspi segments, and initialize them
	create_globalRecvBuffer(nb_recv, nb_recv_sz);
	create_remoteBufferIndexes(recvnode, recvnode_sz, nb_recv);
	create_globalSendBuffer(nb_send, nb_send_sz);
	init_sendBufferIndexes(sendnode, sendnode_sz, nb_send);
	init_expectPerSrcAndLevel();
	
	_nbQueues = 1;
}

/* -------------------------------------------------------------------- */
/*               CREATION DE SEGMENTS  et INITIALISATIONS               */
/* -------------------------------------------------------------------- */
void construct_m2l_communicator(i64 * nb_send, int nb_send_sz, i64 * nb_recv, int nb_recv_sz,
	i64 * sendnode, int sendnode_sz, i64 * recvnode, int recvnode_sz, int nivterm, int levcom,
	i64 * fsend, i64 * send, i64 * frecv, i64 * recv, i64 * nst, i64 * nsp, i64 * fniv, i64 * endlev, 
	i64 * codech, int includeLevcom, Gaspi_FF_communicator *& gCommFF)
{
    // Class Constructor
    gCommFF = new Gaspi_FF_communicator(nb_send, nb_send_sz, nb_recv, nb_recv_sz,
		sendnode, sendnode_sz, recvnode, recvnode_sz, nivterm, levcom,
		fsend, send, frecv, recv, nst, nsp, fniv, endlev, codech, includeLevcom);
}

void Gaspi_FF_communicator::create_globalRecvBuffer(i64 * nb_recv, int nb_recv_sz)
{
	int globalRecvBufferSize = 0;
	for (int i=0; i<nb_recv_sz; i++)
	    globalRecvBufferSize = globalRecvBufferSize + nb_recv[i];   

 	// create segment
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_globalRecvBuffer_id,
			globalRecvBufferSize * sizeof(complex),
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

void Gaspi_FF_communicator::create_remoteBufferIndexes(i64 * recvnode, int recvnode_sz, i64 * nb_recv)
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
    
    // initialize segment with values -> put in rank order
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
				cerr << "[Create_remoteBufferIndexes] Fortran recvnode array should not contain the current rank !"; 
				exit(0);
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
    
    flush_queues(_nbQueues);
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

void Gaspi_FF_communicator::create_globalSendBuffer(i64 * nb_send, int nb_send_sz)
{
	// compute size
	int globalSendBufferSize = 0;
	for (int i=0; i<nb_send_sz; i++)
	    globalSendBufferSize = globalSendBufferSize + nb_send[i];
  
 	// update segment size
 	_seg_globalSendBuffer_size = globalSendBufferSize * sizeof(complex);

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

void Gaspi_FF_communicator::init_sendBufferIndexes(i64 * sendnode, int sendnode_sz, i64 * nb_send)
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

void Gaspi_FF_communicator::init_expectPerSrcAndLevel()
{
	int indexToC = -1;
	
	// alloc
	_expectPerSrcAndLevel = new int* [_wsize]();
	for (int i=0; i<(_wsize); i++)
		_expectPerSrcAndLevel[i] = new int[_nivterm]();
	
	// fill
	for (int src=0; src<_wsize; src++)
	{
		if (src != _rank)
		{
			// box range to recv, FROM RECV ARRAY
			int firstBoxToRecvIdx= _frecv[src] + indexToC;
			int lastBoxToRecvIdx = _frecv[src+1]-1 + indexToC;
							 
			// for each box, find the corresponding level
			for (int k = firstBoxToRecvIdx; k<=lastBoxToRecvIdx; k++)
			{
				int cellID = _recv[k] + indexToC;					
				
				// go through octree levels
				int found = 0;
				int level;
				if(_incLevcom)
					level = _levcom + indexToC;
				else
					level = _levcom + 1 + indexToC;
				
				while ( (level <= _nivterm + indexToC) && (!found))
				{
					// test if box belongs to level 
					if ( (cellID > _endlev[level-1]+indexToC) && (cellID <= _endlev[level]+indexToC) )
					{
						found = 1;
						_expectPerSrcAndLevel[src][level]++;
					}
					level++;
				}
				if (!found)
				{
					cout << "ERROR [Gaspi_FF_communicator::init_expectPerSrcAndLevel]: cell level has not been identified !" << endl;
					exit(-1);
				}
			}
		}
	}
}

/* ---------------------------------------------- */
/*               GASPI BULK VERSION               */
/* ---------------------------------------------- */
void Gaspi_FF_communicator::exchangeFFBulk(complex * bufsave, complex * ff)
{
	double t_begin, t_end, accumul;
	accumul = 0;

	// init global send buffer
	t_begin = MPI_Wtime();
	initGlobalSendSegment(bufsave, ff);
	t_end = MPI_Wtime();
	add_time_sec("FF_write_to_buffer", t_end - t_begin);

	
	gaspi_offset_t local_offset;
    gaspi_offset_t remote_offset;
   	gaspi_notification_id_t notif_offset = _wsize * 2;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_queue_id_t queue=0;
    gaspi_size_t qty;
    
	/** avec gestion de la queue **/
	t_begin = MPI_Wtime();
	flush_queues(_nbQueues);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_flush_queue", t_end - t_begin);
	accumul += (t_end - t_begin);
	
	t_begin = MPI_Wtime();
    // send infos    
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank ) // if not current rank
		{
			if (_nb_send[i] > 0) // if something to send
			{
				local_offset  = _sendBufferIndexes[i] * sizeof(complex);
				remote_offset = _remoteBufferIndexes[i] * sizeof(complex);
				qty = _nb_send[i] * sizeof(complex);
				
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
						0,									// queue
						GASPI_BLOCK							// Gaspi block
					)
				);
			}
			else // no data to send, notify anyway
			{
				SUCCESS_OR_DIE(
					gaspi_notify(
						_seg_globalRecvBuffer_id,			// seg
						i,		 							// receiver rank
						notifyID,							// remote notif ID
						NO_DATA,							// value of the notif to write
						0, 									// queue
						GASPI_BLOCK
					)
				);
			}
		}
	}
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_write_notify", t_end - t_begin);
	accumul += (t_end - t_begin);	
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
			accumul += (t_end_loop - t_begin_loop);

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
			accumul += (t_end_loop - t_begin_loop);			
			
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
			updateFarFields(sender, ff);
			t_end_loop = MPI_Wtime();
			add_time_sec("FF_read_from_buffer", t_end_loop - t_begin_loop);
		}
	}
	add_time_sec("GASPI_FF_sendrecv", accumul);
}

void Gaspi_FF_communicator::initGlobalSendSegment(complex * bufsave, complex * ff)
{
	int indexToC = -1;
	static int first = 1;
	
	// Global send buffer indexes
	int * globalSendBufferIndex = nullptr;
	globalSendBufferIndex = new int[_sendnode_sz]();
	
	for (int i=1; i<_sendnode_sz; i++)
	{
		if (_sendnode[i-1] == 0)
			globalSendBufferIndex[i] = globalSendBufferIndex[i-1];
	    else
			globalSendBufferIndex[i] = globalSendBufferIndex[i-1] + _nb_send[_sendnode[i-1]+indexToC];
	}
	
	// Fill Global send Buffer
	for (int i=0; i<_sendnode_sz; i++)
	{
		int iDest = _sendnode[i];
		int k = _nivterm;
		int q = globalSendBufferIndex[i]-1;
		
		if (iDest > 0)
		{
			// pour toutes les cellules à échanger avec iDest
			// de la dernière cellule à envoyer -> à la première
			for (int j=_fsend[iDest+1+indexToC]-1; j>=_fsend[iDest+indexToC]; j--)
			{
				int jc = j + indexToC;
								
				bool todo = true;
				int p;
				
				// l = numero de la cellule à envoyer
				int l = _send[jc];
				if ((_levcom<=3 || (_levcom == 4 && _nivterm>8))  && (l < _endlev[_levcom + indexToC]))
					todo = false;
				if ((_levcom>4  || (_levcom == 4 && _nivterm<=8)) && (l <= _endlev[_levcom-1 + indexToC]))
					todo = false;
					
				if (todo)
				{
					// place k pour que k = niveau de la cellule à envoyer
					while(l <= _endlev[k-1 + indexToC])
						k--;
					
					// cas 1 : 
					// si l est commun à +ieurs noeuds
					// communique les données sauvegardées dans BUFSAVE
					if (_codech[l + indexToC] > 1)
					{
						p=_codech[l + indexToC]-2;
						for (int st=0; st<_nst[k + indexToC]; st++)
						{
							for (int sp=0; sp<_nsp[k + indexToC]; sp++)
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
						// calcule le bon index pour p (= offset)
						p=_fniv[k+1 + indexToC]+(_endlev[k + indexToC]-l)*_nst[k + indexToC]*_nsp[k + indexToC];
						
						// communique  nst(niveau) x nsp(niveau) complexes depuis FF
						for (int st=0; st<_nst[k + indexToC]; st++)
						{
							for (int sp=0; sp<_nsp[k + indexToC];  sp++)
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
	first = 0;
}

void Gaspi_FF_communicator::updateFarFields(int src, complex * ff)
{		
	int indexToC = -1;
	int k = _nivterm;
	int srcF=src+1;
	
	int qp = _globalRecvBufIdxPerRank[src]-1; // se mettre 1 case avant l'index à lire
	
	for (int j=_frecv[srcF+1+indexToC]-1; j>=_frecv[srcF+indexToC]; j--)
	{
		int jc = j + indexToC;
		int p = _recv[jc];
		bool todo = true;
		
		if ((_levcom<=3 || (_levcom == 4 && _nivterm>8))  && (p < _endlev[_levcom + indexToC]))
			todo = false;
		if ((_levcom>4  || (_levcom == 4 && _nivterm<=8)) && (p <= _endlev[_levcom-1 + indexToC]))
			todo = false;
		
		if (todo)
		{
			// go to the right octree level
			while(p <= _endlev[k-1 + indexToC])
				k--;
			
			// update far fields
			int pp = _fniv[k+1 + indexToC]+(_endlev[k + indexToC]-p)*_nst[k + indexToC]*_nsp[k + indexToC];
			
			for (int st=0; st<_nst[k + indexToC]; st++)
			{
				for (int sp=0; sp<_nsp[k + indexToC]; sp++)
				{
					pp++;
					qp++;
					ff[pp + indexToC] = ff[pp + indexToC] + _globalRecvBuffer[qp];
				}
			}
		}
	}
}

/* ------------------------------------------------- */
/*               GASPI OVERLAP VERSION               */
/* ------------------------------------------------- */

void Gaspi_FF_communicator::send_ff_level(int level, complex * ff)
{
	int indexToC = -1;
	double t_begin, t_end, accumul;
	accumul = 0;
	
	// flush queue
	t_begin = MPI_Wtime();
	flush_queues(_nbQueues);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_wait_queue", t_end - t_begin);
	accumul = accumul + (t_end - t_begin);	
	
	// prepare _globalSendBuffer
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{
			t_begin = MPI_Wtime();
			// box range to send, FROM SEND ARRAY
			int firstBoxToSendIDX= _fsend[dest] + indexToC;
			int lastBoxToSendIDX = _fsend[dest+1]-1 + indexToC;
			
			// cell range from the current level
			int beginlevel = _endlev[level-1]+1 + indexToC; 
			int endlevel   = _endlev[level] + indexToC;
			
			int count = 0;
			int q = _sendBufferIndexes[dest] - 1 + _offsetKeeper[dest];
			
			// remplissage du buffer d'envoi, par ordre décroissant
			for (int k = lastBoxToSendIDX; k>=firstBoxToSendIDX; k--)
			{
				int cellID = _send[k] + indexToC;
				
				if ( (cellID >= beginlevel) && (cellID <= endlevel))	/* todo : stop earlier, cells are ranged in ascending order*/
				{
					count++;
					
					// calcule le bon index pour p (= offset, dans le tableau ff)
					int p=_fniv[level+1]+(_endlev[level]+indexToC-cellID)*_nst[level]*_nsp[level];
					p = p + indexToC;
					
					// ecrit dans globalSendBuffer les nst(niveau) x nsp(niveau) complexes, de cellID, depuis FF
					for (int st=0; st<_nst[level]; st++)
					{
						for (int sp=0; sp<_nsp[level];  sp++)
						{
							p++;
							q++;
							_globalSendBuffer[q]=ff[p];
						}
					}
				}
			}
			
			t_end = MPI_Wtime();
			add_time_sec("FF_write_to_buffer", t_end - t_begin);
			t_begin = MPI_Wtime();
			
			// if something to send
			if (count > 0)
			{
				// local offset
				gaspi_offset_t local_dest_offset = _sendBufferIndexes[dest] * sizeof(complex);
				gaspi_offset_t level_offset = _offsetKeeper[dest] * sizeof(complex);
				gaspi_offset_t local_offset = local_dest_offset + level_offset;
				gaspi_queue_id_t queue=0;
				
				// remote offset
				gaspi_offset_t remote_sender_offset = _remoteBufferIndexes[dest] * sizeof(complex);
				gaspi_offset_t remote_offset = remote_sender_offset + level_offset;


				// update offset, per destinatary
				_offsetKeeper[dest] += count * _nst[level] * _nsp[level];

				int rankMultiple = level;
				gaspi_notification_id_t notif_offset = _wsize * rankMultiple;
				gaspi_notification_id_t notifyID = notif_offset + _rank;
				gaspi_size_t qty= count * _nst[level] * _nsp[level] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_seg_globalSendBuffer_id,			// local seg ID
						local_offset,						// local offset
						dest,								// receiver rank
						_seg_globalRecvBuffer_id,			// remote seg ID
						remote_offset,						// remote offset
						qty,								// size of data to write
						notifyID,							// remote notif ID
						level,								// value of the notif to write
						queue,								// queue
						GASPI_BLOCK							// Gaspi block
					)
				);
			}
			t_end = MPI_Wtime();
			add_time_sec("GASPI_SEND_write_notify", t_end - t_begin);
			accumul = accumul + (t_end - t_begin);
		}

		// cas include levcom
		if(_incLevcom)
		{
			if (level == _levcom + indexToC)
				_offsetKeeper[dest] = 0;
		}
		else // cas allreduce sur levcom
		{
			if (level == _levcom + 1 + indexToC)
				_offsetKeeper[dest] = 0;
		}
	}

	add_time_sec("GASPI_FF_sendrecv", accumul);
}

void Gaspi_FF_communicator::recv_ff_level(int level, complex * ff)
{
	// wait to receive all infos
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
		nbRecvExpected += (_expectPerSrcAndLevel[i][level]>0);
	
	int rankMultiple = level;
	gaspi_notification_id_t notif_offset = _wsize * rankMultiple;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int recvCpt = 0;
	int sender;
	double t_begin, t_end;

	while (recvCpt < nbRecvExpected)
	{
		t_begin = MPI_Wtime();
		//methode 1 - avec GASPI_BLOCK
		// ne surveille que les notifs du level
		while(1)
		{
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_seg_globalRecvBuffer_id,
					notif_offset,				// surveille les notifications depuis 0
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_seg_globalRecvBuffer_id, 
					new_notif_id, 
					&new_notif_val
				)
			);
			
			if (new_notif_val) 
				break;
		}
		
		t_end = MPI_Wtime();
		add_time_sec("GASPI_RECV_notify_waitsome", t_end - t_begin);
		add_time_sec("GASPI_FF_sendrecv", t_end - t_begin);
		

		// test the notification value and compare with array (I know, useless)
		t_begin = MPI_Wtime();
		if (new_notif_val)
		{
			recvCpt++;
			sender = new_notif_id - notif_offset;
			updateFarFields(sender, level, ff);
		}
		t_end = MPI_Wtime();
		add_time_sec("FF_read_from_buffer", t_end - t_begin);
	}
}


/* For the Overlapping, level by level version */
void Gaspi_FF_communicator::updateFarFields(int src, int level, complex * ff)
{		
	int indexToC = -1;
	int k = level + 1;
	
	// calcul de l'offset	
	int levelOffset = 0;
	for (int i=_nivterm+indexToC; i>level; i--)
	{
		levelOffset += _nst[i]*_nsp[i]*_expectPerSrcAndLevel[src][i];
	}
	int q = _globalRecvBufIdxPerRank[src]-1; // se mettre 1 case avant l'index à lire
	q += levelOffset; 	

	// box range to recv, FROM RECV ARRAY
	int firstBoxToRecvIdx= _frecv[src] + indexToC;
	int lastBoxToRecvIdx = _frecv[src+1]-1 + indexToC;

    for (int j=lastBoxToRecvIdx; j>=firstBoxToRecvIdx; j--)
	{
		int cellID = _recv[j] + indexToC;
	
		// test if box belongs to level 
		if ( (cellID > _endlev[level-1]+indexToC) && (cellID <= _endlev[level]+indexToC) )
		{
			// update far fields
			int p = _fniv[level+1]+(_endlev[level]+indexToC-cellID)*_nst[level]*_nsp[level];
			p = p + indexToC;
			
			for (int st=0; st<_nst[level]; st++)
			{
				for (int sp=0; sp<_nsp[level]; sp++)
				{
					p++;
					q++;
					ff[p] = ff[p] + _globalRecvBuffer[q];
				}
			}
		}
	}
}




