#include "Gaspi_FF_communicator.hpp"
using namespace std;

/* ---------------------------------------------- */
/*               GASPI BULK VERSION               */
/* ---------------------------------------------- */

void Gaspi_FF_communicator::exchangeFFBulk(complex * bufsave, complex * ff, int iOct)
{
	//cout << "GASPI BULK" << endl;
	int indexToC = -1;
	double t_begin, t_end, accumul;
	accumul = 0;
	int octree_offset = iOct * _wsize;
	
	// flush queue
	t_begin = MPI_Wtime();
	flush_queues(_nbQueues);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_wait_queue", t_end - t_begin);
	accumul = accumul + (t_end - t_begin);	
	int p2;
	// PREPARATION _FF_sendBuffer
	for (int dest=0; dest<_wsize; dest++)
	{
		int niv = _nivterm[iOct];
		if (dest != _rank)
		{
			t_begin = MPI_Wtime();
			// box range to send, FROM SEND ARRAY
			int firstBoxToSendIDX= _fsend[iOct][dest] + indexToC;
			int lastBoxToSendIDX = _fsend[iOct][dest+1]-1 + indexToC;
			
			int q = _FF_sendLocalOffsets_byOctDest[octree_offset + dest] - 1 /*+ _FF_sendLocalOffsets_keeper[iOct][dest]*/;
			cout << std::scientific;
			// remplissage du buffer d'envoi, par ordre décroissant
			for (int k = lastBoxToSendIDX; k>=firstBoxToSendIDX; k--)
			{
				bool todo = true;
				int p;
				int cellID = _send[iOct][k] + indexToC;

				if ((_levcom[iOct]<=3 || (_levcom[iOct] == 4 && _nivterm[iOct]>8))  && (cellID < _endlev[iOct][_levcom[iOct] + indexToC]))
					todo = false;
				if ((_levcom[iOct]>4  || (_levcom[iOct] == 4 && _nivterm[iOct]<=8)) && (cellID <= _endlev[iOct][_levcom[iOct]-1 + indexToC]))
					todo = false;

				if (todo)
				{
					// place k pour que k = niveau de la cellule à envoyer
					while(cellID <= _endlev[iOct][niv-1 + indexToC])
						niv--;
					
					// cas 1 : 
					// si l est commun à +ieurs noeuds
					// communique les données sauvegardées dans BUFSAVE
					if (_codech[iOct][cellID + indexToC] > 1)
					{
						p=_codech[iOct][cellID + indexToC]-2;
						for (int st=0; st<_nst[iOct][niv + indexToC]; st++)
						{
							for (int sp=0; sp<_nsp[iOct][niv + indexToC]; sp++)
							{
								p++;
								q++;
								_FF_sendBuffer[q]=bufsave[p + indexToC];
							}
						}
					}
					else
					{
						// calcule le bon index pour p (= offset)
						p=_fniv[iOct][niv+1 + indexToC]+(_endlev[iOct][niv + indexToC]-cellID + indexToC)*_nst[iOct][niv + indexToC]*_nsp[iOct][niv + indexToC];						

						// communique  nst(niveau) x nsp(niveau) complexes depuis FF
						for (int st=0; st<_nst[iOct][niv + indexToC]; st++)
						{
							for (int sp=0; sp<_nsp[iOct][niv + indexToC];  sp++)
							{
								p++;
								q++;
								_FF_sendBuffer[q]=ff[p + indexToC];
							}
						}
					}
				}
			}
		}
	}
	t_end = MPI_Wtime();
	add_time_sec("FF_write_to_buffer", t_end - t_begin);
	t_begin = MPI_Wtime();
			
	// ENVOIS
	gaspi_offset_t local_offset;
    gaspi_offset_t remote_offset;
   	gaspi_notification_id_t notif_offset = _wsize * 2;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
    gaspi_queue_id_t queue=0;
    gaspi_size_t qty;
    
	t_begin = MPI_Wtime();
	flush_queues(_nbQueues);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_flush_queue", t_end - t_begin);
	accumul += (t_end - t_begin);
	
	t_begin = MPI_Wtime();   
    for (int dest=0; dest<_wsize; dest++)
    {
		if ( dest != _rank ) // if not current rank
		{
			if (_nb_send[iOct][dest] > 0) // if something to send
			{
				local_offset  = _FF_sendLocalOffsets_byOctDest[octree_offset + dest] * sizeof(complex);
				remote_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
				qty = _nb_send[iOct][dest] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_FF_sendBuf_seg_id,					// local seg ID
						local_offset,						// local offset
						dest,									// receiver rank
						_FF_recvBuf_seg_id,					// remote seg ID
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
						_FF_recvBuf_seg_id,			// seg
						dest,		 							// receiver rank
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
	
	// RECEIVE
	int recvCpt = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int sender;
	double t_begin_loop, t_end_loop;
	
	while (recvCpt < (_wsize-1))
	{
		while(1)
		{
			t_begin_loop = MPI_Wtime();
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_FF_recvBuf_seg_id,
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
					_FF_recvBuf_seg_id, 
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
			updateFarFields(sender, ff, iOct);
			t_end_loop = MPI_Wtime();
			add_time_sec("FF_read_from_buffer", t_end_loop - t_begin_loop);
		}
	}
}

void Gaspi_FF_communicator::updateFarFields(int src, complex * ff, int iOct)
{		
	int indexToC = -1;
	int k = _nivterm[iOct];
	int srcF=src+1;
	
	int octree_offset = iOct * _wsize;
	int qp =_FF_recvOffsets[octree_offset + src]-1; // se mettre 1 case avant l'index à lire
		
	for (int j=_frecv[iOct][srcF+1+indexToC]-1; j>=_frecv[iOct][srcF+indexToC]; j--)
	{
		int jc = j + indexToC;
		int p = _recv[iOct][jc];
		bool todo = true;
		
		if ((_levcom[iOct]<=3 || (_levcom[iOct] == 4 && _nivterm[iOct]>8))  && (p < _endlev[iOct][_levcom[iOct] + indexToC]))
			todo = false;
		if ((_levcom[iOct]>4  || (_levcom[iOct] == 4 && _nivterm[iOct]<=8)) && (p <= _endlev[iOct][_levcom[iOct]-1 + indexToC]))
			todo = false;
		
		if (todo)
		{
			// go to the right octree level
			while(p <= _endlev[iOct][k-1 + indexToC])
				k--;
			
			// update far fields
			int pp = _fniv[iOct][k+1 + indexToC]+(_endlev[iOct][k + indexToC]-p)*_nst[iOct][k + indexToC]*_nsp[iOct][k + indexToC];
			
			for (int st=0; st<_nst[iOct][k + indexToC]; st++)
			{
				for (int sp=0; sp<_nsp[iOct][k + indexToC]; sp++)
				{
					pp++;
					qp++;
					ff[pp + indexToC] = ff[pp + indexToC] + _FF_recvBuffer[qp];
				}
			}
		}
	}
}

/* ------------------------------------------------- */
/*               GASPI MULTIMAT VERSION              */
/* ------------------------------------------------- */
Gaspi_FF_communicator::Gaspi_FF_communicator(int max_send_terms, int max_recv_terms, int max_send_nodes, int max_recv_nodes, int incLevcom, int nbOct)
{
	// update class attributes
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);

	// update segments ids
	_FF_recvBuf_seg_id = 7; 	// Receive Buffer 
	_FF_sendBuf_seg_id = 8; 	// Send Buffer 
    _FF_recvOffsets_seg_id = 9; 	// Index where to write on other ranks	
	_FF_sendRemoteOffsets_seg_id = 10; // Index in global recv buffer, per RANK
	_Infos_sendbuf_seg_id = 11;			// 11
    _Infos_recvbuf_seg_id = 12;
    _Infos_recvOffsets_seg_id = 13; 		// 13 - local @
    _Infos_sendRemoteOffsets_seg_id = 14; 
    
	// update scalar class attributes
	_nbQueues = 1;
	_incLevcom = incLevcom;
	_nbOct = nbOct;
		
	// allocate array class attributes
	alloc_attributes();
	
	// create gaspi segments, and initialize them
	create_segments(max_send_terms, max_recv_terms, max_send_nodes, max_recv_nodes);
}

void init_gaspi_ff(int max_send_terms, int max_recv_terms, int max_send_nodes, int max_recv_nodes, int nbMat, int incLevcom, Gaspi_FF_communicator *& gCommFF)
{
    gCommFF = new Gaspi_FF_communicator(max_send_terms, max_recv_terms, max_send_nodes, max_recv_nodes, incLevcom, nbMat);
}

void Gaspi_FF_communicator::alloc_attributes()
{
	_FF_sendLocalOffsets_keeper = new int*[_nbOct]();			//write into segment
	for (int i=0; i<_nbOct; i++)
		_FF_sendLocalOffsets_keeper[i] = new int[_wsize]();
		
	_FF_sendRemoteOffsets_keeper = new int*[_nbOct]();		//send from segment
	for (int i=0; i<_nbOct; i++)
		_FF_sendRemoteOffsets_keeper[i] = new int[_wsize]();		

	_FF_sendLocalOffsets_byOctDest = new int[_wsize * _nbOct](); // Local SEND offsets x nb of octrees (=nb materials)

	_Expect = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_Expect[i] = new int*[_wsize]();

	_start_send = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_start_send[i] = new int*[_wsize]();

	_stop_send = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_stop_send[i] = new int*[_wsize]();

	_count_send = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_count_send[i] = new int*[_wsize]();

	_start_recv = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_start_recv[i] = new int*[_wsize]();

	_stop_recv = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_stop_recv[i] = new int*[_wsize]();
	
/*	_send_infos_ptr = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_send_infos_ptr[i] = new int[_wsize]();	*/

/*	_recv_infos_ptr = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_recv_infos_ptr[i] = new int[_wsize]();	*/

	_Infos_sendLocalOffsets_keeper = new int*[_nbOct]();			// write into segment
	for (int i=0; i<_nbOct; i++)
		_Infos_sendLocalOffsets_keeper[i] = new int[_wsize]();
		
	_Infos_sendRemoteOffsets_keeper = new int*[_nbOct]();			// send from segment
	for (int i=0; i<_nbOct; i++)
		_Infos_sendRemoteOffsets_keeper[i] = new int[_wsize]();		
	
	_mutexArray = new pthread_mutex_t [_nbOct*_wsize];
	for (int i=0; i<_nbOct*_wsize; i++)
		_mutexArray[i] = PTHREAD_MUTEX_INITIALIZER;

	// Class attributes
	_nivterm = new int[_nbOct]();				
    _levcom = new int[_nbOct]();	
    _fniv = new i64*[_nbOct]();	
    _fsend = new i64*[_nbOct]();	
    _send = new i64*[_nbOct]();	
    _frecv = new i64*[_nbOct]();
    _recv = new i64*[_nbOct]();
    _nst = new i64*[_nbOct]();			
    _nsp = new i64*[_nbOct]();		
    _endlev = new i64*[_nbOct]();
	_codech = new i64*[_nbOct]();
	_nb_send = new i64*[_nbOct]();
	_nb_recv = new i64*[_nbOct]();
	_sendnode = new i64*[_nbOct]();
	_recvnode = new i64*[_nbOct]();
	_nb_send_sz = new int[_nbOct]();
	_nb_recv_sz = new int[_nbOct]();
	_sendnode_sz = new int[_nbOct]();
	_recvnode_sz = new int[_nbOct]();
	//printf("[%d] alloc_attributes for nbOct : %d --> OK\n", _rank, _nbOct);
}

void Gaspi_FF_communicator::create_segments(int max_send_terms, int max_recv_terms, int max_send_nodes, int max_recv_nodes)
{
	/***********/
	/* FF DATA */ 
	/***********/
	// segments reutilises  pour chaque octree, car volumineux
	
	// send 
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_FF_sendBuf_seg_id,
			max_send_terms * sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);	
	gaspi_segment_ptr(_FF_sendBuf_seg_id, &_FF_sendBuf_seg_ptr);
	_FF_sendBuffer = (complex *)_FF_sendBuf_seg_ptr;
 	
 	// recv 
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_FF_recvBuf_seg_id,
			max_recv_terms * sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_FF_recvBuf_seg_id, &_FF_recvBuf_seg_ptr);
	_FF_recvBuffer = (complex *)_FF_recvBuf_seg_ptr;

	/************/
	/* FF INFOS */
	/************/
	// segments reutilises  pour chaque octree
	// parallelisable par octree

	// send
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_Infos_sendbuf_seg_id,
			max_send_nodes * sizeof(int),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_Infos_sendbuf_seg_id, &_Infos_sendbuf_seg_ptr);
	_Infos_sendbuffer = (int *)_Infos_sendbuf_seg_ptr;
	cout << "SEG SEND INFOS SIZE (bytes) : " << max_send_nodes * sizeof(int) << endl;
	
	// recv
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_Infos_recvbuf_seg_id,
			max_recv_nodes * sizeof(int),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_Infos_recvbuf_seg_id, &_Infos_recvbuf_seg_ptr);
	_Infos_recvbuffer = (int *)_Infos_recvbuf_seg_ptr;
	cout << "SEG RECV INFOS SIZE (bytes) : " << max_recv_nodes * sizeof(int) << endl;


	/*****************/
	/* FF   OFFSETS  */
	/*****************/

	// Local DATA Recv Offsets
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_FF_recvOffsets_seg_id,
			_nbOct * _wsize * sizeof(int),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_FF_recvOffsets_seg_id, &_Infos_recvOffsets_seg_ptr);
	_FF_recvOffsets = (int *)_Infos_recvOffsets_seg_ptr;
	
	// Remote DATA Send offsets
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_FF_sendRemoteOffsets_seg_id, 
			_nbOct * _wsize * sizeof(int),
			GASPI_GROUP_ALL, 
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
    );	
	gaspi_segment_ptr(_FF_sendRemoteOffsets_seg_id, &_FF_sendRemoteOffsets_seg_ptr);
	_FF_sendRemoteOffsets = (int *)_FF_sendRemoteOffsets_seg_ptr;	

	/*****************/
	/* INFO OFFSETS  */
	/*****************/
	
	// Local INFOS Recv offsets
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_Infos_recvOffsets_seg_id,
			_nbOct * _wsize * sizeof(int),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_Infos_recvOffsets_seg_id, &_FF_recvOffsets_seg_ptr);
	_Infos_recvOffsets = (int *)_FF_recvOffsets_seg_ptr;	
		
	// Remote INFOS Send offsets
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_Infos_sendRemoteOffsets_seg_id,
			_nbOct * _wsize * sizeof(int),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_Infos_sendRemoteOffsets_seg_id, &_Infos_sendRemoteOffsets_seg_ptr);
	_Infos_sendRemoteOffsets = (int *)_Infos_sendRemoteOffsets_seg_ptr;
}

void Gaspi_FF_communicator::init_gaspi_offsets(i64 * recvnode, int recvnode_sz, i64 * sendnode, int sendnode_sz, 
		i64 * nb_recv, int nb_recv_sz, i64 * nb_send, int nb_send_sz, int iOct, int nDom,
		int nivterm, i64 * frecv, i64 * recv, int levcom, i64 * endlev, i64 * fniv, i64 * fsend, i64 * send, 
		i64 * nst, i64 * nsp, i64 * codech)		
{	
	// attributes
	fill_attributes(iOct, nivterm, levcom, fniv, fsend, send, frecv, recv, nst, nsp, 
		endlev, codech, nb_send, nb_recv, sendnode, recvnode, nb_send_sz, nb_recv_sz, sendnode_sz, recvnode_sz);
	
	// offsets
	fill_remote_send_offsets(recvnode, recvnode_sz, nb_recv, nb_recv_sz, iOct);
	fill_local_send_offsets(sendnode, sendnode_sz, nb_send, nb_send_sz, iOct);
	
	// expected data
	fill_expectations(iOct);
	fill_send_start_stop_count_send(iOct);
	fill_recv_start_stop(iOct);
	
	// send infos pointers
	//fill_send_info_ptrs(iOct);
	fill_recv_info_ptrs(iOct);
	fill_remote_info_ptrs(iOct);
}

void Gaspi_FF_communicator::fill_attributes(int iOct, int nivterm, int levcom, i64 * fniv, i64 * fsend, i64 * send, i64 * frecv, i64 * recv, i64 * nst, i64 * nsp,
	i64 * endlev, i64 * codech, i64 * nb_send, i64 * nb_recv, i64 * sendnode, i64 * recvnode, int nb_send_sz, int nb_recv_sz, int sendnode_sz, int recvnode_sz)
{
	// Class attributes
	 _nivterm[iOct] = nivterm;				
     _levcom[iOct] = levcom;	
    _fniv[iOct] = fniv;	
    _fsend[iOct] = fsend;	
    _send[iOct] = send;	
    _frecv[iOct] = frecv;
    _recv[iOct] = recv;
    _nst[iOct] = nst;			
    _nsp[iOct] = nsp;		
    _endlev[iOct] = endlev;
	_codech[iOct] = codech;
	_nb_send[iOct] = nb_send;
	_nb_recv[iOct] = nb_recv;
	_sendnode[iOct] = sendnode;
	_recvnode[iOct] = recvnode;
	_nb_send_sz[iOct] = nb_send_sz;
	_nb_recv_sz[iOct] = nb_recv_sz;
	_sendnode_sz[iOct] = sendnode_sz;
	_recvnode_sz[iOct] = recvnode_sz;
}

void Gaspi_FF_communicator::fill_remote_send_offsets(i64 * recvnode, int recvnode_sz, i64 * nb_recv, int nb_recv_sz, int iOct)
{	
	// From Fortran to C/C++
	int indexToC = -1;
	
    /* From per Round to per Rank */     
	int * RecvBufferIndexPerROUND = new int[recvnode_sz]();
	for (int i=1; i<recvnode_sz; i++)
	{
		if (recvnode[i-1] == 0)
			RecvBufferIndexPerROUND[i] = RecvBufferIndexPerROUND[i-1];
	    else
			RecvBufferIndexPerROUND[i] = RecvBufferIndexPerROUND[i-1] + nb_recv[recvnode[i-1]+indexToC]; 
    }
     
    // initialize_FF_recvoffsets[octree] segment with values -> put in rank order
    int octree_offset = iOct * _wsize;
    for (int k=0; k<recvnode_sz; k++)
    {
		if (recvnode[k] > 0)
		{
			int from = recvnode[k]-1;
			if (from != _rank)
			{
				_FF_recvOffsets[octree_offset + from] = RecvBufferIndexPerROUND[k];
			}
			else
			{
				cerr << "[Create_remoteBufferIndexes] Fortran recvnode array should not contain the current rank !"; 
				exit(0);
			}
		}
    }
	_FF_recvOffsets[octree_offset + _rank] = -1;
	
	delete [] RecvBufferIndexPerROUND;
	     
    // Data Exchange to Update the Remote Buffer Indexes segment       
    flush_queues(_nbQueues);
    int local_offset;
    int remote_offset = (octree_offset + _rank) * sizeof(int);
    
    // send infos
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank)
		{
			local_offset = (octree_offset + i) * sizeof(int);
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					_FF_recvOffsets_seg_id,					// local seg ID
					local_offset,						// local offset
					i,									// receiver rank
					_FF_sendRemoteOffsets_seg_id,				// remote seg ID
					remote_offset,						// remote offset
					sizeof(int),						// size of data to write
					octree_offset + _rank,				// remote notif ID
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
			_FF_sendRemoteOffsets_seg_id,
			octree_offset,				            // surveille les notifications depuis le domaine
			_wsize,						            // en surveille wsize
			&new_notif_id,
			GASPI_BLOCK) == GASPI_SUCCESS)
		{

			// get notification value, and reset
			gaspi_notify_reset(_FF_sendRemoteOffsets_seg_id, new_notif_id, &new_notif_val);
		
			// test notification and increase counter
			if (new_notif_val == REMOTE_ADDRESS)
			{
				cpt++;
			}
			else
			{
				cerr << "[Create_remoteBufferIndexes] Unexpected gaspi msg."; 
				exit(0);
			}
		}
	}
	
	_FF_sendRemoteOffsets[octree_offset + _rank] = -1;
}

void Gaspi_FF_communicator::fill_remote_info_ptrs(int iOct)
{	
	// From Fortran to C/C++
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
    
    // copy array into segment
    /*for (int i=0; i<_wsize; i++)
		cout << "send : " << _Infos_recvOffsets[octree_offset + i] << endl;*/
    
    // Data Exchange to Update the Remote Buffer Indexes segment       
    flush_queues(_nbQueues);
    int local_offset;
    int remote_offset = (octree_offset + _rank) * sizeof(int);
    
    // send 
    for (int i=0; i<_wsize; i++)
    {
		if ( i != _rank)
		{
			local_offset = (octree_offset + i) * sizeof(int);
			SUCCESS_OR_DIE(
				gaspi_write_notify( 
					_Infos_recvOffsets_seg_id,			// local seg ID
					local_offset,						// local offset
					i,									// receiver rank
					_Infos_sendRemoteOffsets_seg_id,			// remote seg ID
					remote_offset,						// remote offset
					sizeof(int),						// size of data to write
					octree_offset + _rank,				// remote notif ID
					REMOTE_INFO_ADDRESS,						// value of the notif to write
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
			_Infos_sendRemoteOffsets_seg_id,
			octree_offset,				            // surveille les notifications depuis le domaine
			_wsize,						            // en surveille wsize
			&new_notif_id,
			GASPI_BLOCK) == GASPI_SUCCESS)
		{

			// get notification value, and reset
			gaspi_notify_reset(_Infos_sendRemoteOffsets_seg_id, new_notif_id, &new_notif_val);
		
			// test notification and increase counter
			if (new_notif_val == REMOTE_INFO_ADDRESS)
			{
				cpt++;
			}
			else
			{
				cerr << "[Create_remoteInfoAdresses] Unexpected gaspi msg."; 
				exit(0);
			}
		}
	}
	
	_Infos_sendRemoteOffsets[octree_offset + _rank] = -1;
	/*for (int i=0; i<_wsize; i++)
		cout << _rank << " dest " << i << " octre " << iOct << " remote info pointers ---> : " << _Infos_sendRemoteOffsets[octree_offset + i]<< endl;*/
}


void Gaspi_FF_communicator::fill_local_send_offsets(i64 * sendnode, int sendnode_sz, i64 * nb_send, int nb_send_sz, int iOct)
{
	// From Fortran to C/C++
	int indexToC = -1;
	
    // SendBufferIndexes organised by communication ROUND since sendnode is organised by communication round
	int * SendBufferIndexPerROUND = new int[sendnode_sz]();
	for (int i=1; i<sendnode_sz; i++)
	{
		if (sendnode[i-1] == 0)
			SendBufferIndexPerROUND[i] = SendBufferIndexPerROUND[i-1];
	    else
			SendBufferIndexPerROUND[i] = SendBufferIndexPerROUND[i-1] + nb_send[sendnode[i-1]+indexToC]; 
    }
    
    // Prepare _FF_sendBufferIndexes organised by process RANK, and MATERIAL DOMAIN
    int octree_offset = iOct * _wsize;
    for (int k=0; k<sendnode_sz; k++)
    {
		if(sendnode[k] > 0)
		{
			int to = sendnode[k]-1;
			if (to != _rank)
			{
				_FF_sendLocalOffsets_byOctDest[octree_offset + to] = SendBufferIndexPerROUND[k];
			}
			else
			{
				cerr << "Error in Gaspi_m2l_communicator::init_dataToSendIndexes.\nFortran sendnode array should not contain the current rank !";
				exit(0);
			}
		}
    }  
	_FF_sendLocalOffsets_byOctDest[octree_offset + _rank] = -1;
		
	// delete
	delete [] SendBufferIndexPerROUND;
}

void Gaspi_FF_communicator::fill_expectations(int iOct)
{
	int indexToC = -1;

	// init last dim only
	for (int j=0; j<_wsize; j++)
	{
		_Expect[iOct][j] = new int[_nivterm[iOct]]();
	}
	
	// fill
	for (int src=0; src<_wsize; src++)
	{
		if (src != _rank)
		{
			// box range to recv, FROM RECV ARRAY
			int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
			int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;
							 
			// for each box, find the corresponding level
			for (int k = firstBoxToRecvIdx; k<=lastBoxToRecvIdx; k++)
			{
				int cellID = _recv[iOct][k] + indexToC;					
				
				// go through octree levels
				int found = 0;
				int level;
				if(_incLevcom)
					level = _levcom[iOct] + indexToC;
				else
					level = _levcom[iOct] + 1 + indexToC;
				
				while ( (level <= _nivterm[iOct] + indexToC) && (!found))
				{
					// test if box belongs to level 
					if ( (cellID > _endlev[iOct][level-1]+indexToC) && (cellID <= _endlev[iOct][level]+indexToC) )
					{
						found = 1;
						_Expect[iOct][src][level]++;
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

void Gaspi_FF_communicator::fill_send_start_stop_count_send(int iOct)
{
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	int indexToC = -1;

	// init last dim only
	for (int j=0; j<_wsize; j++)
	{
		_start_send[iOct][j] = new int[_nivterm[iOct]]();
		_stop_send[iOct][j]  = new int[_nivterm[iOct]]();
		_count_send[iOct][j] = new int[_nivterm[iOct]]();
	}
	
	// for each rank
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{
			// for each level
			for (int level=_nivterm[iOct]+indexToC; level>=_levcom[iOct]+indexToC; level--)
			{
				// box range to send, FROM SEND ARRAY
				int firstBoxToSendIDX= _fsend[iOct][dest] + indexToC;
				int lastBoxToSendIDX = _fsend[iOct][dest+1]-1 + indexToC;
				
				// cell range from the current level
				int beginlevel = _endlev[iOct][level-1]+1 + indexToC; 
				int endlevel   = _endlev[iOct][level] + indexToC;
				
				int count = 0;
				
				// --- 
				int __endlev = _endlev[iOct][level]+indexToC;
				int __endlevprec = _endlev[iOct][level-1]+indexToC;
				long * __send = _send[iOct];
				bool start = false;
				bool stop = false;
				int firstBox, lastBox;
				
				for (int j=lastBoxToSendIDX; j>=firstBoxToSendIDX; j--)
				{
					int cellID = __send[j] + indexToC;

					if (!start)
					{
						if ( (cellID > __endlevprec) && (cellID <= __endlev) )
						{
							start = true;
							lastBox = j;
							_start_send[iOct][dest][level]=j;
						}
					}
					else
					{
						if (!stop)
						{
							if ( cellID <= __endlevprec) 
							{
								stop = true;
								firstBox = j+1;
								_stop_send[iOct][dest][level] = j+1;
							}
						}
					}
				}

				if (!start)
				{
					lastBox =  -1;
					firstBox = -1;
					_start_send[iOct][dest][level]= -1;
					_stop_send[iOct][dest][level] = -1;
					_count_send[iOct][dest][level] = 0;				
				}
				
				if (start)
				{
					if(!stop)
					{
						firstBox = firstBoxToSendIDX;
						_stop_send[iOct][dest][level] = firstBoxToSendIDX;
					}
					_count_send[iOct][dest][level] = lastBox-firstBox+1;
				}
			}
		}
	}
	t_end = MPI_Wtime();
	add_time_sec("fill_send_start_stop_count", t_end - t_begin);
}

/*void Gaspi_FF_communicator::fill_send_info_ptrs(int iOct)
{
	int indexToC = -1;

	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest == 0)
		{
			_send_infos_ptr[iOct][dest] = 0;
		}
		else
		{
			int sum = 0;			
			for (int level= 0; level<_nivterm[iOct]+indexToC; level++)
			{
				sum += _count_send[iOct][dest][level];
			}	
			_send_infos_ptr[iOct][dest] = _send_infos_ptr[iOct][dest] + sum;
		}
	}
	_send_infos_ptr[iOct][_rank] = -1;
}*/

void Gaspi_FF_communicator::fill_recv_info_ptrs(int iOct)
{
	int indexToC = -1;
	int octree_offset = iOct * _wsize;

	for (int src=0; src<_wsize; src++)
	{
		if (src == 0)
		{
			//_recv_infos_ptr[iOct][src] = 0;
			_Infos_recvOffsets[octree_offset + src] = 0;
		}
		else
		{
			int sum = 0;			
			for (int level= 0; level<_nivterm[iOct]+indexToC; level++)
			{
				sum += _Expect[iOct][src][level];
				//cout << "add : " << _Expect[iOct][src][level] << endl;
			}	
			//_recv_infos_ptr[iOct][src] = _recv_infos_ptr[iOct][src] + sum;
			//cout << "writing sum into array : " << sum << endl;
			_Infos_recvOffsets[octree_offset + src] += sum;
		}
	}
	//_recv_infos_ptr[iOct][_rank] = -1;
	_Infos_recvOffsets[octree_offset + _rank] = -1;
	
	//for (int i=0; i<_wsize; i++)
		//cout << _rank << " octree " << iOct << " from " << i << " ptr : " << _Infos_recvOffsets[octree_offset + i] << endl;
}

void Gaspi_FF_communicator::fill_recv_start_stop(int iOct)
{
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	int indexToC = -1;

	// init last dim only
	for (int j=0; j<_wsize; j++)
	{
		_start_recv[iOct][j] = new int[_nivterm[iOct]]();
		_stop_recv[iOct][j] = new int[_nivterm[iOct]]();
	}
	

	
	// for each rank
	for (int src=0; src<_wsize; src++)
	{
		if (src != _rank)
		{
			// for each level
			for (int level=_nivterm[iOct]+indexToC; level>=_levcom[iOct]+indexToC; level--)
			{
				
				int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
				int __endlev = _endlev[iOct][level]+indexToC;
				int __endlevprec = _endlev[iOct][level-1]+indexToC;
				long * __recv = _recv[iOct];
				int __fnivnextlev = _fniv[iOct][level+1];
				int __nst = _nst[iOct][level];
				int __nsp = _nsp[iOct][level];
	
				// box range to send, FROM SEND ARRAY
				int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
				int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;
				
				
				// find box range in the targeted level
				bool start = false;
				bool stop = false;
				int firstBox, lastBox;
				
				for (int j=lastBoxToRecvIdx; j>=firstBoxToRecvIdx; j--)
				{
					int cellID = __recv[j] + indexToC;

					if (!start)
					{
						if ( (cellID > __endlevprec) && (cellID <= __endlev) )
						{
							start = true;
							lastBox = j;
							_start_recv[iOct][src][level]=j;
						}
					}
					else
					{
						if (!stop)
						{
							if ( cellID <= __endlevprec) 
							{
								stop = true;
								firstBox = j+1;
								_stop_recv[iOct][src][level] = j+1;
							}
						}
					}
				}
				/*if (!stop)
				{
					firstBox = firstBoxToRecvIdx;
					stop = true;
				}*/

				if (!start)
				{
					lastBox =  -1;
					firstBox = -1;
					_start_recv[iOct][src][level]= -1;
					_stop_recv[iOct][src][level] = -1;
				}
				
				if (start)
				{
					if(!stop)
					{
						firstBox = firstBoxToRecvIdx;
						_stop_recv[iOct][src][level] = firstBoxToRecvIdx;
					}
				}
			}
		}
	}
	t_end = MPI_Wtime();
	add_time_sec("fill_send_start_stop_count", t_end - t_begin);
}


void Gaspi_FF_communicator::send_ff_level(int level, complex * ff, int iOct)
{
	int indexToC = -1;
	double t_begin, t_end, accumul;
	accumul = 0;
	int octree_offset = iOct * _wsize;
	
	// flush queue
	t_begin = MPI_Wtime();
	flush_queues(_nbQueues);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_SEND_wait_queue", t_end - t_begin);
	accumul = accumul + (t_end - t_begin);	
	
	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	int * counters = new int[_wsize]();
	
	// prepare _FF_sendBuffer
	t_begin = MPI_Wtime();

// TODO mettre un parallel for ici + scheduler dynamic
// tester à + grande échelle
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{
			int q = _FF_sendLocalOffsets_byOctDest[octree_offset + dest] - 1 + _FF_sendLocalOffsets_keeper[iOct][dest];
			int q0 = q + 1; 

			int counter =  _count_send[iOct][dest][level];
			int lastBox =  _start_send[iOct][dest][level];
			int firstBox = _stop_send[iOct][dest][level];
			
			if (counter>0)
			{
				#pragma omp parallel for private(q)
				for (int k=lastBox; k>=firstBox; k--)
				{
					int cellID = _send[iOct][k] + indexToC;
					int p=_fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];
					q = q0 + (lastBox-k)*__nstnsp;

					/*for (int i=0; i<__nstnsp; i++) 
						_FF_sendBuffer[q+i] = ff[p+i];*/
					_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp]; 
				}
			}

		}
	}
	t_end = MPI_Wtime();
	add_time_sec("FF_write_to_buffer", t_end - t_begin);
	add_time_sec("FF_buffering", t_end - t_begin);	

	// send !
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{
			t_begin = MPI_Wtime();
			
			// if something to send
			int counter = _count_send[iOct][dest][level];
			if (counter > 0)
			{
				// local offset
				gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets_byOctDest[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t level_offset = _FF_sendLocalOffsets_keeper[iOct][dest] * sizeof(complex);
				gaspi_offset_t local_offset = local_dest_offset + level_offset;
				gaspi_queue_id_t queue=0;
				
				// remote offset
				gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t remote_offset = remote_sender_offset + level_offset;


				// update offset, per destinatary
				_FF_sendLocalOffsets_keeper[iOct][dest] += counter * _nst[iOct][level] * _nsp[iOct][level];

				//int rankMultiple = level;
				//TODO : CORRIGER CAR ARBRES PEUVENT AVOIR DES HAUTEURS DIFFERENTES
				gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
				gaspi_size_t qty= counter * _nst[iOct][level] * _nsp[iOct][level] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_FF_sendBuf_seg_id,			// local seg ID
						local_offset,						// local offset
						dest,								// receiver rank
						_FF_recvBuf_seg_id,			// remote seg ID
						remote_offset,						// remote offset
						qty,								// size of data to write
						notifyID,							// remote notif ID
						level,								// value of the notif to write
						queue,								// queue
						GASPI_BLOCK							// Gaspi block
					)
				);
			//if (_rank ==1 && level ==2 && dest == 0)
	//			printf("\n[%d] level %d, notifyID %d notifyValue %d\n", _rank, level, notifyID, level); fflush(stdout);
				printf("\n[%d] sent level %d, notifyID %d remote_offset (bytes) %d, qty (bytes) %d\n", _rank, level, notifyID, (int)remote_offset, (int)qty); fflush(stdout);

			}
			t_end = MPI_Wtime();
			add_time_sec("GASPI_SEND_write_notify", t_end - t_begin);
			accumul = accumul + (t_end - t_begin);
		}
		
		// Une fois tout fini, reset ou non le tableau d'offsets
		if(_incLevcom)
		{
			if (level == _levcom[iOct] + indexToC)
				_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
		}
		else // cas allreduce sur levcom
		{
			if (level == _levcom[iOct] + 1 + indexToC)
				_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
		}
	}
	
	add_time_sec("GASPI_FF_sendrecv", accumul);
}

void Gaspi_FF_communicator::send_task_ff_level(int level, complex * ff, int iOct, int start, int stop)
{
	//printf("Here I am rank %d thread %d , octree %d, level %d, start %d, stop %d\n",_rank, omp_get_thread_num(), iOct, level, start, stop);

	int indexToC = -1;

	int octree_offset = iOct * _wsize;
	
	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	bool levelComplete = false;
	
	// remplissage des buffers
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{
			int counter =  _count_send[iOct][dest][level];
			//cout << "before if, counter : " << counter << endl;
			
			if (counter>0)
			{
				int lastBox =  _start_send[iOct][dest][level];
				int firstBox = _stop_send[iOct][dest][level];
					
				for (int i= firstBox; i <= lastBox; i++)
				{
					int cellIDf = _send[iOct][i];
					if (cellIDf >= start && cellIDf <= stop)
					{
						//printf(" i : %d  cellID to SEND : %d, offset : %d \n", i, cellID, i-firstBox);
						
						int cellID = cellIDf + indexToC;
						
						// PROTECTION DES COMPTEURS D'OFFSETS de DATA et INFOS 
						pthread_mutex_lock(&_mutexArray[octree_offset + dest]);
									int q = _FF_sendLocalOffsets_byOctDest[octree_offset + dest] + _FF_sendLocalOffsets_keeper[iOct][dest];
									_FF_sendLocalOffsets_keeper[iOct][dest] += _nst[iOct][level]*_nsp[iOct][level];
									
									int idxInfos = _Infos_sendLocalOffsets_keeper[iOct][dest];
									_Infos_sendLocalOffsets_keeper[iOct][dest] += 1;
								
									// Ecriture des DATA et des INFOS		
									int p=_fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];
									_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp]; /*for (int i=0; i<__nstnsp; i++) _FF_sendBuffer[q+i] = ff[p+i];*/
									_Infos_sendbuffer[idxInfos] = i-firstBox;	
						pthread_mutex_unlock(&_mutexArray[octree_offset + dest]);
						// SI LE NIVEAU EST FINI : ENVOI
						if (_Infos_sendLocalOffsets_keeper[iOct][dest] == _count_send[iOct][dest][level])
						{
							levelComplete = true;
						}					
					}
				}
			}
			
			// ENVOI d'un niveau complet dans le desordre
			// TODO decouper l'envoi en petits bouts de taille reglable
			if (levelComplete)
			{
				/*printf("\n[%d:%d] Chunk start %d, Chunk stop %d, offsetKeeper %d data to send %d --> I will SEND octree %d level %d\n",
				   _rank, omp_get_thread_num(), start, stop, _Infos_sendLocalOffsets_keeper[iOct][dest], _count_send[iOct][dest][level], iOct, level);*/
				//printf("\n[%d:%d] prepare send iOct %d\n", _rank, omp_get_thread_num(), iOct); fflush(stdout);
				
				flush_queues(_nbQueues);
				
				// local FF send offset
				gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets_byOctDest[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t level_offset = _FF_sendRemoteOffsets_keeper[iOct][dest] * sizeof(complex);/*_FF_sendLocalOffsets_keeper[iOct][dest] * sizeof(complex);*/
				gaspi_offset_t local_offset = local_dest_offset + level_offset;
				gaspi_queue_id_t queue=0;
				
				// remote offset
				//printf("\n[%d:%d] _FF_sendRemoteOffsets contains : %d, level_offset : %d \n", _rank, omp_get_thread_num(), _FF_sendRemoteOffsets[octree_offset + dest], int(level_offset)); fflush(stdout);
				gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t remote_offset = remote_sender_offset + level_offset;

				/*cout << "nst * nsp = " << __nstnsp << endl;
				cout << "counter = " << counter << endl;
				cout << "qty = " << counter * __nstnsp * sizeof(complex);*/

				_FF_sendRemoteOffsets_keeper[iOct][dest] += counter * __nstnsp;
				gaspi_size_t qty= counter * __nstnsp * sizeof(complex);
				
				//cout << "qty : " << qty <<  ", counter : " << counter << endl;
				
				// envoi des DATA sans notification
				SUCCESS_OR_DIE(
					gaspi_write( 
						_FF_sendBuf_seg_id,					// local seg ID
						local_offset,						// local offset
						dest,								// receiver rank
						_FF_recvBuf_seg_id,					// remote seg ID
						remote_offset,						// remote offset
						qty,								// size of data to write
						queue,								// queue
						GASPI_BLOCK							// Gaspi block
					)
				);
				//printf("\n[%d:%d] oct %d sent level %d, remote_offset (bytes) %d, qty (bytes) %d\n", _rank, omp_get_thread_num(), iOct, level, (int)remote_offset, (int)qty); fflush(stdout);
				
				//printf("\n[%d:%d] sent to %d from offset %d at offset %d qty %d\n",
				  // _rank, omp_get_thread_num(), dest, local_offset, remote_offset, qty);
				//int rankMultiple = level;
				
				//TODO : CORRIGER CAR ARBRES PEUVENT AVOIR DES HAUTEURS DIFFERENTES
				gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
				local_offset = _Infos_sendRemoteOffsets_keeper[iOct][dest] * sizeof(int);
				remote_offset = _Infos_sendRemoteOffsets[octree_offset + dest] * sizeof(int);
				_Infos_sendRemoteOffsets_keeper[iOct][dest] += counter;
				qty = counter * sizeof(int);
				
				/*if (_rank == 1)
				{
					cout << "local_offset = " << local_offset << endl;
					cout << "remote_offset = " <<remote_offset << endl;
					cout << "quantity : " << counter << endl;
					for (int i=local_offset; i<local_offset+counter; i++)
						cout << i << " send info : " << _Infos_sendbuffer[i] << endl;					
				}*/
				//cout << "qty = " << qty;
								
				// envoi des INFO et NOTIFICATION
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_Infos_sendbuf_seg_id,					// local seg ID
						local_offset,						// local offset
						dest,								// receiver rank
						_Infos_recvbuf_seg_id,					// remote seg ID
						remote_offset,						// remote offset
						qty,								// size of data to write
						notifyID,							// remote notif ID
						counter,							// value of the notif to write
						queue,								// queue
						GASPI_BLOCK							// Gaspi block
					)
				);
				/*if (_rank == 1){
				printf("\n[%d:%d] Chunk start %d, Chunk stop %d, offsetKeeper %d data to send %d -->Sent data octree %d level %d, remote_offset %d, qty %d\n",
				   _rank, omp_get_thread_num(), start, stop, _Infos_sendLocalOffsets_keeper[iOct][dest], _count_send[iOct][dest][level], iOct, level, remote_offset, counter);
				}*/
				printf(	"\n[%d:%d] oct %d --> gaspi_write_notify "
						"local_offset %d "
						"dest %d "
						"remote_offset %d "
						"qty %d "
						"\n local area : %d to %d (from sendInfos)"
						"\n remote area : %d to %d (to recvInfos)"
						"\n notifyID %d "
						"counter %d\n", 
						_rank, omp_get_thread_num(), iOct, 
						(int) local_offset,
						dest, 
						(int)remote_offset, 
						(int)qty,
						(int) local_offset, (int) local_offset + int(qty),
						(int) remote_offset, (int) remote_offset + int(qty), 						
						(int)notifyID, 
						(int)counter); 
				fflush(stdout);
				
				// Une fois tout fini, reset ou non le tableau d'offsets
				/*if(_incLevcom)
				{
					if (level == _levcom[iOct] + indexToC)
					{
						_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
						_Infos_sendLocalOffsets_keeper[iOct][dest] = 0;
					}
				}
				else // cas allreduce sur levcom
				{
					if (level == _levcom[iOct] + 1 + indexToC)
					{
						_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
						_Infos_sendLocalOffsets_keeper[iOct][dest] = 0;
					}
				}*/
				
				// temporairement reset toujours le tableau d'offset, car n'est utilisé que pour le dernier niveau
				_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
				_Infos_sendLocalOffsets_keeper[iOct][dest] = 0;
				_FF_sendRemoteOffsets_keeper[iOct][dest] = 0;
				_Infos_sendRemoteOffsets_keeper[iOct][dest] = 0;
				
			}
		}
	}
}

void Gaspi_FF_communicator::recv_task_ff_level(int level, complex * ff, int iOct)
{
	//cout << "toto" << endl;
	// wait to receive all infos
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected++;
		}
	}

	gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int recvCpt = 0;
	int sender;
	int counter;
	
	while (recvCpt < nbRecvExpected)
	{
		while(1)
		{
			//if (_rank ==0 && level ==4)
			printf("\n[%d] level %d, waiting from offset : %d\n",_rank, level, notif_offset); fflush(stdout);
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_Infos_recvbuf_seg_id,
					notif_offset,				// surveille les notifications depuis notif_offset
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			//if (_rank ==0 && level ==4)
			printf("\n[%d] level %d, received \n",_rank, level); fflush(stdout);
						
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_Infos_recvbuf_seg_id, 
					new_notif_id, 
					&new_notif_val
				)
			);
			
			if (new_notif_val) 
				break;
		}

		// test the notification value and compare with array (I know, useless)
		if (new_notif_val)
		{
			recvCpt++;
			sender = new_notif_id - notif_offset;
			counter = new_notif_val;
			/*if (_rank ==0 && level ==2)
				printf("\n[%d] level %d, before update far fields : %d\n",_rank, level, notif_offset); fflush(stdout);			*/
			//updateFarFields(sender, level, ff, iOct);
			updateFarFieldsFromInfos(sender, level, ff, counter, iOct);
		}
	}
}


void Gaspi_FF_communicator::recv_ff_level(int level, complex * ff, int iOct)
{

	// wait to receive all infos
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected++;
		}
	}

	gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int recvCpt = 0;
	int sender;
	double t_begin, t_end;
	
	while (recvCpt < nbRecvExpected)
	{
		t_begin = MPI_Wtime();
		while(1)
		{
			//if (_rank ==0 && level ==2)
			printf("\n[%d] level %d, waiting from notif offset : %d\n",_rank, level, notif_offset); fflush(stdout);
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_FF_recvBuf_seg_id,
					notif_offset,				// surveille les notifications depuis notif_offset
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
			//if (_rank ==0 && level ==2)
			printf("\n[%d] level %d received\n",_rank, level); fflush(stdout);
						
			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_FF_recvBuf_seg_id, 
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
			/*if (_rank ==0 && level ==2)
				printf("\n[%d] level %d, before update far fields : %d\n",_rank, level, notif_offset); fflush(stdout);			*/
			updateFarFields(sender, level, ff, iOct);
		}
		t_end = MPI_Wtime();
		add_time_sec("FF_read_from_buffer", t_end - t_begin);
		add_time_sec("FF_buffering", t_end - t_begin);
	}
}

void Gaspi_FF_communicator::updateFarFields(int src, int level, complex * ff, int iOct)
{			
	int indexToC = -1;
	int k = level + 1;
	int octree_offset = iOct * _wsize;
	double t_begin, t_end;
	
	// calcul de l'offset	
	int levelOffset = 0;
	for (int i=_nivterm[iOct] + indexToC; i>level; i--)
	{
		levelOffset += _nst[iOct][i] * _nsp[iOct][i] * _Expect[iOct][src][i];
	}
	int q =_FF_recvOffsets[octree_offset + src]-1; // se mettre 1 case avant l'index à lire
	q += levelOffset;
	int q0 = q + 1; 	
	
	// variables precalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	int __endlev = _endlev[iOct][level]+indexToC;
	long * __recv = _recv[iOct];
	int __fnivnextlev = _fniv[iOct][level+1];

	// box range to recv, FROM RECV ARRAY
	int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
	int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;

	t_begin = MPI_Wtime();
	
	bool start    = (_start_recv[iOct][src][level]>=0);
	bool stop     = ( _stop_recv[iOct][src][level]>=0);
	int  firstBox =   _stop_recv[iOct][src][level]; 
	int  lastBox  =  _start_recv[iOct][src][level];
	
	if (start && stop)
	{
		#pragma omp parallel for private (q)
		for (int j=lastBox; j>=firstBox; j--)
		{
			q = q0 + (lastBox-j)*__nstnsp;		
			int cellID = __recv[j] + indexToC; 
			int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);			
			//for (int i=0; i<__nstnsp; i++)
			/*{
				ff[p0+i].re = ff[p0+i].re + _FF_recvBuffer[q+i].re;
				ff[p0+i].im = ff[p0+i].im + _FF_recvBuffer[q+i].im;
			}*/			
			ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[q:__nstnsp];
		}
	}

	t_end = MPI_Wtime();
	add_time_sec("Boucle_read_from_buffer", t_end - t_begin);
}

void Gaspi_FF_communicator::updateFarFieldsFromInfos(int src, int level, complex * ff, int counter, int iOct)
{			
	int indexToC = -1;
	int k = level + 1;
	int octree_offset = iOct * _wsize;

	int offset = _Infos_recvOffsets[octree_offset + src];
	int firstBox = _stop_recv[iOct][src][level];
	int index;
	int recvDataIdx;

	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	int __endlev = _endlev[iOct][level]+indexToC;
	int __fnivnextlev = _fniv[iOct][level+1];
		
	/*if (_rank == 0)
	{
		cout << "offset in infos : " << offset << " counter : " << counter << " source : " << src << endl;
		*/
		for (int i=offset; i<offset+counter; i++)
		{
		//	cout << i << " should read index : " << _Infos_recvbuffer[i] << endl;
			index = _Infos_recvbuffer[i];
			int cellIDf = _recv[iOct][firstBox+index];
			int cellID = cellIDf + indexToC;
			//cout << "cell idF : " << cellIDf << endl;
			
			// ranger la boîte
			//q = q0 + (lastBox-j)*__nstnsp;		
			//recvDataIdx = 
			
			int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);			
		
			ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[i:__nstnsp];
		}
	//}
	cout << "euh ... updated without segfault !" << endl; 
}
