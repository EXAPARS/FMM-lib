#include "Gaspi_FF_communicator.hpp"
using namespace std;

/* ---------------------------------------------- */
/*               GASPI BULK VERSION               */
/* ---------------------------------------------- */

void Gaspi_FF_communicator::exchangeFFBulk(complex * bufsave, complex * ff, int iOct)
{
	//dumpBuffer(_rank, _codech[iOct], 8000, "toto", "");
	//void dumpBuffer(int rank, T * buffer, int size, string fileName, string message)

	//cout << "GASPI BULK - exchange begin" << endl;
	//fflush(stdout);
	
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
			
			int q = _FF_sendLocalOffsets[iOct][dest] - 1;
			cout << std::scientific;
			// remplissage du buffer d'envoi, par ordre décroissant
			for (int k = lastBoxToSendIDX; k>=firstBoxToSendIDX; k--)
			{
				bool todo = true;
				int p;
				int cellID = _send[iOct][k] + indexToC;

				if ((_levcom[iOct]<=3 || (_levcom[iOct] == 4 && _nivterm[iOct]>8))  && (cellID <  _endlev[iOct][_levcom[iOct] + indexToC] + indexToC))
					todo = false;
				if ((_levcom[iOct]>4  || (_levcom[iOct] == 4 && _nivterm[iOct]<=8)) && (cellID <= _endlev[iOct][_levcom[iOct]-1 + indexToC] + indexToC))
					todo = false;

				if (todo)
				{
					// place k pour que k = niveau de la cellule à envoyer
					while(cellID <= _endlev[iOct][niv-1 + indexToC] + indexToC)
						niv--;
					
					// cas 1 : 
					// si l est commun à +ieurs noeuds
					// communique les données sauvegardées dans BUFSAVE
					//debug("to_" + convert(dest), convert(cellID));
					if (_codech[iOct][cellID] > 1)
					{
						p = _codech[iOct][cellID]-2;
						//debug("to_" + convert(dest), convert(cellID) + " sendBuf : " + convert(q) +  " bufsave " + convert(p+indexToC) + " p : " + convert(p));	
						
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
						p=_fniv[iOct][niv+1 + indexToC]+(_endlev[iOct][niv + indexToC]+indexToC-cellID)*_nst[iOct][niv + indexToC]*_nsp[iOct][niv + indexToC];
						//debug("to_" + convert(dest), convert(cellID) + " sendBuf : " + convert(q) +  "ff " + convert(p+indexToC));	

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
				local_offset  = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
				remote_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
				qty = _nb_send[iOct][dest] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_FF_sendBuf_seg_id,					// local seg ID
						local_offset,						// local offset
						dest,								// receiver rank
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
						_FF_recvBuf_seg_id,					// seg
						dest,		 						// receiver rank
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

	//cout << "GASPI BULK - exchange out" << endl; 
	//fflush(stdout);
}

void Gaspi_FF_communicator::updateFarFields(int src, complex * ff, int iOct)
{		
	int indexToC = -1;
	int k = _nivterm[iOct];
	int srcF = src+1;
	
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
			int pp = _fniv[iOct][k+1 + indexToC]+(_endlev[iOct][k + indexToC] - p )*_nst[iOct][k + indexToC]*_nsp[iOct][k + indexToC];
			
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
	_msgID = 0;
	_max_comm = 65536/_wsize;
		
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
	// FF
	_FF_sendLocalOffsets_counter = new int*[_nbOct]();			//write into segment
	for (int i=0; i<_nbOct; i++)
		_FF_sendLocalOffsets_counter[i] = new int[_wsize]();
	
	_FF_sendLocalOffsets = new int*[_nbOct]();
		for (int i=0; i<_nbOct; i++)
			_FF_sendLocalOffsets[i] =  new int[_wsize]();
		
	_FF_sendRemoteOffsets_counter = new int*[_nbOct]();		//send from segment
	for (int i=0; i<_nbOct; i++)
		_FF_sendRemoteOffsets_counter[i] = new int[_wsize]();		
	
	// Infos
	_Infos_sendLocalOffsets = new int*[_nbOct](); 		//	 Infos Local send Offsets
	for (int i=0; i<_nbOct; i++)
		_Infos_sendLocalOffsets[i] = new int[_wsize]();	
		
	_Infos_sendLocalOffsets_counter = new int*[_nbOct]();			// write into segment
	for (int i=0; i<_nbOct; i++)
		_Infos_sendLocalOffsets_counter[i] = new int[_wsize]();
		
	_Infos_sendRemoteOffsets_counter = new int*[_nbOct]();			// send from segment
	for (int i=0; i<_nbOct; i++)
		_Infos_sendRemoteOffsets_counter[i] = new int[_wsize]();		
	
	// Keepers
	_FF_sendLocalOffsets_keeper = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_FF_sendLocalOffsets_keeper[i] = new int[_wsize]();
	
	_Infos_sendLocalOffsets_keeper = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_Infos_sendLocalOffsets_keeper[i] = new int[_wsize]();
	
	// Quantities
	_Expect = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_Expect[i] = new int*[_wsize]();
    
    _Received = new int**[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_Received[i] = new int*[_wsize]();    

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
	
	_accumul = new int[_wsize]();
	
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

	//cout << _rank << "max send nodes : " << max_send_nodes << " max recv nodes : " << max_recv_nodes << endl;

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
	gaspi_segment_ptr(_FF_recvOffsets_seg_id, &_FF_recvOffsets_seg_ptr);
	_FF_recvOffsets = (int *)_FF_recvOffsets_seg_ptr;
	
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
	gaspi_segment_ptr(_Infos_recvOffsets_seg_id, &_Infos_recvOffsets_seg_ptr);
	_Infos_recvOffsets = (int *)_Infos_recvOffsets_seg_ptr;	
		
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
	
	// FF offsets
	fill_FF_sendRemoteOffsets(recvnode, recvnode_sz, nb_recv, nb_recv_sz, iOct);
	fill_FF_sendLocalOffsets(sendnode, sendnode_sz, nb_send, nb_send_sz, iOct);
	
	// expected data
	fill_expectations(iOct);
	fill_send_start_stop_count(iOct);
	fill_recv_start_stop(iOct);
	
	// send infos pointers
	fill_Infos_recvOffsets(iOct);
	fill_Infos_sendRemoteOffsets(iOct);
	
	// Info offsets
	fill_Info_sendLocalOffsets(iOct);
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

void Gaspi_FF_communicator::fill_FF_sendRemoteOffsets(i64 * recvnode, int recvnode_sz, i64 * nb_recv, int nb_recv_sz, int iOct)
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

void Gaspi_FF_communicator::fill_Infos_sendRemoteOffsets(int iOct)
{	
	// From Fortran to C/C++
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
        
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
}

void Gaspi_FF_communicator::fill_FF_sendLocalOffsets(i64 * sendnode, int sendnode_sz, i64 * nb_send, int nb_send_sz, int iOct)
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
				_FF_sendLocalOffsets[iOct][to] = SendBufferIndexPerROUND[k];
			}
			else
			{
				cerr << "Error in Gaspi_m2l_communicator::init_dataToSendIndexes.\nFortran sendnode array should not contain the current rank !";
				exit(0);
			}
		}
    }  
	_FF_sendLocalOffsets[iOct][_rank] = -1;
		
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
		_Received[iOct][j] = new int[_nivterm[iOct]]();
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

void Gaspi_FF_communicator::fill_send_start_stop_count(int iOct)
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
	
	/*for (int dest=0; dest<_wsize; dest++)
	{
		for (int level=_nivterm[iOct]+indexToC; level>=_levcom[iOct]+indexToC; level--)
		{
			debug("send_indexes", "oct : " + itoa(iOct) + " dest : " + itoa(dest) + " level : " + 
				itoa(level) + " stop : " + itoa(_stop_send[iOct][dest][level]) + " start : " + itoa(_start_send[iOct][dest][level]) + " qty : " + itoa(_count_send[iOct][dest][level]));
		}
	}*/
}

void Gaspi_FF_communicator::fill_Infos_recvOffsets(int iOct)
{
	int indexToC = -1;
	int octree_offset = iOct * _wsize;

	for (int src=0; src<_wsize; src++)
	{
		// zero
		if (src == 0)
		{
			_Infos_recvOffsets[octree_offset + src] = 0;
		}
		// any other
		else
		{
			int sum = 0;
			// pour chaque rank avant
			for (int i=0; i<src; i++)
			{
				// somme de tout ce qui a ete recu
				for (int level= 0; level<=_nivterm[iOct]+indexToC; level++)
				{
					sum += _Expect[iOct][i][level];
				}
			}
			_Infos_recvOffsets[octree_offset + src] = sum;
		}
	}
	// me
	_Infos_recvOffsets[octree_offset + _rank] = -1;
}

void Gaspi_FF_communicator::fill_recv_start_stop(int iOct)
{
	double t_begin, t_end;
	t_begin = MPI_Wtime();
	int indexToC = -1;

	// init last dim only
	for (int src=0; src<_wsize; src++)
	{
		_start_recv[iOct][src] = new int[_nivterm[iOct]];
		_stop_recv[iOct][src] = new int[_nivterm[iOct]];
		
		// attention, ne pas initialiser à zero !
		for (int level=0; level<_nivterm[iOct]; level++)
		{
			_start_recv[iOct][src][level]=-1;
			_stop_recv[iOct][src][level]=-1;
		}
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
	
	/*for (int src=0; src<_wsize; src++)
	{
		for (int level=_nivterm[iOct]+indexToC; level>=_levcom[iOct]+indexToC; level--)
		{
			debug("recv_indexes", "oct : " + itoa(iOct) + " src : " + itoa(src) + " level : " + 
				itoa(level) + " stop : " + itoa(_stop_recv[iOct][src][level]) + " start : " + itoa(_start_recv[iOct][src][level]));
				
				if (src==1 && level ==4)
				{
					if (_start_recv[iOct][src][level] != _stop_recv[iOct][src][level])
					{
						for (int idx=_start_recv[iOct][src][level]; idx>=_stop_recv[iOct][src][level]; idx--)
						{
							debug ("cells_from_1", itoa(_recv[iOct][idx]-1) + " l = " + itoa(level) + " iOct = " + itoa(iOct));
						}
					}
				}
		}
	}*/
}

void Gaspi_FF_communicator::fill_Info_sendLocalOffsets(int iOct)
{
	int indexToC = -1;

	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest == 0)
		{
			_Infos_sendLocalOffsets[iOct][dest] = 0;
		}
		else
		{
			int sum = 0;
			// pour chaque rank avant
			for (int i=0; i<dest; i++)
			{	
				// somme de tout ce qui a été envoyé		
				for (int level= 0; level<=_nivterm[iOct]+indexToC; level++)
				{
					sum += _count_send[iOct][i][level];
				}
			}		
			_Infos_sendLocalOffsets[iOct][dest] = sum;
		}
	}
	// me
	_Infos_sendLocalOffsets[iOct][_rank] = -1;
}

/*********************************************************************************
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
			int q = _FF_sendLocalOffsets[iOct][dest] - 1 + _FF_sendLocalOffsets_counter[iOct][dest];
			int q0 = q + 1; 

			int counter =  _count_send[iOct][dest][level];
			int lastBox =  _start_send[iOct][dest][level];
			int firstBox = _stop_send[iOct][dest][level];
			
			if (counter>0)
			{
				//#pragma omp parallel for private(q)
				for (int k=lastBox; k>=firstBox; k--)
				{
					int cellID = _send[iOct][k] + indexToC;
					//if (level == 4)
						//debug("send", to_string((long long)(cellID)));
					int p=_fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];
					q = q0 + (lastBox-k)*__nstnsp;
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
				gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
				gaspi_offset_t level_offset = _FF_sendLocalOffsets_counter[iOct][dest] * sizeof(complex);
				gaspi_offset_t local_offset = local_dest_offset + level_offset;
				gaspi_queue_id_t queue=0;
				
				// remote offset
				gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t remote_offset = remote_sender_offset + level_offset;


				// update offset, per destinatary
				_FF_sendLocalOffsets_counter[iOct][dest] += counter * _nst[iOct][level] * _nsp[iOct][level];

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
			}
			t_end = MPI_Wtime();
			add_time_sec("GASPI_SEND_write_notify", t_end - t_begin);
			accumul = accumul + (t_end - t_begin);
		}
		
		// Une fois tout fini, reset ou non le tableau d'offsets
		if(_incLevcom)
		{
			if (level == _levcom[iOct] + indexToC)// +1 if invalidated comm at the top
			{
				_FF_sendLocalOffsets_counter[iOct][dest] = 0;
				
				// pour la compatibilité avec les tasks, en cas parallélisation extrapomp
				// RAZ des compteurs
				_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
				_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
				_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
			}
		}
		else // cas allreduce sur levcom
		{
			if (level == _levcom[iOct] + 1 + indexToC) // +1 if invalidated at the top
				_FF_sendLocalOffsets_counter[iOct][dest] = 0;
		}
	}
	add_time_sec("GASPI_FF_sendrecv", accumul);
}
*********************************************************************************/

void Gaspi_FF_communicator::send_task_ff_level(int level, complex * ff, int iOct, int start, int stop)
{
	//pthread_mutex_lock(&mutex);
	int indexToC = -1;
	
	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
	// for each dest
	for (int dest=0; dest<_wsize; dest++)
	{								
		if (dest != _rank)
		{
			int counter =  _count_send[iOct][dest][level];
			// printf("Src %d Dest %d Ioct %d level %d counter %d\n", _rank, dest, iOct, level, counter); fflush(stdout);
			
			// if something to send to this dest
			if (counter>0)
			{
				int lastBox =  _start_send[iOct][dest][level];	// the big one
				int firstBox = _stop_send[iOct][dest][level];	// the small one

				// traverses the cells to send
				for (int i= firstBox; i <= lastBox; i++)
				{
					int cellIDf = _send[iOct][i];
					
					// if the cell is in the task's nodes range
					if (cellIDf >= start && cellIDf <= stop)
					{						
						int cellID = cellIDf + indexToC;
						
						// get indexes 
						int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];						// in FF senbuffer
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];	// cell's corresponding ff in Fortran
						int idxInfos =_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];
	
						// write FF and Infos		
						_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp]; 
						_Infos_sendbuffer[idxInfos] = i-firstBox;
						
						/*if (level == 4)
						{
							debug("infos_to_" + convert(dest), convert(i-firstBox) + " value : " + to_string( _FF_sendBuffer[q].re) + " " + to_string(_FF_sendBuffer[q].im));
							debug("infos_to_" + convert(dest), to_string( _FF_sendBuffer[q].re) + " " + to_string(_FF_sendBuffer[q].im));

						}*/
											
						// update indexes
						_FF_sendLocalOffsets_counter[iOct][dest] += _nst[iOct][level]*_nsp[iOct][level];
						_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
					
						// if chunk is full  --> send !
						if (_Infos_sendLocalOffsets_counter[iOct][dest] == _count_send[iOct][dest][level])
						{	
							_Infos_sendLocalOffsets_counter[iOct][dest]=0;
							send_chunk(iOct, dest, level);
						} 
					}
				}	
			}
		}
	}
	//pthread_mutex_unlock(&mutex); 
}

void Gaspi_FF_communicator::send_chunk(int iOct, int dest, int level)
{
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	int counter =  _count_send[iOct][dest][level];
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];

	flush_queues(_nbQueues);
			
	// local FF send offset
	gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
	gaspi_offset_t level_offset = _FF_sendRemoteOffsets_counter[iOct][dest] * sizeof(complex);
	gaspi_offset_t local_offset_ff = local_dest_offset + level_offset;
	gaspi_queue_id_t queue=0;
	
	// remote offset
	gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
	gaspi_offset_t remote_offset_ff = remote_sender_offset + level_offset;
	_FF_sendRemoteOffsets_counter[iOct][dest] += counter * __nstnsp;
	gaspi_size_t qty= counter * __nstnsp * sizeof(complex);
	
	
	// envoi des DATA sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_FF_sendBuf_seg_id,					// local seg ID
			local_offset_ff,						// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,						// remote offset
			qty,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	//TODO : CORRIGER CAR ARBRES PEUVENT AVOIR DES HAUTEURS DIFFERENTES
	gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
	gaspi_offset_t remote_offset = (_Infos_sendRemoteOffsets[octree_offset + dest] +_Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	_Infos_sendRemoteOffsets_counter[iOct][dest] += counter;

	gaspi_offset_t local_offset = _Infos_sendLocalOffsets[iOct][dest] * sizeof(int);
	qty = counter * sizeof(int);
					
	// envoi des INFO et NOTIFICATION
	SUCCESS_OR_DIE(
		gaspi_write_notify( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_offset,						// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_offset,						// remote offset
			qty,								// size of data to write
			notifyID,							// remote notif ID
			counter,							// value of the notif to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
	/*if (level == 4)
	{
		debug("infos_bis_to_" + convert(dest), convert(_Infos_sendbuffer[local_offset/sizeof(int)]) + " at " + convert(remote_offset/sizeof(int)) +
			" value : " + to_string( _FF_sendBuffer[local_offset_ff/sizeof(complex)].re) + " " + to_string(_FF_sendBuffer[local_offset_ff/sizeof(complex)].im) +
			" ff_off : " + convert(remote_offset_ff/sizeof(complex)) );				
	}*/
	
	// Une fois tout fini, reset ou non le tableau d'offsets
	if(_incLevcom)
	{
		if (level == _levcom[iOct] + indexToC ) 
		{
			_FF_sendLocalOffsets_counter[iOct][dest] = 0;
			_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
			_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
			_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
		}
	}
	else // cas allreduce sur levcom
	{
		if (level == _levcom[iOct] + 1 + indexToC)
		{
			_FF_sendLocalOffsets_counter[iOct][dest] = 0;
			_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
			_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
			_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
		}
	}
	
	//printf("Src %d Dest %d Ioct %d level %d start %d stop %d\n", _rank, dest, iOct, level, start, stop);
	
	/*SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));*/
							
	/*if(level==4 && _rank==1 && dest==0 && iOct==2)
	{
		for (int i=0; i<counter; i++)
			debug("indexes_sent_to_0", itoa(_Infos_sendbuffer[_Infos_sendLocalOffsets[iOct][dest]]+i));
	}*/
	/*for (int i=0; i<counter; i++)
		accumulMSG("iOct" + itoa(iOct) + " to " + itoa(dest) + " " + itoa(_Infos_sendbuffer[_Infos_sendLocalOffsets[iOct][dest]]+i));*/
}

void Gaspi_FF_communicator::recv_task_ff_level(int level, complex * ff, int iOct)
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
	int counter;
	
	while (recvCpt < nbRecvExpected)
	{
		while(1)
		{
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_Infos_recvbuf_seg_id,
					notif_offset,				// surveille les notifications depuis notif_offset
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
						
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
			updateFarFieldsFromInfos(sender, level, ff, counter, iOct);
		}
	}
	/*if (level + 1 == _levcom[iOct])
	{
		// print nb received
		for (int i=0; i<_wsize; i++)
		{
			printf("%d received %d nodes from %d at level %d\n",_rank, _Expect[iOct][i][level], i, level);
		}
		
		//print ff
		switch(iOct) 
		{
			case 0 : dumpBuffer(_rank, ff, 84992, "ff", ""); break;
			case 1 : dumpBuffer(_rank, ff, 129720, "ff", ""); break;
			case 2 : dumpBuffer(_rank, ff, 129720, "ff", ""); break;
			default : cout <<"not the right octree id"; break;
		}
	}*/
}

void Gaspi_FF_communicator::recv_task_ff(int level, complex * ff, int iOct)
{
	//printf("%d enter recv level %d \n", _rank, level); 	fflush(stdout);


	/*debug(
		"recv_task_ff_", 
		"level : " + to_string(level) + ", ioct : " + to_string(iOct) 
	);*/
	
	// compute nb of nodes to receive at current level
	int expect = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			expect += _Expect[iOct][i][level]; // tout le level !
		}
	}
	//printf("%d updated expect level %d \n", _rank, level); 	fflush(stdout);

	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int src;
	int nbNodes;
	int niv;
	
	// initialize cpt
	int cpt = 0;
	for (int i=0; i<_wsize; i++)
	{
		cpt += _Received[iOct][i][level];
	}

//printf("%d updated cpt level %d - cpt %d expect %d\n", _rank, level, cpt, expect); fflush(stdout);

	
	while (cpt < expect)
	{
		// METHODE 1
		/*while(1)
		{
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_Infos_recvbuf_seg_id,
					0,
					65536,
					&new_notif_id,
					GASPI_BLOCK
				)
			);

			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_Infos_recvbuf_seg_id, 
					new_notif_id, 
					&new_notif_val
				)
			);
			
			if (new_notif_val) 
				break;
		}*/

		/*if (new_notif_val)
		{
			src = int(new_notif_id / _max_comm);
			niv = int(new_notif_val >> 16);
			nbNodes = int(0x00001111) & int(new_notif_val);
			printf("[%d] received %d nodes from src %d level %d\n", _rank, nbNodes, src, level); fflush(stdout);

			updateFarFieldChunksFromInfos(src, niv, ff, nbNodes, iOct);
			_Received[iOct][src][niv] += nbNodes;
			cpt += nbNodes;
		}*/
		
		//printf("%d waiting cpt %d expect %d level %d\n",_rank, cpt, expect, level); fflush(stdout);

		
		// METHODE 2
		if(gaspi_notify_waitsome(
			_Infos_recvbuf_seg_id,
			0,
			65536,
			&new_notif_id,
			GASPI_TEST) == GASPI_SUCCESS)
		{
			src = (uint16_t)new_notif_id / _max_comm;
			gaspi_notify_reset(_Infos_recvbuf_seg_id, new_notif_id, &new_notif_val);

			niv = (uint32_t)(new_notif_val) >> 16;
			nbNodes = (uint32_t)(0x0000FFFF) & (uint32_t)(new_notif_val);
			/*debug("received", 
				"src : " + to_string(src) + 
				", niv : " + to_string(niv) + 
				", nbNodes : " + to_string(nbNodes) +
				", iOct : " + to_string(iOct)
				);*/
			//printf("%d received nbNodes %d level %d notif %u size %d\n",_rank, nbNodes, niv, (unsigned int)new_notif_val, sizeof(new_notif_val)); fflush(stdout);

			updateFarFieldChunksFromInfos(src, niv, ff, nbNodes, iOct);
			_Received[iOct][src][niv] += nbNodes;
			if (niv == level)
			{
				cpt += nbNodes;
				//printf("%d received nbNodes %d cpt %d expect %d level %d notif %u\n",_rank, nbNodes, cpt, expect, level, (unsigned int)new_notif_val); fflush(stdout);

			}
		}
	}
	//debug if finished
	/*if (level + 1 == _levcom[iOct])
	{
		// print ff 
		switch(iOct) 
		{
			case 0 : dumpBuffer(_rank, ff, 84992, "ff", ""); break;
			case 1 : dumpBuffer(_rank, ff, 129720, "ff", ""); break;
			case 2 : dumpBuffer(_rank, ff, 129720, "ff", ""); break;
			default : cout <<"not the right octree id"; break;
		}
	}*/
	
	// raz for next time
	for (int i=0; i<_wsize; i++)
	{
		_Received[iOct][i][level] = 0;
	}
	//printf("%d exit recv level %d - cpt %d expect %d\n", _rank, level, cpt, expect);
	
}


void Gaspi_FF_communicator::send_task_ff(int level, complex * ff, int iOct, int start, int stop)
{
	//printf("send task ff\n");

	//pthread_mutex_lock(&mutex);
	int indexToC = -1;
	int gaspiChunkSize = 4000; // in bytes
	bool levelIsFinished = false;
	int termsSize = 0;
	int nbTerms = 0;
	int nbInfos = 0;
	
	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
	// for each dest
	for (int dest=0; dest<_wsize; dest++)
	{
		levelIsFinished = false; 
		if (dest != _rank)
		{
			int counter =  _count_send[iOct][dest][level];
			//printf("Src %d Dest %d Ioct %d level %d counter %d\n", _rank, dest, iOct, level, counter); fflush(stdout);

			// if something to send to this dest ///  S IL EN RESTE ENCORE A ENVOYER
			if (counter>0)
			{
				int lastBox =  _start_send[iOct][dest][level];	// the big one
				int firstBox = _stop_send[iOct][dest][level];	// the small one

				// traverses the cells to send
				for (int i= firstBox; i <= lastBox; i++)
				{
					int cellIDf = _send[iOct][i];
					//printf("test cellIDf %d\n", cellIDf);
					
					// if the cell is in the task's nodes range
					if (cellIDf >= start && cellIDf <= stop)
					{
						int cellID = cellIDf + indexToC;
						
						// get indexes 
						int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];						// in FF senbuffer
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];	// cell's corresponding ff in Fortran
						int idxInfos =_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];
	
						// write FF and Infos		
						_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp]; 
						_Infos_sendbuffer[idxInfos] = i-firstBox;
											
						// update counters
						_FF_sendLocalOffsets_counter[iOct][dest] += _nst[iOct][level]*_nsp[iOct][level];
						_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
						_accumul[dest] += 1;
						
						nbTerms = _FF_sendLocalOffsets_counter[iOct][dest]-_FF_sendLocalOffsets_keeper[iOct][dest];
						termsSize = nbTerms * sizeof(complex);
						nbInfos = _Infos_sendLocalOffsets_counter[iOct][dest]-_Infos_sendLocalOffsets_keeper[iOct][dest];

						//printf("Src %d Dest %d Ioct %d level %d counter %d - accumul %d\n", _rank, dest, iOct, level, counter, _accumul[dest]); fflush(stdout);

						if (_accumul[dest] == counter)
						{	
							levelIsFinished = true;
							_accumul[dest] = 0;
						}
												
						if (termsSize >= gaspiChunkSize)	// if chunk is full, send !
						{
							// send
							send_chunk(iOct, dest, level, nbTerms, nbInfos, levelIsFinished);
							//printf("[%d] send nbinfo = %d - size : %d - level : %d\n", _rank, nbInfos, termsSize, level);
							
							// update keepers
							_FF_sendLocalOffsets_keeper[iOct][dest] = _FF_sendLocalOffsets_counter[iOct][dest];
							_Infos_sendLocalOffsets_keeper[iOct][dest] = _Infos_sendLocalOffsets_counter[iOct][dest];
						}
						else if (levelIsFinished) // if chunk is not full, but level is finished, send !
						{
							// send
							send_chunk(iOct, dest, level, nbTerms, nbInfos, levelIsFinished);
							// printf("[%d] send info = %d\n", _rank, i-firstBox);

							// update keepers
							_FF_sendLocalOffsets_keeper[iOct][dest] = _FF_sendLocalOffsets_counter[iOct][dest];
							_Infos_sendLocalOffsets_keeper[iOct][dest] = _Infos_sendLocalOffsets_counter[iOct][dest];
						} 
					}
				}
			}
		}
	}	
	//pthread_mutex_unlock(&mutex); 
}

void Gaspi_FF_communicator::send_chunk(int iOct, int dest, int level, int nbTerms, int nbInfos, bool levelIsFinished)
{
	int indexToC = -1;
	int octree_offset = iOct * _wsize;

	flush_queues(_nbQueues);
			
	// local FF send offset
	gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
	gaspi_offset_t accumul_offset = _FF_sendRemoteOffsets_counter[iOct][dest] * sizeof(complex);
	gaspi_offset_t local_offset_ff = local_dest_offset + accumul_offset;
	gaspi_queue_id_t queue=0;
	
	// remote offset
	gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
	gaspi_offset_t remote_offset_ff = remote_sender_offset + accumul_offset;
	_FF_sendRemoteOffsets_counter[iOct][dest] += nbTerms;
	gaspi_size_t qty= nbTerms * sizeof(complex);
	
	
	// envoi des DATA sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_FF_sendBuf_seg_id,					// local seg ID
			local_offset_ff,						// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,						// remote offset
			qty,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
		
	// notify ID and VALUE
	// cout <<"max comms : " << _max_comm << endl;
	// cout << "msgID : " << _msgID << endl;
	// gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
	
	gaspi_notification_id_t notifyID = (_rank * _max_comm) + _msgID;
	_msgID = (_msgID + 1) % _max_comm;
	gaspi_notification_t notifyValue = (level << 16) | nbInfos;
	
	gaspi_offset_t remote_offset = (_Infos_sendRemoteOffsets[octree_offset + dest] +_Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	gaspi_offset_t local_offset = (_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int);
	_Infos_sendRemoteOffsets_counter[iOct][dest] += nbInfos;
	qty = nbInfos * sizeof(int);
	
	/*printf("Src %d Dest %d Ioct %d level %d nbInfos %d notifyId %d info %d local offset %d remote offset %d\n", 
		_rank, dest, iOct, level, nbInfos, notifyID, _Infos_sendbuffer[local_offset/sizeof(int)],
		local_offset/sizeof(int), remote_offset/sizeof(int)); fflush(stdout);*/
	//if (level == 4)
	{
		/*debug("infos_to_" + convert(dest), convert(_Infos_sendbuffer[local_offset/sizeof(int)]) + " at " + convert(remote_offset/sizeof(int)) + " : " + 
			convert(_Infos_sendRemoteOffsets[octree_offset + dest] ) + " " + convert(_Infos_sendRemoteOffsets_counter[iOct][dest]) + " ff : " + convert(remote_offset_ff/sizeof(complex))
				+	" : " + convert(_FF_sendRemoteOffsets_counter[iOct][dest]) + " value : " + to_string( _FF_sendBuffer[local_offset_ff/sizeof(complex)].re) +
				" " + to_string(_FF_sendBuffer[local_offset_ff/sizeof(complex)].im) );				*/
		
		/*debug("infos_to_" + convert(dest), convert(_Infos_sendbuffer[local_offset/sizeof(int)]) + " value : " + to_string( _FF_sendBuffer[local_offset_ff/sizeof(complex)].re) +
				" " + to_string(_FF_sendBuffer[local_offset_ff/sizeof(complex)].im));*/
		
		/*debug("infos_to_" + convert(dest), to_string( _FF_sendBuffer[local_offset_ff/sizeof(complex)].re) +
				" " + to_string(_FF_sendBuffer[local_offset_ff/sizeof(complex)].im));*/
	}
	
	// envoi des INFO et NOTIFICATION
	SUCCESS_OR_DIE(
		gaspi_write_notify( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_offset,						// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_offset,						// remote offset
			qty,								// size of data to write
			notifyID,							// remote notif ID
			notifyValue,						// value of the notif to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
	/*debug(
		"gaspi_to_" + convert(dest), 
		"notifyID : " + to_string((int)notifyID) + ", src : " + to_string((uint32_t)notifyID/_max_comm) +  
		", notifyValue : " + to_string((int)notifyValue) + ", niv : " + to_string((int)notifyValue >> 16) + ", nbNodes : " + to_string((int)notifyValue & 0x00001111)
	);*/
	
	// Une fois tout fini, reset ou non le tableau d'offsets
	
	//cout << "LEVEL : " << level << " IncLevcom : " << _incLevcom << " levelIsFinished : " << levelIsFinished << " levcom : " << _levcom[iOct] + indexToC << endl;
	if (levelIsFinished)
	{
		if(_incLevcom)
		{
			if (level == _levcom[iOct] + indexToC ) 
			{
				_FF_sendLocalOffsets_counter[iOct][dest] = 0;
				_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
				_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
				_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
				_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
				_Infos_sendLocalOffsets_keeper[iOct][dest] = 0;
			}
		}
		else // cas allreduce sur levcom
		{
			if (level == _levcom[iOct] + 1 + indexToC)
			{
				_FF_sendLocalOffsets_counter[iOct][dest] = 0;
				_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
				_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
				_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
				_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
				_Infos_sendLocalOffsets_keeper[iOct][dest] = 0;
			}
		}
	}
	
	/*SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));*/
							
	/*if(level==4 && _rank==1 && dest==0 && iOct==2)
	{
		for (int i=0; i<counter; i++)
			debug("indexes_sent_to_0", itoa(_Infos_sendbuffer[_Infos_sendLocalOffsets[iOct][dest]]+i));
	}*/
	/*for (int i=0; i<counter; i++)
		accumulMSG("iOct" + itoa(iOct) + " to " + itoa(dest) + " " + itoa(_Infos_sendbuffer[_Infos_sendLocalOffsets[iOct][dest]]+i));*/
}

/*********************************************************************************
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
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_FF_recvBuf_seg_id,
					notif_offset,				// surveille les notifications depuis notif_offset
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);
						
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
			updateFarFields(sender, level, ff, iOct);
		}
		t_end = MPI_Wtime();
		add_time_sec("FF_read_from_buffer", t_end - t_begin);
		add_time_sec("FF_buffering", t_end - t_begin);
	}
}
*********************************************************************************/

/*********************************************************************************/
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
			ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[q:__nstnsp];
		}
	}

	t_end = MPI_Wtime();
	add_time_sec("Boucle_read_from_buffer", t_end - t_begin);
}
/*********************************************************************************/

void Gaspi_FF_communicator::updateFarFieldsFromInfos(int src, int level, complex * ff, int counter, int iOct)
{
	int indexToC = -1;
	int octree_offset = iOct * _wsize;

	// calcul du level_offset !
	int info_level_offset = 0;
	int ff_level_offset = 0;
	for (int i=_nivterm[iOct]-1; i>level; i--)
	{
		info_level_offset += _Expect[iOct][src][i];
		ff_level_offset += _Expect[iOct][src][i]*_nst[iOct][i]*_nsp[iOct][i];
	}
	
	int info_offset = _Infos_recvOffsets[octree_offset + src] + info_level_offset; // attention offset du level !!!
	int ff_offset = _FF_recvOffsets[octree_offset + src] + ff_level_offset;

	int firstBox = _stop_recv[iOct][src][level];
	
	int index;
	int recvDataIdx;

	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	int __endlev = _endlev[iOct][level]+indexToC;
	int __fnivnextlev = _fniv[iOct][level+1];
	
	int nbCellsToRecv = _frecv[iOct][src+1] - _frecv[iOct][src] + 1;

	int ff_index;
	for (int i=info_offset; i<info_offset+counter; i++)
	{
		index = _Infos_recvbuffer[i];
		int cellIDf = _recv[iOct][firstBox+index];
		int cellID = cellIDf + indexToC;
			/*if(level==4)
				debug("recv", "iOct : " + itoa(iOct) + " from " + itoa(src) + " level " + itoa(level) + " " + itoa(cellID)
				+ " index " + itoa(index) + " i du for : " + itoa(i) + " starting at offset : " + itoa(info_offset));*/

		// ranger la boîte
		int ff_index = (i-info_offset) * __nstnsp;
		int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);
		
		ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+ff_index:__nstnsp];
			
	/*	if(level==4)
		{
			debug("indexes_received_from_" + convert(src), convert(index) + " at offset " + convert (i) + 
				" value : " + to_string(_FF_recvBuffer[ff_offset+ff_index].re) + " " + to_string(_FF_recvBuffer[ff_offset+ff_index].im) + 
				" ff__off : " + convert((ff_offset+ff_index)/sizeof(complex)));
		}*/
	}	
}

void Gaspi_FF_communicator::updateFarFieldChunksFromInfos(int src, int level, complex * ff, int nbNodes, int iOct)
{

	/*debug("enter_update_ff", 
		"src : " + to_string(src) + 
		", niv : " + to_string(level) + 
		", nbNodes : " + to_string(nbNodes) +
		", iOct : " + to_string(iOct)
		);*/

	int indexToC = -1;
	int octree_offset = iOct * _wsize;

	// calcul du level_offset !
	int info_level_offset = 0;
	int ff_level_offset = 0;
	for (int i=_nivterm[iOct]-1; i>level; i--)
	{
		info_level_offset += _Expect[iOct][src][i];
		ff_level_offset += _Expect[iOct][src][i]*_nst[iOct][i]*_nsp[iOct][i];
	}
	
	//printf("[%d] updateFF -> IN, already received %d\n", _rank, _Received[iOct][level]); fflush(stdout);	
	//printf("[%d] update, already received %d nodes from %d level %d\n", _rank, _Received[iOct][level], src, level);	fflush(stdout);

	int info_offset = _Infos_recvOffsets[octree_offset + src] + info_level_offset + _Received[iOct][src][level]; // attention offset du level !!!
	int ff_offset = _FF_recvOffsets[octree_offset + src] + ff_level_offset + _Received[iOct][src][level]*_nst[iOct][level]*_nsp[iOct][level];

	int firstBox = _stop_recv[iOct][src][level];
	
	int index;
	int recvDataIdx;

	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	int __endlev = _endlev[iOct][level]+indexToC;
	int __fnivnextlev = _fniv[iOct][level+1];

	int ff_index;
	for (int i=info_offset; i<info_offset+nbNodes; i++)
	{
		index = _Infos_recvbuffer[i];
		int cellIDf = _recv[iOct][firstBox+index];
		int cellID = cellIDf + indexToC;

		// ranger la boîte
		int ff_index = (i-info_offset) * __nstnsp;
		int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);
		
		/*debug("going_to_update_ff", 
		"src : " + to_string(src) + 
		", niv : " + to_string(level) + 
		", nbNodes : " + to_string(nbNodes) +
		", iOct : " + to_string(iOct)
		);*/
		ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+ff_index:__nstnsp];
		/*debug("updated_ff", 
		"src : " + to_string(src) + 
		", niv : " + to_string(level) + 
		", nbNodes : " + to_string(nbNodes) +
		", iOct : " + to_string(iOct)
		);*/
		
		/*if(level==4)
		{
			debug("indexes_received_from_" +convert(src), convert(index) + " at offset " + convert (i) + 
				" value : " + to_string(_FF_recvBuffer[ff_offset+ff_index].re) + " " + to_string(_FF_recvBuffer[ff_offset+ff_index].im) + 
				" ff_off : " + convert((ff_offset+ff_index)/sizeof(complex)));
		}*/
	}
	//printf("[%d] updateFF -> OUT, already received %d\n", _rank, _Received[iOct][level]); fflush(stdout);
	/*debug("exit_update_ff", 
		"src : " + to_string(src) + 
		", niv : " + to_string(level) + 
		", nbNodes : " + to_string(nbNodes) +
		", iOct : " + to_string(iOct)
		);*/
}

/*******************************
 * 		TENTATIVE RECUP PERF
 *******************************/


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
	
	// prepare _SendBuffer
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{
			t_begin = MPI_Wtime();
			// box range to send, FROM SEND ARRAY
			/*int firstBoxToSendIDX= _fsend[iOct][dest] + indexToC;
			int lastBoxToSendIDX = _fsend[iOct][dest+1]-1 + indexToC;
			
			// cell range from the current level
			int beginlevel = _endlev[iOct][level-1]+1 + indexToC; 
			int endlevel   = _endlev[iOct][level] + indexToC;
			
			//int count = 0;
			int count = _count_send[iOct][dest][level];
			//int q = _LocalSendOffsets[octree_offset + dest] - 1 + _offsetKeeper[iOct][dest];
			int q = _FF_sendLocalOffsets[iOct][dest] - 1 + _FF_sendLocalOffsets_counter[iOct][dest];

			
			// remplissage du buffer d'envoi, par ordre décroissant
			//#pragma omp parallel for private(q)
			for (int k = lastBoxToSendIDX; k>=firstBoxToSendIDX; k--)
			{
				int cellID = _send[iOct][k] + indexToC;
				
				if ( (cellID >= beginlevel) && (cellID <= endlevel))	// todo : stop earlier, cells are ranged in ascending order
				{
					// calcule le bon index pour p (= offset, dans le tableau ff)
					int p=_fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];
					p = p + indexToC;
					
					// ecrit dans SendBuffer les nst(niveau) x nsp(niveau) complexes, de cellID, depuis FF
					for (int st=0; st<_nst[iOct][level]; st++)
					{
						for (int sp=0; sp<_nsp[iOct][level];  sp++)
						{
							p++;
							q++;
							_FF_sendBuffer[q]=ff[p];
						}
					}
				}
			}*/

			int q = _FF_sendLocalOffsets[iOct][dest] - 1 + _FF_sendLocalOffsets_counter[iOct][dest];
			int q0 = q + 1; 

			int count =  _count_send[iOct][dest][level];
			int lastBox =  _start_send[iOct][dest][level];
			int firstBox = _stop_send[iOct][dest][level];
			int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
			
			if (count>0)
			{
				//#pragma omp parallel for private(q)
				for (int k=lastBox; k>=firstBox; k--)
				{
					int cellID = _send[iOct][k] + indexToC;
					int p=_fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];
					q = q0 + (lastBox-k)*__nstnsp;
					_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp];
				}
			}
			
			t_end = MPI_Wtime();
			add_time_sec("FF_write_to_buffer", t_end - t_begin);
			t_begin = MPI_Wtime();
			
			// if something to send
			if (count > 0)
			{
				// local offset
				gaspi_offset_t local_dest_offset = /*_LocalSendOffsets[octree_offset + dest]*/_FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
				gaspi_offset_t level_offset = _FF_sendLocalOffsets_counter[iOct][dest]/*_offsetKeeper[iOct][dest]*/ * sizeof(complex);
				gaspi_offset_t local_offset = local_dest_offset + level_offset;
				gaspi_queue_id_t queue=0;
				
				// remote offset
				gaspi_offset_t remote_sender_offset = /*_RemoteSendOffsets*/_FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t remote_offset = remote_sender_offset + level_offset;

				// update offset, per destinatary
				/*_offsetKeeper[iOct][dest]*/_FF_sendLocalOffsets_counter[iOct][dest] += count * _nst[iOct][level] * _nsp[iOct][level];

				/*int rankMultiple = level;
				gaspi_notification_id_t notif_offset = _wsize * rankMultiple;
				gaspi_notification_id_t notifyID = notif_offset + _rank;
				TODO : CORRIGER CAR ARBRES PEUVENT AVOIR DES HAUTEURS DIFFERENTES*/
				gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
				gaspi_size_t qty= count * _nst[iOct][level] * _nsp[iOct][level] * sizeof(complex);
				
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
			}
			
			t_end = MPI_Wtime();
			add_time_sec("GASPI_SEND_write_notify", t_end - t_begin);
			accumul = accumul + (t_end - t_begin);
		}

		// cas include levcom
		if(_incLevcom)
		{
			if (level == _levcom[iOct] + indexToC)
			{
				/*_offsetKeeper[iOct][dest]*/ _FF_sendLocalOffsets_counter[iOct][dest] = 0;
				// RAZ des compteurs
				_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
				_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
				_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
			}
		}
		else // cas allreduce sur levcom
		{
			if (level == _levcom[iOct] + 1 + indexToC)
				/*_offsetKeeper[iOct][dest]*/ _FF_sendLocalOffsets_counter[iOct][dest] = 0;
		}	
	}
	add_time_sec("GASPI_FF_sendrecv", accumul);
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

	/*int rankMultiple = level;
	gaspi_notification_id_t notif_offset = _wsize * rankMultiple;*/
	gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
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
					_FF_recvBuf_seg_id,
					notif_offset,				// surveille les notifications depuis 0
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);

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
			updateFarFields(sender, level, ff, iOct);
		}
		t_end = MPI_Wtime();
		add_time_sec("FF_read_from_buffer", t_end - t_begin);
	}
}
/*
void Gaspi_FF_communicator::updateFarFields(int src, int level, complex * ff, int iOct)
{		
	int indexToC = -1;
	int k = level + 1;
	int octree_offset = iOct * _wsize;

	// calcul de l'offset	
	int levelOffset = 0;
	for (int i=_nivterm[iOct] + indexToC; i>level; i--)
	{
		levelOffset += _nst[iOct][i] * _nsp[iOct][i] * _Expect[iOct][src][i];
	}
	int q = _FF_recvOffsets[octree_offset + src]-1; // se mettre 1 case avant l'index à lire
	q += levelOffset; 	

	// box range to recv, FROM RECV ARRAY
	int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
	int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;

	//#pragma omp parallel for private (q) REACTIVER, car TEMPORAIREMENT DESACTIVE
    for (int j=lastBoxToRecvIdx; j>=firstBoxToRecvIdx; j--)
	{
		int cellID = _recv[iOct][j] + indexToC;
	
		// test if box belongs to level 
		if ( (cellID > _endlev[iOct][level-1]+indexToC) && (cellID <= _endlev[iOct][level]+indexToC) )
		{
			// update far fields
			int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];
			p = p + indexToC;
			
			for (int st=0; st<_nst[iOct][level]; st++)
			{
				for (int sp=0; sp<_nsp[iOct][level]; sp++)
				{
					p++;
					q++;
					ff[p] = ff[p] + _FF_recvBuffer[q];
				}
			}			
		}
	}
}*/
