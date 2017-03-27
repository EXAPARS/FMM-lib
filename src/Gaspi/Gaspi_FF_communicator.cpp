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
	// PREPARATION _SendBuffer
	for (int dest=0; dest<_wsize; dest++)
	{
		int niv = _nivterm[iOct];
		if (dest != _rank)
		{
			t_begin = MPI_Wtime();
			// box range to send, FROM SEND ARRAY
			int firstBoxToSendIDX= _fsend[iOct][dest] + indexToC;
			int lastBoxToSendIDX = _fsend[iOct][dest+1]-1 + indexToC;
			
			int q = _LocalSendOffsets[octree_offset + dest] - 1 /*+ _offsetKeeper[iOct][dest]*/;
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
								_SendBuffer[q]=bufsave[p + indexToC];
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
								_SendBuffer[q]=ff[p + indexToC];
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
				local_offset  = _LocalSendOffsets[octree_offset + dest] * sizeof(complex);
				remote_offset = _RemoteSendOffsets[octree_offset + dest] * sizeof(complex);
				qty = _nb_send[iOct][dest] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_seg_SendBuffer_id,					// local seg ID
						local_offset,						// local offset
						dest,									// receiver rank
						_seg_RecvBuffer_id,					// remote seg ID
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
						_seg_RecvBuffer_id,			// seg
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
					_seg_RecvBuffer_id,
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
					_seg_RecvBuffer_id, 
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
	int qp = _RecvOffsets[octree_offset + src]-1; // se mettre 1 case avant l'index à lire
		
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
					ff[pp + indexToC] = ff[pp + indexToC] + _RecvBuffer[qp];
				}
			}
		}
	}
}

/* ------------------------------------------------- */
/*               GASPI MULTIMAT VERSION              */
/* ------------------------------------------------- */
Gaspi_FF_communicator::Gaspi_FF_communicator(int max_send, int max_recv, int incLevcom, int nbOct)
{
	// update class attributes
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);

	// update segments ids
	_seg_RecvBuffer_id = 7; 	// Receive Buffer 
	_seg_SendBuffer_id = 8; 	// Send Buffer 
    _seg_RecvOffsets_id = 9; 	// Index where to write on other ranks	
	_seg_RemoteSendOffsets_id = 10; // Index in global recv buffer, per RANK
	
	// update scalar class attributes
	_nbQueues = 1;
	_incLevcom = incLevcom;
	_nbOct = nbOct;
		
	// allocate array class attributes
	alloc_attributes();
	
	// create gaspi segments, and initialize them
	create_segments(max_send, max_recv);
	
}

void init_gaspi_ff(int max_send, int max_recv, int nbMat, int incLevcom, Gaspi_FF_communicator *& gCommFF)
{
    gCommFF = new Gaspi_FF_communicator(max_send, max_recv, incLevcom, nbMat);
}

void Gaspi_FF_communicator::alloc_attributes()
{
	_offsetKeeper = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_offsetKeeper[i] = new int[_wsize]();

	_LocalSendOffsets = new int[_wsize * _nbOct](); // Local SEND offsets x nb of octrees (=nb materials)

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

void Gaspi_FF_communicator::create_segments(int max_send, int max_recv)
{
 	// RECV segment
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_RecvBuffer_id,
			max_recv * sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_seg_RecvBuffer_id, &_ptr_seg_RecvBuffer);
	_RecvBuffer = (complex *)_ptr_seg_RecvBuffer;

	// SEND segment 	
 	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_SendBuffer_id,
			max_send * sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);	
	gaspi_segment_ptr(_seg_SendBuffer_id, &_ptr_seg_SendBuffer);
	_SendBuffer = (complex *)_ptr_seg_SendBuffer;
	
	// Reception Offsets
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_RecvOffsets_id,
			_nbOct * _wsize * sizeof(int),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
	);
	gaspi_segment_ptr(_seg_RecvOffsets_id, &_ptr_seg_RecvOffsets);
	_RecvOffsets = (int *)_ptr_seg_RecvOffsets;
	
	// Remote SEND offsets x nb of octrees (=nb materials)
	SUCCESS_OR_DIE(
		gaspi_segment_create(
			_seg_RemoteSendOffsets_id, 
			_nbOct * _wsize * sizeof(int),
			GASPI_GROUP_ALL, 
			GASPI_BLOCK, 
			GASPI_ALLOC_DEFAULT
		)
    );	
	gaspi_segment_ptr(_seg_RemoteSendOffsets_id, &_ptr_seg_RemoteSendOffsets);
	_RemoteSendOffsets = (int *)_ptr_seg_RemoteSendOffsets;	
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
     
    // initialize _RecvOffsets[octree] segment with values -> put in rank order
    int octree_offset = iOct * _wsize;
    for (int k=0; k<recvnode_sz; k++)
    {
		if (recvnode[k] > 0)
		{
			int from = recvnode[k]-1;
			if (from != _rank)
			{
				_RecvOffsets[octree_offset + from] = RecvBufferIndexPerROUND[k];
			}
			else
			{
				cerr << "[Create_remoteBufferIndexes] Fortran recvnode array should not contain the current rank !"; 
				exit(0);
			}
		}
    }
	_RecvOffsets[octree_offset + _rank] = -1;
	
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
					_seg_RecvOffsets_id,					// local seg ID
					local_offset,						// local offset
					i,									// receiver rank
					_seg_RemoteSendOffsets_id,				// remote seg ID
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
			_seg_RemoteSendOffsets_id,
			octree_offset,				            // surveille les notifications depuis le domaine
			_wsize,						            // en surveille wsize
			&new_notif_id,
			GASPI_BLOCK) == GASPI_SUCCESS)
		{

			// get notification value, and reset
			gaspi_notify_reset(_seg_RemoteSendOffsets_id, new_notif_id, &new_notif_val);
		
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
	
	_RemoteSendOffsets[octree_offset + _rank] = -1;
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
    
    // Prepare _sendBufferIndexes organised by process RANK, and MATERIAL DOMAIN
    int octree_offset = iOct * _wsize;
    for (int k=0; k<sendnode_sz; k++)
    {
		if(sendnode[k] > 0)
		{
			int to = sendnode[k]-1;
			if (to != _rank)
			{
				_LocalSendOffsets[octree_offset + to] = SendBufferIndexPerROUND[k];
			}
			else
			{
				cerr << "Error in Gaspi_m2l_communicator::init_dataToSendIndexes.\nFortran sendnode array should not contain the current rank !";
				exit(0);
			}
		}
    }  
	_LocalSendOffsets[octree_offset + _rank] = -1;
		
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
		_stop_send[iOct][j] = new int[_nivterm[iOct]]();
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
	
	// prepare _SendBuffer
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{
			t_begin = MPI_Wtime();
			int q = _LocalSendOffsets[octree_offset + dest] - 1 + _offsetKeeper[iOct][dest];
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
					/*
					int i;
					#pragma simd
					for (i=0; i<__nstnsp; i++)
						_SendBuffer[q+i] = ff[p+i];*/
					_SendBuffer[q:__nstnsp] = ff[p:__nstnsp]; 
				}
				t_end = MPI_Wtime();
			}
			add_time_sec("FF_write_to_buffer", t_end - t_begin);
			add_time_sec("FF_buffering", t_end - t_begin);
		}
	}

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
				gaspi_offset_t local_dest_offset = _LocalSendOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t level_offset = _offsetKeeper[iOct][dest] * sizeof(complex);
				gaspi_offset_t local_offset = local_dest_offset + level_offset;
				gaspi_queue_id_t queue=0;
				
				// remote offset
				gaspi_offset_t remote_sender_offset = _RemoteSendOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t remote_offset = remote_sender_offset + level_offset;


				// update offset, per destinatary
				_offsetKeeper[iOct][dest] += counter * _nst[iOct][level] * _nsp[iOct][level];

				//int rankMultiple = level;
				//TODO : CORRIGER CAR ARBRES PEUVENT AVOIR DES HAUTEURS DIFFERENTES
				gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
				gaspi_size_t qty= counter * _nst[iOct][level] * _nsp[iOct][level] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_seg_SendBuffer_id,			// local seg ID
						local_offset,						// local offset
						dest,								// receiver rank
						_seg_RecvBuffer_id,			// remote seg ID
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
			if (level == _levcom[iOct] + indexToC)
				_offsetKeeper[iOct][dest] = 0;
		}
		else // cas allreduce sur levcom
		{
			if (level == _levcom[iOct] + 1 + indexToC)
				_offsetKeeper[iOct][dest] = 0;
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
					_seg_RecvBuffer_id,
					notif_offset,				// surveille les notifications depuis 0
					_wsize,						// en surveille wsize
					&new_notif_id,
					GASPI_BLOCK
				)
			);

			SUCCESS_OR_DIE(
				gaspi_notify_reset(
					_seg_RecvBuffer_id, 
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
	int q = _RecvOffsets[octree_offset + src]-1; // se mettre 1 case avant l'index à lire
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
	
	bool start = (_start_recv[iOct][src][level]>=0);
	bool stop  = (_stop_recv[iOct][src][level]>=0);
	int firstBox = _stop_recv[iOct][src][level]; 
	int lastBox  = _start_recv[iOct][src][level];
	
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
				ff[p0+i].re = ff[p0+i].re + _RecvBuffer[q+i].re;
				ff[p0+i].im = ff[p0+i].im + _RecvBuffer[q+i].im;
			}*/
			ff[p0:__nstnsp] = ff[p0:__nstnsp] + _RecvBuffer[q:__nstnsp];
		}
	}

	t_end = MPI_Wtime();
	add_time_sec("Boucle_read_from_buffer", t_end - t_begin);
}

/* void Gaspi_FF_communicator::updateFarFields(int src, int level, complex * ff, int iOct) // optimized OUT
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
	int q = _RecvOffsets[octree_offset + src]-1; // se mettre 1 case avant l'index à lire
	q += levelOffset; 	

	// box range to recv, FROM RECV ARRAY
	int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
	int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;

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
					ff[p] = ff[p] + _RecvBuffer[q];
				}
			}
		}
	}
}*/

/*
void Gaspi_FF_communicator::send_ff_level(int level, complex * ff, int iOct) // optimized OUT
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
			int firstBoxToSendIDX= _fsend[iOct][dest] + indexToC;
			int lastBoxToSendIDX = _fsend[iOct][dest+1]-1 + indexToC;
			
			// cell range from the current level
			int beginlevel = _endlev[iOct][level-1]+1 + indexToC; 
			int endlevel   = _endlev[iOct][level] + indexToC;
			
			int count = 0;
			int q = _LocalSendOffsets[octree_offset + dest] - 1 + _offsetKeeper[iOct][dest];
			
			// remplissage du buffer d'envoi, par ordre décroissant
			for (int k = lastBoxToSendIDX; k>=firstBoxToSendIDX; k--)
			{
				int cellID = _send[iOct][k] + indexToC;
				
				if ( (cellID >= beginlevel) && (cellID <= endlevel))	// todo : stop earlier, cells are ranged in ascending order
				{
					count++;
					
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
							_SendBuffer[q]=ff[p];
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
				gaspi_offset_t local_dest_offset = _LocalSendOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t level_offset = _offsetKeeper[iOct][dest] * sizeof(complex);
				gaspi_offset_t local_offset = local_dest_offset + level_offset;
				gaspi_queue_id_t queue=0;
				
				// remote offset
				gaspi_offset_t remote_sender_offset = _RemoteSendOffsets[octree_offset + dest] * sizeof(complex);
				gaspi_offset_t remote_offset = remote_sender_offset + level_offset;


				// update offset, per destinatary
				_offsetKeeper[iOct][dest] += count * _nst[iOct][level] * _nsp[iOct][level];

				//int rankMultiple = level;
				//gaspi_notification_id_t notif_offset = _wsize * rankMultiple;
				//gaspi_notification_id_t notifyID = notif_offset + _rank;
				//TODO : CORRIGER CAR ARBRES PEUVENT AVOIR DES HAUTEURS DIFFERENTES
				gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
				gaspi_size_t qty= count * _nst[iOct][level] * _nsp[iOct][level] * sizeof(complex);
				
				SUCCESS_OR_DIE(
					gaspi_write_notify( 
						_seg_SendBuffer_id,			// local seg ID
						local_offset,						// local offset
						dest,								// receiver rank
						_seg_RecvBuffer_id,			// remote seg ID
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
				_offsetKeeper[iOct][dest] = 0;
		}
		else // cas allreduce sur levcom
		{
			if (level == _levcom[iOct] + 1 + indexToC)
				_offsetKeeper[iOct][dest] = 0;
		}
	}
	add_time_sec("GASPI_FF_sendrecv", accumul);
}*/
