#include "Gaspi_FF_communicator.hpp"
#include "../Tools/fmm_tools.hpp"
#include <unistd.h>
using namespace std;

static pthread_mutex_t mutex;
/* ---------------------------------------------- */
/*               GASPI BULK VERSION               */
/* ---------------------------------------------- */

void Gaspi_FF_communicator::exchangeFFBulk(complex * bufsave, complex * ff, int iOct)
{

	int indexToC = -1;
	//double t_begin_comm, t_end_comm, t_begin, t_end, accumul_comm;
	//t_begin = MPI_Wtime();
	//accumul_comm = 0;
		
	int octree_offset = iOct * _wsize;
	

	// PREPARATION _FF_sendBuffer
	for (int dest=0; dest<_wsize; dest++)
	{
		int niv = _nivterm[iOct];
		if (dest != _rank)
		{
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
				
				if (_incLevcom)
				{
					if ((_levcom[iOct]<=3 || (_levcom[iOct] == 4 && _nivterm[iOct]>8))  && (cellID <  _endlev[iOct][_levcom[iOct]-1 + indexToC] + indexToC))
						todo = false;
					if ((_levcom[iOct]>4  || (_levcom[iOct] == 4 && _nivterm[iOct]<=8)) && (cellID <= _endlev[iOct][_levcom[iOct]-1 + indexToC] + indexToC))
						todo = false;
				}
				else
				{
					if ((_levcom[iOct]<=3 || (_levcom[iOct] == 4 && _nivterm[iOct]>8))  && (cellID <  _endlev[iOct][_levcom[iOct] + indexToC] + indexToC))
						todo = false;
					if ((_levcom[iOct]>4  || (_levcom[iOct] == 4 && _nivterm[iOct]<=8)) && (cellID <= _endlev[iOct][_levcom[iOct]-1 + indexToC] + indexToC))
						todo = false;
				}

				if (todo)
				{
					// place k pour que k = niveau de la cellule à envoyer
					while(cellID <= _endlev[iOct][niv-1 + indexToC] + indexToC)
						niv--;
					
					// cas 1 : 
					// si l est commun à +ieurs noeuds
					// communique les données sauvegardées dans BUFSAVE
					if (_codech[iOct][cellID] > 1)
					{
						p = _codech[iOct][cellID]-2;
						
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
			
	// ENVOIS
	gaspi_offset_t local_offset;
    gaspi_offset_t remote_offset;
   	gaspi_notification_id_t notif_offset = _wsize * 2;    
    gaspi_notification_id_t notifyID = notif_offset + _rank;
   // gaspi_queue_id_t queue=0;
    gaspi_size_t qty;
    
	flush_queues(_nbQueues);

	
	  
    for (int dest=0; dest<_wsize; dest++)
    {
		if ( dest != _rank ) // if not current rank
		{
			if (_nb_send[iOct][dest] > 0) // if something to send
			{
				local_offset  = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
				remote_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
				qty = _nb_send[iOct][dest] * sizeof(complex);
				
				
				//t_begin_comm = MPI_Wtime(); 
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
				//t_end_comm = MPI_Wtime();
				//accumul_comm = accumul_comm + (t_end_comm - t_begin_comm);
			}
			else // no data to send, notify anyway
			{
				//t_begin_comm = MPI_Wtime();
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
				//t_end_comm = MPI_Wtime();
				//accumul_comm = accumul_comm + (t_end_comm - t_begin_comm);
			}
		}
	}
	
	
	// RECEIVE
	int recvCpt = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int sender;
	//double t_begin_loop, t_end_loop;
	
	while (recvCpt < (_wsize-1))
	{
		//t_begin_comm = MPI_Wtime();
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
		//t_end_comm = MPI_Wtime();
		//accumul_comm = accumul_comm + (t_end_comm - t_begin_comm); 
		
		
		// test the notification value and update counter
		if (new_notif_val == SEND_DATA || new_notif_val == NO_DATA)
			recvCpt++;

		// update the far field array
		if (new_notif_val == SEND_DATA)
		{
			sender = new_notif_id - notif_offset;
			updateFarFields(sender, ff, iOct);
		}
	}
	
//	t_end = MPI_Wtime();
//~ #ifdef TIMING
	//~ add_time_sec("FF_sendrecv_comm", accumul_comm);
	//~ add_time_sec("FF_sendrecv_copy", t_end - t_begin - accumul_comm);
//~ #endif
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
		
		if (_incLevcom)
		{
			if ((_levcom[iOct]<=3 || (_levcom[iOct] == 4 && _nivterm[iOct]>8))  && (p < _endlev[iOct][_levcom[iOct]-1 + indexToC]))
				todo = false;
			if ((_levcom[iOct]>4  || (_levcom[iOct] == 4 && _nivterm[iOct]<=8)) && (p <= _endlev[iOct][_levcom[iOct]-1 + indexToC]))
				todo = false;
		}
		else
		{
			if ((_levcom[iOct]<=3 || (_levcom[iOct] == 4 && _nivterm[iOct]>8))  && (p < _endlev[iOct][_levcom[iOct] + indexToC]))
				todo = false;
			if ((_levcom[iOct]>4  || (_levcom[iOct] == 4 && _nivterm[iOct]<=8)) && (p <= _endlev[iOct][_levcom[iOct]-1 + indexToC]))
				todo = false;
		}
		
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
	
	// mutex
	pthread_mutex_init(&_mutex, NULL);
	pthread_mutex_init(&_mutex1, NULL);
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

/*** ---- ***/
	_ReceivedInfosCpt = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_ReceivedInfosCpt[i] = new int[_wsize](); 

	_ReceivedFFTermsCpt = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_ReceivedFFTermsCpt[i] = new int[_wsize](); 
		
	_SendInfosCpt = new int*[_nbOct]();
	for (int i=0; i<_nbOct; i++)
		_SendInfosCpt[i] = new int[_wsize](); 



	_ExpectSrc = new int*[_nbOct]();
	/*for (int i=0; i<_nbOct; i++)
		_ExpectSrc[i] = new int[_wsize]();*/

	_CountDest = new int*[_nbOct]();
	_CountDestTerms = new int*[_nbOct]();
	/*for (int i=0; i<_nbOct; i++)
		_CountDest[i] = new int[_wsize]();*/
		
/*** ---- ***/

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
	
	_max_packets = new int*[_nbOct](); // for each octree
	for (int i=0; i<_nbOct; i++)
		_max_packets[i] = new int[_wsize](); // for each dest
	
	_write_counters_ptr = new int*[_nbOct](); // for each octree
	for (int i=0; i<_nbOct; i++) 
		_write_counters_ptr[i] = new int[_wsize](); // for each dest
	
	_write_counters = new int**[_nbOct](); // for each octree
	for (int i=0; i<_nbOct; i++) 
		_write_counters[i] = new int*[_wsize](); // for each dest
	
#pragma omp parallel
{
	_nbThreads = omp_get_num_threads();
}

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

	//~ debug("gaspi_buffers", "recv_size : " + itoa(max_recv_terms * sizeof(complex)));
	//~ debug("gaspi_buffers", "send_size : " + itoa(max_send_terms * sizeof(complex)));

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
	
	/* --- tasks --- */
	fill_max_packets(iOct);
	alloc_and_fill_write_counters(iOct);
	
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
	// int indexToC = -1;
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
   // int octree_offset = iOct * _wsize;
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
	
	for (int i=0; i<_nbOct; i++)
		_ExpectSrc[iOct]   = new int[_wsize]();
	
	// fill _expectSrc, no level distinction
	for (int src=0; src<_wsize; src++)
	{
		for (int level=0; level <_nivterm[iOct]; level++)
		{
			_ExpectSrc[iOct][src] += _Expect[iOct][src][level];
		}
	}
}

void Gaspi_FF_communicator::fill_send_start_stop_count(int iOct)
{
	//double t_begin, t_end;
	//t_begin = MPI_Wtime();
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
				//int beginlevel = _endlev[iOct][level-1]+1 + indexToC; 
				//int endlevel   = _endlev[iOct][level] + indexToC;
				
				//int count = 0;
				
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
	//t_end = MPI_Wtime();
	//add_time_sec("fill_send_start_stop_count", t_end - t_begin);

/*	for (int i=0; i<_nbOct; i++)
		_ExpectSrc[iOct]   = new int[_wsize]();*/
	
	for (int i=0; i<_nbOct; i++)
		_CountDest[i] = new int[_wsize]();

	for (int i=0; i<_nbOct; i++)
		_CountDestTerms[i] = new int[_wsize]();
			
	// fill _CountDest, no level distinction
	for (int dest=0; dest<_wsize; dest++)
	{
		for (int level=0; level <_nivterm[iOct]; level++)
		{
			_CountDest[iOct][dest] += _count_send[iOct][dest][level];
			_CountDestTerms[iOct][dest] += _count_send[iOct][dest][level] * _nsp[iOct][level] * _nst[iOct][level];
		}
	}

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
	//double t_begin, t_end;
	//t_begin = MPI_Wtime();
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
				
				//int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
				int __endlev = _endlev[iOct][level]+indexToC;
				int __endlevprec = _endlev[iOct][level-1]+indexToC;
				long * __recv = _recv[iOct];
				//int __fnivnextlev = _fniv[iOct][level+1];
				//int __nst = _nst[iOct][level];
				//int __nsp = _nsp[iOct][level];
	
				// box range to send, FROM SEND ARRAY
				int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
				int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;
				
				
				// find box range in the targeted level
				bool start = false;
				bool stop = false;
				//int firstBox, lastBox;
				
				for (int j=lastBoxToRecvIdx; j>=firstBoxToRecvIdx; j--)
				{
					int cellID = __recv[j] + indexToC;

					if (!start)
					{
						if ( (cellID > __endlevprec) && (cellID <= __endlev) )
						{
							start = true;
							//lastBox = j;
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
								//firstBox = j+1;
								_stop_recv[iOct][src][level] = j+1;
							}
						}
					}
				}
				if (!start)
				{
					//lastBox =  -1;
					//firstBox = -1;
					_start_recv[iOct][src][level]= -1;
					_stop_recv[iOct][src][level] = -1;
				}
				
				if (start)
				{
					if(!stop)
					{
						//firstBox = firstBoxToRecvIdx;
						_stop_recv[iOct][src][level] = firstBoxToRecvIdx;
					}
				}
			}
		}
	}
	//t_end = MPI_Wtime();
	//add_time_sec("fill_send_start_stop_count", t_end - t_begin);
	
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

void Gaspi_FF_communicator::fill_max_packets(int iOct)
{
	//int indexToC = -1;
		
	for (int dest=0; dest<_wsize; dest++) // for each dest
	{
		//_max_packets[iOct][dest] = new int[_nivterm[iOct]](); 
		for (int level= 0; level<_nivterm[iOct]; level++)// for each level
		{
			//int nbCells = _count_send[iOct][dest][level];
			//int nbTerms = nbCells * _nst[iOct][level] * _nsp[iOct][level] * sizeof(complex);
			int nbTerms = _CountDestTerms[iOct][dest] * sizeof(complex);
			_max_packets[iOct][dest]/*[level]*/ = nbTerms / _gaspiChunkSize + 1;
		}
	}
}

void Gaspi_FF_communicator::alloc_and_fill_write_counters(int iOct)
{
	//int indexToC = -1;
	
	// alloc the pointers
	/**for (int dest=0; dest<_wsize; dest++) // for each dest
	{
		_write_counters_ptr[iOct][dest] = new int[_nivterm[iOct]]();
	}**/
	
	// alloc the counters
	for (int dest=0; dest<_wsize; dest++) // for each dest
	{
		///_write_counters[iOct][dest] = new int*[_nivterm[iOct]](); // for each level
		
		///for (int level= 0; level<_nivterm[iOct]; level++) 
		///{
			int sz = _max_packets[iOct][dest];
			_write_counters[iOct][dest] = new int[sz]();
		///}
	}
}

/*
 * Cette fonction remplit les buffers d'envoi de FF et Infos correspondantes 
 * de manière parallèlle à l'intérieur de tâches
 * L'envoi Gaspi est effectué lorsque le buffer est prêt, et par une seule tâche
 */
void Gaspi_FF_communicator::send_task_ff_level(int level, complex * ff, int iOct, int start, int stop)
{
	int indexToC = -1;
	
	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
	// for each dest
	for (int dest=0; dest<_wsize; dest++)
	{								
		if (dest != _rank)
		{
			//nb of cells to be sent
			int counter =  _count_send[iOct][dest][level];
			
			// if something to send to this dest
			if (counter>0)
			{
				int lastBox  = _start_send[iOct][dest][level];	// the big one
				int firstBox = _stop_send[iOct][dest][level];	// the small one

				// traverses the cells to send
				for (int i= firstBox; i <= lastBox; i++)
				{					
					int cellIDf = _send[iOct][i];
					
					// if the cell is in the task's nodes range
					if (cellIDf >= start && cellIDf <= stop)
					{							
						int cellID = cellIDf + indexToC;
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];	// cell's corresponding ff in Fortran						
						int chunkIsFull = 0;
						
						// prend le VERROU
						pthread_mutex_lock(&mutex); 
						
							// lit les offsets
							int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];						// in FF senbuffer
							int idxInfos =_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];
							// update les offsets
							_FF_sendLocalOffsets_counter[iOct][dest] += _nst[iOct][level]*_nsp[iOct][level];
							_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
							// vérifie si le chunk est full
							if (_Infos_sendLocalOffsets_counter[iOct][dest] == _count_send[iOct][dest][level])  chunkIsFull = 1;	// dernier node à envoyer du level
						


						// write FF and Infos
						///_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp]; 
						for (int i=0; i<__nstnsp; i++)
							_FF_sendBuffer[q+i] = ff[p+i];												
						
						_Infos_sendbuffer[idxInfos] = i-firstBox;	//!\ même danger, ts les threads ont-ils fini d'écrire ?

						// rend le VERROU
						pthread_mutex_unlock(&mutex); /// il faudrait s'assurer que tous les threads aient fini d'écrire

						if (chunkIsFull)
						{	
							_Infos_sendLocalOffsets_counter[iOct][dest]=0;
							send_chunk(iOct, dest, level);
						}
					}
				}
			}
		}
	}
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
			local_offset_ff,					// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,					// remote offset
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
}

void Gaspi_FF_communicator::recv_task_ff_level(int level, complex * ff, int iOct)
{
	double t_begin, t_end, accumul;
	accumul = 0;
	
	// compute nb of messages to receive at this level
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected++;
		}
	}

	// 
	gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int recvCpt = 0;
	int sender;
	int counter;
	
	while (recvCpt < nbRecvExpected)
	{
		t_begin = MPI_Wtime();
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

		t_end= MPI_Wtime();
		accumul = accumul + (t_end - t_begin);

		// test the notification value and update
		{
			recvCpt++;
			sender = new_notif_id - notif_offset;
			counter = new_notif_val;
			
			//debug("infoBase", "---------------- received : " + itoa(counter) + " from : " + itoa(sender) + " level : " + itoa(level));
			//debug("infoBase", "expected : " + itoa(_Expect[iOct][sender][level]));										
			// new decoding
			//debug("infoBase", "notifyID  : " + itoa(new_notif_id));
			
			updateFarFieldsFromInfos(sender, level, ff, counter, iOct);
		}
	}
	//add_time_sec("FF_recv_wait", accumul);
}

void Gaspi_FF_communicator::recv_task_ff_level_rest(int level, complex * ff, int iOct)
{
	double t_begin, t_end, accumul;
	accumul = 0;
	
	// compute nb of nodes to receive at this level
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			//nbRecvExpected++;
			nbRecvExpected += _Expect[iOct][i][level];
		}
	}
	
	// compute nb of already received nodes
	int nbReceived = 0;
	for (int src=0; src<_wsize; src++)
	{
		if (_Received[iOct][src][level]>0)
		{
			nbReceived += _Received[iOct][src][level];
		}
	}	

	debug("wait", "level : " + itoa(level) + ", received : " + itoa(nbReceived) + " expected : " + itoa(nbRecvExpected));
	// 
	gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	
	//int recvCpt = 0;
	int sender;
	int counter;
	
	while (nbReceived < nbRecvExpected)
	{
		t_begin = MPI_Wtime();
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

		t_end= MPI_Wtime();
		accumul = accumul + (t_end - t_begin);

		// test the notification value and update
		{
			//recvCpt++;
			nbReceived += counter;
			sender = new_notif_id - notif_offset;
			counter = new_notif_val;
			
			debug("recv", "---------------- received : " + itoa(counter) + " from : " + itoa(sender) + " level : " + itoa(level));
			//debug("infoBase", "expected : " + itoa(_Expect[iOct][sender][level]));
			//new decoding
			//debug("infoBase", "notifyID  : " + itoa(new_notif_id));
			
			updateFarFieldsFromInfos(sender, level, ff, counter, iOct);
			_Received[iOct][sender][level] += counter;			
							
			//~ /// RAZ compteur pour iteration suivantes
			//~ if (_Received[iOct][sender][level] == _Expect[iOct][sender][level])
				//~ _Received[iOct][sender][level] = 0;
		}
	}
	//add_time_sec("FF_recv_wait", accumul);
}

/* reçoit n'importe quoi, 
 * de n'importe qui, 
 * jusqu'à terminaison de la réception du niveau en cours*/
void Gaspi_FF_communicator::recv_task_ff_any_but_complete(int level, complex * ff, int iOct)
{
	double t_begin, t_end, accumul;
	accumul = 0;
	
	cout << "enter recv_task_ff_any_but_complete" << endl;
	
	// compute nb of nodes to receive at this level
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected += _Expect[iOct][i][level];
		}
	}
	
	cout << "nbRecvExpected" << nbRecvExpected << endl;
	// compute nb of already received nodes
	int nbReceived = 0;
	for (int src=0; src<_wsize; src++)
	{
		if (_Received[iOct][src][level]>0)
		{
			nbReceived += _Received[iOct][src][level];
		}
	}	 

	cout << "nbReceived" << nbReceived << endl;

	// reçoit tout et n'importe quoi :)
	// pas de sélection sur les réceptions
	gaspi_notification_id_t notif_offset = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	//int recvCpt = 0;
	int src;
	int counter;
	int recvLevel;
	
	while (nbReceived < nbRecvExpected)
	{
		t_begin = MPI_Wtime();
		while(1)
		{
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_Infos_recvbuf_seg_id,
					notif_offset,				// surveille les notifications depuis 0
					65536,						// surveille le segment complet
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

		t_end= MPI_Wtime();
		accumul = accumul + (t_end - t_begin);

		// test the notification value and update
		{
			recvLevel= new_notif_id / _wsize;
			src = new_notif_id % _wsize;
			counter = new_notif_val;
			
			//cout << "received : " << counter << " from : " << src << " level : " << recvLevel << endl;
			debug("info", "---------------- received : " + itoa(counter) + " from : " + itoa(src) + " level : " + itoa(recvLevel));
			debug("info", "expected : " + itoa(_Expect[iOct][src][recvLevel]));
										
			// new decoding
			debug("info", "notifyID  : " + itoa(new_notif_id));
							
										
			updateFarFieldChunksFromInfos(src, recvLevel, ff, counter, iOct);			
			//updateFarFieldsFromInfos   (src, recvLevel, ff, counter, iOct);
			
			_Received[iOct][src][recvLevel] += counter;			
			if (recvLevel == level)
				nbReceived += counter;
			
			//~ /// RAZ compteur pour iteration suivantes
			//~ if (_Received[iOct][src][recvLevel] == _Expect[iOct][src][recvLevel])
				//~ _Received[iOct][src][recvLevel] = 0;
		}
	}
	//add_time_sec("FF_recv_wait", accumul);
	
	cout << "exit recv_task_ff_any_but_complete" << endl;

	/// attention à l'update des compteurs et risques race condition
	/// penser au RAZ de _received
}


void Gaspi_FF_communicator::send_task_ff(int level, complex * ff, int iOct, int start, int stop)
{
	pthread_mutex_lock(&mutex);
	
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

			// if something to send to this dest ///  S IL EN RESTE ENCORE A ENVOYER
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
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];	// cell's corresponding ff in Fortran					
						
						// get indexes 
						int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];						// in FF senbuffer
						int idxInfos =_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];
	
						// write FF and Infos		
						///_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp]; 
						for (int i=0; i<__nstnsp; i++)
							_FF_sendBuffer[q+i] = ff[p+i];	

						_Infos_sendbuffer[idxInfos] = i-firstBox;
						
						
						
											
						// update counters
						_FF_sendLocalOffsets_counter[iOct][dest] += _nst[iOct][level]*_nsp[iOct][level];
						_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
						_accumul[dest] += 1;

						nbTerms = _FF_sendLocalOffsets_counter[iOct][dest]-_FF_sendLocalOffsets_keeper[iOct][dest];
						termsSize = nbTerms * sizeof(complex);
						nbInfos = _Infos_sendLocalOffsets_counter[iOct][dest]-_Infos_sendLocalOffsets_keeper[iOct][dest];

						if (_accumul[dest] == counter)
						{	
							levelIsFinished = true;
							_accumul[dest] = 0;
						}
												
						if (termsSize >= gaspiChunkSize)	// if chunk is full, send !
						{
							// send
							send_chunk(iOct, dest, level, nbTerms, nbInfos, levelIsFinished);
							
							// update keepers
							_FF_sendLocalOffsets_keeper[iOct][dest] = _FF_sendLocalOffsets_counter[iOct][dest];
							_Infos_sendLocalOffsets_keeper[iOct][dest] = _Infos_sendLocalOffsets_counter[iOct][dest];
						}
						else if (levelIsFinished) // if chunk is not full, but level is finished, send !
						{
							// send
							send_chunk(iOct, dest, level, nbTerms, nbInfos, levelIsFinished);
							
							// update keepers
							_FF_sendLocalOffsets_keeper[iOct][dest] = _FF_sendLocalOffsets_counter[iOct][dest];
							_Infos_sendLocalOffsets_keeper[iOct][dest] = _Infos_sendLocalOffsets_counter[iOct][dest];
						}
					}
				}
			}
		}
	}
	pthread_mutex_unlock(&mutex);
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
			local_offset_ff,					// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,						// remote offset
			qty,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	//double t_begin, t_end;
	//t_begin = MPI_Wtime();
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	//t_end = MPI_Wtime();
//#ifdef TIMING
	//add_time_sec("GASPI_SEND_wait_queue", t_end - t_begin);
//#endif

	gaspi_offset_t remote_offset = (_Infos_sendRemoteOffsets[octree_offset + dest] +_Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	gaspi_offset_t local_offset = (_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int);
	_Infos_sendRemoteOffsets_counter[iOct][dest] += nbInfos;
	qty = nbInfos * sizeof(int);
		
	// envoi des INFOS sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_offset,						// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_offset,						// remote offset
			qty,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
//	t_begin = MPI_Wtime();
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	//t_end = MPI_Wtime();
//~ #ifdef TIMING	
	//~ add_time_sec("GASPI_SEND_wait_queue", t_end - t_begin);
//~ #endif
	
	if (levelIsFinished)
	{
		// ENVOI DE LA NOTIFICATION
		gaspi_notification_id_t notifyID = (_rank * 65536/_wsize) + level;
		gaspi_notification_t notifyValue = level;
	
		SUCCESS_OR_DIE(
			gaspi_notify( 
				_Infos_recvbuf_seg_id,				// remote seg id
				dest,								// receiver rank
				notifyID,
				notifyValue,
				queue,								// queue
				GASPI_BLOCK							// Gaspi block
			)
		);
		
		// RAZ DES COMPTEURS
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
}



/*********************************************************************************/
void Gaspi_FF_communicator::updateFarFields(int src, int level, complex * ff, int iOct)
{	
	//double t_begin, t_end;
	
	//t_begin = MPI_Wtime();
	
	int indexToC = -1;
	//int k = level + 1;
	int octree_offset = iOct * _wsize;

	
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
	//int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
	//int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;
	
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
			//ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[q:__nstnsp];
			for (int i=0; i<__nstnsp; i++)
				ff[p0+i] = ff[p0+i] + _FF_recvBuffer[q+i];
		}
	}

	//t_end = MPI_Wtime();
//~ #ifdef TIMING
	//~ add_time_sec("updateFF", t_end - t_begin); 
//~ #endif
}

void Gaspi_FF_communicator::updateFarFields_multi(int src, int level, complex * ff, int iOct)
{	
	
	int indexToC = -1;
	//int k = level + 1;
	int octree_offset = iOct * _wsize;
	
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
	//int firstBoxToRecvIdx= _frecv[iOct][src] + indexToC;
	//int lastBoxToRecvIdx = _frecv[iOct][src+1]-1 + indexToC;
	
	bool start    = (_start_recv[iOct][src][level]>=0);
	bool stop     = ( _stop_recv[iOct][src][level]>=0);
	int  firstBox =   _stop_recv[iOct][src][level]; 
	int  lastBox  =  _start_recv[iOct][src][level];
		
	if (start && stop)
	{
		for (int j=lastBox; j>=firstBox; j--)
		{
			q = q0 + (lastBox-j)*__nstnsp;		
			int cellID = __recv[j] + indexToC; 
			int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);
//			ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[q:__nstnsp];
			for (int i=0; i<__nstnsp; i++)
				ff[p0+i] = ff[p0+i] + _FF_recvBuffer[q+i];
		}
	}
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
	//int recvDataIdx;

	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	int __endlev = _endlev[iOct][level]+indexToC;
	int __fnivnextlev = _fniv[iOct][level+1];
	
	//int nbCellsToRecv = _frecv[iOct][src+1] - _frecv[iOct][src] + 1;

	//int ff_index;
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
		
		//ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+ff_index:__nstnsp];
		for (int i=0; i<__nstnsp; i++)
			ff[p0+i] = ff[p0+i] + _FF_recvBuffer[ff_offset+ff_index+i];
		
		//debug("write", "node : " + itoa(cellID) + " from : " + itoa(ff_offset+ff_index) + " to : " + itoa(p0));		
			
	/*	if(level==4)
		{
			debug("indexes_received_from_" + convert(src), convert(index) + " at offset " + convert (i) + 
				" value : " + to_string(_FF_recvBuffer[ff_offset+ff_index].re) + " " + to_string(_FF_recvBuffer[ff_offset+ff_index].im) + 
				" ff__off : " + convert((ff_offset+ff_index)/sizeof(complex)));
		}*/
	}	
}



/*******************************
 * GASPI - Final Version
 *******************************/

/*
 * Cette fonction remplit les buffers d'envoi de FF et Infos correspondantes 
 * de manière parallèlle à l'intérieur de tâches
 * L'envoi Gaspi est effectué lorsque le buffer est prêt, et par une seule tâche
 * pour le moment prêt = level complete 
 */
void Gaspi_FF_communicator::send_task_ff_last(int level, complex * ff, int iOct, int start, int stop)
{
	int indexToC = -1;
	
	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
	// for each dest
	for (int dest=0; dest<_wsize; dest++)
	{								
		if (dest != _rank)
		{
			// if something to send to this dest			
			int nbCellsToSend =  _count_send[iOct][dest][level];		
			if (nbCellsToSend>0)
			{
				int lastBox  = _start_send[iOct][dest][level];	// the big one
				int firstBox = _stop_send[iOct][dest][level];	// the small one

				// traverses the cells to send
				for (int i= firstBox; i <= lastBox; i++)
				{					
					int cellIDf = _send[iOct][i];
					
					// if the cell is in the task's nodes range
					if (cellIDf >= start && cellIDf <= stop)
					{							
						int cellID = cellIDf + indexToC;																			// Id du node en C
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];	// @ in ff						
						int chunkIsFull = 0;
						
						// prend le VERROU
						pthread_mutex_lock(&mutex); 
						
							// lit les offsets
							int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];					// @ in FF 	  SendBuffer
							int idxInfos = _Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];		// @ in Infos SendBuffer
							
							// update les compteurs
							_FF_sendLocalOffsets_counter[iOct][dest] += _nst[iOct][level]*_nsp[iOct][level];
							_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
							
							// vérifie si le chunk est full
							if (_Infos_sendLocalOffsets_counter[iOct][dest] == _count_send[iOct][dest][level])  chunkIsFull = 1;	// dernier node à envoyer du level						


						// write FF and Infos
						//_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp]; 
						for (int i=0; i<__nstnsp; i++)
							_FF_sendBuffer[q+i] = ff[p+i]; 

						// encode level + index
						int code = level << 26 | (i-firstBox);
						//int decodeLevel = code >> 26;
//						int mask = 0x03FFFFFF;
						//int decodeIndex = code & mask;
						//int index = i-firstBox;
						
						/*
						debug("encode", 
						"index : " 			+ itoa(index)
						+ " level : " 		+ itoa(level)
						//+ " code : "		+ itoa();
						+ " decodeLevel : "	+ itoa(decodeLevel)
						+ " decodeIndex : " + itoa(decodeIndex));
						*/
						_Infos_sendbuffer[idxInfos] = code /*i-firstBox*/;	//!\ même danger, ts les threads ont-ils fini d'écrire ?
					

						// rend le VERROU
						pthread_mutex_unlock(&mutex); /// il faudrait s'assurer que tous les threads aient fini d'écrire

						if (chunkIsFull)
						{	
							_Infos_sendLocalOffsets_counter[iOct][dest]=0;
							send_chunk_encoded(iOct, dest, level);
						}
					}
				}
			}
		}
	}
}

void Gaspi_FF_communicator::send_chunk_encoded(int iOct, int dest, int level)
{
	/**** ----------- envoi des FF ----------- **/
	
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	int nbNodes =  _count_send[iOct][dest][level];
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];

	flush_queues(_nbQueues);
			
	// local FF send offset (@ + counter)
	gaspi_offset_t local_dest_offset		= _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
	gaspi_offset_t counter_offser			= _FF_sendRemoteOffsets_counter[iOct][dest] * sizeof(complex);
	gaspi_offset_t local_offset_ff			= local_dest_offset + counter_offser;
	gaspi_queue_id_t queue=0;
	
	// remote FF offset (@ + counter)
	gaspi_offset_t remote_sender_offset 	= _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
	gaspi_offset_t remote_offset_ff 		= remote_sender_offset + counter_offser;


//ATTENTION RISQUE DE RACE CONDITION ICI :
	// update remote offset counter
	_FF_sendRemoteOffsets_counter[iOct][dest] += nbNodes * __nstnsp;
	
	// compute size
	gaspi_size_t qty= nbNodes * __nstnsp * sizeof(complex);
	
	// envoi des DATA sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_FF_sendBuf_seg_id,					// local seg ID
			local_offset_ff,					// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,					// remote offset
			qty,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	/**** ----------- envoi des Infos ----------- **/

	//!\ ne fonctionnera pas si arbres de hauteurs différentes
	gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
	
	/// prêt à changer le notif ID ici
	/*int myMin = _rank * 65536/_wsize;
	_notifCpt = (_notifCpt + 1) % (65536/_wsize);
	gaspi_notification_id_t notifyID2 = myMin + _notifCpt;*/
	
	// Infos local offset (@ + counter) + update
	// Infos remote offset (@ + counter) + update
	gaspi_offset_t local_offset 	= (_Infos_sendLocalOffsets[iOct][dest] /* +_Infos_sendRemoteOffsets_counter[iOct][dest]*/)* sizeof(int);
	gaspi_offset_t remote_offset 	= (_Infos_sendRemoteOffsets[octree_offset + dest] +_Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	_Infos_sendRemoteOffsets_counter[iOct][dest] += nbNodes;

	int qtyInfos = nbNodes * sizeof(int);

	/*gaspi_offset_t remote_offset_infos = (_Infos_sendRemoteOffsets[octree_offset + dest] + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	gaspi_offset_t local_offset_infos =    (_Infos_sendLocalOffsets[iOct][dest]            + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int);
	_Infos_sendRemoteOffsets_counter[iOct][dest] += nbInfos;//
	
	gaspi_size_t qtyFF = nbTerms * sizeof(complex);
	gaspi_size_t qtyInfos = nbInfos * sizeof(int);*/

	// envoi des INFO et NOTIFICATION
	SUCCESS_OR_DIE(
		gaspi_write_notify( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_offset,						// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_offset,						// remote offset
			qtyInfos,								// size of data to write
			notifyID,							// remote notif ID
			nbNodes,							// value of the notif to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
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
}

void Gaspi_FF_communicator::recv_task_ff_last_old(int level, complex * ff, int iOct)
{
	double t_begin, t_end, accumul;
	accumul = 0;
	
	// compute nb of messages to receive at this level
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected++;
		}
	}

	// 
	gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int recvCpt = 0;
	int sender;
	int counter;
	
	while (recvCpt < nbRecvExpected)
	{
		t_begin = MPI_Wtime();
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

		t_end= MPI_Wtime();
		accumul = accumul + (t_end - t_begin);

		// test the notification value and update
		{
			recvCpt++;
			sender = new_notif_id - notif_offset;
			counter = new_notif_val;
			
			//debug("infoBase", "---------------- received : " + itoa(counter) + " from : " + itoa(sender) + " level : " + itoa(level));
			//debug("infoBase", "expected : " + itoa(_Expect[iOct][sender][level]));										
			// new decoding
			//debug("infoBase", "notifyID  : " + itoa(new_notif_id));
			
			updateFarFieldsFromInfosEncoded(sender, level, ff, counter, iOct);
		}
	}
	//~ add_time_sec("FF_recv_wait", accumul);
}




/*******************************
 * 		GASPI
 *******************************/


void Gaspi_FF_communicator::send_ff_level(int level, complex * ff, int iOct)
{
	//double t_begin, t_end, t_begin_comm, t_end_comm, accumul_comm;
	
	//t_begin = MPI_Wtime();
	//accumul_comm = 0;
	
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	
	// flush queue
	flush_queues(_nbQueues);
	
	// prepare _SendBuffer
	for (int dest=0; dest<_wsize; dest++)
	{
		if (dest != _rank)
		{

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
					//_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp];
					for (int i=0; i<__nstnsp; i++)		
						_FF_sendBuffer[q+i] = ff[p+i];
				}
			}
			
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
				
				//t_begin_comm = MPI_Wtime();
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
				//t_end_comm = MPI_Wtime();
				//accumul_comm = accumul_comm + (t_end_comm - t_begin_comm);
			}
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
	//t_end = MPI_Wtime();
//~ #ifdef TIMING	
	//~ add_time_sec("FF_sendrecv_send_comm", accumul_comm);
	//~ add_time_sec("FF_sendrecv_send_copy", t_end - t_begin - accumul_comm);
	//~ add_time_sec("FF_sendrecv_send_total", t_end - t_begin);
//~ #endif
}

void Gaspi_FF_communicator::recv_ff_level(int level, complex * ff, int iOct)
{
	//int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
	//double t_begin, t_end, t_begin_comm, t_end_comm, accumul_comm;
	
//	t_begin = MPI_Wtime();
	//accumul_comm = 0;
	
	// wait to receive all infos
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected++;
		}
	}
	
	// create timestamps array
	//vector<uint64_t> timestamps(nbRecvExpected,0);

	/*int rankMultiple = level;
	gaspi_notification_id_t notif_offset = _wsize * rankMultiple;*/
	gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	int recvCpt = 0;
	int sender;
	
	
	while (recvCpt < nbRecvExpected)
	{
		//t_begin_comm = MPI_Wtime();
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
		
		//t_end_comm = MPI_Wtime();
		//accumul_comm = accumul_comm + (t_end_comm - t_begin_comm);		

		// update timestamps array
		//timestamps[recvCpt] = rdtsc();
		
		if (new_notif_val)
		{
			recvCpt++;
			sender = new_notif_id - notif_offset;
			updateFarFields(sender, level, ff, iOct);
		}
	}
	
	//t_end = MPI_Wtime();
//~ #ifdef TIMING	
	//~ add_time_sec("FF_recv_wait", accumul_comm);
	//~ add_time_sec("FF_sendrecv_recv_comm", accumul_comm);
	//~ add_time_sec("FF_sendrecv_recv_copy", t_end - t_begin - accumul_comm);
	//~ add_time_sec("FF_sendrecv_recv_total", t_end - t_begin);
//~ #endif
	// dump timestamps array
	/*int wsize; MPI_Comm_size(MPI_COMM_WORLD,&wsize);
	string file = "Gaspi_async/timestamps/" + to_string((unsigned long long)wsize) +"/recv_timestamps";
	dumpBuffer(rank, timestamps.data(), nbRecvExpected, file,"level " + to_string((unsigned long long)level));*/
}


void Gaspi_FF_communicator::send_ff_level_multithreaded(int level, complex * ff, int iOct)
{
	// timers
	//double t_begin, t_end;
	//double mean_time = 0.0;
	
	// init threadsafe timers
	//int num_threads;
/*
#pragma omp parallel
{
	num_threads = omp_get_num_threads();
}
	double * t_begin_comm = new double[num_threads]();
	double * t_end_comm = new double[num_threads]();
	double * accumul_comm = new double[num_threads]();
*/
		
	//t_begin = MPI_Wtime();
	
	// begin comm prep
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	flush_queues(_nbQueues);

#pragma omp parallel
{
	//int tid = omp_get_thread_num();
	
	// prepare _SendBuffer
	#pragma omp single
	for (int dest=0; dest<_wsize; dest++)
	{
		#pragma omp task
		{
			if (dest != _rank)
			{

				int q = _FF_sendLocalOffsets[iOct][dest] - 1 + _FF_sendLocalOffsets_counter[iOct][dest];
				int q0 = q + 1; 

				int count =  _count_send[iOct][dest][level];
				int lastBox =  _start_send[iOct][dest][level];
				int firstBox = _stop_send[iOct][dest][level];
				int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
				
				if (count>0)
				{
					// copy data into sendBuffer
					for (int k=lastBox; k>=firstBox; k--)
					{
						int cellID = _send[iOct][k] + indexToC;
						int p=_fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*_nst[iOct][level]*_nsp[iOct][level];
						q = q0 + (lastBox-k)*__nstnsp;
						
						//_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp];						
						for (int i=0; i<__nstnsp; i++)		
							_FF_sendBuffer[q+i] = ff[p+i];
					}

					// local offset
					gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
					gaspi_offset_t level_offset = _FF_sendLocalOffsets_counter[iOct][dest] * sizeof(complex);
					gaspi_offset_t local_offset = local_dest_offset + level_offset;
					gaspi_queue_id_t queue=0;
					
					// remote offset
					gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
					gaspi_offset_t remote_offset = remote_sender_offset + level_offset;

					// update offset, per destinatary
					_FF_sendLocalOffsets_counter[iOct][dest] += count * _nst[iOct][level] * _nsp[iOct][level];

					gaspi_notification_id_t notifyID = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize) + _rank;
					gaspi_size_t qty= count * _nst[iOct][level] * _nsp[iOct][level] * sizeof(complex);
					
					//t_begin_comm[tid] = MPI_Wtime();
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
					/*t_end_comm[tid] = MPI_Wtime();
					accumul_comm[tid] = accumul_comm[tid] + (t_end_comm[tid] - t_begin_comm[tid]);*/
				}
			}

			/* tous les tableaux sont privés, car propres au dest*/
			// cas include levcom
			if(_incLevcom)
			{
				if (level == _levcom[iOct] + indexToC)
				{
					 _FF_sendLocalOffsets_counter[iOct][dest] = 0; 
					// RAZ des compteurs
					_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
					_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
					_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
				}
			}
			else // cas allreduce sur levcom
			{
				if (level == _levcom[iOct] + 1 + indexToC)
					 _FF_sendLocalOffsets_counter[iOct][dest] = 0;
			}
		} // task	
	} // omp single
} // omp parallel

	
	/*delete [] accumul_comm;
	delete [] t_begin_comm;
	delete [] t_end_comm;*/
	
	//t_end = MPI_Wtime();
	
	/*
	for (int i=0; i<num_threads; i++)
		mean_time += accumul_comm[i];
	mean_time = mean_time * 1.0 / num_threads;
	*/
	
	//add_time_sec("FF_sendrecv_comm", mean_time);
//~ #ifdef TIMING
	//~ add_time_sec("FF_sendrecv_send_total", t_end - t_begin);
//~ #endif
}

void Gaspi_FF_communicator::recv_ff_level_multithreaded(int level, complex * ff, int iOct)
{
	// timers
	//double t_begin, t_end;
	/*double mean_time = 0.0;
	double mean_time_copy = 0.0;*/
	
	//t_begin = MPI_Wtime();	
	

	// compute nb infos to receive
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected++;
		}
	}


	
#pragma omp parallel
{
	//int tid = omp_get_thread_num();
	
	#pragma omp single
	for (int i=0; i<nbRecvExpected; i++)
	{
		#pragma omp task
		{	
			gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
			gaspi_notification_id_t new_notif_id;
			gaspi_notification_t new_notif_val;
			int sender;

//			t_begin_comm[tid] = MPI_Wtime();			
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
//			t_end_comm[tid] = MPI_Wtime();
//			accumul_comm[tid] = accumul_comm[tid] + (t_end_comm[tid] - t_begin_comm[tid]);
			
			
//			t_begin_copy[tid] = MPI_Wtime();
			if (new_notif_val)
			{
				//recvCpt++;
				sender = new_notif_id - notif_offset;
				
				
				updateFarFields_multi(sender, level, ff, iOct);
			}
//			t_end_copy[tid] = MPI_Wtime();
//			accumul_copy[tid] = accumul_copy[tid] + (t_end_copy[tid] - t_begin_copy[tid]);
		} // task
	} // single
} // parallel

	//t_end = MPI_Wtime();
	
//~ #ifdef TIMING
	//~ add_time_sec("FF_sendrecv_recv_total", t_end - t_begin);
//~ #endif
}

void Gaspi_FF_communicator::recv_ff_level_multithreaded_2(int level, complex * ff, int iOct)
{
	//double t_begin, t_end, t_begin_comm, t_end_comm, t_begin_copy, t_end_copy;
	//t_begin = MPI_Wtime();
	
	// compute nb infos to receive
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected++;
		}
	}

	//t_begin_comm = MPI_Wtime();
#pragma omp parallel for
//{	
	//#pragma omp single
	for (int i=0; i<nbRecvExpected; i++)
	{
		#pragma omp task
		{	
			gaspi_notification_id_t notif_offset = (iOct * _nivterm[iOct] * _wsize) + (level * _wsize);
			gaspi_notification_id_t new_notif_id;
			gaspi_notification_t new_notif_val;
			//int sender;

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
		} // task
	} // single
//} // parallel
	//t_end_comm = MPI_Wtime();

	//t_begin_copy = MPI_Wtime();
	/* update FF en dehors des tasks */
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			updateFarFields(i, level, ff, iOct);
		}
	}
	//t_end_copy = MPI_Wtime();

	//t_end = MPI_Wtime();
//~ #ifdef TIMING	
	//~ add_time_sec("FF_sendrecv_recv_total", t_end - t_begin);
	//~ add_time_sec("FF_sendrecv_recv_comm", t_end_comm - t_begin_comm);
	//~ add_time_sec("FF_sendrecv_recv_copy", t_end_copy - t_begin_copy);	
//~ #endif
}

//#include "Gaspi_FF_communicator.cpp_suite"----------------------------


/*
 * Version Gaspi Tasks Chunks 
 */
void Gaspi_FF_communicator::send_task_ff_dbg2(int level, complex * ff, int iOct, int start, int stop)
{
	int indexToC = -1;
	int termsSize = 0;
	int nbTerms = 0;
	int nbInfos = 0;
	bool levelIsFinished = false;
	bool readyToSend = false;
	//int tid = omp_get_thread_num();

	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
	// for each dest
	for (int dest=0; dest<_wsize; dest++)
	{
		levelIsFinished = false;

		if (dest != _rank)
		{
			// nb of cells to be sent
			int counter =  _count_send[iOct][dest][level];

			// if something to send to this dest 
			if (counter>0)
			{
				int lastBox =  _start_send[iOct][dest][level];	// the big one
				int firstBox = _stop_send[iOct][dest][level];	// the small one

				// traverses the cells to send
				for (int i= firstBox; i <= lastBox; i++)
				{

					int cellIDf = _send[iOct][i];
					readyToSend = false;
					
					// if the cell is in the task's nodes range
					if (cellIDf >= start && cellIDf <= stop)
					{
						int cellID = cellIDf + indexToC;
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*__nstnsp;	// cell's corresponding ff index in Fortran						

/* --- LOCK MUTEX --- */  
pthread_mutex_lock(&mutex);
						// read and update adresses 
						int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];						// in FF sendbuffer
						int idxInfos =_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];			// in Infos sendBuffer
						_FF_sendLocalOffsets_counter[iOct][dest] += __nstnsp;
						_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
						_accumul[dest] += 1;
						
						// update send threashold infos
						nbTerms = _FF_sendLocalOffsets_counter[iOct][dest]-_FF_sendLocalOffsets_keeper[iOct][dest]; 
						termsSize = nbTerms * sizeof(complex);
						nbInfos = _Infos_sendLocalOffsets_counter[iOct][dest]-_Infos_sendLocalOffsets_keeper[iOct][dest]; 

						// Tasks version counters
						/*int cptIndex = _write_counters_ptr[iOct][dest][level];
						_write_counters[iOct][dest][level][cptIndex]--; */

						// test if ready to Send
						if (_accumul[dest] == counter) // tout le level est écrit
						{	
							levelIsFinished = true;
							readyToSend = true;
							_accumul[dest] = 0;
							//_write_counters_ptr[iOct][dest][level] = 0; // raz for next iteration
						}
						else if (termsSize >= _gaspiChunkSize) // comm_size atteint
						{
							readyToSend = true;
							//_write_counters_ptr[iOct][dest][level]++; // update for next packet
						}
						
						// if send --> update keepers and counters pointer
						if (readyToSend)
						{
							_FF_sendLocalOffsets_keeper[iOct][dest] = _FF_sendLocalOffsets_counter[iOct][dest];
							_Infos_sendLocalOffsets_keeper[iOct][dest] = _Infos_sendLocalOffsets_counter[iOct][dest];
						}
//pthread_mutex_unlock(&mutex);
						
						// write FF and Infos
						//_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp];
						for (int i=0; i<__nstnsp; i++)		
							_FF_sendBuffer[q+i] = ff[p+i];
						
						_Infos_sendbuffer[idxInfos] = i-firstBox;
						
//pthread_mutex_lock(&mutex);
						// update le compteur écriture en cours
						//printf("tid : %i, cpt :%i \n", tid, cptIndex);
						//debug("debug","tid: " + itoa(tid) + " cptIndex: " + itoa(cptIndex) + " cell: " + itoa (i) + " dest: " + itoa(dest));
						
						//_write_counters[iOct][dest][level][cptIndex]++;
						if (readyToSend)
						{
							/*while ( _write_counters[iOct][dest][level][cptIndex] != 0 )
							{ 
								sleep(1); printf("sleep");
							}*/
							send_chunk_dbg(iOct, dest, level, nbTerms, nbInfos, levelIsFinished);
						}
/* --- LOCK MUTEX --- */ 
pthread_mutex_unlock(&mutex);
					} // if start < cell < stop
				} // for cells
			} // if counter >0
		} // if dest != rank
	} // for dest
}

void Gaspi_FF_communicator::send_chunk_dbg(int iOct, int dest, int level, int nbTerms, int nbInfos, bool levelIsFinished)
{
//printf("call send chunk dest %i\t level %i\t nbTerms %i\t nbInfos%i\t levelIsFinished %i\n",dest, level, nbTerms, nbInfos, levelIsFinished); 
/* --- LOCK MUTEX --- */ //pthread_mutex_lock(&mutex);

	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	flush_queues(_nbQueues);
			
	// local FF send offset
	gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
	gaspi_offset_t accumul_offset = _FF_sendRemoteOffsets_counter[iOct][dest] * sizeof(complex);
	gaspi_offset_t local_offset_ff = local_dest_offset + accumul_offset;
	gaspi_queue_id_t queue=0;
	
	// remote offsets FF et Infos
	gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
	gaspi_offset_t remote_offset_ff = remote_sender_offset + accumul_offset;
	_FF_sendRemoteOffsets_counter[iOct][dest] += nbTerms;//
	
	gaspi_offset_t remote_offset_infos = (_Infos_sendRemoteOffsets[octree_offset + dest] +_Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	gaspi_offset_t local_offset_infos = (_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int);
	_Infos_sendRemoteOffsets_counter[iOct][dest] += nbInfos;//
	
	gaspi_size_t qtyFF = nbTerms * sizeof(complex);
	gaspi_size_t qtyInfos = nbInfos * sizeof(int);

	
	// envoi des DATA sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_FF_sendBuf_seg_id,					// local seg ID
			local_offset_ff,					// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,						// remote offset
			qtyFF,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
		
	// envoi des INFOS sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_offset_infos,					// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_offset_infos,				// remote offset
			qtyInfos,							// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
	if (levelIsFinished)
	{
		// ENVOI DE LA NOTIFICATION
		gaspi_notification_id_t notifyID = (_rank * 65536/_wsize) + level;
		gaspi_notification_t notifyValue = level;
	
		SUCCESS_OR_DIE(
			gaspi_notify( 
				_Infos_recvbuf_seg_id,				// remote seg id
				dest,								// receiver rank
				notifyID,
				notifyValue,
				queue,								// queue
				GASPI_BLOCK							// Gaspi block
			)
		);
		
		// RAZ DES COMPTEURS lorsqu'on a fini d'envoyer
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
//pthread_mutex_unlock(&mutex);
}

/* attend 1 message par source, pour 1 level donné */
void Gaspi_FF_communicator::recv_task_ff(int level, complex * ff, int iOct)
{
	//cout << "ENTER recv task ff\n";
	//double t_begin, t_end, t_begin_r, t_end_r, accumul_r;
	//accumul_r = 0;
	//t_begin = MPI_Wtime();
	
	// compute nb of nodes to receive at current level
	int nbSrcs = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbSrcs++; // ajoute 1 source
		}		
	}

	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	//int src;
	//int nbNodes;
	int nbMsg = 0;	

	while (nbMsg != nbSrcs)
	{
		//printf ("while recv \n");
		for (int i=0; i<_wsize; i++)
		{
			if(gaspi_notify_waitsome( // test level L for src I
				_Infos_recvbuf_seg_id,
				(65536/_wsize * i) + level,
				1,
				&new_notif_id,
				GASPI_TEST) == GASPI_SUCCESS)
			{
				gaspi_notify_reset(_Infos_recvbuf_seg_id, new_notif_id, &new_notif_val);
				//t_begin_r = MPI_Wtime();
				updateFarFieldChunksFromInfos(i, level, ff, _Expect[iOct][i][level], iOct);
				//t_end_r = MPI_Wtime();
				//accumul_r = accumul_r + t_end_r - t_begin_r;
				nbMsg++;
			}
		}		
	}
	
	// raz for next time
	for (int i=0; i<_wsize; i++)
	{
		_Received[iOct][i][level] = 0;
	}
	
	//t_end = MPI_Wtime();
//~ #ifdef TIMING
	//~ add_time_sec("FF_read_from_buffer", accumul_r);
	//~ add_time_sec("FF_buffering", accumul_r);
	//~ add_time_sec("GASPI_FF_sendrecv", t_end - t_begin - accumul_r);
//~ #endif
	//cout << "EXIT recv task ff\n";
}

void Gaspi_FF_communicator::updateFarFieldChunksFromInfos(int src, int level, complex * ff, int nbNodes, int iOct)
{
	//debug("info", "in --> update FF, received : " + itoa(nbNodes) + " from : " + itoa(src) + " level : " + itoa(level));
	
	int indexToC = -1;
	int octree_offset = iOct * _wsize;


	// calcul de l'offset du niveau
	int info_level_offset = 0;
	int ff_level_offset = 0;
	for (int i=_nivterm[iOct]-1; i>level; i--)
	{
		info_level_offset += _Expect[iOct][src][i];
		ff_level_offset += _Expect[iOct][src][i]*_nst[iOct][i]*_nsp[iOct][i];

	}
	
	int info_offset = _Infos_recvOffsets[octree_offset + src]/*src*/ 
					+ info_level_offset/*lvl*/ 
					+ _Received[iOct][src][level] /*déjà reçus*/; 
	
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];	
	int ff_offset = _FF_recvOffsets[octree_offset + src] 
					+ ff_level_offset 
					+ _Received[iOct][src][level]*__nstnsp;

	/*debug("info", "src offst : " + itoa(_Infos_recvOffsets[octree_offset + src]) 
								 + " lvl offset : " + itoa(info_level_offset) 
								 + " already received offset : " + itoa(_Received[iOct][src][level]));*/

	int firstBox = _stop_recv[iOct][src][level];
	
	int index;


	int __endlev = _endlev[iOct][level]+indexToC;
	int __fnivnextlev = _fniv[iOct][level+1];
	//int ff_index;
	
	for (int i=info_offset; i<info_offset+nbNodes; i++)
	{
		index = _Infos_recvbuffer[i];
		int cellIDf = _recv[iOct][firstBox+index];
		int cellID = cellIDf + indexToC;
		//printf ("received cellID : %i\n", cellID); fflush(stdout);
		//debug("received_cellIDs", "cellID : " + itoa(cellID) + " i/nbNodes : " + itoa(i) + "/" + itoa(nbNodes));


		// ranger la boîte
		int ff_index = (i-info_offset) * __nstnsp;
		int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);		
		//ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+ff_index:__nstnsp];
		for (int i=0; i<__nstnsp; i++)		
			ff[p0+i] = ff[p0+i] + _FF_recvBuffer[ff_offset+ff_index+i];

		
		//debug("write", "node : " + itoa(cellID) + " from : " + itoa(ff_offset+ff_index) + " to : " + itoa(p0));
	}	
}

void Gaspi_FF_communicator::updateFarFieldsFromInfosEncoded(int src, int level, complex * ff, int counter, int iOct)
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
	//int recvDataIdx;

	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	int __endlev = _endlev[iOct][level]+indexToC;
	int __fnivnextlev = _fniv[iOct][level+1];
	
	//int nbCellsToRecv = _frecv[iOct][src+1] - _frecv[iOct][src] + 1;

	//int ff_index;
	for (int i=info_offset; i<info_offset+counter; i++)
	{
		//index = _Infos_recvbuffer[i];

		///int code = level << 26 | (i-firstBox);
		///int decodeLevel = code >> 26;
		///int mask = 0x03FFFFFF;
		///int decodeIndex = code & mask;
		///int index = i-firstBox;
		
		int code = _Infos_recvbuffer[i];
		int mask = 0x03FFFFFF;
		index = code & mask;
		
		int cellIDf = _recv[iOct][firstBox+index];
		int cellID = cellIDf + indexToC;

		// ranger la boîte
		int ff_index = (i-info_offset) * __nstnsp;
		int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);
		
		//ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+ff_index:__nstnsp];		
		for (int i=0; i<__nstnsp; i++)
			ff[p0+i] = ff[p0+i] + _FF_recvBuffer[ff_offset+ff_index+i];
				
	}
}


/*
 * Version Gaspi Tasks Chunks --- DEV --- RECV ANY
 */
 

/// NEW VERSION ******************************
// --------------------------------- SEND ----------------------------------------------------------------------------------- 

void Gaspi_FF_communicator::send_task_ff_chunk(int level, complex * ff, int iOct, int start, int stop)
{
	int indexToC = -1;
	int termsSize = 0;
	int nbInfos = 0;
	bool destIsFinished = false;
	bool readyToSend = false;
	
	// simplifier l'ecriture
	/**int * local_FF_offset =   ;
	int * local_Infos_offset = ;
	int * local_FF_cpt =		 ;
	int * local_Infos_cpt = 	 ;

	int * sendFFcpt = ;
	int * sendInfoscpt = ;*/
	int nbTerms  = _nst[iOct][level]*_nsp[iOct][level];
	int nbTermsToSend = 0;
	
		
	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
	// for each dest
	for (int dest=0; dest<_wsize; dest++)
	{
		destIsFinished = false;

		if (dest != _rank)
		{
			// nb of cells to be sent
			int counter =  _count_send[iOct][dest][level];
			int totalCellsDest = _CountDest[iOct][dest];

			// if something to send to this dest 
			if (counter>0)
			{
				int lastBox =  _start_send[iOct][dest][level];	// the big one
				int firstBox = _stop_send[iOct][dest][level];	// the small one

				// traverses the cells to send
				for (int i= firstBox; i <= lastBox; i++)
				{

					int cellIDf = _send[iOct][i];
					readyToSend = false;
					
					// if the cell is in the task's nodes range
					if (cellIDf >= start && cellIDf <= stop)
					{
						int cellID = cellIDf + indexToC;
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*__nstnsp;	// cell's corresponding ff index in Fortran

/* --- LOCK MUTEX --- */  
pthread_mutex_lock(&mutex);
						// read address
						int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];				// @ in FF sendbuffer
						int idxInfos = _Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];	// @ in Infos sendBuffer

						// update cpt
						_FF_sendLocalOffsets_counter[iOct][dest] += nbTerms;
						_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
						
						// update next comm send cpt
						_accumul[dest] += 1;
					
						// update send threashold infos
						nbTermsToSend = _FF_sendLocalOffsets_counter[iOct][dest]-_FF_sendLocalOffsets_keeper[iOct][dest]; 
						termsSize = nbTermsToSend * sizeof(complex);
						nbInfos = _Infos_sendLocalOffsets_counter[iOct][dest]-_Infos_sendLocalOffsets_keeper[iOct][dest]; 

						/// prend un compteur et update (va écrire dans packet)
						int cptIndex = _write_counters_ptr[iOct][dest];
						_write_counters[iOct][dest][cptIndex]--;

						// test if ready to Send
						if (_accumul[dest] == totalCellsDest) // tout le level est écrit /// changer par : c'est la tte dernière boite
						{	
							destIsFinished = true;
							readyToSend = true;
							_accumul[dest] = 0;
							_write_counters_ptr[iOct][dest] = 0; // raz for next iteration
						}
						else if (termsSize >= _gaspiChunkSize) // comm_size atteint
						{
							readyToSend = true;
							_write_counters_ptr[iOct][dest]++; // update for next packet
						}
						
						// if send --> update keepers and counters pointer
						if (readyToSend)
						{
							_FF_sendLocalOffsets_keeper[iOct][dest] = _FF_sendLocalOffsets_counter[iOct][dest];
							_Infos_sendLocalOffsets_keeper[iOct][dest] = _Infos_sendLocalOffsets_counter[iOct][dest];
						}

						// write FF and Infos
						//_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp];
						for (int i=0; i<__nstnsp; i++)
							_FF_sendBuffer[q+i] = ff[p+i];
						
						
						// encode level + index
						int code = level << 26 | (i-firstBox);
						//int decodeLevel = code >> 26;
//						int mask = 0x03FFFFFF;
						//int decodeIndex = code & mask;
						//int index = i-firstBox;
						_Infos_sendbuffer[idxInfos] = code;
						_write_counters[iOct][dest][cptIndex]++;
						
						if (readyToSend)
						{
							while ( _write_counters[iOct][dest][cptIndex] != 0 ) // tt le monde n'a pas fini d'écrire
							{ 
								sleep(1); 
								printf("sleep");
							}
							send_chunk_dbg_last_notif_each(iOct, dest, level, nbTermsToSend, nbInfos, destIsFinished);
						}
						else
						{
pthread_mutex_unlock(&mutex);
						}
/* --- LOCK MUTEX --- */ 


					} // if start < cell < stop
				} // for cells
			} // if counter >0
		} // if dest != rank
	} // for dest
}

void Gaspi_FF_communicator::send_chunk_dbg_last_notif_each(int iOct, int dest, int level, int nbTerms, int nbInfos, bool destIsFinished)
{
//pthread_mutex_lock(&mutex1);
	//	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	flush_queues(_nbQueues);
			
	// local FF send offset
	gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);		/// non modifié
	gaspi_offset_t accumul_offset = _FF_sendRemoteOffsets_counter[iOct][dest] * sizeof(complex);// modif ici (cette fx) slmt
	gaspi_offset_t local_offset_ff = local_dest_offset + accumul_offset;
	gaspi_queue_id_t queue=0;
	
	// remote offsets FF et Infos
	gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
	gaspi_offset_t remote_offset_ff = remote_sender_offset + accumul_offset;
	_FF_sendRemoteOffsets_counter[iOct][dest] += nbTerms;// modif ici (cette fx) slmt
	
	gaspi_offset_t remote_offset_infos = (_Infos_sendRemoteOffsets[octree_offset + dest] +_Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	gaspi_offset_t local_offset_infos = (_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int);
	_Infos_sendRemoteOffsets_counter[iOct][dest] += nbInfos;// modif ici (cette fx) slmt
	
	gaspi_size_t qtyFF = nbTerms * sizeof(complex);
	gaspi_size_t qtyInfos = nbInfos * sizeof(int);
	
	_notifCpt = (_notifCpt + 1) % (65536/_wsize);		// modif ici slmt
	int myMin = _rank * 65536/_wsize;	
	gaspi_notification_id_t notifyID2 = myMin + _notifCpt;	



	if (destIsFinished) /// ne devrait pas y avoir de race condition possible (2 iters ne se mélangent pas)
	{
		_FF_sendLocalOffsets_counter[iOct][dest] = 0;
		_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
		_Infos_sendLocalOffsets_counter[iOct][dest] = 0;
		_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
		_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
		_Infos_sendLocalOffsets_keeper[iOct][dest] = 0;
	}

// RELEASE MUTEX HERE
pthread_mutex_unlock(&mutex);

	
	
	
	// envoi des DATA sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_FF_sendBuf_seg_id,					// local seg ID
			local_offset_ff,					// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,						// remote offset
			qtyFF,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
	// envoi des INFOS sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_offset_infos,					// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_offset_infos,				// remote offset
			qtyInfos,							// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
	// envoi de la NOTIFICATION
	gaspi_notification_t notifyValue2 = nbInfos;

	SUCCESS_OR_DIE(
		gaspi_notify( 
			_Infos_recvbuf_seg_id,				// remote seg id
			dest,								// receiver rank
			notifyID2,
			notifyValue2,
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
}

// --------------------------------- RECV ----------------------------------------------------------------------------------- 

/* reçoit n'importe quoi, 
 * de n'importe qui, 
 * jusqu'à terminaison de la réception du niveau en cours*/
void Gaspi_FF_communicator::recv_task_ff_any_but_complete_level(int level, complex * ff, int iOct)
{
	//double t_begin, t_end, accumul;
	//accumul = 0;
	
	//cout << "enter recv_task_ff_any_but_complete, wait on level : " << level << endl;
	
	// compute nb of nodes to receive, from all sources, at this level
	int nbRecvExpected = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbRecvExpected += _Expect[iOct][i][level];
		}
	}
	
	//cout << "nbRecvExpected" << nbRecvExpected << endl;
	
	// compute nb of already received nodes
	int nbReceived = 0;
	for (int src=0; src<_wsize; src++)
	{
		if (_Received[iOct][src][level]>0)
		{
			nbReceived += _Received[iOct][src][level];
		}
	}	 

	//cout << "nbReceived" << nbReceived << endl;

	// reçoit tout et n'importe quoi :)
	// pas de sélection sur les réceptions
	//gaspi_notification_id_t notif_offset = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	
	while (nbReceived < nbRecvExpected)
	{
		//t_begin = MPI_Wtime();
		while(1)
		{
			SUCCESS_OR_DIE(
				gaspi_notify_waitsome(
					_Infos_recvbuf_seg_id,
					0,				// surveille les notifications depuis 0
					65536,			// surveille le segment complet
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

		//t_end= MPI_Wtime();
		//accumul = accumul + (t_end - t_begin);

		// test the notification value and update
		{
			int src = new_notif_id/(65536/_wsize);
			int qty = new_notif_val;
			
			//debug("info", "---------------- received : " + itoa(qty) + " from : " + itoa(src));													
			updateFarFieldChunksFromInfos_last_withoutLevel(src, ff, qty, iOct, nbReceived, level);
		}
	}
//	add_time_sec("FF_recv_wait", accumul);
	
	//cout << "exit recv_task_ff_any_but_complete" << endl;

	/// attention à l'update des compteurs et risques race condition
	/// penser au RAZ de _received
}

void Gaspi_FF_communicator::updateFarFieldChunksFromInfos_last_withoutLevel(int src, complex * ff, int nbNodes, int iOct, int &cpt, int expectedLevel)
{
	//debug("info", "in --> update FF, received : " + itoa(nbNodes) + " nodes, from : " + itoa(src));
	
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	
	// lire les infos
	int * recvInfosCpt   = _ReceivedInfosCpt[iOct];	//_ReceivedInfosBytes[iOct];
	int * recvFFTermsCpt = _ReceivedFFTermsCpt[iOct];	//_ReceivedFFOffsets[iOct];
	//debug("info", "nb Infos already received from : " + itoa(src) + " = " + itoa(recvInfosCpt[src]));
	//debug("info", "nb Terms already received from : " + itoa(src) + " = " + itoa(recvFFTermsCpt[src]));

	// nouveau calcul des offsets du niveau 
	int info_offset = _Infos_recvOffsets[octree_offset + src]/*src*/ 
					 + recvInfosCpt[src]; /*déjà reçus*/;
					
	int ff_offset	 = _FF_recvOffsets[octree_offset + src]; 
					/// + recvFFTermsCpt[src];

	// update des compteurs infos
	recvInfosCpt[src] += nbNodes;

	// lecture 
	for (int i=info_offset; i<(info_offset + nbNodes) ; i++)
	{
		// decodage, à reconstituer pour chaque node
		int code = _Infos_recvbuffer[i];
		int mask = 0x03FFFFFF;
		
		int index = code & mask;
		int level = code >> 26;
		
		//debug("info", "index : " + itoa(index) + " level " + itoa(level));
		
		// valeurs qui dépendent de level
		int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
		int firstBox = _stop_recv[iOct][src][level];
		int __endlev = _endlev[iOct][level]+indexToC;
		int __fnivnextlev = _fniv[iOct][level+1];
		
		// identification du node
		int cellIDf = _recv[iOct][firstBox+index];
		int cellID = cellIDf + indexToC;

		// ranger la boîte
		///int ff_index = (i-info_offset) * __nstnsp; /// en prévision : cas mélange des niveaux
		int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);		
		//ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+ff_index:__nstnsp];
		
		
		/// attention : pour le mélange des niveaux
		//ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+recvFFTermsCpt[src]/*ff_index*/:__nstnsp];
		for (int i=0; i<__nstnsp; i++)
			ff[p0+i] = ff[p0+i] + _FF_recvBuffer[ff_offset+recvFFTermsCpt[src]+i];




		//debug("sort", "cellID : " + itoa(cellID) + " level " + itoa(level) + ", ff@ : " + itoa(p0));

		
		// update du compteur de termes (pour chaque node car level peut varier)
		recvFFTermsCpt[src] += __nstnsp;
		
		// update pour la boucle de reception Gaspi
		_Received[iOct][src][level]++;
		if (level == expectedLevel)
		{
			cpt++;
			/// RAZ compteur pour iteration suivantes
			if (_Received[iOct][src][expectedLevel] == _Expect[iOct][src][expectedLevel])
				_Received[iOct][src][expectedLevel] = 0;
		}

	}

	// RAZ si fini
	// calcul nb total de nodes à recevoir pour src // => à sortir pour ne pas refaire à chaque fois
	int totalNodesExpected = _ExpectSrc[iOct][src];
	/*for (int lvl=0; lvl<_nivterm[iOct]; lvl++)
		totalNodesExpected += _Expect[iOct][src][lvl];*/

	//debug("info", "total nodes expected from : " + itoa(src) + " = " + itoa(totalNodesExpected));
	//debug("info", "check total nodes expected from : " + itoa(src) + " = " + itoa(_ExpectSrc[iOct][src]));
	
	
	if (totalNodesExpected  == recvInfosCpt[src])
	{
		recvFFTermsCpt[src]	= 0;
		recvInfosCpt[src]	= 0;
		//debug("info", "raz recv counters from : " + itoa(src));
	}
}

/// OLD VERSIONS ******************************
void Gaspi_FF_communicator::send_task_ff_dbg2_last(int level, complex * ff, int iOct, int start, int stop)
{
	int indexToC = -1;
	int termsSize = 0;
	int nbTerms = 0;
	int nbInfos = 0;
	bool levelIsFinished = false;
	bool readyToSend = false;
	//int tid = omp_get_thread_num();

	// variables précalculables
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
	// for each dest
	for (int dest=0; dest<_wsize; dest++)
	{
		levelIsFinished = false;

		if (dest != _rank)
		{
			// nb of cells to be sent
			int counter =  _count_send[iOct][dest][level];

			// if something to send to this dest 
			if (counter>0)
			{
				int lastBox =  _start_send[iOct][dest][level];	// the big one
				int firstBox = _stop_send[iOct][dest][level];	// the small one

				// traverses the cells to send
				for (int i= firstBox; i <= lastBox; i++)
				{

					int cellIDf = _send[iOct][i];
					readyToSend = false;
					
					// if the cell is in the task's nodes range
					if (cellIDf >= start && cellIDf <= stop)
					{
						int cellID = cellIDf + indexToC;
						int p = _fniv[iOct][level+1]+(_endlev[iOct][level]+indexToC-cellID)*__nstnsp;	// cell's corresponding ff index in Fortran

/* --- LOCK MUTEX --- */  
pthread_mutex_lock(&mutex);
						// read and update adresses 
						int q = _FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_counter[iOct][dest];						// in FF sendbuffer
						int idxInfos =_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_counter[iOct][dest];			// in Infos sendBuffer
						_FF_sendLocalOffsets_counter[iOct][dest] += __nstnsp;
						_Infos_sendLocalOffsets_counter[iOct][dest] += 1;
						_accumul[dest] += 1;
						
						// update send threashold infos
						nbTerms = _FF_sendLocalOffsets_counter[iOct][dest]-_FF_sendLocalOffsets_keeper[iOct][dest]; 
						termsSize = nbTerms * sizeof(complex);
						nbInfos = _Infos_sendLocalOffsets_counter[iOct][dest]-_Infos_sendLocalOffsets_keeper[iOct][dest]; 

						// test if ready to Send
						if (_accumul[dest] == counter) // tout le level est écrit
						{	
							levelIsFinished = true;
							readyToSend = true;
							_accumul[dest] = 0;
							//_write_counters_ptr[iOct][dest][level] = 0; // raz for next iteration
						}
						else if (termsSize >= _gaspiChunkSize) // comm_size atteint
						{
							readyToSend = true;
							//_write_counters_ptr[iOct][dest][level]++; // update for next packet
						}
						
						// if send --> update keepers and counters pointer
						if (readyToSend)
						{
							_FF_sendLocalOffsets_keeper[iOct][dest] = _FF_sendLocalOffsets_counter[iOct][dest];
							_Infos_sendLocalOffsets_keeper[iOct][dest] = _Infos_sendLocalOffsets_counter[iOct][dest];
						}
						
						// write FF and Infos
//						_FF_sendBuffer[q:__nstnsp] = ff[p:__nstnsp];
						for (int i=0; i<__nstnsp; i++)
							_FF_sendBuffer[q+i] = ff[p+i];


						
						// encode level + index
						int code = level << 26 | (i-firstBox);
						//int decodeLevel = code >> 26;
//						int mask = 0x03FFFFFF;
						//int decodeIndex = code & mask;
						//int index = i-firstBox;
						
						/*	debug("encode", 
						"index : " 			+ itoa(index)
						+ " level : " 		+ itoa(level)
						//+ " code : "		+ itoa();
						+ " decodeLevel : "	+ itoa(decodeLevel)
						+ " decodeIndex : " + itoa(decodeIndex));*/

						_Infos_sendbuffer[idxInfos] = code /*i-firstBox*/;	
						
						//_write_counters[iOct][dest][level][cptIndex]++;
						if (readyToSend)
						{
							_SendInfosCpt[iOct][dest] += nbInfos;							
							
							///send_chunk_dbg(iOct, dest, level, nbTerms, nbInfos, levelIsFinished);
							/// idee : faire le même envoi, n'ajouter que encodage/decodage
							/// validé, ==> checkpoint B
							
							send_chunk_dbg_last(iOct, dest, level, nbTerms, nbInfos, levelIsFinished);
							/// Envoi toujours la notif quand le level est fini
							/// Modif ID : src
							/// Value : qty
							/// checkpoint C : inifinity 456
							
						}
/* --- LOCK MUTEX --- */ 
pthread_mutex_unlock(&mutex);
					} // if start < cell < stop
				} // for cells
			} // if counter >0
		} // if dest != rank
	} // for dest
}

void Gaspi_FF_communicator::send_chunk_dbg_last(int iOct, int dest, int level, int nbTerms, int nbInfos, bool levelIsFinished)
{
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	flush_queues(_nbQueues);
			
	// local FF send offset
	gaspi_offset_t local_dest_offset = _FF_sendLocalOffsets[iOct][dest] * sizeof(complex);
	gaspi_offset_t accumul_offset = _FF_sendRemoteOffsets_counter[iOct][dest] * sizeof(complex);
	gaspi_offset_t local_offset_ff = local_dest_offset + accumul_offset;
	gaspi_queue_id_t queue=0;
	
	// remote offsets FF et Infos
	gaspi_offset_t remote_sender_offset = _FF_sendRemoteOffsets[octree_offset + dest] * sizeof(complex);
	gaspi_offset_t remote_offset_ff = remote_sender_offset + accumul_offset;
	_FF_sendRemoteOffsets_counter[iOct][dest] += nbTerms;//
	
	gaspi_offset_t remote_offset_infos = (_Infos_sendRemoteOffsets[octree_offset + dest] +_Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int) ;
	gaspi_offset_t local_offset_infos = (_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(int);
	_Infos_sendRemoteOffsets_counter[iOct][dest] += nbInfos;//
	
	gaspi_size_t qtyFF = nbTerms * sizeof(complex);
	gaspi_size_t qtyInfos = nbInfos * sizeof(int);
	
	// envoi des DATA sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_FF_sendBuf_seg_id,					// local seg ID
			local_offset_ff,					// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_offset_ff,						// remote offset
			qtyFF,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
	// envoi des INFOS sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_offset_infos,					// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_offset_infos,				// remote offset
			qtyInfos,							// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	
	if (levelIsFinished)
	{
		// ENVOI DE LA NOTIFICATION
		//gaspi_notification_id_t notifyID = (_rank * 65536/_wsize) + level;
		//gaspi_notification_t notifyValue = level;
	
		/// prêt à changer le notif ID ici
		int myMin = _rank * 65536/_wsize;
		_notifCpt = (_notifCpt + 1) % (65536/_wsize);
		gaspi_notification_id_t notifyID2 = myMin + _notifCpt;
		gaspi_notification_t notifyValue2 = _SendInfosCpt[iOct][dest];

		//~ debug("write_verif",
			  //~ "write -----------------------------------\nreal src : " 
			//~ + itoa(_rank) 
			//~ + ", decoded src : " + itoa(notifyID2/(65536/_wsize))
			//~ + ", dest : " + itoa(dest)
			//~ + "\nlevel : " + itoa(level)); 
			//~ //+ "\nnbInfos : " + itoa(notifyValue2));
		//~ debug("write_verif", "_Infos_sendRemoteOffsets_counter : " + itoa(_Infos_sendRemoteOffsets_counter[iOct][dest]));
		//~ debug("write_verif", "_SendInfosCpt : " + itoa(_SendInfosCpt[iOct][dest]));


	
		SUCCESS_OR_DIE(
			gaspi_notify( 
				_Infos_recvbuf_seg_id,				// remote seg id
				dest,								// receiver rank
				notifyID2,
				notifyValue2,
				queue,								// queue
				GASPI_BLOCK							// Gaspi block
			)
		);
		
		_SendInfosCpt[iOct][dest] = 0;
		
		// RAZ DES COMPTEURS lorsqu'on a fini d'envoyer
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
}

void Gaspi_FF_communicator::send_chunk_dbg_last_v2(int iOct, int dest, int nbTerms, int nbInfos, bool destIsFinished)
{
	//int indexToC = -1;
	int octree_offset = iOct * _wsize;
	flush_queues(_nbQueues);
			
	// FF 
	gaspi_queue_id_t queue=0;
	gaspi_offset_t local_ff_offset = (_FF_sendLocalOffsets[iOct][dest] + _FF_sendLocalOffsets_keeper[iOct][dest]) * sizeof(complex);
	gaspi_offset_t remote_ff_offset = (_FF_sendRemoteOffsets[octree_offset + dest] + _FF_sendRemoteOffsets_counter[iOct][dest]) * sizeof(complex);
	gaspi_size_t qtyFF = nbTerms * sizeof(complex);

	debug("gaspi_write",
		  "write -----------------------------------\ndest : " 
		+ itoa(dest)
		+ "\nlocal ff offset : "
		+ itoa(local_ff_offset)
		+ "\nlocal keeper : "
		+ itoa( _FF_sendLocalOffsets_keeper[iOct][dest])
		+ "\nremote ff offset : "
		+ itoa(remote_ff_offset)
		+ "\nremote counter : "
		+ itoa(_FF_sendRemoteOffsets_counter[iOct][dest])
		+ "\nnbterms : "
		+ itoa(nbTerms)
		+ "\nqtyFF : "
		+ itoa(qtyFF)
		);

	
	// envoi des FF sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_FF_sendBuf_seg_id,					// local seg ID
			local_ff_offset,					// local offset
			dest,								// receiver rank
			_FF_recvBuf_seg_id,					// remote seg ID
			remote_ff_offset,						// remote offset
			qtyFF,								// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
	debug("gaspi_write",
		  "write -----------------------------------dest : " 
		+ itoa(dest)
		+ "ok"		
		);

	// Infos
	gaspi_offset_t local_Infos_offset = (_Infos_sendLocalOffsets[iOct][dest] + _Infos_sendLocalOffsets_keeper[iOct][dest]) * sizeof(complex);
	gaspi_offset_t remote_Infos_offset = (_Infos_sendRemoteOffsets[octree_offset + dest] + _Infos_sendRemoteOffsets_counter[iOct][dest]) * sizeof(complex);	
	gaspi_size_t qtyInfos = nbInfos * sizeof(int);
	
	// envoi des INFOS sans notification
	SUCCESS_OR_DIE(
		gaspi_write( 
			_Infos_sendbuf_seg_id,				// local seg ID
			local_Infos_offset,					// local offset
			dest,								// receiver rank
			_Infos_recvbuf_seg_id,				// remote seg ID
			remote_Infos_offset,				// remote offset
			qtyInfos,							// size of data to write
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));


	// envoi de la notification	
	int myMin = _rank * 65536/_wsize;
	_notifCpt = (_notifCpt + 1) % (65536/_wsize);
	gaspi_notification_id_t notifyID = myMin + _notifCpt;
	gaspi_notification_t notifyValue = nbInfos;

	//~ /**debug("write_verif",
		  //~ "write -----------------------------------\nreal src : " 
		//~ + itoa(_rank) 
		//~ + ", decoded src : " + itoa(notifyID2/(65536/_wsize))
		//~ + ", dest : " + itoa(dest)
		//~ + "\nlevel : " + itoa(level)); 
	//~ debug("write_verif", "_Infos_sendRemoteOffsets_counter : " + itoa(_Infos_sendRemoteOffsets_counter[iOct][dest]));
	//~ debug("write_verif", "_SendInfosCpt : " + itoa(nbInfos]));**/

	SUCCESS_OR_DIE(
		gaspi_notify( 
			_Infos_recvbuf_seg_id,				// remote seg id
			dest,								// receiver rank
			notifyID,							// src rank
			notifyValue,						// qty
			queue,								// queue
			GASPI_BLOCK							// Gaspi block
		)
	);
	SUCCESS_OR_DIE (gaspi_wait (queue, GASPI_BLOCK));
		
	// update des compteurs	
	
	/// send counters
	_SendInfosCpt[iOct][dest] = 0;
	
	if (destIsFinished)
	{
		/// send address keepers
		_FF_sendLocalOffsets_keeper[iOct][dest] = 0;
		_Infos_sendLocalOffsets_keeper[iOct][dest] = 0;
		
		/// remote counters
		_FF_sendRemoteOffsets_counter[iOct][dest] = 0;
		_Infos_sendRemoteOffsets_counter[iOct][dest] = 0;
	}
	else
	{
		/// send address keepers	
		_FF_sendLocalOffsets_keeper[iOct][dest] += nbTerms;
		_Infos_sendLocalOffsets_keeper[iOct][dest] += nbInfos;	
		
		/// remote counters
		_FF_sendRemoteOffsets_counter[iOct][dest] += nbTerms;
		_Infos_sendRemoteOffsets_counter[iOct][dest] += nbInfos;		
	}
}

/* attend 1 message par source, pour 1 level donné */
void Gaspi_FF_communicator::recv_task_ff_last(int level, complex * ff, int iOct)
{
	//double t_begin, t_end;
	//t_begin = MPI_Wtime();
	
	// compute nb of nodes to receive at current level
	int nbSrcs = 0;
	for (int i=0; i<_wsize; i++)
	{
		if (_Expect[iOct][i][level]>0)
		{
			nbSrcs++; // ajoute 1 source
		}		
	}

	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;
	//int src;
	//int nbNodes;
	int nbMsg = 0;	
	

	while (nbMsg != nbSrcs)
	{
		for (int src=0; src<_wsize; src++)
		{
			if(gaspi_notify_waitsome( // test level L for src I
				_Infos_recvbuf_seg_id,
				(65536/_wsize * src) + level,
				1,
				&new_notif_id,
				GASPI_TEST) == GASPI_SUCCESS)
			{
				gaspi_notify_reset(
					_Infos_recvbuf_seg_id, 
					new_notif_id, 
					&new_notif_val
				);				
				
				updateFarFieldChunksFromInfos_last(src, level, ff, _Expect[iOct][src][level], iOct);
				/// cette version est ok ---> checkpoint B 
				
				///updateFarFieldsFromInfosEncoded(src, level, ff, _Expect[iOct][src][level], iOct);
				/// cette version, plus ancienne, est ok aussi
				nbMsg++;
			}
		}		
	}
	
	// raz for next time
	/*for (int i=0; i<_wsize; i++)
	{
		_Received[iOct][i][level] = 0;
	}*/
	
	//t_end = MPI_Wtime();
//~ #ifdef TIMING
	//~ add_time_sec("GASPI_FF_sendrecv", t_end - t_begin - accumul_r);
//~ #endif
}

void Gaspi_FF_communicator::updateFarFieldChunksFromInfos_last(int src, int level, complex * ff, int nbNodes, int iOct)
{
	//debug("info", "in --> update FF, received : " + itoa(nbNodes) + " nodes, from : " + itoa(src) + " level : " + itoa(level));
	
	int indexToC = -1;
	int octree_offset = iOct * _wsize;
	int __nstnsp = _nst[iOct][level]*_nsp[iOct][level];
	
/* ----- AJOUTS ----- */
	int * recvInfosCpt   = _ReceivedInfosCpt[iOct];	//_ReceivedInfosBytes[iOct];
	int * recvFFTermsCpt = _ReceivedFFTermsCpt[iOct];	//_ReceivedFFOffsets[iOct];
	
	//debug("info", "nb Infos already received from : " + itoa(src) + " = " + itoa(recvInfosCpt[src]));
	
	// nouveau calcul des offsets du niveau 
	//int info_address = _Infos_recvOffsets[octree_offset + src]/*src*/ 
	//				 + recvInfosCpt[src]; /*déjà reçus*/;
					
	//int ff_address	 = _FF_recvOffsets[octree_offset + src] 
		//			 + recvFFTermsCpt[src];
	
	// update les compteurs d'offset
	recvInfosCpt[src] 	+= nbNodes;
	recvFFTermsCpt[src] += nbNodes * __nstnsp;

	// RAZ si fini
	// calcul nb total de nodes à recevoir pour src
	int totalNodesExpected = 0;
	for (int lvl=0; lvl<_nivterm[iOct]; lvl++)
		totalNodesExpected += _Expect[iOct][src][lvl];

	//debug("info", "total nodes expected from : " + itoa(src) + " = " + itoa(totalNodesExpected));
	
	if (totalNodesExpected  == recvInfosCpt[src])
	{
		recvFFTermsCpt[src]	= 0;
		recvInfosCpt[src]	= 0;
		//debug("info", "raz recv counters from : " + itoa(src));
	}

	//debug("info", "Infos, read @ address : " + itoa(info_address));

/* ----- fin AJOUTS ----- */

	// calcul de l'offset du niveau
	int info_level_offset = 0;
	int ff_level_offset = 0;
	for (int i=_nivterm[iOct]-1; i>level; i--)			// somme tout ce qui est attendu avant
	{
		info_level_offset += _Expect[iOct][src][i];
		ff_level_offset += _Expect[iOct][src][i]*_nst[iOct][i]*_nsp[iOct][i];

	}
	
	int info_offset = _Infos_recvOffsets[octree_offset + src]/*src*/ 
					+ info_level_offset/*lvl*/ 
					+ _Received[iOct][src][level] /*déjà reçus*/; 
	

	int ff_offset = _FF_recvOffsets[octree_offset + src] 
					+ ff_level_offset 
					+ _Received[iOct][src][level]*__nstnsp;

	//~ debug("info", "src offst : " + itoa(_Infos_recvOffsets[octree_offset + src]) 
								 //~ + " lvl offset : " + itoa(info_level_offset) 
								 //~ + " already received offset : " + itoa(_Received[iOct][src][level]));

	int firstBox = _stop_recv[iOct][src][level];
	
	int index;

	int __endlev = _endlev[iOct][level]+indexToC;
	int __fnivnextlev = _fniv[iOct][level+1];
	//int ff_index;
	
	for (int i=info_offset; i<info_offset+nbNodes; i++)
	{
		//index = _Infos_recvbuffer[i];		
		int code = _Infos_recvbuffer[i];
		int mask = 0x03FFFFFF;
		index = code & mask;
		int cellIDf = _recv[iOct][firstBox+index];
		int cellID = cellIDf + indexToC;

		// ranger la boîte
		int ff_index = (i-info_offset) * __nstnsp;
		int p0 = __fnivnextlev + ((__endlev-cellID)*__nstnsp);		
		
		//ff[p0:__nstnsp] = ff[p0:__nstnsp] + _FF_recvBuffer[ff_offset+ff_index:__nstnsp];
		for (int i=0; i<__nstnsp; i++)	
			ff[p0+i] = ff[p0+i] + _FF_recvBuffer[ff_offset+ff_index+i];


		debug("sort", "cellID : " + itoa(cellID) + " level " + itoa(level) + ", ff@ : " + itoa(p0));
		
	}
	
	
				/// update recv counters
			/// RAZ counters, when ?
			
			
			//~ _Received[iOct][src][recvLevel] += counter;
			//~ if (recvLevel == level)
				//~ nbReceived += counter;
			
			//~ /// RAZ compteur pour iteration suivantes
			//~ if (_Received[iOct][src][recvLevel] == _Expect[iOct][src][recvLevel])
				//~ _Received[iOct][src][recvLevel] = 0;
}

