#include "Gaspi_M2L_communicator.hpp"
using namespace std;


Gaspi_m2l_communicator::Gaspi_m2l_communicator(i64 * nb_recv, int nb_recv_sz, i64 * recvnode, int recvnode_sz)	
{
	// update class attributes
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);
		
	// update segments ids
	_seg_globalRecvBuffer_id 	    = 7; // Receive Buffer 
    _seg_globalRecvBufIdxPerRank_id = 8; // Index in global recv buffer, per RANK
    _seg_remoteBufferIndexes_id     = 9; // Index where to write on other ranks
	
	// create gaspi segments, and initialize them
	init_globalRecvBuffer(nb_recv, nb_recv_sz);
	init_remoteBufferIndexes(recvnode, recvnode_sz, nb_recv);
}

void Gaspi_m2l_communicator::init_globalRecvBuffer(i64 * nb_recv, int nb_recv_sz)
{
	// compute size
	int globalRecvBufferSize = 0;
	for (int i=0; i<nb_recv_sz; i++)
 	    globalRecvBufferSize = globalRecvBufferSize + nb_recv[i];
 	    
 	// update segment size
 	_seg_globalRecvBuffer_size = globalRecvBufferSize * sizeof(complex);
 	
 	// create segment
 	gaspi_segment_create(
		_seg_globalRecvBuffer_id,
		_seg_globalRecvBuffer_size,
		GASPI_GROUP_ALL,
		GASPI_BLOCK, 
		GASPI_ALLOC_DEFAULT
    );

	// update gaspi segpointer
	gaspi_segment_ptr(_seg_globalRecvBuffer_id, &_ptr_seg_globalRecvBuffer);
	
	// update user pointer
	_globalRecvBuffer = (complex *)_ptr_seg_globalRecvBuffer;
    
    cout << "--- Global Recv Segment has been created. ---" << endl;
}

void Gaspi_m2l_communicator::init_remoteBufferIndexes(i64 * recvnode, int recvnode_sz, i64 * nb_recv)
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
	gaspi_segment_create(
		_seg_globalRecvBufIdxPerRank_id,
		_seg_globalRecvBufIdxPerRank_size,
		GASPI_GROUP_ALL,
		GASPI_BLOCK,
		GASPI_ALLOC_DEFAULT
    );
    
    // update gaspi segpointers
    gaspi_segment_ptr( _seg_globalRecvBufIdxPerRank_id, &_ptr_seg_globalRecvBufIdxPerRank);
	
	// update user pointers
    _globalRecvBufIdxPerRank = (int *)_ptr_seg_globalRecvBufIdxPerRank;    
    
    // initialize segment with values
    for (int k=0; k<recvnode_sz; k++)
    {
        int from = recvnode[k]-1;
        if (from != _rank)
        {
			_globalRecvBufIdxPerRank[from] = globalRecvBufferIndexPerROUND[k];
		}
		else
		{
			cerr << "Error in Gaspi_m2l_communicator::init_remoteBufferIndexes.\nFortran recvnode array should not contain the current rank !";
			exit(0);
		}
			
		_globalRecvBufIdxPerRank[_rank] = -1;
    }

    cout << _rank << " receptions in RANK order." << endl;
    for (int i=0; i<_wsize; i++)
    {
		cout << "From : " << i << " at index : " << _globalRecvBufIdxPerRank[i] << ", qtty : " << nb_recv[i] << endl;
	}
    /*
    cout << "receptions in ROUND order." << endl;
    for (int k=0; k<recvnode_sz; k++)
    {
        int from = recvnode[k]-1;   
        cout << "n° " << k << " - from : " << from << " at index : " << globalRecvBufferIndexPerROUND[k] << " - qtty : " << nb_recv[from] << endl;
    }*/
    
    /* 
     * Prepare Remote Buffer Indexes SEGMENT 
     */
     
    // update segment size
    _seg_remoteBufferIndexes_size = _wsize * sizeof(int);
	
	// create gaspi segment
    gaspi_segment_create(
		_seg_remoteBufferIndexes_id, 
		_seg_remoteBufferIndexes_size,
		GASPI_GROUP_ALL, 
		GASPI_BLOCK, 
		GASPI_ALLOC_DEFAULT
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
			);
		}
	}
	
	cout << "sent all the write notify" << endl;
	
	// wait to receive all infos
	int cpt = 0;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;

	while (cpt < (_wsize-1))
	{
		// wait for notification
		gaspi_notify_waitsome(
			_seg_remoteBufferIndexes_id,
			0,							// surveille les notifications depuis 0
			_wsize,						// en surveille wsize
			&new_notif_id,
			GASPI_TEST
		);

		// get notification value, and reset
		gaspi_notify_reset(
			_seg_remoteBufferIndexes_id,
			new_notif_id,
			&new_notif_val
		);
		
		if (new_notif_val == REMOTE_ADDRESS)
			cpt++;
	}
    
	cout <<"received all notifications." << endl;
	
	displayTab<int>("remote addresses", _remoteBufferIndexes, _wsize);
}

void init_globalSendBuffer()
{
}



void init_dataToSendIndexes()
{
}

void init_gaspi_m2l_segments(i64 * nb_recv, int nb_recv_sz, i64 * recvnode, int recvnode_sz)
{
    cout <<"----------------- computing data to prepare Gaspi ARRAYS --------" << endl;
    int indexToC = -1;
        
    // Global Receive Buffer
    // init_globalRecvBuffer(nb_recv, nb_recv_sz);

    // Global Send Buffer
    init_globalSendBuffer();

    // Remote Buffer Indexes


    // Adresses in send buffer, for sending
    init_dataToSendIndexes();
    
    // update distant buffer indexes
/*
    gaspi_write_notify( // local seg
                        // local offset
                        // receiver rank
                        // remote seg
                        // remote offset
                        // size of data to write
                        // remote notif ID
                        // value of the notif to write
                        // queue
                        // Gaspi block
    );
*/	
}


void create_gaspi_m2l_segments(i64 * nb_send, int nb_send_sz, 
							 i64 * nb_recv, int nb_recv_sz,
							 i64 * sendnode, int sendnode_sz,
							 i64 * recvnode, int recvnode_sz,
							 int levcom, int nivterm,
							 i64 * fsend, i64 * send, i64 * endlev, i64 * codech,
							 i64 * nst, i64 * nsp, complex * bufsave, i64 * fniv,
							 complex * ff,
                             Gaspi_m2l_communicator *& gCommM2L)

{
	// Passage en Gaspi    
    MPI_Barrier(MPI_COMM_WORLD);
    gaspi_rank_t rank, wsize;
    SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
    SUCCESS_OR_DIE(gaspi_proc_num(&wsize));

    // initialize gaspi segments
    gCommM2L = new Gaspi_m2l_communicator(nb_recv, nb_recv_sz, recvnode, recvnode_sz);
    //displayTab<int>("buffer_indexes", (*gCommM2L)._remoteBufferIndexes, wsize);
    //init_gaspi_m2l_segments(nb_recv, nb_recv_sz, recvnode, recvnode_sz);

/*	int indexToC = -1;	
	
	// Alloc global send Buffer
	int globalSendBufferSize = 0;
	for (int i=0; i<nb_send_sz; i++)
	    globalSendBufferSize = globalSendBufferSize + nb_send[i];
	complex * globalSendBuffer = nullptr;
	globalSendBuffer = new complex [globalSendBufferSize]();

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
		//cout << rank <<  " Remplissage du buffer d'envoi, i = " << i << " / " << sendnode_sz << endl;
		int iDest = sendnode[i];
		int k = nivterm;
		int q = globalSendBufferIndex[i]-1;
		//printf("[%d] - idest: %d, k: %d, q: %d\n",rank, iDest, k, q);
		
		if (iDest > 0)
		{
			//printf("[%d] - bornes de j : [ %ld - %ld ]\n",rank, fsend[iDest+1+indexToC]-1, fsend[iDest+indexToC]);
			for (int j=fsend[iDest+1+indexToC]-1; j>=fsend[iDest+indexToC]; j--)
			{
				int jc = j + indexToC;// /!\ j est un numero de cellule -> -1 pour passer au C
				
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
								globalSendBuffer[q]=bufsave[p + indexToC];
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
								globalSendBuffer[q]=ff[p + indexToC];
							}
						}
					}
				}
			}
		}
	}
	
	// compteurs de messages
	int nbMsgToRecv = 0;
	for (int ip=0; ip<recvnode_sz; ip++)
	{
	    int i = recvnode[ip];
	    if (i > 0)
			nbMsgToRecv++;
	}
	
	int nbMsgToSend = 0;
	for (int ip=0; ip<sendnode_sz; ip++)
	{
	    int i = sendnode[ip];
	    if (i > 0) 
			nbMsgToSend++;
	}

	//displayTab<i64>("nb_send", nb_send, nb_send_sz);


	// dealloc
	delete [] globalSendBuffer;
	delete [] globalSendBufferIndex;
*/
    // Rend la main au MPI
    SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
}
