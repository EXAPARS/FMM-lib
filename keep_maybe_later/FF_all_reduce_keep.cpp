
/*    GASPI ALLREDUCE FF    */
void Gaspi_m2l_communicator::create_allReduceBuffers(complex * ff, complex * ne, int nbEltsToReduce)
{
	// use Fortran arrays for FF and NE
	SUCCESS_OR_DIE(
		gaspi_segment_use(
			_seg_ff_allreduce_id,
			ff,
			nbEltsToReduce*sizeof(complex),
			//_ff_sz*sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			0
		)
	);
	
	SUCCESS_OR_DIE(
		gaspi_segment_use(
			_seg_ne_allreduce_id,
			ne,
			nbEltsToReduce*sizeof(complex),
			GASPI_GROUP_ALL,
			GASPI_BLOCK, 
			0
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
	_seg_globalRecvBuffer_size = nbEltsToReduce * sizeof(complex) * _wsize;
}


void Gaspi_m2l_communicator::initAllReduceBuffers(complex * ff, complex * ne)
{
	//int nbElts = _seg_reduce_size / sizeof(complex);
	int nbElts = _nbEltsToReduce;
	
	// ff
	for (int i=0; i<nbElts; i++)
		_reduceFF[i] = ff[i];
	
	// ne
	for (int i=0; i<nbElts; i++)
		_reduceNE[i] = 0;
}


void Gaspi_m2l_communicator::runM2LallReduce(complex * ff, complex * ne)
{
	// test à l'arrache réutilisation de code
	int nbQueues = 1;
	int localOffset = 0;
	int offsetMultiple = 1;
	int nbElts = _nbEltsToReduce;
	
	// SEND
	broadcast_to_global_buffer(
		nbQueues, localOffset, offsetMultiple, nbElts, sizeof(complex), 
		_rank, _wsize, _seg_ff_allreduce_id, _seg_globalRecvBuffer_id, ALLREDUCE, "GASPI_REDUCE_FF_write_notify");
	
	// RECV and reduce on _unknowns array
	copy_local_data<complex>(_reduceNE, _reduceFF, nbElts, "GASPI_REDUCE_FF");
	
	receive_allReduce(offsetMultiple, "GASPI_REDUCE_FF", nbElts, _wsize, _seg_globalRecvBuffer_id, ALLREDUCE, _reduceNE, _globalRecvBuffer);
}

/* Wrapper */
void fmm_handle_allreduce_gaspi_(complex * ff, complex * ne, i64 * size, 
							i64 * recvnode, i64 * recvnode_sz, 
							i64 * sendnode, i64 * sendnode_sz,
							i64 * nb_recv, 	i64 * nb_recv_sz,
							i64 * nb_send, 	i64 * nb_send_sz)
{
	cout << "HANDLE ALLREDUCE GASPI"<< endl;
	
	// Switch to gaspi
	double t_begin, t_end;
	t_begin = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_switch_interop", t_end - t_begin);


	// DEBUG
	gaspi_rank_t _wsize;
	gaspi_rank_t _rank;
	gaspi_proc_rank(&_rank);
	gaspi_proc_num(&_wsize);
	
	// First call : Class instantiation, allocations and Gaspi segment creation
	if (! gCommM2L)
	{	
		cout << "ATTENTION !!! ---> Appel au constructeur du communicator M2L" << endl;
		t_begin = MPI_Wtime();
    	construct_m2l_communicator(
							nb_send, (int)(*nb_send_sz), 
							nb_recv, (int)(*nb_recv_sz),
							sendnode, (int)(*sendnode_sz),
							recvnode, (int)(*recvnode_sz),
							0, 0,
							ff, ne, (int)(*size),
							nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0,
							gCommM2L);
		t_end = MPI_Wtime();
		add_time_sec("GASPI_M2L_init_class_and_create_segments", t_end - t_begin);
	}
	
	// run allReduce
	t_begin = MPI_Wtime();
	gCommM2L->runM2LallReduce(ff, ne);
	t_end = MPI_Wtime();
	add_time_sec("GASPI_allReduce_ff", t_end - t_begin);
}
