/*
  Copyright 2015 - UVSQ
  Authors list: Nathalie M�ller, Eric Petit

  This file is part of the FMM-lib.

  FMM-lib is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later version.

  FMM-lib is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License along with
  the FMM-lib. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FMM_HPP
#define FMM_HPP

#include "mpi.h"
#include "GASPI.h"
#include "../Tools/Complex.hpp"
#include "../Gaspi/Gaspi_FF_communicator.hpp"
//#include "measure.hpp"

extern "C"
{
	// start stop, switch 
	void fmm_gaspi_init_();
	void fmm_gaspi_finalize_();
	void fmm_switch_to_gaspi_();
	void fmm_switch_to_mpi_();
	
	// initialize communicator and create segments
	void init_gaspi_ff_communicator_(i64 * recvnode,	i64 * recvnode_sz,
		i64 * sendnode, i64 * sendnode_sz, i64 * nb_recv, 	i64 * nb_recv_sz,
		i64 * nb_send, i64 * nb_send_sz, i64 * nivterm,i64 * levcom,
		i64 * fniv, i64 * nst, i64 * nsp, i64 * fsend, i64 * send,
		i64 * frecv, i64 * recv, i64 * endlev, i64 * codech, i64 * includeLevcom);

	// ff 
	void fmm_handle_ff_gaspi_bulk_(complex * ff, complex * bufsave, i64 * idom);
	void gaspi_send_ff_(i64 * niv, complex * ff, i64 * idom);
	void gaspi_recv_ff_(i64 * niv, complex * ff, i64 * idom);
	void gaspi_send_ff_multi_(i64 * niv, complex * ff, i64 * idom);
	void gaspi_recv_ff_multi_(i64 * niv, complex * ff, i64 * idom);	
	void gaspi_task_send_ff_(i64 * start, i64 * stop, i64 * niv, complex * ff, i64 * idom);
	void gaspi_task_recv_ff_(i64 * niv, complex * ff, i64 * idom);
	
	// chunks
	void gaspi_task_chunk_send_(i64 * start, i64 * stop, i64 * niv, complex * ff, i64 * idom);
	void gaspi_task_chunk_recv_(i64 * niv, complex * ff, i64 * idom);

	// total tasks
	void gaspi_send_chunk(i64 * start, i64 * stop, i64 * niv, complex * ff, i64 * idom);

	// unknowns
	void fmm_handle_unknowns_broadcast_(complex * xtmp, complex * xtmp2, i64 * size);
	void fmm_handle_unknowns_allreduce_();
	
	// debug tools
	void fmm_dump_cplx_(complex * tab, i64 * size, i64 * fileNum);
	void fmm_dump_i8_(i64 * tab, i64* size, i64 * fileNum);
	void fmm_dump_2ble_(double * tab, i64 * size, i64 * fileNum);
	//
	void fmm_raz_i8_(i64 * tab, i64* size);


	// multimat
	void gaspi_init_ff_(i64 * max_send_terms, i64 * max_recv_terms, i64 * nbMat, i64 * max_send_nodes, i64 * max_recv_nodes, i64 * includeLevcom);
	void gaspi_init_offsets_(i64 * recvnode, i64 * recvnode_sz, i64 * sendnode, i64 * sendnode_sz, 
		i64 * nb_recv, i64 * nb_recv_sz, i64 * nb_send, i64 * nb_send_sz, i64 * idom, i64 * ndom,
		i64 * nivterm, i64 * frecv, i64 * recv, i64 * levcom, i64 * endlev, i64 * fniv, i64 * fsend, i64 * send,
		i64 * nst, i64 * nsp, i64 * codech);
	void fmm_finalize_dump_vector_();
	
	// random MPI � sortir dans lib � part
	void randomizempi_(i64 * wsize, i64 * procIDs);
}

#endif
