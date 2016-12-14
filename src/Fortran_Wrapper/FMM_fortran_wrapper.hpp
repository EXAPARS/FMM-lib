/*
  Copyright 2015 - UVSQ
  Authors list: Nathalie Möller, Eric Petit

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
#include "Gaspi_tools.hpp"

#include "types.hpp"
#include "fmm_tools.hpp"
#include "Complex.hpp"
#include "vec3D.hpp"
#include "Particles.hpp"
#include "Node.hpp"
#include "Decomposition.hpp"
#include "LoadBalancerBase.hpp"
#include "LoadBalancer.hpp"
#include "LBMortonSyncMPI.hpp"
#include "LBHistApprox.hpp"
#include "Gaspi_M2L_communicator.hpp"
#include "Gaspi_UNK_communicator.hpp"


#include "/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/measure.hpp"



extern "C"
{
	// Load Balancing
	void fmm_load_balance_(	i64 * nbElemPerNode, i64 * firstElemAdress, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nodeOwners, 
		double * nodeCenters, i64 * endlev, i64 * nbLevels, double * maxEdge, int * LBstrategy);
	
	void fmm_get_elem_coors_ (int * elemToNode, i64 * nbElem, double * nodesXcoords, double * nodesYcoords, double * nodesZcoords);
	
	// Gaspi
	void fmm_gaspi_init_();
	void fmm_gaspi_finalize_();
	
	// unknowns
	void fmm_handle_unknowns_broadcast_(complex * xtmp, complex * xtmp2, i64 * size);
	void fmm_handle_unknowns_allreduce_();
	
	// ff
	void init_gaspi_ff_communicator_(i64 * recvnode,	i64 * recvnode_sz,
										  i64 * sendnode, 	i64 * sendnode_sz,
										  i64 * nb_recv, 	i64 * nb_recv_sz,
										  i64 * nb_send, 	i64 * nb_send_sz,
										  i64 * nivterm,	i64 * levcom,
										  i64 * fniv,		i64 * nst,			i64 * nsp,
										  complex * ff, complex * ne, i64 * allreduce_sz,
										  i64 * fsend, 		i64 * send,
										  i64 * frecv,		i64 * recv,
										  i64 * endlev,		i64 * codech);
										    
	void fmm_handle_allreduce_gaspi_(complex * ff, complex * ne, i64 * size, 
							i64 * recvnode, i64 * recvnode_sz, 
							i64 * sendnode, i64 * sendnode_sz,
							i64 * nb_recv, 	i64 * nb_recv_sz,
							i64 * nb_send, 	i64 * nb_send_sz);
							
							
	void fmm_handle_comms_gaspi_(i64 * recvnode, 	i64 * recvnode_sz, 
								 i64 * sendnode, 	i64 * sendnode_sz,
								 i64 * nb_recv, 	i64 * nb_recv_sz,
								 i64 * nb_send, 	i64 * nb_send_sz,
								 i64 * nivterm, 
								 i64 * levcom, 
								 i64 * fniv, 
								 i64 * nst, 
								 i64 * nsp,
								 complex * ff,
								 i64 * fsend, 		i64 * send,
								 i64 * frecv,		i64 * recv,
								 i64 * endlev,		i64 * codech,
								 complex * bufsave);
	
	
	void fmm_gaspi_init_handler_();
	void fmm_gaspi_destroy_handler_();
	
	// debug tools
	void fmm_dump_(complex * tab);

}

#endif
