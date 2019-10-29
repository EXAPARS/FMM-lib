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

#include "../Tools/types.hpp"
#include "../Tools/fmm_tools.hpp"
#include "../LB/vec3D.hpp"
#include "../LB/Particles.hpp"
#include "../LB/Node.hpp"
#include "../LB/Decomposition.hpp"
#include "../LB/LoadBalancerBase.hpp"
#include "../LB/LoadBalancer.hpp"
#include "../LB/LBMortonSyncMPI.hpp"
#include "../LB/LBHistApprox.hpp"


//#include "measure.hpp"



extern "C"
{
	// Load Balancing
	void fmm_load_balance_(	i64 * nbElemPerNode, i64 * firstElemAdress, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nodeOwners, 
		double * nodeCenters, i64 * endlev, i64 * nbLevels, double * maxEdge, int * LBstrategy);
	
	void fmm_get_elem_coors_ (int * elemToNode, i64 * nbElem, double * nodesXcoords, double * nodesYcoords, double * nodesZcoords);
}

#endif
