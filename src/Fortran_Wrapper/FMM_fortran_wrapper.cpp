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

#include <string>
#include <iostream>

#include "FMM_fortran_wrapper.hpp"
using namespace std;

// Load Balancing strategy
#define HIST_EXACT 0
#define HIST_APPROX 1
#define	MORTON_MPI_SYNC 2
#define MORTON_GASPI_ASYNC 3

void fmm_load_balance_(i64 * nbElemPerNode, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nbNodes, i64 * nodeOwner)
{
	// Dump octree info with a recursive DFS traversal
	dfs_dump_spectre_octree(nbElemPerNode, nbSonsPerNode, firstSonId, nbNodes, nodeOwner, 0);

}

