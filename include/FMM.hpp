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

#include <string>
using namespace std;

// Load Balancing strategy
#define HIST_EXACT 0
#define HIST_APPROX 1
#define	MORTON_MPI_SYNC 2
#define MORTON_GASPI_ASYNC 3

void FMM_load_balance(string file, int nbParticles, double dist, double tolerance, int LBStrategy);

#endif
