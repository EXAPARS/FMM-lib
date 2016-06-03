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

#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstdint>

/**
 * Enums
 */

enum LB				// Load Balancing 
{
	hist, 			// 0 : exact separators computation with histograms
	hist_approx, 	// 1 : approximate histograms separators on octree grid
	morton_dfs, 	// 2 : octree grid with morton loop-based sorting
	morton_dfs_async // 3 : octree grid with morton loop-based sorting, with asynchronous comms
};

/** 
 * Typedefs
 */

typedef uint32_t ui32;
typedef uint64_t ui64;
typedef int64_t i64;

/**
 * Constants
 */
  
#define COORDMAX 2147483648  // 2 ^ 31, max value on 32 bits is 2^32-1
#define H16_SIZE 65536
#define H12_SIZE 4096
#define H4_SIZE 16

 
/**
 * Macros
 **/
#define nullValue 0


/**
 * MPI TAGS
 **/

#define LEAVES 0
#define NB_ITEMS 1
#define SEPS 2
#define NB_UNTIL 3
#define BUFFER_REQUEST 4
#define SEPS_UPDATE 5
#define COUNTER_UPDATE 6
#define BUFFER_ANSWER 7
#define PARTICLES 8
#define NO_PARTICLES 9

/**
 * MACRO GASPI
 **/
 
#define SUCCESS_OR_DIE(f...)\
do\
{\
	const gaspi_return_t r = f;\
	if (r != GASPI_SUCCESS)\
	{\
        printf("SUCCESS_OR_DIE failed at [%s:%i]\n",__FILE__,__LINE__);\
        gaspi_rank_t rank;\
        gaspi_proc_rank(&rank);\
		printf ("Error rank %d: '%s' [%s:%i]: %i\n",rank, #f, __FILE__, __LINE__, r);\
		exit (EXIT_FAILURE);\
	}\
} while (0)

#define ASSERT(x...)                                                   		 \
do{\
if (!(x))                                                          		   \
{                                                                 		    \
	fprintf (stderr, "Error rank: '%s' [%s:%i]\n", #x, __FILE__, __LINE__);  \
	exit (EXIT_FAILURE);                                              			  \
}\
} while(0)

#endif

