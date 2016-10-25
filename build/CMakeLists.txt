#  Copyright 2015 - UVSQ
#  Authors list: Nathalie Möller, Eric Petit
#
#  This file is part of the FMM-lib.
#
#  FMM-lib is free software: you can redistribute it and/or modify it under the
#  terms of the GNU Lesser General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later version.
#
#  FMM-lib is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with
#  the FMM-lib. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.8.10)

PROJECT (FMM)
set (CMAKE_CXX_COMPILER icpc)

# variables
set (LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})
set (GASPI_VERSION "GPI-2")
set (GASPI_DIR "/tmp_cci/EOD/CODE_OPTIM/d101219/perf_tools/${GASPI_VERSION}")
set (GASPI_INCLUDE "${GASPI_DIR}/include")
set (FMM_DIR "/da/soc/groupes/csc/projet.h4h/d101219/FMM/FMM-lib")
set (FMM_INCLUDE "${FMM_DIR}/include")
set (FMM_HEADERS "${FMM_DIR}/src")

MESSAGE( STATUS "GASPI_DIR:				" ${GASPI_DIR} )
MESSAGE( STATUS "GASPI_INCLUDE:			" ${GASPI_INCLUDE} )
MESSAGE( STATUS "FMM_DIR:				" ${FMM_DIR} )
MESSAGE( STATUS "FMM_INCLUDE:			" ${FMM_INCLUDE} )
MESSAGE( STATUS "FMM_HEADERS:			" ${FMM_HEADERS} )

# sources
file (GLOB_RECURSE FMM_src_list ../src/*)
include_directories(${FMM_INCLUDE} ${FMM_HEADERS})

# NM tools
set (NMTOOLS_LIB 		"/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT/libNMToolkit.so")
set (NMTOOLS_INCLUDE 	"/da/soc/groupes/csc/projet.h4h/d101219/NM_TOOLKIT")

# Boost
set (BOOST_MPI_LIB 		"/tmp_cci/EOD/CODE_OPTIM/d101219/perf_tools/boost/lib/libboost_mpi.so")
set (BOOST_SERIALIZATION_LIB "/tmp_cci/EOD/CODE_OPTIM/d101219/perf_tools/boost/lib/libboost_serialization.so")
set (BOOST_LIB 
		${BOOST_MPI_LIB} 
		${BOOST_SERIALIZATION_LIB}
	)
set (BOOST_INCLUDE 		"/tmp_cci/EOD/CODE_OPTIM/d101219/perf_tools/boost/include")

# set compile FLAGS
set (CMAKE_CXX_FLAGS "-g -std=c++11 -Wall")
set (MPI_COMPILE_FLAGS "-g -std=c++11 -Wall")

# set library
add_library (
    fmm
    SHARED
    ${FMM_src_list}
)

# verbosity
set (CMAKE_VERBOSE_MAKEFILE 1)

# MPI
find_package (MPI)
include_directories(${MPI_INCLUDE_PATH})

# GASPI
include_directories (${GASPI_INCLUDE})

# LINK
target_link_libraries(fmm
	${MPI_LIBRARIES} 
	${GASPI_LIBRARIES}
        ${NMTOOLS_LIB}
        ${BOOST_LIB}
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(fmm
  PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(fmm
  PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
