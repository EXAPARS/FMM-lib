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

cmake_minimum_required (VERSION 2.6)

PROJECT (FMM)
set (CMAKE_CXX_COMPILER g++)

# variables
set (LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

set (GASPI_VERSION "gpi")
set (GASPI_DIR "/home/nath/opt/${GASPI_VERSION}")
set (GASPI_INCLUDE "${GASPI_DIR}/include")

set (FMM_DIR "/home/nath/Documents/Dassault-UVSQ/FMM/FMM-lib")
set (FMM_INCLUDE "${FMM_DIR}/include")
set (FMM_HEADERS "${FMM_DIR}/src")

MESSAGE( STATUS "GASPI_DIR:				" ${GASPI_DIR} )
MESSAGE( STATUS "GASPI_INCLUDE:			" ${GASPI_INCLUDE} )
MESSAGE( STATUS "FMM_DIR:				" ${FMM_DIR} )
MESSAGE( STATUS "FMM_INCLUDE:				" ${FMM_INCLUDE} )
MESSAGE( STATUS "FMM_HEADERS:				" ${FMM_HEADERS} )

# sources
file (GLOB_RECURSE FMM_src_list ../src/*)
include_directories(${FMM_INCLUDE} ${FMM_HEADERS})


# set compile FLAGS
set (CMAKE_CXX_FLAGS "-std=c++11 -Wall")
set (MPI_COMPILE_FLAGS "-std=c++11 -Wall")

# set library
add_library (
	fmm
    SHARED
	${FMM_src_list}
)

# verbosity
set (CMAKE_VERBOSE_MAKEFILE 0)

# MPI
find_package (MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# GASPI
include_directories (${GASPI_INCLUDE})

target_link_libraries(fmm
	${MPI_LIBRARIES} 
	${GASPI_LIBRARIES}
)

if(MPI_COMPILE_FLAGS)
  set_target_properties(fmm
  PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(fmm
  PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
