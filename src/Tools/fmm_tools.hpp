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

#ifndef FMM_TOOLS_HPP
#define FMM_TOOLS_HPP

#include "mpi.h"
#include "types.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream> 
#include <cmath>
#include <algorithm>
#include <fstream>
#include <vector>
#include "Complex.hpp"


using namespace std;

void announce_axis(string axis, int rank);
void displayHexa2Dim (string info, ui64 ** tab, int dim1, int dim2);
void displayDiff (string info, int * tab, int size);
void displayMpiMSG (int source, int tag);
void dfs_dump_spectre_octree(string prefix, i64 * nbElemPerNode, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nbNodes, i64 * nodeOwner, i64 nodeID, double * centers);
void dfs_dump_centers(string prefix, i64 * nbSonsPerNode, i64 * firstSonId, i64 nodeID, double * centers);


void bfs_dump_centers_level_by_level(string prefix, i64 * nbElemPerNode, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nbNodes, i64 * nodeOwner, i64 nodeID, double * centers, 
	i64 * endlev, i64 * nbLevels);


void debug(string prefix, string message);

string convert(int a);
template <typename T>
string to_string(T const& value)
{
	stringstream sstr;
	sstr << value;
	return sstr.str();
}


void accumulMSG(string message);

void dumpMSG(string prefix);


/** templates **/
template<typename T>
void displayTab(string info, T * tab, int size, ostream & out=cout)
{
	out << info << endl;
	for (int i=0; i<size; i++)
		out << "[" << i << "] " << tab[i] << "\n";
	out << endl;
}


/** DUMP BUFFERS **/
template<typename T>
void dumpBuffer(int rank, T * buffer, int size, int fileNum, string fileName, string message)
{
	string file = fileName + "_" + to_string((unsigned long long)rank) + "_file_" + to_string((unsigned long long)fileNum) + ".txt";
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
	
	if (buffer == nullptr)
	{
		cerr << "[dumpBuffer] Buffer is not allocated."; exit(-1);
	}
	out << message << endl;
	out << "size : " << size << endl;
	for (int i=0; i<size; i++)
		out << buffer[i] << "\n";
	out <<  endl;
	
	out.close();
}

/** DUMP TREE **/
void dump_tree_init(int rank);
void dump_tree_add_child(int rank, int64_t parent, int64_t child, int nbParticles);
void dump_tree_close_file(int rank);


/** DIFF DATA **/
void loadAndDiffData(const string & file1, const string & file2);

string itoa(int a);


/** RDTSC **/
uint64_t rdtsc(void);

#endif
