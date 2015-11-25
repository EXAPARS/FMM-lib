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

#include "Histogram.hpp"

BFHistogram::BFHistogram(const char & B, const int & dim, vec3D * coordinates, 
	const int & first, const int & last, const ui64 & tmpMedian, 
	const int & prefixSize, const int & chunkSize)
	: _nbParticles(last-first+1)
{	
	_axis = dim;
	switch (B)
	{
		case 'E' : constructHistE(coordinates, first, last);  break;
		case 'M' : constructHistM(coordinates, prefixSize, chunkSize, first, last, tmpMedian); break;
		default : cerr << "ERROR! BFHistogram::BFHistogram -> Neither Lower nor Heigher Bytes" << endl; exit(0);
	}	
}

void BFHistogram::constructHistE(vec3D * coordinates, const int & first, const int & last)
{	
	_histo.resize(H12_SIZE); // 2¹²	
	ui64 *db_as_int = NULL;
	ui64 s_and_exp;
	
	for (int i=first; i<=last; ++i)
	{			
		db_as_int = (ui64 *)&(coordinates[i][_axis]);		
		s_and_exp = (*db_as_int >> 52) & 0xFFF;		
		_histo[s_and_exp]++;
	}
}


void BFHistogram::constructHistM(vec3D * coordinates, 
	const int & prefixSize, const int & chunkSize,
	const int & first, const int & last, const ui64 & tmpMedian)
{	
	u_int16_t mask;
	ui64 m;
	
	switch (chunkSize)
	{
		case 4  : _histo.resize(16); 	mask = 0xF;   break;		
		case 16 : _histo.resize(65536); mask = 0xFFFF; break;
		default : cerr << "ERROR! BFHistogram::BFHistogram -> Unexpected chunkSize : "<< chunkSize << endl; exit(0);
	}
	
	ui64 *db_as_int = NULL;
	
	for (int i=first; i<=last; ++i)
	{			
		db_as_int = (ui64 *)&(coordinates[i][_axis]);		
		
		if ( ((*db_as_int ^ tmpMedian) >> (64-prefixSize) )== 0) //same prefix
		{
			m = (*db_as_int >> (64 - prefixSize - chunkSize)) & mask;		
			_histo[m]++;
		}
	}
}

void BFHistogram::display()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	cout << endl << "*** Histogram display" << endl;
	int sum=0;
	
	for (decltype(_histo.size()) i=0; i<_histo.size(); ++i)
	{
		if (_histo[i])
		{
			cout << "[" << i << "] = " << _histo[i] << endl;
			sum += _histo[i];
		}
	}
	cout << std::dec << "Nb Particles : " << sum << endl;
}


//---------------------------- fonctions de classe
void displayGlobal(int * v, int size)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	cout << endl << "*** Global histogram [" <<  rank << "]" << endl;		
	int sum=0;

	for (int i=0; i<size; ++i)
	{
		if (v[i])
		{
			sum += v[i];
			cout << "[" << i << "] = " << v[i] << endl;
		}
	}
	cout << std::dec << "Nb Particles : " << sum << endl;
}
