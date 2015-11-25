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

#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "mpi.h"
#include <vector>
#include <iostream>

#include "types.hpp"
#include "fmm_tools.hpp"
#include "vec3D.hpp"

using namespace std;


/**
 *  Bitfield version
 **/

class BFHistogram{	
private:
	int _nbParticles;
	int _axis;
	vector<int> _histo;
public:
	// constructor
	BFHistogram(const char & B, const int & axis, vec3D * coordinates, 
		const int & first, const int & last, const ui64 & HMedian,
		const int & prefixSize, const int & chunkSize);	
	
	// getter
	vector<int> getHistogram() const { return _histo; }
	int getNbParticles() const { return _nbParticles; }	
	
	// utils for contructor
	void constructHistE(vec3D * coordinates, const int & first, const int & last);		
	void constructHistM(vec3D * coordinates, 
		const int & prefixSize, const int & chunkSize,
		const int & first, const int & last, const ui64 & tmpMedian);

	// display
	void display();
};

#endif

	void displayGlobal(int * v, int size);
