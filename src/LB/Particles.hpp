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

#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include "mpi.h"

#include <fstream>
#include <iostream>
#include <iomanip>

#include "Decomposition.hpp"
#include "Histogram.hpp"
#include "Gaspi_communicator.hpp"

using namespace std;

class Particles
{
public:
	static vec3D * _coordinates;
	static double _coeff;
	static double _translate;

protected:
	int _first;
	int _last;
	int _nbParticles;
	
	// bb
	ui32 _edge;
	vec3D _origin;

public:
	// constructors and destructors
	Particles();
	Particles(const Particles & p);
	Particles(const int & nbParticles, const string & file);
	Particles(int nbParticles, string file, Gaspi_communicator * gComm);
	~Particles(){}
	
	// getters
	int getNbParticles() const { return _nbParticles; }
	int getFirstIndex() const { return _first; }
	int getLastIndex() const { return _last; }
	ui32 getEdge() const { return _edge; }
	vec3D getOrigin() const { return _origin; }
	vector<vec3D> getGlobalCoordsV() const;
	vec3D * getGlobalCoords() const { return _coordinates; }
	double * getGlobalCoordsD() const 
	{
		vec3D * coordinates = getGlobalCoords();
		double * coordsD = reinterpret_cast<double*>(coordinates);
		return coordsD;
	}
	double getCoeff() const {return _coeff;}
	double getTranslate() const {return _translate; }

	// setters
	void loadCoordinates(const int & nbParticles, const string & file);		
	void loadCoordinates(const string & file);
	void loadCoordinatesASCII(const string & file);
	void loadCoordinatesASCII(const int & nbParticles, const string & file);
	void loadCoordinatesASCII(const int & nbParticles, const string & file, Gaspi_communicator * gComm);
	void loadCoordinatesASCIIWithoutQuantity(const string & file);
		
	void initScalingParameters();
	void scale();
	void scaleWithParams(double translate, double coeff);
	
	void setAttributes(const int & index, const int & nbParticles, const int & edge, const vec3D & o);
	void setAttrBounds(const int & firstIndex, const int & lastIndex);
	void setNewCoordinates(vec3D * newCoords, const int & nbParticles) 
	{ 	
		_coordinates = newCoords;
		_nbParticles = nbParticles;
	}
	void copyNewCoordinates(vec3D * newCoords, const int & nbParticles)
	{
		if (_coordinates)
			delete _coordinates;
		
		_coordinates = new vec3D[nbParticles];
		
		for (int i=0; i<nbParticles; i++)
		{
			_coordinates[i].x = newCoords[i].x;
			_coordinates[i].y = newCoords[i].y;
			_coordinates[i].z = newCoords[i].z;
		}
		
		_first = 0;
		_last = nbParticles - 1;
		_nbParticles = nbParticles;
	}

	// operators
	Particles &operator=(const Particles & p);	
	bool operator==(const Particles & p);
	bool operator!=(const Particles & p);	

	// divide Octree
	Particles * divideOctree();
	void dispatchOctree(int * quantities, int * destinations);
	void swapOctree(int * qty, int * tabDest);
	
	/** Load Balancing	**/
	
	/// Exact Histograms
	void compSepHistExact(const int & depth, const decompo & decomp,
		const int & nbWorkers, Particles **& p, int *& flatIdxes, int & flatIdxSize);		
	void compSepExpExact(int & sumNbItems, const char & histType, const int & dim, 
		const int & nbSeps, int * nbUnderSep, ui64 ** separators,
		const int & nbWorkers, int *flatIdxes);	
	void compSepMantExact(const int & sumNbItems, const char & histType, const int & dim, 
		const int & prefixSize, const int & chunkSize, const int & nbSeps, int * nbUnderSep, 
		ui64 ** separators, const int & nbWorkers, int *flatIdxes);
		
	/// Approx Histograms
	void compSepHistApprox(const int & depth, const decompo & decomp,
		const int & nbWorkers, Particles **& p, int *& flatIdxes, int & flatIdxSize,
		const ui32 edge, const ui32 height, double ** grid, int nbGridAxis, double * nodeCenters, i64 * nodeOwners, int nbLeaves);
	void compSepMantApprox(const int & sumNbItems, const char & histType, const int & dim, 
		const int & prefixSize, const int & chunkSize, const int & nbSeps, int * nbUnderSep, 
		ui64 ** separators, const int & nbWorkers, int *flatIdxes, ui32 c, ui32 h, double ** grid, int nbGridAxis);		
	void updateBoxOwners(i64 * nodeOwners, int nbLeaves, int dim, double * nodeCenters, ui64 ** separators, int nbParts);

	/// Swaps
	void swap(const int & dim, const int & nbSeps, ui64 ** separators, const int & nbWorkers, 
		int * flatIdxes, int ** SepIdx);
	void swapInterval(const int & first, const int & last, const int & dim, 
		const ui64 * separators, const int nbSeps, int * sepIdx);	

	// Particle Exchange
	void exchangeMPI(const int * flatIdxes);	

	// display
	void display(const int & nb) const;
	void displayInfo() const;
	void displayInfoShort() const;

	// OpenGL tools
	vector<double> compBB();	
};

// Separators
void compSepSE(const int * globalHist, const int & nbSeps, const int & sumNbItems, 
	int * nbUnderSep, ui64 * separators);	
void compSepM(const int * globalHist, const int & nbSeps, const int & sumNbItems,
	int & localSep, const int & numSep, int & nbUnderSep, int chunkSize);

// FlatIdexes
void updateFlatten(int *& flatIdxes, int & flatIdxSize, int nbWorkers, int nbSeps, int ** sepIdx);

// Particles tab
void fillParticles(Particles **& p, int nbWorkers, int nbParts, int * flatIdxes);
		 
// Adjust Separators on Octree grid
void adjustOnGrid(ui64 * separators, int nbWorkers, int nbSeps, int nbBits, ui32 c, ui32 h, int k, double * grid, int nbGridAxis);

// Scale any array of coordinates
void scale(vec3D & coords);
void scaleArray(double * coords, int nbcoords);
void copyAndScaleArray(double * vIN, double * vOUT, int nbCoords);
vec3D scaleBack(vec3D &coords);
double scaleBackDB(double coord);


#endif
