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

#include "Particles.hpp"


// Class variables initialization
vec3D * Particles::_coordinates = NULL;

/*----------------------------------------------------------------------
							Constructors
----------------------------------------------------------------------*/

/**
* Default constructor.
*/
Particles::Particles()
	: _first(0)
	, _last(0)
	, _nbParticles(0)
	, _edge(COORDMAX)
	, _origin(vec3D(0,0,0))	
{}

/**
* Copy constructor.\n
* Constructs an instance of Class Particles by copying another instance of Particles.\n
* No need to modify the Class Variables.\n
* @param particles : an instance of Class Particles.
*/
Particles::Particles(const Particles & particles)
	: _first(particles._first)
	, _last(particles._last)
	, _nbParticles(particles._nbParticles)
	, _edge(particles._edge)
	, _origin(particles._origin)
{}

/**
* Constructor.\n
* Constructs an instance of Class Particles by reading nbParticles from an input file.\n
* Initializes Class Variables table _coordinates which holds all the coordinates,
* and reads the input.
* @param nbParticles : Number of particles
* @param file : File containing the coordinates
*/
Particles::Particles(const int & nbParticles, const string & file)
	: _first(0)
	, _nbParticles(nbParticles)
	, _edge(COORDMAX)
	, _origin(vec3D(0,0,0))
{
	 _last = _first + nbParticles - 1;
	_coordinates = new vec3D[nbParticles];
	loadCoordinatesASCII(file);
}

Particles::Particles(int nbParticles, string file, Gaspi_communicator * gComm)
	: _first(0)
	, _nbParticles(nbParticles)
	, _edge(COORDMAX)
	, _origin(vec3D(0,0,0))
{
	_last = _first + nbParticles - 1;
	_coordinates = gComm->_initCoords;
	loadCoordinatesASCII(file);
}


				
/*----------------------------------------------------------------------
							Getters
----------------------------------------------------------------------*/
/**
* This method is for OpenGL vizualization which expects vectors of doubles.
* @return Returns a vector containing the global array of coordinates
*/
vector<vec3D> Particles::getGlobalCoordsV() const
{
	vector<vec3D> v;
	v.resize(_nbParticles);
	for (int i=0; i<_nbParticles; ++i)
		v[i] = _coordinates[i];		
	return v;
}

/*----------------------------------------------------------------------
							Setters
----------------------------------------------------------------------*/
/**
* Loads the coordinates into the Class Variable _coordinates +by reading from an input file.\n
* @param file : File containing the coordinates.
* \throws Throws an exception if read fails. 
* \warning Class variable must be allocated, and _nbParticles must contain the number of particles to read.
*/
void Particles::loadCoordinates(const string & file)
{	
	if (!_coordinates)
	{
		cerr << "ERROR! : _coordinates was not allocated." << endl;
		exit(0);
	}

	cout << file << endl;
	ifstream in;
	in.exceptions(ifstream::failbit|ifstream::badbit);
	try
	{
		in.open(file, ifstream::in | ifstream::binary);
		in.read(reinterpret_cast<char*> (Particles::_coordinates), 3 * _nbParticles * sizeof(double));		
		in.close();
	}
	catch(ifstream::failure e)
	{
		cerr << "[erreur !] " << e.what() << endl;
		exit(0);
	}
}

/**
* Loads the coordinates into the Class Variable _coordinates +by reading from an --ASCII-- input file.\n
* @param file : File containing the coordinates.
* \throws Throws an exception if read fails. 
* \warning Class variable must be allocated, and _nbParticles must contain the number of particles to read.
*/
void Particles::loadCoordinatesASCII(const string & file)
{	
	if (!_coordinates)
	{
		cerr << "ERROR! : _coordinates was not allocated." << endl;
		exit(0);
	}

	fstream in;
	in.exceptions(ifstream::failbit|ifstream::badbit);
	
	try
	{
		in.open(file, ifstream::in);
		int index = 0;
		while (index < _nbParticles)
		{
			in >> _coordinates[index].x;
			in >> _coordinates[index].y;
			in >> _coordinates[index].z;			
			index++;
		}
		in.close();		
	}
	catch(ifstream::failure e)
	{
		cerr << "[erreur] " << e.what() << endl;
		exit(0);
	}
}

void Particles::scale()
{	
	// compute scaling parameters
	double min, max;
	min = max = 0;	
	for (int i=0; i<_nbParticles; i++)
	{
		for (int j=0; j<3; j++)
		{
			if (_coordinates[i][j] > max)
				max = _coordinates[i][j];
			else if (_coordinates[i][j] < min)
				min = _coordinates[i][j];
		}			
	}
	
	double global_max, global_min;
	MPI_Allreduce(&min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	/** F7X **/
	/** ATTENTION SCALING EN DUR CALCULE PAR RAPPORT A 1 MPI**/	
/**	double min = -11578.99023;
	double max = 21199.99609;
**/
	
	/** DRONERA **/
	//double min = -1436.6006;
	//double max = 1913.9934;	
	
	/** SPHERE 3GHZ **/
	/*double min = -100.0;
	double max = 100.0;*/
	cout << "---------------- global_min " << global_min << endl;
	cout << "---------------- global_max " << global_max << endl;

/*cout << setprecision(8);
cout << min << endl;
cout << max << endl;
*/
	double coeff = COORDMAX / (global_max - global_min) ;
	double translate = global_min * -1.0;
/*
cout << coeff << endl;
cout << translate << endl;
*/	
	// scale the coordinates
	for (int i=0; i<_nbParticles; i++)
		for (int j=0; j<3; j++)
		{
			_coordinates[i][j] = ( _coordinates[i][j] + translate) * coeff;	
		}

	/// ATTENTION VERIF A L'ARRACHE
/**	for (int i=0; i<_nbParticles; i++)
		for (int j=0; j<3; j++)
			if (_coordinates[i][j] < 0)
				cout << i << " " << j << " coordonnée négative : "<< _coordinates[i][j] << endl;
	cout << "Verification des coordonnées terminée, pas de négatif." << endl;
**/	
}



/**
* Loads the coordinates into the Class Variable _coordinates .\n
* Updates the instance attributes, allocates the Class variable _coordinates and reads the coordinates from file.
* @param file : File containing the coordinates.
* @param nbParticles : Number of particles to load.
*/
void Particles::loadCoordinates(const int & nbParticles, const string & file)
{
	_last = _first + nbParticles - 1;
	_nbParticles = nbParticles;
	_coordinates = new vec3D[nbParticles];	
	loadCoordinates(file);
}

void Particles::loadCoordinatesASCII(const int & nbParticles, const string & file, Gaspi_communicator * gComm)
{	
	_last = _first + nbParticles - 1;
	_nbParticles = nbParticles;
	_coordinates = gComm->_initCoords;	
	loadCoordinatesASCII(file);
}


/**
* Loads the coordinates into the Class Variable _coordinates .\n
* Updates the instance attributes, allocates the Class variable _coordinates and reads the coordinates from file.
* @param file : File containing the coordinates.
* @param nbParticles : Number of particles to load.
*/
void Particles::loadCoordinatesASCII(const int & nbParticles, const string & file)
{
	_last = _first + nbParticles - 1;	
	_nbParticles = nbParticles;
	_coordinates = new vec3D[nbParticles];	
	loadCoordinatesASCII(file);
}

/**
* Loads the coordinates into the Class Variable _coordinates - Without knowing the number of coordinates to read .\n
* Updates the instance attributes, allocates the Class variable _coordinates and reads the coordinates from file.
* @param file : File containing the coordinates.
*/
void Particles::loadCoordinatesASCIIWithoutQuantity(const string & file)
{
	fstream in;
	in.exceptions(ifstream::failbit|ifstream::badbit);
	try
	{
		in.open(file, ifstream::in);
		
		// get number of nodes
		int nbParticles = std::count(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>(), '\n');			
		_last = _first + nbParticles - 1;	
		_nbParticles = nbParticles;
		_coordinates = new vec3D[nbParticles];
		
		// load into _coordinates
		in.seekg(0, ios::beg);
		int index = 0;
		while (index < _nbParticles)
		{
			in >> _coordinates[index].x;
			in >> _coordinates[index].y;
			in >> _coordinates[index].z;			
			index++;
		}
		in.close();		
	}
	catch(ifstream::failure e)
	{
		cerr << "[erreur] " << e.what() << endl;
		exit(0);
	}
}



/**
* Sets the attributes _nbParticles, _first, _edge and _origin.
*/
void Particles::setAttributes(const int & index, const int & nbParticles,  const int & edge, const vec3D & o)
{
	_nbParticles = nbParticles;
	if (nbParticles)
	{
		_first = index;	
		_last = index + nbParticles - 1;
	}
	else
	{
		_first = _last = -1;	
	}
	_edge = edge;
	_origin = o;
}

void Particles::setAttrBounds(const int & firstIndex, const int & lastIndex)
{
	_first = firstIndex;
	_last = lastIndex;
	_nbParticles = lastIndex - firstIndex + 1; 
}

/*----------------------------------------------------------------------
							Operators
----------------------------------------------------------------------*/
/**
* Affectation operator.
*/
Particles& Particles::operator=(const Particles & p)
{
	_nbParticles = p.getNbParticles();
	_first = p.getFirstIndex();
	_last = p.getLastIndex();
	_origin = p.getOrigin();
	_edge = p.getEdge();
	return *this;
}

/**
* Equality test operator.
* @return : Boolean true | false.
*/
bool Particles::operator==(const Particles & p)
{
	if ( _nbParticles != p.getNbParticles() ) return false;
	if ( _first != p.getFirstIndex() ) return false;
	if ( _origin != p.getOrigin() ) return false;
	if ( _edge != p.getEdge() ) return false;
	return true;
}

/**
* Inequality test operator.
* @return : Boolean true | false.
*/
bool Particles::operator!=(const Particles & p)
{
	return !(this->operator ==(p));
}


/*----------------------------------------------------------------------
							DIVIDE OCTREE
----------------------------------------------------------------------*/

/**
* Method that divides a node into 8 smaller nodes. \n
* dispatchOctreees the particles. \n
* Swaps the global coordinates. \n
* @return Returns an array of new Particles.
*/
Particles * Particles::divideOctree()
{
	/**
	* qty : Array of integers, holds the number of particles per octree box.
	* boxDest : Array of _nbParticles integers : holds each particle box destination.
	* sums : Array containing for each boxIndex the number of particles from box 0 to the box Index.
	**/
	
	int qty[8] = {0};
	int * destinations = new int[_nbParticles]();

	// Fill the qty and boxDest destinations for the dispatch
	dispatchOctree(qty, destinations);
	
	// intermediate sums computations
	int sums[8] = {0};	
	sums[0] = qty[0];
	sums[1] = sums[0] + qty[1];
	sums[2] = sums[1] + qty[2];
	sums[3] = sums[2] + qty[3];
	sums[4] = sums[3] + qty[4];
	sums[5] = sums[4] + qty[5];
	sums[6] = sums[5] + qty[6];
	sums[7] = sums[6] + qty[7];
	
	// swap the coordinates to fit to the octree
	swapOctree(sums, destinations);
	delete [] destinations;

	// return an array of particles	
	int halfEdge = _edge/2;
	double xMIN = _origin.x;
	double xMID = _origin.x + halfEdge;
	double yMIN = _origin.y;
	double yMID = _origin.y + halfEdge;
	double zMIN = _origin.z;
	double zMID = _origin.z + halfEdge;
	
	// create an array of 8 children, fill it
	Particles * p = new Particles[8];
	
	/// MORTON ORDER	
	p[0].setAttributes( _first,   	      qty[0], halfEdge, vec3D(xMIN, yMID, zMIN) );	//2
	p[1].setAttributes( _first + sums[0], qty[1], halfEdge, vec3D(xMID, yMID, zMIN) );	//3
	p[2].setAttributes( _first + sums[1], qty[2], halfEdge, vec3D(xMIN, yMIN, zMIN) );	//1
	p[3].setAttributes( _first + sums[2], qty[3], halfEdge, vec3D(xMID, yMIN, zMIN) );	//0
	p[4].setAttributes( _first + sums[3], qty[4], halfEdge, vec3D(xMIN, yMID, zMID) );	//6
	p[5].setAttributes( _first + sums[4], qty[5], halfEdge, vec3D(xMID, yMID, zMID) );	//7
	p[6].setAttributes( _first + sums[5], qty[6], halfEdge, vec3D(xMIN, yMIN, zMID) );	//5
	p[7].setAttributes( _first + sums[6], qty[7], halfEdge, vec3D(xMID, yMIN, zMID) );	//4


	return p;
}

/**
* Method that dispatches the node's particles into 8 children nodes.\n
* Fills the \b quantities and \b destinations arrays.\n
* @param quantities : Array of 8 integers, holding the number of particles per children node box.
* @param destinations : Array of nbParticles, holding the id of the destination octree box.
*/
void Particles::dispatchOctree(int * quantities, int * destinations)
{
	double x, y, z, xMIN, xMID, xMAX, yMIN, yMID, yMAX, zMIN, zMID, zMAX;
	ui32 edge, halfEdge;
	
	// Test all particles
	for (int Index=_first ; Index < (_first + _nbParticles) ; ++Index)
	{
		// index between 0 and nbParticles
		int i = Index - _first;
		
		// particle coordinates
		x = _coordinates[Index].x;
		y = _coordinates[Index].y;
		z = _coordinates[Index].z;
		
		// edge, half-edge
		edge = _edge;
		halfEdge = _edge/2;
		
		// min, mid and max coordinates
		xMIN = _origin.x;
		xMID = _origin.x + halfEdge;
		xMAX = _origin.x + edge;
		yMIN = _origin.y;
		yMID = _origin.y + halfEdge;
		yMAX = _origin.y + edge;
		zMIN = _origin.z;
		zMID = _origin.z + halfEdge;
		zMAX = _origin.z + edge;
		
		// Identifies the destination box by testing the particle coordinates
		// increases quantity tab and updates destination tab
		if 		( x>=xMIN   &&  x<xMID  &&  y>=yMIN   &&  y<yMID  &&  z>=zMIN   &&  z<zMID ) 
		{ 
			quantities[2]++; 	// was 0
			destinations[i] = 2;
		}
		else if ( x>=xMID  &&  x<=xMAX  &&  y>=yMIN   &&  y<yMID  &&  z>=zMIN   &&  z<zMID ) 
		{
			quantities[3]++;	// was 1
			destinations[i] = 3;
		}			
		else if ( x>=xMID  	&&  x<=xMAX  &&  y>=yMID	&&  y<=yMAX  &&  z>=zMIN   &&  z<zMID ) 
		{
			quantities[1]++;	// was 2
			destinations[i] = 1;
		}
		else if ( x>=xMIN   &&  x<xMID  &&  y>=yMID  &&  y<=yMAX  &&  z>=zMIN   &&  z<zMID ) 
		{
			quantities[0]++;	// was 3
			destinations[i] = 0;
		}
		else if ( x>=xMIN   &&  x<xMID  &&  y>=yMIN   &&  y<yMID  &&  z>=zMID  &&  z<=zMAX )
		{
			quantities[6]++;	// was 4
			destinations[i] = 6;
		}
		else if ( x>=xMID  &&  x<=xMAX  &&  y>=yMIN   &&  y<yMID  &&  z>=zMID  &&  z<=zMAX )
		{
			quantities[7]++;	// was 5
			destinations[i] = 7;
		}
		else if ( x>=xMID  &&  x<=xMAX  &&  y>=yMID  &&  y<=yMAX  &&  z>=zMID  &&  z<=zMAX )
		{
			quantities[5]++;	// was 6
			destinations[i] = 5;
		}
		else if ( x>=xMIN   &&  x<xMID  &&  y>=yMID  &&  y<=yMAX  &&  z>=zMID  &&  z<=zMAX )
		{
			quantities[4]++;	// was 7
			destinations[i] = 4;
		}
		else 
		{
			cerr << "ERROR! dispatchOctree : Unexpected coordinates." << endl;
			cerr << "index = " << i << endl;
			cerr << "x = "<< x << ", y = " << y << ", z = " << z << endl;
			cerr << "xMIN = "<< xMIN << ", xMax =" << xMAX << " | x>=xMIN : " << (x>=xMIN) << " | x<=xMax : " << (x<=xMAX) << endl;
			cerr << "yMIN = "<< yMIN << ", yMax =" << yMAX << " | y>=yMIN : " << (y>=yMIN) << " | y<=yMax : " << (y<=yMAX) << endl;			
			cerr << "zMIN = "<< zMIN << ", zMax =" << zMAX << " | z>=xMIN : " << (z>=zMIN) << " | z<=zMax : " << (z<=zMAX) << endl;			
			exit(0);
		}
	}
}

/**
* Method that swaps the particle coordinates.\n
* Traverses a table of tags. \n
* Terminates when all cells are tagged. \n
* @param sums : Array containing for each boxIndex the number of particles from 0 to the box Index.
* @param boxDest : Array of nbParticles, holding the id of the destination octree box.
*/
void Particles::swapOctree(int * sums, int * boxDest)
{
	/**
	* tag : Array of nbParticles integers. Values are 0/1/-1 for todo/done/startPoint. 
	* ctr : Array of 8 counters. 
	* round : Boolean, holds the round status. 
	* nextOne : Next coordinates to be written into th global array of coordinates. 
	* tmp : Temporary variable to hold coordinates. 
	* boxIndex : Index of the Octree destination box. 
	* globalIndex : Index of the coordinates in the global array of coordinates. 
	**/ 
	
	int * tag = new int[_nbParticles]();
	
	int ctr[8] = {0};
	bool round = false;
	
	vec3D nextOne, tmp;
	int boxIndex;
	int globalIndex;
	
	
	/**
	Traverses the table of tags.
	Terminates when all tag cells are tagged with 1.
	**/
	int i = 0;	
	while (i < _nbParticles)
	{
		// Go to the first element to treat
		while ( tag[i] == 1  &&  i < _nbParticles)
		{
			++i;
		}
		// Detect swap termination
		if ( i == _nbParticles )
		{
			break;		
		}
		
		// Update tag, round and nextOne
		tag[i] = -1;
		round = true;
		nextOne = Particles::_coordinates[_first + i];
		boxIndex = boxDest[i];		
		
		// Go through a round
		while (round)
		{
			// Compute globalIndex
			switch (boxIndex)
			{
				case 0 : globalIndex = ctr[0]; break;
				case 1 : globalIndex = sums[0] + ctr[1]; break;
				case 2 : globalIndex = sums[1] + ctr[2]; break;
				case 3 : globalIndex = sums[2] + ctr[3]; break;				
				case 4 : globalIndex = sums[3] + ctr[4]; break;
				case 5 : globalIndex = sums[4] + ctr[5]; break;
				case 6 : globalIndex = sums[5] + ctr[6]; break;
				case 7 : globalIndex = sums[6] + ctr[7]; break;
				default : cerr << "ERROR! swapOctree : Unknown box index : " << boxIndex << endl; exit(0);			
			}
			
			// Detect round termination
			if ( tag[globalIndex] == -1 )
				round = false;
			
			// Add the offset
			globalIndex = globalIndex + _first;
			
			// Swap the particles
			tmp = Particles::_coordinates[globalIndex];
			Particles::_coordinates[globalIndex] = nextOne;
			
			
			// substract the offset
			globalIndex = globalIndex - _first;
			
			// Update tag and counters
			tag[globalIndex] = 1;
			ctr[boxIndex]++;
			
			// Prepare next iteration
			if (round)
			{
				nextOne = tmp;			
				boxIndex = boxDest[globalIndex];
			}				
		}		
		// Increase i
		++i;
	}
	delete[] tag;
}

/*----------------------------------------------------------------------
				EXACT HISTOGRAMS
----------------------------------------------------------------------*/

/**
* This method computes the exact separators thanks to histograms.
* The particles are swapped to respect the new separators.
* The flatIdxes array is updated with the new separators.
* P is filled with th new particles in order to update the level of nodes. 
* @param depth : Current depth.
* @param decomp : Prime numbers decomposition list.
* @param level : Complete level of nodes on the current depth
* @param p : Array of particles, which will be used to update the children nodes
* @param flatIdxes is an array of separators, which will be updated with the last recursion.
* @param flatIdxSize is the size of the flatIdxes array of separators
**/
void Particles::compSepHistExact(const int & depth, const decompo & decomp,
	const int & nbWorkers, Particles **& p, int *& flatIdxes, int & flatIdxSize)
{	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// get parts, seps, dim and workers
	int nbParts = decomp._list[depth]; 	 				 
	int nbSeps = nbParts-1;
 	int dim = depth % 3;

	// initialize arrays : counters, separators values and indexes						 		
	int * nbUnderSep = new int [nbSeps](); 

	ui64 ** separators = new ui64* [nbWorkers];	
	for(int i=0; i<nbWorkers; i++)
		separators[i] = new ui64[nbSeps]();	 	

	int ** SepIdx = new int* [nbWorkers];	
	for(int i=0; i<nbWorkers; i++)
		SepIdx[i] = new int[nbSeps]();
	
	// sign + exponent 
 	int sumNbItems = 0;	
	compSepExpExact(sumNbItems, 'E', dim, nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes);	
			
	// 16 , 16, 16, 4 bits of mantissa
	compSepMantExact(sumNbItems, 'M', dim, 12, 16,  nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes);		
	compSepMantExact(sumNbItems, 'M', dim, 28, 16,  nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes);
	compSepMantExact(sumNbItems, 'M', dim, 44, 16,  nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes);
	compSepMantExact(sumNbItems, 'M', dim, 60,  4,  nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes);
	
	// swap		
	swap(dim, nbSeps, separators, nbWorkers, flatIdxes, SepIdx);				

 	// Update flatten Indexes and size
	updateFlatten(flatIdxes, flatIdxSize, nbWorkers, nbSeps, SepIdx);
 
	// update the array of particles
	fillParticles(p, nbWorkers, nbParts, flatIdxes);

	// dealloc
	delete [] nbUnderSep;
	for(int i=0; i<nbWorkers; i++)
	{
		delete [] separators[i];
		delete [] SepIdx[i];
	}
	delete [] separators;
	delete [] SepIdx;
}

/**
* This method computes the 12 first bits of the separators : sign + exponent.
* The Separators array is updated and broadcasted. 
* @param sumNbItems : Sum of particles on all processes in the current part to process.
* @param histType : Exponent or Mantissa, parameter for the histogram computation.
* @param dim : 1, 2 or 3 for the coordinate dimension selection.
* @param nbSeps : Number of separators to update.
* @param nbUnderSep : Number of Particles under the separators.
* @param separators : Array of separators.
* @param nbWorkers : Number of processes taking part to the current computation.
* @param flatIdxes : Array of separators, is updated here.
**/
void Particles::compSepExpExact(int & sumNbItems, const char & histType, const int & dim, 
	const int & nbSeps, int * nbUnderSep, ui64 ** separators,
	const int & nbWorkers, int *flatIdxes)
{		
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// alloc
	int *globalHist = NULL;
	if (rank < nbWorkers)
		globalHist = new int [H12_SIZE]();

	// Histograms and reductions
	int index=0;
	for (int rank=0; rank < nbWorkers; rank++)
	{
		// Everybody computes the histogram		
		BFHistogram H(histType, dim, getGlobalCoords(), flatIdxes[index]+1, flatIdxes[index+1], 0, 0, 0);
		
		// Reduction on responsible worker 		
		MPI_Reduce(H.getHistogram().data(), globalHist, H12_SIZE, MPI_INT, MPI_SUM, rank, MPI_COMM_WORLD);
		index++;
	}

	// compute SumNbItems
	if (rank < nbWorkers)
	{
		int cpt = 0;
		for (int i=0; i<H12_SIZE; i++)
			cpt += globalHist[i];		
		sumNbItems = cpt ;
	}
 	
 	// separators computation
	if (rank < nbWorkers){
		compSepSE(globalHist, nbSeps, sumNbItems, nbUnderSep, separators[rank]);
	}

	// Broadcast
	for(int rank=0; rank < nbWorkers; rank++)
		MPI_Bcast(separators[rank], nbSeps, MPI_UNSIGNED_LONG, rank, MPI_COMM_WORLD);
	
	// dealloc
	if (rank < nbWorkers)
		delete [] globalHist;
}

/**
* This method computes chunkSize bits of mantissa separators.
* The Separators array is updated and broadcasted. 
* @param sumNbItems : Sum of particles on all processes in the current part to process.
* @param histType : Exponent or Mantissa, parameter for the histogram computation.
* @param dim : 1, 2 or 3 for the coordinate dimension selection.
* @param prefixSize : Already computed bits.
* @param chunkSize : Number of bits to compute.
* @param nbSeps : Number of separators to update.
* @param nbUnderSep : Number of Particles under the separators.
* @param separators : Array of separators.
* @param nbWorkers : Number of processes taking part to the current computation.
* @param flatIdxes : Array of separators, is updated here.
**/
void Particles::compSepMantExact(const int & sumNbItems, const char & histType, const int & dim, 
		const int & prefixSize, const int & chunkSize, const int & nbSeps, int * nbUnderSep, 
		ui64 ** separators, const int & nbWorkers, int *flatIdxes)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
			
	// Buffer allocation
	int *globalHist = NULL;
	int histSize = (1 << chunkSize); //pow(2, chunkSize)
	int localSep = 0;
	if (rank < nbWorkers)
		globalHist = new int [histSize]();		

	// Histograms and reductions
	int index=0;
	for (int i=0; i<nbWorkers; i++)
	{
		for (int k=0; k<nbSeps; k++)
		{
			// Everybody computes the histogram
			BFHistogram H(histType, dim, getGlobalCoords(), flatIdxes[index]+1, flatIdxes[index+1], 
				separators[i][k], prefixSize, chunkSize);
			
			// Reduction on responsible worker 
			MPI_Reduce(H.getHistogram().data(), globalHist, histSize, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
			
			// Separator update by responsible worker
			if (rank == i)
			{
				compSepM(globalHist, nbSeps, sumNbItems, localSep, k, nbUnderSep[k], chunkSize);				
				separators[rank][k] += (((ui64)localSep) << (64 - prefixSize - chunkSize));
			}
		}
		index++;
	}
	
	// Broadcast
	for(int rank=0; rank<nbWorkers; rank++)
		MPI_Bcast(separators[rank], nbSeps, MPI_UNSIGNED_LONG, rank, MPI_COMM_WORLD);
	
	// Deallocation
	if (rank < nbWorkers)
		delete [] globalHist;
}


/*----------------------------------------------------------------------
				APPROXIMATED HISTOGRAMS ON OCTREE GRID
----------------------------------------------------------------------*/

/**
* This method computes approximates separators. 
* The separators are modified in order to stick to the octree grid.
* The tree height is computed in order to know the separators coordinates.
* @param depth : Current depth.
* @param decomp : Prime numbers decomposition list.
* @param level : Complete level of nodes on the current depth
* @param p : Array of particles, which will be used to update the children nodes
* @param flatIdxes is an array of separators, which will be updated with the last recursion.
* @param flatIdxSize is the size of the flatIdxes array of separators

**/
void Particles::compSepHistApprox(const int & depth, const decompo & decomp,
	const int & nbWorkers, Particles **& p, int *& flatIdxes, int & flatIdxSize,
	const ui32 edge, const ui32 height)
{	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// get parts, seps, dim and workers
	int nbParts = decomp._list[depth]; 	 				 
	int nbSeps = nbParts-1;
 	int dim = depth % 3;

	// initialize arrays : counters, separators values and indexes						 		
	int * nbUnderSep = new int [nbSeps](); 
	
	ui64 ** separators = new ui64* [nbWorkers];	
	for(int i=0; i<nbWorkers; i++)
		separators[i] = new ui64[nbSeps+1](); // +1 is an ui64 to tag finished values.

	int ** SepIdx = new int* [nbWorkers];	
	for(int i=0; i<nbWorkers; i++)
		SepIdx[i] = new int[nbSeps]();	
	
	// sign + exponent 
 	int sumNbItems = 0;	
	compSepExpExact(sumNbItems, 'E', dim, nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes);	
						
	// 16 , 16, 16, 4 bits of mantissa
	compSepMantApprox(sumNbItems, 'M', dim, 12, 16, nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes, edge, height);	
	compSepMantApprox(sumNbItems, 'M', dim, 28, 16, nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes, edge, height);
	compSepMantApprox(sumNbItems, 'M', dim, 44, 16, nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes, edge, height);	
	compSepMantApprox(sumNbItems, 'M', dim, 60,  4, nbSeps, nbUnderSep, separators, nbWorkers, flatIdxes, edge, height);
	
	// swap		
	swap(dim, nbSeps, separators, nbWorkers, flatIdxes, SepIdx);				

 	// Update flatten Indexes and size
	updateFlatten(flatIdxes, flatIdxSize, nbWorkers, nbSeps, SepIdx);

	// fill the array of particles
	fillParticles(p, nbWorkers, nbParts, flatIdxes);

	// dealloc
	delete [] nbUnderSep;
	for(int i=0; i<nbWorkers; i++)
	{
		delete [] separators[i];
		delete [] SepIdx[i];
	}
	delete [] separators;
	delete [] SepIdx;
}



/**
* This method computes chunkSize bits of mantissa separators and tests if an octree separator can be identified.
* If yes : the separator is updated and tagged. If not, the computation will continue with the next bits.
* The Separators array is updated and broadcasted. 
* @param sumNbItems : Sum of particles on all processes in the current part to process.
* @param histType : Exponent or Mantissa, parameter for the histogram computation.
* @param dim : 1, 2 or 3 for the coordinate dimension selection.
* @param prefixSize : Already computed bits.
* @param chunkSize : Number of bits to compute.
* @param nbSeps : Number of separators to update.
* @param nbUnderSep : Number of Particles under the separators.
* @param separators : Array of separators.
* @param nbWorkers : Number of processes taking part to the current computation.
* @param flatIdxes : Array of separators, is updated here.
* @param h : height of the octree.
**/
void Particles::compSepMantApprox(const int & sumNbItems, const char & histType, const int & dim, 
		const int & prefixSize, const int & chunkSize, const int & nbSeps, int * nbUnderSep, 
		ui64 ** separators, const int & nbWorkers, int *flatIdxes, ui32 c, ui32 h)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
			
	// buffer allocation
	int *globalHist = NULL;
	int histSize = (1 << chunkSize); //pow(2, chunkSize)	
	int localSep = 0;
	if (rank < nbWorkers)
		globalHist = new int [histSize]();		


	// Histograms and reductions
	int index=0;
	for (int i=0; i<nbWorkers; i++)
	{
		for (int k=0; k<nbSeps; k++)
		{
			if( !((1<<k) & separators[i][nbSeps]) ) // si ce séparateur n'est pas déjà traité
			{				
				// Everybody computes the histogram
				BFHistogram H(histType, dim, getGlobalCoords(), flatIdxes[index]+1, flatIdxes[index+1], 
					separators[i][k], prefixSize, chunkSize);
				
				// Reduction on responsible worker 
				MPI_Reduce(H.getHistogram().data(), globalHist, histSize, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
				
				// Separator update by responsible worker
				if (rank == i)
				{
					compSepM(globalHist, nbSeps, sumNbItems, localSep, k, nbUnderSep[k], chunkSize);				
					separators[rank][k] += (((ui64)localSep) << (64 - prefixSize - chunkSize));
					
					// Adjust on the Octree Grid if it's possible
					adjustOnGrid(separators[rank], nbWorkers, nbSeps, (prefixSize + chunkSize), c, h, k);
				}
			}
		}
		index++;
	}
	
	// Broadcast
	for(int rank=0; rank<nbWorkers; rank++)
		MPI_Bcast(separators[rank], nbSeps+1, MPI_UNSIGNED_LONG, rank, MPI_COMM_WORLD);
	
	// dealloc
	if (rank < nbWorkers)
		delete [] globalHist;
}


/*----------------------------------------------------------------------
				SWAP METHODS
----------------------------------------------------------------------*/

/**
* This function swaps an interval of particles around an array of separators
* It takes the separators values as input, and updates the particle coordinates array
* and the separator index array.
* @param first : index of the first particle of the selection
* @param last : index of the last particle of the selection
* @param dim : 1, 2 or 3 for the coordinate dimension selection.
* @param separators : Array of separators.
* @param nbSeps : Number of separators to update.
* @param SepIdx : Array of separators idexes, is updated here.
**/
void Particles::swapInterval(const int & first, const int & last, const int & dim, 
	const ui64 * separators, const int nbSeps, int * sepIdx)
{	
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	// Convert separators array from ui64 to double
	double * seps = new double[nbSeps];
	for(int i=0; i<nbSeps; i++)
		seps[i] = *((double *)&(separators[i]));	
	
	// Counters per part
	int nbCounters = nbSeps+1;
	int * counters = new int[nbCounters]();
	
	// Initialize the counters
	for (int i=first; i<=last; i++)
		for (int k=0; k<nbSeps; k++)
			if(_coordinates[i][dim] <= seps[k])
			{
				counters[k]++;
				break;			// next i
			}	

	// Update separator indexes
	sepIdx[0] = first + counters[0] - 1;
	for (int i=1; i<nbSeps; i++)
		sepIdx[i] = sepIdx[i-1] + counters[i];

	// Reset counters
	for (int i=0; i<nbCounters; i++)
		counters[i] = 0;
	
	// Swap variables
	int index = first;	
	int mark = first;
	int targetInterval = -1;
	int targetIndex = 1;
	vec3D aux;	
	bool found;
	
	// For each separator
	for(int i=0; i<nbSeps; i++)
	{
		// Until the separator is reached
		while(index <= sepIdx[i])
		{	
			// if the particle is well placed, increment the index
			if(_coordinates[index][dim] <= seps[i])
				index++;
			else // -> swap
			{				
				// store & mark
				aux = _coordinates[index];
				mark = index;			

				// compute target interval
				found = false;					
				for (int k=0; k<nbSeps; k++)
					if(_coordinates[index][dim] <= seps[k])
					{
						targetInterval = k;
						found = true;
						break;
					}

				// If the target interval has not been identified, it is the last one
				if (!found)
					targetInterval = nbCounters-1;

				// Compute target index
				targetIndex = sepIdx[targetInterval-1] + counters[targetInterval] + 1;
				counters[targetInterval]++;
				
				// swap
				_coordinates[mark] = _coordinates[targetIndex];
				_coordinates[targetIndex] = aux;
			}
		}
	}
		
	// Dealloc
	delete [] seps; 
	delete [] counters;	
}	

/**
* This function calls the swapInterval function on each interval of the flatIdxes array.
* @param dim : 1, 2 or 3 for the coordinate dimension selection.
* @param nbSeps : Number of separators to update.
* @param separators : Array of separators.
* @param nbWorkers : Number of processes taking part to the current computation.
* @param flatIdxes : Array of separators, is updated here.
* @param SepIdx : Array of separator indexes.
**/
void Particles::swap(const int & dim, const int & nbSeps, ui64 ** separators, const int & nbWorkers, 
	int * flatIdxes, int ** SepIdx)
{	
	int index=0;
	for (int rank=0; rank<nbWorkers; rank++)
	{
		swapInterval(flatIdxes[index]+1, flatIdxes[index+1], dim, separators[rank], nbSeps, SepIdx[rank]); 
		index++;			
	}
}


/*----------------------------------------------------------------------
							Particles Exchange
----------------------------------------------------------------------*/

/**
* This function handles the global MPI particles exchange.
* I also updates the particles attributes.
* @param flatIdxes : Array of separators, is updated here.
**/
void Particles::exchangeMPI(const int * flatIdxes)
{
	cout << "--- [ MPI Exchange ";

	/*cout << flatIdxes[0] << " " <<
			flatIdxes[1] << " " <<
			flatIdxes[2] << " " <<
			flatIdxes[3] << " " <<
			flatIdxes[4] << "\n";
*/
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// sendCount, 3 coordonnées par particule
	int * sendCount = new int[size];		
	for (int i=1; i<(size+1); i++)
		sendCount[i-1] = (flatIdxes[i] - flatIdxes[i-1] ) * 3;

/*	cout << rank << endl;
	cout << sendCount[0] << " " <<
			sendCount[1] << " " <<
			sendCount[2] << " " <<
			sendCount[3] << endl;
*/
	// sendOffsets
	int * sOffset = new int[size];	
	for (int i=0; i<size; i++)
		sOffset[i] = (flatIdxes[i] + 1) * 3;
	
	/*cout << "sendOffsets : " << endl;
	cout << sOffset[0] << " " <<
			sOffset[1] << " " <<
			sOffset[2] << " " <<
			sOffset[3] << endl;	
*/

	for (int i=0; i<size; i++)
	{
		cout << rank << " sends to " << i  << "\t" << sendCount[i] << " coordinates, starting at offset : " << sOffset[i] << endl;
	}
	// recvCount 
	int * recvCount = new int[size];	
	MPI_Alltoall(sendCount, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD);	

	int nbRecvCoords = 0;
	for (int i=0; i<size; i++)
		nbRecvCoords += recvCount[i];
	
	// recvOffsets
	int * rOffset = new int[size];
	rOffset[0] = 0;
	for(int i=1; i<size; i++)
		rOffset[i] = rOffset[i-1] + recvCount[i-1];
			
	// AlltoAllv	
	double * recvbuf = new double[nbRecvCoords]();	
		
	MPI_Alltoallv(reinterpret_cast<double*>(_coordinates), sendCount, sOffset, MPI_DOUBLE,
				  recvbuf, 		recvCount, rOffset, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// Update _nbParticles
	if((nbRecvCoords/3) != _nbParticles)
	{
		// resize _coordinates array and update attributes
		delete [] _coordinates;
		_nbParticles = nbRecvCoords/3;
		_coordinates = new vec3D[_nbParticles];
		_last = _nbParticles -1;
	}

/// TODO : voir s'il est possible d'éviter la recopie par un swap de pointeurs !	
	// Copy recvbuf into _coordinates
	for (int i=0; i<_nbParticles; i++)
		for (int j=0; j<3; j++)
			_coordinates[i][j] = recvbuf[(i*3) + j];	

	delete [] sendCount;
	delete [] sOffset;
	delete [] recvCount;
	delete [] rOffset;
	delete [] recvbuf;
	
	cout << "] --- >> OK" << endl;

}

/*----------------------------------------------------------------------
							Display
----------------------------------------------------------------------*/

/**
* This method displays the particles coordinates.
* @param nb : number of particles to display.
*/ 
void Particles::display(const int & nb) const
{
	if (nb > _nbParticles)
	{
		cerr << "ERROR! Particles::display() => incorrect parameter nb";		
		exit(0);
	}
	for (int i=0; i<nb; i++)
	{
		printf("%f  %f  %f\n", 
			_coordinates[_first + i].x,
			_coordinates[_first + i].y,
			_coordinates[_first + i].z);
	}
}

/**
* This method displays the particles attributes : _origin, _edge, _nbParticles and _first.
*/ 
void Particles::displayInfo() const
{
	cout << "--- Particle Display Info ---" << endl;
	cout << "nbParticles : " << _nbParticles << endl;	
	cout << "FirstIndex : " << _first << endl;
	cout << "LastIndex : " << _last << endl;
}

void Particles::displayInfoShort() const
{
	cout << "| " << _nbParticles << " | \t["<< _first << " - " << _last << "]" << endl;	
}
/*----------------------------------------------------------------------
							OpenGL tools
----------------------------------------------------------------------*/
/**
* This methods computes the bounding box coordinates. \n
* It is useful for the OpenGL vizualization. \n
* @return : Returns a vector of 24 doubles, representing the 8 3-dimensional coordinates.
*/
vector<double> Particles::compBB()
{
	double xMIN = _origin.x;
	double xMAX = _origin.x + _edge;
	double yMIN = _origin.y;
	double yMAX = _origin.y + _edge;
	double zMIN = _origin.z;
	double zMAX = _origin.z + _edge;
	
	vector<double> bb;
	bb.resize(24); // 8 x 3D coordinates
	
	bb[0]  = xMIN;	bb[1]  = yMIN; 	bb[2]  = zMIN;
	bb[3]  = xMAX;	bb[4]  = yMIN; 	bb[5]  = zMIN;
	bb[6]  = xMAX;	bb[7]  = yMAX; 	bb[8]  = zMIN;
	bb[9]  = xMIN;	bb[10] = yMAX; 	bb[11] = zMIN;
	bb[12] = xMIN;	bb[13] = yMIN; 	bb[14] = zMAX;
	bb[15] = xMAX;	bb[16] = yMIN; 	bb[17] = zMAX;
	bb[18] = xMAX;	bb[19] = yMAX; 	bb[20] = zMAX;
	bb[21] = xMIN;	bb[22] = yMAX; 	bb[23] = zMAX;
	
	return bb;
}

/*----------------------------------------------------------------------
							Out of Class functions
----------------------------------------------------------------------*/

/**
* This function computes the 12 first bits of each separator in the separators array.
* It takes a global histogram as input, and updates the separators array.
* @param globalHist : Global histogram.
* @param nbSeps : Number of separators to compute.
* @param sumNbItems : Global number of particles in the global histogram.
* @param nbUnderSep : Array of number of particles under each separator.
* @param separators : Already of separators
**/
void compSepSE(const int * globalHist, const int & nbSeps, const int & sumNbItems, 
	int * nbUnderSep, ui64 * separators)
{
	// Counters
	int * tabCpt = new int[nbSeps];			// Nb of items from beginning to Sep
	int * tabMax = new int[nbSeps];			// Max nb of items from beginning to Sep
	int * tabDone = new int[nbSeps]();		// Flag indicating the separator's state		
	
	// initialize tabMax and tabCpt
	int nbParts = nbSeps+1;	
	for (int k=0; k<nbSeps; k++)
	{
		tabMax[k] = (sumNbItems/nbParts)*(k+1);
		tabCpt[k] = nbUnderSep[k];
	}
	
	// Find Separators
	// Histogram's positive's values only	
	for(int i=0; i<=2046; i++)
		
		// For each separator		
		for (int k=0; k<nbSeps; k++)
		
			// If not already treated
			if(!tabDone[k])
			{
				// Update number of particles
				tabCpt[k] += globalHist[i];
				
				// If the number of particles exceeds the max
				if (tabCpt[k] >= tabMax[k])
				{	
					// Update the separators tab
					separators[k] = ((ui64)i << 52);
					
					// Update the nbUnderSep counter by substracting the last quantity
					// (next step will refine the histogram)
					nbUnderSep[k] = tabCpt[k]-globalHist[i];
					
					// and mark the separator as treated
					tabDone[k]=1;
				}
			}

	// Dealloc
	delete [] tabCpt;
	delete [] tabMax;
	delete [] tabDone;
}	

/**
* This function computes the chunkSize next bits of one separator.
* It takes a global histogram as input, and updates the separator.
* @param globalHist : Global histogram.
* @param nbSeps : Number of separators to compute.
* @param sumNbItems : Global number of particles in the global histogram.
* @param localSep : The separator to update.
* @param numSep : The index of the separator to compute.
* @param nbUnderSep : The number of particles under the separator.
* @param chunkSize : number of bits
**/
void compSepM(const int * globalHist, const int & nbSeps, const int & sumNbItems,
	int & localSep, const int & numSep, int & nbUnderSep, int chunkSize)
{	

	int nbParts = nbSeps+1;		
	int cpt = nbUnderSep;
	int max = (sumNbItems/nbParts)*(numSep+1);
	int histSize = (1 << chunkSize); // = 2^chunkSize;

	// for each histogram bin
	for (int i=0; i<histSize; i++)
	{
		// increase the counter
		cpt += globalHist[i];			
		
		// if the counter exceeds the max		
		if (cpt >= max )
		{	
			// update localSep, counters and exit
			localSep = i; 
			cpt -= globalHist[i]; 
			nbUnderSep = cpt;
			break;
		}
	}
}

/**
* This function updates the flatten array of separators with the new separators stored in 
* the sepIdx array.
* @param flatIdxes : Array of separators, is updated here.
* @param flatIdxSize : Size of the array of separators.
* @param nbWorkers : Last Number of parts.
* @param nbSeps : Number of newly computed separators.
* @param sepIdx : Array of newly computed separators.
**/
void updateFlatten(int *& flatIdxes, int & flatIdxSize, int nbWorkers, int nbSeps, int ** sepIdx)
{
	int nbParts = nbSeps + 1;
  	flatIdxSize = (nbWorkers * nbParts) +1;
	
	// compute new flatIdxes tab
 	int * newFlatIdxes = new int[flatIdxSize];				
	int j=0;
	
	// update the new flatIdxes by inserting the Separator indexes
	for (int rank=0; rank<nbWorkers; rank++)
	{
		newFlatIdxes[j] = flatIdxes[rank];
		j++;
		for (int k=0; k<nbSeps; k++)
		{
			newFlatIdxes[j] = sepIdx[rank][k];
			j++;
		}
	}
	newFlatIdxes[j] = flatIdxes[nbWorkers];
	
	// swap the pointers
	int * aux;
	aux = flatIdxes;
	flatIdxes = newFlatIdxes;		
	delete [] aux;	
}


/**
* This function updates the array of particles that will be used to update the new nodes of particles.
* @param p : Array of particles, is updated here.
* @param nbWorkers : Last Number of workers.
* @param nbParts : Number of parts.
* @param flatIdxes : Array of separators.
**/
void fillParticles(Particles **& p, int nbWorkers, int nbParts, int * flatIdxes)
{
	// allocate the array
	p = new Particles* [nbWorkers]();
	for (int i=0; i<nbWorkers; i++)
		p[i] = new Particles[nbParts];	
	
	// update the particles using the flatten indexes tab
	int index = 0;
	for (int rank=0; rank<nbWorkers; rank++)
		for (int k=0; k<nbParts; k++)
		{
			p[rank][k].setAttrBounds(flatIdxes[index]+1, flatIdxes[index+1]);	
			index++;
		}
}

/**
* This function tests if the computed bits are enough to determine an octree grid separator.
* If yes, update the separators array and the corresponding tag.
* @param separators : Array of separators.
* @param nbWorkers : Number of processes taking part to the current computation.
* @param nbSeps : Number of separators to update.
* @param nbBits : Already computed bits.
* @param c : octree grid square side.
* @param h : octree height.
* @param k : separator index
**/
void adjustOnGrid(ui64 * separators, int nbWorkers, int nbSeps, int nbBits, ui32 c, ui32 h, int k)
{
	// Right (64 - nbBits)  bits to 1
	ui64 mask = (1l<<(64-nbBits)) - 1l;
	
	// Value of the smallest Separator
	int singleSep = COORDMAX/(1<<h); // 1<<h = 2^h

	double min, max, diff = 0;
	ui64 sepMin, sepMax = 0;

	// compute min and max, by filling with Zeros and Ones
	min = *((double *)&(separators[k]));		
	ui64 valMax = separators[k] | mask;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
	max = *((double *)&(valMax));	
#pragma GCC diagnostic pop	

	diff = max - min;
	
	// if min and max are separated by less than a square side
	if ( diff < c ) 
	{
		// fit the min and max seps to the nearest octree grid separators
		sepMin = round (min / singleSep) * singleSep;
		sepMax = round (max / singleSep) * singleSep;
		
		// test if they are the same
		if (sepMin == sepMax)
		{
			// update tag
			separators[nbSeps] |= (1 << k);
			
			// cast the separator to double and update separators array
			double sep = static_cast<double>(sepMin);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
			separators[k] = *(ui64 *)&(sep);
#pragma GCC diagnostic pop	
		}
	}			
}
