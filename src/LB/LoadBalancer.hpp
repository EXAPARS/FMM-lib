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

#ifndef LOAD_BALANCER_HPP
#define LOAD_BALANCER_HPP

#include <iostream>
using namespace std;

#include "Node.hpp"
#include "Decomposition.hpp"
#include "LoadBalancerBase.hpp"
#include "Gaspi_communicator.hpp"


template <typename T, typename LBPolicy>
class LoadBalancer : public LB_Base, private LBPolicy
{
	using LBPolicy::loadBalance;
	
private:
	Node<T> * _tree;
	const decompo _nb1ers;
	const double _dist;
	const double _tol;
	const int _first;
	const int _last;
	const double _maxEdge;
	const vec3D _center;
	Gaspi_communicator * _gComm;
	int64_t * _nodeOwners;
	double * _nodeCenters;
	const int _nbLeaves;

public:
	// constructeur pour utilisation DA
	LoadBalancer<T,LBPolicy>(Node<T> * n, const decompo & nb1ers, const double & dist, double tol,
		const int & first, const int & last, const double & maxEdge, const vec3D & center, Gaspi_communicator * gComm, double * nodeCenters, int64_t * nodeOwners, 
		const int & nbLeaves)
		: _tree(n)
		, _nb1ers(nb1ers)
		, _dist(dist)
		, _tol(tol)
		, _first(first)
		, _last(last)
		, _gComm(gComm)
		, _nodeOwners(nodeOwners)
		, _maxEdge(maxEdge)
		, _center(center)
		, _nodeCenters(nodeCenters)
		, _nbLeaves(nbLeaves)
	{}
	
	// constructeur pour utilisation FMM-viz 	
	LoadBalancer<T,LBPolicy>(Node<T> * n, const decompo & nb1ers, const double & dist, double tol,
		const int & first, const int & last, Gaspi_communicator * gComm)
		: _tree(n)
		, _nb1ers(nb1ers)
		, _dist(dist)
		, _tol(tol)
		, _first(first)
		, _last(last)
		, _gComm(gComm)
		/* non renseignés*/
		, _nodeOwners(nullptr)
		, _maxEdge(0)
		, _center(vec3D(0,0,0))
		, _nodeCenters(nullptr)
		, _nbLeaves(0)
	{}
	
	virtual void run() const { loadBalance(_tree, _nb1ers, _dist, _tol, _first, _last, _maxEdge, _center, *_gComm, _nodeCenters, _nodeOwners, _nbLeaves); }
	
	

};

#endif
