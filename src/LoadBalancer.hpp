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
	Gaspi_communicator * _gComm;

public:
	LoadBalancer<T,LBPolicy>(Node<T> * n, const decompo & nb1ers, const double & dist, double tol,
		const int & first, const int & last, Gaspi_communicator * gComm)
		: _tree(n)
		, _nb1ers(nb1ers)
		, _dist(dist)
		, _tol(tol)
		, _first(first)
		, _last(last)
		, _gComm(gComm)
	{
		cout << "Load Balancer constructor " << endl;
	}
	
	virtual void run() const { loadBalance(_tree, _nb1ers, _dist, _tol, _first, _last, *_gComm); }
};

#endif
