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

#ifndef DECOMPOSITION_HPP
#define DECOMPOSITION_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

struct decompo
{	
	// attributes
	vector<int> _list;
	
	// methods
	decompo(): _list() {}		
	decompo(int nb);	
	decompo(const decompo & d)
	{
		_list.resize(d._list.size());
		copy(d._list.begin(), d._list.end(), _list.data());
	}		
	~decompo(){}


	void create(int nb);	
	int operator[](const int & index) const { return _list[index];}
	void display() const
	{
	    copy(_list.begin(), _list.end(), ostream_iterator<int>(cout," "));
	}
};

#endif
