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

#ifndef VEC3D_HPP
#define VEC3D_HPP

#include <iostream>
using namespace std;

struct vec3D{	

	union
	{
		struct
		{
			double x;
			double y;
			double z;
		};
		double v[3];
	};

	vec3D(double a=0, double b=0, double c=0): x(a), y(b), z(c){}
	bool operator ==(const vec3D & v) { return x==v.x && y==v.y && z==v.z; }
	bool operator !=(const vec3D & v) { return !(this->operator ==(v)); } 
	double & operator[](int i){ return v[i]; }
	double operator[](int i) const { return v[i]; }  
};

ostream& operator << (ostream& out, const vec3D & v);

#endif
