/*
  Copyright 2016 - UVSQ
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

#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <iostream>

using namespace std;

struct complex
{
	double re;
	double im;
	complex(double a=0, double b=0): re(a), im(b){}
	bool operator ==(const complex & c) { return re==c.re && im==c.im; }
	bool operator !=(const complex & c) { return !(this->operator ==(c)); }   
};

ostream& operator << (ostream& out, const complex & c);

#endif
