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

#include "Decomposition.hpp"

decompo::decompo(int nb)
{
	create(nb);
}

void decompo::create(int nb)
{
	int p = 2;
	
	while (p <= nb)
	{
		if (nb%p == 0)
		{
			_list.push_back(p);			
			nb = nb/p;
		}
		else if (p>61)
		{
			cerr <<"The decomposition of the number of processes has a prime number > 61." << endl;
			cerr <<"Please modifiy the number of processes." << endl;

			exit(0);
		}
		else
			p++;
	}
	reverse(_list.begin(), _list.end());
}
