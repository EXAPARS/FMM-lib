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

#include "fmm_tools.hpp"

void announce_axis(string axis, int rank)
{
cout << "======================= " << endl;
cout << "=        "<< axis <<"        | "<< rank << " |" << endl;
cout << "======================= " << endl;	
	
}

void displayHexa2Dim(string info, ui64 ** tab, int dim1, int dim2)
{
	cout << info << endl;
	for (int i=0; i<dim1; i++)
		for (int j=0; j<dim2; j++)
			cout << std::hex << tab[i][j] << endl;
}

void displayDiff (string info, int * tab, int size)
{
	cout << info << endl;	
	for (int i=1; i<size; i++)
		cout << tab[i] - tab[i-1] << endl;
}


void displayMpiMSG (int source, int tag)
{
	switch (tag)
	{
		case NB_ITEMS :
			cout << " sumNbItems\n";
			break;
			
		case SEPS :
			cout << " separators\n";
			break;
			
		case NB_UNTIL :
			cout << " nbUnder\n";
			break;
			
		case BUFFER_REQUEST :
			cout << " buffer request\n";
			break;

		case SEPS_UPDATE :
			cout << " separator update\n";
			break;
			
		case BUFFER_ANSWER :
			cout << " buffer answer\n";
			break;
		
		case PARTICLES :
			cout << " new particles\n";
			break;
			
		default :
			cout <<" !!!! unhandled message : " << tag << endl;
	}
}

void verbose(int rank, string message)
{
	string file = "output/output_" + to_string(rank);
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
	out << message << endl;
	out.close();
}


void dump_tree_init(int rank)
{
	string file = "output/dump_tree_" + to_string (rank);
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
    out << "digraph tree{\n\t";	
	out.close();
}

void dump_tree_add_child(int rank, int64_t parent, int64_t child, int nbParticles)
{
	string file = "output/dump_tree_" + to_string (rank);
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
	
	out << parent << "->" << child << "[label=" << nbParticles << "]\n\t";
	
	out.close();
}

void dump_tree_close_file(int rank)
{
	string file = "output/dump_tree_" + to_string (rank);
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);

    out << "}";	
	out.close();
}
