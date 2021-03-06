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

vector<string> msg;

void announce_axis(string axis, int rank)
{
	cout << " ======================= " << endl;
	cout << " =        "<< axis <<"   | "<< rank << " |" << endl;
	cout << " ======================= " << endl;
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

void debug(string prefix, string message)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	string file = prefix + "_" + to_string((unsigned long long)rank) + ".txt";
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
	out << message << endl;
	out.close();
}

void accumulMSG(string message)
{
	msg.push_back(message);
}

void dumpMSG(string prefix)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	string file = prefix + "_" + to_string((unsigned long long)rank) + ".txt";
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
	for (auto i:msg)
		out << i << " ";
	out << endl;
	out.close();
	msg.clear();	
}


string itoa(int a)
{
	 return to_string((long long)(a));
}

void dump_tree_init(int rank)
{
	string file = "output/dump_tree_" + to_string ((unsigned long long)rank);
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
    out << "digraph tree{\n\t";	
	out.close();
}

void dump_tree_add_child(int rank, int64_t parent, int64_t child, int nbParticles)
{
	string file = "output/dump_tree_" + to_string ((unsigned long long)rank);
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);
	
	out << parent << "->" << child << "[label=" << nbParticles << "]\n\t";
	
	out.close();
}

void dump_tree_close_file(int rank)
{
	string file = "output/dump_tree_" + to_string ((unsigned long long)rank);
	ofstream out;
	out.open (file, std::ofstream::out | std::ofstream::app);

    out << "}";	
	out.close();
}

void dfs_dump_spectre_octree(string prefix, i64 * nbElemPerNode, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nbNodes, i64 * nodeOwner, i64 nodeID, double * centers)
{
	// Open the output file
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) 
	{	
		string file = prefix + "cpp_dfs_dump_octree_" + to_string ((unsigned long long)rank) + ".txt";
		ofstream out;
		out.open (file, std::ofstream::out | std::ofstream::app);
		
		if (nodeID == 0) 
		{
			out << "digraph G{" << endl;
			out <<  nodeID << "[label=" << '"' << "[" <<  nodeID << "] " << "\\n " << nbElemPerNode[0] << '"' << "]"  << endl;
		}

		i64 nbSons = nbSonsPerNode[nodeID];
		i64 firstSonID = firstSonId[nodeID]-1;

		// for each son, write and call
		i64 sonID;
		i64 owner;
		//i64 nbElem;
		string style;
		string color;
		for (int i=0; i<nbSons; i++)
		{
			sonID = firstSonID + i;
			owner = nodeOwner[sonID];
			//nbElem = nbElemPerNode[sonID];
			out << nodeID << " -> " << sonID << ";" << endl;
			out << sonID << "[label=" << '"' << "[" << sonID << "] " << "\\n " /*<< nbElem << "-" << owner */<< "C(" << int(centers[sonID*3]) << "," << int(centers[(sonID*3)+1]) << ",\\n" << int(centers[(sonID*3)+2]) << ")" << '"' << "];"  << endl;
			/* TODO : Au secours que c'est moche --> faire un IO manip*/
			
			style = "filled";
			if (owner == 1) 
				color="red";
			else if (owner == 2) 
				color="orange";
			else if (owner == 3) 
				color="yellow";
			else if (owner == 4) 
				color="green";
			else if (owner == 5) 
				color="blue";
			else if (owner == 6) 
				color="pink";
			else if (owner == 7) 
				color="Coral";
			else if (owner == 8) 
				color="Sienna";
			else
			{
				color="black";
				style="solid";
			}
			out << sonID << "[color = " << color << ", style = " << style << "];" << endl;
			dfs_dump_spectre_octree(prefix, nbElemPerNode, nbSonsPerNode, firstSonId, nbNodes, nodeOwner, sonID, centers);
		}

		if (nodeID == 0)
		{
			out << "}";	
			out.close();
		}
	}
}

void dfs_dump_centers(string prefix, i64 * nbSonsPerNode, i64 * firstSonId, i64 nodeID, double * centers)
{
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) 
	{
		// Open the output file
		string file = prefix + to_string ((unsigned long long)rank) + ".txt";
		ofstream out;
		out.open (file, std::ofstream::out | std::ofstream::app);
		
		// Get sons Infos
		i64 nbSons = nbSonsPerNode[nodeID];
		i64 firstSonID = firstSonId[nodeID]-1;

		out << "nbSons = " << nbSons << endl;
		
		// for each son, write and call
		i64 sonID;
		for (int i=0; i<nbSons; i++)
		{
			sonID = firstSonID + i;
			out <</* "[" << sonID << "]" << */ centers[sonID*3] << " "<< centers[(sonID*3)+1] << " " << centers[(sonID*3)+2] << endl;
			dfs_dump_centers(prefix, nbSonsPerNode, firstSonId, sonID, centers);
		}

		if (nodeID == 0)
		{
			out.close();
		}
	}
}

void bfs_dump_centers_level_by_level(string prefix, i64 * nbElemPerNode, i64 * nbSonsPerNode, i64 * firstSonId, i64 * nbNodes, i64 * nodeOwner, 
	i64 nodeID, double * centers, i64 * endlev, i64 * nbLevels)
{
int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) 
	{
		// BFS
		for (int i=0; i< (*nbLevels); i++)
		{
			// Open the output file
			string file = prefix + "_rank_" + to_string ((unsigned long long)rank) + "_level_" + to_string((unsigned long long)i)+ ".txt";
			ofstream out;
			out.open (file, std::ofstream::out | std::ofstream::app);

			// write the centers
			if (i==0)
			{
				i64 boxID = 0;
				out << centers[boxID*3] << " "<< centers[(boxID*3)+1] << " " << centers[(boxID*3)+2] << endl;
				
			}
			else
			{
				// nbBoxes
				int nbBoxes = endlev[i] - endlev[i-1];
				
				// Box ID and center
				i64 firstBoxID = endlev[i-1];
				for (int j=0; j<nbBoxes; j++)
				{
					i64 boxID = firstBoxID + j;
					out << centers[boxID*3] << " "<< centers[(boxID*3)+1] << " " << centers[(boxID*3)+2] << endl;
				}
			}
			
			// close the output file
			out.close();
		} 
	} 
}

string convert(int a)
{
	return to_string((long long)(a));
}


void loadAndDiffData(const string & file1, const string & file2)
{
	double diffMax = 0.1;
	fstream in1, in2;
	in1.exceptions(ifstream::failbit|ifstream::badbit);
	in2.exceptions(ifstream::failbit|ifstream::badbit);
	try
	{
		in1.open(file1, ifstream::in);
		in2.open(file2, ifstream::in);

		// get number of lines
		int nb1 = std::count(std::istreambuf_iterator<char>(in1), std::istreambuf_iterator<char>(), '\n');
		int nb2 = std::count(std::istreambuf_iterator<char>(in2), std::istreambuf_iterator<char>(), '\n');
		
		nb1--;
		nb2--;
		
		if (nb1 == nb2)
		{
			// alloc arrays
			double * data1 = new double[nb1];
			double * data2 = new double[nb2];			
			
			// load data into arrays
			in1.seekg(0, ios::beg);
			in2.seekg(0, ios::beg);
			double tmp1, tmp2;

			int index = 0;
			while (index < nb1)
			{

				in1 >> tmp1;
				in2 >> tmp2;
				data1[index] = tmp1;
				data2[index] = tmp2;
				index++;
			}
			
			// check
			cout << "Verification at "<< diffMax << endl;
			int cpt = 0;
			vector<int> memo;
			for (int i=0; i<nb1; i++)
			{
				if (abs(data1[i] - data2[i]) > diffMax)
				{
					cpt++;
					memo.push_back(i);
				}
			}
			
			if (cpt == 0)
				cout << "Verification OK" << endl;
			else
			{
				cout << cpt << " values did not pass the test at : " << diffMax << endl;
			}
		}
		else
		{
			cout << "Verification aborted : Files don't have the same size !" << endl;
		}
		in1.close();
		in2.close();
	}
	catch(ifstream::failure &e)
	{
		cerr << "[erreur] " << e.what() << endl;
		exit(0);
	}
	
}


uint64_t rdtsc(void) {
	uint64_t a, d;
	__asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
	return (d<<32) | a;
}
