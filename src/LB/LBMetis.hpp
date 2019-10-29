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

#ifndef LB_METIS
#define LB_METIS

#include "Node.hpp"

using namespace std;

class MetisLB
{
public:
	template <typename T>
	void loadBalance(Node<T> * octree, const decompo & nb1ers, const double & dist, double tol,
		const int & first, const int & last, const double & maxEdge, const vec3D & center, Gaspi_communicator & gComm, 
		double * nodeCenters, i64 * nodeOwners, int nbLeaves, i64* fnear, i64* nnear, i64* near, i64* fson, i64* nson) const; 
	
    template <typename T>
    void computeIS(Node<T> * octree, i64* fnear, i64* nnear, i64* near, i64* fson, i64* nson) const;

};

template <typename T>
void MetisLB::loadBalance(Node<T> * octree, const decompo & nb1ers, const double & dist, double tol,
    const int & first, const int & last, const double & maxEdge, const vec3D & center, Gaspi_communicator & gComm, 
    double * nodeCenters, i64 * nodeOwners, int nbLeaves, i64* fnear, i64* nnear, i64* near, i64* fson, i64* nson) const
{ 
    cout << "--> MetisLB load balancing strategy" << endl;
 
    /**
     * 1 - Compute the interaction set for each node
     **/
    computeIS(octree, fnear, nnear, near, fson, nson);
   
    /**
    * 2 - Prepare Metis' Data Structures
    **/
            

    /**
    * 3 - Call to Metis
    **/

    /**
    * 4 -  C to F +1 pour les ranks MPI au niveau des feuilles
    **/
    for (int i=0; i<nbLeaves; i++)
        nodeOwners[i] = i%4 + 1;        // juste pour éviter la segfault ...
}

template <typename T>
void MetisLB::computeIS(Node<T> * octree, i64* fnear, i64* nnear, i64* near, i64* fson, i64* nson) const
{
    ui32 height = octree->getHeight() + 1; // au sens spectre
    ui32 level = octree->getDepth() + 1; // spectre
    ui32 nbElem = octree->getContent().getNbParticles();
    ui32 coeff =  pow(2,height - (level));
          
    
    cout << "FMMLIB : height : " << height + 1 << endl;
    
    //string prefix, int nodeID_f, int level, i64 height, i64 * fathers, i64* nelem, i64 * nson, i64 * fson, i64 * fnear, i64 * nnear, i64 * near
    
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nodeID_c = octree->getId();
    int nodeID_f = nodeID_c + 1;
    
    string prefix = "dump/octree/LB_oct";
     // Open the output file
    if (rank == 0) 
    {
        string file = prefix + "_" + to_string ((unsigned long long)rank) + ".txt";
        ofstream out;
        out.open (file, std::ofstream::out | std::ofstream::app);

        if (nodeID_c == 0) 
        {
                out << "---------------- ID (Fortran) : " << nodeID_f <<  " ---------------- "  << endl;
        }
        else
        {
            // ID, nbElem, coeff, father
            int64_t father_idxC = octree->getParent()->getId(); 
            out << "---------------- ID (Fortran) : " << nodeID_f <<  " ---------------- "  << endl;
            out << "current level : "  << level << endl;
            out << "oct height  : "  << height << endl;
            out << "nbElem : " << nbElem << endl;
            out << "coeff : "  << coeff << endl;
            out << "father : " << father_idxC +1 << endl;

            // sons
            int nbSons = octree->getNbChildren();
            //int firstSon_IdxF = fson[nodeID_c];            
            if (nbSons)
                int firstSon_IdxC = octree->getChildren()[0]->getId();//firstSon_IdxF - 1;
            
            out << "nb sons : " <<  nbSons << endl;
            for (int i=0; i< nbSons; i++)
            {
                int son_IdxC = octree->getChildren()[i]->getId();
                int son_idxF = son_IdxC + 1;              
                out << "\t" << i <<" th son : " << son_idxF << endl;
            }
            
            // Interaction Set : -> get Father(nodeIndex), -> get near(Father)
            int fatherFirstNearCell_idxF = fnear[father_idxC];
            int fatherFirstNearCell_idxC = fatherFirstNearCell_idxF - 1;
            int fatherNbNearCells = nnear[father_idxC];
            out << "nb cells in father's near zone : " << fatherNbNearCells << endl;            

            // Parcours de tous les nodes en zone near du father
            for (int inear=0; inear<fatherNbNearCells; inear++)
            {
                int currentFatherNearCell_IdxF = near[fatherFirstNearCell_idxC + inear];
                int currentFatherNearCell_IdxC = currentFatherNearCell_IdxF - 1;
                out << "\t" << inear <<" th near cell (F) : " << currentFatherNearCell_IdxF << endl;               
                
                // Récupère chacun des enfants du #inear
                int currentfatherNearCell_nbSons = nson[currentFatherNearCell_IdxC]; /*octree->getNodePtr(currentFatherNearCell_IdxC)->getNbChildren();*/

                int currentFatherNearCell_fson_idxF = fson[currentFatherNearCell_IdxC];
                int currentFatherNearCell_fson_idxC = currentFatherNearCell_fson_idxF - 1;
                out << inear <<" which has nb sons : " << currentfatherNearCell_nbSons << endl;
                
                for (int i=0; i<currentfatherNearCell_nbSons; i++)
                {
                    // teste --> currentFatherNearCell_son_IdxC
                    int currentFatherNearCell_son_IdxC = currentFatherNearCell_fson_idxC + i; /*octree->getNodePtr(currentFatherNearCell_IdxC)->getChildren()[i]->getId();*/
                    
                    // exclusion, si appartient à la zone near du noeud
                    int nodeFirstNearCell_idxF = fnear[nodeID_c];
                    int nodeFirstNearCell_IdxC = nodeFirstNearCell_idxF - 1;                    
                    int nodeNbNearCells = nnear[nodeID_c];
                    
                    int found = 0;
                    // parcours tous les near du node
                    for (int jnear=0; jnear < nodeNbNearCells; jnear++)
                    {
                        int nodeNearCell_idxF = near[nodeFirstNearCell_IdxC + jnear];
                        int nodeNearCell_idxC = nodeNearCell_idxF -1;
                        if (nodeNearCell_idxC == currentFatherNearCell_son_IdxC)
                            found = 1;
                    }
                    if (found == 0)  
                        out << "  Node : " << nodeID_c << "\t" << currentFatherNearCell_son_IdxC + 1 << " Interaction Set" << endl;
                    else
                        out << "  Node : " << nodeID_c << "\t" << currentFatherNearCell_son_IdxC + 1 << " excluded" << endl;
                }
            }
        }
        
        // for each son, recursive call
        int nbSons  = octree->getNbChildren();
        for (int i=0; i<nbSons; i++)
        {
            computeIS(octree->getChildren()[i], fnear, nnear, near, fson ,nson);
        }

        // close the output file
        if (nodeID_c == 0) 
        {
            out << "--- END of FILE ---" << endl;
            out.close();
        }
    }
}

#endif
