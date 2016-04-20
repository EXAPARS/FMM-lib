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
  the FMM-lib. If not, see <http://www.gnu.ogit rg/licenses/>.
*/

#ifndef LB_MORTON_ASYNC_GASPI_HPP
#define LB_MORTON_ASYNC_GASPI_HPP

#include "Node.hpp"
#include "LBMortonBase.hpp"

using namespace std;

class MortonAsyncGASPI : public LBMortonBase
{
public:
	template <typename T>
	void loadBalance(Node<T> * n, const decompo & nb1ers, const double & dist, double tolerance, 
		const int & first, const int & last, Gaspi_communicator & gComm) const;	
	
	template<typename T>
	void WaitForInitBuffers(Node<T> * n, const int & nbLeaves, int *& globalBuffer, 
		Gaspi_communicator & gComm) const;

	template<typename T>
	void initializeSepNodes(Node<T> * n, int * globalBuffer, const int & nbLeaves, const int & nbSeps, 
		int * targets, const int & divRoot, const int & divHeight, Gaspi_communicator & gComm) const;
		
	template<typename T>
	void updateDiffAndSendFirstRequests(Node<T> * n, const int & nbSeps, 
		double * diff, int * targets, const double & tolerance, const int & meanNbItems, int & nbExpectedUpdates,
		int * myBuffer, int * sepMarkers, Gaspi_communicator & gComm) const;
	
	template<typename T>
	void handleBufferRequest(Node<T> * n, int owner, int nbLeaves, Gaspi_communicator & gComm) const;

	template<typename T>
	void StateMachine(Node<T> * n, const int & nbSeps, int & sumNbItems, int & meanNbItems, 
		int * targets, double * diff, double tolerance, int & nbExpectedUpdates, int * myBuffer, 
		int * sepMarkers, const int & first, const int & last, Gaspi_communicator & gComm) const;
	
	template<typename T>	
	void testNbUntilNode(Node<T> * n, const bool & receivedInitSepNodes, bool & receivedInitNbUntil,
		int & meanNbItems, int * targets, const int & nbSeps, double * diff, const double & tolerance, 
		int & nbExpectedUpdates, int * myBuffer, int * sepMarkers, Gaspi_communicator & gComm) const;
	
	template<typename T>	
	void testSepNodes(Node<T> * n, bool & receivedInitSepNodes, const bool & receivedInitNbUntil,
		int & meanNbItems, int * targets, const int & nbSeps, double * diff, const double & tolerance, 
		int & nbExpectedUpdates, int * myBuffer, int * sepMarkers, int & nbSepUpdates, const int & first, 
		const int & last, int & nbNewParticles, int & nbCompletedCoordsUpdates,	int nbLeaves,
		Gaspi_communicator & gComm) const;

	template<typename T>	
	void handleSepUpdate(Node<T> *n, const int & source, int & nbSepUpdates, int * sepMarkers,
		const int & first, const int & last, int & nbNewParticles, int & nbCompletedCoordsUpdates,
		Gaspi_communicator & gComm) const;

	template<typename T>	
	void testMarkers(Node<T> * n, const int & sepID, int * sepMarkers, const int & first, 
		const int & last, int & nbNewParticles, int & nbCompletedCoordsUpdates, 
		Gaspi_communicator & gComm) const;
	
	template<typename T>	
	void sendCoords( Node<T> * n, const int & chunk, const int & first, const int & last, 
		const int & fromSepID, const int & toSepID, const int & dest, int & nbNewParticles, 
		int & nbCompletedCoordsUpdates, Gaspi_communicator & gComm) const;
	
	template<typename T>	
	void computeBounds(Node<T> * n, const int & chunk, const int & first, const int & last, 
		const int & fromSepID, const int & toSepID, int & nbParticlesToSend, int & beginCoord, 
		int & firstParticle, int & firstIndex, int & lastIndex, Gaspi_communicator & gComm) const;
	
	template<typename T>		
	void testReceiveBuffer(Node<T> * n, int * myBuffer, int & nbReceivedBuffers, int * targets, 
		int & nbSepUpdates, const int & meanNbItems, const double & tolerance, int * sepMarkers,
		const int & first, const int & last, vec3D * newCoords, int & nbNewParticles, 
		int & nbCompletedCoordsUpdates, int nbSeps, int nbLeaves,
		Gaspi_communicator & gComm) const;
	
	template<typename T>	
	void handleBufferAnswer(Node<T> * n, int * myBuffer, int & nbReceivedBuffers, 
		int * targets, int & nbSepUpdates, const int & meanNbItems, const double & tolerance, 
		int * sepMarkers, const int & first, const int & last, vec3D * newCoords, int & nbNewParticles, 
		int & nbCompletedCoordsUpdates, int nbSeps, int nbLeaves,
		Gaspi_communicator & gComm) const;
	
	template<typename T>
	void sendBufferRequest(Node<T> * n, int * myBuffer, Gaspi_communicator & gComm) const;
	
	template<typename T>
	void SendInitBuffers(Node<T> * n, const int & nbLeaves, const int & destination,
		Gaspi_communicator & gComm) const;
	

	void updateTargets(const int & sumNbItems, int & meanNbItems, int * targets,
		const double & tolerance, const int nbSeps) const;
			
	void writeCoords(int & nbParticlesToSend, int & beginCoord, int & firstParticle,
		const int & dest, int & nbNewParticles, int & nbCompletedCoordsUpdates, 
		Gaspi_communicator & gComm) const;

	void sendSepUpdate(Gaspi_communicator & gComm) const;
	void testNewCoords(int & nbCompletedCoordsUpdates, Gaspi_communicator & gComm) const;
	void testCommunicationInfos(int & nbNewParticles, Gaspi_communicator & gComm) const;		
};



template <typename T>
void MortonAsyncGASPI::loadBalance(Node<T> * n, const decompo & nb1ers, const double & dist, 
	double tolerance, const int & first, const int & last, Gaspi_communicator & gComm) const
{ 
	// MPI Barrier
	MPI_Barrier(MPI_COMM_WORLD);	
	
	// Gaspi rank and size
	gaspi_rank_t rank, wsize;
	SUCCESS_OR_DIE(gaspi_proc_rank(&rank));
	SUCCESS_OR_DIE(gaspi_proc_num(&wsize));

	// LB init
	int nbParts = wsize;
	int nbSeps = nbParts - 1;
			
	// Octree construction with height = 3	
	int divRoot = 0;	
	int divHeight = 3;
	n->recDivideOctreeH(divHeight);
	
	// diff
	double * diff = new double [nbSeps]();
		
	// variables
	int nbLeaves = 512;
	int * myBuffer = new int [nbLeaves]();
	int * targets = new int [nbSeps]();
	int meanNbItems;
	int nbExpectedUpdates = 0;
	int * sepMarkers = new int [nbSeps];
	for (int i=0; i<nbSeps; i++)
		sepMarkers[i] = 1;

	// SumNbItems
	int sumNbItems = 0;                                                 	
	int nbItems = n->getContent().getNbParticles();	
	gaspi_allreduce(&nbItems, &sumNbItems, 1, GASPI_OP_SUM, GASPI_TYPE_INT, GASPI_GROUP_ALL, GASPI_BLOCK);

	// new coordinates array
	int maxNbItems = ceil((sumNbItems/wsize) * (1 + 2*tolerance));
	gComm.initNewCoords(maxNbItems);	

	/** TODO : virer la barriere et faire le segment create en 2 fonctions **/				
	SUCCESS_OR_DIE(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));	
						
	int root = 0;
	
	if (rank == root)
	{
		//Phase 1
		int * globalBuffer;
		WaitForInitBuffers(n, nbLeaves, globalBuffer, gComm);
		updateTargets(sumNbItems, meanNbItems, targets, tolerance, nbSeps);
			
		initializeSepNodes(n, globalBuffer, nbLeaves, nbSeps, targets, divRoot, divHeight, gComm);
				
		updateDiffAndSendFirstRequests(n, nbSeps, diff, targets, tolerance, meanNbItems, 
			nbExpectedUpdates, myBuffer, sepMarkers, gComm);	
	
		StateMachine(n, nbSeps, sumNbItems, meanNbItems, targets, diff, tolerance, 
			nbExpectedUpdates, myBuffer, sepMarkers, first, last, gComm);
	
	}
	else
	{
		SendInitBuffers(n, nbLeaves, root, gComm);
		updateTargets(sumNbItems, meanNbItems, targets, tolerance, nbSeps);
		StateMachine(n, nbSeps, sumNbItems, meanNbItems, targets, diff, tolerance,
			nbExpectedUpdates, myBuffer,sepMarkers, first, last, gComm);	
	}
}

				/* *****************************************
				 * 	   RANK 0 INITIALIZATIONS              *	
				 ******************************************/

/**
 *  Rank 0 waits for all buffers
 **/
template<typename T>
void MortonAsyncGASPI::WaitForInitBuffers(Node<T> * n, const int & nbLeaves, int *& globalBuffer, 
	Gaspi_communicator & gComm) const
{		
	// receives buffers from the others
	int cpt = 0;
	int nbBuffers = gComm._wsize - 1;
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val;

	while (cpt < nbBuffers)
	{
		// wait for notification
		gaspi_notify_waitsome(
			gComm._seg_RecvBuffer_id,
			0,							// surveille les notifications depuis 0
			gComm._wsize,				// en surveille wsize
			&new_notif_id,
			GASPI_TEST
		);

		// get notification value, and reset
		gaspi_notify_reset(
			gComm._seg_RecvBuffer_id,
			new_notif_id,
			&new_notif_val
		);
		
		if (new_notif_val == INIT_LEAVES_BUFFER)
			cpt++;
	}
	
	// 0 writes its own leaves info into globalBuffer
	globalBuffer = new int [nbLeaves]();
	n->FillSendBuffer(globalBuffer, 3);
			
	// add into global buffer	
	for (int i=0; i<gComm._wsize; i++)	
		for (int j=0; j<nbLeaves; j++)
			globalBuffer[j] += gComm._recvBuffer[i*nbLeaves +j];			
}



/**
 *  Rank 0 initalizes the separators computation :
 * 	SepNodes and nbUntilNode
 *  And sends them to all others
 **/
template<typename T> 
void MortonAsyncGASPI::initializeSepNodes(Node<T> * n, int * globalBuffer, const int & nbLeaves, const int & nbSeps, 
	int * targets, const int & divRoot, const int & divHeight, Gaspi_communicator & gComm) const
{				
	// Update sepNodes and nbUntilNode
	computeMortonSeps(n, globalBuffer, nbLeaves, nbSeps, targets, gComm._nbUntilNode, gComm._sepNodes, 
		divRoot, divHeight);
		
	delete globalBuffer;
	
	int notify_ID = gComm._rank;
	int notify_VAL = INIT_SEP_NODES;	
	for (int i=1; i<gComm._wsize; i++)
	{
		// sepNodes
		gaspi_write_notify(	gComm._seg_SepNodes_id,		// local seg
							0, 						// local offset
							i,			 			// receiver rank
							gComm._seg_SepNodes_id,		// remote seg 
							0, 						// remote offset
							nbSeps * sizeof(int64_t), // size of data to write
							notify_ID,				// remote notif ID
							notify_VAL,				// value of the notif to write
							0, 						// queue
							GASPI_BLOCK
		);		
	}
	
	// Send sepNodes and nbUntilNode with notifications
	notify_ID = gComm._rank;
	notify_VAL = INIT_NB_UNTIL_NODE;	
	for (int i=1; i<gComm._wsize; i++)
	{
		// nbUntilNode
		gaspi_write_notify(	gComm._seg_NbUntilNode_id,		// local seg
							0, 						// local offset
							i,			 			// receiver rank
							gComm._seg_NbUntilNode_id,		// remote seg 
							0, 						// remote offset
							nbSeps * sizeof(int), 	// size of data to write
							notify_ID,				// remote notif ID
							notify_VAL,				// value of the notif to write
							0, 						// queue
							GASPI_BLOCK
		);		
	}
}

/**
 * update the diff array, and send a buffer request if needed.
 */
template<typename T>
void MortonAsyncGASPI::updateDiffAndSendFirstRequests(Node<T> * n, const int & nbSeps, 
	double * diff, int * targets, const double & tolerance, const int & meanNbItems, int & nbExpectedUpdates,
	int * myBuffer, int * sepMarkers, Gaspi_communicator & gComm) const
{
	// update diff array - (diff without box)
	for (int i=0; i<nbSeps; i++)
		diff[i] = abs((targets[i] - gComm._nbUntilNode[i]))/(meanNbItems*1.0);			
	
	// count how many seps there are to refine
	for (int i=0; i<nbSeps; i++)
		if (diff[i] > tolerance)
		{
			nbExpectedUpdates++;
			sepMarkers[i] = 0;
		}
		
	if (gComm._rank)
	{
		int mySep = gComm._rank - 1;
		if (diff[mySep] > tolerance)
			sendBufferRequest(n, myBuffer, gComm);	
	}
}


				/* *****************************************
				 * 	   OTHER RANKS  INITIALIZATIONS        *	
				 ******************************************/
/**
 * This method is called by all processes except rank 0.
 * 
 **/
template<typename T>
void MortonAsyncGASPI::SendInitBuffers(Node<T> * n, const int & nbLeaves, const int & destination,
	Gaspi_communicator & gComm) const
{
	// fill localBuffer at recipient offset
	n->FillSendBuffer(&(gComm._localBuffer[destination * nbLeaves]), 3);
	
	// Gaspi write and notify the receiver	
	int remote_offset = gComm._rank * nbLeaves * sizeof(int);
	int local_offset = destination * nbLeaves * sizeof(int);
	int notify_ID = gComm._rank;
	int notify_VAL = INIT_LEAVES_BUFFER;
	
	gaspi_write_notify(	gComm._seg_LocalBuffer_id,
						local_offset,						// local offset
						destination,		 				// dreceiver rank
						gComm._seg_RecvBuffer_id,			// remote seg 
						remote_offset, 						// remote offset
						nbLeaves * sizeof(int), 			// size of data to write
						notify_ID,							// remote notif ID
						notify_VAL,							// value of the notif to write
						0, 									// queue
						GASPI_BLOCK
	);
}


			/* *****************************************
			 * 				STATE MACHINE              *	
			 ******************************************/

template<typename T>
void MortonAsyncGASPI ::StateMachine(Node<T> * n, 
	const int & nbSeps, 
	int & sumNbItems, 
	int & meanNbItems, 
	int * targets, 
	double * diff, 
	double tolerance, 
	int & nbExpectedUpdates,
	int * myBuffer, 
	int * sepMarkers, 
	const int & first, 
	const int & last,
	Gaspi_communicator & gComm) const
{

	int nbSepUpdates = 0;
	int nbReceivedBuffers = 0;			
	int nbNewParticles = 0;
	int nbCompletedCoordsUpdates = 0;
	
	bool receivedInitSepNodes = false;
	bool receivedInitNbUntil = false;		
	bool receivedAllSepUpdates = false;
	bool receivedAllCoordinates = false;
	bool finished = false;

	int nbLeaves = 512;	
			
	while ( !(finished))	
	{	
		// test NB_UNTIL_NODE
		if( !receivedInitNbUntil )
			testNbUntilNode(n, receivedInitSepNodes, receivedInitNbUntil,
				meanNbItems, targets, nbSeps,
				diff, tolerance, nbExpectedUpdates, myBuffer, sepMarkers, gComm);
		
		// test SEP_NODES : INIT, UPDATE or REQUEST
		testSepNodes(n, receivedInitSepNodes, receivedInitNbUntil, meanNbItems, targets, nbSeps, 
			diff, tolerance, nbExpectedUpdates, myBuffer, sepMarkers, nbSepUpdates, first, last, 
			nbNewParticles, nbCompletedCoordsUpdates, nbLeaves, gComm);
		
		// test RECV_BUFFER : LEAVES_BUFFER_ANSWER		
		testReceiveBuffer(n, myBuffer, nbReceivedBuffers, targets, 
			nbSepUpdates, meanNbItems, tolerance, sepMarkers, first, last, gComm._newCoords, 
			nbNewParticles,	nbCompletedCoordsUpdates, nbSeps, nbLeaves, gComm);		
				
		// test NEWCOORDS
		testNewCoords(nbCompletedCoordsUpdates, gComm);
				
		// test COMMINFOS
		testCommunicationInfos(nbNewParticles, gComm);
		
		// state update
		if ( nbExpectedUpdates && nbSepUpdates && (nbExpectedUpdates == nbSepUpdates) )
			receivedAllSepUpdates = true;

		if ( nbCompletedCoordsUpdates == gComm._wsize )
			receivedAllCoordinates = true;
		finished = receivedAllCoordinates && receivedAllSepUpdates;	
	}
	
	// modifier le tableau de particules
	n->getContent().setNewCoordinates(gComm._newCoords, nbNewParticles);
}


/* ***********************************
 * STATE MACHINE Gaspi test Messages *
 * ***********************************/

template<typename T>
void MortonAsyncGASPI::testReceiveBuffer(Node<T> * n,
	int * myBuffer, int & nbReceivedBuffers,
	int * targets, int & nbSepUpdates,
	const int & meanNbItems, const double & tolerance, int * sepMarkers,
	const int & first, const int & last, vec3D * newCoords, int & nbNewParticles, 
	int & nbCompletedCoordsUpdates,
	int nbSeps, 
	int nbLeaves,
	Gaspi_communicator & gComm) const
{
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val = 0;
	
	if(gaspi_notify_waitsome(
		gComm._seg_RecvBuffer_id,
		0,				
		gComm._wsize,			
		&new_notif_id,
		GASPI_TEST) == GASPI_SUCCESS)		
	{
		gaspi_notify_reset(
			gComm._seg_RecvBuffer_id,
			new_notif_id,
			&new_notif_val);			
		
		if ( new_notif_val == LEAVES_BUFFER_ANSWER)
		{
			int nbLeaves = 512;
			handleBufferAnswer(n, myBuffer, nbReceivedBuffers, targets, 
				nbSepUpdates, meanNbItems, tolerance, sepMarkers, first, last, gComm._newCoords, 
				nbNewParticles,	nbCompletedCoordsUpdates, nbSeps, nbLeaves, gComm);
		}						
	}
}

template<typename T>
void MortonAsyncGASPI::testNbUntilNode(Node<T> * n, const bool & receivedInitSepNodes, bool & receivedInitNbUntil,
	int & meanNbItems, int * targets, const int & nbSeps, 
	double * diff, const double & tolerance, int & nbExpectedUpdates, 
	int * myBuffer, int * sepMarkers,
	Gaspi_communicator & gComm) const
{
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val = 0;

	// test if there is a message
	if(gaspi_notify_waitsome(
		gComm._seg_NbUntilNode_id,
		0,						// surveille les notifications depuis 0
		gComm._wsize,			// en surveille wsize
		&new_notif_id,
		GASPI_TEST) == GASPI_SUCCESS)
	{
		// get notification value, and reset
		gaspi_notify_reset(
			gComm._seg_NbUntilNode_id,
			new_notif_id,
			&new_notif_val);
		
		// update marker & handle
		receivedInitNbUntil = true;				
		
		// if possible update and send requests
		if (receivedInitSepNodes && receivedInitNbUntil)
			updateDiffAndSendFirstRequests(n, nbSeps, diff, targets, tolerance, meanNbItems, nbExpectedUpdates, 
				myBuffer, sepMarkers, gComm);
	}
}

template <typename T>
void MortonAsyncGASPI::testSepNodes(Node<T> * n, bool & receivedInitSepNodes, const bool & receivedInitNbUntil,
	int & meanNbItems, int * targets, const int & nbSeps, double * diff, const double & tolerance, 
	int & nbExpectedUpdates, int * myBuffer, int * sepMarkers, int & nbSepUpdates, const int & first, 
	const int & last, int & nbNewParticles, int & nbCompletedCoordsUpdates,	int nbLeaves,
	Gaspi_communicator & gComm) const
{	
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val = 0;	
	
	if(gaspi_notify_waitsome(
		gComm._seg_SepNodes_id,
		0,				
		gComm._wsize,			
		&new_notif_id,
		GASPI_TEST) == GASPI_SUCCESS)
	{
		int source = (int)new_notif_id;
		gaspi_notify_reset(gComm._seg_SepNodes_id, new_notif_id, &new_notif_val);

		if ( (new_notif_val == INIT_SEP_NODES) && (!receivedInitSepNodes) )
		{
			receivedInitSepNodes = true;
			if (receivedInitSepNodes && receivedInitNbUntil)
				updateDiffAndSendFirstRequests(n, nbSeps, diff, targets, 
					tolerance, meanNbItems, nbExpectedUpdates, myBuffer, sepMarkers, gComm);					
		}
		else if ( new_notif_val == UPDATE_SEP_NODES )
		{
			handleSepUpdate(n, source, nbSepUpdates, sepMarkers, first, last, nbNewParticles, 
				nbCompletedCoordsUpdates, gComm);					
		}
		else if ( new_notif_val == LEAVES_BUFFER_REQUEST)
		{
			handleBufferRequest(n, source, nbLeaves, gComm);
		}
	}
}


/* *************************************
 * STATE MACHINE Gaspi HANDLE Messages *
 * *************************************/

/**
 * Reception of a buffer request, from someone needing a refinement
 * Send the buffer answer **/ 
template<typename T>
void MortonAsyncGASPI::handleBufferRequest(Node<T> * n, int owner, int nbLeaves, 
	Gaspi_communicator & gComm)	const
{
	
	// get nodeID to refine		
	int ownerSepID = owner - 1;
	int64_t sepNodeID = gComm._sepNodes[ownerSepID];
		
	// Get pointer on node
	Node<T> * node = n->getNodePtr(sepNodeID);
		
	// Refine and fill buffer
	// Fill it with leaves information, at the recipient offset
	if ( !node->getChildren() )
	{
		node->divideOctreeNTimes(3);
		node->FillSendBuffer(&(gComm._localBuffer[owner * nbLeaves]), node->getDepth() + 3);
	}
	else
	{
		node->FillSendBuffer(&(gComm._localBuffer[owner * nbLeaves]), node->getDepth() + 3);
	}
		
	// Gaspi write and notify the receiver	
	int remote_offset = gComm._rank * nbLeaves * sizeof(int);
	int local_offset = owner * nbLeaves * sizeof(int);
	int notify_ID = gComm._rank;
	int notify_VAL = LEAVES_BUFFER_ANSWER;
	
	gaspi_write_notify(	gComm._seg_LocalBuffer_id,			// local seg
						local_offset,						// local offset
						owner,		 						// dreceiver rank
						gComm._seg_RecvBuffer_id,			// remote seg 
						remote_offset, 						// remote offset
						nbLeaves * sizeof(int), 			// size of data to write
						notify_ID,							// remote notif ID
						notify_VAL,							// value of the notif to write
						0, 									// queue
						GASPI_BLOCK
	);

}

/**
 * Reception of a buffer answer, consequently to a buffer request
 * Register the new information
 * If all buffers have been received, compute the new separator
 * Test if it is ok with respect to th tolerance criteria
 * If yes, send a separtor update
 * If not, send a buffer request for refinement **/
template<typename T>
void MortonAsyncGASPI::handleBufferAnswer(Node<T> * n, int * myBuffer, int & nbReceivedBuffers, 
	int * targets, int & nbSepUpdates, 
	const int & meanNbItems, const double & tolerance, int * sepMarkers,
	const int & first, const int & last, vec3D * newCoords, int & nbNewParticles, 
	int & nbCompletedCoordsUpdates,
	int nbSeps, 
	int nbLeaves,
	Gaspi_communicator & gComm) const
{		
	// test if has received all needed buffers
	nbReceivedBuffers++;	
	
	// update the separator
	if ( nbReceivedBuffers == (gComm._wsize-1) )
	{
		// sum into myBuffer		
		for (int i=0; i<gComm._wsize; i++)
			for (int j=0; j<nbLeaves; j++)
				myBuffer[j] += gComm._recvBuffer[i*nbLeaves + j];

		// refine
		int mySep = gComm._rank-1;
		computeMortonOneSep(n, myBuffer, 512, targets[mySep], gComm._nbUntilNode[mySep], 
			gComm._sepNodes[mySep], 3); 

		// test if ok with tolerance
		double diff = abs((targets[mySep] - gComm._nbUntilNode[mySep]))/(meanNbItems*1.0);
		
		if (diff > tolerance)
		{
			// re-init for the next step
			nbReceivedBuffers = 0;
			for (int i=0; i<512; i++)
				myBuffer[i] = 0;

			// Ask for the buffers
			sendBufferRequest(n, myBuffer, gComm);
		}
		else
		{
			// send and update counter
			sendSepUpdate(gComm);
			nbSepUpdates++;
			sepMarkers[mySep] = 1;
			testMarkers(n, mySep, sepMarkers, first, last, nbNewParticles, nbCompletedCoordsUpdates, gComm);
		}
	}
}

/**
 * Reception of a separator update.
 * Register the new information.
 * Test if there are particles ready to send.**/
template<typename T>
void MortonAsyncGASPI::handleSepUpdate(Node<T> *n, const int & source, int & nbSepUpdates, int * sepMarkers,
	const int & first, const int & last, int & nbNewParticles, int & nbCompletedCoordsUpdates,
	Gaspi_communicator & gComm) const
{
	nbSepUpdates++;
	sepMarkers[source-1] = 1;
	testMarkers(n, source-1, sepMarkers, first, last, nbNewParticles, nbCompletedCoordsUpdates, gComm);		
}

/** chunk position in the array
 * 0 = start, from the first index
 * 1 = middle
 * 2 = end, to the last index
 **/ 
template<typename T>
void MortonAsyncGASPI::testMarkers(Node<T> * n, const int & sepID, int * sepMarkers, 
	const int & first, const int & last, int & nbNewParticles,
	int & nbCompletedCoordsUpdates,
	Gaspi_communicator & gComm) const
{
	// size
	int nbSeps = gComm._wsize - 1;
	int lastSep = nbSeps - 1;
				
	// starting with the first particle, send to 0
	if (sepID == 0)
		sendCoords(n, 0, first, 0, 0, sepID, 0, nbNewParticles, nbCompletedCoordsUpdates, gComm);
	
	// ending with the last particle send to size-1
	if (sepID == lastSep)
		sendCoords(n, 2, 0, last, sepID, 0, gComm._wsize - 1, nbNewParticles, nbCompletedCoordsUpdates, gComm);
	
	// In the Middle, between sepID and sepID+1
	if (sepID < lastSep)
		if (sepMarkers[sepID + 1])
			sendCoords(n, 1, 0, 0, sepID, sepID+1, sepID+1, nbNewParticles, nbCompletedCoordsUpdates, gComm);
	
	// In the Middle, between sepID-1 and sepID
	if (sepID > 0)
		if (sepMarkers[sepID - 1]) 
			sendCoords(n, 1, 0, 0, sepID-1, sepID, sepID, nbNewParticles, nbCompletedCoordsUpdates, gComm);
}

/* *************************************
 * STATE MACHINE Gaspi SEND Messages *
 * *************************************/

/*** Send a buffer request to all other processes, and write the rank's data into its buffer.*/
template<typename T>
void MortonAsyncGASPI::sendBufferRequest(Node<T> * n, int * myBuffer, Gaspi_communicator & gComm) const

{		
	// get mySep and sepNode
	int mySep = gComm._rank - 1;
	int64_t sepNode = gComm._sepNodes[mySep];	

	// write parameters
	int remote_offset(mySep* sizeof(int64_t));
	int local_offest = remote_offset;
	int notify_ID = gComm._rank;
	int notify_VAL = LEAVES_BUFFER_REQUEST;
	
	for (int i=0; i<gComm._wsize; i++)
	{
		if (i != gComm._rank)
		{
			gaspi_write_notify(	gComm._seg_SepNodes_id,		// local seg
								local_offest,					// local offset
								i,			 					// dreceiver rank
								gComm._seg_SepNodes_id,		// remote seg 
								remote_offset, 					// remote offset
								sizeof(int64_t), 				// size of data to write
								notify_ID,						// remote notif ID
								notify_VAL,						// value of the notif to write
								0, 								// queue
								GASPI_BLOCK
			);
		}	
	}

	// And write my own data into myBuffer
	Node<T> * node = n->getNodePtr(sepNode);
	if ( !node->getChildren() )
	{
		node->divideOctreeNTimes(3);
		node->FillSendBuffer(myBuffer, node->getDepth() + 3);
	}
	else
	{
		node->FillSendBuffer(myBuffer, node->getDepth() + 3);
	}		
/** Todo : à nettoyer, tests à encapsuler dans la fonction
	n->divideOctreeNTimes(3);
	n->FillSendBuffer(myBuffer, n->getDepth() + 3);		
**/	
}




			/* ****************************************************
			*			PHASE 3 - ENVOI DE PARTICULES			  *	
			*******************************************************/


template<typename T>
void MortonAsyncGASPI::sendCoords( Node<T> * n, 
	const int & chunk, const int & first, const int & last, 
	const int & fromSepID, const int & toSepID, const int & dest, 
	int & nbNewParticles, int & nbCompletedCoordsUpdates,
	Gaspi_communicator & gComm) const
{
	// Send parameter
	int nbParticlesToSend = -1;
	int beginCoord = -1;
	int firstParticle = -1;		
	int firstIndex = -1;
	int lastIndex = -1;
	
	// Update the parameters
	computeBounds(n, chunk, first, last, fromSepID, toSepID, 
		nbParticlesToSend, beginCoord, firstParticle, firstIndex, lastIndex, gComm);
	
	// Write the coordinates
	writeCoords(nbParticlesToSend, beginCoord, firstParticle,
		dest, nbNewParticles, nbCompletedCoordsUpdates, gComm);
}

template<typename T>
void MortonAsyncGASPI::computeBounds(Node<T> * n, 
	const int & chunk, const int & first, const int & last, const int & fromSepID, const int & toSepID,
	int & nbParticlesToSend, int & beginCoord, int & firstParticle, int & firstIndex, int & lastIndex,
	Gaspi_communicator & gComm) const
{
	// starting with the first particle
	if (chunk == 0)
	{
		// update last
		lastIndex = n->getNodePtr(gComm._sepNodes[toSepID])->getContent().getLastIndex();	
		if (lastIndex == -1)
			n->findLastParticleIndex(gComm._sepNodes[toSepID], lastIndex);
		
		firstIndex = first;

		// parameters update
		nbParticlesToSend = lastIndex - firstIndex + 1;
		beginCoord = firstIndex * 3;
		firstParticle = firstIndex;		
	}
	
	// middle
	if (chunk == 1)
	{	
		// Update First
		firstIndex = n->getNodePtr(gComm._sepNodes[fromSepID])->getContent().getLastIndex();	
		if (firstIndex == -1)
			n->findLastParticleIndex(gComm._sepNodes[fromSepID], firstIndex);

		// Update Last
		lastIndex = n->getNodePtr(gComm._sepNodes[toSepID])->getContent().getLastIndex();	
		if (lastIndex == -1)
			n->findLastParticleIndex(gComm._sepNodes[toSepID], lastIndex);

		// parameters update
		nbParticlesToSend = lastIndex - firstIndex;
		beginCoord = (firstIndex + 1) * 3;
		firstParticle = firstIndex + 1;
	}

	// ending with last particles
	if (chunk == 2)
	{			
		// update First	
		firstIndex = n->getNodePtr(gComm._sepNodes[fromSepID])->getContent().getLastIndex();	
		if (firstIndex == -1)
			n->findLastParticleIndex(gComm._sepNodes[fromSepID], firstIndex);

		lastIndex = last;
		
		// parameters update
		nbParticlesToSend = lastIndex - firstIndex;
		beginCoord = (firstIndex + 1)* 3;
		firstParticle = firstIndex + 1;
	}
}

#endif
