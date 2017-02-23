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

#include "LBMortonAsyncGASPI.hpp"


/**
 * update the targets array and allocates the new coordinates array
 */
void MortonAsyncGASPI::updateTargets(const int & sumNbItems, int & meanNbItems, int * targets,
	const double & tolerance, const int nbSeps) const
{
	// compute mean number opf particles and update targets
	int wsize = nbSeps + 1;	
	meanNbItems = sumNbItems / wsize;
	for(int i=0; i<nbSeps; i++)
		targets[i] = meanNbItems*(i+1);
}


void MortonAsyncGASPI::writeCoords(int & nbParticlesToSend, int & beginCoord, int & firstParticle,
	const int & dest, int & nbNewParticles, int & nbCompletedCoordsUpdates,
	Gaspi_communicator & gComm) const
{	
	/// test if there are no particles to send	
	if (nbParticlesToSend == 0)
	{
		if (dest != gComm._rank)
		{	
			gaspi_notify(	gComm._seg_NewCoords_id,			// local seg
							dest,		 						// receiver rank
							gComm._rank,						// remote notif ID
							COORDS_EMPTY,						// value of the notif to write
							0, 									// queue
							GASPI_BLOCK
			);						
		}
		else
		{
			nbCompletedCoordsUpdates++;
		}
	}
	
	/// send the particles 
	else
	{			
		// parameters update
		int nbCoordsToSend = nbParticlesToSend * 3;		
		gaspi_atomic_value_t writeCoordIndex = 0;
			
		// Send the particles
		if (dest != gComm._rank)
		{	
			// update commInfos 
			/** ceci est sûr d'être terminé car si la notif arrive, 
			tout ce qui était avant dans la queue est traité **/
			
			gaspi_notify(	gComm._seg_CommInfos_id,			// local seg
							dest,		 						// receiver rank
							gComm._rank,						// remote notif ID
							nbParticlesToSend,					// value of the notif to write
							0, 									// queue
							GASPI_BLOCK
			);	
						
			// get and increase index
			gaspi_atomic_fetch_add (
					gComm._seg_CommInfos_id,
					0,												// offset
					dest,											// rank
					nbCoordsToSend,									// add
					&writeCoordIndex,									// old value
					GASPI_BLOCK
			);
			
			// write data and notify
			int local_offset = beginCoord * sizeof(double);	
			int remote_offset = writeCoordIndex * sizeof(double);
			
			gaspi_write_notify(	gComm._seg_InitCoords_id,			// local seg
								local_offset,						// local offset
								dest,		 						// receiver rank
								gComm._seg_NewCoords_id,			// remote seg 
								remote_offset,						// remote offset
								nbCoordsToSend * sizeof(double),	// size of data to write
								gComm._rank,						// remote notif ID
								COORDS_COMPLETED,					// value of the notif to write
								0, 									// queue
								GASPI_BLOCK
			);
		}

		else
		{		
			// get and increase index
			gaspi_atomic_fetch_add (
					gComm._seg_CommInfos_id,
					0,												// offset
					gComm._rank,									// rank
					nbCoordsToSend,									// add
					&writeCoordIndex,									// old value
					GASPI_BLOCK
			);
				
			// write the data into the buffer
			int writeParticleIndex = writeCoordIndex / 3;
			for (int i=0; i<nbParticlesToSend; i++)
			{
				gComm._newCoords[writeParticleIndex + i].x = gComm._initCoords[firstParticle + i].x;
				gComm._newCoords[writeParticleIndex + i].y = gComm._initCoords[firstParticle + i].y;
				gComm._newCoords[writeParticleIndex + i].z = gComm._initCoords[firstParticle + i].z;			
			}
			
			// update terminatedUpdates counter
			nbCompletedCoordsUpdates++;
			
			// update new particles counter
			nbNewParticles += nbParticlesToSend;
		}
	}		
}

/** Send a separator update to all other processes */	
void MortonAsyncGASPI::sendSepUpdate(Gaspi_communicator & gComm) const
{		 
	// write parameters
	int mySep = gComm._rank - 1;
	int remote_offset(mySep* sizeof(int64_t));
	int local_offest = remote_offset;
	int notify_ID = gComm._rank;
	int notify_VAL = UPDATE_SEP_NODES;
	
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
}

void MortonAsyncGASPI::testNewCoords(int & nbCompletedCoordsUpdates, Gaspi_communicator & gComm) const
{
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val = 0;	
	
	if(gaspi_notify_waitsome(
		gComm._seg_NewCoords_id,
		0,				
		gComm._wsize,			
		&new_notif_id,
		GASPI_TEST) == GASPI_SUCCESS)		
	{		
		// get notification value, and reset
		gaspi_notify_reset(gComm._seg_NewCoords_id, new_notif_id, &new_notif_val);			
		
		// update counter
		if ( new_notif_val == COORDS_COMPLETED || new_notif_val == COORDS_EMPTY) // prevoir les cas : COORDS_TO_BE_CONTINUED || EMPTY
			nbCompletedCoordsUpdates++;
	}
}

void MortonAsyncGASPI::testCommunicationInfos(int & nbNewParticles, Gaspi_communicator & gComm) const
{
	gaspi_notification_id_t new_notif_id;
	gaspi_notification_t new_notif_val = 0;	

	// from 0 to wsize-1 : update new number of particles
	if(gaspi_notify_waitsome(
		gComm._seg_CommInfos_id,
		0,				
		gComm._wsize,			
		&new_notif_id,
		GASPI_TEST) == GASPI_SUCCESS)		
	{		
		// get notification value, and reset
		gaspi_notify_reset(gComm._seg_CommInfos_id, new_notif_id, &new_notif_val);			
			
		// update counter of new particles
		nbNewParticles += new_notif_val;			
	}
}
