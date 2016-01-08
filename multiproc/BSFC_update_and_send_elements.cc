
/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: BSFC_update_and_send_elements.C,v 1.1.1.1 2003/08/13 19:26:11 sorokine Exp $ 
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cassert>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <particle.h>

#include "repartition_BSFC.h"
#include "exvar.h"
#include "pack_data.h"
#include "multiproc.h"
//This function will only move non-brief buckets
//I will have another function for communication of only brief buckets.
void
BSFC_update_and_send_elements (int myid, int numprocs,
                               THashTable * P_table, HashTable * BG_mesh)
{
  int i, j, k;
  int dir1[DIMENSION], dir2[DIMENSION];
  Key mykey, neighkey;

#ifdef DEBUG2
  char filename[20];
  sprintf (filename, "debug%02d.txt", myid);
  FILE * fp = fopen (filename, "a+");
#endif

  int *send_info = new int[numprocs * 3];
  for (i = 0; i < 3 * numprocs; i++)
    send_info[i] = 0;

  // now figure out what processors need to have neighbor info updated 
  // and what processors to send objects to
  Bucket *buck;
  BriefBucket *breif_buck = NULL;
  void * tempptr =NULL;

  HTIterator *itr = new HTIterator (BG_mesh);

  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ()) // this element will get moved to a new processor
		  {
		     // neigh info
		     const int *neigh_proc = buck->get_neigh_proc ();

		     for (j = 0; j < NEIGH_SIZE; j++)
		        if ((neigh_proc[j] != myid) && (neigh_proc[j] > -1))
		          send_info[3 * neigh_proc[j]] += 1;

		     // bucket info
		     send_info[3 * (buck->get_myprocess ()) + 1] += 1;
		     send_info[3 * (buck->get_myprocess ()) + 2] +=
		        (int) buck->get_plist ().size ();
		  }

	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // copy send_info to temp_info
  int *temp_info = new int[3 * numprocs];
  for (i = 0; i < 3 * numprocs; i++)
    temp_info[i] = send_info[i];

  // send info to everyone, so everyone know how much 
  // they are receiving from whom
  int *recv_info = new int[3 * numprocs];
  MPI_Alltoall (temp_info, 3, MPI_INT, recv_info, 3, MPI_INT, MPI_COMM_WORLD);
  delete []temp_info;

  /***********************************/
  /*  first post all of the receives */
  /***********************************/
  MPI_Request *recv_request = new MPI_Request[3 * numprocs];

  // recv_count is total data current proc will receive
  int recv_count[3] = { 0, 0, 0 };
  for (i = 0; i < numprocs; i++)
  {
    recv_count[0] += recv_info[3 * i];
    recv_count[1] += recv_info[3 * i + 1];
    recv_count[2] += recv_info[3 * i + 2];
  }

  int counter_recv[3] = { 0, 0, 0 };
  int blksize = 2 * KEYLENGTH + 1;      // 2 keys + 1 proc
  unsigned *recv_neigh_array = new unsigned[recv_count[0] * blksize];
  BucketPack *recv_buck_array = new BucketPack[recv_count[1]];
  ParticlePack *recv_part_array = new ParticlePack[recv_count[2]];
  int neigh_tag = 99565, buck_tag = 88476, part_tag = 35262;    // random tag numbers

  for (i = 0; i < numprocs; i++)
  {
    if (recv_info[3 * i] != 0)  // receive neighbor info
    {
      j =
        MPI_Irecv ((recv_neigh_array + counter_recv[0]),
                   recv_info[3 * i] * blksize, MPI_UNSIGNED, i, neigh_tag,
                   MPI_COMM_WORLD, (recv_request + 3 * i));
      counter_recv[0] += recv_info[3 * i] * blksize;
    }
    if (recv_info[3 * i + 1] != 0)      // receive buckets
    {
      j =
        MPI_Irecv ((recv_buck_array + counter_recv[1]),
                   recv_info[3 * i + 1], BUCKET_TYPE, i, buck_tag,
                   MPI_COMM_WORLD, (recv_request + 3 * i + 1));
      counter_recv[1] += recv_info[3 * i + 1];
    }
    if (recv_info[3 * i + 2] != 0)      // receive particles
    {
      j =
        MPI_Irecv ((recv_part_array + counter_recv[2]),
                   recv_info[3 * i + 2], PARTICLE_TYPE, i, part_tag,
                   MPI_COMM_WORLD, (recv_request + 3 * i + 2));
      counter_recv[2] += recv_info[3 * i + 2];
    }
  } /* done posting the receives */

  int send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[3 * i];

  //While sending out bucket to other processes,
  //mykey is send out together with the key of the bucket that is being sent
  unsigned *send_neigh_array = new unsigned[send_count * blksize];
  for (i = 0; i < send_count * blksize; i++)
    send_neigh_array[i] = 5;

  int *counter_send_proc = new int[numprocs];
  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] = counter_send_proc[i - 1] + send_info[3 * (i - 1)];

  // pack neighbor information to send out
  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ()) // this element will get moved to a new processor
		  {
		      mykey = buck->getKey ();
		      const int * neigh_proc = buck->get_neigh_proc ();
		      Key * neighbors = buck->get_neighbors ();

		      for (j = 0; j < NEIGH_SIZE; j++)  //go through all neighbors one by one.
		        if ((neigh_proc[j] != myid) && (neigh_proc[j] > -1))
		        {
		          for (k = 0; k < KEYLENGTH; k++)
		          {
		            send_neigh_array[counter_send_proc[neigh_proc[j]] *
		                             blksize + k] = neighbors[j].key[k];
		            send_neigh_array[counter_send_proc[neigh_proc[j]] *
		                             blksize + KEYLENGTH + k] = mykey.key[k];
		          }
		          send_neigh_array[counter_send_proc[neigh_proc[j]] * blksize +
		                           2 * KEYLENGTH] = (unsigned) buck->get_myprocess ();
		          counter_send_proc[neigh_proc[j]] += 1;
		        }
		  }
	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // send out neighbor information
  int counter = 0;
  MPI_Request *send_request = new MPI_Request[numprocs];

  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i] != 0)
    {
      j = MPI_Isend ((send_neigh_array + counter * blksize),
                     send_info[3 * i] * blksize, MPI_UNSIGNED, i,
                     neigh_tag, MPI_COMM_WORLD, (send_request + i));
      counter += send_info[3 * i];
    }

  // update neighbor information on this processor
  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  // this element will get moved to a new processor
		  if (myid != buck->get_myprocess ())
		  {
		      // neigh info
		      for (i = 0; i < 3; i++)
		        for (j = 0; j < 3; j++)
		          for (k = 0; k < 3; k++)
		          {
		            // skip thyself
		            if ((i == 0) && (j == 0) && (k == 0))
		              continue;

		            // direction indices of neighbors
		            int dir1[DIMENSION] = { i, j, k };

		            // mirror the direction indices
		            // indices can be 0, 1 or 2
		            //  1^3 = 2  (bitwise XOR)
		            //  2^3 = 1
		            dir2[0] = (i > 0 ? 3 ^ i : i);
		            dir2[1] = (j > 0 ? 3 ^ j : j);
		            dir2[2] = (k > 0 ? 3 ^ k : k);

		            // update neigh-procs
		            int neigh_proc = buck->which_neigh_proc (dir1);
		            if ( neigh_proc == myid )
		            {
		              Bucket *BkNeigh =
		                  (Bucket *) BG_mesh->lookup (buck->which_neigh (dir1));
		              buck->put_neigh_proc (dir1, BkNeigh->get_myprocess ());
		              BkNeigh->put_neigh_proc (dir2, buck->get_myprocess ());
		            }
		          }
		   }
	  }// end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // wait to receive new neighbor info and then update neighbor info
  // from neighbor information received from other processors
  MPI_Status status;

  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[3 * i] != 0)
    {
      j = MPI_Wait ((recv_request + 3 * i), &status);
      for (j = counter; j < counter + recv_info[3 * i]; j++)
      {
        for (k = 0; k < KEYLENGTH; k++)
        {
          mykey.key[k] = recv_neigh_array[j * blksize + k];
          neighkey.key[k] = recv_neigh_array[j * blksize + KEYLENGTH + k];
        }
        int neigh_proc = (int) recv_neigh_array[j * blksize + 2 * KEYLENGTH];

        buck = (Bucket *) BG_mesh->lookup (mykey);
        assert (buck);

        if (buck->find_neigh_dir (neighkey, dir1))
        {
          buck->put_neigh_proc (dir1, neigh_proc);
        }
        else
        {
          fprintf (stderr, "ERROR: can't find neigh_dir in %s, at %d\n",
                   __FILE__, __LINE__);
          MPI_Abort (MPI_COMM_WORLD, -1);
        }
      }
      counter += recv_info[3 * i];
    }
  delete []recv_neigh_array;

  // wait and then delete sent info that was already received
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_neigh_array;

/****************************************************
 *     migrate particles within the buckets
 ****************************************************/
  send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[3 * i + 2];

  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] =
      counter_send_proc[i - 1] + send_info[3 * (i - 1) + 2];

  // pack particles
  ParticlePack *send_part_array = new ParticlePack[send_count];

  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ())
		  {
		      int myprocess = buck->get_myprocess ();

		      vector < TKey > plist = buck->get_plist ();
		      if (plist.size () > 0)
		      {
		        vector < TKey >::iterator ip;
		        for (ip = plist.begin (); ip != plist.end (); ip++)
		        {
		          Particle *psend = (Particle *) P_table->lookup (*ip);

		          pack_particles (psend,
		                          (send_part_array + counter_send_proc[myprocess]));
		          counter_send_proc[myprocess]++;
		        }
		      }
		    }
	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // now send out packed particles arrays
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 2] != 0)
    {
      j =
        MPI_Isend ((send_part_array + counter), send_info[3 * i + 2],
                   PARTICLE_TYPE, i, part_tag, MPI_COMM_WORLD,
                   (send_request + i));
      counter += send_info[3 * i + 2];
    }

  // Wait for recieves and unpack particles
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[3 * i + 2] != 0)
    {
      j = MPI_Wait ((recv_request + 3 * i + 2), &status);
      for (j = 0; j < recv_info[3 * i + 2]; j++)
      {
        // unpack particles here
        Particle *part = new Particle ();

        unpack_particle ((recv_part_array + counter), part);
        P_table->add (part->getKey (), part);
        counter++;
      }
    }

  delete []recv_part_array;

  // wait and delete send info
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 2] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_part_array;

  ////////////////////////////////////////////////////////////
  // can now migrate any of the buckets which need to be moved
  ////////////////////////////////////////////////////////////

  // first pack buckets in an array
  send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[3 * i + 1];

  BucketPack *send_buck_array = new BucketPack[send_count];

  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] =
      counter_send_proc[i - 1] + send_info[3 * (i - 1) + 1];

  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ())
		  {
		      int myprocess = buck->get_myprocess ();

		      assert ((myprocess > -1) && (myprocess < numprocs));
		      pack_bucket ((send_buck_array + counter_send_proc[myprocess]), buck,
		                   myprocess);
		      counter_send_proc[myprocess] += 1;

		      vector < TKey > plist = buck->get_plist ();
		      for (j = 0; j < (int) plist.size (); j++)
		      {
		        Particle *p2rm = (Particle *) P_table->lookup (plist[j]);

		        P_table->remove (plist[j]);
		        delete p2rm;
		      }
		      BG_mesh->remove (buck->getKey ());
		      delete buck;
		  }
	  }//end of if bucket is not a brief bucket
  }//end of while loop go through all buckets

  delete []counter_send_proc;

  // now send out packed elements
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 1] != 0)
    {
      j =
        MPI_Isend ((send_buck_array + counter), send_info[3 * i + 1],
                   BUCKET_TYPE, i, buck_tag, MPI_COMM_WORLD,
                   (send_request + i));
      counter += send_info[3 * i + 1];
    }

  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[3 * i + 1] != 0)
    {
      j = MPI_Wait ((recv_request + 3 * i + 1), &status);
      for (j = 0; j < recv_info[3 * i + 1]; j++)
      {
        // upack bucket data
        Bucket *buck = new Bucket ();

        unpack_bucket ((recv_buck_array + counter), buck, myid);
        BG_mesh->add (buck->getKey (), buck);
        counter++;
      }
    }

#ifdef DEBUG2
  fclose (fp);
#endif

  delete []recv_buck_array;
  delete []recv_request;
  delete []recv_info;

  // wait for send info before deleting
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 1] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_buck_array;
  delete []send_request;
  delete []send_info;
  delete itr;

  return;
}

//This function will only move brief buckets
//I will have another function for communication of only brief buckets.
void
BSFC_update_and_send_elements (int myid, int numprocs,
                               THashTable * P_table, HashTable * BG_mesh)
{
  int i, j, k;
  int dir1[DIMENSION], dir2[DIMENSION];
  Key mykey, neighkey;

#ifdef DEBUG2
  char filename[20];
  sprintf (filename, "debug%02d.txt", myid);
  FILE * fp = fopen (filename, "a+");
#endif

  int *send_info = new int[numprocs * 3];
  for (i = 0; i < 3 * numprocs; i++)
    send_info[i] = 0;

  // now figure out what processors need to have neighbor info updated
  // and what processors to send objects to
  Bucket *buck;
  BriefBucket *breif_buck = NULL;
  void * tempptr =NULL;

  HTIterator *itr = new HTIterator (BG_mesh);

  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ()) // this element will get moved to a new processor
		  {
		     // neigh info
		     const int *neigh_proc = buck->get_neigh_proc ();

		     for (j = 0; j < NEIGH_SIZE; j++)
		        if ((neigh_proc[j] != myid) && (neigh_proc[j] > -1))
		          send_info[3 * neigh_proc[j]] += 1;

		     // bucket info
		     send_info[3 * (buck->get_myprocess ()) + 1] += 1;
		     send_info[3 * (buck->get_myprocess ()) + 2] +=
		        (int) buck->get_plist ().size ();
		  }

	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // copy send_info to temp_info
  int *temp_info = new int[3 * numprocs];
  for (i = 0; i < 3 * numprocs; i++)
    temp_info[i] = send_info[i];

  // send info to everyone, so everyone know how much
  // they are receiving from whom
  int *recv_info = new int[3 * numprocs];
  MPI_Alltoall (temp_info, 3, MPI_INT, recv_info, 3, MPI_INT, MPI_COMM_WORLD);
  delete []temp_info;

  /***********************************/
  /*  first post all of the receives */
  /***********************************/
  MPI_Request *recv_request = new MPI_Request[3 * numprocs];

  // recv_count is total data current proc will receive
  int recv_count[3] = { 0, 0, 0 };
  for (i = 0; i < numprocs; i++)
  {
    recv_count[0] += recv_info[3 * i];
    recv_count[1] += recv_info[3 * i + 1];
    recv_count[2] += recv_info[3 * i + 2];
  }

  int counter_recv[3] = { 0, 0, 0 };
  int blksize = 2 * KEYLENGTH + 1;      // 2 keys + 1 proc
  unsigned *recv_neigh_array = new unsigned[recv_count[0] * blksize];
  BucketPack *recv_buck_array = new BucketPack[recv_count[1]];
  ParticlePack *recv_part_array = new ParticlePack[recv_count[2]];
  int neigh_tag = 99565, buck_tag = 88476, part_tag = 35262;    // random tag numbers

  for (i = 0; i < numprocs; i++)
  {
    if (recv_info[3 * i] != 0)  // receive neighbor info
    {
      j =
        MPI_Irecv ((recv_neigh_array + counter_recv[0]),
                   recv_info[3 * i] * blksize, MPI_UNSIGNED, i, neigh_tag,
                   MPI_COMM_WORLD, (recv_request + 3 * i));
      counter_recv[0] += recv_info[3 * i] * blksize;
    }
    if (recv_info[3 * i + 1] != 0)      // receive buckets
    {
      j =
        MPI_Irecv ((recv_buck_array + counter_recv[1]),
                   recv_info[3 * i + 1], BUCKET_TYPE, i, buck_tag,
                   MPI_COMM_WORLD, (recv_request + 3 * i + 1));
      counter_recv[1] += recv_info[3 * i + 1];
    }
    if (recv_info[3 * i + 2] != 0)      // receive particles
    {
      j =
        MPI_Irecv ((recv_part_array + counter_recv[2]),
                   recv_info[3 * i + 2], PARTICLE_TYPE, i, part_tag,
                   MPI_COMM_WORLD, (recv_request + 3 * i + 2));
      counter_recv[2] += recv_info[3 * i + 2];
    }
  } /* done posting the receives */

  int send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[3 * i];

  //While sending out bucket to other processes,
  //mykey is send out together with the key of the bucket that is being sent
  unsigned *send_neigh_array = new unsigned[send_count * blksize];
  for (i = 0; i < send_count * blksize; i++)
    send_neigh_array[i] = 5;

  int *counter_send_proc = new int[numprocs];
  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] = counter_send_proc[i - 1] + send_info[3 * (i - 1)];

  // pack neighbor information to send out
  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ()) // this element will get moved to a new processor
		  {
		      mykey = buck->getKey ();
		      const int * neigh_proc = buck->get_neigh_proc ();
		      Key * neighbors = buck->get_neighbors ();

		      for (j = 0; j < NEIGH_SIZE; j++)  //go through all neighbors one by one.
		        if ((neigh_proc[j] != myid) && (neigh_proc[j] > -1))
		        {
		          for (k = 0; k < KEYLENGTH; k++)
		          {
		            send_neigh_array[counter_send_proc[neigh_proc[j]] *
		                             blksize + k] = neighbors[j].key[k];
		            send_neigh_array[counter_send_proc[neigh_proc[j]] *
		                             blksize + KEYLENGTH + k] = mykey.key[k];
		          }
		          send_neigh_array[counter_send_proc[neigh_proc[j]] * blksize +
		                           2 * KEYLENGTH] = (unsigned) buck->get_myprocess ();
		          counter_send_proc[neigh_proc[j]] += 1;
		        }
		  }
	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // send out neighbor information
  int counter = 0;
  MPI_Request *send_request = new MPI_Request[numprocs];

  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i] != 0)
    {
      j = MPI_Isend ((send_neigh_array + counter * blksize),
                     send_info[3 * i] * blksize, MPI_UNSIGNED, i,
                     neigh_tag, MPI_COMM_WORLD, (send_request + i));
      counter += send_info[3 * i];
    }

  // update neighbor information on this processor
  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  // this element will get moved to a new processor
		  if (myid != buck->get_myprocess ())
		  {
		      // neigh info
		      for (i = 0; i < 3; i++)
		        for (j = 0; j < 3; j++)
		          for (k = 0; k < 3; k++)
		          {
		            // skip thyself
		            if ((i == 0) && (j == 0) && (k == 0))
		              continue;

		            // direction indices of neighbors
		            int dir1[DIMENSION] = { i, j, k };

		            // mirror the direction indices
		            // indices can be 0, 1 or 2
		            //  1^3 = 2  (bitwise XOR)
		            //  2^3 = 1
		            dir2[0] = (i > 0 ? 3 ^ i : i);
		            dir2[1] = (j > 0 ? 3 ^ j : j);
		            dir2[2] = (k > 0 ? 3 ^ k : k);

		            // update neigh-procs
		            int neigh_proc = buck->which_neigh_proc (dir1);
		            if ( neigh_proc == myid )
		            {
		              Bucket *BkNeigh =
		                  (Bucket *) BG_mesh->lookup (buck->which_neigh (dir1));
		              buck->put_neigh_proc (dir1, BkNeigh->get_myprocess ());
		              BkNeigh->put_neigh_proc (dir2, buck->get_myprocess ());
		            }
		          }
		   }
	  }// end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // wait to receive new neighbor info and then update neighbor info
  // from neighbor information received from other processors
  MPI_Status status;

  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[3 * i] != 0)
    {
      j = MPI_Wait ((recv_request + 3 * i), &status);
      for (j = counter; j < counter + recv_info[3 * i]; j++)
      {
        for (k = 0; k < KEYLENGTH; k++)
        {
          mykey.key[k] = recv_neigh_array[j * blksize + k];
          neighkey.key[k] = recv_neigh_array[j * blksize + KEYLENGTH + k];
        }
        int neigh_proc = (int) recv_neigh_array[j * blksize + 2 * KEYLENGTH];

        buck = (Bucket *) BG_mesh->lookup (mykey);
        assert (buck);

        if (buck->find_neigh_dir (neighkey, dir1))
        {
          buck->put_neigh_proc (dir1, neigh_proc);
        }
        else
        {
          fprintf (stderr, "ERROR: can't find neigh_dir in %s, at %d\n",
                   __FILE__, __LINE__);
          MPI_Abort (MPI_COMM_WORLD, -1);
        }
      }
      counter += recv_info[3 * i];
    }
  delete []recv_neigh_array;

  // wait and then delete sent info that was already received
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_neigh_array;

/****************************************************
 *     migrate particles within the buckets
 ****************************************************/
  send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[3 * i + 2];

  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] =
      counter_send_proc[i - 1] + send_info[3 * (i - 1) + 2];

  // pack particles
  ParticlePack *send_part_array = new ParticlePack[send_count];

  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ())
		  {
		      int myprocess = buck->get_myprocess ();

		      vector < TKey > plist = buck->get_plist ();
		      if (plist.size () > 0)
		      {
		        vector < TKey >::iterator ip;
		        for (ip = plist.begin (); ip != plist.end (); ip++)
		        {
		          Particle *psend = (Particle *) P_table->lookup (*ip);

		          pack_particles (psend,
		                          (send_part_array + counter_send_proc[myprocess]));
		          counter_send_proc[myprocess]++;
		        }
		      }
		    }
	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // now send out packed particles arrays
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 2] != 0)
    {
      j =
        MPI_Isend ((send_part_array + counter), send_info[3 * i + 2],
                   PARTICLE_TYPE, i, part_tag, MPI_COMM_WORLD,
                   (send_request + i));
      counter += send_info[3 * i + 2];
    }

  // Wait for recieves and unpack particles
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[3 * i + 2] != 0)
    {
      j = MPI_Wait ((recv_request + 3 * i + 2), &status);
      for (j = 0; j < recv_info[3 * i + 2]; j++)
      {
        // unpack particles here
        Particle *part = new Particle ();

        unpack_particle ((recv_part_array + counter), part);
        P_table->add (part->getKey (), part);
        counter++;
      }
    }

  delete []recv_part_array;

  // wait and delete send info
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 2] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_part_array;

  ////////////////////////////////////////////////////////////
  // can now migrate any of the buckets which need to be moved
  ////////////////////////////////////////////////////////////

  // first pack buckets in an array
  send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[3 * i + 1];

  BucketPack *send_buck_array = new BucketPack[send_count];

  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] =
      counter_send_proc[i - 1] + send_info[3 * (i - 1) + 1];

  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ())
		  {
		      int myprocess = buck->get_myprocess ();

		      assert ((myprocess > -1) && (myprocess < numprocs));
		      pack_bucket ((send_buck_array + counter_send_proc[myprocess]), buck,
		                   myprocess);
		      counter_send_proc[myprocess] += 1;

		      vector < TKey > plist = buck->get_plist ();
		      for (j = 0; j < (int) plist.size (); j++)
		      {
		        Particle *p2rm = (Particle *) P_table->lookup (plist[j]);

		        P_table->remove (plist[j]);
		        delete p2rm;
		      }
		      BG_mesh->remove (buck->getKey ());
		      delete buck;
		  }
	  }//end of if bucket is not a brief bucket
  }//end of while loop go through all buckets

  delete []counter_send_proc;

  // now send out packed elements
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 1] != 0)
    {
      j =
        MPI_Isend ((send_buck_array + counter), send_info[3 * i + 1],
                   BUCKET_TYPE, i, buck_tag, MPI_COMM_WORLD,
                   (send_request + i));
      counter += send_info[3 * i + 1];
    }

  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[3 * i + 1] != 0)
    {
      j = MPI_Wait ((recv_request + 3 * i + 1), &status);
      for (j = 0; j < recv_info[3 * i + 1]; j++)
      {
        // upack bucket data
        Bucket *buck = new Bucket ();

        unpack_bucket ((recv_buck_array + counter), buck, myid);
        BG_mesh->add (buck->getKey (), buck);
        counter++;
      }
    }

#ifdef DEBUG2
  fclose (fp);
#endif

  delete []recv_buck_array;
  delete []recv_request;
  delete []recv_info;

  // wait for send info before deleting
  for (i = 0; i < numprocs; i++)
    if (send_info[3 * i + 1] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_buck_array;
  delete []send_request;
  delete []send_info;
  delete itr;

  return;
}
