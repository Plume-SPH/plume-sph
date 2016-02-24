
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
#include <hilbert.h>

#include "repartition_BSFC.h"
#include "exvar.h"
#include "pack_data.h"
#include "multiproc.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

//The following is exactly the same as the one defined in bucket.cc
//To avoid name conflict,
int ineigh_new[3][3][3] = {{{13, 12, 14}, {10,  9, 11}, {16, 15, 17}},
                       {{4,   3,  5}, {1,   0,  2}, {7,   6,  8}},
                       {{22, 21, 23}, {19, 18, 20}, {25, 24, 26}}
                      };
//function that will determine bucket key based on bucket coordinate --> Here is the expense that we pay for using of brief bucket to minimized memory allocation.
void compute_key(double* mindom, double * maxdom, double *crd, unsigned key[])
{
	double normc[DIMENSION];
	for (int l=0; l<DIMENSION; l++)
	    normc[l]=(*(crd+l)- *(mindom+l))/(*(maxdom+l)- *(mindom+l));

	unsigned keylen = KEYLENGTH;

	HSFC3d (normc, & keylen, key);

	return;
}

//function that used to determine neighbor info of one
void compute_neighbor_keys(double*  mindom, double*  maxdom, double*  min, double bucket_size, Key * neighbors)
{
	int ii, jj, kk;
	int i, j;
	double bksize_half = 0.5 * bucket_size;
	double new_dom_max[DIMENSION],new_dom_min[DIMENSION], crd[DIMENSION];
	double normc[DIMENSION];

    unsigned keylen = KEYLENGTH;
	Key key;
	unsigned u_key [KEYLENGTH];

    for (i=0;i<DIMENSION; i++)
    {
    	new_dom_max[i]= *(maxdom+i) - bksize_half;
    	new_dom_min[i]= *(mindom+i) + bksize_half;
    }

	j=0;
	for (ii=0; ii<DIMENSION; ii++)//The loop can guarantee the order of neighbor keys are consistent
	{
		 crd[0]=*min + ii*bucket_size-bksize_half;
		 for (jj=0; jj<DIMENSION; jj++)
		 {
			 crd[1]=*(min+1)+jj*bucket_size-bksize_half;
			 for (kk=0; kk<DIMENSION; kk++)
			 {
				 crd[2]=*(min+2)+kk*bucket_size-bksize_half;
			         for ( int i=0; i<DIMENSION; i++)
			         	  normc[i]=(crd[i]- *(mindom+i))/(*(maxdom+i) - *(mindom+i));

			         HSFC3d (normc, & keylen, u_key);

		         for (i=0; i<KEYLENGTH; i++)
		         	  key.key[i]=u_key[i];

		         *(neighbors+j)=key;

		         j++;

			 }
		 }
	}

	return;
}

//Some comments about this function:
/*
 * Several factors should be considered while design this data movement function
 * 1) Overlap communication and computation --> To minimize communication overhead
 * 2) Usage of memory: Memory should be released as soon as possible
 *                     Avoid allocate too much space for intermediate variable
 *                     The reason why you need to pay attention to memory usage is not only due to out of memory danger
 *                     Allocation of memory might also affect the efficiency of memory access
 * 3) Pay attention to the order of sending information, who comes first, who comes second
 *                                                       Sometimes, Data A should be received before going to the next step to send data B
 */

//Based on Dinesh's code, Zhixuan made the following modification:
/*
 * 1) Separate the migration of brief bucket and non-brief_bucket
 * 2) Add migration of brief bucket
 * 3) Add code for updating neigh info, if the the neigh is a brief bucket
 * Note: Brief bucket does not have neigh info, neigh info will be added at the time when switch brief bucket to bucket
 *       Neigh proc info of brief on a non-brief bucket should be correctly communicated
 */
void
BSFC_update_and_send_elements (int myid, int numprocs, double bucket_size,
                               THashTable * P_table, HashTable * BG_mesh)
{
  int i, j, k;
  int ii, jj, kk;
  int dir1[DIMENSION], dir2[DIMENSION], dir3[DIMENSION];
  int info_num = 4; //number of info that need to communicate: bucket, brief bucket, neigh. particle
  Key mykey, neighkey;
  unsigned neigh_key[KEYLENGTH];

#ifdef DEBUG2
  char filename[20];
  sprintf (filename, "debug%02d.txt", myid);
  FILE * fp = fopen (filename, "a+");
#endif

#ifdef DEBUG
   bool do_check = false;
#endif

  int *send_info = new int[numprocs * info_num]; //Bucket is categorized into two types: brief bucket and bucket: 4: particle, neigh, bucket, brief bucket
  for (i = 0; i < info_num * numprocs; i++)
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
	  {
		  if (myid != breif_buck->get_myprocess ()) // this element will get moved to a new processor
		  {
			// neigh info
			const int *neigh_proc = breif_buck->get_neigh_proc ();

			for (j = 0; j < NEIGH_SIZE; j++)
			   if ((neigh_proc[j] != myid) && (neigh_proc[j] > -1))
			       send_info[info_num * neigh_proc[j]] += 1;

		     // brief buck info
		     send_info[info_num * (breif_buck->get_myprocess ()) + 3] += 1;
		  }
	  }
	  else
	  {
		  buck = (Bucket*) tempptr;
		  if (myid != buck->get_myprocess ()) // this element will get moved to a new processor
		  {
		     // neigh info
		     const int *neigh_proc = buck->get_neigh_proc ();

		     for (j = 0; j < NEIGH_SIZE; j++)
		        if ((neigh_proc[j] != myid) && (neigh_proc[j] > -1))
		          send_info[info_num * neigh_proc[j]] += 1;

		     // bucket info
		     send_info[info_num * (buck->get_myprocess ()) + 1] += 1;
		     send_info[info_num * (buck->get_myprocess ()) + 2] +=
		        (int) buck->get_plist ().size ();
		  }

	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // copy send_info to temp_info
  int *temp_info = new int[info_num * numprocs];
  for (i = 0; i < info_num * numprocs; i++)
    temp_info[i] = send_info[i];

  // send info to everyone, so everyone know how much
  // they are receiving from whom
  int *recv_info = new int[info_num * numprocs];
  MPI_Alltoall (temp_info, info_num, MPI_INT, recv_info, info_num, MPI_INT, MPI_COMM_WORLD);
  delete []temp_info;

  /***********************************/
  /*  first post all of the receives */
  /***********************************/
  MPI_Request *recv_request = new MPI_Request[info_num * numprocs];

  // recv_count is total data current proc will receive
  int recv_count[4] = {0, 0, 0, 0};
  for (i = 0; i < numprocs; i++)
  {
    recv_count[0] += recv_info[info_num * i];
    recv_count[1] += recv_info[info_num * i + 1];
    recv_count[2] += recv_info[info_num * i + 2];
    recv_count[3] += recv_info[info_num * i + 3];
  }

  int counter_recv[4] = { 0, 0, 0, 0 };
  int blksize = 2 * KEYLENGTH + 1;      // 2 keys + 1 proc
  unsigned *recv_neigh_array = new unsigned[recv_count[0] * blksize];
  BucketPack *recv_buck_array = new BucketPack[recv_count[1]];
  ParticlePack *recv_part_array = new ParticlePack[recv_count[2]];
  BriefBucketPack *recv_briefbuck_array = new BriefBucketPack[recv_count[3]];
  int neigh_tag = 99565, buck_tag = 88476, part_tag = 35262, briefbuck_tag=12347;    // random tag numbers

  for (i = 0; i < numprocs; i++)
  {
    if (recv_info[info_num * i] != 0)  // receive neighbor info
    {
      j =
        MPI_Irecv ((recv_neigh_array + counter_recv[0]),
                   recv_info[info_num * i] * blksize, MPI_UNSIGNED, i, neigh_tag,
                   MPI_COMM_WORLD, (recv_request + info_num * i));
      counter_recv[0] += recv_info[info_num * i] * blksize;
    }
    if (recv_info[info_num * i + 1] != 0)      // receive buckets
    {
      j =
        MPI_Irecv ((recv_buck_array + counter_recv[1]),
                   recv_info[info_num * i + 1], BUCKET_TYPE, i, buck_tag,
                   MPI_COMM_WORLD, (recv_request + info_num * i + 1));
      counter_recv[1] += recv_info[info_num * i + 1];
    }
    if (recv_info[info_num * i + 2] != 0)      // receive particles
    {
      j =
        MPI_Irecv ((recv_part_array + counter_recv[2]),
                   recv_info[info_num * i + 2], PARTICLE_TYPE, i, part_tag,
                   MPI_COMM_WORLD, (recv_request + info_num * i + 2));
      counter_recv[2] += recv_info[info_num * i + 2];
    }
    if (recv_info[info_num * i + 3] != 0)      // receive brief buckets
    {
      j =
        MPI_Irecv ((recv_briefbuck_array + counter_recv[3]),
                   recv_info[info_num * i + 3], BRIEF_BUCKET_TYPE, i, briefbuck_tag,
                   MPI_COMM_WORLD, (recv_request + info_num * i + 3));
      counter_recv[3] += recv_info[info_num * i + 3];
    }
  } /* done posting the receives */

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

  int send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[info_num * i]; //send count of neighbor info

  //While sending out bucket to other processes,
  //mykey is send out together with the key of the bucket that is being sent
  unsigned *send_neigh_array = new unsigned[send_count * blksize];
  for (i = 0; i < send_count * blksize; i++)
    send_neigh_array[i] = 5; //5=blksize

  int *counter_send_proc = new int[numprocs];
  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] = counter_send_proc[i - 1] + send_info[info_num * (i - 1)];

  // pack neighbor information to send out --> The purpose is essentially to update neigh_proc info
  double mindom[DIMENSION], maxdom[DIMENSION];
  for (j = 0; j < DIMENSION; j++)
  		mindom[j]= *(BG_mesh->get_minDom () + j);
  for (j = 0; j < DIMENSION; j++)
  	  	maxdom[j]= *(BG_mesh->get_maxDom () + j);
  double crd[DIMENSION], min[DIMENSION];
  double bksize_half = 0.5*bucket_size;

  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
	  {
		  if (myid != breif_buck->get_myprocess ()) // this element will get moved to a new processor
		  {
		      mykey = breif_buck->getKey ();
		      const int * neigh_proc = breif_buck->get_neigh_proc ();
//		      Key * neighbors = compute_neighbor(mindom, maxdom);
		      min[0]=*breif_buck->get_mincrd();
		      min[1]=*(breif_buck->get_mincrd()+1);
		      min[2]=*(breif_buck->get_mincrd()+2);
              //It really matters that the order of neighbor_proc should be consistent with the order of neighbors
		      // --->which is given in ineigh
	    	  for (ii=0; ii<DIMENSION; ii++)
	    	  {
	    		  crd[0]=min[0]+ii*bucket_size-bksize_half;
	    		  dir3[0] = (ii==1 ? 0 : (ii == 0 ? 1 : 2));
	    		  for (jj=0; jj<DIMENSION; jj++)
	    		  {
	    			  crd[1]=min[1]+jj*bucket_size-bksize_half;
	    			  dir3[1] = (jj==1 ? 0 : (jj == 0 ? 1 : 2));
	    			  for (kk=0; kk<DIMENSION; kk++)
	    			  {
	    				  dir3[2] = (kk==1 ? 0 : (kk == 0 ? 1 : 2));
	    				  j=ineigh_new[dir3[0]][dir3[1]][dir3[2]];
	    			      if ((neigh_proc[j] != myid) && (neigh_proc[j] > -1))  //This will guarantee that the neighbor is not invalid buckets (such as bucket that locate outside of the domain)
	    			      {
	    			    	  crd[2]=min[2]+kk*bucket_size-bksize_half;
	    			    	  compute_key(mindom, maxdom, crd, neigh_key);
	    			          for (k = 0; k < KEYLENGTH; k++)
	    			          {
	    			            send_neigh_array[counter_send_proc[neigh_proc[j]] *
	    			                             blksize + k] = neigh_key[k];
	    			            send_neigh_array[counter_send_proc[neigh_proc[j]] *
	    			                             blksize + KEYLENGTH + k] = mykey.key[k];
	    			          }
	    			          send_neigh_array[counter_send_proc[neigh_proc[j]] * blksize +
	    			                           2 * KEYLENGTH] = (unsigned) breif_buck->get_myprocess ();
	    			          counter_send_proc[neigh_proc[j]] += 1;
	    			      }
	    			  }
	    		  }
	    	  }

		  }//end of if this bucket need to be moved to other processes
	  }
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

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

  // send out neighbor information
  int counter = 0;
  MPI_Request *send_request = new MPI_Request[numprocs];

  for (i = 0; i < numprocs; i++)
    if (send_info[info_num * i] != 0)
    {
      j = MPI_Isend ((send_neigh_array + counter * blksize),
                     send_info[info_num * i] * blksize, MPI_UNSIGNED, i,
                     neigh_tag, MPI_COMM_WORLD, (send_request + i));
      counter += send_info[info_num * i];
    }

  // update neighbor information on this processor---> while communication is going on
  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, the neighbor info  is not available --> and not necessary to update it here!
	  {
		  // this element will get moved to a new processor
		  if (myid != breif_buck->get_myprocess ())
		  {
		      min[0]=*breif_buck->get_mincrd();
		      min[1]=*(breif_buck->get_mincrd()+1);
		      min[2]=*(breif_buck->get_mincrd()+2);
		      // neigh info
		      for (i = 0; i < DIMENSION; i++)
		      {
		    	crd[0]=min[0]+i*bucket_size-bksize_half;
		    	dir3[0] = (i==1 ? 0 : (i == 0 ? 1 : 2));
		        for (j = 0; j < DIMENSION; j++)
		        {
		          crd[1]=min[1]+j*bucket_size-bksize_half;
		          dir3[1] = (j==1 ? 0 : (j == 0 ? 1 : 2));
		          for (k = 0; k < DIMENSION; k++)
		          {
		        	dir3[2] = (k==1 ? 0 : (k == 0 ? 1 : 2));
		            // skip itself
		            if ((dir3[0] == 0) && (dir3[1] == 0) && (dir3[2] == 0))
		              continue;

		            // direction indices of neighbors
		            int dir1[DIMENSION] = { dir3[0], dir3[1], dir3[2] };

		            // mirror the direction indices
		            // indices can be 0, 1 or 2
		            //  1^3 = 2  (bitwise XOR)
		            //  2^3 = 1
		            dir2[0] = (dir3[0] > 0 ? 3 ^ dir3[0] : dir3[0]);
		            dir2[1] = (dir3[1] > 0 ? 3 ^ dir3[1] : dir3[1]);
		            dir2[2] = (dir3[2] > 0 ? 3 ^ dir3[2] : dir3[2]);

		            // update neigh-procs
		            int neigh_proc = breif_buck->which_neigh_proc (dir1);
		            if ( neigh_proc == myid ) //which means that this neigh bucket will stay on its old owner process, this bucket should be available on current process
		            {
		              crd[2]=min[2]+k*bucket_size-bksize_half;
		              compute_key(mindom, maxdom, crd, neigh_key);
		              BriefBucket *BkNeigh =
		                  (BriefBucket *) BG_mesh->lookup (neigh_key);//We do not distinguish brief bucket and bucket here, hopefully this will work;
		              breif_buck->put_neigh_proc (dir1, BkNeigh->get_myprocess ());
		              BkNeigh->put_neigh_proc (dir2, breif_buck->get_myprocess ());
		            }
		          }
		        }
		      }
		   }
	  }//end of if bucket is brief bucket
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
		            if ((i == 0) && (j == 0) && (k == 0)) //Does this means the bucket at the center?
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
		              BriefBucket *BkNeigh =
		            		(BriefBucket *) BG_mesh->lookup (buck->which_neigh (dir1)); //We do not distinguish brief bucket and bucket here, hopefully this will work;
		              buck->put_neigh_proc (dir1, BkNeigh->get_myprocess ());
		              BkNeigh->put_neigh_proc (dir2, buck->get_myprocess ());
		            }
		          }
		   }
	  }// end of if bucket is not brief bucket
  }//end of while loop go through all buckets

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

  // wait to receive new neighbor info and then update neighbor info
  // from neighbor information received from other processors
  MPI_Status status;

  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[info_num * i] != 0)
    {
      j = MPI_Wait ((recv_request + info_num * i), &status); //Make sure data is received before start making use of it
      for (j = counter; j < counter + recv_info[info_num * i]; j++)
      {
        for (k = 0; k < KEYLENGTH; k++)
        {
          mykey.key[k] = recv_neigh_array[j * blksize + k];
          neighkey.key[k] = recv_neigh_array[j * blksize + KEYLENGTH + k];
        }
        int neigh_proc = (int) recv_neigh_array[j * blksize + 2 * KEYLENGTH];

        tempptr	=  BG_mesh->lookup (mykey);
        assert(tempptr);
        breif_buck = (BriefBucket *) tempptr;
        if (breif_buck->check_brief()) //if is brief bucket, the neighbor is not necessary to update here
        {
		    min[0]=*breif_buck->get_mincrd();
		    min[1]=*(breif_buck->get_mincrd()+1);
		    min[2]=*(breif_buck->get_mincrd()+2);

		    Key neighbors[NEIGH_SIZE];
            compute_neighbor_keys(mindom, maxdom, min, bucket_size, neighbors);
        	if (breif_buck->find_neigh_dir (neighkey, neighbors, dir1))
            {
        		breif_buck->put_neigh_proc (dir1, neigh_proc);
            }
            else
            {
              fprintf (stderr, "ERROR: can't find neigh_dir in %s, at %d\n",
                       __FILE__, __LINE__);
              MPI_Abort (MPI_COMM_WORLD, -1);
            }
        }//end of if bucket is brief bucket
        else
        {
            buck = (Bucket *) tempptr;
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
        }//end of if bucket is not neighbor info

      }//end of for loop go through received info
      counter += recv_info[info_num * i];
    }
  delete []recv_neigh_array;

  // wait and then delete sent info that was already received
  for (i = 0; i < numprocs; i++)
    if (send_info[info_num * i] != 0)
      j = MPI_Wait ((send_request + i), &status); //make sure all info that was send out is received before go to next step

  delete []send_neigh_array;

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

/****************************************************
 *     migrate particles within the buckets
 ****************************************************/
  send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[info_num * i + 2];

  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] =
      counter_send_proc[i - 1] + send_info[info_num * (i - 1) + 2];

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
    if (send_info[info_num * i + 2] != 0)
    {
      j =
        MPI_Isend ((send_part_array + counter), send_info[info_num * i + 2],
                   PARTICLE_TYPE, i, part_tag, MPI_COMM_WORLD,
                   (send_request + i));
      counter += send_info[info_num * i + 2];
    }

  //pack and send brief bucket here--> not necessary to wait for neighbor info and particles to finish communicate.
  //But allocate too many send buffer might take up too much memory.

  // Wait for recieves and unpack particles
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[info_num * i + 2] != 0)
    {
      j = MPI_Wait ((recv_request + info_num * i + 2), &status);
      for (j = 0; j < recv_info[info_num * i + 2]; j++)
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
    if (send_info[info_num * i + 2] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_part_array;

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif
  ////////////////////////////////////////////////////////////
  // can now migrate any of the buckets which need to be moved
  ////////////////////////////////////////////////////////////

  // first pack buckets in an array
  send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[info_num * i + 1];

  BucketPack *send_buck_array = new BucketPack[send_count];

  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] =
      counter_send_proc[i - 1] + send_info[info_num * (i - 1) + 1];

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

  // now send out packed elements
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (send_info[info_num * i + 1] != 0)
    {
      j =
        MPI_Isend ((send_buck_array + counter), send_info[info_num * i + 1],
                   BUCKET_TYPE, i, buck_tag, MPI_COMM_WORLD,
                   (send_request + i));
      counter += send_info[info_num * i + 1];
    }

  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[info_num * i + 1] != 0)
    {
      j = MPI_Wait ((recv_request + info_num * i + 1), &status);
      for (j = 0; j < recv_info[info_num * i + 1]; j++)
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
//  delete []recv_request;
//  delete []recv_info;

  // wait for send info before deleting
  for (i = 0; i < numprocs; i++)
    if (send_info[info_num * i + 1] != 0)
      j = MPI_Wait ((send_request + i), &status);

  delete []send_buck_array;

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

  ////////////////////////////////////////////////////////////
  // can now migrate any of the brief buckets which need to be moved
  ////////////////////////////////////////////////////////////

  // first pack buckets in an array
  send_count = 0;
  for (i = 0; i < numprocs; i++)
    send_count += send_info[info_num * i + 3];

  BriefBucketPack *send_briefbuck_array = new BriefBucketPack[send_count];

  counter_send_proc[0] = 0;
  for (i = 1; i < numprocs; i++)
    counter_send_proc[i] =
      counter_send_proc[i - 1] + send_info[info_num * (i - 1) + 3];

  itr->reset ();
  while ((tempptr=itr->next ()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
	  {
		  if (myid != breif_buck->get_myprocess ()) //This bucket need to be send to other processes
		  {
		      int myprocess = breif_buck->get_myprocess ();

		      assert ((myprocess > -1) && (myprocess < numprocs));
		      pack_bucket ((send_briefbuck_array + counter_send_proc[myprocess]), breif_buck,
		                   myprocess);
		      counter_send_proc[myprocess] += 1;

		      BG_mesh->remove (breif_buck->getKey ());
		      delete breif_buck;
		  }
	  } //end of if bucket is brief bucket
	  else
		  continue;
  }//end of while loop go through all buckets

  delete []counter_send_proc; //After all info is packed

  // now send out packed elements
  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (send_info[info_num * i + 3] != 0)
    {
      j =
        MPI_Isend ((send_briefbuck_array + counter), send_info[info_num * i + 3],
        		    BRIEF_BUCKET_TYPE, i, briefbuck_tag, MPI_COMM_WORLD,
                   (send_request + i));
      counter += send_info[info_num * i + 3];
    }

  counter = 0;
  for (i = 0; i < numprocs; i++)
    if (recv_info[info_num * i + 3] != 0)
    {
      j = MPI_Wait ((recv_request + info_num * i + 3), &status);
      for (j = 0; j < recv_info[info_num * i + 3]; j++)
      {
        // upack bucket data
    	BriefBucket *bfbuck = new BriefBucket ();

        unpack_bucket ((recv_briefbuck_array + counter), bfbuck, myid);
        BG_mesh->add (bfbuck->getKey (), bfbuck);
        counter++;
      }
    }

#ifdef DEBUG2
  fclose (fp);
#endif

  delete []recv_briefbuck_array;

  delete []recv_request;
  delete []recv_info;

  // wait for send info before deleting
  for (i = 0; i < numprocs; i++)
    if (send_info[info_num * i + 3] != 0)
      j = MPI_Wait ((send_request + i), &status);

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

  delete []send_briefbuck_array;

  delete []send_request;
  delete []send_info;
  delete itr;

  return;
}
