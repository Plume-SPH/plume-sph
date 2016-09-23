
/*
 * repartition_BSFC.cc
 *
 *  Created on: Oct 18, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cassert>
#include <vector>
#include <algorithm>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <hilbert.h>
#include <constant.h>

#include "exvar.h"
#include "buckhead.h"
#include "repartition_BSFC.h"

#include <properties.h>

#ifdef DEBUG
#  include <debug_header.h>
#endif

/* decent values for this partitioning scheme, the user can
 * set them in the as parameters to tune for better performance 
 */
/* minimum amount of coarse bins on each processor */
#define BINS_PER_PROC 25

/* amount of subbins a bin is divided into */
#define SUBBINS_PER_BIN 25 

/* amount of refinement of the bins */
#define MAX_REFINEMENT_LEVEL 10

/* if you're going to send fewer than this number of 
 * buckets to the processor before or after you on the space filling curve
 * don't bother because the communication isn't worth it, this number is 
 * machine dependent, should probably be much larger than 5
 */
#define MIN_NUM_2_SEND 5

int
repartition (vector < BucketHead > & PartitionTable, THashTable * P_table,
             HashTable * BG_mesh, MatProps * matprops, int * my_comm)
{
  int i, j, k;                  /* local variables */
  int num_local_objects;        /* the number of objects this processor owns */
  float total_weight = 0, wght;      /* the amount of work for all objects */
  const double TINY = 0.1e-4;
  int balanced_flag;
  int number_of_cuts = 0;
  int bits_used = 0;
  int size_of_unsigned = sizeof (unsigned);
  unsigned sfc_key[2];
  unsigned buck_key[KEYLENGTH];
  double xx[2];

  double bucket_size = PARTICLE_DENSITY * (matprops->smoothing_length);

  // direction vectors for neighbors
  int Up[DIMENSION] = { 0, 0, 2 };
  int Down[DIMENSION] = { 0, 0, 1 };

  int myid, numprocs;
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  if (numprocs < 2)
    return 0;


#ifdef DEBUG2
  char filename[20];
  sprintf (filename, "debug%02d.txt", myid);
  FILE * fp = fopen (filename, "w");
#endif

#ifdef DEBUG
   bool do_check = false;
#endif

  Bucket * buck = NULL;
  Bucket *curr_buck;
  BriefBucket *breif_buck = NULL;
  void * tempptr =NULL;

  /* get application data (number of objects, ids, weights, and coords */
  vector < BucketHead >::iterator ibuck;
  vector < float > weights;
  for (ibuck = PartitionTable.begin (); ibuck != PartitionTable.end (); ibuck++)
  {
    for (i = 0; i < KEYLENGTH; i++)
      buck_key[i] = ibuck->get_buck_head ()[i];

    tempptr = BG_mesh->lookup (buck_key);
	breif_buck = (BriefBucket *) tempptr;
	if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
	{
		wght = 0.;
	}
	else
	{
	    curr_buck = (Bucket *) tempptr;
	    assert (curr_buck);
	    wght = 0.;

	    vector < TKey > particles = curr_buck->get_plist ();
	    vector < TKey >::iterator p_itr;
	    for (p_itr = particles.begin (); p_itr != particles.end (); p_itr++)
	    {
	        Particle *p_curr = (Particle *) P_table->lookup (*p_itr);
	        if ((!p_curr->is_guest ()))
	        	switch (p_curr->get_bc_type ())
	        	{
	        	case 100:
	        		wght +=realp_load;
	        		break;
	        	case 2:
	        		wght +=wallp_load;
	        		break;
	        	case 0:
	        		wght +=eruptp_load;
	        		break;
	        	case 1:
	        		wght +=pressp_load;
	        		break;
	        	}
	     }//go through all particles in the bucket
	}//end of if bucket is not brief bucket
    weights.push_back (wght);
    total_weight += wght;
  }//end of while loop go through all buckets in the partition table

  double * work_per_proc = new double [numprocs];
  MPI_Allgather (&total_weight, 1, MPI_FLOAT,
                 work_per_proc, 1, MPI_FLOAT, MPI_COMM_WORLD); //Work load is averaged among several processes

  double global_weight = 0;
  for (i = 0; i < numprocs; i++)
    global_weight += work_per_proc[i];

  // Allocate memeory for SFC_VERTICES array
  num_local_objects = (int) PartitionTable.size ();
  sfc_vertex_t * sfc_vertices = new sfc_vertex_t [num_local_objects];

  // fill up the sfc_vertices array which stores 
  // all the necessary info about the load-balancing objects
  for (j = 0; j < num_local_objects; j++)
  {
    sfc_vertices[j].lb_weight =  weights[j];
    for (k = 0; k < 2; k++)
      sfc_vertices[j].sfc_key[k] = PartitionTable[j].getKey ()[k];
    for (k = 0; k < KEYLENGTH; k++)
      sfc_vertices[j].obj_key[k] = PartitionTable[j].get_buck_head()[k];
  }

  // allocate space for atucal work allocated
  float * global_actual_work_allocated = new float [numprocs+1];
  float * work_percent_array = new float [numprocs];
  unstructured_communication verts_in_cut_info;
  verts_in_cut_info.used_flag = 0;

  /* create bins, fill global weight vector and perform initial partition --->Based on SFC only*/
  BSFC_create_bins (num_local_objects, sfc_vertices, & bits_used,
                    size_of_unsigned, global_actual_work_allocated, 
                    work_percent_array, & total_weight, & balanced_flag,
                    & verts_in_cut_info, & number_of_cuts, BINS_PER_PROC, myid,
                    numprocs);

  if ( balanced_flag != BSFC_BALANCED)
  {
    
    float * work_prev_allocated = NULL;

    if (verts_in_cut_info.recv_count == 0 || myid == 0)
      balanced_flag = BSFC_BALANCED;

    BSFC_create_refinement_info (& number_of_cuts, global_actual_work_allocated, 
                                 total_weight, work_percent_array, 
                                 verts_in_cut_info, & work_prev_allocated, myid,
                                 numprocs);

    int * ll_bins_head = new int [number_of_cuts + 1];

    if (number_of_cuts == 0)
      balanced_flag = BSFC_BALANCED;

    if (ll_bins_head != NULL)
      ll_bins_head[number_of_cuts] = 0;

    for (i = 0; i < number_of_cuts; i++)
      ll_bins_head[i] = -1;

    int * local_balanced_flag_array = new int [number_of_cuts + 1];
    for (i = 0; i < number_of_cuts; i++)
      local_balanced_flag_array[i] = BSFC_BALANCED;
    local_balanced_flag_array[number_of_cuts] = balanced_flag;

    /* refine bins until a satisfactory partition tolerance in obtained */
    int refinement_level_counter = 0;
    while (balanced_flag != BSFC_BALANCED &&
           refinement_level_counter < MAX_REFINEMENT_LEVEL)
    {
      BSFC_refine_partition (& balanced_flag, & bits_used, 
                             verts_in_cut_info.recv_count,
                             verts_in_cut_info.recv_sfc_vert,
                             work_percent_array, total_weight,
                             global_actual_work_allocated, number_of_cuts,
                             ll_bins_head, work_prev_allocated, SUBBINS_PER_BIN,
                             local_balanced_flag_array, myid, numprocs);
      refinement_level_counter++;
    }
    delete [] local_balanced_flag_array;
    delete [] ll_bins_head;
    if ( work_prev_allocated )
      delete [] work_prev_allocated;
  }

  /* if the load-balancing objects were moved to different processors,
     we need to move them back now */
  if (verts_in_cut_info.used_flag != 0)
  {
    int recv_count = 0;
    sfc_vertex_t * send_sfc_vert = 
      new sfc_vertex_t [verts_in_cut_info.send_count];

    // fill up the send array
    int * proc_counter = new int [numprocs];
    proc_counter[0] = 0;
    for (i = 1; i < numprocs; i++)
      proc_counter[i] = proc_counter[i-1] + 
                        verts_in_cut_info.send_procs_ptr[i-1];
    int tag = 21054;
    recv_count = 0;
    MPI_Request * send_request = new MPI_Request [numprocs];
    MPI_Request * recv_request = new MPI_Request [numprocs];
    for (i = 0; i < numprocs; i++)
    {
      if ( i != myid )
      {
        if (verts_in_cut_info.send_procs_ptr[i] != 0)
          j = MPI_Irecv ((send_sfc_vert + proc_counter[i]), 
                         verts_in_cut_info.send_procs_ptr[i], LB_VERT_TYPE, i,
                         tag, MPI_COMM_WORLD, (send_request + i));
 
        if (verts_in_cut_info.recv_procs_ptr[i] != 0)
          j = MPI_Isend (& verts_in_cut_info.recv_sfc_vert[recv_count],
                         verts_in_cut_info.recv_procs_ptr[i], LB_VERT_TYPE,
                         i, tag, MPI_COMM_WORLD, (recv_request + i));
      }
      recv_count += verts_in_cut_info.recv_procs_ptr[i];
    }
    // wait until the info is send and received
    MPI_Status status;
    for (i = 0;  i < numprocs; i++)
      if ( i != myid )
      {
        if (verts_in_cut_info.send_procs_ptr[i] != 0)
          j = MPI_Wait ((send_request + i), & status);

        if (verts_in_cut_info.recv_procs_ptr[i] != 0)
          j = MPI_Wait ((recv_request + i), & status);
      }
    recv_count = 0;
    for (i = 0; i < myid; i++)
      recv_count += verts_in_cut_info.recv_procs_ptr[i];

    for (i = 0; i < num_local_objects; i++)
      if ( sfc_vertices[i].cut_bin_flag == BSFC_CUT)
      {
        if (sfc_vertices[i].destination_proc != myid)
        {
          j = sfc_vertices[i].destination_proc;
          sfc_vertices[i].destination_proc =
            send_sfc_vert[proc_counter[j]].destination_proc;
          proc_counter[j]++;
        }
        else
        {
          sfc_vertices[i].destination_proc = 
            verts_in_cut_info.recv_sfc_vert[recv_count].destination_proc;
          recv_count++;
        }
      }
    delete [] proc_counter;
    delete [] send_request;
    delete [] recv_request;
    delete [] send_sfc_vert;
  }
  delete [] global_actual_work_allocated;
  delete [] work_percent_array;

  // update proc-info
  for (i = 0; i < num_local_objects; i++)
  {
	tempptr=BG_mesh->lookup (sfc_vertices[i].obj_key);
	breif_buck = (BriefBucket*) tempptr; //Whether this works? ---> need test!
	breif_buck->put_myprocess (sfc_vertices[i].destination_proc);
//    do
//    {
//      buck = (Bucket *) BG_mesh->lookup (buck->which_neigh (Up));
//      buck->put_myprocess (sfc_vertices[i].destination_proc);
//    }
//    while (buck->which_neigh_proc (Up) != -1);
  }
 
#ifdef DEBUG2
  // check final load-distribution
  float * local_load = new float [numprocs];
  float * global_load = new float [numprocs];
  for (i = 0; i < numprocs; i++)
    local_load[i] = global_load[i] = 0.;

  for (i = 0; i < num_local_objects; i++)
    local_load[sfc_vertices[i].destination_proc] +=
      sfc_vertices[i].lb_weight;

  MPI_Reduce (local_load, global_load, numprocs, 
              MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid == 0)
    for (i = 0; i < numprocs; i++)
    {
      fprintf (stdout, "%d -> %f\n", i, global_load[i]);
      fflush (stdout);
    }
  delete [] local_load;
  delete [] global_load;
#endif

#ifdef DEBUG
    bool check_sndspd =false;
    bool ng_sndspd =false ;
    if (check_sndspd)
    {
        ng_sndspd = check_particles_sndspd (P_table);
        if (ng_sndspd)
  	      cout << "negative sound speed shows up before this point!" << endl;
    }
#endif

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

  // a name can't be more self-explanatory
  MPI_Barrier (MPI_COMM_WORLD);
  BSFC_update_and_send_elements (myid, numprocs, bucket_size, P_table, BG_mesh);

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

#ifdef DEBUG
    if (check_sndspd)
    {
        ng_sndspd = check_particles_sndspd (P_table);
        if (ng_sndspd)
  	      cout << "negative sound speed shows up before this point!" << endl;
    }
#endif

  // reset communication flags 
  for (i = 0; i < numprocs; i++)
    my_comm [i] = 0;

  // update the repartion info
  PartitionTable.clear();
  HTIterator *itr = new HTIterator (BG_mesh);
  const double * mindom = BG_mesh->get_minDom();
  const double * maxdom = BG_mesh->get_maxDom();
  unsigned keylen = KEYLENGTH;
  double normc[2];

//  while ((tempptr=itr->next ()))
//  {
//	  breif_buck = (BriefBucket *) tempptr;
//	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
//	  {
//		  const int * neigh_proc = breif_buck->get_neigh_proc ();
//		  for (i = 0; i < NEIGH_SIZE; i++)
//		  if (neigh_proc[i] > -1)
//		      my_comm[neigh_proc[i]] = 1;
//
//		  Key bkey = breif_buck->getKey ();
//		  BucketHead bhead (bkey.key, bkey.key);
//		  PartitionTable.push_back (bhead);
//	  }
//	  else
//	  {
//		  buck = (Bucket*) tempptr;
//		  const int * neigh_proc = buck->get_neigh_proc ();
//		  for (i = 0; i < NEIGH_SIZE; i++)
//		  if (neigh_proc[i] > -1)
//		      my_comm[neigh_proc[i]] = 1;
//
//		  Key bkey = buck->getKey ();
//		  BucketHead bhead (bkey.key, bkey.key);
//		  PartitionTable.push_back (bhead);
//	  }
//  }

#ifdef DEBUG
	if (do_check)
	   BG_mesh_check(BG_mesh);
#endif

  itr->reset ();
  while ((breif_buck = (BriefBucket *) itr->next ())) //Whether this works? ---> need test!
  {
	assert(breif_buck);
    // find out communication buddies of this prcess
    const int * neigh_proc = breif_buck->get_neigh_proc ();
    for (i = 0; i < NEIGH_SIZE; i++)
      if (neigh_proc[i] > -1)
        my_comm[neigh_proc[i]] = 1;

//    if (buck->which_neigh_proc (Down) == -1)
//    {
      Key bkey = breif_buck->getKey ();
//      for (i = 0; i < DIMENSION; i++)
//      {
//        xx[i] = (*(breif_buck->get_mincrd () + i) + *(breif_buck->get_maxcrd () + i)) * 0.5;
//        normc[i] = (xx[i] - mindom[i]) / (maxdom[i] - mindom[i]);
//      }
//      HSFC2d (normc , & keylen, sfc_key);
//      BucketHead bhead (sfc_key, bkey.key);

//      Key sfc_bkey = breif_buck->getKey ();
      BucketHead bhead (bkey.key, bkey.key);
      PartitionTable.push_back (bhead);
//    }
  }//end of while loop

  // no self communication
  my_comm[myid] = 0;

#ifdef DEBUG2
      fclose (fp);
#endif

  // sort bucket-heads
  sort (PartitionTable.begin (), PartitionTable.end ());

  // clean up
  delete[]work_per_proc;
  delete[]sfc_vertices;
  delete itr;

  return 1;
}
