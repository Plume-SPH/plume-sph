/*
 * shift_brief_buck.cc
 *
 *  Created on: Aug 5, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

#include <mpi.h>
#include <hashtab.h>
#include <hilbert.h>
#include <bucket.h>
#include <particler.h>
#include <multiproc.h>
#include <properties.h>

#include "sph_header.h"
#include "constant.h"
#include "parameters.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

//function that shif the layer outside the pressure ghost boundary from brief bucket to bucket
void
shift_brief_buck (HashTable * BG_mesh, MatProps * matprops, TimeProps *timeprops, int myid)
{
	  int i, j;
	  double mindom[DIMENSION], maxdom[DIMENSION];
	  double mindom_o[DIMENSION], maxdom_o[DIMENSION];

#ifdef DEBUG
	  bool check_buck = false;
      bool check_mesh_err = false;
#endif

	  // get min-max domain from hashtable, for key generation
	  for (i = 0; i < DIMENSION; i++)
	  {
	    mindom[i] = *(BG_mesh->get_minDom() + i);
	    maxdom[i] = *(BG_mesh->get_maxDom() + i);
	  }

	// add particles
	  HTIterator * itr = new HTIterator (BG_mesh);
	  Bucket * Bnd_buck;
	  BriefBucket *breif_buck = NULL;
	  void * tempptr =NULL;
	  BriefBucket * breif_neigh =NULL;

 	  mindom_o[0]= Lx_P[0];  //original domain
 	  mindom_o[1]= Ly_P[0];  //original domain
 	  mindom_o[2]= Lz_P[0];  //original domain

      maxdom_o[0]= Lx_P[1];  //original domain
      maxdom_o[1]= Ly_P[1];  //original domain
 	  maxdom_o[2]= Lz_P[1];  //original domain

 	  Key btkey;
 	  double len_scale = matprops->LENGTH_SCALE;
 	  double bucket_size = PARTICLE_DENSITY * (matprops->smoothing_length);

      // will not add pressure ghost on UNDERGROUND bucket
	  while ((tempptr=itr->next ()))
	  {
		  breif_buck = (BriefBucket *) tempptr;
		  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
			  continue;
		  else
		  {
			  Bnd_buck = (Bucket*) tempptr;
			  // 1) Bnd_buck has involved flag = 0
			  // 2) has pressure ghost particle
	          // I did not get rid of guest bucket??
	  	      /*
	  	     * Assume Bnd_buck is not guest, but one of its neighbor is guest, this guest bucket will be switched to a bucket without any issue
	  	     * But the corresponding no-guest bucket on another process might not be switched to bucket as Band_buck is now a guest bucket on this process
	  	     * Then we need a inverse communicate : send info from guest to no-guest
	  	     * To avoid this, we added this section.
	  	     */
		      if (!(Bnd_buck->get_has_involved()) && Bnd_buck->has_pressure_ghost_particles() && Bnd_buck->get_erupt_flag ()==0)
			  {
		        const int * neigh_proc = Bnd_buck->get_neigh_proc ();
		        Key * neighbors = Bnd_buck->get_neighbors ();

		        for (i = 0; i < NEIGH_SIZE; i++)
		        {
		           if (neigh_proc[i] > -1)
		           {
						// some neighs may not of available on current process --> This is true for guest buckets, so we have to syn data after domain adapt
						breif_neigh = (BriefBucket *) BG_mesh->lookup (neighbors[i]); //Even though the bucket might be no brief bucket, but it still works: cast the pointer points to a Bucket to a pointer points to a BriefBucket
						if ( breif_neigh && breif_neigh->check_brief()) //If the neighbor exist and is brief bucket
						{
							Bucket * buck = NULL;
							switch_brief(breif_neigh, mindom, maxdom, mindom_o, maxdom_o, bucket_size,len_scale, &buck); //creation of new bucket is within this function
							buck->put_guest_flag(breif_neigh->is_guest());
							btkey = breif_neigh->getKey();
							BG_mesh->remove(btkey);
							delete breif_neigh;
							btkey = buck->getKey(); // should not necessary as key for briefbucket and bucket are the same.
							BG_mesh->add(btkey, buck);
						}// end of if the neighbor exist and is brief bucket

#ifdef DEBUG
  if (check_mesh_err)
	  BG_mesh_err_check (BG_mesh, myid);
#endif

		           }//make sure that current buck is not on the most out layer---> process id is valid

		        }//end of go through all neighbors
			 }//end of if
		}//end if else
	}//end of go through all buckets

    delete itr;

	return;
}
