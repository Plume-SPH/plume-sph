/*
 * adapt_domain.cc
 *
 *  Created on: Aug 7, 2015
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

#include <hashtab.h>
#include <thashtab.h>
#include <hilbert.h>
#include <bucket.h>
#include <particler.h>
#include <multiproc.h>
#include <properties.h>

#include "constant.h"
#include "parameters.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

/*
 * What happened before this function calling is: scan the most outside potential involved bucket, if any of particles becomes involved (two ways), turn the bucket to be has_involved (has_involved = 3).
 * function which scans these buckets which were originally pressure ghost buckets, if any of its neighbor buckets is has_involved>0, 1) turn it to be a potential_involved_bucket, 2) turn the particle in the bucket to be real particle, 3)particle involved =1 (potential involved);
 * What will happen after this function calling is: add a new layer of pressure ghost bucket as the old pressure ghost bucket become potential_involved.
 */
void adapt_domain(THashTable * P_table, HashTable * BG_mesh, int numproc, int myid)
{
	  int i;
	  int PRESSURE_GHOST = 1;
	  int REAL = 100;
	  HTIterator * itr = new HTIterator (BG_mesh);
	  Bucket * Bnd_buck = NULL;
	  Bucket * Curr_buck = NULL;
	  Bucket * neigh = NULL;

	  vector < TKey > plist;
	  vector < TKey >::iterator p_itr;

//	  Particle * p_curr;

	 while ((Bnd_buck = (Bucket *) itr->next ()))
	 {
		 //make sure that 1)the bucket is not empty, 2)has_involved = 0, 3) is not underground bucket
	      if (!(Bnd_buck->get_has_involved()) && ((Bnd_buck->get_plist ()).size()) && (Bnd_buck->get_bucket_index())[5] > -1)
		  {
	    	  const int * neigh_proc = Bnd_buck->get_neigh_proc ();
	    	  Key * neighbors = Bnd_buck->get_neighbors ();

	    	  for (i = 0; i < NEIGH_SIZE; i++)
	    	  {
	    	     if (neigh_proc[i] > -1)
	    	     {
	                // some neighs may not of available on current process
	            	neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);
	                if ( neigh && neigh->is_has_involved ())
	                {
	                  Bnd_buck->set_has_involved (false);
	                  Bnd_buck->set_has_potential_involved (true);
	                  Bnd_buck->set_pressure_ghost_particles (false);
	                  plist = Bnd_buck->get_plist();

	                  if ((Bnd_buck->get_bucket_index())[4] > -1) //if not MIXED --> is OVERGROUND
	                  {
	                	for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
	                	{
	                	    Particle *p_curr = (Particle *) P_table->lookup(*p_itr);
	                	    assert(p_curr);

	                	    p_curr->put_bc_type (REAL); //turn pressure ghost into real
	                	    p_curr->set_involved_flag (POTENTIAL_INVOLVED); //involved to be 1 (from non-involved to potential involved!)
	                	}
	                  }// end of if bucket is not MIXED
	                  else
	                  {
		                	for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
		                	{
		                	    Particle *p_curr = (Particle *) P_table->lookup(*p_itr);
		                	    assert(p_curr);
		                	    // turn involved = 0 to be involved = 1;
		                	    if (p_curr->is_press_ghost ()) //p_curr might be wall ghost in a MIXED bucket, for OVERGROUND bucket, all particle is pressure ghost!
		                	    {
		                	       p_curr->put_bc_type (REAL); //turn pressure ghost into real
		                	       p_curr->set_involved_flag (POTENTIAL_INVOLVED); //involved to be 1 (from non-involved to potential involved!)
		                	    }
		                	}
	                  }// end of bucket is MIXED

	                	goto next_bnd_buck;
	                }// if any neigh is has_involved (has_involved = 2)
	    	     }
	    	  } // end of loop go through all buckets

		  }
	      next_bnd_buck:
	      1; //This line is useless actually, but without this one, the code would be invalid ----> lable : }  is invalid!
	 }// end of go through all buckets

	return;
}
