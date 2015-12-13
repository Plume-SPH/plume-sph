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
//	  int PRESSURE_GHOST = 1;
	  int REAL = 100;
	  double pos[DIMENSION];
	  HTIterator * itr = new HTIterator (BG_mesh);
	  Bucket * Bnd_buck = NULL;
//	  Bucket * Curr_buck = NULL;
	  Bucket * neigh = NULL;

	  vector < TKey > plist;
	  vector < TKey >::iterator p_itr;

//	  Particle * p_curr;

	 while ((Bnd_buck = (Bucket *) itr->next ()))
	 {
		 //make sure that 1)the bucket is not empty, 2)has_involved = 0, 3) is not underground buckets 4) is not PRESS_BC.
		 //Note: in the old version of the code, the PRESS_BC is not excluded, actually, it should be, because "extension of the domain should stop when it reached to the largest domain" ---> This should not influence adding of pressure ghost because when adding pressure ghost, we do not care about whether one of our neigh has_involved or has_potential_involved, we only care about whether has_inolved==0
		 //The new adjustment for 3D domain decomposition: we only care about no-guest bucket, data for guest bucket will be syn after ths step.
		 //---> As a results we will only consider the no-guest bucket
	      if (!(Bnd_buck->get_has_involved()) && ((Bnd_buck->get_plist ()).size()) && Bnd_buck->get_bucket_type()!=UNDERGROUND && Bnd_buck->get_bucket_type()!=PRESS_BC && !Bnd_buck->is_guest())
		  {
	    	  const int * neigh_proc = Bnd_buck->get_neigh_proc ();
	    	  Key * neighbors = Bnd_buck->get_neighbors ();

	    	  for (i = 0; i < NEIGH_SIZE; i++)
	    	  {
	    	     if (neigh_proc[i] > -1)
	    	     {
	                // some neighs may not of available on current process --> This is true for guest buckets, so we have to syn data after domain adapt
	            	neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);
	                if ( neigh && neigh->is_has_involved ())
	                {
	                  Bnd_buck->set_has_involved (false);
	                  Bnd_buck->set_has_potential_involved (true);
	                  plist = Bnd_buck->get_plist();

	                  if (Bnd_buck->get_bucket_type()==MIXED && (Bnd_buck->get_bucket_index())[4] ==-1) //bucket is on ground mixed
	                  {
	                	    Bnd_buck->set_pressure_ghost_particles (false); //But if Bnd_buck is MIXED on top or side, it should still have some pressure ghost particles
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
	                  }// end of bucket is on ground MIXED

	                  else if(Bnd_buck->get_bucket_type()==MIXED && (Bnd_buck->get_bucket_index())[4] >-1) //bucket is mixed except for on ground
	                  {
	                	    for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
		                	{
		                	    Particle *p_curr = (Particle *) P_table->lookup(*p_itr);
		                	    assert(p_curr);
		                	    //get pos of particle
		                	    for (i = 0; i < DIMENSION; i++)
		                	        pos[i] = *(p_curr->get_coords() + i);
		                	    // turn involved = 0 to be involved = 1;
		                	    if (!Bnd_buck->determine_escape(pos)) //If particle is not locate at outside of the domain. ---> These particles which locate outside of the domain will never be transfered into real
		                	    {
		                	       p_curr->put_bc_type (REAL); //turn pressure ghost into real
		                	       p_curr->set_involved_flag (POTENTIAL_INVOLVED); //involved to be 1 (from non-involved to potential involved!)
		                	    }
		                	}
	                  }// end of bucket is MIXED except for on ground
	                  else
	                  {
	                  	  Bnd_buck->set_pressure_ghost_particles (false); //But if Bnd_buck is MIXED on top or side, it should still have some pressure ghost particles
	                  	  for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
	                  	  {
	                  	      Particle *p_curr = (Particle *) P_table->lookup(*p_itr);
	                  	      assert(p_curr);

	                  	      p_curr->put_bc_type (REAL); //turn pressure ghost into real
	                  	      p_curr->set_involved_flag (POTENTIAL_INVOLVED); //involved to be 1 (from non-involved to potential involved!)
	                  	   }
	                  	}// end of if bucket is not MIXED
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
