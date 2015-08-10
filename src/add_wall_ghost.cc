/*
 * add_wall_ghost.cc
 *
 *  Created on: Aug 5, 2015
 *      Author: zhixuanc
 */

/*
 * add_pressure_ghost.cc
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
#include <thashtab.h>
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

void
add_wall_ghost (THashTable * P_table, HashTable * BG_mesh,
        MatProps * matprops, TimeProps *timeprops, int numproc, int myid)
{
	  int i, j, k, ii;
	  int tempid;
	  unsigned pkey[TKEYLENGTH];
	  unsigned tkeylen = TKEYLENGTH;
	  double mincrd[DIMENSION], maxcrd[DIMENSION];
	  double mindom[DIMENSION], maxdom[DIMENSION];
	  double normc[DIMENSION];
	  double pcrd[DIMENSION];
	  double poly[DIMENSION + 1];
	  double bnd[2*DIMENSION], flag[DIMENSION],index[2*DIMENSION];

	  Bucket *Curr_buck = NULL;
	  int num_particle = 0;

	  // start putting piles
	  double mass = matprops->particle_mass;
	  double smlen = matprops->smoothing_length;
	  double dx = smlen;
	  double dx2 = 0.5 * dx;

	  // direction vectors for neighbors
	  int Up[DIMENSION] = { 0, 0, 2 };
	  int Down[DIMENSION] = { 0, 0, 1 };

#ifdef DEBUG
     bool do_search = false;
     double check[DIMENSION] = {-250, -250, 4250};
//     double temp[DIMENSION];
#endif

	  // get min-max domain from hashtable, for key generation
	  for (i = 0; i < DIMENSION; i++)
	  {
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	  }

	  //Initial particle property
	  int bctp =2; //wall condition ghost particle type.
	  int ptype = 0; //real particle
	  double prss = 0.;
	  double masfrc = 0.;
	  double gmm = 1.4;
	  double sndspd = 340.;
	  int phs_num = 1; //phase 1, air particle
	  unsigned add_step = timeprops->step; // eruption particle adding step

      int Nx, Ny, Nz;
	// add particles
	  HTIterator * itr = new HTIterator (BG_mesh);
	  Bucket * Bnd_buck;

	  /*temporarily use the following way to get bnd, Finally, I will put bnd either
	   * in SimulProps or parameters.h
	   * will modify this part later
	  */
       bnd[0]=Lx_P[0];
       bnd[1]=Lx_P[1];
       bnd[2]=Ly_P[0];
       bnd[3]=Ly_P[1];
       bnd[4]=Lz_P[0];
       bnd[5]=Lz_P[1];

	  while ((Bnd_buck = (Bucket *) itr->next ()))
	 //--->The reason why need to make sure the bucket is active is that for MIXED bucket that need to add wall ghost, it should already have either real particles or pressure ghost particle on the above boundary portion of the bucket.
	 //But I need to make sure that wall ghost particles will not be added repeatedly
     //Why do I need to make sure that bucket is active? ---->because MIXED in which wall ghost needed should be active, as some pressure ghost has already been added and the bucket was marked as active.
	  if ((Bnd_buck->get_bucket_index())[4] == -1 && Bnd_buck->get_bucket_type () == MIXED && Bnd_buck->is_active() && !Bnd_buck->has_wall_ghost_particles ()) //make sure bucket is on-ground mixed, go through all on ground MiXED and then go down to underground buckets
	  {
//      	    Bnd_buck->mark_active ();

    	    for (i = 0; i < DIMENSION; i++)
    	    {
    	    	mincrd[i] = *(Bnd_buck->get_mincrd () + i);
    	    	maxcrd[i] = *(Bnd_buck->get_maxcrd () + i);
    	        index[2*i] = *(Bnd_buck->get_bucket_index () + 2*i);
    	        index[2*i+1] = *(Bnd_buck->get_bucket_index () + 2*i+1);
    	    }
    	    //generate air particles
    	    Nx = (int) round((maxcrd[0] - mincrd[0]) / dx);
    	    Ny = (int) round((maxcrd[1] - mincrd[1]) / dx);
    	    Nz = (int) round((maxcrd[2] - mincrd[2]) / dx);

    	    Curr_buck = Bnd_buck;
    	    for (i = 0; i < Nx; i++)
    	    	for (j = 0; j < Ny; j++)
    	    		for (k = 0; k < Nz; k++)
    	    		{
    	    				pcrd[0] = mincrd[0] + dx2 + i * dx;
    	    				pcrd[1] = mincrd[1] + dx2 + j * dx;
    	    				pcrd[2] = mincrd[2] + dx2 + k * dx;

    	    				//I need to figure out one way to make sure that particles will not be added repeatedly in MIXED bucket. ---> should keep in mind that efficiency is very important!
    	    				if (pcrd[2]<bnd[4])
    	    				{
    	    					for (ii = 0; ii < DIMENSION; ii++)
    	    					    normc[ii] = (pcrd[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    	    					THSFC3d (normc, add_step, &tkeylen, pkey);
    		    				// check for duplicates
    		    				if (P_table->lookup(pkey))
    		    				{
    		    				    fprintf(stderr, "ERROR: Trying to add particle "
    		    				            "twice on same location.\n");
    		    				    exit(1);
    		    				 }
    		    				Particle * pnew = new Particle(pkey, pcrd, mass, smlen, prss, masfrc, gmm, sndspd, phs_num, myid, bctp);

    		    				//default involved is zero
 					           if (Bnd_buck->is_guest())
 					           {
 					               pnew->put_guest_flag(true);
 					               tempid = Bnd_buck->get_myprocess ();
 					               pnew->put_my_processor(tempid);
 					           }

    		    				// add to hash-table
    		    				P_table->add(pkey, pnew);
    		    				num_particle++;
    		    				TKey tmpkey(pkey);

    		    				Curr_buck->add_wall_ghost_particle(tmpkey);
    	    				}
    	    		}

    	    // add wall ghost on under ground bucket
	        Bucket *Down_buck = (Bucket *) BG_mesh->lookup (Bnd_buck->which_neigh (Down));

	        Down_buck->mark_active ();// the down bucket is not marked as active yet! so need to mark it here

    	    for (i = 0; i < DIMENSION; i++)
    		{
    		   mincrd[i] = *(Down_buck->get_mincrd () + i);
    		   maxcrd[i] = *(Down_buck->get_maxcrd () + i);
    		}
    		//generate air particles
    		Nx = (int) round((maxcrd[0] - mincrd[0]) / dx);
    		Ny = (int) round((maxcrd[1] - mincrd[1]) / dx);
    		Nz = (int) round((maxcrd[2] - mincrd[2]) / dx);

    		Curr_buck = Down_buck;
    		for (i = 0; i < Nx; i++)
    		    for (j = 0; j < Ny; j++)
    		        for (k = 0; k < Nz; k++)
    		        {
    		        		pcrd[0] = mincrd[0] + dx2 + i * dx;
    			        	pcrd[1] = mincrd[1] + dx2 + j * dx;
    			        	pcrd[2] = mincrd[2] + dx2 + k * dx;

    			            for (ii = 0; ii < DIMENSION; ii++)
    			               normc[ii] = (pcrd[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

 	    					THSFC3d (normc, add_step, &tkeylen, pkey);
 		    				// check for duplicates
 		    				if (P_table->lookup(pkey))
 		    				{
 		    				    fprintf(stderr, "ERROR: Trying to add particle "
 		    				           "twice on same location.\n");
 		    				    exit(1);
 		    				}
 		    				Particle * pnew = new Particle(pkey, pcrd, mass, smlen, prss, masfrc, gmm, sndspd, phs_num, myid, bctp);

 		    				//default involved is zero
 					        if (Down_buck->is_guest())
					        {
					            pnew->put_guest_flag(true);
					            tempid = Down_buck->get_myprocess();
					            pnew->put_my_processor(tempid);
					        }



 		    				// add to hash-table
 		    				P_table->add(pkey, pnew);
 		    				num_particle++;
 		    				TKey tmpkey(pkey);

 		    				Curr_buck->add_pressure_ghost_particle(tmpkey);
    		         }

	  }//finish go though all buckets


	  //delete wall ghost in erupt vent

	  //determine the rough range of eruption duck
		double range_x[2];
		double range_y[2];
		double range_z[2];
	  	range_x[0] = -rv_P;
	  	range_x[1] = rv_P;
	  	range_y[0] = -rv_P;
	  	range_y[1] = rv_P;
	    range_z[1] = 0.; // not exact, should use ground height
	    range_z[0] = range_z[1]-(matprops->smoothing_length)*1.5*PARTICLE_DENSITY;
	    unsigned tempkey[TKEYLENGTH];
	    double dist;
	    double rvsq=rv_P*rv_P;

	    vector < TKey > pnew;//particle key
	    vector < TKey > plist;
	    vector < TKey >::iterator p_itr;
	    Particle *pj;
	    itr->reset();
	    //go through all buckets
	    while ((Bnd_buck = (Bucket *) itr->next ()))
	      if (Bnd_buck->is_erupt () && Bnd_buck->has_wall_ghost_particles ())
	      {
	    	  //check all particle in the erupt bucket and remove them when necessary!
	    	  pnew.clear();
	    	  plist = Bnd_buck->get_particle_list ();
	    	  for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
	    	  {
	    		pj = (Particle *) P_table->lookup(*p_itr);
	    		assert (pj);
	    		if (!pj->is_erupt_ghost())
	    		{
	    		  for (k=0; k<DIMENSION; k++)
	    			  pcrd[k] = *(pj->get_coords() + k);

	    		  dist = 0;
	    		  for (j=0; j<2; j++)
	    		      dist += (pcrd[j]*pcrd[j]);

	    		  if ((dist <= rvsq) && pcrd[2]<=range_z[1] && pcrd[2]>=range_z[0] && pj->is_wall_ghost() ) // if particle is in the duct and particle is wall ghost
	    		  {
	    			  for (k = 0; k<TKEYLENGTH; k++)
	    			      tempkey[k] = pj->getKey().key[k];

	    			  /*
	    			   * Here what I did is only remove them from P_table and bucket particle list!
	    			   * But the particle in other particles neighbour list is not deleted!
	    			   * And the particle, if it has image, the image should also be removed!
	    			   */
	    			  P_table->remove(tempkey); //remove from the hashtable
	    		  }
	    		  else
	    			  pnew.push_back(*p_itr);//remove from the from particle list!
	    		}//end of only check non_erupted particles
	    	  }

	    	  //update particle list
	    	  Bnd_buck->put_new_plist (pnew);
	    	  Bnd_buck->update_particles();
	      }// end of loop go through all buckets

	  delete itr;

	  return;
}
