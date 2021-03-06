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

#if CODE_DIMENSION==3 || CODE_DIMENSION==2
void
add_wall_ghost (THashTable * P_table, HashTable * BG_mesh, SimProps* simprops,
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
//	  double poly[DIMENSION + 1];
	  double bnd[2*DIMENSION], index[2*DIMENSION];

#ifdef DEBUG
      bool do_search = false;
      unsigned keycheck[TKEYLENGTH] = {71898628, 982285492, 0};
      unsigned keytemp[TKEYLENGTH] ;
#endif

	  Bucket *Curr_buck = NULL;
	  int num_particle = 0;

	  // start putting piles
	  double mass = matprops->particle_mass;
	  double smlen = matprops->smoothing_length;
	  double dx = smlen;
	  double dx2 = 0.5 * dx;

	  // direction vectors for neighbors
#if CODE_DIMENSION==3
	  int Down[DIMENSION] = { 0, 0, 1 };
#elif CODE_DIMENSION==2
	  int Down[DIMENSION] = { 0, 1 };   // in x direction, in the same column, in y direction, the down direction --> need to double check here!
#endif

	  // get min-max domain from hashtable, for key generation
	  for (i = 0; i < DIMENSION; i++)
	  {
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	  }

	  //Initial particle property
	  int bctp =2; //wall condition ghost particle type.
	  double prss = 0.;
	  double masfrc = 0.;

#if FLUID_COMPRESSIBILITY==0 //using EOS of ideal gas
	  double sndspd = 340.;
	  double gmm = 1.4;
#elif FLUID_COMPRESSIBILITY==1
	  double sndspd = 1482.;
	  double gmm = 7.0;
#endif

	  int phs_num = 1; //phase 1, air particle
	  unsigned add_step = timeprops->step; // eruption particle adding step

      int Nx, Ny, Nz;
	// add particles
	  HTIterator * itr = new HTIterator (BG_mesh);
	  Bucket * Bnd_buck;
	  BriefBucket *breif_buck = NULL;
	  void * tempptr =NULL;

	  /*temporarily use the following way to get bnd, Finally, I will put bnd either
	   * in SimulProps or parameters.h
	   * will modify this part later
	  */
	  for (i = 0; i < DIMENSION; i++)
	  {
	     bnd[i*2]=Ll_P[i];
	     bnd[i*2+1]=Lu_P[i];
	  }

 	  while ((tempptr=itr->next ()))
 	  {
 		  breif_buck = (BriefBucket *) tempptr;
 		  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
 			  continue;
 		  else
 		  {
 			  Bnd_buck = (Bucket*) tempptr;
	          //--->The reason why need to make sure the bucket is not empty is that for MIXED bucket that need to add wall ghost, it should already have either real particles or pressure ghost particle on up part of buckets which above boundary
	          //But I need to make sure that wall ghost particles will not be added repeatedly ---> !Bnd_buck->has_wall_ghost_particles ()
              //Why do I need to make sure that bucket is non_empty? ---->because MIXED in which wall ghost needed should be active, as some pressure ghost has already been added and the bucket was marked as active.
 			  //Please also notice that I need to make sure the bucket is not eruption bucket
 			  //To be consistent, should also make sure that the Bnd_buck is not guest bucket. ----> Should always followed by syn
 			  //Or you can add particle for guest buckets here ---> and not necessary to syn after adding!---> But as we are combining undergeround buckets and ongreound MIXED together, this routine does not work as guest bucket does not have information about its neighbors.
 			  //So I need to make sure that Bnd_buck is not guest bucket here.
 			  if (!Bnd_buck->is_guest() && Bnd_buck->is_onground_or_underground() && Bnd_buck->get_bucket_type () == MIXED && (Bnd_buck->get_plist ()).size() && !Bnd_buck->has_wall_ghost_particles () && Bnd_buck->get_erupt_flag ()==0 ) //make sure bucket is on-ground mixed, go through all on ground MiXED and then go down to underground buckets
 			  {
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
 		    	    				if (pcrd[2]<bnd[4]) // make sure newly added particle is underground
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

#ifdef DEBUG
 		    		    				if (do_search)
 		    		    				{
 		    		    					for (ii = 0; ii < TKEYLENGTH; ii++)
 		    		    						keytemp[ii] = pnew->getKey ().key[ii];

 		    		    					if (find_particle (keytemp, keycheck))
 		    		    						cout << "The particle found!" << endl;
 		    		    				}
#endif

 		    		    				//After creating these new particles, gave each particle a initial status
 		    		    				//The purpose is if the state of these wall ghost particles did not get updated, a "static" wall ghost particle will be added instead
 		    		    				//Which is not accurate, but will help make the code more stable
 		    		    				//In addition, in current version, only particle in the "On-ground MIXED" bucket will be updated based on image
 		    		    				//These particles in underground buckets will not been updated!
 		 					           initial_air (pnew, simprops);

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
 		    	    				}//end if the particle is below the ground
 		    	    		}

 		    	    // add wall ghost on under ground bucket
 		    	    //--> For 3D decomposition, if the Mixed bucket is guest bucket, it DOWN bucket will not be able to find, SO
 		    	    //1) make sure the current bucket is not guest bucket ---> This is not necessary, because I added particles for guest MIXED bucket, I should keep consistent!
 		    	    //2) On the process on which the current bucekt is host, the DOWN bucket (no matter the DOWN bucket is guest or not) can be found and, wall ghost particles can be added there.
 		    	    //3) I just need to point out that the idea that adding particles for ghost bucket to avoid syn still works for 3D decomposition. ---> but need to be more careful to avoid mistake like looking for neighbor for a guest bucket.
 		    	    if (!Bnd_buck->is_guest()) //I need this because the guest bucket does not have neighbor information and will not be able to find its down bucket. ---> So the idea that adding guest particle to avoid syn is not practical for adding wall ghost particles.
 		    	    {
 		    	        Bucket *Down_buck = (Bucket *) BG_mesh->lookup (Bnd_buck->which_neigh (Down));
 		                assert(Down_buck); //make sure Down_buck exist.

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

 		    		    				//After creating these new particles, gave each particle a initial status
 		    		    				//The purpose is if the state of these wall ghost particles did not get updated, a "static" wall ghost particle will be added instead
 		    		    				//Which is not accurate, but will help make the code more stable
 		    		    				//In addition, in current version, only particle in the "On-ground MIXED" bucket will be updated based on image
 		    		    				//These particles in underground buckets will not been updated!
 		 					            initial_air (pnew, simprops);

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

 		     		    				Curr_buck->add_wall_ghost_particle(tmpkey);
 		        		         }
 		    	    }//end of if bnd_bucket is guest
 			  }//end of if bucket is .....
 		  }//end of if bucket is not brief

	  }//finish go though all buckets


	  //delete wall ghost in erupt vent -------> To be honest, I do not think this is necessary!
	  //because
	  //1) wall ghost will be deleted while set up eruption condition
	  //2) all of the erupted bucket should be within initial domian, otherwise, the initial domain is too small.

	  delete itr;

	  return;
}

#endif //CODE_DIMENSION==2 or CODE_DIMENSION=3
