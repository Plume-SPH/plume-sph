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
add_pressure_ghost (THashTable * P_table, HashTable * BG_mesh,
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

	  Bucket *Curr_buck = NULL;

	  // direction indices on upper bucket
//	  int Up[DIMENSION] = { 0, 0, 2 };
	  int num_particle = 0;

	  // start putting piles
	  double mass = matprops->particle_mass;
	  double smlen = matprops->smoothing_length;
	  double dx = smlen;
	  double dx2 = 0.5 * dx;

//#ifdef DEBUG
////     bool do_search = false;
////     double check[DIMENSION] = {-250, -250, 4250};
////     double temp[DIMENSION];
//#endif

	  // get min-max domain from hashtable, for key generation
	  for (i = 0; i < DIMENSION; i++)
	  {
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	  }

	  //Initial particle property
	  int bctp =1; //pressure condition ghost particle type.
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
	  Bucket * neigh = NULL;

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
      // will not add pressure ghost on UNDERGROUND bucket
	  while ((Bnd_buck = (Bucket *) itr->next ()))
	  // 1) Bnd_buck has involved flag = 0
	  // 2) pressure ghost has not been added in Bnd_buck yet
	  // 3) Bnd_buck is not erupt ghost bucket
	  // 4) Because the overall work flow is changed slightly, we will not try add pressure ghost for guest bucket, so also need to make sure that : Bnd_buck is not guest.
      if (!(Bnd_buck->get_has_involved()) && !Bnd_buck->has_pressure_ghost_particles() && !Bnd_buck->is_erupt () && !Bnd_buck->is_guest())
	  {
        // if any neighbor has_potential_involved or has_involved
    	//1) I need to make sure that Bnd_buck is not guest bucket, otherwise, no neighbor info will be available --> So finally, pressure ghost for guest bucket will not be added at this step.
    	//2) need to syn data after this function
        const int * neigh_proc = Bnd_buck->get_neigh_proc ();
        Key * neighbors = Bnd_buck->get_neighbors ();

        for (i = 0; i < NEIGH_SIZE; i++)
        {

              if (neigh_proc[i] > -1)
              {
                  // some neighs may not of available on current process
            	  neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);
                  if ( neigh && neigh->get_has_involved())
                	  //sometimes, neigh might be on other processors
                	  //--> actually, for guest buckets, it is possible that some of its neighbors can not be found on current process!
                	  //--> So syn is necessary after adding of pressure ghost, because, it is possible that some guest buckets that should contain pressure ghost particles will not add pressure ghost particles at this step.
                	  //--> Finally, I decided that I will not try to add pressure ghost particles within this function, instead, I will syn data after this function.
                  {
                	Bnd_buck->mark_active ();
              	    if (Bnd_buck->get_bucket_type () == MIXED)
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
              	    				 ptype = 0; // is real
              	    				 pcrd[0] = mincrd[0] + dx2 + i * dx;
              	    				 pcrd[1] = mincrd[1] + dx2 + j * dx;
              	    				 pcrd[2] = mincrd[2] + dx2 + k * dx;

              	    				 if ((Bnd_buck->get_bucket_index())[4] == -1) //on ground MIXED
              	    				 {
              	    				   if (pcrd[2]>bnd[4]) //This is only for on-ground MIXED, for Mixed in the air, the bnd[4]=0.1, which does not make any sense. for over-ground MIXED and PRESSURE_BC, pressure ghost need to be full of the bucket.
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

             					           //determine initial parameter of air
             					          initial_air (pnew);

              		    				   // add to hash-table
              		    				   P_table->add(pkey, pnew);
              		    				   num_particle++;
              		    				   TKey tmpkey(pkey);

              		    				   Curr_buck->add_pressure_ghost_particle(tmpkey);
              	    				   }// make sure it is above ground.
              	    				 }// end of if bucket is on ground MIXED
              	    				 else //over ground MIXED ---> MIXED on the side of the domain and on the top of the domain
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

            		    				//default involved is zero ---> in the constructor
          					           if (Bnd_buck->is_guest())
          					           {
          					               pnew->put_guest_flag(true);
          					               tempid = Bnd_buck->get_myprocess ();
          					               pnew->put_my_processor(tempid);
          					           }

         					          //determine initial parameter of air
         					          initial_air (pnew);

            		    				// add to hash-table
            		    				P_table->add(pkey, pnew);
            		    				num_particle++;
            		    				TKey tmpkey(pkey);

            		    				Curr_buck->add_pressure_ghost_particle(tmpkey);
              	    				 }// end of else -->is over ground MIXED

              	    			 }

              	    }//finish bucket is Mixed
              	    else if (Bnd_buck->get_bucket_type () == OVERGROUND || (Bnd_buck->get_bucket_type () == PRESS_BC))
              	    {
              	    	for (i = 0; i < DIMENSION; i++)
              		     {
              		        mincrd[i] = *(Bnd_buck->get_mincrd () + i);
              		        maxcrd[i] = *(Bnd_buck->get_maxcrd () + i);
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

         					           //determine initial parameter of air
         					           initial_air (pnew);

           		    				    // add to hash-table
           		    				    P_table->add(pkey, pnew);
           		    				    num_particle++;
           		    				    TKey tmpkey(pkey);

           		    				    Curr_buck->add_pressure_ghost_particle(tmpkey);

              		        		}

              	      }//finish bucket is overground or pressure_BC

              	  goto next_bnd_buck;// no need to look any further
                  }//end of if neigh satisfies has_involved > 0

                }//make sure that current buck is not on the most out layer

        }//go through all neighbors

        next_bnd_buck:
        1;

	  }//finish go though all buckets

	  delete itr;

	  return;
}
