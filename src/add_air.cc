/*
 * add_air.cc
 *
 *  Created on: Mar 8, 2015
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
add_air (THashTable * P_table, HashTable * BG_mesh,
        MatProps * matprops, int numproc, int myid)
{
	  int i, j, k, ii;
	  unsigned pkey[TKEYLENGTH];
	  unsigned tkeylen = TKEYLENGTH;
	  double mincrd[DIMENSION], maxcrd[DIMENSION];
	  double mindom[DIMENSION], maxdom[DIMENSION];
	  double normc[DIMENSION];
	  double pcrd[DIMENSION];
	  double poly[DIMENSION + 1];
	  double bnd[2*DIMENSION], flag[DIMENSION],index[2*DIMENSION];
	  Bucket *Curr_buck = NULL;

	  // direction indices on upper bucket
	  int Up[DIMENSION] = { 0, 0, 2 };
	  int num_particle = 0;

	  // start putting piles
	  double mass = matprops->particle_mass;
	  double smlen = matprops->smoothing_length;
	  double dx = smlen;
	  double dx2 = 0.5 * dx;

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
	  int bctp; //boundary condition ghost particle type.
	  int ptype = 0; //real particle
	  double prss = 0.;
	  double masfrc = 0.;
	  double gmm = 1.4;
	  double sndspd = 340.;
	  int phs_num = 1; //phase 1, air particle
	  unsigned add_step=0; // eruption particle adding step

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
      if (!Bnd_buck->is_guest())
	  {
	    if (Bnd_buck->get_bucket_type () == MIXED )
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

#ifdef DEBUG
		                   if (do_search)
		                   {
		                        if (find_particle (pcrd, check))
			                    cout << "The particle found!" << endl;
		                    }
#endif

	    				   for (ii = 0; ii < DIMENSION; ii++)
	    					   flag[ii]=index[2*ii]+index[2*ii+1];

	    				   for (ii = 0; ii < DIMENSION; ii++)
	    					   if ((pcrd[ii]<bnd[2*ii]) || (pcrd[ii]>bnd[2*ii+1]))
	    						   ptype = 1;//ptype is now an index for determine bc_type

	    				   for (ii = 0; ii < DIMENSION; ii++)
	    				       normc[ii] = (pcrd[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

	    				   //determine bc_type
	    				   if (ptype ==0)
	    					   bctp = 100;
	    				   else
	    				   {
	    				       if ( pcrd[2]<bnd[4] )
	    					       bctp = 2; //wall ghost particle
	    				       else
	    					       bctp = 1 ;// pressure bc ghost particle
	    				   }

	    				   THSFC3d (normc, add_step, &tkeylen, pkey);

	    				   // check for duplicates
	    				   if (P_table->lookup(pkey))
	    				   {
	    				       fprintf(stderr, "ERROR: Trying to add particle "
	    				               "twice on same location.\n");
	    				       exit(1);
	    				   }
	    				   Particle * pnew = new Particle(pkey, pcrd, mass, smlen, prss, masfrc, gmm, sndspd, phs_num, myid, bctp);

	    				   // add to hash-table
	    				   P_table->add(pkey, pnew);
	    				   num_particle++;
	    				   TKey tmpkey(pkey);

	    				   switch (bctp)
	    				   {
	    				   case 1 :
	    					   Curr_buck->add_pressure_ghost_particle(tmpkey);
	    					   break;
	    				   case 2 :
	    					   Curr_buck->add_wall_ghost_particle(tmpkey);
	    					   break;
	    				   case 100 :
	    					   Curr_buck->add_real_particle(tmpkey);
	    					   break;
	    				   default:
	    					   cout << "bctp incorrect in function add_air!\n" << endl;
	    			      }
	    			 }

	    }//finish bucket is on ground and Mixed
	    else if (Bnd_buck->get_bucket_type () == OVERGROUND)
	    {
	    	bctp = 100;
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
#ifdef DEBUG
		                   if (do_search)
		                   {
		                        if (find_particle (pcrd, check))
			                    cout << "The particle found!" << endl;
		                    }
#endif
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

			                // add to hash-table
			                P_table->add(pkey, pnew);
			                num_particle++;
			                TKey tmpkey(pkey);
			                Curr_buck->add_real_particle(tmpkey);
		        		}

	    }//finish bucket is overground

	    else if ((Bnd_buck->get_bucket_type () == UNDERGROUND) || (Bnd_buck->get_bucket_type () == PRESS_BC) )
	    {

	    	if (Bnd_buck->get_bucket_type () == UNDERGROUND)
	    		bctp = 2;  //wall bc
	    	else
	    		bctp = 1;   //pressure bc

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
#ifdef DEBUG
		                 if (do_search)
		                 {
		                     if (find_particle (pcrd, check))
			                     cout << "The particle found!" << endl;
		                 }
#endif
	    				 for (ii = 0; ii < DIMENSION; ii++)
	    				      normc[ii] = (pcrd[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

//	    				 HSFC3d (normc, &keylen, key);
	    				 THSFC3d (normc, add_step, &tkeylen, pkey);

	    				 // check for duplicates
	    				 if (P_table->lookup(pkey))
	    				 {
	    				     fprintf(stderr, "ERROR: Trying to add particle "
	    				                     "twice on same location.\n");
	    				     exit(1);
	    				 }
	    				 Particle * pnew = new Particle(pkey, pcrd, mass, smlen, prss, masfrc, gmm, sndspd, phs_num, myid, bctp);

	    				 // add to hash-table
	    				 P_table->add(pkey, pnew);
	    				 num_particle++;
	    				 TKey tmpkey(pkey);
	    				 switch (bctp)
  					     {
  					     case 1 :
  					         Curr_buck->add_pressure_ghost_particle(tmpkey);
  					         break;
  					     case 2 :
  						     Curr_buck->add_wall_ghost_particle(tmpkey);
  						     break;
  					     default:
  						     cout << "bctp incorrect in function add_air!\n" << endl;
  					     }
	    			}

	    }//finish bucket is pressure bc or underground
	    else
	    	cout << "Invalid bucket type! The bucket is none of OVERGROUND, PRESSURE_BC, UNDERGROUND, MIXED" <<endl;
	  }//finish go though all buckets

	  int ntotal = 0;

	#ifdef MULTI_PROC
	  MPI_Allreduce (&num_particle, &ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  if ( ntotal < 100 )
	  {
	    if ( myid == 0 )
	    {
	      fprintf (stderr, "No. of particles = %d\n", ntotal);
	      fprintf (stderr, "Not enough particles.\n");
	    }
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	#else
	  if ( num_particle < 100 )
	  {
	    fprintf (stderr, "No. of particles = %d\n", ntotal);
	    fprintf (stderr, "Not enough particles.\n");
	    exit (1);
	  }
	#endif

	  //  mark all neighbors with real particles active
	  itr->reset();
	  while ((Curr_buck = (Bucket *) itr->next()))
	  {
	    Curr_buck->mark_inactive ();
	    Key *nkey = Curr_buck->get_neighbors ();
	    for (i = 0; i < NEIGH_SIZE; i++)
	      // the guest buckets are not moved yet,
	      // hence only check buckets which belong to current proc
	      if (*(Curr_buck->get_neigh_proc () + i) == myid )
	      {
	        Bucket *neigh = (Bucket *) BG_mesh->lookup (nkey[i]);
	        if (neigh->has_real_particles ())
	        {
	          Curr_buck->mark_active();
	          break;
	         }
	       }
	    }

	  delete itr;

	  return;
}
