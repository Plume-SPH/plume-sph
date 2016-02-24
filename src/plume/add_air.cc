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
	  int tempid;
	  unsigned pkey[TKEYLENGTH];
	  unsigned tkeylen = TKEYLENGTH;
	  double mincrd[DIMENSION], maxcrd[DIMENSION];
	  double mindom[DIMENSION], maxdom[DIMENSION];
	  double normc[DIMENSION];
	  double pcrd[DIMENSION];
	  double bnd[2*DIMENSION], index[2*DIMENSION];

	  Bucket *Curr_buck = NULL;

	  // direction indices on upper bucket
	  int num_particle = 0;

	  // start putting piles
	  double mass = matprops->particle_mass;
	  double smlen = matprops->smoothing_length;
	  double dx = smlen;
	  double dx2 = 0.5 * dx;

#ifdef DEBUG
     bool do_search = false;
     bool search_by_coor=false;
     double check[DIMENSION] = {-250, -250, 4250};
     unsigned keycheck[TKEYLENGTH] = {70540343, 676583194, 0};
     unsigned keytemp[TKEYLENGTH];
#endif

	  // get min-max domain from hashtable, for key generation
	  for (i = 0; i < DIMENSION; i++)
	  {
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	  }

	  //Initial particle property
	  int bctp = 100; //boundary condition ghost particle type.
	  int ptype = 0; //real particle
	  double prss = 0.;
	  double masfrc = 0.;
	  double gmm = 1.4;
	  double sndspd = 340.;
	  int phs_num = 1; //phase 1, air particle
	  unsigned add_step = 0; // eruption particle adding step
	  int involved = 1; //potential involved

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
       bnd[0]=Lx_P[0];
       bnd[1]=Lx_P[1];
       bnd[2]=Ly_P[0];
       bnd[3]=Ly_P[1];
       bnd[4]=Lz_P[0];
       bnd[5]=Lz_P[1];

	  while ((tempptr=itr->next ()))
	  {
		  breif_buck = (BriefBucket *) tempptr;
		  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
			  continue;
		  else
		  {
			  Bnd_buck = (Bucket*) tempptr;
		      if (Bnd_buck->get_has_involved() && !((Bnd_buck->get_plist ()).size()))  // make sure bucket is has_involved >0 ; is not added yet;
			  {
		    	Bnd_buck->mark_active ();
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

				                   if ((Bnd_buck->get_bucket_index())[4] == -1) //on ground MIXED
				                   {
				                	 if (pcrd[2]>bnd[4]) //This is only for on-ground MIXED, for Mixed in the air, the bnd[4]=0.1, which does not make any sense. for over-ground MIXED and PRESSURE_BC, pressure ghost need to be full of the bucket.
				                	 {
			    				       for (ii = 0; ii < DIMENSION; ii++)
			    				           normc[ii] = (pcrd[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

			    					   THSFC3d (normc, add_step, &tkeylen, pkey);

		#ifdef DEBUG
				if (do_search)
				{
				    for (ii = 0; ii< TKEYLENGTH; ii++)
					    keytemp[ii] = pkey[ii];

				    if (find_particle (keytemp, keycheck))
					    cout << "The particle found!" << endl;
				}
		#endif

				    				   // check for duplicates
				    				   if (P_table->lookup(pkey))
				    				   {
				    				       fprintf(stderr, "ERROR: Trying to add particle "
				    				               "twice on same location.\n");
				    				       exit(1);
				    				   }
				    				   Particle * pnew = new Particle(pkey, pcrd, mass, smlen, prss, masfrc, gmm, sndspd, phs_num, myid, bctp);
							           pnew->set_involved_flag(involved);

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
				                   }//end of if bucket is on ground MIXED
				                   else // bucket is over ground MIXED
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
							           pnew->set_involved_flag(involved);

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
			    			 }

			    }//finish bucket is Mixed
			    else if (Bnd_buck->get_bucket_type () == OVERGROUND)
			    {
		//	    	bctp = 100;
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
				                   if (search_by_coor)
				                   {
				                        if (find_particle (pcrd, check))
					                    cout << "The particle found!" << endl;
				                    }
		#endif
					                for (ii = 0; ii < DIMENSION; ii++)
					                  normc[ii] = (pcrd[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

					                THSFC3d (normc, add_step, &tkeylen, pkey);
		#ifdef DEBUG
				if (do_search)
				{
				    for (ii = 0; ii< TKEYLENGTH; ii++)
					    keytemp[ii] = pkey[ii];

				    if (find_particle (keytemp, keycheck))
					    cout << "The particle found!" << endl;
				}
		#endif
					                // check for duplicates
					                if (P_table->lookup(pkey))
					                {
					                  fprintf(stderr, "ERROR: Trying to add particle "
					                                  "twice on same location.\n");
					                  exit(1);
					                }
					                Particle * pnew = new Particle(pkey, pcrd, mass, smlen, prss, masfrc, gmm, sndspd, phs_num, myid, bctp);
						            pnew->set_involved_flag(involved);

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
					                Curr_buck->add_real_particle(tmpkey);
				        		}

			    }//finish bucket is overground
			    else if (Bnd_buck->get_bucket_type () == PRESS_BC)
			    	cout << "The initial domain is larger than the simulation domain, please double check the input domain parameters!" << endl;
			    else if (Bnd_buck->get_bucket_type () == UNDERGROUND)
			    	cout << "Wow! you marked underground bucket as has_potential_involved, something was wrong!" << endl;
			  }//finish if
		  }//end if bucket is no brief bucket
	  }//finish  go through all buckets

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

	  //The following is actually not necessary ---> we have already done this in the read_data
//	  //  mark all neighbors with real particles active
//	  Bucket * neigh = NULL;
//	  itr->reset();
//	  while ((Curr_buck = (Bucket *) itr->next()))
//	  {
//	    Curr_buck->mark_inactive ();
//	    Key *nkey = Curr_buck->get_neighbors ();
//	    for (i = 0; i < NEIGH_SIZE; i++)
//	      // Particle is also added on guest bucket
//	      // So it is not necessary to check bucket only on current process
//	      if (*(Curr_buck->get_neigh_proc () + i) == myid )
//	      {
//	        neigh = (Bucket *) BG_mesh->lookup (nkey[i]);
//	        // newly added, make sure neigh exist!
//	        assert(neigh);
//
//	        if (neigh->has_real_particles ())
//	        {
//	          Curr_buck->mark_active();
//	          break;
//	        }
//	      }
//	    }

	  delete itr;

	  return;
}
