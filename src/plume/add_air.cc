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
#include <cmath> //sin
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

#if CODE_DIMENSION==1
void
set_up_shock_tube (THashTable * P_table, MatProps * matprops, SimProps* simprops, int numproc, int myid)
{
	  int i, j, k, ii;
	  int tempid;
	  unsigned pkey[TKEYLENGTH];
	  unsigned tkeylen = TKEYLENGTH;
	  double pcrd[DIMENSION], bnd[DIMENSION*2];
	  double mass, des;
	  vector < TKey > neighs;

	  // direction indices on upper bucket
	  int num_particle = 0;

	  // start putting pile

	  double smlen = matprops->smoothing_length;

	  //Initial particle property
	  int bctp_real = 100; //boundary condition ghost particle type.
	  int bctp_prss = 0; //boundary condition ghost particle type.
	  int bctp_nsrp = 99; //no_save real particles

	  class Shock_Inputs
	  {
	  private:
		  double prss_l;  //Parameter
		  double prss_r ;  //Parameter
		  double des_l ;  //Parameter
		  double des_r ;  //Parameter
		  double vel_l ;  //Parameter
		  double vel_r ;  //Parameter
		  double middle_point;
		  double A;  //amplitude for density of Shu-Osher
		  double f; //frequency for density of Shu-Osher
	  public:
		  Shock_Inputs(); //The default contructor ---> Use default value
		  //basic constructor for constant input
		  Shock_Inputs(double dl, double pl, double ul, double dr, double pr, double ur, double mid)
		  {
			 prss_l = pl;  //Parameter
			 prss_r = pr;  //Parameter
			 des_l = dl;  //Parameter
			 des_r = dr;  //Parameter
			 vel_l = ul;  //Parameter
			 vel_r = ur;  //Parameter
			 middle_point=mid;
			 A=0.0;  //amplitude for density of Shu-Osher
			 f=1.0; //frequency for density of Shu-Osher
		  };
		  //constructor for Shu-Osher ---> Density is a function of location
		  Shock_Inputs(double dl, double pl, double ul, double dr, double pr, double ur, double mid, double amp, double fre)
		  {
			  prss_l = pl;  //Parameter
			  prss_r = pr;  //Parameter
			  des_l = dl;  //Parameter
			  des_r = dr;  //Parameter
			  vel_l = ul;  //Parameter
			  vel_r = ur;  //Parameter
			  middle_point=mid;
			  A=amp; //amplitude for density of Shu-Osher
			  f=fre; //frequency for density of Shu-Osher
		  };
		  double get_rho(double x)
		  {
			  double dl=des_l;
			  double dr=des_r+A*sin(f*x);
			  return (x > middle_point ? dr : dl);
		  };

		  double get_press (double x)
		  {
			  return (x > middle_point ? prss_r : prss_l);
		  };

		  double get_vel (double x)
		  {
			  return (x > middle_point ? vel_r : vel_l);
		  };

	  };

#if SHOCK_TUBE_TESTS==0
	  //Sod shock tube input ---Initial GSPH test parameters
	  double prss_l = 1.0;  //Parameter
	  double prss_r = 0.2;  //Parameter
	  double des_l = 1.0;  //Parameter
	  double des_r = 0.5;  //Parameter
	  double vel_l = 0.0;  //Parameter
	  double vel_r = 0.0;  //Parameter
	  double middle_point=0.;
	  Shock_Inputs SIPT (des_l, prss_l, vel_l, des_r, prss_r, vel_r, middle_point);
#elif SHOCK_TUBE_TESTS==1
	  //Sod shock tube input ---Famous test case in FV
	  double prss_l = 1.0;  //Parameter
	  double prss_r = 0.1;  //Parameter
	  double des_l = 1.0;  //Parameter
	  double des_r = 0.125;  //Parameter
	  double vel_l = 0.0;  //Parameter
	  double vel_r = 0.0;  //Parameter
	  double middle_point=0.;
	  Shock_Inputs SIPT (des_l, prss_l, vel_l, des_r, prss_r, vel_r, middle_point);
#elif SHOCK_TUBE_TESTS==2
	  //Shu-Osher problem ---Famous test case in FV
	  double prss_l = 10.33333;  //Parameter
	  double prss_r = 1.0;  //Parameter
	  double des_l = 3.857143;  //Parameter
	  double des_r = 1.0;  //---> Is a function of
	  double vel_l = 2.629369;  //Parameter
	  double vel_r = 0.0;  //Parameter
	  double middle_point=-4.0;
	  double amp = 0.2;
	  double fre = 5.0;
	  Shock_Inputs SIPT (des_l, prss_l, vel_l, des_r, prss_r, vel_r, middle_point, amp, fre);
#endif

	  double dx_l=smlen;

#if EQUAL_PART_MASS==1
	  double dx_r = dx_l*des_l/des_r;
#else
	  double dx_r = smlen;
#endif

	  //When there is no adaptive smoothing length, it it necessary to set the smoothing length to a larger value so that each kernal has enough particles within its support.
	  smlen=max(dx_r, dx_l);

	  double dx2_l = 0.5 * dx_l;
	  double dx2_r = 0.5 * dx_r;

#if FLUID_COMPRESSIBILITY==0 //using EOS of ideal gas
	  double sndspd = 340.;
	  double gmm = 1.4;
#elif FLUID_COMPRESSIBILITY==1
	  double sndspd = 1482.;
	  double gmm = 7.0;
#endif

	  unsigned add_step = 0; // eruption particle adding step
	  int involved = 2; //involved
	  int not_involved = 0; //involved


	  /*temporarily use the following way to get bnd, Finally, I will put bnd either
	   * in SimulProps or parameters.h
	   * will modify this part later
	  */
	  for (i = 0; i < DIMENSION; i++)
	  {
	     bnd[i*2]=Ll_P[i];
	     bnd[i*2+1]=Lu_P[i];
	  }

	for (k = 0; k < Nnsrp_l_P; k++)
	{
		pcrd[0] = bnd[0] - ( dx2_l + k*dx_l);
		pkey[0]=num_particle;
		pkey[1]=0;
		pkey[2]=add_step;

		des=SIPT.get_rho(pcrd[0]);
	    mass=des* dx_l;
		// check for duplicates
		if (P_table->lookup(pkey))
		{
		  fprintf(stderr, "ERROR: Trying to add particle "
						  "twice on same location.\n");
		  exit(1);
		}

		Particle * pnew = new Particle(pkey, pcrd, mass, smlen , des, vel_l, prss_l, gmm, sndspd, bctp_nsrp, not_involved);
		// add to hash-table
		P_table->add(pkey, pnew);
		num_particle++;
		neighs.push_back(pkey);
#if DENSITY_UPDATE_SML==0
		pnew->put_smlen_original(dx_l);
#endif
	}

	for (k = 0; k < Nb_P; k++)
	{
		pcrd[0] -= dx_l;
		pkey[0]=num_particle;
		pkey[1]=0;
		pkey[2]=add_step;

		des=SIPT.get_rho(pcrd[0]);
	    mass=des* dx_l;
		// check for duplicates
		if (P_table->lookup(pkey))
		{
		  fprintf(stderr, "ERROR: Trying to add particle "
						  "twice on same location.\n");
		  exit(1);
		}

		Particle * pnew = new Particle(pkey, pcrd, mass, smlen , des, vel_l, prss_l, gmm, sndspd, bctp_prss, not_involved);
		// add to hash-table
		P_table->add(pkey, pnew);
		num_particle++;
		neighs.push_back(pkey);
#if DENSITY_UPDATE_SML==0
		pnew->put_smlen_original(dx_l);
#endif
	}


    pcrd[0]=bnd[0]+dx2_l;
	while (pcrd[0]<=middle_point)
	{
		pkey[0]=num_particle;
		pkey[1]=0;
		pkey[2]=add_step;

		des=SIPT.get_rho(pcrd[0]);
	    mass=des* dx_l;
		// check for duplicates
		if (P_table->lookup(pkey))
		{
		  fprintf(stderr, "ERROR: Trying to add particle "
						  "twice on same location.\n");
		  exit(1);
		}

		Particle * pnew = new Particle(pkey, pcrd, mass, smlen , des, vel_l, prss_l, gmm, sndspd, bctp_real, involved);
		// add to hash-table
		P_table->add(pkey, pnew);
		num_particle++;
		pcrd[0] += dx_l;
		neighs.push_back(pkey);
#if DENSITY_UPDATE_SML==0
		pnew->put_smlen_original(dx_l);
#endif
	}

	pcrd[0] = pcrd[0] - dx_l+ dx_r;
	while (pcrd[0]<=bnd[1])
	{
		pkey[0]=num_particle;
		pkey[1]=0;
		pkey[2]=add_step;

		des=SIPT.get_rho(pcrd[0]);
	    mass=des* dx_r;
		// check for duplicates
		if (P_table->lookup(pkey))
		{
		  fprintf(stderr, "ERROR: Trying to add particle "
						  "twice on same location.\n");
		  exit(1);
		}

		Particle * pnew = new Particle(pkey, pcrd, mass, smlen, des, vel_r, prss_r, gmm, sndspd, bctp_real, involved);
		// add to hash-table
		P_table->add(pkey, pnew);
		num_particle++;
		pcrd[0] += dx_r;
		neighs.push_back(pkey);
#if DENSITY_UPDATE_SML==0
		pnew->put_smlen_original(dx_r);
#endif
	}

	for (k = 0; k < Nnsrp_r_P; k++)
	{
		pkey[0]=num_particle;
		pkey[1]=0;
		pkey[2]=add_step;

		des=SIPT.get_rho(pcrd[0]);
	    mass=des* dx_r;
		// check for duplicates
		if (P_table->lookup(pkey))
		{
		  fprintf(stderr, "ERROR: Trying to add particle "
						  "twice on same location.\n");
		  exit(1);
		}

		Particle * pnew = new Particle(pkey, pcrd, mass, smlen , des, vel_r, prss_r, gmm, sndspd, bctp_nsrp, involved);
		// add to hash-table
		P_table->add(pkey, pnew);
		num_particle++;
		pcrd[0] += dx_r;
		neighs.push_back(pkey);
#if DENSITY_UPDATE_SML==0
		pnew->put_smlen_original(dx_r);
#endif
	}

	for (k = 0; k < Nb_P; k++)
	{
		pkey[0]=num_particle;
		pkey[1]=0;
		pkey[2]=add_step;

		des=SIPT.get_rho(pcrd[0]);
	    mass=des* dx_r;
		// check for duplicates
		if (P_table->lookup(pkey))
		{
		  fprintf(stderr, "ERROR: Trying to add particle "
						  "twice on same location.\n");
		  exit(1);
		}

		Particle * pnew = new Particle(pkey, pcrd, mass, smlen , des, vel_r, prss_r, gmm, sndspd, bctp_prss, involved);
		// add to hash-table
		P_table->add(pkey, pnew);
		num_particle++;
		pcrd[0] += dx_r;
		neighs.push_back(pkey);
#if DENSITY_UPDATE_SML==0
		pnew->put_smlen_original(dx_r);
#endif
	}

	cout << "The total number of real particles is:" << num_particle-1-2*Nb_P << endl;

	// updating neighbor list
	THTIterator * itr = new THTIterator (P_table);
	Particle * pi = NULL;
	while ((pi = (Particle *) itr->next ()))
	{
	    if (pi->need_neigh())
	    {
	    	pi->put_neighs(neighs);
	    }
	}

	//clean up
	delete itr;

	return;
}

//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------

#elif CODE_DIMENSION==3
void
add_air (THashTable * P_table, HashTable * BG_mesh,
        MatProps * matprops, SimProps* simprops, int numproc, int myid)
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
     unsigned keycheck[TKEYLENGTH] = {268485399, 978591, 0};
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
		 					           initial_air (pnew, simprops);

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
		 					           initial_air (pnew, simprops);

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
							         initial_air (pnew, simprops);

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

	  delete itr;

	  return;
}
#endif  //CODE_DIMENSION==3
