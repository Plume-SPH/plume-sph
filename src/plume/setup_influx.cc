/*
 * setup_influx.cc
 *
 *  Created on: Oct 27, 2016
 *      Author: zhixuanc
 */

//For setting up influx boundary condition for ash dispersion and transportation simulation
#include <math.h>  //sqrt
#include <cassert>

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <mpi.h>
#include <hilbert.h>
#include <particle.h>
#include <multiproc.h>
#include <properties.h>

#include <Involved_header.h>

#include "constant.h"
#include "sph_header.h"
#include "parameters.h"

using namespace std;

#if CODE_DIMENSION==3

#ifdef SIMULATE_ASH

int
setup_influx(int myid, THashTable * P_table, HashTable * BG_mesh,
        TimeProps * timeprops, MatProps* matprops, SimProps *simprops, int numprocs )
{
	double crd_p[DIMENSION];
	double normc[DIMENSION];
	unsigned key[TKEYLENGTH];
	unsigned tempkey[TKEYLENGTH];
	unsigned buckkey[KEYLENGTH];
	bool erpt;
	int tempid;
	int involved;
    int i, j, k, l;
    unsigned tkeylen = TKEYLENGTH;

#if FLUID_COMPRESSIBILITY==0 //using EOS of ideal gas
	  double sndspd = 340.;
#elif FLUID_COMPRESSIBILITY==1
	  double sndspd = 1482.;
#endif

    double dist, dist_sqrt;
    double velx, vely;
    double rinsq=r_in_P*r_in_P;
    double routsq=r_out_P*r_out_P;

    unsigned add_step = 0;
    int bctp ;

	//determine the mass  and sml of influx particle
	double mss, sml;
	double des = rhov_P;
	sml=matprops->smoothing_length;
	mss=Compute_mass (sml, des);

	simprops->mass_of_phase2 = mss;
	simprops->sml_of_phase2 = sml;

	//determine the rough range of influx buck
	double range_x[2];
	double range_y[2];
	double range_z[2];

	range_x[0] = -r_out_P;
	range_x[1] = r_out_P;
	range_y[0] = -r_out_P;
	range_y[1] = r_out_P;
    range_z[1] = h_top_P;
    range_z[0] = Ll_P[2];  //+ PARTICLE_DENSITY*0.5*sml;


    //Find all buckets as influx buckets
    HTIterator * itr = new HTIterator (BG_mesh);
    Bucket * Curr_buck;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

    double mincrd[DIMENSION], maxcrd[DIMENSION];

    /*go through all erupt buckets and
     * 1) delete all non-erupt ghost particles that within r_out
     * 2) mark bucket as influx
     */
    vector < TKey > pnew; //particle key
    vector < TKey > plist;
    vector < TKey >::iterator p_itr;
    Particle *pj;

    while ((tempptr=itr->next ()))
    {
    	breif_buck = (BriefBucket *) tempptr;
    	if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
    		continue;
    	else
    	{
    		Curr_buck = (Bucket*) tempptr;
        	/* Again --->even guest bucket need to be checked,
        	 * ---> so that the synchronization is not needed
        	 * */
    	    for (i = 0; i < DIMENSION; i++)
    	    {
    	        mincrd[i] = *(Curr_buck->get_mincrd () + i);
    	        maxcrd[i] = *(Curr_buck->get_maxcrd () + i);
    	    }

    	    //determine whether some portion of the bucket included in the influx particle range
    	    erpt = determine_erupt_buket (mincrd, maxcrd, range_x, range_y, range_z); // any bucket within r_out will be marked as influx buckets

    	    if (erpt)
    	    {
    	    	  //Mark particles
    	    	  Curr_buck->mark_influx();

    	    	  //set particles type to be 0 ---> Need double check, for these buckets that has real particle, is this correct? ---> Yes, because lately, we re-set the particles_type based on plist
				  Curr_buck->put_particles_type (0);

				  //set has no involved
				  Curr_buck->set_has_no_involved();

				  //check all particle in the erupt bucket and turn them to pressure ghost necessary!
				  pnew.clear();
				  plist = Curr_buck->get_particle_list ();
				  for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
				  {
					  pj = (Particle *) P_table->lookup(*p_itr);
					  assert (pj);
					  for (k=0; k<DIMENSION; k++)
						  crd_p[k] = *(pj->get_coords() + k);

					  dist = 0;
					  for (j=0; j<2; j++)
						  dist += (crd_p[j]*crd_p[j]);

					  if (dist <= routsq && crd_p[2] < h_top_P && crd_p[2] > range_z[0])
					  {
						  if (dist > rinsq && crd_p[2] > h_bot_P)
						  {
							  for (k = 0; k<TKEYLENGTH; k++)
								  tempkey[k] = pj->getKey().key[k];

							  /*
							   * Here what I did is only remove them from P_table and bucket particle list on this processor.
							   * But the particle in other particles neighbour list is not updated!---> will be done in search neighbor (please double check)
							   * And the particle as guest on other processes is not removed!--->well, all particle on other processes will be syn as other processes will execute the same command. ----> you should have already been aware of the fact that no "non_guest" constrain is imposed here.
							   * And the particle, if it has image, the image should also be removed! ---> will be done in search_image and imposing BC (please double check----> yes, I checked that, no problem)
							   */

							  P_table->remove(tempkey); //remove from the hashtable
						  }
						  //maybe syn after adding particles will always a better idea.
						  else
						  {
							  pj->put_bc_type(1); //switch the particle type to be pressure ghost
							  Curr_buck->set_pressure_ghost_particles(true);
							  pj->set_involved_flag (0); //switch the particle type to be not involved
							  double temvel[DIMENSION]={0.0, 0.0, 0.0};
							  pj->put_vel(temvel); //make sure they are stationary

							  pnew.push_back(*p_itr);
							  Curr_buck->set_pressure_ghost_particles(true);

							  //change the phase num to 2, ---> this will be used to make sure that these particles in the eruption area will not be turn to involved and real.
							  pj->set_phs_num(2);
						  }

					  }

					  else
					  {
						  pnew.push_back(*p_itr);//remove from the from particle list!
						  bctp = pj->get_bc_type ();
						  switch (bctp)
						  {
							  case 0 :
								  Curr_buck->set_erupt_ghost_particles(true);
								  break;
							  case 2 :
								  Curr_buck->set_wall_ghost_particles(true);
								  break;
							  case 100 :
								  Curr_buck->set_real_particles(true);
								  break;
							  case 1 :
								  Curr_buck->set_pressure_ghost_particles(true);
								  break;
							  default:
								  cout << "bctp incorrect in function set_up_influx!\n" << endl;
						  }

						  involved = pj->get_involved ();
						  switch (involved)
						  {
							  case 2 :
								  Curr_buck->set_has_involved (true);
								  break;
							  case 1 :
								  Curr_buck->set_has_potential_involved(true);
								  break;
							  case 0 : //0 means that the particle is not involved. So will not change the has_involved flag
								  break;
							  default:
								  cout << "involved incorrect in function set_up_influx!\n" << endl;
						  }
					  }// end of else particle is out of the indlux domain
				  }//end of loop go through all particles in the buckets
				  Curr_buck-> put_new_plist(pnew);
				  Curr_buck-> update_particles ();
    	    }//end of if bucket is erupt
    	}//end of if bucket is not brief
    }// end of loop go through all buckets


	// get min-max domain from hashtable, for key generation
    double mindom[DIMENSION], maxdom[DIMENSION];
	for (i = 0; i < DIMENSION; i++)
	{
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	}

	double r_add = r_in_P + 0.99 * sml; //The larger radius of particle adding ring
	double raddsq=r_add*r_add;
    //adding influx particles
	double sml2 = 0.5*sml;

    itr->reset();
    while ((tempptr=itr->next ()))
    {
    	breif_buck = (BriefBucket *) tempptr;
    	if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
    		continue;
    	else
    	{
    		  Curr_buck = (Bucket*) tempptr;

    	      if (Curr_buck->is_influx ())
    	      {
    	    	  for (l = 0; l < DIMENSION; l++)
    	    	  {
    	    	      mincrd[l] = *(Curr_buck->get_mincrd () + l);
    	    	      maxcrd[l] = *(Curr_buck->get_maxcrd () + l);
    	    	  }
    	    	  for (i = 0; i<PARTICLE_DENSITY; i++)
    	    	  {
    	    		  crd_p[0]=mincrd[0] + sml2 + i*sml;
    	    		  for (j = 0; j<PARTICLE_DENSITY; j++)
    	    		  {
    	    			  crd_p[1]=mincrd[1]+sml2+j*sml;

    	    			  dist =0.0;
    	    			  for (l = 0; l<2; l++)
    	    				  dist += crd_p[l]*crd_p[l];

						  if (dist<routsq && dist>rinsq)
						  {
							  dist_sqrt = sqrt(dist);

							  velx=crd_p[0]/dist_sqrt*vel0_P;
							  vely=crd_p[1]/dist_sqrt*vel0_P;

							  for (k=0; k<PARTICLE_DENSITY; k++)
							  {
								  crd_p[2] = mincrd[2] + sml2 + k*sml;

								  if (crd_p[2] < h_top_P && crd_p[2] > h_bot_P)
								  {
									  for (l = 0; l < DIMENSION; l++)
									       normc[l] = (crd_p[l] - mindom[l]) /(maxdom[l] - mindom[l]);
									  THSFC3d (normc, add_step, &tkeylen, key);
									  Particle * pnew = new Particle (key, crd_p, mss, sml, myid, msfc0_P,
									      		velx, vely, ev0_P, rhov_P,  pv0_P,  gamma_v_P,  sndspd,
									      		ng0_P,  Cvs_P,  Cvg_P,  Cva_P,  Rg_P,  Ra_P);

			    	    			  /*
			    	    			   * Here what I did is only add them into P_table and bucket particle list!
			    	    			   * But the particle in other particles neighbour list is not added!---> will be done in search neighbor (please double check)
			    	    			   * And the particle as guest on other processes is not removed!--->well, all particle on other processes will be added as other processes will execute the same command. ----> you should have already been aware of the fact that no "non_guest" constrain is imposed here.
			    	    			   * And the particle, if it has image, the image should also be added! ---> will be done in search_image and imposing BC (please double check----> yes, I checked that, no problem)
			    	    			   */
									  // add to hash-table
									  P_table->add(key, pnew);
									  TKey pkey (key);
									  Curr_buck->add_erupt_ghost_particle (pkey);

								      //default involved is zero---> need to double check what involved should be allocated here-->No problem, should be not involved.
								      if (Curr_buck->is_guest())
								      {
								    	  pnew->put_guest_flag(true);
										  tempid = Curr_buck->get_myprocess ();
										  pnew->put_my_processor(tempid);
								      }

									  //The next step is for find the particle adding position.
									  if (dist < raddsq)
									  {
										  InfluxAddingPos addpos (crd_p);
										  Curr_buck->add_adding_pos(addpos);
									  }
								   }
							    }

						     }//end of if dist<routsq && dist>rinsq

    	    		  }//end of loop in y dir

    	    	  }//end of loop in x dir

    	      }//end of if bucket is influx

    	}//end of if bucket is not brief
    }//end of while


    //clear up:
    delete itr;
}


int
add_new_influx(int myid, THashTable * P_table, HashTable * BG_mesh,
        TimeProps * timeprops, MatProps* matprops, SimProps *simprops, int numprocs)
{
	int i, j, k, l;
	bool need_add;
	double crd_p[DIMENSION], crd_n[DIMENSION];

	double normc[DIMENSION];
	unsigned key[TKEYLENGTH];
    unsigned tkeylen = TKEYLENGTH;
    double dist, dist_sqrt, hsq;
    double sml = simprops->sml_of_phase2; //influx particles are phase 2
    hsq = sml*sml*0.999; //0.999 is used to avoid repeat adding

    Key *neighbors;//bucket key

	unsigned add_step;
    //the more robust way: use time step as add_step
	add_step = timeprops->step;
	double mss = simprops->mass_of_phase2;

	// get min-max domain from hashtable, for key generation
    double mindom[DIMENSION], maxdom[DIMENSION];
	for (i = 0; i < DIMENSION; i++)
	{
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	}

#if FLUID_COMPRESSIBILITY==0 //using EOS of ideal gas
	 double sndspd = 340.;
#elif FLUID_COMPRESSIBILITY==1
	 double sndspd = 1482.;
#endif

    double velx, vely;

	vector <InfluxAddingPos>::iterator psit;

	vector < TKey > plist;
	vector < TKey >::iterator p_itr;
	Particle *pj;

    //Find all buckets as influx buckets
    HTIterator * itr = new HTIterator (BG_mesh);
    Bucket * Curr_buck;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

	//We do not add influx particles for guest bucket ---> will syn after adding!
    itr->reset();
    while ((tempptr=itr->next ()))
    {
    	breif_buck = (BriefBucket *) tempptr;
    	if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
    		continue;
    	else
    	{
    		  Curr_buck = (Bucket*) tempptr;
    		  if (Curr_buck->has_adding_pos () && ! Curr_buck->is_guest ())
    		  {
    			  plist = Curr_buck->get_particle_list ();
			      neighbors = Curr_buck->get_neighbors();
			      const int *neigh_proc = Curr_buck->get_neigh_proc();

    			  vector <InfluxAddingPos> Alist =  Curr_buck->get_adding_list ();
    			  for (psit = Alist.begin(); psit != Alist.end(); psit++)
    			  {
    				    need_add =true;

    					for (l=0; l<DIMENSION; l++)
    						crd_p[l]=psit->crd[l];


    					for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
    					{
    						  pj = (Particle *) P_table->lookup(*p_itr);
    			              if ( ! pj )
    			              {
    			                  fprintf (stderr, "Particle: {%u, %u, %u} missing on %d.\n",
    			                           p_itr->key[0], p_itr->key[1], p_itr->key[2], myid);
    			                  fflush (stderr);
    			              }
    						  assert (pj);

    						  if (pj->is_erupt_ghost ()) // only the erupt ghost particles should be considered
    						  {
    							  for (k=0; k<DIMENSION; k++)
    								  crd_n[k] = *(pj->get_coords() + k);

    							  if (crd_p[2] == crd_n[2]) //Because these non involved influx particles move only horizontally, so they are well-layered
    							  {
    								  dist = 0.0;
    								  for (j=0; j<2; j++)
    									  dist += ( (crd_n[j]-crd_p[j])*(crd_n[j]-crd_p[j]) );

    								  if (dist < hsq)
    									  need_add = false;
    							  }
    						  }

    					}// end of go through all particles in the buckets.

    					//Need also go through neighbor bucket if the adding position are at the most outside layer of Curr_buck

    			        for (i = 0; i < NEIGH_SIZE; i++)
    			            if (neigh_proc[i] > -1)
    			            {
    			               Bucket * neigh_buck = (Bucket *) BG_mesh->lookup (neighbors[i]);

    			               // if neighbor is not found and it belongs to
    			               // different process, it is not needed as it
    			               // doesn't have any particles.
    			               if ((! neigh_buck) && (neigh_proc[i] != myid))
    			                   continue;
    			               else
    			                   assert(neigh_buck);

    			               vector < TKey > plist2 = neigh_buck->get_plist();
    			               vector < TKey >::iterator it2;
    			               for (it2 = plist2.begin(); it2 != plist2.end(); it2++)
    			               {
    			                    pj = (Particle *) P_table->lookup(*it2);
    			                    if ( ! pj )
    			                    {
    			                        fprintf (stderr, "Particle: {%u, %u, %u} missing on %d.\n",
    			                        it2->key[0], it2->key[1], it2->key[2],myid);
    			                        fflush (stderr);
    			                    }
    			                    assert(pj);
    			      			  if (pj->is_erupt_ghost ()) // only the erupt ghost particles should be considered
    			      			  {
    			      			      for (k=0; k<DIMENSION; k++)
    			      				      crd_n[k] = *(pj->get_coords() + k);

    			    				  if (crd_p[2] == crd_n[2])
    			    				  {
    			    					  dist = 0.0;
    			    					  for (j=0; j<2; j++)
    			    						  dist += ( (crd_n[j]-crd_p[j])*(crd_n[j]-crd_p[j]) );

    			    					  if (dist < hsq)
    			    						  need_add = false;
    			    				  }
    			      			  }
    			               }
    			            }//end of if bucket is not itself

    					if (need_add)
    					{
    						  dist_sqrt = sqrt(crd_p[0]*crd_p[0]+crd_p[1]*crd_p[1]);
    						  velx=crd_p[0]/dist_sqrt*vel0_P;
    						  vely=crd_p[1]/dist_sqrt*vel0_P;
    						  for (l = 0; l < DIMENSION; l++)
    							  normc[l] = (crd_p[l] - mindom[l]) /(maxdom[l] - mindom[l]);

    						  THSFC3d (normc, add_step, &tkeylen, key);
    						  Particle * pnew = new Particle (key, crd_p, mss, sml, myid, msfc0_P,
    								  velx, vely, ev0_P, rhov_P,  pv0_P,  gamma_v_P,  sndspd,
    								  ng0_P,  Cvs_P,  Cvg_P,  Cva_P,  Rg_P,  Ra_P);

    						  /*
    						   * Here what I did is only add them into P_table and bucket particle list!
    						   * But the particle in other particles neighbour list is not added!---> will be done in search neighbor (please double check)
    						   * And the particle as guest on other processes is not removed!--->well, all particle on other processes will be added as other processes will execute the same command. ----> you should have already been aware of the fact that no "non_guest" constrain is imposed here.
    						   * And the particle, if it has image, the image should also be added! ---> will be done in search_image and imposing BC (please double check----> yes, I checked that, no problem)
    						   */
    						  // add to hash-table
    						  P_table->add(key, pnew);
    						  TKey pkey (key);
    						  Curr_buck->add_erupt_ghost_particle (pkey);
    					}
    			  }// end of go through the addinglist
    		  }
    	}// end of if bucket is not brief
    }// end of go through all bucket


    //clear up:
    delete itr;
}


#endif //SIMULATE_ASH

#endif  //#if CODE_DIMENSION==3
