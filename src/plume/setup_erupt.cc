/*
 * setup_erupt.cc
 *
 *  Created on: Apr 27, 2015
 *      Author: zhixuanc
 */

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <mpi.h>
#include <hilbert.h>
#include <particle.h>
#include <multiproc.h>
#include <properties.h>

#include "particler.h"
#include "constant.h"
#include "sph_header.h"
#include "parameters.h"

//function that used to determine the type of bucket
bool determine_erupt_buket (double *mincrd, double *maxcrd, double *xcrd, double *ycrd, double *zcrd)
{
	int flag[DIMENSION];
//	int sum = 0;
	int k;
	bool erpt_flag = true;
	int bt[6];

    bt[0]=determine_face_type(xcrd[0],maxcrd[0],mincrd[0]);
    bt[1]=determine_face_type(xcrd[1],maxcrd[0],mincrd[0]);
    bt[2]=determine_face_type(ycrd[0],maxcrd[1],mincrd[1]);
    bt[3]=determine_face_type(ycrd[1],maxcrd[1],mincrd[1]);
    bt[4]=determine_face_type(zcrd[0],maxcrd[2],mincrd[2]);
    bt[5]=determine_face_type(zcrd[1],maxcrd[2],mincrd[2]);

    for (k=0; k<DIMENSION; k++)
    	flag[k]=abs(bt[2*k] + bt[2*k+1]);

    for (k=0; k<DIMENSION; k++)
    	if (flag[k] == 2)
    		erpt_flag = false;

    return erpt_flag;
}

/*
 * The strategy for setting up eruption boundary condition and adding new erupted particles is
 * 1) adding guest particle in guest buckets ---> well, actually, not really sure whether this is a good idea or not,
 *  The good thing is synchronization is not needed, In pressure ghost adding and wall ghost adding, this strategy is adopted
 *  Another good thing is that this might help to improve efficiency --> does not require additional communication
 *  The bad thing is that this will make the code difficult to maintain, as in one function, I am trying to do too many things.
 *
 *  2) The second strategy is about how to make sure that newly added particles should belong to "this" process
 *  simple: first put all particle in a temp hashtable, then go through all bucket on "this" process and determine whether particles belong to "this" process
 */
int
setup_erupt(int myid, THashTable * P_table, HashTable * BG_mesh,
        TimeProps * timeprops, MatProps* matprops, SimProps *simprops, int numprocs)
{
	double range_x[2];
	double range_y[2];
	double range_z[2];
	double crd_p[DIMENSION];
	unsigned key[TKEYLENGTH];

	unsigned tempkey[TKEYLENGTH];
	bool erpt;
	//unsigned num_particle = 0;
    int i, j, k;
    int tempid;
    unsigned tkeylen = TKEYLENGTH;
   
    double dist;
    double vel;
    double rvsq=rv_P*rv_P;
    unsigned add_step = 0;

    double mindom[DIMENSION], maxdom[DIMENSION];

	// get min-max domain from hashtable, for key generation
	for (i = 0; i < DIMENSION; i++)
	{
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	}

	THashTable * P_temp = new THashTable(ERUPT_TABLE_SIZE, 2017, mindom, maxdom);; //particle hash table
	//We do not need range_z actually, for simple model with flat ground,
	//the max for z is zero and the min for z is determined by number of ghost particle layers.
    //double range_z[2];

	//determine the rough range of eruption duck
	range_x[0] = -rv_P;
	range_x[1] = rv_P;
	range_y[0] = -rv_P;
	range_y[1] = rv_P;
    range_z[1] = 0.; // not exact, should use ground height
    range_z[0] = range_z[1]-(matprops->smoothing_length)*1.5*PARTICLE_DENSITY;


	//determine the mass and smooth length of the particle
	double mss, sml;
    double des = rhov_P;
	int np = num_erupt_particles; // total number of erupt ghost particle
	Compute_mass (np, range_x, range_y, range_z, des, &mss, &sml);

	simprops->mass_of_phase2 = mss;
	simprops->sml_of_phase2 = sml;

    double sml2 = 0.5*sml;
    double bot = range_z[0];
    simprops->update_bot(bot);

	//create new particles and put them into a temporary hash table: coordinate -> key -> new particle
    double normc[DIMENSION];
    crd_p[2]=range_z[1]-sml2;
    while ( crd_p[2] >= range_z[0] )
    {
    	crd_p[0]=Pos_v_P[0]+sml2;
    	while ( crd_p[0] <= range_x[1] )
    	{
    		//y+ direction
    		crd_p[1]=Pos_v_P[1]+sml2;
    		while ( crd_p[1] <= range_y[1])
    		{
    			dist = 0;
    			for (j=0; j<2; j++)
    				dist += (crd_p[j]*crd_p[j]);
    			if (dist <= rvsq) // if particle is in the duct
    			{
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
    			    vel = vel_prof (rv_P, dist, Vv0_P);
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    			}//end of if
    			crd_p[1] += sml;
    		}//end of while y+

    		//y- direction
    		crd_p[1]=Pos_v_P[1]-sml2;
    		while ( crd_p[1] >= range_y[0])
    		{
    			dist = 0;
    			for (j=0; j<2; j++)
    				dist += (crd_p[j]*crd_p[j]);
    			if (dist <= rvsq) // if particle is in the duct
    			{
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
    			    vel = vel_prof (rv_P, dist, Vv0_P);
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    			}//end of if
    			crd_p[1] -= sml;
    		}//end of while y-

    		crd_p[0] += sml;
    	}//end of while x+

    	//x- direction
    	crd_p[0]=Pos_v_P[0]-sml2;
    	while ( crd_p[0] >= range_x[0] )
    	{
    		//y+ direction
    		crd_p[1]=Pos_v_P[1]+sml2;
    		while ( crd_p[1] <= range_y[1])
    		{
    			dist = 0;
    			for (j=0; j<2; j++)
    				dist += (crd_p[j]*crd_p[j]);
    			if (dist <= rvsq) // if particle is in the duct
    			{
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
    			    vel = vel_prof (rv_P, dist, Vv0_P);
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    			}//end of if
    			crd_p[1] += sml;
    		}//end of while y+

    		//y- direction
    		crd_p[1]=Pos_v_P[1]-sml2;
    		while ( crd_p[1] >= range_y[0])
    		{
    			dist = 0;
    			for (j=0; j<2; j++)
    				dist += (crd_p[j]*crd_p[j]);
    			if (dist <= rvsq) // if particle is in the duct
    			{
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
    			    vel = vel_prof (rv_P, dist, Vv0_P);
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    			}//end of if
    			crd_p[1] -= sml;
    		}//end of while y-

    		crd_p[0] -= sml;
    	}//end of while x-


    	crd_p[2] -= sml;
    }//end of while z

    //mark bucket as erupt_source
    HTIterator * itr = new HTIterator (BG_mesh);
    Bucket * Curr_buck;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

    double coordtmp[DIMENSION];
    double mincrd[DIMENSION], maxcrd[DIMENSION];
    THTIterator * itr2 = new THTIterator (P_temp);
    Particle * Curr_part;
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

    	    //determine whether some portion of the bucket include
    	    erpt = determine_erupt_buket (mincrd, maxcrd, range_x, range_y, range_z);

    	    if (erpt)
    	    	 Curr_buck->mark_erupt();
    	}//end of if bucket is not brief
    }// end of loop go through all buckets

    /*go through all erupt buckets and
     * 1) delete all non-erupt ghost particles that within the range of the duct
     * 2) add newly generated erupt particles into hash table of particles.
     * 3) put newly generated particles into bucket
    */
    vector < TKey > pnew;//particle key
    vector < TKey > plist;
    vector < TKey >::iterator p_itr;
    Particle *pj;
    itr->reset();
    int bctp ;
    //go through all buckets
    while ((tempptr=itr->next ()))
    {
    	breif_buck = (BriefBucket *) tempptr;
    	if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
    		continue;
    	else
    	{
    		Curr_buck = (Bucket*) tempptr;
    	      if (Curr_buck->is_erupt ())
    	      {
    	    	  //set particles type to be 0
    	    	  Curr_buck->put_particles_type (0);

    	    	  //check all particle in the erupt bucket and remove them when necessary!
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

    	    		  if ((dist <= rvsq) && crd_p[2]<=range_z[1] && crd_p[2]>=range_z[0] &&  (!pj->is_erupt_ghost())) // if particle is in the duct
    	    		  {
    	    			  for (k = 0; k<TKEYLENGTH; k++)
    	    			      tempkey[k] = pj->getKey().key[k];

    	    			  /*
    	    			   * Here what I did is only remove them from P_table and bucket particle list!
    	    			   * But the particle in other particles neighbour list is not deleted!---> will be done in search neighbor (please double check)
    	    			   * And the particle as guest on other processes is not removed!--->well, all particle on other processes will be removed as other processes will execute the same command. ----> you should have already been aware of the fact that no "non_guest" constrain is imposed here.
    	    			   * And the particle, if it has image, the image should also be removed! ---> will be done in search_image and imposing BC (please double check----> yes, I checked that, no problem)
    	    			   */
    	    			  P_table->remove(tempkey); //remove from the hashtable

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
    						      cout << "bctp incorrect in function set_up_erupt!\n" << endl;
    					   }
    	    		  }
    	    	  }

    	    	  //update particle list
    	    	  Curr_buck->put_new_plist (pnew);
    	    	  Curr_buck->update_particles();

    	    	  //go through all temporarily added particles and added it into

    	    	  /*
    	    	  * Here what I did is only add them to P_table and bucket particle list!
    	    	  * But the particle in other particles neighbour list is not deleted!
    	    	  * And the particle as guest on other processes is not removed!
    	    	  */
    	  	      for (i = 0; i < DIMENSION; i++)
    	  	      {
    	  	          mincrd[i] = *(Curr_buck->get_mincrd () + i);
    	  	          maxcrd[i] = *(Curr_buck->get_maxcrd () + i);
    	  	      }
    	    	  itr2->reset();
    	    	  while ((Curr_part = (Particle *) itr2->next() ))
    	    	  {
    	    	 	  for (i=0; i < DIMENSION; i++)
    	    	 	      coordtmp[i]= *(Curr_part->get_coords ()+i);

    	    	 	  if (in_bucket(maxcrd, mincrd, coordtmp))
    	    	 	  {
    	    	           for (k=0; k<TKEYLENGTH; k++)
    	    	                key[k]=Curr_part->getKey().key[k];

    	    	           //
    	                   //default involved is zero---> need to double check what involved should be allocated here-->No problem, should be not involved.
    	                   if (Curr_buck->is_guest())
    	                   {
    	                	      Curr_part->put_guest_flag(true);
    	                          tempid = Curr_buck->get_myprocess ();
    	                          Curr_part->put_my_processor(tempid);
    	                   }

    	    	 	       TKey tmpkey(key);
    	    	 	       Curr_buck->add_erupt_ghost_particle(tmpkey);
    	    	    	   P_table->add(key, Curr_part);
    	    	 	   }
    	    	  }

    	      }// end of loop go through all buckets
    	}
    }//end of go through all buckets


    //clear up:
    delete itr;
    delete itr2;
    delete P_temp;

return 0;
}

//function for adding new ghost erupt particles at the bottom of the duck
void
add_new_erupt(int myid, THashTable * P_table, HashTable * BG_mesh,
        TimeProps * timeprops, MatProps* matprops, SimProps *simprops, double dt)
{
    double t_add, t_each;
    int n;
    int tempid;
    double bot;
    double crd_p[DIMENSION];
    double range_x[2];
    double range_y[2];
	double range_z[2];
    double normc[3];
    double dist;
    double vel;
    double rvsq=rv_P*rv_P;
	unsigned key[TKEYLENGTH]; //should be time depend key for particle
    int i, j, k;
//    int ii;
    int num_particle=0;
    unsigned tkeylen = TKEYLENGTH;
    double mindom[DIMENSION], maxdom[DIMENSION];

	// get min-max domain from hashtable, for key generation
	for (i = 0; i < DIMENSION; i++)
	{
	    mindom[i] = *(P_table->get_minDom() + i);
	    maxdom[i] = *(P_table->get_maxDom() + i);
	}

	THashTable * P_temp = new THashTable(ERUPT_TABLE_SIZE, 2017, mindom, maxdom);; //particle hash table

    bot = simprops->get_bot();
    double sml = simprops->sml_of_phase2;
    double sml2 =0.5*sml;
    double t_total=timeprops->time;
    double mss = simprops->mass_of_phase2;
	unsigned add_step;
    //the more robust way: use time step as add_step
	add_step = timeprops->step;

    //determine the rough range of eruption duck
    range_x[0] = -rv_P;
    range_x[1] = rv_P;
    range_y[0] = -rv_P;
    range_y[1] = rv_P;
    range_z[1] = 0.; // not exact, should use ground height
    range_z[0] = range_z[1]-(matprops->smoothing_length)*1.5*PARTICLE_DENSITY;

    t_total -= dt; //get t_total of the previous time step

	//create new particles: coordinate -> key -> new particle
    //x + direction
    crd_p[0]=Pos_v_P[0]+sml2;
    while ( crd_p[0] <= range_x[1] )
    {
    	//y+ direction
    	crd_p[1]=Pos_v_P[1]+sml2;
    	while ( crd_p[1] <= range_y[1])
    	{
    		dist = 0;
    		for (j=0; j<2; j++)
    		    dist += (crd_p[j]*crd_p[j]);
    		if (dist <= rvsq) // if particle is in the duct
    		{
    		   vel = vel_prof (rv_P, dist, Vv0_P); //Velocity depends on dist
    		   t_each = sml/vel;
    		   t_add = fmod (t_total, t_each) + dt; // fmod (t_total, t_each) is the balance from previous time_step, the balance should be smaller than t_each
    		   n = floor(t_add/t_each);
    		   t_add = t_add - n*t_each;  //t_add here is balance that will have after adding new particles
    		   for (i=0; i<n; i++)
    		   {
    			    crd_p[2] = bot + (t_add + i * t_each) * vel;
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
		    		if (P_table->lookup(key))
		    		{
		    			fprintf(stderr, "ERROR: Trying to add particle "
		    				    	 	"twice on same location.\n");
		    			exit(1);
		    		}
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    		    }//end of for
    		}
    		crd_p[1] += sml;
    	}//end of while y+

    	//y- direction
    	crd_p[1]=Pos_v_P[1]-sml2;
    	while ( crd_p[1] >= range_y[0])
    	{
    		dist = 0;
    		for (j=0; j<2; j++)
    		    dist += (crd_p[j]*crd_p[j]);
    		if (dist <= rvsq) // if particle is in the duct
    		{
    		   vel = vel_prof (rv_P, dist, Vv0_P);
    		   t_each = sml/vel;
    		   t_add = fmod (t_total, t_each) + dt;
    		   n = floor(t_add/t_each);
    		   t_add = t_add - n*t_each;
    		   for (i=0; i<n; i++)
    		   {
    			    crd_p[2] = bot + (t_add + i * t_each) * vel;
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
		    		if (P_table->lookup(key))
		    		{
		    			fprintf(stderr, "ERROR: Trying to add particle "
		    				    	 	"twice on same location.\n");
		    			exit(1);
		    		}
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    		    }//end of for
    		}
    		crd_p[1] -= sml;
    	}//end of while y-

    	crd_p[0] += sml;
    }//end of while x+

    //x - direction
    crd_p[0]=Pos_v_P[0]-sml2;
    while ( crd_p[0] >= range_x[0] )
    {
    	//y+ direction
    	crd_p[1]=Pos_v_P[1]+sml2;
    	while ( crd_p[1] <= range_y[1])
    	{
    		dist = 0;
    		for (j=0; j<2; j++)
    		    dist += (crd_p[j]*crd_p[j]);
    		if (dist <= rvsq) // if particle is in the duct
    		{
    		   vel = vel_prof (rv_P, dist, Vv0_P);
    		   t_each = sml/vel;
    		   t_add = fmod (t_total, t_each) + dt;
    		   n = floor(t_add/t_each);
    		   t_add = t_add - n*t_each;
    		   for (i=0; i<n; i++)
    		   {
    			    crd_p[2] = bot + (t_add + i * t_each) * vel;
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
		    		if (P_table->lookup(key))
		    		{
		    			fprintf(stderr, "ERROR: Trying to add particle "
		    				    	 	"twice on same location.\n");
		    			exit(1);
		    		}
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    		    }//end of for
    		}
    		crd_p[1] += sml;
    	}//end of while y+

    	//y- direction
    	crd_p[1]=Pos_v_P[1]-sml2;
    	while ( crd_p[1] >= range_y[0])
    	{
    		dist = 0;
    		for (j=0; j<2; j++)
    		    dist += (crd_p[j]*crd_p[j]);
    		if (dist <= rvsq) // if particle is in the duct
    		{
    		   vel = vel_prof (rv_P, dist, Vv0_P);
    		   t_each = sml/vel;
    		   t_add = fmod (t_total, t_each) + dt;
    		   n = floor(t_add/t_each);
    		   t_add = t_add - n*t_each;
    		   for (i=0; i<n; i++)
    		   {
    			    crd_p[2] = bot + (t_add + i * t_each) * vel;
    			    for (int ii = 0; ii < DIMENSION; ii++)
    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

    			    THSFC3d (normc, add_step, &tkeylen, key);
		    		if (P_table->lookup(key))
		    		{
		    			fprintf(stderr, "ERROR: Trying to add particle "
		    				    	 	"twice on same location.\n");
		    			exit(1);
		    		}
    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
    			    // add to hash-table
    			    P_temp->add(key, pnew);
    		    }//end of for
    		}
    		crd_p[1] -= sml;
    	}//end of while y-

    	crd_p[0] -= sml;
    }//end of while x-

	// put newly generated particles into the bucket, particles has already been added into P_table
	    HTIterator * itr = new HTIterator (BG_mesh);
	    Bucket * Curr_buck;
		BriefBucket *breif_buck = NULL;
		void * tempptr =NULL;

	    double coordtmp[DIMENSION];
	    double mincrd[DIMENSION], maxcrd[DIMENSION];
	    THTIterator * itr2 = new THTIterator (P_temp);
	    Particle * Curr_part;

		while ((tempptr=itr->next ()))
		{
			breif_buck = (BriefBucket *) tempptr;
			if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
				continue;
			else
			{
				Curr_buck = (Bucket*) tempptr;
			      if (Curr_buck->is_erupt())
			      {
				    for (i = 0; i < DIMENSION; i++)
				    {
				        mincrd[i] = *(Curr_buck->get_mincrd () + i);
				        maxcrd[i] = *(Curr_buck->get_maxcrd () + i);
				     }

				     itr2->reset();
				     while ((Curr_part = (Particle *) itr2->next() ))
				     {
				    	for (i=0; i < DIMENSION; i++)
				    		coordtmp[i]= *(Curr_part->get_coords ()+i);
				    	if (in_bucket(maxcrd, mincrd, coordtmp))
				    	{
			                for (k=0; k<TKEYLENGTH; k++)
			                    key[k]=Curr_part->getKey().key[k];
			    	           //
			                //default involved is zero---> need to double check what involved should be allocated here
			                if (Curr_buck->is_guest())
			                {
			                  Curr_part->put_guest_flag(true);
			                  tempid = Curr_buck->get_myprocess ();
			                  Curr_part->put_my_processor(tempid);
			                }
				    	   TKey tmpkey(key);
				    	   Curr_buck->add_erupt_ghost_particle(tmpkey);
		                   P_table->add(key, Curr_part);
		                   num_particle++; //will be used to update THASHTAB
				        }
				     }
			      }// end of if
			}//end of if bucket is not brief
		}//end of while loop go through all buckets

	//clean up
    delete itr;
    delete itr2;
    delete P_temp;

    return;
}
