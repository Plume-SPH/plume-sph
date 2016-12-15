/*
 * update_pos.cc
 *
 *  Created on: Feb 25, 2015
 *      Author: zhixuanc
 */

#include <vector>
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <particler.h>
#include <properties.h>

#include "constant.h"
#include "sph_header.h"
#include "parameters.h"

#ifdef DEBUG
#  include <debug_header.h>
#  include <iostream>
#endif

void
update_pos(int myid, THashTable * P_table, HashTable * BG_mesh,
            TimeProps * timeprops, MatProps * matprops, int *lost)
{
  int dir[DIMENSION];
  double pos[DIMENSION], coord[DIMENSION], vel[DIMENSION];
  double mincrd[DIMENSION], maxcrd[DIMENSION];

  double sml_of_phase2 = matprops->smoothing_length;

  int i;

  double dt = timeprops->dtime;
  bool los;

#ifdef DEBUG
   bool do_search = false;
   bool check_contain = true;
   unsigned keycheck[TKEYLENGTH] = {69562537, 292385725, 0};
   unsigned keytemp[TKEYLENGTH] ;
#endif

  // update particle positions
  THTIterator *itr = new THTIterator(P_table);
  Particle *p;

  while ((p = (Particle *) itr->next()))
  {
    /*guest particles are moved as well as some of them
	* The reason why I also move guest particles is guest particles might move in to native doamin.
	* And it is necessary to monitor it.
    * will move to different partitions */
    if (p->is_real ())
    {

#ifdef DEBUG
		if (do_search)
		{
		    for (i = 0; i < TKEYLENGTH; i++)
			    keytemp[i] = p->getKey ().key[i];

		    if (find_particle (keytemp, keycheck))
			    cout << "The particle found!" << endl;
		}
#endif

      // velocity and coordinates
      for (i = 0; i < DIMENSION; i++)
      {
#ifdef HAVE_TURBULENCE_LANS
         vel[i] = *(p->get_smoothed_velocity() + i);
#else
         vel[i] = *(p->get_vel() + i);
#endif
         coord[i] = *(p->get_coords() + i);
      }

      /* update particle positions
       * The way I update position of erupted particles is not the perfect way:
       * all particles move with the average upward velocity
       * but the velocity of each particle as a primitive variable is obtained from a parabolic profile
       */
//      if (p->is_real ())
         for (i = 0; i < DIMENSION; i++)
             pos[i] = coord[i] + dt*vel[i];

       p->put_coords(pos);

    }//end of if particle is real

    //Because updating of erupted particle is based on unsmoothed velocity, ---> Why?
    //SO it is better to do real and erupted separately
    else if (p->is_erupt_ghost()  && timeprops->iferupt() )
    {

#ifdef DEBUG
		if (do_search)
		{
		    for (i = 0; i < TKEYLENGTH; i++)
			    keytemp[i] = p->getKey ().key[i];

		    if (find_particle (keytemp, keycheck))
			    cout << "The particle found!" << endl;
		}
#endif

      // velocity and coordinates
      for (i = 0; i < DIMENSION; i++)
      {
         vel[i] = *(p->get_vel() + i);
         coord[i] = *(p->get_coords() + i);
      }

      /* update particle positions
       * The way I update position of erupted particles is not the perfect way:
       * all particles move with the average upward velocity
       * but the velocity of each particle as a primitive variable is obtained from a parabolic profile
       */
//      if (p->is_real ())
         for (i = 0; i < DIMENSION; i++)
             pos[i] = coord[i] + dt*vel[i];

      p->put_coords(pos);

#ifndef SIMULATE_ASH
      //Change ghost to real; change erupt to false
      if ((pos[2]>=Lz_P[0]) ) //need to make it more general!
      {
    	  p->erupt_turn_real();
    	  p->put_smlen(sml_of_phase2);
    	  p->set_involved_flag(INVOLVED); //equivalent to set involved to be true
      }
#else
      //Change ghost to real; change erupt to false
      double routsq=r_out_P*r_out_P;
      double dist;
      dist = 0.0;
      for (i=0; i<2; i++)
    	  dist += pos[i]*pos[i];
      if ((dist>=routsq)) //need to make it more general!
      {
    	  p->erupt_turn_real();
    	  p->put_smlen(sml_of_phase2);
    	  p->set_involved_flag(INVOLVED); //equivalent to set involved to be true
      }
#endif
    } //end of if particle is erupted

  }//end of go through all particles

  // move-in and move out particles from buckets
  vector < TKey > plist;
  vector < TKey >::iterator p_itr;
  vector < TKey > my_realp, my_particles;

  HTIterator *it2 = new HTIterator(BG_mesh);
  Bucket *curr_bucket;
  BriefBucket *breif_buck = NULL;
  void * tempptr =NULL;

  while ((tempptr= it2->next()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  curr_bucket = (Bucket*) tempptr;
		  if (curr_bucket->is_active())
		  {
		      my_realp.clear();
		      my_particles.clear();

		      // get bucket limits
		      for (i = 0; i < DIMENSION; i++)
		      {
		        mincrd[i] = *(curr_bucket->get_mincrd() + i);
		        maxcrd[i] = *(curr_bucket->get_maxcrd() + i);
		      }

		      // get list of particles in the bucket
		      plist = curr_bucket->get_plist();
		      for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
		      {
		        Particle *p_curr = (Particle *) P_table->lookup(*p_itr);

		        assert(p_curr);
		        if (p_curr->is_real ()  || p_curr->is_erupt_ghost())//if particle is real (no matter it is guest or not.)
		        {
		    		for (i = 0; i < DIMENSION; i++)
		    		     pos[i] = *(p_curr->get_coords() + i);

		    		/*
		    		 * need to re-think about this part. --->MIXED bucket should be only on ground.
		    		 * For the bucket over ground, the boundary is not fixed and not necessary to be MIXED
		    		 */
		    		los = false;
		    		//only real particle will lost, erupt ghost will not lost.
		    		if (p_curr->is_real ())
		        	{
		        		/* there is potential problem that the particles cross the boundary and goes to other other buckets which is next to the MIXED buckets.
		        		 * ---> Hopefully this will not happen!
		        		 */
		            //if (curr_bucket->get_bucket_type() == MIXED && (curr_bucket->get_bucket_index())[4] == -1) // only happens for underground MIXED bucket---> The old code, should also happens at the top and side of the domain.
		    			if (curr_bucket->get_bucket_type() == MIXED)
		    			   los = curr_bucket->determine_escape(pos);
		             }// end of if particles is real ---> to delete particles which crossed the boundary!

		    		// if particle crosses boundary -- remove it
		    		if (los)
		    		{
		    		    P_table->remove(*p_itr);
		    		    delete p_curr;

		    		    lost++;
		    		    continue;
		    		}

		          // check where particle is going,
		          // ... if going anywhere at all
		          bool left_curr_bucket = false;

		          for (i = 0; i < DIMENSION; i++)
		          {
		            if (pos[i] < mincrd[i])
		            {
		              dir[i] = 1;
		              left_curr_bucket = true;
		            }
		            else if (pos[i] >= maxcrd[i])
		            {
		              dir[i] = 2;
		              left_curr_bucket = true;
		            }
		            else
		              dir[i] = 0;
		          }

		          // determine, to which bucket particle migrating
		          // there are 3 special cases, that are worth considering, i.e.
		          //    1) a native particle enters a guest bucket
		          //    2) a guest particle  enters native bucket
		          //    3) leaves a guest bucket and leaves domain

		          if (left_curr_bucket)
		          {
		            if (!p_curr->is_guest())    // should always have neighbors
		            {
		              // 1) a native particle enters a guest bucket
		              Key neigh_key = curr_bucket->which_neigh(dir); //bucket key
		              Bucket *neigh = (Bucket *) BG_mesh->lookup(neigh_key);

#ifdef DEBUG
		              if (check_contain)
		                  if (!neigh->contains(pos))
		                      cout << "no contain happens" << endl;
#endif

		              assert(neigh->contains(pos));

		              // if real particle moves into guest bucket
		              // delete it
		              if (neigh->is_guest())
		              {
		                P_table->remove(*p_itr);
		                delete p_curr;

		                continue;
		              }
		              else  // if neighbor bucket is a native ...
		              {
		                //if real particle moves into a pressure bucket,
		                //remove the particle and los++;
		                if (!neigh->get_has_involved())//if has_involved = 0
		                {
		                	  P_table->remove(*p_itr);
		                	  delete p_curr;

		                	  lost++;// real particle lost on local proc
		                	  continue;
		                }
		                else //if the neigh is not pressure bucket
		                {
		                   // if an empty bucket got particle ...
		                   //   1 -> turn on have_realp flag
		                   //   2 -> turn on adapt flag
		                   if ((! neigh->has_real_particles()) && (p_curr->is_real()))
		                   {
		                     neigh->set_real_particles(true);
		                   }

		                   if ((! neigh->has_erupt_ghost_particles()) && (p_curr->is_erupt_ghost ()))
		                   {
		                        neigh->set_erupt_ghost_particles(true);
		                   }
		                   neigh->add_particle(*p_itr); // 1) If the particle list in neigh has not been updated, after while, when update the particle list in that bucket, the particle that was added here will just been kept.
		                                                // 2) If the particle list in neigh has already been updated, it does not hurt to add a new particle into it.
		                 }
		              }
		            }//end of if  (!p_curr->is_guest())

		            else if (curr_bucket->which_neigh_proc(dir) == myid)
		            {
		              //case 2 : guest particle enter native bucket
		              Key neigh_key = curr_bucket->which_neigh(dir);
		              Bucket *neigh = (Bucket *) BG_mesh->lookup(neigh_key);

		#ifdef DEBUG
		              if (check_contain)
		                  if (!neigh->contains(pos))
		                      cout << "no contain happens" << endl;
		#endif
		              assert(neigh->contains(pos));

		              //if real particle moves into a pressure guest bucket,
		              //remove the particle and los++;
		              if (!neigh->get_has_involved())//if has_involved = 0
		              {
		            	  P_table->remove(*p_itr);
		            	  delete p_curr;

		//            	  lost++;// guest particle lost on local proc
		            	  continue;
		              }

		              // the particle is coming from partition not
		              // on current proc
		              p_curr->put_guest_flag(false); //turn the particle to be no-guest particle --->but need update neigh info
		              neigh->add_particle(*p_itr);
		              if ((! neigh->has_real_particles()) && (p_curr->is_real()))
		              {
		                   neigh->set_real_particles(true);
		              }

		              if ((! neigh->has_erupt_ghost_particles()) && (p_curr->is_erupt_ghost ()))
		              {
		                   neigh->set_erupt_ghost_particles(true);
		              }
		            }
		            else                // a guest particles leaves curent domain
		            {
		              // case 3: guest particle leave the current domain
		              P_table->remove(*p_itr);
		              delete p_curr;

		              continue;
		            }
		          }// end of if left_curr_bucket is true.
		          else
		            my_realp.push_back(*p_itr);

		        }//end of if (if particle is either real or erupted ghost)

		        else //keep these ghost particles in the bucket as they did not move
		          my_particles.push_back(*p_itr);
		      }//end of go through all particles in the bucket


		      if (!my_realp.empty())
		        my_particles.insert(my_particles.end(), my_realp.begin(),
		                            my_realp.end());
		      else
		      {
		        curr_bucket->set_real_particles(false);
		        curr_bucket->set_erupt_ghost_particles(false);
		      }

		      // update new list
		      curr_bucket->put_new_plist(my_particles);
		  }//end of if ...
	  }//end of if bucket is not brief bucket
  }//end of while loop go through all buckets

  // Update particle lists in the buckets--->go through all buckets again!
  it2->reset();

  while ((tempptr= it2->next()))
  {
	  breif_buck = (BriefBucket *) tempptr;
	  if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
		  continue;
	  else
	  {
		  curr_bucket = (Bucket*) tempptr;
		  if ((curr_bucket->is_active()) &&
		      //  (/*curr_bucket->get_bucket_type() != UNDERGROUND &&*/ !curr_bucket->is_erupt ()) &&
		        (curr_bucket->get_bucket_type() != PRESS_BC))
		  {
		      curr_bucket->update_particles();
		      int numofp = curr_bucket->get_plist().size();

		      if (numofp > MAX_PARTICLES_PER_BUCKET)
		      {
		        int nfl = 0, ngh = 0;

		        vector < TKey > pl = curr_bucket->get_plist();
		        for (i = 0; i < pl.size(); i++)
		        {
		          Particle *p = (Particle *) P_table->lookup(pl[i]);

		          fprintf(stderr, "%f, %f, %f, %d\n", *(p->get_coords()),
		                  *(p->get_coords() + 1), *(p->get_coords() + 2), p->need_neigh());
		          if (p->is_real())
		            nfl++;
		          else
		            ngh++;
		        }
		        Key ck = curr_bucket->getKey();

		        fprintf(stderr, "Problem in (%u, %u, %u) bucket.\n \
		                        No. of particles exceeded threshold: %d\n", ck.key[0], ck.key[1], ck.key[2], numofp);
		        fprintf(stderr, "real particles: %d, ghost particles: %d\n", nfl, ngh);
		        fprintf(stderr, " Min coords: %e, %e, %e\n", *(curr_bucket->get_mincrd()),
		                *(curr_bucket->get_mincrd() + 1),*(curr_bucket->get_mincrd() + 2));
		        fprintf(stderr, " Max coords: %e, %e, %e\n", *(curr_bucket->get_maxcrd()),
		                *(curr_bucket->get_maxcrd() + 1), *(curr_bucket->get_maxcrd() + 2));
		        fprintf(stderr, " ... at %s: %d\n\n", __FILE__, __LINE__);
		        return;
		      }//end of if (if number of particles in the bucket excess the maximum number that allowed
		  }//end if bucket is not pressure_bc, for the underground bucket, I should consider them here because the erupted particles is from underground.
	  }//end of if curr_bucket is not brief bucket
  }//end of while loop go through all buckets
  // clean up
  delete itr;
  delete it2;

//  return adapt;
  return;
}
