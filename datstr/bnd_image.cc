/*
 * bnd_image.cc
 *
 *  Created on: Mar 13, 2015
 *      Author: zhixuanc
 */


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <vector>
#include <list>
#include <iostream>
#include <cassert>

using namespace std;

#include <constant.h>
#include <hashtab.h>
#include <thashtab.h>

#include "bucket.h"

#include <particle.h>
#include <bnd_image.h>
#include <exvar.h>


void
search_bnd_images (int myid, THashTable * P_table, HashTable * BG_mesh,
                   list < BndImage > & Image_table, int reset)
{
  bool found_bucket;
  int i, j;
  int dir[DIMENSION];
  double coord[DIMENSION], tmpc[DIMENSION];
  double refc[DIMENSION], intsct[DIMENSION];
  unsigned ghst_key[TKEYLENGTH], buck_key[KEYLENGTH];

  // get a domain length
  /*
   * There is some problem with the following way to determine the domain length, as now the domain is adaptive
   * There should be another way to determine the domain_length --->go through all bucket -->found the most left non-empty bucket and most right bucket-->determine the domain length in that direction
   * However, as the shape of the domain is not regular cube, the most reliable way is go through all bucket and then determine whether the image is belong to certain bucket or not.
   *
   * Another problem that the following method has is: the domain size is different in different direction
   *
   * **Based on assumption that the domain is flat, we only need to check the z direction.
   * */


  //This is the old code:
//  double len_domain = *(BG_mesh->get_maxDom()) - *(BG_mesh->get_minDom());

  //The new code only check the whether image go out of the domain from the top.
  //even if image is in a bucket which is not actved yet, it does not hurt-->that wall ghost corresponding to that image should be far away from the wall boundary ---> so it does not matter.
  double len_domain = *(BG_mesh->get_maxDom()+2) - *(BG_mesh->get_minDom()+2);

  // always reset after repatition. If there is no repartitioning, there is no necessary to reset Image_table
  /*
   * Both wall ghost particles and their images keep unchanged.
   * There are two cases that need pay attention to:
   * 1) new wall ghost particles is added while domain expanding--> In this case, the old wall ghosts and their images keep unchanged
   * what need to do is search for wall ghost for the newly added wall ghost, of course, the image of newly generated wall ghost need to be found and added in the Image_table
   * --->But the old Image_table is not necessary to be changed.---> the has_reflection will indicate whether the reflection is been searched or not.
   * 2) repartitioning  --> In this case, the image_table need to be reset-->The reason is:
   * after repartition, some wall ghost originally belong to "this process" might go to other processes and some other particles originally from other processes might now belong this process
   * ---> As a result, the old table might have some images that belongs to wall ghost particles that in other processes
   * */
  if (reset)
  {
    Image_table.clear ();
    THTIterator *p_itr = new THTIterator (P_table);
    Particle *p_temp = NULL;

    while ((p_temp = (Particle *) p_itr->next ()))
    	if (p_temp->is_wall_ghost())
            p_temp->set_reflection(false);//set reflection to false.
    delete p_itr;
  }

  Bucket *buck = NULL, *buck1 = NULL, *buck2 = NULL;
  Bucket *buck_neigh;
  HTIterator *itr = new HTIterator (BG_mesh);

  while ((buck = (Bucket *) itr->next ()))//go through all buckets
    if ((buck->has_wall_ghost_particles ()) &&  //bucket has wall ghost particle and is not guest
        (!buck->is_guest ()))
    {
      vector < TKey > plist = buck->get_plist ();
      vector < TKey >::iterator p_itr;
      const int * neigh_proc = buck->get_neigh_proc ();
      const Key * neighbors = buck->get_neighbors ();

      if (plist.size () > 0)
        for (p_itr = plist.begin (); p_itr != plist.end (); p_itr++)
        {
          Particle *p_ghost = (Particle *) P_table->lookup (*p_itr);
          if (p_ghost->is_wall_ghost())//wall_bc_ghost
        	  if (! p_ghost->has_reflection ())//do not have reflection yet ---> because of this, it is not necessary to reset Image_table if new wall ghost is added.
        	  {
                  {
                    // coordinate of ghost particle
                    for (i = 0; i < DIMENSION; i++)
                      coord[i] = *(p_ghost->get_coords () + i);

                    // initalize distance and state-variables
                    double bnddist = 1.0E+10;     // very large number
                    int img_proc = myid;

                    // first look in current bucket,
                    // ... if current bucket has ground boundary information
                    if ((buck->get_bucket_type () == MIXED) && (buck->is_onground ()))
                    {
                      buck2 = buck;
                      bnddist = buck->get_bnddist (coord, intsct);
                    }

                    // then check neighboring buckets, and compare
                    // ... perpendicular distance
                    for (i = 0; i < NEIGH_SIZE; i++)
                      if (neigh_proc[i] > -1)
                      {
                        buck_neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);

                        // skip neighbor if it belongs to foreign process
                        // but is not synchronized----->Yes, this is a problem.
                        if (!(buck_neigh) && (neigh_proc[i] != myid))
                          continue;

                        // involve only if neighbor bucket is a on ground boundary bucket
                        if ((buck_neigh->get_bucket_type () == MIXED) && (buck_neigh->is_onground ()))
                        {
                          double tmpdist = buck_neigh->get_bnddist (coord, tmpc);
                          if (tmpdist < bnddist)//What needed is the shortest distance.
                          {
                            bnddist = tmpdist;
                            buck2 = buck_neigh;
                            img_proc = neigh_proc[i]; //I have question here: how could you make sure that the image is in the neighbor bucket? It is also possible that the image is not in the MIXED bucket.
                            for (j = 0; j < DIMENSION; j++)
                              intsct[j] = tmpc[j];
                          }
                        }
                      }//finish check all neighbor buckets, if the process id of neighbor bucket is valid

                    // if bnd-distace is too large, then we were unable to
                    // find the reflection, hopefully we'll be fine
                    if (bnddist > len_domain) //check whether the particle is inside the domain or not.
                    {
                      fprintf(stderr, "image of wall_ghost particle with coordinate: %f, %f, %f is out of the domain! \n", coord[0], coord[1], coord[2]);
                      continue;
                    }

                    /*
                     *  Now, we have found the samllest perpendicular distance
                     *  and bucket which contains the nearest boundary,
                     *  we can find the mirror image, i.e. a point
                     *  along normal at 2*bnddist from current particle coords
                     */
                    for (j = 0; j < DIMENSION; j++)
                      refc[j] = 2 * intsct[j] - coord[j];

                    // find the bucket containing reflection
                    /*
                     * The problem here is that:
                     * if the image is not contained by any bucket on current processor, what should I do?
                     *
                     * There are three possible situations:
                     * 1) image is in a local bucket
                     * 2) image is in a guest bucket
                     * 3) image is not in neith a local bucket nor a guest bucket on current processors.
                     * The third situation will cause some trouble for me, the only thing that I can is to pray that this will happen.
                     * If this happens, what should I do? ---> there are some ways to overcome this issue, but up to now, I did not encountered any situation likes this---> it is not necessary to worry about this.
                     *
                     * */
                    do
                    {
                      found_bucket = true;
                      for (j = 0; j < DIMENSION; j++)
                      {
                        if (refc[j] < *(buck2->get_mincrd () + j))//buck2 is the last neighbor
                        {
                          dir[j] = 1;//dir[j] = 1 means that the target bucket should be at the "left" of current bucket. This is a very good idea to search for buckets.
                          found_bucket = false;
                        }
                        else if (refc[j] > *(buck2->get_maxcrd () + j))
                        {
                          dir[j] = 2;//dir[j] = 1 means that the target bucket should be at the "right" of current bucket. This is a very good idea to search for buckets.
                          found_bucket = false;
                        }
                        else
                        {
                          dir[j] = 0;
                        }
                      }
                      // if can't find the bucket, search in "dir" direction
                      if (!found_bucket)//if bucket is not found
                      {
                        buck1 = buck2; //buck1
                        img_proc = buck2->which_neigh_proc (dir);//process id of the bucket which locates in the correct search direction
                        buck2 = (Bucket *) BG_mesh->lookup (buck2->which_neigh (dir));//find pointer to the bucket that the image belongs to.
                        if ((! buck2) && //if the buck can not been found--> the bucket is not on current process or the bucket does not exist at all (exceed the domain).
                            (buck1->which_neigh_proc (dir) != myid)) // Not sure why he need this?
                        {
                          fprintf (stderr, "leaving bucket  behind myid = %d, neigh_proc = %d\n",
                                   myid, img_proc);
                          fprintf (stderr, "%f, %f, %f\n", *(buck->get_mincrd()),
                                  *(buck->get_mincrd()+1), *(buck->get_mincrd()+2));
                          fprintf (stderr, "coord : = {%f, %f, %f}\n", coord[0], coord[1], coord[2]);
                          fprintf (stderr, "intsct : = {%f, %f, %f}\n", intsct[0], intsct[1], intsct[2]);
                          fprintf (stderr, "refc : = {%f, %f, %f}\n", refc[0], refc[1], refc[2]);
                          fprintf (stderr, "key := {%u, %u, %u}\n", p_ghost->getKey().key[0],
                                            p_ghost->getKey().key[1], p_ghost->getKey().key[2]);
                          fflush (stderr);
                          break;
                        }
                      }//end of if statement
                    }
                    while (!found_bucket);//found the bucket to which the image belongs.

                    // if reflection is found, create an image data-structure
                    if ( found_bucket )
                    {
                      for (j = 0; j < KEYLENGTH; j++)
                        buck_key[j] = buck2->getKey ().key[j];

                      for (j = 0; j < TKEYLENGTH; j++)
                        ghst_key[j] = p_ghost->getKey ().key[j];

                      BndImage bnd_img (img_proc, myid, buck_key, ghst_key, refc);
                      Image_table.push_back (bnd_img);
                      p_ghost->set_reflection (true);
                    }

                  }//end if: if particle reflection is not found yet.
        	  }//end if: if particle is wall ghost

        } // end of loop: the loop go through all partiles in certain bucket
    } //end of go through all buckets belong to certain processor: and the bucket has to be MIXED.-----> I am wondering whether MIXED bucket should be considered? What about UNDERGROUND?

  delete itr;
}
