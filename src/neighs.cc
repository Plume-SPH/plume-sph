/*
 * neighs.cc
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */


#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <particle.h>
#include <particler.h>

#include "constant.h"
#include "sph_header.h"
#include "parameters.h"

/*
 * First: search for neighbors --->to guarantee that smooth_update is using the most recent information;
 * Second: update smooth length;
 * Third: search for neighbors;
 * This will be called in the main loop of SPH update in time;
 */
int
search_neighs (int myid, THashTable * P_table, HashTable * BG_mesh)
{
  int i, j, k;

  double H_star, H_old, H_new;
  double rho_star;
  double supp ;
  double mj;
  Particle * pi = NULL ;
  Particle * pj = NULL ;
  double dx[DIMENSION], xi[DIMENSION], si[DIMENSION];
  double wght;
  THTIterator *itr2 = new THTIterator(P_table);

  double hi;
  Key *neighbors;//bucket key
  // create a Hash Table iterators instance
  HTIterator *igrd = new HTIterator(BG_mesh);
  // iterate over the bucket table
  Bucket *curr_bucket = NULL;

  while ((curr_bucket = (Bucket *) igrd->next()))
  {
    assert(curr_bucket);
    if (!curr_bucket->is_guest() && (curr_bucket->get_plist().size() > 0))
    {
      vector < TKey > plist = curr_bucket->get_plist();
      vector < TKey >::iterator itr;
      for (itr = plist.begin(); itr != plist.end(); itr++)
      {
        pi = (Particle *) P_table->lookup(*itr);
        assert(pi);

        // 3*sqrt(2) < 4.25
        hi = 4.25 * pi->get_smlen();
        for (int k = 0; k < DIMENSION; k++)
          xi[k] = *(pi->get_coords() + k);

        // all particles in current bucket are neighbors
		vector < TKey > pneighs = pi->get_neighs();
		vector < TKey >::iterator p_itr;
        pneighs.clear();
        pneighs.insert(pneighs.begin(), plist.begin(), plist.end());

        // search the neighboring buckets for neighbors
        neighbors = curr_bucket->get_neighbors();
        const int *neigh_proc = curr_bucket->get_neigh_proc();

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

              // get difference of position vectors between i,j
              double ds[DIMENSION];
              for (j = 0; j < DIMENSION; j++)
                ds[j] = xi[j] - *(pj->get_coords() + j);

              // if within support, add to neigh list
              if (in_support(ds, hi))
                pneighs.push_back(*it2);
            }
          }
        pi->put_neighs(pneighs);
      }

    }// end of if :bucket is not guest and has particles
  } //end of while


  /*
   * before searching for neighbors, adjust smooth length
   * The algorithm is a simplified version of the one introduced in Shu-ichiro Inutsuka's paper:
   * 1) The algorithm is excuted only for one loop instead of a iteration process
   * 2) I did not search for neighbors for the new smooth length h*, just used the old neighbor information.
   * 3) Consider all phases while computing rho_star;
   * 4) Consider even ghost particles while sum up ---> Ghost particles also contribute to momentum and energy update.
   * 5) Only the smooth length of particles who need neighbor information will be updated!
   */

  while ((pi = (Particle *) itr2->next ()))
   {
 	  if (pi->need_neigh ())
 	  {
		  for (i = 0; i < DIMENSION; i++)
			  xi[i] = *(pi->get_coords() + i);

		  // expanded smoothing length for Momentum equation
		  H_old = pi->get_smlen();
		  H_star = H_old *  C_smooth_P;
		  supp = 3 * H_star;

		  vector < TKey > pneighs = pi->get_neighs();
		  vector < TKey >::iterator p_itr;
 		  rho_star = 0.;
 		  for (p_itr = pneighs.begin(); p_itr != pneighs.end(); p_itr++)
 		  {
 		      pj = (Particle *) P_table->lookup(*p_itr);
 		      assert (pj);

 		      for (i = 0; i < DIMENSION; i++)
 		      {
 		           dx[i] = xi[i] - *(pj->get_coords() + i);
 		           si[i] = dx[i] / H_star;
 		      }

 		      if (in_support(dx, supp))
 		      {
 		          mj = pj->get_mass();
 		          wght = weight(si, H_star);
 		          rho_star += wght * mj;
 		      }
 		  }// end loop over neighs

 		   switch (DIMENSION)
 		   {
 		   case 3:
 			  H_new = eta_smooth_P * cbrt(pi->get_mass()/rho_star);
 			  break;
 		   case 2:
 			  H_new = eta_smooth_P * sqrt(pi->get_mass()/rho_star);
 			  break;
 		   case 1:
 			  H_new = eta_smooth_P * (pi->get_mass()/rho_star);
 			  break;
 		   }
 		   pi->put_smlen(H_new);
 	  }
   }

  igrd->reset();
  while ((curr_bucket = (Bucket *) igrd->next()))
  {
    assert(curr_bucket);
    if (!curr_bucket->is_guest() && (curr_bucket->get_plist().size() > 0))
    {
      vector < TKey > plist = curr_bucket->get_plist();
      vector < TKey >::iterator itr;
      for (itr = plist.begin(); itr != plist.end(); itr++)
      {
        pi = (Particle *) P_table->lookup(*itr);
        assert(pi);

        // 3*sqrt(2) < 4.25
        hi = 4.25 * pi->get_smlen();
        for (int k = 0; k < DIMENSION; k++)
          xi[k] = *(pi->get_coords() + k);

        // all particles in current bucket are neighbors
		vector < TKey > pneighs = pi->get_neighs();
		vector < TKey >::iterator p_itr;
        pneighs.clear();
        pneighs.insert(pneighs.begin(), plist.begin(), plist.end());

        // search the neighboring buckets for neighbors
        neighbors = curr_bucket->get_neighbors();
        const int *neigh_proc = curr_bucket->get_neigh_proc();

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

              // get difference of position vectors between i,j
              double ds[DIMENSION];
              for (j = 0; j < DIMENSION; j++)
                ds[j] = xi[j] - *(pj->get_coords() + j);

              // if within support, add to neigh list
              if (in_support(ds, hi))
                pneighs.push_back(*it2);
            }
          }
        pi->put_neighs(pneighs);
      }

    }// end of if :bucket is not guest and has particles
  } //end of while

  // delete iterator
  delete igrd;
//  delete itr;
//  delete p_itr;
  delete itr2;
//  delete it2;

  return 0;
}

// update neighbor information without updating smooth length --> Will be called at the time of particle set up;
int
search_neighs_consth (int myid, THashTable * P_table, HashTable * BG_mesh)
{
  int i, j, k;

  Particle * pi = NULL ;
  Particle * pj = NULL ;
  double dx[DIMENSION], xi[DIMENSION], si[DIMENSION];
  double wght;

  double hi;
  Key *neighbors;//bucket key
  // create a Hash Table iterators instance
  HTIterator *igrd = new HTIterator(BG_mesh);

  // iterate over the bucket table
  Bucket *curr_bucket = NULL;

  while ((curr_bucket = (Bucket *) igrd->next()))
  {
    assert(curr_bucket);
    if (!curr_bucket->is_guest() && (curr_bucket->get_plist().size() > 0))
    {
      vector < TKey > plist = curr_bucket->get_plist();
      vector < TKey >::iterator itr;
      for (itr = plist.begin(); itr != plist.end(); itr++)
      {
        pi = (Particle *) P_table->lookup(*itr);
        assert(pi);

        // 3*sqrt(2) < 4.25
        hi = 4.25 * pi->get_smlen();
        for (int k = 0; k < DIMENSION; k++)
          xi[k] = *(pi->get_coords() + k);

        // all particles in current bucket are neighbors
		vector < TKey > pneighs = pi->get_neighs();
		vector < TKey >::iterator p_itr;
        pneighs.clear();
        pneighs.insert(pneighs.begin(), plist.begin(), plist.end());

        // search the neighboring buckets for neighbors
        neighbors = curr_bucket->get_neighbors();
        const int *neigh_proc = curr_bucket->get_neigh_proc();

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

              // get difference of position vectors between i,j
              double ds[DIMENSION];
              for (j = 0; j < DIMENSION; j++)
                ds[j] = xi[j] - *(pj->get_coords() + j);

              // if within support, add to neigh list
              if (in_support(ds, hi))
                pneighs.push_back(*it2);
            }
          }
        pi->put_neighs(pneighs);
      }

    }// end of if :bucket is not guest and has particles
  } //end of while

  // delete iterator
  delete igrd;
//  delete itr;
//  delete p_itr;
//  delete it2;

  return 0;
}
