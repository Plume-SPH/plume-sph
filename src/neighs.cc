/*
 * neighs.cc
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */


#include <vector>
#include <iostream>
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <particler.h>

#include "constant.h"
#include "sph_header.h"

int
search_neighs (int myid, THashTable * P_table, HashTable * BG_mesh)
{
  int i, j;
  //unsigned tempkey[KEYLENGTH];
  double xi[DIMENSION];
  double hi;

  vector < TKey > pneighs;//particle key
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
        Particle * pi = (Particle *) P_table->lookup(*itr);
        assert(pi);

        // 3*sqrt(2) < 4.25
        hi = 4.25 * pi->get_smlen();
        for (int k = 0; k < DIMENSION; k++)
          xi[k] = *(pi->get_coords() + k);

        // all particles in current bucket are neighbors
        pneighs.clear();

        //  all particles of current bucket are neighbors
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
              Particle * pj = (Particle *) P_table->lookup(*it2);
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
  }

  // delete iterator
  delete igrd;
  //delete it2, itr;

  return 0;
}
