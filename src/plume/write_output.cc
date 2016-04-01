/*
 * write_output.cc
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <cstdio>
#include <cassert>
#include <vector>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>

#include <properties.h>
#include <particler.h>
#include <outforms.h>
#include "constant.h"
#include "sph_header.h"

void
write_output(int myid, int numprocs,
             THashTable * P_table, HashTable * BG_mesh,
             vector <BucketHead> & partition_table,
             TimeProps * timeprops, int format)
{

  if (format & 1)
    write_h5part(myid, numprocs, P_table, timeprops);

  if (format & 2)
    write_matlab(myid, P_table, BG_mesh, timeprops, partition_table);

  return;
}

void
write_output_show(int myid, int numprocs,
             THashTable * P_table, HashTable * BG_mesh,
             vector <BucketHead> & partition_table,
             TimeProps * timeprops, int format)
{

  if (format & 1)
    write_h5part_show(myid, numprocs, P_table, timeprops);

  if (format & 2)
    write_matlab(myid, P_table, BG_mesh, timeprops, partition_table);

  return;
}

void
write_debug_info(int myid, THashTable * P_table, int index)
{

  char fname[20];

  sprintf(fname, "Step%02d%06d.dat", myid, index);
  THTIterator *itr = new THTIterator(P_table);
  Particle *p_curr = NULL;

  while ((p_curr = (Particle *) itr->next()))
    if ((p_curr->getKey().key[0] == 1748747027) && (p_curr->need_neigh()))
    {
      FILE *fp = fopen(fname, "w");

      vector < TKey > neighs = p_curr->get_neighs();
      vector < TKey >::iterator itr_p;
      fprintf(fp, "%e, %e, %e, %e, %e, %e, %e, %e\n",
              *(p_curr->get_coords()),
              *(p_curr->get_coords() + 1),
              *(p_curr->get_coords() + 2),
              *(p_curr->get_state_vars()),
              *(p_curr->get_state_vars() + 1),
              *(p_curr->get_state_vars() + 2),
              *(p_curr->get_state_vars() + 3),
              *(p_curr->get_state_vars() + 4));
      for (itr_p = neighs.begin(); itr_p != neighs.end(); itr_p++)
      {
        Particle *p_neigh = (Particle *) P_table->lookup(*itr_p);

        assert(p_neigh);
        if (*p_curr == *p_neigh)
          continue;
        fprintf(fp, "%e, %e, %e, %e, %e, %e, %e, %e \n",
                *(p_neigh->get_coords()),
                *(p_neigh->get_coords() + 1),
                *(p_neigh->get_coords() + 2),
                *(p_neigh->get_state_vars()),
                *(p_neigh->get_state_vars() + 1),
                *(p_neigh->get_state_vars() + 2),
                *(p_neigh->get_state_vars() + 3),
                *(p_neigh->get_state_vars() + 4)
                );
      }
      fclose(fp);
    }
}
