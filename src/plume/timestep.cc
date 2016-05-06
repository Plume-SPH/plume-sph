/*
 * timestep.cc
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#ifdef DEBUG
#  include <cstdio>
#  include <cstdlib>
using namespace std;
#endif

#include <algorithm>
#include <cassert>

#include <hashtab.h>
#include <particler.h>
#include <particle.h>
#include <properties.h>

#include "constant.h"
#include "parameters.h"

double
timestep(THashTable * P_table, TimeProps* timeprops)
{
//  int i, j, k;
  double dt, temp;

  dt = 1.0E+10;                 // Initialize to very high value
  THTIterator *itr = new THTIterator(P_table);
  Particle *p_curr = NULL;

//  //before computing time step, second variables need to be updated, because of repartition.
//  while ((p_curr = (Particle *) itr->next ()))
//  {
//        if (p_curr->need_neigh())
//        {
//        	p_curr->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P);
//        }
//  }//end of go through all particles
//
//  itr->reset();
  while ((p_curr = (Particle *) itr->next()))
    if (p_curr->need_neigh())
    {
      // calc speed of sound through the medium
      double c = p_curr->get_sound_speed();

#ifdef DEBUG
      bool check_sndspd = true;
      if (check_sndspd)
    	  if (!(c>0))
    		  cout << "sound speed = " << c << endl;
#endif
      assert (c>0);

      temp = p_curr->get_smlen() / c;

      if (temp < dt)
        dt = temp;
    }

  // delete HT Iterator
  delete itr;

  //set up constrain on time step, such that at most one layer of eruption particle will be added
  double t;
//  t_each = (timeprops->TIME_SCALE)*(timeprops->t_each);
//  double t_each = timeprops->t_each;
//  t= std::min (t_each, 0.2 * dt);
  t=CFL_P*dt;
  return t;
}
