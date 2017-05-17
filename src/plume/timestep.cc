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
#include <cmath>

#include <hashtab.h>
#include <particler.h>
#include <particle.h>
#include <properties.h>

#include "constant.h"
#include "parameters.h"

using namespace std;

double
pressure_bc_step(HashTable * BG_mesh, THashTable * P_table)
{
	double dt, temp;
    dt = 1.0E+10;

	double p_spd_u, p_spd_v, p_spd_w, p_spd;

	HTIterator * itr = new HTIterator (BG_mesh);
	Bucket * Bnd_buck;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

	vector < TKey > plist;
	vector < TKey >::iterator p_itr;

	while ((tempptr=itr->next ()))
	{
		breif_buck = (BriefBucket *) tempptr;
		if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
			  continue;
		else
		{
			Bnd_buck = (Bucket*) tempptr;
			// 1) Bnd_buck has involved flag = 1 --->Only has potential involved particles
			// 2) is not guest
		    if ((Bnd_buck->get_has_involved()>=1) && !Bnd_buck->is_guest()) // should be (Bnd_buck->get_has_involved()>=1), previously it was (Bnd_buck->get_has_involved()==1)
			{
		    	plist = Bnd_buck->get_plist();
                //assert(plist.size()>0); //make sure Bnd_buck is not empty ---> not necessary! Some bucket might be empty due to particle movement
		    	if (plist.size()>0)
			    	for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
			    	{
			    		 Particle *p_curr = (Particle *) P_table->lookup(*p_itr);
			    		 assert(p_curr);
#if HAVE_TURBULENCE_LANS !=0
			    		 p_spd_u = *(p_curr->get_smoothed_velocity());
			    		 p_spd_v = *(p_curr->get_smoothed_velocity()+1);
			    		 p_spd_w = *(p_curr->get_smoothed_velocity()+2);
#else
			    		 p_spd_u = *(p_curr->get_vel());
			    		 p_spd_v = *(p_curr->get_vel()+1);
			    		 p_spd_w = *(p_curr->get_vel()+2);
#endif
			    		 p_spd = max((abs(p_spd_u)>abs(p_spd_v)) ? abs(p_spd_u):abs(p_spd_v), abs(p_spd_w));

			    		 temp = p_curr->get_smlen()/p_spd;

			    		 if (temp < dt)
			    		    dt = temp;
			    	}// loop go through all particles in the bucket

			}//end of if bucket only contains potential involved particles
		}
	}//end of go through all buckets

	delete itr;

	return dt;
}

//main function determines the time step
double
timestep(HashTable * BG_mesh, THashTable * P_table, TimeProps* timeprops)
{
  double dt, temp;

  dt = 1.0E+10;                 // Initialize to very high value
  THTIterator *itr = new THTIterator(P_table);
  Particle *p_curr = NULL;

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
#if TIME_STEP_CONSTRAIN ==1
  double dt_bc = pressure_bc_step (BG_mesh, P_table);
  t = min (CFL_P*dt, dt_bc*CFL_BC_P);
#elif TIME_STEP_CONSTRAIN ==0
  t = CFL_P*dt;
#endif

  return t;
}

#if CODE_DIMENSION==1
//main function determines the time step
double
timestep(THashTable * P_table, TimeProps* timeprops)
{
  double dt, temp;

  dt = 1.0E+10;                 // Initialize to very high value
  THTIterator *itr = new THTIterator(P_table);
  Particle *p_curr = NULL;

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
  t = CFL_P*dt;


  return t;
}
#endif //CODE_DIMENSION==1
