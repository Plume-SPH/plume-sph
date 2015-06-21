/*
 * mom_engr_update.cc
 *
 *  Created on: Feb 24, 2015
 *      Author: zhixuanc
 */

#include <vector>
#include <cmath>
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>

#include "particler.h"
#include "constant.h"
#include "sph_header.h"
#include "parameters.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

const double sqrt2 = 1.41421356237310;
const int NUM_RHS=DIMENSION+1;

int
mom_engr_update(int myid, THashTable * P_table, HashTable * BG_mesh,
                TimeProps * timeprops)
{
  int i, k;

  vector < TKey > pneighs;
  double dx[DIMENSION], xi[DIMENSION], si[DIMENSION];
  double rhs_v[DIMENSION], unew[NO_OF_EQNS], rhs_e;
  double dwdx[DIMENSION];
  double veli[DIMENSION], velj[DIMENSION], velij[DIMENSION];
  double vis;
  double deltae;
  Particle *pi=NULL;

  // three point gauss-quadrature points
  double gravity[3] = { 0., 0., -9.81};

  // time-step
  double dt = timeprops->dtime;
  //create a hash table iterator instance
  THTIterator *itr = new THTIterator(P_table);


#ifdef DEBUG
   bool do_search = false;
   bool find;
   unsigned keycheck[TKEYLENGTH] = {92878763, 2157710878, 0};
   unsigned keytemp[TKEYLENGTH] ;

   bool search_bypos = false;
   double range[2*DIMENSION] = {-5500., -5000., -5500., -5000, 6000, 6500};
#endif
//
//  //before moment and energy update, update secondary variables for guest particles
//  while ((pi = (Particle *) itr->next ()))
//  {
//      if ((pi->is_guest()))
//      {
//          pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P);
//      }
//  }//end of go through all particles

  // go through particle table again
//  itr->reset();

  while ((pi = (Particle *) itr->next ()))
  {
	  if (pi->need_neigh ())
	  {
		  for (i = 0; i < DIMENSION; i++)
			  xi[i] = *(pi->get_coords() + i);

#ifdef DEBUG
		  if (do_search)
		  {
		      for (i = 0; i < TKEYLENGTH; i++)
			      keytemp[i] = pi->getKey ().key[i];

		      if (find_particle (keytemp, keycheck))
			      cout << "The particle found!" << endl;
		  }

		  if (search_bypos)
		  {
			  if (find_particle_pos_range(xi, range))
				  cout << "The particle found!" << endl;
		  }
#endif

		  // expanded smoothing length for Momentum equation
		  double hi = pi->get_smlen();
		  double supp = 3 * hi;
		  double pressi = (pi->get_pressure());

		  const double *uvec = pi->get_state_vars();
		  // density must always be positive
		  assert(uvec[0] > 0);
		  double Vi = 1.0 / uvec[0];
          double pvsqi=pressi*Vi*Vi;//p*v^2

          // velocity veli
          for (k = 0; k < DIMENSION; k++)
               veli[k] = uvec[k+1];

          double mpvsqij= 0.;   //mpvsqij=mj*(pi/rhoi^2+pj/rhoj^2);
		  // reset rhs to zero, for current particle
	//	  double wnorm = 0.;

		  for (i = 0; i < DIMENSION; i++)
		      rhs_v[i] = 0;
		  rhs_e=0;

		  // list of neighbors
		  vector < TKey > pneighs = pi->get_neighs();
		  vector < TKey >::iterator p_itr;
		  for (p_itr = pneighs.begin(); p_itr != pneighs.end(); p_itr++)
		  {

		      Particle *pj = (Particle *) P_table->lookup(*p_itr);
		      assert (pj);

		      // self contribution is zero as dw(0)=0
		      if (*pi == *pj)
		         continue;

		      double dist_sq = 0;
//		      double dist = 0;

		      for (i = 0; i < DIMENSION; i++)
		      {
		           dx[i] = xi[i] - *(pj->get_coords() + i);
		           si[i] = dx[i] / hi;
		           dist_sq += dx[i] * dx[i];
		      }
//		      dist = sqrt(dist_sq);
		      // if dx < 3 * sqrt(2) * h
		      if (in_support(dx, supp))
		      {
		          double mj = pj->get_mass();
		          double pressj = (pj->get_pressure());
		          double Vj = 1.0 / ((pj->get_density()));

				  const double *uvecj = pj->get_state_vars();
				  // density must always be positive
				  assert(uvecj[0] > 0);
		          // velocity veli
		          for (k = 0; k < DIMENSION; k++)
		               velj[k] = uvecj[k+1];

		          // pre-computing, velocity difference veli-velj
		          for (k = 0; k < DIMENSION; k++)
		               velij[k] = veli[k]- velj[k];

                  //pre-compute, artificial viscosity
		          double rhoab= 0.5 * ((pi->get_density()) + (pj->get_density()));
		          double sndspdab = 0.5 * ((pi->get_sound_speed()) + (pj->get_sound_speed()));
		          vis= art_vis (rhoab, sndspdab, dx, velij, dist_sq, hi);

		          //pre-computing
		          mpvsqij = mj*(pressj * Vj * Vj + pvsqi + vis);
		          // pre-compute weight function derivatives
		          for (k = 0; k < DIMENSION; k++)
		              dwdx[k] = d_weight (si, hi, k);
		          // Velocity rhs
		          for (k = 0; k < DIMENSION; k++)
		              rhs_v[k] -= mpvsqij* dwdx[k];
		          // Energy rhs
		          deltae = 0.;
		          for (k = 0; k < DIMENSION; k++)
		              deltae += 0.5* (mpvsqij* dwdx[k])* velij[k];
		          rhs_e += deltae;

		      }
		  } // end loop over neighs

		   //In Dinesh's code, density is also updated here. But it is not necessary in my code, I will update density somewhere else!
		   // density keep unchange, will be updated later!
		   unew[0] = uvec[0];

           //  x-velocity
           unew[1] = uvec[1] + dt * (rhs_v[0]);

		   // y-velocity
		   unew[2] = uvec[2] + dt * (rhs_v[1]);

		   // z-velocity
		   unew[3] = uvec[3] + dt * (rhs_v[2] + gravity[2]);

		   // energy
		   unew[4] = uvec[4] + dt * (rhs_e + veli[2] * gravity[2]);

		   pi->put_new_state_vars(unew);
	  }
  }

  //pi = (Particle *) P_table->lookup(*ip); // seem useless pointer

  // iterate over hashtable to update state_variables,
  //what need to note is that density was not updated up to this step.
  THTIterator *it2 = new THTIterator(P_table);
  while ((pi = (Particle *) it2->next()))
    if (pi->need_neigh())
    {
      pi->update_state_vars();// Maybe I only need to update it for once
      pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P);
    }

  // clean up
  delete itr, it2;

  return 0;
}
