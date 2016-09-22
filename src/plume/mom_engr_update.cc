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
//#include <bucket.h>

#include "particler.h"
#include "constant.h"
#include "sph_header.h"
#include "parameters.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

const double sqrt2 = 1.41421356237310;
const double two_k = 2*lamda_P; //this is the coefficient that in front of

int
mom_engr_update(int myid, THashTable * P_table, HashTable * BG_mesh,
                TimeProps * timeprops, SimProps *simprops)
{
  int i, k;

  vector < TKey > pneighs;
  double dx[DIMENSION], xi[DIMENSION], si[DIMENSION], s_heat[DIMENSION];
  double rhs_v[DIMENSION], unew[NO_OF_EQNS], rhs_e;
  double dwdx[DIMENSION], dwdx_heat[DIMENSION];
  double veli[DIMENSION], velj[DIMENSION], velij[DIMENSION];
  double vis;
  double Fij, Cp_ij, kij, Cp_i;
  double deltae;
  double turb_stress; //turbulent viscosity term

  Particle *pi=NULL;

  // three point gauss-quadrature points
  double gravity[3] = { 0., 0., -g_P};

  // time-step
  double dt = timeprops->dtime;

  //create a hash table iterator instance
  THTIterator *itr = new THTIterator(P_table);


#ifdef DEBUG
   bool do_search = false;
   unsigned keycheck[TKEYLENGTH] = {71888048, 253371, 0};
   unsigned keytemp[TKEYLENGTH] ;

   bool search_bypos = false;
   double range[2*DIMENSION] = {-5500., -5000., -5500., -5000, 6000, 6500};

   bool check_vis = false;
   double vis_2d;
   int phasej;

   bool check_engr = false;
   double engr_thresh = 400000;
#endif

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
		      {
			      cout << "The particle found!" << endl;
			      cout << "energy: "<< pi->get_energy() <<endl;
		      }
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
		  double tempi= pi->get_temperature ();
		  Cp_i = pi->get_specific_heat_p();

		  const double *uvec = pi->get_state_vars();
		  // density must always be positive
		  assert(uvec[0] > 0);
		  double Vi = 1.0 / uvec[0];

#if MOMENTUM_DISCRETIZE ==1
		  double p_external;
          if (pi->which_phase() == 1) //Because when I added the external pressure in the momentum equation, the fluctuation in atmosphere become larger, so use two different forms for two phases here.
        	  p_external=0.0;
          else
        	  p_external=determine_pressure(simprops, *(pi->get_coords()+2));

		  double pvsqi=(pressi-p_external)*Vi*Vi;//p*v^2
#else
          double pvsqi=pressi*Vi*Vi;//p*v^2
#endif

          // velocity veli
          for (k = 0; k < DIMENSION; k++)
#ifdef HAVE_TURBULENCE_LANS
		  for (k = 0; k < DIMENSION; k++)
		      veli[k] = *(pi->get_smoothed_velocity()+k);
#else
          for (k = 0; k < DIMENSION; k++)
               veli[k] = uvec[k+1];
#endif

          double mpvsqij= 0.;   //mpvsqij=mj*(pi/rhoi^2+pj/rhoj^2);


          double heat_tran = 0.; //heat transfer term heat_tran=mj/rho_j (Ti-Tj) F_ij

          // reset rhs to zero, for current particle
		  for (i = 0; i < DIMENSION; i++)
		      rhs_v[i] = 0;

		  rhs_e=0;

		  // list of neighbors
		  vector < TKey > pneighs = pi->get_neighs();
		  vector < TKey >::iterator p_itr;
		  for (p_itr = pneighs.begin(); p_itr != pneighs.end(); p_itr++)
		  {

		      Particle *pj = (Particle *) P_table->lookup(*p_itr);

		      if (!pj)
		    	  cout << "here it is, the missed particles" << endl;
		      assert (pj);

		      // self contribution is zero as dw(0)=0
		      if (*pi == *pj)
		         continue;

		      double dist_sq = 0;

		      for (i = 0; i < DIMENSION; i++)
		      {
		           dx[i] = xi[i] - *(pj->get_coords() + i);
		           si[i] = dx[i] / hi;
		           dist_sq += dx[i] * dx[i];
		      }
		      if (in_support(dx, supp))
		      {
		          double mj = pj->get_mass();
		          double pressj = (pj->get_pressure());
		          double Vj = 1.0 / ((pj->get_density()));

				  const double *uvecj = pj->get_state_vars();
				  // density must always be positive
				  assert(uvecj[0] > 0);

#ifdef HAVE_TURBULENCE_LANS
		          for (k = 0; k < DIMENSION; k++)
		               velj[k] = *(pj->get_smoothed_velocity()+k);
#else
		          for (k = 0; k < DIMENSION; k++)
		          	   velj[k] = uvecj[k+1];
#endif

		          // pre-computing, velocity difference veli-velj
		          for (k = 0; k < DIMENSION; k++)
		               velij[k] = veli[k]- velj[k];

                  //pre-compute, artificial viscosity
		          double rhoab= 0.5 * ((pi->get_density()) + (pj->get_density()));
		          double sndspdab = 0.5 * ((pi->get_sound_speed()) + (pj->get_sound_speed()));
		          vis= art_vis (rhoab, sndspdab, dx, velij, dist_sq, hi);

		          //compute turbulent stress
#ifdef HAVE_TURBULENCE_LANS
		          turb_stress=SPH_epsilon_mom(velij, Vj);
		          vis -=turb_stress; //The minus comes from the fact that there is a minus sign at the front of discretized equation
#endif

#ifdef DEBUG
		          if (check_vis)
		          {
		        	  vis_2d = art_vis_2d (rhoab, sndspdab, &(dx[1]), &(velij[1]), dist_sq, hi);
		        	  phasej = pj->which_phase();
		          }
#endif

#if MOMENTUM_DISCRETIZE ==1
		          mpvsqij = mj*((pressj-p_external) * Vj * Vj + pvsqi + vis);
#else
		          mpvsqij = mj*(pressj * Vj * Vj + pvsqi + vis);
#endif
		          mpvsqij = mj*(pressj * Vj * Vj + pvsqi + vis);

	              for (int ii=0; ii<DIMENSION; ii++)
	            	  s_heat[ii]=si[ii]*HEAT_TRANS_SCALE_RATIO;
		          // pre-compute weight function derivatives
		          for (k = 0; k < DIMENSION; k++)
		          {
		              dwdx[k] = d_weight (si, hi, k);
		              dwdx_heat[k] = d_weight (s_heat, hi/HEAT_TRANS_SCALE_RATIO, k);
		          }

		          // Velocity rhs
		          for (k = 0; k < DIMENSION; k++)
		              rhs_v[k] -= mpvsqij* dwdx[k];

		          // Energy rhs
		          //turbulent heat transfer term
		          heat_tran = 0.; //This line is for the situation when we do not have turbulence in our model
#ifdef HAVE_TURBULENCE_LANS
		          Cp_ij = 0.5* (pj->get_specific_heat_p() + Cp_i);
//		          kij=SPH_epsilon_heat_conductivity(Cp_ij, dx, velij);
		          double hab=0.5*(hi+pj->get_smlen());
		          kij=SPH_epsilon_heat_conductivity(Cp_ij, dx, velij, rhoab, hab, sndspdab);
		          Fij=compute_F(dwdx_heat,dx);
		          heat_tran =mj*Vj*Vi* kij*(tempi- pj->get_temperature ())* Fij;
#endif
		          deltae = 0.;
		          for (k = 0; k < DIMENSION; k++)
		              deltae += 0.5* (mpvsqij* dwdx[k])* velij[k];

		          rhs_e += (deltae + Vi*heat_tran);
		      }
		  } // end loop over neighs

		   //In Dinesh's code, density is also updated here. But it is not necessary in my code, I will update density somewhere else!
		   // density keep unchange, will be updated later!
		   unew[0] = uvec[0];

           // x-velocity
           unew[1] = uvec[1] + dt * (rhs_v[0]);

		   // y-velocity
		   unew[2] = uvec[2] + dt * (rhs_v[1]);

		   // z-velocity
		   unew[3] = uvec[3] + dt * (rhs_v[2] + gravity[2]);

		   // energy
//		   unew[4] = uvec[4] + dt * (rhs_e + veli[2] * gravity[2]); //the variable is only internal energy, it has nnothing to do with mechanical energy, so gravity should not appear here!
		   unew[4] = uvec[4] + dt * rhs_e;
#if HAVE_ENERGY_CUT==1
		   if (unew[4]<=0)
			   unew[4] = ENERGY_CUT;
#endif

#ifdef DEBUG
		  if (check_engr)
		  {
		      if (unew[4]>= engr_thresh) //find the particle that has large energy
			      cout << "The particle found!" << endl;
		  }

#endif


		   pi->put_new_state_vars(unew);
	  } //end of if particle need to up date particle momentum and energy
  }//end of go through all particles

  //pi = (Particle *) P_table->lookup(*ip); // seem useless pointer

  // iterate over hashtable to update state_variables,
  //what need to note is that density was not updated up to this step.
  THTIterator *it2 = new THTIterator(P_table);
  while ((pi = (Particle *) it2->next()))
    if (pi->need_neigh())
    {
      pi->update_state_vars();// Maybe I only need to update it for once
      pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P);
    }

  // clean up
  delete itr;
  delete it2;

  return 0;
}
