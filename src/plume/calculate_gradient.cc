/*
 * calculate_gradient.cc
 *
 *  Created on: Feb 14, 2017
 *      Author: zhixuanc
 */
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include <vector>
#include <cassert>
#include <cstdio>
using namespace std;

#include <hashtab.h>

#include "particler.h"
#include "constant.h"
#include "sph_header.h"

//Function for computing gradient of all state variables.
void calc_gradients(THashTable *P_table)
{
	  int i, j, k;
	  double xi[DIMENSION], dx[DIMENSION], s[DIMENSION], hi;
	  double mj, rhoj;

	  double dwdx[DIMENSION], dwm_rho[DIMENSION];
      double AA[DIMENSION][DIMENSION], FF[DIMENSION][NO_OF_EQNS], ff[DIMENSION][NO_OF_EQNS];
      double gradU[NOEQxDIM];

      vector<Key> neighs;

	  // create a Hash Table iterator instance
	  THTIterator * itr = new THTIterator (P_table);
	  Particle * pi = NULL;

#ifdef DEBUG
	   bool do_search = true;
	   unsigned keycheck[TKEYLENGTH] = {268315389, 3363091095, 0};
	   unsigned keytemp[TKEYLENGTH];

	   bool search_by_pos =false;
	   double coord_check[DIMENSION]={7350, 8850, 8550};
	   double coord_temp[DIMENSION];

#endif

	  // iterate over the table
	  while ((pi = (Particle *) itr->next ()))
	  {
		//All no-guest particles will need gradient.
	    if (!pi->is_guest ())
	    {

#ifdef DEBUG
			if (do_search)
			{
			    for (i = 0; i < TKEYLENGTH; i++)
				    keytemp[i] = pi->getKey ().key[i];

			    if (find_particle (keytemp, keycheck))
				    cout << "The particle found!" << endl;
			}

			if (search_by_pos)
			{
			    for (i = 0; i < DIMENSION; i++)
			    	coord_temp[i] = *(pi->get_coords ()+i);

			    if (find_particle (coord_temp, coord_check))
				    cout << "The particle found!" << endl;
			}
#endif

	      for (i = 0; i < DIMENSION; i++)
	          xi[i] = (*(pi->get_coords () + i));

#if DENSITY_UPDATE_SML==0
	      hi = pi->get_original_smlen ();
#elif DENSITY_UPDATE_SML==1
	      hi = pi->get_smlen ();
#endif

	      double supp = 3.0 * hi;

	      int nreal=0;

	      // Initialize summation variables to ZERO
	      for (i = 0; i < DIMENSION; i++)
	      {
	        for (j = 0; j < DIMENSION; j++)
	          AA[i][j] = 0.0;

	        for (j = 0; j < NO_OF_EQNS; j++)
	        {
	          FF[i][j] = 0.0;
	          ff[i][j] = 0.0; //It is necessary to make initial value of ff to be zero
	        }
	      }

	      vector <TKey> neighs = pi->get_neighs ();
	      vector <TKey> :: iterator p_itr;

	      //go through all neighbors
	      for (p_itr = neighs.begin (); p_itr != neighs.end (); p_itr++)
	      {
	          Particle *pj = (Particle *) P_table->lookup (*p_itr);
	          assert (pj);

		      // self contribution is zero as dw(0)=0
		      if (*pi == *pj)
		         continue;

	          // All particles should be included for updating gradient ..
	          //if (pj->contr_dens())
	          //{

			  // the boundary deficiency and interface deficiency is handled by wnorm!
			  for (i = 0; i < DIMENSION; i++)
				  dx[i] = xi[i] - *(pj->get_coords() + i);

			  if (in_support(dx, supp))
			  {

				  nreal++;
				  mj=pj->get_mass();
				  rhoj=pj->get_density();

				  for (k = 0; k < DIMENSION; k++)
					  s[k] = dx[k] / hi;

				  for (k = 0; k < DIMENSION; k++)
				  {
					  dwdx[k] = d_weight (s, hi, k);
					  dwm_rho[k] = dwdx[k] * mj/rhoj; //dwdx * mj * rhoj
				  }

				  // Get matrix A
				  for (i = 0; i < DIMENSION; i++)
					for (j = 0; j < DIMENSION; j++)
						AA[i][j] -= dx[i]*dwm_rho[j];

				  for (i = 0; i < DIMENSION; i++)
				  {
					  FF[i][0]-=((pi->get_density() - pj->get_density())*dwm_rho[i]);
					  FF[i][1]-=((*(pi->get_vel()) - *(pj->get_vel()))*dwm_rho[i]);
					  FF[i][2]-=((*(pi->get_vel()+1) - *(pj->get_vel()+1))*dwm_rho[i]);
					  FF[i][3]-=((*(pi->get_vel()+2) - *(pj->get_vel()+2))*dwm_rho[i]);
					  FF[i][4]-=((pi->get_pressure() - pj->get_pressure())*dwm_rho[i]);
				  }
			   }
	          //}//end of if particle will contribute to the density.

	      }//end of go through all neighbors

	      if (nreal <= 20)
	         for (i = 0; i < NOEQxDIM; i++)
	        	 **(ff+i) = 0.;
	      else
	    	 linsolve(&AA[0][0], &FF[0][0], NO_OF_EQNS, &ff[0][0]);


	      //put result in a one dimensional array for the convenience of updating!
	      // check and save derivatives
	      for (i = 0; i < NO_OF_EQNS; i++)
	        for (j = 0; j < DIMENSION; j++)
	          gradU[i * DIMENSION + j] = ff[j][i];


	      //  if there is a problem ... then die a violent death
//	      for (i = 0; i < NOEQxDIM; i++)
//	        if (isnan(gradU[i]) && (nreal > 20))
//	        {
//	          fprintf(stderr,"FATAL ERROR: calc_gradients() failed\n");
//	          fprintf(stderr,".. at (%f, %f, %f), no of neighbors: %d \n",
//	                          xi[0], xi[1], xi[2], nreal);
//
//	          exit (EXIT_FAILURE);
//	        }

		  pi->put_density_d(gradU);
		  pi->put_velocity_u_d(gradU+DIMENSION);
		  pi->put_velocity_v_d(gradU+2*DIMENSION);
		  pi->put_velocity_w_d(gradU+3*DIMENSION);
		  pi->put_pressure_d(gradU+4*DIMENSION);

		}//end of if --> if the density of that particle need to be updated based summation
	  }//end of loop -->go through all particle


	  // clean up stuff
	  delete itr;

	  return;
}
