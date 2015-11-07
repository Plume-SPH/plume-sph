/*
 * smooth_velocity.cc
 *
 *  Created on: Nov 5, 2015
 *      Author: zhixuanc
 */

//function that will smooth velocity ---> in SPH-epsilon model

#include <vector>
#include <cassert>
#include <iostream>
using namespace std;

#include <thashtab.h>

#include "particler.h"
#include "constant.h"
#include "sph_header.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

void
smooth_velocity(THashTable * P_table)
{
  int i, j, k, no_of_neighs;
  //unsigned jkey[TKEYLENGTH];
  double xi[DIMENSION], ds[DIMENSION], s[DIMENSION];
  double wght, temp_dv[DIMENSION];
  double vi[DIMENSION], vj[DIMENSION];
  TKey tmpkey;
  double smoothed_v[DIMENSION];
  double wnorm, mvw;
  // create a Hash Table iterator instance
  THTIterator * itr = new THTIterator (P_table);
  Particle * pi = NULL;

#ifdef DEBUG
   bool check_den = false;
   bool do_search = false;
   bool check_mssfrac = false;
   bool find;
   unsigned keycheck[TKEYLENGTH] =  {70615689, 3732101160, 0};
   unsigned keytemp[TKEYLENGTH] ;
#endif

  // iterate over the table
  while ((pi = (Particle *) itr->next ()))
  {
    if (pi->need_neigh()) //Only real && no-guest particle need their velocity to be smoothed
    {

#ifdef DEBUG
		if (do_search)
		{
		    for (i = 0; i < TKEYLENGTH; i++)
			    keytemp[i] = pi->getKey ().key[i];

		    if (find_particle (keytemp, keycheck))
			    cout << "The particle found!" << endl;
		}
#endif

      for (i = 0; i < DIMENSION; i++)
        xi[i] = (*(pi->get_coords () + i));

      for (i = 0; i < DIMENSION; i++)
        vi[i] = (*(pi->get_vel() + i));

      double hi = pi->get_smlen ();
      double supp = 3.0 * hi;
      TKey ki = pi->getKey ();

      for (i = 0; i < DIMENSION; i++)
    	  temp_dv[i] = 0.;

      wnorm = 0.;

      vector <TKey> neighs = pi->get_neighs ();
      vector <TKey> :: iterator p_itr;

      //go through all neighbors
      for (p_itr = neighs.begin (); p_itr != neighs.end (); p_itr++)
      {
          // Transfer key into particle #
          Particle *pj = (Particle *) P_table->lookup (*p_itr);
          assert (pj);

          //update velocity
      	  // guests are included ghosts are not --> the wall effect is handled by wnorm!
          // It is necessary to double think whether normalization is necessary for velocity smoothing
      	  if (pj->contr_vel_smooth()) //Double think is necessary about which ghost particles should be considered here
      	  {
      		 for (i = 0; i < DIMENSION; i++)
      		      ds[i] = xi[i] - *(pj->get_coords() + i);

      		 if (in_support(ds, supp))
      		 {
      		     for (k = 0; k < DIMENSION; k++)
      		        s[k] = ds[k] / hi;
      		      wght = weight(s, hi);
      		      mvw =  wght * (pj->get_mass()) / pj->get_density();

      		      for (i = 0; i < DIMENSION; i++)
      		            vj[i] = (*(pj->get_vel() + i));
      		      for (i = 0; i < DIMENSION; i++)
      		         temp_dv[i] += mvw * (vj[i] - vi[i]);
    		      wnorm += mvw;
      		  }
      	    }
      }//end of go through all neighbors

      assert(wnorm > 0);
      for (i = 0; i < DIMENSION; i++)
    	  temp_dv[i]  /= wnorm;
      for (i = 0; i < DIMENSION; i++)
    	  smoothed_v[i] = vi[i] + temp_dv[i];

      pi->put_smoothed_velocity(smoothed_v);

    }//end of if --> if the velocity of that particle need to be updated based summation
  }//end of loop -->go through all particle

  // clean up stuff
  delete itr;

  return;
}


