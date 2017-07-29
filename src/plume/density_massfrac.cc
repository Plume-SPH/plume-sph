/*
 * density_massfrac.cc
 *
 *  Created on: Feb 24, 2015
 *      Author: zhixuanc
 */

// Function for density and mass fraction updating

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
smooth_density(THashTable * P_table)
{
  int i, k;
  double xi[DIMENSION], ds[DIMENSION], s[DIMENSION], hi;
  double wght, density=0.0;
  TKey tmpkey;
  int phs_i;
  double tmprho[PHASE_NUM]={0.0}, phaserho[PHASE_NUM]={0.0}, mssfrc;
  double wnorm[PHASE_NUM], norm;
  double wm, mj, rhoj;

#if DENSITY_UPDATE_SPH==13
   double mi;
#endif

#if HAVE_ENERGY_SMOOTH ==1
  double engr, wm_e, wnorm_e, s_e[DIMENSION];
#elif HAVE_ENERGY_SMOOTH ==2
  double engr, wm_e, s_e[DIMENSION];
#endif

#if ADAPTIVE_SML==2
   double   hi_ori, tmprho_ori[PHASE_NUM], s_ori_i[DIMENSION], s_ori_j[DIMENSION], wght_ori, density_ori;
#endif

  // create a Hash Table iterator instance
  THTIterator * itr = new THTIterator (P_table);
  Particle * pi = NULL;

#ifdef DEBUG
   bool check_den = true;
   bool do_search = false;
   bool check_mssfrac = false;
   unsigned keycheck[TKEYLENGTH] =  {224, 0, 0};
   unsigned keytemp[TKEYLENGTH] ;
#endif

  // iterate over the table
  while ((pi = (Particle *) itr->next ()))
  {
    if (pi->need_neigh())
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

#if DENSITY_UPDATE_SML==0
      hi = pi->get_original_smlen ();
#elif (DENSITY_UPDATE_SML==1) || (DENSITY_UPDATE_SML==11) || (DENSITY_UPDATE_SML==12) || (DENSITY_UPDATE_SML==13)
      hi = pi->get_smlen ();
#endif

#if ADAPTIVE_SML==2
      hi_ori=pi->get_original_smlen ();
#endif

      double supp = CUTOFF2 * hi;   // usually cut_off at 3, considering we might use a different h in place of hi, use a larger sml to guarantee number of particles are enough ---> Probably, I should use a even larger number: CUTOFF

      for (i = 0; i < PHASE_NUM; i++)
    	  tmprho[i] = 0.;

#if ADAPTIVE_SML==2
     for (i = 0; i < PHASE_NUM; i++)
		  tmprho_ori[i] = 0.;
#endif

      for (i = 0; i < PHASE_NUM; i++)
    	  phaserho[i] = 0.;

//#if ADAPTIVE_SML==2
//     for (i = 0; i < PHASE_NUM; i++)
//    	 phaserho_ori[i] = 0.;
//#endif

      for (i = 0; i < PHASE_NUM; i++)
    	  wnorm [i] = 0.;

#if HAVE_ENERGY_SMOOTH == 1
      wnorm_e = 0.0;
      engr = 0.0;
#elif HAVE_ENERGY_SMOOTH == 2
      engr = 0.0;
#endif

#if DENSITY_UPDATE_SPH==13
     mi=pi->get_mass();
#endif

      vector <TKey> neighs = pi->get_neighs ();
      vector <TKey> :: iterator p_itr;

      //go through all neighbors
      for (p_itr = neighs.begin (); p_itr != neighs.end (); p_itr++)
      {
          Particle *pj = (Particle *) P_table->lookup (*p_itr);
          assert (pj);

          // guests are included ghosts are not,
          if (pj->contr_dens()) //In my old code, I did not have this so I will get very large density at the very beginning of the eruption.
          //update density for each phase
          //Probably, a switch case is better than for loop in terms of efficiency!
          {
        	  mj=pj->get_mass();
        	  rhoj=pj->get_density();
			  for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
			  {
				  //for each phase
				  if (pj->which_phase()==phs_i)
				  {

					// the boundary deficiency and interface deficiency is handled by wnorm!
					for (i = 0; i < DIMENSION; i++)
						ds[i] = xi[i] - *(pj->get_coords() + i);

					if (in_support(ds, supp))
					{

#if DENSITY_UPDATE_SML==11
						hi=0.5*(hi+pj->get_smlen ()); //Overwrite original hi with 0.5*(hi+hj)
#elif DENSITY_UPDATE_SML==12
						hi=pj->get_smlen ();          //Overwrite original hi with hj
#endif

						for (k = 0; k < DIMENSION; k++)
						    s[k] = ds[k] / hi;

#if HAVE_ENERGY_SMOOTH==1
						for (k = 0; k < DIMENSION; k++)
							s_e[k] = E_SMOOTH_RATIO*ds[k] / hi;
						wm_e=weight(s_e, hi/E_SMOOTH_RATIO) * mj;
						wnorm_e +=wm_e/rhoj;
						engr +=wm_e*pj->get_energy()/rhoj;
#elif HAVE_ENERGY_SMOOTH==2
						for (k = 0; k < DIMENSION; k++)
							s_e[k] = E_SMOOTH_RATIO*ds[k] / hi;
						wm_e=weight(s_e, hi/E_SMOOTH_RATIO) * mj;
						engr +=wm_e*pj->get_energy()/rhoj;
#endif

#if DENSITY_UPDATE_SML==13
						wght = weight(s, hi);
						double hj=pj->get_smlen ();
						for (k = 0; k < DIMENSION; k++)
							s[k] = ds[k] / hj;
						wght += weight(s, hj);
						wght = 0.5*wght;
#else
						wght = weight(s, hi);
#endif // DENSITY_UPDATE_SML==13

#if DENSITY_UPDATE_SPH==13
						wm = wght * mi;
#else
						wm = wght * mj;
#endif
						tmprho[phs_i-1] += wm;

#if ADAPTIVE_SML==2
						for (k = 0; k < DIMENSION; k++)
						    s_ori_i[k] = ds[k]/hi_ori;
						for (k = 0; k < DIMENSION; k++)
						    s_ori_j[k] = ds[k]/(pj->get_original_smlen ());

						wght_ori = 0.5*(weight(s_ori_i,hi_ori)+weight(s_ori_j,pj->get_original_smlen ()));
						tmprho_ori[phs_i-1] += wght_ori*mj;
#endif

#if (DENSITY_UPDATE_SPH == 1) || (DENSITY_UPDATE_SPH == 22)
						//View multiple phase SPH as one set of discretized point
						wnorm[phs_i-1] += wm / rhoj;
#elif (DENSITY_UPDATE_SPH == 2) || (DENSITY_UPDATE_SPH == 12)
						//View multiple phase SPH as two independent sets of discretized point
						double phase_des=*(pj->get_phase_density ()+ phs_i -1);
						if (phase_des>0)
						    wnorm[phs_i-1] += wm /phase_des;
#endif
					 }
				  }//end of if pj->which_phase()==phs_i

			  }//end of loop for phase
          }//end of if particle will contribute to the density.

      }//end of go through all neighbors

      density=0.0;
      norm = 0.0;

#if ADAPTIVE_SML==2
      density_ori=0.0;
#endif

#if HAVE_ENERGY_SMOOTH==1 //Also do normalization
      assert (wnorm_e > 0);
      engr = engr/wnorm_e;
      pi->put_new_energy(engr);
#elif HAVE_ENERGY_SMOOTH==2 //Do not normalize it
      pi->put_new_energy(engr);
#endif

#if (DENSITY_UPDATE_SPH == 1) || (DENSITY_UPDATE_SPH == 12) || (DENSITY_UPDATE_SPH == 10)
	  for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++) //add wnorm up, so start from second phase.
	       norm += wnorm [phs_i-1];

      for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
      {
           assert (norm > 0);
    	   phaserho[phs_i-1]= tmprho[phs_i-1]/norm;
    	   density +=phaserho[phs_i-1];
      }
#elif (DENSITY_UPDATE_SPH == 2) || (DENSITY_UPDATE_SPH == 22)
	 for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
	 {
	    if (wnorm[phs_i-1] > 0)
		   phaserho[phs_i-1]= tmprho[phs_i-1] / wnorm[phs_i-1];
	    else
		   phaserho[phs_i-1]=0;

	    density +=phaserho[phs_i-1];
	 }
//#elif DENSITY_UPDATE_SPH == 10
//	  for (phs_i = 2; phs_i<= PHASE_NUM; phs_i++) //add wnorm up, so start from second phase.
//	  {
//	       wnorm[0] += wnorm [phs_i-1];
//	  }
//
//     for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
//     {
//       assert (wnorm[0] > 0);
//   	   phaserho[phs_i-1]= tmprho[phs_i-1] / wnorm[0];
//   	   density +=phaserho[phs_i-1];
//     }
#elif (DENSITY_UPDATE_SPH == 11) || (DENSITY_UPDATE_SPH==13)
     for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
     {
       phaserho[phs_i-1]=tmprho[phs_i-1];
   	   density +=phaserho[phs_i-1];
     }
#endif

#if ADAPTIVE_SML==2
     for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
     {
//       phaserho_ori[phs_i-1]=tmprho_ori[phs_i-1];
   	   density_ori +=tmprho_ori[phs_i-1];
     }
#endif

#ifndef SIMULATE_ASH
     mssfrc=phaserho[1]/density;
#else
     mssfrc=phaserho[1]/density*msfc0_P; //mass fraction of erupted material in the mixture. Phase on is mixture of erupted material with initial mass fraction of erupted material of mssfrac0_P;
#endif

#if DENSITY_UPDATE_SPH == 10 //for this case, density will be updated based on equation according to Suzuki's 2005 paper
     pi->calc_density_suzuki(mssfrc);
#else
     assert(density > 0);
//     assert(!isnan(density));
     pi->put_new_density(density);
#endif //DENSITY_UPDATE_SPH

#if ADAPTIVE_SML==2
     pi->put_dx_density(density_ori);
#endif

#ifdef DEBUG
      if (check_den)
      {
         if (!(density > 0))
        	 cout << "density is invalid" << endl;
      }
#endif

#ifdef DEBUG
      if (check_mssfrac)
      {
         if (!(mssfrc <=1))
        	 cout << "mass fraction larger than 1" << endl;
      }
#endif
      assert(mssfrc <= 1);
      pi->put_new_mass_frac(mssfrc);

//      /*
//       * the bad thing here is: particle do not have information about its bgmesh, the better way is set adapt +=1 only if the bucket where the particle belong to originally does not have involved particles.
//       * Well why not add a new member Key mybucket in particle ---> this is another story, keep track bgmesh for each particle requires additional effort
//       *
//       * Any way, what is the solution? The solution is give up determine the value of adapt(which will determine whether updating of background mesh is necessary!) within density updating function.
//       * update_adpt_domain will be called at more frequent to determine whether updating of background mesh is necessary.
//       */
//
//      if (mssfrc >= MSFRC_THRESH && !pi->is_involved())
//          pi->set_involved_flag(INVOLVED);


    }//end of if --> if the density of that particle need to be updated based summation
  }//end of loop -->go through all particle


  //In the old code, mass fraction is updated while computing, put it separately after density updating done
   THTIterator *it2 = new THTIterator(P_table);
   while ((pi = (Particle *) it2->next()))
     if (pi->need_neigh())
     {
       pi->update_mass_frac();
       mssfrc = pi->get_mass_frac();
       if ( mssfrc>= MSFRC_THRESH && !pi->is_involved())
          pi->set_involved_flag(INVOLVED);

       pi->update_density();
#if (HAVE_ENERGY_SMOOTH==1 || HAVE_ENERGY_SMOOTH==2)
       pi->update_energy();
#endif

       pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P); //The secondary varible was not updated after updating density and mass fraction.
     }

  // clean up stuff
  delete itr;
  delete it2;

  return;
}

