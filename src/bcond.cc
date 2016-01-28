/*
 * bcond.cc
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */

#include <vector>
#include <list>
#include <cassert>
#include <iostream>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <particler.h>
#include <bnd_image.h>

#include "constant.h"
#include "particler.h"
#include "sph_header.h"

int
apply_bcond(int myid, THashTable * P_table, HashTable * BG_mesh,
            MatProps * matprops, list < BndImage > & Image_table)
{

  int i, j;
  double refc[DIMENSION];
  double normal[DIMENSION], velrot[DIMENSION];
  double uvec[NO_OF_EQNS], state_vars[NO_OF_EQNS];
  double dx[DIMENSION], s[DIMENSION], pcoord[DIMENSION];
  double supp, bnddist, wnorm;
//  double engr;

  Key * neighbors;//bucket key
  Particle * p_ghost = NULL;
  vector < TKey > plist;
  vector < TKey >::iterator p_itr;
  list < BndImage >::iterator i_img;

  for (i_img = Image_table.begin(); i_img != Image_table.end(); i_img++)
    if (i_img->buckproc == myid)
    {
      // reflection coordinates
      for (i = 0; i < DIMENSION; i++)
        refc[i] = i_img->coord[i];

      // grab to particle from image-reflection
      if (i_img->partproc == myid )
      {
        p_ghost = (Particle *) P_table->lookup(i_img->ghost_key);
        if ( ! p_ghost )
        {
          fprintf (stderr, "Error: ghost particle missing on proc %d\n", myid);
          return 3;
        }

        double d = 0.;
        for (i = 0; i < DIMENSION; i++)
        {
          pcoord[i]= *(p_ghost->get_coords() + i);
          normal[i] = refc[i] - pcoord[i];
          d += normal[i] * normal[i];
        }

        for (i = 0; i < DIMENSION; i++)
          normal[i] /= sqrt(d);
      }

      // reset variables
      bnddist = 1.0E10;
      wnorm = 0.;
      for (i = 0; i < NO_OF_EQNS; i++)
        uvec[i] = 0.;

      // get hold of bucket containing the image --> it should be a no brief bucket, so no worries about whether the bucket is brief bucket or not!
      Bucket *buck = (Bucket *) BG_mesh->lookup(i_img->bucket_key);
      assert(buck);

      // go through every particle in the bucket to calculate its
      // contribution at reflection
      plist = buck->get_plist();
      for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
      {
        Particle *pj = (Particle *) P_table->lookup (*p_itr);

        if (pj->contr_image())
        {
          double h = pj->get_smlen();

          supp = 3 * h;
          for (j = 0; j < DIMENSION; j++)
            dx[j] = refc[j] - *(pj->get_coords() + j);

          // if Particle(j) is in 3-h of reflection ...
          if (in_support(dx, supp))
          {
            for (j = 0; j < DIMENSION; j++)
              s[j] = dx[j] / h;
            for (j = 0; j < NO_OF_EQNS; j++)
              state_vars[j] = *(pj->get_state_vars() + j);
            double w = weight(s, h);
            double mj = pj->get_mass();
            uvec[0] += mj * w; //rho
            uvec[1] += mj * w * state_vars[1] / state_vars[0];//velocity u
            uvec[2] += mj * w * state_vars[2] / state_vars[0];//velocity v
            uvec[3] += mj * w * state_vars[3] / state_vars[0];//velocity w
            uvec[4] += mj * w * state_vars[4] / state_vars[0];//internal energy
            wnorm += mj * w / state_vars[0];
          }
        }
      }//end of loop go through all particles of local bucket

      // now search neighboring buckets ...
      neighbors = buck->get_neighbors();
      for (i = 0; i < NEIGH_SIZE; i++)
        if (*(buck->get_neigh_proc() + i) > -1)
        {
          Bucket *buck_neigh = (Bucket *) BG_mesh->lookup(neighbors[i]);

          // if the neighbor is not found, and it belongs to different
          // to different proc, skip the iteration
          if ((!buck_neigh) && (*(buck->get_neigh_proc() + i) != myid))
            continue;

          // otherwise the neighbor must exist, or it is a bug
          else
            assert(buck_neigh);

          if (buck_neigh->check_brief())
        	  continue;
          else
          {
              // search buckets for real particles in 3h neighborhood
              if (buck_neigh->has_real_particles()  || buck_neigh->has_pressure_ghost_particles())
              {
                plist = buck_neigh->get_plist();
                for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
                {
                  Particle *pj = (Particle *) P_table->lookup(*p_itr);
                  assert(pj);

                  if (pj->contr_image())
                  {
                    double h = pj->get_smlen();

                    supp = 3 * h;
                    for (j = 0; j < DIMENSION; j++)
                      dx[j] = refc[j] - *(pj->get_coords() + j);

                    // if Particle(j) is in 3-h of reflection ...
                    if (in_support(dx, supp))
                    {
                      for (j = 0; j < DIMENSION; j++)
                        s[j] = dx[j] / h;
                      for (j = 0; j < NO_OF_EQNS; j++)
                        state_vars[j] = *(pj->get_state_vars() + j);
                      double w = weight(s, h);
                      double mj = pj->get_mass();

                      uvec[0] += mj * w;
                      uvec[1] += mj * w * state_vars[1] / state_vars[0];
                      uvec[2] += mj * w * state_vars[2] / state_vars[0];
                      uvec[3] += mj * w * state_vars[3] / state_vars[0];
                      uvec[4] += mj * w * state_vars[4] / state_vars[0];//internal energy
                      wnorm += mj * w / state_vars[0];
                    }
                  }
                }//go through all particles in that bucket
              }//if gas real particles or has pressure ghost particles
          }//end of if neigh bucket is not brief
        }// end of search all neighbor buckets

      // renormalize
      if ( wnorm > 0. )
        for (i = 0; i < NO_OF_EQNS; i++)
          uvec[i] /= wnorm;
      else //wnorm might equal to zero only when there is no other particles within kernel support.
      {
        uvec[0] = 1.;
        for (i = 1; i < NO_OF_EQNS; i++)
          uvec[i] = 0.;
      }

//      /*In old code, the following part is inside if statement, actually, they should be at outside of if*/
      reflect (&uvec[1], velrot, normal);//change the direction of image velocity according to norm of wall.
      for (i = 0; i < DIMENSION; i++)
           uvec[i+1] = velrot[i];

      /*In old code, internal energy is not given at the wall, internal energy given---> to force a essentially boundary condition
       *But actually, the way that I am adopting here to imposing essential boundary is not proper---> to imposing proper essential boundary condition requires solving of system of equations implicitly.
       */
//      engr=air_engr_hydro(pcoord);
//      uvec[NO_OF_EQNS-1]=engr;
      // if ghost belongs to my proc, update it
      if (i_img->partproc == myid)
      {
    	  /*In old code, the following part is inside if, actually, they should be at outside of if*/
//        reflect (&uvec[1], velrot, normal);//change the direction of image velocity according to norm of wall.
//        for (i = 0; i < DIMENSION; i++)
//          uvec[i+1] = velrot[i];
        p_ghost->put_state_vars(uvec);
        p_ghost->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P); //This is added later, there was a bug in old code as the secondary variable is not updated after imposing of boundary condition.
        p_ghost->put_update_delayed(false);
      }
      // else store the values to snyc at appropriate time--> The late snyc will be done in move_bnd_img.cc
      else
      {
        for (i = 0; i < NO_OF_EQNS; i++)
          i_img->state_vars[i] = uvec[i];
      }
    }//end of go through all wall ghost particles
  return 0;
}
