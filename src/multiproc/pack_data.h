
/*
 * =====================================================================================
 *
 *       Filename:  pack_data.h
 *
 *    Description:  
 *
 *        Created:  01/31/2011 05:28:38 PM EST
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *        License:  GNU General Public License
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */

#ifndef PACK_DATA_H
#  define PACK_DATA_H

#  include <constant.h>

#  include <Involved_header.h>

struct ParticlePack
{
  // Integers
  int bc_type; //100: real, 0:erupt bc, 1: pressure_bc, 2: wall bc
  int phase_num;
  int update_delayed;
  int involved;
  int myprocess; //originally, myprocess is not updated correctly add this into ParticlePack to make sure that process id for guest particles is correct.

  //unsigned
  unsigned key[TKEYLENGTH];

  //doubles
  double mass;
  double smlen;
#if DENSITY_UPDATE_SML==0 || ADAPTIVE_SML==2
  double smlen_original;
#endif
#if ADAPTIVE_SML==2
  double rho_based_on_dx;
#endif
  double mass_frac;
  double smoothed_e;
#if ADAPTIVE_SML==31
  double m_ind;
#endif
  double sound_speed; //actually sound speed is not primary variable, but is needed when computing pressure (When the fluid is weakly compressible flow)
  double coords[DIMENSION];
#if ADAPTIVE_SML==3 || ADAPTIVE_SML==31
  double dm[DIMENSION];
#endif
  double smoothed_v[DIMENSION];
  double state_vars[NO_OF_EQNS];// rho, v , e

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //derivatives
  double d_rho[DIMENSION];
  double d_u[DIMENSION];
  double d_v[DIMENSION];
  double d_w[DIMENSION];
  double d_p[DIMENSION];
#endif

};

struct BriefBucketPack
{
	int myprocess;
	int is_brief;
	int neigh_proc[NEIGH_SIZE];
	unsigned key[KEYLENGTH];
	double mincrd[DIMENSION];
};

//It is not a good idea to derive BucketPack from BriefBucketPack
//Because that will change the order of members --> as results, corresponding MPI data structure also need to be changed
//It will definitely cause trouble in the future
struct BucketPack
{
  // Integers
  int myprocess;
  int activeflag; //1 means active, 0 means not active
  int erupt_flag;/*flag that used to indicate the bucket is source bucket or not
                    * if erupt_flag = 1, it is eruption bucket for plume
                    * if erupt_flag = 2, it is influx bucket for umbrella
                    * if erupt_flag = 0, it is not eruption bucket
                     * */
  int bucket_type;//0, 1, 2, 3: Mixed, ect
  int particles_type; //have real particle or not? have ghost particles or not...
  int has_involved;//1 means have involved particles, 0 means does not
  int NumParticles;
  int neigh_proc[NEIGH_SIZE];
  int bucket_index[2*DIMENSION];

  // Unsigned integers
  unsigned key[KEYLENGTH];
  unsigned neighs[NEIGH_SIZE * KEYLENGTH];
  unsigned particles[MAX_PARTICLES_PER_BUCKET * TKEYLENGTH];

  // Doubles
  double mincrd[DIMENSION];
  double maxcrd[DIMENSION];
  double poly[4];
  double bnd[2*DIMENSION];

};
//It is not a good idea to derive BucketPackAdd from BucketPack
//Because that will change the order of members --> as results, corresponding MPI data structure also need to be changed
//It will definitely cause trouble in the future
struct BucketPackAdd
{
  // Integers
  int myprocess;
  int activeflag; //1 means active, 0 means not active
  int erupt_flag;/*flag that used to indicate the bucket is source bucket or not
					* if erupt_flag = 1, it is eruption bucket for plume
					* if erupt_flag = 2, it is influx bucket for umbrella
					* if erupt_flag = 0, it is not eruption bucket
					 * */
  int bucket_type;//0, 1, 2, 3: Mixed, ect
  int particles_type; //have real particle or not? have ghost particles or not...
  int has_involved;//1 means have involved particles, 0 means does not
  int NumParticles;
  int NumAdd;
  int neigh_proc[NEIGH_SIZE];
  int bucket_index[2*DIMENSION];

  // Unsigned integers
  unsigned key[KEYLENGTH];
  unsigned neighs[NEIGH_SIZE * KEYLENGTH];
  unsigned particles[MAX_PARTICLES_PER_BUCKET * TKEYLENGTH];

  // Doubles
  double mincrd[DIMENSION];
  double maxcrd[DIMENSION];
  double poly[4];
  double bnd[2*DIMENSION];
  double velo_inlet_add[ADDING_NUM*DIMENSION];
};
#endif // PACK_DATA__H
