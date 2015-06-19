
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

struct ParticlePack
{
  // Integers
  int bc_type; //100: real, 0:erupt bc, 1: pressure_bc, 2: wall bc
  int phase_num;
  int update_delayed;

  //unsigned
  unsigned key[TKEYLENGTH];

  //doubles
  double mass;
  double smlen;
  double coords[DIMENSION];
  double state_vars[NO_OF_EQNS];// rho, v , e
  double mass_frac;
};

struct BucketPack
{
  // Integers
  int myprocess;
  int activeflag; //1 means active, 0 means not active
  int bucket_type;//0, 1, 2, 3: Mixed, ect
  int particles_type; //have real particle or not? have ghost particles or not...
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

#endif // PACK_DATA__H
