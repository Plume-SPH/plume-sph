/*
 * =====================================================================================
 *
 *       Filename:  pack_data.cc
 *
 *    Description:  
 *
 *        Created:  02/03/2011 11:55:04 AM EST
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
#include <iostream>
#include <vector>
using namespace std;

#include <constant.h>
#include <parameters.h>
#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <particle.h>
#include "pack_data.h"

// pack brief buck
void pack_bucket (BriefBucketPack *buckpack, BriefBucket *sendbuck, int process)
{
  int i;
  buckpack->myprocess = process;
  buckpack->is_brief = (int) sendbuck->is_brief;
  for (i=0; i<KEYLENGTH; i++)
    buckpack->key[i] = sendbuck->key.key[i];

  for (i = 0; i < NEIGH_SIZE; i++)
    buckpack->neigh_proc[i] = sendbuck->neigh_proc[i];

  for (i=0; i < DIMENSION; i++)
    buckpack->mincrd[i] = sendbuck->mincrd[i];

  return;
}

// pack buck
void pack_bucket (BucketPack *buckpack, Bucket *sendbuck, int process)
{
  int j, i;
  buckpack->myprocess = process;
  buckpack->bucket_type = sendbuck->bucket_type;
  buckpack->particles_type = sendbuck->particles_type;
  buckpack->erupt_flag = (int) sendbuck->erupt_flag;
  buckpack->activeflag = (int) sendbuck->active;
  buckpack->has_involved = sendbuck->has_involved;
  for (j=0; j<KEYLENGTH; j++)
    buckpack->key[j] = sendbuck->key.key[j];

  for (i = 0; i < NEIGH_SIZE; i++)
  {
    buckpack->neigh_proc[i] = sendbuck->neigh_proc[i];
    for ( j=0; j<KEYLENGTH; j++)
      buckpack->neighs[i*KEYLENGTH+j] = sendbuck->neighbors[i].key[j];
  }

  vector<TKey> particles = sendbuck->particles;
  int psize = (int) particles.size();
  if ( psize > MAX_PARTICLES_PER_BUCKET)
  {
    cerr << "Number of particles exceed Maximum allowable limit." << endl;
    cerr << "Error at line " << __LINE__ <<" in file " << __FILE__ << endl;
    exit (1);
  }

  // pack particle keys within the bucket
  buckpack->NumParticles = psize;
  for ( i=0; i < psize; i++ )
    for ( j=0; j < TKEYLENGTH; j++ )
      buckpack->particles[i*TKEYLENGTH+j] = particles[i].key[j];

  // bucket upper and lower limits
  for (i=0; i < DIMENSION; i++)
  {
    buckpack->mincrd[i] = sendbuck->mincrd[i];
    buckpack->maxcrd[i] = sendbuck->maxcrd[i];
  }

  // boundary fucntion
  for (i=0; i<4; i++)
    buckpack->poly[i] = sendbuck->poly[i];

  for (i=0; i < 2*DIMENSION; i++)
  {
    buckpack->bucket_index[i] = sendbuck->bucket_index[i];
    buckpack->bnd[i] = sendbuck->bnd[i];
  }

  return;
}


void pack_particles (Particle *psend, ParticlePack *pack_array)
{
  int j;
  pack_array->update_delayed = (int) psend->update_delayed;
  pack_array->involved = psend->involved;
  pack_array->bc_type  = psend->bc_type;
  pack_array->phase_num  = psend->phase_num;

  pack_array->mass  = psend->mass;
  pack_array->smlen = psend->smlen;
  pack_array->mass_frac = psend->mass_frac;
  pack_array->myprocess = psend->myprocess;
   
  for (j=0; j < TKEYLENGTH; j++)
    pack_array->key[j] = psend->key.key[j];

  for (j=0; j < DIMENSION; j++)
    pack_array->coords[j] = psend->coord[j];

  for (j=0; j < DIMENSION; j++)
      pack_array->smoothed_v[j] = psend->smoothed_v[j];

  for (j=0; j < NO_OF_EQNS; j++)
    pack_array->state_vars[j] = psend->state_vars[j];

  return;
}

//unpack brief bucket
void unpack_bucket (BriefBucketPack *recvdBuck, BriefBucket *buck, int myid)
{
  int i;

  buck->myprocess = myid;//Why do not use recvdBuck->myprocess? Need double check!
  buck->is_brief = (bool) recvdBuck->is_brief;

  for ( i=0; i < KEYLENGTH; i++ )
    buck->key.key[i] = recvdBuck->key[i];

  for ( i=0; i < NEIGH_SIZE; i++ )
    buck->neigh_proc[i] = recvdBuck->neigh_proc[i];

  for ( i=0; i < DIMENSION; i++ )
    buck->mincrd[i] = recvdBuck->mincrd[i];

  return;
}

//unpack bucket
void unpack_bucket (BucketPack *recvdBuck, Bucket *buck, int myid)
{

  int i, j;
  vector<TKey> plist;

  buck->myprocess = myid;//Why do not use recvdBuck->myprocess? Need double check!
  buck->bucket_type = recvdBuck->bucket_type;
  buck->particles_type = recvdBuck->particles_type;
  buck->active = (bool) recvdBuck->activeflag;
  buck->erupt_flag = (bool) recvdBuck->erupt_flag;
  buck->has_involved = recvdBuck->has_involved;

  for ( i=0; i < KEYLENGTH; i++ )
    buck->key.key[i] = recvdBuck->key[i];

  for ( i=0; i < NEIGH_SIZE; i++ )
  {
    buck->neigh_proc[i] = recvdBuck->neigh_proc[i];
    for ( j=0; j < KEYLENGTH; j++ )
      buck->neighbors[i].key[j]  = recvdBuck->neighs[i*KEYLENGTH+j];
  }

  for ( i=0; i < DIMENSION; i++ )
  {
    buck->mincrd[i] = recvdBuck->mincrd[i];
    buck->maxcrd[i] = recvdBuck->maxcrd[i];
  }

  for ( i=0; i < 4; i++ )
  {
    buck->poly[i] = recvdBuck->poly[i];
  }

  // unpack particle keys
  buck->particles.clear();
  int psize = recvdBuck->NumParticles;
  TKey tmpkey;
  for ( i=0; i < psize; i++ )
  {
    for ( j=0; j < TKEYLENGTH; j++ )
      tmpkey.key[j] = recvdBuck->particles[i*TKEYLENGTH+j];
    buck->particles.push_back(tmpkey);
  }

  for (i=0; i < 2*DIMENSION; i++)
  {
	  buck->bucket_index[i] = recvdBuck->bucket_index[i];
	  buck->bnd[i] = recvdBuck->bnd[i];
  }

  return;
}

void unpack_particle (ParticlePack *packet, Particle *part)
{
  int i;
  part->update_delayed = (bool) packet->update_delayed;
  part->involved = packet->involved;
  part->bc_type = packet->bc_type;
  part->phase_num = packet->phase_num;

  for ( i=0; i < TKEYLENGTH; i++ )
    part->key.key[i] = packet->key[i];
  
  part->mass  = packet->mass;
  part->smlen = packet->smlen;
  part->myprocess = packet->myprocess;
  
  for ( i=0; i < DIMENSION; i++ )
    part->coord[i] = packet->coords[i];

  for ( i=0; i < DIMENSION; i++ )
      part->smoothed_v[i] = packet->smoothed_v[i];

  for ( i=0; i < NO_OF_EQNS; i++ )
    part->state_vars[i] = packet->state_vars[i];

  part->mass_frac = packet->mass_frac;

  //update second variable every time after communication as secondary variable is not been send and receive while communication.
  part->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P);

  return;
}
