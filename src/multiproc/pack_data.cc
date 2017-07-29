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

// pack buck ---> overload
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

// pack buck -->Overload
void pack_bucket (BucketPackAdd *buckpack, Bucket *sendbuck, int process)
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

//----------------------------------------------------------------------------------------
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
//----------------------------------------------------------------------------------------

  vector<InfluxAddingPos> add_list = sendbuck->velo_inlet_add;
  int addsize = (int) add_list.size();
  if (addsize > ADDING_NUM)
  {
    cerr << "Number of adding position exceed Maximum allowable limit." << endl;
    cerr << "Error at line " << __LINE__ <<" in file " << __FILE__ << endl;
    exit (1);
  }

  // pack adding position within the bucket
  buckpack->NumAdd = addsize;
  for ( i=0; i < addsize; i++ )
    for ( j=0; j < DIMENSION; j++ )
      buckpack->velo_inlet_add[i*DIMENSION+j] = add_list[i].crd[j];
//--------------------------------------------------------------------------------------------

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
#if DENSITY_UPDATE_SML==0 || ADAPTIVE_SML==2
  pack_array->smlen_original = psend->smlen_original;
#endif

#if ADAPTIVE_SML==2
  pack_array->rho_based_on_dx = psend->rho_based_on_dx;
#endif

#if ADAPTIVE_SML==31
  pack_array->m_ind = psend->m_ind;
#endif

  pack_array->mass_frac = psend->mass_frac;
  pack_array->smoothed_e = psend->smoothed_e;
  pack_array->sound_speed = psend->sound_speed;
  pack_array->myprocess = psend->myprocess;
   
  for (j=0; j < TKEYLENGTH; j++)
    pack_array->key[j] = psend->key.key[j];

  for (j=0; j < DIMENSION; j++)
    pack_array->coords[j] = psend->coord[j];

#if ADAPTIVE_SML==3 || ADAPTIVE_SML==31
  for (j=0; j < DIMENSION; j++)
    pack_array->dm[j] = psend->dm[j];
#endif

  for (j=0; j < DIMENSION; j++)
      pack_array->smoothed_v[j] = psend->smoothed_v[j];

  for (j=0; j < NO_OF_EQNS; j++)
    pack_array->state_vars[j] = psend->state_vars[j];

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  for (j = 0; j < DIMENSION; j++)
	  pack_array->d_rho[j]=psend->d_rho[j];

  for (j = 0; j < DIMENSION; j++)
	  pack_array->d_u[j]=psend->d_u[j];

  for (j = 0; j < DIMENSION; j++)
	  pack_array->d_v[j]=psend->d_v[j];

  for (j = 0; j < DIMENSION; j++)
	  pack_array->d_w[j]=psend->d_w[j];

  for (j = 0; j < DIMENSION; j++)
	  pack_array->d_p[j]=psend->d_p[j];

#endif

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
  buck->erupt_flag = (unsigned) recvdBuck->erupt_flag;
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

//unpack bucket
void unpack_bucket ( BucketPackAdd *recvdBuck, Bucket *buck, int myid)
{

  int i, j;
  vector<TKey> plist;

  buck->myprocess = myid;//Why do not use recvdBuck->myprocess? Need double check!
  buck->bucket_type = recvdBuck->bucket_type;
  buck->particles_type = recvdBuck->particles_type;
  buck->active = (bool) recvdBuck->activeflag;
  buck->erupt_flag = (unsigned) recvdBuck->erupt_flag;
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

  //unpack adding position
  buck->velo_inlet_add.clear();
  int addsize = recvdBuck->NumAdd;
  InfluxAddingPos tmpAdd;
  for ( i=0; i < addsize; i++ )
  {
    for ( j=0; j < DIMENSION; j++ )
    	tmpAdd.crd[j] = recvdBuck->velo_inlet_add[i*DIMENSION+j];
    buck->velo_inlet_add.push_back(tmpAdd);
  }

  //
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
#if DENSITY_UPDATE_SML==0 || ADAPTIVE_SML==2
  part->smlen_original = packet->smlen_original;
#endif
#if ADAPTIVE_SML==2
  part->rho_based_on_dx = packet->rho_based_on_dx;
#endif

#if ADAPTIVE_SML==31
  part->m_ind = packet->m_ind;
#endif
  part->myprocess = packet->myprocess;
  
  for ( i=0; i < DIMENSION; i++ )
    part->coord[i] = packet->coords[i];

#if ADAPTIVE_SML==3 || ADAPTIVE_SML==31
  for ( i=0; i < DIMENSION; i++ )
    part->dm[i] = packet->dm[i];
#endif

  for ( i=0; i < DIMENSION; i++ )
      part->smoothed_v[i] = packet->smoothed_v[i];

  for ( i=0; i < NO_OF_EQNS; i++ )
    part->state_vars[i] = packet->state_vars[i];

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //derivatives
  for ( i=0; i < DIMENSION; i++ )
    part->d_rho[i] = packet->d_rho[i];

  for ( i=0; i < DIMENSION; i++ )
    part->d_u[i] = packet->d_u[i];

  for ( i=0; i < DIMENSION; i++ )
    part->d_v[i] = packet->d_v[i];

  for ( i=0; i < DIMENSION; i++ )
    part->d_w[i] = packet->d_w[i];

  for ( i=0; i < DIMENSION; i++ )
    part->d_p[i] = packet->d_p[i];
#endif

  part->mass_frac = packet->mass_frac;

  part->smoothed_e = packet->smoothed_e;

  part->sound_speed = packet->sound_speed;

  //update second variable every time after communication as secondary variable is not been send and receive while communication.
  part->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P);

  return;
}
