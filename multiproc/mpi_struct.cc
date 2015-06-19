
/*
 * =====================================================================================
 *
 *       Filename:  mpi_struct.cc
 *
 *    Description:  
 *
 *        Created:  01/31/2011 05:05:46 PM EST
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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <mpi.h>
#include <constant.h>
#include <bnd_image.h>

#include "pack_data.h"
#include "repartition_BSFC.h"

MPI_Datatype BUCKET_TYPE;
MPI_Datatype PARTICLE_TYPE;
MPI_Datatype BND_IMAGE_TYPE;
MPI_Datatype LB_VERT_TYPE;

void
GMFG_new_MPI_Datatype ()
{

  int blockcounts[3];
  MPI_Datatype types[3];
  MPI_Aint displs[3];
  int d;
  BucketPack * buck = new BucketPack;

  blockcounts[0] = 5 + NEIGH_SIZE + 2*DIMENSION;
  blockcounts[1] = KEYLENGTH * (1 + NEIGH_SIZE) + MAX_PARTICLES_PER_BUCKET * TKEYLENGTH ;
  blockcounts[2] = (4 * DIMENSION) + 4;
  MPI_Address (&(buck->myprocess), &displs[0]);
  MPI_Address (&(buck->key[0]), &displs[1]);
  MPI_Address (&(buck->mincrd[0]), &displs[2]);

  types[0] = MPI_INT;
  types[1] = MPI_UNSIGNED;
  types[2] = MPI_DOUBLE;

  for (d = 2; d >= 0; d--)// disp is the relative displacement, what MPI_Address return is absolute disp. That is why we need this step here!
    displs[d] -= displs[0];

  MPI_Type_struct (3, blockcounts, displs, types, &BUCKET_TYPE);
  MPI_Type_commit (&BUCKET_TYPE);

  //create the 2nd new d_type
  int blockcounts2[3];
  MPI_Datatype types2[3];
  MPI_Aint displs2[3];

  ParticlePack * particlePack = new ParticlePack;

  // 1 int , 2 unsigned, bunch of doubles
  blockcounts2[0] = 3;
  blockcounts2[1] = TKEYLENGTH;
  blockcounts2[2] = 3 + DIMENSION + NO_OF_EQNS;

  // get adresses
  MPI_Address (&(particlePack->bc_type), &displs2[0]);
  MPI_Address (&(particlePack->key), &displs2[1]);
  MPI_Address (&(particlePack->mass), &displs2[2]);

  types2[0] = MPI_INT;
  types2[1] = MPI_UNSIGNED;
  types2[2] = MPI_DOUBLE;

  for (d = 2; d >= 0; d--)
    displs2[d] -= displs2[0];

  MPI_Type_struct (3, blockcounts2, displs2, types2, &PARTICLE_TYPE);
  MPI_Type_commit (&PARTICLE_TYPE);

  // create 3rd datatype for Boundary Images
  int blockcounts3[3] = { 2, TKEYLENGTH + KEYLENGTH, DIMENSION + NO_OF_EQNS };
  MPI_Datatype type3[3] = { MPI_INT, MPI_UNSIGNED, MPI_DOUBLE };
  MPI_Aint displs3[3];

  // get adresses
  BndImage * bndimage = new BndImage ();

  MPI_Address (&(bndimage->buckproc), &(displs3[0]));
  MPI_Address (&(bndimage->bucket_key), &(displs3[1]));
  MPI_Address (&(bndimage->coord), &(displs3[2]));

  for (d = 2; d >= 0; d--)
    displs3[d] -= displs3[0];

  MPI_Type_struct (3, blockcounts3, displs3, type3, &BND_IMAGE_TYPE);
  MPI_Type_commit (&BND_IMAGE_TYPE);

  //This is not necessary to change!
  // create 4th datatype for SFC verticies
  int blockcounts6[3] = {3,KEYLENGTH+1,1};
  MPI_Datatype types6[3] = {MPI_INT, MPI_UNSIGNED, MPI_FLOAT};
  MPI_Aint displs6[3];

  sfc_vertex_t * sfc_vert_ptr = new sfc_vertex_t;

  MPI_Address(&(sfc_vert_ptr->destination_proc), &displs6[0]);
  MPI_Address(&(sfc_vert_ptr->sfc_key[0]), &displs6[1]);
  MPI_Address(&(sfc_vert_ptr->lb_weight), &displs6[2]);

  for(d=2; d>=0; d--)
    displs6[d]-=displs6[0]; 

  MPI_Type_struct(3, blockcounts6, displs6, types6, &LB_VERT_TYPE);
  MPI_Type_commit(&LB_VERT_TYPE);
  
}//New data types are created at this point
