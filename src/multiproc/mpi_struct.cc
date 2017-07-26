
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
#include <Involved_header.h>

#include "pack_data.h"
#include "repartition_BSFC.h"

MPI_Datatype BRIEF_BUCKET_TYPE;
MPI_Datatype BUCKET_ADD_TYPE;
MPI_Datatype BUCKET_TYPE;
MPI_Datatype PARTICLE_TYPE;
MPI_Datatype BND_IMAGE_TYPE;
MPI_Datatype LB_VERT_TYPE;
MPI_Datatype INVOLVED_HEADER_TYPE;

void
GMFG_new_MPI_Datatype ()
{

//MPI data structure for brief bucket
  int blockcounts0[3];
  MPI_Datatype types0[3];
  MPI_Aint displs0[3];
  int d;
  BriefBucketPack * brief_buck = new BriefBucketPack;

  blockcounts0[0] = 2+NEIGH_SIZE;
  blockcounts0[1] = KEYLENGTH;
  blockcounts0[2] = DIMENSION;
  MPI_Address (&(brief_buck->myprocess), &displs0[0]);
  MPI_Address (&(brief_buck->key[0]), &displs0[1]);
  MPI_Address (&(brief_buck->mincrd[0]), &displs0[2]);

  types0[0] = MPI_INT;
  types0[1] = MPI_UNSIGNED;
  types0[2] = MPI_DOUBLE;

  for (d = 2; d >= 0; d--)// disp is the relative displacement, what MPI_Address return is absolute disp. That is why we need this step here!
    displs0[d] -= displs0[0];

  MPI_Type_struct (3, blockcounts0, displs0, types0, &BRIEF_BUCKET_TYPE);
  MPI_Type_commit (&BRIEF_BUCKET_TYPE);

  //MPI data strucutre for bucketAdd
    int blockcounts1[3];
    MPI_Datatype types1[3];
    MPI_Aint displs1[3];
    BucketPackAdd * buck_add = new BucketPackAdd;

    blockcounts1[0] = 8 + NEIGH_SIZE + 2*DIMENSION;
    blockcounts1[1] = KEYLENGTH * (1 + NEIGH_SIZE) + MAX_PARTICLES_PER_BUCKET * TKEYLENGTH ;
    blockcounts1[2] = (4 * DIMENSION) + 4 + ADDING_NUM*DIMENSION;
    MPI_Address (&(buck_add->myprocess), &displs1[0]);
    MPI_Address (&(buck_add->key[0]), &displs1[1]);
    MPI_Address (&(buck_add->mincrd[0]), &displs1[2]);

    types1[0] = MPI_INT;
    types1[1] = MPI_UNSIGNED;
    types1[2] = MPI_DOUBLE;

    for (d = 2; d >= 0; d--)// disp is the relative displacement, what MPI_Address return is absolute disp. That is why we need this step here!
      displs1[d] -= displs1[0];

    MPI_Type_struct (3, blockcounts1, displs1, types1, &BUCKET_ADD_TYPE);
    MPI_Type_commit (&BUCKET_ADD_TYPE);

   //MPI data structure for bucket
    int blockcounts[3];
    MPI_Datatype types[3];
    MPI_Aint displs[3];
    BucketPack * buck = new BucketPack;

    blockcounts[0] = 7 + NEIGH_SIZE + 2*DIMENSION;
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


  //MPI data structure for particle
  int blockcounts2[3];
  MPI_Datatype types2[3];
  MPI_Aint displs2[3];

  ParticlePack * particlePack = new ParticlePack;
#if USE_GSPH==0
  blockcounts2[0] = 5;
  blockcounts2[1] = TKEYLENGTH;
  blockcounts2[2] = 5 + 2*DIMENSION + NO_OF_EQNS;
#elif (USE_GSPH==1 || USE_GSPH==2)
  blockcounts2[0] = 5;
  blockcounts2[1] = TKEYLENGTH;
  blockcounts2[2] = 5 + 7*DIMENSION + NO_OF_EQNS;  // For GSPH there will be 5 additional variable for gradient
#endif

#if DENSITY_UPDATE_SML==0 || ADAPTIVE_SML==2  //This is for original smoothing length
  blockcounts2[2]++;
#endif

#if ADAPTIVE_SML==2
  blockcounts2[2]++;
#endif

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

  // create datatype for Boundary Images
#if BC_FOR_KX==1
  int blockcounts3[3] = { 2, TKEYLENGTH + KEYLENGTH, DIMENSION + NO_OF_EQNS +1};  //if impose natural BC for KX, need additional variables in image
#else
  int blockcounts3[3] = { 2, TKEYLENGTH + KEYLENGTH, DIMENSION + NO_OF_EQNS };
#endif
  MPI_Datatype type3[3] = { MPI_INT, MPI_UNSIGNED, MPI_DOUBLE };
  MPI_Aint displs3[3];

  BndImage * bndimage = new BndImage ();

  MPI_Address (&(bndimage->buckproc), &(displs3[0]));
  MPI_Address (&(bndimage->bucket_key), &(displs3[1]));
  MPI_Address (&(bndimage->coord), &(displs3[2]));

  for (d = 2; d >= 0; d--)
    displs3[d] -= displs3[0];

  MPI_Type_struct (3, blockcounts3, displs3, type3, &BND_IMAGE_TYPE);
  MPI_Type_commit (&BND_IMAGE_TYPE);

  // create 6th datatype for SFC verticies
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
  
  // create 7th datatype for Boundary Images
  int blockcounts5[2] = { 1, KEYLENGTH };
  MPI_Datatype type5[2] = { MPI_INT, MPI_UNSIGNED};
  MPI_Aint displs5[2];

  InvolvedHead * involvedhead = new InvolvedHead();

  MPI_Address (&(involvedhead->has_involved), &(displs5[0]));
  MPI_Address (&(involvedhead->bucket_key), &(displs5[1]));

  for (d = 1; d >= 0; d--)
    displs5[d] -= displs5[0];

  MPI_Type_struct (2, blockcounts5, displs5, type5, &INVOLVED_HEADER_TYPE);
  MPI_Type_commit (&INVOLVED_HEADER_TYPE);

}//New data types are created at this point
