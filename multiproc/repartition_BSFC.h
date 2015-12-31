#ifndef REPARTATION_BSFC_H
#  define REPARTATION_BSFC_H

/* define some constants that are used in multiple files */
#  define BSFC_NO_CUT 0
#  define BSFC_CUT 1
#  define BSFC_NOT_BALANCED 1
#  define BSFC_BALANCED 0
#  define BSFC_COARSE_LEVEL_FLAG 2

#  include <vector>
using namespace std;

#  include "pack_data.h"
#  include <bucket.h>
#  include <hashtab.h>
#  include <thashtab.h>
#  include <particle.h>

/* BSFC vertex  */
struct sfc_vertex_t
{
  int destination_proc;
  /*!
   * cut_bin_flag = BSFC_NO_CUT = 0 if this object does 
   * not belong to a bin with a cut in it
   * = BSFC_CUT = 1 if this object belongs to a bin with
   * a cut in it 
   */
  int cut_bin_flag;
  /*! 
   * used for creating a linked list of vertices 
   *   (linklist is fortran style) 
   */
  int next_sfc_vert_index;
  /*! 
   * space-filling curve key 
   */
  unsigned sfc_key[KEYLENGTH];
  unsigned obj_key[KEYLENGTH];
  unsigned my_bin;
  float lb_weight;
};

struct unstructured_communication
{
  int send_count;
  int recv_count;
  int used_flag;
  int * send_procs_ptr;
  int * recv_procs_ptr;
  sfc_vertex_t * recv_sfc_vert;
};

// get power of two
inline int
BSFC_pow2(int intexp)
{
  return (1 << intexp);
}

/* declare functions in the sfc routines */
//This is a useless declaration
void BSFC_update_element_proc(int, int, HashTable *, sfc_vertex_t *);

void BSFC_update_and_send_elements(int, int, THashTable *, HashTable *);

int BSFC_find_imbalance
  (float *work_percent,
   float cumulative_work, float total_work, int which_proc, int numprocs);

int BSFC_get_array_location(
                             //! Number of bins
                             int,
                             //! Number of bits
                             int,
                             //! Number of previously used bits 
                             int,
                             //! Array of Vertices
                             sfc_vertex_t *);

void BSFC_refine_partition(
                            //! Array of local imbalanced flags
                            int *,
                            //! Amount of bits used
                            int *,
                            //! Number of vetices in cut 
                            int,
                            //! Array of Vertices
                            sfc_vertex_t *,
                            //! Array of work percentage
                            float *,
                            //! total weight 
                            float,
                            //! Global actual work allocated
                            float *,
                            //! Number of cuts
                            int,
                            //! Linked list of bin heads
                            int *,
                            //! Previously allocated work
                            float *,
                            //! Sub-bins per bin
                            int,
                            //! Array of local balanced flags
                            int *,
                            //! My process ID
                            int myid,
                            //! Total number of processes
                            int numprocs);

void BSFC_create_refinement_info(
                                  //! Number of Cuts in SFC
                                  int *,
                                  //! Global work allocated
                                  float *,
                                  //! Total weight
                                  float,
                                  //! Work percent array
                                  float *,
                                  //! Vertices cut info
                                  unstructured_communication,
                                  //! Array to store priviously allocated work
                                  float **,
                                  //! My process ID
                                  int myid,
                                  //! Total number of processes
                                  int numprocs);

//! Create Bins from Space filling curve
void BSFC_create_bins(
                       //! Number of local buckets
                       int,
                       //! BSFC_Vertex pointer
                       sfc_vertex_t *,
                       //! Amount of bits used 
                       int *,
                       //! size of unsigned
                       int ,
                       //! Global actucal work allocated
                       float *,
                       //! Work percent array
                       float *,
                       //! Total weight pointer
                       float *,
                       //! Balanced flag
                       int *,
                       //! Vertcies cut info 
                       unstructured_communication *,
                       //! Number of cuts
                       int *,
                       //! Bins per process
                       int,
                       //! My process id
                       int,
                       //! Number of processes 
                       int numprocs);

//! Packs bucket data for migration
//for brief bucket
void pack_bucket(
                  //! packet for bucket data
		          BriefBucketPack *,
                  //! pointer to the bucket being packed
				  BriefBucket *,
                  //! destination process id
                  int);
//for bucket
void pack_bucket(
                  //! packet for bucket data
                  BucketPack *,
                  //! pointer to the bucket being packed
                  Bucket *,
                  //! destination process id
                  int);

//! Unpacks bucket data from recieved packets
//for brief bucket
void unpack_bucket(
                    //! pointer to recieved bucket packet
		            BucketPack *,
                    //! pointer to bucket to be filled
                    Bucket *,
                    //! current process id
                    int);
//for bucket
void unpack_bucket(
                    //! pointer to recieved bucket packet
		            BriefBucketPack *,
                    //! pointer to bucket to be filled
					BriefBucket *,
                    //! current process id
                    int);

//! Packs particles in a bucket for migration
void pack_particles(
                     //!  Particle to be packed
                     Particle *,
                     //! array of particle packets
                     ParticlePack *);

//! Unpacks data to paticle class from recieved particle packets
void unpack_particle(
                      //! pointer to recieved particle packet
                      ParticlePack *,
                      //! pointer to particle, to which data will be filled
                      Particle *);

#endif // REPARTATION_BSFC__H
