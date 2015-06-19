/*
 * bucket.h
 *
 *  Created on: Mar 5, 2015
 *      Author: zhixuanc
 */

#ifndef BUCKET_H
#define BUCKET_H

#  include <vector>
using namespace std;

#  include <hashtab.h>
#  include <thashtab.h>
#  include <constant.h>
#  include <properties.h>
#  include <pack_data.h>

const int FIRST_BIT_UP = 0x1;
const int SECND_BIT_UP = 0x2;
const int THIRD_BIT_UP = 0x4;
const int FOURTH_BIT_UP = 0x8;

const int EXPT_FIRST_BIT_UP = 0xE;
const int EXPT_SECND_BIT_UP = 0xD;
const int EXPT_THIRD_BIT_UP = 0xB;
const int EXPT_FOURTH_BIT_UP = 0x7;

// Bucket is a unit of background mesh
class Bucket
{
  // friends in repartition
  friend void pack_bucket (BucketPack *, Bucket *, int);
  friend void unpack_bucket (BucketPack *, Bucket *, int);

private:
    Key key;
  double lb_weight;
  double mincrd[DIMENSION];
  double maxcrd[DIMENSION];//domain information
  double bnd[2*DIMENSION];
  double poly[4];
  bool active;
  bool guest_flag;
  bool erupt_flag;/*flag that used to indicate the bucket is source bucket or not
                   * if erupt_flag = true, it is eruption bucket
                   * if erupt_flag = false, it is not eruption bucket
                   * */
//  int newold;
  int myprocess;
  int particles_type; // used to determine whether bucket has ghost or not.
  int bucket_type; //mixed, pressure_bc, underground, overground... 0: invalid bucket
  int bucket_index[2*DIMENSION];//used to determine bucket type.
  int neigh_proc[NEIGH_SIZE];  //neighbor processes
  Key neighbors[NEIGH_SIZE];   //neighbor buckets
  vector < TKey > particles;
  vector < TKey > new_plist;

public:
  //! Contructors
  //! base constructor
    Bucket ();

    // constructors with bucket_index as input without
    // bnd is determined.
    // bnd_xcrd might be useful when ground is not flat
    Bucket (
    		//! Hashtable key of the Bucket
    		unsigned *,
    		//! Minimum coordinates
    		double *,
    		 //! Maximum coordinates
    		double* ,
    		 //! Bucket type ( over-ground, underground or boundary)
    		int,
    		//! boundary elevations in order [i,j; i+1,j; i+1,j+1; i,j+1]
            double *,
            //! my procees id
            int ,
            //! Neighbor process info
            int *,
            //! Array of neighboring buckets HT keys
            Key *,
            //bucket type index
            int*,
            //array bnd
            double*
            );

  //! change process id (only called from repartition)
  void put_myprocess (int myid)
  {
    myprocess = myid;
  }

  //! put repartition weights
  void put_lb_weight (double wght)
  {
    lb_weight = wght;
  }

  //! mark bucket as guest, i.e. belongs to different proc
  void put_guest_flag (int fl)
  {
    guest_flag = fl;
  }

  //! add a particle to the bucket
  void add_particle (TKey pk)
  {
    new_plist.push_back (pk);
  }

  //! put a list of particles
  void put_plist (vector < TKey > pl)
  {
    particles = pl;
  }

  //! put a list of new particles
  void put_new_plist (vector < TKey > pl)
  {
    new_plist.insert (new_plist.end (), pl.begin (), pl.end ());
  }

  //! update particles to new particles
  void update_particles ()
  {
    particles = new_plist;
    new_plist.clear ();
  }

  //! empty particle vector
  void empty_plist ()
  {
    particles.clear ();
    new_plist.clear ();
  }

  //! puts ghosts in underground buckets
  void add_ghosts (vector < TKey > ghosts)
  {
    particles.insert (particles.end (), ghosts.begin (), ghosts.end ());
  }

  //! mark the bucket active
  void mark_active ()
  {
    active = true;
  }

  //! mark bucket inactive
  void mark_inactive ()
  {
    active = false;
  }

  //! mark bucket has eruption
  void mark_erupt ()
  {
    erupt_flag = true;
  }

  //! check if point is contained in bucket or not
  bool contains (double pnt[]) const
  {
    for (int i = 0; i < DIMENSION; i++)
      if ((pnt[i] < mincrd[i]) || (pnt[i] >= maxcrd[i]))
        return false;
    return true;
  }

  //! put have real particles flags
  void set_real_particles (bool flag)
  {
    if (flag)
      particles_type |= FIRST_BIT_UP;
    else
      particles_type &= EXPT_FIRST_BIT_UP;
  }

  //! put have wall ghosts
  void set_wall_ghost_particles (bool flag)
  {
    if (flag)
      particles_type |= SECND_BIT_UP;
    else
      particles_type &= EXPT_SECND_BIT_UP;
  }

  //! put have pressure ghosts
  void set_pressure_ghost_particles (bool flag)
  {
    if (flag)
      particles_type |= THIRD_BIT_UP;
    else
      particles_type &= EXPT_THIRD_BIT_UP;
  }

  //! put have erupt ghosts
  void set_erupt_ghost_particles (bool flag)
  {
    if (flag)
      particles_type |= FOURTH_BIT_UP;
    else
      particles_type &= EXPT_FOURTH_BIT_UP;
  }
//  //! put newold information
//  void put_new_old (int info)
//  {
//    newold = info;
//  }

  // access methods
  //! Access HT key of current bucket
  Key getKey () const
  {
    return key;
  }

  //! get my process id
  int get_myprocess ()
  {
    return myprocess;
  }

  //! Access minimum coordinates of current bucket
  const double *get_mincrd () const
  {
    return mincrd;
  }

  //! Access maximum coordinates of current bucket
  const double *get_maxcrd () const
  {
    return maxcrd;
  }

  //! Access bucket-type
  int get_bucket_type ()
  {
    return bucket_type;
  }

  //! check if bucket belongs to different proc
  int is_guest ()
  {
    return guest_flag;
  }

  //!get particle list in the buckets
  vector <TKey> get_particle_list () const
  {
	  return particles;
  }

//  //! get newold info
//  int get_new_old ()
//  {
//    return newold;
//  }

  //! get repartition weights
  double get_lb_weight () const
  {
    return lb_weight;
  }

  //! Compare hash-table keys for equality
  bool compare_keys (Key k1, Key k2) const
  {
    for (int i = 0; i < KEYLENGTH; i++)
      if (k1.key[i] != k2.key[i])
        return false;
    return true;
  }

  //! Value of elevation z(x,y) using linear interpolation, it will be only called by Mixed Buckets
 double get_bndZ (double coord[]) const
 {
    int i;
    double x[DIMENSION-1];
    for (i = 0; i < DIMENSION-1; i++)
      x[i] = coord[i] - mincrd[i];

    register double pl[4];

    for (i = 0; i < 4; i++)
        pl[i] = poly[i];

    double ground_type;
    for (i=0; i<3; i++) //first 3 cof is zero, means the ground is flat.
    	 ground_type += abs(pl[i]);

    if (ground_type == 0)//ground_type == 0 means that the ground BC is a flat ground.
       return bnd[4];
    else
       return (poly[0] * x[0] + poly[1] * x[1] +
              poly[2] * x[0] * x[1] + poly[3] + mincrd[2]);
 }

  // get bnd
  const double *get_bnd () const
  {
    return bnd;
  }

  const int *get_bucket_index () const
  {
    return bucket_index;
  }

  // update neigh_proc
  void put_neigh_proc (int *, int);

  /*!
   * distance of point from the boundary,
   * also find point of intersection with perpendicular
   */

  //! Check is the bucket has some portion underground or not
  bool is_onground () const
  {
	  return (bucket_index[4] == -1 && bucket_index[5] == 0);
  }

  //! Check is the bucket is active/inactive
  bool is_active () const
  {
    return active;
  }

  //! Check is the eruption source or not
  bool is_erupt () const
  {
    return erupt_flag;
  }

  //! check if bucket has any real particles
  bool has_real_particles () const
  {
    return (particles_type & FIRST_BIT_UP);
  }

  //! check if bucket has wall ghost particles
  bool has_wall_ghost_particles () const
  {
    return (particles_type & SECND_BIT_UP);
  }

  //! check if bucket has pressure ghost particles
  bool has_pressure_ghost_particles () const
  {
    return (particles_type & THIRD_BIT_UP);
  }

  //! check if bucket has erupt ghost particles
  bool has_erupt_ghost_particles () const
  {
    return (particles_type & FOURTH_BIT_UP);
  }

  //! get list of particles in the current bucket
  vector < TKey > get_plist ()
  {
    return particles;
  }

  //!  array of neighboring buckets
  Key *get_neighbors ()
  {
    return neighbors;
  }

  //! get neigh_proc info. neigh_proc also tells if neigh is boundary
  const int *get_neigh_proc () const
  {
    return neigh_proc;
  }

  //! overload get neigh_proc info. neigh_proc also tells if neigh is boundary
  void get_neigh_proc (int * neigh_pc, int *np, int np_total)
  {
    int i, k;
    k=0;
    for (i=0; i<NEIGH_SIZE; i++)
    {
    	if (neigh_proc[i]>0 && neigh_proc[i]<np_total) //make sure that the id of processor is valid
    	{
    		neigh_pc[k]=neigh_proc[i];
    		k++;
    	}
    }

    *np = k;
  }

  //! get HT key of the neighbor is  up,down etc direction
  Key which_neigh (int dir[]) const;

  //! get neigh_proc info in up, down etc direction
  int which_neigh_proc (int dir[]) const;

  //! find out direction of the neighbor
  bool find_neigh_dir (Key k, int dir[]) const;

  /*!
   *  get boundary normal and return 0 if the bucket is a boundary bucket
   *  else return 1
   *  @param pnt : coordinates of any point, ususally a SPH particle
   *  @param normal : cosines of normal to the boundary
   *  @return
   *  0 if current bucket is a boundary bucket, 1 otherwise
   */
  int get_bnd_normal (double * pnt, double * normal) const;

  //!  get distance of point x fom the boundary,
  double get_bnddist (double * pnt, double * intsct) const;

  /*! calculate intersection of line and the ground boundary, such
   * that the line is normal to boundary at pt. of intersection
   */
  int calc_intersection_ground (double * pnt, double * intsct) const;

  /*function that used to determine whether the real particle will escape the domain and lost
   * Only consider Mixed Bucket,
   * The input pos should be coordinates of an real particle
   * */
  bool determine_escape (
		  double * // position of particles
		  );

  //*! Add particles for initial piles.
  void add_real_particle (TKey k)
  {
    active = true;
    set_real_particles (true);
    particles.push_back (k);
  }

  //*! Add ghost particles
  void add_wall_ghost_particle (TKey k)
  {
    active = true;
    set_wall_ghost_particles (true);
    particles.push_back (k);
  }

  //*! Add ghost particles
  void add_pressure_ghost_particle (TKey k)
  {
    active = true;
    set_pressure_ghost_particles (true);
    particles.push_back (k);
  }

  //*! Add ghost particles
  void add_erupt_ghost_particle (TKey k)
  {
    active = true;
    set_erupt_ghost_particles (true);
    particles.push_back (k);
  }

  //! Add ghost particles to the bucket
  int put_ghost_particles (
    //! Particle HashTable
    THashTable *,
    //! Background Mesh
    HashTable *,
    //! Material properties
    MatProps *,
    //! Time properties
    TimeProps*
    );
  //! overloading Add ghost particles to the bucket
  int put_ghost_particles ();

  void Copy_data (
    //! void pointer to datastream
    void *,
    //! HashTable of particles
    HashTable *,
    //! current process id
    int);
};

#endif /* BUCKET_H_ */
