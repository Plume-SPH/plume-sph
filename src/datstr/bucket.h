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
#  include <Involved_header.h>

const int FIRST_BIT_UP = 0x1;
const int SECND_BIT_UP = 0x2;
const int THIRD_BIT_UP = 0x4;
const int FOURTH_BIT_UP = 0x8;

const int EXPT_FIRST_BIT_UP = 0xE;
const int EXPT_SECND_BIT_UP = 0xD;
const int EXPT_THIRD_BIT_UP = 0xB;
const int EXPT_FOURTH_BIT_UP = 0x7;

//brief summary of bucket
//This is used to store empty bucket
//other class will be derived based on this most basic class:
class BriefBucket
{
	  // friends in repartition
	  friend void pack_bucket (BriefBucketPack *, BriefBucket *, int);
	  friend void unpack_bucket (BriefBucketPack *, BriefBucket *, int);
protected:
	  int myprocess;   //the myprocess is the id of process where the bucket originally belong to--> that is to say, for guest buckets, myprocess should be the id of its own home.
	  int neigh_proc[NEIGH_SIZE];  //neighbor processes
	  Key key;
	  double mincrd[DIMENSION];
	  bool is_brief; //Flag that used to indicate the basic type of bucket:
	                 //True: the bucket is brief bucket
	                 //False: the bucket is not brief bucket
	  bool guest_flag; //if true, bucket is guest
	                   //if false, bucket is not guest
public:
	  BriefBucket();
	  BriefBucket(
	    		//!key
	    		unsigned *,
	    		//! Minimum coordinates
	    		double *,
	            //! my procees id
	            int ,
	            //! Neighbor process info
	            int *
	            );
	  //! change process id (only called from repartition)
	  void put_myprocess (int myid)
	  {
	    myprocess = myid;
	  }

	  //! get my process id
	  int get_myprocess ()
	  {
	    return myprocess;
	  }

	  // update neigh_proc
	  void put_neigh_proc (int *, int);

	  //! get neigh_proc info. neigh_proc also tells if neigh is boundary
	  const int *get_neigh_proc () const
	  {
	    return neigh_proc;
	  }

	  //! get neigh_proc info in up, down etc direction
	  int which_neigh_proc (int dir[]) const;

	  //! Access HT key of current bucket
	  Key getKey () const
	  {
	    return key;
	  }

	  //! Compare hash-table keys for equality
	  bool compare_keys (Key k1, Key k2) const
	  {
	    for (int i = 0; i < KEYLENGTH; i++)
	      if (k1.key[i] != k2.key[i])
	        return false;
	    return true;
	  }

	  //! Access minimum coordinates of current bucket
	  const double *get_mincrd () const
	  {
	    return mincrd;
	  }
	  //! mark bucket as guest, i.e. belongs to different proc
	  void put_guest_flag (bool fl)
	  {
	    guest_flag = fl;
	  }

	  //! check if bucket belongs to different proc
	  bool is_guest () const
	  {
	    return guest_flag;
	  }

	  //! figure out whether the bucket is brief or not
	  bool check_brief() const
	  {
		return is_brief;
	  }

	  //The following are virtual member functions
	  //! Access has_involved flag
	  virtual int get_has_involved ()
	  {
	    return HAS_NON_INVOLVED;
	  }

	  //! get list of particles in the current bucket -->return nothing
	  virtual vector < TKey > get_plist ()
	  {
		  return vector<TKey>();
	  }

	  //! Access bucket-type
	  virtual int get_bucket_type ()
	  {
	    return BREIF;
	  }

//	  //! check if bucket belongs to different proc
//	  virtual int is_guest ()
//	  {
//	    return false;
//	  }

	  //! Check is the bucket is active/inactive
	  virtual bool is_active ()
	   {
	     return false;
	   }

	   //! Check is the eruption source or not
	  virtual bool is_erupt ()
	   {
	     return false;
	   }

	   //! check if bucket has any real particles
	  virtual bool has_real_particles ()
	   {
	     return false;
	   }

	   //! check if bucket has wall ghost particles
	  virtual bool has_wall_ghost_particles ()
	   {
	     return false;
	   }

	   //! check if bucket has pressure ghost particles
	  virtual bool has_pressure_ghost_particles ()
	   {
	     return false;
	   }

	   //! check if bucket has erupt ghost particles
	  virtual bool has_erupt_ghost_particles ()
	   {
	     return false;
	   }


	   //! check if bucket has potential involved particle
	  virtual bool is_has_potential_involved ()
	   {
	     return false;
	   }

	   //! check if bucket has potential involved particle
	  virtual bool is_has_involved ()
	   {
	     return false;
	   }

	  //! Check if the bucket is on ground or under ground
	  virtual  bool is_onground_or_underground ()
	  {
		  return false;
	  }

	  virtual bool find_neigh_dir (
			  Key ,  //keyin,
			  Key* , //neighbor_key,
			  int dir[] //direction
				);

};

// Bucket is a unit of background mesh
class Bucket: public BriefBucket
{
  // friends in repartition
  friend void pack_bucket (BucketPack *, Bucket *, int);
  friend void unpack_bucket (BucketPack *, Bucket *, int);

  friend void pack_bucket (BucketPackAdd *, Bucket *, int);
  friend void unpack_bucket (BucketPackAdd *, Bucket *, int);

protected:
  bool active; //active and inactive flag make sense for non brief buckets
               //Active flag should be always false for brief bucket --> not necessary to define it as a member in the class.
  unsigned erupt_flag;/*flag that used to indicate the bucket is inlet source bucket or not
                    * if erupt_flag = 1, it is eruption bucket for plume
                    * if erupt_flag = 2, it is influx bucket for umbrella
                    * if erupt_flag = 0, it is not eruption bucket
                    * */

  int has_involved; /*flag that used to determine what kind of particles does the bucekt contain
                      * 0 : no involved particles
                      * 1 : has potential involved particles
                      * 2 : only has involved particles
                      * 3 : has both involved particles and potential involved particles
                     */
 //  int newold;
   int particles_type; // used to determine whether bucket has ghost or not.
                        /*
                         * 1 bit: real
                         * 2 bit: wall
                         * 3 bit: pressure
                         * 4 bit: erupt/influx
                         */
   int bucket_type; //mixed, pressure_bc, underground, overground... 0: invalid bucket

  double lb_weight;

  int bucket_index[2*DIMENSION];//used to determine bucket type.

  Key neighbors[NEIGH_SIZE];   //neighbor buckets

  double maxcrd[DIMENSION];//domain information
  double poly[4];
  double bnd[2*DIMENSION];

  vector < TKey > particles;
  vector < TKey > new_plist;

  vector <InfluxAddingPos> velo_inlet_add;

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

  //! put repartition weights
  void put_lb_weight (double wght)
  {
    lb_weight = wght;
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
    erupt_flag = 1;
  }

  //! mark bucket influx bucket
  void mark_influx ()
  {
    erupt_flag = 2;
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

  // set has_potential involved
  void set_has_potential_involved (bool flag)
  {
    if (flag)
      has_involved |= FIRST_BIT_UP;
    else
      has_involved &= SECND_BIT_UP;
  }
  // set has_involved
  void set_has_involved (bool flag)
  {
    if (flag)
        has_involved |= SECND_BIT_UP;
    else
    	has_involved &= FIRST_BIT_UP;
  }

  //void set has no involved
  void set_has_no_involved ()
  {
	  has_involved = 0;
  }

  //put particles_type
  void put_particles_type (int in)
   {
	  particles_type = in;
   }

  // access methods

  //! Access maximum coordinates of current bucket
  const double *get_maxcrd () const
  {
    return maxcrd;
  }

  //! Access bucket-type
  int get_bucket_type () const
  {
    return bucket_type;
  }

  //! Access has_involved flag
  int get_has_involved () const
  {
    return has_involved;
  }

  //! get is erupt flag
  unsigned get_erupt_flag () const
  {
	  return erupt_flag;
  }
  //!get particle list in the buckets
  vector <TKey> get_particle_list () const
  {
	  return particles;
  }

  //!get inlet velocity bc particle adding list
  vector <InfluxAddingPos> get_adding_list () const
  {
	  return velo_inlet_add;
  }

  //! check whethter the bucket has adding position
  bool has_adding_pos ()
  {
	  return velo_inlet_add.size();
  }

  //! get repartition weights
  double get_lb_weight () const
  {
    return lb_weight;
  }

  //! get particles type:
  int get_particles_type() const
  {
	  return particles_type;
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

  /*!
   * distance of point from the boundary,
   * also find point of intersection with perpendicular
   */

  //! Check is the bucket has some portion underground or not
  bool is_onground () const
  {
	  return (bucket_index[4] == -1 && bucket_index[5] == 0);
  }

  //! Check if the bucket is on ground or under ground
  bool is_onground_or_underground () const
  {
	  return (bucket_index[4] == -1);
  }

  //! Check is the bucket is active/inactive
  bool is_active () const
  {
    return active;
  }

  //! Check is the eruption source or not
  bool is_erupt () const
  {
    return erupt_flag==1;
  }

  //! Check is the influx source or not
  bool is_influx () const
  {
    return erupt_flag==2;
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


  //! check if bucket has potential involved particle
  bool is_has_potential_involved () const
  {
    return (has_involved & FIRST_BIT_UP);
  }

  //! check if bucket has potential involved particle
  bool is_has_involved () const
  {
    return (has_involved & SECND_BIT_UP);
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

  //! get HT key of the neighbor is  up,down etc direction
  Key which_neigh (int dir[]) const;

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

  //function that used to determine the type of bucket
  //--->This is actually exactly the same as the bucket type determine function from preprocess
  //--->Remember to update this section when any modification is made in preprocess
  void determine_bucket_type (
		  double *,
          double *);

  //! add a particle to the bucket
  void add_particle (TKey pk)
  {
    new_plist.push_back (pk);
  }

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

  //*! Add velocity inlet bc particle adding position
  void add_adding_pos (InfluxAddingPos add)
  {
	  velo_inlet_add.push_back (add);
  }

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
