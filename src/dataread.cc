/*
 * dataread.cc
 *
 *  Created on: Mar 5, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#define deg2rad(A)  ((A)*(0.01745329252))

#include <iostream>
#include <string>
#include <limits>
#include <fstream>

using namespace std;

#include <hdf5.h>
#include "hdf5calls.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

#ifdef MULTI_PROC
#  include <mpi.h>
#endif

#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <buckstr.h>
#include <buckhead.h>
#include <particler.h>
#include <properties.h>

#include "sph_header.h"
#include "parameters.h"
#include "constant.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif


//function that used to determine the involved flag
//can also be used for determine initial active flag
int determine_bucket_involved_flag (double *mincrd, double *maxcrd, double *left, double *right)
{
	int flag[DIMENSION];
	int bt[2*DIMENSION];
	int k,sum;

    bt[0]=determine_face_type(left[0],maxcrd[0],mincrd[0]);
    bt[1]=determine_face_type(right[0],maxcrd[0],mincrd[0]);
    bt[2]=determine_face_type(left[1],maxcrd[1],mincrd[1]);
    bt[3]=determine_face_type(right[1],maxcrd[1],mincrd[1]);
    bt[4]=determine_face_type(left[2],maxcrd[2],mincrd[2]);
    bt[5]=determine_face_type(right[2],maxcrd[2],mincrd[2]);

    for (k=0; k<DIMENSION; k++)
    	flag[k]=abs(bt[2*k]+bt[2*k+1]);

    for (k=0; k<DIMENSION; k++)
    	sum += flag[k];

    if ( (flag[0]==2) || (flag[1]==2) || (flag[2]==2))
        return 0; //has no involved
    else
    	return 1; //has potential involved
}

int
Read_Data (MatProps * matprops, TimeProps * timeprops, SimProps * simprops, int *format)
{

  ifstream inD1 ("scale.data", ios::in);
  if (inD1.good())
  {
    inD1 >> matprops->LENGTH_SCALE;
    inD1 >> matprops->GRAVITY_SCALE;
  }
  else
  {
    matprops->LENGTH_SCALE = 1.;
    matprops->GRAVITY_SCALE = 1.; //not sure why I need scale? It seems that it is not necessary at all!
  }

  ifstream inD2 ("simulation.data", ios::in);
  if (inD2.fail())
  {
    cerr << "ERROR: Can't find \"simulation.data\" input file." << endl;
    exit(1);
  }

  double temp;
  double len_scale = matprops->LENGTH_SCALE;

  double time_scale = 1.;
  // simulation time properties
  inD2 >> timeprops->max_time;
  inD2 >> timeprops->max_steps;
  inD2 >> timeprops->timeoutput;
  inD2 >> timeprops->stat_erupt;
  inD2 >> timeprops->end_erupt;
  inD2 >> *(format);
  timeprops->TIME_SCALE = time_scale;
  timeprops->ndtimeoutput = timeprops->timeoutput / time_scale;
  timeprops->ndmax_time = timeprops->max_time / time_scale;

  // material properties
  inD2 >> matprops->P_CONSTANT;
  inD2 >> matprops->GAMMA; //Gamma should be changeable in my problem.
  inD2 >> temp;
  matprops->smoothing_length = temp / len_scale;
  matprops->particle_mass = rhoa0_P*pow(matprops->smoothing_length, DIMENSION);//This is incorrect

  //simulation properties
  inD2 >>simprops->Idom_x_min;
  inD2 >>simprops->Idom_x_max;
  inD2 >>simprops->Idom_y_min;
  inD2 >>simprops->Idom_y_max;
  inD2 >>simprops->Idom_z_min;
  inD2 >>simprops->Idom_z_max;


  inD2.close();
 
  return 0;
}

int
Read_Grid (THashTable ** P_table, HashTable ** BG_mesh,
          vector < BucketHead > & partition_tab, MatProps * matprops,
          SimProps* simprops, int myid, int numprocs, int * my_comm)
{
  int No_of_Buckets;
  int BG_TABLE_SIZE = 400000;
  int P_TABLE_SIZE = 4000000;
  double mindom[DIMENSION], maxdom[DIMENSION];
  double mindom_i[DIMENSION], maxdom_i[DIMENSION];
  double mindom_new[DIMENSION], maxdom_new[DIMENSION];
  unsigned btkey[KEYLENGTH];
  Key tempbtkey;
  double hvars[6], min_crd[DIMENSION], max_crd[DIMENSION];
  char filename[14];

  int Down[DIMENSION] = { 0, 0, 1 };

  BucketStruct *bucket;
  Key neigh_btkeys[NEIGH_SIZE];
  int neigh_proc[NEIGH_SIZE];
  double elev[4];
  int i, j, k;
  int btflag, btype, has_involved, active;

  double len_scale = matprops->LENGTH_SCALE;
  double bucket_size = PARTICLE_DENSITY * (matprops->smoothing_length);

#ifdef DEBUG
   bool do_search = false;
   unsigned keycheck[KEYLENGTH] = {205661440, 304357968};
//   unsigned keytemp[KEYLENGTH] ;
#endif

  // set all flags to down
  for (i = 0; i < numprocs; i++)
    my_comm[i] = 0;

  //The following, it seems that they are useless
  // infinity
  double infty;
  if (numeric_limits < double >::has_infinity)
    infty = numeric_limits < double >::infinity();
  else
    infty = HUGE_VAL;

  // Read Hash-table related constants
  sprintf(filename, "funky%04d.h5", myid);

  // compare time-stamps of simulation.data and funkyxxxx.h5
  struct stat fstat1, fstat2;
  int file1 = open (filename, O_RDONLY);
  int file2 = open ("simulation.data", O_RDONLY);
  fstat (file1, & fstat1);
  fstat (file2, & fstat2);

  if ( fstat2.st_mtime > fstat1.st_mtime )
    fprintf (stderr,"WARNING: \"simulation.data\" is more recent than "
                    " \"%s\"\n", filename);
  close (file1);
  close (file2);

  hid_t fp = GH5_fopen_serial(filename, 'r');


  // Read Hash table constants
  GH5_readdata(fp, "/hashtable_constants", hvars);
  mindom[0] = hvars[0];         // min x
  maxdom[0] = hvars[1];         // max x
  mindom[1] = hvars[2];         // min y
  maxdom[1] = hvars[3];         // max y
  mindom[2] = hvars[4];         // min z
  maxdom[2] = hvars[5];         // max z

  // Create two new Hash-tables
  for (i = 0; i < DIMENSION; i++)
  {
    mindom[i] /= matprops->LENGTH_SCALE;
    maxdom[i] /= matprops->LENGTH_SCALE;
  }

  //initial domain info
  mindom_i[0]=simprops->Idom_x_min;
  mindom_i[1]=simprops->Idom_y_min;
  mindom_i[2]=simprops->Idom_z_min;

  maxdom_i[0]=simprops->Idom_x_max;
  maxdom_i[1]=simprops->Idom_y_max;
  maxdom_i[2]=simprops->Idom_z_max;

  //The initial domain is extended to have a layer of empty non-brief buckets ---> for adding of
  for (i = 0; i < DIMENSION; i++)
  {
    mindom_new[i] = mindom_i[i] - bucket_size;
    maxdom_new[i] = maxdom_i[i] + bucket_size;
  }

  // create hash-table for back-ground mesh
  *BG_mesh = new HashTable (BG_TABLE_SIZE, 2017, mindom, maxdom);

  // get the size of BG Mesh
  hsize_t dims[2];
  int num = GH5_getsize(fp, "/Buckets", dims);

  No_of_Buckets = (int) dims[0];

  // allocate memory for buckets
  bucket = new BucketStruct[No_of_Buckets];

  // read BG Mesh data
  GH5_read_grid_data(fp, "/Buckets", bucket);

  double flat[2*DIMENSION];
  flat[0]=Lx_P[0];
  flat[1]=Lx_P[1];
  flat[2]=Ly_P[0];
  flat[3]=Ly_P[1];
  flat[4]=Lz_P[0];
  flat[5]=Lz_P[1];
  for (i = 0; i < No_of_Buckets; i++)
  {
	// hash-table keys
    for (j = 0; j < KEYLENGTH; j++)
    	btkey[j] = bucket[i].key[j];

#ifdef DEBUG
    if (do_search)
    {
    	if (find_bucket (btkey, keycheck))
    		cout << "The bucket is found!" << endl;
    }
#endif

    // min coordinates
    min_crd[0] = bucket[i].xcoord[0] / len_scale;
    min_crd[1] = bucket[i].ycoord[0] / len_scale;
    min_crd[2] = bucket[i].zcoord[0] / len_scale;

    // max coordinates
    max_crd[0] = bucket[i].xcoord[1] / len_scale;
    max_crd[1] = bucket[i].ycoord[1] / len_scale;
    max_crd[2] = bucket[i].zcoord[1] / len_scale;

    if (bucket[i].myproc != myid)
    {
      fprintf(stderr, "ERROR: Input data is not correct. Aborting.\n");
      fprintf(stderr, "myid = %d, data_proc = %d\n", myid, bucket[i].myproc);
#ifdef MULTI_PROC
      MPI_Abort(MPI_COMM_WORLD, myid);
#else
      exit(1);
#endif
    }

    for (j = 0; j < NEIGH_SIZE; j++)
    {
      for (k = 0; k < KEYLENGTH; k++)
        tempbtkey.key[k] = bucket[i].neighs[j * KEYLENGTH + k];

      neigh_btkeys[j] = tempbtkey;
      neigh_proc[j] = bucket[i].neigh_proc[j];

      // turn on all my procs that communicate with this proc
      if ( neigh_proc[j] > -1 )
        my_comm[neigh_proc[j]] = 1;
    }

    active=0;
    active=determine_bucket_involved_flag(mindom_new, maxdom_new, min_crd, max_crd);

    if(active)
    {
        // buckettype flag
    	int bt_index[2*DIMENSION];
    	for (k=0; k<2*DIMENSION; k++)
    	    bt_index[k]= bucket[i].bucket_index[k];

    	for (j = 0; j < 4; j++)
    	    elev[j] = bucket[i].elev[j];

    	btype = bucket[i].buckettype;
    	switch (btype)
    	{
    	  case 1:
    	    btflag = UNDERGROUND;
    	    break;
    	  case 2:
    	    btflag = MIXED;
    	    if (bt_index[4] == -1)
    	       for (j = 0; j < 4; j++)
    	           elev[j] = bucket[i].elev[j] / len_scale;
    	    break;
    	  case 3:
    	    btflag = OVERGROUND;
    	    break;
    	  case 4:
    	    btflag = PRESS_BC;
    	    break;
    	  default:
    	    fprintf(stderr,"ERROR: Unknown buckettype flag.\nCheck the preoprocessor\n");
    	    exit(1);
    	}

    	// create a new bucket
    	Bucket * buck = new Bucket(btkey, min_crd, max_crd, btflag, elev,
                myid, neigh_proc, neigh_btkeys, bt_index, flat);

    	//determine has_involved flag for the bucket
    	has_involved=0;
    	has_involved=determine_bucket_involved_flag(mindom_i, maxdom_i, min_crd, max_crd);
    	// set has_potential involved
    	buck->set_has_potential_involved ((bool) has_involved); //Initially, all buckets "contains" the initial domain should be has_potential_involved.
    	(*BG_mesh)->add(btkey, buck);
    } //end of if bucket is active
    else
    {
    	// create a new bucket
    	BriefBucket * briefbuck = new BriefBucket(btkey, min_crd, myid, neigh_proc);
 	    (*BG_mesh)->add(btkey, briefbuck);
    }
  } //end of go through all buckets, loop index: i

  // please don't talk to yourself
  my_comm [myid] = 0;

  // free memory
  delete [] bucket;

  // clear out the partition table
  partition_tab.clear();

  // get the size of BG Mesh
  GH5_getsize (fp, "/partition_table", dims);
  if ((dims[0] <= 0))
  {
    fprintf (stderr,"Error! unable to read partition table.\n");
    exit (1);
  }

  int keylength = KEYLENGTH;
  unsigned * part_keys = new unsigned [dims[0]];//dims[0] should be the size of partition_table.
  int incr = 2 + keylength; //The incr is very important 2 is size of sfc-key in BucketHead, and length is keylength of bucket key
  GH5_readdata (fp, "/partition_table", part_keys);
  for (i = 0; i < dims[0]; i += incr)
  {
    //Bucket * buck = (Bucket *) (*BG_mesh)->lookup (part_keys + i + 2);
    partition_tab.push_back (BucketHead (part_keys + i, part_keys + i + 2));
  }

  // Create hash-table for particles
  *P_table = new THashTable(P_TABLE_SIZE, 2017, mindom, maxdom);
  return 0;
}
