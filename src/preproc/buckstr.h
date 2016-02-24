/*
 * buckstr.h
 *
 *  Created on: Mar 7, 2015
 *      Author: zhixuanc
 */

#ifndef BUCKSTR_H
#define BUCKSTR_H


#include <constant.h>

struct BucketStruct
{
  unsigned key[KEYLENGTH];
  unsigned neighs[NEIGH_SIZE*KEYLENGTH]; //neighbor buckets
  int      buckettype; // bucket type: Mixed, over ground ect....
  int      myproc;
  int      neigh_proc[NEIGH_SIZE];//corresponding process id of neighbor buckets
  int      bucket_index[2*DIMENSION];
  double   xcoord[2];
  double   ycoord[2];
  double   zcoord[2];
  double   elev[4]; //elev only needed for onground MIXED bucket.
};

#endif /* BUCKSTR_H_ */
