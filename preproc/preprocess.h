/*
 * preprocess.h
 *
 *  Created on: Mar 7, 2015
 *      Author: zhixuanc
 */

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include <vector>
#include <iostream>
using namespace std;

#include <buckstr.h>

struct ColumnHead
{
  int xind, yind, proc;
  unsigned key[KEYLENGTH];

  // constructor
  ColumnHead (int i, int j, unsigned keyi[])
  {
    xind = i;
    yind = j;
    proc = 0;
    for (i = 0; i < KEYLENGTH; i++)
      key[i] = keyi[i];
  }


    bool operator < (const ColumnHead & rhs) const
    {
      if ( key[0] < rhs.key[0] )
        return true;
      else if ( key[0] > rhs.key[0] )
        return false;
      else if ( key[1] < rhs.key[1] )
        return true;
      else
        return false;
    }
};

//! Write Background mesh and particle data to HDF5 file
void createfunky(
                 //! Numboer of processes in a multiproc run
                 int ,
                 //! Size of Hash Table constants array
                 int,
                 //! Hash Table Constants
                 double *,
                 //! Background grid data
                 vector <BucketStruct> &,
                 //! partition table keys
                 vector <unsigned> &
                );

//! Generate key based on location
void determine_the_key (
                //! Normalized coordinates
                double ,
                //! Keylength
                unsigned ,
                //! key to be generated
                unsigned ,
                //! Max key
                unsigned ,
                //! Min key
                unsigned
                );
//function that used to determine the type of bucket
void determine_bucket_type
               (
               //min of domain
               double *,
               //max of domain
               double *,
               // [xmin, xmax] of bucket
               double *,
               // [ymin, ymax] of bucket
               double *,
               // [zmin, zmax] of bucket
               double *,
               //return buckete type
               int *,
               //return bucket index
               int *
               );


#endif /* PREPROCESS_H_ */
