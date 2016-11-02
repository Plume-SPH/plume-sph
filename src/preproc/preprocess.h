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

struct PartiHead
{
  int xind, yind, zind, proc;
  unsigned key[KEYLENGTH];

  bool is_underground;
  unsigned top_key[KEYLENGTH]; //key of the bucket above it, newly added, will be used only when the bucket is underground bucket.

  // constructor
  PartiHead (int i, int j, int k, unsigned keyi[])
  {
    xind = i;
    yind = j;
    zind = k; //This is newly added, as a new member int k is added
    proc = 0;
    for (i = 0; i < KEYLENGTH; i++)
      key[i] = keyi[i];

    is_underground = false;

    for (i = 0; i < KEYLENGTH; i++)
      top_key[i] = 0;

  }

  // constructor overload --> constructor for underground bucket
  PartiHead (int i, int j, int k, unsigned keyi[], unsigned top_keyi[])
  {
    xind = i;
    yind = j;
    zind = k; //This is newly added, as a new member int k is added
    proc = 0;
    for (i = 0; i < KEYLENGTH; i++)
      key[i] = keyi[i];

    is_underground = true;

    for (i = 0; i < KEYLENGTH; i++)
      top_key[i] = top_keyi[i];

  }

  //The basic idea of this operator overload:
  /*
   * For underground bucket, use the top_key for comparison
   * 1) 0.1 is smaller than 1, which is the minimum unit for keys
   * 2) By adding 0.1, the  underground bucket will always be the one after its top bucket in a sorted list
   */
      bool operator < (const PartiHead & rhs) const
      {
    	  if (is_underground)
    	  {
        	  if (rhs.is_underground)
        	  {
         	     if ( top_key[0] < rhs.top_key[0] )
         	       return true;
         	     else if ( top_key[0] > rhs.top_key[0] )
         	       return false;
         	     else if ( (top_key[1] + 0.1)< (rhs.top_key[1] + 0.1)) //set the key of underground bucket to be a little bit larger than its above.
         	       return true;
         	     else
         	       return false;
        	  } //end of is rhs is underground
        	  else
        	  {
        	     if ( top_key[0] < rhs.key[0] )
        	       return true;
        	     else if ( top_key[0] > rhs.key[0] )
        	       return false;
        	     else if ( (top_key[1] + 0.1) < rhs.key[1] )
        	       return true;
        	     else
        	       return false;
        	  }//end of is rhs is not underground
    	  }// end of if bucket is underground
    	  else
    	  {
        	  if (rhs.is_underground)
        	  {
         	     if ( key[0] < rhs.top_key[0] )
         	       return true;
         	     else if ( key[0] > rhs.top_key[0] )
         	       return false;
         	     else if ( key[1] < (rhs.top_key[1] + 0.1)) //set the key of underground bucket to be a little bit larger than its above.
         	       return true;
         	     else
         	       return false;
        	  }//end of is rhs is underground
        	  else
        	  {
        	     if ( key[0] < rhs.key[0] )
        	       return true;
        	     else if ( key[0] > rhs.key[0] )
        	       return false;
        	     else if ( key[1] < rhs.key[1] )
        	       return true;
        	     else
        	       return false;
        	  }//end of is rhs is not underground
    	  }// end of if bucket is not underground
      }

//    bool operator < (const PartiHead & rhs) const
//    {
//      if ( key[0] < rhs.key[0] )
//        return true;
//      else if ( key[0] > rhs.key[0] )
//        return false;
//      else if ( key[1] < rhs.key[1] )
//        return true;
//      else
//        return false;
//    }
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
