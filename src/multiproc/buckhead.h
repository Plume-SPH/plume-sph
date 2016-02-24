
/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id:$ 
 */

#ifndef BUCKET_HEAD_H
#  define BUCKET_HEAD_H

#  include <hashtab.h>

class BucketHead
{
protected:
  unsigned sfc_key[2];
  unsigned buck_key[KEYLENGTH];

public:
  //! constructor
  BucketHead(unsigned * skey, unsigned * bkey)
  {
    int i;
    for (i = 0; i < 2; i++)
      sfc_key[i] = skey[i];

    for (i = 0; i < KEYLENGTH; i++)
      buck_key[i] = bkey[i];
  }

  //! get key of first bucket
  const unsigned * get_buck_head () const
  {
    return buck_key;
  }

  //! get key for the current object
  const unsigned * getKey () const
  {
    return sfc_key;
  }

  //! comparison operator for sorting-->comparison is based on SFC-key
  bool operator < (const BucketHead & rhs) const
  {
    if ( sfc_key[0] < rhs.getKey ()[0] )
      return true;
    else if ( sfc_key[0] > rhs.getKey ()[0] )
      return false;
    else if ( sfc_key[1] < rhs.getKey ()[1] )
      return true;
    else
      return false;
  }
};

#endif //BUCKET_HEAD__H
