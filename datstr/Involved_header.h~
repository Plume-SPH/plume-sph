/*
 * Involved_header.h
 *
 *  Created on: Jul 24, 2015
 *      Author: zhixuanc
 */

#ifndef INVOLVED_HEADER_H_
#define INVOLVED_HEADER_H_
#  include <constant.h>

struct InvolvedHead
{
  //!  process id of bucket containing image
  int has_involved;
  //!  key of bucket that contains image

  unsigned bucket_key[KEYLENGTH];

  // constructor 1
  InvolvedHead ()
  {
	has_involved = 0;
    int i;
    for (i = 0; i < KEYLENGTH; i++)
      bucket_key[i] = 0;
  }

  // constructor 2
  InvolvedHead (int inv, unsigned *bkey)
  {
	has_involved = inv;
    int i;
    for (i = 0; i < KEYLENGTH; i++)
      bucket_key[i] = *(bkey + i);

  }
};

#endif /* INVOLVED_HEADER_H_ */
