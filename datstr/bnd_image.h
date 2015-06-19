/*
 * bnd_image.h
 *
 *  Created on: Mar 13, 2015
 *      Author: zhixuanc
 */

#ifndef BND_IMAGE_H
#define BND_IMAGE_H
#  include <constant.h>
struct BndImage
{
  //!  process id of bucket containing image
  int buckproc;
  //!  process id of ghost particle
  int partproc;
  //!  key of bucket that contains image
  unsigned bucket_key[KEYLENGTH];
  //!  key of ghost whose image it is
  unsigned ghost_key[TKEYLENGTH];
  //!  coordinate of image
  double coord[DIMENSION];
  //!  interpolate values of state_vars
  double state_vars[NO_OF_EQNS];

  // constructor 1
  BndImage ()
  {
	int i ;
    buckproc = -1;
    partproc = -1;
    state_vars[0] = 1.0;
    for (i = 1; i < NO_OF_EQNS; i++)
      state_vars[i] = 0;
  }

  // constructor 2
  BndImage (int p1, int p2, unsigned *bkey, unsigned *pkey, double *crd)
  {
    int i;

    buckproc = p1;
    partproc = p2;
    // keys
    for (i = 0; i < KEYLENGTH; i++)
      bucket_key[i] = *(bkey + i);
    for (i = 0; i < TKEYLENGTH; i++)
      ghost_key[i] = *(pkey + i);
    // coords
    for (i = 0; i < DIMENSION; i++)
      coord[i] = *(crd + i);

    // inter polated values
    state_vars[0] = 1.;
    for (i = 1; i < NO_OF_EQNS; i++)
      state_vars[i] = 0;
  }
};



#endif /* BND_IMAGE_H_ */
