/*
 * hilbert.h
 *
 *  Created on: Mar 5, 2015
 *      Author: zhixuanc
 */

#ifndef HILBERT_H
#define HILBERT_H

/* 2D Hilbert Space-filling curve */

#  ifdef __cplusplus
extern "C"
{
#  endif

  void hsfc2d (unsigned coord[],        /* IN: Normalized integer coordinates */
               unsigned *nkey,  /* IN: Word length of key */
               unsigned key[]   /* OUT: space-filling curve key */
    );

/* 3D Hilbert Space-filling curve */

  void hsfc3d (unsigned coord[],        /* IN: Normalized integer coordinates */
               unsigned *nkey,  /* IN: Word length of 'key' */
               unsigned key[]   /* OUT: space-filling curve key */
    );

/* API for 2-D Hilbert Space Filling Curve */
  void HSFC2d (double coord[],  /* IN: Normalized floating point coordinates */
               unsigned *nkey,  /* IN: Word length of key */
               unsigned key[]   /* OUT: space-filling curve key */
    );

/* API for 3-D Hilbert Space Filling Curve */
  void HSFC3d (double coord[],  /* IN: Normalized floating point coordinates */
               unsigned *nkey,  /* IN: Word length of key */
               unsigned key[]   /* OUT: space-filling curve key */
    );

/* API for 2-D Hilbert Space Filling Curve with time indicator*/
    void THSFC2d (double coord[],  /* IN: Normalized floating point coordinates */
    		     unsigned add_step,     /* IN: Normalized time*/
                 unsigned *nkey,  /* IN: Word length of key */
                 unsigned key[]   /* OUT: space-filling curve key */
      );

/* API for 3-D Hilbert Space Filling Curve with time indicator*/
    void THSFC3d (double coord[],  /* IN: Normalized floating point coordinates */
    		     unsigned add_step,     /* IN: Normalized time*/
                 unsigned *nkey,  /* IN: Word length of key */
                 unsigned key[]   /* OUT: space-filling curve key */
      );
#  ifdef __cplusplus
}
#  endif

#endif /* HILBERT_H_ */
