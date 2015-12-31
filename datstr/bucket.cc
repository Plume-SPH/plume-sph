/*
 * bucket.cc
 *
 *  Created on: Mar 5, 2015
 *      Author: zhixuanc
 */

#define PARABOLIC_PILE

#include <cstdio>
#include <cassert>
#include <cmath>
using namespace std;

#include <particle.h>
#include <hilbert.h>
#include <hashtab.h>
#include <thashtab.h>
#include <sph_header.h>
#include <constant.h>
#include "bucket.h"

int ineigh[3][3][3] = {{{13, 12, 14}, {10,  9, 11}, {16, 15, 17}},
                       {{4,   3,  5}, {1,   0,  2}, {7,   6,  8}},
                       {{22, 21, 23}, {19, 18, 20}, {25, 24, 26}}
                      };
//default constructor
BriefBucket::BriefBucket ()
{
  int i;

  myprocess = -1;
  for (i = 0; i < DIMENSION; i++)
    mincrd[i] = 0.;
  for (i = 0; i < NEIGH_SIZE; i++)
    neigh_proc[i] = -1;
  for (i = 0; i < KEYLENGTH; i++)
	  key.key[i] = 0;

  return;
}

BriefBucket::BriefBucket (unsigned *keyi, double *minx, int myid, int *nproc)
{
  int i;
  myprocess = myid;
  for (i = 0; i < KEYLENGTH; i++)
    key.key[i] = keyi[i];
  for (i = 0; i < DIMENSION; i++)
    mincrd[i] = minx[i];
  for (i = 0; i < NEIGH_SIZE; i++)
    neigh_proc[i] = nproc[i];
} // end constuctor

void
BriefBucket::put_neigh_proc (int dir[], int proc)
{
  int i = dir[0], j = dir[1], k = dir[2];

  neigh_proc[ineigh[i][j][k]] = proc;
  return;
}

int
BriefBucket::which_neigh_proc (int dir[]) const
{
  int i, j, k;

  i = dir[0];
  j = dir[1];
  k = dir[2];
  return neigh_proc[ineigh[i][j][k]];
}

// constructors with bucket_index as input
// bnd_xcrd ect does not matter
// bnd is determined.
// bnd_xcrd might be useful when ground is not flat
Bucket::Bucket (unsigned *keyi, double *minx, double *maxx, int buck_type,
                double *elev, int myid, int *nproc, Key * neigh, int* bt,
				double* flat) : BriefBucket (keyi, minx, myid, nproc)
{
  erupt_flag = false;
  int i;

  bucket_type = buck_type;
//  myprocess = myid;
  guest_flag = 0;
  active = false;
  particles_type = 0;
  has_involved = 0;

//  for (i = 0; i < KEYLENGTH; i++)
//    key.key[i] = keyi[i];

  for (i = 0; i < DIMENSION; i++)
  {
//    mincrd[i] = minx[i];
    maxcrd[i] = maxx[i];
  }
  for (i = 0; i < NEIGH_SIZE; i++)
  {
    neighbors[i] = neigh[i];
//    neigh_proc[i] = nproc[i];
  }

  for (i=0; i<2*DIMENSION; i++)
	  bucket_index[i]=bt[i];

  if (bucket_type == MIXED)// only mixed bucket has all boundary information
  {
	 for (i=0; i<2*DIMENSION; i++)
	 {
	     if (bt[i] != 0)
		    bnd[i]= flat[i];
	     else
	    	bnd[i]= 0.1;
	 }

  }
  else
  {
	for (i=0; i<2*DIMENSION; i++)
		bnd[i] = 0.1; //What else value could I use?
  }

  // list of particles in the bucket
  particles.clear ();
  new_plist.clear ();

//(i,j+1)_________(i+1,j+1)
//    |           |
//    |           |
//    |           |
//    |           |
//    |           |
//(i,j)-----------(i+1,j)
//

  //usually, poly is useless and set all value to zero
  for (i = 0; i < 4; i++)//This need to be modified in the future.
    poly[i] = 0.2; //0.2 for poly means useless for this kind of buckets.

  //only for on ground mixed buckets, it is necessary to use poly.
  if ((bucket_type == MIXED) && (bt[4] == -1))
  {

    // transform to local coordinate sys
    double xcrd[2], ycrd[2];
    xcrd[0] = 0;    //x(i)
    xcrd[1] = maxcrd[0] - mincrd[0];    //x(i+1)
    ycrd[0] = 0;    //y(j)
    ycrd[1] = maxcrd[1] - mincrd[1];    //y(j+1)

    // and transform elevs to local coordinates too
    for (i = 0; i < 4; i++)
      elev[i] -= mincrd[2];

    // fit the surface // 4 pts - 4 constants
    poly_surf (xcrd, ycrd, elev, poly);//How could I make sure that this ploy_surf also works for flat ground.

    for (i=0; i<4; i++)
      if (isnan (poly[i]))
        exit (51);
  }
} // end constuctor

Bucket::Bucket () : BriefBucket()
{
  erupt_flag = false;

  int i, j;

//  myprocess = -1;
  active = false;
  guest_flag = 0;
  bucket_type = 0;
  particles_type = 0;
  has_involved = 0;

  for (i = 0; i < DIMENSION; i++)
  {
//    mincrd[i] = 0.;
    maxcrd[i] = 0.;
  }

  for (i = 0; i < 4; i++)
    poly[i] = 0.2; //0.2 for poly means useless for this kind of buckets.


  for (i = 0; i < NEIGH_SIZE; i++)
  {
//    neigh_proc[i] = -1;
    for (j = 0; j < KEYLENGTH; j++)
      neighbors[i].key[j] = 0;
  }
  particles.clear ();
  new_plist.clear ();
  return;
}

/*! boundary normal from point \f$\mathbf{x}\f$, provided
 * /f$ \mathbf{x} \in \Gamma^i\f$.
 */
int
Bucket::get_bnd_normal (double point[], double normal[]) const
{
  int i;
  // if this isn't a boundary bucket, it shouldn't
  if ( bucket_type != MIXED )
    return 1;

  double pnt[DIMENSION-1];
  for (i = 0; i < DIMENSION-1; i++)
	  pnt[i] = point[i] - bnd[2*i];

  //The following expression also works for flat ground:
  //when the ground is flat the norm of the ground would be (0 0 1)
  normal[0] = -(poly[0] + poly[2] * pnt[1]);    // P1 + P3*x
  normal[1] = -(poly[1] + poly[2] * pnt[0]);    // P2 + P3*y
  normal[2] = 1.;

  // normalize
  double d = 0.;
  for (i = 0; i < DIMENSION; i++)
    d += normal[i] * normal[i];

  for (i = 0; i < DIMENSION; i++)
    normal[i] /= sqrt(d);
  return 0;
}

Key
Bucket::which_neigh (int dir[]) const
{
  int i, j, k;

  i = dir[0];
  j = dir[1];
  k = dir[2];
  return neighbors[ineigh[i][j][k]];
}

bool
Bucket::find_neigh_dir (Key keyin, int dir[]) const
{
  int i, j, k;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
      {
        if (neighbors[ineigh[i][j][k]] == keyin)
        {
          dir[0] = i;
          dir[1] = j;
          dir[2] = k;
          return true;
        }
      }
  return false;
}

/*
 * Normal Distance between point and boundary:
 * finite value is returned only if the
 * intersection between normal from point
 * and boundary is within current bucket,
 * else infinity is returned
 */
double
Bucket::get_bnddist (double pnt[], double intsct[]) const
{
  if (!calc_intersection_ground (pnt , intsct))//means that calculation of intersect point is correctly done.
  {
    double d = 0;
    for (int i=0; i<DIMENSION; i++)
      d += (pnt[i] - intsct[i]) * (pnt[i] - intsct[i]);
    return sqrt (d);
  }
    else
      return 1.E+10;
}

// overload of the put_ghost_particle with doing nothing.
//Because in my code, I never need put_new ghost particles, what I need is just update state of existing ghost particles.
int
Bucket::put_ghost_particles ()
{
	return 0;
}
/*! iterates for point \f$\mathbf{x}_b\f$ on boundary from point
 *  \f$\mathbf{x}\f$, outside the boundary, such that
 *  \f$ \mathbf{x}_b - \mathbf{x} \f$ is normal to the boundary,
 *  given it is in the form : \f$ a x + b y + c x y + d - z = 0 \f$
 *  Newton's method, upto 5 cycles is used.
 */
int
Bucket::calc_intersection_ground (double * point, double * xnew) const
{
  register int i, j;
  register double xold[DIMENSION];
  register double pl[4];
  double pt[DIMENSION], tmp[DIMENSION];
  double tol = 1.0E-5;
  double err = 0.;

  //There is two situations: flat ground and curve ground
  /*
   * pl is only useful for the situation when bucket is "ground MIXED"
   * for any other situations, the default value for ploy is 0.2 (an useless value)
   * When the bucket is "ground MIXED", value of poly will depend on the ground
   * If ground is flat, ploy[0-2]=0, poly[3]= bnd[4] ---that is incorrect! poly[3] = bnd[4] - mincrd[2], it is based on local coordinates;
   * If ground is not flat, at least one of poly[0-2] != 0;
   * */
  // hopefully it will stay in registers
  for (i = 0; i < 4; i++)
    pl[i] = poly[i];
  double ground_type;
  for (i=0; i<3; i++) //first 3 cof is zero, means the ground is flat.
	  ground_type += abs(pl[i]);

  //A more robust way:
  if ((bucket_type == MIXED) && (bucket_index[4] == -1))//to calculate the intersection, the bucket has to be on ground MIXED first
  {
	  if (ground_type == 0)//ground_type == 0 means that the ground BC is a flat ground.
	    {
            //I am assuming that the z direction is the elevation direction
		    xnew[0] = *point;
	        xnew[1] = *(point+1);
	        xnew[2] = bnd[4];
	    }
	  else
	    {
		  // transform point to local coordinates
		  for (i = 0; i < DIMENSION; i++)
			  pt[i] = point[i] - mincrd[i];

	  	  // inital guess
	  	  xnew[0] = 0.5 * (maxcrd[0] - mincrd[0]);
	  	  xnew[1] = 0.5 * (maxcrd[1] - mincrd[1]);
	  	  xnew[2] = pl[0] * xnew[0] + pl[1] * xnew[1] +
	  	            pl[2] * xnew[0] * xnew[1] + pl[3];

	  	  // Newton's method
	  	  for (i=0; i<8; i++)
	  	  {
	  	    // copy x(n+1) to x(n)
	  	    for (j = 0; j < DIMENSION; j++)
	  	      xold[j] = xnew[j];

	  	    // function at xold
	  	    double fn[3] = { (pl[0] * xold[0]) + (pl[1] * xold[1]) +
	  	                     (pl[2] * xold[0] * xold[1]) + pl[3] - xold[2]     ,
	  	                     (xold[0] - pt[0]) +
	  	                     ((pl[0] + (pl[2] * xold[0])) * (xold[2] - pt[2])) ,
	  	                     (xold[1] - pt[1]) +
	  	                     ((pl[1] + (pl[2] * xold[1])) * (xold[2] - pt[2]))
	  	                   };

	  	    // jacobian at xold
	  	    double a1 = pl[0] + pl[2] * xold[1];
	  	    double a2 = pl[1] + pl[2] * xold[0];
	  	    double a3 = pl[2] * (xold[2] - pt[2]);

	  	    /* jacobian = {{ a1, a2, -1},
	  	                   { 1,  a3, a1},
	  	                   { a3, 1,  a2}
	  	                  }
	  	     */
	  	    // --- begin --- maple generated code
	  	    double t1 = (-a3 *  a2 + a1);
	  	    double t2 = ( a2 *  a2);
	  	    double t3 = ( a1 *  a1);
	  	    double t4 = ( a2 *  a1);
	  	    double t5 = t3 + t2 + (-2 * t4 - a3) * a3 + 1;
	  	    t4 = t4 + a3;
	  	    t5 = 1. / t5;
	  	    double t6 = -a2 + a1 * a3;
	  	    tmp[0] = (t1 * fn[0] + (t2 + 1.) * fn[1] - t4 * fn[2]) * t5;
	  	    tmp[1] = (-t6 * fn[0] - t4 * fn[1] + (t3 + 1.0) * fn[2]) * t5;
	  	    tmp[2] = ((-1 + a3 * a3) * fn[0] + t1 * fn[1] - t6 * fn[2]) * t5;
	  	    // --- end ----- maple generated code

	  	    // upate x(n+1)
	  	    for (j = 0; j < DIMENSION; j++)
	  	      xnew[j] = xold[j] - tmp[j];

	  	    // check error
	  	    err = 0;
	  	    for (j = 0; j < DIMENSION; j++)
	  	      err += (xnew[j] - xold[j])*(xnew[j] - xold[j]);
	  	    if (sqrt (err) < tol ) break;
	  	  }
	  	  if ( sqrt (err) > tol ) return 1;

	  	  // transform back to global coordinates
	  	  for (i = 0; i < DIMENSION; i++)
	  	    xnew[i] += mincrd[i];
	    }
  }
  else
	 std::cout << "You are trying to find a on ground intersection in a bucket which is not a on ground MIXED! \n" << endl;

  return 0;
}

/*function that used to determine whether the real particle will escape the domain and lost
 * Only consider Mixed Bucket ---> This function should only be called by MIXED bucket which has correct bnd info.
 * The input pos should be coordinates of an real particle
 * */
bool
Bucket::determine_escape (double * pos)
{
	int i;
	for (i = 0; i < DIMENSION; i++)
	{
		if (bucket_index[i*2] == -1) //If the Mixed bucket is at "left side" of the domain
		{
			if ((pos[i] - bnd[i*2]) < -TINY)
				return true;
		}
		else if (bucket_index[i*2+1] == 1) //If the Mixed bucket is at "right side" of the domain
		{
			if ((pos[i] - bnd[i*2+1]) > TINY)
				return true;
		}
	}
	return false;
}
