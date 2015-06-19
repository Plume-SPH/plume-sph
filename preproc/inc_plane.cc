/*
 * =====================================================================================
 *
 *       Filename:  inc_plane.cc
 *
 *    Description:  
 *
 *        Created:  07/26/2010 03:15:45 PM EDT
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *        License:  GNU General Public License
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */

#include <cmath>
#include <constant.h>
#include "gisapi.h"

double abs ( double a)
{
  return (a > 0 ? a : -a);
}

double GIS_get_elevation (double x, int *polytype, double smlen)
{
  double dx = 5*smlen;
  double x1 = 1.0;
  double z1 = x1*Slope + Intcpt;
  double z = Slope * x + Intcpt;

  if ( x >= x1 )
    return z;

  double hbox = 0.09;
  double xbox = hbox*sin(inc_angle);

  // make sure the step exits
  if ( xbox < dx )
    xbox = dx;

  double x2   = x1-xbox;
  double z2  = z1 + hbox*cos(inc_angle);
  if ( (x >= x2) && (x < x1 ) )
    return (z1 - (x - x1)/Slope);

  double lbox = 2*hbox;
  double x3 = x2 - lbox*cos(inc_angle);
  double z3 = z2 + Slope*(x3 - x2);
  if ((x >= x3) && (x < x2))
    return (z2 + Slope*(x - x2));

  double x4 = x3 - dx;
  double z4 = Slope*x4 + Intcpt;
  
  if ((x >= x4) && ( x < x3 ))
    return (z4 + ((z3-z4)/(x3-x4))*(x - x4));

  return z;
}

int Get_xmin(double *crd)
{
  *crd = 0.5;
  return 0;
}

int Get_ymin(double *crd)
{
  *crd = 0.4;
  return 0;
}

int Get_xmax(double *crd)
{
  *crd = 3.5;
  return 0;
}

int Get_ymax (double *crd)
{
  *crd = 2.0;
  return 0;
}
