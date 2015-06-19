/*
 * =====================================================================================
 *
 *       Filename:  gisapi.h
 *
 *    Description:  
 *
 *        Created:  01/17/2011 11:19:26 AM EST
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

#ifndef GISAPI__H
#define GISAPI__H

const int HYBRID = 0x3;

const double inc_angle = 0.41888;     // 24 deg
const double Slope  = 0.4452; // tan(24)
const double Intcpt = 0.24;   // wess's exp

const double Radius = 0.08;
const double xCen   = 0.8;

double GIS_get_elevation (double , int *, double );
int    Get_xmin(double *);
int    Get_xmax(double *);
int    Get_ymin(double *);
int    Get_ymax(double *);

#endif // GISAPI__H
