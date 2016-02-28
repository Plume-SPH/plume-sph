/*
 * meteo.h
 *
 *  Created on: Feb 27, 2016
 *      Author: zhixuanc
 */

#ifndef METEO_H
#define METEO_H

#include <iostream>
#include <vector>

using namespace std;

/*
 * class for meteo data
 * will do linear interpolation for points within data range
 * will do zeros order exterpolation for points out of data range
 */
class Meteo
{
protected:
   int number_of_data;
   int number_of_props;
   vector<double> z;
   vector<double> *data;
public:
   Meteo();
   Meteo(double*, int, int);
   ~Meteo();
   //return interpolation value of all properties
   void interpolate (double, double*);
   //overload of interpolation -->only return specific property
   double interpolate (double, int);
   int get_number_of_props();
};

#endif /*METEO_H*/
