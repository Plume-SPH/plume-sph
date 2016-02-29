/*
 * meteo.cc
 *
 *  Created on: Feb 28, 2016
 *      Author: zhixuanc
 */

#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <algorithm>   // std::lower_bound, std::upper_bound, std::find
#include "meteo.h"
using namespace std;

//construction
Meteo::Meteo()
{
   number_of_data = 0;
   number_of_props = 0;
   data = NULL;
}

//The *input_mat should be pointer to a pointer
Meteo::Meteo(double* input_mat, int nr, int np)
{
//check data
  if(nr<=0)
  {
    cout<<"input matrix is empty or invalid row number!" <<endl;
//    return;
  }
  if(np<=0)
  {
    cout<<"input matrix is empty or invalid row number!" <<endl;
//    return;
  }

  number_of_data = nr;
  number_of_props = np;

  data = new vector<double>[number_of_props];
//put data
  int i, j;
  int nc = np + 1; //number of columes in mat is number of columns in data + 1
  for (i=0; i<number_of_data; i++)
  {
	  z.push_back(*(input_mat+i*nc));
	  for (j=0; j<number_of_props; j++)
	  {
		  data[j].push_back(*(input_mat+i*nc+j+1)); //j+1 because the first column is put into z
	  }
  }
//make sure data is sorted
  if (!(is_sorted(z.begin(),z.end())))
  {
 	 cout <<"input data is not sorted based on height" << endl;
 	 exit (EXIT_FAILURE);
  }
}

Meteo::~Meteo()
{
	delete data;
}
//
void Meteo::interpolate(double z0, double* value)
{
	int i;
	int index, li, ui;
	vector<double>::iterator low, up, fd;

	//if the data is sample data, output the sample data and no interpolation is needed
	fd = find (z.begin(), z.end(), z0);
	if (fd != z.end()) //means do not need interpolation
	{
		index = fd - z.begin();
		for (i = 0; i< number_of_props; i++)
			*(value+i)= data[i].at(index);
		return;
	}
	else
	{
		//before doing interpolation, first  find out whether z0 is out of range
		if (z0 < *(z.begin()))
			for (i = 0; i< number_of_props; i++)
				*(value+i) = *(data[i].begin());
		else if (z0 > *(z.end()-1))
			for (i = 0; i< number_of_props; i++)
				*(value+i) = *(data[i].end() - 1);
		else //if inside the range
		{
			// search for the data points that will be interpolate between.
			low = lower_bound (z.begin(), z.end(), z0);
			up = upper_bound (z.begin(), z.end(), z0);
			li = low - z.begin() - 1;
			ui = up - z.begin();
			//then do linear iterpolation
			for (i = 0; i< number_of_props; i++)
				*(value+i)= data[i].at(li) + (z0- z.at(li))/(z.at(ui)- z.at(li)) * (data[i].at(ui) - data[i].at(li));
		}
		return;
	}
}

//overloading of interpolate, only return specific property
/*
 * 0 density
 * 1 pressure
 * 2 temperature
 * 3 specific humidity
 * 4 wind velocity West->East
 * 5 wind velocity North->South
 */
double Meteo::interpolate(double z0, int prop_index)
{
	int i;
	int index, li, ui;
	vector<double>::iterator low, up, fd;

	//if the data is sample data, output the sample data and no interpolation is needed
	fd = find (z.begin(), z.end(), z0);
	if (fd != z.end()) //means do not need interpolation
	{
		index = fd - z.begin();
	    return data[prop_index].at(index);
	}
	else
	{
		//before doing interpolation, first  find out whether z0 is out of range
		if (z0 < *(z.begin()))
			return *(data[prop_index].begin());
		else if (z0 > *(z.end()-1))
			return *(data[prop_index].end() - 1);
		else //if inside the range
		{
			// search for the data points that will be interpolate between.
			low = lower_bound (z.begin(), z.end(), z0);
			up = upper_bound (z.begin(), z.end(), z0);
			li = low - z.begin() - 1;
			ui = up - z.begin();
			//then do linear iterpolation
			return (data[prop_index].at(li) + (z0- z.at(li))/(z.at(ui)- z.at(li)) * (data[prop_index].at(ui) - data[prop_index].at(li)));
		}
	}
}
int Meteo::get_number_of_props()
{
	return number_of_props;
}
