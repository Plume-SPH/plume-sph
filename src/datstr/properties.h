/*
 * properties.h
 *
 *  Created on: Feb 17, 2015
 *      Author: zhixuanc
 */

#ifndef PROPERTIES_H
#define PROPERTIES_H
#  include <list>
#  include <cmath>
using namespace std;

#  include <hashtab.h>
#  include <thashtab.h>

#  include "meteo.h"
#  include  "Involved_header.h"

struct TimeProps
{
  //! current non-dimensional time
  double time;

  //! current non-dimensional time-step
  double dtime;

  //! maximum simulation time in sec
  double max_time;

  //! non-dimensional max simulation time
  double ndmax_time;

  //! time ineterval after which output is to be written
  double timeoutput;

  //! non-dimensional timeoutput
  double ndtimeoutput;

  //! Scale to Normalize the time
  double TIME_SCALE;\

  //! start erupt time --> in seconds
  double stat_erupt;

  //! end erupt time --> in seconds
  double end_erupt;

  //! current time-step
  int step;

  //! count of previous outputs
  int ioutput;

  //! maximum time-steps allowed
  int max_steps;

  //! constructor to allocate default values
    TimeProps ()
  {
    max_time = 180.0;
    timeoutput = 1.0;
    TIME_SCALE = 1.0;
    ndmax_time = max_time = 180.;
    ndtimeoutput = 0.0;
    ioutput = 0;
    step = 0;
    max_steps = 10;
  }

  void incrtime (double *dt)
  {
    // first reduce dt to hit output or end time "exactly"
    if ((time + (*dt)) > ndtimeoutput)
      (*dt) = ndtimeoutput - time;
    if ((time + (*dt)) > ndmax_time)
      (*dt) = ndmax_time - time;
    // then increment time
    time += (*dt);
    dtime = (*dt);
    step++;
  }

  bool ifstart ()
  {
    return (step == 0);
  }                             //before first time step
  bool iffirst ()
  {
    return (step == 1);
  }                             //at first time step
  bool ifend ()
  {
    return ((time >= ndmax_time) || (step > max_steps));
  }

  bool is_int_time () //function that used to determine whether time is a integer or not, different from ifoutput, it does not change values of member data
  {
    if (time >= ndtimeoutput)
      return true;
    else
      return false;
  }

  bool ifoutput ()
  {
    if (time >= ndtimeoutput)
    {
      ioutput++;                //using ioutput eliminates roundoff
      ndtimeoutput = ((ioutput + 1) * timeoutput) / TIME_SCALE;
      return true;
    }
    else
      return false;
  }

  //function to determine whether the volcano is eruptting or not.
  //Return TRUE if it is erupt
  //Return FALSE if it is not erupt
  bool iferupt()
  {
	  if((time * TIME_SCALE) < stat_erupt)
		 return false;
	  else if ((time * TIME_SCALE) > end_erupt)
		 return false;
	  else
		 return true;
  }

  void chunktime (int *hours, int *minutes, double *seconds)
  {
    double dimtime = time * TIME_SCALE;

    *hours = ((int) dimtime) / 3600;
    *minutes = (((int) dimtime) % 3600) / 60;
    *seconds = dimtime - (double) (*hours * 3600 + *minutes * 60);
  }

  double timesec ()
  {
    return (time * TIME_SCALE);
  }

};


// This class is not necessary in my code
struct MatProps
{
  //! constant in equation of state
  double P_CONSTANT;

  //! slope limiting stuff
  double GAMMA;

  //! SPH smoothing length---> for phase1, sml in TimeProps is for phase2
  double smoothing_length;

  //! mass of individual particle ---> for phase1, sml in TimeProps is for phase2
  double particle_mass;

  //! length scaling factor
  double LENGTH_SCALE;

  //! gravity scaling factor
  double GRAVITY_SCALE;

  //! constructor allocates default properties
    MatProps ()
  {
    P_CONSTANT = 1.;
    GAMMA = 1.4;
    LENGTH_SCALE = 1.;
    GRAVITY_SCALE = 9.81;
//    TINY = 1.0E-06;
  }
};

struct SimProps
{
   double Idom_x_min;
   double Idom_x_max;
   double Idom_y_min;
   double Idom_y_max;
   double Idom_z_min;
   double Idom_z_max;
   double bot;
   double mass_of_phase2;
   // smooth length of erupted particle---> sml in Matprops is for air particles
   double sml_of_phase2;
   //Atmosphere data (including wind field)
   Meteo * meteo_data;

   SimProps()
   {
	   Idom_x_min = -2000;
	   Idom_x_max = 2000;
	   Idom_y_min = -2000;
	   Idom_y_max = 2000;
	   Idom_z_min = 0;
	   Idom_z_max = 7000;

	   bot = 0.;
	   mass_of_phase2 =0.;
	   sml_of_phase2 =0.;
	   meteo_data  = NULL;

   }

   ~SimProps()
   {
	   delete meteo_data;
   }

   // bot is the z coordinate of the lowest layer in the duct
   double get_bot () const
   {
 	  return bot;
   }

    // bot is the z coordinate of the lowest layer in the duct
    void update_bot (double val)
    {
   	   bot = val;
    }
};
#endif /* PROPERTIES_H_ */
