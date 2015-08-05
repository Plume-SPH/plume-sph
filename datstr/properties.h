/*
 * properties.h
 *
 *  Created on: Feb 17, 2015
 *      Author: zhixuanc
 */

#ifndef PROPERTIES_H
#define PROPERTIES_H

#  include <cmath>
using namespace std;

#  include <hashtab.h>
#  include <thashtab.h>

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

  //! current time-step
  int step;

  //! count of previous outputs
  int ioutput;

  //! maximum time-steps allowed
  int max_steps;

  //! non-dimensional t_add, used to determine new layers of particles should be added or not.
  //! It is essentially the time from previous adding of new layers
  double t_add;

  //! non-dimensional t_each, the time that one erupt ghost particle need to move through a length of h (smooth length)--->computed based on average velocity
  double t_each;

  //! the z coordinate of the lowest layer of erupt ghost particles.
  double bot;
  //! the number of particles per layer in duct
  double cof;
  //! Scale to Normalize the time
  double TIME_SCALE;
  //mass of erupted particle---> mass in Matprops is for air particles
  double mass;
  // smooth length of erupted particle---> sml in Matprops is for air particles
  double sml;
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

    t_add=0.;
    t_each=1.;
    bot = 0;
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

  void update_tadd (double dt)
  {
	  t_add += dt;
  }

  void set_back_tadd (double value)
  {
	  t_add = value;
  }

  void update_teach (double value)
  {
	  t_each = value;
  }
  // bot is the z coordinate of the lowest layer in the duct
  void update_bot (double val)
  {
	  bot = val;
  }

  double get_tadd()
  {
	  return t_add;
  }

//  //get dimensional dt
//  double get_dt()
//  {
//	  return (dtime * TIME_SCALE);
//  }


//time for adding a new layer, t_each = dist/vel;
  //where dist is the distance between eahch layer
  //and vel is the eruption velocity
  double get_teach()
  {
	  return t_each;
  }
  // bot is the z coordinate of the lowest layer in the duct
  double get_bot ()
  {
	  return bot;
  }
//  //normalized time
//  double get_ntime ()
//  {
//	  return ndmax_time;
//  }

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

  //! SPH smoothing length
  double smoothing_length;

  //! mass of individual particle
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

   SimProps()
   {
	   Idom_x_min = -2000;
	   Idom_x_max = 2000;
	   Idom_y_min = -2000;
	   Idom_y_max = 2000;
	   Idom_z_min = 0;
	   Idom_z_max = 7000;
   }
};
#endif /* PROPERTIES_H_ */
