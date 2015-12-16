/*
 * options.h
 *
 *  Created on: Oct 18, 2015
 *      Author: zhixuanc
 */

#ifndef OPTIONS_H
#define OPTIONS_H

////using Gaussian Kernel
#ifndef USE_GAUSSIAN
#define USE_GAUSSIAN
#endif

//define use summation method to update density
#ifndef USE_SUMMATION
#define USE_SUMMATION
#endif

//Define have LANS turbulent model in the code
//--> it is a stupid idea to do module management in C++ in this way, I should make use of the template, inherit, overloading as much as possible
#ifndef HAVE_TURBULENCE_LANS
#define HAVE_TURBULENCE_LANS
#endif

////Define have physics viscosity
//#ifndef USE_PHYSICS_VIS
//#define USE_PHYSICS_VIS
//#endif

//Have heat transfer
//#ifndef HAVE_HEAT_TRANSFER
//#define HAVE_HEAT_TRANSFER
//#endif

////Define eruption time period, so the eruption will only happen during a period of the whole simulation
//#ifndef HAVE_ERUPT_PERIOD
//#define HAVE_ERUPT_PERIOD

//No eruption
//#ifndef NO_ERUPTION
//#define NO_ERUPTION 1
//#endif

////Define a situation where erupted material is pure air
//#ifndef ERUPT_PURE_AIR
//#define ERUPT_PURE_AIR 1
//#endif

////Define a situation where erupted material is the at the same temperature as air
//#ifndef ERUPT_COOL_MATERIAL
//#define ERUPT_COOL_MATERIAL 1
//#endif

////Define a situation where atmosphere is uniform property ---> also no gravity
//#ifndef UNIFORM_ATMOSPHERE
//#define UNIFORM_ATMOSPHERE
//#endif

//---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//--------------------------------OPTIONS FOR DEBUG--------------------------------------------
#ifdef DEBUG

#ifndef OUT_PUT_EXCUT_TIME
#define OUT_PUT_EXCUT_TIME
#endif

//Output ghost particles
//#ifndef WRITE_GHOSTS
//#define WRITE_GHOSTS
//#endif

//output PID
#ifndef WRITE_PID
#define WRITE_PID
#endif

#endif
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
#endif /* OPTIONS_H */
