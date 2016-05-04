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

////Define a situation where erupted material is pure air
//#ifndef ERUPT_PURE_AIR
//#define ERUPT_PURE_AIR 1
//#endif

////Define a situation where erupted material is the at the same temperature as air
//#ifndef ERUPT_COOL_MATERIAL
//#define ERUPT_COOL_MATERIAL 1
//#endif

//Define the atmosphere type
/*
 * 0: realistic -->Based on equations
 * 1: hydro-static
 * 2: uniform
 * 3: uniform-temperature, atmosphere stratefied due to gravity
 * 4: realistic interpolation --->read realistic atmosphere data and do interpolation to determine temperature, pressure, density
 */
#ifndef ATMOSPHERE_TYPE
#define ATMOSPHERE_TYPE 2 //The default value represents hydro-static atmosphere
#endif


//Define the erupt velocity profile type
/*
 * 0: uniform
 * 1: parabolic
 */
#ifndef ERUPT_VELOCITY_PROF
#define ERUPT_VELOCITY_PROF 1 //The default value represents hydro-static atmosphere
#endif

//Turn on and turn off domain adjusting ---> might be useful for performance analysis
// if define ADJUST_DOMAIN, will enable the feature of domain adjusting
// if domain adjusting feature is turned off, the initial domain (parameters giving in input file) should be consistent with computational domain
#ifndef ADJUST_DOMAIN
#define ADJUST_DOMAIN
#endif

//define different EOS
/*
 * 0: ideal gas
 * 1: weakly compressible ---> based on equation in "Weakly compressible SPH for free surface flows" by Markus Becker
 *                             The sound speed should be computed accordingly
 *                             In this case, energy equation is essentially decoupled from the other equations
 *                             If the flow is multiple phase flow, should be careful while use this EOS
 *              CAUTIOn: Atmosphere should be uniform when choose FLUID_COMPRESSIBILITY=1
 *                       Because when the fluid is water, the gradient of temperature and density is almost vanish
 */
#ifndef FLUID_COMPRESSIBILITY
#define FLUID_COMPRESSIBILITY 1
#endif

//---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//--------------------------------OPTIONS FOR DEBUG--------------------------------------------
#ifdef DEBUG

//output accumulative execute time, this might be helpful for performance analysis
#ifndef OUT_PUT_EXCUT_TIME
#define OUT_PUT_EXCUT_TIME
#endif

//Output ghost particles
#ifndef WRITE_GHOSTS
#define WRITE_GHOSTS
#endif

//output PID
#ifndef WRITE_PID
#define WRITE_PID
#endif

//output particle mass ---> particle mass is necessary in post process
#ifndef WRITE_PMASS
#define WRITE_PMASS
#endif

#endif //DEBUG
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
#endif /* OPTIONS_H */
