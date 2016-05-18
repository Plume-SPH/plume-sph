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

/*
 * Define the smoothing length that will be used in density update: ---> only need when SPH summation formulism is used!
 * 0 : use original smoothing length
 * 1 : use current smoothing length
 * Note: to make sure conservation of momentum and energy conservation, smoothing length might be changed so that sml for two phases are equal in mixing region
 * Any way, SPH will have some trouble if the smoothing length is different for two different phases.
 */
#ifndef DENSITY_UPDATE_SML
#define DENSITY_UPDATE_SML 1
#endif

/*
 * There are two ways to view particles of two phases in multiphase SPH method:
 * 0: update of density will not based on SPH, instead, it will based on an equation (Can be found in Suzuki's 2005 paper)
 * 1: They are nothing but discretized points, phase num of each particle is just a flag of that particle (or properties of particle)
 * 2: Two different sets of discretized points. essentially independent and only interact with each other by the interact terms (include explicit terms like drag force or implicit terms such as the pressure force term.) in the governing equation.
 */
#ifndef DENSITY_UPDATE_SPH
#define DENSITY_UPDATE_SPH 1  //It seems that only DENSITY_UPDATE_SPH=1 works well
#endif

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


//Define the atmosphere type
/*
 * 0: realistic -->Based on equations
 * 1: hydro-static --->With temperature gradient (Described in Suzuki's 2005 plume modeling paper), based on pure hydro-static relation
 * 2: uniform    ---> No gravity, no temperature gradient, no density gradient
 * 3: uniform-temperature, atmosphere stratefied due to gravity
 * 4: realistic interpolation --->read realistic atmosphere data and do interpolation to determine temperature, pressure, density
 */
#ifndef ATMOSPHERE_TYPE
#define ATMOSPHERE_TYPE 4
//The default value represents hydro-static atmosphere
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
#define FLUID_COMPRESSIBILITY 0
#endif

//Define the time step constrain
/*
 * 0: Only CFL
 * 1: CFL + boundary stablization
 */
#ifndef TIME_STEP_CONSTRAIN
#define TIME_STEP_CONSTRAIN 1 //The default value represents hydro-static atmosphere
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
