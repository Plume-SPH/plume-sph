/*
 * options.h
 *
 *  Created on: Oct 18, 2015
 *      Author: zhixuanc
 */

#ifndef OPTIONS_H
#define OPTIONS_H

////using Gaussian Kernel currently, only one kind of kernel is available
#ifndef USE_GAUSSIAN
#define USE_GAUSSIAN
#endif

//define use summation method to update density  ---> Another way is based on discretized mass conservation
#ifndef USE_SUMMATION
#define USE_SUMMATION

//define to use GSPH or not ---> I need move this option to configure after test is done
/*
 * USE_GSPH 0 : SPH
 * USE_GSPH 1 : GSPH
 */

#ifndef USE_GSPH
#define USE_GSPH 1
#endif

//Defines the way to compute derivative. ---> This option only works when USE_GSPH = 1
/*
 * NORM_DERIVATIVE 0: derivative is not normalized
 * NORM_DERIVATIVE 1: derivative is normalized
 */
#ifndef NORM_DERIVATIVE
#define NORM_DERIVATIVE 1
#endif


//Define whether use adaptive smoothing length or not ---> adaptively adjust sml at a given interval, to avoid sml update at every time step
/*
 * 0: Not adaptive  -->DENSITY_UPDATE_SML an be 1 or 0
 * 1: adaptive      -->DENSITY_UPDATE_SML shold always be 1
 */
#ifndef ADAPTIVE_SML
#define ADAPTIVE_SML 0
#endif

// only when no sml adaptive is used, it will be necessary to decide to use either original sml or current sml. If sml is adaptive, always use current sml
/*
 * Define the smoothing length that will be used in density update: ---> only need when SPH summation formulism is used!
 * 0 : use original smoothing length   ---> It is easier to get negative energy
 * 1 : use current smoothing length    ----> More stable choice
 * Note: to make sure conservation of momentum and energy conservation, smoothing length might be changed so that sml for two phases are equal in mixing region
 * Any way, SPH will have some trouble if the smoothing length is different for two different phases.
 */
#ifndef DENSITY_UPDATE_SML
#define DENSITY_UPDATE_SML 1
#endif


/*
 * Based on assumption of immediate thermodynamics equilibrium, a internal energy smooth might be necessary to make this assumption to be true
 * 0: do not use energy smooth
 * 1: use energy smooth with normalization
 * 2: use energy smooth without normalization
 *
 * Note: this option is not available for GSPH yet.
 */
#ifndef HAVE_ENERGY_SMOOTH
#define HAVE_ENERGY_SMOOTH 0
#endif

/*
 * There are two ways to view particles of two phases in multiphase SPH method:
 * 1: They are nothing but discretized points, phase num of each particle is just a flag of that particle (or properties of particle)
// * 10: The same as option 1 update of density will not based on SPH, instead, it will based on an equation (Can be found in Suzuki's 2005 paper)
// * 11: The same as option 1 , except that does not do normalization. ---> Original smoothing length should be used.
 * 12 : The same as option 1, the only difference is that use phase density instead of total density in computing normalization.
 * 2: Two different sets of discretized points. essentially independent and only interact with each other by the interact terms (include explicit terms like drag force or implicit terms such as the pressure force term.) in the governing equation.
 * 22 : The same as option 2, the only difference is that use the total density instead of phase density in computing normalization.
 */
#ifndef DENSITY_UPDATE_SPH
#define DENSITY_UPDATE_SPH 1 //DENSITY_UPDATE_SPH=2, view particle as two different phases will lead to very large density
                           //When DENSITY_UPDATE_SPH=0, DENSITY_UPDATE_SML should be set to 1;
#endif

#endif // end of USE_SUMMATION



//Define which format of to use for discretized momentum equation
/*
 * 0: The basic symmetric format which can conserve momentum, without any further modification
 * 1: Based on 0, an external pressure is deduce by every pressure, the external pressure the pressure of atmosphere at corresponding height of particle a. ---> The purpose of this is to make sure when pressure gradient vanish the acceleration will be zero.
 *    Note: Please note that for GSPH, such tricky is not necessary, as GSPH can guarantee zero RHS of momentum equation when pressure gradient is zero.
 *    So you should not use this option together with
 */
#ifndef MOMENTUM_DISCRETIZE
#define MOMENTUM_DISCRETIZE 0
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


//Define the erupt velocity profile type --> for plume modeling
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
#define TIME_STEP_CONSTRAIN 0 //The default value represents hydro-static atmosphere
#endif

/*
 * Define whether use cut method to avoid negative energy
 * The purpose of this is to avoid negative energy which is not physical
 * Why I will have negative energy is not clear yet, ---> Need more analysis on SPH method
 * The cut method is a very natural way to avoid negative energy. ---> At least one paper reports they are using this method.
 * "An SPH model for multiphase flows with complex interface and large density difference" by Z. Chen
 *
 * 0: No, do not use energy cut
 * 1: yes, use energy cut
 */
#ifndef HAVE_ENERGY_CUT
#define HAVE_ENERGY_CUT 0
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

//output particle smoothing length ---> useful for debug
#ifndef WRITE_SML
#define WRITE_SML
#endif

#endif //DEBUG
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
#endif /* OPTIONS_H */
