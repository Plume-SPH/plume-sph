/*
 * options.h
 *
 *  Created on: Oct 18, 2015
 *      Author: zhixuanc
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

//##################################################################################################################################################################
#if CODE_DIMENSION==1
////using Gaussian Kernel currently, only one kind of kernel is available
#ifndef USE_GAUSSIAN
#define USE_GAUSSIAN
#endif

//define use summation method to update density  ---> Another way is based on discretized mass conservation
/*
 * USE_SUMMATION=1: Use weighted summation formulation for density update
 * USE_SUMMATION=0: Use discretized mass conservation equation ---> Currently, did not consider updating for mass fraction yet. ---> Is this one works well with shock tube, try it for Plume model
 */
#ifndef USE_SUMMATION
#define USE_SUMMATION 1

// only when no sml adaptive is used, it will be necessary to decide to use either original sml or current sml. If sml is adaptive, always use current sml
// NOTE: currently, all of these schemes are only applied to density updating. Actually, some of these ideas can be used in momentum and energy updating.
/*
 * Define the smoothing length that will be used in density update: ---> only need when SPH summation formulism is used!
 * 0 : use original smoothing length   ---> It is easier to get negative energy
 * 1 : use current smoothing length hi  ----> More stable choice
 * 11: same as option1, the only difference is use (hi+hj)*0.5 instead of hi
 * 12: same as option1, the only difference is use hj instead of hi
 * 13: This one is using similar idea of 11, but use [w(hi)+w(hj)]*0.5 instead of
 * Note: to make sure conservation of momentum and energy conservation, smoothing length might be changed so that sml for two phases are equal in mixing region
 * Any way, SPH will have some trouble if the smoothing length is different for two different phases.
 */
#ifndef DENSITY_UPDATE_SML
#define DENSITY_UPDATE_SML 1
#endif

//Define the smoothing length used in momentum and energy update
//Note: To be safe, it would be better to have more neighbors than needed by support.--->That's why I set CUTOFF2 = 2*CUTOFF
/*
 * 0 : use hi for updating of energy and momentum of particle i --> The default manner  		-->Only apply for SPH
 * 1 : use hj for updating of energy and momentum of particle i                         		-->Only apply for SPH
 * 2 : use 0.5(hi+hj) for updating of energy and momentum of particle i                 		-->apply for SPH and GSPH
 * 3 : use 0.5[w(hi)+w(hj)] for updating of energy and momentum of particle i           		-->apply for SPH and GSPH
 * 4 : use (w(hi)/rho_i^2+w(hj)/rho_j^2) --->This is actually the default formulation for GSPH  -->Only apply for GSPH
 */
#ifndef ME_UPDATE_SML
#define ME_UPDATE_SML 2
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

//This option only make sense when USE_SUMMATION==1
/*
 * There are two ways to view particles of two phases in multiphase SPH method.
 * In addition, there are several different formulations for updating density by weighted summation, For example, "A multi-phase SPH method for macroscopic and mesoscopic flows" by Xu
 * 1  : They are nothing but discretized points, phase num of each particle is just a flag of that particle (or properties of particle)
 * 10 : The same as option 1 update of density will not based on SPH, instead, it will based on an equation (Can be found in Suzuki's 2005 paper)
 * 11 : The same as option 1 except that does not do normalization.
 * 12 : The same as option 1, the only difference is that use phase density instead of total density in computing normalization.
 * 13 : The same as option 1, except that adopt the weighted summation formulation in the paper: "A multi-phase SPH method for macroscopic and mesoscopic flows". ---> And do not do normalization  ---> Probably DENSITY_UPDATE_SPH=12 should be used when use this form.
 * 2  : Two different sets of discretized points. essentially independent and only interact with each other by the interact terms (include explicit terms like drag force or implicit terms such as the pressure force term.) in the governing equation.
 * 22 : The same as option 2, the only difference is that use the total density instead of phase density in computing normalization.
 */
#ifndef DENSITY_UPDATE_SPH
#define DENSITY_UPDATE_SPH 11 //DENSITY_UPDATE_SPH=2, view particle as two different phases will lead to very large density
                              //When DENSITY_UPDATE_SPH=0, DENSITY_UPDATE_SML should be set to 1;
#endif

//Define whether should pressure ghost particles should be took into account for density update
//PGHOST_CONTRIBUTE_DES == 2: yes, take the eruption ghost particles and pressure ghost particles into account
//PGHOST_CONTRIBUTE_DES == 1: yes, take the pressure ghost particles into account
//PGHOST_CONTRIBUTE_DES == 0: No,  do not take the pressure ghost particles into account
#ifndef PGHOST_CONTRIBUTE_DES
#define PGHOST_CONTRIBUTE_DES 2
#endif

#endif // end of USE_SUMMATION
//-----------------------------------------------------------------------------------------------------------------------

//Different cases for shock tube test:
/*
 * 0: Sod shock in first GSPH
 * 1: another Sod shock test: see paper: assessment of localized artificial diffusive scheme for large-eddy simulation of compressible turbulent flow.
 * 2: shu-Osher problem, see paper: assessment of localized artificial diffusive scheme for large-eddy simulation of compressible turbulent flow.
 * 3: Sjogreen test: See paper Approximate Riemann Solvers for the Godunov SPH
 * 4: strong blast test: See Toro's book "Riemann Solvers and numerical method for fluid dynamics"
 * 5: Double expansion
 * 6: Double shock
 */
#ifndef SHOCK_TUBE_TESTS
#define SHOCK_TUBE_TESTS 4
#endif

//Whether smooth the initial distribution for Shock tube problem or not?
/*
 * 0: No
 * 1: Yes
 */
#ifndef SHOCK_TUBE_SMOOTH
#define SHOCK_TUBE_SMOOTH 0
#endif

//if defined, use equal particle mass and different sml
//Otherwise, use different particle mass to guarantee equal sml
//This option is only for 1D problem
/*
 * EQUAL_PART_MASS 1: yes equal mass
 * EQUAL_PART_MASS 0: No, different mass
 */
#if (SHOCK_TUBE_TESTS==0) || (SHOCK_TUBE_TESTS==1)
#ifndef EQUAL_PART_MASS
#define EQUAL_PART_MASS 0
#endif
#endif

//define to use GSPH or not ---> I need move this option to configure after test is done
/*
 * USE_GSPH 0 : SPH
 * USE_GSPH 1 : GSPH
 * USE_GSPH 2 : RCMSPH
 */
#ifndef USE_GSPH
#define USE_GSPH 2
#endif


//define the order that used to approximate specific volume between two particles for GSPH
/* GSPH_SPECIFIC_VOL_APP 0: zeroth order --> means use fi as fr and fj as fl directly --> It is actually piece wise constant
 * GSPH_SPECIFIC_VOL_APP 1: linear --> See Shu-ichiro's 2002 paper  --> It is actually piece wise linear
 * GSPH_SPECIFIC_VOL_APP 3: Cubic --> See Shu-ichiro's 2002 paper   --> It is actually piece wise cubic
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef GSPH_SPECIFIC_VOL_APP
#define GSPH_SPECIFIC_VOL_APP 0
#endif
#endif

//define wether to apply the minimum derivative condition and how to apply it
/*
 * MINI_DERIVATIVE_COND 0: Do not apply
 * MINI_DERIVATIVE_COND 1: apply to variables one by one and side by side, that is to say, if negative pressure on left side shows up, only apply the minimum derivative condition for pressure on the left side.
 * MINI_DERIVATIVE_COND 2: apply to variables one by one for both side, that is to say, if negative pressure on left side shows up, only apply the minimum derivative condition for pressure on both sides.
 * MINI_DERIVATIVE_COND 3: apply to all variables on both sides, that is to say, if negative pressure on one side shows up, apply the minimum derivative condition for all variables (pressure, density and velocity) on both sides.
 * MINI_DERIVATIVE_COND 10: based on 1, use piece-wise linear constant if C_SHOCK*abs(ul-ur)>min(CSi, CSj)
 * MINI_DERIVATIVE_COND 4: without check physics value, directly use zeroth order construction of RP
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef  MINI_DERIVATIVE_COND
#define  MINI_DERIVATIVE_COND 4
#endif
#endif

//define whether to use original monotonicity condition or modified monotonicity condition.
/* GSPH_MODIFIED_MONOTONICITY 0: original --> See Shu-ichiro's 2002 paper
 * GSPH_MODIFIED_MONOTONICITY 1: modified --> Based on Shu-ichiro's 2002 paper, the second monotonicity condition is change be to use absolute value of (ul-ur) There is no abs in Shu-ichiro's 2002 paper
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef GSPH_MODIFIED_MONOTONICITY
#define GSPH_MODIFIED_MONOTONICITY 0
#endif
#endif

//define which kind of Riemann Solver
/* 0: --> Roe
 * 1: --> HLLC
 * #2: --> HLLC , based on paper:A robust HLLC-type Riemann solver for strong shock  --> No difference between this one and previous one.
 * 3: --> An iterative solver based by Van Leer
 */
//For RCM SPH, only use HLLC Riemann Solver or HLLC Riemann Solver.
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef RIEMANN_SOLVER
#define RIEMANN_SOLVER 0
#endif
#endif

//define Sl and Sr evaluation in HLL and HLLC Riemann Solvers
/* 0: --> Kunal Puri, Approximate Riemann solvers for the Godunov SPH (GSPH) For HLLC Riemann Solver in Puri's paper --> should always use this one.
 * ###1--> B. Einfeldt On Godunov-type methods near low density ---> Are actually the same as 0, Notice that 0 is actually the for a moving coordinate system.
 * 2: --> S.F. Davis Simplified second-order Godunov-type methods
 * 3: --> E.F. Toro Riemann Solvers and Numerical Methods for Fluid Dynamics
 * ###4: --> P. Batten Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows ---> The same as 0
 * ###5: --> A modified version of 0 : us ul and ur instead of vl (=ul-ulr) and vr  ---> Actually 1, but incorrect if directly use, should be converted to a moving coordinate system which moves together with interface.
 *
 */
#if (RIEMANN_SOLVER==1) || (RIEMANN_SOLVER==2) //Use HLL type of Riemann Solver
#ifndef HLL_WAVE_SPEED_EVA
#define HLL_WAVE_SPEED_EVA 0
#endif
#endif

//Defines the way to compute derivative. ---> This option will be needed only when USE_GSPH = 1
/*
 * NORM_DERIVATIVE 0: derivative is not normalized  --> Currently not available
 * NORM_DERIVATIVE 1: derivative is normalized
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef NORM_DERIVATIVE
#define NORM_DERIVATIVE 1
#endif
#endif

//Defines the way to approximate shear velocity in GSPH. ---> This option will be needed only when USE_GSPH == 1
//It is a common practice to approximate the shear velocity to use (weighted) average of shear velocity on both sides.
//---> It has been shown that distance based average introduces more dissipation than Roe average (Which is essentially a square root of density weighted average).
//---> The finally solution to this problem is solving shear velocity wave from Riemann Solver.

/*
 * SHEAR_VEL_APP 0: arithmatical mean
 * SHEAR_VEL_APP 1: distanced weighted average. --->The same as what did in Shu-ichiro Inutsuka's paper.
 * SHEAR_VEL_APP 2: Roe average ---> Weighted by square root of density
 * SHEAR_VEL_APP 10: From solving a Riemann problem taking the shear velocity wave into account.
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef SHEAR_VEL_APP
#define SHEAR_VEL_APP 2
#endif
#endif

//The GSPH will also introduce artificial viscosity for expansion ---> To turn this off, one possible way is to use an if-else condition similar to what used in SPH
/*
 * SWITCH_OFF_AV_FOR_EXPAN 0: Do not switch off artificial viscosity for expansion
 * SWITCH_OFF_AV_FOR_EXPAN 1: switch off artificial viscosity for expansion
 * SWITCH_OFF_AV_FOR_EXPAN 2: use an switch coefficient that is proportional to the density gradient to determine the amount of artificial viscosity: Use the approximate numerical density gradient (AROUND 25) at the location where there is a shock as the reference to determine the value of this switch coefficient
 */
#if  (USE_GSPH==1 || USE_GSPH==2)
#ifndef SWITCH_OFF_AV_FOR_EXPAN
#define SWITCH_OFF_AV_FOR_EXPAN 0
#endif
#endif

//For particle mass with Non-uniform particle mass, use particle mass weighted average when solving Riemann Problems.
/*
 * RP_MASS_WEIGHTED==0: do not use mass weighted
 * RP_MASS_WEIGHTED==1: use m/rho as weight
 * RP_MASS_WEIGHTED==2: use m as weight  ---> Use mass weighted project only when project back
 * RP_MASS_WEIGHTED==21:use m as weight  ---> Use mass weighted project for both construct interface and project back.
 * RP_MASS_WEIGHTED==3: replace Roe average, which is weighted by sqrt(rho) with sqrt(m/rho)
 */
#if  (USE_GSPH==1 || USE_GSPH==2) && EQUAL_PART_MASS==0
#ifndef RP_MASS_WEIGHTED
#define RP_MASS_WEIGHTED 0
#endif
#endif

//Define whether use nature boundary condition for ks or essentiall boundary for ks
/*
 * 0: essential
 * 1: nature
 */
#ifndef BC_FOR_KX
#define BC_FOR_KX 1
#endif


//Define whether use adaptive smoothing length or not ---> adaptively adjust sml at a given interval, to avoid sml update at every time step
/*
 * 0: Not adaptive  -->DENSITY_UPDATE_SML an be 1 or 0
 * 1: adaptive scheme one: based on the assumption that smoothing length should be proportional to ratio between particle mass and density. This is the traditional way of adaptive smoothing length -->DENSITY_UPDATE_SML should always be 1
 * 11: based on 1, set an up and down limit on smoothing length ---> Due to current data structure, the current code does not support smoothing length of h > (ADDING_NUM/CUTOFF)*dx ---> currently, set the up limit as 2 times of dx, and the lower bound as the 2/3 of dx
 * 2: adaptive scheme two: based ADKE: see paper "Adaptive kernel estimation and SPH tensile instability", In this case, need original sml, even though updating of density is based on current sml, we still have original sml.
 * 3: adaptive scheme three: my own scheme
 * 31: adaptive scheme three: my own scheme, use the idea of mass gradient to adaptively change smoothing length
 * 32: based on 31, set an up and down limit on smoothing length ---> Due to current data structure, the current code does not support smoothing length of h > (ADDING_NUM/CUTOFF)*dx ---> currently, set the up limit as 2 times of dx, and the lower bound as the 2/3 of dx
 */
#ifndef ADAPTIVE_SML
#define ADAPTIVE_SML 32
#endif

//Whether apply adaptive smooth length to ghost particle or not
/*
 * 0: No
 * 1: Yes  ---> Currently not able to apply, need to 1) let ghost particles have neighbour information
 *                                                   2) Will have difficult to deal with boundary deficiency. ---> Normalized it? Not sure
 * 11: Instead of adaptively change smoothing length of ghost particles, set initial sml satisfy the common requirement on sml.
 */
#ifndef ADAPTIVE_SML_GHOST
#define ADAPTIVE_SML_GHOST 11
#endif

//Define which format of to use for discretized momentum equation
/*
 * 0: The basic symmetric format which can conserve momentum, without any further modification
 * 1: Based on 0, an external pressure is deduce by every pressure, the external pressure the pressure of atmosphere at corresponding height of particle a. ---> The purpose of this is to make sure when pressure gradient vanish the acceleration will be zero.
 *    Note: Please note that for GSPH, such trick is not necessary, as GSPH can guarantee zero RHS of momentum equation when pressure gradient is zero.
 *    So you should not use this option together with GSPH
 */
#ifndef MOMENTUM_DISCRETIZE
#define MOMENTUM_DISCRETIZE 0
#endif

//Define have LANS turbulent model in the code
//--> it is a stupid idea to do module management in C++ in this way, I should make use of the template, inherit, overloading as much as possible
/*
 * 0 : No HAVE_TURBULENCE_LANS
 * 1 : only filter velocity
 * 2 : filter both velocity and energy  ---> In which case, it is not necessary to smooth energy. ---> For energy smooth, it is OK to use a different filter scale length.
 */
#ifndef HAVE_TURBULENCE_LANS
#define HAVE_TURBULENCE_LANS 0
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
#define ATMOSPHERE_TYPE 2
//The default value represents hydro-static atmosphere
#endif

//Define whether use variable gravity or not, the gravity is a function of height:
//   g = 9.80665 * (6400/(6400+h(km)))^2
/*
 * 0: use constant gravity
 * 1: use
 */
#if (ATMOSPHERE_TYPE==0) || (ATMOSPHERE_TYPE==1) || (ATMOSPHERE_TYPE==4)
#ifndef VARIABLE_GRAVITY
#define VARIABLE_GRAVITY 0
#endif
#endif


//Define the erupt velocity profile type --> for plume modeling
/*
 * 0: uniform
 * 1: parabolic
 */
#ifndef ERUPT_VELOCITY_PROF
#define ERUPT_VELOCITY_PROF 0 //The default value represents hydro-static atmosphere
#endif

//Turn on and turn off domain adjusting ---> might be useful for performance analysis
// if define ADJUST_DOMAIN, will enable the feature of domain adjusting
// if domain adjusting feature is turned off, the initial domain (parameters giving in input file) should be consistent with computational domain
//#ifndef ADJUST_DOMAIN
//#define ADJUST_DOMAIN
//#endif

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
/*
 * Not defined: do not output
 * If defined, output to h5part file for sure
 * If  WRITE_GHOSTS==2: besides output to h5part file, also output to csv file  ---> Only apply to 1D code.
 */
#ifndef WRITE_GHOSTS
#define WRITE_GHOSTS 0
#endif

//output PID
#ifndef WRITE_PID
#define WRITE_PID
#endif

//output particle mass ---> particle mass is necessary in post process
#ifndef WRITE_PMASS
#define WRITE_PMASS
#endif
//
////output particle smoothing length ---> useful for debug
#ifndef WRITE_SML
#define WRITE_SML
#endif

#endif //DEBUG

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//##################################################################################################################################################################

#else //CODE_DIMENSION==2 or 3
////using Gaussian Kernel currently, only one kind of kernel is available
#ifndef USE_GAUSSIAN
#define USE_GAUSSIAN
#endif

//define use summation method to update density  ---> Another way is based on discretized mass conservation
/*
 * USE_SUMMATION=1: Use weighted summation formulation for density update
 * USE_SUMMATION=0: Use discretized mass conservation equation ---> Currently, did not consider updating for mass fraction yet. ---> Is this one works well with shock tube, try it for Plume model
 */
#ifndef USE_SUMMATION
#define USE_SUMMATION 1


// only when no sml adaptive is used, it will be necessary to decide to use either original sml or current sml. If sml is adaptive, always use current sml
// NOTE: currently, all of these schemes are only applied to density updating. Actually, some of these ideas can be used in momentum and energy updating.
/*
 * Define the smoothing length that will be used in density update: ---> only need when SPH summation formulism is used!
 * 0 : use original smoothing length   ---> It is easier to get negative energy
 * 1 : use current smoothing length hi  ----> More stable choice
 * 11: same as option1, the only difference is use (hi+hj)*0.5 instead of hi
 * 12: same as option1, the only difference is use hj instead of hi
 * 13: This one is using similar idea of 11, but use [w(hi)+w(hj)]*0.5 instead of
 * Note: to make sure conservation of momentum and energy conservation, smoothing length might be changed so that sml for two phases are equal in mixing region
 * Any way, SPH will have some trouble if the smoothing length is different for two different phases.
 */
#ifndef DENSITY_UPDATE_SML
#define DENSITY_UPDATE_SML 1
#endif


//Define the smoothing length used in momentum and energy update
//Note: To be safe, it would be better to have more neighbors than needed by support.--->That's why I set CUTOFF2 = 2*CUTOFF
/*
 * 0 : use hi for updating of energy and momentum of particle i --> The default manner                  -->Only apply for SPH
 * 1 : use hj for updating of energy and momentum of particle i                                         -->Only apply for SPH
 * 2 : use 0.5(hi+hj) for updating of energy and momentum of particle i                                 -->apply for SPH and GSPH
 * 3 : use 0.5[w(hi)+w(hj)] for updating of energy and momentum of particle i                           -->apply for SPH and GSPH
 * 4 : use (w(hi)/rho_i^2+w(hj)/rho_j^2) --->This is actually the default formulation for GSPH  -->Only apply for GSPH
 */
#ifndef ME_UPDATE_SML
#define ME_UPDATE_SML 0
#endif


/*
 * Based on assumption of immediate thermodynamics equilibrium, a internal energy smooth might be necessary to make this assumption to be true  ---> is necessary when Turbulence model is included. ---> If the turbulence has an energy smoothing process, it is not necessary to do smooth again here!
 * 0: do not use energy smooth
 * 1: use energy smooth with normalization
 * 2: use energy smooth without normalization
 *
 * Note: this option is not available for GSPH yet.
 */
#ifndef HAVE_ENERGY_SMOOTH
#define HAVE_ENERGY_SMOOTH 1
#endif

/*
 * There are two ways to view particles of two phases in multiphase SPH method.
 * In addition, there are several different formulations for updating density by weighted summation, For example, "A multi-phase SPH method for macroscopic and mesoscopic flows" by Xu
 * 1  : They are nothing but discretized points, phase num of each particle is just a flag of that particle (or properties of particle)
 * 10 : The same as option 1 update of density will not based on SPH, instead, it will based on an equation (Can be found in Suzuki's 2005 paper)
 * 11 : The same as option 1 , except that does not do normalization.
 * 12 : The same as option 1, the only difference is that use phase density instead of total density in computing normalization.
 * 13 : The same as option 1, except that adopt the weighted summation formulation in the paper: "A multi-phase SPH method for macroscopic and mesoscopic flows". ---> And do not do normalization
 * 2  : Two different sets of discretized points. essentially independent and only interact with each other by the interact terms (include explicit terms like drag force or implicit terms such as the pressure force term.) in the governing equation.
 * 22 : The same as option 2, the only difference is that use the total density instead of phase density in computing normalization.
 */
#ifndef DENSITY_UPDATE_SPH
#define DENSITY_UPDATE_SPH 11 //DENSITY_UPDATE_SPH=2, view particle as two different phases will lead to very large density
                           //When DENSITY_UPDATE_SPH=0, DENSITY_UPDATE_SML should be set to 1;
#endif

//Define whether should pressure ghost particles should be took into account for density update
//PGHOST_CONTRIBUTE_DES == 2: yes, take the eruption ghost particles and pressure ghost particles into account
//PGHOST_CONTRIBUTE_DES == 1: yes, take the pressure ghost particles into account
//PGHOST_CONTRIBUTE_DES == 0: No,  do not take the pressure ghost particles into account
#ifndef PGHOST_CONTRIBUTE_DES
#define PGHOST_CONTRIBUTE_DES 1
#endif

#endif // end of USE_SUMMATION


//define to use GSPH or not ---> I need move this option to configure after test is done
/*
 * USE_GSPH 0 : SPH
 * USE_GSPH 1 : GSPH
 * USE_GSPH 2 : RCMSPH
 */

#ifndef USE_GSPH
#define USE_GSPH 0
#endif


//define the order that used to approximate specific volume between two particles for GSPH
/* GSPH_SPECIFIC_VOL_APP 0: zeroth order --> means use fi as fr and fj as fl directly
 * GSPH_SPECIFIC_VOL_APP 1: linear --> See Shu-ichiro's 2002 paper
 * GSPH_SPECIFIC_VOL_APP 3: Cubic --> See Shu-ichiro's 2002 paper
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef GSPH_SPECIFIC_VOL_APP
#define  GSPH_SPECIFIC_VOL_APP 3
#endif
#endif

//define wether to apply the minimum derivative condition and how to apply it
/*
 * MINI_DERIVATIVE_COND 0: Do not apply
 * MINI_DERIVATIVE_COND 1: apply to variables one by one and side by side, that is to say, if negative pressure on left side shows up, only apply the minimum derivative condition for pressure on the left side.
 * MINI_DERIVATIVE_COND 2: apply to variables one by one for both side, that is to say, if negative pressure on left side shows up, only apply the minimum derivative condition for pressure on both sides.
 * MINI_DERIVATIVE_COND 3: apply to all variables on both sides, that is to say, if negative pressure on one side shows up, apply the minimum derivative condition for all variables (pressure, density and velocity) on both sides.
 * MINI_DERIVATIVE_COND 10: based on 1, use piece-wise linear constant if C_SHOCK*abs(ul-ur)>min(CSi, CSj)
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef  MINI_DERIVATIVE_COND
#define  MINI_DERIVATIVE_COND 10
#endif
#endif

//define whether to use original monotonicity condition or modified monotonicity condition.
/* GSPH_MODIFIED_MONOTONICITY 0: original --> See Shu-ichiro's 2002 paper
 * GSPH_MODIFIED_MONOTONICITY 1: modified --> Based on Shu-ichiro's 2002 paper, the second monotonicity condition is change be to use absolute value of (ul-ur) There is no abs in Shu-ichiro's 2002 paper
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef GSPH_MODIFIED_MONOTONICITY
#define GSPH_MODIFIED_MONOTONICITY 1
#endif
#endif

//define which kind of Riemann Solver
/* 0: --> Roe
 * 1: --> HLLC
 * #2: --> HLLC , based on paper:A robust HLLC-type Riemann solver for strong shock  --> No difference between this one and previous one.
 */
//For RCM SPH, only use HLLC Riemann Solver or HLLC Riemann Solver.
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef RIEMANN_SOLVER
#define RIEMANN_SOLVER 0
#endif
#endif

//define Sl and Sr evaluation in HLL and HLLC Riemann Solvers
/* 0: --> Kunal Puri, Approximate Riemann solvers for the Godunov SPH (GSPH) For HLLC Riemann Solver in Puri's paper --> should always use this one.
 * ###1--> B. Einfeldt On Godunov-type methods near low density ---> Are actually the same as 0, Notice that 0 is actually the for a moving coordinate system.
 * 2: --> S.F. Davis Simplified second-order Godunov-type methods
 * 3: --> E.F. Toro Riemann Solvers and Numerical Methods for Fluid Dynamics
 * ###4: --> P. Batten Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows ---> The same as 0
 * ###5: --> A modified version of 0 : us ul and ur instead of vl (=ul-ulr) and vr  ---> Actually 1, but incorrect if directly use, should be converted to a moving coordinate system which moves together with interface.
 *
 */
#if (RIEMANN_SOLVER==1) || (RIEMANN_SOLVER==2) //Use HLL type of Riemann Solver
#ifndef HLL_WAVE_SPEED_EVA
#define HLL_WAVE_SPEED_EVA 0
#endif
#endif

//Defines the way to compute derivative. ---> This option will be needed only when USE_GSPH = 1
/*
 * NORM_DERIVATIVE 0: derivative is not normalized
 * NORM_DERIVATIVE 1: derivative is normalized
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef NORM_DERIVATIVE
#define NORM_DERIVATIVE 1
#endif
#endif

//Defines the way to approximate shear velocity in GSPH. ---> This option will be needed only when USE_GSPH == 1
//It is a common practice to approximate the shear velocity to use (weighted) average of shear velocity on both sides.
//---> It has been shown that distance based average introduces more dissipation than Roe average (Which is essentially a square root of density weighted average).
//---> The finally solution to this problem is solving shear velocity wave from Riemann Solver.

/*
 * SHEAR_VEL_APP 0: arithmatical mean
 * SHEAR_VEL_APP 1: distanced weighted average. --->The same as what did in Shu-ichiro Inutsuka's paper.
 * SHEAR_VEL_APP 2: Roe average ---> Weighted by square root of density
 * SHEAR_VEL_APP 10: From solving a Riemann problem taking the shear velocity wave into account.
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef SHEAR_VEL_APP
#define SHEAR_VEL_APP 2
#endif
#endif

//The GSPH will also introduce artificial viscosity for expansion ---> To turn this off, one possible way is to use an if-else condition similar to what used in SPH
/*
 * SWITCH_OFF_AV_FOR_EXPAN 0: Do not switch off artificial viscosity for expansion
 * SWITCH_OFF_AV_FOR_EXPAN 1: switch off artificial viscosity for expansion
 * SWITCH_OFF_AV_FOR_EXPAN 2: use an switch coefficient that is proportional to the density gradient to determine the amount of artificial viscosity: Use the approximate numerical density gradient (AROUND 25) at the location where there is a shock as the reference to determine the value of this switch coefficient
 */
#if (USE_GSPH==1 || USE_GSPH==2)
#ifndef SWITCH_OFF_AV_FOR_EXPAN
#define SWITCH_OFF_AV_FOR_EXPAN 0
#endif
#endif

//For particle mass with Non-uniform particle mass, use particle mass weighted average when solving Riemann Problems.
/*
 * RP_MASS_WEIGHTED==0: do not use mass weighted
 * RP_MASS_WEIGHTED==1: use m/rho as weight
 * RP_MASS_WEIGHTED==2: use m as weight  ---> Use mass weighted project only when project back
 * RP_MASS_WEIGHTED==21:use m as weight  ---> Use mass weighted project for both construct interface and project back.
 * RP_MASS_WEIGHTED==3: replace Roe average, which is weighted by sqrt(rho) with sqrt(m/rho)
 */
#if  (USE_GSPH==1 || USE_GSPH==2) && EQUAL_PART_MASS==0
#ifndef RP_MASS_WEIGHTED
#define RP_MASS_WEIGHTED 0
#endif
#endif


//Define whether use nature boundary condition for ks or essentiall boundary ks
/*
 * 0: essential
 * 1: nature
 */
#ifndef BC_FOR_KX
#define BC_FOR_KX 1
#endif


//Define whether use adaptive smoothing length or not ---> adaptively adjust sml at a given interval, to avoid sml update at every time step
/*
 * 0: Not adaptive  -->DENSITY_UPDATE_SML an be 1 or 0
 * 1: adaptive scheme one: based on the assumption that smoothing length should be proportional to ratio between particle mass and density. This is the traditional way of adaptive smoothing length -->DENSITY_UPDATE_SML should always be 1
 * 11: based on 1, set an up and down limit on smoothing length ---> Due to current data structure, the current code does not support smoothing length of h > (ADDING_NUM/CUTOFF)*dx ---> currently, set the up limit as 2 times of dx, and the lower bound as the 2/3 of dx
 * 2: adaptive scheme two: based ADKE: see paper "Adaptive kernel estimation and SPH tensile instability", In this case, need original sml, even though updating of density is based on current sml, we still have original sml.
 * 3: adaptive scheme three: my own scheme
 * 31: adaptive scheme three: my own scheme, use the idea of mass gradient to adaptively change smoothing length
 * 32: based on 31, set an up and down limit on smoothing length ---> Due to current data structure, the current code does not support smoothing length of h > (ADDING_NUM/CUTOFF)*dx ---> currently, set the up limit as 2 times of dx, and the lower bound as the 2/3 of dx
 */
#ifndef ADAPTIVE_SML
#define ADAPTIVE_SML 0
#endif


//Define which format of to use for discretized momentum equation
/*
 * 0: The basic symmetric format which can conserve momentum, without any further modification
 * 1: Based on 0, an external pressure is deduce by every pressure, the external pressure the pressure of atmosphere at corresponding height of particle a. ---> The purpose of this is to make sure when pressure gradient vanish the acceleration will be zero.
 *    Note: Please note that for GSPH, such tricky is not necessary, as GSPH can guarantee zero RHS of momentum equation when pressure gradient is zero.
 *    So you should not use this option together with GSPH
 */
#ifndef MOMENTUM_DISCRETIZE
#define MOMENTUM_DISCRETIZE 0
#endif

//Define have LANS turbulent model in the code
//--> it is a stupid idea to do module management in C++ in this way, I should make use of the template, inherit, overloading as much as possible
/*
 * 0 : No HAVE_TURBULENCE_LANS
 * 1 : only filter velocity
 * 2 : filter both velocity and energy  ---> In which case, it is not necessary to smooth energy. ---> For energy smooth, it is OK to use a different filter scale length.
 */
#ifndef HAVE_TURBULENCE_LANS
#define HAVE_TURBULENCE_LANS 1
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

//Define whether use variable gravity or not, the gravity is a function of height:
//   g = 9.80665 * (6400/(6400+h(km)))^2
/*
 * 0: use constant gravity
 * 1: use variable gravity as a function of height
 */
#if (ATMOSPHERE_TYPE==0) || (ATMOSPHERE_TYPE==1) || (ATMOSPHERE_TYPE==4)
#ifndef VARIABLE_GRAVITY
#define VARIABLE_GRAVITY 0
#endif
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
 * To avoid creating new energy, we actually convert kinetic energy into internal energy to avoid negative pressure (or in other words negative internal energy)
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
//
////output particle smoothing length ---> useful for debug
//#ifndef WRITE_SML
//#define WRITE_SML
//#endif

#endif //DEBUG
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//##################################################################################################################################################################
#endif //CODE_DIMENSION

#endif /* OPTIONS_H */
