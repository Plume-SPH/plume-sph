/*
 * constant.h
 *
 *  Created on: Feb 16, 2015
 *      Author: zhixuanc
 */

#ifndef CONSTANT_H
#define CONSTANT_H
#include "options.h"

#if CODE_DIMENSION==1

const int DIMENSION=1;
const int PHASE_NUM=2;
const int KEYLENGTH=2;
const int TKEYLENGTH=KEYLENGTH+1;
const int NO_OF_EQNS=DIMENSION+2;/*For my current model which assumed immediate thermal and dynamic equilibrium!*/
const int NOEQxDIM=NO_OF_EQNS*DIMENSION; //total number of gradients
const int DIMSQRD=1;
const int NEIGH_SIZE=3;
const int MAX_PARTICLES_PER_BUCKET=6000; // usually 6000 for JPUE, 4500 for plume modelling
//const int MAX_ADD_STEPS = 1000; //maximum number of steps for eruption material adding

// Directions
const int XDIR=0;
const int YDIR=1;
const int ZDIR=2;

// Boundary Conditions
const int NONE    =0;
const int DIRCHLET=1;
const int NEUMANN =2;

// Ghost particles
const double EXT_DOM_COF = 1.5;//Coefficient that used to extend the domain to have some buckets for ghost particles.
const double EXT_DOM_COF_BOT = 1.5;//Coefficient that used to extend  -z direction of the domain to have some buckets for ghost particles.
                                   //This number can be different from EXT_DOM_COF, usually take 1.6 for PARTICLE_DENSITY=5
// Number of particles per cell per dimension EXT_DOM_COF, to avoid placing particles on the boundary-->placing of particles on the boundary will cause instability of simulation.
const int PARTICLE_DENSITY=6;
const int ADDING_NUM= 6; //For 1D is 6
const double CUTOFF = 3.0; // CUTOFF is very critical
const double CUTOFF2 = 12.0; // CUTOFF is very critical --> Need to be at least CUTOFF*sqrt(2)
const double CUTOFF3 = 12.0;

// Bucket TYPES
const int BREIF = 0x0;  //Brief bucket which is totally empty
const int UNDERGROUND = 0xA;  //10
const int MIXED       = 0xB;  //11  -->MIXED is important, it contains boundary info
const int OVERGROUND  = 0xC;  //12
const int PRESS_BC    = 0xD;  //13

const int INVOLVED = 2; //involved for particle
const int POTENTIAL_INVOLVED = 1; //potential involved for particle
const int NON_INVOLVED = 0; //not involved for particle

const int HAS_INVOLVED = 2; //has only involved for bucket
const int HAS_BOTH = 3; //has both involved and potential involved particles for bucket
const int HAS_POTENTIAL_INVOLVED = 1; // has only potential involved for bucket
const int HAS_NON_INVOLVED = 0; //has not involved for bucket

//Load balance
const float  LOAD_BALANCE_TOLERANCE = 1.001;
const double PI = 3.14159265358979;
const double TINY = 1.0E-08;


const double MSFRC_THRESH = 1.0E-06; //Threshold for determine whether air particle is involved in or not.

//For adding newly erupted/influx particles, that hash table is a temporiry hash table
const int ERUPT_TABLE_SIZE = 40000;

//For turbulence modeling with SPH-epsilon method
const double EPSILON = 0.8;
const double EPSILON_HALF = 0.4; //For the efficiency of computation EPSILON_HALF = EPSILON/2
const double PRANDTL_NUM = 0.85; //Prandtl number

#if HAVE_ENERGY_CUT==1
 //Constant for cut energy method ---> a simple way to avoid negative pressure (or negative energy)
const double ENERGY_CUT = 60.0;
#endif

//heat transfer spatial ratio to momentum exchange
const double HEAT_TRANS_SCALE_RATIO = 1.0;
const double E_SMOOTH_RATIO = 4.0;  //Ratio of energy smoothing length to smoothing length of particle (which is used in velocity filtering for LANS turbulence model)
const double E_FILTER_HRATIO = 2.0;  //Ratio of energy filtering smoothing length to smoothing length of particle (which is used in velocity filtering for LANS turbulence model)

//Adaptive sml
const int SML_UPDATE_INT = 0.1;
const double RATIO_SML_DX = 1.0; //When use "conservation forms" Has to be consistent with "eta_smooth_P" --> Otherwise, it will cause numerical perturbation
								 //The ratio between sml and dx, the maximum value allowed by current code design is 2.0 (as number of particles in each direction in each bucket is 6, the kernal cut off is 3).
                                 /*
                                  * This constant is used within particle class (which might not be a good design),
                                  * 1) For whatever sml given in constructors of particle, the sml will be enlarged by a ratio of RATIO_SML_DX
                                  * 2) For whatever sml given to method "put_smlen (double h)", the sml will be enlarged by a ratio of RATIO_SML_DX
                                  *
                                  * By doing this, what we actually do is having enough particles within kernel without using adaptive sml (which will cause losing of conservation of energy and momentum)
                                  * This is especially necessary for simulation with relatively large density ratio.
                                  */
const double MASS_IND_THRESH=0.;

//Parameters for GSPH
const double C_SHOCK = 3.0; //corresponding to number of particles in shock surface.

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
#else

const int NUM_NODES=4;
const int DIMENSION=3;
const int PHASE_NUM=2;
const int KEYLENGTH=2;
const int TKEYLENGTH=KEYLENGTH+1;
const int NO_OF_EQNS=DIMENSION+2;/*For my current model which assumed immediate thermal and dynamic equilibrium!*/
const int NOEQxDIM=NO_OF_EQNS*DIMENSION; //total number of gradients
const int DIMSQRD=9;
const int NEIGH_SIZE=27;
const int MAX_PARTICLES_PER_BUCKET=6000; // usually 6000 for JPUE, 4500 for plume modelling
//const int MAX_ADD_STEPS = 1000; //maximum number of steps for eruption material adding

// Directions
const int XDIR=0;
const int YDIR=1;
const int ZDIR=2;

// Boundary Conditions
const int NONE    =0;
const int DIRCHLET=1;
const int NEUMANN =2;

// Ghost particles
const double EXT_DOM_COF = 1.5;//Coefficient that used to extend the domain to have some buckets for ghost particles.
const double EXT_DOM_COF_BOT = 1.5;//Coefficient that used to extend  -z direction of the domain to have some buckets for ghost particles.
                                   //This number can be different from EXT_DOM_COF, usually take 1.6 for PARTICLE_DENSITY=5
// Number of particles per cell per dimension EXT_DOM_COF, to avoid placing particles on the boundary-->placing of particles on the boundary will cause instability of simulation.
const int PARTICLE_DENSITY=6;
const int ADDING_NUM= 256; //For 3D is 256
const double CUTOFF = 3.0; // CUTOFF is very critical
const double CUTOFF2 = 6.0; // CUTOFF is very critical --> Need to be at least CUTOFF*sqrt(2)

// Bucket TYPES
const int BREIF = 0x0;  //Brief bucket which is totally empty
const int UNDERGROUND = 0xA;  //10
const int MIXED       = 0xB;  //11  -->MIXED is important, it contains boundary info
const int OVERGROUND  = 0xC;  //12
const int PRESS_BC    = 0xD;  //13

const int INVOLVED = 2; //involved for particle
const int POTENTIAL_INVOLVED = 1; //potential involved for particle
const int NON_INVOLVED = 0; //not involved for particle

const int HAS_INVOLVED = 2; //has only involved for bucket
const int HAS_BOTH = 3; //has both involved and potential involved particles for bucket
const int HAS_POTENTIAL_INVOLVED = 1; // has only potential involved for bucket
const int HAS_NON_INVOLVED = 0; //has not involved for bucket

//Load balance
const float  LOAD_BALANCE_TOLERANCE = 1.001;
const double PI = 3.14159265358979;
const double TINY = 1.0E-08;

//For domain adjusting
#ifdef SIMULATE_ASH
const double MSFRC_THRESH = 1.0E-08; //Threshold for determine whether air particle is involved in or not --> use smaller value for ash transportation.
#else
const double MSFRC_THRESH = 1.0E-06; //Threshold for determine whether air particle is involved in or not.
#endif

//For adding newly erupted/influx particles, that hash table is a temporiry hash table
const int ERUPT_TABLE_SIZE = 40000;

//For turbulence modeling with SPH-epsilon method
const double EPSILON = 0.8;
const double EPSILON_HALF = 0.4; //For the efficiency of computation EPSILON_HALF = EPSILON/2
const double PRANDTL_NUM = 0.85; //Prandtl number

#if HAVE_ENERGY_CUT==1
 //Constant for cut energy method ---> a simple way to avoid negative pressure (or negative energy)
const double ENERGY_CUT = 60.0;
#endif

//heat transfer spatial ratio to momentum exchange
const double HEAT_TRANS_SCALE_RATIO = 1.0;
const double E_SMOOTH_RATIO = 4.0;  //Ratio of energy smoothing length to smoothing length of particle (which is used in velocity filtering for LANS turbulence model)
const double E_FILTER_HRATIO = 2.0;  //Ratio of energy filtering smoothing length to smoothing length of particle (which is used in velocity filtering for LANS turbulence model)  --> will be used when HAVE_TURBULENCE_LANS == 2

//Adaptive sml
const int SML_UPDATE_INT = 3.;
const double RATIO_SML_DX = 1.0; //The ratio between sml and dx, the maximum value allowed by current code design is 2.0 (as number of particles in each direction in each bucket is 6, the kernal cut off is 3).
                                 /*
                                  * This constant is used within particle class (which might not be a good design),
                                  * 1) For whatever sml given in constructors of particle, the sml will be enlarged by a ratio of RATIO_SML_DX
                                  * 2) For whatever sml given to method "put_smlen (double h)", the sml will be enlarged by a ratio of RATIO_SML_DX
                                  *
                                  * By doing this, what we actually do is having enough particles within kernel without using adaptive sml (which will cause losing of conservation of energy and momentum)
                                  * This is especially necessary for simulation with relatively large density ratio.
                                  */

const double MASS_IND_THRESH=0.;

//Parameters for GSPH
const double C_SHOCK = 3.0; //corresponding to number of particles in shock surface.
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

#endif //CODE_DIMENSION

#endif /* CONSTANT_H_ */
