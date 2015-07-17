/*
 * constant.h
 *
 *  Created on: Feb 16, 2015
 *      Author: zhixuanc
 */

#ifndef CONSTANT_H
#define CONSTANT_H

const int NUM_NODES=4;
const int DIMENSION=3;
const int PHASE_NUM=2;
const int KEYLENGTH=2;
const int TKEYLENGTH=KEYLENGTH+1;
const int NO_OF_EQNS=DIMENSION+2;/*For my current model which assumed immediate thermal and dynamic equilibrium!*/
const int DIMSQRD=9;
const int NEIGH_SIZE=27;
const int MAX_PARTICLES_PER_BUCKET=3000;
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
//const int NUM_GHOST_ROWS=12;//Was increased from 6 to 12 so that a reliable static pressure BC can be imposed.
const double EXT_DOM_COF = 1.5;//Coefficient that used to extend the domain to have some buckets for ghost particles.
// Number of particles per cell per dimension
const int PARTICLE_DENSITY=6;

// Bucket TYPES
const int UNDERGROUND = 0xA;  //10
const int MIXED       = 0xB;  //11
const int OVERGROUND  = 0xC;  //12
const int PRESS_BC    = 0xD;  //13
const float  LOAD_BALANCE_TOLERANCE = 1.001;
const double PI = 3.14159265358979;
const double TINY = 1.0E-08;

//For adding newly erupted particles, that hash table is a temporiry hash table
const int ERUPT_TABLE_SIZE = 40000;
#endif /* CONSTANT_H_ */
