/*
 * particle.cc
 *
 *  Created on: Feb 17, 2015
 *      Author: zhixuanc
 */

#include <cmath>
#include <iostream>
using namespace std;

#include "particle.h"
#include "constant.h"

const double PHASE_DENS=1.0/PHASE_NUM; //This is just used for constructor of particle

//constructor from no input parameter--> the basic constructor
Particle::Particle ()
{
  int i;

  for (i = 0; i < TKEYLENGTH; i++)
    key.key[i] = 0;

  mass = 0.;
  smlen = 0.;

  specific_heat_p = 0.;

  update_delayed = false;
  guest = false;
  reflection = false;
  new_old = 0; //default new_old for non-guest particles is 0.
                 //for guest particles: -1: old, 1: new
  bc_type = 100;
  involved = 0;

  for (i = 0; i < DIMENSION; i++)
    coord[i] = 0.;

  for (i = 0; i < DIMENSION; i++)
	  smoothed_v[i] = 0.;

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0.;

  for (i = 0; i < PHASE_NUM; i++)
	  phase_density[i] = 0.;

  pressure = 0.;

  mass_frac = 0.;
  gamma = 0;
  sound_speed = 0.;
  temperature = 0.;
  phase_num = 1;
  myprocess = 10000; //The default process id

  return;
}

/*
 * No matter which kind of particle I am creating, the second variable should be a function of primitive variable
 * But in this constructor, secondary variable is given
 * This is not good, I should double check and try to make the sure that secondary variable computed is consistent with these given.
 *
 * So the following constructor should be replaced by some function in the future
 *
 * Well, because the density is determine later by SPH summation, so computing of secondary variables from primitive is not feasible!!!
 * Actually, energy, pressure, density and mass will be determined in initial_air(), secondary variable updating will also be called there.
 */
// constructor for initial air
Particle::Particle (unsigned *keyin, double *crd, double m, double h, double prss, double masfrc, double gmm, double sndspd, int phs_num, int id, int bc)
{
  int i;

  mass = m;
  smlen = h;

  specific_heat_p = 0.;

  update_delayed = false;
  guest = false;
  reflection = false; //The newly added wall ghost will not have image untill search, so need to be false.
  new_old = 0; //default new_old for non-guest particles is 0.
               //for guest particles: -1: old, 1: new
  bc_type = bc;
  involved = 0;

  for (i = 0; i < TKEYLENGTH; i++)
    key.key[i] = *(keyin + i);

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = *(crd + i);
  }

  for (i = 0; i < DIMENSION; i++)
	  smoothed_v[i] = 0.;

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0;
//
//  for (i = 0; i < PHASE_NUM; i++)
//  	  phase_density[i] = PHASE_DENS; //divided by PHASE_NUM to make sure sum of density of all phase will be 1.
//                                        //Density will be updated in updating of secondary variable

  state_vars[0] = 1.0;//default density is 1.0;

  pressure = prss;

  mass_frac = masfrc;
  gamma = gmm;
  sound_speed = sndspd;
  phase_num = phs_num;
  myprocess = id;

  return;
}

//overloading
//Constructor for adding eruption initial condition: all parameters come from input
Particle::Particle (unsigned *keyin, double *crd, double m, double h, int id,
		double Vv0, double ev0, double rhov, double pv0, double gamma_v,
		double ng0_P, double Cvs_P, double Cvg_P, double Cva_P, double Rg_P, double Ra_P )
{
  int i;
  myprocess = id;
  mass = m;
  smlen = h;

  specific_heat_p = 0.;

  update_delayed = false;
  guest = false;
  reflection = false;
  new_old = 0; //default new_old for non-guest particles is 0.
                 //for guest particles: -1: old, 1: new

  bc_type = 0;
  involved = 0;

  for (i = 0; i < TKEYLENGTH; i++)
    key.key[i] = *(keyin + i);

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = *(crd + i);
//    bedfrict[i] = 0;
  }

  for (i = 0; i < DIMENSION; i++)
	  smoothed_v[i] = 0.;

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0;//velocity is set to zero

  for (i = 0; i < PHASE_NUM; i++)
       phase_density[i] = PHASE_DENS; //divided by PHASE_NUM to make sure sum of density of all phase will be 1.
                                         //Density will be updated in updating of secondary variable

    state_vars[3] = Vv0;//velocity in z direction is set to its initial value

// The following is newly added
    state_vars[NO_OF_EQNS - 1]= ev0;
    phase_num= 2;
    mass_frac= 1.; //default is 1

    double des, engr;
    des = rhov;
    state_vars[0] = des;
    engr = ev0;
    state_vars[NO_OF_EQNS-1] = engr;

    update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P);
    pressure = pv0; //Why I have this? --> I need double check this!

//    gamma=gamma_v;
//    sound_speed=sndspd;

  return;

}

//overload the operator ==
bool Particle::operator== (const Particle & rhs) const
{
  for (int i = 0; i < TKEYLENGTH; i++)
    if (key.key[i] != (rhs.getKey ()).key[i])
      return false;
  return true;
}

//member function that used to update second variables based on give primitive variables.
void Particle::update_second_var(double ng0_P, double Cvs_P, double Cvg_P, double Cva_P, double Rg_P, double Ra_P)
{
	//desm should be density of mixture of erupt material and air
	double desm = state_vars[0];
	double engr = state_vars[NO_OF_EQNS-1];

	double ng=ng0_P*mass_frac;
	double ns=mass_frac-ng;
	double na=1-mass_frac;

	double Cvm=ns*Cvs_P+ng*Cvg_P+na*Cva_P;
	//In the old code, I am using EOS for idea gas for mixture of gas and solid
	//Now, I am using EOS for gas mixture
	double Rm=ng*Rg_P+na*Ra_P;
	double Rmg=Rm/(ng+na); //gas constant for gas mixture

	gamma=1+Rm/Cvm; // gamma for mixture

//	double lmd=desm*Rm/Cvm; not necessary
	// for gas mixture
//	//The following way of density computing is not correct -->The reason is na and ng is density with respect to mixture of solid and gas, In third equation, I was assuming that na is mass fraction with respect to phase1 and ng is with respect to pahse2
//	double des1=(1-mass_frac)*desm;//density of phase1 : air
//	double des2=mass_frac*desm;    //density of phase2 : erupted material
//	double desmg=des1*na+des2*ng;     //density of mixture of erupted gas and air
	//The correct way:
	// double desmg=des1+des2*ng/mass_frac;
   //The computational efficient way is:
	double desmg=desm*(na+ng);

	//implement EOS for idea gas on gas phase (mixture of air and gas)
	pressure = Rmg*desmg*engr/Cvm; //attention: engr/Cvm is temperature of mixture, also temperature of gas mixture

	//This is the old way of computing Cvmg --> The mistake that I made here is again using mass fraction with respect to mixture of gas and solid to compute property for gas mixture.
//	double Cvmg=ng*Cvg_P+na*Cva_P;

	//The correct way of computing Cvmg
	double Cvmg=(ng*Cvg_P+na*Cva_P)/(ng+na);

//	//This is the old way
//	double engrg=engr*Cvmg/Cvm; //specific internal energy for air+gas from erupted material
//	sound_speed=pow((Rm*(press/desmg+engrg)/Cvmg),0.5);  //sound speed is assumed to only depends on gas phase---> definitely not exact, but in SPH, sound speed is only useful when determine the time steps.
//	sound_speed=pow((Rmg*(pressure/desmg+engrg)/Cvmg),0.5);  //sound speed is assumed to only depends on gas phase---> definitely not exact, but in SPH, sound speed is only useful when determine the time steps.

//	//The new way
	sound_speed=pow((1+Rmg/Cvmg)*pressure/desmg, 0.5); //use sound speed of only gas phase, that will make sound speed larger than the real value, but safe in determine time steps.

//Updating single phase density
	phase_density[0]=state_vars[0]*na;
	phase_density[1]=state_vars[0]*mass_frac;
//update temperature
	temperature = engr/Cvm;

//update specific heat
#ifdef HAVE_TURBULENCE_LANS
	specific_heat_p=Cvm;
#endif
}


