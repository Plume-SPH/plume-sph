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

// constructor for initial air
Particle::Particle (unsigned *keyin, double *crd, double m, double h, double prss, double masfrc, double gmm, double sndspd, int phs_num, int id, int bc)
{
  int i;

  mass = m;
  smlen = h;

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

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0;

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

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0;//velocity is set to zero

  state_vars[3] = Vv0;//velocity in z direction is set to its initial value

// The following is newly added
    state_vars[NO_OF_EQNS - 1]= ev0;
    phase_num= 2;
    mass_frac= 1.;

    double des, engr;
    des = rhov;
    state_vars[0] = des;
    engr = ev0;
    state_vars[NO_OF_EQNS-1] = engr;

    update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P);
    pressure = pv0;

//    gamma=gamma_v;
//    sound_speed=sndspd;

  return;

}
//constructor from no input parameter
Particle::Particle ()
{
  int i;

  for (i = 0; i < TKEYLENGTH; i++)
    key.key[i] = 0;

  mass = 0.;
  smlen = 0.;
  guest = false;
  reflection = false;
  bc_type = 100;
  involved = 0;

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = 0.;
  }

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0.;

}

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
	double mass_frac_in =  mass_frac;

	double gmm, press, sndspd;

	double ng=ng0_P*mass_frac_in;
	double ns=mass_frac_in-ng;
	double na=1-mass_frac_in;

	double Cvm=ns*Cvs_P+ng*Cvg_P+na*Cva_P;
	double Rm=ng*Rg_P+na*Ra_P;

	gmm=1+Rm/Cvm;

	double lmd=desm*Rm/Cvm;
	press=lmd*engr;

	// for gas mixture
	double des1=(1-mass_frac_in)*desm;//density of phase1 : air
	double des2=mass_frac_in*desm;    //density of phase2 : erupted material
	double desmg=des1*na+des2*ng;     //density of mixture of erupted gas and air
	double Cvmg=ng*Cvg_P+na*Cva_P;
	double engrg=engr*Cvmg/Cvm;

	sndspd=pow((Rm*(press/desmg+engrg)/Cvmg),0.5);  //sound speed is assumed to be only depends on gas phase.

	pressure = press ;
	sound_speed = sndspd ;
	gamma = gmm ;
}


