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

#if DENSITY_UPDATE_SML==0 || ADAPTIVE_SML==2
  smlen_original = 0.;
#endif


  specific_heat_p = 0.;

  update_delayed = false;
  guest = false;
  reflection = false;
  new_old = 0; //default new_old for non-guest particles is 0.
                 //for guest particles: -1: old, 1: new
  bc_type = 100;
  involved = 0; //default value

  for (i = 0; i < DIMENSION; i++)
    coord[i] = 0.;

  for (i = 0; i < DIMENSION; i++)
	  smoothed_v[i] = 0.;

  smoothed_e =0.0;

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0.;

  for (i = 0; i < PHASE_NUM; i++)
	  phase_density[i] = 0.;

#if ADAPTIVE_SML==2
  rho_based_on_dx=0.0;
#endif

  pressure = 0.;

  mass_frac = 0.;
  new_mass_frac = 0.;
  gamma = 0;
  sound_speed = 0.;
  temperature = 0.;
  phase_num = 1;
  myprocess = 10000; //The default process id

#if (USE_GSPH==1 || USE_GSPH==2)//Assume 3D
  //derivatives
  for (i = 0; i < DIMENSION; i++)
	  d_rho[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_u[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_v[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_w[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_p[i]=0.0;
#endif

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
  smlen = RATIO_SML_DX* h;

#if DENSITY_UPDATE_SML==0
  smlen_original = h;
#elif ADAPTIVE_SML==2
  smlen_original = h*ADKE_sml_ratio_P;
#endif

  specific_heat_p = 0.;

  update_delayed = false;
  guest = false;
  reflection = false; //The newly added wall ghost will not have image untill search, so need to be false.
  new_old = 0; //default new_old for non-guest particles is 0.
               //for guest particles: -1: old, 1: new
  bc_type = bc;
  involved = 0; //The default involved for all particles are not involved. When adding air particle, the value will be reset to 1 (potential involved.)

  for (i = 0; i < TKEYLENGTH; i++)
    key.key[i] = *(keyin + i);

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = *(crd + i);
  }

  for (i = 0; i < DIMENSION; i++)
	  smoothed_v[i] = 0.;

  smoothed_e =0.0;

  for (i = 0; i < NO_OF_EQNS; i++)
      state_vars[i] = 0.;
//
//  for (i = 0; i < PHASE_NUM; i++)
//  	  phase_density[i] = PHASE_DENS; //divided by PHASE_NUM to make sure sum of density of all phase will be 1.
//                                        //Density will be updated in updating of secondary variable

  state_vars[0] = 1.0;//default density is 1.0;

#if ADAPTIVE_SML==2
  rho_based_on_dx=1.0;
#endif

  pressure = prss;

  mass_frac = masfrc;
  new_mass_frac = 0.;
  gamma = gmm;
  sound_speed = sndspd;
  phase_num = phs_num;
  myprocess = id;

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //derivatives
  for (i = 0; i < DIMENSION; i++)
	  d_rho[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_u[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_v[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_w[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_p[i]=0.0;
#endif

  return;
}

//overloading
//Constructor for adding eruption initial condition: all parameters come from input
Particle::Particle (unsigned *keyin, double *crd, double m, double h, int id,
		double Vv0, double ev0, double rhov, double pv0, double gamma_v, double sndspd,
		double ng0_P, double Cvs_P, double Cvg_P, double Cva_P, double Rg_P, double Ra_P )
{
  int i;
  myprocess = id;
  mass = m;
  smlen = RATIO_SML_DX*h;

#if DENSITY_UPDATE_SML==0
  smlen_original = h;
#elif ADAPTIVE_SML==2
  smlen_original = h*ADKE_sml_ratio_P;
#endif


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
	  smoothed_v[i] = 0.0;

  smoothed_e = ev0;

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

    double des;
    des = rhov;
    state_vars[0] = des;

#if ADAPTIVE_SML==2
  rho_based_on_dx=des;
#endif

    update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P);
    pressure = pv0; //Why I have this? --> I need double check this!

//    gamma=gamma_v;
    sound_speed=sndspd;

    gamma=gamma_v;

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //derivatives
  for (i = 0; i < DIMENSION; i++)
	  d_rho[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_u[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_v[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_w[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_p[i]=0.0;
#endif

  return;

}


//overloading
//Constructor for adding influx particles
Particle::Particle (unsigned *keyin, double *crd, double m, double h, int id, double msfc0,
		double Vx, double Vy, double ev0, double rhov, double pv0, double gamma_v, double sndspd,
		double ng0_P, double Cvs_P, double Cvg_P, double Cva_P, double Rg_P, double Ra_P )
{
  int i;
  myprocess = id;
  mass = m;
  smlen = RATIO_SML_DX*h;

#if DENSITY_UPDATE_SML==0
  smlen_original = h;
#elif ADAPTIVE_SML==2
  smlen_original = h*ADKE_sml_ratio_P;
#endif

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
  }

  for (i = 0; i < DIMENSION; i++)
	  smoothed_v[i] = 0.;

  smoothed_e = ev0;

  for (i = 0; i < NO_OF_EQNS; i++)
      state_vars[i] = 0;//velocity is set to zero

  for (i = 0; i < PHASE_NUM; i++)
       phase_density[i] = PHASE_DENS; //divided by PHASE_NUM to make sure sum of density of all phase will be 1.
                                         //Density will be updated in updating of secondary variable

    state_vars[1] = Vx;//velocity in x direction
    state_vars[2] = Vy;//velocity in y direction

// The following is newly added
    state_vars[NO_OF_EQNS - 1]= ev0;
    phase_num= 2;
    mass_frac= msfc0;

    double des;
    des = rhov;
    state_vars[0] = des;

#if ADAPTIVE_SML==2
    rho_based_on_dx=des;
#endif

    update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P);
    pressure = pv0; //Why I have this? --> I need double check this!

//    gamma=gamma_v;
    sound_speed=sndspd;

    gamma=gamma_v;

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //derivatives
  for (i = 0; i < DIMENSION; i++)
	  d_rho[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_u[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_v[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_w[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_p[i]=0.0;
#endif

  return;

}

#if CODE_DIMENSION==1
// constructor for initial condition of Shocktube problem
Particle::Particle (unsigned *keyin, double *crd, double m, double h , double des, double vel, double prss, double gmm, double sndspd, int bc, int invloved)
{
  int i;

  mass = m;
  smlen = RATIO_SML_DX* h;

#if DENSITY_UPDATE_SML==0
  smlen_original = h;
#elif ADAPTIVE_SML==2
  smlen_original = h*ADKE_sml_ratio_P;
#endif

  specific_heat_p = 0.;

  update_delayed = false;
  guest = false;
  reflection = false; //The newly added wall ghost will not have image untill search, so need to be false.
  new_old = 0; //default new_old for non-guest particles is 0.
               //for guest particles: -1: old, 1: new
  bc_type = bc;
  involved = invloved;

  for (i = 0; i < TKEYLENGTH; i++)
    key.key[i] = *(keyin + i);

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = *(crd + i);
  }

  for (i = 0; i < DIMENSION; i++)
	  smoothed_v[i] = 0.;

  smoothed_e = prss/((gmm-1)*des);
//
//  for (i = 0; i < PHASE_NUM; i++)
//  	  phase_density[i] = PHASE_DENS; //divided by PHASE_NUM to make sure sum of density of all phase will be 1.
//                                        //Density will be updated in updating of secondary variable

  state_vars[0] = des;//default density is 1.0;
  state_vars[1] = vel;
  state_vars[2] = prss/((gmm-1)*des);

#if ADAPTIVE_SML==2
  rho_based_on_dx=des;
#endif

  pressure = prss;

  mass_frac = 0.;
  new_mass_frac = 0.;
  gamma = gmm;
  sound_speed = sndspd;
  phase_num = 1;
  myprocess = 1;

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //derivatives
  for (i = 0; i < DIMENSION; i++)
	  d_rho[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_u[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_v[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_w[i]=0.0;

  for (i = 0; i < DIMENSION; i++)
	  d_p[i]=0.0;
#endif

  return;
}

#endif

//overload the operator ==
bool Particle::operator== (const Particle & rhs) const
{
  for (int i = 0; i < TKEYLENGTH; i++)
    if (key.key[i] != (rhs.getKey ()).key[i])
      return false;
  return true;
}

//member function that used to update second variables based on give primitive variables.
void Particle::update_second_var(double ng0_P, double Cvs_P, double Cvg_P, double Cva_P, double Rg_P, double Ra_P, double rho_0)
{
	//desm should be density of mixture of erupt material and air
	double desm = state_vars[0];

	double na=1-mass_frac;

#if HAVE_TURBULENCE_LANS==2
	double engr = smoothed_e;
#else
	double engr = state_vars[NO_OF_EQNS-1];
#endif

	double ng=ng0_P*mass_frac;
	double ns=mass_frac-ng;

	double Cvm=ns*Cvs_P+ng*Cvg_P+na*Cva_P;

	//implement EOS to get pressure
#if FLUID_COMPRESSIBILITY==0 //using EOS of ideal gas
	//In the old code, I am using EOS for idea gas for mixture of gas and solid
	//Now, I am using EOS for gas mixture
	double Rm=ng*Rg_P+na*Ra_P;
	double Rmg=Rm/(ng+na); //gas constant for gas mixture

	gamma=1+Rm/Cvm; // gamma for mixture
	//gamma=1.4;

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
	pressure = Rmg*desmg*engr/Cvm; //attention: engr/Cvm is temperature of mixture, also temperature of gas mixture
#elif FLUID_COMPRESSIBILITY==1
	double gama=7.0;
	double B;
	B=rho_0*sound_speed*sound_speed/gama;  //sound speed used here is sound speed from previous time step.
    pressure=B*(pow((desm/rho_0), gama)-1) + pa0_P; //rho_0 should be the density when pressure is zero, here, the pressure refers to pressure due to hydrostatic. Atmosphere pressure should be added
#endif

#if FLUID_COMPRESSIBILITY==0 //using EOS of ideal gas
//This is the old way of computing Cvmg --> The mistake that I made here is again using mass fraction with respect to mixture of gas and solid to compute property for gas mixture.
//	double Cvmg=ng*Cvg_P+na*Cva_P;

	//The correct way of computing Cvmg
	double Cvmg=(ng*Cvg_P+na*Cva_P)/(ng+na);

//	//This is the old way
//	double engrg=engr*Cvmg/Cvm; //specific internal energy for air+gas from erupted material
//	sound_speed=pow((Rm*(press/desmg+engrg)/Cvmg),0.5);  //sound speed is assumed to only depends on gas phase---> definitely not exact, but in SPH, sound speed is only useful when determine the time steps.
//	sound_speed=pow((Rmg*(pressure/desmg+engrg)/Cvmg),0.5);  //sound speed is assumed to only depends on gas phase---> definitely not exact, but in SPH, sound speed is only useful when determine the time steps.

    //The new way
	sound_speed=pow(gamma*pressure/desmg, 0.5); //use sound speed of only gas phase, that will make sound speed larger than the real value, but safe in determine time steps.
	//sound_speed=pow(1.4*pressure/desmg, 0.5); //use sound speed of only gas phase, that will make sound speed larger than the real value, but safe in determine time steps.
#elif FLUID_COMPRESSIBILITY==1
	sound_speed=pow(gama*(pressure+B-pa0_P)/desm, 0.5); //The equation of sound_speed comes from paper: "Historical Review of Real-Fluid Isentropio Flow Models" by D. A. Sullivan
#endif //FLUID_COMPRESSIBILIT

//Updating single phase density
	phase_density[0]=state_vars[0]*na;
	phase_density[1]=state_vars[0]*mass_frac;

//update temperature
	temperature = engr/Cvm;

//update specific heat
#if HAVE_TURBULENCE_LANS != 0
	specific_heat_p=Cvm;
#endif
}

//function that calculate the density of mixture base on equation in suzuki's 2005 paper
void Particle::calc_density_suzuki(double msf)
{
	double msf_air = 1 - msf;
	double Cp_eupt = msf*Cpv_P;
	double Cp_air = msf_air*Cpa_P;
    //update density_new  ----> density will be updated later
	new_state_vars[0] =  rhoa0_P*Ra_P/(msf*ng0_P*Rg_P+msf_air*Ra_P)*(Cp_eupt+Cp_air)/(Cp_eupt*Tv0_P+Cp_air*Ta0_P)*Ta0_P;
}
