/*
 * sph_lib.cc
 *
 *  Created on: Feb 18, 2015
 *      Author: zhixuanc
 */

//! The following two is for quick search function
#include <cstdio>
#include <cstdlib>

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <iostream> //for readFile
#include <fstream>  //for readFile
#include <string>   //for readFile
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <cmath>

//#include <petsc.h>
using namespace std;

#include <particler.h>
#include <sph_header.h>
#include <constant.h>
#  include <IndMap.h>
#include <parameters.h>

#include <properties.h> //for atmosphere interpolation
#include "meteo.h" //for atmosphere interpolation
#include "options.h"

MPI_Datatype MPI_TKEY, MPI_IndMap;

//! calculate dot-product
double
dot(double vec1[], double vec2[])
{
  double sum = 0.;
  for (int i = 0; i < DIMENSION; i++)
    sum += vec1[i] * vec2[i];
  return sum;
}


/*! Rotate vector of state vars {\it i.e.}
 * $\mathbf{u} = (\mathbf{u} \cdot \hat{n}) \hat{n}$
 */
void
rotate (double * u, double * cos)
{
  // u_r = dot (u, n_r) * n_r
  int i;
  double temp = dot (u, cos);

  for (i = 0; i < DIMENSION; i++)
    u[i] = temp * cos[i];
  return;
}

/*! reflect vector $\mathbf{u}_i$ in plane $\hat{n}$
 * $\mathbf{u}_r = \mathbf{u}^i - 2 (\mathbf{u}^i \cdot \hat{n})\hat{n}$
 */
void
reflect (double * ui, double * ur, double * n)
{
  int i;
  double vdotn = dot (ui, n);
  for (i = 0; i < DIMENSION; i++)
    ur[i] = ui[i] - 2 * vdotn * n[i];
  return;
}

/****************************************************
 * h is the smoothing length
 * the weight extends up to 3*h in all directions
 ***************************************************/
const double t1 = 0.564189583547756;    // 1/sqrt(pi)

// weight function
double
weight(double *s, double h)
{
  int i;
  double norm, expt, w;


  for (i = 0; i < DIMENSION; i++)
    if (abs(s[i]) > CUTOFF)
      return 0;

  norm = 1;
  expt = 0;
  for (i = 0; i < DIMENSION; i++)
  {
    norm *= t1 / h;
    expt += s[i]*s[i];
  }
  w = norm * exp(-expt);
  return w;
}

/****************************************************
 * d_weight is derivative of weight function in
 * any given direction
 * x-dir = 0
 * y-dir = 1
 * z-dir = 2
 ****************************************************/
// d(weight function)
double
d_weight (double *s, double h, int dir)
{
  int i;
  double dw, tmp;

  for (i = 0; i < DIMENSION; i++)
    if (abs (s[i]) > CUTOFF)
      return 0;

  tmp = -2 * s[dir] / h;
  dw = tmp * weight (s, h);
  return dw;
}

//Function to compute F form r_ab and W_ab
//The relationship is: F r_ab = dW_ab
double compute_F(double* dwdx, double* dx)
{
	int i;
	double r_sq=0.;
	double W_r=0.;
	double result;
	for (i=0; i<DIMENSION; i++)
	{
		r_sq += (*(dx+i)) * (*(dx+i));
		W_r += (*(dwdx+i)) * (*(dx+i));
	}

	result=W_r/r_sq;

#ifdef DEBUG
	bool print = true;
	if (print && result>=0)
        cout << "You got a non-negative Fab, something is wrong!" << endl;
#endif

	return result;
}

/*
 * z = p1 x + p2 y + p3 x y + p4
 */
void poly_surf (double * x, double * y, double * z, double * poly)
{

  // pt0 (x=0, y=0)
  poly[3] = z[0];

  // pt1 (x = dx, y = 0)
  poly[0] = (z[1] - poly[3]) / x[1];

  // pt3 (x = 0, y = dy)
  poly[1] = (z[3] - poly[3]) / y[1];

  // pt2 (x = dx, y = dx)
  poly[2] = (z[2] - (poly[0] * x[1] + poly[1] * y[1] + poly[3])) /
            (x[1] * y[1]);
  return;
}


/*
 * fit least square surface with 3 constants
 * z = p1 x + p2 y + p3
 */
void
lsq_surf3 (double * x, double * y, double * z, double * poly)
{
  double t1 = y[0] * y[0];
  double t2 = 2.;
  double t3 = y[1] * y[1];
  double t4 = x[1] * x[1];
  double t5 = x[2] * x[2];
  double t6 = x[0] * x[0];
  double t7 = y[2] * y[2];
  double t8 = x[0] * x[2];
  double t9 = x[1] * x[2];
  double t10 = x[0] * x[1];
  double t11 = (x[0] + x[1] - x[2]) * x[2];
  double t12 = (t3 + t1) * t5;
  double t13 = (t1 + t7) * t4 + t12 + t6 * (t7 + t3) +
               (-t9 * t1 + (y[0] * (-t10 + t11) - t8 * y[1]) * y[1] +
               ((-t8 + (x[2] + x[0] - x[1]) * x[1]) * y[0] +
                (-t9 + (x[2] + x[1] - x[0]) * x[0]) * y[1] -
                t10 * y[2]) * y[2]) * t2;
  double t14 = x[0] * z[0] + x[1] * z[1] + x[2] * z[2];
  double t15 = t2 * y[1];
  double t16 = (t15 - y[0] - y[2]) * x[1] +
               (t2 * y[2] - y[0] - y[1]) * x[2] +
               (t2 * y[0] - y[1] - y[2]) * x[0];
  double t17 = y[0] * z[0] + y[1] * z[1] + y[2] * z[2];
  double t18 = y[0] * y[1];
  double t19 = y[1] * y[2];
  double t20 = y[0] * y[2];
  t7 = (t18 + t19 - t1 - t7) * x[1] +
       (t20 + t19 - t1 - t3) * x[2] +
       (t18 + t20 - t3 - t7) * x[0];
  t18 = z[0] + z[1] + z[2];
  t13 = 1. / t13;
  t5 = (t4 + t5 - t10 - t8) * y[0] +
       (t6 + t4 - t8 - t9) * y[2] +
       (t6 + t5 - t10 - t9) * y[1];

  // solution
  poly[0] = (t2 * (t1 + (-y[0] + y[1]) * y[1] +
            (-y[0] - y[1] + y[2]) * y[2]) * t14 - t16 * t17 + t7 * t18) * t13;
  poly[1] = (-t16 * t14 + t2 * (t6 + (-x[0] + x[1]) * x[1] - t11) *
             t17 - t5 * t18) * t13;
  poly[2] = (t7 * t14 - t5 * t17 +
             (t4 * t1 - t15 * t10 * y[0] + t12 +
             t6 * t3 + (-t2 * x[2] * (x[0] * y[0] + x[1] * y[1])
                        + (t6 + t4) * y[2]) * y[2]) * t18) * t13;
  return;
}


//general function for artificial viscosity
double art_vis ( Particle * pi, Particle * pj)
{
	int    k;
	double rhoab = 0.5 * ((pi->get_density () ) + (pj->get_density () ));
	double sndspdab = 0.5 * ((pi->get_sound_speed () ) + (pj->get_sound_speed () ));
	double rab [DIMENSION];
	double vab [DIMENSION];
	double miuab;
	double vrab = 0.;
	double rsqab = 0.;
	double vis = 0.;
	double h = (pi->get_smlen());

	for (k = 0; k < DIMENSION; k++)
		rab [k]= *(pi->get_coords() + k) - *(pj->get_coords() + k);

	for (k = 0; k < DIMENSION; k++)
		vab [k]= *(pi->get_vel() + k) - *(pj->get_vel() + k);

	for (k = 0; k < DIMENSION; k++)
		vrab += rab[k] * vab[k];

	for (k = 0; k < DIMENSION; k++)
		rsqab += (rab[k] * rab[k]);

	miuab = h * vrab / (rsqab + 0.01 * h * h);

	switch ( vrab >= 0 )
	{
	case 1:
		vis = 0.;
	    break;
	case 0:
		vis = (2. * miuab - sndspdab) * miuab / rhoab;
	    break;
	}

	return vis;
}

//overloading function for artificial viscosity
double art_vis ( double rhoab, double sndspdab, double rab[DIMENSION], double vab[DIMENSION], double rsqab, double h)
{
	int    k;
	double miuab;
	double vrab = 0.;
	double vis = 0.;

	for (k = 0; k < DIMENSION; k++)
	     vrab += (rab[k] * vab[k]);

	miuab = h * vrab / (rsqab + ata_P * h * h);

#ifdef 	USE_PHYSICS_VIS
	vis =  (- alf_P*sndspdab) * miuab / rhoab; //beta will become zero, and it is not necessary to turn off viscosity for receding
#else
    if (vrab > 0)
		vis = 0.;
    else if (vrab < 0)
		vis = (beta_P * miuab - alf_P*sndspdab) * miuab / rhoab;
#endif

	return vis;
}

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) and internal energy
//Based on realistic atmosphere model
void air_prop_realistic (double *coord, double * energy, double *pressure, double * density)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop! \n" <<endl;

	double T;
	T = Ta0_P *(h < 0) + (Ta0_P-miu1_P*h)*((h>=0)&&(h<H1_P))+(Ta0_P-miu1_P*H1_P)*((h>=H1_P)&&(h<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h-H2_P))*((h>=H2_P)&&(h<H3_P));
	double C0 = -0.034193145144839; //coefficient in expression of pressure: C0=-28.97*g/(6.02*1000*1.3806448)
	*pressure = pa0_P*exp(C0*h/T)*(h>=0) +  pa0_P*(h<0);
	*density = (*pressure) /(Ra_P*T) ;
	*energy = (*pressure) /(*density * (gamma_P-1));

#ifdef DEBUG
	bool print = false;
	if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
}

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Based on realistic atmosphere model
void air_prop_realistic (double *coord, double *range, double * energy, double *pressure, double * density, double * mass)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop! \n" <<endl;

	double T;
	T = Ta0_P *(h < 0) + (Ta0_P-miu1_P*h)*((h>=0)&&(h<H1_P))+(Ta0_P-miu1_P*H1_P)*((h>=H1_P)&&(h<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h-H2_P))*((h>=H2_P)&&(h<H3_P));
	double C0 = -0.034193145144839; //coefficient in expression of pressure: C0=-28.97*g/(6.02*1000*1.3806448)
	*pressure = pa0_P*exp(C0*h/T)*(h>=0)+pa0_P*(h<0);
	*density = (*pressure) /(Ra_P*T) ;
	*energy = (*pressure) /(*density * (gamma_P-1));

	//integration is based on interpolation of order 2 (three points)
	double h1=range[4];
	double h2=range[5];
	double T1, T2;
	T1 = Ta0_P *(h1 < 0) + (Ta0_P-miu1_P*h1)*((h1>=0)&&(h1<H1_P))+(Ta0_P-miu1_P*H1_P)*((h1>=H1_P)&&(h1<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h1-H2_P))*((h1>=H2_P)&&(h1<H3_P));
	T2 = Ta0_P *(h2 < 0) + (Ta0_P-miu1_P*h2)*((h2>=0)&&(h2<H1_P))+(Ta0_P-miu1_P*H1_P)*((h2>=H1_P)&&(h2<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h2-H2_P))*((h2>=H2_P)&&(h2<H3_P));
	double p1 = pa0_P*exp(C0*h1/T1)*(h1>=0)+pa0_P*(h1<0);
	double d1 = p1/(Ra_P*T1);
	double p2 = pa0_P*exp(C0*h2/T2)*(h2>=0)+pa0_P*(h2<0);
	double d2 = p2/(Ra_P*T2);
	double d = *density;

	//The following code is based on numerical integration, coefficient are 1/6, 4/6 1/6
	*mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])*(0.1666667*d1+0.6666666*d+0.1666667*d2);

#ifdef DEBUG
	bool print = false;
	if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
}

//function that will determine pressure based hydro-static equation
double determine_pressure_hydro (double h)
{
	if (h < 0)
		return pa0_P;
	else if ((h>=0)&&(h<H1_P))
		return (Ata1_p*pow(Ta0_P-miu1_P*h, Ate1_p));
	else if ((h>=H1_P)&&(h<H2_P))
		return (Ata2_p*exp(Atb2_p*h));
	else if ((h>=H2_P)&&(h<H3_P))
		return (AtC3_p*pow(Atb3_p*h+Ata3_p,Ate3_p));
	else
		cout << "input height makes nonsense!" <<endl;
}

//function that will determine temperature based hydro-static equation
double determine_temperature_hydro (double h)
{
	if (h < 0)
		return Ta0_P;
	else if ((h>=0)&&(h<H1_P))
		return (Ta0_P-miu1_P*h);
	else if ((h>=H1_P)&&(h<H2_P))
		return (Ta0_P-miu1_P*H1_P);
	else if ((h>=H2_P)&&(h<H3_P))
		return (Ta0_P-miu1_P*H1_P+miu2_P*(h-H2_P));
	else
		cout << "input height makes nonsense!" <<endl;
}
//function that used to determine the property of air: density, pressure, (temperature not explicitly output) and internal energy
//Based on hydro-static equation dp/dz=-rho*g ---> This will give a less realistic initial atmosphere, but more consistent with current model.
void air_prop_hydro (double *coord, double * energy, double *pressure, double * density)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop! \n" <<endl;

	double T;
	//T = Ta0_P *(h < 0) + (Ta0_P-miu1_P*h)*((h>=0)&&(h<H1_P))+(Ta0_P-miu1_P*H1_P)*((h>=H1_P)&&(h<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h-H2_P))*((h>=H2_P)&&(h<H3_P));
	T= determine_temperature_hydro (h);
	//*pressure = pa0_P *(h < 0) + (Ata1_p*pow(Ta0_P-miu1_P*h, Ate1_p))*((h>=0)&&(h<H1_P))+(Ata2_p*exp(Atb2_p*h))*((h>=H1_P)&&(h<H2_P))+(AtC3_p*pow(Atb3_p*h+Ata3_p,Ate3_p))*((h>=H2_P)&&(h<H3_P));
	*pressure = determine_pressure_hydro(h);

	*density = (*pressure) /(Ra_P*T) ;
	*energy = (*pressure) /(*density * (gamma_P-1));

#ifdef DEBUG
	bool print = false;
	if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
}

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Based on hydro-static equation dp/dz=-rho*g ---> This will give a less realistic initial atmosphere, but more consistent with current model.
void air_prop_hydro (double *coord, double *range, double * energy, double *pressure, double * density, double * mass)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop! \n" <<endl;

	double T;
	//T = Ta0_P *(h < 0) + (Ta0_P-miu1_P*h)*((h>=0)&&(h<H1_P))+(Ta0_P-miu1_P*H1_P)*((h>=H1_P)&&(h<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h-H2_P))*((h>=H2_P)&&(h<H3_P));
	T= determine_temperature_hydro (h);

	//*pressure = pa0_P *(h < 0) + (Ata1_p*pow(Ta0_P-miu1_P*h, Ate1_p))*((h>=0)&&(h<H1_P))+(Ata2_p*exp(Atb2_p*h))*((h>=H1_P)&&(h<H2_P))+(AtC3_p*pow(Atb3_p*h+Ata3_p,Ate3_p))*((h>=H2_P)&&(h<H3_P));
	*pressure = determine_pressure_hydro(h);

	*density = (*pressure) /(Ra_P*T) ;
	*energy = (*pressure) /(*density * (gamma_P-1));

	//integration is based on interpolation of order 2 (three points)
	double h1=range[4];
	double h2=range[5];
	double T1, T2;
//	T1 = Ta0_P *(h1 < 0) + (Ta0_P-miu1_P*h1)*((h1>=0)&&(h1<H1_P))+(Ta0_P-miu1_P*H1_P)*((h1>=H1_P)&&(h1<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h1-H2_P))*((h1>=H2_P)&&(h1<H3_P));
//	T2 = Ta0_P *(h2 < 0) + (Ta0_P-miu1_P*h2)*((h2>=0)&&(h2<H1_P))+(Ta0_P-miu1_P*H1_P)*((h2>=H1_P)&&(h2<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h2-H2_P))*((h2>=H2_P)&&(h2<H3_P));
	T1= determine_temperature_hydro (h1);
	T2= determine_temperature_hydro (h2);
	//double p1 = pa0_P *(h1 < 0) + (Ata1_p*pow(Ta0_P-miu1_P*h1, Ate1_p))*((h1>=0)&&(h1<H1_P))+(Ata2_p*exp(Atb2_p*h1))*((h1>=H1_P)&&(h1<H2_P))+(AtC3_p*pow(Atb3_p*h1+Ata3_p,Ate3_p))*((h1>=H2_P)&&(h1<H3_P));
	double p1 = determine_pressure_hydro(h1);
	double d1 = p1/(Ra_P*T1);
	//There was a bug here, h1 was used for computing of p2----> a stupid mistake
	//double p2 = pa0_P *(h2 < 0) + (Ata1_p*pow(Ta0_P-miu1_P*h2, Ate1_p))*((h2>=0)&&(h2<H1_P))+(Ata2_p*exp(Atb2_p*h2))*((h2>=H1_P)&&(h2<H2_P))+(AtC3_p*pow(Atb3_p*h2+Ata3_p,Ate3_p))*((h2>=H2_P)&&(h2<H3_P));
	double p2 = determine_pressure_hydro(h2);
	double d2 = p2/(Ra_P*T2);
	double d = *density;

	//The following code is based on numerical integration, coefficient are 1/6, 4/6 1/6
	*mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])*(0.1666667*d1+0.6666666*d+0.1666667*d2);

#ifdef DEBUG
	bool print = false;
	if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
}

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy
//Function that will give a "uniform temperature" atmosphere ---> which is not realistic for sure
//Used for code testing
//Gravity coefficient will be set to 9.81 -->This is different from uniform environment.
void air_prop_uniformT (double *coord, double * energy, double *pressure, double * density)
{
        double h=coord[2];

        if (h>H3_P)
                cout << "height of domain exceeds the maximum height of atmosphere, in air_prop! \n" <<endl;

        double T;
        T = Ta0_P;
        *pressure = pa0_P*exp(Atf_P * h);
        *density = (*pressure) /(Ra_P*T) ;
        *energy = (*pressure) /(*density * (gamma_P-1));

#ifdef DEBUG
        bool print = false;
        if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
}
//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Function that will give a "uniform temperature" atmosphere ---> which is not realistic for sure
//Used for code testing
//Gravity coefficient will be set to 9.81 -->This is different from uniform environment.
void air_prop_uniformT (double *coord, double *range, double * energy, double *pressure, double * density, double * mass)
{
        double h=coord[2];

        if (h>H3_P)
                cout << "height of domain exceeds the maximum height of atmosphere, in air_prop! \n" <<endl;

        double T;
        T = Ta0_P;
        *pressure = pa0_P*exp(Atf_P * h);
        *density = (*pressure) /(Ra_P*T) ;
        *energy = (*pressure) /(*density * (gamma_P-1));

#ifdef DEBUG
        bool print = false;
        if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
        //integration is based on interpolation of order 2 (three points)
        double h1=range[4];
        double h2=range[5];
        double p1 = pa0_P * exp(Atf_P * h1);
        double d1 = p1/(Ra_P*Ta0_P);
        double p2 = pa0_P * exp(Atf_P * h2);
        double d2 = p2/(Ra_P*Ta0_P);
        double d = *density;
        //The following code is based on numerical integration, coefficient are 1/6, 4/6 1/6
        *mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])*(0.1666667*d1+0.6666666*d+0.1666667*d2);        
}

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy
//Function that will give a "uniform" atmosphere ---> which is not realistic for sure
//Used for code testing
//Gravity coefficient will be set to zero -->This is real uniform environment.
void air_prop_uniform (double *coord, double * energy, double *pressure, double * density)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop_uniform ! \n" <<endl;

	double T;
	T = Ta0_P;
	*pressure = pa0_P;

#if FLUID_COMPRESSIBILITY==0
	*density = (*pressure) /(Ra_P*T) ;
	*energy = (*pressure) /(*density * (gamma_P-1));
#elif FLUID_COMPRESSIBILITY==1
	//Here I am using a constant for water sound speed, need to figure out a more genereal way
	double gama=7.0;
	double water_sndspd = 1482;
	double B=rhoa0_P*water_sndspd*water_sndspd/gama;
	*density = rhoa0_P*pow((B/(*pressure)+1),1/gama);
	*energy = T*Cva_P;
#endif

#ifdef DEBUG
	bool print = false;
	if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
}

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Function that will give a "uniform" atmosphere ---> which is not realistic for sure
//Used for code testing
////Gravity coefficient will be set to zero -->This is real uniform environment.
void air_prop_uniform (double *coord, double *range, double * energy, double *pressure, double * density, double * mass)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop_uniform ! \n" <<endl;

	double T;
	T = Ta0_P;
	*pressure = pa0_P;

#if FLUID_COMPRESSIBILITY==0
	*density = (*pressure) /(Ra_P*T) ;
	*energy = (*pressure) /(*density * (gamma_P-1));

	double d = *density;

	//The following code is based on numerical integration, coefficient are 1/6, 4/6 1/6
    //	*mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])*(0.1666667*d1+0.6666666*d+0.1666667*d2);

	*mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])*d;
#elif FLUID_COMPRESSIBILITY==1
	//Here I am using a constant for water sound speed, need to figure out a more genereal way
	double gama=7.0;
	double water_sndspd = 1482;
	double B=rhoa0_P*water_sndspd*water_sndspd/gama;
	*density = rhoa0_P*pow(((*pressure - pa0_P)/B+1),1/gama);
	*energy = T*Cva_P;

	*mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])* (*density);
#endif

#ifdef DEBUG
	bool print = false;
	if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
}

#if ATMOSPHERE_TYPE==4
//function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy
//Based on input meteo data, using interpolation
void air_prop_meteo_based (SimProps * simprops, double *coord, double * energy, double *pressure, double * density)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop_uniform ! \n" <<endl;

	int np = (simprops->meteo_data)->get_number_of_props();
	double prop[np];

	(simprops->meteo_data)->interpolate (h*0.001, prop); //need to convert meters to kilometers
	*pressure = prop[1]*100; //because the metric is hPa
	*density = prop[0] ;
	*energy = (*pressure) /(*density * (gamma_P-1));
}

//Overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy
//Based on input meteo data, using interpolation
void air_prop_meteo_based (SimProps * simprops, double *coord, double *range, double * energy, double *pressure, double * density, double * mass)
{
	double h=coord[2];

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop_uniform ! \n" <<endl;

	int np = (simprops->meteo_data)->get_number_of_props();
	double prop[np];

	(simprops->meteo_data)->interpolate (h*0.001, prop);
	*pressure = prop[1]*100; //because the metric is hPa
	*density = prop[0] ;
	*energy = (*pressure) /(*density * (gamma_P-1));

	//integration is based on interpolation of order 2 (three points)
	double h1=range[4];
	double h2=range[5];
	double d1 = (simprops->meteo_data)->interpolate(h1*0.001, 0);
	double d2 = (simprops->meteo_data)->interpolate(h2*0.001, 0);
	double d = *density;

	//The following code is based on numerical integration, coefficient are 1/6, 4/6 1/6
	*mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])*(0.1666667*d1+0.6666666*d+0.1666667*d2);
}
#endif

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) and internal energy
void air_prop(SimProps * simprops, double *coord, double * energy, double *pressure, double * density)
{
#if ATMOSPHERE_TYPE==0
	air_prop_realistic(coord, energy, pressure, density);
#elif ATMOSPHERE_TYPE==1
	air_prop_hydro(coord, energy, pressure, density);
#elif ATMOSPHERE_TYPE==2
	air_prop_uniform(coord, energy, pressure, density);
#elif ATMOSPHERE_TYPE==3
    air_prop_uniformT(coord, energy, pressure, density);
#elif ATMOSPHERE_TYPE==4
    air_prop_meteo_based(simprops, coord, energy, pressure, density);
#endif
}

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
void air_prop(SimProps * simprops, double *coord, double *range, double * energy, double *pressure, double * density, double * mass)
{
#if ATMOSPHERE_TYPE==0
	air_prop_realistic(coord, range, energy, pressure, density, mass);
#elif ATMOSPHERE_TYPE==1
	air_prop_hydro(coord, range, energy, pressure, density, mass);
#elif ATMOSPHERE_TYPE==2
	air_prop_uniform(coord, range, energy, pressure, density, mass);
#elif ATMOSPHERE_TYPE==3
    air_prop_uniformT(coord, range, energy, pressure, density, mass);
#elif ATMOSPHERE_TYPE==4
    air_prop_meteo_based(simprops, coord, range, energy, pressure, density, mass);
#endif
}

//function that used to determine the pressure of atmosphere
double determine_pressure(SimProps * simprops, double h)
{
    if (h>H3_P)
    	cout << "height of domain exceeds the maximum height of atmosphere, in air_prop_uniform ! \n" <<endl;

	double pressure;

#if ATMOSPHERE_TYPE==0
	double C0 = -0.034193145144839; //coefficient in expression of pressure: C0=-28.97*g/(6.02*1000*1.3806448)
	double T = Ta0_P *(h < 0) + (Ta0_P-miu1_P*h)*((h>=0)&&(h<H1_P))+(Ta0_P-miu1_P*H1_P)*((h>=H1_P)&&(h<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h-H2_P))*((h>=H2_P)&&(h<H3_P));
	pressure = pa0_P*exp(C0*h/T)*(h>=0)+pa0_P*(h<0);
#elif ATMOSPHERE_TYPE==1
	pressure = determine_pressure_hydro(h);
#elif ATMOSPHERE_TYPE==2
	pressure = pa0_P;
#elif ATMOSPHERE_TYPE==3
	pressure = pa0_P*exp(Atf_P * h);
#elif ATMOSPHERE_TYPE==4

    int np = (simprops->meteo_data)->get_number_of_props();
    double prop[np];

    (simprops->meteo_data)->interpolate (h*0.001, prop);
    pressure = prop[1]*100; //because the metric is hPa
#endif
    return pressure;
}

//function that determines parameters of certain particle
void initial_air (Particle * pi, SimProps * simprops)
{
	  int i;

	  double xi[DIMENSION];
	  double vel[DIMENSION] = {0., 0., 0.};  //initial velocity need to be set to zero
	  double prss, erg, dens, mss;
	  double range[6];
	  double sml2;

	  	  for (i = 0; i < DIMENSION; i++)
	  		  xi[i] = *(pi->get_coords() + i);

	  	  sml2 = 0.5*(pi->get_smlen ());
	  	  range[0]=xi[0]-sml2;
	  	  range[1]=xi[0]+sml2;
	  	  range[2]=xi[1]-sml2;
	  	  range[3]=xi[1]+sml2;
	  	  range[4]=xi[2]-sml2;
	  	  range[5]=xi[2]+sml2;

	  	  air_prop(simprops, xi, range, &erg, &prss, &dens, &mss);

	      //put data back into particle:
	      pi->put_density(dens);
	      pi->put_energy(erg);
	      pi->put_pressure(prss);
	      pi->put_velocity(vel);
	      pi->put_mass(mss);

	      //the second variable need to be updated.
	      pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P);

	return;
}

//function that used to determine only internal energy based on altitude.
//This function is based on a less realistic model: hydrostatic model
//This function will be used while imposing wall boundary condition
//actually, the way to imposing essential boundary condition for internal energy is not a proper way in bcond,cc
//But temporarily, I just use this not proper way to impose boundary condition.---> need to read more papers on how to imposing essential boundary condition in SPH method.
double air_engr_hydro (double *coord)
{
	double h=coord[2];
	double density, pressure, energy;
	double T;

	if (h>H3_P)
		cout << "height of domain exceeds the maximum height of atmosphere, in air_prop! \n" <<endl;

	T = Ta0_P *(h < 0) + (Ta0_P-miu1_P*h)*((h>=0)&&(h<H1_P))+(Ta0_P-miu1_P*H1_P)*((h>=H1_P)&&(h<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h-H2_P))*((h>=H2_P)&&(h<H3_P));
	pressure = pa0_P *(h < 0) + (Ata1_p*pow(Ta0_P-miu1_P*h, Ate1_p))*((h>=0)&&(h<H1_P))+(Ata2_p*exp(Atb2_p*h))*((h>=H1_P)&&(h<H2_P))+(AtC3_p*pow(Atb3_p*h+Ata3_p,Ate3_p))*((h>=H2_P)&&(h<H3_P));
	density = (pressure) /(Ra_P*T) ;
	energy = (pressure) /(density * (gamma_P-1));

	return energy;
}

////Function that used Petsc GESRM to solve system of equations to return mass of each particles
//PETSC_EXTERN PetscErrorCode PCCreate_Jacobi(PC);
//void solve_mass(double **A, double *b, double * x, int N, int NP, int id)
///*note:
// *  *     dimension of A is nx * ny;
// *   *     dimension of b is ny * 1
// *    *     */
////The solver has already been tested. it works very well! I can directly use it in my code
//{
//  int argc=0;
//  char **argv=NULL;
//  int i, j;
//  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
//
//  //PetscErrorCode ierr;
//  Mat A_petsc;
//  Vec b_petsc, x_petsc;//, u;
//
//  VecCreate(PETSC_COMM_WORLD,&b_petsc);
//  VecSetSizes(b_petsc,N,PETSC_DECIDE);
//  VecSetFromOptions(b_petsc);
//
//  MatCreate(PETSC_COMM_WORLD,&A_petsc);
//  MatSetSizes(A_petsc,N,N,PETSC_DECIDE,PETSC_DECIDE);
//  MatSetFromOptions(A_petsc);
//  MatSetUp(A_petsc);
//
//  PetscInt Low, Up;
//  VecGetOwnershipRange(b_petsc, &Low, &Up);
//  for (i=Low; i<Up; i++)
//     VecSetValues(b_petsc,1,&i, &(b[i-Low]),INSERT_VALUES);
//
//  VecAssemblyBegin(b_petsc);
//  VecAssemblyEnd (b_petsc);
//
//  for (i=Low; i<Up; i++)
//     for (j=0; j<NP; j++)
//         MatSetValues(A_petsc, 1, &i, 1, &j, &(A[i-Low][j]), INSERT_VALUES);
//
//  MatAssemblyBegin(A_petsc, MAT_FINAL_ASSEMBLY);
//  MatAssemblyEnd (A_petsc, MAT_FINAL_ASSEMBLY);
///*
//   PetscViewer viewer;
//   VecView(b_petsc, PETSC_VIEWER_STDOUT_WORLD);
//   MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD);
//*/
// //solve_try( argc, argv, A_petsc, x_petsc, u, N, N, id, &its);
// // /*
// //  *  *      Create linear solver context
// //   *   *        */
//
//  KSP ksp;
//  PC  pc;
//  //PetscScalar one = 1.0, neg_one = -1.0;
//  //PetscReal      norm;
//
//  //MatMult(A_petsc,b_petsc,x_petsc);
//
//  KSPCreate(PETSC_COMM_WORLD,&ksp);
//  //KSPSetOperators(ksp, A_petsc, A_petsc, DIFFERENT_NONZERO_PATTERN);
//  KSPSetOperators(ksp, A_petsc, A_petsc);
//
//  PCRegister("ourjacobi",PCCreate_Jacobi);
//  KSPGetPC(ksp,&pc);
//  PCSetType(pc,"ourjacobi");
//  KSPSetFromOptions(ksp);
//
//  VecDuplicate(b_petsc,&x_petsc);
//  KSPSolve(ksp,b_petsc,x_petsc);
//
//  //VecAXPY(b_petsc,neg_one,u);
//  //VecNorm(b_petsc,NORM_2,&norm);
//  PetscInt its;
//  KSPGetIterationNumber(ksp,&its);
//
////  PetscPrintf(PETSC_COMM_WORLD,"number of iterations %D\n",*its);
//  //VecView(x_petsc, PETSC_VIEWER_STDOUT_WORLD);
//  //VecView(b_petsc, PETSC_VIEWER_STDOUT_WORLD);
//  //VecView(u, PETSC_VIEWER_STDOUT_WORLD);
//
//  PetscScalar hh;
//  VecGetOwnershipRange(b_petsc, &Low, &Up);
//  for (i=Low; i<Up; i++)
//  {
//      VecGetValues(x_petsc, 1, &i, &hh);
//     x[i-Low]=hh;
//  }
//  KSPDestroy(&ksp);
//
//  MatDestroy (&A_petsc);
//  VecDestroy (&b_petsc);
//  VecDestroy (&x_petsc);
//  //VecDestroy (&u);
//
//  PetscFinalize();
//}

/*function that get mass of and smooth length
 * it is suitable for space where density is uniform
 * Domain is 2.5D, dimension of the third direction is 1
 *
 * input: np number of particles per layer
 * dom, length of the domain within which np particles will be added, the third dimension of the domain should be 1
 * density of material within the domain: assuming uniform density here
 * output:mss and sml
*/
void Compute_mass (int np, double * dom, double des, double *mss, double *sml)
{
	double volume;
	int i;
	volume = 1;
	for (i=0; i<2; i++) //assuming unit thickness in z direction
		volume *= dom[i];

	*sml = sqrt(volume/np);
	*mss = *sml * des*volume/np;
}

/*OVERLOAD FUNCTION
 * function that get mass of and smooth length
 * it is suitable for space where density is uniform
 * Domain is 3D
 *
 * input: np: total number of particles
 * dom, length of the domain within which np particles will be added
 * density of material within the domain: assuming uniform density here
 * output:mss and sml
*/
void Compute_mass (int np, double *range_x, double *range_y, double *range_z, double des, double *mss, double *sml)
{
	int i;
	double dom[DIMENSION];
	    dom[0]=range_x[1]-range_x[0];
	    dom[1]=range_y[1]-range_y[0];
	    dom[2]=range_z[1]-range_z[0];

	double volume;

	volume = 1.;
	for (i=0; i<DIMENSION; i++)
		volume *= dom[i];

	*sml = cbrt(volume/np);
	*mss = des*volume/np;
}

//function that determine whether particle is within bucket or nor
bool in_bucket(double *bkmax, double *bkmin, double *pcrd)
{
	bool temp = true;
	int i;
	for (i=0; i< DIMENSION; i++)
		temp &=((pcrd[i]<=bkmax[i])&&(pcrd[i]>=bkmin[i]));

	return temp;
}

//function that used to sort process id to form comm
void id_sort(int arr[],int n)
{
   int temp;
   int i, j;
   for (i=0; i< n; i++)
   {
       for (j=0; j< n; j++)
           {
              if (arr[j]<arr[i])
                {
                   temp = arr[i];
                   arr[i] = arr[j];
                   arr[j] = temp;
                }
           }
   }
}

//Define a new MPI data type corresponding to IndMap
void create_indmpi_struct()
{
    int one=1;
    MPI_Aint zero=0;
    int tkeylength = TKEYLENGTH;
    MPI_Datatype unsignd = MPI_UNSIGNED;
    MPI_Type_struct( one, &tkeylength, &zero, &unsignd, &MPI_TKEY );
    MPI_Datatype mpitype[2] = {MPI_TKEY, MPI_INT};
    int block_sz[2]={1,3};
    MPI_Aint initial_bt[2]={0,12};//initial byte displacement
    MPI_Type_struct( one, block_sz, initial_bt, mpitype, &MPI_IndMap);
    MPI_Type_commit(&MPI_IndMap);
}

//Send and receive IndexMap among neighbors
void
exchange_indmap (int nump, int myid, int * my_comm,
		IndMap *array_local, IndMap **array, int nrp, int * nrp_all)
{
    // if number of procs == 1 , don't waste time here
    if (nump < 2)
        return;

    int i, j, ierr;
    MPI_Status status1;
    MPI_Status status;

    // send_info array
    int *send_info = new int[nump];
    int *recv_info = new int[nump];

    for (i = 0; i < nump; i++)
    {
        send_info[i] = 0;
        recv_info[i] = 0;
    }

    MPI_Request * reqinfo = new MPI_Request [2*nump];
    int tag1 = 454004;

    // post recveives for size info
    for (i = 0; i < nump; i++)
           if ( my_comm[i] )
                ierr = MPI_Irecv ((recv_info + i), 1, MPI_INT, i, tag1,
                                  MPI_COMM_WORLD, (reqinfo + i));

    //determine how many IndMaps should I send to each process
    for (i=0; i<nump; i++)
    {
    	if (my_comm[i])
    		send_info[i]=nrp;
    }

    /* send out size information */
    for (i = 0; i < nump; i++)
        if ( my_comm[i] )
            ierr = MPI_Isend ((send_info + i), 1, MPI_INT, i, tag1,
                              MPI_COMM_WORLD, (reqinfo + i + nump));

    // Wait to receive size info from others
    for (i = 0; i < nump; i++)
        if ( my_comm[i] )
            MPI_Wait ((reqinfo + i), & status1);

    /* post receives */
    // size of data to be received
    int recv_count= 0;
    for (i = 0; i < nump; i++)
        recv_count += recv_info[i];    //number of Indmap

    *nrp_all = recv_count + nrp;      //total number of Indmap should also include local

    // allocate space for data to be received
    MPI_Request *recv_req = new MPI_Request[nump];
    *array = new IndMap[*nrp_all];
    IndMap *rec_buf = *array;
    int counter_recv=0;
    int tag2 = 487358;

    for (i = 0; i < nump; i++)
    {
        if (recv_info[i] != 0)
        {
            j = MPI_Irecv ((rec_buf + counter_recv), recv_info[i],
            		MPI_IndMap, i, tag2, MPI_COMM_WORLD,
                   (recv_req +i));
            counter_recv += recv_info[i];
        }
    }    /* done with receives */

    //send Indmap
    MPI_Request *send_req = new MPI_Request[nump];
    int counter = 0;
    for (i = 0; i < nump; i++)
    {
        if (send_info[i] != 0)
        {
            j = MPI_Isend ((array_local + counter), send_info[i],
            		MPI_IndMap, i, tag2, MPI_COMM_WORLD,
                           (send_req + i));
            counter += send_info[i];
        }
    }

    //Make sure that receive is completed!
    for (i = 0; i < nump; i++)
           if (recv_info[i] != 0)
        	    MPI_Wait ((recv_req + i), &status);

    // make sure all the sends are completed
    for (i = 0; i < nump; i++)
        if ( my_comm[i] )
            ierr = MPI_Wait ((reqinfo + nump + i), & status1);

    // now check data sends
    for (i = 0; i < nump; i++)
        if (send_info[i] != 0)
            ierr = MPI_Wait ((send_req + i), &status);

    //Finally, I need add all of my own particles in the global array.
    for (i=0; i<nrp; i++)
    	rec_buf[recv_count+i]= array_local [i];

    // clean up
    delete [] reqinfo;
    delete [] recv_info;
    delete [] recv_req;

    delete [] send_req;
    delete [] send_info;

    return;
}

// function that used to determine the velocity based on a parabolic velocity profile

double parabolic_vel(double R, double rsq, double u_avg)
{
	double umax=u_avg*2;
	return umax*(1-(rsq/(R*R)));
}

// function that used to determine the velocity based on a uniform velocity profile
double uniform_vel(double u_avg)
{
	return u_avg;
}

//function that determine the velocity profile
double vel_prof(double R, double rsq, double u_avg)
{
#if ERUPT_VELOCITY_PROF==0
	return uniform_vel(u_avg);
#elif ERUPT_VELOCITY_PROF==1
	return parabolic_vel(R, rsq, u_avg);
#endif
}

//function that used to determine the value of face by the face's index;
//--->This is actually exactly the same as the face type determine function from preprocess
//--->Remember to update this section when any modification is made in preprocess
int determine_face_type (double crd, double max, double min)
{
	int flag;

	if (crd<min)
		flag = -1;
	else if (crd > max)
		flag =1;
	else
		flag =0;

	return flag;
}

//function that used to determine the type of bucket
//--->This is actually exactly the same as the face type determine function from preprocess
//--->Remember to update this section when any modification is made in preprocess
void determine_bucket_type (double *mincrd, double *maxcrd, double *xcrd, double *ycrd, double *zcrd, int* type, int* bt)
{
	int flag[DIMENSION];
	int sum = 0;
	int k;

    bt[0]=determine_face_type(xcrd[0],maxcrd[0],mincrd[0]);
    bt[1]=determine_face_type(xcrd[1],maxcrd[0],mincrd[0]);
    bt[2]=determine_face_type(ycrd[0],maxcrd[1],mincrd[1]);
    bt[3]=determine_face_type(ycrd[1],maxcrd[1],mincrd[1]);
    bt[4]=determine_face_type(zcrd[0],maxcrd[2],mincrd[2]);
    bt[5]=determine_face_type(zcrd[1],maxcrd[2],mincrd[2]);

    for (k=0; k<DIMENSION; k++)
    	flag[k]=abs(bt[2*k]+bt[2*k+1]);

    for (k=0; k<DIMENSION; k++)
    	sum += flag[k];

    if ( (flag[0]==2) || (flag[1]==2) || (flag[2]==2))
    {
    	if ((bt[4] == -1) && (bt[5] == -1)) //only bt[5] == -1 is enough!
    		*type = 1; //underground
    	else
    		*type = 4; //pressure bc
    }

    else if ( (flag[0]==1) || (flag[1]==1) || (flag[2]==1) )
    	*type = 2;     //Mixed
    else if(sum == 0)
    	*type = 3;     //overground
    else
    {
    	*type = 0;     //invalid
    	cout << "Invalid input bucket_index!!!" << endl;
    }

    /*Some buckets which are "underground" but also has ground boundary information are shift to MIXED

        |
        |                                        ORIGINAL DOMAIN
        |
        |                 Left boundary
        |______________|_______!_______
        |Orig:UNDG     |       !       |
        |._._._._._._._|._._._.!_._._._|._._._.  Bottom Boundary
        |Changed to:   |       !       |
        |Mixed___ _____|_______!_______|_______
        |              |       !       |
        |              |       !       |
        |              |       !       |
        |______________|_______!_______|________________________

         The bucket at the left
         corner should still
         be UNDERGROUND!
    */

    if (bt[4]==-1 && bt[5]==0)
    	*type = 2;

    return ;
}


//function that used to compute the additional term in momentum equation if SPH_epsilon turbulence model is adopted
double SPH_epsilon_mom(double* vab, double V_b) //V_b is the specific volume
{
   double dotv=0.;
   for (int i=0; i<DIMENSION; i++)
	   dotv += (*(vab+i)) * (*(vab+i)); //v dot v

   return (EPSILON_HALF*dotv*V_b);
}

//function that used to compute turbulent heat conductivity in energy equation if SPH_epsilon turbulence model is adopted
double SPH_epsilon_heat_conductivity(double Cp_ab, double * ds, double *vab)
{
	int i;
	double dotvv=0., dotvr=0., dotrr=0.;
	double kab;

	/*
	 * when ds=[-577.7861  577.7746  -74.8549], vab=[-5.8261    5.9188   90.6049]  ds*vab is very small, this is not reasonable
	 */

	for (i=0; i<DIMENSION; i++)
		dotvr += abs((*(vab+i)) * (*(ds+i))); //add abs here to avoid superficial large heat conductivity
                                              //Actually, I believe this kind of treatment is reasonable, at least will not be worse than original equation.
                                              //The derivation of original expression is based on some viscous equation which comes from one paper in 1980s and no one explain how the equation obtained.
                                              //Any way, theories on SPH is far away from complete and there are so many abitary treatment.....

	if(dotvr==0.) //To be consistent with artificial viscosity --> But in my opinion, it is not really reasonable!
		return 0.;
	else
	{

	  for (i=0; i<DIMENSION; i++)
		dotrr += (*(ds+i)) * (*(ds+i));

	  for (i=0; i<DIMENSION; i++)
	  	dotvv += (*(vab+i)) * (*(vab+i));

      if (DIMENSION==3)
	    kab=6.0*EPSILON_HALF*Cp_ab*dotrr*dotvv/(PRANDTL_NUM*dotvr); //For 3D shear viscosity is 1/6 of h*alf*soundspeed
      else if (DIMENSION==2)
        kab=4.8*EPSILON_HALF*Cp_ab*dotrr*dotvv/(PRANDTL_NUM*dotvr); //For 2D shear viscosity is 5/24 of h*alf*soundspeed
      else
    	  cout<< "Dimension is neither 2 or 3, this program is not supposed to handle it!" <<endl;

#ifdef DEBUG
	bool print = true;
	if (print && kab<0)
        cout << "You got a negative heat conductivity, something is wrong!" << endl;
#endif

	  return 0.02*kab; //This 0.02 is not justified, but simulation experiment shows that heat transfer is too faster and I just want to make it slower ----> This is old way, not correct! ---> See the following for the correct way
	                  //Some thing that not sure currently is
	                  //1) JJ Monaghan's SPH discretize of second order derivative for 2D/3D is based on 1D Talor's series expansion...
	                  //2) He transfer from using grad (w_ab) to F_ab is not justified, maybe he is correct, but I am still suspecting...
	                  //3) viscosity in JJ Monaghan's equation depends on h? physically, it should be independent of numerical simulations
	                  //4) He turned off viscosity for departing and turned on viscosity for approaching
	                  //5) there is an additional coefficient, which has no physical significance and just numerical make up
	                  //6) The turbulence turn (turbulence stress), behaves not like a shear stress term---> more like a attractive force...
	                  //7) I might have negative heat conductivity term or superficial large heat conductivity, In my opinion, this comes from problem 1) and 2)
//	  return 0;
	}
}

//function that used to compute turbulent heat conductivity in energy equation if SPH_epsilon turbulence model is adopted
// In JJ Monagha's paper, viscosity is defined as alf*h*c
// h and c should be divided here ---> As the physics viscosity should be independent of h and c
double SPH_epsilon_heat_conductivity(double Cp_ab, double * ds, double *vab, double dab, double hab, double cab)
{
	int i;
	double dotvv=0., dotvr=0., dotrr=0.;
	double kab;

	/*
	 * when ds=[-577.7861  577.7746  -74.8549], vab=[-5.8261    5.9188   90.6049]  ds*vab is very small, this is not reasonable
	 */

	for (i=0; i<DIMENSION; i++)
		dotvr += abs((*(vab+i)) * (*(ds+i)));
	if(dotvr==0.) //To be consistent with artificial viscosity --> But in my opinion, it is not really reasonable!
		return 0.;
	else
	{

	  for (i=0; i<DIMENSION; i++)
		dotrr += (*(ds+i)) * (*(ds+i));

	  for (i=0; i<DIMENSION; i++)
	  	dotvv += (*(vab+i)) * (*(vab+i));

      if (DIMENSION==3)
	    kab=6.0*EPSILON*Cp_ab*dab*dotrr*dotvv/(PRANDTL_NUM*dotvr); //For 3D shear viscosity is 1/6 of h*alf*soundspeed
      else if (DIMENSION==2)
        kab=4.8*EPSILON*Cp_ab*dab*dotrr*dotvv/(PRANDTL_NUM*dotvr); //For 2D shear viscosity is 5/24 of h*alf*soundspeed
      else
    	  cout<< "Dimension is neither 2 or 3, this program is not supposed to handle it!" <<endl;

#ifdef DEBUG
	bool print = true;
	if (print && kab<0)
        cout << "You got a negative heat conductivity, something is wrong!" << endl;
#endif

	 return kab/(hab*cab); //viscosity should be independent of smoothing length and sound speed, this one works well when use classical momentum discretization.

	}//end of else
}

//function that switch brief bucket to a bucket
void switch_brief(BriefBucket * breif_neigh, double * mindom, double * maxdom, double * mindom_o, double * maxdom_o, double bucket_size, double len_scale, Bucket ** buck)
{
	unsigned btkey[KEYLENGTH];
	Key tempbtkey;
	Key neigh_btkeys[NEIGH_SIZE];
	int neigh_proc[NEIGH_SIZE];
	double min_crd[DIMENSION], max_crd[DIMENSION], cent_crd[DIMENSION], normc[DIMENSION];
	double xcrd[2], ycrd[2], zcrd[2];
	double neigh_crd[DIMENSION];
	double elev[4];
	int i, j, k, l;
	unsigned keylen= KEYLENGTH;
	int btflag, has_involved=0, type;
	int bk_index[2*DIMENSION];

	double flat[2*DIMENSION];
	flat[0]=Lx_P[0];
	flat[1]=Lx_P[1];
	flat[2]=Ly_P[0];
	flat[3]=Ly_P[1];
	flat[4]=Lz_P[0];
	flat[5]=Lz_P[1];

	for (j = 0; j < DIMENSION; j++)
		min_crd[j]= *(breif_neigh->get_mincrd () + j);
	for (j = 0; j < DIMENSION; j++)
		max_crd[j]= min_crd[j] + bucket_size;
	for (j = 0; j < DIMENSION; j++)
		cent_crd[j]= 0.5 * (min_crd[j] + max_crd[j]);

	for (j = 0; j <  NEIGH_SIZE; j++)
	    neigh_proc[j]= *(breif_neigh->get_neigh_proc () + j);

	//go through all neighbors and determine their keys based on their location, then add the key into neigh_btkeys
	double bucket_size_half = bucket_size*0.5;
	int count = 0;
	for (i=-1; i<2; i++)
	{
		neigh_crd[0] = cent_crd[0] + i * bucket_size;
		for (j=-1; j<2; j++)
		{
			neigh_crd[1] = cent_crd[1] + j * bucket_size;
			for (k=-1; k<2; k++)
			{
				neigh_crd[2] = cent_crd[2] + k * bucket_size;
				xcrd[0]=neigh_crd[0]- bucket_size_half;
				xcrd[1]=neigh_crd[0]+ bucket_size_half;
				ycrd[0]=neigh_crd[1]- bucket_size_half;
				ycrd[1]=neigh_crd[1]+ bucket_size_half;
				zcrd[0]=neigh_crd[2]- bucket_size_half;
				zcrd[1]=neigh_crd[2]+ bucket_size_half;
				determine_bucket_type (mindom, maxdom, xcrd, ycrd, zcrd, &type, bk_index);
				if (type == 3)//if the bucket is within the domain ---> the so called OVERGROUND buckets
				{
				    for ( l=0; l<DIMENSION; l++)
				        normc[l]=(neigh_crd[l]- *(mindom+l))/(*(maxdom+l)- *(mindom+l));
				    HSFC3d (normc, &keylen, btkey);

				    for (l=0; l<KEYLENGTH; l++)
				    	neigh_btkeys[count].key[l] = btkey[l];

				    count ++;
				}
				else
				{
				    for (l=0; l<KEYLENGTH; l++)
						neigh_btkeys[count].key[l] = 0;
				    count ++;
				}
			}
		}
	}

	for (i=0; i<4; i++)
	   elev[i]= mindom_o[2]; //the boundary on the ground --> This is only for flat ground.

	xcrd[0]=min_crd[0];
	ycrd[0]=min_crd[1];
	zcrd[0]=min_crd[2];
	xcrd[1]=max_crd[0];
	ycrd[1]=max_crd[1];
	zcrd[1]=max_crd[2];
	determine_bucket_type (mindom_o, maxdom_o, xcrd, ycrd, zcrd, &type, bk_index);
	switch (type)
	{
	  case 1:
	    btflag = UNDERGROUND;
	    break;
	  case 2:
	    btflag = MIXED;
	    if (bk_index[4] == -1)
	       for (j = 0; j < 4; j++)
	           elev[j] = elev[j] / len_scale;
	    break;
	  case 3:
	    btflag = OVERGROUND;
	    break;
	  case 4:
	    btflag = PRESS_BC;
	    break;
	  default:
	    fprintf(stderr,"ERROR: Unknown buckettype flag.\nCheck the preoprocessor\n");
	    exit(1);
	}

	tempbtkey=breif_neigh->getKey ();
	for (j = 0; j < KEYLENGTH; j++)
	    btkey[j] = tempbtkey.key[j];

	int myid = breif_neigh->get_myprocess();
	*buck = new Bucket(btkey, min_crd, max_crd, btflag, elev,
	                myid, neigh_proc, neigh_btkeys, bk_index, flat);
	(*buck)->mark_active();
	(*buck)->set_has_potential_involved ((bool) has_involved); //Initially, all buckets "contains" the initial domain should be has_potential_involved.
}

// Read file to matrix
// mat is m row n col
void readFile(string fileName, double ** mat, int m, int n)
{
    // Create streamobject
    ifstream infile;
    infile.open(fileName);

    // Exit if file opening failed
    if (!infile.is_open())
    {
        cerr<<"Opening failed"<<endl;
        exit(1);
    }

    string line;
    *mat = new double[m*n];
    string temp, val;
    int find;
    int i=0;
    while( getline (infile,line) && i<m)
    {
        temp = line;
        for(int j=0; j<n; j++)
        {
            find= temp.find_first_of(" ");
            val = temp.substr(0, find);
            temp = temp.substr(find+1);
            (*mat)[i*n+j] = stof(val);
        }
        i++;
    }
    infile.close();
}

//function that used to determine the type of bucket
/*
 * In this function, a face based manner is used ---> This works well for the case where the domain is a box
 * For arbitrary domain, nodes based manner is better.
 */
bool determine_erupt_buket (double *mincrd, double *maxcrd, double *xcrd, double *ycrd, double *zcrd)
{
	int flag[DIMENSION];
//	int sum = 0;
	int k;
	bool erpt_flag = true;
	int bt[6];

    bt[0]=determine_face_type(xcrd[0],maxcrd[0],mincrd[0]);
    bt[1]=determine_face_type(xcrd[1],maxcrd[0],mincrd[0]);
    bt[2]=determine_face_type(ycrd[0],maxcrd[1],mincrd[1]);
    bt[3]=determine_face_type(ycrd[1],maxcrd[1],mincrd[1]);
    bt[4]=determine_face_type(zcrd[0],maxcrd[2],mincrd[2]);
    bt[5]=determine_face_type(zcrd[1],maxcrd[2],mincrd[2]);

    for (k=0; k<DIMENSION; k++)
    	flag[k]=abs(bt[2*k] + bt[2*k+1]);

    for (k=0; k<DIMENSION; k++)
    	if (flag[k] == 2)
    		erpt_flag = false;

    return erpt_flag;
}


