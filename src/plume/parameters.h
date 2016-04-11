/*
 * parameters.h
 *
 *  Created on: Feb 17, 2015
 *      Author: zhixuanc
 */

// These are input parameters

#include <cmath>
#include <constant.h>

#ifndef PARAMETERS_H
#define PARAMETERS_H

/*
 * These parameters are for running jobs in the 2005 paper.
 */
//----------------------------------------------------------------------------------------
//% parameter for domain definition

const double Lx_P[2]={-12000,12000};
const double Ly_P[2]={-12000,12000};
const double Lz_P[2]={0, 26000};

const int N_total_P=2700;  // This one is actually not used in current code!

//----------------------------------------------------------------------------------------
//% parameter for phase1 (air), use  to generated initial atmosphere condition
#if ATMOSPHERE_TYPE==2
const double g_P=0.;
#else
const double g_P=9.81;
#endif

//the following 4 parameters are not independent!, they should satisfy the EOS
const double Ta0_P=273;
const double pa0_P=1.013e5;
const double rhoa0_P=1.29;
const double Ra_P=287;

const double Cvs_P=1617;/*specific heat*/
const double Cvg_P=1155;
const double Cva_P=717;
const double gamma_P=1.4; /*specific heat ratio for ideal gas */
const double H1_P=11000; /*height of tropopause*/
const double H2_P=20000;/*height of straightpause*/
const double H3_P=100000;/*height of atmosphere*/
const double miu1_P=6.5/1000;/* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert*/
const double miu2_P= 2./1000;
const double m_ratio_P=1.3;/*ratio between particles mass of phase2 to that of phase1*/

//parameter for pressure approximation--->these parameters are obtained by solving hydro static equation togehter with EOS and atmosphere temperature
//#if ATMOSPHERE_TYPE==1
const double Ata1_p =1.561024741e-8;
const double Ate1_p =5.258643795;
const double Ata2_p = 1.321758528e+5;
const double Atb2_p = -0.16963e-3;
const double Ata3_p = 4.7248e+03;
const double Atb3_p = 58.5117;
const double AtC3_p = 5642.519912;
const double Ate3_p = -0.1709059897e-1;
//#elif ATMOSPHERE_TYPE==3
const double Atf_P = -g_P/(Ra_P*Ta0_P);
//#endif
//----------------------------------------------------------------------------------------
// %parameter at vent, used to impose eruption condition
// %v: vent
const double Rg_P=462.; /*%gas constant for volcanic gases*/
const double ng0_P=0.05; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
const double Uv0_P=0.0; /* velocity in horizontal direction*/
const double Vv0_P=150; /* velocity in verticle direction*/
//const double Vv0_P=0.00015;  /* velocity in verticle direction, this velocity is used for testing non-eruption*/
const double Tv0_P=1000;
const double pv0_P=pa0_P;

const double Rv_P=ng0_P*Rg_P;  /*as constant for erupted material*/
const double Cvv_P=ng0_P*Cvg_P+(1-ng0_P)*Cvs_P; /*specific heat of erupted material*/
const double gamma_v_P=1+Rv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper*/

//const double rhov_P=pv0_P/(Rv_P*Tv0_P); //This is incorrect, I am using EOS of idea gas for mixture of solid and gas with gas mass fraction is only 0.05
//The correct way
//const double Rv_g_P=Rg_P;
//const double pv0_g_P=pv0_P; //I am assuming that solid did not influence the pressure of gas, to check whether this is reasonable: Vs=0.95/1000=0.00095, Vg=0.05/1=0.05, Vg/Vs=0.05/0.001=50, so even for ng0=0.05, this kind of assumption is still not bad
//const double rhov_g_P=pv0_g_P/(Rv_g_P*Tv0_P);

//The numerically efficient way of computing
const double rhov_g_P=pv0_P/(Rg_P*Tv0_P);
const double rhov_P = rhov_g_P/ng0_P; //just notice that this equation is exactly the same as the "incorrect equation"

const double lamda_v_P=rhov_P*Rv_P/Cvv_P; /*A new defined parameter which lamda_P=rho_m*Rm/Cm, where, rho_m: is the density of the mixture, while Rm is gas constant for the mixture, Cm specific heat of at constatnt volume.*/
const double ev0_P=ng0_P*Cvg_P*Tv0_P+(1-ng0_P)*Cvs_P*Tv0_P; /*energy of erupted material*/


//This is the old way to give parameters
//const double Mv_P=3.9811e+07;  /*mass flow rate, it is not directly used in simulation.*/
////const double Mv_P=1.0;   //for test non-erupt condition
//
//const double rv_P= pow ((Mv_P/(rhov_P*Vv0_P*3.1415926)),0.5); /*radius of vent*/

//Here is the new way to give parameters --->So it is also OK if sml2 is given as a parameter
const double rv_P = 140;
const double Mv_P = rv_P*rv_P*rhov_P*Vv0_P*3.1415926;

const double Pos_v_P[DIMENSION]={0., 0., 0. };  /*position of vent at the origin*/

const int num_erupt = 2; /*this parameter should be used to determine total number */

const int num_erupt_perlayer = 10;
const int num_erupt_particles = 1646; //number of particle in the initial erupt duct

//----------------------------------------------------------------------------------------
// For artificial viscosity
const double alf_P=0.3;
const double beta_P=0.6;
const double ata_P = 0.01;
//----------------------------------------------------------------------------------------
// for variable smooth length
const int num_loop_P=2;
const double thresh_P=1e-5;
const double C_smooth_P = 2.0;
const double eta_smooth_P = 1.2;

//----------------------------------------------------------------------------------------
//CFL coefficient for time step update
const double CFL_P=0.25;

//----------------------------------------------------------------------------------------
//Heat transfer coefficient
const double lamda_P = 0.03; //using a constant heat conduction coefficient... That for air is 0.024


//----------------------------------------------------------------------------------------
//Work load for each type of particles ----> the data is obtained by profiling
const double realp_load = 2.0;
const double wallp_load = 0.75;
const double pressp_load = 0.;
const double eruptp_load = 0.02;


///*
// * These parameters are for repeat the most recent excise!
// *
// * 1) Consider the weak plume
// */
////----------------------------------------------------------------------------------------
////% parameter for domain definition
//
////for weak
////const double Lx_P[2]={-2430,2430};
////const double Ly_P[2]={-2430,2430};
////const double Lz_P[2]={1500, 14000};
//
////for strong
//const double Lx_P[2]={-28800,28800};
//const double Ly_P[2]={-28800,28800};
//const double Lz_P[2]={1500,50000};
//
//const int N_total_P=2700;  // This one is actually not used in current code!
//
////----------------------------------------------------------------------------------------
////% parameter for phase1 (air), use  to generated initial atmosphere condition
//#if ATMOSPHERE_TYPE==2
//const double g_P=0.;
//#else
//const double g_P=9.80665;
//#endif
//
////the following 4 parameters are not independent!, they should satisfy the EOS
//const double Ta0_P=273; //this will not been used if temperature
//const double pa0_P=1.013e5; //This will not been used also
//const double rhoa0_P=1.29;  //This will not been used also
//const double Ra_P=287;
//
//const double Cvs_P=1100;/*specific heat*/
//const double Cvg_P=1348;
//const double Cva_P=717;
////Additional properties
//const double Cpg_P=1810;
//const double Cpa_P=1000;
//const double Cps_P=Cvs_P;
////Cps_P should be the same as Cvs_P
//
//const double gamma_g_P=Cpg_P/Cvg_P;
//const double gamma_a_P=Cpa_P/Cva_P;
//
//const double gamma_P=1.4; /*specific heat ratio for ideal gas ----->I need to double check this one and make sure that I am using the correct value in computation */
//const double H1_P=11000; /*height of tropopause -->will not been used!*/
//const double H2_P=20000;/*height of straightpause-->will not been used!*/
//const double H3_P=100000;/*height of atmosphere-->will not been used!*/
//const double miu1_P=6.5/1000; /* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert -->will not been used!*/
//const double miu2_P= 2./1000; //-->will not been used!*/
//const double m_ratio_P=1.3;/*ratio between particles mass of phase2 to that of phase1*/
//
////parameter for pressure approximation--->these parameters are obtained by solving hydro static equation togehter with EOS and atmosphere temperature
////#if ATMOSPHERE_TYPE==1
//const double Ata1_p =1.561024741e-8;
//const double Ate1_p =5.258643795;
//const double Ata2_p = 1.321758528e+5;
//const double Atb2_p = -0.16963e-3;
//const double Ata3_p = 4.7248e+03;
//const double Atb3_p = 58.5117;
//const double AtC3_p = 5642.519912;
//const double Ate3_p = -0.1709059897e-1;
////#elif ATMOSPHERE_TYPE==3
//const double Atf_P = -g_P/(Ra_P*Ta0_P);
////#endif
////----------------------------------------------------------------------------------------
//// %parameter at vent, used to impose eruption condition
//// %v: vent
//const double Rg_P=462.; /*%gas constant for volcanic gases*/
//
////for strong
//const double ng0_P=0.05; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
//////for weak
////const double ng0_P=0.03; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
//
//const double Uv0_P=0.0; /* velocity in horizontal direction*/
//
////for strong
//const double Vv0_P=275; /* velocity in verticle direction*/
////for weak
////const double Vv0_P=135; /* velocity in verticle direction*/
//
////for strong
//const double Tv0_P=1053;
//////for weak
////const double Tv0_P=1273;
//
//const double pv0_P=84363.4; //assuming pressure-balanced jet
//
//const double Rv_P=ng0_P*Rg_P;  /*as constant for erupted material*/
//const double Cvv_P=ng0_P*Cvg_P+(1-ng0_P)*Cvs_P; /*specific heat of erupted material*/
//const double gamma_v_P=1+Rv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper   ----> double check this to make sure they are consistent*/
//
////const double rhov_P=pv0_P/(Rv_P*Tv0_P); //This is incorrect, I am using EOS of idea gas for mixture of solid and gas with gas mass fraction is only 0.05
////The correct way
////const double Rv_g_P=Rg_P;
////const double pv0_g_P=pv0_P; //I am assuming that solid did not influence the pressure of gas, to check whether this is reasonable: Vs=0.95/1000=0.00095, Vg=0.05/1=0.05, Vg/Vs=0.05/0.001=50, so even for ng0=0.05, this kind of assumption is still not bad
////const double rhov_g_P=pv0_g_P/(Rv_g_P*Tv0_P);
//
////The numerically efficient way of computing
//const double rhov_g_P=pv0_P/(Rg_P*Tv0_P);
//const double rhov_P = rhov_g_P/ng0_P; //just notice that this equation is exactly the same as the "incorrect equation"
//
//const double lamda_v_P=rhov_P*Rv_P/Cvv_P; /*A new defined parameter which lamda_P=rho_m*Rm/Cm, where, rho_m: is the density of the mixture, while Rm is gas constant for the mixture, Cm specific heat of at constatnt volume.*/
//const double ev0_P=ng0_P*Cvg_P*Tv0_P+(1-ng0_P)*Cvs_P*Tv0_P; /*energy of erupted material*/
//
//
////This is the old way to give parameters
////for strong
//const double Mv_P=1.5e+09;  /*mass flow rate, it is not directly used in simulation.*/
//////for weak
////const double Mv_P=1.5e+06;  /*mass flow rate, it is not directly used in simulation.*/
//
//
//const double rv_P= pow ((Mv_P/(rhov_P*Vv0_P*3.1415926)),0.5); /*radius of vent*/
//
//////Here is the new way to give parameters --->So it is also OK if sml2 is given as a parameter
////const double rv_P = 140;
////const double Mv_P = rv_P*rv_P*rhov_P*Vv0_P*3.1415926;
//
//const double Pos_v_P[DIMENSION]={0., 0., Lz_P[0] };  /*position of vent at the origin*/
//
//const int num_erupt = 2; /*this parameter should be used to determine total number */
//
//const int num_erupt_perlayer = 10;
//
////for strong
//const int num_erupt_particles = 5238; //number of particle in the initial erupt duct
//////for weak
////const int num_erupt_particles = 2541; //4963 for sml2=5.44//number of particle in the initial erupt duct
//
////----------------------------------------------------------------------------------------
//// For artificial viscosity
//const double alf_P=0.3;
//const double beta_P=0.6;
//const double ata_P = 0.01;
////----------------------------------------------------------------------------------------
//// for variable smooth length
//const int num_loop_P=2;
//const double thresh_P=1e-5;
//const double C_smooth_P = 2.0;
//const double eta_smooth_P = 1.2;
//
////----------------------------------------------------------------------------------------
////CFL coefficient for time step update
//const double CFL_P=0.25;
//
////----------------------------------------------------------------------------------------
////Heat transfer coefficient
//const double lamda_P = 0.03; //using a constant heat conduction coefficient... That for air is 0.024
//
////----------------------------------------------------------------------------------------
////Work load for each type of particles ----> the data is obtained by profiling
//const double realp_load = 2.0;
//const double wallp_load = 0.75;
//const double pressp_load = 0.;
//const double eruptp_load = 0.02;
//
#endif /* PARAMETERS_H_ */
