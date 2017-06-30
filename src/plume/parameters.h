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

//##################################################################################################################################################################################
#if CODE_DIMENSION==3

#ifndef SIMULATE_ASH
/*
 * These parameters are for running jobs in the 2005 paper of SK-3D.
 */
//----------------------------------------------------------------------------------------
////% parameter for domain definition
//
////const double Lx_P[2]={-12000,12000};
////const double Ly_P[2]={-12000,12000};
////const double Lz_P[2]={0, 26000};
//const double Ll_P[DIMENSION]={-12000,-12000, 0};
//const double Lu_P[DIMENSION]={12000,12000, 26000};
//const int N_total_P=2700;  // This one is actually not used in current code!
//
////----------------------------------------------------------------------------------------
////% parameter for phase1 (air), use  to generated initial atmosphere condition
//#if ATMOSPHERE_TYPE==2
//const double g_P=0.;
//#else
//const double g_P=9.81;
//#endif
//
////the following 4 parameters are not independent!, they should satisfy the EOS
//const double Ta0_P=273;
//const double pa0_P=1.013e5;
//const double rhoa0_P=1.29;
//const double Ra_P=287;
//
//const double Cvs_P=1617;/*specific heat*/
//const double Cvg_P=1155;
//const double Cva_P=717;
//const double gamma_P=1.4; /*specific heat ratio for ideal gas */
//const double H1_P=11000; /*height of tropopause*/
//const double H2_P=20000;/*height of straightpause*/
//const double H3_P=100000;/*height of atmosphere*/
//const double miu1_P=6.5/1000;/* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert*/
//const double miu2_P= 2./1000;
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
//const double ng0_P=0.05; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
//const double Uv0_P=0.0; /* velocity in horizontal direction*/
//const double Vv0_P=150; /* velocity in verticle direction*/
////const double Vv0_P=0.00015;  /* velocity in verticle direction, this velocity is used for testing non-eruption*/
//const double Tv0_P=1000;
//const double pv0_P=pa0_P;
//
////Additional properties
//const double Cpg_P=Cpg_P+Rg_P;
//const double Cpa_P=Ra_P+Cva_P;
//const double Cps_P=Cvs_P;

//// The old equation: --> ng0_P is the mass fraction of water (vapor)
////const double Cpv_P=Cps_P*ng0_P + (1-ng0_P)*Cpg_P; /*specific heat of erupted material under constant pressure*/
// The new equation:
//const double Cpv_P=Cps_P*(1-ng0_P) + ng0_P*Cpg_P; /*specific heat of erupted material under constant pressure*/
//
//const double Rv_P=ng0_P*Rg_P;  /*as constant for erupted material*/
//const double Cvv_P=ng0_P*Cvg_P+(1-ng0_P)*Cvs_P; /*specific heat of erupted material*/
//const double gamma_v_P=1+Rv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper*/
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
////const double Mv_P=3.9811e+07;  /*mass flow rate, it is not directly used in simulation.*/
//////const double Mv_P=1.0;   //for test non-erupt condition
////
////const double rv_P= pow ((Mv_P/(rhov_P*Vv0_P*3.1415926)),0.5); /*radius of vent*/
//
////Here is the new way to give parameters --->So it is also OK if sml2 is given as a parameter
//const double rv_P = 140;
//const double Mv_P = rv_P*rv_P*rhov_P*Vv0_P*3.1415926;
//
//const double Pos_v_P[DIMENSION]={0., 0., Ll_P[2] };  /*position of vent at the origin*/
//
//const int num_erupt = 2; /*this parameter should be used to determine total number */
//
//const int num_erupt_perlayer = 10;
//const int num_erupt_particles = 1646; //number of particle in the initial erupt duct
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
//const double CFL_BC_P=0.08; //CLF number to stable fluctuation near the boundary, 1/6.
//
////----------------------------------------------------------------------------------------
////Heat transfer coefficient
//const double lamda_P = 0.03; //using a constant heat conduction coefficient... That for air is 0.024
//
//
////----------------------------------------------------------------------------------------
////Work load for each type of particles ----> the data is obtained by profiling
//const double realp_load = 2.0;
//const double wallp_load = 0.75;
//const double pressp_load = 0.;
//const double eruptp_load = 0.02;

////work load check interval, this value should be optimized to get good balance, for problems with different time scale, this value should be different
//const int balancing_check_int_P = 3;


/*
 * These parameters are for repeat the most recent excise!
 *
 * 1) Consider the weak plume
 */
////----------------------------------------------------------------------------------------
////% parameter for domain definition: Lx_P the lower and up boundary in x direction
////                                   Ly_p the lower and up boundary in y direction
////                                   Lz_p the lower and up boundary in z direction
//
////for weak sml1=15m or sml1=30
//const double Lx_P[2]={-2430,2430};
//const double Ly_P[2]={-2430,2430};
//const double Lz_P[2]={1500, 14000};

////for strong sml1=170 smaller domain
//const double Lx_P[2]={-29580,29580};
//const double Ly_P[2]={-29580,29580};
//const double Lz_P[2]={1500,50000};

////for strong sml1=170 larger domain
//const double Lx_P[2]={-39780,39780};
//const double Ly_P[2]={-39780,39780};
//const double Lz_P[2]={1500,55000};
//const double Ll_P[DIMENSION]={-39780,-39780, 1500};
//const double Lu_P[DIMENSION]={39780,39780, 55000};

//for strong sml1=400
const double Ll_P[DIMENSION]={-40800, -40800, 1500};
const double Lu_P[DIMENSION]={40800, 40800, 55000};

//for strong sml1=300 or sml1=200, sml1=150
//const double Lx_P[2]={-28800,28800};
//const double Ly_P[2]={-28800,28800};
//const double Lz_P[2]={1500,50000};

////larger domain sml1=300 or sml1=200, sml1=150
//const double Lx_P[2]={-39600,39600};
//const double Ly_P[2]={-39600,39600};
//const double Lz_P[2]={1500,55000};

////larger domain sml1=300 or sml1=200, sml1=150
//const double Lx_P[2]={-100800,100800};
//const double Ly_P[2]={-100800,100800};
//const double Lz_P[2]={1500,42000};

////for strong coarse resolution sml1=400
//const double Lx_P[2]={-31200,31200};
//const double Ly_P[2]={-31200,31200};
//const double Lz_P[2]={1500,50000};

////for strong coarse resolution sml1=500
//const double Lx_P[2]={-30000,30000};
//const double Ly_P[2]={-30000,30000};
//const double Lz_P[2]={1500,50000};

//////for strong coarse resolution sml1=600
//const double Lx_P[2]={-28800,28800};
//const double Ly_P[2]={-28800,28800};
//const double Lz_P[2]={1500,50000};


//const int N_total_P=2700;  // This one is actually not used in current code!

//----------------------------------------------------------------------------------------
//% parameter for phase1 (air), use  to generated initial atmosphere condition
#if ATMOSPHERE_TYPE==2
const double g_P=0.;
#else
const double g_P=9.80665;
#endif

//the following 4 parameters are not independent!, they should satisfy the EOS
const double Ta0_P=268.;
const double pa0_P= 0.852321e5;
const double rhoa0_P=1.104;
const double Ra_P=287;

const double Cvs_P=1100;/*specific heat*/
const double Cvg_P=1348;
const double Cva_P=717;
//Additional properties
const double Cpg_P=1810;
const double Cpa_P=1000;
const double Cps_P=Cvs_P;
//Cps_P should be the same as Cvs_P

const double gamma_g_P=Cpg_P/Cvg_P;
const double gamma_a_P=Cpa_P/Cva_P;

const double gamma_P=1.4; /*specific heat ratio for ideal gas ----->I need to double check this one and make sure that I am using the correct value in computation */
const double H1_P=11000; /*height of tropopause -->will not been used!*/
const double H2_P=20000;/*height of straightpause-->will not been used!*/
const double H3_P=100000;/*height of atmosphere-->will not been used!*/
const double miu1_P=6.5/1000; /* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert -->will not been used!*/
const double miu2_P= 2./1000; //-->will not been used!*/
//const double m_ratio_P=1.3;/*ratio between particles mass of phase2 to that of phase1*/

//parameter for pressure approximation--->these parameters are obtained by solving hydro static equation togehter with EOS and atmosphere temperature
const double Ata1_p =1.561024741e-8;
const double Ate1_p =5.258643795;
const double Ata2_p = 1.321758528e+5;
const double Atb2_p = -0.16963e-3;
const double Ata3_p = 4.7248e+03;
const double Atb3_p = 58.5117;
const double AtC3_p = 5642.519912;
const double Ate3_p = -0.1709059897e-1;
const double Atf_P = -g_P/(Ra_P*Ta0_P);
//#endif
//----------------------------------------------------------------------------------------
// %parameter at vent, used to impose eruption condition
// %v: vent
const double Rg_P=462.; /*%gas constant for volcanic gases*/

//for stron
const double ng0_P=0.05; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
//////for weak
//const double ng0_P=0.03; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/

const double Uv0_P=0.0; /* velocity in horizontal direction*/

//for strong
const double Vv0_P=275; /* velocity in verticle direction*/
////for weak
//const double Vv0_P=135; /* velocity in verticle direction*/

//for strong
const double Tv0_P=1053;
////for weak
//const double Tv0_P=1273;

const double pv0_P=84363.4; //assuming pressure-balanced jet

const double Rv_P=ng0_P*Rg_P;  /*as constant for erupted material*/
const double Cvv_P=ng0_P*Cvg_P+(1-ng0_P)*Cvs_P; /*specific heat of erupted material*/

// The old equation: --> ng0_P is the mass fraction of water (vapor)
//const double Cpv_P=Cps_P*ng0_P + (1-ng0_P)*Cpg_P; /*specific heat of erupted material under constant pressure*/
// The new equation:
const double Cpv_P=Cps_P*(1-ng0_P) + ng0_P*Cpg_P; /*specific heat of erupted material under constant pressure*/
const double gamma_v_P=1+Rv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper   ----> double check this to make sure they are consistent*/

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
//for strong
const double Mv_P=1.5e+09;  /*mass flow rate, it is not directly used in simulation.*/
//for weak
//const double Mv_P=1.5e+06;  /*mass flow rate, it is not directly used in simulation.*/


const double rv_P= pow ((Mv_P/(rhov_P*Vv0_P*3.1415926)),0.5); /*radius of vent*/

////Here is the new way to give parameters --->So it is also OK if sml2 is given as a parameter
//const double rv_P = 140;
//const double Mv_P = rv_P*rv_P*rhov_P*Vv0_P*3.1415926;

const double Pos_v_P[DIMENSION]={0., 0., Ll_P[2] };  /*position of vent at the origin*/

const int num_erupt = 2; /*this parameter should be used to determine total number */

const int num_erupt_perlayer = 10;

////for strong
//for strong coarse resolution
//--> 2544 for 10 each direction sml1=400, 825 for 6 each direction, sml1=600, 1628 for 8 each direction sml1=500, 1908 for sml1=300, 10 each direction (sml2=141.5),  900 for sml1=141.5=sml2, 1764 for sml1 =100, sml2 =101.4
//----> 7056 for sml1=50, sml2=50.5, 3493 for sml1=200, sml2=101, 954 for sml1=150; sml2=141.5, 1337 for sml1=150, sml2=101; 1515 for sm1=170 sml2=101; 549 for sml1=400, sml2=236
const int num_erupt_particles = 549; //number of particle in the initial erupt duct
////for weak
//4963 for sml1=30, sml2=5.44 ; 1271 for sml1=15, sml2=6.8
//const int num_erupt_particles = 1271; //number of particle in the initial erupt duct

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
const double eta_smooth_P = 1.2; //We can try different value to get best results

//----------------------------------------------------------------------------------------
//CFL coefficient for time step update
const double CFL_P=0.40;
const double CFL_BC_P=0.08; //CLF number to stable fluctuation near the boundary, 1/6.

//----------------------------------------------------------------------------------------
//Heat transfer coefficient
const double lamda_P = 0.03; //using a constant heat conduction coefficient... That for air is 0.024

//----------------------------------------------------------------------------------------
//Work load for each type of particles ----> the data is obtained by profiling
const double realp_load = 2.0;
const double wallp_load = 0.75;
const double pressp_load = 0.;
const double eruptp_load = 0.02;


//work load check interval, this value should be optimized to get good balance, for problems with different time scale, this value should be different
//For strong
const int balancing_check_int_P = 3;
//for weak
//const int balancing_check_int_P = 1;


/*
 * These parameters are for JPUE of incompressible flow
 *
 */
// parameter for domain definition

////const double Lx_P[2]={-528,528};
////const double Ly_P[2]={-528,528};
////const double Lz_P[2]={0, 6000};
//const double Ll_P[DIMENSION]={-528,-528,0};
//const double Lu_P[DIMENSION]={528,528,6000};
//
////----------------------------------------------------------------------------------------
////% parameter for phase1 (air), use  to generated initial atmosphere condition
//#if ATMOSPHERE_TYPE==2
//const double g_P=0.;
//#else
//const double g_P=9.81;
//#endif
//
////the following 4 parameters are not independent!, they should satisfy the EOS
//const double Ta0_P=273;
//const double pa0_P=1.013e5;
//const double rhoa0_P= 1000;
////const double Ra_P=287;
//const double Ra_P=25182.; //4197*6
//
//const double Cvs_P=1617;/*specific heat*/
//const double Cvg_P=4197;
//const double Cva_P=4197;
////const double gamma_P=1.4; /*specific heat ratio for ideal gas */
//const double gamma_P=7.; /*specific heat ratio for ideal gas */
//const double H1_P=11000; /*height of tropopause ---> will not used */
//const double H2_P=20000;/*height of straightpause  ---> will not used */
//const double H3_P=100000;/*height of atmosphere  ---> will not used */
//const double miu1_P=6.5/1000;/* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert*/
//const double miu2_P= 2./1000;
//const double m_ratio_P=1.3;/*ratio between particles mass of phase2 to that of phase1*/
//
////parameter for pressure approximation--->these parameters are obtained by solving hydro static equation togehter with EOS and atmosphere temperature
//const double Ata1_p =1.561024741e-8;
//const double Ate1_p =5.258643795;
//const double Ata2_p = 1.321758528e+5;
//const double Atb2_p = -0.16963e-3;
//const double Ata3_p = 4.7248e+03;
//const double Atb3_p = 58.5117;
//const double AtC3_p = 5642.519912;
//const double Ate3_p = -0.1709059897e-1;
//
//const double Atf_P = -g_P/(Ra_P*Ta0_P);
//
////----------------------------------------------------------------------------------------
//// %parameter at vent, used to impose eruption condition
//// %v: vent
////const double Rg_P=287.; /*%gas constant for volcanic gases*/
//const double Rg_P=25182.; //4197*6
//const double ng0_P=1.; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
//const double Uv0_P=0.0; /* velocity in horizontal direction*/
//const double Vv0_P=500; /* velocity in verticle direction*/
//const double Tv0_P=273;
//const double pv0_P=pa0_P;
//
////Additional properties
//const double Cpg_P=Cpg_P+Rg_P;
//const double Cpa_P=Ra_P+Cva_P;
//const double Cps_P=Cvs_P;
//
//// The old equation: --> ng0_P is the mass fraction of water (vapor)
////const double Cpv_P=Cps_P*ng0_P + (1-ng0_P)*Cpg_P; /*specific heat of erupted material under constant pressure*/
////The new equation:
//const double Cpv_P=Cps_P*(1-ng0_P) + ng0_P*Cpg_P; /*specific heat of erupted material under constant pressure*/
//
//const double Rv_P=ng0_P*Rg_P;  /*as constant for erupted material*/
//const double Cvv_P=ng0_P*Cvg_P+(1-ng0_P)*Cvs_P; /*specific heat of erupted material*/
//const double gamma_v_P=1+Rv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper*/
//
////const double rhov_P=pv0_P/(Rv_P*Tv0_P); //This is incorrect, I am using EOS of idea gas for mixture of solid and gas with gas mass fraction is only 0.05
////The correct way
////const double Rv_g_P=Rg_P;
////const double pv0_g_P=pv0_P; //I am assuming that solid did not influence the pressure of gas, to check whether this is reasonable: Vs=0.95/1000=0.00095, Vg=0.05/1=0.05, Vg/Vs=0.05/0.001=50, so even for ng0=0.05, this kind of assumption is still not bad
////const double rhov_g_P=pv0_g_P/(Rv_g_P*Tv0_P);
//
////The numerically efficient way of computing
////const double rhov_g_P=pv0_P/(Rg_P*Tv0_P);
////const double rhov_P = rhov_g_P/ng0_P; //just notice that this equation is exactly the same as the "incorrect equation"
//const double rhov_P = rhoa0_P;
//const double lamda_v_P=rhov_P*Rv_P/Cvv_P; /*A new defined parameter which lamda_P=rho_m*Rm/Cm, where, rho_m: is the density of the mixture, while Rm is gas constant for the mixture, Cm specific heat of at constatnt volume.*/
//const double ev0_P=ng0_P*Cvg_P*Tv0_P+(1-ng0_P)*Cvs_P*Tv0_P; /*energy of erupted material*/
//
////--->old code:
////const double Mv_P=5.471084e+07;  /*mass flow rate, it is not directly used in simulation.*/
////const double rv_P= pow ((Mv_P/(rhov_P*Vv0_P*3.1415926)),0.5); /*radius of vent*/
//
////--->new code:
//const double rv_P= 20; /*radius of vent*/
//const double Mv_P=(rv_P*rv_P*3.1415926)*(rhov_P*Vv0_P); /*mass flow rate, it is not directly used in simulation.*/
//
//const double Pos_v_P[DIMENSION]={0., 0., Ll_P[2] };  /*position of vent at the origin*/
//
//const int num_erupt = 2; /*this parameter should be used to determine total number */
//
//const int num_erupt_perlayer = 10;
//const int num_erupt_particles = 1800; //number of particle in the initial erupt duct
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
//const double CFL_P=0.45;
//const double CFL_BC_P=0.08; //CLF number to stable fluctuation near the boundary, 1/6.
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
////work load check interval, this value should be optimized to get good balance, for problems with different time scale, this value should be different
//const int balancing_check_int_P = 1;

/*
 * Parameters for running RP1 with sml1=200
 * In this simulation, there was fluctuation near the boundary
 *
 */

//----------------------------------------------------------------------------------------
//% parameter for domain definition

//const double Lx_P[2]={-12000,12000};
//const double Ly_P[2]={-12000,12000};
//const double Lz_P[2]={0, 23000};

//const double Ll_P[DIMENSION]={-12000,-12000, 0};
//const double Lu_P[DIMENSION]={12000,12000, 23000};
//
//const int N_total_P=2700;  // This one is actually not used in current code!
//
////----------------------------------------------------------------------------------------
////% parameter for phase1 (air), use  to generated initial atmosphere condition
//#if ATMOSPHERE_TYPE==2
//const double g_P=0.;
//#else
//const double g_P=9.81;
//#endif
//
////the following 4 parameters are not independent!, they should satisfy the EOS
//const double Ta0_P=273;
//const double pa0_P=1.013e5;
//const double rhoa0_P=1.29;
//const double Ra_P=287;
//
//const double Cvs_P=1617;/*specific heat*/
//const double Cvg_P=1155;
//const double Cva_P=717;
//const double gamma_P=1.4; /*specific heat ratio for ideal gas */
//const double H1_P=11000; /*height of tropopause*/
//const double H2_P=20000;/*height of straightpause*/
//const double H3_P=100000;/*height of atmosphere*/
//const double miu1_P=6.5/1000;/* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert*/
//const double miu2_P= 2./1000;
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
//const double ng0_P=0.05; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
//const double Uv0_P=0.0; /* velocity in horizontal direction*/
//const double Vv0_P=150; /* velocity in verticle direction*/
////const double Vv0_P=0.00015;  /* velocity in verticle direction, this velocity is used for testing non-eruption*/
//const double Tv0_P=1000;
//const double pv0_P=pa0_P;
//
////Additional properties
//const double Cpg_P=Cpg_P+Rg_P;
//const double Cpa_P=Ra_P+Cva_P;
//const double Cps_P=Cvs_P;

//// The old equation: --> ng0_P is the mass fraction of water (vapor)
////const double Cpv_P=Cps_P*ng0_P + (1-ng0_P)*Cpg_P; /*specific heat of erupted material under constant pressure*/
// The new equation:
//const double Cpv_P=Cps_P*(1-ng0_P) + ng0_P*Cpg_P; /*specific heat of erupted material under constant pressure*/
//
//const double Rv_P=ng0_P*Rg_P;  /*as constant for erupted material*/
//const double Cvv_P=ng0_P*Cvg_P+(1-ng0_P)*Cvs_P; /*specific heat of erupted material*/
//const double gamma_v_P=1+Rv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper*/
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
//const double Mv_P=3.9811e+07;  /*mass flow rate, it is not directly used in simulation.*/
////const double Mv_P=1.0;   //for test non-erupt condition
//
//const double rv_P= pow ((Mv_P/(rhov_P*Vv0_P*3.1415926)),0.5); /*radius of vent*/
//
//const double Pos_v_P[DIMENSION]={0., 0., Ll_P[2] };  /*position of vent at the origin*/
//
//const int num_erupt = 2; /*this parameter should be used to determine total number */
//
//const int num_erupt_perlayer = 10;
//const int num_erupt_particles = 1800; // 1800 for sml1 =200, 3600 for sml1=400 number of particle in the initial erupt duct
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
//const double CFL_P=0.20;
//const double CFL_BC_P=0.08; //CLF number to stable fluctuation near the boundary, 1/6.
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

////work load check interval, this value should be optimized to get good balance, for problems with different time scale, this value should be different
//const int balancing_check_int_P = 3;



///*
// * Parameters for running plume simulation: Grimsvotn
// *
// */
//
////----------------------------------------------------------------------------------------
////% parameter for domain definition
//
//////const double Lx_P[2]={-28800,28800};
//////const double Ly_P[2]={-28800,28800};
//////const double Lz_P[2]={54, 20000};
//const double Ll_P[DIMENSION]={-28800,-28800, 54};
//const double Lu_P[DIMENSION]={28800,28800, 20000};
////
////const double Lx_P[2]={-8800,8800};
////const double Ly_P[2]={-8800,8800};
////const double Lz_P[2]={54, 20000};
//
//
//const int N_total_P=2700;  // This one is actually not used in current code!
//
////----------------------------------------------------------------------------------------
////% parameter for phase1 (air), use  to generated initial atmosphere condition
//#if ATMOSPHERE_TYPE==2
//const double g_P=0.;
//#else
//const double g_P=9.81;
//#endif
//
////the following 4 parameters are not independent!, they should satisfy the EOS
//const double Ta0_P=277.95;
//const double pa0_P=100100;
//const double rhoa0_P=1.26921570 ;
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
//const double Ata1_p =1.561024741e-8;
//const double Ate1_p =5.258643795;
//const double Ata2_p = 1.321758528e+5;
//const double Atb2_p = -0.16963e-3;
//const double Ata3_p = 4.7248e+03;
//const double Atb3_p = 58.5117;
//const double AtC3_p = 5642.519912;
//const double Ate3_p = -0.1709059897e-1;
//const double Atf_P = -g_P/(Ra_P*Ta0_P);
//
//
////----------------------------------------------------------------------------------------
//// %parameter at vent, used to impose eruption condition
//// %v: vent
//const double Rg_P=462.; /*%gas constant for volcanic gases*/
//const double ng0_P=0.01; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/
//const double Uv0_P=0.0; /* velocity in horizontal direction*/
//const double Vv0_P=50; /* velocity in verticle direction*/
//const double Tv0_P=1200;
//const double pv0_P=pa0_P;
//
////Additional properties
////const double Cpg_P=Cpg_P+Rg_P;
////const double Cpa_P=Ra_P+Cva_P;
////const double Cps_P=Cvs_P;

//// The old equation: --> ng0_P is the mass fraction of water (vapor)
//// const double Cpv_P=Cps_P*ng0_P + (1-ng0_P)*Cpg_P; /*specific heat of erupted material under constant pressure*/
// The new equation:
//const double Cpv_P=Cps_P*(1-ng0_P) + ng0_P*Cpg_P; /*specific heat of erupted material under constant pressure*/
//
//const double Rv_P=ng0_P*Rg_P;  /*as constant for erupted material*/
//const double Cvv_P=ng0_P*Cvg_P+(1-ng0_P)*Cvs_P; /*specific heat of erupted material*/
//const double gamma_v_P=1+Rv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper*/
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
//const double Mv_P=1772187.67219867;  /*mass flow rate, it is not directly used in simulation.*/
//
//const double rv_P= pow ((Mv_P/(rhov_P*Vv0_P*3.1415926)),0.5); /*radius of vent*/
//
//const double Pos_v_P[DIMENSION]={0., 0., Ll_P[2]};  /*position of vent at the origin*/
//
//const int num_erupt = 2; /*this parameter should be used to determine total number */
//
//const int num_erupt_perlayer = 10;
//const int num_erupt_particles = 778;
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
//const double CFL_P=0.20;
//const double CFL_BC_P=0.08; //CLF number to stable fluctuation near the boundary, 1/6.
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
////work load check interval, this value should be optimized to get good balance, for problems with different time scale, this value should be different
//const int balancing_check_int_P = 1;

#else //ASH_SIMULATION
/*
 * Parameters for umbrella simulation   ---> based on Pinatubo eruption 1990
 * These parameter are used for initial test.
 */

//larger domain sml1=300 or sml1=200, sml1=150
//const double Lx_P[2]={-109600,109600};
//const double Ly_P[2]={-109600,109600};
//const double Lz_P[2]={1500,50000};

const double Ll_P[DIMENSION]={-109600,-109600, 1500};
const double Lu_P[DIMENSION]={109600,109600, 50000};

//----------------------------------------------------------------------------------------
//% parameter for phase1 (air), use  to generated initial atmosphere condition
#if ATMOSPHERE_TYPE==2
const double g_P=0.;
#else
const double g_P=9.80665;
#endif

//the following 4 parameters are not independent!, they should satisfy the EOS
const double Ta0_P=268.; //Should be temperature on the ground
const double pa0_P= 0.852321e5; //Should be pressure on the ground
const double rhoa0_P=1.104;  //density of air on the ground
const double Ra_P=287;   //gas constant of pure air --> Gas constant of influx mixture will based on influx mass fraction and initial mass fraction of vapor in erupted material (ng0)
                                                  //--> gas constant of any mxiture will be computed based on properties of pure material

const double Cvs_P=1100;/*specific heat*/
const double Cvg_P=1348;
const double Cva_P=717;
//Additional properties
const double Cpg_P=1810;
const double Cpa_P=1000;
const double Cps_P=Cvs_P;//Cps_P should be the same as Cvs_P

const double gamma_g_P=Cpg_P/Cvg_P;
const double gamma_a_P=Cpa_P/Cva_P;

const double gamma_P=1.4; /*specific heat ratio for ideal gas ----->I need to double check this one and make sure that I am using the correct value in computation */
const double H1_P=11000; /*height of tropopause -->will not been used!*/
const double H2_P=20000;/*height of straightpause-->will not been used!*/
const double H3_P=100000;/*height of atmosphere-->will not been used!*/
const double miu1_P=6.5/1000; /* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert -->will not been used!*/
const double miu2_P= 2./1000; //-->will not been used!*/
//const double m_ratio_P=1.3;/*ratio between particles mass of phase2 to that of phase1*/

//parameter for pressure approximation--->these parameters are obtained by solving hydro static equation togehter with EOS and atmosphere temperature
const double Ata1_p =1.561024741e-8;
const double Ate1_p =5.258643795;
const double Ata2_p = 1.321758528e+5;
const double Atb2_p = -0.16963e-3;
const double Ata3_p = 4.7248e+03;
const double Atb3_p = 58.5117;
const double AtC3_p = 5642.519912;
const double Ate3_p = -0.1709059897e-1;
const double Atf_P = -g_P/(Ra_P*Ta0_P);
//#endif
//----------------------------------------------------------------------------------------
// %parameter at vent, used to impose eruption condition
// %v: vent
const double Rg_P=462.; /*%gas constant for volcanic gases*/

//for stron
const double ng0_P=0.05; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)   ---> Please be aware of that currently, phase 2 is a mixture erupted material and air, another initial parameter: mssfrc_0 will needed also.*/
//////for weak
//const double ng0_P=0.03; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/

//radial influx velocity
 const double vel0_P=50.0;

//temperature of influx material
const double Tv0_P=243; //unit K

//Initial pressure ---> assume pressure balance
const double pv0_P=4363.4;

//Initial mass fraction of erupted material, for umbrella simulation, the mass fraction of erupted material should be smaller than 1.
const double msfc0_P=0.08;

const double Rv_P=ng0_P*Rg_P*msfc0_P + (1-msfc0_P)*Ra_P;  //Not used in current code
const double Cvv_P=msfc0_P*ng0_P*Cvg_P+msfc0_P*(1-ng0_P)*Cvs_P + (1-msfc0_P)* Cva_P; /*specific heat of erupted material*/

const double Cpv_P=Cps_P*(1-ng0_P)*msfc0_P* + msfc0_P*ng0_P*Cpg_P + (1-msfc0_P)*Cpa_P; /*specific heat of erupted material under constant pressure*/
const double gamma_v_P=Cpv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper   ----> double check this to make sure they are consistent*/

//The numerically efficient way of computing
const double rhov_g_P=pv0_P/(Rg_P*Tv0_P); //density of mixture of air and vapor
const double rhov_P = rhov_g_P/(ng0_P*msfc0_P + 1-msfc0_P); //just notice that this equation is exactly the same as the "incorrect equation"

const double lamda_v_P=rhov_P*Rv_P/Cvv_P; /*A new defined parameter which lamda_P=rho_m*Rm/Cm, where, rho_m: is the density of the mixture, while Rm is gas constant for the mixture, Cm specific heat of at constatnt volume.*/ //----> Not used here
const double ev0_P=msfc0_P*ng0_P*Cvg_P*Tv0_P+(1-ng0_P)*Cvs_P*Tv0_P*msfc0_P + (1-msfc0_P)*Cva_P*Tv0_P; /*energy of erupted material*/

//This is the old way to give parameters
//for strong
const double Mv_P=1.5e+09;  /*mass flow rate, it is not directly used in simulation.*/
//for weak
//const double Mv_P=1.5e+06;  /*mass flow rate, it is not directly used in simulation.*/

const double r_out_P= 10000.0; /*outside radius of influx particle adding */
const double r_in_P= 5000.0; /*outside radius of influx particle adding */
const double h_top_P= 22000.0; /*Top height of the plume */
const double h_bot_P= 17000.0; /*bottom height of the plume */

////Here is the new way to give parameters --->So it is also OK if sml2 is given as a parameter
//const double rv_P = 140;
//const double Mv_P = rv_P*rv_P*rhov_P*Vv0_P*3.1415926;

const double Pos_v_P[DIMENSION]={0., 0., Ll_P[2] };  /*position of vent at the origin*/

//const int num_erupt = 2; /*this parameter should be used to determine total number */
//
//const int num_erupt_perlayer = 10;
//
//////for strong
////for strong coarse resolution
////--> 2544 for 10 each direction sml1=400, 825 for 6 each direction, sml1=600, 1628 for 8 each direction sml1=500, 1908 for sml1=300, 10 each direction (sml2=141.5),  900 for sml1=141.5=sml2, 1764 for sml1 =100, sml2 =101.4
////----> 7056 for sml1=50, sml2=50.5, 3493 for sml1=200, sml2=101, 954 for sml1=150; sml2=141.5, 1337 for sml1=150, sml2=101; 1515 for sm1=170 sml2=101;
//const int num_erupt_particles = 3493; //number of particle in the initial erupt duct
//////for weak
////4963 for sml1=30, sml2=5.44 ; 1271 for sml1=15, sml2=6.8
////const int num_erupt_particles = 1271; //number of particle in the initial erupt duct

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
const double eta_smooth_P = 1.2; //We can try different value to get best results

//----------------------------------------------------------------------------------------
//CFL coefficient for time step update
const double CFL_P=0.4;
const double CFL_BC_P=0.08; //CLF number to stable fluctuation near the boundary, 1/6.

//----------------------------------------------------------------------------------------
//Heat transfer coefficient
const double lamda_P = 0.03; //using a constant heat conduction coefficient... That for air is 0.024

//----------------------------------------------------------------------------------------
//Work load for each type of particles ----> the data is obtained by profiling
const double realp_load = 2.0;
const double wallp_load = 0.75;
const double pressp_load = 0.;
const double eruptp_load = 0.02;


//work load check interval, this value should be optimized to get good balance, for problems with different time scale, this value should be different
//For strong
const int balancing_check_int_P = 30.0;
//for weak
//const int balancing_check_int_P = 1;
#endif /*SIMULATE_ASH*/
//CODE_DIMENSION==3
//##################################################################################################################################################################################
#elif CODE_DIMENSION==1

/*
 * These parameters are for 1D shock tube simulation
 *
 */
const int Nb_P=15;
#if (SHOCK_TUBE_TESTS==0 || SHOCK_TUBE_TESTS==1)
//For Sod shock tube
const double Ll_P[DIMENSION]={-0.4};
const double Lu_P[DIMENSION]={0.4};
const int Nnsrp_P=20; //No saving real particles
#elif SHOCK_TUBE_TESTS==2
//for shu-Osher problem, see paper: assessment of localized artificial diffusive scheme for large-eddy simulation of compressible turbulent flow.
const double Ll_P[DIMENSION]={-5.0};
const double Lu_P[DIMENSION]={5.0};
const int Nnsrp_P=50;
#endif

//----------------------------------------------------------------------------------------
//% parameter for phase1 (air), use  to generated initial atmosphere condition
#if ATMOSPHERE_TYPE==2
const double g_P=0.;
#else
const double g_P=9.80665;
#endif

//the following 4 parameters are not independent!, they should satisfy the EOS
const double Ta0_P=268.; //Should be temperature on the ground
const double pa0_P= 1.0125e5; //Should be pressure on the ground
const double rhoa0_P=1.104;  //density of air on the ground
const double Ra_P=287;   //gas constant of pure air --> Gas constant of influx mixture will based on influx mass fraction and initial mass fraction of vapor in erupted material (ng0)
                                                  //--> gas constant of any mxiture will be computed based on properties of pure material

const double Cvs_P=1100;/*specific heat*/
const double Cvg_P=1348;
const double Cva_P=717;
//Additional properties
const double Cpg_P=1810;
const double Cpa_P=1000;
const double Cps_P=Cvs_P;//Cps_P should be the same as Cvs_P

const double gamma_g_P=Cpg_P/Cvg_P;
const double gamma_a_P=Cpa_P/Cva_P;

const double gamma_P=1.4; /*specific heat ratio for ideal gas ----->I need to double check this one and make sure that I am using the correct value in computation */
const double H1_P=11000; /*height of tropopause -->will not been used!*/
const double H2_P=20000;/*height of straightpause-->will not been used!*/
const double H3_P=100000;/*height of atmosphere-->will not been used!*/
const double miu1_P=6.5/1000; /* per meter, in original paper, it is per kilometer, any way, divided by 1000 is for unit convert -->will not been used!*/
const double miu2_P= 2./1000; //-->will not been used!*/
//const double m_ratio_P=1.3;/*ratio between particles mass of phase2 to that of phase1*/

//parameter for pressure approximation--->these parameters are obtained by solving hydro static equation togehter with EOS and atmosphere temperature
const double Ata1_p =1.561024741e-8;
const double Ate1_p =5.258643795;
const double Ata2_p = 1.321758528e+5;
const double Atb2_p = -0.16963e-3;
const double Ata3_p = 4.7248e+03;
const double Atb3_p = 58.5117;
const double AtC3_p = 5642.519912;
const double Ate3_p = -0.1709059897e-1;
const double Atf_P = -g_P/(Ra_P*Ta0_P);
//#endif
//----------------------------------------------------------------------------------------
// %parameter at vent, used to impose eruption condition
// %v: vent
const double Rg_P=462.; /*%gas constant for volcanic gases*/

//for stron
const double ng0_P=0.05; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)   ---> Please be aware of that currently, phase 2 is a mixture erupted material and air, another initial parameter: mssfrc_0 will needed also.*/
//////for weak
//const double ng0_P=0.03; /* initial mass fraction of volcanic gas: (mass of volcanic gas)/(total mass of erupted material)*/

//radial influx velocity
const double vel0_P=50.0;

//temperature of influx material
const double Tv0_P=243; //unit K

//Initial pressure ---> assume pressure balance
const double pv0_P=4363.4;

//Initial mass fraction of erupted material, for umbrella simulation, the mass fraction of erupted material should be smaller than 1.
const double msfc0_P=0.1;

const double Rv_P=ng0_P*Rg_P*msfc0_P + (1-msfc0_P)*Ra_P;  //Not used in current code
const double Cvv_P=msfc0_P*ng0_P*Cvg_P+msfc0_P*(1-ng0_P)*Cvs_P + (1-msfc0_P)* Cva_P; /*specific heat of erupted material*/

const double Cpv_P=Cps_P*(1-ng0_P)*msfc0_P* + msfc0_P*ng0_P*Cpg_P + (1-msfc0_P)*Cpa_P; /*specific heat of erupted material under constant pressure*/
const double gamma_v_P=Cpv_P/Cvv_P; /*gamma of the mixture, using equation (6) in main reference paper   ----> double check this to make sure they are consistent*/

//The numerically efficient way of computing
const double rhov_g_P=pv0_P/(Rg_P*Tv0_P); //density of mixture of air and vapor
const double rhov_P = rhov_g_P/(ng0_P*msfc0_P + 1-msfc0_P); //just notice that this equation is exactly the same as the "incorrect equation"

const double lamda_v_P=rhov_P*Rv_P/Cvv_P; /*A new defined parameter which lamda_P=rho_m*Rm/Cm, where, rho_m: is the density of the mixture, while Rm is gas constant for the mixture, Cm specific heat of at constatnt volume.*/ //----> Not used here
const double ev0_P=msfc0_P*ng0_P*Cvg_P*Tv0_P+(1-ng0_P)*Cvs_P*Tv0_P*msfc0_P + (1-msfc0_P)*Cva_P*Tv0_P; /*energy of erupted material*/

//This is the old way to give parameters
//for strong
const double Mv_P=1.5e+09;  /*mass flow rate, it is not directly used in simulation.*/
//for weak
//const double Mv_P=1.5e+06;  /*mass flow rate, it is not directly used in simulation.*/

const double r_out_P= 10000.0; /*outside radius of influx particle adding */
const double r_in_P= 5000.0; /*outside radius of influx particle adding */
const double h_top_P= 25000.0; /*Top height of the plume */
const double h_bot_P= 20000.0; /*bottom height of the plume */

////Here is the new way to give parameters --->So it is also OK if sml2 is given as a parameter
//const double rv_P = 140;
//const double Mv_P = rv_P*rv_P*rhov_P*Vv0_P*3.1415926;

//const double Pos_v_P[DIMENSION]={0., 0., Ll_P[2] };  /*position of vent at the origin*/

//const int num_erupt = 2; /*this parameter should be used to determine total number */
//
//const int num_erupt_perlayer = 10;
//
//////for strong
////for strong coarse resolution
////--> 2544 for 10 each direction sml1=400, 825 for 6 each direction, sml1=600, 1628 for 8 each direction sml1=500, 1908 for sml1=300, 10 each direction (sml2=141.5),  900 for sml1=141.5=sml2, 1764 for sml1 =100, sml2 =101.4
////----> 7056 for sml1=50, sml2=50.5, 3493 for sml1=200, sml2=101, 954 for sml1=150; sml2=141.5, 1337 for sml1=150, sml2=101; 1515 for sm1=170 sml2=101;
//const int num_erupt_particles = 3493; //number of particle in the initial erupt duct
//////for weak
////4963 for sml1=30, sml2=5.44 ; 1271 for sml1=15, sml2=6.8
////const int num_erupt_particles = 1271; //number of particle in the initial erupt duct

//----------------------------------------------------------------------------------------
// For artificial viscosity
const double alf_P=1.0;
const double beta_P=2.0;
const double ata_P = 0.01;
//----------------------------------------------------------------------------------------
// for variable smooth length
const int num_loop_P=2;
const double thresh_P=1e-5;
const double C_smooth_P = 2.0;
const double eta_smooth_P = 1.2; //We can try different value to get best results

//----------------------------------------------------------------------------------------
//CFL coefficient for time step update
const double CFL_P=0.1;
const double CFL_BC_P=0.08; //CLF number to stable fluctuation near the boundary, 1/6.

//----------------------------------------------------------------------------------------
//Heat transfer coefficient
const double lamda_P = 0.03; //using a constant heat conduction coefficient... That for air is 0.024

//----------------------------------------------------------------------------------------
//Work load for each type of particles ----> the data is obtained by profiling
const double realp_load = 2.0;
const double wallp_load = 0.75;
const double pressp_load = 0.;
const double eruptp_load = 0.02;


//work load check interval, this value should be optimized to get good balance, for problems with different time scale, this value should be different
//For strong
const int balancing_check_int_P = 3;

#endif //CODE_DIMENSION
//##################################################################################################################################################################################


#endif /* PARAMETERS_H_ */
