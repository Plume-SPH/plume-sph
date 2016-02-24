/*
 * sph_header.h
 *
 *  Created on: Feb 17, 2015
 *      Author: zhixuanc
 */

#ifndef SPH_HEADER_H
#define SPH_HEADER_H
#  include <cmath>
#  include <ctime>
#  include <cstdlib>
#  include <iostream>
using namespace std;

#  include <hashtab.h>
#  include <properties.h>
#  include <particle.h>
#  include <IndMap.h>
#  include "constant.h"

/*
 *  Walltime
 */
inline double
walltime()
{
  return (double) clock() / (double) CLOCKS_PER_SEC;
}

/*
 *  sgn
 */
inline double
sgn(double x)
{
  if (x < 0)
    return -1.0;
  else
    return 1.0;
}

/*!
 * in_support() checks if the neigboring particles in within support or not.
 */
inline bool
in_support(
            //! cartesian projections of inter particle distance
            const double *dx,
            //! support
            double h)
{
  for (int i = 0; i < DIMENSION; i++)
    if (abs(*(dx + i)) > h)
      return false;
  return true;
};

/*!
 * weight() computes and returns the weight of a point, based on its
 * location. Gussian weight function is used.
 */
double weight(
               //!  s = (x_i-x_j) / h
               double *s,
               //! smoothing length
               double h);

/*!
 * d_weight computes and returns the derivative of Gaussian weight
 * function and returns its value. Similar to weight(), function is
 * overloaded.
 */
double d_weight(
                 //! sx = (x-xj)/h
                 double * s,
                 //! smoothing length
                 double h,
                 //! { 0, 1, 2 }: direction of differentiation
                 int dir);
//Function to compute F form r_ab and W_ab
//The relationship is: F r_ab = W_ab
double compute_F(
		double* ,  //dwdx
		double*    //dx
		);

//These basic algebra code is not used in SPH code, but will be used in GSPH code
//The recommended option is to use libraries for such basic things, such as blas, Lapack, Petsc

/*!
 * rotates the vector \f$\mathbf{u}\f$ along unit vector \em cosines.
 * It involves rotation of velocity vector and rotation of stress-tensor.
 * It is used to rotate velocity etc. to inter-particle local coordinate system
 */
void rotate(
             //! the vector of state variables
             double * u,
             //! consines of rotation
             double * cosines);

/*!
 * Reflects vector \f$\overline{u}\f$ in plane with \f$\hat{n}\f$
 * normal
 */
void reflect(
              //! incident vector
              double *,
              //! reflection
              double *,
              //! normal to plane of reflection
              double *);

/*!
 * dot product of vector A and vector B
 */
double dot(
            //! Vector A
            double *,
            //! Vecotr B
            double *);

/*!
 *  Rotation of Vector.
 *  For any higher order multiplication BLAS etc is preferred
 */
void matrix_vec_mult(
                      //! Matix A ( assumed square matrix )
                      double *A,
                      // matirx leading dim, i.e. N
                      int N,
                      //! Vector b
                      double *b,
                      //! c = A*b
                      double *c);

void lsq_surf3 (
                //! [x(i), x(i+1)]
                double * x,
                //! [y(j), y(j+1)]
                double * y,
                //! [z(i,j), ..., z(i,j+1)]
                double * z,
                //! polynomial constants
                double * poly);


void poly_surf (
                //! [x(i), x(i+1)]
                double * x,
                //! [y(j), y(j+1)]
                double * y,
                //! [z(i,j), ..., z(i,j+1)]
                double * z,
                //! polynomial constants
                double * poly);

//general file for artificial viscosity
double art_vis (
		// pointer to particles i
		Particle *,
		// pointer to particles j
		Particle *
		);

//overloading file for artificial viscosity
double art_vis (
		// rhoab
		double ,
		// sndspdab
		double ,
		// rab
		double [3],
		//vab
		double [3],
		//rsqab
		double ,
		//h
		double
		);

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) and internal energy
void air_prop (
		double *, //coordinate of particle
		double *, //internal energy of particle
		double *, //pressure
		double *  //density
        );

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
void air_prop (
		double *, //coordinate of particle
		double *, //range of the space which is occupied by the particle [xmin, xmax, ymin, ymax, zmin, zmax]
		double *, //internal energy of particle
		double *, //pressure
		double *, //density
		double *  //particle mass
		);

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) and internal energy
//Based on realistic atmosphere model
void air_prop_realistic (
		double *, //coordinate of particle
		double *, //internal energy of particle
		double *, //pressure
		double *  //density
        );

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Based on realistic atmosphere model
void air_prop_realistic (
		double *, //coordinate of particle
		double *, //range of the space which is occupied by the particle [xmin, xmax, ymin, ymax, zmin, zmax]
		double *, //internal energy of particle
		double *, //pressure
		double *, //density
		double *  //particle mass
		);

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) and internal energy
//Based on hydro-static equation dp/dz=-rho*g ---> This will give a less realistic initial atmosphere, but more consistent with current model.
void air_prop_hydro (
		double *, //coordinate of particle
		double *, //internal energy of particle
		double *, //pressure
		double *  //density
        );

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Based on hydro-static equation dp/dz=-rho*g ---> This will give a less realistic initial atmosphere, but more consistent with current model.
void air_prop_hydro (
		double *, //coordinate of particle
		double *, //range of the space which is occupied by the particle [xmin, xmax, ymin, ymax, zmin, zmax]
		double *, //internal energy of particle
		double *, //pressure
		double *, //density
		double *  //particle mass
		);

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy
//Function that will give a "uniform" atmosphere ---> which is not realistic for sure
//Used for code testing
void air_prop_uniform (
		double *, //coordinate of particle
		double *, //internal energy of particle
		double *, //pressure
		double *  //density
        );

//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Function that will give a "uniform" atmosphere ---> which is not realistic for sure
//Used for code testing
void air_prop_uniform (
		double *, //coordinate of particle
		double *, //range of the space which is occupied by the particle [xmin, xmax, ymin, ymax, zmin, zmax]
		double *, //internal energy of particle
		double *, //pressure
		double *, //density
		double *  //particle mass
		);
//function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy
//Function that will give a "uniform temperature" atmosphere ---> which is not realistic for sure
////Used for code testing
//Gravity coefficient will be set to 9.81 -->This is different from uniform environment.
void air_prop_uniformT (
                double *, //coordinate of particle
                double *, //internal energy of particle
                double *, //pressure
                double *  //density
        );
//overloading of function that used to determine the property of air: density, pressure, (temperature not explicitly output) , internal energy and mass of particles
//Function that will give a "uniform temperature" atmosphere ---> which is not realistic for sure
//Used for code testing
//Gravity coefficient will be set to 9.81 -->This is different from uniform environment.
void air_prop_uniformT (
                double *, //coordinate of particle
                double *, //range of the space which is occupied by the particle [xmin, xmax, ymin, ymax, zmin, zmax]
                double *, //internal energy of particle
                double *, //pressure
                double *, //density
                double *  //particle mass
                );
//function that used to determine only internal energy based on altitude.
//This function is based on a less realistic model: hydrostatic model
//This function will be used while imposing wall boundary condition
//actually, the way to imposing essential boundary condition for internal energy is not a proper way in bcond,cc
//But temporarily, I just use this not proper way to impose boundary condition.---> need to read more papers on how to imposing essential boundary condition in SPH method.
double air_engr_hydro (
		double * //coordinate
		);

////Function that used Petsc GESRM to solve system of equations to return mass of each particles
//void solve_mass(
//		double **, //matrix A      |
//		double *,  //vector b      | the function is get x by solving Ax=b;
//		double * , //solution x    |
//		int ,      //number of rows, matrices are decomposed by rows
//		int ,     //number of columns, each matrix has the same number of columns
//		int       //process id
//		);

//function determine mass for uniform density field discretized by uniformally distributed particles of only one layer.
void Compute_mass (
		//! number of particles
		int np,
		//! Domain information
		double * dom,
		//density
		double des,
		//mass, out-put
		double *mss,
		//smooth length: out_put
		double *sml
		);

//function determine mass for uniform density field discretized by uniformally distributed particles.
void Compute_mass (
		//! number of particles
		int np,
		//! Domain information x direction
		double * range_x,
		//! Domain information x direction
		double * range_y,
		//! Domain information x direction
		double * range_z,
		//density
		double des,
		//mass, out-put
		double *mss,
		//smooth length: out_put
		double *sml
		);

//function that determine whether particle is within bucket or nor
bool in_bucket(
		double *, //bucket max coordinate
		double *, //bucket min coordinate
		double *  //particle coordinate
		);

//function that used to sort process id to form comm
void id_sort(
		//input 1D array
		int *,
		//number of elements in the array.
		int
		);

//Define a new MPI data type corresponding to IndMap
void create_indmpi_struct();

//Send and receive IndexMap among neighbors
void
exchange_indmap (
		//number of processes
		int,
		//myid
		int,
		//my_comm
		int * ,
		//local array of IndMap, will send out to neighbors
		IndMap *,
		//global array of IndMap, receive data from all neighbors
		IndMap **,
		//number of non-guest particles in local process
		int,
		//number of non-guest particles in whole neighbor comm
		int *
		);

// function that used to determine the velocity based on a parabolic velocity profile

double parabolic_vel(
		double , //R radius of the pipe
		double , //rsq square of distance from given point to center
		double   //maximum velocity
		);

// function that used to determine the velocity based on a uniform velocity profile
double uniform_vel(
		double  //average velocity
		);

//function that determine the velocity profile
double vel_prof(
		double , //R radius of the pipe
		double , //rsq square of distance from given point to center
		double   //average velocity
		);

//function that used to determine the value of face by the face's index;
int determine_face_type (
		double, // coordinate (x, y, z)
		double, // max, max domain
		double // min, min domain
		);

//function that used to determine the type of bucket
void determine_bucket_type (
		double *,  //mincrd,
		double *,  //maxcrd,
		double *,  //xcrd,
		double *,  //ycrd,
		double *,  //zcrd,
		int*,      // type,
		int*       //bt
		);


//function that determines parameters of certain particle
void initial_air (
		Particle * //pi
		);

//function that used to compute the additional term in momentum equation if SPH_epsilon turbulence model is adopted
double SPH_epsilon_mom(
		double*, //vab,
		double   //V_b
		);

//function that used to compute turbulent heat conductivity in energy equation if SPH_epsilon turbulence model is adopted
double SPH_epsilon_heat_conductivity(
		double ,  //Cp_ab,
		double*,  //ds,
		double *  //vab
		);
//function that switch brief bucket to a bucket
void switch_brief(
		BriefBucket *,   // breif_neigh,
		double * ,       //mindom,
		double * ,       //maxdom,
		double * ,       //mindom_o,
		double * ,       //maxdom_o,
		double ,         //bucket_size,
		double ,         //len_scale
		Bucket **         //buck
		);
#endif  //SPH_HEADER_H_
