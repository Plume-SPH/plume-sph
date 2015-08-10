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

#include <cmath>
//#include <petsc.h>
using namespace std;

#include <particler.h>
#include <sph_header.h>
#include <constant.h>
#  include <IndMap.h>
#include <parameters.h>

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
const double cutoff = 3.0;
const double t1 = 0.564189583547756;    // 1/sqrt(pi)

// weight function
double
weight(double *s, double h)
{
  int i;
  double norm, expt, w;

  for (i = 0; i < DIMENSION; i++)
    if (abs(s[i]) > cutoff)
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
    if (abs (s[i]) > cutoff)
      return 0;

  tmp = -2 * s[dir] / h;
  dw = tmp * weight (s, h);
  return dw;
}

void
matrix_vec_mult(double *A, int N, double *b, double *c)
{
  register int i, j;
  for (i = 0; i < N; i++)
    c[i] = 0;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      c[i] += A[i * N + j] * b[j];
  return;
}

void
matrix_matrix_mult (double * A,   //! input A : N x P matrix
                    int N,        //! no. of rows of LHS
                    int P,        //! no. of columns of LHS
                    double * B,   //! input B : P x M matrix
                    int M,        //! no. of columns of RHS matrix
                    double * C    //! output C: N x M matrix
                   )
{
  register int i, j, k;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      C[i * M + j] = 0;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      for (k = 0; k < P; k++)
        C[i * M + j] += A[i * P + k] * B[k * M + j];

  return;
}

// code generated by Maple ( assummed correct )
void
inv3 (double * A, //! input: 3 x 3 Matrix
      double * Ap //! output: 3 x 3 Matrix inverse
     )
{
  double t4 = A[2] * A[3];
  double t6 = A[2] * A[6];
  double t8 = A[1] * A[3];
  double t10 = A[1] * A[6];
  double t12 = A[0] * A[4];
  double t14 = A[0] * A[7];
  double t17 =
    1.0 / (t4 * A[7] - t6 * A[4] - t8 * A[8] + t10 * A[5] + t12 * A[8] -
           t14 * A[5]);
  Ap[0] = (A[4] * A[8] - A[7] * A[5]) * t17;
  Ap[1] = -(A[3] * A[8] - A[6] * A[5]) * t17;
  Ap[2] = (A[3] * A[7] - A[6] * A[4]) * t17;
  Ap[3] = -(-A[2] * A[7] + A[1] * A[8]) * t17;
  Ap[4] = (-t6 + A[0] * A[8]) * t17;
  Ap[5] = -(-t10 + t14) * t17;
  Ap[6] = (-A[2] * A[4] + A[1] * A[5]) * t17;
  Ap[7] = -(-t4 + A[0] * A[5]) * t17;
  Ap[8] = (-t8 + t12) * t17;

}

// solve linear equations for LHS: 3 x 3 and RHS: 3 x N
void
linsolve (double * A,  //! input: LHS 3 x 3 Matrix
          double * d,  //! input: RHS 3 x N Matix
          int N,       //! number of columns in RHS
          double * x   //! output  3 x N Matrix
         )
{
  double Ainv[9];
  inv3 (A, Ainv);
  matrix_matrix_mult (Ainv, 3, 3, d, N, x);
  return;
}


// linear solver for 3 x 3 matrix and 3 x 1 vector
void
linsolve31 (double * A, double * d, double * sol )
{
  double t1 = A[4] * A[8] - A[5] * A[7];
  double t2 = A[1] * A[5] - A[2] * A[4];
  double t3 = -A[1] * A[8] + A[2] * A[7];
  double t4 = t3 * A[3] + t2 * A[6] + A[0] * t1;
  t4 = 1. / t4;
  sol[0] = (t1 * d[0] + t3 * d[1] + t2 * d[2]) * t4;
  sol[1] = (-(A[3] * A[8] - A[5] * A[6]) * d[0] +
           (A[0] * A[8] - A[6] * A[2]) * d[1] -
           (A[0] * A[5] - A[2] * A[3]) * d[2]) * t4;
  sol[2] = ((A[3] * A[7] - A[4] * A[6]) * d[0] -
            (A[0] * A[7] - A[6] * A[1]) * d[1] +
            (A[0] * A[4] - A[3] * A[1]) * d[2]) * t4;
  return;
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


//general file for artificial viscosity
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

//overloading file for artificial viscosity
double art_vis ( double rhoab, double sndspdab, double rab[3], double vab[3], double rsqab, double h)
{
	int    k;
	double miuab;
	double vrab = 0.;
	double vis = 0.;

	for (k = 0; k < DIMENSION; k++)
	     vrab += (rab[k] * vab[k]);

	miuab = h * vrab / (rsqab + ata_P * h * h);

    if (vrab > 0)
		vis = 0.;
    else if (vrab < 0)
		vis = (beta_P * miuab - alf_P*sndspdab) * miuab / rhoab;

	return vis;
}

//function that used to determine the property of air: density, pressure, (temperature not explicitly output) and internal energy
void air_prop (double *coord, double * energy, double *pressure, double * density)
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
void air_prop (double *coord, double *range, double * energy, double *pressure, double * density, double * mass)
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

	//integration is based on interpolation of order 2 (three points)
	double h1=range[4];
	double h2=range[5];
	double T1, T2;
	T1 = Ta0_P *(h1 < 0) + (Ta0_P-miu1_P*h1)*((h1>=0)&&(h1<H1_P))+(Ta0_P-miu1_P*H1_P)*((h1>=H1_P)&&(h1<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h1-H2_P))*((h1>=H2_P)&&(h1<H3_P));
	T2 = Ta0_P *(h2 < 0) + (Ta0_P-miu1_P*h2)*((h2>=0)&&(h2<H1_P))+(Ta0_P-miu1_P*H1_P)*((h2>=H1_P)&&(h2<H2_P))+(Ta0_P-miu1_P*H1_P+miu2_P*(h2-H2_P))*((h2>=H2_P)&&(h2<H3_P));
	double p1 = pa0_P*exp(C0*h1/T1)*(h1>=0) +  pa0_P*(h1<0);
	double d1 = p1/(Ra_P*T1);
	double p2 = pa0_P*exp(C0*h2/T2)*(h2>=0) +  pa0_P*(h2<0);
	double d2 = p2/(Ra_P*T2);
	double d = *density;

	//The following code is based in numerical integration, coefficient are 1/6, 4/6 1/6
	*mass = (range[1]-range[0])*(range[3]-range[2])*(range[5]-range[4])*(0.1666667*d1+0.6666666*d+0.1666667*d2);

#ifdef DEBUG
	bool print = false;
	if (print)
        cout << "h=" << h<< "\t T=" << T << "\t p="  << *pressure << "\t rho=" << *density << "\t e=" << *energy << "\n" << endl;
#endif
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

double parabolic_vel(double R, double rsq, double umax)
{
	return umax*(1-(rsq/(R*R)));
}

//function that used to determine the value of face by the face's index;
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

//function that determines parameters of certain particle
void initial_air (Particle * pi)
{
	  int i;

	  double xi[DIMENSION];
	  double vel[DIMENSION] = {0., 0., 0.};  //initial velocity need to be set to zero
	  double prss, erg, dens, mss;
	  double range[6];
	  double sml2;

	  // go through particle table
//	  if (!pi->is_guest())//count number of non-guest particles.
//	     {

	  	  for (i = 0; i < DIMENSION; i++)
	  		  xi[i] = *(pi->get_coords() + i);

	  	  sml2 = 0.5*(pi->get_smlen ());
	  	  range[0]=xi[0]-sml2;
	  	  range[1]=xi[0]+sml2;
	  	  range[2]=xi[1]-sml2;
	  	  range[3]=xi[1]+sml2;
	  	  range[4]=xi[2]-sml2;
	  	  range[5]=xi[2]+sml2;

	  	  air_prop (xi, range, &erg, &prss, &dens, &mss);

	      //put data back into particle:
	      pi->put_density(dens);
	      pi->put_energy(erg);
	      pi->put_pressure(prss);
	      pi->put_velocity(vel);
	      pi->put_mass(mss);

	      //the second variable need to be updated.
	      pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P);

//	     }//end of go through all particles

	return;
}

#ifdef DEBUG
      //function to check where does the negative sound speed comes from
      bool check_particles_sndspd (THashTable * P_table)
      {
    	  bool ng_sndspd = false;

    	  THTIterator *itr = new THTIterator(P_table);
    	  Particle *p_curr = NULL;

    	  while ((p_curr = (Particle *) itr->next()))
    	      if (p_curr->need_neigh())
    	      {
    	        // calc speed of sound through the medium
    	        double c = p_curr->get_sound_speed();
    	        if (c < 0)
    	        {
    	        	ng_sndspd =true;
    	        	cout << "c = " << c << endl;
//    	        	return ng_sndspd;
    	        }
    	      }

    	  return ng_sndspd;
      }

#endif


