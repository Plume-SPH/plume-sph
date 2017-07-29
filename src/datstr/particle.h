/*
 * particle.h
 *
 *  Created on: Feb 16, 2015
 *      Author: zhixuanc
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#ifdef DEBUG
#   include <cstdio>
#   include <cstdlib>
#   include <debug_header.h>
#endif

#include <vector>

using namespace std;

#  include <constant.h>
#  include <thashtab.h>
#  include <pack_data.h>

class Particle
{
  // friend functions in repartition
  friend void pack_particles (Particle *, ParticlePack *);
  friend void unpack_particle (ParticlePack *, Particle *);
#ifdef DEBUG
  friend void particle_deb (int myid);
#endif

protected:

  //! Data Structure properties
  TKey key;

  //! mass of the each particle
  double mass;

  //! support of each particle
  double smlen;

  //! original smoothing length
  /*
   * When dealing with multiple phase, each phase might have different smoothing length for computational efficiency purpose
   * However, to conserve momentum and energy, the same smoothing length should be used ---> So smlen of phase2 will be changed to equal to that of phase1 as soon as it erupted out.
   * Then in density updating, too large smoothing length will smear out the density and will lead to non-physics results --> for example, density of the plume near the vent should be almost equals to that in the vent, but when very large smoothing length is used, density near the vent will be much smaller than that in the vent
   *
   * The strategy to get rid of this issue is:
   * 1) updating of momentum and energy use larger smoothing length. ---> conserve energy and momentum
   * 2) updating of density based on summation expression and use its original smoothing length --> The normalized summation expression will always guarantee conservation of mass
   */
#if DENSITY_UPDATE_SML==0 || ADAPTIVE_SML==2
  double smlen_original;  //Original sml will not change during whole simulation ---> if sml is adaptive, we will not need sml_original
#endif

  //! pressure
  double pressure;

  //! mass fraction
  double mass_frac;

  //! new mass fraction
  double new_mass_frac;

  //! gamma
  double gamma;

  //! sound speed
  double sound_speed;

  //Temperature
  double temperature;

  double specific_heat_p; //specific heat under constant pressure

#if ADAPTIVE_SML==2
  double rho_based_on_dx; //This is an intermediate density calculated based on original smoothing length. ---> Needed for ADKE
#endif

#if ADAPTIVE_SML==31
  //! An indicate used to show whether the particle mass is relatively larger than surrounding or smaller on average.
  /*
   * m_ind > 0 particle mass is smaller than surrounding on average
   * m_ind <0  particle mass is larger than surrounding on average
   * m_ind = 0, no need to consider particle mass non-equilibrium
   */
  double m_ind;
#endif

  //density of each phase
  double phase_density[PHASE_NUM]; //phase1 = air, phase2 =erupted material
                               //According to the simplest Japanese Model, d1=d*(1-ks), d2=d*ks
                               //For more complicated model the above equation does not hold
                               //What also need to be mentioned here is that for a real multiple phase model, the density of each phase is supposed to be primitive variable.

  //density of each phase new variables ---> will need if I want to update density based on mass conservation instead of SPH summation
  double new_phase_density[PHASE_NUM];

  //Only particle need_neigh will have velocity smoothed -->guest particle, smoothed velocity will be updated in syn process
  double smoothed_v[DIMENSION]; //smoothed velocity

  //! current (normalized) postion of each particle
  double coord[DIMENSION];

#if ADAPTIVE_SML==3 || ADAPTIVE_SML==31
  //! gradient of mass
  double dm[DIMENSION];
#endif

  //! state_vars are normalized density, velocities;  stresses is not included in my model
  double state_vars[NO_OF_EQNS];

  //! new_state_vars stores updated state variables
  double new_state_vars[NO_OF_EQNS];

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //derivatives
  double d_rho[DIMENSION];
  double d_u[DIMENSION];
  double d_v[DIMENSION];
  double d_w[DIMENSION];
  double d_p[DIMENSION];
#endif

  double smoothed_e; //smoothed energy

  //!to indicate whether the particle is eruption ghost particle or not
  int bc_type; //bc_type=0: eruption ghost
               //bc_type=1: pressure ghost
               //bc_type=2: wall ghost
               //bc_type=99; No_save_real: treated like real particle, the only different is do not save it when output: This is useful for 1D shock tube problem test, we can use such particles for boundary conditions
               //bc_type=100: corresponding to real partilce in the old data
               //First of all, the particle has to be ghost particles.
               //The tricky here is that as long as bc_type is not zero,
               //it will not move up while imposing eruption boundary condition.

  //! new_old for wall guest updates across process boundaries
  int new_old;     //new_old = 1 : new
                   //new_old = -1 : old
                   //new_old = 0 :default

  /*
   * involved flag: based on velocity and mass fraction
   * only for bctp = 100 or bctp = 2
   * involve = 0 if not involved
   * involve = 1 if potential involved
   * involve = 2 if involved
   */
  int involved;

  //! indicate which phase does the particles belong to
  int phase_num; //1 for air, 2 for erupt material

  //indicate which processor does the particles belong to
  int myprocess; //For guest particles, myprocess should be the "home" process id, not the current process on which the guest particle locate on

  //! if the particle is real, guest or a ghost
  bool guest;

  //! update delayed
  bool update_delayed; //if my image goes to other process, the delay will be delayed!

  //! aleardy searched reflection, if flag is up
  bool reflection; //is true if the ghost particle has already been reflected-->its imaged has already be found

  //!  neighbors
  vector < TKey > neighs;

public:
  //! Constructors
  //! Base constructor
    Particle ();

    //! Constructor for initial air
     Particle (
                 //! hash-table key
                 unsigned *,
                 //! particle coodinates
                 double *,
                 //! particle mass
                 double,
                 //! smoothing length
                 double,
                 //! pressure
                 double,
                 //! mass_fraction
                 double,
                 //! gamma
                 double,
                 //! sound speed
                 double,
                 //! phase id
                 int,
                 //! process id
                 int,
                 //!bc_type
                 int
      );
  //Constructor for adding eruption initial condition
   Particle (
		      //! hash-table key
		      unsigned*,
		      //! particle coodinates
		      double *,
		      //! particle mass
		      double ,
		      //! smoothing length
		      double,
		      //! myid
		      int,
		      // Verticle velocity
		      double,
		      // internal energy
		      double,
		      // density
		      double,
		      // pressure
		      double,
		      // gamma
		      double,
			  //sound speed
			  double,
			  //Initial mass fraction of solid in erupted material ng0_P
			  double,
			  //Cvs_P
			  double,
			  //Cvg_P
			  double,
			  //Cva_P
			  double,
			  //Rg_P
			  double,
			  //Ra_P
			  double
		      );

   //overloading
   //Constructor for adding influx particles
   Particle (
				    unsigned *, //keyin,
				    double *,   //crd,
				    double ,	//m,
				    double ,	//h,
				    int , 		//id,
				    double,		// msfc0,
					double ,	//Vx,
					double , 	//Vy,
					double ,	//ev0,
					double ,	//rhov,
					double ,	//pv0,
					double ,	//gamma_v,
					double ,	//sndspd,
					double ,	//ng0_P,
					double ,	//Cvs_P,
					double ,	//Cvg_P,
					double ,	//Cva_P,
					double ,	//Rg_P,
					double   	//Ra_P
		            );

#if CODE_DIMENSION==1
// constructor for initial condition of Shocktube problem
  Particle (
		  unsigned*,  // *keyin,
		  double *,  //crd,
		  double,    //m,
		  double ,   //h ,
		  double ,   //des,
		  double ,   //vel,
		  double ,   //prss,
		  double ,   //gmm,
		  double ,   //sndspd,
		  int ,      //bc,
		  int        //invloved
		  );
#endif

  //overload update secondary variables
  void update_second_var(
		  //Initial mass fraction of solid in erupted material ng0_P
		  double,
		  //Cvs_P
		  double,
		  //Cvg_P
		  double,
		  //Cva_P
		  double,
		  //Rg_P
		  double,
		  //Ra_P
		  double,
		  //rho_0
		  double
		                 );

  //function that calculate the density of mixture base on equation in suzuki's 2005 paper
  void calc_density_suzuki(
		  double // msf
		  );

  //! get hash-table key
  TKey getKey () const
  {
    return key;
  }

  //! get particle mass
  double get_mass () const
  {
    return mass;
  }

  //! get smoothing length of current paricle
  double get_smlen () const
  {
    return smlen;
  }

#if DENSITY_UPDATE_SML==0 || ADAPTIVE_SML==2
  //! get original smoothing length of current paricle
  double get_original_smlen () const
  {
    return smlen_original;
  }

  //! update original smoothing length
  void put_smlen_original (double h)
  {
    smlen_original = h;
  }
#endif


  //!get coordinates
  const double *get_coords () const
  {
    return coord;
  }

#if ADAPTIVE_SML==3 || ADAPTIVE_SML==31
  //!get mass gradient
  const double *get_mass_grad () const
  {
    return dm;
  }
#endif

#if ADAPTIVE_SML==31
  //!get mass indicator
  const double get_mass_indicator () const
  {
    return m_ind;
  }

  //return a inter value 1, -1 or 0 for adaptive sml
  int which_mass_ind ()
  {
    if (m_ind < -MASS_IND_THRESH)
    	return -1 ; //Should reduce equivalent particle mass in adaptive sml--> Means this particle is larger than surrounding particles on average
    else if (m_ind > MASS_IND_THRESH)
    	return 1;  ///Should increase equivalent particle mass in adaptive sml--> Means this particle is smaller than surrounding particles on average
    else
    	return 0;
  }
#endif

  //! get velocity
  const double *get_vel () const
  {
    return (state_vars + 1);
  }

  //!get state variable
  const double *get_state_vars () const
  {
    return state_vars;
  }

  //!
  const double *get_new_state_vars () const
  {
    return new_state_vars;
  }

  //get phase density
  const double *get_phase_density () const
  {
    return phase_density;
  }

  //get new phase density
  const double *get_new_phase_density () const
  {
    return new_phase_density;
  }

  //get new_mass_frac
  const double get_new_mass_frac () const
  {
    return new_mass_frac;
  }

  //! get density
  const double get_density () const
  {
    return state_vars[0];/*In C++, array index starts from zero*/;
  }

#if ADAPTIVE_SML==2
  //! get dx_density
  const double get_dx_density () const
  {
    return rho_based_on_dx;
  }
#endif

  //! get new density
  const double get_new_density () const
  {
    return new_state_vars[0];/*In C++, array index starts from zero*/
  }
  //! get energy
  const double get_energy () const
  {
    return state_vars[NO_OF_EQNS-1];/*In C++, array index starts from zero*/;
  }

  //! get new energy
  const double get_new_energy () const
  {
    return new_state_vars[NO_OF_EQNS-1];/*In C++, array index starts from zero*/
  }
  //! get pressure
  const double get_pressure () const
  {
    return pressure;
  }

  //! get mass_frac
  const double get_mass_frac () const
  {
    return mass_frac;
  }

  //! get gamma
  const double get_gamma () const
  {
    return gamma;
  }

  //! get sound_speed
  const double get_sound_speed () const
  {
    return sound_speed;
  }

  //! get temperature
  const double get_temperature () const
  {
    return temperature;
  }

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //! get density derivative
  const double * get_density_derivative () const
  {
    return d_rho;
  }

  //! get velocity_x derivative
  const double * get_velocity_u_derivative () const
  {
    return d_u;
  }

  //! get velocity_y derivative
  const double * get_velocity_v_derivative () const
  {
    return d_v;
  }

  //! get velocity_z derivative
  const double * get_velocity_w_derivative () const
  {
    return d_w;
  }

  //! get velocity_x derivative
  const double * get_pressure_derivative () const
  {
    return d_p;
  }
#endif


  //! get my processor
  const int get_my_processor () const
  {
    return myprocess;
  }

  //! get bc_type
  const int get_bc_type () const
  {
    return bc_type;
  }

  //! get involved
  const int get_involved () const
  {
    return involved;
  }

  //! get table of keys of all neighbors
  vector < TKey > get_neighs () const
  {
    return neighs;
  }

  //! get new_old info
  int get_new_old ()
  {
#ifdef DEBUG
	  if (!guest)
		  cout << "You are trying to get new_old for non-guest particles! \n" <<endl;
#endif

	  return new_old;
  }

  //! check if particle is need neighbor information: only these real and non-guest particles need neighbor
  bool need_neigh () const
  {
    return ((bc_type == 100  || bc_type == 99) && (!guest));
  }

  //! check if particle will contribute to the image of wall ghost which reflected into the domain
  // The contribution of pressure ghost particles is also considered.
  bool contr_image() const
  {
	  return ((bc_type == 100) || (bc_type == 99) || (bc_type == 1));//real and pressure ghost will contribute to image
  }

  //! check if particle will contribute to the density estimation?
  // only "real" particles considered while updating density.
  //contr_dens is equivalent to (is_real)
  bool contr_dens() const
  {
#if PGHOST_CONTRIBUTE_DES==1
	  return (bc_type == 100 || bc_type == 99 || bc_type == 1); //real+pressure
#elif PGHOST_CONTRIBUTE_DES==2
	  return (bc_type == 100 || bc_type == 99 || bc_type == 1 || bc_type == 0); ////real+pressure+erupt
#else
	  return (bc_type == 100 || bc_type == 99);
#endif //PGHOST_CONTRIBUTE_DES

  }

  //! check if particle will contribute to the velocity smoothing
  // only "real" particles should be considered for sure
  // In my opinion, wall ghost should also contribute to the velocity smoothing ---> to get real no-slip bc ---> But this will also cause rotation (not sure whether it is physical or no physical)
  // What about pressure ghost? ---> not sure --> no, because that will cause additional turbulence near the boundary
  // erupted ghost should not be considered for sure
  bool contr_vel_smooth() const
  {
	  return (bc_type == 100 || bc_type == 99 || bc_type == 2);
  }

  //! check if particle is erupt ghost particle or not
  bool is_erupt_ghost () const
  {
    return (bc_type == 0);
  }

  //! check if particle is pressure ghost particle or not
  bool is_press_ghost () const
  {
    return (bc_type == 1);
  }

  //! check if particle is wall ghost particle or not
  bool is_wall_ghost () const
  {
    return (bc_type == 2);
  }

  //! check if particle is real particle or not
  bool is_real () const
  {
    return (bc_type == 100 || bc_type == 99);
  }

  //! check if particle is guest particle or not
  bool is_guest () const
  {
    return guest;
  }

  //! return 0 if particle is guest, if not return 1
  double get_guest() const
  {
    return (double) (guest ?  0 : 1 );
  }

  //! check whether wall-ghost particles is not_updated
  bool is_not_updated () const
  {
#ifdef DEBUG
	  if (bc_type!=2)
		  cout << "You are trying to get update_delayed for particle which is not a wall-bc-ghost! \n" <<endl;
#endif
      return update_delayed;
  }

  //! check if wall-ghost particles already searched the ghost reflection
  bool has_reflection ()
  {
#ifdef DEBUG
	  if (bc_type!=2)
		  cout << "You are trying to get reflection for particle which is not a wall-bc-ghost! \n" <<endl;
#endif
      return reflection;
  }

  //check whether certain particle is involved or not
   bool is_involved() const
   {
       return involved == 2;
   }

   //check whether certain particle is potential involved or not
   bool is_potential_involved() const
   {
       return involved == 1;
   }

  //!check phase number
  int which_phase()
  {
	return phase_num;
  }

  //! update particles mass
  void put_mass (double pmss)
  {
    mass = pmss;
  }

  //! update particles velocity
  void put_vel (double vel[])
  {
    for (int i=0; i<DIMENSION; i++)
    	state_vars[i+1]=vel[i];

  }

  //! update density value
  void put_density (double den)
  {
    state_vars[0] = den;
  }

#if ADAPTIVE_SML==2
  //! update dx_density value
  void put_dx_density (double den)
  {
	  rho_based_on_dx = den;
  }
#endif

  //! update new_density value
  void put_new_density (double den)
  {
    new_state_vars[0] = den;
  }

  //! update smoothing length
  void put_smlen (double h)
  {
    smlen = RATIO_SML_DX*h;
  }

  //! update neighbor information
  void put_neighs (vector < TKey > n)
  {
    neighs = n;
  }

  //! get delayed update flag
  void put_update_delayed (bool val)
  {
#ifdef DEBUG
	  if (bc_type!=2)
		  cout << "You are trying to put update_delayed for particle which is not a wall-bc-ghost! \n" <<endl;
#endif

    update_delayed = val;
  }

  //! update new state variables
  void put_new_state_vars (double u[])
  {
	  int i;
	  for (i = 0; i < NO_OF_EQNS; i++)
           new_state_vars[i] = u[i];
  }

  //! update new phase density
  void put_new_phase_density (double phsd[])
  {
	  int i;
	  for (i = 0; i < PHASE_NUM; i++)
		  new_phase_density[i] = phsd[i];
  }

  //! put new_mass_frac
  void put_new_mass_frac (double mssfrc)
  {
	  new_mass_frac=mssfrc;
  }

  //! put state_vars. This should only be uesed for ghost particles
  void put_state_vars (double u[])
  {
    for (int i = 0; i < NO_OF_EQNS; i++)
      state_vars[i] = u[i];
  }

  //! update phase density
  void put_phase_density (double phsd[])
  {
	  int i;
	  for (i = 0; i < PHASE_NUM; i++)
		  phase_density[i] = phsd[i];
  }

  //put energy
  void put_energy (double engr)
  {
	  state_vars[ NO_OF_EQNS-1]=engr;
  }

  //put pressure
  void put_pressure (double prss)
  {
	  pressure = prss;
  }

  //put velocity
  void put_velocity (double v[])
  {
	  int i;
	  for ( i = 1; i < NO_OF_EQNS -1 ; i++)
	      new_state_vars[i] = v[i];
  }
  void put_new_energy (double new_engr)
  {
      new_state_vars[ NO_OF_EQNS-1]=new_engr;
  }

#if (USE_GSPH==1 || USE_GSPH==2)  //Assume 3D
  //! put density derivative
  void put_density_d (double d[])
  {
      for (int i=0; i<DIMENSION; i++)
    	  d_rho[i]=d[i];
  }

  //! put velocity u derivative
  void put_velocity_u_d (double d[])
  {
      for (int i=0; i<DIMENSION; i++)
    	  d_u[i]=d[i];
  }

  //! put velocity v derivative
  void put_velocity_v_d (double d[])
  {
      for (int i=0; i<DIMENSION; i++)
    	  d_v[i]=d[i];
  }

  //! put velocity w derivative
  void put_velocity_w_d (double d[])
  {
      for (int i=0; i<DIMENSION; i++)
    	  d_w[i]=d[i];
  }

  //! put pressure derivative
  void put_pressure_d (double d[])
  {
      for (int i=0; i<DIMENSION; i++)
    	  d_p[i]=d[i];
  }
#endif

  //! update state_vars
  void update_state_vars ()
  {
    for (int i = 0; i < NO_OF_EQNS; i++)
      state_vars[i] = new_state_vars[i];
  }

  //! update state_vars
  void update_phase_density ()
  {
    for (int i = 0; i < PHASE_NUM; i++)
    	phase_density[i] = new_phase_density[i];
  }

  //! update mass_frac
  void update_mass_frac ()
  {
	  mass_frac=new_mass_frac;
  }

  //! update density
  void update_density ()
  {
	  state_vars[0]=new_state_vars[0];
  }

  //! update density
  void update_energy ()
  {
	  state_vars[NO_OF_EQNS-1]=new_state_vars[NO_OF_EQNS-1];
  }

  //! update mass_frac
  void put_mass_frac (double mssfrc)
  {
      mass_frac=mssfrc;
  }

  //! update gamma
  void put_gamma (double gmm)
  {
      gamma=gmm;
  }

  //! put sound_speed
  void put_sound_speed (double sndspd)
  {
      sound_speed=sndspd;
  }

  //! put temperature
  void put_temperature (double temp)
  {
	  temperature=temp;
  }

  //! put my_processor
  void put_my_processor (int id)
  {
	  myprocess = id;
  }

  // update partilce positions
  void put_coords (double *x)
  {
    for (int i = 0; i < DIMENSION; i++)
      coord[i] = *(x + i);
  }

#if ADAPTIVE_SML==3 || ADAPTIVE_SML==31
  // update mass gradient
  void put_mass_grad (double *m_grad)
  {
    for (int i = 0; i < DIMENSION; i++)
      dm[i] = *(m_grad + i);
  }
#endif

#if ADAPTIVE_SML==31
  // update mass indicator
  void put_mass_indicator (double ind)
  {
	  m_ind = ind;
  }
#endif

  //! update guest info
  void put_guest_flag (bool val)
  {
    guest = val;
  }

  //! update new_old info
  void put_new_old (int info)
  {
#ifdef DEBUG
	  if (!guest)
		  cout << "You are trying to put new_old for non-guest particles! \n" <<endl;
#endif

      new_old = info;
  }

  //put bc type
  void put_bc_type (int in)
  {
	  bc_type = in;
  }

  //! change erupt to real particle as the erupted particle come out from duct
  void erupt_turn_real ()
  {
      bc_type = 100;
  }

  //! set reflection flag
  void set_reflection (bool val)
  {
#ifdef DEBUG
	  if (bc_type!=2)
		  cout << "You are trying to get reflection for particle which is not a wall-bc-ghost! \n" <<endl;
#endif
      reflection = val;
  }

  void set_involved_flag( int inv)
  {
      involved = inv;
  }

  //! set phase number
  void set_phs_num (int phs_nmb)
  {
	 phase_num = phs_nmb;
  }

  // Operator overloads
  bool operator== (const Particle & rhs) const;

    //! get specific heat under constant pressure
    const double get_specific_heat_p() const
    {
      return specific_heat_p;
    }

    //! get smoothed velocity
    const double *get_smoothed_velocity() const
    {
      return smoothed_v;
    }

    //! get smoothed energy
    const double get_smoothed_energy() const
    {
      return smoothed_e;
    }

    //put smoothed velocity
    void put_smoothed_velocity (double smed_u[])
    {
  	  for(int i=0; i<DIMENSION; i++)
  		 smoothed_v[i]=smed_u[i];
    }

    //put smoothed energy
    void put_smoothed_energy (double smed_e)
    {
    	smoothed_e = smed_e;
    }

    //put smoothed specific heat
    void put_specific_heat_p (double Cp)
    {
  	  specific_heat_p = Cp;
    }
};

#endif /* PARTICLE_H */
