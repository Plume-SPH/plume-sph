/*
 * debug_lib.cc
 *
 *  Created on: Jun 16, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

#include <mpi.h>
#include <hashtab.h>
#include <thashtab.h>
#include <hilbert.h>
#include <bucket.h>
#include <particler.h>
#include <multiproc.h>
#include <properties.h>
#include <outforms.h>

#include "sph_header.h"
#include "constant.h"
#include "parameters.h"

#include <hdf5.h>
#include "hdf5calls.h"

/*find particle by the key
 * keyin is input particle key
 * keycheck is the key that is given and all keyin will be compared with keycheck
  */
bool find_particle (unsigned* keyin, unsigned* keycheck)
{
	int i;
	for (i=0; i<TKEYLENGTH; i++)
		if (keyin[i] != keycheck[i])
			return false;

    return true;
}

/*find particle by the pos
 *
 */
bool find_particle (double* in, double* check)
{
	int i;
	for (i=0; i<DIMENSION; i++)
		if (in[i] != check[i])
			return false;

    return true;
}

/*find particle by the range of particle position
 *
 */
bool find_particle_pos_range  (double *in, double* check)
{
	int i;
	for (i=0; i<DIMENSION; i++)
		if (!((in[i] >= check[2*i]) && (in[i]<= check[2*i+1])))
			return false;

    return true;
}

//function to find particle with given key, will be useful in debugging.
void check_particle_bykey (THashTable * P_table)
{

    bool do_search = true;
    bool find;
    unsigned keycheck[TKEYLENGTH] = {71901794, 4252428424, 0}; //key of its neighbor which is missing
    //{93400792, 2478187757, 0}; //key of the particle itself is:
    //
    unsigned keytemp[TKEYLENGTH] ;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	for (i = 0; i < TKEYLENGTH; i++)
		  		keytemp[i] = p_curr->getKey ().key[i];

		  	if (find_particle (keytemp, keycheck))
		  		cout << "The particle found!" << endl;
		 }
	 }

}

//overload function to find particle with given key, will be useful in debugging. --> output my process
bool check_particle_bykey (THashTable * P_table, int* id)
{

    bool do_search = true;
    bool find;
    unsigned keycheck[TKEYLENGTH] = {71606976, 8257279, 0}; //key of its neighbor which is missing
    //{93400792, 2478187757, 0}; //key of the particle itself is:
    //
    unsigned keytemp[TKEYLENGTH] ;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	for (i = 0; i < TKEYLENGTH; i++)
		  		keytemp[i] = p_curr->getKey ().key[i];

		  	if (find_particle (keytemp, keycheck))
		  	{
		  		cout << "The particle found!" << endl;
//		  		*id = p_curr->get_my_processor (); //when I want to track myprocess
		  		*id = p_curr->get_new_old ();
		  		return true;
		  	}

		 }
	 }
	*id = 100000;
	return false;
}
//go through all particle and check their type!
void check_particle_all_type (THashTable * P_table)
{
	bool do_search = true;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	if (p_curr->is_real())
		  		cout << "Type: real; \t pos:" << *p_curr->get_coords () <<*(p_curr->get_coords ()+1) << *(p_curr->get_coords ()+2) << endl;
		  	else if (p_curr->is_wall_ghost())
                cout << "Type: wall ghost; \t pos:" << *p_curr->get_coords () <<*(p_curr->get_coords ()+1) << *(p_curr->get_coords ()+2) << endl;
		  	else if (p_curr->is_press_ghost())
		  	    cout << "Type: pressure ghost; \t pos:" << *p_curr->get_coords () <<*(p_curr->get_coords ()+1) << *(p_curr->get_coords ()+2) << endl;
		  	else if (p_curr->is_erupt_ghost())
		  		cout << "Type: eruption ghost; \t pos:" << *p_curr->get_coords () <<*(p_curr->get_coords ()+1) << *(p_curr->get_coords ()+2) << endl;
		  	else
		  		cout << "Qu ni mei de, che jiba dan!" <<endl;
		 }
	 }

}

//check particles in certain region
#if CODE_DIMENSION==3
void check_particle_bypos (THashTable * P_table)
{

    bool do_search = true;
    double range_x[2]={-300,700};
    double range_y[2]={-300,700};
    double range_z[2]={0, 700};
    double pcrd[DIMENSION];
    int bctp;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	for (i = 0; i < DIMENSION; i++)
		  		pcrd[i] = *(p_curr->get_coords ()+i);

		  	if ((pcrd[0]>range_x[0] && pcrd[0]<range_x[1]) &&
		  		(pcrd[1]>range_y[0] && pcrd[1]<range_y[1]) &&
		  		(pcrd[2]>range_z[0] && pcrd[2]<range_z[1]))
		  	{
		  		bctp = p_curr->get_bc_type();
		  		cout << "The particle found!" << endl;
		  		cout << "done!" << endl;
		  	}
		 }
	 }

}
#elif CODE_DIMENSION==2
void check_particle_bypos (THashTable * P_table)
{

    bool do_search = true;
    double range_x[2]={-300,700};
    double range_y[2]={-300,700};
    double pcrd[DIMENSION];
    int bctp;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	for (i = 0; i < DIMENSION; i++)
		  		pcrd[i] = *(p_curr->get_coords ()+i);

		  	if ((pcrd[0]>range_x[0] && pcrd[0]<range_x[1]) &&
		  		(pcrd[1]>range_y[0] && pcrd[1]<range_y[1]) )
		  	{
		  		bctp = p_curr->get_bc_type();
		  		cout << "The particle found!" << endl;
		  		cout << "done!" << endl;
		  	}
		 }
	 }

}
#elif CODE_DIMENSION==1
void check_particle_bypos (THashTable * P_table)
{

    bool do_search = true;
    double range_x[2]={20,1e123};
    double pcrd[DIMENSION];
    int bctp;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	for (i = 0; i < DIMENSION; i++)
		  		pcrd[i] = *(p_curr->get_coords ()+i);

		  	if ((pcrd[0]>range_x[0] && pcrd[0]<range_x[1]))
		  	{
		  		bctp = p_curr->get_bc_type();
		  		cout << "The particle found!" << endl;
		  		cout << "done!" << endl;
		  	}
		 }
	 }

}
#endif

//find the particle with non-physical density
void find_large_density_particle (THashTable * P_table)
{

    bool do_search = true;
    double threshold = 20.;
    double density;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	density=p_curr->get_density();

		  	if (density >=threshold)
		  	{
		  		cout << "The particle found!" << endl;

		  	}
		 }
	 }//end of go through all particles

}

/*find all particle of given pahse
 * This will be used to find particle of phase 2 ----> get its key and then track the movement of that particle
 */
void find_particle_by_phase (THashTable * P_table)
{

    bool do_search = true;
    int phase_need = 2;
    int phase;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;

	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search && p_curr->is_real()) //only look for real particles
		{
		  	phase=p_curr->which_phase();

		  	if (phase ==phase_need)
		  	{
		  		cout << "The particle found!" << endl;
		  	}
		 }
	 }//end of go through all particles
}

/*
 * Find the highest position of all phase2 particles
 */
void find_highest_z (THashTable * P_table)
{

    bool do_search = true;
    int phase_need = 2;
    double zmax=0., z;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;
	Particle *p_max = NULL;

	if (do_search )
	   while ((p_curr = (Particle *) itr->next()))
		   if ((p_curr->is_real()) && (p_curr->which_phase()==phase_need))
	       {
		      z = *(p_curr->get_coords()+2);

		  	  if (z > zmax)
		  	  {
		  		 zmax = z;
		  		 p_max = p_curr;
		  	  }
	        }//end of go through all particles

	cout << "maximum height of phase2 particle is: " << zmax << endl;
	cout << "Up speed of that particle is: " <<  *(p_max->get_vel()+2) <<endl;
}

//output certain type of particle
void
write_h5part_bctp(int myid, int numproc, THashTable * P_table, TimeProps * timepros, int bctp, char *pre)
{
  int i, j;
  vector < double >x, y, z, Vx, Vy, Vz, rho, engr, mssfrc;
  char filename[18];
  static int step = 0;
  hid_t fp;
  herr_t ierr;

#ifdef PARALLEL_IO
  sprintf(filename, "%01c__bctp%03d_pvplot_out.h5part", *pre, prefix);
#else
  sprintf(filename, "%01c__bctp%03d_pvplot%03d.h5part", *pre, bctp, myid);
#endif

  if (timepros->ifstart())
    fp = GH5_fopen(filename, 'w');
  else
    fp = GH5_fopen(filename, 'a');

  THTIterator *itr = new THTIterator(P_table);
  Particle *pi = NULL;
  int my_count = 0;

  while ((pi = (Particle *) itr->next()))
  {
#ifndef WRITE_GHOSTS
//    if (!(pi->is_guest()) && (pi->get_bc_type () == bctp))//non guest and is real
//    {
#endif
      rho.push_back(pi->get_density());
      x.push_back(*(pi->get_coords()));
      Vx.push_back(*(pi->get_vel()));
      y.push_back(*(pi->get_coords() + 1));
      Vy.push_back(*(pi->get_vel() + 1));
      my_count++;
      z.push_back(*(pi->get_coords() + 2));
      Vz.push_back(*(pi->get_vel() + 2));
      engr.push_back(pi->get_energy ());
      mssfrc.push_back(pi->get_mass_frac());
#ifndef WRITE_GHOSTS
//    }
#endif
  }
  char group[10];

  sprintf(group, "Step#%d", step);
  step++;

  hid_t gid = GH5_gopen(fp, group, 'w');

  int start;
  int size = 0;
  int *id_lims = new int[numproc + 1];

  for (i = 0; i < numproc + 1; i++)
    id_lims[i] = 0;

#ifdef PARALLEL_IO
  // make space to recv num of particles on all proc
  int *size_arr = new int[numproc];

  for (i = 0; i < numproc; i++)
    size_arr[i] = 0;

  // get sizes from all around the world
  MPI_Allgather(&my_count, 1, MPI_INT, size_arr, 1, MPI_INT, MPI_COMM_WORLD);

  // create array index partitons
  id_lims[0] = 0;
  for (i = 1; i < numproc + 1; i++)
    id_lims[i] = id_lims[i - 1] + size_arr[i - 1];
  start = id_lims[myid];

  // get problem size by adding up
  for (i = 0; i < numproc; i++)
    size += size_arr[i];

  delete[]size_arr;
#else
  size = my_count;
  id_lims[myid] = 0;
  id_lims[myid + 1] = size;
  start = 0;
#endif
  int dims[2] = { size, 0 };
  double *buf = new double[my_count];
  int *ibuf = new int[my_count];

  j = 0;
  for (i = id_lims[myid]; i < id_lims[myid + 1]; i++)
  {
    ibuf[j] = i;
    j++;
  }

  // particle ids
  ierr = GH5_Write(gid, "ID", dims, (void *) ibuf, start, my_count, INTTYPE);

  // x-coordinates
  copy(x.begin(), x.end(), buf);
  ierr = GH5_Write(gid, "x", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // y-coordinates
  copy(y.begin(), y.end(), buf);
  ierr = GH5_Write(gid, "y", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // z-coordinates
  copy(z.begin(), z.end(), buf);
  ierr = GH5_Write(gid, "z", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // density of particles
  copy(rho.begin(), rho.end(), buf);
  ierr = GH5_Write(gid, "Rho", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // Velocity in x-directions
  copy(Vx.begin(), Vx.end(), buf);
  ierr = GH5_Write(gid, "Vx", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // Velocity in y-directions
  copy(Vy.begin(), Vy.end(), buf);
  ierr = GH5_Write(gid, "Vy", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // Velocity in z-directions
  copy(Vz.begin(), Vz.end(), buf);
  ierr = GH5_Write(gid, "Vz", dims, (void *) buf, start, my_count, DOUBLETYPE);

  //energy
  copy(engr.begin(), engr.end(), buf);
  ierr = GH5_Write(gid, "engr", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // mass fraction
  copy(mssfrc.begin(), mssfrc.end(), buf);
  ierr = GH5_Write(gid, "mssfrc", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // close file and group
  ierr = GH5_gclose(gid);

  // flush output
  ierr = H5Fflush(fp, H5F_SCOPE_LOCAL);
  ierr = H5Fclose(fp);

  // free memory
  delete[]id_lims;
  delete[]buf;
  delete[]ibuf;
  delete itr;

  return;
}


//function that call output sub_functions
void
write_particles_debug(int myid, int numprocs,
             THashTable * P_table, HashTable * BG_mesh,
             TimeProps * timeprops, int format, char *prefix)
{
  int bctp;
  if (format & 1)
  {
	  bctp = 0;
	  write_h5part_bctp(myid, numprocs, P_table, timeprops, bctp, prefix);

	  bctp = 1;
	  write_h5part_bctp(myid, numprocs, P_table, timeprops, bctp, prefix);

	  bctp = 2;
	  write_h5part_bctp(myid, numprocs, P_table, timeprops, bctp, prefix);

	  bctp = 100;
	  write_h5part_bctp(myid, numprocs, P_table, timeprops, bctp, prefix);
  }

//  if (format & 2)
//    write_matlab(myid, P_table, BG_mesh, timeprops, partition_table);

  return;
}

/*
 * function to debug the artificial viscosity code:
 */
void debug_vis ()
{
	const int DIMENSION2 = 2;
	const int nneigh=2;
	double rhoa = 5.83839246693786;
	double rhob[nneigh] = {4.6913, 8.7833};
	double Xa = 100.548573323395;
	double Xb[nneigh]= {69.3988, 69.3988};
	double Ya = 100.309592291433;
	double Yb[nneigh]= {146.056669986189, 7.25897673849247};
	double Ua = 11.3988024296422;
	double Ub[nneigh]= {8.40225864033622, 12.4756819006781};
	double Va = 23.9486170594628;
	double Vb[nneigh]= {164.466807083100, 140.394739074049};
	double CSa = 257.577586032606;
	double CSb[nneigh]= {732.888508148973, 715.288482521613};

	double rhoab;
	double sndspdab;
	double rab [DIMENSION2];
	double vab [DIMENSION2];
	double dist_sq;
	double h=107.675357293927;
    double vis;

	int i, j;

	for (i = 0; i < nneigh; i++)
	{
       rhoab =0.5*(rhoa + rhob[i]);
       sndspdab = 0.5 * (CSa + CSb[i]);

	   rab[0] = Xa-Xb[i];
	   rab[1] = Ya-Yb[i];

	   vab[0] = Ua-Ub[i];
	   vab[1] = Va-Vb[i];

	   dist_sq = 0;
       for (j = 0; j<DIMENSION2; j++)
    	   dist_sq += (rab[j]*rab[j]);

       vis= art_vis_2d (rhoab, sndspdab, rab, vab, dist_sq, h);
       cout << "artificial viscosity is :" << vis << endl;
	}
}

//overloading file for artificial viscosity
double art_vis_2d ( double rhoab, double sndspdab, double rab[2], double vab[2], double rsqab, double h)
{
	const int DIMENSION2 = 2;
	int    k;
	double miuab;
	double vrab = 0.;
	double vis = 0.;

	for (k = 0; k < DIMENSION2; k++)
	     vrab += (rab[k] * vab[k]);

	miuab = h * vrab / (rsqab + ata_P * h * h);

    if (vrab > 0)
		vis = 0.;
    else if (vrab < 0)
		vis = (beta_P * miuab - alf_P*sndspdab) * miuab / rhoab;

	return vis;
}

/*
 * function that used to check neighors of certain particles
 */
void check_neigh_part(THashTable * P_table)
{

    bool do_search = true;
    bool find;
    unsigned keycheck1[TKEYLENGTH] = {69674717, 3383906357, 0};
    unsigned keycheck2[TKEYLENGTH] = {69562537, 292385725, 0};
    unsigned keytemp[TKEYLENGTH] ;

    int i;
	THTIterator *itr = new THTIterator(P_table);
	Particle *p_curr = NULL;
    int phase;
    int num_diff_phase, num_same_phase;
	while ((p_curr = (Particle *) itr->next()))
	{
		if (do_search)
		{
		  	for (i = 0; i < TKEYLENGTH; i++)
		  		keytemp[i] = p_curr->getKey ().key[i];

		  	if (find_particle (keytemp, keycheck2))
		  	{
		  		 // list of neighbors
		  		vector < TKey > pneighs = p_curr->get_neighs();
		  		vector < TKey >::iterator p_itr;
		  		phase = p_curr->which_phase();
		  		num_diff_phase = 0;
		  		num_same_phase = 0;
		  		for (p_itr = pneighs.begin(); p_itr != pneighs.end(); p_itr++)
		  		{

		  			Particle *pj = (Particle *) P_table->lookup(*p_itr);
		  			assert (pj);

		  			if (phase == pj->which_phase())
		  				num_same_phase ++;
		  			else
		  				num_diff_phase ++;
		  		}
		  		cout << "----------particle of phase2-----:" << endl;
		  		cout << "number of same phase: " << num_same_phase << "number of different phase: " << num_diff_phase<< endl;
		  	}

		  	if (find_particle (keytemp, keycheck1))
		  	{
		  		 // list of neighbors
		  		vector < TKey > pneighs = p_curr->get_neighs();
		  		vector < TKey >::iterator p_itr;
		  		phase = p_curr->which_phase();
		  		num_diff_phase = 0;
		  		num_same_phase = 0;
		  		for (p_itr = pneighs.begin(); p_itr != pneighs.end(); p_itr++)
		  		{

		  			Particle *pj = (Particle *) P_table->lookup(*p_itr);
		  			assert (pj);

		  			if (phase == pj->which_phase())
		  				num_same_phase ++;
		  			else
		  				num_diff_phase ++;
		  		}
		  		cout << "----------particle of phase1-----:" << endl;
		  		cout << "number of same phase: " << num_same_phase << "number of different phase: " << num_diff_phase<< endl;
		  	}
		 }
	 }//end of go through all particles

}

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
