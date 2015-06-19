/*
 * particle_deb.cc
 *
 *  Created on: Jun 1, 2015
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

#include "sph_header.h"
#include "constant.h"
#include "parameters.h"

#include <hdf5.h>
#include "hdf5calls.h"

inline void print_particle(Particle * pnew)
{
	int i;
	cout << "key = {" << endl;
	for (i=0; i<TKEYLENGTH; i++ )
	   cout << pnew->getKey().key[i] << "  " <<endl;
	cout << "} \n" << endl;

	cout << "is_erupt_ghost=" << pnew->is_erupt_ghost () <<   endl;
	cout << "is_press_ghost=" << pnew->is_press_ghost () <<   endl;
	cout << "is_wall_ghost=" << pnew->is_wall_ghost () <<   endl;
	cout << "is_real=" << pnew->is_real () <<   endl;

	cout << "is_not_updated=" << pnew->is_not_updated () <<   endl;
	cout << "has_reflection=" << pnew->has_reflection () <<   endl;

	cout << "is_guest=" << pnew->is_guest () <<   endl;
	cout << "new_old=" << pnew->get_new_old () <<   endl;

	cout << "mass=" << pnew->get_mass () <<   endl;

	cout << "need_neigh=" << pnew->need_neigh () <<   endl;
	cout << "contr_image=" << pnew->contr_image() <<   endl;
	cout << "contr_dens=" << pnew->contr_dens () <<   endl;

}

inline void create_part (double *pcrd, int myid, double prss, int phs_num,
                        double masfrc, double gmm, double sndspd,
                        unsigned add_step, double *mindom, double *maxdom,
                        double mass, double smlen, THashTable *P_table, Particle * pnew)
{
	int i, ii;
	unsigned tkeylen = TKEYLENGTH;
	double bnd[2*DIMENSION], flag[DIMENSION],index[2*DIMENSION];
	double normc[DIMENSION];
	unsigned pkey[TKEYLENGTH];
	int bctp;
	int ptype;

	  /*temporarily use the following way to get bnd, Finally, I will put bnd either
	   * in SimulProps or parameters.h
	   * will modify this part later
	  */
     bnd[0]=Lx_P[0];
     bnd[1]=Lx_P[1];
     bnd[2]=Ly_P[0];
     bnd[3]=Ly_P[1];
     bnd[4]=Lz_P[0];
     bnd[5]=Lz_P[1];

	ptype = 0; //is real
	for (i = 0; i < DIMENSION; i++)
		 cout << *(pcrd+i) << "\t" << endl;
	cout << "\n" << endl;

	for (ii = 0; ii < DIMENSION; ii++)
		flag[ii]=index[2*ii]+index[2*ii+1];

	for (ii = 0; ii < DIMENSION; ii++)
	     if ((pcrd[ii]<bnd[2*ii]) || (pcrd[ii]>bnd[2*ii+1]))
	          ptype = 1;//ptype is now an index for determine bc_type

	//determine bc_type
	if (ptype ==0)
	   	bctp = 100;
	else
	{
	   	if ( pcrd[2]<bnd[4] )
	   	     bctp = 2; //wall ghost particle
	    else
	       	 bctp = 1 ;// pressure bc ghost particle
	 }
	for (ii = 0; ii < DIMENSION; ii++)
		normc[ii] = (pcrd[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);

	THSFC3d (normc, add_step, &tkeylen, pkey);
	// check for duplicates
	if (P_table->lookup(pkey))
	{
		fprintf(stderr, "ERROR: Trying to add particle "
			  	    	"twice on same location.\n");
//		exit(1);
	}
    pnew = new Particle(pkey, pcrd, mass, smlen, prss, masfrc, gmm, sndspd, phs_num, myid, bctp);

	// add to hash-table
	P_table->add(pkey, pnew);

	//
	//out put results
	print_particle(pnew);
}

void
particle_deb (int myid)
{

	  bool debug_particles = false;
	  if (debug_particles)
	  {
		  // Create hash-table for particles
		  int P_TABLE_SIZE = 400;
		  double mindom[DIMENSION], maxdom[DIMENSION];

		  // Read Hash table constants
		  // Read Hash-table related constants
		  char filename[14];
		  double hvars[6];
		  sprintf(filename, "funky%04d.h5", myid);
		  hid_t fp = GH5_fopen_serial(filename, 'r');
		  GH5_readdata(fp, "/hashtable_constants", hvars);
		  mindom[0] = hvars[0];         // min x
		  maxdom[0] = hvars[1];         // max x
		  mindom[1] = hvars[2];         // min y
		  maxdom[1] = hvars[3];         // max y
		  mindom[2] = hvars[4];         // min z
		  maxdom[2] = hvars[5];         // max z
		  THashTable *P_table = new THashTable(P_TABLE_SIZE, 2017, mindom, maxdom);

		  int i, j, k, ii;
		  unsigned pkey[TKEYLENGTH];
		  unsigned tkeylen = TKEYLENGTH;
		  double normc[DIMENSION];
		  double pcrd[DIMENSION];
		  double bnd[2*DIMENSION], flag[DIMENSION],index[2*DIMENSION];
		  // start putting piles
		  double mass = 5400;
		  double smlen = 500;
		  double dx = smlen;
		  double dx2 = 0.5 * dx;

		  //Initial particle property
		  int bctp; //boundary condition ghost particle type.
		  int ptype = 0; //real particle
		  double prss = 0.;
		  double masfrc = 0.;
		  double gmm = 1.4;
		  double sndspd = 340.;
		  int phs_num = 1; //phase 1, air particle
		  unsigned add_step=0; // eruption particle adding step
		  Particle * pnew;

		  cout << "*************************test determine particle type******************** \n" << endl;
		  cout << "------case1: real\n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = 500;
		  pcrd[2] = 500;

          //create particles
		  create_part (pcrd, myid, prss, phs_num, masfrc, gmm, sndspd, add_step,
                       mindom, maxdom, mass, smlen, P_table, pnew);

		  cout << "------case2: wall ghost \n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = 500;
		  pcrd[2] = -500;
          //create particles
		  create_part (pcrd, myid, prss, phs_num, masfrc, gmm, sndspd, add_step,
		                         mindom, maxdom, mass, smlen, P_table, pnew);
		  //out put results
		  print_particle(pnew);

		  cout << "------case3: wall ghost\n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = -7500;
		  pcrd[2] = -2500;

          //create particles
		  create_part (pcrd, myid, prss, phs_num, masfrc, gmm, sndspd, add_step,
		                         mindom, maxdom, mass, smlen, P_table, pnew);

		  cout << "------case4: pressure bc\n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = -7500;
		  pcrd[2] = 500;

          //create particles
		  create_part (pcrd, myid, prss, phs_num, masfrc, gmm, sndspd, add_step,
		                         mindom, maxdom, mass, smlen, P_table, pnew);

		  cout << "------case5: pressure bc \n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = 500;
		  pcrd[2] = 3800;

          //create particles
		  create_part (pcrd, myid, prss, phs_num, masfrc, gmm, sndspd, add_step,
		                         mindom, maxdom, mass, smlen, P_table, pnew);
		  ;
		  cout << "------case6: Add particle the same position\n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = 500;
		  pcrd[2] = 500;

          //create particles
		  create_part (pcrd, myid, prss, phs_num, masfrc, gmm, sndspd, add_step,
		                         mindom, maxdom, mass, smlen, P_table, pnew);


		  cout << "*************************test update second variable******************** \n" << endl;
		  cout << "------case1: real\n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = 500;
		  pcrd[2] = 500;
		  //create particles

		  cout << "*************************test erupt constructor******************** \n" << endl;
		  cout << "------case1: real\n" << endl;
		  ptype = 0; // is real
		  pcrd[0] = 400;
		  pcrd[1] = 500;
		  pcrd[2] = 500;
		  //create particles
	  }
	  return;
}
