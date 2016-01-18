/*
 * debug_header.h
 *
 *  Created on: Jun 1, 2015
 *      Author: zhixuanc
 */
#ifndef DEBUG_HEADER_H
#define DEBUG_HEADER_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <cmath>

using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bgmesh.h>
#include <bnd_image.h>
#include <properties.h>
#include <buckhead.h>
#include <hilbert.h>
#include <bucket.h>
#include <particle.h>
#include <multiproc.h>
#include <properties.h>

#include <constant.h>
#include <IndMap.h>
#include <parameters.h>

void
particle_deb (
		//myid
		int
		);

/*find particle by the key
 * keyin is input particle key
 * keycheck is the key that is given and all keyin will be compared with keycheck
  */
bool find_particle (
		unsigned*, // keyin
		unsigned* //keycheck
		);

/*find particle by the pos
 *
 */
bool find_particle (
		double*, // kin
		double* //check
		);

/*find bucket by the key
 * keyin is input particle key
 * keycheck is the key that is given and all keyin will be compared with keycheck
  */
bool find_bucket (
		unsigned*, // keyin
		unsigned* //keycheck
		);

/*find bucket by the pos
 *
 */
bool find_bucket (
		double*, // kin
		double* //check
		);


//function that used to go through all buckets to check whether some buckets in the BG_mesh table is wired or not.
void BG_mesh_check (
		HashTable *  //BG_mesh
		);

//function that used to check mesh error
void BG_mesh_err_check(
		HashTable * , // BG_mesh
		int           // myid
		);


/*find particle by the range of particle position
 *
 */
bool find_particle_pos_range  (double* in, double* check);


//function to find particle with given key, will be useful in debugging.
void check_particle_bykey (
		THashTable *
        );

//overload function to find particle with given key, will be useful in debugging. --> output my process
bool check_particle_bykey (
		THashTable *,// P_table
		int* //id
		);


//function to find bucket with given key, will be useful in debugging.
void check_bucket_bykey (
		HashTable *
        );

//find the particle with non-physical density
void find_large_density_particle (
		THashTable *
        );

//function that call output sub_functions
void
write_particles_debug(
		     int , //myid
		     int , //numprocs
             THashTable * , //P_table
             HashTable * ,  //BG_mesh
             TimeProps * , //timeprops
             int,  //format
             char *//prefix
             );

//output certain type of particle
void
write_h5part_bctp(
		int , //myid
		int , //numproc
		THashTable * , //P_table
		TimeProps * ,  //timepros
		int ,//bctp
		char *//prefix
		);

//go through all particle and check their type!
void check_particle_all_type (
		THashTable * //P_table
		);

//check particles in certain region
void check_particle_bypos (
		THashTable * //P_table
		);

/*find all particle of given pahse
 * This will be used to find particle of phase 2 --> get its key and then track the movement of that particle
 */
void find_particle_by_phase (
		THashTable *  //P_table
        );

/*
 * Find the highest position of all phase2 particles
 */
void find_highest_z (
		THashTable * //P_table
		);


//overloading file for artificial viscosity
double art_vis_2d (
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


/*
 * function to debug the artificial viscosity code:
 */
void debug_vis ();


/*
 * function that used to check neighors of certain particles
 */
void check_neigh_part(
		THashTable *// P_table
		);


//function to find particle with given key, will be useful in debugging.
void check_bucket_guest (
		HashTable * //BG_mesh
		);

//function to check where does the negative sound speed comes from
bool check_particles_sndspd (
    	       THashTable * /*particle hash table*/
               );

#endif /* DEBUG_HEADER_H */
