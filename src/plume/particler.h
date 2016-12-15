/*
 * particler.h
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */

#ifndef PARTICLER_H
#define PARTICLER_H

#include <vector>
#include <list>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bnd_image.h>
#include <multiproc.h>
#include <properties.h>
#include "options.h"

//add air particles
void
add_air (
		//P_table
		THashTable *,
		//BG_table
		HashTable * ,
        //Mat property
        MatProps *,
		// simprops
		SimProps*,
        //number of processor
        int ,
        //my ID
        int
       );

//add pressure ghost particles
void
add_pressure_ghost(
		//P_table
		THashTable *,
		//BG_table
		HashTable * ,
		//simprops
		SimProps *,
        //Mat property
        MatProps *,
        //TimeProps
        TimeProps *,
        //number of processor
        int ,
        //my ID
        int
       );

//add wall ghost particles
void
add_wall_ghost(
		//P_table
		THashTable *,
		//BG_table
		HashTable * ,
		//simprops,
		 SimProps* ,
        //Mat property
        MatProps *,
        //TimeProps
        TimeProps *,
        //number of processor
        int ,
        //my ID
        int
       );

/*
 * What happened before this function calling is: scan the most outside potential involved bucket, if any of particles becomes involved (two ways), turn the bucket to be has_involved.
 * function which scans these buckets which were originally pressure ghost buckets, if any of its neighbor buckets is has_involved, turn it to be a potential_involved_bucket.
 * What will happen after this function calling is: add a new layer of pressure ghost bucket as the old pressure ghost bucket become potential_involved.
 */
void adapt_domain(
		THashTable * , // P_table
		HashTable * ,  // BG_mesh
		MatProps * ,   //matprops,
		int ,          // numproc
		int            // myid
		);

//! Apply boundary conditions
int apply_bcond (
                 int        ,  //! Process ID
                 THashTable *,  //! Hashtable of SPH partilces
                 HashTable *,  //! Background mesh
                 MatProps *,    //! Structure of material properties
                 list <BndImage> & //! table of ghost images
                );

//! Momentum update function. This where most of the work is done
int mom_engr_update(
               int ,           //! my Proc ID
               THashTable *,    //! HashTable of SPH partilces
               HashTable *,    //! HashTable of Mesh elements
               TimeProps *,    //! struct for simulation props
			   SimProps *      //simprops
              );

//! Read simulation data \f${\it i.e.}\f$ simulation time, output format, pile pros etc
int Read_Data (
               MatProps *,     //! Structure containg material mroperties
               TimeProps *,    //! Structure containg time mroperties
               SimProps *,     //simulation properties
               int *           //! output file format
              );

//! Read background grid created by preprocessor
int Read_Grid (
               THashTable **,   //! Pointer to Hashtable for background mesh
               HashTable **,   //! Pointer to HashTable for partilces
               vector<BucketHead> &,  //! Vector of sorted partition table keys
               MatProps  *,    //! Structure containing material properties
               SimProps *,     //simulation properties
               int ,           //! my process id
               int ,           //! number of total processes
               int *           //! array of flags, for cummnication with other procs
              );

// update neighbor information without updating smooth length --> Will be called at the time of particle set up;
int  search_neighs_consth(
                         int myid   , //! My processor ID
                         THashTable *, //! HashTable of SPH partilces
                         HashTable *  //! HashTable of cells of background mesh
                         );

/*
 * First: search for neighbors --->to guarantee that smooth_update is using the most recent information;
 * Second: update smooth length --> without loop
 * Third: search for neighbors;
 * This will be called in the main loop of SPH update in time;
 */
int  search_neighs(
                         int myid   , //! My processor ID
                         THashTable *, //! HashTable of SPH partilces
                         HashTable *  //! HashTable of cells of background mesh
                         );

//function that used to set up initial atmosphere and determine the mass of air particles
int setup_ini(
		int ,          //myid
		THashTable * , //P_table
		HashTable * ,  //BG_mesh
        TimeProps * ,  //timeprops
        int ,          //numprocs
        int*           //my_comm
        );

//! Smooth density (low-pass filter)
void smooth_density (
                     THashTable * //! HashTable of SPH partilces
                    );

//! Smooth velocity --> for SPH-epsilon
void smooth_velocity (
                     THashTable * //! HashTable of SPH partilces
                    );

//! Update particle positions and their relationship with background Mesh
void update_pos(
               int ,         //! my proc id
               THashTable *, //! HashTable of SPH particles
               HashTable *,  //! HashTable of cells of background mesh
               TimeProps *,  //! Time properties struct
               MatProps * ,  //! Mat properties struct
               int *         //! Pointer to number of particles removed
              );

//! Write output, in specified format
void write_output(
                  int ,        //! my process id
                  int ,        //! total number of procs
                  THashTable *, //! HashTable of SPH particles
                  HashTable *, //! HashTable of Background Cells
                  vector<BucketHead> &, //! table of partitioned keys
                  TimeProps *, //! Structure time properties
                  int     //! file format, {\it i.e.} hdf5, tecplot etc
                 );

//! Write output, in specified format ---> output data without ghost or guest particles, clean for show the results.
void write_output_show(
                  int ,        //! my process id
                  int ,        //! total number of procs
                  THashTable *, //! HashTable of SPH particles
                  HashTable *, //! HashTable of Background Cells
                  vector<BucketHead> &, //! table of partitioned keys
                  TimeProps *, //! Structure time properties
                  int     //! file format, {\it i.e.} hdf5, tecplot etc
                 );

#ifndef SIMULATE_ASH
//set up eruption particles
int setup_erupt(
		        int, //myid
		        THashTable * ,  //P_table
		        HashTable * ,  //BG_mesh
                TimeProps * ,  //timeProps
                MatProps * ,   //matprops
				SimProps *,    //simprops
                int            //Number of processor
                );

//function for adding new ghost erupt particles at the bottom of the duck
void add_new_erupt(
		          int ,          //myid
		          THashTable *,  //particle hash table
		          HashTable *,   //background mesh table
                  TimeProps *,   //Time propos
				  MatProps *,    //matprops
				  SimProps *,    //simprops
                  double         //Time steps
                  );

#else
//set up influx particles
int setup_influx(
		        int, //myid
		        THashTable * ,  //P_table
		        HashTable * ,  //BG_mesh
                TimeProps * ,  //timeProps
                MatProps * ,   //matprops
				SimProps *,    //simprops
                int            //Number of processor
                );

//add influx particles
int add_new_influx(
		        int, //myid
		        THashTable * ,  //P_table
		        HashTable * ,  //BG_mesh
                TimeProps * ,  //timeProps
                MatProps * ,   //matprops
				SimProps *,    //simprops
                int            //Number of processor
                );

#endif

//set up initial air
int
setup_ini(
		  int,         //myid
		  THashTable *, //P_table
		  HashTable *, //BG_mesh
          TimeProps *, //timeProps
		  SimProps * , // simpros,
          int          //Number of processor
          );

//determine time step
double
timestep (

		  HashTable *, //BG_mesh
		  THashTable *, //particle table
		  TimeProps *    //timeprops
         );

//determine time step that can stable fluctuation near the boundary
double
pressure_bc_step(
		HashTable *, // BG_mesh,
		THashTable * //P_table
		);

//scan the most outside layer of buckets satisfy has_involved>0,
//if they have any involved particles (involved = 2), then make the most_out_side layer to be has_involved.
//What will happen next is the domain will be adjust such that at least one layer of potential involved bucket is added outside the domain.
int
scan_outside_layer
        (
		THashTable * , // P_table,
		HashTable *,   // BG_mesh,
		int ,          //numproc,
		int            //myid
		);

//function that shif the layer outside the pressure ghost boundary from brief bucket to bucket
void
shift_brief_buck (
		HashTable *, // BG_mesh,
		MatProps *,  // matprops,
		TimeProps *, //timeprops,
		int          //myid
		);


////This function is useless should be removed.
//int
//put_ghost_particles (
//		// P_table
//		THashTable * ,
//		//BG_mesh
//		HashTable * ,
//		// timeprops
//		TimeProps *,
//		// partition_table
//        vector<BucketHead> & ,
//        // matprops
//        MatProps * ,
//        // my_comm
//        int *,
//        //myid
//        int ,
//        //numproc
//        int ,
//        //added_ghosts
//        int *
//        );

#endif /* PARTICLER_H_ */
