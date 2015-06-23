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


//add air particles
void
add_air (
		//P_table
		THashTable *,
		//BG_table
		HashTable * ,
        //Mat property
        MatProps *,
        //number of processor
        int ,
        //my ID
        int
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
               TimeProps *     //! struct for simulation props
              );

//! Read simulation data \f${\it i.e.}\f$ simulation time, output format, pile pros etc
int Read_Data (
               MatProps *,     //! Structure containg material mroperties
               TimeProps *,    //! Structure containg time mroperties
               int *           //! output file format
              );

//! Read background grid created by preprocessor
int Read_Grid (
               THashTable **,   //! Pointer to Hashtable for background mesh
               HashTable **,   //! Pointer to HashTable for partilces
               vector<BucketHead> &,  //! Vector of sorted partition table keys
               MatProps  *,    //! Structure containing material properties
               int ,           //! my process id
               int ,           //! number of total processes
               int *           //! array of flags, for cummnication with other procs
              );

//! Search neighbors for their proximity
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

//! Calculate time increment, depending upon CFL condition
double timestep(
                THashTable *,    //! HashTable of SPH partilces
                MatProps *      //! Structure of material properties
               );

//! Update particle positions and their relationship with background Mesh
int update_pos(
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

//set up eruption particles
int setup_erupt(
		        int, //myid
		        THashTable * ,  //P_table
		        HashTable * ,  //BG_mesh
                TimeProps * ,  //timeProps
                MatProps * ,   //matprops
                int            //Number of processor
                );

//function for adding new ghost erupt particles at the bottom of the duck
void add_new_erupt(
		          int ,          //myid
		          THashTable *,  //particle hash table
		          HashTable *,   //background mesh table
                  TimeProps *,   //Time propos
				  MatProps *,    //matprops
                  double         //Time steps
                  );
//set up initial air
int
setup_ini(
		  int,         //myid
		  THashTable *, //P_table
		  HashTable *, //BG_mesh
          TimeProps *, //timeProps
          int          //Number of processor
          );

//determine time step
double
timestep (
		  THashTable *, //particle table
		  TimeProps *    //timeprops
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
