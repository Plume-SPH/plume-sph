/*
 * main.cc
 *
 *  Created on: Mar 5, 2015
 *      Author: zhixuanc
 */

#include <iostream>
#include <vector>
#include <list>
using namespace std;

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef MULTI_PROC
#  include <mpi.h>
#  include <multiproc.h>
#endif

#ifdef DEBUG
#  include <debug_header.h>
#endif

#include <hashtab.h>
#include <thashtab.h>
#include <bgmesh.h>
#include <bnd_image.h>
#include <properties.h>
#include <buckhead.h>

#include "sph_header.h"
#include "particler.h"


int
main(int argc, char **argv)
{
  int i, j, ierr = 0;
  double dt;
  MatProps *matprops = new MatProps();
  TimeProps *timeprops = new TimeProps();
  SimProps *simprops = new SimProps();
  THashTable *P_table;
  HashTable *BG_mesh;

  vector < BucketHead > partition_table;
//  vector < OutSideBucket > outside_table;
  list < BndImage > Image_table;
  int format = 0;
  int adapt = 0;
  int lost = 0;
  int lostsum = 0;
  double local_data[3], global_data[3];
  int myid, numprocs;

#ifdef MULTI_PROC
  double start, finish;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  GMFG_new_MPI_Datatype();
  start = MPI_Wtime();
#else
  numprocs = 1;
  myid = 0;
#endif

#ifdef DEBUG
  bool check_part = true;
  int  id;
  bool find = false;
  bool check_buck = false;
  bool check_part_tp =false;
  bool check_bypos = false;
  bool find_large_density = false;
  bool search_byphase = false;
  bool find_maxz = false;
  bool vis_flag = false;
  bool check_neigh =false;
  bool check_bucket = false;
#endif

  // allocate communcation array
  int * my_comm = new int [numprocs];

  //Read input data
  if (Read_Data(matprops, timeprops, simprops, &format) != 0)
  {
    cerr << "ERROR: Can't read input data\n";
    exit(1);
  }

  //read initial particle distribution
  if (Read_Grid (&P_table, &BG_mesh, partition_table, matprops,
		  simprops, myid, numprocs, my_comm) != 0)
  {
    cerr << "ERROR: Can't read Initial grid\n";
    exit(1);
  }

  //add air particles and put particles into bucket, bc_type is determined in this process
  add_air (P_table, BG_mesh, matprops, numprocs, myid);

#ifdef DEBUG
  if (check_bucket)
	  check_bucket_guest (BG_mesh);
#endif

  // scan mesh and mark buckets active/inactive
  update_bgmesh (BG_mesh, myid, numprocs, my_comm);

#ifdef DEBUG
  if (check_part)
  {
	  find = false;
	  find = check_particle_bykey (P_table, &id);
	  if (find)
	     cout <<"its myid is : " << id << endl;
  }
#endif
  // sync data
  move_data (numprocs, myid, my_comm, P_table, BG_mesh);

#ifdef DEBUG
  if (check_part)
  {
	  find = false;
	  find = check_particle_bykey (P_table, &id);
	  if (find)
	     cout <<"its myid is : " << id << endl;
  }
#endif

//The reason why I put the following functions late is that adding of pressure_ghost and adding of wall ghost need neighbor info.
  //add pressure ghost
   add_pressure_ghost(P_table, BG_mesh, matprops, timeprops, numprocs, myid);

   //add wall ghost
   add_wall_ghost(P_table, BG_mesh, matprops, timeprops, numprocs, myid);

   //initialized the mass of all particles---> actually, only initial air need to set up initial! wall ghost and pressure ghost will set up parameter while adding particles
//   setup_ini(myid,  P_table,  BG_mesh, timeprops, numprocs, my_comm);

   //Adding eruption boundary condition
   setup_erupt(myid, P_table, BG_mesh, timeprops, matprops, numprocs);


#ifdef MULTI_PROC
  // wait till initialization has finished
  MPI_Barrier (MPI_COMM_WORLD);

  // remove guest buckets and particles
  delete_guest_buckets (BG_mesh, P_table);

  // Initial repartition
  repartition (partition_table, P_table, BG_mesh, my_comm);

  // move data
  move_data (numprocs, myid, my_comm, P_table, BG_mesh);
#endif

  // search and update neighbors ---> can I use a send neighbor here? --> not, I do not think it is necessary!
  search_neighs_consth (myid, P_table, BG_mesh);

  // search mirror imgaes of ghost particles into boundary
  search_bnd_images(myid, P_table, BG_mesh, Image_table, 1);

//  //set up initial outside bucket layer list
//  update_out_layer (P_table,BG_mesh, outside_table, numproc, myid);

#ifdef DEBUG
  if (check_part)
  {
	  find = false;
	  find = check_particle_bykey (P_table, &id);
	  if (find)
	     cout <<"its myid is : " << id << endl;
  }
#endif

#ifdef DEBUG
  if (check_bucket)
	  check_bucket_guest (BG_mesh);
#endif


#ifdef MULTI_PROC
  // send reflections that belong to other procs
  send_foreign_images (myid, numprocs, BG_mesh, Image_table, my_comm);
#endif

  // apply boundary conditions
  apply_bcond(myid, P_table, BG_mesh, matprops, Image_table);

  // Write inital configuration
  write_output (myid, numprocs, P_table, BG_mesh,
                partition_table, timeprops, format);

  /*
   *
   * The time-stepping loop
   *
   */
  while (!timeprops->ifend())
  {
    // calculate time-step
    dt = timestep(P_table, timeprops);

#ifdef MULTI_PROC
    // get-minimum time-step for multiproc run
    local_data[0] = dt;
    local_data[1] = (double) -ierr;
    local_data[2] = (double) -adapt;
    MPI_Allreduce(local_data, global_data, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = global_data[0];

    // increment time
    timeprops->incrtime(&dt);
    if (myid == 0)
      cout << "Time-step: " << timeprops->step << " dt=" << dt
        << " time=" << timeprops->timesec() << endl;

    ierr = (int) round (global_data[1]);
    if ( ierr )
      MPI_Abort(MPI_COMM_WORLD, ierr);

    // update mesh
    adapt = (int) round (global_data[2]);
    if (adapt)
    {
      //sync data, in outside layer scanning, guest bucket is not able to updated!
      move_data (numprocs, myid, my_comm, P_table, BG_mesh);

      // adapt_domain -> What did here is make the domain larger--->because the involved particle enter the most outside bucket layer (has_potential_involved = 1;)
      adapt_domain (P_table, BG_mesh, numprocs, myid);

      //sync data, as some ghost_pressure bucket becomes has_potential_involved (also for the particles inside)-->these information is needed when adding pressure ghost and as well as wall ghost.
      move_data (numprocs, myid, my_comm, P_table, BG_mesh);

      //add pressure ghost --> will not delete the old pressure ghost, just add pressure ghost where computational domain is "exposed"
      add_pressure_ghost (P_table, BG_mesh, matprops, timeprops, numprocs, myid);

      //add wall ghost
      add_wall_ghost(P_table, BG_mesh, matprops, timeprops, numprocs, myid);

      // sync data again
      move_data (numprocs, myid, my_comm, P_table, BG_mesh);

#ifdef DEBUG
  if (check_buck)
	  check_bucket_bykey  (BG_mesh);
#endif

#ifdef DEBUG
  if (check_part)
	  check_particle_bykey (P_table);
#endif

      // scan buckets and make them active / inactive
      update_bgmesh (BG_mesh, myid, numprocs, my_comm);

      // sync data again
      move_data (numprocs, myid, my_comm, P_table, BG_mesh);

      // search reflections for new ghost particles
      search_bnd_images (myid, P_table, BG_mesh, Image_table, 0);

      // send the images that belong on other subdomains
      send_foreign_images (myid, numprocs, BG_mesh, Image_table, my_comm);

      // apply boundary conditions--> as some new wall ghost added, the boundary condition need to be updated!
      ierr += apply_bcond (myid, P_table, BG_mesh, matprops, Image_table);

//      //set up initial outside bucket layer list
//      update_out_layer (P_table,BG_mesh, outside_table, numproc, myid);
    }
#endif
    ierr = 0;  // reset error code
    adapt = 0; // and adapt flag

    // search and update neighbors
    search_neighs_consth (myid, P_table, BG_mesh);

#ifdef MULTI_PROC
    // ghost need to be updated only before updating momentum
    move_bnd_images (myid, numprocs, P_table, BG_mesh, Image_table);
#endif

    // update momentum and energy
    int err2 = mom_engr_update (myid, P_table, BG_mesh, timeprops);
    if ( err2 )
    {
      cerr << "Momentum update failed on proc" << myid <<
        " at time-step : " << timeprops->step << endl;
      cerr << "Check outfile for proc " << myid << " for errors." << endl;
      ierr += err2;
    }

#ifdef MULTI_PROC
    // update guests on all procs
    move_data(numprocs, myid, my_comm, P_table, BG_mesh);
#endif

    // smooth out density oscillations (if any)
    smooth_density(P_table);

#ifdef MULTI_PROC
    // update guests on all procs
    move_data(numprocs, myid, my_comm, P_table, BG_mesh);
#endif

#ifdef DEBUG
  if (check_part)
  {
	  find = false;
	  find = check_particle_bykey (P_table, &id);
	  if (find)
	     cout <<"its myid is : " << id << endl;
  }
#endif

    // update particle positions
    update_pos (myid, P_table, BG_mesh, timeprops, matprops, &lost);

    // add new layers of particle in the duct
    add_new_erupt(myid, P_table, BG_mesh, timeprops, matprops, dt);

#ifdef MULTI_PROC
    // update guests as density has changed since last update
    move_data(numprocs, myid, my_comm, P_table, BG_mesh);

    //scan most outside layer of has_potential_involved buckets
    adapt = scan_outside_layer (P_table, BG_mesh, numprocs, myid);

    /* DYNAMIC LOAD BALANCING */
    if ((numprocs > 1) && (timeprops->step % 350 == 0))
    {
      // remove guest buckets and particles
      delete_guest_buckets (BG_mesh, P_table);

      // repartition the domain
      i = repartition (partition_table, P_table, BG_mesh, my_comm);

      // send new guests
      move_data (numprocs, myid, my_comm, P_table, BG_mesh);

      // make buckets active / inactive, if needed
      update_bgmesh (BG_mesh, myid, numprocs, my_comm);

      // send new guests
      move_data (numprocs, myid, my_comm, P_table, BG_mesh);

      // search ghost refections
      search_bnd_images (myid, P_table, BG_mesh, Image_table, 1);

      // send reflections that belong to other procs
      send_foreign_images (myid, numprocs, BG_mesh, Image_table, my_comm);

//      // remove ghosts not needed anymore
//      delete_unused_ghosts (P_table, BG_mesh, myid);
    }
#endif

    // apply boundary conditions
    ierr += apply_bcond (myid, P_table, BG_mesh, matprops, Image_table);

    // write output if needed
    if (timeprops->ifoutput())
      write_output (myid, numprocs, P_table, BG_mesh,
                    partition_table, timeprops, format);
  }

#ifdef MULTI_PROC
  move_data (numprocs, myid, my_comm, P_table, BG_mesh);
#endif

#ifdef MULTI_PROC
  finish = MPI_Wtime();
  MPI_Reduce (&lost, &lostsum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Finalize ();
#endif

  // just for the sake of good practice
  delete [] my_comm;

  if (myid == 0)
  {
    printf("A total of %d SPH Particles were lost.\n", lostsum);
    double walltime = finish - start;
    int hours = (int) (walltime / 3600.);
    int mins = (int) ((walltime - hours * 3600) / 60);
    double secs = walltime - (hours * 3600) - (mins * 60);

    printf ("Computation time for a %d proc run was %d:%02d:%f\n",
           numprocs, hours, mins, secs);
  }

  return 0;
}
