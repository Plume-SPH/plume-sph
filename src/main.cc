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
  THashTable *P_table;
  HashTable *BG_mesh;

  vector < BucketHead > partition_table;
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
  bool check_part = false;
  bool check_part_tp =false;
  bool check_bypos = false;
  bool find_large_density = false;
#endif

  // allocate communcation array
  int * my_comm = new int [numprocs];

  //Read input data
  if (Read_Data(matprops, timeprops, &format) != 0)
  {
    cerr << "ERROR: Can't read input data\n";
    exit(1);
  }

  //read initial particle distribution
  if (Read_Grid (&P_table, &BG_mesh, partition_table, matprops,
		  myid, numprocs, my_comm) != 0)
  {
    cerr << "ERROR: Can't read Initial grid\n";
    exit(1);
  }

  //add air particles and put particles into bucket, bc_type is determined in this process
  add_air (P_table, BG_mesh, matprops, numprocs, myid);

#ifdef DEBUG
  if (check_part)
	  // find certain particle and check its values
      check_particle_bykey (P_table);
#endif

  // scan mesh and mark buckets active/inactive
  update_bgmesh (BG_mesh, myid, numprocs, my_comm);

  // sync data
  move_data (numprocs, myid, my_comm, P_table, BG_mesh);

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

  // search and update neighbors
  search_neighs (myid, P_table, BG_mesh);

  //initialized the mass of all particles
  setup_ini(myid,  P_table,  BG_mesh, timeprops, numprocs, my_comm);

  //Adding eruption boundary condition
  setup_erupt(myid, P_table, BG_mesh, timeprops, matprops, numprocs);

#ifdef DEBUG
  if (check_part_tp)
	  check_particle_all_type (P_table);
#endif

#ifdef DEBUG
  if (check_bypos)
	  check_particle_bypos (P_table);
#endif

  // scan mesh and mark buckets active/inactive
  update_bgmesh (BG_mesh, myid, numprocs, my_comm);

  // sync data again
  move_data (numprocs, myid, my_comm, P_table, BG_mesh);

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

  // search mirror imgaes of ghost particles into boundary
  search_bnd_images(myid, P_table, BG_mesh, Image_table, 1);

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
      // scan buckets and make them active / inactive
      update_bgmesh (BG_mesh, myid, numprocs, my_comm);

      // sync data again
      move_data (numprocs, myid, my_comm, P_table, BG_mesh);

      // search reflections for new ghost particles

      search_bnd_images(myid, P_table, BG_mesh, Image_table, 0);

      // send the images that belong on other subdomains
      send_foreign_images (myid, numprocs, BG_mesh, Image_table, my_comm);
    }
#endif
    ierr = 0;  // reset error code
    adapt = 0; // and adapt flag

    // search and update neighbors
    search_neighs (myid, P_table, BG_mesh);

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

#ifdef DEBUG
  if (find_large_density)
	  // find certain particle and check its values
	  find_large_density_particle (P_table);
#endif

#ifdef MULTI_PROC
    // update guests on all procs
    move_data(numprocs, myid, my_comm, P_table, BG_mesh);
#endif


    // update particle positions
    adapt = update_pos (myid, P_table, BG_mesh, timeprops, &lost);

    // add new layers of particle in the duct
    add_new_erupt(myid, P_table, BG_mesh, timeprops, matprops, dt);

#ifdef MULTI_PROC
    // update guests as density has changed since last update
    move_data(numprocs, myid, my_comm, P_table, BG_mesh);

    /* DYNAMIC LOAD BALANCING */
    if ((numprocs > 1) && (timeprops->step % 2500 == 0))
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
