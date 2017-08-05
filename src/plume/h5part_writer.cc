/*
 * h5part_writer.cc
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <hdf5.h>

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <cstdio>
#include <vector>
using namespace std;

#include <thashtab.h>
#include "constant.h"
#include "particler.h"
#include "sph_header.h"
#include "hdf5calls.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

#if CODE_DIMENSION==1
#include <iostream>
#include <fstream>
#include <string>
#endif

//#define WRITE_GHOSTS

void
write_h5part(int myid, int numproc, THashTable * P_table, TimeProps * timepros)
{
  int i, j;
  vector < double >x, y, z, Vx, Vy, Vz, rho, engr, mssfrc, prss;
  vector < int >  phase, bctp, guest;


#ifdef DEBUG
  vector < int > involved;
#endif

#ifdef WRITE_PID
  vector < int > myprocess;
#endif

#ifdef WRITE_PMASS
  vector < double > mymass;
#endif

#if defined (WRITE_SML) && (ADAPTIVE_SML >=1)
  vector < double > mysml;
#endif

  char filename[18];
  static int step = 0;
  hid_t fp;
  herr_t ierr;

#ifdef PARALLEL_IO
  sprintf(filename, "pvplot_out.h5part");
#else
  sprintf(filename, "pvplot%03d.h5part", myid);
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
    if ((!pi->is_guest ()) && (pi->get_bc_type ()==100) )//non guest and is real
    {
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
      prss.push_back(pi->get_pressure());
      phase.push_back(pi->which_phase());
      bctp.push_back(pi->get_bc_type ());
      guest.push_back(pi->get_guest());

#ifdef DEBUG
      involved.push_back(pi->get_involved ());
#endif

#ifdef WRITE_PID
      myprocess.push_back(myid);
#endif

#ifdef WRITE_PMASS
      mymass.push_back(pi->get_mass ());
#endif

#if defined (WRITE_SML) && (ADAPTIVE_SML >=1)
      mysml.push_back(pi->get_smlen ());
#endif

#ifndef WRITE_GHOSTS
    }
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
  for (i = id_lims[myid]; i < id_lims[myid + 1]; i++) //for ID
  {
    ibuf[j] = i;
    j++;
  }

  // particle ids
  ierr = GH5_Write(gid, "ID", dims, (void *) ibuf, start, my_count, INTTYPE);

  //involved
#ifdef DEBUG
  copy(involved.begin(), involved.end(), ibuf);
  ierr = GH5_Write(gid, "Involved", dims, (void *) ibuf, start, my_count, INTTYPE);
#endif

  //myprocess
#ifdef WRITE_PID
  copy(myprocess.begin(), myprocess.end(), ibuf);
  ierr = GH5_Write(gid, "myprocess", dims, (void *) ibuf, start, my_count, INTTYPE);
#endif

  //mymass
#ifdef WRITE_PMASS
  copy(mymass.begin(), mymass.end(), buf);
  ierr = GH5_Write(gid, "mymass", dims, (void *) buf, start, my_count,  DOUBLETYPE);
#endif

  //mysml
#if defined (WRITE_SML) && (ADAPTIVE_SML >=1)
  copy(mysml.begin(), mysml.end(), buf);
  ierr = GH5_Write(gid, "sml", dims, (void *) buf, start, my_count,  DOUBLETYPE);
#endif

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

  //pressure
  copy(prss.begin(), prss.end(), buf);
  ierr = GH5_Write(gid, "prss", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // phase
  copy(phase.begin(), phase.end(), ibuf);
  ierr = GH5_Write(gid, "phase", dims, (void *) ibuf, start, my_count,INTTYPE);

  // bctp
  copy(bctp.begin(), bctp.end(), ibuf);
  ierr = GH5_Write(gid, "bctp", dims, (void *) ibuf, start, my_count, INTTYPE);

  // guest flag
  copy(guest.begin(), guest.end(), ibuf);
  ierr = GH5_Write(gid, "guest", dims, (void *) ibuf, start, my_count, INTTYPE);

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


void
write_h5part_show(int myid, int numproc, THashTable * P_table, TimeProps * timepros)
{
  int i, j;
  vector < double >x, y, z, Vx, Vy, Vz, rho, engr, mssfrc, prss;
  vector < int >  phase, bctp, guest;


  char filename[18];
  static int step = 0;
  hid_t fp;
  herr_t ierr;

#ifdef PARALLEL_IO
  sprintf(filename, "pvplot_show_out.h5part");
#else
  sprintf(filename, "pvplot_show%03d.h5part", myid);
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
    if (pi->need_neigh())//non guest and is real
    {
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
      prss.push_back(pi->get_pressure());
      phase.push_back(pi->which_phase());
      bctp.push_back(pi->get_bc_type ());
      guest.push_back(pi->get_guest());
    }
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
  for (i = id_lims[myid]; i < id_lims[myid + 1]; i++) //for ID
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

  //pressure
  copy(prss.begin(), prss.end(), buf);
  ierr = GH5_Write(gid, "prss", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // phase
  copy(phase.begin(), phase.end(), ibuf);
  ierr = GH5_Write(gid, "phase", dims, (void *) ibuf, start, my_count,INTTYPE);

  // bctp
  copy(bctp.begin(), bctp.end(), ibuf);
  ierr = GH5_Write(gid, "bctp", dims, (void *) ibuf, start, my_count, INTTYPE);

  // guest flag
  copy(guest.begin(), guest.end(), ibuf);
  ierr = GH5_Write(gid, "guest", dims, (void *) ibuf, start, my_count, INTTYPE);

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

#if CODE_DIMENSION==1
//function that used to write output into a cvs file ---> useful for small size output, such as output for shock tube problem.
void
write_csv(int myid, int numproc, THashTable * P_table, TimeProps * timepros)
{

	  char filename[18];
	  static int step = 0;

	  ofstream outputFile;

//	  if (timepros->ifstart())
		  outputFile.open("pvplot_out.csv", std::fstream::out | std::fstream::app);
//	  else
//		  outputFile.open("pvplot_out.csv", std::fstream::app);

	  // write the file headers
	  //if (timepros->ifstart())
		  outputFile << "RHO" << ","<<"X"<< "," << "U"<< ","<<"E"<< "," <<"KS"<< "," << "P"<< "," << "PHASE"<< "," << "BCTP"<< ","<< "GUEST_FLAG"<< ","<< "INVOLVE_FLAG" << "smlen" << std::endl;

	  THTIterator *itr = new THTIterator(P_table);
	  Particle *pi = NULL;

	  while ((pi = (Particle *) itr->next()))
	  {
#if WRITE_GHOSTS!=2
	    if (pi->get_bc_type()==100) //only output real
	    {
#endif
	        outputFile << pi->get_density() << "," << *(pi->get_coords()) << ","<< *(pi->get_vel())<< ","<<pi->get_energy () << ","<<pi->get_mass_frac() << ","<< pi->get_pressure() << "," << pi->which_phase()<< "," <<pi->get_bc_type ()<< ","<< pi->get_guest()<< ","<< pi->get_involved () << ","<< pi->get_involved ()<< std::endl;
#if WRITE_GHOSTS!=2
	    }
#endif
	  }

	  // close the output file
	  outputFile.close();
}
#endif  //CODE_DIMENSION==1
