/*
 * =====================================================================================
 *
 *       Filename:  multiproc.h
 *
 *    Description:  
 *
 *        Created:  03/14/2011 08:55:06 PM EDT
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *        License:  GNU General Public License
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */

#ifndef MULTIPROC_H
#define MULTIPROC_H

#include <vector>
#include <list>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bnd_image.h>

#include "buckhead.h"

//! register MPI_structs 
void GMFG_new_MPI_Datatype();

void send_foreign_images(
    //! My process-id
    int,
    //! Total size of MPI pool
    int,
    //! Hash-Table of Buckets
    HashTable *,
    //! STL list of ghost reflections
    list < BndImage > &,
    //! flags for communication
    int *);

//! Update data across, inter-proc boundaries
void move_data(
    //! Number of total processes
    int,
    //! process's MPI rank 
    int,
    //! array of flags for communication with other procs
    int *,
    //! Hash-Table of SPH particles
    THashTable *,
    //! Hash-Table of Buckets 
    HashTable *);

//! Update ghost reflections from other processors
int move_bnd_images(
    //! My process ID
    int,
    //! Total number of process
    int,
    //! HashTable of particles
    THashTable *,
    //! HashTable of buckets
    HashTable *,
    //! STL Vector of Images
    list < BndImage > &);

//! repartion the domain if load-balance has changed
int repartition(
    //! STL vector of Partition Table Keys
    vector < BucketHead > &,
    //! Hash-Table of SPH particles
    THashTable *,
    //! Hash-Table of Buckets
    HashTable *,
    //! communication yes/ no flag
    int *);

//! delete guest buckets and particles
void delete_guest_buckets(
    //! Hash-Table of Buckets
    HashTable *,
    //! Hash-Table of particles
    THashTable *);

#endif // MULTIPROC__H
