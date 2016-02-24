
/*
 * =====================================================================================
 *
 *       Filename:  bgmesh.h
 *
 *    Description:  
 *
 *        Created:  08/31/2010 02:46:57 PM EDT
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

#ifndef BGMESH_H
#  define BGMESH_H

# include <vector>
# include <list>
using namespace std;

#include <hashtab.h>
#include <bnd_image.h>
#include <properties.h>
#include <buckhead.h>

int update_bgmesh (
    //! Bucket Hash-table 
    HashTable *, 
    //! Process ID 
    int ,
    //! Total number of procs
    int ,
    //! array of switches to communicate with other procs
    int *);

//void
//delete_unused_ghosts (
//		THashTable * , //P_table,
//		HashTable * ,  //BG_mesh,
//		int            //myid
//		);


void search_bnd_images (
    //! ProcessID
    int,
    //! Hash-table of particles
    THashTable *,
    //! Hash-table of buckets
    HashTable *,
    //! Vector of Boundary reflections
    list < BndImage > &,
    //! flag to reset image-table
    int);
#endif
