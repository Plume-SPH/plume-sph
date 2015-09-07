/*
 * scan_outside_layer.cc
 *
 *  Created on: Aug 8, 2015
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

#include <hashtab.h>
#include <thashtab.h>
#include <hilbert.h>
#include <bucket.h>
//#include <particler.h>
#include <multiproc.h>
#include <properties.h>

#include "constant.h"
#include "parameters.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

//scan the most outside layer of buckets satisfy has_involved>0,
//if they have any involved particles (involved = 2), then make the bucket to be has_involved.
//What will happen next is the domain will be adjust such that at least one layer of potential involved bucket is added outside the domain.

/*
 * Modification made on Sep 7th
 * Originally, only scan the most outside layer.
 * After modification, all has_potential_involved bucket will be scanned.
 * The only difference is that only the most outside layer will trigger the adapt++
 */
int scan_outside_layer (THashTable * P_table, HashTable * BG_mesh, int numproc, int myid)
{
	int adapt = 0;
	int i;
	int REAL = 100;
	HTIterator * itr = new HTIterator (BG_mesh);
	Bucket * Bnd_buck = NULL;
	Bucket * Curr_buck = NULL;
	Bucket * neigh = NULL;

	vector < TKey > plist;
	vector < TKey >::iterator p_itr;

//	Particle * p_curr;

	while ((Bnd_buck = (Bucket *) itr->next ()))
		/*
		 * 1) has_involved = 1
		 * 2) is not guest --> not necessary to do extra repeated work ---> communication will be executed after scan of the most outer layer
		 * --> more important, we need neighbor info to the following work,
		 * --> But some neighbors of guest buckets might be on other procs
		 */
//    if (Bnd_buck->is_has_potential_involved () && !Bnd_buck->is_guest() ) //for Bnd_buck->get_bucket_type () != UNDERGROUND, has_involved is always zero --> the old code, modified, because it is not necessary to check buckets that already has involved particles
	if (Bnd_buck->get_has_involved() ==1  && !Bnd_buck->is_guest() ) //for Bnd_buck->get_bucket_type () != UNDERGROUND, has_involved is always zero
	{
    	plist = Bnd_buck->get_plist();
    	for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
    	{
    	    Particle *p_curr = (Particle *) P_table->lookup(*p_itr);
    	    assert(p_curr);
    	    if (p_curr->is_involved())
    	    {
    	    	Bnd_buck->set_has_involved (HAS_BOTH);


    	    	//if the bucket is the most outside layer, will trigger adapt++
    	    	//check whether bucket is the most outside layer or not
    	    	// if any neighbor is pressure ghost
    	    	const int * neigh_proc = Bnd_buck->get_neigh_proc ();
    	    	Key * neighbors = Bnd_buck->get_neighbors ();
    	    	for (i = 0; i < NEIGH_SIZE; i++)
    	    	if (neigh_proc[i] > -1)
    	    	{
    	    	    neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);
    	    	    /*
    	    	    * 1) is not empty
    	    	    * 2) has_involved = 0 --> up to now, The only possible are pressure ghost buckets and underground buckets
    	    	    * 3) so, need to make sure that bucket is not UNDERGROUND bucket
    	    	    */
    	    	    if (((neigh->get_plist ()).size()) && !neigh->get_has_involved() && neigh->get_bucket_type () != UNDERGROUND)
    	    	    {
    	    	    	adapt++;
    	    	        break; //jump out from the for loop of going through all neighbor --> have already confirmed that Bnd_buck is the most out layer of the
    	    	    }//end of if one neighbor is

    	    	}// end of go through all neighbor of bucket Bnd_buck

    	    	goto next_bucket; // if found involved particle in current bucket --> then mark the current bucket as has_involved (actually, has both involved and potential involved) and go and check next bucket
    	    }
    	}// loop: go through all particles in the bucket

    	next_bucket:
    	1;
	}//end of go through all buckets

	return adapt;

}
