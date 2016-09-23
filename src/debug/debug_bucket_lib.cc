/*
 * debug_bucket_lib.cc
 *
 *  Created on: Jul 27, 2015
 *      Author: zhixuanc
 */

/*
 * debug_lib.cc
 *
 *  Created on: Jun 16, 2015
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
#include <outforms.h>

#include "sph_header.h"
#include "constant.h"
#include "parameters.h"

#include <hdf5.h>
#include "hdf5calls.h"


/*find bucket by the key
 * keyin is input particle key
 * keycheck is the key that is given and all keyin will be compared with keycheck
  */
bool find_bucket (unsigned* keyin, unsigned* keycheck)
{
	int i;
	for (i=0; i<KEYLENGTH; i++)
		if (keyin[i] != keycheck[i])
			return false;

    return true;
}

/*find bucket by the pos
 *
 */
bool find_bucket (double* in, double* check)
{
	int i;
	for (i=0; i<DIMENSION; i++)
		if (in[i] != check[i])
			return false;

    return true;
}

//function to find particle with given key, will be useful in debugging.
void check_bucket_bykey (HashTable * BG_mesh)
{

    bool do_search = true;
    bool find;
    unsigned keycheck[KEYLENGTH] = {279233414, 853448706} ; //the bucket contain the missing particle
    //The above bucket is the neighbor of  {202798183, 343262493}
    unsigned keytemp[KEYLENGTH] ;

    int i;
	HTIterator *itr = new HTIterator(BG_mesh);
	Bucket *b_curr = NULL;

	while ((b_curr = (Bucket *) itr->next()))
	{
		if (do_search)
		{
		  	for (i = 0; i < KEYLENGTH; i++)
		  		keytemp[i] = b_curr->getKey().key[i];

		  	if (find_bucket (keytemp, keycheck))
		  	{
		  		cout << "The bucket found!" << endl;
		  		cout << "particles type in this bucket is :" << b_curr->get_particles_type() <<endl;
		  		cout << "has involved of this bucket is :" << b_curr->get_has_involved() <<endl;
		  	}
		 }
	 }

}

//go though all buckets and check if any of them are guest!

//function to find particle with given key, will be useful in debugging.
void check_bucket_guest (HashTable * BG_mesh)
{

    bool do_search = true;
    bool find;
    unsigned keycheck[KEYLENGTH] = {125131437, 1454069162};
    //The above bucket is the neighbor of  {202798183, 343262493}
    unsigned keytemp[KEYLENGTH] ;

    int i;
	HTIterator *itr = new HTIterator(BG_mesh);
	Bucket *b_curr = NULL;

	while ((b_curr = (Bucket *) itr->next()))
	{
		if (do_search)
		{
		  	if (b_curr->is_guest())
		  		cout << "The guest bucket found!" << endl;
		 }
	 }

}

//function that used to go through all buckets to check whether some buckets in the BG_mesh table is wired or not.
void BG_mesh_check (HashTable * BG_mesh)
{
	HTIterator *itr = new HTIterator(BG_mesh);
	Bucket *b_curr = NULL;

	while ((b_curr = (Bucket *) itr->next()))
	{
		if ((b_curr->getKey()).key[0] == 0)
		{
		  	cout << "One wired bucket found!" << endl;
		}
	}

}

//function that used to check mesh error
void BG_mesh_err_check(HashTable * BG_mesh, int myid)
{
    int i;

    // visit every bucket
    HTIterator * itr = new HTIterator (BG_mesh);
    Bucket *curr_bucket = NULL, *neigh = NULL;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

    while (tempptr=itr->next ())
    {
    	breif_buck = (BriefBucket *) tempptr;
    	if(breif_buck->check_brief()) //Currently, I do not care about the mesh error in brief bucket. But actually I should also check that for brief bucket
    		continue;
    	else
    	{
    		curr_bucket = (Bucket*) tempptr;

            // if any neighbor has any particle, mark current bucket active
            const int * neigh_proc = curr_bucket->get_neigh_proc ();
            Key * neighbors = curr_bucket->get_neighbors ();

            for (i = 0; i < NEIGH_SIZE; i++)
                if (neigh_proc[i] > -1)
                {
                    // some neighs may not of available on current process

                    neigh = (Bucket *) BG_mesh->lookup (neighbors[i]); //The neighbor bucket might be brief bucket. but we do not care about that, we can set a pointer for bucket point to brief bucket
                    // check for mesh errors
                    if ((! neigh) && (neigh_proc[i] == myid))
                    {
                        fprintf (stderr,"Error: neigh-bucket missing on %d\n",
                                 myid);
                        fprintf (stderr,"trace: %s:%d\n", __FILE__, __LINE__);
                        return ;
                    }
                }

    	} //end of if bucket is not brief bucket
    }//end of while loop go through all buckets
}
