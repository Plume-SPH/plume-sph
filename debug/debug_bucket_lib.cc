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



///*find particle by the key
// * keyin is input particle key
// * keycheck is the key that is given and all keyin will be compared with keycheck
// * pi is the pointer points to the particle corresponding to keyin
// */
//Particle* find_particle (unsigned* keyin, unsigned* keycheck, Particle* pi)
//{
//	int i;
//	bool find = true;
//	for (i=0; i<TKEYLENGTH; i++)
//		if (keyin[i] != keycheck[i])
//		{
//			find = false;
//			return NULL;
//		}
//
//    return find;
//}


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
		  	for (i = 0; i < KEYLENGTH; i++)
		  		keytemp[i] = b_curr->getKey().key[i];

		  	if (find_bucket (keytemp, keycheck))
		  		cout << "The bucket found!" << endl;
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
