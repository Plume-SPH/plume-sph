/*
 * update_out_layer.cc
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
#include <bnd_image.h>

#include "constant.h"
#include "parameters.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

/*
 * to improve efficiency, use a vector to store the buckets which are the most outside layer of the domain.
 * But I did not make decision yet, I have several options to handle this problem
 * 1) without using any thing sepcial
 * 2) using a new hash table to handle the active buckets
 * 3) use the old hash table, but add more indicators about bucket types like; most outside layer bucket ect.
 *
 */
void update_out_layer (THashTable * P_table, HashTable * BG_mesh, vector < OutSideBucket > outside_table, int numproc, int myid)
{

	HTIterator * itr = new HTIterator (BG_mesh);
	Bucket * Bnd_buck = NULL;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

	vector < TKey > plist;
	vector < TKey >::iterator p_itr;

	while (tempptr=itr->next ())
	{
		breif_buck = (BriefBucket *) tempptr;
		if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
			continue;
		else
		{
			Bnd_buck = (Bucket*) tempptr;
			/*
			 * 1) is not empty
			 * 2) is not guest -->not necessary to do extra work to check guest bucket
			 * 3) has_involved = 0 --> up to now, what only left is pressure ghost buckets and underground buckets
			 * 4) so, need to make sure that bucket is not UNDERGROUND bucket
			 *
			 */
			if (((Bnd_buck->get_plist ()).size()) && !Bnd_buck->is_guest() && !Bnd_buck->get_has_involved() && Bnd_buck->get_bucket_type () != UNDERGROUND)
			{

			}
		}
	}
	return;

}
