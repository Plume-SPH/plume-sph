/*
 * matlab_writer.cc
 *
 *  Created on: Apr 30, 2015
 *      Author: zhixuanc
 */

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstdio>
#include <vector>
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <constant.h>
#include <particler.h>
#include <buckhead.h>
#include <bucket.h>
#include <outforms.h>
#include <hdf5calls.h>


//Function for writing buckets decomposition in a hdf5 file
void
write_matlab(int myid, THashTable * P_table, HashTable * BG_mesh,
             TimeProps * timeprops, vector <BucketHead>  & partition_table)
{
    int i, j;
    char file1[25];
    static int icount = 0;
    double mincrd[DIMENSION], maxcrd[DIMENSION];
    Bucket * buck2 = NULL;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

    sprintf(file1, "buckets%03d%06d.h5", myid, timeprops->step);
    icount++;

    hid_t fp = GH5_fopen_serial (file1, 'w');

    int size = (int) partition_table.size ();
    double * xcoord = new double [size * 8];
    double * ycoord = new double [size * 8];
    double * zcoord = new double [size * 8];
    int * active = new int [size];

    vector <BucketHead>::iterator itr;
    j = 0;
    for (itr = partition_table.begin(); itr != partition_table.end(); itr++)
    {
    	tempptr = BG_mesh->lookup ( (unsigned *)itr->get_buck_head ());
    	breif_buck = (BriefBucket *) tempptr;
    	if (breif_buck->check_brief())
    		continue;
    	else
    	{
            Bucket * buck = (Bucket *) tempptr;
            assert (buck);

            active[j] = 0;
            buck2 = buck;

            /*I am not sure about the reason he did this, but in my opinion, it is not necessary.
             * what's worse, his way of search for the overground bucket is not suitable for me,
             * in my code, the situation is more complicated and I need to use the bucket index
             * to determine the direction of searching!
            */
    //        while (buck2->get_bucket_type () != OVERGROUND)
    //        {
    //            buck2 = (Bucket *) BG_mesh->lookup (buck2->which_neigh (Up));
    //
    //            if ( buck2->has_real_particles ())
    //            {
    //                active[j] += 4;
    //                break;
    //            }
    //        }

            for (i = 0; i < DIMENSION; i++)
            {
                mincrd[i] = *(buck->get_mincrd()+i);
                maxcrd[i] = *(buck->get_maxcrd()+i);
            }

            /* corner 0*/
            xcoord[8*j] = mincrd[0];
            ycoord[8*j] = mincrd[1];
            zcoord[8*j] = mincrd[2];

            /* corner 1*/
            xcoord[8*j +1] = maxcrd[0];
            ycoord[8*j +1] = mincrd[1];
            zcoord[8*j +1] = mincrd[2];

            /* corner 2*/
            xcoord[8*j + 2] = maxcrd[0];
            ycoord[8*j + 2] = maxcrd[1];
            zcoord[8*j + 2] = mincrd[2];

            /* corner 3*/
            xcoord[8*j + 3] = mincrd[0];
            ycoord[8*j + 3] = maxcrd[1];
            zcoord[8*j + 3] = mincrd[2];

            /* corner 4*/
            xcoord[8*j + 4] = mincrd[0];
            ycoord[8*j + 4] = mincrd[1];
            zcoord[8*j + 4] = maxcrd[2];

            /* corner 5*/
            xcoord[8*j +5] = maxcrd[0];
            ycoord[8*j +5] = mincrd[1];
            zcoord[8*j +5] = maxcrd[2];

            /* corner 6*/
            xcoord[8*j + 6] = maxcrd[0];
            ycoord[8*j + 6] = maxcrd[1];
            zcoord[8*j + 6] = maxcrd[2];

            /* corner 7*/
            xcoord[8*j + 7] = mincrd[0];
            ycoord[8*j + 7] = maxcrd[1];
            zcoord[8*j + 7] = maxcrd[2];

            if ( buck->is_active () )
                active[j] += 1;

            // add 1 for paritcles
            if ( buck->has_wall_ghost_particles ()) //not sure whether this is correct or not. Not sure about the meaning of active.
                active[j]++;

            // increment the counter
            j++;
    	}//end of if bucket is not brief bucket
    }//end of for loop go through partition table

    int dims1[2] = {size, 8};
    int dims2[2] = {size, 1};
    GH5_WriteS (fp, "/xcoord", dims1, (void *) xcoord, 0, 0, DOUBLETYPE);
    GH5_WriteS (fp, "/ycoord", dims1, (void *) ycoord, 0, 0, DOUBLETYPE);
    GH5_WriteS (fp, "/zcoord", dims1, (void *) zcoord, 0, 0, DOUBLETYPE);
    GH5_WriteS (fp, "/active", dims2, (void *) active, 0, 0, INTTYPE);
    GH5_fclose (fp);

    delete [] xcoord;
    delete [] ycoord;
    delete [] zcoord;
    delete [] active;
    return;
}
