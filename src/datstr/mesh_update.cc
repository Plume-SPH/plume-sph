/*
 * mesh_update.cc
 *
 *  Created on: Mar 12, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <vector>
using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <properties.h>
#include <bucket.h>
#include <particle.h>
#include <buckhead.h>

#ifdef MULTI_PROC
#include <mpi.h>
#endif

/*
 * Scan mesh and mark buckets active and inactive
 */
//It would be better if the inactive bucket will be turn to brief bucket. ---> For more complicated application, currently the domain is expanding, so this should not be a problem.
//This will help to save memory!
int
update_bgmesh (HashTable * BG_mesh, int myid, int numproc, int * my_comm)
{
    int i, j;
    unsigned btkey[KEYLENGTH];

    // visit every bucket
    HTIterator * itr = new HTIterator (BG_mesh);
    Bucket *curr_bucket = NULL, *neigh = NULL;
	BriefBucket *breif_buck = NULL;
	void * tempptr =NULL;

    while (tempptr=itr->next ())
    {
    	breif_buck = (BriefBucket *) tempptr;
    	if(breif_buck->check_brief()) //Brief buckets are default inactive
    		continue;
    	else
    	{
    		curr_bucket = (Bucket*) tempptr;

            // mark it inactive to start with
            curr_bucket->mark_inactive ();

            // if any neighbor has any particle, mark current bucket active
            const int * neigh_proc = curr_bucket->get_neigh_proc ();
            Key * neighbors = curr_bucket->get_neighbors ();

            for (i = 0; i < NEIGH_SIZE; i++)
                if (neigh_proc[i] > -1)
                {
                    // some neighs may not of available on current process

                    neigh = (Bucket *) BG_mesh->lookup (neighbors[i]); //The neighbor bucket might be brief bucket. but we do not care about that, we can set a pointer for bucket point to brief bucket
                    if ( neigh && (neigh->get_plist ().size () > 0) && !neigh->check_brief()) //It is not necessary to check whether neighbor bucket is brief bucket or not, if it is brief bucket, it will return a empty particle list
                    {
                        curr_bucket->mark_active ();
                        break; // no need to look any further
                    }

                    // check for mesh errors
                    if ((! neigh) && (neigh_proc[i] == myid))
                    {
                        fprintf (stderr,"Error: neigh-bucket missing on %d\n",
                                 myid);
                        fprintf (stderr,"trace: %s:%d\n", __FILE__, __LINE__);
                        return -1;
                    }
                }

    	} //end of if bucket is not brief bucket
    }//end of while loop go through all buckets

    /*
     * To switch buckets on foreign procs ---> because for guest buckets, it is very possible that its non-empty neighbor is not available locally
     * do some communications
     */
    int * recv_info = new int [numproc];
    int * send_info = new int [numproc];
    for (i = 0; i < numproc; i++)
    {
        send_info[i] = 0;
        recv_info[i] = 0;
    }

    // make calls to get recv_info from other procs
    int tag1 = 122311;
    MPI_Request * request = new int [2 * numproc];
    for (i = 0; i < numproc; i++)
        if ( my_comm[i] )
            MPI_Irecv ((recv_info + i), 1, MPI_INT, i, tag1,
                       MPI_COMM_WORLD, (request + i));

    // calculate send_info
    itr->reset ();
    while (tempptr=itr->next ())
    {
		breif_buck = (BriefBucket *) tempptr;
		if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
			continue;
		else
		{
			curr_bucket = (Bucket*) tempptr;
	        if ((! curr_bucket->is_guest ()) &&
	            (! curr_bucket->get_plist ().empty ())) //bucket is not empty, this is important: In the first loop go through all buckets, active inactive is marked based info only available locally.
	        	                                        //This part will exchange information which are not available locally, to mark active and inactive.
	        {
	            const int * neigh_proc = curr_bucket->get_neigh_proc ();
	            for (i = 0; i < NEIGH_SIZE; i++)
	                if ((neigh_proc[i] > -1) &&
	                    (neigh_proc[i] != myid))
	                    send_info[neigh_proc[i]] += KEYLENGTH;
	        }
		}
    }//end of while loop

    // make send calls
    for (i = 0; i < numproc; i++)
        if ( my_comm[i] )
            MPI_Isend ((send_info + i), 1, MPI_INT, i, tag1,
                       MPI_COMM_WORLD, (request + numproc + i));

    // allocate memory to sendbuf
    unsigned ** sendbuf = new unsigned * [numproc];
    for (i = 0; i < numproc; i++)
        if ( send_info[i] > 0 )
            sendbuf[i] = new unsigned [send_info[i]];

    // check if receives are finished
    MPI_Status status;
    for (i = 0; i < numproc; i++)
        if ( my_comm[i] )
            MPI_Wait ((request + i), & status);

    // allocate memory to recvbuf
    unsigned ** recvbuf = new unsigned * [numproc];
    for (i = 0; i < numproc; i++)
        if ( recv_info[i] > 0 )
            recvbuf[i] = new unsigned [recv_info[i]];

    // need to re-use MPI_Request array, check here if sends
    // are done, or wait. Theoretically they should be done.
    for (i = 0; i < numproc; i++)
        if ( my_comm[i] )
            MPI_Wait ((request + numproc + i), & status);

    // post receives
    recv_info[myid] = 0;
    int tag = 812012;
    for (i = 0; i < numproc; i++)
        if ( recv_info[i] > 0 )
            MPI_Irecv (recvbuf[i], recv_info[i], MPI_UNSIGNED, i, tag,
                       MPI_COMM_WORLD, (request + i));

    // pack data for sending
    int * count = new int [numproc];
    for (i = 0; i <numproc; i++)
        count[i] = 0;

    itr->reset ();
    while (tempptr=itr->next ())
    {
		breif_buck = (BriefBucket *) tempptr;
		if (breif_buck->check_brief()) //if is brief bucket, this bucket contains nothing!
			continue;
		else
		{
			curr_bucket = (Bucket*) tempptr;
	        if ((! curr_bucket->is_guest ()) &&
	            (! curr_bucket->get_plist ().empty ()))
	        {
	            const int * neigh_proc = curr_bucket->get_neigh_proc ();
	            Key * neighbors = curr_bucket->get_neighbors ();
	            for ( i = 0; i < NEIGH_SIZE; i++)
	                if ((neigh_proc[i] > -1) &&
	                    (neigh_proc[i] != myid))
	                {
	                    int k = neigh_proc[i];
	                    for (j = 0; j < KEYLENGTH; j++)
	                        sendbuf[k][count[k]++] = neighbors[i].key[j];
	                }
	        }
		}//end of if bucket is not brief buckets
    }//end of while loop go through all buckets


    // post sends
    send_info[myid] = 0;
    for (i = 0; i < numproc; i++)
        if ( send_info[i] > 0 )
            MPI_Isend (sendbuf[i], send_info[i], MPI_UNSIGNED, i, tag,
                       MPI_COMM_WORLD, (request + numproc + i));

    // wait for receive, and turn buckets ON
    for (i = 0; i < numproc; i++)
        if ( recv_info[i] > 0 )
        {
            MPI_Wait ((request + i), & status);
            for (j = 0; j < recv_info[i]; j += KEYLENGTH)
            {
                for (int k = 0; k < KEYLENGTH; k++)
                	btkey[k] = recvbuf[i][j + k];
                curr_bucket = (Bucket *) BG_mesh->lookup (btkey); //This btkey should always be the key of a bucket. it is impossible that the bucket is a brief bucket
                assert (curr_bucket);

                // regardless of its current status / make it active
                // overwriting is faster than checking
                curr_bucket->mark_active ();
            }
        }

    // wait for all the sends to finish
    for (i = 0; i < numproc; i++)
        if (send_info[i] > 0)
            MPI_Wait ((request + numproc + i), & status);

    // clean up
    for (i = 0; i < numproc; i++)
    {
        if (send_info[i] > 0)
            delete [] sendbuf[i];
        if (recv_info[i] > 0)
            delete [] recvbuf[i];
    }

    delete [] sendbuf;
    delete [] recvbuf;
    delete [] send_info;
    delete [] recv_info;
    delete [] request;
    delete [] count;
    delete itr;
    return 0;
}

//void
//delete_unused_ghosts (THashTable * P_table, HashTable * BG_mesh, int myid)
//{
//    int i;
//    int Up[3] = {0, 0, 2};
//    Bucket *buck = NULL, *neigh = NULL;
//    HTIterator * itr = new HTIterator (BG_mesh);
//    while ((buck = (Bucket *) itr->next ()))
//        if ((! buck->is_guest ()) &&
//            (buck->has_ghost_particles ()))
//        {
//
//            // some underground buckets have ghost-particle flag up,
//            // without having any particles,
//            // if such pretender is encountered leave it alone
//            if ((buck->get_bucket_type () == UNDERGROUND) &&
//                (buck->get_plist ().empty ()))
//                continue;
//
//            // set flag to true to begin with
//            bool delete_ghosts = true;
//
//            // check if the bucket is on partition-edge
//            // if yes, don't remove ghost yet.
//            const int * neigh_proc = buck->get_neigh_proc ();
//            Key * neighbors = buck->get_neighbors ();
//            for (i = 0; i < NEIGH_SIZE; i++)
//                if ((neigh_proc[i] != myid) &&
//                    (neigh_proc[i] > -1))
//                {
//                    delete_ghosts = false;
//                    break;
//                }
//                // if the bucket is not on partions-edge
//                //  check if any neighbor has real particles
//                else if ( neigh_proc[i] == myid )
//                {
//                    neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);
//                    do
//                    {
//                        if ( neigh->has_real_particles () )
//                        {
//                            delete_ghosts = false;
//                            break;
//                        }
//                        else
//                            neigh = (Bucket *)
//                                 BG_mesh->lookup (neigh->which_neigh (Up));
//                    }
//                    while (neigh->get_bucket_type () != OVERGROUND);
//                    if ( ! delete_ghosts )
//                        break;
//                }
//            // if both tests cleared, go ahead and remove the ghosts
//            if ( delete_ghosts )
//            {
//                vector<TKey> particles = buck->get_plist ();
//                vector<TKey>::iterator p_itr;
//                for (p_itr = particles.begin (); p_itr != particles.end();
//                     p_itr++);
//                {
//                    Particle * p_ghost = (Particle *) P_table->lookup (*p_itr);
//                    P_table->remove (*p_itr);
//                    delete p_ghost;
//                }
//            }
//        }
//
//    delete itr;
//    return;
//}
