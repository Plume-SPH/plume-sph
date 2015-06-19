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


//int
//put_ghost_particles (THashTable * P_table, HashTable * BG_mesh,TimeProps *timeprops,
//                     vector<BucketHead> & partition_table, MatProps * matprops,
//                     int * my_comm, int myid, int numproc, int * added_ghosts)
//{
//    int i, j;
//    int Up[3]   = {0, 0, 2}; // pos-z search direction
//    int Down[3] = {0, 0, 1}; // neg-z search direction
//    unsigned btkey[KEYLENGTH];
//    Bucket *buck0 = NULL, *buck1=NULL, *buck2=NULL;
//
//#ifdef DEBUG
//    char filename[20];
//    sprintf (filename, "put_ghost%03d.log\0", myid);
//    FILE * fp = fopen (filename, "a+");
//#endif
//
//    /*****************************************
//     * put ghost particle on current process
//     ****************************************/
//    // initiall set added_ghosts to 0
//    *added_ghosts = 0;
//
//    // scan all the boundary buckets to see if it, or any of
//    // its neighbors have real particles, if yes, add ghost particles
//    vector<BucketHead>::iterator b_itr;
//    for (b_itr = partition_table.begin ();
//         b_itr != partition_table.end (); b_itr++)
//    {
//        // get first two bucket from the bottom
//        for (i = 0; i < KEYLENGTH; i++)
//        	btkey[i] = *(b_itr->get_buck_head () + i);
//
//        // declare important bucket pointers
//        buck0 = (Bucket *) BG_mesh->lookup (btkey);
//        buck1 = (Bucket *) BG_mesh->lookup (buck0->which_neigh (Up));
//        buck2 = NULL;
//        bool put_ghosts = false;
//
//        // if currect column has ghost particles, go to next one
//        if ( buck1->has_ghost_particles () )
//            continue;
//
//        // go through all the boundary buckets in the column
//        do
//        {
//            const int * neigh_proc = buck1->get_neigh_proc ();
//            Key * neighbors = buck1->get_neighbors ();
//            for (i = 0; i < NEIGH_SIZE; i++)
//                if ( neigh_proc[i] > -1 )
//                {
//                    Bucket * neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);
//                    // the bucket may or may not be present if it doesn't
//                    // belong to current process
//                    if ( (neigh) && (neigh->has_real_particles ()))
//                        put_ghosts = true;
//
//                    // check for mesh errors
//                    if ((! neigh) && (neigh_proc[i] == myid))
//                    {
//                        fprintf (stderr,"Error: neigh-bucket missing on %d\n",
//                                 myid);
//                        fprintf (stderr,"trace: %s:%d\n", __FILE__, __LINE__);
//                        return -1;
//                    }
//                }
//            buck2 = buck1;
//            buck1 = (Bucket *) BG_mesh->lookup (buck1->which_neigh (Up));
//        }
//        while (buck1->get_bucket_type () == MIXED);
//
//        if ( put_ghosts )
//        {
////            buck2->put_ghost_particles (P_table, BG_mesh, matprops, timeprops);
//        	 buck2->put_ghost_particles ();
//            *added_ghosts = 1;
//
//            // this should be already done at this point
//            // just to be safe
//            buck0->mark_active ();
//            buck0->set_wall_ghost_particles (true);
//        }
//    }
//
//    /*****************************************
//     * Now send info about real particles,
//     * for the ghosts on partition boundaries
//     ****************************************/
//    int * send_info = new int [numproc];
//    int * recv_info = new int [numproc];
//
//    for (i = 0; i < numproc; i++)
//    {
//        send_info[i] = 0;
//        recv_info[i] = 0;
//    }
//
//    // make calls to get recv_info from other procs
//    int tag1 = 122312;
//    MPI_Request * request = new int [2 * numproc];
//    for (i = 0; i < numproc; i++)
//        if ( my_comm[i] )
//            MPI_Irecv ((recv_info + i), 1, MPI_INT, i, tag1,
//                       MPI_COMM_WORLD, (request + i));
//
//   // calculate send_info
//    HTIterator * itr = new HTIterator (BG_mesh);
//    while ((buck0 = (Bucket *) itr->next ()))
//        if ((! buck0->is_guest ()) &&
//            (buck0->has_real_particles ()) &&
//            (buck0->get_bucket_type () == MIXED))
//        {
//            const int * neigh_proc = buck0->get_neigh_proc ();
//            for (i = 0; i < NEIGH_SIZE; i++)
//                if ((neigh_proc[i] > -1) &&
//                    (neigh_proc[i] != myid))
//                    send_info[neigh_proc[i]] += KEYLENGTH;
//        }
//
//    // make send calls
//    for (i = 0; i < numproc; i++)
//        if ( my_comm[i] )
//            MPI_Isend ((send_info +i), 1, MPI_INT, i, tag1,
//                       MPI_COMM_WORLD, (request + numproc + i));
//
//    // allocate memory to sendbuf
//    unsigned ** sendbuf = new unsigned * [numproc];
//    for (i = 0; i < numproc; i++)
//        if ( send_info[i] > 0 )
//            sendbuf[i] = new unsigned [send_info[i]];
//
//    // check if receives are finished
//    MPI_Status recv_st;
//    for (i = 0; i < numproc; i++)
//        if ( my_comm[i] )
//            MPI_Wait ((request + i), & recv_st);
//
//    // allocate memory to recvbuf
//    unsigned ** recvbuf = new unsigned * [numproc];
//    for (i = 0; i < numproc; i++)
//        if ( recv_info[i] > 0 )
//            recvbuf[i] = new unsigned [recv_info[i]];
//
//    // need to re-use MPI_Request array, check here if sends
//    // are done, or wait. Theoretically they should be done.
//    MPI_Status send_st;
//    for (i = 0; i < numproc; i++)
//        if ( my_comm[i] )
//            MPI_Wait ((request + numproc + i), & send_st);
//
//    // post receives
//    recv_info[myid] = 0;
//    int tag = 812012;
//    for (i = 0; i < numproc; i++)
//        if ( recv_info[i] > 0 )
//            MPI_Irecv (recvbuf[i], recv_info[i], MPI_UNSIGNED, i, tag,
//                       MPI_COMM_WORLD, (request + i));
//
//    // pack data for sending
//    int * count = new int [numproc];
//    for (i = 0; i <numproc; i++)
//        count[i] = 0;
//
//    itr->reset ();
//    while ((buck0 = (Bucket *) itr->next ()))
//        if ((! buck0->is_guest ()) &&
//            (buck0->has_real_particles ()) &&
//            (buck0->get_bucket_type () == MIXED))
//        {
//            const int * neigh_proc = buck0->get_neigh_proc ();
//            Key * neighbors = buck0->get_neighbors ();
//            for (i = 0; i < NEIGH_SIZE; i++)
//                if ((neigh_proc[i] > -1) &&
//                    (neigh_proc[i] != myid))
//                {
//                    int k = neigh_proc[i];
//                    for (j = 0; j < KEYLENGTH; j++)
//                        sendbuf[k][count[k]++] = neighbors[i].key[j];
//                }
//        }
//
//    // post sends
//    send_info[myid] = 0;
//    for (i = 0; i < numproc; i++)
//        if ( send_info[i] > 0 )
//            MPI_Isend (sendbuf[i], send_info[i], MPI_UNSIGNED, i, tag,
//                       MPI_COMM_WORLD, (request + numproc + i));
//
//    // wait for receive and add ghosts if needed
//    MPI_Status status;
//    for (i = 0; i < numproc; i++)
//        if ( recv_info[i] > 0 )
//        {
//            MPI_Wait ((request + i), &status);
//            for (j = 0; j < recv_info[i]; j += KEYLENGTH)
//            {
//                for (int k = 0; k < KEYLENGTH; k++)
//                	btkey[k] = recvbuf[i][j + k];
//                buck0 = (Bucket *) BG_mesh->lookup (btkey);
//                assert (buck0);
//                // if bucket has ghosts; skip
//                if ((buck0->get_bucket_type () == MIXED) &&
//                    (! buck0->has_ghost_particles ()))
//                {
//                    // search the top most boundary bucket
//                    buck1 = buck0;
//                    while ( buck1->get_bucket_type () == MIXED )
//                    {
//                        buck2 = buck1;
//                        buck1 = (Bucket *)
//                            BG_mesh->lookup (buck1->which_neigh (Up));
//                    }
////                    buck2->put_ghost_particles (P_table, BG_mesh, matprops, timeprops);
//                    buck2->put_ghost_particles ();
//                    *added_ghosts = 1;
//                }
//            }
//        }
//
//    // wait for all the sends to finish
//    for (i = 0; i < numproc; i++)
//        if (send_info[i] > 0)
//            MPI_Wait ((request + numproc + i), & status);
//
//    // clean up
//    for (i = 0; i < numproc; i++)
//    {
//        if (send_info[i] > 0)
//            delete [] sendbuf[i];
//        if (recv_info[i] > 0)
//            delete [] recvbuf[i];
//    }
//
//    delete [] sendbuf;
//    delete [] recvbuf;
//    delete [] send_info;
//    delete [] recv_info;
//    delete [] request;
//    delete [] count;
//    delete itr;
//
//    return 0;
//}

/*
 * Scan mesh and mark buckets active and inactive
 */
int
update_bgmesh (HashTable * BG_mesh, int myid, int numproc, int * my_comm)
{
    int i, j;
    unsigned btkey[KEYLENGTH];

    // visit every bucket
    HTIterator * itr = new HTIterator (BG_mesh);
    Bucket *curr_bucket = NULL, *neigh = NULL;
    while ((curr_bucket = (Bucket *) itr->next ()))
    {
        // mark it inactive to start with
        curr_bucket->mark_inactive ();

        // if any neighbor as any particle, mark current bucket active
        const int * neigh_proc = curr_bucket->get_neigh_proc ();
        Key * neighbors = curr_bucket->get_neighbors ();

        for (i = 0; i < NEIGH_SIZE; i++)
            if (neigh_proc[i] > -1)
            {
                // some neighs may not of available on current process
                neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);
                if ( neigh && (neigh->get_plist ().size () > 0))
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
    }

    /* To switch buckets on foreign procs,
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
    while ((curr_bucket = (Bucket *) itr->next ()))
        if ((! curr_bucket->is_guest ()) &&
            (! curr_bucket->get_plist ().empty ()))
        {
            const int * neigh_proc = curr_bucket->get_neigh_proc ();
            for (i = 0; i < NEIGH_SIZE; i++)
                if ((neigh_proc[i] > -1) &&
                    (neigh_proc[i] != myid))
                    send_info[neigh_proc[i]] += KEYLENGTH;
        }

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
    while ((curr_bucket = (Bucket *) itr->next ()))
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
                curr_bucket = (Bucket *) BG_mesh->lookup (btkey);
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
