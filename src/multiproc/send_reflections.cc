/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author:      dkumar
 * Description: send over the reflections that belong to other prcs
 *
 *******************************************************************
 */

#include <list>
using namespace std;

#include <mpi.h>

#include <hashtab.h>
#include <bucket.h>
#include <bnd_image.h>
#include <exvar.h>


// send reflctions that belong to neighboring partitions
/*
 * It happens that the bucket which contains the images is a guest bucket, in this case, the image need to be send out to the process where the guest bucket "originally" belong to
 * At the time of "apply bc", these image will not be able to find their corresponding wall ghost particle and imposing of wall bc will be "delayed"
 *
 * */
void
send_foreign_images (int myid, int numprocs, HashTable * BG_mesh, 
                     list < BndImage > & Image_table, int * my_comm)
{

  if (numprocs < 2)
    return;

  int i, j;
//  int dir[DIMENSION];
//  double coord[DIMENSION], tmpc[DIMENSION];
  double intsct[DIMENSION];
//  unsigned ghst_key[KEYLENGTH], buck_key[KEYLENGTH];

  int *send_info = new int [numprocs];//number of particles need to be send corresponding to each process. If the value is zero, then no communication needed with that process
  int *recv_info = new int [numprocs];//number of particles will receive corresponding to each process. if the value is zero, no communication with that process
  for (i = 0; i < numprocs; i++)
  {
    send_info[i] = 0;
    recv_info[i] = 0;
  }

  MPI_Request * info_request = new MPI_Request [2 * numprocs];
  int tag0 = 43451, tag = 13249;

  // post receives
  for (i = 0; i < numprocs; i++)
    if ( my_comm[i] )
        j = MPI_Irecv ((recv_info + i), 1, MPI_INT, i, tag0,
                       MPI_COMM_WORLD, (info_request + i));

  // count number of BoundaryImages to send to other procs
  list < BndImage >::iterator i_img;
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->buckproc != myid)// if image of ghost particles does not belong to current process
      send_info[i_img->buckproc]++;

  
  for (i = 0; i < numprocs; i++)
    if ( my_comm[i] )
      j = MPI_Isend ((send_info + i), 1, MPI_INT, i, tag0,
                     MPI_COMM_WORLD, (info_request + numprocs + i));

  // allocate send buffer
  BndImage **sendbuf = new BndImage *[numprocs];//the data that will be send.
  for (i = 0; i < numprocs; i++)
    if (send_info[i] > 0)
      sendbuf[i] = new BndImage[send_info[i]];

  // wait for receives to finish
  MPI_Status status;
  for (i = 0; i < numprocs; i++)
    if ( my_comm[i] )
        j = MPI_Wait ((info_request + i), & status);

  // allocate receive buffer
  BndImage **recvbuf = new BndImage *[numprocs];
  for (i = 0; i < numprocs; i++)
    if (recv_info[i] > 0)
      recvbuf[i] = new BndImage[recv_info[i]];

  // post receive 
  MPI_Request *recv_req = new MPI_Request[numprocs];

  for (i = 0; i < numprocs; i++)
    if (recv_info[i] > 0)
      MPI_Irecv (recvbuf[i], recv_info[i], BND_IMAGE_TYPE,
                 i, tag, MPI_COMM_WORLD, &recv_req[i]);

  // pack and send images
  int *count = new int[numprocs];

  for (i = 0; i < numprocs; i++)
    count[i] = 0;

  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->buckproc != myid)// not communicate with yourself
    {
      j = i_img->buckproc;
      sendbuf[j][count[j]] = *i_img;
      count[j]++;
    }
  // post sends
  MPI_Request *send_req = new MPI_Request[numprocs];
  for (i = 0; i < numprocs; i++)
    if (send_info[i] > 0)
      MPI_Isend (sendbuf[i], send_info[i], BND_IMAGE_TYPE,
                 i, tag, MPI_COMM_WORLD, &send_req[i]);

  // take some time are remove old foreign images
  // so that there are no duplicates
  i_img = Image_table.begin ();
  while ( i_img !=  Image_table.end ())
    if ( i_img->partproc != myid )
      i_img = Image_table.erase (i_img);
    else
      i_img++;
        
  // add to local Image_table after receives
  for (i = 0; i < numprocs; i++)
    if (recv_info[i] != 0)
    {
      j = MPI_Wait (&recv_req[i], &status);
      for (j = 0; j < recv_info[i]; j++)
        Image_table.push_back (recvbuf[i][j]);
    }

  // clean up
  for (i = 0; i < numprocs; i++)
    if (recv_info[i] > 0)
      delete []recvbuf[i];
  delete [] recvbuf;
  delete [] recv_req;
  delete [] recv_info;

  // Wait for sends to finish
  for (i = 0; i < numprocs; i++)
    if ( my_comm[i] )
      j = MPI_Wait ((info_request + numprocs + i), & status);

  for (i = 0; i < numprocs; i++)
    if (send_info[i] != 0)
      j = MPI_Wait ((send_req + i), &status);

  for (i = 0; i < numprocs; i++)
    if (send_info[i] > 0)
      delete [] sendbuf[i];

  delete [] sendbuf;
  delete [] send_info;
  delete [] send_req;
  delete [] info_request;
  delete [] count;
}
