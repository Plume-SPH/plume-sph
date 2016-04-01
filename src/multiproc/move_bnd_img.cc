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
 * Description: update ghost particles, whose reflections are on
 *              other procs
 *
 *******************************************************************
 * $Id: $
 */
#include <list>
#include <cassert>
using namespace std;

#include <mpi.h>

#include <hashtab.h>
#include <thashtab.h>
#include <particle.h>
#include <bucket.h>
#include <bnd_image.h>
#include <parameters.h> //This will be needed while updating secondary variables
#include <exvar.h>

int
move_bnd_images (int myid, int nump, THashTable * P_table, HashTable * BG_mesh,
                 list < BndImage > & Image_table)
{

  if (nump < 2)
    return 0;

  int i, j, k;
  double uvec[NO_OF_EQNS];

  int img_tag = 3152;
  int *send_info = new int[nump];
  int *recv_info = new int[nump];
  int *recv_info2 = new int[nump];

  for (i = 0; i < nump; i++)
  {
    send_info[i] = 0;
    recv_info[i] = 0;
    recv_info2[i] = 0;
  }

  /* if the bucket that contains reflected image belongs to 
   * foreing process, current proc needs to receive data from
   * that proc
   */
  
  list < BndImage >::iterator i_img;
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->buckproc != myid)
      recv_info2[i_img->buckproc]++;

  /* calulate the amount of data that needs to be sent */
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->partproc != myid)
      send_info[i_img->partproc]++;

  // communicate send_info to update recv_info
  MPI_Request *req = new MPI_Request[2 * nump];
  int tag0 = 5132;

  /* start receiving data from foreign procs */
  for (i = 0; i < nump; i++)
    if (recv_info2[i] > 0)
      MPI_Irecv ((recv_info + i), 1, MPI_INT, i, tag0, MPI_COMM_WORLD,
                 (req + i));

  /* start sending data to foreign procs */
  for (i = 0; i < nump; i++)
    if (send_info[i] > 0)
      MPI_Isend ((send_info + i), 1, MPI_INT, i, tag0, MPI_COMM_WORLD,
                 (req + nump + i));

  MPI_Status status;

  // wait for sends to finish
  for (i = 0; i < nump; i++)
    if (recv_info2[i] > 0)
      k = MPI_Wait ((req + i), &status);

  // wait for recvs to finish
  for (i = 0; i < nump; i++)
    if (send_info[i] > 0)
      k = MPI_Wait ((req + nump + i), &status);

  delete [] req;

  // check data integrity
  for (i = 0;  i < nump; i++)
    if ( recv_info[i] != recv_info2[i] )
    {
      fprintf(stderr,"data sizes don't match %s : %d\n", __FILE__, __LINE__);
      exit (1);
    }

  //  allocate recv_buffer
  BndImage **recv_buf = new BndImage *[nump];

  for (j = 0; j < nump; j++)
    if (recv_info[j] > 0)
      recv_buf[j] = new BndImage[recv_info[j]];

  // post reveive calls
  MPI_Request *img_recv_req = new MPI_Request[nump];
  for (j = 0; j < nump; j++)
    if (recv_info[j] > 0)
      MPI_Irecv (recv_buf[j], recv_info[j], BND_IMAGE_TYPE, j,
                 img_tag, MPI_COMM_WORLD, (img_recv_req + j));

  // allocate memory
  int *img_counter = new int[nump];
  BndImage **send_buf = new BndImage *[nump];

  for (j = 0; j < nump; j++)
  {
    if (send_info[j] > 0)
      send_buf[j] = new BndImage[send_info[j]];
    img_counter[j] = 0;
  }

  // pack Images
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->partproc != myid)
    {
      j = i_img->partproc;
      send_buf[j][img_counter[j]] = *i_img;
      img_counter[j]++;
    }

  // send buffer
  MPI_Request *img_send_req = new MPI_Request[nump];

  for (j = 0; j < nump; j++)
    if (send_info[j] > 0)
      MPI_Isend (send_buf[j], send_info[j], BND_IMAGE_TYPE, j,
                 img_tag, MPI_COMM_WORLD, (img_send_req + j));

  // Wait and update ghost particles
  for (j = 0; j < nump; j++)
    if (recv_info[j] > 0)
    {
      MPI_Wait ((img_recv_req + j), &status);
      for (k = 0; k < recv_info[j]; k++)
      {
        Particle *pghost =
          (Particle *) P_table->lookup (recv_buf[j][k].ghost_key);
        assert (pghost);
        for (int i2 = 0; i2 < NO_OF_EQNS; i2++)
          uvec[i2] = recv_buf[j][k].state_vars[i2];
        pghost->put_state_vars (uvec);
        pghost->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P); //This is added later, there was a bug in old code as the secondary variable is not updated after imposing of boundary condition.
        pghost->put_update_delayed(false); //This also added later, but I do not think this mistake is a big deal in the old code.
      }
    }

  // Wait for sends to finish
  for (j = 0; j < nump; j++)
    if (send_info[j] > 0)
      MPI_Wait ((img_send_req + j), &status);

  delete [] img_send_req;
  delete [] img_recv_req;
  delete [] img_counter;
  for (j = 0; j < nump; j++)
  {
    if (send_info[j] > 0)
      delete [] send_buf[j];
    if (recv_info[j] > 0)
      delete [] recv_buf[j];
  }
  delete [] send_buf;
  delete [] recv_buf;
  delete [] send_info;
  delete [] recv_info;
  delete [] recv_info2;

  return 0;
}


