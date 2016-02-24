
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
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: BSFC_create_bins.C,v 1.2 2003/11/25 22:13:04 kdalbey Exp $ 
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <mpi.h>
#include <constant.h>

#include "exvar.h"
#include "repartition_BSFC.h"

/* 
 * global_actual_work_allocated is used to make sure that each processor
 * knows how much extra work every processor has allocated
 * extra work[proc] =
 *     global_actual_work_allocated[proc] - work_percent_array[proc]*total_work 
 */

void
BSFC_create_bins (int num_local_objects, sfc_vertex_t * sfc_vert_ptr,
                  int * amount_of_bits_used, int size_of_unsigned,
                  float * global_actual_work_allocated,
                  float *work_percent_array, float *total_weight_ptr,
                  int *balanced_flag,
                  unstructured_communication * verts_in_cut_info,
                  int * number_of_cuts, int bins_per_proc, int myid,
                  int numprocs)
{
  int i, j, number_of_bins, ierr = 0;

  float *binned_weight_array;
  float my_work_percent;
  int *bin_proc_array;
  float scanned_work_prev_allocated;    /*scanned_work_prev_allocated is the 
                                           amount of work allocated to higher
                                           ranked procs */
  int amount_of_bits;
  int current_proc;
  int local_balanced_flag;

  /*assume initially that each processor has the same amount of bins */
  number_of_bins = numprocs * bins_per_proc;
  i = 0;
  while (number_of_bins > BSFC_pow2 (i))
    i++;
  amount_of_bits = i;
  /* check to see that we have not used up all of the bits */
  if (amount_of_bits > 8 * size_of_unsigned * KEYLENGTH)
    amount_of_bits = 8 * size_of_unsigned * KEYLENGTH;
  number_of_bins = BSFC_pow2 (i);
  *amount_of_bits_used = amount_of_bits;

  float *tmp_float_array = new float[number_of_bins + 1];

  for (i = 0; i < number_of_bins; i++)
    tmp_float_array[i] = 0;

  for (i = 0; i < num_local_objects; i++)
  {
    sfc_vert_ptr[i].my_bin =
      BSFC_get_array_location (number_of_bins, amount_of_bits, 0,
                               (sfc_vert_ptr + i));
    tmp_float_array[sfc_vert_ptr[i].my_bin] += sfc_vert_ptr[i].lb_weight;
  }
  tmp_float_array[number_of_bins] = *total_weight_ptr;
  binned_weight_array = new float[number_of_bins + 1];

  i = MPI_Allreduce (tmp_float_array, binned_weight_array,
                     number_of_bins + 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  delete[]tmp_float_array;
  *total_weight_ptr = binned_weight_array[number_of_bins];
  /* global weight array has been created,
     now perform the scan operation on it */

  /* put in desired amount of work here, needs to
     be changed for varying workloads */

  // currently we assume that each processor should get the same amount of work
  my_work_percent = 1.0 / ((float) numprocs);
  work_percent_array[0] = 2.0;
  for (i = 1; i < numprocs; i++)
    work_percent_array[i] = ((float) (numprocs - i)) * my_work_percent;

  /* each processor needs to know which bins get partitioned
     into which processor, bin_proc_array lists max bin that
     a processor should get */
  for (i = 0; i < numprocs; i++)
    global_actual_work_allocated[i] = 0;

  bin_proc_array = new int[number_of_bins];

  current_proc = numprocs - 1;
  i = number_of_bins - 1;
  scanned_work_prev_allocated = 0;
  while (i >= 0 && current_proc > -1)
  {
    scanned_work_prev_allocated += binned_weight_array[i];
    bin_proc_array[i] = current_proc;
    if (scanned_work_prev_allocated >
        work_percent_array[current_proc] * (*total_weight_ptr))
    {
      global_actual_work_allocated[current_proc] = scanned_work_prev_allocated;
      bin_proc_array[i] = -current_proc;
      current_proc--;
      // the while statement is if there is more than 1 cut in a bin
      while (current_proc >= 0 &&
             scanned_work_prev_allocated >
             work_percent_array[current_proc] * *total_weight_ptr)
      {
        global_actual_work_allocated[current_proc] =
          scanned_work_prev_allocated;
        current_proc--;
      }
    }
    i--;
  }

  /* make sure that the last bin does not have a cut in it */
  if (bin_proc_array[0] > 0)
    bin_proc_array[0] = 0;
  global_actual_work_allocated[0] = *total_weight_ptr;

  /* specify which processor an object belongs to,
     we will know this because we know what bin an object 
     belongs to and we know what processor a bin belongs to */
  for (i = 0; i < num_local_objects; i++)
  {
    if (bin_proc_array[sfc_vert_ptr[i].my_bin] >= 0)
    {
      sfc_vert_ptr[i].cut_bin_flag = BSFC_NO_CUT;
      sfc_vert_ptr[i].destination_proc = bin_proc_array[sfc_vert_ptr[i].my_bin];
    }
    else
    {
      sfc_vert_ptr[i].cut_bin_flag = BSFC_CUT;
      sfc_vert_ptr[i].destination_proc =
        -bin_proc_array[sfc_vert_ptr[i].my_bin];
    }
  }

  /* check to see if any cut-bin has too many objects in it and refine it. 
     the problem with too many objects in a cut-bin is that the processor 
     that gets assigned to that bin will get swamped with communication
     and might not have enough memory to hold all of the information.  we
     detect the cut-bins with too many objects by how many cuts are in that 
     bin.  numprocs - 1 is the amount of cuts and if this is less than
     or equal to max_cuts_in_bin then there is no possibility that a
     coarse bin is overloaded */

  local_balanced_flag =
    BSFC_find_imbalance (work_percent_array, global_actual_work_allocated[myid],
                         *total_weight_ptr, myid, numprocs);

  ierr = MPI_Allreduce (&local_balanced_flag, balanced_flag, 1,
                        MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  delete[]bin_proc_array;
  delete[]binned_weight_array;

  /* if the current partitioning is acceptable, the algorithm is finished */
  if (*balanced_flag == BSFC_BALANCED)
  {
    //printf("first level is good enough on proc %d \n", myid);
    return;
  }

  //printf("#############    need to refine some bins....   ##############\n");

  /* move the sfc objects that belong to any bin that contains a cut
     to the proper processor */
  verts_in_cut_info->used_flag = 1;
  verts_in_cut_info->send_procs_ptr = new int[numprocs];
  verts_in_cut_info->recv_procs_ptr = new int[numprocs];

  for (i = 0; i < numprocs; i++)
    verts_in_cut_info->send_procs_ptr[i] = 0;
  for (i = 0; i < num_local_objects; i++)
    if (sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT)       // this vert is in a cut... 
      verts_in_cut_info->send_procs_ptr[sfc_vert_ptr[i].destination_proc] += 1;

  i = MPI_Alltoall (verts_in_cut_info->send_procs_ptr, 1, MPI_INT,
                    verts_in_cut_info->recv_procs_ptr,
                    1, MPI_INT, MPI_COMM_WORLD);

  //recalculate send_procs because it probably got changed
  for (i = 0; i < numprocs; i++)
    verts_in_cut_info->send_procs_ptr[i] = 0;
  for (i = 0; i < num_local_objects; i++)
    if (sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT)       // this vert is in a cut...
      verts_in_cut_info->send_procs_ptr[sfc_vert_ptr[i].destination_proc] += 1;

  MPI_Request *send_request = new MPI_Request[numprocs];
  MPI_Request *recv_request = new MPI_Request[numprocs];
  int send_count = 0, recv_count = 0;

  for (i = 0; i < numprocs; i++)
    send_count += verts_in_cut_info->send_procs_ptr[i];
  for (i = 0; i < numprocs; i++)
    recv_count += verts_in_cut_info->recv_procs_ptr[i];
  verts_in_cut_info->send_count = send_count;
  verts_in_cut_info->recv_count = recv_count;

  //temp storage for objects that get sent out...
  sfc_vertex_t * send_sfc_vert = new sfc_vertex_t [send_count];

  verts_in_cut_info->recv_sfc_vert = new sfc_vertex_t [recv_count];

  // fill up the send array...
  int *proc_counter = new int[numprocs];

  proc_counter[0] = 0;
  for (i = 1; i < numprocs; i++)
    proc_counter[i] =
      proc_counter[i - 1] + verts_in_cut_info->send_procs_ptr[i - 1];

  recv_count = 0;
  for (i = 0; i < myid; i++)
    recv_count += verts_in_cut_info->recv_procs_ptr[i];
  //printf("proc %d has recv_count of %d \n", myid, recv_count);
  for (i = 0; i < num_local_objects; i++)
    if (sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT)
    {
      if (sfc_vert_ptr[i].destination_proc != myid)
      {
        send_sfc_vert[proc_counter[sfc_vert_ptr[i].destination_proc]] =
          sfc_vert_ptr[i];
        proc_counter[sfc_vert_ptr[i].destination_proc] += 1;
      }
      else
      {                         // if i need to send to myself... 
        verts_in_cut_info->recv_sfc_vert[recv_count] = sfc_vert_ptr[i];
        recv_count++;
      }
    }
  //done filling up the send array
  int tag = 21503;

  proc_counter[0] = 0;
  for (i = 1; i < numprocs; i++)
    proc_counter[i] =
      proc_counter[i - 1] + verts_in_cut_info->send_procs_ptr[i - 1];

  recv_count = 0;
  for (i = 0; i < numprocs; i++)
  {
    if (i != myid)
    {
      // send out necessary info here...
      if (verts_in_cut_info->send_procs_ptr[i] != 0)
      {
        j = MPI_Isend ((send_sfc_vert + proc_counter[i]),
                       verts_in_cut_info->send_procs_ptr[i], LB_VERT_TYPE,
                       i, tag, MPI_COMM_WORLD, (send_request + i));
      }
      // receive necessary info here...
      if (verts_in_cut_info->recv_procs_ptr[i] != 0)
      {
        j = MPI_Irecv (&(verts_in_cut_info->recv_sfc_vert[recv_count]),
                       verts_in_cut_info->recv_procs_ptr[i], LB_VERT_TYPE,
                       i, tag, MPI_COMM_WORLD, (recv_request + i));
      }
    }
    recv_count += verts_in_cut_info->recv_procs_ptr[i];
  }
  delete[]proc_counter;

  // wait until the info is sent and received...
  for (i = 0; i < numprocs; i++)
    if (i != myid)
    {
      if (verts_in_cut_info->send_procs_ptr[i] != 0)
      {
        MPI_Status status;

        j = MPI_Wait ((send_request + i), &status);
      }
      if (verts_in_cut_info->recv_procs_ptr[i] != 0)
      {
        MPI_Status status;

        j = MPI_Wait ((recv_request + i), &status);
      }
    }

  delete[]send_request;
  delete[]recv_request;
  delete[]send_sfc_vert;
  return;

}                               /*  done sfc_create_bins routine */

/* create info before starting the multi-level refinement of the bins 
 */
void
BSFC_create_refinement_info (int * number_of_cuts, float * gl_work,
                             float total_weight, float * work_percent,
                             unstructured_communication verts_in_cut_info,
                             float ** work_prev_allocated, int myid,
                             int numprocs)
{
  float my_work, work2;
  int i = 0;

  //printf("the total weight is %e\n",total_weight);

  /* find out how many cuts are in this bin. */
  my_work = gl_work[myid];
  while (my_work > total_weight * work_percent[myid - i])
    i++;

  *number_of_cuts = i;
  if (verts_in_cut_info.recv_count == 0 || myid == 0)
  {
    *number_of_cuts = 0;
    return;
  }

  /* create link list for objects in the array.  link list
     is set up so that the objects in the array are
     traversed consecutively */
  for (i = 0; i < (verts_in_cut_info.recv_count - 1); i++)
    verts_in_cut_info.recv_sfc_vert[i].next_sfc_vert_index = i + 1;
  verts_in_cut_info.recv_sfc_vert[verts_in_cut_info.recv_count -
                                  1].next_sfc_vert_index = -1;

  /* update work previously allocated to include work in all bins
     with higher keys than this bin */
  work2 = 0;
  for (i = 0; i < verts_in_cut_info.recv_count; i++)
    work2 += verts_in_cut_info.recv_sfc_vert[i].lb_weight;

  work2 = gl_work[myid] - work2;

  // note that work_prev_allocated data will only be correct for the
  // actual processors that are getting examined
  *work_prev_allocated = new float [numprocs];
  int ncuts = *(number_of_cuts);

  for (i = myid - ncuts; i <= myid; i++)
    (*work_prev_allocated)[i] = work2;

  return;
}

/* 
   routine calculates the new bin number of an object using its
   sfc_key.  prev_used_bits is how many bits have already been
   used and number_of_bits is how many bits to use to calculate
   the key.  the output is the bin number and will be a value
   between 0 and (2^number_of_bits - 1)
*/

int
BSFC_get_array_location (int number_of_bins, int number_of_bits,
                         int prev_used_bits, sfc_vertex_t * sfc_vert_ptr)
{
  unsigned ilocation = 0;
  unsigned ilocation2 = 0;
  //int pub = prev_used_bits;
  int size_of_unsigned = sizeof (unsigned);
  int counter = 0;

  if (prev_used_bits == 0)
    ilocation =
      (sfc_vert_ptr->sfc_key[0]) >> (size_of_unsigned * 8 - number_of_bits);
  else
  {
    /* in case prev_used_bits is larger than an unsigned integer */
    while ((counter + 1) * size_of_unsigned * 8 < prev_used_bits)
      counter++;
    prev_used_bits = prev_used_bits - counter * size_of_unsigned * 8;

    ilocation2 = (sfc_vert_ptr->sfc_key[counter]) << prev_used_bits;
    ilocation = ilocation2 >> (size_of_unsigned * 8 - number_of_bits);
    /* if some of the bits that we want are in the next array value
       this might not be correct!!! */
    if (prev_used_bits + number_of_bits > size_of_unsigned * 8)
      ilocation += ((sfc_vert_ptr->sfc_key[counter + 1]) >>
                    (2 * size_of_unsigned * 8 - prev_used_bits -
                     number_of_bits));
  }
  if (ilocation >= (unsigned) number_of_bins)
  {
    //int myid;

    ilocation = number_of_bins - 1;
  }

  return (ilocation);
}

/* routine finds the imbalance for processor which_proc given the
   cumulative work assigned to processors which_proc through the
   last processor (numprocs-1).  routine returns value which 
   indicates whether this imbalance is beyond a specified 
   tolerance.   */
int
BSFC_find_imbalance (float *work_percent_array,
                     float cumulative_work,
                     float total_work, int which_proc, int numprocs)

/* NOTE:  which_proc is not required to be this procs rank */
{
  int balanced_flag;
  float my_extra_work;
  float my_ideal_work;

  my_extra_work = cumulative_work - work_percent_array[which_proc] * total_work;

  if (which_proc != numprocs - 1)
    my_ideal_work = (work_percent_array[which_proc] -
                     work_percent_array[which_proc + 1]) * total_work;
  else
    my_ideal_work = work_percent_array[which_proc] * total_work;

  /* if processor which_proc is not supposed to have any work,
     it is imbalanced if it has any amount of work greater
     than 0 */
  if (my_ideal_work == 0)
  {
    if (my_extra_work != 0)
      balanced_flag = BSFC_NOT_BALANCED;
    else
      balanced_flag = BSFC_BALANCED;
  }
  else
  {
    if (1 + my_extra_work / my_ideal_work > LOAD_BALANCE_TOLERANCE)
      balanced_flag = BSFC_NOT_BALANCED;
    else
      balanced_flag = BSFC_BALANCED;
  }

  return (balanced_flag);
}
