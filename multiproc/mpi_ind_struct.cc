/*
 * mpi_ind_struct.cc
 *
 *  Created on: May 26, 2015
 *      Author: zhixuanc
 */

#include <mpi.h>
#include <constant.h>

void create_indmpi_struct()
{
 //Define a new MPI data type corresponding to IndMap
    int one=1;
    MPI_Aint zero=0;
    int tkeylength = TKEYLENGTH;
    MPI_Datatype MPI_TKEY, MPI_IndMap;
    MPI_Datatype unsignd = MPI_UNSIGNED;
    MPI_Type_struct( one, &tkeylength, &zero, &unsignd, &MPI_TKEY );
    MPI_Datatype mpitype[2] = {MPI_TKEY, MPI_INT};
    int block_sz[2]={1,3};
    MPI_Aint initial_bt[2]={0,12};//initial byte displacement
    MPI_Type_struct( one, block_sz, initial_bt, mpitype, &MPI_IndMap);
    MPI_Type_commit(&MPI_IndMap);
}
