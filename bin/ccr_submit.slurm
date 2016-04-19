#!/bin/sh  
#SBATCH --partition=general-compute
#SBATCH --time=71:30:00 
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=12
#SBATCH --mem=24000 
##SBATCH --account=mae610s15
##Memory per node specification is in MB. It is optional. # The default limit is 3GB per core. 
#SBATCH --mail-user=zhixuanc@buffalo.edu 
#SBATCH --mail-type=ALL 
#SBATCH --constraint=IB&CPU-E5645
##SBATCH --exclusive
#SBATCH --requeue 
#SBATCH --job-name="sml100" 
#SBATCH --output=sml100.out
#SBATCH --error=sml100.err
#SBATCH --qos=supporters
#
echo "SLURM_JOBID="$SLURM_JOBID 
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST 
echo "SLURM_NNODES"=$SLURM_NNODES 
echo "SLURMTMPDIR="$SLURMTMPDIR 
cd $SLURM_SUBMIT_DIR echo "working directory = "$SLURM_SUBMIT_DIR  
#
module load intel-mpi/5.0.2
module load mkl/11.2
module load hdf5/1.8.15p1
module list

#whether the following is still necessary or not need double check.
export HDF5_DISABLE_VERSION_CHECK=1

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

#for MPI debug, the larger the value, the more details will given
export I_MPI_DEBUG=1

echo "srun with $SLURM_NNODES nodes, each node has 8 processors"
srun ./particler 1>out 2>err
echo "All Done!"

