#!/bin/sh  
#SBATCH --partition=general-compute
#SBATCH --time=71:30:00 
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=1
##SBATCH --mem=12000 
##SBATCH --account=mae610s15
##Memory per node specification is in MB. It is optional. # The default limit is 3GB per core. 
#SBATCH --mail-user=zhixuanc@buffalo.edu 
#SBATCH --mail-type=ALL 
##SBATCH --constraint=CPU-E5645
#SBATCH --requeue 
#SBATCH --job-name="sml130" 
#SBATCH --output=sml130.out
#SBATCH --error=sml130.err
#
echo "SLURM_JOBID="$SLURM_JOBID 
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST 
echo "SLURM_NNODES"=$SLURM_NNODES 
echo "SLURMTMPDIR="$SLURMTMPDIR 
cd $SLURM_SUBMIT_DIR echo "working directory = "$SLURM_SUBMIT_DIR  
#
module load intel-mpi/5.0.2
module load mkl/11.2
module load hdf
module list

#echo "Run with 32 threads."
#mpirun -n 32 ./my_4000 > only_mpi_slurm.txt
#

#echo "Run with 64 threads."
#mpirun -n 64 ./my_4000 >> only_mpi_slurm.txt
#

#echo "N,	t,	flops,	   w"
#./try 10 >try.txt
#for ((np=30; np<=10000000000; np *= 2)) do
#    ./try $np >>try.txt
#done

echo "run with 32 processes"
mpirun -n 32 ./particler >out
echo "All Done!"

