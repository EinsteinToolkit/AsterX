#!/bin/bash

#SBATCH -p rtx-dev         # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 40
#SBATCH -n 4               # Total # of mpi tasks 280
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH -A PHY20010        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --gpu-bind=map_gpu:0,1,2,3
##SBATCH --mail-type=all    # Send email at begin and end of job
##SBATCH --mail-user=@rit.edu

# Any other commands must follow all #SBATCH directives...
#export OMP_NUM_THREADS=4
ml

source $SPACK_DIR/share/spack/setup-env.sh
spack load gcc@11.2.0
spack load cuda@11.5.2

ibrun ./cactus_CarpetX-cuda qc0.par
