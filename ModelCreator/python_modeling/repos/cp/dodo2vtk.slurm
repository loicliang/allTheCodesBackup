#! /bin/bash

#SBATCH --nodes 1
#SBATCH -c 1
#SBATCH -p work
#SBATCH --mem=10G
source /opt/share/intel/bin/ifortvars.sh intel64
module load mpi/mpich-intel
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
srun --mpi=pmi2 /home/xliang/INSTALL6/bin/dodo2vtk tt

